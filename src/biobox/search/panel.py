"""One-pass multi-HMM panel search, per-ORF hit table, and co-occurring system caller (card W-39).

Searches a panel of K target HMMs (+ optional outgroup/composition decoy HMMs) against a release's protein
pack (or a `--stages` subset) in a single streaming pass per `(shard_slice, model_batch)` sub-step, keeping
every sub-step under 200 CPU-minutes and < 2 GB RSS (AGENTS §4.1, §7.2).

Outputs are written outside the data root (`--out DIR`) with explicit PyArrow schemas defined here (never
modifying `src/biobox/registry/schemas.py` or schema locks):
  models.parquet / models.tsv
  protein_hits.parquet
  orf_hits.parquet / orf_hits.tsv
  loci_summary.parquet / loci_summary.tsv
  loci/locus_*.tsv
  controls.tsv, coverage.json, unannotated_proteins.txt, RUN, README.txt
plus `<runs>/run=<run_id>/run.parquet` (`recipe = "panel.v1"`) and `panel_hits.parquet`.
"""

from __future__ import annotations

import concurrent.futures
import contextlib
import hashlib
import heapq
import io
import json
import math
import os
import platform
import re
import resource
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyhmmer
from blake3 import blake3

from biobox import catalog, ids, provenance, runstep
from biobox.l1 import domains
from biobox.l1.domains import DOM_Z, Z_PROT, iter_hmms
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.registry.schemas import TABLES, UID, table
from biobox.search import arrays, engine, export, fastz, nbr, pack, runner

PANEL_RECIPE = "panel.v1"
# <= 10M residues per streaming block keeps two prefetched blocks well under 2 GB RSS
PANEL_BLOCK_RESIDUES = 10_000_000
# Conservative CPU-min per (residue * model_len) derived from GUIDE.md §6 (9.5 CPU-min / 3.4e9 aa / ~200 aa)
EST_CPU_MIN_PER_RES_M = 1.5e-11
DEFAULT_MAX_SUBSTEP_CPU_MIN = 200.0
ROW_GROUP_ROWS = 122_880
SMORF_MAX_AA = 100
PANEL_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
CAS_DOMAIN_RE = re.compile(
    r"(?:^|[^A-Za-z0-9])(cas[0-9]|cas_|crispr_cas|csy[0-9]|cse[0-9]|csm[0-9]|cmr[0-9]|csn[0-9]|csx[0-9]|csd[0-9]|csc[0-9]|cpf1)",
    re.IGNORECASE,
)


class PanelError(ValueError):
    """Refusal before running a panel search (exit 2)."""


# --- Explicit PyArrow schemas (never modifying schemas.py) ----------------------------------------------

MODELS_SCHEMA = pa.schema([
    pa.field("model_key", pa.string(), nullable=False),
    pa.field("namespace", pa.string(), nullable=False),
    pa.field("name", pa.string(), nullable=False),
    pa.field("acc", pa.string(), nullable=False),
    pa.field("model_len", pa.uint32(), nullable=False),
    pa.field("hmm_sha256", pa.string(), nullable=False),
    pa.field("is_decoy", pa.bool_(), nullable=False),
    pa.field("cutoff_mode", pa.string(), nullable=False),
    pa.field("seq_bits", pa.float32(), nullable=False),
    pa.field("dom_bits", pa.float32(), nullable=False),
    pa.field("i_evalue_max", pa.float64(), nullable=True),
    pa.field("hmm_cov_min", pa.float32(), nullable=True),
    pa.field("min_dom_frac", pa.float32(), nullable=True),
    pa.field("n_hit_proteins", pa.uint64(), nullable=False),
    pa.field("n_best_proteins", pa.uint64(), nullable=False),
    pa.field("n_hit_orfs", pa.uint64(), nullable=False),
    pa.field("pos_controls_passed", pa.uint32(), nullable=False),
    pa.field("pos_controls_total", pa.uint32(), nullable=False),
    pa.field("neg_control_hits", pa.uint32(), nullable=False),
    pa.field("neg_shuffled_hits", pa.int64(), nullable=True),
    pa.field("neg_shuffled_pass", pa.bool_(), nullable=True),
])

PROTEIN_HITS_SCHEMA = pa.schema([
    pa.field("protein_uid", UID, nullable=False),
    pa.field("model_key", pa.string(), nullable=False),
    pa.field("namespace", pa.string(), nullable=False),
    pa.field("model_name", pa.string(), nullable=False),
    pa.field("model_acc", pa.string(), nullable=False),
    pa.field("model_len", pa.uint32(), nullable=False),
    pa.field("is_decoy", pa.bool_(), nullable=False),
    pa.field("seq_bitscore", pa.float32(), nullable=False),
    pa.field("seq_evalue_zref", pa.float64(), nullable=False),
    pa.field("bitscore", pa.float32(), nullable=False),
    pa.field("i_evalue_zref", pa.float64(), nullable=False),
    pa.field("env_start", pa.uint32(), nullable=False),
    pa.field("env_end", pa.uint32(), nullable=False),
    pa.field("ali_start", pa.uint32(), nullable=False),
    pa.field("ali_end", pa.uint32(), nullable=False),
    pa.field("hmm_from", pa.uint32(), nullable=False),
    pa.field("hmm_to", pa.uint32(), nullable=False),
    pa.field("hmm_cov", pa.float32(), nullable=False),
    pa.field("n_model_domains", pa.uint16(), nullable=False),
    pa.field("dom_seq_frac", pa.float32(), nullable=False),
    pa.field("best_on_overlap", pa.bool_(), nullable=False),
    pa.field("best_in_namespace", pa.bool_(), nullable=False),
    pa.field("vetoed_by_decoy", pa.bool_(), nullable=False),
])

ORF_HITS_SCHEMA = pa.schema([
    pa.field("contig_uid", UID, nullable=False),
    pa.field("start", pa.int64(), nullable=False),
    pa.field("end", pa.int64(), nullable=False),
    pa.field("strand", pa.int8(), nullable=False),
    pa.field("feature_uid", UID, nullable=False),
    pa.field("ordinal", pa.uint32(), nullable=False),
    pa.field("tape_version", pa.string(), nullable=False),
    pa.field("contig_accession", pa.string(), nullable=True),
    pa.field("contig_len", pa.int64(), nullable=False),
    pa.field("is_circular", pa.bool_(), nullable=False),
    pa.field("genetic_code", pa.uint8(), nullable=True),
    pa.field("stage", pa.string(), nullable=False),
    pa.field("genome_kind", pa.string(), nullable=True),
    pa.field("taxonomy", pa.string(), nullable=True),
    pa.field("protein_uid", UID, nullable=False),
    pa.field("nr90_uid", UID, nullable=True),
    pa.field("nr50_uid", UID, nullable=True),
    pa.field("domain_state", pa.string(), nullable=False),
    pa.field("pfam_architecture", pa.string(), nullable=True),
    pa.field("model_key", pa.string(), nullable=False),
    pa.field("namespace", pa.string(), nullable=False),
    pa.field("model_name", pa.string(), nullable=False),
    pa.field("model_acc", pa.string(), nullable=False),
    pa.field("model_len", pa.uint32(), nullable=False),
    pa.field("is_decoy", pa.bool_(), nullable=False),
    pa.field("seq_bitscore", pa.float32(), nullable=False),
    pa.field("seq_evalue_zref", pa.float64(), nullable=False),
    pa.field("bitscore", pa.float32(), nullable=False),
    pa.field("i_evalue_zref", pa.float64(), nullable=False),
    pa.field("env_start", pa.uint32(), nullable=False),
    pa.field("env_end", pa.uint32(), nullable=False),
    pa.field("ali_start", pa.uint32(), nullable=False),
    pa.field("ali_end", pa.uint32(), nullable=False),
    pa.field("hmm_from", pa.uint32(), nullable=False),
    pa.field("hmm_to", pa.uint32(), nullable=False),
    pa.field("hmm_cov", pa.float32(), nullable=False),
    pa.field("n_model_domains", pa.uint16(), nullable=False),
    pa.field("dom_seq_frac", pa.float32(), nullable=False),
    pa.field("best_on_overlap", pa.bool_(), nullable=False),
    pa.field("best_in_namespace", pa.bool_(), nullable=False),
    pa.field("vetoed_by_decoy", pa.bool_(), nullable=False),
])

LOCI_SUMMARY_SCHEMA = pa.schema([
    pa.field("contig_uid", UID, nullable=False),
    pa.field("anchor_feature_uid", UID, nullable=False),
    pa.field("contig_accession", pa.string(), nullable=True),
    pa.field("stage", pa.string(), nullable=False),
    pa.field("genome_kind", pa.string(), nullable=True),
    pa.field("is_circular", pa.bool_(), nullable=False),
    pa.field("contig_len", pa.int64(), nullable=False),
    pa.field("cluster_start", pa.int64(), nullable=False),
    pa.field("cluster_end", pa.int64(), nullable=False),
    pa.field("span_bp", pa.int64(), nullable=False),
    pa.field("censored", pa.bool_(), nullable=False),
    pa.field("models_present", pa.string(), nullable=False),
    pa.field("models_ordered", pa.string(), nullable=False),
    pa.field("n_models", pa.uint32(), nullable=False),
    pa.field("n_panel_orfs", pa.uint32(), nullable=False),
    pa.field("n_unmatched_internal_orfs", pa.uint32(), nullable=False),
    pa.field("fusion_orfs", pa.uint32(), nullable=False),
    pa.field("n_crispr_arrays", pa.uint32(), nullable=False),
    pa.field("max_adjacent_igr_bp", pa.int64(), nullable=False),
    pa.field("adjacent_igr_uids", pa.string(), nullable=False),
    pa.field("bridged_noncoding_bp", pa.int64(), nullable=False),
    pa.field("cas_less_crispr", pa.bool_(), nullable=False),
    pa.field("split_model_pair", pa.bool_(), nullable=False),
])


@dataclass
class ModelSpec:
    model_key: str
    namespace: str
    name: str
    acc: str
    model_len: int
    hmm_sha256: str
    is_decoy: bool
    cutoff_mode: str
    seq_bits: float
    dom_bits: float
    i_evalue_max: float | None
    hmm_cov_min: float | None
    min_dom_frac: float | None
    hmm: object  # pyhmmer.plan7.HMM

    def signature(self) -> tuple:
        """Exact carry-forward identity tuple (MR-12)."""
        return (
            self.model_key,
            self.hmm_sha256,
            self.is_decoy,
            self.cutoff_mode,
            round(float(self.seq_bits), 6),
            round(float(self.dom_bits), 6),
            None if self.i_evalue_max is None else float(self.i_evalue_max),
            None if self.hmm_cov_min is None else round(float(self.hmm_cov_min), 6),
            None if self.min_dom_frac is None else round(float(self.min_dom_frac), 6),
        )


# --- Panel loading & cutoffs ----------------------------------------------------------------------------


def _split_ns_file(spec: str) -> tuple[str, Path]:
    """Parse `[NAMESPACE:]FILE`. Only treat leading token before `:` as namespace when `spec` is not an existing path."""
    p = Path(spec)
    if p.exists():
        return "", p
    if ":" in spec:
        ns, rest = spec.split(":", 1)
        if ns and "/" not in ns and rest:
            return ns, Path(rest)
    return "", p


def _split_lib_spec(spec: str) -> tuple[str, str, list[str] | None]:
    """Parse `[NAMESPACE@]LIB[:ACC1,ACC2,...]`."""
    ns = ""
    rest = spec
    if "@" in rest:
        ns, rest = rest.split("@", 1)
    if ":" in rest:
        lib, acc_str = rest.split(":", 1)
        accs = [x.strip() for x in acc_str.split(",") if x.strip()]
        if not accs:
            raise PanelError(f"library spec {spec!r} has empty accession list after ':'")
        return ns, lib, accs
    return ns, rest, None


_VERSION_TAIL = re.compile(r"-(\d+(?:\.\d+)*)$")
_DOTTED = re.compile(r"\d+(?:\.\d+)*$")


def _version_key(path: Path, stop: Path | None = None) -> tuple:
    """Sort key of a library file: the trailing `-<dotted version>` of its name as integers (38.10 is newer than
    38.2), else the nearest ancestor directory below `stop` whose whole name is a dotted version (the refdb layout
    `pfam/38.2/Pfam-A.hmm.gz`). A file with neither sorts before every versioned one, so it never wins over one."""
    stem = path.name
    for ext in (".gz", ".tgz", ".tar", ".hmm"):
        stem = stem.removesuffix(ext)
    m = _VERSION_TAIL.search(stem)
    if m is not None:
        return (1, tuple(int(x) for x in m.group(1).split(".")))
    if stop is not None:
        for anc in path.parents:
            if anc == stop or stop not in anc.parents:
                break
            if _DOTTED.fullmatch(anc.name):
                return (1, tuple(int(x) for x in anc.name.split(".")))
    return (0, ())


def _newest(paths: list[Path], stop: Path | None = None) -> Path:
    """The path with the highest version key; among equal keys the last in sorted path order (as before R-06)."""
    return sorted(sorted(paths), key=lambda p: _version_key(p, stop))[-1]


def _resolve_library_path(lib: str, pack_dir: Path | None, root: Path | None, refdb_dirs: Iterable[Path]) -> Path:
    """Resolve a library name (e.g. 'pfam', 'pfam-38.2', 'ncbifam') to an HMM archive path."""
    candidates: list[Path] = []
    if pack_dir is not None:
        mdir = Path(pack_dir) / "models"
        if mdir.is_dir():
            for p in sorted(mdir.iterdir()):
                stem = p.name.removesuffix(".hmm.gz").removesuffix(".hmm")
                if stem == lib or stem.startswith(f"{lib}-"):
                    candidates.append(p)
    if candidates:
        return _newest(candidates)
    search_bases = list(refdb_dirs)
    if root is not None:
        search_bases.append(Path(root) / "l0" / "refdb")
    for base in search_bases:
        d = Path(base) / lib
        if d.is_dir():
            for ext in ("*.hmm", "*.hmm.gz", "*.tgz", "*.tar.gz"):
                found = sorted(d.rglob(ext))
                if found:
                    return _newest(found, d)
    raise PanelError(f"library {lib!r} not found in pack models or refdb directories")


def parse_cutoffs_tsv(path: Path) -> dict[str, dict]:
    """Parse `model_key<TAB>seq_bits[<TAB>dom_bits[<TAB>i_evalue_max[<TAB>hmm_cov_min[<TAB>min_dom_frac]]]]]`."""
    out: dict[str, dict] = {}
    for lineno, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if lineno == 1 and cols[0] == "model_key":
            continue
        if len(cols) < 2:
            raise PanelError(f"{path}:{lineno}: expected at least model_key<TAB>seq_bits")
        key = cols[0].strip()
        if not key:
            raise PanelError(f"{path}:{lineno}: empty model_key")

        def _f(idx: int) -> float | None:
            if idx >= len(cols) or cols[idx].strip() in ("", "-", "none", "None"):
                return None
            return float(cols[idx].strip())

        seq_b = _f(1)
        if seq_b is None:
            raise PanelError(f"{path}:{lineno}: seq_bits is required for {key}")
        out[key] = {
            "seq_bits": seq_b,
            "dom_bits": _f(2),
            "i_evalue_max": _f(3),
            "hmm_cov_min": _f(4),
            "min_dom_frac": _f(5),
        }
    return out


def _hmm_cutoff(hmm, mode: str) -> tuple[float, float] | None:
    if mode == "gathering" and hmm.cutoffs.gathering_available():
        g = hmm.cutoffs.gathering
        return float(g[0]), float(g[1])
    if mode == "trusted" and hmm.cutoffs.trusted_available():
        t = hmm.cutoffs.trusted
        return float(t[0]), float(t[1])
    if mode == "noise" and hmm.cutoffs.noise_available():
        n = hmm.cutoffs.noise
        return float(n[0]), float(n[1])
    return None


def load_panel_specs(
    hmm_specs: list[str],
    library_specs: list[str],
    decoy_hmm_specs: list[str],
    decoy_library_specs: list[str],
    *,
    cutoff: str | None = None,
    fallback_bits: float | None = None,
    bits: float | None = None,
    dom_bits: float | None = None,
    cutoffs_tsv: Path | None = None,
    pack_dir: Path | None = None,
    root: Path | None = None,
    refdb_dirs: Iterable[Path] = (),
) -> list[ModelSpec]:
    """Load target and decoy HMMs, assign deterministic `model_key`s, and resolve heterogeneous cutoffs."""
    if not hmm_specs and not library_specs:
        raise PanelError("panel needs at least one target --hmm or --library")
    tsv_cuts = parse_cutoffs_tsv(cutoffs_tsv) if cutoffs_tsv is not None else {}

    raw_items: list[tuple[str, bool, object]] = []

    def _add_file(spec: str, is_decoy: bool) -> None:
        ns, path = _split_ns_file(spec)
        if not path.is_file():
            raise PanelError(f"HMM file not found: {path}")
        loaded = list(iter_hmms(path))
        if not loaded:
            raise PanelError(f"HMM file {path} contains no profiles")
        for h in loaded:
            if not h.alphabet.is_amino():
                raise PanelError(f"model {engine.text(h.name)} in {path} is not an amino-acid profile")
            raw_items.append((ns, is_decoy, h))

    def _add_lib(spec: str, is_decoy: bool) -> None:
        ns, lib, want_accs = _split_lib_spec(spec)
        path = _resolve_library_path(lib, pack_dir, root, refdb_dirs)
        want_set = set(want_accs) if want_accs is not None else None
        matched_want: set[str] = set()
        for h in iter_hmms(path):
            if not h.alphabet.is_amino():
                continue
            if want_set is not None:
                acc = engine.text(h.accession)
                name = engine.text(h.name)
                keys = {k for k in (acc, acc.split(".")[0], name) if k}
                hit = keys & want_set
                if not hit:
                    continue
                matched_want |= hit
            raw_items.append((ns or lib, is_decoy, h))
        if want_set is not None:
            missing = sorted(want_set - matched_want)
            if missing:
                raise PanelError(f"library {lib!r} ({path}) does not contain model(s): {', '.join(missing)}")

    for s in hmm_specs:
        _add_file(s, False)
    for s in library_specs:
        _add_lib(s, False)
    for s in decoy_hmm_specs:
        _add_file(s, True)
    for s in decoy_library_specs:
        _add_lib(s, True)

    # Check accession uniqueness per namespace so model_key is [NS:]ACC iff ACC is present & unique
    acc_counts: dict[tuple[str, str], int] = defaultdict(int)
    for ns, _, h in raw_items:
        acc = engine.text(h.accession)
        if acc:
            acc_counts[(ns, acc)] += 1

    specs: list[ModelSpec] = []
    seen_keys: set[str] = set()
    for ns, is_decoy, h in raw_items:
        name = engine.text(h.name)
        acc = engine.text(h.accession)
        base_id = acc if (acc and acc_counts[(ns, acc)] == 1) else name
        model_key = f"{ns}:{base_id}" if ns else base_id
        if model_key in seen_keys:
            raise PanelError(f"duplicate model_key {model_key!r} in panel; use distinct NAMESPACE: prefixes")
        seen_keys.add(model_key)

        hmm_sha = hashlib.sha256(engine._hmm_text(h)).hexdigest()

        # Resolve cutoff: cutoffs-tsv > --cutoff (GA/TC/NC) > --fallback-bits > --bits
        tc = tsv_cuts.get(model_key)
        if tc is None and acc:
            tc = tsv_cuts.get(acc) or tsv_cuts.get(acc.split(".")[0])
        if tc is None:
            tc = tsv_cuts.get(name)

        if tc is not None:
            cmode = "tsv"
            sb = float(tc["seq_bits"])
            db = float(tc["dom_bits"]) if tc["dom_bits"] is not None else 0.0
            ie_max = tc["i_evalue_max"]
            cov_min = tc["hmm_cov_min"]
            mdf = tc["min_dom_frac"]
        elif cutoff is not None:
            got = _hmm_cutoff(h, cutoff)
            if got is not None:
                cmode = cutoff
                sb, db = got
            elif fallback_bits is not None:
                cmode = f"fallback_bits"
                sb = float(fallback_bits)
                db = float(dom_bits) if dom_bits is not None else 0.0
            elif bits is not None:
                cmode = "bits"
                sb = float(bits)
                db = float(dom_bits) if dom_bits is not None else 0.0
            else:
                raise PanelError(
                    f"model {model_key!r} has no {cutoff} cutoff header; supply --fallback-bits or --cutoffs-tsv"
                )
            ie_max, cov_min, mdf = None, None, None
        elif bits is not None:
            cmode = "bits"
            sb = float(bits)
            db = float(dom_bits) if dom_bits is not None else 0.0
            ie_max, cov_min, mdf = None, None, None
        else:
            raise PanelError(f"no cutoff specified for model {model_key!r}; give --cutoff, --bits, or --cutoffs-tsv")

        specs.append(ModelSpec(
            model_key=model_key, namespace=ns, name=name, acc=acc, model_len=int(h.M),
            hmm_sha256=hmm_sha, is_decoy=is_decoy, cutoff_mode=cmode, seq_bits=sb, dom_bits=db,
            i_evalue_max=ie_max, hmm_cov_min=cov_min, min_dom_frac=mdf, hmm=h,
        ))

    if not any(not m.is_decoy for m in specs):
        raise PanelError("panel has no non-decoy target models")
    return specs


# --- Filtering & overlap / decoy competition ------------------------------------------------------------


def passes_min_dom_frac(max_dom_bits: float, seq_bitscore: float, min_dom_frac: float | None) -> bool:
    """Reject multi-repeat bitscore stacking where max_dom_bits / seq_bitscore < min_dom_frac."""
    if min_dom_frac is None or min_dom_frac <= 0:
        return True
    if seq_bitscore <= 0:
        return False
    return (max_dom_bits / seq_bitscore) >= min_dom_frac


def _overlap_ge_half(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    ov = min(a_end, b_end) - max(a_start, b_start)
    if ov <= 0:
        return False
    min_len = min(a_end - a_start, b_end - b_start)
    return min_len > 0 and ov >= 0.5 * min_len


def annotate_overlap_competition(rows: list[dict]) -> None:
    """Compute `best_on_overlap`, `best_in_namespace`, and `vetoed_by_decoy` in-place per `protein_uid`."""
    by_prot: dict[bytes, list[dict]] = defaultdict(list)
    for r in rows:
        by_prot[r["protein_uid"]].append(r)

    for grp in by_prot.values():
        grp.sort(key=lambda r: (-float(r["bitscore"]), float(r["i_evalue_zref"]), r["model_key"], int(r["env_start"])))
        accepted: list[dict] = []
        accepted_ns: dict[str, list[dict]] = defaultdict(list)
        decoys_seen: list[dict] = []

        for r in grp:
            es, ee = int(r["env_start"]), int(r["env_end"])
            # best_on_overlap across both target and decoy envelopes
            if any(_overlap_ge_half(es, ee, int(a["env_start"]), int(a["env_end"])) for a in accepted):
                r["best_on_overlap"] = False
            else:
                r["best_on_overlap"] = True
                accepted.append(r)

            # best_in_namespace within same namespace
            ns = r["namespace"]
            if any(_overlap_ge_half(es, ee, int(a["env_start"]), int(a["env_end"])) for a in accepted_ns[ns]):
                r["best_in_namespace"] = False
            else:
                r["best_in_namespace"] = True
                accepted_ns[ns].append(r)

            # vetoed_by_decoy: higher-scoring decoy overlapping >= 50%
            if r["is_decoy"]:
                r["vetoed_by_decoy"] = False
                decoys_seen.append(r)
            else:
                r["vetoed_by_decoy"] = any(
                    _overlap_ge_half(es, ee, int(d["env_start"]), int(d["env_end"])) for d in decoys_seen
                )


# --- Sub-step planner (< 200 CPU-min) -------------------------------------------------------------------


def plan_tasks(
    shards: list[dict],
    models: list[dict] | list[ModelSpec],
    *,
    max_cpu_min: float = DEFAULT_MAX_SUBSTEP_CPU_MIN,
) -> list[dict]:
    """Partition `(shard_slice, model_batch)` deterministically so every sub-step estimates <= `max_cpu_min`."""
    if not shards or not models:
        return []
    m_lens = [int(m.model_len if isinstance(m, ModelSpec) else m["model_len"]) for m in models]
    s_res = [int(s.get("n_residues", 0)) for s in shards]

    # First group models into batches that fit even the largest single shard
    max_shard_res = max(s_res) if s_res else 0
    model_batches: list[list[int]] = []
    cur_m: list[int] = []
    cur_mlen = 0
    for idx, ml in enumerate(m_lens):
        est_if_add = max_shard_res * (cur_mlen + ml) * EST_CPU_MIN_PER_RES_M
        if cur_m and est_if_add > max_cpu_min:
            model_batches.append(cur_m)
            cur_m, cur_mlen = [idx], ml
        else:
            cur_m.append(idx)
            cur_mlen += ml
    if cur_m:
        model_batches.append(cur_m)

    # For each model batch, pack consecutive shards while staying <= max_cpu_min
    tasks: list[dict] = []
    for mb in model_batches:
        mb_len = sum(m_lens[i] for i in mb)
        cur_s: list[int] = []
        cur_res = 0
        for s_idx, r in enumerate(s_res):
            est_if_add = (cur_res + r) * mb_len * EST_CPU_MIN_PER_RES_M
            if cur_s and est_if_add > max_cpu_min:
                tasks.append({
                    "task_id": len(tasks),
                    "shard_indices": cur_s,
                    "model_indices": mb,
                    "est_cpu_min": cur_res * mb_len * EST_CPU_MIN_PER_RES_M,
                })
                cur_s, cur_res = [s_idx], r
            else:
                cur_s.append(s_idx)
                cur_res += r
        if cur_s:
            tasks.append({
                "task_id": len(tasks),
                "shard_indices": cur_s,
                "model_indices": mb,
                "est_cpu_min": cur_res * mb_len * EST_CPU_MIN_PER_RES_M,
            })
    return tasks


# --- Streaming multi-HMM scanning -----------------------------------------------------------------------


def _sorted_uid_stream(protein_seq_files: list[Path]) -> Iterator[bytes]:
    """Stream distinct sorted `protein_uid`s from `protein_seq` Parquet files in bounded memory (one row group at a time)."""
    def _file_iter(path: Path) -> Iterator[bytes]:
        pf = pq.ParquetFile(path)
        for rg in range(pf.metadata.num_row_groups):
            col = pf.read_row_group(rg, columns=["protein_uid"], use_threads=False).column(0).to_pylist()
            yield from col

    prev = None
    for uid in heapq.merge(*(_file_iter(p) for p in sorted(protein_seq_files, key=str))):
        if uid != prev:
            prev = uid
            yield uid


def _iter_protein_seq_blocks(
    protein_seq_files: list[Path],
    block_residues: int,
    skip_uids_iter: Iterator[bytes] | None = None,
) -> Iterator[pyhmmer.easel.DigitalSequenceBlock]:
    """Stream distinct `(protein_uid, seq)` from `protein_seq` Parquet files into digitized blocks, skipping `skip_uids_iter`."""
    def _file_iter(path: Path) -> Iterator[tuple[bytes, str]]:
        pf = pq.ParquetFile(path)
        for rg in range(pf.metadata.num_row_groups):
            t = pf.read_row_group(rg, columns=["protein_uid", "seq"], use_threads=False)
            yield from zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist())

    skip_next = next(skip_uids_iter, None) if skip_uids_iter is not None else None
    seqs: list[pyhmmer.easel.TextSequence] = []
    n_res = 0
    prev_uid = None

    for uid, seq in heapq.merge(*(_file_iter(p) for p in sorted(protein_seq_files, key=str)), key=lambda x: x[0]):
        if uid == prev_uid:
            continue
        prev_uid = uid
        while skip_next is not None and skip_next < uid:
            skip_next = next(skip_uids_iter, None)
        if skip_next is not None and skip_next == uid:
            continue
        clean = seq.replace("*", "")
        seqs.append(pyhmmer.easel.TextSequence(name=uid.hex().encode(), sequence=clean))
        n_res += len(clean)
        if n_res >= block_residues:
            yield pyhmmer.easel.TextSequenceBlock(seqs).digitize(engine.AMINO)
            seqs, n_res = [], 0
    if seqs:
        yield pyhmmer.easel.TextSequenceBlock(seqs).digitize(engine.AMINO)


def _iter_fastz_blocks_filtered(
    shard_paths: list[Path],
    block_residues: int,
    skip_uids_iter: Iterator[bytes] | None = None,
) -> Iterator[pyhmmer.easel.DigitalSequenceBlock]:
    """Stream digitized blocks from `.faa.zst` shards, optionally skipping sorted `skip_uids_iter`."""
    if skip_uids_iter is None:
        for p in shard_paths:
            yield from fastz.blocks(p, engine.AMINO, block_residues)
        return

    skip_next = next(skip_uids_iter, None)
    seqs: list[pyhmmer.easel.TextSequence] = []
    n_res = 0
    for p in shard_paths:
        for name, seq in fastz.read_records(p):
            uid = bytes.fromhex(name)
            while skip_next is not None and skip_next < uid:
                skip_next = next(skip_uids_iter, None)
            if skip_next is not None and skip_next == uid:
                continue
            seqs.append(pyhmmer.easel.TextSequence(name=name.encode(), sequence=seq))
            n_res += len(seq)
            if n_res >= block_residues:
                yield pyhmmer.easel.TextSequenceBlock(seqs).digitize(engine.AMINO)
                seqs, n_res = [], 0
    if seqs:
        yield pyhmmer.easel.TextSequenceBlock(seqs).digitize(engine.AMINO)


def _scan_block_multi(
    block: pyhmmer.easel.DigitalSequenceBlock,
    models: list[ModelSpec],
    cpus: int,
    out_rows: list[dict],
    *,
    neg_hmms: list[tuple[str, float, object]] | None = None,
    neg_counts: dict[str, int] | None = None,
) -> None:
    """Run one `pyhmmer.hmmsearch` pass of `models` (+ optional shuffled negative models) on `block`."""
    if not models:
        return
    all_hmms = [m.hmm for m in models]
    min_seq = min(m.seq_bits for m in models)
    min_dom = min(m.dom_bits for m in models)
    if neg_hmms:
        all_hmms.extend(h for _, _, h in neg_hmms)
        min_seq = min(min_seq, *(sb for _, sb, _ in neg_hmms))

    # Pass T/domT floor to pyhmmer so it prunes noise early while retaining every model's passing hits
    kw: dict = {"cpus": cpus, "Z": Z_PROT, "domZ": DOM_Z}
    if min_seq > 0:
        kw["T"] = min_seq
        kw["incT"] = min_seq
    if min_dom > 0:
        kw["domT"] = min_dom
        kw["incdomT"] = min_dom

    n_models = len(models)
    for idx, top_hits in enumerate(pyhmmer.hmmsearch(all_hmms, block, **kw)):
        if idx >= n_models:
            # Shuffled negative control model
            mk, sb, _ = neg_hmms[idx - n_models]
            cnt = sum(1 for h in top_hits if h.score >= sb)
            if neg_counts is not None:
                neg_counts[mk] = neg_counts.get(mk, 0) + cnt
            continue
        m = models[idx]
        for hit in top_hits:
            if hit.score < m.seq_bits:
                continue
            doms = list(hit.domains)
            if not doms:
                continue
            n_mod_doms = len(doms)
            max_dom_b = max(d.score for d in doms)
            dom_seq_frac = float(max_dom_b / hit.score) if hit.score > 0 else 0.0
            if not passes_min_dom_frac(max_dom_b, hit.score, m.min_dom_frac):
                continue
            uid = bytes.fromhex(engine.text(hit.name))
            for d in doms:
                if d.score < m.dom_bits:
                    continue
                if m.i_evalue_max is not None and d.i_evalue > m.i_evalue_max:
                    continue
                a = d.alignment
                hf, ht = a.hmm_from - 1, a.hmm_to
                cov = (ht - hf) / m.model_len if m.model_len else 0.0
                if m.hmm_cov_min is not None and cov < m.hmm_cov_min:
                    continue
                out_rows.append({
                    "protein_uid": uid,
                    "model_key": m.model_key,
                    "namespace": m.namespace,
                    "model_name": m.name,
                    "model_acc": m.acc,
                    "model_len": m.model_len,
                    "is_decoy": m.is_decoy,
                    "seq_bitscore": float(hit.score),
                    "seq_evalue_zref": float(hit.evalue),
                    "bitscore": float(d.score),
                    "i_evalue_zref": float(d.i_evalue),
                    "env_start": int(d.env_from - 1),
                    "env_end": int(d.env_to),
                    "ali_start": int(a.target_from - 1),
                    "ali_end": int(a.target_to),
                    "hmm_from": int(hf),
                    "hmm_to": int(ht),
                    "hmm_cov": float(cov),
                    "n_model_domains": min(65535, int(n_mod_doms)),
                    "dom_seq_frac": dom_seq_frac,
                    "best_on_overlap": False,
                    "best_in_namespace": False,
                    "vetoed_by_decoy": False,
                })


def scan_fastz_shards(
    shard_paths: list[Path],
    models: list[ModelSpec],
    *,
    cpus: int = 1,
    block_residues: int = PANEL_BLOCK_RESIDUES,
    skip_uids_iter: Iterator[bytes] | None = None,
    neg_hmms: list[tuple[str, float, object]] | None = None,
    neg_counts: dict[str, int] | None = None,
) -> list[dict]:
    """One-pass streaming multi-HMM search over `.faa.zst` shards."""
    rows: list[dict] = []
    with contextlib.closing(engine.prefetch(
        lambda: _iter_fastz_blocks_filtered(shard_paths, block_residues, skip_uids_iter)
    )) as blocks:
        for block in blocks:
            _scan_block_multi(block, models, cpus, rows, neg_hmms=neg_hmms, neg_counts=neg_counts)
            del block
    return rows


def scan_protein_seq_files(
    protein_seq_files: list[Path],
    models: list[ModelSpec],
    *,
    cpus: int = 1,
    block_residues: int = PANEL_BLOCK_RESIDUES,
    skip_uids_iter: Iterator[bytes] | None = None,
    neg_hmms: list[tuple[str, float, object]] | None = None,
    neg_counts: dict[str, int] | None = None,
) -> list[dict]:
    """One-pass streaming multi-HMM search over `protein_seq` Parquet files (for `--stages` subsets)."""
    rows: list[dict] = []
    with contextlib.closing(engine.prefetch(
        lambda: _iter_protein_seq_blocks(protein_seq_files, block_residues, skip_uids_iter)
    )) as blocks:
        for block in blocks:
            _scan_block_multi(block, models, cpus, rows, neg_hmms=neg_hmms, neg_counts=neg_counts)
            del block
    return rows


# --- Circular-aware co-occurring system / operon caller + 4 diagnostics ---------------------------------


def _is_cas_domain(name: str) -> bool:
    return bool(CAS_DOMAIN_RE.search(name or ""))


def _has_split_model_pair(cluster_slots: list[int], hits_by_slot: dict[int, list[dict]], tape: dict) -> bool:
    """True iff same-strand adjacent hit ORFs hit non-overlapping N- and C-terminal halves of the same model family."""
    strands = tape["strands"]
    n_slots = int(tape["n_slots"])
    is_circ = bool(tape["is_circular"])
    for idx in range(len(cluster_slots) - 1):
        s_a, s_b = cluster_slots[idx], cluster_slots[idx + 1]
        # Adjacent or at most 1 intervening ORF (to catch short frameshift/intron bridge)
        gap = (s_b - s_a - 1) if s_b > s_a else ((s_b + n_slots - s_a - 1) if is_circ else 999)
        if gap > 1:
            continue
        if strands[s_a] != strands[s_b]:
            continue
        strand = strands[s_a]
        # Transcription order: on +1, s_a is upstream (N-term), s_b is downstream (C-term); on -1 reversed
        up_hits, dn_hits = (hits_by_slot[s_a], hits_by_slot[s_b]) if strand == 1 else (hits_by_slot[s_b], hits_by_slot[s_a])
        for ha in up_hits:
            for hb in dn_hits:
                same_family = (
                    ha["model_key"] == hb["model_key"]
                    or (ha["model_name"] and ha["model_name"] == hb["model_name"])
                )
                if not same_family:
                    continue
                af, at = int(ha["hmm_from"]), int(ha["hmm_to"])
                bf, bt = int(hb["hmm_from"]), int(hb["hmm_to"])
                if af >= bf:
                    continue
                ov = max(0, min(at, bt) - max(af, bf))
                min_span = min(at - af, bt - bf)
                if min_span > 0 and ov < 0.2 * min_span:
                    return True
    return False


def _igr_and_bridged_diagnostics(
    tape: dict,
    cluster_slots: list[int],
    igrs: list[dict],
    protein_lengths: dict[bytes, int],
    protein_domains: dict[bytes, set[str]],
    hit_slot_set: set[int],
) -> tuple[int, str, int]:
    """Compute (max_adjacent_igr_bp, adjacent_igr_uids, bridged_noncoding_bp) within +-2 slots of cluster."""
    n_slots = int(tape["n_slots"])
    is_circ = bool(tape["is_circular"])
    f_uids = tape["feature_uids"]
    p_uids = tape["protein_uids"]
    starts = tape["start_bps"]
    ends = tape["end_bps"]
    contig_len = int(tape["contig_len"])

    # Collect feature_uids within +-2 slots of any slot in the cluster
    near_slots: set[int] = set()
    for s in cluster_slots:
        for d in range(-2, 3):
            if is_circ:
                near_slots.add((s + d) % n_slots)
            elif 0 <= s + d < n_slots:
                near_slots.add(s + d)
    near_fus = {f_uids[s] for s in near_slots}

    adj_igrs = [g for g in igrs if g.get("left_feature_uid") in near_fus or g.get("right_feature_uid") in near_fus]
    max_igr_bp = 0
    adj_uids: list[str] = []
    for g in adj_igrs:
        length = int(g["end"]) - int(g["start"])
        if length > max_igr_bp:
            max_igr_bp = length
        ruid = ids.region_uid(g["contig_uid"], int(g["start"]), int(g["end"]), int(g["strand"]))
        adj_uids.append(ruid.hex())
    adj_uids = sorted(set(adj_uids))

    # Bridged noncoding bp: bridge across adjacent unannotated <100 aa smORFs starting from cluster boundaries
    def _is_unannotated_smorf(slot: int) -> bool:
        if slot in hit_slot_set:
            return False
        pu = p_uids[slot]
        if pu is None:
            return True
        aa_len = protein_lengths.get(pu)
        if aa_len is None:
            aa_len = max(0, (int(ends[slot]) - int(starts[slot])) // 3)
        if aa_len >= SMORF_MAX_AA:
            return False
        return not bool(protein_domains.get(pu))

    bridged_max = max_igr_bp

    # Walk outward from cluster's first slot (left) and last slot (right) across consecutive unannotated smORFs
    if cluster_slots:
        # Rightward from cluster_slots[-1]
        anchor_r = cluster_slots[-1]
        cur = anchor_r
        anchor_end = int(ends[anchor_r])
        visited = {anchor_r}
        for _ in range(min(n_slots, 20)):
            nxt = (cur + 1) % n_slots if is_circ else cur + 1
            if nxt < 0 or nxt >= n_slots or nxt in visited:
                break
            visited.add(nxt)
            nxt_start = int(starts[nxt])
            if is_circ and nxt_start < anchor_end:
                nxt_start += contig_len
            span = nxt_start - anchor_end
            if span > bridged_max:
                bridged_max = span
            if not _is_unannotated_smorf(nxt):
                break
            cur = nxt

        # Leftward from cluster_slots[0]
        anchor_l = cluster_slots[0]
        cur = anchor_l
        anchor_start = int(starts[anchor_l])
        visited = {anchor_l}
        for _ in range(min(n_slots, 20)):
            prv = (cur - 1) % n_slots if is_circ else cur - 1
            if prv < 0 or prv >= n_slots or prv in visited:
                break
            visited.add(prv)
            prv_end = int(ends[prv])
            if is_circ and prv_end > anchor_start:
                prv_end -= contig_len
            span = anchor_start - prv_end
            if span > bridged_max:
                bridged_max = span
            if not _is_unannotated_smorf(prv):
                break
            cur = prv

    return max_igr_bp, ",".join(adj_uids), max(0, bridged_max)


def call_loci_for_contig(
    tape: dict,
    orf_hits: list[dict],
    *,
    max_gap_genes: int = 5,
    max_gap_bp: int = 10000,
    min_models: int = 2,
    distinct_orfs: bool = False,
    require: set[str] | None = None,
    forbid: set[str] | None = None,
    all_overlaps: bool = False,
    neighbors: int = 5,
    igrs: list[dict] | None = None,
    protein_lengths: dict[bytes, int] | None = None,
    protein_domains: dict[bytes, set[str]] | None = None,
    crispr_arrays: list[dict] | None = None,
) -> list[dict]:
    """Call circular-aware co-occurring system clusters on one contig tape."""
    igrs = igrs or []
    protein_lengths = protein_lengths or {}
    protein_domains = protein_domains or {}
    crispr_arrays = crispr_arrays or []
    req_set = set(require or ())
    forbid_set = set(forbid or ())

    # Filter hits to non-decoy (+ best_on_overlap unless all_overlaps)
    usable = [h for h in orf_hits if not h["is_decoy"] and (all_overlaps or h["best_on_overlap"])]
    if not usable:
        return []

    hits_by_slot: dict[int, list[dict]] = defaultdict(list)
    for h in usable:
        hits_by_slot[int(h["ordinal"])].append(h)

    hit_slots = sorted(hits_by_slot)
    n_slots = int(tape["n_slots"])
    contig_len = int(tape["contig_len"])
    is_circ = bool(tape["is_circular"])
    starts = tape["start_bps"]
    ends = tape["end_bps"]
    strands = tape["strands"]
    f_uids = tape["feature_uids"]
    p_uids = tape["protein_uids"]

    # Linear sweep to form initial clusters
    clusters: list[list[int]] = [[hit_slots[0]]]
    for s in hit_slots[1:]:
        prev = clusters[-1][-1]
        g_gap = s - prev - 1
        b_gap = int(starts[s]) - int(ends[prev])
        if g_gap <= max_gap_genes and b_gap <= max_gap_bp:
            clusters[-1].append(s)
        else:
            clusters.append([s])

    # Circular origin wrap-around merge
    wraps_origin = False
    if is_circ and len(clusters) >= 2:
        last_s = clusters[-1][-1]
        first_s = clusters[0][0]
        wrap_g = (first_s + n_slots - last_s) - 1
        wrap_b = (int(starts[first_s]) + contig_len) - int(ends[last_s])
        if wrap_g <= max_gap_genes and wrap_b <= max_gap_bp:
            merged = clusters[-1] + clusters[0]
            clusters = [merged] + clusters[1:-1]
            wraps_origin = True

    out: list[dict] = []
    for c_idx, c_slots in enumerate(clusters):
        c_wraps = wraps_origin and c_idx == 0
        c_hits = [h for s in c_slots for h in hits_by_slot[s]]
        models_set = {h["model_key"] for h in c_hits}
        if len(models_set) < min_models:
            continue
        if distinct_orfs and len(c_slots) < min_models:
            continue
        if req_set and not (req_set <= models_set):
            continue

        # Build window slot set (+-neighbors) for forbid / Cas-less CRISPR checks
        win_slots: set[int] = set()
        for s in c_slots:
            for d in range(-neighbors, neighbors + 1):
                if is_circ:
                    win_slots.add((s + d) % n_slots)
                elif 0 <= s + d < n_slots:
                    win_slots.add(s + d)

        if forbid_set:
            # Check panel hits and Pfam/domain hits in cluster + window
            win_models = set(models_set)
            for s in win_slots:
                for h in hits_by_slot.get(s, ()):
                    win_models.add(h["model_key"])
                    win_models.add(h["model_name"])
                pu = p_uids[s]
                if pu is not None:
                    win_models |= protein_domains.get(pu, set())
            if forbid_set & win_models:
                continue

        # Anchor: highest bitscore hit (preferring first required model if --require given)
        cand_anchor_hits = [h for h in c_hits if h["model_key"] in req_set] if req_set else c_hits
        anchor_hit = min(cand_anchor_hits, key=lambda h: (-float(h["bitscore"]), h["feature_uid"]))
        anchor_fu = anchor_hit["feature_uid"]

        # Ordered models string e.g. TniQ(+)>Cas7f(+)
        slot_labels = []
        fusion_orfs = 0
        for s in c_slots:
            sh = sorted(hits_by_slot[s], key=lambda h: (int(h.get("env_start", 0)), h["model_key"]))
            uniq_m = list(dict.fromkeys(h["model_key"] for h in sh))
            if len(uniq_m) > 1:
                fusion_orfs += 1
            sign = "+" if strands[s] == 1 else "-"
            slot_labels.append(f"{'+'.join(uniq_m)}({sign})")
        models_ordered = ">".join(slot_labels)

        # Span and unmatched internal ORFs
        first_s, last_s = c_slots[0], c_slots[-1]
        c_start = int(starts[first_s])
        c_end = int(ends[last_s])
        if c_wraps and c_end <= c_start:
            span_bp = (c_end + contig_len) - c_start
            total_span_slots = (last_s + n_slots - first_s) + 1
        else:
            span_bp = c_end - c_start
            total_span_slots = (last_s - first_s) + 1
        n_unmatched_internal = max(0, total_span_slots - len(c_slots))

        # Censoring: linear contig whose +-neighbors window hits a contig end
        if is_circ:
            censored = False
        else:
            censored = (min(c_slots) < neighbors) or (n_slots - 1 - max(c_slots) < neighbors)

        # CRISPR arrays overlapping cluster +-neighbors window
        w_min = min(win_slots)
        w_max = max(win_slots)
        w_lo = int(starts[w_min])
        w_hi = int(ends[w_max])
        n_arrays = 0
        for arr in crispr_arrays:
            a_s, a_e = int(arr["start"]), int(arr["end"])
            if is_circ and c_wraps:
                # Any array in [c_start, contig_len) or [0, c_end)
                if a_e > c_start or a_s < c_end:
                    n_arrays += 1
            elif min(a_e, w_hi) > max(a_s, w_lo):
                n_arrays += 1

        # Cas-less CRISPR: n_crispr_arrays > 0 and zero Cas domains within +-10 genes
        cas_less = False
        if n_arrays > 0:
            cas_slots: set[int] = set()
            for s in c_slots:
                for d in range(-10, 11):
                    if is_circ:
                        cas_slots.add((s + d) % n_slots)
                    elif 0 <= s + d < n_slots:
                        cas_slots.add(s + d)
            has_cas = False
            for s in cas_slots:
                for h in hits_by_slot.get(s, ()):
                    if _is_cas_domain(h["model_name"]) or _is_cas_domain(h["model_key"]):
                        has_cas = True
                        break
                if has_cas:
                    break
                pu = p_uids[s]
                if pu is not None and any(_is_cas_domain(d) for d in protein_domains.get(pu, ())):
                    has_cas = True
                    break
            cas_less = not has_cas

        max_igr_bp, adj_igr_uids, bridged_bp = _igr_and_bridged_diagnostics(
            tape, c_slots, igrs, protein_lengths, protein_domains, set(c_slots),
        )
        split_pair = _has_split_model_pair(c_slots, hits_by_slot, tape)

        out.append({
            "contig_uid": tape["contig_uid"],
            "anchor_feature_uid": anchor_fu,
            "contig_accession": anchor_hit.get("contig_accession"),
            "stage": anchor_hit.get("stage", ""),
            "genome_kind": anchor_hit.get("genome_kind"),
            "is_circular": is_circ,
            "contig_len": contig_len,
            "cluster_start": c_start,
            "cluster_end": c_end,
            "span_bp": int(span_bp),
            "censored": censored,
            "models_present": ",".join(sorted(models_set)),
            "models_ordered": models_ordered,
            "n_models": len(models_set),
            "n_panel_orfs": len(c_slots),
            "n_unmatched_internal_orfs": n_unmatched_internal,
            "fusion_orfs": fusion_orfs,
            "n_crispr_arrays": n_arrays,
            "max_adjacent_igr_bp": max_igr_bp,
            "adjacent_igr_uids": adj_igr_uids,
            "bridged_noncoding_bp": bridged_bp,
            "cas_less_crispr": cas_less,
            "split_model_pair": split_pair,
        })

    out.sort(key=lambda r: (r["contig_uid"], r["anchor_feature_uid"]))
    return out


# --- Parquet & TSV writers (deterministic, row-group stats) ---------------------------------------------


def _write_parquet_sorted(rows: list[dict], schema: pa.Schema, path: Path) -> None:
    """Write rows atomically with `ROW_GROUP_ROWS` row groups and min/max statistics enabled."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with pq.ParquetWriter(tmp, schema, compression="zstd", write_statistics=True) as w:
            if not rows:
                w.write_table(pa.Table.from_pylist([], schema=schema))
            else:
                for i in range(0, len(rows), ROW_GROUP_ROWS):
                    w.write_table(pa.Table.from_pylist(rows[i : i + ROW_GROUP_ROWS], schema=schema))
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _write_tsv(rows: list[dict], columns: list[str], path: Path) -> None:
    lines = ["\t".join(columns)]
    for r in rows:
        lines.append("\t".join(export.cell(c, r.get(c)) for c in columns))
    export.write_lines(lines, path)


# --- SQL panel mounting (`bio-box sql --panel NAME=DIR`) ------------------------------------------------


def validate_panel_name(name: str) -> str:
    """Validate `--panel NAME=DIR` identifier and refuse collision with `schemas.TABLES`."""
    if not PANEL_NAME_RE.match(name):
        raise ValueError(f"invalid panel name {name!r}: must match ^[A-Za-z_][A-Za-z0-9_]*$")
    if name in TABLES:
        raise ValueError(f"panel name {name!r} collides with a release table in schemas.TABLES")
    for suffix in ("models", "protein_hits", "orf_hits", "loci"):
        if f"{name}_{suffix}" in TABLES:
            raise ValueError(f"panel view {name}_{suffix!r} collides with a release table in schemas.TABLES")
    return name


PANEL_VIEW_FILES = (
    ("models", "models.parquet"),
    ("protein_hits", "protein_hits.parquet"),
    ("orf_hits", "orf_hits.parquet"),
    ("loci", "loci_summary.parquet"),
)


def resolve_panel_mount(spec: str) -> tuple[str, Path]:
    """Parse and validate `--panel NAME=DIR`."""
    if "=" not in spec:
        raise ValueError(f"--panel {spec!r}: expected NAME=DIR")
    name, d = spec.split("=", 1)
    validate_panel_name(name)
    pdir = Path(d)
    if not pdir.is_dir():
        raise ValueError(f"--panel {name}: directory {pdir} does not exist")
    for _, fname in PANEL_VIEW_FILES:
        if not (pdir / fname).is_file():
            raise ValueError(f"--panel {name}: {pdir} is missing {fname}")
    return name, pdir


def panel_parquet_paths(panels: list[tuple[str, Path]]) -> list[Path]:
    out: list[Path] = []
    for _, pdir in panels:
        for _, fname in PANEL_VIEW_FILES:
            out.append(pdir / fname)
    return out


def register_panel_views(con, panels: list[tuple[str, Path]], release: str) -> None:
    """Register read-only DuckDB views `<NAME>_models`, `<NAME>_protein_hits`, `<NAME>_orf_hits`, `<NAME>_loci`."""
    for name, pdir in panels:
        cov_path = pdir / "coverage.json"
        if cov_path.is_file():
            try:
                cov = json.loads(cov_path.read_text(encoding="utf-8"))
                p_rel = cov.get("release")
                if p_rel and p_rel != release:
                    print(
                        f"note: panel {name} ({pdir}) was built on release {p_rel}, not {release}",
                        file=sys.stderr,
                    )
            except (OSError, ValueError):
                pass
        for suffix, fname in PANEL_VIEW_FILES:
            fpath = str((pdir / fname).resolve())
            con.execute(
                f'CREATE TEMP VIEW "{name}_{suffix}" AS '
                f"SELECT * FROM read_parquet([{export._q(fpath)}], hive_partitioning = false)"
            )


# --- Carry-forward (`--carry-from PREV_RUN`) ------------------------------------------------------------


def _load_carry_from(prev_dir: Path, models: list[ModelSpec], root: Path) -> tuple[set[str], list[Path], list[dict]]:
    """Inspect `PREV_RUN` and return (reusable_model_keys, covered_protein_seq_paths, carried_hit_rows)."""
    cov_path = prev_dir / "coverage.json"
    mod_path = prev_dir / "models.parquet"
    hits_path = prev_dir / "protein_hits.parquet"
    if not (cov_path.is_file() and mod_path.is_file() and hits_path.is_file()):
        raise PanelError(f"--carry-from {prev_dir}: missing coverage.json, models.parquet, or protein_hits.parquet")

    cov = json.loads(cov_path.read_text(encoding="utf-8"))
    cur_tools = provenance.tool_versions("pyhmmer", "pyarrow")
    if (
        float(cov.get("Z_PROT", 0)) != float(Z_PROT)
        or float(cov.get("DOM_Z", 0)) != float(DOM_Z)
        or cov.get("tool_versions") != cur_tools
    ):
        return set(), [], []

    # Verify all covered protein_seq files still exist on disk
    covered_paths: list[Path] = []
    for entry in cov.get("protein_seq_files", []):
        p = Path(root) / entry["path"]
        if not p.is_file():
            return set(), [], []
        covered_paths.append(p)

    prev_models = {r["model_key"]: r for r in pq.read_table(mod_path).to_pylist()}
    reusable: set[str] = set()
    for m in models:
        pm = prev_models.get(m.model_key)
        if pm is None:
            continue
        prev_sig = (
            pm["model_key"],
            pm["hmm_sha256"],
            bool(pm["is_decoy"]),
            pm["cutoff_mode"],
            round(float(pm["seq_bits"]), 6),
            round(float(pm["dom_bits"]), 6),
            None if pm["i_evalue_max"] is None else float(pm["i_evalue_max"]),
            None if pm["hmm_cov_min"] is None else round(float(pm["hmm_cov_min"]), 6),
            None if pm["min_dom_frac"] is None else round(float(pm["min_dom_frac"]), 6),
        )
        if prev_sig == m.signature():
            reusable.add(m.model_key)

    if not reusable:
        return set(), [], []

    # Stream-copy hits for reusable models from PREV_RUN/protein_hits.parquet
    carried: list[dict] = []
    pf = pq.ParquetFile(hits_path)
    want_arr = pa.array(sorted(reusable), pa.string())
    for rg in range(pf.metadata.num_row_groups):
        t = pf.read_row_group(rg, use_threads=False)
        t = t.filter(pc.is_in(t.column("model_key"), value_set=want_arr))
        if t.num_rows:
            carried.extend(t.to_pylist())
    return reusable, covered_paths, carried


# --- Sub-step execution helper (bounded CPU-min, STOP check, atomic partials) ---------------------------


def _same_fs(a: Path, b: Path) -> bool:
    try:
        return os.stat(a).st_dev == os.stat(b).st_dev
    except OSError:
        return False


def _run_substep(
    task: dict,
    *,
    shards_meta: list[dict],
    pack_dir: Path | None,
    protein_seq_files: list[Path] | None,
    models: list[ModelSpec],
    work_dir: Path,
    partial_path: Path,
    cpus: int,
    skip_covered_paths: list[Path] | None,
    neg_hmms: list[tuple[str, float, object]] | None,
) -> tuple[list[dict], dict[str, int]]:
    """Execute one `(shard_slice, model_batch)` sub-step and write its atomic partial Parquet."""
    sub_models = [models[i] for i in task["model_indices"]]
    sub_keys = {m.model_key for m in sub_models}
    sub_neg = [x for x in (neg_hmms or []) if x[0] in sub_keys]
    neg_counts: dict[str, int] = {}

    if partial_path.is_file() and not sub_neg:
        return pq.read_table(partial_path).to_pylist(), neg_counts

    skip_iter = _sorted_uid_stream(skip_covered_paths) if skip_covered_paths else None

    if protein_seq_files is not None:
        # Stage-subset streaming directly from l0seq/protein_seq
        sub_files = [protein_seq_files[i] for i in task["shard_indices"]]
        rows = scan_protein_seq_files(
            sub_files, sub_models, cpus=cpus, skip_uids_iter=skip_iter,
            neg_hmms=sub_neg or None, neg_counts=neg_counts,
        )
    else:
        assert pack_dir is not None
        sub_shards = [shards_meta[i] for i in task["shard_indices"]]
        same_fs = _same_fs(pack_dir, work_dir)
        staged_tmp: list[Path] = []
        scan_paths: list[Path] = []
        try:
            for s in sub_shards:
                src = pack_dir / s["path"]
                if same_fs:
                    if ids.sha256_file(src).hex() != s["sha256"]:
                        raise PanelError(f"pack shard {src} sha256 differs from pack_manifest.json")
                    scan_paths.append(src)
                else:
                    dst = work_dir / f"stage_t{task['task_id']}" / Path(s["path"]).name
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(src, dst)
                    if ids.sha256_file(dst).hex() != s["sha256"]:
                        raise PanelError(f"staged shard {dst} sha256 differs from pack_manifest.json")
                    staged_tmp.append(dst)
                    scan_paths.append(dst)
            rows = scan_fastz_shards(
                scan_paths, sub_models, cpus=cpus, skip_uids_iter=skip_iter,
                neg_hmms=sub_neg or None, neg_counts=neg_counts,
            )
        finally:
            for p in staged_tmp:
                p.unlink(missing_ok=True)

    rows.sort(key=lambda r: (r["protein_uid"], r["model_key"], r["env_start"]))
    _write_parquet_sorted(rows, PROTEIN_HITS_SCHEMA, partial_path)
    return rows, neg_counts


# --- ORF expansion (bounded DuckDB row-group-pruning joins) ---------------------------------------------


def _map_domain_state(raw: str | None) -> str:
    if raw in ("measured", "inferred"):
        return "annotated"
    if raw == "none":
        return "annotated_empty"
    return "unscanned"


def _expand_to_orf_hits(
    root: Path,
    release: str,
    files,
    protein_hits: list[dict],
    *,
    gene_call_version: str,
    tape_version: str,
    stages_filter: list[str] | None,
    include_decoys: bool,
    library: str | None,
    threads: int | None,
    tmp: Path | None,
) -> list[dict]:
    """Expand `protein_hits` to `orf_hits` via bounded, row-group-pruning DuckDB joins."""
    kept_ph = [r for r in protein_hits if include_decoys or not r["is_decoy"]]
    if not kept_ph:
        return []

    hit_uids = sorted({r["protein_uid"] for r in kept_ph})
    hits_by_prot: dict[bytes, list[dict]] = defaultdict(list)
    for r in kept_ph:
        hits_by_prot[r["protein_uid"]].append(r)

    with export.open_store(
        root, release=release, gene_call_version=gene_call_version, tape_version=tape_version,
        manifest=files, threads=threads, tmp=tmp,
    ) as st:
        con = st.con
        lib = st.library(library)
        con.register("_p_in", pa.table({"protein_uid": pa.array(hit_uids, UID)}))
        con.execute("CREATE TEMP TABLE _p AS SELECT protein_uid FROM _p_in")
        con.unregister("_p_in")

        # Add feature_protein view for genetic_code
        fp_files = nbr._store_files(root, "feature_protein", {"gene_call_version": gene_call_version}, files)
        export._view(con, "fp", "feature_protein", fp_files, ("feature_uid", "genetic_code"), ())

        # Add contig_dim tax_order view
        cd_files = nbr._store_files(root, "contig_dim", {"release": release}, files)
        export._view(con, "cdim_tax", "contig_dim", cd_files, ("contig_uid", "is_circular", "tax_order"), ())

        stage_clause = ""
        if stages_filter:
            in_stages = ", ".join(export._q(s) for s in stages_filter)
            stage_clause = f"AND o.stage IN ({in_stages})"

        con.execute(f"""
            CREATE TEMP TABLE _occ AS
            SELECT o.protein_uid, o.contig_uid, o.feature_uid, o.start, o."end", o.strand, o.stage
            FROM occ o
            WHERE o.protein_uid IN (SELECT protein_uid FROM _p) {stage_clause}
        """)
        con.execute("CREATE TEMP TABLE _c AS SELECT DISTINCT contig_uid FROM _occ")
        export._annotate_contigs(con, with_ordinals=True, genomes=False)
        export._annotate_proteins(con, lib, st.stages_of(lib))

        con.execute("""
            CREATE TEMP TABLE _fp_code AS
            SELECT feature_uid, min(genetic_code) AS genetic_code
            FROM fp WHERE feature_uid IN (SELECT feature_uid FROM _occ)
            GROUP BY feature_uid
        """)
        con.execute("""
            CREATE TEMP TABLE _kind AS
            SELECT contig_uid, stage, min(kind) AS kind
            FROM gin WHERE contig_uid IN (SELECT contig_uid FROM _c)
            GROUP BY contig_uid, stage
        """)
        con.execute("""
            CREATE TEMP TABLE _ctax AS
            SELECT contig_uid, bool_or(is_circular) AS is_circular, min(tax_order) AS tax_order
            FROM cdim_tax WHERE contig_uid IN (SELECT contig_uid FROM _c)
            GROUP BY contig_uid
        """)
        # Fallback is_circular from contig_tapes when contig_dim is absent
        tp_files = nbr._store_files(root, "contig_tapes", {"tape_version": tape_version}, files)
        export._view(con, "tp_circ", "contig_tapes", tp_files, ("contig_uid", "contig_len", "is_circular"), ())
        con.execute("""
            CREATE TEMP TABLE _tc AS
            SELECT contig_uid, max(contig_len) AS contig_len, bool_or(is_circular) AS is_circular
            FROM tp_circ WHERE contig_uid IN (SELECT contig_uid FROM _c)
            GROUP BY contig_uid
        """)

        sql = """
            SELECT o.contig_uid, o.start, o."end", o.strand, o.feature_uid,
                   ord.ordinal, l.accession AS contig_accession,
                   coalesce(tc.contig_len, ci.contig_len) AS contig_len,
                   coalesce(ct.is_circular, tc.is_circular, false) AS is_circular,
                   fc.genetic_code, o.stage, k.kind AS genome_kind,
                   CASE WHEN ct.tax_order IS NOT NULL THEN CAST(ct.tax_order AS VARCHAR) ELSE NULL END AS taxonomy,
                   o.protein_uid, cl.nr90 AS nr90_uid, cl.nr50 AS nr50_uid,
                   coalesce(ds.domain_state, 'not scanned') AS raw_domain_state,
                   a.arch AS pfam_architecture
            FROM _occ o
            JOIN _ord ord ON ord.contig_uid = o.contig_uid AND ord.feature_uid = o.feature_uid
            LEFT JOIN _lab l ON l.contig_uid = o.contig_uid
            LEFT JOIN _clen ci ON ci.contig_uid = o.contig_uid
            LEFT JOIN _tc tc ON tc.contig_uid = o.contig_uid
            LEFT JOIN _ctax ct ON ct.contig_uid = o.contig_uid
            LEFT JOIN _fp_code fc ON fc.feature_uid = o.feature_uid
            LEFT JOIN _kind k ON k.contig_uid = o.contig_uid AND k.stage = o.stage
            LEFT JOIN _cl cl ON cl.protein_uid = o.protein_uid
            LEFT JOIN _dstate ds ON ds.protein_uid = o.protein_uid
            LEFT JOIN _arch a ON a.protein_uid = o.protein_uid
            ORDER BY o.contig_uid, o.start, o.feature_uid
        """
        orf_rows: list[dict] = []
        for g in export._stream(con, sql):
            dstate = _map_domain_state(g["raw_domain_state"])
            for h in hits_by_prot[g["protein_uid"]]:
                orf_rows.append({
                    "contig_uid": g["contig_uid"],
                    "start": int(g["start"]),
                    "end": int(g["end"]),
                    "strand": int(g["strand"]),
                    "feature_uid": g["feature_uid"],
                    "ordinal": int(g["ordinal"]),
                    "tape_version": tape_version,
                    "contig_accession": g["contig_accession"],
                    "contig_len": int(g["contig_len"]),
                    "is_circular": bool(g["is_circular"]),
                    "genetic_code": int(g["genetic_code"]) if g["genetic_code"] is not None else None,
                    "stage": g["stage"],
                    "genome_kind": g["genome_kind"],
                    "taxonomy": g["taxonomy"],
                    "protein_uid": g["protein_uid"],
                    "nr90_uid": g["nr90_uid"],
                    "nr50_uid": g["nr50_uid"],
                    "domain_state": dstate,
                    "pfam_architecture": g["pfam_architecture"],
                    "model_key": h["model_key"],
                    "namespace": h["namespace"],
                    "model_name": h["model_name"],
                    "model_acc": h["model_acc"],
                    "model_len": h["model_len"],
                    "is_decoy": h["is_decoy"],
                    "seq_bitscore": h["seq_bitscore"],
                    "seq_evalue_zref": h["seq_evalue_zref"],
                    "bitscore": h["bitscore"],
                    "i_evalue_zref": h["i_evalue_zref"],
                    "env_start": h["env_start"],
                    "env_end": h["env_end"],
                    "ali_start": h["ali_start"],
                    "ali_end": h["ali_end"],
                    "hmm_from": h["hmm_from"],
                    "hmm_to": h["hmm_to"],
                    "hmm_cov": h["hmm_cov"],
                    "n_model_domains": h["n_model_domains"],
                    "dom_seq_frac": h["dom_seq_frac"],
                    "best_on_overlap": h["best_on_overlap"],
                    "best_in_namespace": h["best_in_namespace"],
                    "vetoed_by_decoy": h["vetoed_by_decoy"],
                })
    orf_rows.sort(key=lambda r: (r["contig_uid"], r["start"], r["feature_uid"], r["model_key"], r["env_start"]))
    return orf_rows


# --- Full panel runner ----------------------------------------------------------------------------------


def _load_igrs_for_contigs(files, contigs: set[bytes]) -> dict[bytes, list[dict]]:
    out: dict[bytes, list[dict]] = defaultdict(list)
    if not contigs:
        return out
    sorted_c = sorted(contigs)
    for f in files.table("igr"):
        t = pq.read_table(f, filters=[("contig_uid", "in", sorted_c)])
        for r in t.to_pylist():
            out[r["contig_uid"]].append(r)
    return out


def _load_protein_meta_for_contigs(
    root: Path, release: str, files, tapes_map: dict[bytes, dict],
    hits_by_contig: dict[bytes, list[dict]], neighbors: int,
    gene_call_version: str, threads: int | None, tmp: Path | None,
) -> tuple[dict[bytes, int], dict[bytes, set[str]]]:
    """Load (protein_lengths, protein_domains) for proteins within +-max(10, neighbors+2) slots of hit ORFs
    in bounded chunks using row-group-pruning DuckDB joins (AGENTS §4.1, §4.4)."""
    rad = max(10, int(neighbors) + 2)
    near_puids: set[bytes] = set()
    for cu, chits in hits_by_contig.items():
        t = tapes_map.get(cu)
        if t is None:
            continue
        n_slots = int(t["n_slots"])
        is_circ = bool(t["is_circular"])
        p_uids = t["protein_uids"]
        for h in chits:
            s = int(h["ordinal"])
            for d in range(-rad, rad + 1):
                idx = (s + d) % n_slots if is_circ else (s + d)
                if 0 <= idx < n_slots:
                    pu = p_uids[idx]
                    if pu is not None:
                        near_puids.add(pu)
    all_puids = sorted(near_puids)
    if not all_puids:
        return {}, {}
    lengths: dict[bytes, int] = {}
    doms: dict[bytes, set[str]] = defaultdict(set)
    chunk_size = 50_000
    with export.open_store(
        root, release=release, gene_call_version=gene_call_version, manifest=files,
        threads=threads, tmp=tmp, relations=("prot", "dom"),
    ) as st:
        con = st.con
        for i in range(0, len(all_puids), chunk_size):
            batch = all_puids[i : i + chunk_size]
            con.register("_p_in", pa.table({"protein_uid": pa.array(batch, UID)}))
            con.execute("CREATE OR REPLACE TEMP TABLE _p AS SELECT protein_uid FROM _p_in")
            con.unregister("_p_in")
            for r in export._stream(
                con,
                "SELECT protein_uid, max(length) AS length FROM prot WHERE protein_uid IN (SELECT protein_uid FROM _p) GROUP BY protein_uid",
            ):
                lengths[r["protein_uid"]] = int(r["length"])
            for r in export._stream(
                con,
                "SELECT DISTINCT protein_uid, model_name FROM dom WHERE protein_uid IN (SELECT protein_uid FROM _p)",
            ):
                if r["model_name"]:
                    doms[r["protein_uid"]].add(r["model_name"])
    return lengths, dict(doms)


def _panel_readme(has_arrays: bool) -> str:
    txt = """\
Files written by `bio-box search panel` (README.txt; card W-39)

models.parquet / models.tsv         per-model metadata, cutoffs, hit counts, and control outcomes.
protein_hits.parquet                one row per passing (protein_uid, model_key, domain) with overlap & decoy flags.
orf_hits.parquet / orf_hits.tsv     non-decoy per-ORF hits joined with coordinates, tape ordinal, clusters, and Pfam.
loci_summary.parquet / .tsv         circular-aware co-occurring system / operon clusters + 4 RT/RNAP diagnostics.
loci/                               GUIDE.md-compatible per-locus tables (locus_<accession>_<anchor_feature_uid>.tsv).
controls.tsv                        positive, negative-clade, and shuffled-consensus control results.
coverage.json                       covered release, stages, Z constants, tool versions, and protein_seq file list.
unannotated_proteins.txt            unscanned protein_uids in locus windows (for `bio-box annotate round`).
RUN                                 run id and release (`bio-box results archive --under` finds this directory).
"""
    if has_arrays:
        txt += export.README_ARRAYS
    return txt


def run_panel(
    root: Path,
    release: str,
    *,
    hmm_specs: list[str],
    library_specs: list[str],
    decoy_hmm_specs: list[str],
    decoy_library_specs: list[str],
    out_dir: Path,
    packs_dir: Path,
    runs_dir: Path,
    work_base: Path,
    positives: list[tuple[str | None, bytes, str]],
    neg_controls: list[tuple[str | None, bytes, str]] = (),
    negative_mode: str = "shuffled_seed",
    cutoff: str | None = None,
    fallback_bits: float | None = None,
    bits: float | None = None,
    dom_bits: float | None = None,
    cutoffs_tsv: Path | None = None,
    stages: list[str] | None = None,
    max_gap_genes: int = 5,
    max_gap_bp: int = 10000,
    min_models: int | None = None,
    distinct_orfs: bool = False,
    require: list[str] = (),
    forbid: list[str] = (),
    all_overlaps: bool = False,
    include_decoys: bool = False,
    neighbors: int = 5,
    carry_from: Path | None = None,
    cpus: int = 1,
    jobs: int = 1,
    max_substep_cpu_min: float = DEFAULT_MAX_SUBSTEP_CPU_MIN,
    gene_call_version: str = GENE_CALL_VERSION,
    library: str | None = None,
    threads: int | None = None,
    user_command: list[str] | None = None,
    refdb_dirs: Iterable[Path] = (),
) -> tuple[dict, Path]:
    """Execute a complete `bio-box search panel` run and write all outputs atomically. Returns (run_row, run_dir)."""
    t0 = time.monotonic()
    root = Path(root).resolve()
    out_dir = pack.outside_root(root, Path(out_dir), "the output directory")
    packs_dir = pack.outside_root(root, Path(packs_dir), "the pack directory")
    runs_dir = pack.outside_root(root, Path(runs_dir), "the runs directory")
    work_base = pack.outside_root(root, Path(work_base), "the work directory")

    if out_dir.exists() and (not out_dir.is_dir() or any(out_dir.iterdir())):
        raise PanelError(f"output directory {out_dir} is not empty")
    if (root / "STOP").exists():
        raise PanelError(f"STOP file exists in {root}; refusing to start (AGENTS §7.4)")

    catalog.require_not_archived(root, release)
    files = catalog.ReleaseFiles.of(root, release)
    all_stages = list(files.stages)
    if stages is not None:
        bad_s = [s for s in stages if s not in all_stages]
        if bad_s:
            raise PanelError(f"--stages names stage(s) not in release {release}: {', '.join(bad_s)}")
        target_stages = list(dict.fromkeys(stages))
    else:
        target_stages = all_stages
    use_stage_subset = stages is not None and set(target_stages) != set(all_stages)

    pdir = pack.find(packs_dir, release, root=root)
    pman, pid = pack.read_manifest(pdir)
    tv = pack.manifest_tape_version(pman)
    if not use_stage_subset:
        probs, _ = pack.search_check(pdir, release)
        if probs:
            raise PanelError(f"search pack {pdir} fails check: {probs[:3]}")

    models = load_panel_specs(
        hmm_specs, library_specs, decoy_hmm_specs, decoy_library_specs,
        cutoff=cutoff, fallback_bits=fallback_bits, bits=bits, dom_bits=dom_bits,
        cutoffs_tsv=cutoffs_tsv, pack_dir=pdir, root=root, refdb_dirs=refdb_dirs,
    )
    target_models = [m for m in models if not m.is_decoy]
    target_keys = {m.model_key for m in target_models}
    eff_min_models = min_models if min_models is not None else (2 if len(target_models) >= 2 else 1)

    if not positives:
        raise PanelError("search panel requires at least one --positive [MODEL_KEY=]SPEC (S-04)")
    for mk, _, _ in positives:
        if mk is not None and mk not in target_keys:
            raise PanelError(f"--positive model_key {mk!r} is not a non-decoy target model in the panel")

    # Canonical config & run_id
    submitted = datetime.now(timezone.utc)
    panel_cfg = {
        "recipe": PANEL_RECIPE,
        "release": release,
        "stages": target_stages,
        "models": [
            {"model_key": m.model_key, "hmm_sha256": m.hmm_sha256, "is_decoy": m.is_decoy,
             "cutoff_mode": m.cutoff_mode, "seq_bits": m.seq_bits, "dom_bits": m.dom_bits,
             "i_evalue_max": m.i_evalue_max, "hmm_cov_min": m.hmm_cov_min, "min_dom_frac": m.min_dom_frac}
            for m in models
        ],
        "negative": negative_mode,
        "z_prot": Z_PROT,
        "dom_z": DOM_Z,
    }
    cfg_text = json.dumps(panel_cfg, sort_keys=True, indent=2) + "\n"
    rh = hashlib.sha256(cfg_text.encode("utf-8")).hexdigest()
    run_id = runner.new_run_id(rh, release, pid, submitted)
    rdir = runner.run_dir(runs_dir, run_id)
    shards_out = rdir / "shards"
    shards_out.mkdir(parents=True, exist_ok=True)
    work = work_base / run_id.hex()
    work.mkdir(parents=True, exist_ok=True)

    # Carry-forward check
    reusable_keys: set[str] = set()
    covered_pseq_paths: list[Path] = []
    carried_rows: list[dict] = []
    if carry_from is not None:
        reusable_keys, covered_pseq_paths, carried_rows = _load_carry_from(Path(carry_from), models, root)

    # Shuffled-consensus negative controls (deterministic per model_key)
    neg_hmms: list[tuple[str, float, object]] = []
    neg_exp: dict[str, float] = {}
    neg_seed_base = int.from_bytes(hashlib.sha256(rh.encode()).digest()[:4], "big")
    target_n_seq = int(pman["targets"]["unique"]["n_seq"])
    if negative_mode == "shuffled_seed":
        for i, m in enumerate(target_models):
            cons = m.hmm.consensus.upper()
            shuf = engine.shuffled(cons, neg_seed_base + i)
            nh = engine.build_from_sequence(f"shuf_{m.model_key}", shuf)
            neg_hmms.append((m.model_key, m.seq_bits, nh))
            neg_exp[m.model_key] = engine.tail_expectation(nh, m.seq_bits, max(1, target_n_seq))

    # Determine input streams: pack shards (full release) vs stage-filtered protein_seq files
    if use_stage_subset:
        pseq_files = [
            p for p, parts in files.entries("protein_seq", gene_call_version=gene_call_version)
            if str(parts.get("stage")) in target_stages
        ]
        pseq_files.sort(key=str)
        shards_meta = [{"path": str(p), "n_residues": pq.read_metadata(p).num_rows * 300} for p in pseq_files]
    else:
        pseq_files = None
        shards_meta = [s for s in pman["shards"] if s["target"] == "unique"]

    # Split models into (reusable: scan only new UIDs not in covered_pseq_paths) vs (fresh: scan all UIDs)
    reusable_models = [m for m in models if m.model_key in reusable_keys]
    fresh_models = [m for m in models if m.model_key not in reusable_keys]

    all_hits: list[dict] = list(carried_rows)
    neg_counts: dict[str, int] = defaultdict(int)

    def _run_group(group_models: list[ModelSpec], skip_paths: list[Path] | None, tag: str) -> None:
        if not group_models:
            return
        tasks = plan_tasks(shards_meta, group_models, max_cpu_min=max_substep_cpu_min)
        if jobs <= 1 or len(tasks) <= 1:
            for t in tasks:
                if (root / "STOP").exists():
                    raise PanelError(f"STOP file appeared in {root} between sub-steps (AGENTS §7.4)")
                part = shards_out / f"{tag}_task_{t['task_id']:04d}.parquet"
                rows, nc = _run_substep(
                    t, shards_meta=shards_meta, pack_dir=pdir, protein_seq_files=pseq_files,
                    models=group_models, work_dir=work, partial_path=part, cpus=cpus,
                    skip_covered_paths=skip_paths, neg_hmms=neg_hmms or None,
                )
                all_hits.extend(rows)
                for k, v in nc.items():
                    neg_counts[k] += v
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
                futs = []
                for t in tasks:
                    if (root / "STOP").exists():
                        raise PanelError(f"STOP file appeared in {root} between sub-steps (AGENTS §7.4)")
                    part = shards_out / f"{tag}_task_{t['task_id']:04d}.parquet"
                    futs.append(ex.submit(
                        _run_substep, t, shards_meta=shards_meta, pack_dir=pdir, protein_seq_files=pseq_files,
                        models=group_models, work_dir=work, partial_path=part, cpus=cpus,
                        skip_covered_paths=skip_paths, neg_hmms=neg_hmms or None,
                    ))
                for fut in futs:
                    rows, nc = fut.result()
                    all_hits.extend(rows)
                    for k, v in nc.items():
                        neg_counts[k] += v

    _run_group(reusable_models, covered_pseq_paths if covered_pseq_paths else None, "delta")
    _run_group(fresh_models, None, "full")

    # Sort deterministically & compute overlap/decoy competition across full current panel
    all_hits.sort(key=lambda r: (r["protein_uid"], r["model_key"], r["env_start"]))
    annotate_overlap_competition(all_hits)
    all_hits.sort(key=lambda r: (r["protein_uid"], r["model_key"], r["env_start"]))

    # Evaluate positive & negative controls
    hits_by_prot_model: dict[tuple[bytes, str], list[dict]] = defaultdict(list)
    target_hit_prots: set[bytes] = set()
    for r in all_hits:
        hits_by_prot_model[(r["protein_uid"], r["model_key"])].append(r)
        if not r["is_decoy"]:
            target_hit_prots.add(r["protein_uid"])

    controls_rows: list[dict] = []
    pos_passed_per_model: dict[str, int] = defaultdict(int)
    pos_total_per_model: dict[str, int] = defaultdict(int)
    all_pos_ok = True
    for mk, puid, spec_str in positives:
        if mk is not None:
            pos_total_per_model[mk] += 1
            matched = hits_by_prot_model.get((puid, mk), [])
            ok = bool(matched)
            if ok:
                pos_passed_per_model[mk] += 1
            best_b = max((r["bitscore"] for r in matched), default=0.0)
            best_e = min((r["i_evalue_zref"] for r in matched), default=1.0)
        else:
            ok = puid in target_hit_prots
            matched = [r for r in all_hits if r["protein_uid"] == puid and not r["is_decoy"]]
            best_b = max((r["bitscore"] for r in matched), default=0.0)
            best_e = min((r["i_evalue_zref"] for r in matched), default=1.0)
        if not ok:
            all_pos_ok = False
        controls_rows.append({
            "kind": "positive", "model_key": mk or "*", "spec": spec_str, "protein_uid": puid.hex(),
            "recovered": ok, "best_bitscore": round(best_b, 2), "best_i_evalue": best_e, "passed": ok,
        })

    neg_hits_per_model: dict[str, int] = defaultdict(int)
    for mk, puid, spec_str in neg_controls:
        check_models = [mk] if mk is not None else sorted(target_keys)
        for cmk in check_models:
            matched = hits_by_prot_model.get((puid, cmk), [])
            hit_above = bool(matched)
            if hit_above:
                neg_hits_per_model[cmk] += 1
            best_b = max((r["bitscore"] for r in matched), default=0.0)
            best_e = min((r["i_evalue_zref"] for r in matched), default=1.0)
            controls_rows.append({
                "kind": "negative_control", "model_key": cmk, "spec": spec_str, "protein_uid": puid.hex(),
                "recovered": hit_above, "best_bitscore": round(best_b, 2), "best_i_evalue": best_e,
                "passed": not hit_above,
            })

    neg_shuf_pass: dict[str, bool] = {}
    if negative_mode == "shuffled_seed":
        for m in target_models:
            nh = neg_counts.get(m.model_key, 0)
            exp = neg_exp.get(m.model_key, 0.0)
            ok = nh <= max(1, math.ceil(exp))
            neg_shuf_pass[m.model_key] = ok
            controls_rows.append({
                "kind": "shuffled_seed", "model_key": m.model_key, "spec": "shuffled_consensus",
                "protein_uid": "-", "recovered": nh > 0, "best_bitscore": float(nh),
                "best_i_evalue": exp, "passed": ok,
            })

    status = "ok"
    error = None
    if not all_pos_ok:
        status = "failed"
        error = "positive control not recovered"

    # Expand to ORF hits & call loci
    orf_hits = _expand_to_orf_hits(
        root, release, files, all_hits,
        gene_call_version=gene_call_version, tape_version=tv,
        stages_filter=target_stages if use_stage_subset else None,
        include_decoys=include_decoys, library=library, threads=threads, tmp=work,
    )

    # Call circular-aware loci on hit contigs
    hits_by_contig: dict[bytes, list[dict]] = defaultdict(list)
    for r in orf_hits:
        hits_by_contig[r["contig_uid"]].append(r)

    hit_contigs = set(hits_by_contig)
    tapes_map = nbr._tapes(files, tv, hit_contigs) if hit_contigs else {}
    igrs_map = _load_igrs_for_contigs(files, hit_contigs) if hit_contigs else {}
    prot_lens, prot_doms = _load_protein_meta_for_contigs(
        root, release, files, tapes_map, hits_by_contig, neighbors, gene_call_version, threads, work,
    )
    has_arr = arrays.has_arrays(files)
    arr_by_contig: dict[bytes, list[dict]] = {}
    if has_arr and hit_contigs:
        arr_by_contig, _ = arrays.read(files, hit_contigs, tv, tmp=work, threads=threads)

    loci_rows: list[dict] = []
    for cu in sorted(hit_contigs):
        tape = tapes_map.get(cu)
        if tape is None:
            continue
        loci_rows.extend(call_loci_for_contig(
            tape, hits_by_contig[cu],
            max_gap_genes=max_gap_genes, max_gap_bp=max_gap_bp, min_models=eff_min_models,
            distinct_orfs=distinct_orfs, require=set(require), forbid=set(forbid),
            all_overlaps=all_overlaps, neighbors=neighbors,
            igrs=igrs_map.get(cu, []), protein_lengths=prot_lens, protein_domains=prot_doms,
            crispr_arrays=arr_by_contig.get(cu, []),
        ))
    loci_rows.sort(key=lambda r: (r["contig_uid"], r["anchor_feature_uid"]))

    # Build models summary table
    hit_prots_by_m: dict[str, set[bytes]] = defaultdict(set)
    best_prots_by_m: dict[str, set[bytes]] = defaultdict(set)
    hit_orfs_by_m: dict[str, set[bytes]] = defaultdict(set)
    for r in all_hits:
        hit_prots_by_m[r["model_key"]].add(r["protein_uid"])
        if r["best_on_overlap"]:
            best_prots_by_m[r["model_key"]].add(r["protein_uid"])
    for r in orf_hits:
        hit_orfs_by_m[r["model_key"]].add(r["feature_uid"])

    models_rows: list[dict] = []
    for m in sorted(models, key=lambda x: x.model_key):
        models_rows.append({
            "model_key": m.model_key, "namespace": m.namespace, "name": m.name, "acc": m.acc,
            "model_len": m.model_len, "hmm_sha256": m.hmm_sha256, "is_decoy": m.is_decoy,
            "cutoff_mode": m.cutoff_mode, "seq_bits": float(m.seq_bits), "dom_bits": float(m.dom_bits),
            "i_evalue_max": m.i_evalue_max, "hmm_cov_min": m.hmm_cov_min, "min_dom_frac": m.min_dom_frac,
            "n_hit_proteins": len(hit_prots_by_m[m.model_key]),
            "n_best_proteins": len(best_prots_by_m[m.model_key]),
            "n_hit_orfs": len(hit_orfs_by_m[m.model_key]),
            "pos_controls_passed": pos_passed_per_model[m.model_key],
            "pos_controls_total": pos_total_per_model[m.model_key],
            "neg_control_hits": neg_hits_per_model[m.model_key],
            "neg_shuffled_hits": neg_counts.get(m.model_key, 0) if (negative_mode == "shuffled_seed" and not m.is_decoy) else None,
            "neg_shuffled_pass": neg_shuf_pass.get(m.model_key) if (negative_mode == "shuffled_seed" and not m.is_decoy) else None,
        })

    # Write outputs atomically into temp dir, then rename into out_dir
    tmp_out = out_dir.parent / f".{out_dir.name}.tmp-{os.getpid()}"
    shutil.rmtree(tmp_out, ignore_errors=True)
    tmp_out.mkdir(parents=True, exist_ok=True)
    try:
        _write_parquet_sorted(models_rows, MODELS_SCHEMA, tmp_out / "models.parquet")
        _write_tsv(models_rows, MODELS_SCHEMA.names, tmp_out / "models.tsv")
        _write_parquet_sorted(all_hits, PROTEIN_HITS_SCHEMA, tmp_out / "protein_hits.parquet")
        _write_parquet_sorted(orf_hits, ORF_HITS_SCHEMA, tmp_out / "orf_hits.parquet")
        _write_tsv(orf_hits, ORF_HITS_SCHEMA.names, tmp_out / "orf_hits.tsv")
        _write_parquet_sorted(loci_rows, LOCI_SUMMARY_SCHEMA, tmp_out / "loci_summary.parquet")
        _write_tsv(loci_rows, LOCI_SUMMARY_SCHEMA.names, tmp_out / "loci_summary.tsv")
        _write_tsv(controls_rows, ["kind", "model_key", "spec", "protein_uid", "recovered", "best_bitscore", "best_i_evalue", "passed"],
                   tmp_out / "controls.tsv")

        # Write GUIDE.md-compatible loci/*.tsv + unannotated_proteins.txt
        loci_dir = tmp_out / "loci"
        loci_dir.mkdir(exist_ok=True)
        unannotated_uids: set[bytes] = set()
        if loci_rows:
            anchor_puids = set()
            anchor_fus = {(r["contig_uid"], r["anchor_feature_uid"]) for r in loci_rows}
            for cu, afu in anchor_fus:
                t = tapes_map.get(cu)
                if t is not None:
                    for k, fu in enumerate(t["feature_uids"]):
                        if fu == afu and t["protein_uids"][k] is not None:
                            anchor_puids.add(t["protein_uids"][k])
                            break
            policy = nbr.Policy.parse(f"{neighbors}g", "any", None, "truncate", "auto")
            nbr_rows, _ = nbr.fetch(root, anchor_puids, policy, gene_call_version, tape_version=tv, files=files)
            nbr_rows = [r for r in nbr_rows if (r["contig_uid"], r["anchor_feature_uid"]) in anchor_fus]
            hit_bits = {}
            for h in all_hits:
                if not h["is_decoy"]:
                    hit_bits[h["protein_uid"]] = max(hit_bits.get(h["protein_uid"], 0.0), float(h["seq_bitscore"]))
            near_arr = arrays.near(nbr_rows, policy, files=files, tmp=work, threads=threads) if has_arr and nbr_rows else None
            if nbr_rows:
                for fname, lines in export.locus_tables(
                    root, nbr_rows, release=release, gene_call_version=gene_call_version,
                    hit_bits=hit_bits, manifest=files, library=library, threads=threads, tmp=work, arrays=near_arr,
                ):
                    export.write_lines(lines, loci_dir / fname)
                # Collect unscanned window proteins
                win_puids = {r["protein_uid"] for r in nbr_rows if r.get("protein_uid") is not None}
                with export.open_store(
                    root, release=release, gene_call_version=gene_call_version, manifest=files,
                    threads=threads, tmp=work, relations=("prot", "dom", "clus"),
                ) as st:
                    lib = st.library(library)
                    if lib is None:
                        unannotated_uids = win_puids
                    else:
                        st.con.register("_p_in", pa.table({"protein_uid": pa.array(sorted(win_puids), UID)}))
                        st.con.execute("CREATE TEMP TABLE _p AS SELECT protein_uid FROM _p_in")
                        export._annotate_proteins(st.con, lib, st.stages_of(lib))
                        for (pu,) in st.con.execute("SELECT protein_uid FROM _dstate WHERE domain_state = 'not scanned'").fetchall():
                            unannotated_uids.add(pu)

        export.write_lines([u.hex() for u in sorted(unannotated_uids)], tmp_out / "unannotated_proteins.txt")

        # Write coverage.json for incremental --carry-from
        man_dict = files.manifest or {}
        pseq_manifest = []
        for f in man_dict.get("files", []):
            if (
                f.get("table") == "protein_seq"
                and f.get("partitions", {}).get("gene_call_version") == gene_call_version
                and f.get("partitions", {}).get("stage") in target_stages
            ):
                pseq_manifest.append({
                    "path": f["path"],
                    "sha256": f["sha256"],
                    "stage": f["partitions"]["stage"],
                    "protein_prefix": f["partitions"].get("protein_prefix", ""),
                })
        pseq_manifest.sort(key=lambda x: x["path"])
        cov = {
            "schema_version": 1,
            "release": release,
            "stages": sorted(target_stages),
            "gene_call_version": gene_call_version,
            "tape_version": tv,
            "Z_PROT": Z_PROT,
            "DOM_Z": DOM_Z,
            "tool_versions": provenance.tool_versions("pyhmmer", "pyarrow"),
            "protein_seq_files": pseq_manifest,
        }
        (tmp_out / "coverage.json").write_text(json.dumps(cov, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        (tmp_out / "README.txt").write_text(_panel_readme(has_arr), encoding="utf-8")
        export.write_run_file(tmp_out, run_id.hex(), release)
        os.replace(tmp_out, out_dir)
    finally:
        shutil.rmtree(tmp_out, ignore_errors=True)

    # Write <runs>/run=<run_id>/run.parquet + panel_hits.parquet
    wall_s = time.monotonic() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_bytes = rss if sys.platform == "darwin" else rss * 1024
    run_row = {
        "run_id": run_id,
        "recipe_hash": rh,
        "recipe_yaml": cfg_text,
        "argv": list(user_command or []),
        "tool_versions_json": json.dumps(provenance.tool_versions("bio-box", "pyhmmer", "pyarrow"), sort_keys=True),
        "runtime_image_digest": os.environ.get("BIOBOX_IMAGE_DIGEST", "none"),
        "written_against": release,
        "ran_against": release,
        "pack_id": pid,
        "shard_sha256": [bytes.fromhex(s["sha256"]) for s in shards_meta if "sha256" in s],
        "z_prot": float(Z_PROT),
        "z_nt_mb": float(pack.Z["z_nt_mb"]),
        "query_hash": rh,
        "submitted_at": submitted,
        "wall_s": wall_s,
        "peak_rss_bytes": rss_bytes,
        "cpu_min_max_pid": None,
        "host": socket.gethostname(),
        "exit_code": 0 if status == "ok" else 1,
        "status": status,
        "partial_policy": "flag",
        "max_output_gb": 25.0,
        "target": "unique",
        "target_size": target_n_seq,
        "rounds": 1,
        "final_round": 1,
        "rounds_json": json.dumps([{"round": 1, "n_models": len(models), "n_protein_hits": len(all_hits),
                                    "n_orf_hits": len(orf_hits), "n_loci": len(loci_rows)}]),
        "stop_reason": "panel_complete",
        "controls_json": json.dumps(controls_rows, sort_keys=True),
        "warnings_json": json.dumps([]),
        "limits_json": json.dumps({"command": list(user_command or []), "stages": target_stages,
                                   "cpus": cpus, "jobs": jobs}, sort_keys=True),
        "error": error,
        "seed_json": json.dumps({"recipe": PANEL_RECIPE, "n_models": len(models),
                                 "model_keys": [m.model_key for m in models]}, sort_keys=True),
    }
    _write_parquet_sorted([run_row], table("search_run").schema, rdir / "run.parquet")
    shutil.copyfile(out_dir / "protein_hits.parquet", rdir / "panel_hits.parquet")
    shutil.rmtree(work, ignore_errors=True)
    return run_row, rdir

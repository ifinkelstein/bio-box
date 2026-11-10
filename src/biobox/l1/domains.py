"""Library domain scans with pyhmmer (operation O-SCAN).

Every distinct protein of a gene-calling version is searched with every model of
one HMM library, using the library's own bit-score cutoffs (Pfam GA, NCBIfam TC)
and a fixed search space (Z = domZ = 1e8, SEARCH_PLAN SD-04, DEC-04), so stored
E-values do not change as the corpus grows. Only included domains are kept.

Writes `<root>/l1/domains/library=<name>-<version>/stage=<S>/part-<k>.parquet` through
`biobox.io.write_parts`, cut at `max_file_bytes` of estimated size (card B-19), and the step's
ledger `ledger/stage=<S>/step=domains-<name>-<version>`.

Two ways to scan a stage, with byte-identical output (card B-24):

- `run`: the stage's proteins in one process.
- In shards, each in its own process under the 300 CPU-minute limit (AGENTS §7.2):
  `plan` cuts the stage's proteins, sorted by protein_uid, into shards of at most `shard_residues`
  residues each (cumulative `length`; a protein longer than that is a shard of its own) and writes
  the plan, `l1/domain_scan/library=<L>/stage=<S>/plan.json`; `scan_shard` scans shard k into a
  sorted scratch run, `<tmp>/biobox-domain-scan/<plan>/shard-<k>.parquet` (`<plan>`: the first 16
  hex of the plan's sha256), and writes its record `.../shard-<k>.json` beside the plan (the run's
  sha256 and counts); `merge` checks every run against its record and writes the runs, in shard
  order, as the files and ledger `run` writes. If scratch is cleared between the shards and the
  merge (a reboot clears /tmp), `rescan_shard` (`python -m biobox.l1.domains rescan`, card C-12)
  scans a shard again into scratch without writing its record again: the new run must have the
  sha256 the record holds, so the merge is of the same runs as before.

At a fixed Z a protein's hits do not depend on the other proteins searched (MR-12), the shards are
consecutive ranges of protein_uid, and `scan` sorts each one's rows by the table's sort key, so the
runs in shard order are the rows of `run` in its order, and `io.write_parts` cuts the same parts
from them. The plan and the shard records are the `domain_scan` artifact: small by nature, never in
a release manifest (a release is the same whether its scans were sharded or not), and what lets
`bio-box run-step` publish each step (a step must publish a file the registry claims). The runs stay
in scratch on local disk, never in the data root. The plan names its `protein_seq` files and the
sha256 of the HMM file, which every shard checks. Thread counts are the caller's (`cpus`).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import sys
import tarfile
import tempfile
from collections.abc import Iterator
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyhmmer

from biobox import ids, layout, ledger, provenance
from biobox import io as bio
from biobox.registry.schemas import ARTIFACTS, table

Z_PROT = 1e8
DOM_Z = 1e8
# estimate: zstd bytes per domain row; cuts parts only. Measured 53.8 B/row for Pfam and NCBIfam on the
# golden set (planning/measurements/2026-11-09_design_stress_test.md §3).
EST_ROW_BYTES = 54
PLAN_FORMAT = 1
SCAN_DIR = "biobox-domain-scan"  # under the scratch directory: the shards' runs, one directory per plan
_READ_ROWS = 16_384  # protein_seq rows per read while planning or reading a shard


def iter_hmms(path: Path):
    """Read HMMs from a .hmm, .hmm.gz, or a .tgz of .HMM files (NCBIfam), in file order."""
    path = Path(path)
    name = path.name
    if name.endswith((".tgz", ".tar.gz")):
        with tarfile.open(path, "r:gz") as tar:
            for member in sorted((m for m in tar.getmembers() if m.isfile() and m.name.upper().endswith(".HMM")),
                                 key=lambda m: m.name):
                data = tar.extractfile(member).read()
                with pyhmmer.plan7.HMMFile(io.BytesIO(data)) as hf:
                    yield from hf
        return
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rb") as fh, pyhmmer.plan7.HMMFile(fh) as hf:
        yield from hf


def _protein_files(root: Path, gene_call_version: str, stage: str) -> list[Path]:
    return sorted((Path(root) / "l0seq/protein_seq").glob(
        f"gene_call_version={gene_call_version}/stage={stage}/protein_prefix=*/part-*.parquet"))


def load_proteins(root: Path, gene_call_version: str, stage: str) -> dict[bytes, str]:
    out: dict[bytes, str] = {}
    for f in _protein_files(root, gene_call_version, stage):
        t = pq.read_table(f, columns=["protein_uid", "seq"])
        out.update(zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist()))
    return out


def _decode(x) -> str:
    return x.decode() if isinstance(x, bytes) else (x or "")


def scan(
    proteins: dict[bytes, str],
    hmm_path: Path,
    cutoff: str,
    cpus: int = 1,
) -> list[dict]:
    """Search `proteins` with every model in `hmm_path`; return included domains as rows."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    uids = sorted(proteins)
    block = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=u.hex().encode(), sequence=proteins[u].replace("*", "")) for u in uids
    ).digitize(alphabet)
    kwargs = {"cpus": cpus, "Z": Z_PROT, "domZ": DOM_Z}
    if cutoff != "none":
        kwargs["bit_cutoffs"] = cutoff
    rows = []
    for hits in pyhmmer.hmmsearch(iter_hmms(hmm_path), block, **kwargs):
        q = hits.query
        acc, name, mlen = _decode(q.accession) or _decode(q.name), _decode(q.name), q.M
        for hit in hits:
            if not hit.included:
                continue
            for dom in hit.domains:
                if not dom.included:
                    continue
                a = dom.alignment
                rows.append({
                    "protein_uid": bytes.fromhex(_decode(hit.name)),
                    "model_acc": acc,
                    "model_name": name,
                    "model_len": mlen,
                    "seq_bitscore": hit.score,
                    "seq_evalue_zref": hit.evalue,
                    "dom_bitscore": dom.score,
                    "dom_ievalue_zref": dom.i_evalue,
                    "env_start": dom.env_from - 1,
                    "env_end": dom.env_to,
                    "ali_start": a.target_from - 1,
                    "ali_end": a.target_to,
                    "hmm_from": a.hmm_from - 1,
                    "hmm_to": a.hmm_to,
                    "cutoff": cutoff,
                    "z": Z_PROT,
                    "dom_z": DOM_Z,
                })
    rows.sort(key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
    return rows


def _batches(rows: list[dict], spec) -> list[pa.RecordBatch]:
    return [pa.RecordBatch.from_pylist(rows[i : i + bio.ROW_GROUP_ROWS], schema=spec.schema)
            for i in range(0, len(rows), bio.ROW_GROUP_ROWS)]


def _names(stage: str, library: str, library_version: str) -> tuple[str, str]:
    """The domain files' directory and the step's name."""
    lib = f"{library}-{library_version}"
    return f"l1/domains/library={lib}/stage={stage}", f"domains-{lib}"


def _write(out: Path, stage: str, lib_dir: str, step: str, batches, n_in: int, hit: int,
           max_file_bytes: int) -> None:
    """The stage's domain files from sorted `batches`, and the step's ledger: `n_in` proteins in, `hit`
    with an included domain (the one writer of `run` and `merge`, so their files are the same)."""
    spec = table("l1_domains")
    bio.write_parts(batches, spec, out / lib_dir, EST_ROW_BYTES, max_file_bytes)
    lib = lib_dir.split("/")[2].removeprefix("library=").encode()
    tools = provenance.tool_versions("pyhmmer", "pyarrow")
    ledger.write(out, stage, step,
                 sorted([ledger.row(lib, "library", n_in, hit, tools=tools),
                         ledger.row(lib, "library", n_in, hit, n_rejected=n_in - hit,
                                    reason="no_hit", tools=tools)],
                        key=lambda r: (r["object_key"], r["reason_code"])))


def run(root: Path, stage: str, gene_call_version: str, library: str, library_version: str,
        hmm_path: Path, cutoff: str, cpus: int = 1, max_file_bytes: int = layout.MAX_FILE_BYTES,
        out: Path | None = None) -> dict:
    """Scan the stage's distinct proteins (read from `root`) with one library; write its domains and
    ledger under `out` (default `root`; card C-08). Refuses, before reading a protein, when either
    already holds a file in `root` or `out`."""
    out = Path(root) if out is None else Path(out)
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, out, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    proteins = load_proteins(root, gene_call_version, stage)
    rows = scan(proteins, hmm_path, cutoff, cpus=cpus)
    hit = len({r["protein_uid"] for r in rows})
    _write(out, stage, lib_dir, step, _batches(rows, table("l1_domains")), len(proteins), hit, max_file_bytes)
    return {"proteins": len(proteins), "domains": len(rows),
            "proteins_with_hit": hit,
            "models_with_hit": len({r["model_acc"] for r in rows})}


# --- the scan in shards (card B-24) -----------------------------------------------------------------


def scan_dir(stage: str, library: str, library_version: str) -> str:
    """The root-relative directory of the `domain_scan` artifact for one stage and library."""
    return f"{ARTIFACTS['domain_scan'].location}/library={library}-{library_version}/stage={stage}"


def _plan_rel(stage: str, library: str, library_version: str) -> str:
    return f"{scan_dir(stage, library, library_version)}/plan.json"


def _shard_rel(stage: str, library: str, library_version: str, k: int) -> str:
    return f"{scan_dir(stage, library, library_version)}/shard-{k:05d}.json"


def _json_bytes(doc: dict) -> bytes:
    return json.dumps(doc, sort_keys=True, indent=1).encode() + b"\n"


def _write_once(root: Path, out: Path, rel: str, data: bytes, what: str) -> None:
    for base in dict.fromkeys([Path(root).resolve(), Path(out).resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: {what} is written once; refusing to overwrite")
    path = Path(out) / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _find(root: Path, out: Path | None, rel: str) -> Path | None:
    """`rel` in the root, else in `out` (written earlier by a direct call into another output root)."""
    for base in [Path(root)] + ([] if out is None else [Path(out)]):
        if (base / rel).is_file():
            return base / rel
    return None


def _plan_bytes(root: Path, stage: str, library: str, library_version: str, out: Path | None = None) -> bytes:
    p = _find(root, out, _plan_rel(stage, library, library_version))
    if p is None:
        raise FileNotFoundError(f"{Path(root) / _plan_rel(stage, library, library_version)}: no shard plan for "
                                f"domains-{library}-{library_version} of stage {stage}; write it with domains.plan")
    return p.read_bytes()


def read_plan(root: Path, stage: str, library: str, library_version: str, out: Path | None = None) -> dict:
    """The stage's shard plan for one library (from the root, else from `out`)."""
    return json.loads(_plan_bytes(root, stage, library, library_version, out))


def run_path(tmp: Path | None, root: Path, stage: str, library: str, library_version: str, k: int,
             out: Path | None = None) -> Path:
    """Where `scan_shard` writes shard k's sorted run: under the scratch directory (`tmp`, else the
    system's), in a directory named by the plan's sha256, so runs of another plan never mix in."""
    digest = hashlib.sha256(_plan_bytes(root, stage, library, library_version, out)).hexdigest()[:16]
    return Path(tmp if tmp is not None else tempfile.gettempdir()) / SCAN_DIR / digest / f"shard-{k:05d}.parquet"


def _uid_lengths(files: list[Path]) -> Iterator[tuple[bytes, int]]:
    """(protein_uid, length) of the stage's `protein_seq` files in uid order: in bucket and part order the
    files are one sorted run (the table's sort key, bucketed by uid prefix); a uid that does not follow
    the one before is a RuntimeError."""
    last = None
    for f in files:
        for b in pq.ParquetFile(f).iter_batches(batch_size=_READ_ROWS, columns=["protein_uid", "length"]):
            for u, n in zip(b.column(0).to_pylist(), b.column(1).to_pylist()):
                if last is not None and u <= last:
                    raise RuntimeError(f"{f}: protein_seq rows are not in strictly increasing protein_uid order")
                last = u
                yield u, n


def _cut(rows: Iterator[tuple[bytes, int]], shard_residues: int) -> Iterator[dict]:
    """Consecutive shards of at most `shard_residues` residues; a protein that would not fit closes the
    shard before it (a protein above the limit is a shard of its own)."""
    cur = None
    for u, n in rows:
        if cur is not None and cur["residues"] + n > shard_residues:
            yield cur
            cur = None
        if cur is None:
            cur = {"first": u.hex(), "last": u.hex(), "proteins": 0, "residues": 0}
        cur["last"] = u.hex()
        cur["proteins"] += 1
        cur["residues"] += n
    if cur is not None:
        yield cur


def plan(root: Path, stage: str, gene_call_version: str, library: str, library_version: str, hmm_path: Path,
         cutoff: str, shard_residues: int, out: Path | None = None) -> dict:
    """Write the stage's shard plan for one library under `out` (default `root`): its proteins sorted by
    protein_uid, cut into shards of at most `shard_residues` residues. Refuses when the stage's domains of
    this library are written already, or a plan exists, in `root` or `out`. Returns the shard count and
    the totals."""
    root = Path(root)
    out = root if out is None else Path(out)
    if not isinstance(shard_residues, int) or shard_residues < 1:
        raise ValueError(f"shard_residues {shard_residues!r}: must be a positive integer")
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, out, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    rel = _plan_rel(stage, library, library_version)
    for base in dict.fromkeys([root.resolve(), out.resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: the shard plan is written once; refusing to overwrite")
    files = _protein_files(root, gene_call_version, stage)
    if not files:
        raise FileNotFoundError(f"stage {stage}: no protein_seq files for gene_call_version={gene_call_version}")
    shards = list(_cut(_uid_lengths(files), shard_residues))
    doc = {"format": PLAN_FORMAT, "stage": stage, "gene_call_version": gene_call_version, "library": library,
           "library_version": library_version, "cutoff": cutoff, "hmm_sha256": ids.sha256_file(Path(hmm_path)).hex(),
           "shard_residues": shard_residues,
           "files": [{"path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size} for f in files],
           "proteins": sum(s["proteins"] for s in shards), "residues": sum(s["residues"] for s in shards),
           "shards": shards}
    _write_once(root, out, rel, _json_bytes(doc), "the shard plan")
    return {"shards": len(shards), "proteins": doc["proteins"], "residues": doc["residues"]}


def _shard_proteins(root: Path, p: dict, s: dict) -> dict[bytes, str]:
    """The proteins of shard `s` (uids `first` to `last`) from the plan's files, reading only the row
    groups whose protein_uid range meets the shard's."""
    lo, hi = bytes.fromhex(s["first"]), bytes.fromhex(s["last"])
    out: dict[bytes, str] = {}
    for f in p["files"]:
        path = Path(root) / f["path"]
        if path.stat().st_size != f["bytes"]:
            raise ValueError(f"{path}: {path.stat().st_size} bytes, the shard plan lists {f['bytes']}")
        pf = pq.ParquetFile(path)
        col = pf.schema_arrow.get_field_index("protein_uid")
        for g in range(pf.metadata.num_row_groups):
            st = pf.metadata.row_group(g).column(col).statistics
            if st is not None and st.has_min_max and (bytes(st.max) < lo or bytes(st.min) > hi):
                continue
            t = pf.read_row_group(g, columns=["protein_uid", "seq"])
            for u, seq in zip(t.column(0).to_pylist(), t.column(1).to_pylist()):
                if lo <= u <= hi:
                    out[u] = seq
    return out


def _shard_inputs(root: Path, stage: str, library: str, library_version: str, hmm_path: Path, k: int,
                  out: Path | None = None) -> tuple[bytes, dict, dict, dict[bytes, str]]:
    """(the plan's bytes, the plan, its shard k, the shard's proteins), checked: k is one of the plan's shards,
    `hmm_path` is the plan's HMM file, and the proteins read are the plan's counts."""
    raw = _plan_bytes(root, stage, library, library_version, out)
    p = json.loads(raw)
    if not 0 <= k < len(p["shards"]):
        step = _names(stage, library, library_version)[1]
        raise ValueError(f"shard {k}: the plan of {step} for stage {stage} has {len(p['shards'])} shards")
    have = ids.sha256_file(Path(hmm_path)).hex()
    if have != p["hmm_sha256"]:
        raise ValueError(f"{hmm_path}: sha256 {have} is not the plan's HMM file ({p['hmm_sha256']})")
    s = p["shards"][k]
    proteins = _shard_proteins(root, p, s)
    if len(proteins) != s["proteins"] or sum(len(x) for x in proteins.values()) != s["residues"]:
        raise ValueError(f"shard {k}: {len(proteins)} proteins read, the plan lists {s['proteins']} "
                         f"({s['residues']} residues)")
    return raw, p, s, proteins


def scan_shard(root: Path, stage: str, library: str, library_version: str, hmm_path: Path, k: int, cpus: int,
               tmp: Path | None = None, out: Path | None = None) -> dict:
    """Scan shard `k` of the stage's plan with `cpus` threads: write its sorted run under scratch
    (`run_path`) and its record (the run's sha256 and counts) under `out` (default `root`). Refuses a
    shard whose record exists in `root` or `out`, an HMM file other than the plan's, and proteins that
    do not match the plan's counts."""
    root = Path(root)
    out = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out)
    if not isinstance(cpus, int) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give the thread count explicitly (at least 1)")
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, out, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    rel = _shard_rel(stage, library, library_version, k)
    for base in dict.fromkeys([root.resolve(), out.resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: shard {k} is scanned once; refusing to overwrite")
    raw, p, s, proteins = _shard_inputs(root, stage, library, library_version, hmm_path, k, out)
    rows = scan(proteins, hmm_path, p["cutoff"], cpus=cpus)
    path = run_path(tmp, root, stage, library, library_version, k, out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.unlink(missing_ok=True)  # a run left by an attempt that died before its record: scanned again
    spec = table("l1_domains")
    bio.write_table(_batches(rows, spec), spec, path)
    rec = {"format": PLAN_FORMAT, "plan_sha256": hashlib.sha256(raw).hexdigest(), "shard": k, "first": s["first"],
           "last": s["last"], "proteins": s["proteins"], "residues": s["residues"], "domains": len(rows),
           "proteins_with_hit": len({r["protein_uid"] for r in rows}), "run_sha256": ids.sha256_file(path).hex(),
           "run_bytes": path.stat().st_size}
    _write_once(root, out, rel, _json_bytes(rec), f"the record of shard {k}")
    return {"shard": k, "proteins": s["proteins"], "residues": s["residues"], "domains": len(rows),
            "proteins_with_hit": rec["proteins_with_hit"]}


# --- a shard whose scratch run was lost (card C-12) --------------------------------------------------------


def shard_record(root: Path, stage: str, library: str, library_version: str, k: int) -> dict | None:
    """The published record of shard k (in the root), or None when there is none."""
    p = _find(root, None, _shard_rel(stage, library, library_version, k))
    return None if p is None else json.loads(p.read_bytes())


def run_problem(root: Path, stage: str, library: str, library_version: str, k: int, tmp: Path | None) -> str | None:
    """Why the scratch run of shard k, whose record is published, cannot be merged as it is (`merge` would refuse
    it, or could not check it), or None: "is missing", "differs from the sha256 its record holds", or "cannot be
    checked: its record holds no run sha256" (a record written without one)."""
    rec = shard_record(root, stage, library, library_version, k)
    if rec is None:
        raise FileNotFoundError(f"shard {k} of domains-{library}-{library_version}, stage {stage}: no record")
    path = run_path(tmp, root, stage, library, library_version, k)
    if not path.is_file():
        return "is missing"
    if "run_sha256" not in rec:
        return "cannot be checked: its record holds no run sha256"
    if ids.sha256_file(path).hex() != rec["run_sha256"]:
        return "differs from the sha256 its record holds"
    return None


def rescan_shard(root: Path, stage: str, library: str, library_version: str, hmm_path: Path, k: int, cpus: int,
                 tmp: Path | None = None) -> dict:
    """Scan shard k again into scratch (`run_path`) for a shard whose record is published in `root` but whose run
    is gone or damaged, e.g. after a reboot cleared /tmp (card C-12). The record is not written again: the new run
    must reproduce it. Its sha256 must equal the record's `run_sha256`; a record that holds none (written without
    one) is checked by its domain and hit counts instead, and the result says `"verified": "counts"`. On a mismatch
    the new run is deleted, the run in scratch is left as it was, and a ValueError names the shard. Refuses, like
    `scan_shard`, when the stage's domains are written, the HMM file is not the plan's, or the proteins do not
    match the plan."""
    root = Path(root)
    ledger.check_scratch(tmp, root)
    if not isinstance(cpus, int) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give the thread count explicitly (at least 1)")
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, root, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    rec = shard_record(root, stage, library, library_version, k)
    if rec is None:
        raise FileNotFoundError(f"shard {k} of {step}, stage {stage}: no published record to scan it again for; "
                                "scan it with domains.scan_shard")
    raw, p, s, proteins = _shard_inputs(root, stage, library, library_version, hmm_path, k)
    if rec["plan_sha256"] != hashlib.sha256(raw).hexdigest() or rec["shard"] != k:
        raise ValueError(f"shard {k} of {step}, stage {stage}: its record is of another plan or shard")
    rows = scan(proteins, hmm_path, p["cutoff"], cpus=cpus)
    path = run_path(tmp, root, stage, library, library_version, k)
    path.parent.mkdir(parents=True, exist_ok=True)
    new = path.with_name(f".rescan-{os.getpid()}-{path.name}")
    spec = table("l1_domains")
    bio.write_table(_batches(rows, spec), spec, new)
    have, hits = ids.sha256_file(new).hex(), len({r["protein_uid"] for r in rows})
    want = rec.get("run_sha256")
    if want is not None:
        bad = None if have == want else f"its run has sha256 {have}, the record holds {want}"
    else:
        got, rec_counts = (len(rows), hits), (rec["domains"], rec["proteins_with_hit"])
        bad = None if got == rec_counts else (f"it gives {got[0]} domains in {got[1]} proteins, the record "
                                              f"{rec_counts[0]} in {rec_counts[1]}")
    if bad is not None:
        new.unlink()
        raise ValueError(f"shard {k} of {step}, stage {stage}: scanned again, {bad}; the published record is not "
                         f"the scan of these inputs, so the merge cannot run (the run in scratch is left as it was)")
    os.replace(new, path)
    return {"shard": k, "domains": len(rows), "proteins_with_hit": hits, "run_sha256": have,
            "verified": "sha256" if want is not None else "counts"}


def merge(root: Path, stage: str, library: str, library_version: str, tmp: Path | None = None,
          max_file_bytes: int = layout.MAX_FILE_BYTES, out: Path | None = None) -> dict:
    """Write the shards' runs, in shard order, as the stage's domain files and ledger under `out`
    (default `root`): the files `run` writes. Every shard needs its record (in `root` or `out`) and a run
    whose sha256 is the record's (for a record that holds no run sha256, whose row count is the record's
    domains; build-stage scans such a shard again first, `rescan_shard`); otherwise nothing is written."""
    root = Path(root)
    out_root = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out_root)
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, out_root, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    raw = _plan_bytes(root, stage, library, library_version, out)
    p, digest = json.loads(raw), hashlib.sha256(raw).hexdigest()
    missing = [k for k in range(len(p["shards"])) if _find(root, out, _shard_rel(stage, library, library_version, k)) is None]
    if missing:
        raise FileNotFoundError(f"{step} of stage {stage}: no record of shard {', '.join(map(str, missing))} "
                                f"of {len(p['shards'])}; scan it with domains.scan_shard")
    recs, runs = [], []
    for k in range(len(p["shards"])):
        rec = json.loads(_find(root, out, _shard_rel(stage, library, library_version, k)).read_bytes())
        if rec["plan_sha256"] != digest or rec["shard"] != k:
            raise ValueError(f"shard {k}: its record is of another plan or shard")
        path = run_path(tmp, root, stage, library, library_version, k, out)
        if not path.is_file():
            raise FileNotFoundError(f"{path}: the run of shard {k} is not in the scratch directory; its record "
                                    "is published, so scan the shard again into this scratch directory")
        if "run_sha256" not in rec:  # a record written without one: its row count is checked (card C-12)
            if pq.ParquetFile(path).metadata.num_rows != rec["domains"]:
                raise ValueError(f"{path}: row count differs from the domains in the record of shard {k}")
        elif ids.sha256_file(path).hex() != rec["run_sha256"]:
            raise ValueError(f"{path}: sha256 differs from the record of shard {k}")
        recs.append(rec)
        runs.append(path)
    models: set[str] = set()

    def batches() -> Iterator[pa.RecordBatch]:
        for path in runs:
            for b in pq.ParquetFile(path).iter_batches(batch_size=bio.ROW_GROUP_ROWS):
                models.update(pc.unique(b.column("model_acc")).to_pylist())
                yield b

    hit = sum(r["proteins_with_hit"] for r in recs)
    _write(out_root, stage, lib_dir, step, batches(), p["proteins"], hit, max_file_bytes)
    return {"proteins": p["proteins"], "domains": sum(r["domains"] for r in recs), "proteins_with_hit": hit,
            "models_with_hit": len(models), "shards": len(recs)}


def main(argv: list[str] | None = None) -> int:
    """`python -m biobox.l1.domains rescan --root R --stage S --library L --library-version V --hmm F --shard K
    --cpus N --tmp D`: `rescan_shard` in its own process, as build-stage runs it (card C-12). It writes nothing
    into the root, so it is not a run-step. Prints one JSON line; exit 2 with the reason on a refusal."""
    ap = argparse.ArgumentParser(prog="python -m biobox.l1.domains")
    ap.add_argument("action", choices=["rescan"])
    for name in ("--root", "--stage", "--library", "--library-version", "--hmm", "--tmp"):
        ap.add_argument(name, required=True)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--cpus", type=int, required=True)
    a = ap.parse_args(argv)
    try:
        res = rescan_shard(Path(a.root), a.stage, a.library, a.library_version, Path(a.hmm), a.shard, a.cpus,
                           tmp=Path(a.tmp))
    except (ValueError, FileNotFoundError, FileExistsError) as e:
        print(f"domains rescan: {e}", file=sys.stderr)
        return 2
    print(json.dumps({"action": a.action, "result": res, "stage": a.stage}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

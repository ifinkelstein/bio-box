"""Card C-64 (DEC-69, DEC-71): Targeted Logan SRA contig bait-pull (LGRT1) via logan50
for ART-RT and novel prokaryotic RTs.

Subcommands:
    build-bait       OUT_DIR [--both-hmm PATH] [--toro-hmm PATH] [--rvtall-hmm PATH] [--df-dir PATH]
    scan             BAIT_HMM OUT_DIR [--decoy-hmm PATH] [--cpus 6] [--batch-size 50000]
                     [--max-compressed-bytes BYTES] [--max-wall-s SEC] [--min-bits 32.0]
                     [--min-dom-frac 0.45] [--start-byte 0]
    resolve          HITS_JSONL OUT_DIR [--index-path PATH] [--max-members-per-cluster 50]
                     [--max-frames 200] [--max-runs 500]
    fetch-and-stage  RESOLVED_JSON OUT_META_DIR OUT_FASTA_DIR [--final-fasta-dir DIR]
                     [--fasta-name lgrt1_baited.fna.gz] [--min-len 2000]
                     [--root ROOT] [--root-stages ...] [--owner-inputs NAME=DIR ...]
                     [--skip-ena] [--threads 2]
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import duckdb
import pyhmmer

sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_inputs  # noqa: E402
from biobox import ids  # noqa: E402

Refused = make_inputs.Refused

UA = "curl/8.5.0"
MIN_STOR_FILE_BYTES = 64 * 1024 * 1024  # AGENTS.md §7.3

LOGAN_BASE = "https://logan-pub.s3.amazonaws.com"
LOGAN50_FASTA_KEY = "p/logan50_c1.0_nonhuman_complete.fa.zst"
LOGAN50_C2P_KEY = "p/logan50_c1.0_nonhuman_complete_cluster-to-protein.tsv.zst"
LOGAN50_C2P_INDEX_KEY = "p/logan50_c1.0_nonhuman_complete_cluster-to-protein.index.tsv"
LOGAN50_C2P_TOTAL_BYTES = 864_901_503_638
LOGAN50_FASTA_TOTAL_BYTES = 271_342_209_328

ENA_PORTAL_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

DEFAULT_BOTH_HMM = Path("/stor/scratch/Wilke/if2453/fragdef/retrons-new/task14/s5b/both.hmm")
DEFAULT_TORO_HMM = Path("/stor/scratch/Wilke/if2453/fragdef/retrons-new/task20/toro_types.hmm")
DEFAULT_RVTALL_HMM = Path("/stor/scratch/Wilke/if2453/fragdef/retrons-new/task20/rvtall.hmm")
DEFAULT_DF_PROFILES = Path("/stor/scratch/Wilke/if2453/fragdef/db/defense-finder/defense-finder-models/profiles")

ART_TARGET_MODELS = ("NTERM", "CTERM", "PROXORF", "DISTORF")
ART_DECOY_MODELS = (
    "NTERM_REV", "NTERM_PERM",
    "CTERM_REV", "CTERM_PERM",
    "PROXORF_REV", "PROXORF_PERM",
    "DISTORF_REV", "DISTORF_PERM",
)
TORO_BAIT_MODELS = ("clade11-A", "clade11-B")
RVTALL_BAIT_MODELS = ("RVT-G2L4", "RVT-CRISPR-G1", "RVT-CRISPR-G2", "RVT-CRISPR-G3")
DF_BAIT_FILES = ("DRT9__DRT9.hmm", "UG27__UG27.hmm")

# Per-model bitscore thresholds calibrated from W-39 (DEC-71)
MODEL_MIN_BITS: dict[str, float] = {
    "NTERM": 32.0,
    "CTERM": 32.0,
    "PROXORF": 32.0,
    "DISTORF": 32.0,
    "clade11-A": 40.0,
    "clade11-B": 40.0,
    "DRT9": 45.0,
    "DRT9__DRT9": 45.0,
    "UG28": 45.0,
    "UG27": 45.0,
    "UG27__UG27": 45.0,
    "RVT-G2L4": 45.0,
    "RVT-CRISPR-G1": 45.0,
    "RVT-CRISPR-G2": 45.0,
    "RVT-CRISPR-G3": 45.0,
}

_CL_TAG_RE = re.compile(r"(?:^|\s)CL:([^\s#;]+)")
_CIRCULAR_RE = re.compile(
    r"(?:(?:^|[\s|;,])(?:circular\s*=\s*(?:true|yes|y|1)|circ\s*=\s*(?:true|yes|y|1)|"
    r"is_circular\s*=\s*(?:true|yes|y|1)|topology\s*=\s*circular|circular_contig|dtr\b))",
    re.IGNORECASE,
)


def _as_str(val: str | bytes | None) -> str:
    if val is None:
        return ""
    return val.decode("utf-8", "replace") if isinstance(val, (bytes, bytearray)) else str(val)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def _sql_lit(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _sql_list(files: Iterable[Path]) -> str:
    return "[" + ", ".join(_sql_lit(f) for f in files) + "]"


def parse_protein_id(protein_id: str) -> tuple[str, str, int]:
    """Parse Logan protein ID `<run_acc>_<contig_num>_<gene_num>` -> (run_acc, contig_id, gene_idx)."""
    tok = protein_id.strip()
    if tok.startswith("CL:"):
        tok = tok[3:]
    parts = tok.rsplit("_", 2)
    if len(parts) != 3:
        raise ValueError(f"Unrecognized Logan protein_id format: {protein_id!r}")
    run_acc, contig_num, gene_str = parts
    return run_acc, f"{run_acc}_{contig_num}", int(gene_str)


def parse_logan_protein_header(header: str) -> dict:
    """Parse a Logan logan50 FASTA header (supporting both bare `<run>_<ctg>_<gene>` and `CL:<id>` tags)."""
    raw = header.lstrip(">").strip()
    tokens = raw.split()
    if not tokens:
        raise ValueError(f"Empty FASTA header: {header!r}")

    cl_match = _CL_TAG_RE.search(raw)
    first = tokens[0]
    if first.startswith("CL:"):
        cluster_id = first[3:]
        protein_id = tokens[1] if len(tokens) > 1 and not tokens[1].startswith("#") else cluster_id
    else:
        protein_id = first
        cluster_id = cl_match.group(1) if cl_match else protein_id

    run_acc, contig_id, gene_idx = parse_protein_id(protein_id)
    return {
        "protein_id": protein_id,
        "cluster_id": cluster_id,
        "run_acc": run_acc,
        "contig_id": contig_id,
        "gene_idx": gene_idx,
    }


def load_cluster_index(index_path: Path) -> list[tuple[str, int]]:
    """Load Logan `cluster-to-protein.index.tsv` -> sorted list of (first_cluster_id, byte_offset)."""
    entries: list[tuple[str, int]] = []
    with open(index_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            entries.append((parts[0], int(parts[1])))
    return entries


def locate_cluster_frames(
    index: list[tuple[str, int]],
    cluster_ids: Iterable[str],
    total_bytes: int,
) -> list[tuple[int, int, set[str]]]:
    """Map requested cluster_ids to sorted, deduplicated (start_byte, end_byte, {cluster_ids}) zstd frames."""
    if not index:
        return []
    keys = [k for k, _ in index]
    by_frame: dict[int, set[str]] = defaultdict(set)
    for cid in cluster_ids:
        pos = bisect.bisect_right(keys, cid) - 1
        if pos < 0:
            pos = 0
        by_frame[pos].add(cid)

    frames: list[tuple[int, int, set[str]]] = []
    for pos in sorted(by_frame):
        start_b = index[pos][1]
        end_b = (index[pos + 1][1] - 1) if (pos + 1 < len(index)) else (total_bytes - 1)
        frames.append((start_b, end_b, by_frame[pos]))
    return frames


def resolve_cluster_members_from_lines(
    lines: Iterable[str],
    wanted_clusters: set[str],
    *,
    max_members_per_cluster: int | None = None,
) -> dict[str, list[str]]:
    """Parse `<cluster_rep> <member_protein_id>` lines for wanted_clusters."""
    members: dict[str, set[str]] = defaultdict(set)
    for line in lines:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        cid, mid = parts[0], parts[1]
        if cid in wanted_clusters:
            members[cid].add(mid)

    out: dict[str, list[str]] = {}
    for cid in sorted(members):
        mset = members[cid]
        mset.add(cid)
        ordered = [cid] + sorted(x for x in mset if x != cid)
        if max_members_per_cluster is not None and len(ordered) > max_members_per_cluster:
            ordered = ordered[:max_members_per_cluster]
        out[cid] = ordered
    return out


def detect_exact_dtr(seq: bytes, min_k: int = 20, max_k: int = 1000) -> int:
    """Return length k of exact terminal direct repeat (20..1,000 bp) per DEC-65, or 0 if none/low-complexity."""
    n = len(seq)
    if n < 2 * min_k + 1:
        return 0
    upper = seq.upper()
    hi = min(max_k, n // 2 - 1)
    for k in range(hi, min_k - 1, -1):
        if upper[:k] == upper[n - k :]:
            rep = upper[:k]
            if len(set(rep)) <= 2:
                return 0
            return k
    return 0


def extract_baited_contigs_from_fasta(
    records: Iterable[tuple[str, bytes]],
    wanted_contigs: set[str],
    *,
    min_len: int = 2_000,
) -> tuple[list[dict], dict]:
    """Extract wanted contigs >= min_len bp with DEC-65 exact DTR trimming where circular-flagged."""
    kept: list[dict] = []
    stats = Counter()
    for hdr, raw_seq in records:
        cid = hdr.lstrip(">").split()[0]
        if cid not in wanted_contigs:
            continue
        seq = raw_seq.strip().upper()
        raw_len = len(seq)
        if raw_len < min_len:
            stats["short_dropped"] += 1
            stats["short_bases"] += raw_len
            continue
        is_circ_flag = bool(_CIRCULAR_RE.search(hdr))
        trimmed_seq = seq
        trimmed_k = 0
        is_circ = False
        if is_circ_flag:
            k = detect_exact_dtr(seq)
            if k > 0:
                trimmed_seq = seq[:-k]
                trimmed_k = k
                is_circ = True
                stats["dtr_trimmed"] += 1
            else:
                is_circ = True
        if len(trimmed_seq) < min_len:
            stats["short_dropped"] += 1
            stats["short_bases"] += raw_len
            continue
        run_acc = cid.rsplit("_", 1)[0]
        kept.append({
            "run_acc": run_acc,
            "contig_id": cid,
            "seq": trimmed_seq,
            "circular": is_circ,
            "trimmed_k": trimmed_k,
            "raw_len": raw_len,
        })
        stats["kept"] += 1
        stats["kept_bases"] += len(trimmed_seq)
    return kept, dict(stats)


def build_bait_hmm(
    out_dir: Path,
    *,
    both_hmm: Path = DEFAULT_BOTH_HMM,
    toro_hmm: Path = DEFAULT_TORO_HMM,
    rvtall_hmm: Path = DEFAULT_RVTALL_HMM,
    df_profiles_dir: Path = DEFAULT_DF_PROFILES,
) -> dict:
    """Extract the 12 target rare-RT bait HMMs and 8 inline ART decoys into deterministic files."""
    out_dir.mkdir(parents=True, exist_ok=True)
    target_path = out_dir / "bait_rt12.hmm"
    decoy_path = out_dir / "bait_art_decoys8.hmm"

    target_models: list[dict] = []
    decoy_models: list[dict] = []

    def _extract_named(hmm_path: Path, wanted: set[str], out_fh, record_list: list[dict], pack: str) -> set[str]:
        found: set[str] = set()
        with pyhmmer.plan7.HMMFile(hmm_path) as hf:
            for hmm in hf:
                name = _as_str(hmm.name)
                if name in wanted:
                    hmm.write(out_fh, binary=False)
                    found.add(name)
                    record_list.append({"name": name, "length": hmm.M, "pack": pack, "source": str(hmm_path)})
        return found

    with open(target_path, "wb") as tfh, open(decoy_path, "wb") as dfh:
        f_art = _extract_named(both_hmm, set(ART_TARGET_MODELS), tfh, target_models, "lab_art_rt")
        if f_art != set(ART_TARGET_MODELS):
            raise Refused(f"Missing ART target models: {set(ART_TARGET_MODELS) - f_art}")
        f_dec = _extract_named(both_hmm, set(ART_DECOY_MODELS), dfh, decoy_models, "lab_art_rt_decoy")
        if f_dec != set(ART_DECOY_MODELS):
            raise Refused(f"Missing ART decoy models: {set(ART_DECOY_MODELS) - f_dec}")

        f_toro = _extract_named(toro_hmm, set(TORO_BAIT_MODELS), tfh, target_models, "spire_retrons_2026")
        if f_toro != set(TORO_BAIT_MODELS):
            raise Refused(f"Missing Toro models: {set(TORO_BAIT_MODELS) - f_toro}")

        f_rvt = _extract_named(rvtall_hmm, set(RVTALL_BAIT_MODELS), tfh, target_models, "myrt")
        if f_rvt != set(RVTALL_BAIT_MODELS):
            raise Refused(f"Missing myRT models: {set(RVTALL_BAIT_MODELS) - f_rvt}")

        for fname in DF_BAIT_FILES:
            p = df_profiles_dir / fname
            if not p.is_file():
                raise Refused(f"Missing DefenseFinder HMM: {p}")
            with pyhmmer.plan7.HMMFile(p) as hf:
                for hmm in hf:
                    name = _as_str(hmm.name)
                    hmm.write(tfh, binary=False)
                    target_models.append({"name": name, "length": hmm.M, "pack": "defensefinder_rt", "source": str(p)})

    manifest = {
        "target_hmm": str(target_path),
        "target_hmm_sha256": _sha256_file(target_path),
        "target_count": len(target_models),
        "target_models": target_models,
        "decoy_hmm": str(decoy_path),
        "decoy_hmm_sha256": _sha256_file(decoy_path),
        "decoy_count": len(decoy_models),
        "decoy_models": decoy_models,
    }
    (out_dir / "bait_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def _load_earlier_contig_owners(
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path] | None,
    *,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict[bytes, str]:
    """Map contig_uid (16 bytes) -> earliest owning stage across built root stages and staged owner_inputs."""
    owners: dict[bytes, str] = {}
    if root is not None:
        stages = list(root_stages)
        if not stages:
            ci_dir = root / "l0seq" / "contig_index"
            if ci_dir.is_dir():
                stages = sorted(p.name.split("=", 1)[1] for p in ci_dir.glob("stage=*") if p.is_dir())
        if stages:
            con = duckdb.connect()
            con.execute(f"SET threads={int(threads)}")
            con.execute(f"SET memory_limit='{memory_limit}'")
            for st in stages:
                files = sorted((root / "l0seq" / "contig_index" / f"stage={st}").glob("*.parquet"))
                if not files:
                    continue
                for (uid,) in con.execute(f"SELECT contig_uid FROM read_parquet({_sql_list(files)})").fetchall():
                    b = bytes(uid)
                    if b not in owners:
                        owners[b] = st
            con.close()

    for st, d in (owner_inputs or {}).items():
        gjson = Path(d) / "genomes.json"
        if not gjson.is_file():
            continue
        with open(gjson, encoding="utf-8") as fh:
            for g in json.load(fh):
                for h in g.get("contig_uids", ()):
                    b = bytes.fromhex(h)
                    if b not in owners:
                        owners[b] = st
    return owners


def stage_lgrt1(
    *,
    baited_contigs: list[dict],
    out_meta_dir: Path,
    out_fasta_dir: Path,
    final_fasta_dir: Path | None = None,
    fasta_name: str = "lgrt1_baited.fna.gz",
    root: Path | None = None,
    root_stages: Iterable[str] = (),
    owner_inputs: dict[str, Path] | None = None,
    retrieved_at: str = "2026-11-24T00:00:00Z",
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict:
    """Deduplicate baited Logan contigs by contig_uid and write LGRT1 stage inputs conforming to s2_check_inputs.py."""
    out_meta_dir.mkdir(parents=True, exist_ok=True)
    out_fasta_dir.mkdir(parents=True, exist_ok=True)

    earlier_owners = _load_earlier_contig_owners(
        root, root_stages, owner_inputs, threads=threads, memory_limit=memory_limit
    )

    fasta_path = out_fasta_dir / fasta_name
    tmp_fasta = fasta_path.with_name(fasta_path.name + ".tmp")

    written_uids: set[bytes] = set()
    genomes_out: list[dict] = []
    records_rows: list[list] = []
    repeated_rows: list[list] = []
    trims_rows: list[list] = []
    uid_len: dict[bytes, int] = {}

    counts = Counter()
    biomes = Counter()
    models_counter = Counter()
    runs_kept: set[str] = set()
    total_bp = 0

    ordered = sorted(baited_contigs, key=lambda r: (r["run_acc"], r["contig_id"]))
    with gzip.open(tmp_fasta, "wb", compresslevel=6) as gz_out:
        for rec in ordered:
            seq: bytes = rec["seq"]
            uid = ids.contig_uid(seq)
            cid = rec["contig_id"]
            run_acc = rec["run_acc"]
            is_circ = bool(rec.get("circular", False))
            trimmed_k = int(rec.get("trimmed_k", 0))
            raw_len = int(rec.get("raw_len", len(seq) + trimmed_k))
            bait_models = sorted(set(rec.get("bait_models", [])))
            cluster_ids = sorted(set(rec.get("cluster_ids", [])))
            biome = str(rec.get("biome", "unknown"))

            if uid in written_uids:
                counts["repeated_contigs"] += 1
                repeated_rows.append([cid, run_acc, "duplicate_contig_uid", uid.hex()])
                continue

            written_uids.add(uid)
            uid_len[uid] = len(seq)
            runs_kept.add(run_acc)
            biomes[biome] += 1
            for bm in bait_models:
                models_counter[bm] += 1

            gz_out.write(b">" + uid.hex().encode("ascii") + b"\n")
            for i in range(0, len(seq), 80):
                gz_out.write(seq[i : i + 80] + b"\n")
            total_bp += len(seq)

            if is_circ:
                counts["circular_contigs"] += 1
            if trimmed_k > 0:
                counts["dtr_trimmed"] += 1
                trims_rows.append([run_acc, cid, raw_len, trimmed_k, len(seq)])

            genomes_out.append({
                "contig_uids": [uid.hex()],
                "genetic_code": 11,
                "kind": "prokaryote",
                "circular": [uid.hex()] if is_circ else [],
                "assembly_quality": "contig",
            })
            records_rows.append([
                cid, run_acc, uid.hex(), len(seq), int(is_circ),
                ",".join(bait_models), ",".join(cluster_ids), biome,
            ])

    # Enforce AGENTS.md §7.3: never persist < 64 MB files under /stor/
    fasta_bytes = tmp_fasta.stat().st_size
    if str(fasta_path.resolve()).startswith("/stor/") and fasta_bytes < MIN_STOR_FILE_BYTES:
        tmp_fasta.unlink(missing_ok=True)
        raise Refused(
            f"Refusing to write {fasta_bytes} B (< 64 MB) to /stor path {fasta_path} (AGENTS.md §7.3); "
            f"stage under /tmp/if2453/s2/inputs/LGRT1/ instead."
        )
    tmp_fasta.replace(fasta_path)
    fasta_sha256 = _sha256_file(fasta_path)

    shared_rows: list[list] = []
    shared_by_stage = Counter()
    shared_bases_by_stage = Counter()
    for uid in sorted(written_uids):
        owner = earlier_owners.get(uid)
        if owner is not None:
            ln = uid_len[uid]
            shared_rows.append([uid.hex(), owner, ln])
            shared_by_stage[owner] += 1
            shared_bases_by_stage[owner] += ln

    genomes_path = out_meta_dir / "genomes.json"
    genomes_path.write_text(json.dumps(genomes_out, separators=(",", ":")) + "\n", encoding="utf-8")

    obj_staging = [{
        "path": str(fasta_path),
        "source_name": "logan_contigs_baited",
        "source_release": "logan-v1-c1.0-lgrt1",
        "license": "INSDC / AWS Open Data",
        "redistributable": True,
        "retrieved_at": retrieved_at,
    }]
    obj_staging_path = out_meta_dir / "objects_staging.json"
    obj_staging_path.write_text(json.dumps(obj_staging, indent=2) + "\n", encoding="utf-8")

    final_path = (final_fasta_dir / fasta_name) if final_fasta_dir is not None else fasta_path
    obj_final = [{**obj_staging[0], "path": str(final_path)}]
    obj_final_path = out_meta_dir / "objects.json"
    obj_final_path.write_text(json.dumps(obj_final, indent=2) + "\n", encoding="utf-8")

    def _write_tsv(path: Path, header: list[str], rows: list[list]) -> str:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)
        return _sha256_file(path)

    rec_sha = _write_tsv(
        out_meta_dir / "records.tsv",
        ["contig_id", "run_acc", "contig_uid", "length", "circular", "bait_models", "cluster_ids", "biome"],
        records_rows,
    )
    rep_sha = _write_tsv(
        out_meta_dir / "repeated_records.tsv",
        ["contig_id", "run_acc", "reason", "contig_uid"],
        repeated_rows,
    )
    shr_sha = _write_tsv(
        out_meta_dir / "shared_records.tsv",
        ["contig_uid", "owner_stage", "length"],
        shared_rows,
    )
    trm_sha = _write_tsv(
        out_meta_dir / "trims.tsv",
        ["run_acc", "contig_id", "raw_len", "trimmed_k", "final_len"],
        trims_rows,
    )

    shared_contigs = len(shared_rows)
    shared_bp = sum(shared_bases_by_stage.values())

    report = {
        "stage": "LGRT1",
        "fasta_path": str(fasta_path),
        "fasta_bytes": fasta_bytes,
        "fasta_sha256": fasta_sha256,
        "genomes_json_sha256": _sha256_file(genomes_path),
        "objects_json_sha256": _sha256_file(obj_final_path),
        "objects_staging_json_sha256": _sha256_file(obj_staging_path),
        "records_tsv_sha256": rec_sha,
        "repeated_records_tsv_sha256": rep_sha,
        "shared_records_tsv_sha256": shr_sha,
        "trims_tsv_sha256": trm_sha,
        "distinct_contigs": len(written_uids),
        "distinct_runs": len(runs_kept),
        "total_bp": total_bp,
        "repeated_contigs": counts["repeated_contigs"],
        "circular_contigs": counts["circular_contigs"],
        "dtr_trimmed": counts["dtr_trimmed"],
        "shared_contigs": shared_contigs,
        "shared_bp": shared_bp,
        "net_new_contigs": len(written_uids) - shared_contigs,
        "net_new_bp": total_bp - shared_bp,
        "shared_by_stage": dict(shared_by_stage),
        "shared_bases_by_stage": dict(shared_bases_by_stage),
        "contigs_by_model": dict(models_counter),
        "top_biomes": biomes.most_common(30),
    }
    (out_meta_dir / "stage_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def _iter_wanted_fasta_bytes_stream(
    stream: io.BufferedIOBase,
    wanted_ids: set[str],
) -> Iterable[tuple[str, bytes]]:
    """Fast single-pass FASTA stream reader that only buffers sequences whose ID is in wanted_ids."""
    wanted_b = {w.encode("ascii", "replace") for w in wanted_ids}
    header: str | None = None
    chunks: list[bytes] = []
    keep = False
    for raw in stream:
        if raw.startswith(b">"):
            if keep and header is not None:
                yield header, b"".join(chunks)
                chunks = []
            tok = raw[1:].split(None, 1)[0]
            if tok in wanted_b:
                header = raw[1:].rstrip(b"\r\n").decode("utf-8", "replace")
                chunks = []
                keep = True
            else:
                header = None
                keep = False
        elif keep:
            chunks.append(raw.strip())
    if keep and header is not None:
        yield header, b"".join(chunks)


def _iter_fasta_bytes_stream(stream: io.BufferedIOBase) -> Iterable[tuple[str, bytes]]:
    header: str | None = None
    chunks: list[bytes] = []
    for raw in stream:
        if raw.startswith(b">"):
            if header is not None:
                yield header, b"".join(chunks)
                chunks = []
            header = raw[1:].rstrip(b"\r\n").decode("utf-8", "replace")
        else:
            if header is not None:
                chunks.append(raw.strip())
    if header is not None:
        yield header, b"".join(chunks)


def _scan_batch(
    batch: list[tuple[dict, bytes]],
    target_hmms: list[pyhmmer.plan7.HMM],
    decoy_hmms: list[pyhmmer.plan7.HMM],
    alphabet: pyhmmer.easel.Alphabet,
    *,
    cpus: int,
    default_min_bits: float,
    min_dom_frac: float,
    evalue: float,
) -> list[dict]:
    if not batch:
        return []
    seqs: list[pyhmmer.easel.DigitalSequence] = []
    meta_by_name: dict[str, tuple[dict, bytes]] = {}
    for parsed, aa in batch:
        clean = aa.rstrip(b"*")
        if len(clean) < 60:
            continue
        pid = parsed["protein_id"]
        try:
            tb = pyhmmer.easel.TextSequence(name=pid.encode("ascii", "replace"), sequence=clean.decode("ascii", "replace"))
            seqs.append(tb.digitize(alphabet))
            meta_by_name[pid] = (parsed, clean)
        except Exception:
            continue
    if not seqs:
        return []

    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)
    candidate_hits: list[dict] = []
    art_hit_names: set[str] = set()

    for th in pyhmmer.hmmer.hmmsearch(target_hmms, block, cpus=cpus, E=evalue, domE=evalue):
        model_name = _as_str(th.query.name)
        min_bits = MODEL_MIN_BITS.get(model_name, default_min_bits)
        is_art = model_name in ART_TARGET_MODELS
        for hit in th:
            if not hit.included:
                continue
            if hit.score < min_bits:
                continue
            max_dom = max((d.score for d in hit.domains), default=0.0)
            dom_frac = (max_dom / hit.score) if hit.score > 0 else 0.0
            if is_art and dom_frac < min_dom_frac:
                continue
            pid = _as_str(hit.name)
            parsed, _ = meta_by_name[pid]
            rec = {
                "protein_id": parsed["protein_id"],
                "cluster_id": parsed["cluster_id"],
                "run_acc": parsed["run_acc"],
                "contig_id": parsed["contig_id"],
                "gene_idx": parsed["gene_idx"],
                "model": model_name,
                "score": round(float(hit.score), 2),
                "max_dom_score": round(float(max_dom), 2),
                "dom_frac": round(float(dom_frac), 4),
                "evalue": float(hit.evalue),
                "vetoed_by_decoy": False,
            }
            candidate_hits.append(rec)
            if is_art:
                art_hit_names.add(pid)

    # Decoy veto only on the tiny set of ART hits
    if decoy_hmms and art_hit_names:
        art_seqs = []
        for pid in sorted(art_hit_names):
            _, clean = meta_by_name[pid]
            tb = pyhmmer.easel.TextSequence(name=pid.encode("ascii", "replace"), sequence=clean.decode("ascii", "replace"))
            art_seqs.append(tb.digitize(alphabet))
        art_block = pyhmmer.easel.DigitalSequenceBlock(alphabet, art_seqs)
        decoy_best: dict[str, float] = {}
        for th in pyhmmer.hmmer.hmmsearch(decoy_hmms, art_block, cpus=min(cpus, 2), E=1e-2):
            for hit in th:
                pid = _as_str(hit.name)
                if hit.score > decoy_best.get(pid, 0.0):
                    decoy_best[pid] = float(hit.score)
        for rec in candidate_hits:
            if rec["model"] in ART_TARGET_MODELS:
                d_sc = decoy_best.get(rec["protein_id"], 0.0)
                if d_sc >= rec["score"]:
                    rec["vetoed_by_decoy"] = True
                    rec["decoy_score"] = round(d_sc, 2)

    return [r for r in candidate_hits if not r["vetoed_by_decoy"]]


def cmd_build_bait(args: argparse.Namespace) -> None:
    manifest = build_bait_hmm(
        Path(args.out_dir),
        both_hmm=Path(args.both_hmm),
        toro_hmm=Path(args.toro_hmm),
        rvtall_hmm=Path(args.rvtall_hmm),
        df_profiles_dir=Path(args.df_dir),
    )
    print(json.dumps(manifest, indent=2))


def cmd_scan(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hits_path = out_dir / "hits.jsonl"
    ckpt_path = out_dir / "checkpoint.json"

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(args.bait_hmm) as hf:
        target_hmms = list(hf)
    decoy_hmms: list[pyhmmer.plan7.HMM] = []
    if args.decoy_hmm and Path(args.decoy_hmm).is_file():
        with pyhmmer.plan7.HMMFile(args.decoy_hmm) as hf:
            decoy_hmms = list(hf)

    url = f"{LOGAN_BASE}/{LOGAN50_FASTA_KEY}"
    headers = {"User-Agent": UA}
    if args.start_byte > 0:
        end_byte = (
            args.start_byte + args.max_compressed_bytes - 1
            if args.max_compressed_bytes > 0
            else LOGAN50_FASTA_TOTAL_BYTES - 1
        )
        headers["Range"] = f"bytes={args.start_byte}-{end_byte}"
    elif args.max_compressed_bytes > 0:
        headers["Range"] = f"bytes=0-{args.max_compressed_bytes - 1}"

    t0 = time.monotonic()
    c0 = time.process_time()
    comp_bytes = [0]
    stop_flag = threading.Event()

    proc = subprocess.Popen(
        ["zstd", "-dc", "--long=31"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1 << 20,
    )

    def _feeder() -> None:
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=120) as resp:
                while not stop_flag.is_set():
                    chunk = resp.read(1 << 18)  # 256 KiB
                    if not chunk:
                        break
                    comp_bytes[0] += len(chunk)
                    if proc.stdin is not None:
                        try:
                            proc.stdin.write(chunk)
                        except BrokenPipeError:
                            break
                    if args.max_compressed_bytes > 0 and comp_bytes[0] >= args.max_compressed_bytes:
                        break
                    if args.max_wall_s > 0 and (time.monotonic() - t0) >= args.max_wall_s:
                        break
        finally:
            if proc.stdin is not None:
                try:
                    proc.stdin.close()
                except Exception:
                    pass

    feeder = threading.Thread(target=_feeder, daemon=True)
    feeder.start()

    proteins_scanned = 0
    hits_found = 0
    hits_by_model = Counter()
    last_protein_id = ""
    batch: list[tuple[dict, bytes]] = []

    with open(hits_path, "w", encoding="utf-8") as hfh:
        assert proc.stdout is not None
        for hdr, seq in _iter_fasta_bytes_stream(proc.stdout):
            try:
                parsed = parse_logan_protein_header(hdr)
            except ValueError:
                continue
            last_protein_id = parsed["protein_id"]
            batch.append((parsed, seq))
            proteins_scanned += 1
            if len(batch) >= args.batch_size:
                hits = _scan_batch(
                    batch, target_hmms, decoy_hmms, alphabet,
                    cpus=args.cpus, default_min_bits=args.min_bits,
                    min_dom_frac=args.min_dom_frac, evalue=args.evalue,
                )
                for h in hits:
                    hfh.write(json.dumps(h) + "\n")
                    hits_found += 1
                    hits_by_model[h["model"]] += 1
                hfh.flush()
                batch.clear()

                wall_s = time.monotonic() - t0
                ckpt = {
                    "url": url,
                    "start_byte": args.start_byte,
                    "compressed_bytes_read": comp_bytes[0],
                    "total_archive_bytes": LOGAN50_FASTA_TOTAL_BYTES,
                    "pct_archive_scanned": round(100.0 * comp_bytes[0] / LOGAN50_FASTA_TOTAL_BYTES, 4),
                    "proteins_scanned": proteins_scanned,
                    "hits_found": hits_found,
                    "hits_by_model": dict(hits_by_model),
                    "last_protein_id": last_protein_id,
                    "wall_s": round(wall_s, 2),
                    "proteins_per_s": round(proteins_scanned / wall_s, 1) if wall_s > 0 else 0.0,
                    "comp_mb_per_s": round((comp_bytes[0] / 1e6) / wall_s, 3) if wall_s > 0 else 0.0,
                    "cpus": args.cpus,
                    "status": "running",
                }
                ckpt_path.write_text(json.dumps(ckpt, indent=2) + "\n", encoding="utf-8")
                if args.max_wall_s > 0 and wall_s >= args.max_wall_s:
                    stop_flag.set()
                    break

        if batch:
            hits = _scan_batch(
                batch, target_hmms, decoy_hmms, alphabet,
                cpus=args.cpus, default_min_bits=args.min_bits,
                min_dom_frac=args.min_dom_frac, evalue=args.evalue,
            )
            for h in hits:
                hfh.write(json.dumps(h) + "\n")
                hits_found += 1
                hits_by_model[h["model"]] += 1
            hfh.flush()

    stop_flag.set()
    try:
        proc.terminate()
    except Exception:
        pass
    proc.wait(timeout=15)
    feeder.join(timeout=15)

    wall_s = time.monotonic() - t0
    cpu_s = time.process_time() - c0
    final_ckpt = {
        "url": url,
        "start_byte": args.start_byte,
        "compressed_bytes_read": comp_bytes[0],
        "total_archive_bytes": LOGAN50_FASTA_TOTAL_BYTES,
        "pct_archive_scanned": round(100.0 * comp_bytes[0] / LOGAN50_FASTA_TOTAL_BYTES, 4),
        "proteins_scanned": proteins_scanned,
        "hits_found": hits_found,
        "hits_by_model": dict(hits_by_model),
        "last_protein_id": last_protein_id,
        "wall_s": round(wall_s, 2),
        "cpu_s": round(cpu_s, 2),
        "proteins_per_s": round(proteins_scanned / wall_s, 1) if wall_s > 0 else 0.0,
        "comp_mb_per_s": round((comp_bytes[0] / 1e6) / wall_s, 3) if wall_s > 0 else 0.0,
        "cpus": args.cpus,
        "status": "complete",
    }
    ckpt_path.write_text(json.dumps(final_ckpt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(final_ckpt, indent=2))


def cmd_resolve(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load hit clusters
    cluster_models: dict[str, set[str]] = defaultdict(set)
    cluster_best_score: dict[str, float] = {}
    with open(args.hits_jsonl, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            h = json.loads(line)
            cid = h["cluster_id"]
            cluster_models[cid].add(h["model"])
            sc = float(h["score"])
            if sc > cluster_best_score.get(cid, 0.0):
                cluster_best_score[cid] = sc

    idx_path = Path(args.index_path) if args.index_path else (out_dir / "cluster-to-protein.index.tsv")
    if not idx_path.is_file():
        url = f"{LOGAN_BASE}/{LOGAN50_C2P_INDEX_KEY}"
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as resp:
            idx_path.write_bytes(resp.read())

    index = load_cluster_index(idx_path)
    frames = locate_cluster_frames(index, cluster_models.keys(), LOGAN50_C2P_TOTAL_BYTES)

    # Optionally prioritize frames with the highest-scoring / most clusters if capped
    if args.max_frames > 0 and len(frames) > args.max_frames:
        frames.sort(key=lambda f: (-len(f[2]), -max(cluster_best_score[c] for c in f[2]), f[0]))
        frames = sorted(frames[: args.max_frames], key=lambda f: f[0])

    resolved_members: dict[str, list[str]] = {}
    deferred_rows: list[list] = []
    frames_fetched = 0
    c2p_url = f"{LOGAN_BASE}/{LOGAN50_C2P_KEY}"

    for start_b, end_b, wanted in frames:
        req = urllib.request.Request(c2p_url, headers={"User-Agent": UA, "Range": f"bytes={start_b}-{end_b}"})
        proc = subprocess.Popen(
            ["zstd", "-dc"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )

        def _feed_frame(r=req, p=proc) -> None:
            try:
                with urllib.request.urlopen(r, timeout=180) as resp:
                    while chunk := resp.read(1 << 18):
                        if p.stdin is not None:
                            p.stdin.write(chunk)
            finally:
                if p.stdin is not None:
                    try:
                        p.stdin.close()
                    except Exception:
                        pass

        t = threading.Thread(target=_feed_frame, daemon=True)
        t.start()
        assert proc.stdout is not None
        lines = (raw.decode("utf-8", "replace") for raw in proc.stdout)
        frame_res = resolve_cluster_members_from_lines(lines, wanted, max_members_per_cluster=None)
        proc.wait(timeout=60)
        t.join(timeout=60)
        frames_fetched += 1

        for cid, mlist in frame_res.items():
            if args.max_members_per_cluster > 0 and len(mlist) > args.max_members_per_cluster:
                kept = mlist[: args.max_members_per_cluster]
                for extra in mlist[args.max_members_per_cluster :]:
                    deferred_rows.append([cid, extra, "over_cluster_member_cap"])
                resolved_members[cid] = kept
            else:
                resolved_members[cid] = mlist

    # Ensure every hit cluster rep is included even if its frame wasn't fetched
    for cid in sorted(cluster_models):
        if cid not in resolved_members:
            resolved_members[cid] = [cid]

    # Group by SRA run accession -> contig_id -> {models, clusters}
    by_run: dict[str, dict[str, dict]] = defaultdict(lambda: defaultdict(lambda: {"models": set(), "clusters": set()}))
    for cid, mlist in resolved_members.items():
        mods = cluster_models[cid]
        for mid in mlist:
            try:
                run_acc, contig_id, _ = parse_protein_id(mid)
            except ValueError:
                continue
            by_run[run_acc][contig_id]["models"].update(mods)
            by_run[run_acc][contig_id]["clusters"].add(cid)

    # Order runs deterministically by (# hit contigs desc, run_acc asc)
    ordered_runs = sorted(by_run.keys(), key=lambda r: (-len(by_run[r]), r))
    if args.max_runs > 0 and len(ordered_runs) > args.max_runs:
        for r in ordered_runs[args.max_runs :]:
            for contig_id, info in sorted(by_run[r].items()):
                deferred_rows.append([",".join(sorted(info["clusters"])), contig_id, "over_lgrt1_run_cap"])
        ordered_runs = ordered_runs[: args.max_runs]

    runs_payload = {}
    for r in ordered_runs:
        runs_payload[r] = {
            ctg: {"models": sorted(info["models"]), "clusters": sorted(info["clusters"])}
            for ctg, info in sorted(by_run[r].items())
        }

    def_path = out_dir / "LGRT1_deferred_members.tsv.gz"
    with gzip.open(def_path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster_id", "member_or_contig_id", "reason"])
        w.writerows(deferred_rows)

    out_payload = {
        "hit_clusters": len(cluster_models),
        "frames_fetched": frames_fetched,
        "total_runs_resolved": len(by_run),
        "selected_runs": len(ordered_runs),
        "selected_contigs": sum(len(v) for v in runs_payload.values()),
        "deferred_members_count": len(deferred_rows),
        "deferred_members_path": str(def_path),
        "runs": runs_payload,
    }
    res_path = out_dir / "resolved_runs.json"
    res_path.write_text(json.dumps(out_payload, indent=2) + "\n", encoding="utf-8")
    summary = {k: v for k, v in out_payload.items() if k != "runs"}
    print(json.dumps(summary, indent=2))


def _fetch_ena_biomes(run_accs: list[str], *, batch_size: int = 40) -> dict[str, dict]:
    """Query ENA Portal API for read_run biome / library_source metadata."""
    out: dict[str, dict] = {}
    fields = "run_accession,sample_accession,study_accession,scientific_name,tax_id,library_source,library_strategy,instrument_platform,base_count"
    refusals = 0
    for i in range(0, len(run_accs), batch_size):
        chunk = run_accs[i : i + batch_size]
        query = " OR ".join(f'run_accession="{r}"' for r in chunk)
        params = urllib.parse.urlencode({
            "result": "read_run",
            "query": query,
            "fields": fields,
            "limit": "0",
            "format": "tsv",
        })
        url = f"{ENA_PORTAL_URL}?{params}"
        time.sleep(0.35)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as resp:
                text = resp.read().decode("utf-8", "replace")
            refusals = 0
            reader = csv.DictReader(io.StringIO(text), delimiter="\t")
            for row in reader:
                ra = (row.get("run_accession") or "").strip()
                if ra:
                    out[ra] = dict(row)
        except Exception:
            refusals += 1
            if refusals >= 5:
                break
    return out


def cmd_fetch_and_stage(args: argparse.Namespace) -> None:
    resolved = json.loads(Path(args.resolved_json).read_text(encoding="utf-8"))
    runs_map: dict[str, dict[str, dict]] = resolved["runs"]
    run_list = sorted(runs_map.keys())

    biomes_map: dict[str, dict] = {}
    if not args.skip_ena and run_list:
        biomes_map = _fetch_ena_biomes(run_list)

    baited_contigs: list[dict] = []
    fetch_stats = Counter()
    failed_runs: list[str] = []
    deferred_large_runs: list[list] = []
    out_meta = Path(args.out_meta_dir)
    out_meta.mkdir(parents=True, exist_ok=True)
    prog_path = out_meta / "fetch_progress.json"
    t0 = time.monotonic()
    max_run_bytes = int(args.max_run_bytes)

    def _fetch_one_run(run_acc: str) -> tuple[str, list[dict], dict, str | None, int]:
        wanted_map = runs_map[run_acc]
        wanted_set = set(wanted_map.keys())
        url = f"{LOGAN_BASE}/c/{run_acc}/{run_acc}.contigs.fa.zst"
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            resp = urllib.request.urlopen(req, timeout=120)
            cl_hdr = resp.headers.get("Content-Length")
            cl_bytes = int(cl_hdr) if cl_hdr and cl_hdr.isdigit() else 0
            if max_run_bytes > 0 and cl_bytes > max_run_bytes:
                resp.close()
                return run_acc, [], {}, "deferred_large_run", cl_bytes

            proc = subprocess.Popen(
                ["zstd", "-dc"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
            )
            err_box: list[Exception] = []

            def _feed(r=resp, p=proc, eb=err_box) -> None:
                try:
                    with r:
                        while chunk := r.read(1 << 18):
                            if p.stdin is not None:
                                p.stdin.write(chunk)
                except Exception as exc:
                    eb.append(exc)
                finally:
                    if p.stdin is not None:
                        try:
                            p.stdin.close()
                        except Exception:
                            pass

            t = threading.Thread(target=_feed, daemon=True)
            t.start()
            assert proc.stdout is not None
            kept, st = extract_baited_contigs_from_fasta(
                _iter_wanted_fasta_bytes_stream(proc.stdout, wanted_set),
                wanted_set,
                min_len=args.min_len,
            )
            proc.wait(timeout=60)
            t.join(timeout=60)
            if err_box:
                raise err_box[0]

            ena_row = biomes_map.get(run_acc, {})
            biome = (ena_row.get("scientific_name") or "unknown").strip() or "unknown"
            lib_src = (ena_row.get("library_source") or "").strip()
            if lib_src and lib_src.upper() == "METATRANSCRIPTOMIC":
                biome = f"{biome} [METATRANSCRIPTOMIC]"

            for rec in kept:
                info = wanted_map.get(rec["contig_id"], {})
                rec["bait_models"] = info.get("models", [])
                rec["cluster_ids"] = info.get("clusters", [])
                rec["biome"] = biome
            return run_acc, kept, st, None, cl_bytes
        except Exception as exc:
            return run_acc, [], {}, str(exc), 0

    workers = max(1, min(int(args.workers), 4))
    done_runs = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_fetch_one_run, r) for r in run_list]
        for fut in as_completed(futs):
            run_acc, kept, st, err, cl_bytes = fut.result()
            done_runs += 1
            if err is None:
                fetch_stats["runs_ok"] += 1
                for k, v in st.items():
                    fetch_stats[k] += v
                baited_contigs.extend(kept)
            elif err == "deferred_large_run":
                fetch_stats["runs_deferred_large"] += 1
                deferred_large_runs.append([run_acc, cl_bytes, len(runs_map[run_acc]), "over_max_run_bytes"])
            else:
                failed_runs.append(run_acc)
                fetch_stats["runs_failed"] += 1
            if done_runs % 50 == 0 or done_runs == len(run_list):
                prog_path.write_text(
                    json.dumps({
                        "done_runs": done_runs,
                        "total_runs": len(run_list),
                        "runs_ok": fetch_stats["runs_ok"],
                        "runs_deferred_large": fetch_stats["runs_deferred_large"],
                        "runs_failed": fetch_stats["runs_failed"],
                        "contigs_kept": fetch_stats["kept"],
                        "short_dropped": fetch_stats["short_dropped"],
                        "kept_bases": fetch_stats["kept_bases"],
                        "wall_s": round(time.monotonic() - t0, 1),
                    }, indent=2) + "\n",
                    encoding="utf-8",
                )

    if deferred_large_runs:
        deferred_large_runs.sort(key=lambda r: r[0])
        with open(out_meta / "LGRT1_deferred_large_runs.tsv", "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["run_acc", "compressed_bytes", "wanted_contigs", "reason"])
            w.writerows(deferred_large_runs)

    owner_inputs = {}
    for item in args.owner_inputs or []:
        k, v = item.split("=", 1)
        owner_inputs[k] = Path(v)
    if biomes_map:
        bpath = out_meta / "run_biomes.tsv"
        bcols = [
            "run_accession", "sample_accession", "study_accession",
            "scientific_name", "tax_id", "library_source",
            "library_strategy", "instrument_platform", "base_count",
        ]
        with open(bpath, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=bcols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in sorted(biomes_map):
                w.writerow(biomes_map[r])

    rep = stage_lgrt1(
        baited_contigs=baited_contigs,
        out_meta_dir=out_meta,
        out_fasta_dir=Path(args.out_fasta_dir),
        final_fasta_dir=Path(args.final_fasta_dir) if args.final_fasta_dir else None,
        fasta_name=args.fasta_name,
        root=Path(args.root) if args.root else None,
        root_stages=args.root_stages.split(",") if args.root_stages else (),
        owner_inputs=owner_inputs,
        retrieved_at=args.retrieved_at,
        threads=args.threads,
        memory_limit=args.memory_limit,
    )
    rep["fetch_stats"] = dict(fetch_stats)
    rep["failed_runs_count"] = len(failed_runs)
    (out_meta / "stage_report.json").write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(rep, indent=2))


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Card C-64: Targeted Logan SRA contig bait-pull (LGRT1)")
    sub = p.add_subparsers(dest="cmd", required=True)

    pb = sub.add_parser("build-bait")
    pb.add_argument("out_dir")
    pb.add_argument("--both-hmm", default=str(DEFAULT_BOTH_HMM))
    pb.add_argument("--toro-hmm", default=str(DEFAULT_TORO_HMM))
    pb.add_argument("--rvtall-hmm", default=str(DEFAULT_RVTALL_HMM))
    pb.add_argument("--df-dir", default=str(DEFAULT_DF_PROFILES))
    pb.set_defaults(func=cmd_build_bait)

    ps = sub.add_parser("scan")
    ps.add_argument("bait_hmm")
    ps.add_argument("out_dir")
    ps.add_argument("--decoy-hmm", default=None)
    ps.add_argument("--cpus", type=int, default=6)
    ps.add_argument("--batch-size", type=int, default=50_000)
    ps.add_argument("--max-compressed-bytes", type=int, default=0)
    ps.add_argument("--max-wall-s", type=float, default=0.0)
    ps.add_argument("--min-bits", type=float, default=32.0)
    ps.add_argument("--min-dom-frac", type=float, default=0.45)
    ps.add_argument("--evalue", type=float, default=1e-5)
    ps.add_argument("--start-byte", type=int, default=0)
    ps.set_defaults(func=cmd_scan)

    pr = sub.add_parser("resolve")
    pr.add_argument("hits_jsonl")
    pr.add_argument("out_dir")
    pr.add_argument("--index-path", default=None)
    pr.add_argument("--max-members-per-cluster", type=int, default=50)
    pr.add_argument("--max-frames", type=int, default=200)
    pr.add_argument("--max-runs", type=int, default=500)
    pr.set_defaults(func=cmd_resolve)

    pf = sub.add_parser("fetch-and-stage")
    pf.add_argument("resolved_json")
    pf.add_argument("out_meta_dir")
    pf.add_argument("out_fasta_dir")
    pf.add_argument("--final-fasta-dir", default=None)
    pf.add_argument("--fasta-name", default="lgrt1_baited.fna.gz")
    pf.add_argument("--min-len", type=int, default=2_000)
    pf.add_argument("--root", default=None)
    pf.add_argument("--root-stages", default="")
    pf.add_argument("--owner-inputs", action="append", default=[])
    pf.add_argument("--retrieved-at", default="2026-11-24T00:00:00Z")
    pf.add_argument("--skip-ena", action="store_true")
    pf.add_argument("--max-run-bytes", type=int, default=100_000_000)
    pf.add_argument("--workers", type=int, default=4)
    pf.add_argument("--threads", type=int, default=2)
    pf.add_argument("--memory-limit", default="4GB")
    pf.set_defaults(func=cmd_fetch_and_stage)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

"""Per-release search pack: build and verify (SEARCH_PLAN SD-06, SD-07; data_layer_and_tests.md §2B.2; C-12).

    $ATD_ROOT/search_pack/release=<R>/
      pack_manifest.json         shards, meta files, counts, sources, rules, builder versions
      unique/part-<k>.faa.zst    every exact-unique protein, named by protein_uid hex
      igr/part-<k>.fna.zst       IGRs >= 30 nt (+ strand), named by region_uid hex
      igr/igr_meta.parquet       search_region rows for the IGRs
      dark_orf/part-<k>.faa.zst  dark ORFs (biobox.search.orfs), named by region_uid hex
      dark_orf/dark_orf_meta.parquet
      models/<library>.hmm       the release's HMM libraries, re-written as plain HMMER3 text

S0 targets are `unique`, `igr`, `dark_orf` and `models` (ST0: no reduction). `spacers` waits
for CRISPR calls (O-ARRAY) and is listed under `missing_targets`; `nr50`, `nr90`, `segments`
and `igr90` come from S1, and the DIAMOND index of `unique` when DIAMOND enters the image.

A shard is cut once it holds `shard_text_bytes` of uncompressed FASTA (2 GiB: about 1 GiB of
protein and 0.6 GiB of nucleotide zstd), so every shard but the last of a target is over
256 MB. `pack_id` is the BLAKE3-128 of `pack_manifest.json`; the manifest holds no
timestamps, so the same inputs give the same pack_id (D-01).

The pack sits under the data root but outside the release manifest (catalog.EXCLUDED_TREES):
it has its own manifest and is verified by `verify` (C-12) before the release is signed.
"""

from __future__ import annotations

import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import ids, layout, provenance
from biobox.l0seq.store import SeqStore
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l1.igr import igr_version
from biobox.l3.geometry import junction
from biobox.l3.tapes import tape_version
from biobox.registry.schemas import table
from biobox.search import fastz, orfs

PACK_FORMAT = "biobox-pack-1"
SHARD_TEXT_BYTES = 2 * 2**30
IGR_MIN_NT = 30
Z = {"z_prot": domains.Z_PROT, "dom_z": domains.DOM_Z, "z_nt_mb": 1e5}  # SD-04, DEC-04
TARGET_KIND = {"unique": "protein", "igr": "igr", "dark_orf": "dark_orf"}
MISSING_AT_S0 = {"spacers": "no CRISPR spacer calls in this release (O-ARRAY not built)"}


def pack_dir(root: Path, release: str) -> Path:
    return Path(root) / f"search_pack/release={release}"


def _rel(root: Path, p: Path) -> str:
    return p.relative_to(root).as_posix()


def _source_hash(root: Path, paths: list[Path]) -> dict:
    h = blake3()
    for p in sorted(paths):
        h.update(f"{_rel(root, p)}:{ids.sha256_file(p).hex()}\n".encode())
    return {"files": len(paths), "blake3": h.hexdigest()}


def _protein_files(root: Path, gcv: str) -> list[Path]:
    return sorted((root / "l0seq/protein_seq").glob(f"gene_call_version={gcv}/stage=*/protein_prefix=*/part-*.parquet"))


def _occurrence_files(root: Path, gcv: str) -> list[Path]:
    return sorted((root / "l1/protein_occurrence").glob(f"gene_call_version={gcv}/stage=*/part-*.parquet"))


def _igr_files(root: Path, gcv: str) -> list[Path]:
    return sorted((root / "l0seq/igr").glob(f"igr_version={igr_version([gcv])}/stage=*/part-*.parquet"))


def _tape_files(root: Path, gcv: str) -> list[Path]:
    # the gene-call producer's tape version (card B-07)
    return sorted((root / "l3/tapes").glob(f"tape_version={tape_version([gcv])}/stage=*/part-*.parquet"))


def _contigs(root: Path, gcv: str) -> tuple[dict[bytes, int], set[bytes]]:
    lengths, circular = {}, set()
    for f in _tape_files(root, gcv):
        for r in pq.read_table(f, columns=["contig_uid", "contig_len", "is_circular"]).to_pylist():
            lengths[r["contig_uid"]] = r["contig_len"]
            if r["is_circular"]:
                circular.add(r["contig_uid"])
    return lengths, circular


def _slots(root: Path, gcv: str) -> dict[bytes, list[dict]]:
    by_contig: dict[bytes, list[dict]] = defaultdict(list)
    for f in _occurrence_files(root, gcv):
        for r in pq.read_table(f, columns=["contig_uid", "slot_ordinal", "start", "end", "strand"]).to_pylist():
            by_contig[r["contig_uid"]].append(r)
    for v in by_contig.values():
        v.sort(key=lambda r: r["slot_ordinal"])
    return by_contig


def _codes(root: Path, gcv: str) -> dict[bytes, int]:
    """Genetic code per contig: the commonest code among its slots (11 for a contig with none)."""
    counts: dict[bytes, Counter] = defaultdict(Counter)
    for f in sorted((root / "l1/slot_protein").glob(f"gene_call_version={gcv}/stage=*/part-*.parquet")):
        for r in pq.read_table(f, columns=["contig_uid", "genetic_code"]).to_pylist():
            counts[r["contig_uid"]][r["genetic_code"]] += 1
    return {cu: c.most_common(1)[0][0] for cu, c in counts.items()}


def _gc(seq: bytes) -> float:
    s = seq.upper()
    acgt = sum(s.count(b) for b in (b"A", b"C", b"G", b"T"))
    return (s.count(b"G") + s.count(b"C")) / acgt if acgt else 0.0


def _region(kind: str, cu: bytes, start: int, end: int, strand: int, n_residues: int, **kw) -> dict:
    row = {"region_uid": ids.region_uid(cu, start, end, strand), "target_kind": kind, "contig_uid": cu,
           "start": start, "end": end, "strand": strand, "n_residues": n_residues, "left_slot": None,
           "right_slot": None, "orientation": None, "gc": None, "dist_to_contig_end_bp": None,
           "partial": False, "at_contig_end": False, "gene_overlap_frac": None, "genetic_code": None}
    row.update(kw)
    return row


def _write_meta(rows: list[dict], path: Path) -> None:
    rows.sort(key=lambda r: r["region_uid"])
    pq.write_table(pa.Table.from_pylist(rows, schema=table("search_region").schema), path, compression="zstd",
                   row_group_size=122_880)


def _build_unique(root: Path, out: Path, gcv: str, shard_bytes: int) -> list[dict]:
    prots = {}
    for f in _protein_files(root, gcv):
        t = pq.read_table(f, columns=["protein_uid", "seq"])
        prots.update(zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist()))
    return fastz.write_records(out / "unique", ".faa.zst",
                               ((u.hex(), prots[u].replace("*", "")) for u in sorted(prots)), shard_bytes)


def _build_igr(root: Path, out: Path, gcv: str, store: SeqStore, lengths, circular, slots, shard_bytes):
    meta, records = [], []
    for f in _igr_files(root, gcv):
        for r in pq.read_table(f).to_pylist():
            cu, s, e = r["contig_uid"], r["start"], r["end"]
            if e - s < IGR_MIN_NT:
                continue
            L, circ = lengths[cu], cu in circular
            seq = store.get(cu, s, e, 1, circular=circ)
            ls = slots[cu][r["left_slot"]] if r["left_slot"] is not None else None
            rs = slots[cu][r["right_slot"]] if r["right_slot"] is not None else None
            # flanks of an IGR never nest, so only their strands matter
            orient = junction((0, 1, ls["strand"]), (2, 3, rs["strand"])) if ls and rs else None
            row = _region("igr", cu, s, e, 1, e - s, left_slot=r["left_slot"], right_slot=r["right_slot"],
                          orientation=orient, gc=_gc(seq), dist_to_contig_end_bp=None if circ else min(s, L - e),
                          partial=not circ and (s == 0 or e == L))
            meta.append(row)
            records.append((row["region_uid"].hex(), seq.decode("ascii")))
    records.sort()
    shards = fastz.write_records(out / "igr", ".fna.zst", records, shard_bytes)
    _write_meta(meta, out / "igr/igr_meta.parquet")
    return shards


def _build_dark_orf(root: Path, out: Path, gcv: str, store: SeqStore, lengths, circular, slots, shard_bytes):
    codes = _codes(root, gcv)
    meta, records = [], []
    for cu in sorted(lengths):
        circ = cu in circular
        seq = store.get(cu)
        genes = [(r["start"], r["end"], r["strand"]) for r in slots.get(cu, [])]
        code = codes.get(cu, 11)
        for o, frac, at_end in orfs.dark_orfs(seq, circ, genes, code):
            row = _region("dark_orf", cu, o.start, o.end, o.strand, len(o.aa), partial=o.partial,
                          at_contig_end=at_end, gene_overlap_frac=frac, genetic_code=code)
            meta.append(row)
            records.append((row["region_uid"].hex(), o.aa))
    records.sort()
    shards = fastz.write_records(out / "dark_orf", ".faa.zst", records, shard_bytes)
    _write_meta(meta, out / "dark_orf/dark_orf_meta.parquet")
    return shards


def _build_models(out: Path, models: dict[str, Path]) -> list[dict]:
    shards = []
    (out / "models").mkdir(parents=True, exist_ok=True)
    for lib in sorted(models):
        path = out / f"models/{lib}.hmm"
        n, m = 0, 0
        with open(path, "wb") as fh:
            for hmm in domains.iter_hmms(models[lib]):
                hmm.write(fh)
                n += 1
                m += hmm.M
        shards.append({"path": path, "n_seq": n, "n_residues": m})
    return shards


def build(root: Path, release: str, gene_call_version: str = GENE_CALL_VERSION,
          models: dict[str, Path] | None = None, shard_text_bytes: int = SHARD_TEXT_BYTES) -> dict:
    """Build the release's pack from the data root; returns the manifest (with `pack_id`)."""
    root = Path(root).resolve()
    out = pack_dir(root, release)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    gcv = gene_call_version
    lengths, circular = _contigs(root, gcv)
    slots = _slots(root, gcv)
    # a protein-only root (a synthetic world) has no sequence store: its nucleotide targets are empty
    store = SeqStore(root) if (root / "l0seq/seq_locator.parquet").exists() else None
    if store is None:
        lengths = {}
    try:
        built = {
            "unique": _build_unique(root, out, gcv, shard_text_bytes),
            "igr": _build_igr(root, out, gcv, store, lengths, circular, slots, shard_text_bytes),
            "dark_orf": _build_dark_orf(root, out, gcv, store, lengths, circular, slots, shard_text_bytes),
            "models": _build_models(out, models or {}),
        }
    finally:
        if store is not None:
            store.close()
    seq_files = sorted((root / "l0seq/seq_shards").glob("stage=*/shard-*.fa.gz"))
    sources = {
        "unique": {"protein_seq": _source_hash(root, _protein_files(root, gcv))},
        "igr": {"igr": _source_hash(root, _igr_files(root, gcv)), "seq_shards": _source_hash(root, seq_files)},
        "dark_orf": {"protein_occurrence": _source_hash(root, _occurrence_files(root, gcv)),
                     "contig_tapes": _source_hash(root, _tape_files(root, gcv)),
                     "seq_shards": _source_hash(root, seq_files)},
        "models": {lib: {"files": 1, "sha256": ids.sha256_file(p).hex()} for lib, p in sorted((models or {}).items())},
    }
    manifest = {
        "pack_format": PACK_FORMAT,
        "release_id": release,
        "gene_call_version": gcv,
        "igr_version": igr_version([gcv]),
        "z": Z,
        "rules": {"igr_min_nt": IGR_MIN_NT, "dark_orf_min_aa": orfs.MIN_AA, "dark_orf_end_bp": orfs.END_BP,
                  "dark_orf_max_gene_overlap": orfs.MAX_GENE_OVERLAP, "shard_text_bytes": shard_text_bytes},
        "builder": provenance.tool_versions("bio-box", "pyhmmer", "pyarrow"),
        "targets": {t: {"kind": TARGET_KIND.get(t, t), "n_seq": sum(s["n_seq"] for s in sh),
                        "n_residues": sum(s["n_residues"] for s in sh), "sources": sources[t]}
                    for t, sh in built.items()},
        "missing_targets": MISSING_AT_S0,
        "shards": [{"target": t, "path": _rel(out, s["path"]), "sha256": ids.sha256_file(s["path"]).hex(),
                    "bytes": s["path"].stat().st_size, "n_seq": s["n_seq"], "n_residues": s["n_residues"]}
                   for t, sh in built.items() for s in sh],
        "meta": [{"target": t, "path": f"{t}/{t}_meta.parquet",
                  "sha256": ids.sha256_file(out / f"{t}/{t}_meta.parquet").hex(),
                  "bytes": (out / f"{t}/{t}_meta.parquet").stat().st_size} for t in ("igr", "dark_orf")],
    }
    data = json.dumps(manifest, sort_keys=True, indent=1).encode()
    (out / "pack_manifest.json").write_bytes(data)
    return {**manifest, "pack_id": pack_id_of(data).hex()}


def pack_id_of(manifest_bytes: bytes) -> bytes:
    return blake3(manifest_bytes).digest(length=16)


def read_manifest(root: Path, release: str) -> tuple[dict, bytes]:
    data = (pack_dir(root, release) / "pack_manifest.json").read_bytes()
    return json.loads(data), pack_id_of(data)


def _count_source(root: Path, manifest: dict, target: str, out: Path) -> int | None:
    gcv = manifest["gene_call_version"]
    if target == "unique":
        return len({u for f in _protein_files(root, gcv)
                    for u in pq.read_table(f, columns=["protein_uid"]).column("protein_uid").to_pylist()})
    if target == "igr":
        return sum(1 for f in _igr_files(root, gcv)
                   for r in pq.read_table(f, columns=["start", "end"]).to_pylist()
                   if r["end"] - r["start"] >= manifest["rules"]["igr_min_nt"])
    if target == "dark_orf":
        return pq.read_metadata(out / "dark_orf/dark_orf_meta.parquet").num_rows
    return None


def verify(root: Path, release: str, pack_id: bytes | None = None,
           min_shard_bytes: int = layout.TARGET_FILE_BYTES) -> list[str]:
    """C-12: problems with the release's pack; empty means it verifies."""
    root = Path(root).resolve()
    out = pack_dir(root, release)
    if not (out / "pack_manifest.json").exists():
        return [f"{_rel(root, out)}: no pack_manifest.json"]
    manifest, pid = read_manifest(root, release)
    problems = []
    if pack_id is not None and pack_id != pid:
        problems.append(f"pack_id {pack_id.hex()} differs from the manifest hash {pid.hex()}")
    listed = {"pack_manifest.json"}
    by_target: dict[str, list[dict]] = defaultdict(list)
    for entry in manifest["shards"] + manifest["meta"]:
        listed.add(entry["path"])
        p = out / entry["path"]
        if not p.exists():
            problems.append(f"{entry['path']}: listed but missing")
            continue
        if p.stat().st_size != entry["bytes"]:
            problems.append(f"{entry['path']}: {p.stat().st_size} bytes, manifest says {entry['bytes']}")
        if ids.sha256_file(p).hex() != entry["sha256"]:
            problems.append(f"{entry['path']}: sha256 differs from the manifest")
    for entry in manifest["shards"]:
        by_target[entry["target"]].append(entry)
    for p in sorted(x for x in out.rglob("*") if x.is_file()):
        if _rel(out, p) not in listed:
            problems.append(f"{_rel(out, p)}: on disk but not in the manifest")
    for target, shards in sorted(by_target.items()):
        info = manifest["targets"][target]
        if sum(s["n_seq"] for s in shards) != info["n_seq"]:
            problems.append(f"{target}: shard counts do not sum to the target's n_seq")
        # models are one file per library (small by nature, SEARCH_PLAN §3.2), not size-cut shards
        for s in (sorted(shards, key=lambda s: s["path"])[:-1] if target != "models" else []):
            if s["bytes"] < min_shard_bytes:
                problems.append(f"{s['path']}: {s['bytes']} bytes, under {min_shard_bytes} and not the last shard")
        if target == "models":
            for s in shards:
                p = out / s["path"]
                if p.exists() and sum(1 for _ in domains.iter_hmms(p)) != s["n_seq"]:
                    problems.append(f"{s['path']}: model count differs from the manifest")
            continue
        n_records = 0
        for s in shards:
            p = out / s["path"]
            if p.exists():
                n = sum(1 for _ in fastz.read_records(p))
                n_records += n
                if n != s["n_seq"]:
                    problems.append(f"{s['path']}: {n} records, manifest says {s['n_seq']}")
        expect = _count_source(root, manifest, target, out)
        if expect is not None and expect != info["n_seq"]:
            problems.append(f"{target}: {info['n_seq']} sequences, its source has {expect}")
    return problems

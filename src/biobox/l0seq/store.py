"""Write and read the L0-seq store (operations O-INGEST, O-SLICE, O-EXTRACT).

Layout under `<root>/` (planning/data_layer_and_tests.md §2A.1, §3.1):

    l0/l0_manifest.parquet
    l0seq/store_params.json
    l0seq/seq_chunks/stage=<S>/digest_prefix=<xx>/part-0.parquet
    l0seq/contig_occurrence/stage=<S>/source=<name>/part-0.parquet
    l0seq/contig_index/stage=<S>/part-0.parquet
    l0seq/seq_locator.parquet

A stage never rewrites an earlier stage's files: a sequence already stored adds
only occurrence rows. `seq_locator`, `l0_manifest` and `store_params.json` are
small metadata files rewritten per stage.
"""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox import fasta, ids
from biobox.registry.schemas import table

CHUNK_SIZE = 65_536
ROW_GROUP_CHUNKS = 16  # 1 MiB of sequence per row group; K-13 / LY-06 (2026-10-30): 1,024 decoded 36 MB per 1 kb slice

_COMP = bytes.maketrans(b"ACGTRYKMBVDHNacgtrykmbvdhn", b"TGCAYRMKVBHDNtgcayrmkvbhdn")


def revcomp(seq: bytes) -> bytes:
    return seq.translate(_COMP)[::-1]


@dataclass(frozen=True)
class SourceObject:
    path: Path
    source_name: str
    source_release: str
    license: str = "unknown"
    redistributable: bool = False
    precedence: int = 100  # lower wins when choosing a contig's default variant
    retrieved_at: datetime = field(default_factory=lambda: datetime(1970, 1, 1, tzinfo=timezone.utc))


# Min/max statistics only on key columns: on the 64 KB `seq` strings they bloat the footer
# (4.8 MB for one golden-set file with 4 KB chunks, K-13) and never prune anything.
SEQ_CHUNK_STATS = ["raw_digest", "chunk_ix", "contig_len"]


def _write(tbl: pa.Table, path: Path, row_group_size: int, statistics: list[str] | bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(tbl, path, compression="zstd", row_group_size=row_group_size, write_statistics=statistics)


def _existing_locator(root: Path) -> pa.Table | None:
    p = root / "l0seq" / "seq_locator.parquet"
    return pq.read_table(p) if p.exists() else None


def _existing_index(root: Path) -> pa.Table | None:
    d = root / "l0seq" / "contig_index"
    if not d.exists():
        return None
    return pa.concat_tables(pq.read_table(f) for f in sorted(d.glob("stage=*/part-0.parquet")))


def ingest(
    objects: list[SourceObject],
    root: Path,
    stage: str,
    chunk_size: int = CHUNK_SIZE,
    row_group_chunks: int = ROW_GROUP_CHUNKS,
) -> dict:
    """Ingest FASTA source objects for one stage. Returns counts for the ledger."""
    root = Path(root)
    params_path = root / "l0seq" / "store_params.json"
    if params_path.exists():
        params = json.loads(params_path.read_text())
        if params["chunk_size"] != chunk_size:
            raise ValueError("chunk_size differs from the existing store")
    else:
        params_path.parent.mkdir(parents=True, exist_ok=True)
        params_path.write_text(json.dumps({"chunk_size": chunk_size}, sort_keys=True))

    locator = _existing_locator(root)
    known = set(locator.column("raw_digest").to_pylist()) if locator is not None else set()
    old_index = _existing_index(root)
    old_default = (
        {u for u, d in zip(old_index.column("contig_uid").to_pylist(), old_index.column("is_default_variant").to_pylist()) if d}
        if old_index is not None
        else set()
    )
    old_pairs = (
        set(zip(old_index.column("contig_uid").to_pylist(), old_index.column("raw_digest").to_pylist()))
        if old_index is not None
        else set()
    )

    manifest_rows, occ_rows = [], []
    new_seqs: dict[bytes, bytes] = {}
    variants: dict[tuple[bytes, bytes], tuple] = {}  # (contig_uid, raw_digest) -> best precedence key
    for obj in objects:
        obj_sha = ids.sha256_file(obj.path)
        src_uid = ids.source_uid(obj.source_name, obj.source_release)
        manifest_rows.append({
            "source_uid": src_uid,
            "source_name": obj.source_name,
            "source_release": obj.source_release,
            "object_path": Path(obj.path).name,
            "object_sha256": obj_sha,
            "bytes": Path(obj.path).stat().st_size,
            "license": obj.license,
            "redistributable": obj.redistributable,
            "retrieved_at": obj.retrieved_at,
        })
        for ordinal, rec in enumerate(fasta.read(Path(obj.path))):
            rd = ids.raw_digest(rec.seq)
            cu = ids.contig_uid(rec.seq)
            occ_rows.append({
                "source_uid": src_uid,
                "source_release": obj.source_release,
                "l0_object_sha256": obj_sha,
                "record_ordinal": ordinal,
                "accession": rec.accession,
                "original_header": rec.header,
                "raw_digest": rd,
                "contig_uid": cu,
                "line_width": rec.line_width,
                "_source": obj.source_name,
            })
            if rd not in known:
                new_seqs[rd] = rec.seq
            key = (obj.precedence, obj_sha, ordinal)
            if (cu, rd) not in variants or key < variants[(cu, rd)]:
                variants[(cu, rd)] = key

    # seq_chunks: one file per digest prefix, sorted by (raw_digest, chunk_ix)
    by_prefix: dict[str, list[bytes]] = defaultdict(list)
    for rd in new_seqs:
        by_prefix[rd[:1].hex()].append(rd)
    spec = table("seq_chunks")
    new_loc = []
    for prefix in sorted(by_prefix):
        rows = {"raw_digest": [], "chunk_ix": [], "seq": [], "contig_len": []}
        firsts = []
        for rd in sorted(by_prefix[prefix]):
            seq = new_seqs[rd]
            firsts.append((rd, len(rows["raw_digest"])))
            n = max(1, -(-len(seq) // chunk_size))
            for ix in range(n):
                rows["raw_digest"].append(rd)
                rows["chunk_ix"].append(ix)
                rows["seq"].append(seq[ix * chunk_size : (ix + 1) * chunk_size].decode("ascii"))
                rows["contig_len"].append(len(seq))
        rel = f"l0seq/seq_chunks/stage={stage}/digest_prefix={prefix}/part-0.parquet"
        _write(pa.table(rows, schema=spec.schema), root / rel, row_group_chunks, SEQ_CHUNK_STATS)
        for rd, first in firsts:
            new_loc.append({"raw_digest": rd, "file": rel, "row_group": first // row_group_chunks, "first_chunk_row": first})

    # contig_occurrence: partitioned by source, sorted by (object sha, ordinal)
    occ_spec = table("contig_occurrence")
    by_source: dict[str, list[dict]] = defaultdict(list)
    for r in occ_rows:
        by_source[r["_source"]].append({k: v for k, v in r.items() if k != "_source"})
    for src in sorted(by_source):
        rows = sorted(by_source[src], key=lambda r: (r["l0_object_sha256"], r["record_ordinal"]))
        rel = f"l0seq/contig_occurrence/stage={stage}/source={src}/part-0.parquet"
        _write(pa.Table.from_pylist(rows, schema=occ_spec.schema), root / rel, 122_880)

    # contig_index: new (contig_uid, raw_digest) pairs; one default per contig across stages
    idx_spec = table("contig_index")
    best: dict[bytes, tuple] = {}
    for (cu, rd), key in variants.items():
        if cu not in best or key < best[cu][0]:
            best[cu] = (key, rd)
    lengths = {rd: len(s) for rd, s in new_seqs.items()}
    if locator is not None:
        lengths.update(_lengths_from_store(root, locator, set(rd for _, rd in variants) - set(lengths)))
    idx_rows = [
        {
            "contig_uid": cu,
            "raw_digest": rd,
            "contig_len": lengths[rd],
            "is_default_variant": cu not in old_default and best[cu][1] == rd,
        }
        for (cu, rd) in sorted(variants)
        if (cu, rd) not in old_pairs
    ]
    if idx_rows:
        _write(
            pa.Table.from_pylist(idx_rows, schema=idx_spec.schema),
            root / f"l0seq/contig_index/stage={stage}/part-0.parquet",
            122_880,
        )

    # seq_locator and l0_manifest: small metadata, rewritten with all stages
    loc_spec = table("seq_locator")
    loc = pa.Table.from_pylist(new_loc, schema=loc_spec.schema)
    if locator is not None:
        loc = pa.concat_tables([locator, loc])
    loc = loc.sort_by("raw_digest")
    _write(loc, root / "l0seq/seq_locator.parquet", 122_880)

    man_spec = table("l0_manifest")
    man_path = root / "l0/l0_manifest.parquet"
    man = pa.Table.from_pylist(manifest_rows, schema=man_spec.schema)
    if man_path.exists():
        man = pa.concat_tables([pq.read_table(man_path), man])
    man = man.sort_by([("source_uid", "ascending"), ("object_path", "ascending")])
    _write(man, man_path, 122_880)

    return {"objects": len(objects), "records": len(occ_rows), "new_sequences": len(new_seqs), "index_rows": len(idx_rows)}


def _lengths_from_store(root: Path, locator: pa.Table, digests: set[bytes]) -> dict[bytes, int]:
    store = SeqStore(root)
    return {rd: store.length_of(rd) for rd in digests}


class SeqStore:
    """Reader: slices by (contig_uid, start, end, strand), 0-based half-open."""

    def __init__(self, root: Path):
        self.root = Path(root)
        self.chunk_size = json.loads((self.root / "l0seq/store_params.json").read_text())["chunk_size"]
        loc = pq.read_table(self.root / "l0seq/seq_locator.parquet").to_pylist()
        self._loc = {r["raw_digest"]: r for r in loc}
        idx = _existing_index(self.root)
        self._default: dict[bytes, bytes] = {}
        self._len: dict[bytes, int] = {}
        for r in idx.to_pylist():
            self._len[r["raw_digest"]] = r["contig_len"]
            if r["is_default_variant"]:
                self._default[r["contig_uid"]] = r["raw_digest"]
        self._files: dict[str, pq.ParquetFile] = {}
        self._rg_cache: dict[tuple[str, int], list[str]] = {}  # small LRU of decoded row groups

    _RG_CACHE_MAX = 64

    def _row_group_seqs(self, rel: str, rg: int) -> list[str]:
        key = (rel, rg)
        hit = self._rg_cache.pop(key, None)
        if hit is None:
            hit = self._file(rel).read_row_group(rg, columns=["seq"]).column(0).to_pylist()
            if len(self._rg_cache) >= self._RG_CACHE_MAX:
                self._rg_cache.pop(next(iter(self._rg_cache)))
        self._rg_cache[key] = hit
        return hit

    def _file(self, rel: str) -> pq.ParquetFile:
        if rel not in self._files:
            self._files[rel] = pq.ParquetFile(self.root / rel)
        return self._files[rel]

    def length_of(self, rd: bytes) -> int:
        if rd in self._len:
            return self._len[rd]
        loc = self._loc[rd]
        pf = self._file(loc["file"])
        rg = pf.read_row_group(loc["row_group"], columns=["contig_len"])
        return rg.column(0)[loc["first_chunk_row"] - self._rg_start(pf, loc["row_group"])].as_py()

    @staticmethod
    def _rg_start(pf: pq.ParquetFile, rg: int) -> int:
        return sum(pf.metadata.row_group(i).num_rows for i in range(rg))

    def raw(self, rd: bytes, start: int = 0, end: int | None = None) -> bytes:
        """Forward-strand bytes of one stored variant."""
        loc = self._loc[rd]
        length = self.length_of(rd)
        end = length if end is None else end
        if not 0 <= start <= end <= length:
            raise ValueError(f"slice {start}:{end} outside 0:{length}")
        if start == end:
            return b""
        cs = self.chunk_size
        first_ix, last_ix = start // cs, (end - 1) // cs
        pf = self._file(loc["file"])
        rows_wanted = range(loc["first_chunk_row"] + first_ix, loc["first_chunk_row"] + last_ix + 1)
        parts: list[str] = []
        offset = 0
        for rg in range(pf.metadata.num_row_groups):
            n = pf.metadata.row_group(rg).num_rows
            lo, hi = offset, offset + n
            offset = hi
            if hi <= rows_wanted.start or lo >= rows_wanted.stop:
                continue
            seqs = self._row_group_seqs(loc["file"], rg)
            for row in range(max(lo, rows_wanted.start), min(hi, rows_wanted.stop)):
                parts.append(seqs[row - lo])
        joined = "".join(parts).encode("ascii")
        base = first_ix * cs
        return joined[start - base : end - base]

    def get(
        self, contig_uid: bytes, start: int = 0, end: int | None = None, strand: int = 1, circular: bool = False
    ) -> bytes:
        """Slice the default variant. On a circular contig `end` may exceed the length (origin wrap)."""
        rd = self._default[contig_uid]
        n = self.length_of(rd)
        if circular and end is not None and end > n:
            if not 0 <= start < n < end <= start + n:
                raise ValueError(f"wrapping slice {start}:{end} invalid for length {n}")
            seq = self.raw(rd, start, n) + self.raw(rd, 0, end - n)
        else:
            seq = self.raw(rd, start, end)
        if strand == -1:
            return revcomp(seq)
        if strand != 1:
            raise ValueError("strand must be +1 or -1")
        return seq

    def reconstruct_object(self, l0_object_sha256: bytes) -> list[bytes]:
        """FASTA lines of one source object, from occurrence rows and chunks (R-02)."""
        occ = pa.concat_tables(
            pq.read_table(f) for f in sorted((self.root / "l0seq/contig_occurrence").glob("stage=*/source=*/part-0.parquet"))
        )
        rows = sorted(
            (r for r in occ.to_pylist() if r["l0_object_sha256"] == l0_object_sha256),
            key=lambda r: r["record_ordinal"],
        )
        lines: list[bytes] = []
        for r in rows:
            lines.extend(fasta.format_record(r["original_header"], self.raw(r["raw_digest"]), r["line_width"]))
        return lines

"""Write and read the L0-seq store (operations O-INGEST, O-SLICE, O-EXTRACT).

Layout under `<root>/` (planning/data_layer_and_tests.md §2A.1, §3.1):

    l0/l0_manifest.parquet
    l0seq/store_params.json
    l0seq/seq_shards/stage=<S>/shard-<NNNN>.fa.gz   (+ .fa.gz.fai, .fa.gz.gzi)
    l0seq/contig_occurrence/stage=<S>/source=<name>/part-0.parquet
    l0seq/contig_index/stage=<S>/part-0.parquet
    l0seq/seq_locator.parquet

Sequences are bgzip FASTA shards (DEC-17, after K-13): each stored byte sequence once,
named by its `raw_digest` hex, records sorted by `raw_digest`, a shard closed once it
holds `shard_bases` bases. samtools, pysam and MMseqs2 read the shards as they are.

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

from biobox import fasta, ids, ledger, provenance
from biobox.l0seq.bgzf import BgzfFasta, write_bgzf_fasta
from biobox.registry.schemas import table

STORE_FORMAT = "bgzf-fasta-1"
SHARD_BASES = 2**30  # ~0.33 GB on disk at 0.31 bytes/bp: inside the 256 MB-1 GB band (LY-04)
LINE_WIDTH = 80
BGZF_LEVEL = 6

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
    # only the columns the store uses, so files from before an appended column still concatenate (LY-09)
    cols = ["contig_uid", "raw_digest", "contig_len", "is_default_variant"]
    return pa.concat_tables(pq.read_table(f, columns=cols) for f in sorted(d.glob("stage=*/part-0.parquet")))


def ingest(
    objects: list[SourceObject],
    root: Path,
    stage: str,
    shard_bases: int = SHARD_BASES,
    line_width: int = LINE_WIDTH,
) -> dict:
    """Ingest FASTA source objects for one stage. Returns counts for the ledger."""
    root = Path(root)
    params_path = root / "l0seq" / "store_params.json"
    params = {"format": STORE_FORMAT, "shard_bases": shard_bases, "line_width": line_width, "bgzf_level": BGZF_LEVEL}
    if params_path.exists():
        if json.loads(params_path.read_text()) != params:
            raise ValueError("store parameters differ from the existing store")
    else:
        params_path.parent.mkdir(parents=True, exist_ok=True)
        params_path.write_text(json.dumps(params, sort_keys=True))

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

    tools = provenance.tool_versions("pyarrow")
    manifest_rows, occ_rows, ledger_rows = [], [], []
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
        n_in = n_kept = 0
        for ordinal, rec in enumerate(fasta.read(Path(obj.path))):
            n_in += 1
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
            if rd not in known and rd not in new_seqs:
                new_seqs[rd] = rec.seq
                n_kept += 1
            key = (obj.precedence, obj_sha, ordinal)
            if (cu, rd) not in variants or key < variants[(cu, rd)]:
                variants[(cu, rd)] = key
        ledger_rows.append(ledger.row(obj_sha, "l0_object", n_in, n_kept, n_in - n_kept, tools=tools))

    # seq_shards: new sequences sorted by raw_digest, cut into shards of about shard_bases
    groups: list[list[bytes]] = [[]]
    held = 0
    for rd in sorted(new_seqs):
        n = len(new_seqs[rd])
        if groups[-1] and held + n > shard_bases:
            groups.append([])
            held = 0
        groups[-1].append(rd)
        held += n
    new_loc = []
    for i, group in enumerate(g for g in groups if g):
        rel = f"l0seq/seq_shards/stage={stage}/shard-{i:04d}.fa.gz"
        write_bgzf_fasta(((rd.hex(), new_seqs[rd]) for rd in group), root / rel, line_width, BGZF_LEVEL)
        new_loc += [{"raw_digest": rd, "shard": rel, "contig_len": len(new_seqs[rd])} for rd in group]

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
        lengths.update(zip(locator.column("raw_digest").to_pylist(), locator.column("contig_len").to_pylist()))
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

    ledger.write(root, stage, "ingest", ledger_rows)
    return {"objects": len(objects), "records": len(occ_rows), "new_sequences": len(new_seqs), "index_rows": len(idx_rows)}


class SeqStore:
    """Reader: slices by (contig_uid, start, end, strand), 0-based half-open.

    `seq_locator` maps a raw_digest to its shard; the shard's `.fai` and `.gzi` find the
    blocks. `counters()` sums blocks decompressed and bytes read over all shards (LY-06).
    """

    def __init__(self, root: Path, cache_blocks: int = 64):
        self.root = Path(root)
        loc = pq.read_table(self.root / "l0seq/seq_locator.parquet")
        self._loc: dict[bytes, tuple[str, int]] = dict(zip(
            loc.column("raw_digest").to_pylist(),
            zip(loc.column("shard").to_pylist(), loc.column("contig_len").to_pylist())))
        idx = _existing_index(self.root)
        self._default: dict[bytes, bytes] = {
            r["contig_uid"]: r["raw_digest"] for r in idx.to_pylist() if r["is_default_variant"]}
        self._cache_blocks = cache_blocks
        self._shards: dict[str, BgzfFasta] = {}

    def _shard(self, rel: str) -> BgzfFasta:
        if rel not in self._shards:
            self._shards[rel] = BgzfFasta(self.root / rel, cache_blocks=self._cache_blocks)
        return self._shards[rel]

    def shard_of(self, rd: bytes) -> str:
        return self._loc[rd][0]

    def counters(self) -> dict[str, int]:
        s = self._shards.values()
        return {"blocks_read": sum(b.blocks_read for b in s),
                "compressed_bytes_read": sum(b.compressed_bytes_read for b in s),
                "uncompressed_bytes_decoded": sum(b.uncompressed_bytes_decoded for b in s)}

    def close(self) -> None:
        for b in self._shards.values():
            b.close()

    def length_of(self, rd: bytes) -> int:
        return self._loc[rd][1]

    def raw(self, rd: bytes, start: int = 0, end: int | None = None) -> bytes:
        """Forward-strand bytes of one stored variant."""
        rel, length = self._loc[rd]
        end = length if end is None else end
        if not 0 <= start <= end <= length:
            raise ValueError(f"slice {start}:{end} outside 0:{length}")
        return self._shard(rel).fetch(rd.hex(), start, end)

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

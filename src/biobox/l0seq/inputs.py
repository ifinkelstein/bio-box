"""Persisted genome inputs (card B-04): the `genome_input` table, one row per contig.

`GenomeInput` is the gene caller's description of one genome; it lives here, below L1, and
`biobox.l1.genecall` imports it (card B-05). `write_genome_input` turns `GenomeInput` objects into
rows sorted by `contig_uid` and streams them through `biobox.io.write_table`; `iter_genome_input`
streams the objects back one genome at a time in `genome_uid` order, and `read_genome_input` is
its list. `tier` and `is_species_rep` are not stored (they belong to `contig_dim`) and come back
at the dataclass defaults.
"""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

import duckdb
import pyarrow as pa

from biobox import ids, io
from biobox.registry.schemas import table

_SPEC = table("genome_input")
_READ_BATCH_ROWS = 65_536
_DUCK_MEMORY_LIMIT = "1GB"  # sorts larger than this spill to a private temporary directory


@dataclass(frozen=True)
class GenomeInput:
    contig_uids: tuple[bytes, ...]
    genetic_code: int = 11  # used in single mode; meta modes report the code per gene
    kind: str = "prokaryote"
    circular: tuple[bytes, ...] = ()  # contig_uids that are circular
    assembly_quality: str = "complete"
    tier: int = 1
    is_species_rep: bool = False


def _path(root, stage: str) -> Path:
    return Path(root) / _SPEC.location / f"stage={stage}" / "part-0.parquet"


def _rows(genomes: Iterable[GenomeInput]) -> list[tuple]:
    rows, seen = [], set()
    for g in genomes:
        members = set(g.contig_uids)
        if len(members) != len(g.contig_uids):
            raise ValueError("a genome lists the same contig twice")
        stray = set(g.circular) - members
        if stray:
            raise ValueError(f"circular names contig {sorted(stray)[0].hex()} that the genome does not contain")
        dup = members & seen
        if dup:
            raise ValueError(f"contig {sorted(dup)[0].hex()} appears in more than one genome")
        seen |= members
        guid = ids.genome_uid(list(g.contig_uids))
        circ = set(g.circular)
        rows += [(c, guid, c in circ, g.kind, g.genetic_code, g.assembly_quality) for c in g.contig_uids]
    rows.sort(key=lambda r: r[0])
    return rows


def _batches(rows: list[tuple]) -> Iterator[pa.RecordBatch]:
    for i in range(0, len(rows), io.ROW_GROUP_ROWS):
        cols = list(zip(*rows[i : i + io.ROW_GROUP_ROWS]))
        yield pa.RecordBatch.from_arrays(
            [pa.array(cols[k], type=f.type) for k, f in enumerate(_SPEC.schema)], schema=_SPEC.schema)


def write_genome_input(root, stage: str, genomes: Iterable[GenomeInput]) -> dict:
    """Write `<root>/l0seq/genome_input/stage=<stage>/part-0.parquet`; refuses to overwrite it."""
    rows = _rows(genomes)
    io.write_table(_batches(rows), _SPEC, _path(root, stage))
    return {"genomes": len({r[1] for r in rows}), "contigs": len(rows)}


def _genome(guid: bytes, rs: list[dict]) -> GenomeInput:
    for col in ("kind", "genetic_code", "assembly_quality"):
        if len({r[col] for r in rs}) != 1:
            raise ValueError(f"genome {guid.hex()}: its contigs disagree on {col}")
    uids = tuple(r["contig_uid"] for r in rs)
    if ids.genome_uid(list(uids)) != guid:
        raise ValueError(f"genome {guid.hex()}: stored genome_uid does not match its contigs")
    return GenomeInput(uids, rs[0]["genetic_code"], kind=rs[0]["kind"],
                       circular=tuple(r["contig_uid"] for r in rs if r["is_circular"]),
                       assembly_quality=rs[0]["assembly_quality"])


def iter_genome_input(root, stage: str, tmp=None) -> Iterator[GenomeInput]:
    """Stream the genomes of a stage one at a time, sorted by `genome_uid` (contigs by `contig_uid`).

    DuckDB sorts the table and an Arrow reader hands it over in batches, so only the rows of the
    genome being assembled are held. `tmp` holds DuckDB's spill directory (system temp by default).
    """
    path = _path(root, stage)
    if not path.exists():
        raise FileNotFoundError(f"{path}: no genome_input for stage {stage}")
    scratch = Path(tempfile.mkdtemp(prefix="biobox-genome-input-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": _DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        lit = str(path).replace("'", "''")
        reader = con.sql(f"SELECT * FROM read_parquet('{lit}') ORDER BY genome_uid, contig_uid"
                         ).to_arrow_reader(_READ_BATCH_ROWS)
        guid, rs = None, []
        for b in reader:
            for r in b.to_pylist():
                if r["genome_uid"] != guid:
                    if rs:
                        yield _genome(guid, rs)
                    guid, rs = r["genome_uid"], []
                rs.append(r)
        if rs:
            yield _genome(guid, rs)
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)


def read_genome_input(root, stage: str) -> list[GenomeInput]:
    """Rebuild the genomes of a stage, sorted by `genome_uid`."""
    return list(iter_genome_input(root, stage))

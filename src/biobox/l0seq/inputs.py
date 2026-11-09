"""Persisted genome inputs (card B-04): the `genome_input` table, one row per contig.

`write_genome_input` turns `GenomeInput` objects into rows sorted by `contig_uid` and streams them
through `biobox.io.write_table`; `read_genome_input` rebuilds the objects. `tier` and
`is_species_rep` are not stored (they belong to `contig_dim`) and come back at the dataclass defaults.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox import ids, io
from biobox.l1.genecall import GenomeInput
from biobox.registry.schemas import table

_SPEC = table("genome_input")


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


def read_genome_input(root, stage: str) -> list[GenomeInput]:
    """Rebuild the genomes of a stage, sorted by `genome_uid`."""
    t = pq.read_table(_path(root, stage))
    groups: dict[bytes, list[dict]] = {}
    for r in t.to_pylist():
        groups.setdefault(r["genome_uid"], []).append(r)
    out = []
    for guid in sorted(groups):
        rs = sorted(groups[guid], key=lambda r: r["contig_uid"])
        for col in ("kind", "genetic_code", "assembly_quality"):
            if len({r[col] for r in rs}) != 1:
                raise ValueError(f"genome {guid.hex()}: its contigs disagree on {col}")
        uids = tuple(r["contig_uid"] for r in rs)
        if ids.genome_uid(list(uids)) != guid:
            raise ValueError(f"genome {guid.hex()}: stored genome_uid does not match its contigs")
        out.append(GenomeInput(uids, rs[0]["genetic_code"], kind=rs[0]["kind"],
                               circular=tuple(r["contig_uid"] for r in rs if r["is_circular"]),
                               assembly_quality=rs[0]["assembly_quality"]))
    return out

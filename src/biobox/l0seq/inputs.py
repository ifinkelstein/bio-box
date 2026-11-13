"""Persisted genome inputs (card B-04): the `genome_input` table, one row per contig and genome.

`GenomeInput` is the gene caller's description of one genome; it lives here, below L1, and
`biobox.l1.genecall` imports it (card B-05). `write_genome_input` turns `GenomeInput` objects into
rows sorted by (`contig_uid`, `genome_uid`) and streams them through `biobox.io.write_table`;
`iter_genome_input` streams the objects back one genome at a time in `genome_uid` order, and
`read_genome_input` is its list. `tier` and `is_species_rep` are not stored (they belong to
`contig_dim`) and come back at the dataclass defaults.

A contig may be in several genomes of one stage (DEC-27, card C-21a): it has one row in each. The
genomes must agree on its call inputs (kind, genetic code, topology), so that one set of gene calls
serves them all; the genome with the least `genome_uid` writes them (`shared_contigs`,
`store.ContigOwners`). A genome's own attributes (`assembly_quality`) may differ.
"""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

import duckdb
import pyarrow as pa

from biobox import ids, io, layout, ledger
from biobox.registry.schemas import table

_SPEC = table("genome_input")
_READ_BATCH_ROWS = 65_536
_DUCK_MEMORY_LIMIT = "1GB"  # sorts larger than this spill to a private temporary directory
# estimate: zstd bytes per row (contig_uid and genome_uid, 16 random bytes each, + small columns); cuts
# parts only. No measurement yet (planning/measurements/2026-11-09_design_stress_test.md §3 lists none).
EST_ROW_BYTES = 36


@dataclass(frozen=True)
class GenomeInput:
    contig_uids: tuple[bytes, ...]
    genetic_code: int = 11  # used in single mode; meta modes report the code per gene
    kind: str = "prokaryote"
    circular: tuple[bytes, ...] = ()  # contig_uids that are circular
    assembly_quality: str = "complete"
    tier: int = 1
    is_species_rep: bool = False


def _dir(root, stage: str) -> Path:
    return Path(root) / _SPEC.location / f"stage={stage}"


# the columns of a row that are the contig's call inputs: genomes that share a contig must agree on them (DEC-27)
_CALL_INPUTS = ((2, "is_circular"), (3, "kind"), (4, "genetic_code"))


def _rows(genomes: Iterable[GenomeInput]) -> list[tuple]:
    """The rows of `genomes`, sorted by (contig_uid, genome_uid). Refused: a genome that lists a contig twice or
    names a circular contig it does not hold; one genome given twice (the same contigs: one `genome_uid`); a contig
    whose genomes disagree on its call inputs (`_CALL_INPUTS`), named with both genomes."""
    rows, guids = [], set()
    for g in genomes:
        members = set(g.contig_uids)
        if len(members) != len(g.contig_uids):
            raise ValueError("a genome lists the same contig twice")
        stray = set(g.circular) - members
        if stray:
            raise ValueError(f"circular names contig {sorted(stray)[0].hex()} that the genome does not contain")
        guid = ids.genome_uid(list(g.contig_uids))
        if guid in guids:
            raise ValueError(f"genome {guid.hex()} is listed twice (two entries with the same contigs)")
        guids.add(guid)
        circ = set(g.circular)
        rows += [(c, guid, c in circ, g.kind, g.genetic_code, g.assembly_quality) for c in g.contig_uids]
    rows.sort(key=lambda r: (r[0], r[1]))
    for a, b in zip(rows, rows[1:]):  # agreement of neighbours is agreement of all a contig's rows
        if a[0] != b[0]:
            continue
        for i, col in _CALL_INPUTS:
            if a[i] != b[i]:
                raise ValueError(
                    f"contig {a[0].hex()} is in genomes {a[1].hex()} and {b[1].hex()}, which disagree on its {col} "
                    f"({a[i]!r} against {b[i]!r}): a contig shared by genomes of one stage has one set of gene calls, "
                    "made from one set of call inputs (DEC-27; planning/design/C-21_shared_contigs.md §3)")
    return rows


def _batches(rows: list[tuple]) -> Iterator[pa.RecordBatch]:
    for i in range(0, len(rows), io.ROW_GROUP_ROWS):
        cols = list(zip(*rows[i : i + io.ROW_GROUP_ROWS]))
        yield pa.RecordBatch.from_arrays(
            [pa.array(cols[k], type=f.type) for k, f in enumerate(_SPEC.schema)], schema=_SPEC.schema)


def write_genome_input(root, stage: str, genomes: Iterable[GenomeInput],
                       max_file_bytes: int = layout.MAX_FILE_BYTES, out=None) -> dict:
    """Write `<out>/l0seq/genome_input/stage=<stage>/part-<k>.parquet`, cut at `max_file_bytes` of
    estimated size (`io.write_parts`). `out` is the output root (default `root`; under run-step,
    `$BIOBOX_STAGING`, card C-08); refuses, before reading `genomes`, when the stage's partition
    already holds a file in `root` or in `out`, and refuses the genomes `_rows` refuses. Returns the
    genomes and distinct contigs and, when a contig is in several genomes, the rows and those contigs."""
    out = Path(root) if out is None else Path(out)
    ledger.refuse_written(root, out, [_dir("", stage).as_posix()], f"genome_input of stage {stage}")
    rows = _rows(genomes)
    io.write_parts(_batches(rows), _SPEC, _dir(out, stage), EST_ROW_BYTES, max_file_bytes)
    contigs = len({r[0] for r in rows})
    res = {"genomes": len({r[1] for r in rows}), "contigs": contigs}
    shared = len({a[0] for a, b in zip(rows, rows[1:]) if a[0] == b[0]})
    return {**res, "rows": len(rows), "shared_contigs": shared} if shared else res


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
    files = io.part_files(_dir(root, stage))
    if not files:
        raise FileNotFoundError(f"{_dir(root, stage)}: no genome_input for stage {stage}")
    scratch = Path(tempfile.mkdtemp(prefix="biobox-genome-input-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": _DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in files)
        reader = con.sql(f"SELECT * FROM read_parquet([{lits}], hive_partitioning = false) "
                         "ORDER BY genome_uid, contig_uid").to_arrow_reader(_READ_BATCH_ROWS)
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


def shared_contigs(root, stage: str, tmp=None) -> dict[bytes, bytes]:
    """For each contig in more than one genome of the stage's `genome_input`, the least of their `genome_uid`s: the
    genome that writes the contig's gene calls (DEC-27). One DuckDB aggregate (spilling under `tmp`); only the
    shared contigs are held (4 of 54,265 contigs at 50 Gbp were shared at all, *measured*, card C-20)."""
    files = io.part_files(_dir(root, stage))
    if not files:
        raise FileNotFoundError(f"{_dir(root, stage)}: no genome_input for stage {stage}")
    scratch = Path(tempfile.mkdtemp(prefix="biobox-genome-input-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": _DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in files)
        return dict(con.execute(f"SELECT contig_uid, min(genome_uid) FROM read_parquet([{lits}], hive_partitioning = "
                                "false) GROUP BY contig_uid HAVING count(*) > 1 ORDER BY contig_uid").fetchall())
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)


def stage_genomes(root, stage: str, genomes: Iterable[GenomeInput] | None, tmp=None
                  ) -> tuple[Iterable[GenomeInput], dict[bytes, bytes]]:
    """The genomes a stage's producer reads and `shared_contigs` of them. None reads the stage's `genome_input`
    (streamed, `iter_genome_input`); genomes given explicitly (tests, pilots) are held as a list, which is what
    their caller already holds."""
    if genomes is None:
        return iter_genome_input(root, stage, tmp=tmp), shared_contigs(root, stage, tmp=tmp)
    genomes = list(genomes)
    least: dict[bytes, bytes] = {}
    n: dict[bytes, int] = {}
    for g in genomes:
        guid = ids.genome_uid(list(g.contig_uids))
        for cu in g.contig_uids:
            least[cu] = min(guid, least.get(cu, guid))
            n[cu] = n.get(cu, 0) + 1
    return genomes, {cu: least[cu] for cu in sorted(n) if n[cu] > 1}

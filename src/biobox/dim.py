"""Dimension tables (card B-12): `contig_dim`, one row per default-variant contig for a release.

Attributes that can change without the sequence changing live here and nowhere else (AGENTS §3.6).
This card creates the table and the place where later steps fill it. What it can fill today comes
from `contig_index` (`contig_len`) and `genome_input` (`genome_uid`, `is_circular`,
`assembly_quality`). Taxonomy and tiering steps do not exist yet, so `tax_order` and `tier` are null
and `is_species_rep` is false for every row; those steps will write them.

Writes `<root>/dim/contig_dim/release=<release>/part-<k>.parquet` (cut at `max_file_bytes` of
estimated size, `io.write_parts`, card B-19), once: the writer refuses to overwrite, so a release's
dimension table cannot be rewritten.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import duckdb

from biobox import catalog, io, layout, ledger
from biobox.l1.genecall import DUCK_MEMORY_LIMIT
from biobox.registry.schemas import table

# estimate: zstd bytes per contig_dim row (contig_uid and genome_uid, 16 random bytes each, + small
# columns); cuts parts only. No measurement yet (planning/measurements/2026-11-09_design_stress_test.md).
EST_ROW_BYTES = 44


def _sql_list(files: list[Path]) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def run(root, release: str, stages: list[str], tmp=None, max_file_bytes: int = layout.MAX_FILE_BYTES,
        out=None) -> dict:
    """Build `contig_dim` for `release` over `stages`, the release's stages (card B-25): a stage the
    release does not list is not read, a listed stage without `contig_index` files raises
    FileNotFoundError, and the release's `release_def` is written, or checked to list the same stages
    (`catalog.define_release`).

    Inputs are read from `root`; the table and (when the root has none) the `release_def` are written
    under `out` (default `root`; card C-08). Refuses, before any work, when the release's partition
    already holds a file in `root` or `out`; with `out` given, `tmp` may lie under neither.

    Errors: a default-variant contig without a `genome_input` row, a `genome_input` row whose
    contig is not a default-variant contig of the given stages, and a contig whose `genome_input`
    rows differ between stages. Identical rows from several stages collapse to one.
    Returns `{"contigs": n, "genomes": g, "circular": c}`.
    """
    root = Path(root)
    if out is not None:
        ledger.check_scratch(tmp, root, out)
    out = root if out is None else Path(out)
    stages = catalog.check_stages(stages)
    spec = table("contig_dim")
    rel_base = f"{spec.location}/release={release}"
    catalog.verify_release_def(root, release, stages)  # a disagreeing definition is named before a refusal
    ledger.refuse_written(root, out, [rel_base], f"contig_dim of release {release} (a release's dimension table "
                                                  "is written once)")
    base = out / rel_base
    if not stages:
        raise ValueError(f"release {release}: no stages")
    idx, gin = [], []
    for s in stages:  # in the release's order; the output is sorted by contig_uid and does not depend on it
        found = io.part_files(root / f"l0seq/contig_index/stage={s}")
        if not found:
            raise FileNotFoundError(f"stage {s}: no contig_index files under {root}")
        idx += found
        gin += io.part_files(root / f"l0seq/genome_input/stage={s}")
    catalog.define_release(root, release, stages, out=out)
    scratch = Path(tempfile.mkdtemp(prefix="biobox-dim-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        con.execute("CREATE VIEW ci AS SELECT DISTINCT contig_uid, contig_len "
                    f"FROM read_parquet({_sql_list(idx)}) WHERE is_default_variant")
        if gin:
            con.execute("CREATE VIEW gi AS SELECT contig_uid, genome_uid, is_circular, assembly_quality "
                        f"FROM read_parquet({_sql_list(gin)})")
        else:
            con.execute("CREATE VIEW gi AS SELECT NULL::BLOB AS contig_uid, NULL::BLOB AS genome_uid, "
                        "NULL::BOOLEAN AS is_circular, NULL::VARCHAR AS assembly_quality WHERE false")
        n, first = con.execute("SELECT count(*), min(contig_uid) FROM ci "
                               "WHERE contig_uid NOT IN (SELECT contig_uid FROM gi)").fetchall()[0]
        if n:
            raise ValueError(f"{n} default-variant contigs have no genome_input row; first: {first.hex()}")
        n, first = con.execute("SELECT count(*), min(contig_uid) FROM gi "
                               "WHERE contig_uid NOT IN (SELECT contig_uid FROM ci)").fetchall()[0]
        if n:
            raise ValueError(f"{n} genome_input rows name a contig that contig_index does not hold as a "
                             f"default variant; first: {first.hex()}")
        n, first = con.execute("SELECT count(*), min(contig_uid) FROM (SELECT contig_uid FROM gi "
                               "GROUP BY contig_uid HAVING count(DISTINCT (genome_uid, is_circular, "
                               "assembly_quality)) > 1)").fetchall()[0]
        if n:
            raise ValueError(f"{n} contigs have conflicting genome_input rows across stages; first: {first.hex()}")
        con.execute("CREATE VIEW dim AS SELECT ci.contig_uid, g.genome_uid, ci.contig_len, g.is_circular, "
                    "NULL::UTINYINT AS tier, g.assembly_quality, NULL::UINTEGER AS tax_order, "
                    "false AS is_species_rep FROM ci JOIN (SELECT DISTINCT contig_uid, genome_uid, "
                    "is_circular, assembly_quality FROM gi) g USING (contig_uid)")
        contigs, genomes, circular = con.execute(
            "SELECT count(*), count(DISTINCT genome_uid), count(*) FILTER (WHERE is_circular) FROM dim").fetchall()[0]
        io.write_parts(io.duck_batches(con, "SELECT * FROM dim ORDER BY contig_uid", spec), spec, base,
                       EST_ROW_BYTES, max_file_bytes)
        return {"contigs": contigs, "genomes": genomes, "circular": circular}
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)

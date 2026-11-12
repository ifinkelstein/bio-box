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
from biobox.l0seq.store import stored_by
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

    A contig shared by genomes of two stages (DEC-26, card C-20) belongs to the stage that stored it
    first, the one whose `contig_index` holds its default variant: its row is that stage's genome, and a
    different row of another stage is skipped and listed (`_views`). Identical rows from several stages
    collapse to one. A contig in several genomes of its owning stage (DEC-27, card C-21a) keeps the
    genome with the least `genome_uid`, the one that writes its gene calls; the others are listed as
    skipped too, until card C-21b counts a contig's genomes.

    Errors: a default-variant contig without a `genome_input` row; a `genome_input` row whose contig
    is not a default-variant contig of the given stages (naming the stage on disk that owns it, which
    the release must then list); rows for a contig from the stages that own it that disagree on its
    topology, or two different rows for one genome (two stages both holding its default variant).
    Returns `{"contigs": n, "genomes": g, "circular": c}`, and `"skipped"` (one entry per skipped row:
    contig, stage, genome, owning stage, the genome kept) when there is one.
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
    for s in stages:  # every listed stage has contig_index files
        if not io.part_files(root / f"l0seq/contig_index/stage={s}"):
            raise FileNotFoundError(f"stage {s}: no contig_index files under {root}")
    catalog.define_release(root, release, stages, out=out)
    scratch = Path(tempfile.mkdtemp(prefix="biobox-dim-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        why = _views(con, root, release, stages)
        if why is not None:
            raise ValueError(why)
        contigs, genomes, circular = con.execute(
            "SELECT count(*), count(DISTINCT genome_uid), count(*) FILTER (WHERE is_circular) FROM dim").fetchall()[0]
        skipped = [{"contig_uid": cu.hex(), "stage": s, "genome_uid": g.hex(), "owner_stage": o,
                    "kept_genome_uid": k.hex()} for cu, s, g, o, k in con.execute(
                        "SELECT contig_uid, stage, genome_uid, owner_stage, kept_genome_uid FROM skipped "
                        "ORDER BY contig_uid, stage, genome_uid").fetchall()]
        io.write_parts(io.duck_batches(con, "SELECT * FROM dim ORDER BY contig_uid", spec), spec, base,
                       EST_ROW_BYTES, max_file_bytes)
        res = {"contigs": contigs, "genomes": genomes, "circular": circular}
        return {**res, "skipped": skipped} if skipped else res
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)


def _stage_files(root: Path, rel: str, stages: list[str], cols: str) -> str:
    """The rows of every listed stage's files under `rel`, each with its `stage`, as one SQL query (None: no file)."""
    parts = [f"SELECT '{s}' AS stage, {cols} FROM read_parquet({_sql_list(fs)}, hive_partitioning = false)"
             for s in stages if (fs := io.part_files(root / rel / f"stage={s}"))]
    return " UNION ALL ".join(parts) if parts else None


def _views(con, root: Path, release: str, stages: list[str]) -> str | None:
    """Make the views `dim` (the table's rows) and `skipped` (genome_input rows of a contig that another stage
    owns, DEC-26) over the release's stages, or return why the release is refused.

    A contig's owning stage is the stage whose `contig_index` holds its default variant: ingest stores a
    sequence once, in the first stage that sees it. Its row comes from that stage's `genome_input`, from the
    genome with the least `genome_uid` when the stage lists it in several (DEC-27); every other row that
    differs is skipped (that genome is one contig short here: DEC-26's consequence, until card C-21b). Without
    a row in an owning stage, the row is chosen the same way from all its rows. Refused: a default-variant
    contig without any row; a row naming a contig that no listed stage owns (the stage on disk that owns it is
    named); candidate rows that disagree on `is_circular`; two different rows for the chosen genome (two
    stages that both hold the default variant)."""
    ci = _stage_files(root, "l0seq/contig_index", stages, "contig_uid, contig_len, is_default_variant")
    gi = _stage_files(root, "l0seq/genome_input", stages, "contig_uid, genome_uid, is_circular, assembly_quality")
    con.execute(f"CREATE TEMP TABLE own AS SELECT DISTINCT stage, contig_uid, contig_len FROM ({ci}) "
                "WHERE is_default_variant")
    con.execute("CREATE TEMP TABLE gi AS " + (gi if gi else
                "SELECT NULL::VARCHAR AS stage, NULL::BLOB AS contig_uid, NULL::BLOB AS genome_uid, "
                "NULL::BOOLEAN AS is_circular, NULL::VARCHAR AS assembly_quality WHERE false"))
    con.execute("CREATE TEMP VIEW ci AS SELECT DISTINCT contig_uid, contig_len FROM own")
    n, first = con.execute("SELECT count(*), min(contig_uid) FROM ci "
                           "WHERE contig_uid NOT IN (SELECT contig_uid FROM gi)").fetchall()[0]
    if n:
        return f"{n} default-variant contigs have no genome_input row; first: {first.hex()}"
    n, first = con.execute("SELECT count(*), min(contig_uid) FROM gi "
                           "WHERE contig_uid NOT IN (SELECT contig_uid FROM ci)").fetchall()[0]
    if n:
        msg = (f"{n} genome_input rows name a contig that contig_index does not hold as a default variant; "
               f"first: {first.hex()}")
        owner = stored_by(root, [first], skip=set(stages)).get(first)
        if owner is not None:
            (first_stage,), = con.execute("SELECT min(stage) FROM gi WHERE contig_uid = ?", [first]).fetchall()
            msg += (f" (in stage {first_stage}'s genome_input): stage {owner} stored it first and owns it (DEC-26), "
                    f"and release {release} lists stages {stages} without it; list stage {owner} in the release")
        return msg
    # the rows each contig's row is chosen from: its owning stage's, else all it has
    con.execute("CREATE TEMP TABLE gio AS SELECT gi.*, (o.stage IS NOT NULL) AS from_owner FROM gi "
                "LEFT JOIN own o ON o.stage = gi.stage AND o.contig_uid = gi.contig_uid")
    con.execute("CREATE TEMP TABLE cand AS SELECT contig_uid, genome_uid, is_circular, assembly_quality FROM gio g "
                "WHERE from_owner OR NOT EXISTS (SELECT 1 FROM gio x WHERE x.contig_uid = g.contig_uid AND x.from_owner)")
    n, first = con.execute("SELECT count(*), min(contig_uid) FROM (SELECT contig_uid FROM cand "
                           "GROUP BY contig_uid HAVING count(DISTINCT is_circular) > 1)").fetchall()[0]
    if n:
        return (f"{n} contigs have conflicting genome_input rows in the stages that own them: their genomes disagree "
                f"on is_circular, the contig's own topology; first: {first.hex()}")
    # several genomes of the owning stage (DEC-27): the one with the least genome_uid, which writes the gene calls
    con.execute("CREATE TEMP TABLE chosen AS SELECT DISTINCT c.* FROM cand c JOIN "
                "(SELECT contig_uid, min(genome_uid) AS genome_uid FROM cand GROUP BY contig_uid) USING (contig_uid, "
                "genome_uid)")
    n, first = con.execute("SELECT count(*), min(contig_uid) FROM (SELECT contig_uid FROM chosen "
                           "GROUP BY contig_uid HAVING count(*) > 1)").fetchall()[0]
    if n:
        return (f"{n} contigs have conflicting genome_input rows in the stages that own them (two stages holding the "
                f"default variant, with different rows for one genome); first: {first.hex()}")
    con.execute("CREATE TEMP VIEW skipped AS SELECT g.contig_uid, g.stage, g.genome_uid, "
                "(SELECT min(o.stage) FROM own o WHERE o.contig_uid = g.contig_uid) AS owner_stage, "
                "c.genome_uid AS kept_genome_uid FROM gio g JOIN chosen c USING (contig_uid) "
                "WHERE g.genome_uid IS DISTINCT FROM c.genome_uid "
                "OR g.is_circular IS DISTINCT FROM c.is_circular OR g.assembly_quality IS DISTINCT FROM c.assembly_quality")
    con.execute("CREATE TEMP VIEW dim AS SELECT ci.contig_uid, c.genome_uid, ci.contig_len, c.is_circular, "
                "NULL::UTINYINT AS tier, c.assembly_quality, NULL::UINTEGER AS tax_order, "
                "false AS is_species_rep FROM ci JOIN chosen c USING (contig_uid)")
    return None


def owner_problem(root, release: str, stages: list[str], tmp=None) -> str | None:
    """Why `contig_dim` of `release` over `stages` would be refused (`run`'s checks, DEC-26), or None; reads only.
    `build-stage` calls it before a release's clustering (card C-20), so a release that lists a later stage
    without the stage that owns one of its contigs is refused before hours of work, naming that stage."""
    root = Path(root)
    stages = catalog.check_stages(stages)
    for s in stages:
        if not io.part_files(root / f"l0seq/contig_index/stage={s}"):
            return f"stage {s}: no contig_index files under {root}"
    scratch = Path(tempfile.mkdtemp(prefix="biobox-dim-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        return _views(con, root, release, stages)
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)

"""Dimension tables (card B-12): `contig_dim`, one row per default-variant contig for a release, and
`genome_dim`, one row per genome (card C-21b).

Attributes that can change without the sequence changing live here and nowhere else (AGENTS §3.6).
What can be filled today comes from `contig_index` (`contig_len`) and `genome_input` (`genome_uid`,
`is_circular`, `kind`, `genetic_code`, `assembly_quality`). Taxonomy and tiering steps do not exist
yet, so `tax_order` and `tier` are null and `is_species_rep` is false for every row; those steps will
write them.

A contig may belong to several genomes (DEC-27): the release's membership is the union of its stages'
`genome_input` rows (the catalog reads it through the view `genome_member`). `contig_dim` keeps one
row per contig: its `genome_uid` is the least of the contig's genomes, `n_genomes` their number, and
a genome attribute (`assembly_quality`) is filled only where all those genomes agree. `genome_dim`
holds each genome whole: its distinct contigs and their bases, and `kind`, `genetic_code` and
`assembly_quality` where the genome's rows in the release's stages agree. Both are functions of the
release's set of stages, not of the order its `release_def` lists them nor of the order they were
built in (design note `planning/design/C-21_shared_contigs.md` §2a, §4).

Writes `<root>/dim/contig_dim/release=<release>/part-<k>.parquet` and
`<root>/dim/genome_dim/release=<release>/part-<k>.parquet` (cut at `max_file_bytes` of estimated
size, `io.write_parts`, card B-19), once: the writer refuses to overwrite, so a release's dimension
tables cannot be rewritten.
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
# estimate: zstd bytes per genome_dim row (genome_uid, 16 random bytes, + small columns); cuts parts only.
# No measurement yet
EST_GENOME_ROW_BYTES = 30


def _sql_list(files: list[Path]) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def run(root, release: str, stages: list[str], tmp=None, max_file_bytes: int = layout.MAX_FILE_BYTES,
        out=None) -> dict:
    """Build `contig_dim` and `genome_dim` for `release` over `stages`, the release's stages (card
    B-25): a stage the release does not list is not read, a listed stage without `contig_index` files
    raises FileNotFoundError, and the release's `release_def` is written, or checked to list the same
    stages (`catalog.define_release`).

    Inputs are read from `root`; the tables and (when the root has none) the `release_def` are written
    under `out` (default `root`; card C-08). Refuses, before any work, when the release's partition of
    either table already holds a file in `root` or `out`; with `out` given, `tmp` may lie under neither.

    A contig in several genomes (DEC-27, card C-21b), of one stage or of several: its `contig_dim` row
    names the least of its genomes and counts them (`n_genomes`); every genome is whole in `genome_dim`.
    Its topology (`is_circular`, a fact of the contig) comes from the rows of the stage that stored it
    first, the one whose `contig_index` holds its default variant and whose gene calls it carries
    (DEC-26); identical rows from several stages count once.

    Errors: a default-variant contig without a `genome_input` row; a `genome_input` row whose contig
    is not a default-variant contig of the given stages (naming the stage on disk that owns it, which
    the release must then list); rows for a contig from the stages that own it that disagree on its
    topology, or two different rows for one genome (two stages both holding its default variant).
    Returns `{"contigs": n, "genomes": g, "circular": c}` (`g` the rows of `genome_dim`), and
    `"shared_contigs"`, the contigs in more than one genome, when there is one.
    """
    root = Path(root)
    if out is not None:
        ledger.check_scratch(tmp, root, out)
    out = root if out is None else Path(out)
    stages = catalog.check_stages(stages)
    spec, gspec = table("contig_dim"), table("genome_dim")
    rel_base, g_base = f"{spec.location}/release={release}", f"{gspec.location}/release={release}"
    catalog.verify_release_def(root, release, stages)  # a disagreeing definition is named before a refusal
    ledger.refuse_written(root, out, [rel_base, g_base], f"contig_dim and genome_dim of release {release} "
                                                         "(a release's dimension tables are written once)")
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
        contigs, circular, shared = con.execute(
            "SELECT count(*), count(*) FILTER (WHERE is_circular), count(*) FILTER (WHERE n_genomes > 1) "
            "FROM dim").fetchall()[0]
        (genomes,), = con.execute("SELECT count(*) FROM gdim").fetchall()
        io.write_parts(io.duck_batches(con, "SELECT * FROM dim ORDER BY contig_uid", spec), spec, base,
                       EST_ROW_BYTES, max_file_bytes)
        io.write_parts(io.duck_batches(con, "SELECT * FROM gdim ORDER BY genome_uid", gspec), gspec, out / g_base,
                       EST_GENOME_ROW_BYTES, max_file_bytes)
        res = {"contigs": contigs, "genomes": genomes, "circular": circular}
        return {**res, "shared_contigs": shared} if shared else res
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)


def _stage_files(root: Path, rel: str, stages: list[str], cols: str) -> str:
    """The rows of every listed stage's files under `rel`, each with its `stage`, as one SQL query (None: no file)."""
    parts = [f"SELECT '{s}' AS stage, {cols} FROM read_parquet({_sql_list(fs)}, hive_partitioning = false)"
             for s in stages if (fs := io.part_files(root / rel / f"stage={s}"))]
    return " UNION ALL ".join(parts) if parts else None


def _views(con, root: Path, release: str, stages: list[str]) -> str | None:
    """Make the views `dim` (`contig_dim`'s rows) and `gdim` (`genome_dim`'s) over the release's stages, or
    return why the release is refused.

    Membership is every `genome_input` row of the release's stages (DEC-27): a contig's `genome_uid` is the least
    of its genomes, `n_genomes` their number, and its `assembly_quality` the one value its genomes' rows hold, else
    null. A genome's `n_contigs` and `genome_len` count its distinct contigs; its `kind`, `genetic_code` and
    `assembly_quality` are the one value its rows hold, else null. Aggregates of a set of rows: no stage order enters.

    A contig's topology is the contig's own fact and comes from its owning stage, the stage whose `contig_index`
    holds its default variant (ingest stores a sequence once, in the first stage that sees it, and that stage
    writes its gene calls and tape, DEC-26); without a row there, from all its rows. Refused: a default-variant
    contig without any row; a row naming a contig that no listed stage owns (the stage on disk that owns it is
    named); owning rows that disagree on `is_circular`; two different owning rows for one genome (two stages
    that both hold the default variant)."""
    ci = _stage_files(root, "l0seq/contig_index", stages, "contig_uid, contig_len, is_default_variant")
    gi = _stage_files(root, "l0seq/genome_input", stages,
                      "contig_uid, genome_uid, is_circular, kind, genetic_code, assembly_quality")
    con.execute(f"CREATE TEMP TABLE own AS SELECT DISTINCT stage, contig_uid, contig_len FROM ({ci}) "
                "WHERE is_default_variant")
    con.execute("CREATE TEMP TABLE gi AS " + (gi if gi else
                "SELECT NULL::VARCHAR AS stage, NULL::BLOB AS contig_uid, NULL::BLOB AS genome_uid, "
                "NULL::BOOLEAN AS is_circular, NULL::VARCHAR AS kind, NULL::UTINYINT AS genetic_code, "
                "NULL::VARCHAR AS assembly_quality WHERE false"))
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
    # the rows a contig's topology is read from: its owning stage's, else all it has
    con.execute("CREATE TEMP TABLE gio AS SELECT gi.*, (o.stage IS NOT NULL) AS from_owner FROM gi "
                "LEFT JOIN own o ON o.stage = gi.stage AND o.contig_uid = gi.contig_uid")
    con.execute("CREATE TEMP TABLE cand AS SELECT contig_uid, genome_uid, is_circular, assembly_quality FROM gio g "
                "WHERE from_owner OR NOT EXISTS (SELECT 1 FROM gio x WHERE x.contig_uid = g.contig_uid AND x.from_owner)")
    n, first = con.execute("SELECT count(*), min(contig_uid) FROM (SELECT contig_uid FROM cand "
                           "GROUP BY contig_uid HAVING count(DISTINCT is_circular) > 1)").fetchall()[0]
    if n:
        return (f"{n} contigs have conflicting genome_input rows in the stages that own them: their genomes disagree "
                f"on is_circular, the contig's own topology; first: {first.hex()}")
    n, first = con.execute("SELECT count(*), min(contig_uid) FROM (SELECT DISTINCT contig_uid FROM cand "
                           "GROUP BY contig_uid, genome_uid HAVING count(DISTINCT assembly_quality) > 1)").fetchall()[0]
    if n:
        return (f"{n} contigs have conflicting genome_input rows in the stages that own them (two stages holding the "
                f"default variant, with different rows for one genome); first: {first.hex()}")

    def agreed(col: str) -> str:  # the one value of `col` in a group's rows, else null
        return f"CASE WHEN count(DISTINCT {col}) = 1 THEN min({col}) END"

    con.execute("CREATE TEMP TABLE topo AS SELECT contig_uid, bool_or(is_circular) AS is_circular FROM cand "
                "GROUP BY contig_uid")
    con.execute("CREATE TEMP TABLE memb AS SELECT contig_uid, min(genome_uid) AS genome_uid, "
                f"count(DISTINCT genome_uid)::UINTEGER AS n_genomes, {agreed('assembly_quality')} AS assembly_quality "
                "FROM gi GROUP BY contig_uid")
    con.execute("CREATE TEMP VIEW dim AS SELECT ci.contig_uid, m.genome_uid, ci.contig_len, t.is_circular, "
                "NULL::UTINYINT AS tier, m.assembly_quality, NULL::UINTEGER AS tax_order, "
                "false AS is_species_rep, m.n_genomes FROM ci JOIN memb m USING (contig_uid) "
                "JOIN topo t USING (contig_uid)")
    con.execute("CREATE TEMP TABLE gattr AS SELECT genome_uid, "
                f"{agreed('kind')} AS kind, {agreed('genetic_code')} AS genetic_code, "
                f"{agreed('assembly_quality')} AS assembly_quality FROM gi GROUP BY genome_uid")
    con.execute("CREATE TEMP TABLE gsize AS SELECT genome_uid, count(*)::UINTEGER AS n_contigs, "
                "sum(ci.contig_len)::BIGINT AS genome_len "
                "FROM (SELECT DISTINCT genome_uid, contig_uid FROM gi) JOIN ci USING (contig_uid) GROUP BY genome_uid")
    con.execute("CREATE TEMP VIEW gdim AS SELECT genome_uid, s.n_contigs, s.genome_len, a.kind, a.genetic_code, "
                "a.assembly_quality, NULL::UTINYINT AS tier, NULL::UINTEGER AS tax_order, false AS is_species_rep "
                "FROM gsize s JOIN gattr a USING (genome_uid)")
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

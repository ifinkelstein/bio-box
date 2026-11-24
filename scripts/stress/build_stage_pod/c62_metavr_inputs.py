"""Card C-62 (DEC-71): stage MetaVR vOTU representatives and >= 30 kb jumbo/large UViGs as stage MV1.

    python c62_metavr_inputs.py stage PARQUET FNA_BGZ OUT_DIR [--root ROOT] [--root-stages S1,P1,...]
        [--owner-inputs NAME=DIR ...] [--max-contigs 3000000] [--max-bases 40000000000]
        [--retrieved-at 2026-10-23T15:22:00Z] [--final-dir DIR] [--work WORK_DIR] [--bgzip]

Selection rule (deterministic, independent of row order):
1. Exclude all `is_concatemer == True` rows (deferred with reason `concatemer`).
2. Within each `votu`, choose one best non-concatemer representative ordered by:
   quality tier (`Complete` < `High-quality` < `Medium-quality` < `Low-quality` < `Not-determined`),
   highest `min(completeness, 100.0)` (NULL -> -1.0), lowest `checkv_contamination` (NULL -> 999.0),
   longest `length`, least `uvig`.
3. Priority tiers for `MV1` (non-concatemer):
   - Tier 1 (`jumbo_ge200k`): all `length >= 200_000` UViGs (reps first, then non-reps).
   - Tier 2 (`rep_complete_hq`): `Complete` and `High-quality` vOTU reps (`< 200 kb`, `contamination <= 10.0`).
   - Tier 3 (`rep_mq_ge30k`): `Medium-quality` vOTU reps (`30 kb .. < 200 kb`, `contamination <= 10.0`).
   - Tier 4 (`rep_mq_lt30k`): `Medium-quality` vOTU reps (`< 30 kb`, `contamination <= 10.0`).
   - Tier 5 (`rep_large_ge30k_lq_nd`): `Low-quality` and `Not-determined` vOTU reps (`30 kb .. < 200 kb`,
     `contamination <= 10.0`), `High` viral_confidence first, then longest.
4. Take the longest strict prefix of Priority 1..5 within `max_contigs` (3,000,000) and `max_bases`
   (40,000,000,000 bp). Every remaining priority candidate goes to `MV1_left_for_later.tsv.gz` with
   reason `mv2_spillover` in exact priority order (`mv2_rank`), alongside `concatemer`,
   `high_contamination`, `not_representative`, and `lq_nd_lt30k` so zero rows are dropped silently (DEC-59).
5. Stream `METAVR.fna.bgz` once in bounded memory, apply DEC-65 exact `20..1,000` bp DTR trimming on
   single-contig `Direct terminal repeat` records (recording `trims.tsv`), deduplicate repeated contig sets
   into `repeated_records.tsv` (DEC-27), and check `contig_uid` owners against earlier stages into
   `shared_records.tsv` (DEC-26).
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import shutil
import subprocess
import sys
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parent))

import c48c_virus_inputs as v2mod  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402
from biobox.l0seq import store  # noqa: E402

Refused = make_inputs.Refused

DEFAULT_MAX_CONTIGS = 3_000_000
DEFAULT_MAX_BASES = 40_000_000_000
MAX_CONTAMINATION = 10.0
JUMBO_BP = 200_000
LARGE_BP = 30_000

MV1_FASTA = "metavr_mv1.fna"
LEFT_GZ = "MV1_left_for_later.tsv.gz"
METAVR_LICENSE = (
    "MetaVR public download (https://portal.nersc.gov/cfs/m342/METAVR/, Fiamenghi et al. NAR 2026, "
    "doi:10.1093/nar/gkaf1283); contains JGI IMG-derived data, cite per JGI Data Utilization Policy; "
    "held local, not redistributed (DEC-50, DEC-51, DEC-71)"
)
METAVR_SOURCE = {
    "source_name": "metavr",
    "source_release": "MetaVR (2026-08-26 METAVR.fna.bgz / 2026-06-24 METAVR_main_table.parquet); "
                      "MV1 selection by c62_metavr_inputs.py (DEC-65, DEC-71)",
    "license": METAVR_LICENSE,
    "redistributable": False,
}

TRIM_HEAD = v2mod.TRIM_HEAD
RECORDS_HEAD = [
    "ordinal", "source", "genome", "votu", "contig_uid", "length", "quality", "viral_confidence",
    "completeness", "contamination", "topology", "source_genetic_code", "genetic_code", "circular",
    "tier", "is_rep", "genome_type", "origin_source", "ecosystem", "host_taxonomy",
    "dtr_class", "source_contig_uid", "source_length", "trimmed_bases",
]
LEFT_HEAD = [
    "reason", "mv2_rank", "uvig", "votu", "quality", "viral_confidence", "length",
    "completeness", "checkv_contamination", "uvig_topology", "is_concatemer", "genetic_code",
]


def _sql_lit(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _sql_list(files: Iterable[Path]) -> str:
    return "[" + ", ".join(_sql_lit(f) for f in files) + "]"


def _uvig_from_header(header: str) -> str:
    return header.split("|", 1)[0].split(None, 1)[0] if header else ""


def apply_dec65(topology: str, seq: bytes) -> tuple[bytes, bool, str, int]:
    """DEC-65 DTR trimming for a single-contig MetaVR UViG.

    Returns (staged_seq, circular, dtr_class, trimmed_k).
    """
    if topology != "Direct terminal repeat":
        return seq, False, "", 0
    e = v2mod.ends(seq)
    dtr_cls = v2mod.classify(e)
    if dtr_cls == "exact":
        return v2mod.trim(seq, e.exact), True, dtr_cls, e.exact
    if dtr_cls == "none":
        return seq, True, dtr_cls, 0
    return seq, False, dtr_cls, 0


def _build_selection_tables(
    con: duckdb.DuckDBPyConnection,
    parquet_path: Path,
    *,
    max_contigs: int,
    max_bases: int,
) -> dict:
    """Create `mv_selected` and `mv_deferred` tables in DuckDB and return summary stats."""
    pq_lit = _sql_lit(parquet_path)
    con.execute(f"""
        CREATE TEMP TABLE mv_raw AS
        SELECT
            uvig,
            votu,
            coalesce(quality, '') AS quality,
            coalesce(viral_confidence, '') AS viral_confidence,
            coalesce(uvig_topology, '') AS uvig_topology,
            coalesce(is_concatemer, false) AS is_concatemer,
            completeness,
            checkv_contamination,
            length::BIGINT AS length,
            coalesce(CAST(genetic_code AS VARCHAR), '') AS genetic_code,
            coalesce(genome_type, '') AS genome_type,
            coalesce(source, '') AS source,
            coalesce(Ecosystem, '') AS ecosystem,
            coalesce(host_taxonomy, '') AS host_taxonomy,
            CASE quality
                WHEN 'Complete' THEN 0
                WHEN 'High-quality' THEN 1
                WHEN 'Medium-quality' THEN 2
                WHEN 'Low-quality' THEN 3
                WHEN 'Not-determined' THEN 4
                ELSE 5
            END AS qrank,
            least(coalesce(completeness, -1.0), 100.0) AS comp_cap,
            coalesce(checkv_contamination, 0.0) AS contam_eff,
            coalesce(checkv_contamination, 999.0) AS contam_rep
        FROM read_parquet({pq_lit})
    """)

    total_rows = int(con.execute("SELECT count(*) FROM mv_raw").fetchone()[0])

    con.execute("""
        CREATE TEMP TABLE mv_reps AS
        SELECT uvig FROM (
            SELECT
                uvig,
                row_number() OVER (
                    PARTITION BY votu
                    ORDER BY qrank ASC, comp_cap DESC, contam_rep ASC, length DESC, uvig ASC
                ) AS rn
            FROM mv_raw
            WHERE NOT is_concatemer
        ) WHERE rn = 1
    """)

    con.execute(f"""
        CREATE TEMP TABLE mv_classified AS
        SELECT
            m.*,
            (r.uvig IS NOT NULL) AS is_rep,
            CASE
                WHEN m.is_concatemer THEN 'concatemer'
                WHEN m.length >= {JUMBO_BP} THEN 'jumbo_ge200k'
                WHEN r.uvig IS NULL THEN 'not_representative'
                WHEN m.contam_eff > {MAX_CONTAMINATION} THEN 'high_contamination'
                WHEN m.quality IN ('Complete', 'High-quality') THEN 'rep_complete_hq'
                WHEN m.quality = 'Medium-quality' AND m.length >= {LARGE_BP} THEN 'rep_mq_ge30k'
                WHEN m.quality = 'Medium-quality' THEN 'rep_mq_lt30k'
                WHEN m.length >= {LARGE_BP} THEN 'rep_large_ge30k_lq_nd'
                ELSE 'lq_nd_lt30k'
            END AS cls
        FROM mv_raw m
        LEFT JOIN mv_reps r ON m.uvig = r.uvig
    """)

    con.execute("""
        CREATE TEMP TABLE mv_candidates AS
        SELECT
            *,
            row_number() OVER (
                ORDER BY
                    CASE cls
                        WHEN 'jumbo_ge200k' THEN 1
                        WHEN 'rep_complete_hq' THEN 2
                        WHEN 'rep_mq_ge30k' THEN 3
                        WHEN 'rep_mq_lt30k' THEN 4
                        WHEN 'rep_large_ge30k_lq_nd' THEN 5
                    END ASC,
                    -- Tier-specific tie-breakers matching brute_select oracle exactly
                    CASE cls
                        WHEN 'jumbo_ge200k' THEN (CASE WHEN is_rep THEN 0 ELSE 1 END)
                        WHEN 'rep_complete_hq' THEN qrank
                        WHEN 'rep_large_ge30k_lq_nd' THEN (CASE WHEN viral_confidence = 'High' THEN 0 ELSE 1 END)
                        ELSE 0
                    END ASC,
                    CASE cls
                        WHEN 'jumbo_ge200k' THEN qrank
                        ELSE 0
                    END ASC,
                    CASE cls
                        WHEN 'rep_complete_hq' THEN -comp_cap
                        WHEN 'rep_mq_lt30k' THEN -comp_cap
                        ELSE 0.0
                    END ASC,
                    CASE cls
                        WHEN 'rep_complete_hq' THEN contam_eff
                        ELSE 0.0
                    END ASC,
                    -length ASC,
                    CASE cls
                        WHEN 'rep_mq_ge30k' THEN -comp_cap
                        ELSE 0.0
                    END ASC,
                    contam_eff ASC,
                    uvig ASC
            ) AS prio_ord
        FROM mv_classified
        WHERE cls IN ('jumbo_ge200k', 'rep_complete_hq', 'rep_mq_ge30k', 'rep_mq_lt30k', 'rep_large_ge30k_lq_nd')
    """)

    con.execute("""
        CREATE TEMP TABLE mv_cum AS
        SELECT
            *,
            sum(length) OVER (ORDER BY prio_ord ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)::BIGINT AS cum_bp
        FROM mv_candidates
    """)

    first_over = con.execute(
        f"SELECT min(prio_ord) FROM mv_cum WHERE prio_ord > {int(max_contigs)} OR cum_bp > {int(max_bases)}"
    ).fetchone()[0]

    if first_over is None:
        cutoff_sql = "TRUE"
        spill_sql = "FALSE"
    else:
        cutoff_sql = f"prio_ord < {int(first_over)}"
        spill_sql = f"prio_ord >= {int(first_over)}"

    con.execute(f"""
        CREATE TEMP TABLE mv_selected AS
        SELECT
            uvig, votu, quality, viral_confidence, uvig_topology, completeness, checkv_contamination,
            length, genetic_code, genome_type, source, ecosystem, host_taxonomy, is_rep,
            cls AS tier, prio_ord
        FROM mv_cum
        WHERE {cutoff_sql}
        ORDER BY prio_ord
    """)

    con.execute(f"""
        CREATE TEMP TABLE mv_deferred AS
        SELECT
            'mv2_spillover' AS reason,
            (prio_ord - coalesce({int(first_over) if first_over is not None else 0}, 0) + 1)::BIGINT AS mv2_rank,
            uvig, votu, quality, viral_confidence, length, completeness, checkv_contamination,
            uvig_topology, is_concatemer, genetic_code
        FROM mv_cum
        WHERE {spill_sql}
        UNION ALL
        SELECT
            cls AS reason,
            0::BIGINT AS mv2_rank,
            uvig, votu, quality, viral_confidence, length, completeness, checkv_contamination,
            uvig_topology, is_concatemer, genetic_code
        FROM mv_classified
        WHERE cls NOT IN ('jumbo_ge200k', 'rep_complete_hq', 'rep_mq_ge30k', 'rep_mq_lt30k', 'rep_large_ge30k_lq_nd')
    """)

    sel_cnt, sel_bp = con.execute("SELECT count(*), coalesce(sum(length), 0)::BIGINT FROM mv_selected").fetchone()
    def_cnt, def_bp = con.execute("SELECT count(*), coalesce(sum(length), 0)::BIGINT FROM mv_deferred").fetchone()
    if int(sel_cnt) + int(def_cnt) != total_rows:
        raise Refused(f"row accounting mismatch: {sel_cnt} + {def_cnt} != {total_rows}")

    by_tier = {
        r[0]: {"uvigs": int(r[1]), "bp": int(r[2])}
        for r in con.execute("SELECT tier, count(*), sum(length)::BIGINT FROM mv_selected GROUP BY 1 ORDER BY 1").fetchall()
    }
    by_reason = {
        r[0]: {"uvigs": int(r[1]), "bp": int(r[2])}
        for r in con.execute("SELECT reason, count(*), sum(length)::BIGINT FROM mv_deferred GROUP BY 1 ORDER BY 1").fetchall()
    }
    return {
        "total_rows": total_rows,
        "selected_uvigs": int(sel_cnt),
        "selected_bases": int(sel_bp),
        "deferred_uvigs": int(def_cnt),
        "deferred_bases": int(def_bp),
        "by_tier": by_tier,
        "by_deferred_reason": by_reason,
    }


def _fmt_num(val) -> str:
    if val is None:
        return ""
    return f"{val:.4f}".rstrip("0").rstrip(".")


def _stream_deferred_tsv_gz(con: duckdb.DuckDBPyConnection, out_gz: Path) -> int:
    """Stream `mv_deferred` sorted by (reason, mv2_rank, uvig) to a deterministic gzip (mtime=0)."""
    cur = con.execute("""
        SELECT reason, mv2_rank, uvig, votu, quality, viral_confidence, length,
               completeness, checkv_contamination, uvig_topology, is_concatemer, genetic_code
        FROM mv_deferred
        ORDER BY reason ASC, mv2_rank ASC, uvig ASC
    """)
    n = 0
    with open(out_gz, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(("\t".join(LEFT_HEAD) + "\n").encode("utf-8"))
        while True:
            batch = cur.fetchmany(100_000)
            if not batch:
                break
            lines = []
            for r in batch:
                lines.append(
                    f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t"
                    f"{_fmt_num(r[7])}\t{_fmt_num(r[8])}\t{r[9]}\t{str(bool(r[10])).lower()}\t{r[11]}\n"
                )
            gz.write("".join(lines).encode("utf-8"))
            n += len(batch)
    return n


def select_mv1(
    parquet_path: Path,
    *,
    max_contigs: int = DEFAULT_MAX_CONTIGS,
    max_bases: int = DEFAULT_MAX_BASES,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> tuple[list[dict], list[dict], dict]:
    """Pure selection returning (selected_rows, deferred_rows, stats) — suitable for unit tests and small tables."""
    cfg = {"threads": int(threads), "memory_limit": memory_limit}
    if work is not None:
        duck_dir = Path(work) / "duck"
        duck_dir.mkdir(parents=True, exist_ok=True)
        cfg["temp_directory"] = str(duck_dir)
    con = duckdb.connect(config=cfg)
    try:
        stats = _build_selection_tables(con, Path(parquet_path), max_contigs=max_contigs, max_bases=max_bases)
        sel_cols = [d[0] for d in con.execute("SELECT * FROM mv_selected LIMIT 0").description]
        sel_rows = [dict(zip(sel_cols, r)) for r in con.execute("SELECT * FROM mv_selected ORDER BY prio_ord").fetchall()]
        def_cols = [d[0] for d in con.execute("SELECT * FROM mv_deferred LIMIT 0").description]
        def_rows = [
            dict(zip(def_cols, r))
            for r in con.execute("SELECT * FROM mv_deferred ORDER BY reason ASC, mv2_rank ASC, uvig ASC").fetchall()
        ]
        return sel_rows, def_rows, stats
    finally:
        con.close()


def _owners_from_records(
    records_tsv: Path,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    *,
    work: Path,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict[bytes, str]:
    """Find first owning stage for each distinct contig_uid in records.tsv using DuckDB (bounded memory)."""
    if root is None and not owner_inputs:
        return {}
    duck_dir = work / "duck_own"
    duck_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(config={"threads": int(threads), "memory_limit": memory_limit, "temp_directory": str(duck_dir)})
    try:
        from biobox import io as bbio

        con.execute(
            f"CREATE TEMP TABLE v AS SELECT DISTINCT from_hex(contig_uid) AS contig_uid "
            f"FROM read_csv({_sql_lit(records_tsv)}, delim='\t', header=true, quote='', all_varchar=true)"
        )
        out: dict[bytes, str] = {}
        for s in root_stages:
            files = bbio.part_files(Path(root) / store.INDEX_DIR / f"stage={s}") if root is not None else []
            if root is not None and not files:
                raise Refused(f"stage {s}: no contig_index under {root}")
            if not files:
                continue
            for (cu,) in con.execute(
                f"SELECT DISTINCT contig_uid FROM read_parquet({_sql_list(files)}, hive_partitioning = false) "
                "WHERE is_default_variant AND contig_uid IN (SELECT contig_uid FROM v) ORDER BY 1"
            ).fetchall():
                out.setdefault(bytes(cu), s)
        for name, d in owner_inputs.items():
            gj = Path(d) / "genomes.json"
            if gj.is_file():
                p = _sql_lit(gj)
                for (cu,) in con.execute(
                    f"SELECT DISTINCT from_hex(u) FROM (SELECT unnest(contig_uids) AS u FROM read_json({p}, "
                    "format='array', columns={'contig_uids': 'VARCHAR[]'})) "
                    "WHERE from_hex(u) IN (SELECT contig_uid FROM v) ORDER BY 1"
                ).fetchall():
                    out.setdefault(bytes(cu), name)
                continue
            rec = Path(d) / "records.tsv"
            if rec.is_file():
                for (cu,) in con.execute(
                    f"SELECT DISTINCT from_hex(contig_uid) FROM read_csv({_sql_lit(rec)}, delim='\t', "
                    "header=true, quote='', all_varchar=true) "
                    "WHERE from_hex(contig_uid) IN (SELECT contig_uid FROM v) ORDER BY 1"
                ).fetchall():
                    out.setdefault(bytes(cu), name)
        return out
    finally:
        con.close()


def _open_fasta_stream(fna_path: Path, use_bgzip: bool):
    """Return (stream, cleanup_fn) for plain, .gz, or .bgz FASTA."""
    if use_bgzip and fna_path.suffix in (".bgz", ".gz"):
        cmd = shutil.which("bgzip") or shutil.which("gzip")
        if cmd:
            args = [cmd, "-dc", str(fna_path)]
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, bufsize=1 << 22)
            assert proc.stdout is not None

            def _close():
                proc.stdout.close()
                rc = proc.wait()
                if rc != 0:
                    raise Refused(f"{args[0]} -dc {fna_path} exited {rc}")

            return proc.stdout, _close
    if fna_path.suffix in (".bgz", ".gz"):
        fh = gzip.open(fna_path, "rb")
        return io.BufferedReader(fh, 1 << 22), fh.close
    fh = open(fna_path, "rb", buffering=1 << 22)
    return fh, fh.close


def stage_mv1(
    parquet_path: Path,
    fna_bgz: Path,
    out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    max_contigs: int = DEFAULT_MAX_CONTIGS,
    max_bases: int = DEFAULT_MAX_BASES,
    retrieved_at: str = "2026-10-23T15:22:00Z",
    final_dir: Path | None = None,
    fasta_dir: Path | None = None,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
    use_bgzip: bool = False,
) -> dict:
    """Select MV1 from METAVR_main_table.parquet, stream METAVR.fna.bgz once, apply DEC-65 DTR trimming,
    write staged inputs atomically to OUT_DIR (and large FASTA/deferred files to fasta_dir if given),
    and return a JSON report."""
    out = Path(out)
    if out.exists():
        raise Refused(f"{out} exists")
    stage_tmp = out.with_name(f"{out.name}.tmp")
    shutil.rmtree(stage_tmp, ignore_errors=True)
    stage_tmp.mkdir(parents=True)
    work_dir = Path(work) if work is not None else stage_tmp / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    large_out = Path(fasta_dir) if fasta_dir is not None else stage_tmp
    large_out.mkdir(parents=True, exist_ok=True)
    obj_dir = Path(final_dir) if final_dir is not None else (Path(fasta_dir) if fasta_dir is not None else out)

    try:
        rep = _stage_mv1(
            Path(parquet_path), Path(fna_bgz), obj_dir, stage_tmp, large_out, work_dir,
            root, list(root_stages), owner_inputs,
            max_contigs=max_contigs, max_bases=max_bases,
            retrieved_at=retrieved_at, threads=threads, memory_limit=memory_limit,
            use_bgzip=use_bgzip,
        )
    except BaseException:
        shutil.rmtree(stage_tmp, ignore_errors=True)
        raise

    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    stage_tmp.rename(out)

    meta_names = [
        "genomes.json", "objects.json", "records.tsv", "repeated_records.tsv",
        "shared_records.tsv", "trims.tsv",
    ]
    large_final = Path(fasta_dir) if fasta_dir is not None else out
    rep["sha256"] = {
        **{n: ids.sha256_file(out / n).hex() for n in meta_names},
        LEFT_GZ: ids.sha256_file(large_final / LEFT_GZ).hex(),
        MV1_FASTA: ids.sha256_file(large_final / MV1_FASTA).hex(),
    }
    rep["bytes"] = {
        **{n: (out / n).stat().st_size for n in meta_names},
        LEFT_GZ: (large_final / LEFT_GZ).stat().st_size,
        MV1_FASTA: (large_final / MV1_FASTA).stat().st_size,
    }
    return rep


def _stage_mv1(
    parquet_path: Path,
    fna_bgz: Path,
    obj_dir: Path,
    stage_tmp: Path,
    large_out: Path,
    work_dir: Path,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    *,
    max_contigs: int,
    max_bases: int,
    retrieved_at: str,
    threads: int,
    memory_limit: str,
    use_bgzip: bool,
) -> dict:
    duck_dir = work_dir / "duck_sel"
    duck_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(config={"threads": int(threads), "memory_limit": memory_limit, "temp_directory": str(duck_dir)})
    try:
        stats = _build_selection_tables(con, parquet_path, max_contigs=max_contigs, max_bases=max_bases)
        left_rows = _stream_deferred_tsv_gz(con, large_out / LEFT_GZ)
        # Load only selected UViGs (~1M rows) into a compact dict
        chosen: dict[str, tuple] = {}
        cur = con.execute("""
            SELECT uvig, votu, quality, viral_confidence, uvig_topology, completeness,
                   checkv_contamination, length, genetic_code, genome_type, source,
                   ecosystem, host_taxonomy, is_rep, tier
            FROM mv_selected
        """)
        while True:
            batch = cur.fetchmany(100_000)
            if not batch:
                break
            for r in batch:
                chosen[r[0]] = r[1:]
    finally:
        con.close()

    # Stream METAVR.fna.bgz once in bounded memory
    # uvig -> [(cu, staged_len, src_len, rec_ord, dtr_cls, src_cu, trimmed_k)]
    parts: dict[str, list[tuple[bytes, int, int, int, str, bytes, int]]] = {}
    circ: dict[str, bool] = {}
    trims: list[tuple] = []
    dtr_classes: Counter = Counter()
    n_fasta = n_recs = 0

    stream, close_stream = _open_fasta_stream(fna_bgz, use_bgzip=use_bgzip)
    try:
        with open(large_out / MV1_FASTA, "wb", buffering=1 << 22) as w:
            for header, r in v2mod.iter_fasta_taken(stream, lambda hd: _uvig_from_header(hd) in chosen):
                n_fasta += 1
                if r is None:
                    continue
                uvig = _uvig_from_header(header)
                meta = chosen[uvig]
                topo = meta[3]
                expected_len = int(meta[6])
                src_len = len(r.seq)
                src_cu = ids.contig_uid(r.seq)
                is_single = (uvig not in parts) and (src_len == expected_len)
                if is_single and topo == "Direct terminal repeat":
                    seq, is_c, dtr_cls, trimmed_k = apply_dec65(topo, r.seq)
                    dtr_classes[dtr_cls] += 1
                    circ[uvig] = is_c
                else:
                    seq, dtr_cls, trimmed_k = r.seq, "", 0
                    if uvig not in circ:
                        circ[uvig] = False
                    else:
                        circ[uvig] = False
                w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
                cu = ids.contig_uid(seq) if trimmed_k else src_cu
                if trimmed_k:
                    trims.append((
                        "metavr", uvig, r.accession, src_cu.hex(), cu.hex(),
                        src_len, trimmed_k, trimmed_k, len(seq),
                    ))
                parts.setdefault(uvig, []).append((cu, len(seq), src_len, n_recs, dtr_cls, src_cu, trimmed_k))
                n_recs += 1
    finally:
        close_stream()

    missing = sorted(set(chosen) - set(parts))
    if missing:
        raise Refused(f"{len(missing)} selected UViGs have no FASTA record; first {missing[0]}")
    multi_dtr = [u for u, p in parts.items() if chosen[u][3] == "Direct terminal repeat" and len(p) > 1]
    if multi_dtr:
        raise Refused(f"{len(multi_dtr)} Direct terminal repeat UViGs have >1 FASTA record; first {multi_dtr[0]}")
    length_mismatches = [
        {
            "uvig": u,
            "fasta_bp": sum(sn for _, _, sn, _, _, _, _ in p),
            "table_bp": int(chosen[u][6]),
            "topology": chosen[u][3],
            "records": len(p),
        }
        for u, p in sorted(parts.items())
        if sum(sn for _, _, sn, _, _, _, _ in p) != int(chosen[u][6])
    ]
    trims.sort(key=lambda t: (t[0], t[1], t[2]))

    # Write records.tsv first so DuckDB can query distinct contig_uids for owner lookup without big executemany
    order = sorted(parts, key=lambda u: parts[u][0][3])
    records_path = stage_tmp / "records.tsv"
    with open(records_path, "w", encoding="utf-8", buffering=1 << 20) as rf:
        rf.write("\t".join(RECORDS_HEAD) + "\n")
        ordinal = 0
        for uvig in order:
            votu, qual, conf, topo, comp, contam, _len, src_gc, gtype, orig_src, eco, host, is_rep, tier = chosen[uvig]
            is_c = circ[uvig]
            for cu, n, src_n, _, dtr_cls, src_cu, trimmed_k in parts[uvig]:
                rf.write(
                    f"{ordinal}\tmetavr\t{uvig}\t{votu}\t{cu.hex()}\t{n}\t{qual}\t{conf}\t"
                    f"{_fmt_num(comp)}\t{_fmt_num(contam)}\t{topo}\t{src_gc}\t11\t{str(is_c).lower()}\t"
                    f"{tier}\t{str(bool(is_rep)).lower()}\t{gtype}\t{orig_src}\t{eco}\t{host}\t"
                    f"{dtr_cls}\t{src_cu.hex()}\t{src_n}\t{trimmed_k}\n"
                )
                ordinal += 1

    owners = _owners_from_records(
        records_path, root, root_stages, owner_inputs,
        work=work_dir, threads=threads, memory_limit=memory_limit,
    )

    listing = make_inputs.Listing()
    genomes, repeated, shared = [], [], []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ordinal = 0
    for uvig in order:
        is_rep = bool(chosen[uvig][12])
        is_c = circ[uvig]
        uids = [cu for cu, _, _, _, _, _, _ in parts[uvig]]
        for cu, n, _, _, _, _, _ in parts[uvig]:
            contig_len.setdefault(cu, n)
            if cu in owners:
                shared.append((ordinal, uvig, cu.hex(), owners[cu]))
            ordinal += 1
        aq = "vOTU" if is_rep else "contig"
        g = listing.add(uids, 11, "virus", aq)
        guid = ids.genome_uid(list(dict.fromkeys(uids)))
        if g is None:
            repeated.append((ordinal - len(uids), uvig, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [uids[0].hex()] if is_c else []
        genomes.append(g)
        first_of[guid] = uvig

    make_inputs._write_tsv(stage_tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(stage_tmp / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"], shared)
    make_inputs._write_tsv(stage_tmp / "trims.tsv", TRIM_HEAD, trims)
    objects = [{"path": str(obj_dir / MV1_FASTA), **METAVR_SOURCE, "retrieved_at": retrieved_at}]
    (stage_tmp / "objects.json").write_text(json.dumps(objects, indent=1))
    (stage_tmp / "genomes.json").write_text(json.dumps(genomes))

    owner_contigs: Counter = Counter()
    owner_bases: Counter = Counter()
    for cu, st in owners.items():
        if cu in contig_len:
            owner_contigs[st] += 1
            owner_bases[st] += contig_len[cu]

    return {
        **stats,
        "fasta_records_scanned": n_fasta,
        "fasta_records_written": n_recs,
        "left_for_later_rows": left_rows,
        "multi_record_genomes": sum(1 for p in parts.values() if len(p) > 1),
        "source_bases": sum(sn for p in parts.values() for _, _, sn, _, _, _, _ in p),
        "staged_bases": sum(n for p in parts.values() for _, n, _, _, _, _, _ in p),
        "trims": len(trims),
        "trimmed_bases": sum(t[7] for t in trims),
        "dtr_classes": dict(sorted(dtr_classes.items())),
        "length_mismatches_count": len(length_mismatches),
        "length_mismatches": length_mismatches[:50],
        "genomes": len(genomes),
        "contigs": len(contig_len),
        "contig_bases": sum(contig_len.values()),
        "circular_genomes": sum(1 for g in genomes if g["circular"]),
        "repeated_genomes": len(repeated),
        "shared_records": len(shared),
        "owner_contigs": dict(sorted(owner_contigs.items())),
        "owner_bases": dict(sorted(owner_bases.items())),
        "new_contigs": sum(1 for cu in contig_len if cu not in owners),
        "new_bases": sum(n for cu, n in contig_len.items() if cu not in owners),
        **listing.report(),
        "objects": len(objects),
    }


def _parse_owners(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise Refused(f"--owner-inputs {item!r}: expected NAME=DIR")
        name, path = item.split("=", 1)
        out[name] = Path(path)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage MetaVR MV1 inputs (card C-62, DEC-71)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("stage")
    p.add_argument("parquet", type=Path)
    p.add_argument("fna_bgz", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--root-stages", default="")
    p.add_argument("--owner-inputs", action="append", default=[])
    p.add_argument("--max-contigs", type=int, default=DEFAULT_MAX_CONTIGS)
    p.add_argument("--max-bases", type=int, default=DEFAULT_MAX_BASES)
    p.add_argument("--retrieved-at", default="2026-10-23T15:22:00Z")
    p.add_argument("--final-dir", type=Path, default=None)
    p.add_argument("--fasta-dir", type=Path, default=None)
    p.add_argument("--work", type=Path, default=None)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--memory-limit", default="4GB")
    p.add_argument("--bgzip", action="store_true")
    p.add_argument("--report-out", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        stages = [s for s in a.root_stages.split(",") if s]
        rep = stage_mv1(
            a.parquet, a.fna_bgz, a.out,
            root=a.root, root_stages=stages, owner_inputs=_parse_owners(a.owner_inputs),
            max_contigs=a.max_contigs, max_bases=a.max_bases,
            retrieved_at=a.retrieved_at, final_dir=a.final_dir, fasta_dir=a.fasta_dir, work=a.work,
            threads=a.threads, memory_limit=a.memory_limit, use_bgzip=a.bgzip,
        )
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    txt = json.dumps(rep, indent=2, sort_keys=True)
    if a.report_out is not None:
        a.report_out.parent.mkdir(parents=True, exist_ok=True)
        a.report_out.write_text(txt + "\n")
    print(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

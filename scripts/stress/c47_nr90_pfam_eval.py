"""Card C-47 (Step 2): NR90 representative Pfam inference evaluation on stage C5 (DEC-32, DEC-35).

Implements the exact DEC-35 / B-47 §3 rules:
- Length rule: a non-representative member with `|L_m - L_r| > 0.2 * L_r` is scanned (`role = 'length'`).
- Alignment coverage rule (DEC-35 (4), B-47 §3.6): for every domain row of the representative,
  `f = max(0, min(env_end, tend) - max(env_start, tstart)) / (env_end - env_start)`.
  If ANY representative domain row has `f < 0.8` (or the alignment is missing), the member is
  scanned (`role = 'unmapped'`), never written with a guessed envelope.
- Envelope mapping (B-47 §3.1-3.4): CIGAR walk with query = representative (`M` consumes both,
  `I` representative only, `D` member only); `env_start` maps to the first `M` rep residue >= `env_start`,
  `env_end` maps to one past the last `M` rep residue < `env_end`.
- Stop gate (DEC-35 (6)): stops above 1% wrong pairs (`false / reported > 0.01`) OR above 1% missed
  pairs (`missed / true > 0.01`).
- Cross-stage comparison (DEC-35 (7)): compares stage-local NR90 scanned residue fraction against
  release-level NR90 (`r26.11.1 + C5`), flagging if stage-local scans > 10 percentage points more.

Usage (on wilkcomp02, read-only on /tmp/if2453/r10/root):
    python c47_nr90_pfam_eval.py eval --root ROOT --stage C5 \
        --stage-cm STAGE_CM.parquet --stage-aln STAGE_ALN.parquet \
        --selection C5/selection.tsv --objects C5/objects.json \
        --out OUT.json --tmp TMP_DIR \
        [--release-cm REL_CM.parquet --release-aln REL_ALN.parquet]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
import time
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path

import duckdb
import numpy as np

CIGAR_RE = re.compile(rb"(\d+)([MID])")
LIB = "pfam-38.2"
STOP_THRESHOLD = 0.01
DEC35_7_GAP_PP = 10.0


def length_mismatch(lm: int, lr: int) -> bool:
    """DEC-32 / B-47 §3.7: strict inequality |L_m - L_r| > 0.2 * L_r."""
    return abs(lm - lr) > 0.2 * lr


def envelope_coverage(env_start: int, env_end: int, tstart: int, tend: int) -> float:
    """B-47 §3.6: fraction of representative envelope [env_start, env_end) inside [tstart, tend)."""
    span = env_end - env_start
    if span <= 0:
        return 0.0
    ov = max(0, min(env_end, tend) - max(env_start, tstart))
    return ov / span


def classify_member(
    *,
    is_rep: bool,
    lm: int,
    lr: int,
    rep_envs: Sequence[tuple[int, int]],
    tstart: int | None,
    tend: int | None,
) -> str:
    """Classify a protein into 'rep', 'length', 'unmapped', or 'inferred' under DEC-35."""
    if is_rep:
        return "rep"
    if length_mismatch(lm, lr):
        return "length"
    if tstart is None or tend is None:
        return "unmapped"
    for es, ee in rep_envs:
        if envelope_coverage(es, ee, tstart, tend) < 0.8:
            return "unmapped"
    return "inferred"


def map_envelope(
    cigar: bytes,
    qs: int,
    qe: int,
    ts: int,
    te: int,
    env_start: int,
    env_end: int,
) -> tuple[int, int] | None:
    """B-47 §3.1-3.4: map representative envelope [env_start, env_end) to member coordinates.

    CIGAR has query = representative: M consumes both, I rep only, D member only.
    Returns (m_start, m_end) 0-based half-open, or None if no M column falls in [env_start, env_end).
    """
    rp, mp = ts, qs
    first_m: int | None = None
    last_m: int | None = None
    for num, op in CIGAR_RE.findall(cigar):
        n = int(num)
        if op == b"M":
            r0, r1 = rp, rp + n
            lo = max(r0, env_start)
            hi = min(r1, env_end)
            if lo < hi:
                m_lo = mp + (lo - r0)
                m_hi = mp + (hi - 1 - r0)
                if first_m is None:
                    first_m = m_lo
                last_m = m_hi
            rp += n
            mp += n
        elif op == b"I":
            rp += n
        else:  # b"D"
            mp += n
    if rp != te or mp != qe:
        raise ValueError(f"CIGAR end mismatch: walked ({rp},{mp}) != ({te},{qe})")
    if first_m is None or last_m is None:
        return None
    return first_m, last_m + 1


def f_bin(f: float) -> str:
    if f >= 1.0:
        return "1.0"
    if f >= 0.8:
        return "[0.8,1)"
    if f >= 0.5:
        return "[0.5,0.8)"
    if f > 0.0:
        return "(0,0.5)"
    return "0"


def _summarize_bucket(rows: list[dict]) -> dict:
    n_prot = len(rows)
    tot_res = sum(r["len"] for r in rows)
    reps = sum(1 for r in rows if r["role"] == "rep")
    len_sc = sum(1 for r in rows if r["role"] == "length")
    unm_sc = sum(1 for r in rows if r["role"] == "unmapped")
    inferred = sum(1 for r in rows if r["role"] == "inferred")
    sc_prot = reps + len_sc + unm_sc
    sc_res = sum(r["len"] for r in rows if r["role"] != "inferred")
    copied_prot = sum(1 for r in rows if r["copied"])
    copied_res = sum(r["len"] for r in rows if r["copied"])
    sc_prot_b47b = sum(1 for r in rows if r["role"] != "inferred" and not r["copied"])
    sc_res_b47b = sum(r["len"] for r in rows if r["role"] != "inferred" and not r["copied"])
    true_pairs = sum(r["true_pairs"] for r in rows)
    right = sum(r["right"] for r in rows)
    false = sum(r["false"] for r in rows)
    missed = sum(r["missed"] for r in rows)
    reported = right + false
    wrong_rate = (false / reported) if reported else 0.0
    missed_rate = (missed / true_pairs) if true_pairs else 0.0
    return {
        "proteins": n_prot,
        "total_residues": tot_res,
        "reps": reps,
        "length_scanned": len_sc,
        "unmapped_scanned": unm_sc,
        "inferred_members": inferred,
        "scanned_proteins": sc_prot,
        "scanned_residues": sc_res,
        "scanned_protein_fraction": (sc_prot / n_prot) if n_prot else 0.0,
        "scanned_residue_fraction": (sc_res / tot_res) if tot_res else 0.0,
        "copied_proteins": copied_prot,
        "copied_residues": copied_res,
        "scanned_proteins_after_b47b": sc_prot_b47b,
        "scanned_residues_after_b47b": sc_res_b47b,
        "scanned_residue_fraction_after_b47b": (sc_res_b47b / tot_res) if tot_res else 0.0,
        "true_pairs": true_pairs,
        "reported_pairs": reported,
        "pairs_right": right,
        "pairs_false": false,
        "pairs_missed": missed,
        "wrong_rate": wrong_rate,
        "missed_rate": missed_rate,
        "stop_triggered": (wrong_rate > STOP_THRESHOLD) or (missed_rate > STOP_THRESHOLD),
    }


def evaluate_clusters(
    proteins: dict[str, dict],
    clusters: dict[str, tuple[str, int | None, int | None]],
    domains: dict[str, list[tuple[str, int, int]]],
) -> dict:
    """In-memory DEC-35 evaluator (used by unit tests and small fixtures)."""
    per_row: list[dict] = []
    by_phylum: dict[str, list[dict]] = defaultdict(list)
    fbins: dict[str, dict[str, int]] = {
        b: {"rows": 0, "true_pair": 0} for b in ("1.0", "[0.8,1)", "[0.5,0.8)", "(0,0.5)", "0")
    }
    for uid, pinfo in sorted(proteins.items()):
        rep, ts, te = clusters[uid]
        is_rep = uid == rep
        lm = int(pinfo["len"])
        lr = int(proteins[rep]["len"])
        rep_doms = domains.get(rep, [])
        mem_doms = domains.get(uid, [])
        rep_envs = [(s, e) for _, s, e in rep_doms]
        role = classify_member(is_rep=is_rep, lm=lm, lr=lr, rep_envs=rep_envs, tstart=ts, tend=te)
        mem_set = {acc for acc, _, _ in mem_doms}
        rep_set = {acc for acc, _, _ in rep_doms}
        if not is_rep and not length_mismatch(lm, lr) and ts is not None and te is not None:
            for acc, es, ee in rep_doms:
                b = f_bin(envelope_coverage(es, ee, ts, te))
                fbins[b]["rows"] += 1
                if acc in mem_set:
                    fbins[b]["true_pair"] += 1
        if role == "inferred":
            right = len(mem_set & rep_set)
            false = len(rep_set - mem_set)
            missed = len(mem_set - rep_set)
        else:
            right = len(mem_set)
            false = 0
            missed = 0
        rec = {
            "uid": uid,
            "len": lm,
            "copied": bool(pinfo.get("copied", False)),
            "role": role,
            "true_pairs": len(mem_set),
            "right": right,
            "false": false,
            "missed": missed,
        }
        per_row.append(rec)
        by_phylum[pinfo.get("phylum", "unknown")].append(rec)
    return {
        "overall": _summarize_bucket(per_row),
        "by_phylum": {ph: _summarize_bucket(rs) for ph, rs in sorted(by_phylum.items())},
        "f_bins": fbins,
    }


def _sql_list(paths: Sequence[Path]) -> str:
    return "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for p in paths) + "]"


def _record_files(root: Path, stage: str, step: str, prefix: str) -> list[Path]:
    rec = json.loads((root / "runs" / f"stage={stage}" / f"{step}.json").read_text())
    return [root / o["path"] for o in rec["outputs"] if o["path"].startswith(prefix) and o["path"].endswith(".parquet")]


def _phylum_map_sql(
    con: duckdb.DuckDBPyConnection,
    root: Path,
    stage: str,
    selection_tsv: Path,
    objects_json: Path,
    gcv: str,
) -> None:
    """Create table `prot_phylum(protein_uid BLOB, phylum VARCHAR)` mapping each C5 protein to its phylum.

    Uses `selection.tsv` (accession -> phylum) and the stage's ingest files: `l0/l0_manifest` (`object_path`, the
    object's file name `<accession>.fna.gz`, and `object_sha256`), `l0seq/contig_occurrence` (`l0_object_sha256`,
    `contig_uid`) and the gene calls' `l1/feature_protein` (`contig_uid`, `protein_uid`). `objects.json` is not
    needed for the mapping (the manifest names each object). If a protein appears in multiple phyla, assigns the
    deterministic lexicographically first phylum (and records multi-phylum count separately).
    """
    acc_to_phylum: dict[str, str] = {}
    with open(selection_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row["accession"]
            ph = row.get("phylum") or ""
            if not ph and row.get("gtdb_taxonomy"):
                parts = row["gtdb_taxonomy"].split(";")
                ph = parts[1] if len(parts) > 1 else ""
            acc_to_phylum[acc] = ph or "unknown"

    def accession(name: str) -> str:
        stem = Path(name).name
        for suf in (".fna.gz", ".fna", ".fa.gz", ".fa"):
            if stem.endswith(suf):
                return stem[: -len(suf)]
        return stem

    l0_man = _record_files(root, stage, "ingest", f"l0/l0_manifest/stage={stage}/")
    occ = _record_files(root, stage, "ingest", f"l0seq/contig_occurrence/stage={stage}/")
    fp = _record_files(root, stage, "genecall", f"l1/feature_protein/gene_call_version={gcv}/stage={stage}/")
    obj_to_phylum = [
        (bytes(sha), acc_to_phylum.get(accession(path), "unknown"))
        for path, sha in con.execute(
            f"SELECT object_path, object_sha256 FROM read_parquet({_sql_list(l0_man)}, hive_partitioning = false)"
        ).fetchall()
    ]
    con.execute("CREATE TEMP TABLE obj_phylum (l0_object_sha256 BLOB, phylum VARCHAR)")
    con.executemany("INSERT INTO obj_phylum VALUES (?, ?)", obj_to_phylum)
    con.execute(f"""
        CREATE TABLE prot_phylum AS
        WITH contig_ph AS (
            SELECT DISTINCT o.contig_uid, coalesce(op.phylum, 'unknown') AS phylum
            FROM read_parquet({_sql_list(occ)}, hive_partitioning = false) o
            LEFT JOIN obj_phylum op USING (l0_object_sha256)
        ),
        fp AS (
            SELECT DISTINCT contig_uid, protein_uid FROM read_parquet({_sql_list(fp)}, hive_partitioning = false)
        )
        SELECT fp.protein_uid, min(coalesce(cp.phylum, 'unknown')) AS phylum,
               count(DISTINCT coalesce(cp.phylum, 'unknown')) AS n_phyla
        FROM fp LEFT JOIN contig_ph cp USING (contig_uid)
        GROUP BY fp.protein_uid
    """)


def _eval_duckdb(
    con: duckdb.DuckDBPyConnection,
    cm_files: list[Path],
    aln_files: list[Path],
    tag: str,
    *,
    with_coords: bool = True,
) -> dict:
    """Evaluate NR90 clustering (`cm_files`, `aln_files`) against C5's full Pfam scan in DuckDB.

    Expects tables already created in `con`:
      - `c5_prot(uid BLOB, len BIGINT, copied BOOLEAN)`
      - `prot_phylum(protein_uid BLOB, phylum VARCHAR, n_phyla BIGINT)`
      - `all_prot_len(uid BLOB, len BIGINT)` (for release-level reps outside C5)
      - `dom_rows(p BLOB, acc VARCHAR, env_start INT, env_end INT)` (full scan rows for C5 + any earlier stages)
      - `dom_pair(p BLOB, acc VARCHAR)` (distinct (protein_uid, model_acc))
      - `dom_n(p BLOB, n BIGINT)` (count of distinct models per protein)
    """
    con.execute(f"""
        CREATE OR REPLACE TABLE cm_{tag} AS
        SELECT member_uid AS m, rep_uid AS r
        FROM read_parquet({_sql_list(cm_files)}, hive_partitioning = false)
        WHERE level = 'nr90'
    """)
    con.execute(f"""
        CREATE OR REPLACE TABLE aln_{tag} AS
        SELECT member_uid AS m, rep_uid AS r, qstart::BIGINT AS qs, qend::BIGINT AS qe,
               tstart::BIGINT AS ts, tend::BIGINT AS te, cigar, pident::DOUBLE AS pid
        FROM read_parquet({_sql_list(aln_files)}, hive_partitioning = false)
        WHERE level = 'nr90' AND member_uid <> rep_uid
    """)
    # Filter to C5 proteins, join representative lengths and alignment spans
    con.execute(f"""
        CREATE OR REPLACE TABLE mem_{tag} AS
        SELECT p.uid AS m, c.r, p.len AS lm, pr.len AS lr, p.copied,
               coalesce(pp.phylum, 'unknown') AS phylum,
               (p.uid = c.r) AS is_rep,
               (p.uid <> c.r AND abs(p.len - pr.len) > 0.2 * pr.len) AS len_mis,
               a.ts, a.te, (a.m IS NOT NULL) AS has_aln, a.pid
        FROM c5_prot p
        JOIN cm_{tag} c ON c.m = p.uid
        JOIN all_prot_len pr ON pr.uid = c.r
        LEFT JOIN prot_phylum pp ON pp.protein_uid = p.uid
        LEFT JOIN aln_{tag} a ON a.m = p.uid AND a.r = c.r
    """)
    # f-coverage per representative domain row for non-rep, non-length-mismatch members
    con.execute(f"""
        CREATE OR REPLACE TABLE cov_{tag} AS
        SELECT m.m, m.r, d.acc, d.env_start, d.env_end, m.ts, m.te,
               CASE WHEN m.has_aln AND d.env_end > d.env_start THEN
                    greatest(0, least(d.env_end, m.te) - greatest(d.env_start, m.ts))::DOUBLE
                      / (d.env_end - d.env_start)
                    ELSE 0.0 END AS f
        FROM mem_{tag} m
        JOIN dom_rows d ON d.p = m.r
        WHERE NOT m.is_rep AND NOT m.len_mis
    """)
    con.execute(f"""
        CREATE OR REPLACE TABLE unmapped_{tag} AS
        SELECT m FROM mem_{tag} WHERE NOT is_rep AND NOT len_mis AND NOT has_aln
        UNION
        SELECT DISTINCT m FROM cov_{tag} WHERE f < 0.8
    """)
    con.execute(f"""
        CREATE OR REPLACE TABLE role_{tag} AS
        SELECT m.m, m.r, m.lm, m.lr, m.copied, m.phylum,
               CASE
                 WHEN m.is_rep THEN 'rep'
                 WHEN m.len_mis THEN 'length'
                 WHEN u.m IS NOT NULL THEN 'unmapped'
                 ELSE 'inferred'
               END AS role
        FROM mem_{tag} m
        LEFT JOIN unmapped_{tag} u ON u.m = m.m
    """)
    # Pair intersection for inferred members
    con.execute(f"""
        CREATE OR REPLACE TABLE inter_{tag} AS
        SELECT r.m, count(*) AS ni
        FROM role_{tag} r
        JOIN dom_pair dm ON dm.p = r.m
        JOIN dom_pair dr ON dr.p = r.r AND dr.acc = dm.acc
        WHERE r.role = 'inferred'
        GROUP BY r.m
    """)

    def _agg(where_sql: str) -> dict:
        row = con.execute(f"""
            SELECT
                count(*),
                coalesce(sum(r.lm), 0),
                count(*) FILTER (WHERE r.role = 'rep'),
                count(*) FILTER (WHERE r.role = 'length'),
                count(*) FILTER (WHERE r.role = 'unmapped'),
                count(*) FILTER (WHERE r.role = 'inferred'),
                count(*) FILTER (WHERE r.role <> 'inferred'),
                coalesce(sum(r.lm) FILTER (WHERE r.role <> 'inferred'), 0),
                count(*) FILTER (WHERE r.copied),
                coalesce(sum(r.lm) FILTER (WHERE r.copied), 0),
                count(*) FILTER (WHERE r.role <> 'inferred' AND NOT r.copied),
                coalesce(sum(r.lm) FILTER (WHERE r.role <> 'inferred' AND NOT r.copied), 0),
                coalesce(sum(coalesce(dm.n, 0)), 0),
                coalesce(sum(coalesce(dm.n, 0)) FILTER (WHERE r.role <> 'inferred'), 0)
                  + coalesce(sum(coalesce(ix.ni, 0)) FILTER (WHERE r.role = 'inferred'), 0),
                coalesce(sum(coalesce(dr.n, 0) - coalesce(ix.ni, 0)) FILTER (WHERE r.role = 'inferred'), 0),
                coalesce(sum(coalesce(dm.n, 0) - coalesce(ix.ni, 0)) FILTER (WHERE r.role = 'inferred'), 0)
            FROM role_{tag} r
            LEFT JOIN dom_n dm ON dm.p = r.m
            LEFT JOIN dom_n dr ON dr.p = r.r
            LEFT JOIN inter_{tag} ix ON ix.m = r.m
            WHERE {where_sql}
        """).fetchone()
        (n, res, reps, len_sc, unm_sc, inf, sc_p, sc_r, cp_p, cp_r, sc_p_b, sc_r_b,
         true_p, right_p, false_p, missed_p) = (int(x) for x in row)
        reported = right_p + false_p
        wrong_r = (false_p / reported) if reported else 0.0
        missed_r = (missed_p / true_p) if true_p else 0.0
        return {
            "proteins": n,
            "total_residues": res,
            "reps": reps,
            "length_scanned": len_sc,
            "unmapped_scanned": unm_sc,
            "inferred_members": inf,
            "scanned_proteins": sc_p,
            "scanned_residues": sc_r,
            "scanned_protein_fraction": round(sc_p / n, 6) if n else 0.0,
            "scanned_residue_fraction": round(sc_r / res, 6) if res else 0.0,
            "copied_proteins": cp_p,
            "copied_residues": cp_r,
            "scanned_proteins_after_b47b": sc_p_b,
            "scanned_residues_after_b47b": sc_r_b,
            "scanned_residue_fraction_after_b47b": round(sc_r_b / res, 6) if res else 0.0,
            "true_pairs": true_p,
            "reported_pairs": reported,
            "pairs_right": right_p,
            "pairs_false": false_p,
            "pairs_missed": missed_p,
            "wrong_rate": round(wrong_r, 6),
            "missed_rate": round(missed_r, 6),
            "stop_triggered": (wrong_r > STOP_THRESHOLD) or (missed_r > STOP_THRESHOLD),
        }

    overall = _agg("TRUE")
    phyla = [r[0] for r in con.execute(f"SELECT DISTINCT phylum FROM role_{tag} ORDER BY phylum").fetchall()]
    by_phylum = {ph: _agg(f"r.phylum = '{ph.replace(chr(39), chr(39)*2)}'") for ph in phyla}

    f_rows = con.execute(f"""
        SELECT CASE WHEN f >= 1 THEN '1.0' WHEN f >= 0.8 THEN '[0.8,1)' WHEN f >= 0.5 THEN '[0.5,0.8)'
                    WHEN f > 0 THEN '(0,0.5)' ELSE '0' END AS bin,
               count(*) AS rows, count(z.p) AS true_pair
        FROM cov_{tag} c
        LEFT JOIN dom_pair z ON z.p = c.m AND z.acc = c.acc
        GROUP BY bin ORDER BY bin
    """).fetchall()
    f_bins = {b: {"rows": int(n), "true_pair": int(tp)} for b, n, tp in f_rows}

    coord_summary = None
    if with_coords:
        # Sample up to 200,000 inferred domain pairs where member and rep both have the model
        sample = con.execute(f"""
            SELECT a.qs, a.qe, a.ts, a.te, a.cigar, dr.env_start, dr.env_end, dm.env_start, dm.env_end
            FROM role_{tag} r
            JOIN aln_{tag} a ON a.m = r.m AND a.r = r.r
            JOIN dom_rows dr ON dr.p = r.r
            JOIN dom_rows dm ON dm.p = r.m AND dm.acc = dr.acc
            WHERE r.role = 'inferred'
            ORDER BY r.m, dr.acc, dr.env_start
            LIMIT 200000
        """).fetchall()
        ds_list, de_list = [], []
        for qs, qe, ts, te, cig, r_es, r_ee, m_es, m_ee in sample:
            mapped = map_envelope(bytes(cig), int(qs), int(qe), int(ts), int(te), int(r_es), int(r_ee))
            if mapped is not None:
                ds_list.append(abs(mapped[0] - int(m_es)))
                de_list.append(abs(mapped[1] - int(m_ee)))
        if ds_list:
            ds_arr = np.asarray(ds_list, dtype=float)
            de_arr = np.asarray(de_list, dtype=float)
            both = np.maximum(ds_arr, de_arr)
            coord_summary = {
                "n_pairs": len(ds_list),
                "ds_median": float(np.quantile(ds_arr, 0.5)),
                "ds_p95": float(np.quantile(ds_arr, 0.95)),
                "de_median": float(np.quantile(de_arr, 0.5)),
                "de_p95": float(np.quantile(de_arr, 0.95)),
                "both_within_5": float((both <= 5).mean()),
                "both_within_10": float((both <= 10).mean()),
            }

    return {
        "overall": overall,
        "by_phylum": by_phylum,
        "f_bins": f_bins,
        "coordinates": coord_summary,
    }


# --- the breakdown (card C-47, after DEC-35 (6) stopped C5 on missed pairs, 2026-11-21) --------------------------------
# Measurements of where the missed and wrong pairs come from and of alternative rules, for the PI. Nothing here is a rule
# of the build: B-47a to B-47e are held.

MARGIN_EDGES = (0, 1, 2, 5, 10, 25)
SWEEP_LENGTH = (0.2, 0.15, 0.1, 0.05)
SWEEP_IDENTITY = (None, 0.95, 0.97, 0.99)
SWEEP_F = (0.8, 0.9, 1.0)
CHEAPEST_LEVELS = (0.01, 0.005)


def parse_ga(lines) -> dict[str, tuple[float, float]]:
    """Each model's gathering thresholds (sequence, domain) from a Pfam HMM file's `ACC` and `GA` lines, keyed by the
    accession without its version (`PF00001.26` -> `PF00001`)."""
    out: dict[str, tuple[float, float]] = {}
    acc = None
    for line in lines:
        if line.startswith("ACC "):
            acc = line.split()[1].split(".")[0]
        elif line.startswith("GA ") and acc is not None:
            f = line.split()
            out[acc] = (float(f[1]), float(f[2].rstrip(";")))
        elif line.startswith("//"):
            acc = None
    return out


def margin_bin(x: float | None) -> str:
    """The bin of a bit score's margin over its gathering threshold."""
    if x is None:
        return "no_score"
    if x < 0:
        return "<0"
    for lo, hi in zip(MARGIN_EDGES, MARGIN_EDGES[1:]):
        if x < hi:
            return f"[{lo},{hi})"
    return f"[{MARGIN_EDGES[-1]},inf)"


def _margin_sql(x: str) -> str:
    """margin_bin as a SQL expression over `x` (the same edges)."""
    arms = " ".join(f"WHEN {x} < {hi} THEN '[{lo},{hi})'" for lo, hi in zip(MARGIN_EDGES, MARGIN_EDGES[1:]))
    return f"CASE WHEN {x} IS NULL THEN 'no_score' WHEN {x} < 0 THEN '<0' {arms} ELSE '[{MARGIN_EDGES[-1]},inf)' END"


def match_instance(mapped: tuple[int, int], member_envs) -> tuple[int, int] | None:
    """The member's own envelope of the same model that overlaps the mapped envelope most (the first by start on a
    tie), or None when none overlaps: a model can occur more than once in a protein, and only an overlapping instance
    is the one the inference describes."""
    best, best_ov = None, 0
    for s, e in sorted(member_envs):
        ov = min(e, mapped[1]) - max(s, mapped[0])
        if ov > best_ov:
            best, best_ov = (s, e), ov
    return best


def _scanned_sql(length_cut: float, identity_min: float | None, f_cut: float) -> str:
    """A member is scanned under the rule (DEC-35 with these parameters): it is a representative, its length differs by
    more than `length_cut` of its representative's, it has no alignment, a representative domain lies less than `f_cut`
    inside the alignment, or (when given) its identity to the representative is under `identity_min`."""
    ident = "" if identity_min is None else f" OR coalesce(pid, 0) < {identity_min}"
    return f"(is_rep OR abs(lm - lr) > {length_cut} * lr OR NOT has_aln OR min_f < {f_cut}{ident})"


def _rule_row(con, length_cut: float, identity_min: float | None, f_cut: float) -> dict:
    sc = _scanned_sql(length_cut, identity_min, f_cut)
    (n, res, sp, sr, spb, srb, true_p, right, false_p, missed) = (int(v or 0) for v in con.execute(f"""
        SELECT count(*), sum(lm),
               count(*) FILTER (WHERE sc), coalesce(sum(lm) FILTER (WHERE sc), 0),
               count(*) FILTER (WHERE sc AND NOT copied), coalesce(sum(lm) FILTER (WHERE sc AND NOT copied), 0),
               sum(true_n), sum(CASE WHEN sc THEN true_n ELSE ni END),
               sum(CASE WHEN sc THEN 0 ELSE rep_n - ni END), sum(CASE WHEN sc THEN 0 ELSE true_n - ni END)
        FROM (SELECT *, {sc} AS sc FROM bd_feat)
    """).fetchone())
    reported = right + false_p
    return {
        "length_cut": length_cut, "identity_min": identity_min, "f_cut": f_cut,
        "scanned_proteins": sp, "scanned_residues": sr,
        "scanned_protein_fraction": round(sp / n, 6) if n else 0.0,
        "scanned_residue_fraction": round(sr / res, 6) if res else 0.0,
        "scanned_proteins_after_b47b": spb, "scanned_residues_after_b47b": srb,
        "scanned_residue_fraction_after_b47b": round(srb / res, 6) if res else 0.0,
        "true_pairs": true_p, "reported_pairs": reported, "pairs_right": right, "pairs_false": false_p,
        "pairs_missed": missed,
        "wrong_rate": (false_p / reported) if reported else 0.0,
        "missed_rate": (missed / true_p) if true_p else 0.0,
    }


def _bins(con, sql: str) -> dict[str, int]:
    return {b: int(n) for b, n in con.execute(sql).fetchall()}


def _breakdown(con, tag: str, dom_files: list[Path], ga: dict[str, tuple[float, float]], overall: dict) -> dict:
    """The breakdown of the evaluation's tables for `tag` (after `_eval_duckdb`): causes, margins over GA, identity and
    length bins, models, the sweep of rules and the coordinate block with instances matched by position."""
    con.execute("CREATE OR REPLACE TEMP TABLE ga (acc_base VARCHAR, ga_seq DOUBLE, ga_dom DOUBLE)")
    con.executemany("INSERT INTO ga VALUES (?, ?, ?)", [(k, s, d) for k, (s, d) in sorted(ga.items())])
    con.execute(f"""
        CREATE OR REPLACE TABLE dom_score AS
        SELECT protein_uid AS p, model_acc AS acc, max(seq_bitscore)::DOUBLE AS s, max(dom_bitscore)::DOUBLE AS d
        FROM read_parquet({_sql_list(dom_files)}, hive_partitioning = false) GROUP BY 1, 2
    """)
    con.execute(f"""
        CREATE OR REPLACE TABLE bd_inter AS
        SELECT x.m, count(*) AS ni FROM mem_{tag} x
        JOIN dom_pair dm ON dm.p = x.m JOIN dom_pair dr ON dr.p = x.r AND dr.acc = dm.acc
        WHERE NOT x.is_rep GROUP BY x.m
    """)
    con.execute(f"CREATE OR REPLACE TABLE bd_minf AS SELECT m, min(f) AS min_f FROM cov_{tag} GROUP BY m")
    con.execute(f"""
        CREATE OR REPLACE TABLE bd_feat AS
        SELECT x.m, x.r, x.lm, x.lr, x.copied, x.is_rep, x.has_aln, x.pid, coalesce(f.min_f, 1.0) AS min_f,
               coalesce(dm.n, 0) AS true_n, coalesce(dr.n, 0) AS rep_n, coalesce(i.ni, 0) AS ni
        FROM mem_{tag} x
        LEFT JOIN bd_minf f ON f.m = x.m
        LEFT JOIN dom_n dm ON dm.p = x.m
        LEFT JOIN dom_n dr ON dr.p = x.r
        LEFT JOIN bd_inter i ON i.m = x.m
    """)

    # the rule as built (DEC-35: 20%, f 0.8, no identity cut); it must reproduce the evaluation's totals
    base = _rule_row(con, 0.2, None, 0.8)
    keys = ("scanned_proteins", "scanned_residues", "scanned_residues_after_b47b", "true_pairs", "pairs_right",
            "pairs_false", "pairs_missed")
    consistent = all(base[k] == overall[k] for k in keys)
    base["total_residues"] = int(overall["total_residues"])
    con.execute(f"CREATE OR REPLACE TABLE bd_inf AS SELECT * FROM bd_feat WHERE NOT {_scanned_sql(0.2, None, 0.8)}")
    con.execute("""
        CREATE OR REPLACE TABLE bd_missed AS
        SELECT x.m, x.r, x.rep_n, dm.acc FROM bd_inf x JOIN dom_pair dm ON dm.p = x.m
        WHERE NOT EXISTS (SELECT 1 FROM dom_pair dr WHERE dr.p = x.r AND dr.acc = dm.acc)
    """)
    con.execute("""
        CREATE OR REPLACE TABLE bd_wrong AS
        SELECT x.m, x.r, dr.acc FROM bd_inf x JOIN dom_pair dr ON dr.p = x.r
        WHERE NOT EXISTS (SELECT 1 FROM dom_pair dm WHERE dm.p = x.m AND dm.acc = dr.acc)
    """)
    n_missed = con.execute("SELECT count(*) FROM bd_missed").fetchone()[0]
    n_wrong = con.execute("SELECT count(*) FROM bd_wrong").fetchone()[0]
    consistent = consistent and n_missed == base["pairs_missed"] and n_wrong == base["pairs_false"]
    not_mapped = con.execute(f"""
        SELECT count(*) FROM cov_{tag} c SEMI JOIN bd_inf x ON x.m = c.m WHERE c.f < 0.8
    """).fetchone()[0]
    missed_by_cause = {
        "rep_has_no_row_for_the_model": int(n_missed),
        "rep_row_not_mapped": int(not_mapped),  # 0 by construction: such a member is scanned (DEC-35 (4))
        "of_which_rep_has_no_domain_at_all": int(con.execute(
            "SELECT count(*) FROM bd_missed WHERE rep_n = 0").fetchone()[0]),
    }
    wrong_by_cause = {"member_scan_lacks_the_model": int(n_wrong)}

    def margins(pairs_sql: str) -> dict:
        q = (f"FROM ({pairs_sql}) t LEFT JOIN dom_score s ON s.p = t.p AND s.acc = t.acc "
             f"LEFT JOIN ga g ON g.acc_base = split_part(t.acc, '.', 1)")
        return {"seq": _bins(con, f"SELECT {_margin_sql('s.s - g.ga_seq')} AS b, count(*) {q} GROUP BY b ORDER BY b"),
                "dom": _bins(con, f"SELECT {_margin_sql('s.d - g.ga_dom')} AS b, count(*) {q} GROUP BY b ORDER BY b")}

    margin_doc = {
        "missed": margins("SELECT m AS p, acc FROM bd_missed"),  # the member has the hit
        "wrong": margins("SELECT r AS p, acc FROM bd_wrong"),  # the representative has the hit
        "right_inferred": margins(  # for comparison: the inferred pairs that are right, the member's own hit
            "SELECT x.m AS p, dm.acc FROM bd_inf x JOIN dom_pair dm ON dm.p = x.m "
            "JOIN dom_pair dr ON dr.p = x.r AND dr.acc = dm.acc"),
    }

    def by(expr: str) -> dict:
        rows = con.execute(f"""
            SELECT {expr} AS b, count(*), sum(true_n), sum(rep_n), sum(rep_n - ni), sum(true_n - ni)
            FROM bd_inf GROUP BY b ORDER BY b
        """).fetchall()
        return {b: {"members": int(n), "true_pairs": int(t), "reported_pairs": int(rp), "pairs_false": int(f),
                    "pairs_missed": int(mi)} for b, n, t, rp, f, mi in rows}

    by_identity = by("CASE WHEN pid IS NULL THEN 'none' WHEN pid < 0.90 THEN '<0.90' WHEN pid < 0.95 THEN '[0.90,0.95)' "
                     "WHEN pid < 0.97 THEN '[0.95,0.97)' WHEN pid < 0.99 THEN '[0.97,0.99)' "
                     "WHEN pid < 1.0 THEN '[0.99,1)' ELSE '1.0' END")
    by_length = by("CASE WHEN abs(lm - lr) = 0 THEN '0' WHEN abs(lm - lr) <= 0.05 * lr THEN '(0,0.05]' "
                   "WHEN abs(lm - lr) <= 0.10 * lr THEN '(0.05,0.10]' WHEN abs(lm - lr) <= 0.15 * lr THEN '(0.10,0.15]' "
                   "WHEN abs(lm - lr) <= 0.20 * lr THEN '(0.15,0.20]' ELSE '>0.20' END")

    def top(table: str, total: int) -> list[dict]:
        rows = con.execute(f"""
            SELECT split_part(acc, '.', 1) AS model, count(*) AS n FROM {table} GROUP BY model ORDER BY n DESC, model
            LIMIT 20
        """).fetchall()
        return [{"model": mo, "pairs": int(n), "share": round(n / total, 6) if total else 0.0} for mo, n in rows]

    sweep = [_rule_row(con, lc, im, fc) for lc in SWEEP_LENGTH for im in SWEEP_IDENTITY for fc in SWEEP_F]

    def cheapest(level: float) -> dict | None:
        ok = [r for r in sweep if r["wrong_rate"] <= level and r["missed_rate"] <= level]
        key = lambda r: (r["scanned_residues"], -r["length_cut"],  # noqa: E731
                         -1.0 if r["identity_min"] is None else r["identity_min"], r["f_cut"])
        return min(ok, key=key) if ok else None

    # the coordinate block with instances matched by position (the evaluation's block joined every member row of the
    # same model, so a model occurring twice paired mismatched instances)
    rep_rows = int(con.execute("SELECT count(*) FROM bd_inf x JOIN dom_rows dr ON dr.p = x.r").fetchone()[0])
    sample = con.execute(f"""
        SELECT a.qs, a.qe, a.ts, a.te, a.cigar, dr.env_start, dr.env_end,
               list(struct_pack(s := dm.env_start, e := dm.env_end)) FILTER (WHERE dm.p IS NOT NULL) AS mine
        FROM bd_inf x
        JOIN aln_{tag} a ON a.m = x.m AND a.r = x.r
        JOIN dom_rows dr ON dr.p = x.r
        LEFT JOIN dom_rows dm ON dm.p = x.m AND dm.acc = dr.acc
        GROUP BY x.m, a.qs, a.qe, a.ts, a.te, a.cigar, dr.acc, dr.env_start, dr.env_end
        ORDER BY x.m, dr.acc, dr.env_start
        LIMIT 200000
    """).fetchall()
    ds, de, no_inst, unmappable = [], [], 0, 0
    for qs, qe, ts, te, cig, r_es, r_ee, mine in sample:
        mapped = map_envelope(bytes(cig), int(qs), int(qe), int(ts), int(te), int(r_es), int(r_ee))
        if mapped is None:
            unmappable += 1
            continue
        hit = match_instance(mapped, [(int(x["s"]), int(x["e"])) for x in (mine or [])])
        if hit is None:
            no_inst += 1
            continue
        ds.append(abs(mapped[0] - hit[0]))
        de.append(abs(mapped[1] - hit[1]))
    coords: dict = {"rep_rows": rep_rows, "sampled": len(sample), "matched": len(ds), "no_member_instance": no_inst,
                    "unmappable": unmappable}
    if ds:
        a1, a2 = np.asarray(ds, dtype=float), np.asarray(de, dtype=float)
        both = np.maximum(a1, a2)
        coords.update({"ds_median": float(np.quantile(a1, 0.5)), "ds_p95": float(np.quantile(a1, 0.95)),
                       "de_median": float(np.quantile(a2, 0.5)), "de_p95": float(np.quantile(a2, 0.95)),
                       "both_within_5": float((both <= 5).mean()), "both_within_10": float((both <= 10).mean())})

    return {
        "base": base,
        "consistent_with_eval": bool(consistent),
        "missed_by_cause": missed_by_cause,
        "wrong_by_cause": wrong_by_cause,
        "margins": margin_doc,
        "by_identity": by_identity,
        "by_length": by_length,
        "by_model": {"missed": top("bd_missed", int(n_missed)), "wrong": top("bd_wrong", int(n_wrong))},
        "sweep": sweep,
        "cheapest": {str(lv): cheapest(lv) for lv in CHEAPEST_LEVELS},
        "coordinates": coords,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("eval", "breakdown"))
    ap.add_argument("--hmm", type=Path, default=None, help="breakdown: Pfam-A.hmm[.gz], for each model's GA")
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--stage", default="C5")
    ap.add_argument("--earlier", default="S1,P1,P2,P3,P4,P5,P6,H1,PL1")
    ap.add_argument("--stage-cm", required=True, nargs="+", type=Path)
    ap.add_argument("--stage-aln", required=True, nargs="+", type=Path)
    ap.add_argument("--release-cm", nargs="+", type=Path, default=None)
    ap.add_argument("--release-aln", nargs="+", type=Path, default=None)
    ap.add_argument("--selection", required=True, type=Path)
    ap.add_argument("--objects", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--tmp", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--memory", default="12GB")
    a = ap.parse_args(argv)

    t0 = time.monotonic()
    from biobox.l1.genecall import GENE_CALL_VERSION as GCV

    a.tmp.mkdir(parents=True, exist_ok=True)
    db_path = a.tmp / "c47_eval.duckdb"
    db_path.unlink(missing_ok=True)
    con = duckdb.connect(
        str(db_path),
        config={
            "threads": a.threads,
            "memory_limit": a.memory,
            "temp_directory": str(a.tmp / "spill"),
            "preserve_insertion_order": False,
        },
    )

    root, st = a.root, a.stage
    earlier = [s for s in a.earlier.split(",") if s]

    # Load C5 proteins + B-47b copied status from C5's scan plan
    c5_ps = _record_files(root, st, "genecall", f"l0seq/protein_seq/gene_call_version={GCV}/stage={st}/")
    plan_p = root / f"l1/domain_scan/library={LIB}/stage={st}/plan.json"
    plan = json.loads(plan_p.read_text())
    src_ps = [root / x["path"] for s in (plan.get("copy") or {}).get("sources", []) for x in s["protein_files"]]
    if src_ps:
        con.execute(f"""
            CREATE TABLE c5_prot AS
            WITH src AS (SELECT DISTINCT protein_uid FROM read_parquet({_sql_list(src_ps)}, hive_partitioning = false))
            SELECT p.protein_uid AS uid, p.length::BIGINT AS len, (s.protein_uid IS NOT NULL) AS copied
            FROM read_parquet({_sql_list(c5_ps)}, hive_partitioning = false) p
            LEFT JOIN src s ON s.protein_uid = p.protein_uid
        """)
    else:
        con.execute(f"""
            CREATE TABLE c5_prot AS
            SELECT protein_uid AS uid, length::BIGINT AS len, FALSE AS copied
            FROM read_parquet({_sql_list(c5_ps)}, hive_partitioning = false)
        """)

    all_stages = earlier + [st] if a.release_cm else [st]
    all_ps = [
        p for s in all_stages
        for p in _record_files(root, s, "genecall", f"l0seq/protein_seq/gene_call_version={GCV}/stage={s}/")
    ]
    con.execute(f"""
        CREATE TABLE all_prot_len AS
        SELECT protein_uid AS uid, max(length)::BIGINT AS len
        FROM read_parquet({_sql_list(all_ps)}, hive_partitioning = false)
        GROUP BY protein_uid
    """)

    # Load full-scan domains for C5 (and earlier stages if release-level counterfactual is run)
    all_dom = [
        p for s in all_stages
        for p in _record_files(root, s, f"domains-{LIB}", f"l1/domains/library={LIB}/stage={s}/")
    ]
    con.execute(f"""
        CREATE TABLE dom_rows AS
        SELECT DISTINCT protein_uid AS p, model_acc AS acc, env_start::INT AS env_start, env_end::INT AS env_end
        FROM read_parquet({_sql_list(all_dom)}, hive_partitioning = false)
    """)
    con.execute("CREATE TABLE dom_pair AS SELECT DISTINCT p, acc FROM dom_rows")
    con.execute("CREATE TABLE dom_n AS SELECT p, count(*) AS n FROM dom_pair GROUP BY p")

    _phylum_map_sql(con, root, st, a.selection, a.objects, GCV)

    if a.mode == "breakdown":
        if a.hmm is None:
            raise SystemExit("breakdown needs --hmm (the Pfam-A.hmm[.gz] whose GA lines the scan used)")
        stage_res = _eval_duckdb(con, list(a.stage_cm), list(a.stage_aln), "stage", with_coords=False)
        opener = gzip.open if str(a.hmm).endswith(".gz") else open
        with opener(a.hmm, "rt") as fh:
            ga = parse_ga(fh)
        doc = {"stage": st, "library": LIB, "hmm": str(a.hmm), "ga_models": len(ga),
               **_breakdown(con, "stage", all_dom, ga, stage_res["overall"])}
        doc["elapsed_s"] = round(time.monotonic() - t0, 1)
        con.close()
        db_path.unlink(missing_ok=True)
        a.out.parent.mkdir(parents=True, exist_ok=True)
        tmp_out = a.out.with_name(f".{a.out.name}.tmp")
        tmp_out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
        tmp_out.replace(a.out)
        print(f"CONSISTENT_WITH_EVAL: {doc['consistent_with_eval']}", flush=True)
        return 0 if doc["consistent_with_eval"] else 3

    stage_res = _eval_duckdb(con, list(a.stage_cm), list(a.stage_aln), "stage", with_coords=True)
    doc: dict = {
        "stage": st,
        "library": LIB,
        "stage_local": stage_res,
        "dec35_6_stop_triggered": stage_res["overall"]["stop_triggered"],
    }

    if a.release_cm and a.release_aln:
        rel_res = _eval_duckdb(con, list(a.release_cm), list(a.release_aln), "rel", with_coords=False)
        gap_pp = round(
            (stage_res["overall"]["scanned_residue_fraction"] - rel_res["overall"]["scanned_residue_fraction"]) * 100.0,
            4,
        )
        gap_b47b_pp = round(
            (
                stage_res["overall"]["scanned_residue_fraction_after_b47b"]
                - rel_res["overall"]["scanned_residue_fraction_after_b47b"]
            ) * 100.0,
            4,
        )
        doc["release_level"] = rel_res
        doc["dec35_7_comparison"] = {
            "stage_local_scanned_residue_pct": round(stage_res["overall"]["scanned_residue_fraction"] * 100.0, 4),
            "release_level_scanned_residue_pct": round(rel_res["overall"]["scanned_residue_fraction"] * 100.0, 4),
            "gap_percentage_points": gap_pp,
            "stage_local_after_b47b_pct": round(stage_res["overall"]["scanned_residue_fraction_after_b47b"] * 100.0, 4),
            "release_level_after_b47b_pct": round(rel_res["overall"]["scanned_residue_fraction_after_b47b"] * 100.0, 4),
            "gap_after_b47b_pp": gap_b47b_pp,
            "revisit_triggered_gt_10pp": gap_pp > DEC35_7_GAP_PP,
        }

    doc["elapsed_s"] = round(time.monotonic() - t0, 1)
    con.close()
    db_path.unlink(missing_ok=True)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    tmp_out = a.out.with_name(f".{a.out.name}.tmp")
    tmp_out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    tmp_out.replace(a.out)
    print(f"STOP_TRIGGERED: {doc['dec35_6_stop_triggered']}", flush=True)
    print(json.dumps(doc["stage_local"]["overall"], indent=1, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

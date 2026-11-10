"""The query set of the 10 Gbp rebuild through one release's catalog (wall time and what each read). Read only.

    python r10_queries.py ROOT RELEASE OUT_JSON

Run in the image, in a capped scope. Opens `ROOT/catalog/release=RELEASE.duckdb` with
`catalog.open_catalog` (4 DuckDB threads, 16 GB) and times, with `catalog.profile_query` where it
applies (row groups scanned of the total each scan could read; the point lookups of `bb_slot_tokens`
and `bb_token_of` are scans of their own):

  counts      rows of this release: cluster members and clusters per level, tokens per level,
              postings, tape tokens, contig_dim (plain counts, not profiled)
  seq_get     1,000 random genes +-500 bp, picked as the ladder picked them (`ORDER BY
              hash(feature_uid), feature_uid LIMIT 1000` over feature_protein joined to
              contig_dim_now()), through `seq_get`; then the same again (warm)
  tapes_of    tapes of one contig (the first pick)
  postings    the largest NR50 cluster (most slots in token_dict): `token_rows(tok)` (the B-31 macro)
              and the ladder's literal filter on token_postings
  partners    NR50 partners within 3 slots of the largest NR50 cluster and of the largest prokaryote
              NR50 cluster (most postings on contigs of genome_input kind 'prokaryote'; skipped when
              the release has none): `partners_of(tok, 3)` (card B-31: partner tokens from
              tape_tokens of the anchors' contigs), and the ladder's query (partners from
              token_postings filtered by contig_uid) for comparison on the same store
  membership  cluster_member rows of one protein (the representative of the last of those clusters),
              by member_uid
  label       contig_label of one contig (the first pick) and the reverse, accession -> contig_uid

Writes one JSON object to OUT_JSON and prints it.
"""

import json
import sys
import time
from pathlib import Path

import pyarrow as pa

from biobox import catalog


def slim(p: dict) -> dict:
    return {"scans": [{k: s.get(k) for k in ("row_groups_scanned", "row_groups_total", "files_read", "rows_scanned",
                                                "lookup")} | {"table": _table(s)} for s in p.get("scans", [])],
            "bytes_read": p.get("bytes_read"), "bgzf": p.get("bgzf"), "index": p.get("index")}


def _table(s: dict) -> str | None:
    fs = s.get("files") or []
    if not fs:
        return None
    parts = Path(fs[0]).parts
    for i, x in enumerate(parts):
        if x in ("l0", "l0seq", "l1", "l2", "l3", "dim", "ledger"):
            return parts[i + 1] if i + 1 < len(parts) else x
    return None


def timed(con, sql: str, params=None) -> dict:
    t = time.perf_counter()
    p = catalog.profile_query(con, sql, params)
    return {"wall_s": round(time.perf_counter() - t, 3), "result": [list(r) for r in p["result"]], **slim(p)}


def main(root: Path, rel: str, out: Path) -> dict:
    R: dict = {"release": rel, "root": str(root)}
    def connect():
        return catalog.open_catalog(root / "catalog" / f"release={rel}.duckdb", threads=4, memory_limit="16GB",
                                    temp_directory=out.parent / f"duck-q-{rel}")

    t = time.perf_counter()
    con = connect()
    R["open_catalog_s"] = round(time.perf_counter() - t, 3)

    t = time.perf_counter()
    R["counts"] = {
        "cluster_member": con.execute("SELECT level, count(*), count(DISTINCT cluster_uid) FROM cluster_member "
                                      "WHERE release = ? GROUP BY level ORDER BY level", [rel]).fetchall(),
        "token_dict": con.execute("SELECT level, count(*) FROM token_dict WHERE release = ? GROUP BY level "
                                  "ORDER BY level", [rel]).fetchall(),
        "token_postings": con.execute("SELECT count(*) FROM token_postings WHERE release = ?", [rel]).fetchone()[0],
        "tape_tokens": con.execute("SELECT count(*) FROM tape_tokens WHERE release = ?", [rel]).fetchone()[0],
        "contig_dim": con.execute("SELECT count(*) FROM contig_dim_now()").fetchone()[0],
        "contig_index_default": con.execute("SELECT count(*) FROM contig_index WHERE is_default_variant").fetchone()[0],
        "feature": con.execute("SELECT count(*) FROM feature").fetchone()[0],
        "protein_seq": con.execute("SELECT count(*) FROM protein_seq").fetchone()[0],
        "igr": con.execute("SELECT count(*) FROM igr").fetchone()[0],
    }
    R["counts_wall_s"] = round(time.perf_counter() - t, 2)

    t = time.perf_counter()
    rows = con.execute(
        "SELECT f.contig_uid, f.start, f.\"end\", f.strand, d.contig_len FROM feature_protein f "
        "JOIN contig_dim_now() d USING (contig_uid) ORDER BY hash(f.feature_uid), f.feature_uid LIMIT 1000").fetchall()
    R["pick_genes_s"] = round(time.perf_counter() - t, 2)
    picks = pa.table({"contig_uid": pa.array([r[0] for r in rows], pa.binary()),
                      "s": pa.array([max(0, r[1] - 500) for r in rows], pa.int64()),
                      "e": pa.array([min(r[4], r[2] + 500) for r in rows], pa.int64()),
                      "strand": pa.array([r[3] for r in rows], pa.int64())})
    con.register("picks", picks)
    q = "SELECT count(*), sum(length(seq_get(contig_uid, s, e, strand))) FROM picks"
    R["seq_get_1000"] = timed(con, q) | {
        "expected_bp": int(sum(e - s for s, e in zip(picks.column("s").to_pylist(), picks.column("e").to_pylist()))),
        "distinct_contigs": len(set(picks.column("contig_uid").to_pylist()))}
    R["seq_get_1000_warm"] = timed(con, q)

    cu = rows[0][0]
    R["tapes_of_one"] = {"contig_uid": cu.hex()} | timed(
        con, f"SELECT count(*), sum(n_slots), sum(len(nr50_tokens)) FROM tapes_of([unhex('{cu.hex()}')])")

    tok, key, nsl = con.execute("SELECT token_id, key, n_slots FROM token_dict WHERE release = ? AND level = 'nr50' "
                                "ORDER BY n_slots DESC, token_id LIMIT 1", [rel]).fetchone()
    R["postings_lookup_nr50"] = {
        "token_id": tok, "cluster_uid": key, "n_slots_dict": nsl,
        "token_rows": timed(con, f"SELECT count(*), sum(ordinal) FROM token_rows({tok}::UBIGINT)"),
        "literal": timed(con, f"SELECT count(*), sum(ordinal) FROM token_postings WHERE release = '{rel}' "
                              f"AND token_bucket_ok(token_bucket, {tok}::UBIGINT) AND token_id = {tok}::UBIGINT")}

    lo, hi = 2 << 61, 3 << 61

    def neighbours(ptok: int) -> dict:
        """NR50 partners within 3 slots of the slots of token `ptok`: the B-31 macro and the ladder's query. On a
        new connection, so the tape-token lookup's row-group cache holds nothing from an earlier query."""
        con = connect()
        d = {"partners_of": timed(
                 con, f"SELECT count(*), sum(n_loci), sum(n_pairs) FROM partners_of({ptok}::UBIGINT, 3) "
                      f"WHERE partner_token >= {lo}::UBIGINT AND partner_token < {hi}::UBIGINT"),
             "anchors": timed(con, f"SELECT count(*), count(DISTINCT contig_uid) FROM token_slots({ptok}::UBIGINT)")}
        t = time.perf_counter()
        d["top10"] = [[r[0], r[1], r[2]] for r in con.execute(
            f"SELECT d.key, p.n_loci, p.n_pairs FROM partners_of({ptok}::UBIGINT, 3) p JOIN token_dict d "
            f"ON d.token_id = p.partner_token AND d.release = '{rel}' WHERE p.partner_token >= {lo}::UBIGINT "
            f"AND p.partner_token < {hi}::UBIGINT ORDER BY p.n_loci DESC, d.key LIMIT 10").fetchall()]
        d["top10_wall_s"] = round(time.perf_counter() - t, 3)
        ladder = (f"WITH a AS (SELECT contig_uid, ordinal, n_slots, is_circular FROM token_postings "
                  f"  WHERE release = '{rel}' AND token_bucket_ok(token_bucket, {ptok}::UBIGINT) AND token_id = {ptok}::UBIGINT), "
                  f"p AS (SELECT token_id, contig_uid, ordinal FROM token_postings WHERE release = '{rel}' "
                  f"  AND token_id >= {lo}::UBIGINT AND token_id < {hi}::UBIGINT "
                  f"  AND contig_uid IN (SELECT contig_uid FROM a)), "
                  "nb AS (SELECT a.contig_uid, a.ordinal AS anchor, p.token_id AS partner FROM a JOIN p USING (contig_uid) "
                  "  WHERE p.ordinal <> a.ordinal AND (abs(p.ordinal::BIGINT - a.ordinal::BIGINT) <= 3 "
                  "  OR (a.is_circular AND a.n_slots - abs(p.ordinal::BIGINT - a.ordinal::BIGINT) <= 3))) "
                  "SELECT (SELECT count(*) FROM (SELECT DISTINCT contig_uid, ordinal FROM a)), "
                  "(SELECT count(DISTINCT contig_uid) FROM a), (SELECT count(DISTINCT (contig_uid, anchor)) FROM nb), "
                  "(SELECT count(*) FROM nb)")
        d["ladder_query"] = timed(con, ladder)
        return d

    R["partners_top_nr50"] = {"token_id": tok, "cluster_uid": key} | neighbours(tok)
    t = time.perf_counter()
    pk = con.execute(
        "SELECT p.token_id, d.key, count(*) AS n FROM token_postings p JOIN token_dict d "
        "ON d.token_id = p.token_id AND d.release = p.release "
        "WHERE p.release = ? AND d.level = 'nr50' AND p.contig_uid IN "
        "(SELECT contig_uid FROM genome_input WHERE kind = 'prokaryote') "
        "GROUP BY 1, 2 ORDER BY n DESC, p.token_id LIMIT 1", [rel]).fetchone()
    find_s = round(time.perf_counter() - t, 2)
    mkey = key
    if pk is not None:  # a release without prokaryote contigs has none
        ptok, mkey, pn = pk
        R["partners_prok_top_nr50"] = {"token_id": ptok, "cluster_uid": mkey, "prokaryote_postings": pn,
                                       "find_s": find_s} | neighbours(ptok)

    rep = con.execute("SELECT rep_uid FROM cluster_member WHERE release = ? AND level = 'nr50' AND cluster_uid = unhex(?) "
                      "LIMIT 1", [rel, mkey]).fetchone()
    if rep:
        R["membership_one_protein"] = {"member_uid": rep[0].hex()} | timed(
            con, f"SELECT level, hex(cluster_uid), hex(rep_uid) FROM cluster_member WHERE release = '{rel}' "
                 f"AND member_uid = unhex('{rep[0].hex()}') ORDER BY level")
    R["label_of_contig"] = {"contig_uid": cu.hex()} | timed(
        con, f"SELECT accession, n_occurrences FROM contig_label WHERE contig_uid = unhex('{cu.hex()}')")
    acc = R["label_of_contig"]["result"][0][0] if R["label_of_contig"]["result"] else None
    if acc is not None:
        R["contig_of_accession"] = {"accession": acc} | timed(
            con, "SELECT DISTINCT hex(contig_uid) FROM contig_occurrence WHERE accession = ?", [acc])
    out.write_text(json.dumps(R, indent=1, default=str))
    return R


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)
    print(json.dumps(main(Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])), indent=1, default=str))

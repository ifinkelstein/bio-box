"""Card C-48c (part 2), DEC-59: what the V2 and GV1 selections leave for later, written down so that it can be staged
without rediscovering it.

    python c48c_v2_gv1_left.py v2 HC_TABLE V2_INPUTS OUT_TSV_GZ
    python c48c_v2_gv1_left.py v4-putative FULL_TABLE OUT_TSV_GZ
    python c48c_v2_gv1_left.py gv1 GVMAGS_GENOMES_TSV OUT_TSV_GZ

v2: every row of IMG/VR v4.1's high-confidence table (HC_TABLE, IMGVR_all_Sequence_information-high_confidence.tsv as
`c48c_v2_gv1_look.py imgvr` copied it out, md5 checked) that V2 does not take (V2_INPUTS/records.tsv, source imgvr):
reason `votu_without_hq_member` (a vOTU none of whose rows is MIUViG High-quality) or `not_representative` (a member of a
vOTU V2 takes, other than its representative).
v4-putative: every row of IMG/VR v4.0's full table (FULL_TABLE, IMGVR_all_Sequence_information.tsv, as
`c48c_v2_gv1_est.py imgvr-full` copied it out, md5 checked) whose Confidence is Putative: in no v4.1 high-confidence
file, so not in V2's pool; reason `putative`.
gv1: every genome of GVMAGs V2's all_18kgenomes_fna.tar.gz (GVMAGS_GENOMES_TSV of `c48c_v2_gv1_look.py gvmags`):
GV1 takes none yet (its selection waits on a decision), reason `gv1_selection_undecided`.

Columns: sequence id (UViG or genome), vOTU (empty for GV1), quality tier (MIUViG for IMG/VR, empty for GV1), length (bp),
records (contigs), reason, archive, member. Rows sorted by reason then id; written with gzip mtime 0, level 6, through a
temporary name, so a second run gives the same bytes. DuckDB, one thread, 3 GB. Prints counts and bases by reason
and quality tier, rows and the sha256 of the file.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c48c_v2_gv1_look as look  # noqa: E402

HC_ZIP = "jgi_img_vr_v4.1_hc_20261118.zip"
HC_FNA = "Cus_VR/IMG_VR_2022-12-19_7.1/IMGVR_all_nucleotides-high_confidence.fna.gz"
FULL_ZIP = "jgi_img_vr_v4_full_20261118.zip"
FULL_FNA = "Cus_VR/IMG_VR_2022-12-19_7/IMGVR_all_nucleotides.fna.gz"
GV_TAR = "gvmags_v2_20261118.tar"
GV_MEMBER = "gvmags_v2/all_18kgenomes_fna.tar.gz"
HEAD = ["sequence_id", "votu", "quality", "length", "records", "reason", "archive", "member"]


def _table(con, path: Path) -> None:
    cols = "{" + ", ".join(f"'{c}': 'VARCHAR'" for c in look.IMGVR_COLS) + "}"
    con.execute(f"CREATE VIEW f AS SELECT * FROM read_csv('{path}', delim='\t', header=true, quote='', escape='', "
                f"columns={cols})")


def write(con, sql: str, out: Path) -> dict:
    """Rows of `sql` (the HEAD columns, already sorted) to `out` as deterministic gzip; counts by reason and tier."""
    tmp = out.with_name(f".{out.name}.tmp")
    rows = 0
    with open(tmp, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as gz:
        gz.write(("\t".join(HEAD) + "\n").encode())
        reader = con.sql(sql).to_arrow_reader(100_000)  # as genecall.stage_items reads a query
        for b in reader:
            cols = [b.column(i).to_pylist() for i in range(len(HEAD))]
            gz.write("".join("\t".join("" if v is None else str(v) for v in r) + "\n" for r in zip(*cols)).encode())
            rows += b.num_rows
    tmp.replace(out)
    by = con.execute(f"SELECT reason, quality, count(*), sum(length::BIGINT), sum(records::BIGINT) FROM ({sql}) "
                     "GROUP BY 1, 2 ORDER BY 1, 2").fetchall()
    return {"path": str(out), "rows": rows, "bytes": out.stat().st_size,
            "sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
            "by_reason_and_quality": [{"reason": r, "quality": q, "rows": n, "bases": int(s), "records": int(c)}
                                      for r, q, n, s, c in by]}


def main(argv: list[str]) -> int:
    con = duckdb.connect(config={"threads": 1, "memory_limit": "3GB", "preserve_insertion_order": True})
    if argv[:1] == ["v2"] and len(argv) == 4:
        _table(con, Path(argv[1]))
        con.execute(f"CREATE TABLE taken AS SELECT DISTINCT genome AS uvig, votu FROM read_csv('{Path(argv[2]) / 'records.tsv'}', "
                    "delim='\t', header=true, quote='', all_varchar=true) WHERE source = 'imgvr'")
        sql = (f"SELECT uvig, votu, miuvig AS quality, length::BIGINT AS length, 1 AS records, "
               f"CASE WHEN votu IN (SELECT votu FROM taken) THEN 'not_representative' ELSE 'votu_without_hq_member' END "
               f"AS reason, '{HC_ZIP}' AS archive, '{HC_FNA}' AS member FROM f WHERE uvig NOT IN (SELECT uvig FROM taken) "
               "ORDER BY reason, uvig")
        rep = write(con, sql, Path(argv[3]))
        rep["taken_uvigs"], = con.execute("SELECT count(*) FROM taken").fetchone()
        rep["hc_rows"], = con.execute("SELECT count(*) FROM f").fetchone()
        rep["votus_left_whole"], = con.execute(
            "SELECT count(DISTINCT votu) FROM f WHERE votu NOT IN (SELECT votu FROM taken)").fetchone()
    elif argv[:1] == ["v4-putative"] and len(argv) == 3:
        _table(con, Path(argv[1]))
        sql = (f"SELECT uvig, votu, miuvig AS quality, length::BIGINT AS length, 1 AS records, 'putative' AS reason, "
               f"'{FULL_ZIP}' AS archive, '{FULL_FNA}' AS member FROM f WHERE confidence = 'Putative' ORDER BY uvig")
        rep = write(con, sql, Path(argv[2]))
        rep["putative_votus"], = con.execute("SELECT count(DISTINCT votu) FROM f WHERE confidence = 'Putative'").fetchone()
        rep["putative_only_votus"], = con.execute(
            "SELECT count(*) FROM (SELECT votu FROM f GROUP BY votu HAVING bool_and(confidence = 'Putative'))").fetchone()
    elif argv[:1] == ["gv1"] and len(argv) == 3:
        con.execute(f"CREATE VIEW g AS SELECT * FROM read_csv('{argv[1]}', delim='\t', header=true, quote='', "
                    "all_varchar=true) WHERE \"set\" = 'all_18kgenomes_fna'")
        sql = (f"SELECT genome, '' AS votu, '' AS quality, bases::BIGINT AS length, records::BIGINT AS records, "
               f"'gv1_selection_undecided' AS reason, '{GV_TAR}' AS archive, '{GV_MEMBER}::' || member AS member "
               "FROM g ORDER BY genome")
        rep = write(con, sql, Path(argv[2]))
    else:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

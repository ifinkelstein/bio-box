"""Card C-48c (part 2) step 4: what stages V2 and GV1 share by `contig_uid` with the stages on disk, with V1's inputs
and with each other (read only).

    python c48c_v2_gv1_overlap.py --root ROOT --stages S1,P1,...,C5 --v1 V1_INPUTS --v2 V2_INPUTS \
        --gvmags GVMAGS_TAR --gv-genomes GVMAGS_GENOMES_TSV --out OUT

V2_INPUTS holds what `c48c_virus_inputs.py v2` wrote (records.tsv). GV1 has no inputs (its selection waits on a
decision, see the note), so its candidate contigs are read here from the archive: every FASTA record of every genome
of all_18kgenomes_fna.tar.gz in GVMAGS_TAR (streamed, one genome at a time), its contig uid (`ids.contig_uid`) and
length, written to OUT/gv1_candidate_contigs.tsv; GVMAGS_GENOMES_TSV (`c48c_v2_gv1_look.py gvmags`) says which genomes
are also in fna_files.tar.gz (the provider's 8,508-genome subset). The identity is `ids.contig_uid`, no strand or
rotation canonicalisation. Counted with DuckDB (one thread, 3 GB), for each set (V2; GV1 all 18,726 genomes; GV1's
8,508 named in fna_files, their contigs as the 18k copy has them):

- distinct contigs and bases; per stage of --stages, the contigs (and bases) its contig_index stores as default
  variants (it owns them, DEC-26) and the contigs its genome_input lists (DEC-27);
- the contigs V1's inputs list (records.tsv), and between V2 and GV1;
- what the stage would add: contigs and bases no earlier stage holds (for V2: the stages and V1; for GV1: those and V2),
  and its genomes all of whose contigs an earlier stage holds (second genomes of contigs an earlier stage keeps; DEC-26,
  DEC-27) or some of them.
Writes OUT/summary.json and prints it.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import tarfile
from pathlib import Path

import duckdb

from biobox import io as bbio
from biobox import ids
from biobox.l0seq import store

import c48c_v2_gv1_look as look  # noqa: E402  (beside this script; the run script sets the working directory)
import c48c_virus_inputs as vi  # noqa: E402

GV18K = "gvmags_v2/all_18kgenomes_fna.tar.gz"


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def gv1_contigs(tar: Path, out: Path) -> Path:
    path = out / "gv1_candidate_contigs.tsv"
    tmp = out / ".gv1_candidate_contigs.tsv.tmp"
    with look.open_member(tar, GV18K) as fh, tarfile.open(fileobj=fh, mode="r|gz") as tf, open(tmp, "w") as w:
        w.write("genome\trecord\tcontig_uid\tlength\n")
        for m in tf:
            if not m.isfile() or not m.name.endswith(".fna"):
                continue
            g = m.name.split("/")[-1][:-4]
            for k, r in enumerate(vi.iter_fasta_stream(io.BytesIO(tf.extractfile(m).read()))):
                w.write(f"{g}\t{k}\t{ids.contig_uid(r.seq).hex()}\t{len(r.seq)}\n")
    tmp.replace(path)
    return path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--stages", required=True)
    ap.add_argument("--v1", required=True, type=Path)
    ap.add_argument("--v2", required=True, type=Path)
    ap.add_argument("--gvmags", required=True, type=Path)
    ap.add_argument("--gv-genomes", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    gv_tsv = a.out / "gv1_candidate_contigs.tsv"
    if not gv_tsv.is_file():
        gv1_contigs(a.gvmags, a.out)
    with open(a.gv_genomes, newline="") as fh:
        fna_files = {r["genome"] for r in csv.DictReader(fh, delimiter="\t") if r["set"] == "fna_files"}

    con = duckdb.connect(config={"threads": 1, "memory_limit": "3GB", "temp_directory": str(a.out / "duckdb")})
    q = lambda sql: con.execute(sql).fetchall()  # noqa: E731
    rd = "delim='\t', header=true, quote='', all_varchar=true"
    con.execute(f"CREATE TABLE v2 AS SELECT genome, source, topology, unhex(contig_uid) AS cu, length::BIGINT AS n "
                f"FROM read_csv('{a.v2 / 'records.tsv'}', {rd})")
    con.execute(f"CREATE TABLE gv AS SELECT genome, unhex(contig_uid) AS cu, length::BIGINT AS n "
                f"FROM read_csv('{gv_tsv}', {rd})")
    con.execute("CREATE TABLE fna (genome VARCHAR)")
    con.executemany("INSERT INTO fna VALUES (?)", [(g,) for g in sorted(fna_files)])
    con.execute("CREATE TABLE gvf AS SELECT * FROM gv WHERE genome IN (SELECT genome FROM fna)")
    con.execute(f"CREATE TABLE v1 AS SELECT DISTINCT unhex(contig_uid) AS cu FROM read_csv('{a.v1 / 'records.tsv'}', {rd})")
    con.execute("CREATE TABLE allc AS SELECT DISTINCT cu, n FROM (SELECT cu, n FROM v2 UNION ALL SELECT cu, n FROM gv)")
    stages = a.stages.split(",")
    con.execute("CREATE TABLE own (cu BLOB, stage VARCHAR)")
    con.execute("CREATE TABLE mem (cu BLOB, stage VARCHAR)")
    for s in stages:
        ci = bbio.part_files(a.root / store.INDEX_DIR / f"stage={s}")
        gi = bbio.part_files(a.root / "l0seq/genome_input" / f"stage={s}")
        if not ci or not gi:
            raise SystemExit(f"stage {s}: contig_index {len(ci)} files, genome_input {len(gi)} files")
        con.execute(f"INSERT INTO own SELECT DISTINCT contig_uid, '{s}' FROM read_parquet({sql_list(ci)}, "
                    "hive_partitioning = false) WHERE is_default_variant AND contig_uid IN (SELECT cu FROM allc)")
        con.execute(f"INSERT INTO mem SELECT DISTINCT contig_uid, '{s}' FROM read_parquet({sql_list(gi)}, "
                    "hive_partitioning = false) WHERE contig_uid IN (SELECT cu FROM allc)")

    def one(table: str, earlier: str) -> dict:
        """`earlier`: a table of the contig uids (column cu) an earlier stage holds."""
        con.execute(f"CREATE OR REPLACE TABLE t AS SELECT x.genome, x.cu, x.n, x.cu IN (SELECT cu FROM {earlier}) AS h "
                    f"FROM {table} x")
        table = "t"
        (contigs, bases), = q(f"SELECT count(*), coalesce(sum(n), 0) FROM (SELECT DISTINCT cu, n FROM {table})")
        res: dict = {"genomes": q(f"SELECT count(DISTINCT genome) FROM {table}")[0][0], "contigs": contigs,
                     "bases": int(bases), "stages": {}}
        for s in stages:
            (oc, ob), = q(f"SELECT count(*), coalesce(sum(n), 0) FROM (SELECT DISTINCT cu, n FROM {table}) "
                          f"WHERE cu IN (SELECT cu FROM own WHERE stage = '{s}')")
            (mc,), = q(f"SELECT count(DISTINCT cu) FROM {table} WHERE cu IN (SELECT cu FROM mem WHERE stage = '{s}')")
            res["stages"][s] = {"owns_contigs": oc, "owns_bases": int(ob), "member_contigs": mc}
        (vc, vb), = q(f"SELECT count(*), coalesce(sum(n), 0) FROM (SELECT DISTINCT cu, n FROM {table}) "
                      "WHERE cu IN (SELECT cu FROM v1)")
        res["v1_inputs"] = {"contigs": vc, "bases": int(vb)}
        (sc, sb), = q(f"SELECT count(*), coalesce(sum(n), 0) FROM (SELECT DISTINCT cu, n FROM {table} WHERE h)")
        res["held_earlier"] = {"contigs": sc, "bases": int(sb)}
        res["adds"] = {"contigs": contigs - sc, "bases": int(bases) - int(sb)}
        (allh, someh), = q(f"SELECT count(*) FILTER (WHERE hk = k), count(*) FILTER (WHERE hk > 0 AND hk < k) FROM "
                           f"(SELECT genome, count(*) k, count(*) FILTER (WHERE h) hk FROM {table} GROUP BY genome)")
        res["genomes_all_contigs_held_earlier"] = allh
        res["genomes_some_contigs_held_earlier"] = someh
        return res

    con.execute("CREATE TABLE held_v2 AS SELECT cu FROM own UNION SELECT cu FROM v1")
    con.execute("CREATE TABLE held_gv AS SELECT cu FROM held_v2 UNION SELECT cu FROM v2")
    con.execute("CREATE TABLE v2_imgvr AS SELECT * FROM v2 WHERE source = 'imgvr'")
    con.execute("CREATE TABLE v2_inovirus AS SELECT * FROM v2 WHERE source = 'inovirus'")
    out: dict = {"stages": stages}
    out["V2"] = one("v2", "held_v2")
    out["V2"]["by_source"] = {s: one(f"v2_{s}", "held_v2") for s in ("imgvr", "inovirus")}
    out["V2"]["inovirus_in_imgvr_reps"] = q(
        "SELECT count(DISTINCT cu) FROM v2_inovirus WHERE cu IN (SELECT cu FROM v2_imgvr)")[0][0]
    out["GV1_all_18k"] = one("gv", "held_gv")
    out["GV1_fna_files_8508"] = one("gvf", "held_gv")
    (c, b), = q("SELECT count(*), coalesce(sum(n), 0) FROM (SELECT DISTINCT cu, n FROM gv) "
                "WHERE cu IN (SELECT cu FROM v2)")
    out["V2_and_GV1_all_18k"] = {"contigs": c, "bases": int(b), "v2_by_topology": dict(q(
        "SELECT topology, count(DISTINCT cu) FROM v2 WHERE cu IN (SELECT cu FROM gv) GROUP BY 1 ORDER BY 2 DESC")),
        "v2_gvmag_genomes": q("SELECT count(DISTINCT genome) FROM v2 WHERE topology = 'GVMAG'")[0][0],
        "v2_gvmag_contigs": q("SELECT count(DISTINCT cu) FROM v2 WHERE topology = 'GVMAG'")[0][0],
        "v2_gvmag_contigs_in_gv1": q("SELECT count(DISTINCT cu) FROM v2 WHERE topology = 'GVMAG' AND cu IN "
                                     "(SELECT cu FROM gv)")[0][0]}
    out["GV1_repeats"] = {"contigs_in_several_genomes": q(
        "SELECT count(*) FROM (SELECT cu FROM gv GROUP BY cu HAVING count(DISTINCT genome) > 1)")[0][0]}
    held = sorted(r[0] for r in q("SELECT DISTINCT lower(hex(cu)) FROM gv WHERE cu IN (SELECT cu FROM held_gv)"))
    (a.out / "gv1_held_contigs.txt").write_text("".join(f"{h}\n" for h in held))  # for the estimate step
    con.close()
    (a.out / "summary.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()

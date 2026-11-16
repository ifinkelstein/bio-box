"""Card C-39: what stage V1 (UHGV vOTU representatives) shares by contig_uid with the stages on disk (read only).

    python c39_overlap.py --root ROOT --inputs IN --uhgv-meta UHGV_METADATA_ZST --out OUT

IN holds what `c39_inputs.py uhgv ... --root ROOT` wrote (records.tsv, repeated_records.tsv, shared_records.tsv,
genomes.json). The identity is `ids.contig_uid` (the store's; computed by c39_inputs.py), no strand or rotation
canonicalisation. Writes OUT/summary.json (and prints it):

- for every stage on disk: the V1 contigs it owns (its contig_index holds the default variant, DEC-26), their records
  and bases, and the V1 contigs its genome_input lists (DEC-27), with the owned contigs' UHGV source study
  (`original_study_alias` of the representative in uhgv_metadata.tsv) and CheckV quality;
- repeats inside V1 (records whose bytes repeat an earlier record's);
- UHGV's own figures from uhgv_metadata.tsv: genomes, vOTUs, representatives, members per vOTU, the source studies of
  the representatives and of all genomes, the CheckV tiers of each, and the representatives' contigs and bases by
  CheckV tier and viral confidence and by source study, to say what the representatives stand for.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from compression import zstd
from pathlib import Path

import duckdb

from biobox import io
from biobox.l0seq import store

ZOPT = {zstd.DecompressionParameter.window_log_max: 31}  # UHGV's files use a long window (zstd --long)


def tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tiers(first: list[dict], meta: dict) -> dict:
    """(CheckV quality, viral confidence) -> [contigs, bases] over the first records."""
    out: dict = {}
    for r in first:
        m = meta[r["genome"]]
        t = out.setdefault((m["checkv_quality"], m["viral_confidence"]), [0, 0])
        t[0] += 1
        t[1] += int(r["length"])
    return out


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--inputs", required=True, type=Path)
    ap.add_argument("--uhgv-meta", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    out: dict = {}

    recs = tsv(a.inputs / "records.tsv")
    reps = tsv(a.inputs / "repeated_records.tsv")
    shared = tsv(a.inputs / "shared_records.tsv")
    rep_ord = {r["ordinal"] for r in reps}
    first = [r for r in recs if r["ordinal"] not in rep_ord]  # the first record of each contig
    length = {bytes.fromhex(r["contig_uid"]): int(r["length"]) for r in first}
    by_contig = {bytes.fromhex(r["contig_uid"]): r for r in first}
    records_of = Counter(bytes.fromhex(r["contig_uid"]) for r in recs)
    groups = [n for n in records_of.values() if n > 1]
    out["repeats"] = {"records": len(recs), "contigs": len(first), "records_beyond_the_first": len(reps),
                      "sequences_in_several_records": len(groups), "largest_group": max(groups, default=0),
                      "pairs": [(r["genome"], r["first_genome"]) for r in reps][:50]}

    meta = {}
    with zstd.open(a.uhgv_meta, "rt", newline="", options=ZOPT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["uhgv_genome"]] = r
    reps_meta = [m for m in meta.values() if m["votu_representative"] == "Yes"]
    votu_size = Counter(m["uhgv_votu"] for m in meta.values())
    sizes = sorted(votu_size.values())
    out["uhgv_tables"] = {
        "genomes": len(meta), "votus": len(votu_size), "representatives": len(reps_meta),
        "representatives_in_fasta": sum(1 for r in recs if meta.get(r["genome"], {}).get("votu_representative") == "Yes"),
        "fasta_records_not_representative": [r["genome"] for r in recs
                                             if meta.get(r["genome"], {}).get("votu_representative") != "Yes"][:20],
        "votus_of_one_genome": sum(1 for n in sizes if n == 1),
        "members_per_votu_median": sizes[len(sizes) // 2] if sizes else 0, "members_per_votu_max": max(sizes, default=0),
        "study_of_representatives": dict(Counter(m["original_study_alias"] for m in reps_meta).most_common()),
        "study_of_all_genomes": dict(Counter(m["original_study_alias"] for m in meta.values()).most_common()),
        "checkv_of_representatives": dict(Counter(m["checkv_quality"] for m in reps_meta).most_common()),
        "checkv_of_all_genomes": dict(Counter(m["checkv_quality"] for m in meta.values()).most_common()),
        "viral_confidence_of_representatives": dict(Counter(m["viral_confidence"] for m in reps_meta).most_common()),
        "representatives_by_checkv_and_confidence": {
            f"{q}|{c}": {"n": n, "bases": b} for (q, c), (n, b) in sorted(tiers(first, meta).items())},
        "representative_bases_by_study": dict(Counter(
            {s: sum(int(r["length"]) for r in first if meta[r["genome"]]["original_study_alias"] == s)
             for s in {m["original_study_alias"] for m in reps_meta}}).most_common()),
        "genome_length_of_representatives": {
            "lt_5kb": sum(1 for r in first if int(r["length"]) < 5000),
            "5_to_10kb": sum(1 for r in first if 5000 <= int(r["length"]) < 10000),
            "10_to_20kb": sum(1 for r in first if 10000 <= int(r["length"]) < 20000),
            "20_to_50kb": sum(1 for r in first if 20000 <= int(r["length"]) < 50000),
            "50_to_200kb": sum(1 for r in first if 50000 <= int(r["length"]) < 200000),
            "ge_200kb": sum(1 for r in first if int(r["length"]) >= 200000)},
    }

    con = duckdb.connect(config={"threads": 2, "memory_limit": "4GB", "temp_directory": str(a.out / "duckdb")})
    con.execute("CREATE TABLE v (contig_uid BLOB)")
    con.executemany("INSERT INTO v VALUES (?)", [(c,) for c in sorted(length)])
    owner = {bytes.fromhex(r["contig_uid"]): r["owner_stage"] for r in shared}
    per_stage = {}
    for s in store.stages_on_disk(a.root):
        row = {}
        ci = io.part_files(a.root / store.INDEX_DIR / f"stage={s}")
        gi = io.part_files(a.root / "l0seq/genome_input" / f"stage={s}")
        if ci:
            (n, b), = con.execute(
                f"SELECT count(*), coalesce(sum(contig_len), 0) FROM read_parquet({sql_list(ci)}, "
                "hive_partitioning = false) WHERE is_default_variant AND contig_uid IN (SELECT contig_uid FROM v)"
            ).fetchall()
            row["owns_contigs"], row["owns_bases"] = n, int(b)
        if gi:
            members = {r[0] for r in con.execute(
                f"SELECT DISTINCT contig_uid FROM read_parquet({sql_list(gi)}, hive_partitioning = false) "
                "WHERE contig_uid IN (SELECT contig_uid FROM v)").fetchall()}
            row["member_contigs"] = len(members)
            row["member_not_owner_contigs"] = sum(1 for c in members if owner.get(c) != s)
        owned = [c for c, st in owner.items() if st == s]
        row["owns_records"] = sum(records_of[c] for c in owned)
        row["owned_by_study"] = dict(Counter(meta.get(by_contig[c]["genome"], {}).get("original_study_alias", "?")
                                             for c in owned).most_common())
        row["owned_by_checkv"] = dict(Counter(by_contig[c]["checkv_quality"] for c in owned).most_common())
        per_stage[s] = row
    out["stages"] = per_stage
    out["shared_total"] = {"contigs": len(owner), "records": sum(records_of[c] for c in owner),
                           "bases": sum(length[c] for c in owner)}
    out["new"] = {"contigs": len(length) - len(owner), "bases": sum(n for c, n in length.items() if c not in owner)}
    (a.out / "summary.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()

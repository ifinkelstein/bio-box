"""PLSDB metadata against its FASTA and against the prokaryote stages, read-only (S2 inputs, 2026-11-11).

    python s2_plsdb_meta.py META_DIR SEQSTATS_DIR STAGE_UIDS_TSV ROOT

META_DIR: the unpacked `meta.tar.gz`. SEQSTATS_DIR: `s2_seqstats.py`'s output for the PLSDB FASTA
(`records.tsv`, `overlap.tsv`). Reads `nuccore.csv` (one row per plasmid record), `taxonomy.csv`,
`nucc_identical.csv` with the `csv` module. Reports value counts of the columns that carry topology,
completeness, source database, genome tag and status; whether the FASTA's records are the nuccore
rows (accession and length); taxonomy coverage; and the content overlaps with each stage split by
topology and source. Also an accession check: PLSDB accessions found among the stages' record
accessions (`contig_occurrence` of the stored stages in ROOT, read-only; the accession column of
STAGE_UIDS_TSV for the fetched ones), and how many of those are also content matches.
Prints JSON.
"""

import csv
import json
import sys
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

from biobox import io

csv.field_size_limit(1 << 30)


def rows(path: Path):
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh)


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    meta, ss, stage_tsv, root = Path(argv[0]), Path(argv[1]), Path(argv[2]), Path(argv[3])
    nuc = list(rows(meta / "nuccore.csv"))
    cols = list(nuc[0].keys())
    tax = {r["TAXONOMY_UID"]: r for r in rows(meta / "taxonomy.csv")}
    ident = list(rows(meta / "nucc_identical.csv"))
    by_acc = {r["NUCCORE_ACC"]: r for r in nuc}
    rec = {}
    with open(ss / "records.tsv") as fh:
        fh.readline()
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            rec[f[1]] = int(f[3])
    in_both = set(rec) & set(by_acc)
    len_diff = sum(1 for a in in_both if str(rec[a]) != by_acc[a]["NUCCORE_Length"])

    def counts(col, rs=nuc):
        return dict(Counter(r.get(col, "") for r in rs).most_common())

    ov = {}
    with open(ss / "overlap.tsv") as fh:
        fh.readline()
        for ln in fh:
            _, acc, _, stage, how = ln.rstrip("\n").split("\t")
            if how == "same":
                ov.setdefault(stage, []).append(by_acc.get(acc, {}))
    stage_acc: dict[str, set[str]] = {}
    with open(stage_tsv) as fh:
        fh.readline()
        for ln in fh:
            s, how, _, _, a = ln.rstrip("\n").split("\t")
            if how == "fasta":
                stage_acc.setdefault(f"{s}/fasta", set()).add(a)
    for part in sorted((root / "l0seq" / "contig_occurrence").glob("stage=*")):
        s = part.name.split("=", 1)[1]
        for f in io.part_files(part, "source=*/part-*.parquet"):
            stage_acc.setdefault(f"{s}/store", set()).update(pq.read_table(f, columns=["accession"])
                                                             .column("accession").to_pylist())
    content, uid_of = {}, {}
    with open(ss / "overlap.tsv") as fh:
        fh.readline()
        for ln in fh:
            _, acc, uid, stage, how = ln.rstrip("\n").split("\t")
            if how == "same":
                content.setdefault(stage, set()).add(acc)
                uid_of[acc] = uid
    # one set per stage: the stored copy where the stage is stored (P3 is both; s2_stage_uids checks they agree)
    stage_sets = {}
    for s in sorted(content):
        name, how = s.split("/")
        if how == "store" or f"{name}/store" not in content:
            stage_sets[name] = content[s]
    prok = set().union(*(v for k, v in stage_sets.items() if k.startswith("P")))
    n_stages = Counter(a for k, v in stage_sets.items() if k.startswith("P") for a in v)
    union = {"stages_counted": sorted(stage_sets),
             "plsdb_records_in_any_P_stage": len(prok),
             "distinct_uids_in_any_P_stage": len({uid_of[a] for a in prok}),
             "plsdb_records_in_two_or_more_P_stages": sum(1 for v in n_stages.values() if v > 1),
             "plsdb_records_in_any_stage_incl_S1": len(set().union(*stage_sets.values())),
             "P_records_by_topology": dict(Counter(by_acc.get(a, {}).get("NUCCORE_Topology", "(none)")
                                                   for a in prok).most_common())}
    acc_check = {}
    for s, accs in sorted(stage_acc.items()):
        hit = set(by_acc) & accs
        acc_check[s] = {"stage_accessions": len(accs), "plsdb_accessions_found": len(hit),
                        "of_which_content_match": len(hit & content.get(s, set())),
                        "content_match_without_accession_match": len(content.get(s, set()) - hit)}
    rep = {
        "nuccore_rows": len(nuc), "nuccore_columns": cols, "nuccore_distinct_acc": len(by_acc),
        "nuccore_length_sum": sum(int(r["NUCCORE_Length"]) for r in nuc),
        "topology": counts("NUCCORE_Topology"), "completeness": counts("NUCCORE_Completeness"),
        "genome": counts("NUCCORE_Genome"), "source": counts("NUCCORE_Source"), "status": counts("STATUS"),
        "with_assembly": sum(1 for r in nuc if r["ASSEMBLY_UID"] not in ("", "-1")),
        "with_biosample": sum(1 for r in nuc if r["BIOSAMPLE_UID"] not in ("", "-1")),
        "taxonomy_rows": len(tax), "nuccore_with_taxonomy_row": sum(1 for r in nuc if r["TAXONOMY_UID"] in tax),
        "superkingdom": dict(Counter(tax.get(r["TAXONOMY_UID"], {}).get("TAXONOMY_superkingdom", "(none)")
                                     for r in nuc).most_common()),
        "nucc_identical_rows": len(ident), "nucc_identical_targets": len({r["NUCCORE_Identical"] for r in ident}),
        "fasta_records": len(rec), "fasta_in_nuccore": len(in_both), "fasta_only": len(set(rec) - set(by_acc)),
        "nuccore_only": len(set(by_acc) - set(rec)), "length_differs": len_diff,
        "overlap_same_by_stage": {s: {"records": len(v), "topology": counts("NUCCORE_Topology", v),
                                      "source": counts("NUCCORE_Source", v),
                                      "completeness": counts("NUCCORE_Completeness", v)}
                                  for s, v in sorted(ov.items())},
        "accession_check": acc_check,
        "union": union,
    }
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

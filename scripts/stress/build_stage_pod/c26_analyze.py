"""Card C-26 part 1: what the PLSDB stage PL1 adds, read from its inputs and the data root (read only).

    python c26_analyze.py --root ROOT --inputs IN --meta META_DIR --out OUT

IN holds the inputs `make_inputs.py plasmid ... --root ROOT` wrote (genomes.json, repeated_records.tsv,
shared_records.tsv); META_DIR is PLSDB's unpacked metadata (nuccore.csv, proteins.csv, proteins.fasta). Nothing in
ROOT is written. Writes OUT/summary.json (and prints it):

- records per contig (the first record of a sequence is its genome, repeated_records.tsv the others): sequences
  that several records carry, their records, the largest group, and how often the repeats disagree on topology or
  quality;
- for every stage on disk, the PL1 contigs it owns (its contig_index holds the default variant: DEC-26) and the PL1
  contigs its genome_input lists (its genomes that carry the contig: DEC-27), in contigs, records and bases;
- the density of the plasmids already called: for the PL1 contigs that a stage P1 to P6 owns, that stage's genes
  (feature_protein), their distinct proteins and residues (protein_seq `length`), per Mbp;
- PGAP's annotation of the same records (proteins.csv rows per accession; proteins.fasta sequences, hashed with
  ids.protein_uid after removing a trailing '*'): CDS per Mbp over every record and over the first records of the
  contigs PL1 stores, and the distinct proteins and residues of those;
- the bytes each stage's files hold in ROOT (paths with `stage=<S>/`) beside the bases its contig_index stores.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path

import duckdb

from biobox import fasta, ids, io
from biobox.l0seq import store


def tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--inputs", required=True, type=Path)
    ap.add_argument("--meta", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    out: dict = {}

    genomes = json.loads((a.inputs / "genomes.json").read_text())
    # the gene-calling mode the code gives these genomes (genecall._finder: kind plasmid is meta mode at any length)
    from biobox.cli import _genomes
    from biobox.l1 import genecall

    gs = _genomes(str(a.inputs / "genomes.json"))
    out["kinds"] = dict(Counter(g.kind for g in gs))
    out["gene_call_mode"] = {f"{g.kind} {n} bp": genecall._finder(g, {g.contig_uids[0]: "ACGT" * (n // 4)})[1]
                             for g in gs[:1] for n in (4_000, 400_000)}
    out["genetic_codes"] = dict(Counter(g["genetic_code"] for g in genomes))
    contigs = {bytes.fromhex(g["contig_uids"][0]) for g in genomes}
    reps = tsv(a.inputs / "repeated_records.tsv")
    shared = tsv(a.inputs / "shared_records.tsv")
    records_of = Counter({cu: 1 for cu in contigs})
    for r in reps:
        records_of[bytes.fromhex(r["contig_uid"])] += 1
    groups = [n for n in records_of.values() if n > 1]
    out["repeats"] = {"sequences_in_several_records": len(groups), "records_in_those": sum(groups),
                      "records_beyond_the_first": len(reps), "largest_group": max(groups, default=0),
                      "group_sizes": dict(sorted(Counter(groups).items())),
                      "repeats_other_topology": sum(r["same_topology"] == "false" for r in reps),
                      "repeats_other_quality": sum(r["same_quality"] == "false" for r in reps)}

    # contig lengths of the PL1 contigs: from PLSDB's nuccore.csv (NUCCORE_Length equals the FASTA length for
    # every record, measured 2026-11-11) by the first record of each contig
    length_of_acc = {}
    values = {c: Counter() for c in ("NUCCORE_Source", "NUCCORE_Genome", "VIRAL_category", "VIRAL_putative", "STATUS",
                                     "NUCCORE_Topology", "NUCCORE_Completeness")}
    with open(a.meta / "nuccore.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            length_of_acc[r["NUCCORE_ACC"]] = int(r["NUCCORE_Length"])
            for c, n in values.items():
                n[r.get(c) or ""] += 1
    out["metadata_values_by_record"] = {c: dict(sorted(n.items())) for c, n in values.items()}
    owned_contigs = {r["contig_uid"] for r in tsv(a.inputs / "shared_records.tsv")}
    new_g = [g for g in genomes if g["contig_uids"][0] not in owned_contigs]
    out["genomes"] = {
        "all": {"n": len(genomes), "circular": sum(1 for g in genomes if g["circular"]),
                "quality": dict(Counter(g["assembly_quality"] for g in genomes))},
        "on_contigs_the_stage_stores": {"n": len(new_g), "circular": sum(1 for g in new_g if g["circular"]),
                                        "quality": dict(Counter(g["assembly_quality"] for g in new_g))},
        "on_contigs_an_earlier_stage_owns": {"n": len(genomes) - len(new_g),
                                             "circular_by_plsdb": sum(1 for g in genomes if g["circular"]) -
                                             sum(1 for g in new_g if g["circular"])}}
    rep_accs = {r["accession"] for r in reps}
    owner = {bytes.fromhex(r["contig_uid"]): r["owner_stage"] for r in shared}

    con = duckdb.connect(config={"threads": 2, "memory_limit": "4GB", "temp_directory": str(a.out / "duckdb")})
    con.execute("CREATE TABLE pl (contig_uid BLOB)")
    con.executemany("INSERT INTO pl VALUES (?)", [(c,) for c in sorted(contigs)])

    stages = store.stages_on_disk(a.root)
    per_stage = {}
    for s in stages:
        ci = io.part_files(a.root / store.INDEX_DIR / f"stage={s}")
        gi = io.part_files(a.root / "l0seq/genome_input" / f"stage={s}")
        row = {}
        if ci:
            (n, b), = con.execute(
                f"SELECT count(*), coalesce(sum(contig_len), 0) FROM read_parquet({sql_list(ci)}, "
                "hive_partitioning = false) WHERE is_default_variant AND contig_uid IN (SELECT contig_uid FROM pl)"
            ).fetchall()
            row["owns_contigs"], row["owns_bases"] = n, int(b)
        if gi:
            (n,), = con.execute(
                f"SELECT count(DISTINCT contig_uid) FROM read_parquet({sql_list(gi)}, hive_partitioning = false) "
                "WHERE contig_uid IN (SELECT contig_uid FROM pl)").fetchall()
            members = {r[0] for r in con.execute(
                f"SELECT DISTINCT contig_uid FROM read_parquet({sql_list(gi)}, hive_partitioning = false) "
                "WHERE contig_uid IN (SELECT contig_uid FROM pl)").fetchall()}
            row["member_contigs"] = n
            row["member_records"] = sum(records_of[c] for c in members)
            row["member_not_owner_contigs"] = sum(1 for c in members if owner.get(c) != s)
        row["owns_records"] = sum(1 for r in shared if r["owner_stage"] == s)
        per_stage[s] = row
    out["stages"] = per_stage

    # density of the plasmids the P stages already called (their genes in genome context, single mode)
    gcv = sorted(os.listdir(a.root / "l1" / "feature_protein"))
    if len(gcv) != 1:
        raise SystemExit(f"gene call versions {gcv}: expected one")
    gcv = gcv[0].removeprefix("gene_call_version=")
    dens = {"genes": 0, "bases": 0, "contigs": 0}
    con.execute("CREATE TABLE prot (protein_uid BLOB)")
    for s in [s for s in stages if s.startswith("P")]:
        fp = io.part_files(a.root / "l1/feature_protein" / f"gene_call_version={gcv}" / f"stage={s}")
        ci = io.part_files(a.root / store.INDEX_DIR / f"stage={s}")
        con.execute(f"CREATE OR REPLACE TEMP VIEW own AS SELECT contig_uid, contig_len FROM read_parquet({sql_list(ci)},"
                    " hive_partitioning = false) WHERE is_default_variant AND contig_uid IN (SELECT contig_uid FROM pl)")
        (nc, nb), = con.execute("SELECT count(*), coalesce(sum(contig_len), 0) FROM own").fetchall()
        (ng,), = con.execute(f"SELECT count(*) FROM read_parquet({sql_list(fp)}, hive_partitioning = false) "
                             "WHERE contig_uid IN (SELECT contig_uid FROM own)").fetchall()
        con.execute(f"INSERT INTO prot SELECT DISTINCT protein_uid FROM read_parquet({sql_list(fp)}, "
                    "hive_partitioning = false) WHERE contig_uid IN (SELECT contig_uid FROM own)")
        dens["genes"] += ng
        dens["bases"] += int(nb)
        dens["contigs"] += nc
    ps = [f for s in stages if s.startswith("P") for f in io.part_files(
        a.root / "l0seq/protein_seq" / f"gene_call_version={gcv}" / f"stage={s}", "protein_prefix=*/part-*.parquet")]
    (nd, res), = con.execute(
        f"SELECT count(*), sum(length) FROM (SELECT DISTINCT protein_uid, length FROM read_parquet({sql_list(ps)}, "
        "hive_partitioning = false) WHERE protein_uid IN (SELECT protein_uid FROM prot))").fetchall()
    dens.update({"distinct_proteins": nd, "residues": int(res or 0), "gene_call_version": gcv})
    if dens["bases"]:
        dens["genes_per_mbp"] = round(dens["genes"] / dens["bases"] * 1e6, 3)
        dens["distinct_per_gene"] = round(nd / dens["genes"], 4)
        dens["residues_per_distinct"] = round(res / nd, 2)
    out["p_stage_plasmids"] = dens

    # PGAP's annotation: CDS per record (proteins.csv), proteins of the records (proteins.fasta)
    cds = Counter()
    with open(a.meta / "proteins.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            cds[r["NUCCORE_ACC"]] += 1
    # a record is the first record of its contig exactly when it is not in repeated_records.tsv; the first records of
    # the contigs PL1 stores are those no earlier stage owns (not in shared_records.tsv)
    accs_all = set(length_of_acc)
    first_accs = accs_all - rep_accs
    shared_accs = {r["accession"] for r in shared}
    new_first = first_accs - shared_accs
    out["pgap"] = {
        "records_with_cds": sum(1 for x in accs_all if cds[x]), "records": len(accs_all),
        "cds_all_records": sum(cds.values()), "bases_all_records": sum(length_of_acc.values()),
        "first_records_new": len(new_first), "first_records_new_with_cds": sum(1 for x in new_first if cds[x]),
        "cds_first_records_new": sum(cds[x] for x in new_first),
        "bases_first_records_new": sum(length_of_acc[x] for x in new_first),
        "bases_first_records_new_with_cds": sum(length_of_acc[x] for x in new_first if cds[x]),
    }
    prot: dict[bytes, int] = {}
    n_seq = 0
    for r in fasta.iter_records(a.meta / "proteins.fasta"):
        parts = r.header.split()
        if len(parts) < 2 or parts[1] not in new_first:
            continue
        aa = r.seq.decode().rstrip("*")
        n_seq += 1
        prot.setdefault(ids.protein_uid(aa), len(aa))
    out["pgap"].update({"proteins_first_records_new": n_seq, "distinct_proteins_first_records_new": len(prot),
                        "residues_distinct_first_records_new": sum(prot.values())})

    # bytes of each stage's files in the root, beside the bases the stage stores
    size = Counter()
    for dp, _, fs in os.walk(a.root):
        for f in fs:
            p = os.path.join(dp, f)
            for part in Path(p).relative_to(a.root).parts:
                if part.startswith("stage="):
                    size[part.removeprefix("stage=")] += os.path.getsize(p)
                    break
    stored = {}
    for s in stages:
        ci = io.part_files(a.root / store.INDEX_DIR / f"stage={s}")
        (b,), = con.execute(f"SELECT coalesce(sum(contig_len), 0) FROM read_parquet({sql_list(ci)}, "
                            "hive_partitioning = false) WHERE is_default_variant").fetchall()
        stored[s] = {"bytes": size[s], "stored_bases": int(b)}
    out["root_bytes_by_stage"] = stored
    (a.out / "summary.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()

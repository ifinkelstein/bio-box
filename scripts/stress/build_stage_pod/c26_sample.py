"""Card C-26 part 1: a sample of the PLSDB records for the rehearsal, copied byte for byte from the PLSDB FASTA.

    python3 -I c26_sample.py IN NUCCORE_CSV FASTA OUT_FASTA OUT_TSV [P1,P2]

IN is the full stage's inputs directory (`make_inputs.py plasmid ... --root ROOT`: repeated_records.tsv,
shared_records.tsv). The sample, chosen in FASTA order so that it is the same on every run:

- owned: 40 records whose contig stage P1 owns and 60 whose contig P2 owns (DEC-26: they stay with that stage);
- repeat groups owned by P1 or P2: up to 3 groups, every record of each (several PLSDB records and a P genome);
- repeat groups owned by no stage: up to 15 groups of at most 4 records each, the contig under 60 kbp;
- new: 150 records owned by no stage and repeating no other record, under 60 kbp, every 50th such record, so that the
  sample spans the file; with 10 of them linear and 10 with an empty NUCCORE_Completeness where the file has them.

No record whose contig another stage (S1, P3 to P6, H1) owns is taken: the rehearsal root links only the stages
named (default P1,P2). Writes OUT_FASTA (the records' lines as in FASTA) and OUT_TSV (ordinal, accession, why).
"""
import csv
import sys

inp, meta_csv, fa, out_fa, out_tsv = sys.argv[1:6]
linked = (sys.argv[6] if len(sys.argv) > 6 else "P1,P2").split(",")


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


meta = {}
order = []
with open(meta_csv, newline="") as fh:
    for r in csv.DictReader(fh):
        meta[r["NUCCORE_ACC"]] = r
reps = tsv(f"{inp}/repeated_records.tsv")
shared = tsv(f"{inp}/shared_records.tsv")
owner_of_acc = {r["accession"]: r["owner_stage"] for r in shared}
contig_of_rep = {r["accession"]: r["contig_uid"] for r in reps}
group = {}  # contig -> accessions of its records, the first record first
for r in reps:
    group.setdefault(r["contig_uid"], [r["first_accession"]]).append(r["accession"])
owner_of_contig = {r["contig_uid"]: r["owner_stage"] for r in shared}

# the FASTA order of the accessions (headers only)
with open(fa, "rb") as fh:
    for line in fh:
        if line.startswith(b">"):
            order.append(line[1:].split()[0].decode())
pos = {a: i for i, a in enumerate(order)}
in_group = {a for g in group.values() for a in g}
pick: dict[str, str] = {}

for st, n in zip(linked, (40, 60)):
    for a in order:
        if n == 0:
            break
        if owner_of_acc.get(a) == st and a not in in_group:
            pick[a] = f"owned:{st}"
            n -= 1
k = 0
for cu, accs in sorted(group.items(), key=lambda t: pos[t[1][0]]):
    if k < 3 and owner_of_contig.get(cu) in linked and len(accs) <= 4:
        for a in accs:
            pick[a] = f"repeat-owned:{owner_of_contig[cu]}"
        k += 1
k = 0
for cu, accs in sorted(group.items(), key=lambda t: pos[t[1][0]]):
    if k < 15 and cu not in owner_of_contig and len(accs) <= 4 and int(meta[accs[0]]["NUCCORE_Length"]) < 60_000:
        for a in accs:
            pick[a] = "repeat"
        k += 1
new = [a for a in order if a not in owner_of_acc and a not in in_group and int(meta[a]["NUCCORE_Length"]) < 60_000]
for a in new[::50][:150]:
    pick[a] = "new"
for cond, why in ((lambda r: r["NUCCORE_Topology"] == "linear", "new-linear"),
                  (lambda r: r["NUCCORE_Completeness"] == "", "new-incomplete")):
    n = 10
    for a in new:
        if n == 0:
            break
        if a not in pick and cond(meta[a]):
            pick[a] = why
            n -= 1

k = 0
with open(fa, "rb") as fh, open(out_fa, "wb") as w:
    keep = False
    for line in fh:
        if line.startswith(b">"):
            keep = line[1:].split()[0].decode() in pick
            k += keep
        if keep:
            w.write(line)
with open(out_tsv, "w") as w:
    w.write("ordinal\taccession\twhy\n")
    for a in sorted(pick, key=pos.get):
        w.write(f"{pos[a]}\t{a}\t{pick[a]}\n")
print({"records": k, "picked": len(pick), "by_why": {y: sum(1 for v in pick.values() if v == y)
                                                      for y in sorted(set(pick.values()))}})

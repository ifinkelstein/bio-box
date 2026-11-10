"""Step 6, compare: bio-box vs samtools slices, vs samtools-built indexes, vs Prodigal 2.6.3, vs RefSeq CDS.

    python step6_compare.py STEP6_DIR ROOT
"""

import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

if len(sys.argv) < 3 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)

d, root = Path(sys.argv[1]), Path(sys.argv[2])
out = {}

# 1. slices
for strand in ("plus", "minus"):
    exp = [ln.split("\t") for ln in (d / f"expected_{strand}.tsv").read_text().splitlines()]
    got = [ln.split("\t") for ln in (d / "pod" / f"{strand}.tsv").read_text().splitlines()]
    same = sum(1 for e, g in zip(exp, got) if e == g)
    out[f"samtools_slices_{strand}"] = {"expected": len(exp), "returned": len(got), "identical": same}

# 2. indexes written by the repo vs indexes samtools builds for the same shard
shard = root / "l0seq/seq_shards/stage=S0/shard-0000.fa.gz"
for suf in (".fai", ".gzi"):
    a, b = Path(str(shard) + suf).read_bytes(), (d / "pod" / f"shard-0000.fa.gz{suf}").read_bytes()
    out[f"index{suf}_identical_to_samtools"] = a == b
    if a != b:
        out[f"index{suf}_sizes"] = [len(a), len(b)]


# 3. gene calls
def load_gff(path, opener=open):
    genes = defaultdict(set)
    with opener(path, "rt") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            if "pseudo=true" in f[8] or "partial=true" in f[8]:
                continue
            genes[f[0]].add((int(f[3]) - 1, int(f[4]), 1 if f[6] == "+" else -1))
    return genes


bb = defaultdict(set)
lens = {}
for ln in (d / "genes.tsv").read_text().splitlines():
    acc, s, e, st, n = ln.split("\t")
    bb[acc].add((int(s), int(e), int(st)))
    lens[acc] = int(n)


def stop_key(g, n):
    s, e, st = g
    return (e % n if e > n else e, st) if st == 1 else (s, st)


def compare(a, b, n):
    """a = bio-box, b = reference; exact = same start, stop, strand; stop = same 3' end and strand."""
    exact = a & b
    sa, sb = {stop_key(g, n) for g in a}, {stop_key(g, n) for g in b}
    return {"biobox": len(a), "reference": len(b), "exact": len(exact), "same_stop": len(sa & sb),
            "only_biobox_by_stop": len(sa - sb), "only_reference_by_stop": len(sb - sa)}


prod = {}
tot = defaultdict(int)
for gff in sorted((d / "pod").glob("*.gff")):
    ref = load_gff(gff)
    for acc, genes in ref.items():
        key = acc.split()[0]
        c = compare(bb[key], genes, lens[key])
        prod[f"{gff.stem}:{key}"] = c
        for k, v in c.items():
            tot[k] += v
        if c["exact"] != c["biobox"] or c["exact"] != c["reference"]:
            c["only_biobox"] = sorted(bb[key] - genes)[:6]
            c["only_reference"] = sorted(genes - bb[key])[:6]
out["prodigal_2.6.3_single_mode"] = {"total": dict(tot), "per_contig": prod}

ref = load_gff(d / "GCF_000005845.2.gff.gz", gzip.open)
c = compare(bb["NC_000913.3"], ref["NC_000913.3"], lens["NC_000913.3"])
c["sensitivity_by_stop"] = round(c["same_stop"] / c["reference"], 4)
c["exact_start_given_same_stop"] = round(c["exact"] / c["same_stop"], 4)
out["refseq_cds_ecoli_k12"] = c
print(json.dumps(out, indent=1))

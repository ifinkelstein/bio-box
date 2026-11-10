"""Score distributions of both ssRNAP profiles and the low tail, to choose an inclusion cutoff."""
import math

from common import OUT, by_query, coverage, hist, read_domtbl

pfam = by_query(read_domtbl(f"{OUT}/pfam_RNA_pol.domtbl"))["RNA_pol"]
s8 = by_query(read_domtbl(f"{OUT}/seed8.domtbl"))["ssRNAP_seed8"]
panel = by_query(read_domtbl(f"{OUT}/pfam_panel.domtbl"))


def best_panel(u, skip=()):
    best = None
    for q, hits in panel.items():
        if q in skip or u not in hits:
            continue
        if best is None or hits[u]["seq_score"] > best[1]:
            best = (q, hits[u]["seq_score"])
    return best


edges = [10, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1500, 3000]
for name, hits in (("Pfam RNA_pol (PF00940.25, M=449, GA 22)", pfam), ("ssRNAP_seed8 (M=815)", s8)):
    print(f"== {name}: {len(hits)} proteins reported at >= 10 bits")
    for lo, hi, n in hist([h["seq_score"] for h in hits.values()], edges):
        print(f"   [{lo:>4}, {hi:>4}) bits: {n}")
    print("   low tail (10-80 bits): score, E at Z=1e8, length, coverage, best other panel model")
    for u, h in sorted(hits.items(), key=lambda x: x[1]["seq_score"]):
        if h["seq_score"] >= 80:
            break
        bp = best_panel(u, skip=("RNA_pol",))
        print(f"     {h['seq_score']:6.1f} {h['seq_evalue']:9.2g} {h['tlen']:5d} aa cov {coverage(h):.2f}  {u}  {bp}")

# overlap of the two at a few cutoffs
for c in (25, 30, 40, 50):
    a = {u for u, h in pfam.items() if h["seq_score"] >= c}
    b = {u for u, h in s8.items() if h["seq_score"] >= c}
    print(f"cutoff {c}: Pfam {len(a)}, seed8 {len(b)}, both {len(a & b)}, Pfam only {len(a - b)}, seed8 only {len(b - a)}")
# E-value of a bitscore at Z=1e8 for the record (HMMER: E = Z * exp(-lambda*(S - tau)) for the forward tail)
print("panel models with hits:", {q: len(h) for q, h in panel.items()})

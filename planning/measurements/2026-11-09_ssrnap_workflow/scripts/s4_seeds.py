"""Step 1b: build an ssRNAP profile from seed sequences, as a user with a few seeds would.

1. Find phage T7's RNAP in the store: the T7 contig (AY264774) by accession, its slots, the slot
   protein that the Pfam RNA_pol model hits.
2. Candidates: confident, full-length Pfam RNA_pol hits (seq score >= 300 bits, >= 90% model
   coverage, 750-1200 aa), one per NR50 cluster (the cluster representative's sequence when it is
   itself a candidate, else the best-scoring member).
3. No MSA program in the image (no mafft/muscle/clustalo/famsa): build a profile from T7 alone,
   hmmalign the candidates to it, and pick 8 seeds by greedy max-min identity from T7 over the
   aligned match columns.
4. hmmalign the 8 seeds to the T7 profile and hmmbuild the final model from that alignment.
"""
import io
import json
import subprocess
import time

import pyhmmer

from common import OUT, by_query, catalog, coverage, read_domtbl, uid_table

T0 = time.time()
AMINO = pyhmmer.easel.Alphabet.amino()
HMMD = "/tmp/if2453/ssrnap/hmm"
con = catalog()
pf = by_query(read_domtbl(f"{OUT}/pfam_RNA_pol.domtbl"))["RNA_pol"]

# 1. T7's RNAP
t7 = con.execute("""
  SELECT lower(hex(o.protein_uid)), o.slot_ordinal, o.start, o."end", o.strand, p.length
  FROM contig_occurrence c JOIN protein_occurrence o USING (contig_uid) JOIN protein_seq p USING (protein_uid)
  WHERE c.accession = 'AY264774' ORDER BY o.slot_ordinal""").fetchall()
t7_hits = [r for r in t7 if r[0] in pf]
print(f"T7 contig AY264774: {len(t7)} slots; Pfam RNA_pol hits: "
      + "; ".join(f"slot {r[1]} {r[2]}-{r[3]} ({r[4]:+d}) {r[5]} aa, {pf[r[0]]['seq_score']:.1f} bits" for r in t7_hits))
t7_uid = max(t7_hits, key=lambda r: pf[r[0]]["seq_score"])[0]

# 2. candidates, one per NR50 cluster
print(f"T7 RNAP Pfam coverage {coverage(pf[t7_uid]):.2f}, domains {pf[t7_uid]['domains']}")
conf = [u for u, h in pf.items() if h["seq_score"] >= 300 and coverage(h) >= 0.85 and 750 <= h["tlen"] <= 1200]
conf = sorted(set(conf) | {t7_uid})
uid_table(con, "cand", conf)
rows = con.execute("""
  SELECT lower(hex(c.protein_uid)), lower(hex(m.cluster_uid)), lower(hex(m.rep_uid)), p.seq
  FROM cand c JOIN cluster_member m ON m.member_uid = c.protein_uid AND m.level = 'nr50'
  JOIN protein_seq p ON p.protein_uid = c.protein_uid""").fetchall()
seq = {r[0]: r[3].rstrip("*") for r in rows}
clusters = {}
for u, cl, rep, _ in rows:
    clusters.setdefault(cl, []).append(u)
pick = {}
for cl, members in clusters.items():
    reps = [m for m in members if any(r[0] == m and r[2] == m for r in rows)]
    pick[cl] = reps[0] if reps else max(members, key=lambda m: pf[m]["seq_score"])
t7_cl = next(cl for cl, ms in clusters.items() if t7_uid in ms)
pick[t7_cl] = t7_uid
print(f"confident full-length Pfam hits: {len(conf)} proteins in {len(clusters)} NR50 clusters")

# 3. T7-only profile; align candidates; greedy max-min identity
t7_seq = pyhmmer.easel.TextSequence(name=b"T7_gp1", sequence=seq[t7_uid]).digitize(AMINO)
t7_hmm, _, _ = pyhmmer.plan7.Builder(AMINO, seed=42).build(t7_seq, pyhmmer.plan7.Background(AMINO))
t7_hmm.name = b"T7_gp1_single"
with open(f"{HMMD}/t7_single.hmm", "wb") as fh:
    t7_hmm.write(fh)
cands = sorted(pick.values())
block = pyhmmer.easel.DigitalSequenceBlock(AMINO, [
    pyhmmer.easel.TextSequence(name=u.encode(), sequence=seq[u]).digitize(AMINO) for u in cands])
aligner = pyhmmer.plan7.TraceAligner()
msa = aligner.align_traces(t7_hmm, block, aligner.compute_traces(t7_hmm, block), trim=True)
aln = {n.decode() if isinstance(n, bytes) else n: s for n, s in zip(msa.names, msa.alignment)}
match_cols = [i for i, c in enumerate(aln[t7_uid]) if c.isalpha() and c.isupper() or c == "-"]


def ident(a, b):
    sa, sb = aln[a], aln[b]
    same = tot = 0
    for i in match_cols:
        x, y = sa[i], sb[i]
        if x.isupper() and y.isupper():
            tot += 1
            same += x == y
    return same / tot if tot else 0.0


seeds = [t7_uid]
while len(seeds) < 8:
    best = max((u for u in cands if u not in seeds), key=lambda u: min(1 - ident(u, s) for s in seeds))
    seeds.append(best)
meta = {r[0]: r for r in con.execute("""
  SELECT lower(hex(s.protein_uid)), any_value(c.accession), any_value(c.original_header), count(*)
  FROM cand s JOIN protein_occurrence o ON o.protein_uid = s.protein_uid
  JOIN contig_occurrence c ON c.contig_uid = o.contig_uid GROUP BY 1""").fetchall()}
print("seeds (greedy max-min identity from T7):")
for u in seeds:
    print(f"  {u} {len(seq[u])} aa, Pfam {pf[u]['seq_score']:.0f} bits, id to T7 {ident(u, t7_uid):.2f}, "
          f"max id to another seed {max(ident(u, s) for s in seeds if s != u):.2f}, {meta[u][3]} occ, "
          f"{meta[u][1]} {meta[u][2]}")
with open(f"{HMMD}/seeds.faa", "w") as fh:
    for u in seeds:
        fh.write(f">{u}\n{seq[u]}\n")

# 4. align seeds to the T7 profile (HMMER's hmmalign binary), hmmbuild
subprocess.run(["hmmalign", "--trim", "--amino", "-o", f"{HMMD}/seeds.sto", f"{HMMD}/t7_single.hmm",
                f"{HMMD}/seeds.faa"], check=True)
subprocess.run(["hmmbuild", "--amino", "-n", "ssRNAP_seed8", "-o", f"{HMMD}/hmmbuild.log", f"{HMMD}/ssRNAP_seed8.hmm",
                f"{HMMD}/seeds.sto"], check=True)
print(open(f"{HMMD}/hmmbuild.log").read().split("# idx")[1])
json.dump({"t7_uid": t7_uid, "seeds": seeds}, open(f"{HMMD}/seeds.json", "w"))
print(f"# {time.time() - T0:.1f} s")

"""Biological sanity checks on the hit set: split genes, doubles, fusions, contig ends, giant vRNAPs."""
import csv
from collections import Counter, defaultdict

from common import OUT, by_query, catalog, read_domtbl, uid_table

rows = list(csv.DictReader(open(f"{OUT}/export/hits.tsv"), delimiter="\t"))
pfam = by_query(read_domtbl(f"{OUT}/pfam_RNA_pol.domtbl"))["RNA_pol"]
s8 = by_query(read_domtbl(f"{OUT}/seed8.domtbl"))["ssRNAP_seed8"]
con = catalog()
by_c = defaultdict(list)
for r in rows:
    by_c[r["contig_uid"]].append(r)
kinds = Counter()
pairs = Counter()
examples = defaultdict(list)
for cu, rr in by_c.items():
    if len(rr) < 2:
        continue
    rr.sort(key=lambda r: int(r["slot_ordinal"]))
    full = [r for r in rr if r["klass"].startswith("full")]
    slots = [int(r["slot_ordinal"]) for r in rr]
    adjacent = max(slots) - min(slots) <= len(rr)
    k = ("full+partial" if full and len(full) < len(rr) else "2+ full" if len(full) >= 2 else "all partial") + \
        (", adjacent" if adjacent else ", apart")
    kinds[k] += 1
    if k.startswith("all partial") and adjacent:
        pairs[tuple(sorted(int(r["length"]) for r in rr))] += 1
    examples[k].append(f"{rr[0]['accession']} {rr[0]['header'][:60]} | " +
                       "; ".join(f"slot {r['slot_ordinal']} {r['length']} aa {r['strand']} {r['pfam_architecture']} "
                                 f"(Pfam {r['pfam_bits'] or '-'}, seed8 {r['s8_bits'] or '-'})" for r in rr))
print("contigs with 2+ hits:", dict(kinds))
for k, ex in sorted(examples.items()):
    print(f"  {k}: e.g.")
    for e in ex[:4]:
        print("    ", e)
print("adjacent all-partial length tuples (top 12):", pairs.most_common(12))
for label, pred in (("fusion >1300 aa", lambda r: r["klass"] == "partial(>1300aa)"),
                    ("RNA_pol+GFP", lambda r: "GFP" in r["pfam_architecture"]),
                    ("full but low score", lambda r: r["klass"] == "partial(700-1300,lowaa)")):
    sel = [r for r in rows if pred(r)]
    print(f"{label}: {len(sel)} genes")
    for r in sel[:6]:
        print(f"   {r['accession']} {r['header'][:60]} {r['length']} aa {r['pfam_architecture']} Pfam {r['pfam_bits']} "
              f"seed8 {r['s8_bits']} dist_to_end {r['dist_to_contig_end_bp']}")
end = [r for r in rows if r["dist_to_contig_end_bp"] == "0"]
print("hit genes touching a contig end:", len(end), Counter(r["klass"] for r in end))

# giant virion RNAPs (N4 vRNAP-like) by Pfam in the store, and what our two profiles gave them
g = con.execute("""
  SELECT lower(hex(d.protein_uid)), p.length, string_agg(DISTINCT d.model_name, '+')
  FROM l1_domains d JOIN protein_seq p USING (protein_uid)
  WHERE d.model_name IN ('vRNAP_dom', 'vRNAP_dom_2', 'vRNAP_plug', 'N4_RNAP_helical') GROUP BY 1, 2""").fetchall()
L = [x[1] for x in g]
print(f"proteins with N4 vRNAP Pfam domains: {len(g)}; length quartiles "
      f"{sorted(L)[len(L) // 4] if L else None}/{sorted(L)[len(L) // 2] if L else None}/{sorted(L)[3 * len(L) // 4] if L else None}, "
      f">=2500 aa: {sum(x >= 2500 for x in L)}")
best = sorted(((pfam.get(u, {}).get("seq_score"), s8.get(u, {}).get("seq_score"), n, u) for u, n, _ in g),
              key=lambda x: -(x[1] or 0))[:5]
print("  best seed8 / Pfam RNA_pol scores on them:", [(b, a, n) for a, b, n, _ in best])
n4 = con.execute("""
  SELECT o.slot_ordinal, o.start, o."end", p.length, lower(hex(o.protein_uid))
  FROM contig_occurrence c JOIN protein_occurrence o USING (contig_uid) JOIN protein_seq p USING (protein_uid)
  WHERE c.accession = 'NC_008720' AND p.length > 2500""").fetchall()
for s, a, b, ln, u in n4:
    arch = con.execute("SELECT string_agg(model_name, '+' ORDER BY env_start) FROM l1_domains WHERE protein_uid = ?",
                       [bytes.fromhex(u)]).fetchall()[0][0]
    print(f"  N4 slot {s} {a}-{b} {ln} aa: Pfam {arch}; RNA_pol {pfam.get(u, {}).get('seq_score')}, "
          f"seed8 {s8.get(u, {}).get('seq_score')}")
# phiKMV-like (mid-genome RNAP) vs T7-like (early) position, by genus-like grouping via the neighbour kinase
pos = defaultdict(list)
for r in rows:
    if r["klass"].startswith("full"):
        L_ = int(r["contig_len"])
        x = (int(r["start"]) if r["strand"] == "1" else L_ - int(r["end"])) / L_
        pos["early (<0.2)" if x < 0.2 else "mid (0.2-0.6)" if x < 0.6 else "late (>=0.6)"].append(r["header"])
for k, v in pos.items():
    words = Counter(w for h in v for w in h.split()[1:3])
    print(f"position {k}: {len(v)} genes; commonest host words {words.most_common(6)}")

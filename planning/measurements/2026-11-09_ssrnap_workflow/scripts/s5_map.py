"""Steps 3-5: map hits to genes and genomes, neighborhoods (biobox.search.nbr), exports.

usage: s5_map.py <cutoff_bits>
Hit set: proteins with Pfam RNA_pol OR ssRNAP_seed8 sequence score >= cutoff (both at Z = 1e8).
"""
import json
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow as pa

from biobox.search import export, nbr
from common import KNOWN, OUT, Z, by_query, catalog, coverage, hist, read_domtbl, uid_table

T = {}
t0 = time.time()
CUT = float(sys.argv[1])
ROOT = Path("/stor/scratch/Wilke/if2453/all-the-data/provisional/rS1.pilot1")  # durable copy, read-only
D = Path(OUT) / "export"
D.mkdir(exist_ok=True)

pfam = by_query(read_domtbl(f"{OUT}/pfam_RNA_pol.domtbl"))["RNA_pol"]
s8 = by_query(read_domtbl(f"{OUT}/seed8.domtbl"))["ssRNAP_seed8"]
panel = by_query(read_domtbl(f"{OUT}/pfam_panel.domtbl"))
hits = sorted({u for u, h in pfam.items() if h["seq_score"] >= CUT} | {u for u, h in s8.items() if h["seq_score"] >= CUT})
print(f"hit proteins (Pfam RNA_pol or seed8 >= {CUT} bits): {len(hits)}")


def sc(tab, u, k="seq_score"):
    return tab[u][k] if u in tab else None


def length(u):
    return (pfam.get(u) or s8.get(u))["tlen"]


def klass(u):
    L, p = length(u), max(sc(pfam, u) or 0, 0)
    if L >= 2500:
        return "giant(>=2500aa)"
    if 700 <= L <= 1300 and (p >= 100 or (sc(s8, u) or 0) >= 200):
        return "full(700-1300aa)"
    return f"partial({'<500' if L < 500 else '500-699' if L < 700 else '>1300' if L > 1300 else '700-1300,low'}aa)"


con = catalog()
uid_table(con, "hit", hits, {
    "pfam_bits": pa.array([sc(pfam, u) for u in hits], pa.float64()),
    "pfam_E": pa.array([sc(pfam, u, "seq_evalue") for u in hits], pa.float64()),
    "s8_bits": pa.array([sc(s8, u) for u in hits], pa.float64()),
    "s8_E": pa.array([sc(s8, u, "seq_evalue") for u in hits], pa.float64()),
    "klass": pa.array([klass(u) for u in hits]),
})
t1 = time.time()
occ = con.execute("""
  WITH acc AS (SELECT contig_uid, min(accession) AS accession, min(original_header) AS header,
                      count(*) AS n_records FROM contig_occurrence GROUP BY contig_uid)
  SELECT lower(hex(h.protein_uid)) AS protein_uid, h.klass, h.pfam_bits, h.pfam_E, h.s8_bits, h.s8_E,
         p.length, lower(hex(o.contig_uid)) AS contig_uid, lower(hex(o.genome_uid)) AS genome_uid, o.slot_ordinal,
         o.start, o."end", o.strand, o.dist_to_contig_end_bp, t.contig_len, t.is_circular, t.n_slots,
         a.accession, a.header, a.n_records,
         lower(hex(m90.cluster_uid)) AS nr90, lower(hex(m50.cluster_uid)) AS nr50, lower(hex(m50.rep_uid)) AS nr50_rep
  FROM hit h JOIN protein_occurrence o ON o.protein_uid = h.protein_uid
  JOIN protein_seq p ON p.protein_uid = h.protein_uid
  JOIN contig_tapes t ON t.contig_uid = o.contig_uid
  JOIN acc a ON a.contig_uid = o.contig_uid
  LEFT JOIN cluster_member m90 ON m90.member_uid = h.protein_uid AND m90.level = 'nr90'
  LEFT JOIN cluster_member m50 ON m50.member_uid = h.protein_uid AND m50.level = 'nr50'
  ORDER BY a.accession, o.start""").fetchall()
cols = ["protein_uid", "klass", "pfam_bits", "pfam_E", "s8_bits", "s8_E", "length", "contig_uid", "genome_uid",
        "slot_ordinal", "start", "end", "strand", "dist_to_contig_end_bp", "contig_len", "is_circular", "n_slots",
        "accession", "header", "n_records", "nr90", "nr50", "nr50_rep"]
occ = [dict(zip(cols, r)) for r in occ]
T["map_query_s"] = time.time() - t1
print(f"gene occurrences: {len(occ)}; genomes: {len({o['genome_uid'] for o in occ})}; contigs: "
      f"{len({o['contig_uid'] for o in occ})}; NR90 clusters: {len({o['nr90'] for o in occ})}; "
      f"NR50 clusters: {len({o['nr50'] for o in occ})}  (map query {T['map_query_s']:.1f} s)")

# --- classes, lengths, hits per genome ---
by_cls = defaultdict(set)
for u in hits:
    by_cls[klass(u)].add(u)
print("classes (distinct proteins / gene occurrences / genomes):")
for k in sorted(by_cls):
    oo = [o for o in occ if o["klass"] == k]
    print(f"  {k}: {len(by_cls[k])} / {len(oo)} / {len({o['genome_uid'] for o in oo})}")
L = [length(u) for u in hits]
print("protein length (distinct hit proteins):", hist(L, [0, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
                                                        1300, 2000, 2500, 3500, 5000]))
print(f"  median {statistics.median(L)}, min {min(L)}, max {max(L)}")
per_g = Counter(o["genome_uid"] for o in occ)
print("hits per genome:", sorted(Counter(per_g.values()).items()))
per_g_full = Counter(o["genome_uid"] for o in occ if o["klass"].startswith("full"))
print("full-length hits per genome:", sorted(Counter(per_g_full.values()).items()))
both = Counter((o["genome_uid"]) for o in occ)

# panel profile of each class: which other models hit them
print("panel models on hit proteins, by class (n proteins with model >= 25 bits):")
for k in sorted(by_cls):
    c = Counter(q for q, hh in panel.items() for u in by_cls[k] if u in hh and hh[u]["seq_score"] >= 25)
    print(f"  {k}: {dict(c.most_common())}")
# Pfam-vs-seed8 differences
po = {u for u in hits if (sc(pfam, u) or 0) >= CUT and (sc(s8, u) or 0) < CUT}
so = {u for u in hits if (sc(s8, u) or 0) >= CUT and (sc(pfam, u) or 0) < CUT}
print(f"Pfam-only: {len(po)} {Counter(klass(u) for u in po)}; seed8-only: {len(so)} {Counter(klass(u) for u in so)}")
for name, s in (("Pfam-only", po), ("seed8-only", so)):
    ex = defaultdict(list)
    for o in occ:
        if o["protein_uid"] in s:
            ex[o["protein_uid"]].append(o)
    print(f"  {name} examples:")
    for u in sorted(s, key=lambda u: -(sc(s8, u) or sc(pfam, u) or 0))[:12]:
        o = ex[u][0]
        bp = sorted(((q, hh[u]["seq_score"]) for q, hh in panel.items() if u in hh), key=lambda x: -x[1])[:3]
        cv = f"{coverage(s8[u]):.2f}" if u in s8 else "-"
        cp = f"{coverage(pfam[u]):.2f}" if u in pfam else "-"
        print(f"    {u} {length(u)} aa Pfam {sc(pfam, u)} (cov {cp}) seed8 {sc(s8, u)} (cov {cv}) | "
              f"{o['accession']} {o['header'][:70]} | {bp}")

# --- known positives ---
print("known ssRNAP phages (accession found by header text):")
for lab, acc in KNOWN.items():
    rows = [o for o in occ if o["accession"] == acc]
    allacc = con.execute("SELECT count(*) FROM contig_occurrence WHERE accession = ?", [acc]).fetchall()[0][0]
    if not allacc:
        print(f"  {lab} {acc}: not in INPHARED")
        continue
    if not rows:  # the contig may be listed under another accession of the same sequence
        cu = con.execute("SELECT lower(hex(contig_uid)) FROM contig_occurrence WHERE accession = ?", [acc]).fetchall()[0][0]
        rows = [o for o in occ if o["contig_uid"] == cu]
    print(f"  {lab} {acc}: " + ("; ".join(
        f"slot {o['slot_ordinal']} {o['start']}-{o['end']} ({o['strand']:+d}) {o['length']} aa Pfam {o['pfam_bits']} "
        f"seed8 {o['s8_bits']} [{o['klass']}]" for o in rows) or "NO HIT"))

# --- position in genome ---
pos = [(o["start"] / o["contig_len"]) for o in occ if o["klass"].startswith("full")]
posn = [((o["start"] if o["strand"] == 1 else o["contig_len"] - o["end"]) / o["contig_len"])
        for o in occ if o["klass"].startswith("full")]
fwd = sum(1 for o in occ if o["klass"].startswith("full") and o["strand"] == 1)
print(f"full-length RNAP position, start/contig_len (as deposited): deciles "
      f"{[round(x, 2) for x in statistics.quantiles(pos, n=10)]}; on + strand: {fwd}/{len(pos)}")
print(f"  strand-normalised (distance from the end the gene points away from): deciles "
      f"{[round(x, 2) for x in statistics.quantiles(posn, n=10)]}; frac < 0.2: {sum(x < 0.2 for x in posn) / len(posn):.2f}")
print("  slot ordinal (strand-normalised) median:",
      statistics.median([(o["slot_ordinal"] if o["strand"] == 1 else o["n_slots"] - 1 - o["slot_ordinal"])
                         for o in occ if o["klass"].startswith("full")]))
partial_genes = [o for o in occ if o["dist_to_contig_end_bp"] == 0]
print(f"hit genes at a contig end (dist_to_contig_end_bp = 0, possibly truncated): {len(partial_genes)}")

# --- neighborhoods via the project's code ---
t1 = time.time()
policy = nbr.Policy.parse(window="5g", strand="any", edge="truncate")
rows, counts = nbr.fetch(ROOT, {bytes.fromhex(u) for u in hits}, policy)
T["nbr_fetch_s"] = time.time() - t1
print(f"nbr.fetch window 5g: {counts}, {T['nbr_fetch_s']:.1f} s")
with open(D / "neighborhoods.ndjson", "w") as fh:  # the project's exporter, used outside a run
    for line in export.ndjson(export.nbr_records(rows)):
        fh.write(line + "\n")

# --- neighbour annotation: NR50 cluster and Pfam 38.2 domains from the store (l1_domains) ---
t1 = time.time()
nb_uids = sorted({r["protein_uid"] for r in rows})
uid_table(con, "nb", [u.hex() for u in nb_uids])
cl = {r[0]: r[1:] for r in con.execute("""
  SELECT lower(hex(n.protein_uid)), lower(hex(m.cluster_uid)), lower(hex(m.rep_uid)), p.length
  FROM nb n JOIN cluster_member m ON m.member_uid = n.protein_uid AND m.level = 'nr50'
  JOIN protein_seq p ON p.protein_uid = n.protein_uid""").fetchall()}
rep_len = dict(con.execute("""
  SELECT lower(hex(m.rep_uid)), any_value(p.length) FROM cluster_member m JOIN protein_seq p ON p.protein_uid = m.rep_uid
  WHERE m.level = 'nr50' AND m.rep_uid IN (SELECT m2.rep_uid FROM nb n JOIN cluster_member m2
       ON m2.member_uid = n.protein_uid AND m2.level = 'nr50') GROUP BY 1""").fetchall())
doms = defaultdict(list)
for u, name, acc, s in con.execute("""
  SELECT lower(hex(d.protein_uid)), d.model_name, d.model_acc, d.env_start FROM l1_domains d
  WHERE d.protein_uid IN (SELECT protein_uid FROM nb) ORDER BY 1, 4""").fetchall():
    doms[u].append(name)
lib = con.execute("SELECT DISTINCT library, cutoff FROM l1_domains").fetchall()
T["nbr_annot_s"] = time.time() - t1


def arch(u):
    a = []
    for n in doms.get(u, []):
        if not a or a[-1] != n:
            a.append(n)
    return "+".join(a)


genome_of = {o["contig_uid"]: o["genome_uid"] for o in occ}
hitset = set(hits)
n_nb = {r["protein_uid"].hex() for r in rows if r["offset"] != 0}
print(f"neighbour annotation ({lib}): {sum(1 for u in n_nb if u in doms)}/{len(n_nb)} distinct neighbour proteins "
      f"have a Pfam domain; neighbour slots with one: "
      f"{sum(1 for r in rows if r['offset'] and r['protein_uid'].hex() in doms)}/{sum(1 for r in rows if r['offset'])} "
      f"({T['nbr_annot_s']:.1f} s)")
print("hit proteins' own Pfam architectures (l1_domains, distinct proteins):")
for k in sorted(by_cls):
    c = Counter(arch(u) or "(none)" for u in by_cls[k])
    print(f"  {k}: {c.most_common(8)}")

# cluster-mates' Pfam families, for unannotated NR50 clusters
def summarise(rows, label, n_top=15):
    anchors = {(r["contig_uid"], r["anchor_ordinal"]) for r in rows}
    fam_loci, fam_g, fam_off = defaultdict(set), defaultdict(set), defaultdict(Counter)
    cl_loci, cl_g, cl_off = defaultdict(set), defaultdict(set), defaultdict(Counter)
    for r in rows:
        if r["offset"] == 0:
            continue
        key = (r["contig_uid"], r["anchor_ordinal"])
        u = r["protein_uid"].hex()
        g = genome_of[r["contig_uid"].hex()]
        if u in doms:
            for f in set(doms[u]):
                fam_loci[f].add(key)
                fam_g[f].add(g)
                fam_off[f][r["offset_rel"]] += 1
        else:
            c = cl[u][0]
            cl_loci[c].add(key)
            cl_g[c].add(g)
            cl_off[c][r["offset_rel"]] += 1
    n = len(anchors)
    out = {"n_loci": n, "families": [], "clusters": []}
    print(f"== {label}: {n} anchor loci; offset_rel = slots along the RNAP's strand, > 0 downstream")
    print(f"   top {n_top} Pfam families among neighbours (loci with the family within +-5 slots):")
    for f in sorted(fam_loci, key=lambda f: (-len(fam_loci[f]), f))[:n_top]:
        o = dict(sorted(fam_off[f].items()))
        print(f"     {f}: {len(fam_loci[f])} loci ({len(fam_loci[f]) / n:.0%}), {len(fam_g[f])} genomes, offsets {o}")
        out["families"].append({"family": f, "n_loci": len(fam_loci[f]), "n_genomes": len(fam_g[f]), "offsets": o})
    top = sorted(cl_loci, key=lambda c: (-len(cl_loci[c]), c))[:n_top]
    if top:
        uid_table(con, "topcl", top)
        mates = defaultdict(Counter)
        for c, name in con.execute("""
          SELECT lower(hex(m.cluster_uid)), d.model_name FROM cluster_member m
          JOIN l1_domains d ON d.protein_uid = m.member_uid
          WHERE m.level = 'nr50' AND m.cluster_uid IN (SELECT protein_uid FROM topcl)""").fetchall():
            mates[c][name] += 1
    print(f"   top {n_top} NR50 clusters among neighbours with no Pfam domain:")
    for c in top:
        rep = next(cl[r["protein_uid"].hex()][1] for r in rows if r["offset"] and cl[r["protein_uid"].hex()][0] == c)
        o = dict(sorted(cl_off[c].items()))
        hit_flag = " [cluster holds RNAP hits]" if any(v[0] == c and u in hitset for u, v in cl.items()) else ""
        print(f"     {c[:12]} rep {rep[:12]} {rep_len.get(rep)} aa: {len(cl_loci[c])} loci ({len(cl_loci[c]) / n:.0%}), "
              f"{len(cl_g[c])} genomes, offsets {o}; Pfam on cluster-mates {dict(mates[c].most_common(3))}{hit_flag}")
        out["clusters"].append({"nr50": c, "rep": rep, "rep_len": rep_len.get(rep), "n_loci": len(cl_loci[c]),
                                "n_genomes": len(cl_g[c]), "offsets": o, "mates_pfam": dict(mates[c].most_common(3))})
    return out


full_anchor = {bytes.fromhex(u) for u in hits if klass(u).startswith("full")}
part_anchor = {bytes.fromhex(u) for u in hits if klass(u).startswith("partial")}
summ = {"all": summarise(rows, "all hits"),
        "full": summarise([r for r in rows if r["anchor_protein_uid"] in full_anchor], "full-length hits"),
        "partial": summarise([r for r in rows if r["anchor_protein_uid"] in part_anchor], "partial hits")}
json.dump(summ, open(f"{OUT}/neighbour_summary.json", "w"), indent=1)

# --- example loci ---
hscore = {o["protein_uid"]: (o["pfam_bits"], o["s8_bits"]) for o in occ}
LOCI = {"T7": "AY264774", "phiKMV": "NC_005045", "N4": "NC_008720"}
lcols = ["offset", "offset_rel", "start", "end", "strand", "length_aa", "pfam_architecture", "nr50_cluster",
         "rnap_hit_pfam_bits", "rnap_hit_seed8_bits", "protein_uid"]
for lab, acc in LOCI.items():
    cu = con.execute("SELECT contig_uid FROM contig_occurrence WHERE accession = ?", [acc]).fetchall()[0][0]
    hdr = con.execute("SELECT original_header FROM contig_occurrence WHERE accession = ?", [acc]).fetchall()[0][0]
    rr = sorted((r for r in rows if r["contig_uid"] == cu), key=lambda r: (r["anchor_ordinal"], r["offset"]))
    with open(D / f"locus_{lab}_{acc}.tsv", "w") as fh:
        fh.write(f"# {hdr}\n# anchor(s): RNAP hit gene(s) at offset 0; window +-5 slots (bio-box nbr, policy 5g)\n")
        fh.write("anchor_ordinal\t" + "\t".join(lcols) + "\n")
        print(f"== locus {lab} {acc} ({hdr})")
        for r in rr:
            u = r["protein_uid"].hex()
            hs = hscore.get(u, (None, None))
            vals = [r["offset"], r["offset_rel"], r["start"], r["end"], r["strand"], cl[u][2] if u in cl else "",
                    arch(u) or "-", cl[u][0][:12], hs[0] if hs[0] is not None else "", hs[1] if hs[1] is not None else "", u]
            fh.write(f"{r['anchor_ordinal']}\t" + "\t".join(str(v) for v in vals) + "\n")
            print(f"   anchor {r['anchor_ordinal']:3d} " + "  ".join(str(v) for v in vals[:9]))

# --- exports ---
sel = ["accession", "header", "contig_len", "start", "end", "strand", "length", "pfam_bits", "pfam_E", "s8_bits", "s8_E",
       "klass", "pfam_architecture", "protein_uid", "nr90", "nr50", "nr50_rep", "genome_uid", "contig_uid",
       "slot_ordinal", "dist_to_contig_end_bp", "n_records"]
uid_table(con, "hit2", hits)
hdoms = defaultdict(list)
for u, name in con.execute("""SELECT lower(hex(d.protein_uid)), d.model_name FROM l1_domains d
                              WHERE d.protein_uid IN (SELECT protein_uid FROM hit2) ORDER BY 1, d.env_start""").fetchall():
    if not hdoms[u] or hdoms[u][-1] != name:
        hdoms[u].append(name)
with open(D / "hits.tsv", "w") as fh:
    fh.write("\t".join(sel) + "\n")
    for o in occ:
        o["pfam_architecture"] = "+".join(hdoms.get(o["protein_uid"], [])) or "-"
        fh.write("\t".join("" if o[c] is None else (f"{o[c]:.3g}" if c.endswith("_E") else str(o[c])) for c in sel) + "\n")
first = {}
for o in occ:
    first.setdefault(o["protein_uid"], o)
nocc = Counter(o["protein_uid"] for o in occ)
seq_of = dict(con.execute("SELECT lower(hex(h.protein_uid)), p.seq FROM hit h JOIN protein_seq p USING (protein_uid)"
                          ).fetchall())
with open(D / "hits.faa", "w") as fh:
    for u in hits:
        o = first[u]
        fh.write(f">{u} len={length(u)} pfam_bits={sc(pfam, u)} seed8_bits={sc(s8, u)} class={klass(u)} "
                 f"n_genes={nocc[u]} example={o['accession']}:{o['start']}-{o['end']}({o['strand']:+d}) {o['header']}\n"
                 f"{seq_of[u].rstrip('*')}\n")
# the project's hit exporter, fed with rows shaped like search_hit from the hmmsearch domtbl (outside any run)
fake = [{"target_kind": "protein", "target_uid": bytes.fromhex(u), "seq_bitscore": sc(pfam, u),
         "seq_evalue_zref": sc(pfam, u, "seq_evalue"), "bitscore": max(d["score"] for d in pfam[u]["domains"]),
         "target_size": 1134550, "model_len": 449} for u in hits if u in pfam]
occ_by = defaultdict(list)
for r in nbr.occurrences(ROOT, {bytes.fromhex(u) for u in hits}):
    occ_by[r["protein_uid"]].append(r)
with open(D / "hits_project_exporter.ndjson", "w") as fh:
    for line in export.ndjson(export.hit_records(fake, occ_by, top_n=None)):
        fh.write(line + "\n")
T["total_s"] = time.time() - t0
print("timings:", {k: round(v, 1) for k, v in T.items()})

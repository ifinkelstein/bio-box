"""Card 0-08 design note: the model rules applied to Rfam 15.1, with the 0-09c measurement's per-model CPU.

Two steps. Neither runs a search; both are arithmetic over files that already exist.

    python -I 0-08_model_rules.py join <0-08_rfam15.1_models.tsv> <m09c out dir> <out: 0-08_model_costs.tsv>
    python -I 0-08_model_rules.py rules <0-08_rfam15.1_models.tsv> <0-08_model_costs.tsv>

`join` (needs the measurement's outputs, on the pod under /tmp/if2453/m09c/out): per model, the truth hits
of the whole-genome producer-settings scan of the 8 genomes (`merged/wg_default.hits.jsonl`), the genomes
they lie in, and the CPU seconds of the model's own `# CPU time` line (user + system) in the runs
`wgA_default` + `wgB_default` (whole genome), `igr_m200_notrunc` and `igr_m0_notrunc` (`*.permodel.json`,
written by `scripts/stress/annotation_cost/infernal_igr.py permodel`).

`rules` (needs only the two TSVs beside this script): classifies every model by the rules below and prints
the model count and CPU-hours per 10 Gbp of each rule and of each step of the proposed set. CPU-hours per
10 Gbp is derived as  sum of the set's per-model CPU-s / 18,106,231 genome bases * 1e10 / 3600  (the bases
of the 8 genomes the runs cover; the intergenic runs cover the same genomes, so the margin is inside the
number). The sum of the per-model lines is 0.4% below the processes' measured CPU (12,403.2 against
12,450.1 CPU-s for the whole-genome runs), so full-library figures here read 1,903 where the measurement
reads 1,910.

The class rules read the model's NAME and DESC only. They are a heuristic: Rfam's family type (`TP`) and
clan membership are in `family.txt` and `Rfam.clanin`, which this project has not fetched.
"""

from __future__ import annotations

import collections
import csv
import json
import re
import sys
from pathlib import Path

GENOME_BASES = 18_106_231  # the 8 genomes of the 0-09c measurement


def h10(cpu_s: float) -> float:
    return cpu_s / GENOME_BASES * 1e10 / 3600


def _desc(pattern: str):
    rx = re.compile(pattern, re.I)
    return lambda n, s: bool(rx.search(s))


# (b) included whatever the 8 genomes show. Edit here; order matters only for the report.
INCLUDE = [
    ("b1", "rRNA (5S, 5.8S, SSU, LSU models)", lambda n, s: bool(re.match(r"(5S|5_8S)_rRNA$|(SSU|LSU)_rRNA_", n))),
    ("b2", "tRNA, tRNA-Sec", lambda n, s: n in ("tRNA", "tRNA-Sec")),
    ("b3", "RNase P, bacterial and archaeal", lambda n, s: n in ("RNaseP_bact_a", "RNaseP_bact_b", "RNaseP_arch",
                                                                 "RNaseP-T")),
    ("b4", "tmRNA", lambda n, s: n in ("tmRNA", "alpha_tmRNA", "beta_tmRNA", "cyano_tmRNA", "tmRNA_permuted")),
    ("b5", "SRP, bacterial and archaeal", lambda n, s: n in ("Bacteria_small_SRP", "Bacteria_large_SRP",
                                                             "Archaea_SRP", "Bacteroidales_small_SRP")),
    ("b6", "6S", lambda n, s: n in ("6S", "6S-Flavo", "Bacteroidales-1")),
    ("b7", "CRISPR direct repeats, tracrRNA", lambda n, s: n.startswith("CRISPR-DR") or n == "tracrRNA"),
    ("b8", "group I and II introns", lambda n, s: n in ("Intron_gpI", "Intron_gpII") or n.startswith("group-II-D1D4")),
    ("b9", "riboswitches, T-box", lambda n, s: bool(re.search(r"riboswitch|T-box", s, re.I)) or "T-box" in n),
    ("b10", "sRNAs named for a model prokaryote", lambda n, s: bool(re.search(
        r"\b(E\.? ?coli|Escherichia|Enterobacteri|Salmonella|Shigella|Bacillus|B\. subtilis|Mycobacteri|"
        r"M\. tuberculosis|Staphylococc|S\. aureus|Pseudomonas|Vibrio|Listeria|Streptococc|Synechocystis|"
        r"Caulobacter|Sinorhizobium|Clostridi|Helicobacter|Neisseria)", s + " " + n, re.I))),
    # Rfam describes the archaeal box C/D sRNAs as "Small nucleolar RNA sR<n>"; 10 of them hit the archaeon
    ("b11", "archaeal box C/D sRNAs (described as snoRNA)", lambda n, s: (
        bool(re.match(r"sR\d+$", n)) and bool(re.search(r"nucleolar|snoRNA", s, re.I))) or n == "snoPyro_CD"),
]
# (c) excluded: eukaryote-only classes, read from the description
EXCLUDE = [
    ("c1", "snoRNA, scaRNA", _desc(r"small nucleolar|snoRNA|scaRNA|Cajal")),
    ("c2", "microRNA", _desc(r"microRNA")),
    ("c3", "snRNA (spliceosomal, U7)", _desc(r"spliceosomal|small nuclear")),
    ("c4", "lncRNA (named human and mouse lncRNAs)", lambda n, s: bool(re.search(
        r"lncRNA|long non-?coding|long ncRNA|lincRNA|\bXIST|HOTAIR|MALAT|NEAT1|MEG3|\bH19\b|KCNQ1OT1|TUG1|PVT1|GAS5|"
        r"HOTTIP|MESTIT|CDKN2B|Evf|DLEU|FTX|JPX|Taurine upregulated|Maternally expressed|X-chromosome inactivation|"
        r"Deleted in lymphocytic|Metastasis associated|Pvt1|Hovlinc", s + " " + n, re.I))),
    ("c5", "vault", _desc(r"vault")),
    ("c6", "telomerase", _desc(r"telomerase")),
    ("c7", "other eukaryotic housekeeping (nuclear RNase P, MRP, eukaryotic SRP, 7SK, Y RNA)", lambda n, s: bool(
        re.search(r"Nuclear RNase P|RNase MRP|Plasmodium RNase_P|Metazoan signal|Fungal signal|Plant signal|"
                  r"Protozoan signal|Dictyostelium signal|7SK", s, re.I)) or n == "Y_RNA"),
    ("c8", "elements of eukaryotic viruses", lambda n, s: bool(re.search(r"virus|viral|HIV|HCV|retro", s, re.I))
     and not re.search(r"phage", s, re.I)),
    ("c9", "organellar (mitochondrion, plastid)", lambda n, s: bool(
        re.search(r"mitochondri|plastid|chloroplast", s, re.I) or n.startswith("mt-"))),
]
# a hit model dropped anyway: the eukaryotic rRNA models, whose every hit in the 8 genomes lies on a locus
# that a bacterial or archaeal rRNA model also hits (see the note, section 2)
EUK_RRNA = ("SSU_rRNA_eukarya", "SSU_rRNA_microsporidia", "LSU_rRNA_eukarya", "5_8S_rRNA")
# RNA that lies inside coding sequence: the 11 families every intergenic margin lost entirely in 0-09c
# (margin 200, --notrunc), and models whose name or description says antisense, frameshift or SECIS
IN_CDS_LOST = ("DnaX", "AS1726", "ASdes", "ASpks", "Ms_AS-5", "5_ureB_sRNA", "MicL", "AsdA", "AsrC", "Erse", "Ysr276")
IN_CDS_RX = re.compile(r"antisense|asRNA|frameshift|recoding|SECIS|selenocysteine insertion", re.I)


def classify(r: dict) -> tuple[str, str]:
    n, s = r["name"], r["desc"]
    if int(r["hits"]) > 0:
        return "a", "a"
    for k, _, f in INCLUDE:
        if f(n, s):
            return "b", k
    for k, _, f in EXCLUDE:
        if f(n, s):
            return "c", k
    return "d", "d"


def join(models_tsv: Path, out_dir: Path, out_tsv: Path) -> None:
    hdr = list(csv.DictReader(open(models_tsv), delimiter="\t"))

    def cpu(names):
        acc = collections.defaultdict(float)
        for n in names:
            for m in json.loads((out_dir / f"{n}.permodel.json").read_text()):
                acc[m["acc"]] += m["cpu_s"] or 0.0
        return acc

    wg, m200, m0 = cpu(["wgA_default", "wgB_default"]), cpu(["igr_m200_notrunc"]), cpu(["igr_m0_notrunc"])
    hits = [json.loads(x) for x in (out_dir / "merged/wg_default.hits.jsonl").read_text().splitlines() if x.strip()]
    nh = collections.Counter(h["acc"] for h in hits)
    gen = collections.defaultdict(set)
    for h in hits:
        gen[h["acc"]].add(h["genome"])
    assert {r["acc"] for r in hdr} == set(wg) == set(m200) == set(m0)
    with open(out_tsv, "w") as fo:
        fo.write("acc\tname\thits\tgenomes\tcpu_s_wg\tcpu_s_igr200nt\tcpu_s_igr0nt\n")
        for r in sorted(hdr, key=lambda r: r["acc"]):
            a = r["acc"]
            fo.write(f"{a}\t{r['name']}\t{nh.get(a, 0)}\t{','.join(sorted(gen.get(a, ()))) or '-'}\t"
                     f"{wg[a]:.2f}\t{m200[a]:.2f}\t{m0[a]:.2f}\n")


def rules(models_tsv: Path, costs_tsv: Path) -> None:
    desc = {r["acc"]: r for r in csv.DictReader(open(models_tsv), delimiter="\t")}
    rows = []
    for r in csv.DictReader(open(costs_tsv), delimiter="\t"):
        r.update(desc=desc[r["acc"]]["desc"], clen=desc[r["acc"]]["clen"], bps=desc[r["acc"]]["bps"])
        r["grp"], r["rule"] = classify(r)
        rows.append(r)
    tot = {k: sum(float(r[k]) for r in rows) for k in ("cpu_s_wg", "cpu_s_igr200nt", "cpu_s_igr0nt")}

    def cost(sel):
        return (len(sel),) + tuple(sum(float(r[k]) for r in sel) for k in ("cpu_s_wg", "cpu_s_igr200nt", "cpu_s_igr0nt"))

    def line(label, sel):
        n, w, m, z = cost(sel)
        print(f"| {label} | {n:,} | {w:,.1f} | {w / tot['cpu_s_wg']:.1%} | {h10(w):,.0f} | {h10(m):,.0f} | {h10(z):,.0f} |")

    print("| rule | models | CPU-s, whole genome | share | CPU-h/10 Gbp, whole genome | intergenic 200, --notrunc | intergenic 0, --notrunc |")
    print("|---|---|---|---|---|---|---|")
    line("full library", rows)
    line("(a) a hit in the 8 genomes", [r for r in rows if r["grp"] == "a"])
    for k, label, _ in INCLUDE:
        line(f"({k}) {label}", [r for r in rows if r["rule"] == k])
    for k, label, _ in EXCLUDE:
        line(f"({k}) {label}", [r for r in rows if r["rule"] == k])
    line("(d) long tail", [r for r in rows if r["grp"] == "d"])
    # would an exclusion rule have removed a model that hits prokaryotes or phages?
    for k, label, f in EXCLUDE:
        bad = [r["name"] for r in rows if r["grp"] == "a" and f(r["name"], r["desc"])]
        if bad:
            print(f"rule {k} also matches hit models: {bad}")
    euk = [r for r in rows if r["name"] in EUK_RRNA]
    print("eukaryotic rRNA models:", [(r["name"], r["hits"]) for r in euk])
    print("\n| step | models | CPU-s, whole genome | share | CPU-h/10 Gbp, whole genome | intergenic 200, --notrunc | intergenic 0, --notrunc |")
    print("|---|---|---|---|---|---|---|")
    a = [r for r in rows if r["grp"] == "a"]
    ab = [r for r in rows if r["grp"] in "ab"]
    abx = [r for r in ab if r["name"] not in EUK_RRNA]
    abx_no10 = [r for r in abx if r["rule"] != "b10"]
    line("0. today: full library", rows)
    line("1. (a)", a)
    line("2. (a) + (b)", ab)
    line("2'. (a) + (b) without the 4 eukaryotic rRNA models", abx)
    line("2''. as 2' without (b10)", abx_no10)
    line("3. full library without (c1)-(c6)", [r for r in rows if not (r["grp"] == "c" and r["rule"] <= "c6"
                                                                    and len(r["rule"]) == 2)])
    line("3'. full library without (c1)-(c9)", [r for r in rows if r["grp"] != "c"])
    line("(d) alone", [r for r in rows if r["grp"] == "d"])
    rr = [r for r in rows if re.match(r"(5S|5_8S)_rRNA$|(SSU|LSU)_rRNA_", r["name"])]
    n, w, m, _ = cost(rr)
    print(f"\nrRNA models: {n}, {w:.1f} CPU-s wg ({w / tot['cpu_s_wg']:.1%}), {m:.1f} igr200nt ({m / tot['cpu_s_igr200nt']:.1%})")
    for sel_name, sel in (("2'", abx), ("2", ab)):
        _, sw, sm, _ = cost(sel)
        _, rw, rm, _ = cost([r for r in rr if r in sel])
        print(f"rRNA share of step {sel_name}: wg {rw / sw:.1%}, igr200nt {rm / sm:.1%}")
    for r in rr:
        print(f"  {r['acc']} {r['name']} clen {r['clen']} hits {r['hits']} wg {r['cpu_s_wg']} igr200nt {r['cpu_s_igr200nt']}")
    t = [r for r in rows if r["name"] in ("tRNA", "tRNA-Sec")]
    print("tRNA models:", [(r["name"], r["hits"], r["cpu_s_wg"]) for r in t])
    cds = [r for r in abx if r["name"] in IN_CDS_LOST or IN_CDS_RX.search(r["desc"] + " " + r["name"])]
    _, cw, cm, _ = cost(cds)
    _, sw, sm, _ = cost(abx)
    print(f"\nin-CDS class within step 2': {len(cds)} models, wg {h10(cw):.1f} h, their intergenic share {h10(cm):.1f} h;"
          f" step 2' intergenic {h10(sm):.1f} h, with the class scanned whole genome {h10(sm - cm + cw):.1f} h")
    print("  ", ", ".join(f"{r['acc']} {r['name']}" for r in sorted(cds, key=lambda r: r["acc"])))
    lost = [r for r in rows if r["name"] in IN_CDS_LOST]
    _, lw, lm, _ = cost(lost)
    print(f"the 11 lost families alone, whole genome: {h10(lw):.1f} h; step 2' + them whole genome {h10(sm - lm + lw):.1f} h")
    _, dw, dm, _ = cost([r for r in rows if r["grp"] == "d"])
    print(f"(d) long tail on a 425 Mbp sample, whole genome: {dw / GENOME_BASES * 4.25e8 / 3600:.1f} CPU-h;"
          f" full library on it: {tot['cpu_s_wg'] / GENOME_BASES * 4.25e8 / 3600:.1f} CPU-h")
    # appendix A and B
    print("\n## Appendix A\n\n| accession | name | description | hits | CPU-s | share |\n|---|---|---|---|---|---|")
    for r in sorted(a, key=lambda r: r["acc"]):
        print(f"| {r['acc']} | {r['name']} | {r['desc']} | {r['hits']} | {float(r['cpu_s_wg']):.2f} | "
              f"{float(r['cpu_s_wg']) / tot['cpu_s_wg']:.3%} |")
    print("\n## Appendix B\n")
    for k, label, _ in INCLUDE:
        sel = sorted((r for r in rows if r["rule"] == k), key=lambda r: r["acc"])
        print(f"- ({k}) {label}, {len(sel)}: " + ", ".join(f"{r['acc']} {r['name']}" for r in sel))


if __name__ == "__main__":
    if sys.argv[1] == "join":
        join(Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
    elif sys.argv[1] == "rules":
        rules(Path(sys.argv[2]), Path(sys.argv[3]))
    else:
        raise SystemExit(__doc__)

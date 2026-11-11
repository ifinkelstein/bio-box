"""Card 0-09c: compare every finished run with the whole-genome truth and print the tables (run ON THE POD).

    python infernal_report.py <m09c dir> <out.json>

Truth is the producer-equivalent whole-genome scan with default filters (`wgA_default` + `wgB_default`).
Each variant is compared by `infernal_igr.compare` (same model, same strand, overlap >= half the truth
hit). A run that has not finished is left out and named. CPU-hours per 10 Gbp is derived as
cpu_s / genome_bases * 1e10 / 3600, where cpu_s is the measured CPU of the run(s) covering the genome
bases (user + system, timed.py) and genome_bases the bases of the genomes those runs cover (not the
bases searched, so the intergenic fraction is inside the number).

RNA class is a heuristic of this script, from the model's name and its Rfam description line in the
cmsearch output: rRNA, tRNA, housekeeping (tmRNA, RNase P, SRP, 6S), riboswitch or T-box, cis leader
(leader, attenuator, thermometer, UTR), intron, CRISPR repeat, other (mostly sRNAs and antisense RNAs).
"""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import infernal_igr as ig  # noqa: E402

VARIANTS = {  # option: (runs, fasta maps)
    "wg_default": (["wgA_default", "wgB_default"], ["wg_A", "wg_B"]),
    "wg_rfam": (["wgA_rfam", "wgB_rfam"], ["wg_A", "wg_B"]),
    "wg_subset": (["wgA_subset", "wgB_subset"], ["wg_A", "wg_B"]),
    "wg_panel20": (["wgAll_panel20"], ["wg_all"]),
    "igr_m0_default": (["igr_m0_default"], ["igr_m0"]),
    "igr_m50_default": (["igr_m50_default"], ["igr_m50"]),
    "igr_m200_default": (["igr_m200_default"], ["igr_m200"]),
    "wgB_hmmonly": (["wgB_hmmonly"], ["wg_B"]),  # truth restricted to batch B's genomes below
    "igr_m0_notrunc": (["igr_m0_notrunc"], ["igr_m0"]),
    "igr_m200_notrunc": (["igr_m200_notrunc"], ["igr_m200"]),
}
_RES = re.compile(r"\((\d+) residues")


def residues(out_path: Path) -> dict:
    """Residues cmsearch searched over all models (both strands), and those re-searched for truncated hits,
    summed from its per-model pipeline statistics."""
    main = trunc = 0
    with open(out_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("Target sequences:"):
                m = _RES.search(line)
                main += int(m.group(1)) if m else 0
            elif line.startswith("Target sequences re-searched for truncated hits:"):
                m = _RES.search(line)
                trunc += int(m.group(1)) if m else 0
    return {"residues_searched": main, "residues_trunc_researched": trunc}


def rna_class(name: str, desc: str) -> str:
    n, d = name, (desc or "").lower()
    if "rRNA" in n:
        return "rRNA"
    if n.startswith("tRNA"):
        return "tRNA"
    if n in ("tmRNA", "6S") or "RNaseP" in n or "SRP" in n:
        return "housekeeping"
    if "riboswitch" in d or n == "T-box" or "t-box" in d:
        return "riboswitch/T-box"
    if any(w in d for w in ("leader", "attenuator", "thermometer", "utr", "5' element", "motif")) or "leader" in n.lower():
        return "cis leader/motif"
    if "intron" in d or "Intron" in n or "group-II" in n:
        return "intron"
    if "CRISPR" in n:
        return "CRISPR repeat"
    return "other (sRNA etc.)"


def main(mdir: Path, out: Path) -> None:
    runs = {p.stem: json.loads(p.read_text()) for p in sorted((mdir / "runs").glob("*.json"))}
    root = mdir / "root"
    meta = ig._genomes_json(root)
    genome_bases = sum(v["bases"] for v in meta["genomes"].values())
    desc = {}
    for p in sorted((mdir / "out").glob("*.permodel.json")):
        for m in json.loads(p.read_text()):
            if m.get("desc"):
                desc[m["acc"]] = m["desc"]
    res = {"genome_bases": genome_bases, "runs": {}, "options": {}, "missing": []}
    for name, r in runs.items():
        res["runs"][name] = {k: r[k] for k in ("wall_s", "cpu_min", "peak_rss_mib", "exit_code", "started_at")}
        res["runs"][name]["extra_flags"] = [a for a in r["argv"] if a in ("--rfam", "--hmmonly", "--notrunc")]
        o = mdir / "out" / f"{name}.out"
        if o.exists() and r["argv"][0] == "cmsearch":
            res["runs"][name].update(residues(o))
            fa = Path(r["argv"][-1])
            mp = fa.with_name(fa.name[:-4] + ".map.json")
            if mp.exists():
                d = json.loads(mp.read_text())
                res["runs"][name]["records"] = d.get("records", len(d["map"]))
                res["runs"][name]["bases_searched"] = d.get("bases_searched")
    tdir = mdir / "out" / "merged"
    tdir.mkdir(exist_ok=True)
    for opt, (rnames, maps) in VARIANTS.items():
        if not all(n in runs and runs[n]["exit_code"] == 0 and (mdir / "out" / f"{n}.hits.jsonl").exists()
                   for n in rnames):
            res["missing"].append(opt)
            continue
        merged = tdir / f"{opt}.hits.jsonl"
        merged.write_text("".join((mdir / "out" / f"{n}.hits.jsonl").read_text() for n in rnames))
        mdocs = [json.loads((mdir / "fa" / f"{m}.map.json").read_text()) for m in maps]
        searched = sum(d["bases_searched"] for d in mdocs)
        covered = sum(d["genome_bases"] for d in mdocs)
        cpu_s = sum(runs[n]["cpu_min"] for n in rnames) * 60
        wall = sum(runs[n]["wall_s"] for n in rnames)
        smap = "-"
        if opt.startswith("igr"):
            smap = str(mdir / "fa" / f"{maps[0]}.map.json")
        o = {"runs": rnames, "bases_searched": searched, "genome_bases": covered,
             "fraction_searched": round(searched / covered, 4), "cpu_s": round(cpu_s, 1), "wall_s": round(wall, 1),
             "cpu_per_wall": round(cpu_s / wall, 3), "bases_per_cpu_s_genome": round(covered / cpu_s, 1),
             "bases_per_cpu_s_searched": round(searched / cpu_s, 1),
             "cpu_h_per_10gbp": round(cpu_s / covered * 1e10 / 3600, 1),
             "peak_rss_mib": max(runs[n]["peak_rss_mib"] for n in rnames)}
        if opt != "wg_default" and "wg_default" in res["options"]:
            cmp_path = tdir / f"{opt}.compare.json"
            truth = tdir / "wg_default.hits.jsonl"
            if opt == "wg_panel20":  # compare with the truth restricted to the panel's families
                panel = set(ig.g04_module(mdir / "tree").PANEL)
                tp = tdir / "wg_default_panel20.hits.jsonl"
                tp.write_text("".join(x + "\n" for x in truth.read_text().splitlines()
                                      if json.loads(x)["acc"] in panel))
                truth = tp
            if opt == "wgB_hmmonly":
                tp = tdir / "wg_default_B.hits.jsonl"
                tp.write_text((mdir / "out" / "wgB_default.hits.jsonl").read_text())
                truth = tp
            if opt == "wg_subset":
                keys = set((mdir / "lib" / "subset.keys").read_text().split())
                tp = tdir / "wg_default_subset.hits.jsonl"
                tp.write_text("".join(x + "\n" for x in truth.read_text().splitlines()
                                      if json.loads(x)["acc"] in keys))
                truth = tp
            ig.compare(root, truth, merged, smap, cmp_path)
            c = json.loads(cmp_path.read_text())
            for k in ("truth_hits", "variant_hits", "recovered", "recovered_exact_coords", "recall_hits",
                      "precision_hits", "truth_families", "families_recovered", "recall_families",
                      "precision_families", "families_lost_entirely", "lost_by_reason"):
                o[k] = c[k]
            cls = defaultdict(lambda: defaultdict(int))
            for h in c["lost"]:
                cls[rna_class(h["name"], desc.get(h["acc"], ""))][h["why"]] += 1
            o["lost_by_class"] = {k: dict(v) for k, v in sorted(cls.items())}
            tot = defaultdict(int)
            for h in ig.read_hits(truth):
                tot[rna_class(h["name"], desc.get(h["acc"], ""))] += 1
            o["truth_by_class"] = dict(sorted(tot.items()))
            o["extra_by_family"] = {}
            for h in c["extra"]:
                o["extra_by_family"][h["name"]] = o["extra_by_family"].get(h["name"], 0) + 1
            o["lost_families_names"] = [c["per_family"][a]["name"] for a in c["families_lost_entirely"]]
        res["options"][opt] = o
    # per-model cost of the truth runs: by class, the costliest, HMM-only
    pm = defaultdict(float)
    hmmonly = {}
    for n in ("wgA_default", "wgB_default"):
        p = mdir / "out" / f"{n}.permodel.json"
        if p.exists():
            for m in json.loads(p.read_text()):
                pm[(m["acc"], m["name"], m["clen"])] += m["cpu_s"] or 0
                hmmonly[m["acc"]] = m["hmmonly"]
    if pm:
        tot = sum(pm.values())
        res["per_model"] = {
            "models": len(pm), "cpu_s_sum": round(tot, 1), "hmmonly_models": sum(hmmonly.values()),
            "top20": [{"acc": a, "name": nm, "clen": cl, "cpu_s": round(s, 1), "share": round(s / tot, 4)}
                      for (a, nm, cl), s in sorted(pm.items(), key=lambda kv: -kv[1])[:20]],
        }
        keys = (mdir / "lib" / "subset.keys")
        if keys.exists():
            ks = set(keys.read_text().split())
            res["per_model"]["subset_models"] = len(ks)
            res["per_model"]["subset_cpu_s_in_full_run"] = round(sum(s for (a, _, _), s in pm.items() if a in ks), 1)
        panel = set(ig.g04_module(mdir / "tree").PANEL)
        res["per_model"]["panel20_cpu_s_in_full_run"] = round(sum(s for (a, _, _), s in pm.items() if a in panel), 1)
        bycl = defaultdict(float)
        for (a, nm, _), s in pm.items():
            bycl[rna_class(nm, desc.get(a, ""))] += s
        res["per_model"]["cpu_s_by_class"] = {k: round(v, 1) for k, v in sorted(bycl.items())}
        clen_sum = sum(cl for (_, _, cl) in pm)
        res["per_model"]["clen_sum"] = clen_sum
    # spike-ins
    for v in ("default", "rfam", "hmmonly"):
        p = mdir / "out" / f"g04_{v}.score.json"
        if p.exists() and p.read_text().strip():
            res[f"g04_{v}"] = json.loads(p.read_text())
    # the E-01 producer's own hits against this runner's whole-genome default
    e01 = mdir / "out" / "e01_truth.hits.jsonl"
    if e01.exists() and "wg_default" in res["options"]:
        ig.compare(root, e01, tdir / "wg_default.hits.jsonl", "-", tdir / "e01_vs_wg_default.compare.json")
        c = json.loads((tdir / "e01_vs_wg_default.compare.json").read_text())
        res["e01_check"] = {k: c[k] for k in ("truth_hits", "variant_hits", "recovered", "recovered_exact_coords",
                                              "recall_hits", "precision_hits")}
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("runs",)}, indent=1, sort_keys=True))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))

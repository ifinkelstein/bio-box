"""Card C-39 step 4: a nucleotide comparison at the vOTU threshold on samples, with skani from the image.

    python c39_ani.py --inputs IN --s1 S1_FASTA --h1 H1_FASTA --out OUT [--threads 8] [--every 20]

The vOTU threshold (MIUViG; UHGV's and the soil atlas's own): ANI >= 95% over an aligned fraction >= 85% of the shorter
genome. Samples, deterministic: every EVERY-th contig of V1 in FASTA order (the first record of each sequence,
IN/records.tsv), and every EVERY-th record of S1's and of H1's FASTA (the stages' inputs). Two `skani dist` runs,
each sequence its own genome (--qi, --ri), with -c 30 -m 200 (skani's settings for genomes under 20 kb, without
--faster-small, which drops short genomes from the screen), --min-af 0:

  A. the V1 sample against all of S1 and H1: how many V1 vOTU representatives are vOTUs of a phage the root holds;
  B. the V1 sample, the S1 sample and the H1 sample against all of V1: how many S1 and H1 phages fall in a UHGV vOTU,
     and (the calibration) how many V1 representatives meet the threshold with another V1 representative, which
     for vOTUs made at this threshold should be few.

A pair meets the threshold when ANI >= 95 and the aligned fraction of the shorter of the two (by length) >= 85
(rule `shorter`); the summary also counts rule `both` (ANI >= 95 and both aligned fractions >= 85) and rule
`near_identical` (ANI >= 99, both fractions >= 90), and splits the
V1 queries by CheckV tier (records.tsv). --reuse summarises the tables of an earlier run again without skani.
skani prints a pair only past its screen (about 80% identity on sketches); a pair it does not print does not meet
it. Writes OUT/a.tsv and OUT/b.tsv (skani's tables), OUT/summary.json (counts by query set, by V1 length class, and
the CPU seconds of each run, getrusage of the children) and prints the summary. Sample FASTAs and lists in OUT/work.
"""

from __future__ import annotations

import argparse
import csv
import json
import resource
import subprocess
from collections import Counter
from pathlib import Path

from biobox import fasta

ANI_MIN, AF_MIN = 95.0, 85.0


def child_cpu() -> float:
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def length_class(n: int) -> str:
    return "lt_10kb" if n < 10_000 else "10_to_30kb" if n < 30_000 else "30_to_100kb" if n < 100_000 else "ge_100kb"


def write_sample(src: Path, out: Path, every: int, keep=None) -> dict[str, int]:
    """Every `every`-th record of `src` (or of the records whose accession is in `keep`, in file order) to `out`."""
    lengths, k = {}, 0
    with open(out, "w") as w:
        for r in fasta.iter_records(src):
            if keep is not None and r.accession not in keep:
                continue
            k += 1
            if k % every:
                continue
            w.write(f">{r.accession}\n{r.seq.decode()}\n")
            lengths[r.accession] = len(r.seq)
    return lengths


def all_lengths(src: Path) -> dict[str, int]:
    return {r.accession: len(r.seq) for r in fasta.iter_records(src)}


def skani(threads: int, queries: list[Path], refs: list[Path], out: Path) -> float:
    t0 = child_cpu()
    subprocess.run(["skani", "dist", "-t", str(threads), "-c", "30", "-m", "200", "--min-af", "0", "--qi", "--ri",
                    "--short-header", "-q", *map(str, queries), "-r", *map(str, refs), "-o", str(out)], check=True)
    return child_cpu() - t0


def hits(table: Path, length: dict[str, int]) -> dict[str, dict[str, set[str]]]:
    """rule -> query name -> reference names it meets the threshold with (self pairs dropped). Rules: `shorter`
    (the vOTU threshold: ANI and the shorter genome's aligned fraction) and `both` (ANI and both fractions >= AF_MIN:
    the two genomes are each nearly covered, not one a fragment of the other) and `near_identical` (ANI >= 99 and both
    fractions >= 90: the same virus within skani's error)."""
    out: dict[str, dict[str, set[str]]] = {"shorter": {}, "both": {}, "near_identical": {}}
    with open(table, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, ref = r["Query_name"], r["Ref_name"]
            if q == ref:
                continue
            ani, afr, afq = float(r["ANI"]), float(r["Align_fraction_ref"]), float(r["Align_fraction_query"])
            if ani < ANI_MIN:
                continue
            af_short = afq if length[q] <= length[ref] else afr
            if af_short >= AF_MIN:
                out["shorter"].setdefault(q, set()).add(ref)
            if min(afq, afr) >= AF_MIN:
                out["both"].setdefault(q, set()).add(ref)
            if ani >= 99.0 and min(afq, afr) >= 90.0:
                out["near_identical"].setdefault(q, set()).add(ref)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True, type=Path)
    ap.add_argument("--s1", required=True, type=Path)
    ap.add_argument("--h1", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--every", type=int, default=20)
    ap.add_argument("--reuse", action="store_true", help="summarise a.tsv and b.tsv of an earlier run again")
    a = ap.parse_args()
    work = a.out / "work"
    work.mkdir(parents=True, exist_ok=True)
    v1_fa = a.inputs / "uhgv_v1.0_votus_full.fna"
    with open(a.inputs / "repeated_records.tsv", newline="") as fh:
        repeated = {r["genome"] for r in csv.DictReader(fh, delimiter="\t")}
    v1_len = all_lengths(v1_fa)
    first = {g for g in v1_len if g not in repeated}
    s1_len, h1_len = all_lengths(a.s1), all_lengths(a.h1)
    if set(s1_len) & set(v1_len) or set(h1_len) & set(v1_len) or set(s1_len) & set(h1_len):
        raise SystemExit("record names are shared between the sets; names would be ambiguous")
    length = {**v1_len, **s1_len, **h1_len}
    sv = write_sample(v1_fa, work / "v1_sample.fna", a.every, keep=first)
    ss = write_sample(a.s1, work / "s1_sample.fna", a.every)
    sh = write_sample(a.h1, work / "h1_sample.fna", a.every)
    res: dict = {"every": a.every, "threads": a.threads, "threshold": {"ani_min": ANI_MIN, "af_shorter_min": AF_MIN},
                 "samples": {"v1": len(sv), "s1": len(ss), "h1": len(sh)},
                 "sample_bases": {"v1": sum(sv.values()), "s1": sum(ss.values()), "h1": sum(sh.values())},
                 "reference_records": {"v1": len(v1_len), "s1": len(s1_len), "h1": len(h1_len)}}
    if a.reuse:
        res["cpu_s"] = json.loads((a.out / "summary.json").read_text())["cpu_s"]
    else:
        cpu_a = skani(a.threads, [work / "v1_sample.fna"], [a.s1, a.h1], a.out / "a.tsv")
        cpu_b = skani(a.threads, [work / "v1_sample.fna", work / "s1_sample.fna", work / "h1_sample.fna"], [v1_fa],
                      a.out / "b.tsv")
        res["cpu_s"] = {"a": round(cpu_a, 1), "b": round(cpu_b, 1)}
    ha, hb = hits(a.out / "a.tsv", length), hits(a.out / "b.tsv", length)
    with open(a.inputs / "records.tsv", newline="") as fh:
        quality = {r["genome"]: r["checkv_quality"] for r in csv.DictReader(fh, delimiter="\t")}

    def summary(sample: dict[str, int], hh: dict[str, dict[str, set[str]]], ref_names: set[str]) -> dict:
        res = {"queries": len(sample)}
        by = Counter(length_class(n) for n in sample.values())
        byq = Counter(quality.get(q, "-") for q in sample)
        for rule, h in hh.items():
            met = {q for q in sample if h.get(q, set()) & ref_names}
            by_met = Counter(length_class(sample[q]) for q in met)
            byq_met = Counter(quality.get(q, "-") for q in met)
            res[rule] = {"meet": len(met),
                         "by_length": {c: {"queries": by[c], "meet": by_met[c]} for c in sorted(by)},
                         "by_checkv": {c: {"queries": byq[c], "meet": byq_met[c]} for c in sorted(byq)}}
        return res

    s1n, h1n, v1n = set(s1_len), set(h1_len), set(v1_len)
    res["a_v1_sample_vs_s1"] = summary(sv, ha, s1n)
    res["a_v1_sample_vs_h1"] = summary(sv, ha, h1n)
    res["a_v1_sample_vs_s1_or_h1"] = summary(sv, ha, s1n | h1n)
    res["b_s1_sample_vs_v1"] = summary(ss, hb, v1n)
    res["b_h1_sample_vs_v1"] = summary(sh, hb, v1n)
    res["b_v1_sample_vs_other_v1"] = summary(sv, hb, v1n)
    res["b_v1_examples_both"] = sorted((q, sorted(hb["both"][q])[:3]) for q in sv if q in hb["both"])[:20]
    (a.out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, sort_keys=True))


if __name__ == "__main__":
    main()

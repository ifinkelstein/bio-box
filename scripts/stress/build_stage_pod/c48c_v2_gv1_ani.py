"""Card C-48c (part 2) step 4: stage V2 at the vOTU threshold against S1, H1 and V1, on samples, with skani from the
image (card C-39's method, `c39_ani.py`, whose rules and skani settings are reused).

    python c48c_v2_gv1_ani.py --v2 V2_INPUTS --s1 S1_FASTA --h1 H1_FASTA --v1 V1_FASTA --out OUT [--threads 3]
        [--every 20]

The vOTU threshold (MIUViG; IMG/VR's own): ANI >= 95% over an aligned fraction >= 85% of the shorter genome. Samples,
deterministic: every EVERY-th record of V2's two FASTA objects taken in order (IMG/VR then Inovirus), skipping a record
whose contig an earlier record has (repeated_records.tsv) and the records of multi-record genomes (GVMAG bins, compared
contig by contig they would mislead); every EVERY-th record of S1's, H1's and V1's FASTA. Two `skani dist` runs, each
record its own genome (--qi, --ri), -c 30 -m 200 --min-af 0 (`c39_ani.skani`):

  A. the V2 sample against all of S1, H1 and V1: how many V2 representatives are vOTUs of a virus the root or V1 holds;
  B. the V2, S1, H1 and V1 samples against all of V2: how many S1, H1 and V1 viruses fall in an IMG/VR high-quality
     vOTU, and (the calibration) how many V2 representatives meet the threshold with another V2 representative, which
     for vOTUs made at this threshold should be few.

Rules `shorter`, `both`, `near_identical` as `c39_ani.hits`. Writes OUT/a.tsv, OUT/b.tsv, OUT/summary.json (counts by
query set, by length class and by V2 source; CPU seconds of each run, getrusage of the children).
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import c39_ani  # noqa: E402  (beside this script)
from biobox import fasta


def v2_sample(v2: Path, every: int, out: Path) -> tuple[dict[str, int], dict[str, str], dict[str, int]]:
    """The V2 sample written to `out`; its lengths, each sampled record's source, and every V2 record's length."""
    with open(v2 / "repeated_records.tsv", newline="") as fh:
        repeated = {r["genome"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(v2 / "records.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    per_genome = Counter(r["genome"] for r in rows)
    objects = json.loads((v2 / "objects.json").read_text())
    lengths, source, sample, k = {}, {}, {}, 0
    seen: set[str] = set()
    with open(out, "w") as w:
        for o in objects:
            for r in fasta.iter_records(Path(o["path"])):
                name = r.accession
                lengths[name] = len(r.seq)
                g = name.split("|", 1)[0]
                if g in repeated or per_genome.get(g, 1) > 1 or g in seen:
                    continue
                seen.add(g)
                k += 1
                if k % every:
                    continue
                w.write(f">{name}\n{r.seq.decode()}\n")
                sample[name] = len(r.seq)
                source[name] = o["source_name"]
    return sample, source, lengths


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v2", required=True, type=Path)
    ap.add_argument("--s1", required=True, type=Path)
    ap.add_argument("--h1", required=True, type=Path)
    ap.add_argument("--v1", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=3)
    ap.add_argument("--every", type=int, default=20)
    ap.add_argument("--reuse", action="store_true", help="summarise a.tsv and b.tsv of an earlier run again")
    a = ap.parse_args()
    work = a.out / "work"
    work.mkdir(parents=True, exist_ok=True)
    sv, src, v2_len = v2_sample(a.v2, a.every, work / "v2_sample.fna")
    s1_len, h1_len, v1_len = (c39_ani.all_lengths(p) for p in (a.s1, a.h1, a.v1))
    names = [set(v2_len), set(s1_len), set(h1_len), set(v1_len)]
    if sum(len(n) for n in names) != len(set().union(*names)):
        raise SystemExit("record names are shared between the sets; names would be ambiguous")
    length = {**v2_len, **s1_len, **h1_len, **v1_len}
    ss = c39_ani.write_sample(a.s1, work / "s1_sample.fna", a.every)
    sh = c39_ani.write_sample(a.h1, work / "h1_sample.fna", a.every)
    s1v = c39_ani.write_sample(a.v1, work / "v1_sample.fna", a.every)
    v2_objects = [Path(o["path"]) for o in json.loads((a.v2 / "objects.json").read_text())]
    res: dict = {"every": a.every, "threads": a.threads,
                 "threshold": {"ani_min": c39_ani.ANI_MIN, "af_shorter_min": c39_ani.AF_MIN},
                 "samples": {"v2": len(sv), "s1": len(ss), "h1": len(sh), "v1": len(s1v)},
                 "sample_bases": {"v2": sum(sv.values()), "s1": sum(ss.values()), "h1": sum(sh.values()),
                                  "v1": sum(s1v.values())},
                 "reference_records": {"v2": len(v2_len), "s1": len(s1_len), "h1": len(h1_len), "v1": len(v1_len)}}
    if a.reuse:
        res["cpu_s"] = json.loads((a.out / "summary.json").read_text())["cpu_s"]
    else:
        cpu_a = c39_ani.skani(a.threads, [work / "v2_sample.fna"], [a.s1, a.h1, a.v1], a.out / "a.tsv")
        cpu_b = c39_ani.skani(a.threads, [work / "v2_sample.fna", work / "s1_sample.fna", work / "h1_sample.fna",
                                          work / "v1_sample.fna"], v2_objects, a.out / "b.tsv")
        res["cpu_s"] = {"a": round(cpu_a, 1), "b": round(cpu_b, 1)}
    ha, hb = c39_ani.hits(a.out / "a.tsv", length), c39_ani.hits(a.out / "b.tsv", length)

    def summary(sample: dict[str, int], hh, ref_names: set[str]) -> dict:
        out = {"queries": len(sample)}
        by = Counter(c39_ani.length_class(n) for n in sample.values())
        bys = Counter(src.get(q, "-") for q in sample)
        for rule, h in hh.items():
            met = {q for q in sample if (h.get(q, set()) - {q}) & ref_names}
            out[rule] = {"meet": len(met),
                         "by_length": {c: {"queries": by[c], "meet": sum(1 for q in met if c39_ani.length_class(
                             sample[q]) == c)} for c in sorted(by)},
                         "by_source": {s: {"queries": bys[s], "meet": sum(1 for q in met if src.get(q, "-") == s)}
                                       for s in sorted(bys)}}
        return out

    s1n, h1n, v1n, v2n = set(s1_len), set(h1_len), set(v1_len), set(v2_len)
    res["a_v2_sample_vs_s1"] = summary(sv, ha, s1n)
    res["a_v2_sample_vs_h1"] = summary(sv, ha, h1n)
    res["a_v2_sample_vs_v1"] = summary(sv, ha, v1n)
    res["a_v2_sample_vs_s1_h1_or_v1"] = summary(sv, ha, s1n | h1n | v1n)
    res["b_s1_sample_vs_v2"] = summary(ss, hb, v2n)
    res["b_h1_sample_vs_v2"] = summary(sh, hb, v2n)
    res["b_v1_sample_vs_v2"] = summary(s1v, hb, v2n)
    res["b_v2_sample_vs_other_v2"] = summary(sv, hb, v2n)
    res["b_v2_examples_both"] = sorted((q, sorted(hb["both"][q] - {q})[:3]) for q in sv
                                       if hb["both"].get(q, set()) - {q})[:20]
    (a.out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, sort_keys=True))


if __name__ == "__main__":
    main()

"""Card C-35: the estimates of release r26.11.1's steps before its launch, and the check of each published run record
against them (AGENTS §6: a measurement more than 20% from its estimate stops the card; DEC-43: only steps estimated at
24 CPU-minutes or more are judged).

    python c35_est.py before --root ROOT --release R --stages S1,...,PL1 --earlier rM.2 --calib CALIB.json
                             --parts 24 --profile P [--jobs J --cpus C] [--tmp D]        (runtime image; reads only)
    python3 -I c35_est.py rates --root ROOT --release R --est BEFORE.json                 (host; standard library)

`before` (every number printed is labelled: the code's functions give *estimate*s, counts read from the root are
*measured*, the scaling of the earlier release's records is *derived*):
- the proteins of the release: rows of each stage's `protein_seq` files (the files its gene-call run record lists;
  `stage.release_proteins` sums the same rows) and the distinct protein uids of the release and of the earlier
  release's stages (one DuckDB `count(DISTINCT)` each over those explicit files, under a memory limit, spilling to
  --tmp); the genes (`feature_protein` rows of the same records) of both;
- the code's decisions: the one-process estimate `stage.cluster_cpu_min(rows)` and whether the release clusters in
  parts (`stage._cluster_split`'s rule: over `SHARD_MARGIN` of the budget, or `--cluster-parts`); the plan's cut
  (`l2cluster.cut_of`, as build-stage records it with `--cluster-parts` and the profile's `cluster_cpu_min_max`); for
  every release of the root with cluster plan records, why its clustering is or is not taken (card B-49:
  `l2cluster._source_problem` on this release's plan inputs, `l2cluster.plan_inputs`, which read no file's content);
- the earlier release's clustering call by call (CALIB, from `c35_calib.py`: each call's parts, the cut's estimates
  and the measured run records) scaled to this release: the model's own growth (`l2cluster.call_cpu_min` at both
  protein counts, the call's units taken in proportion to the distinct proteins, and a prefilter part's own table in
  proportion to them too), times the earlier release's measured cost of the call. Per call: the projected total, its
  largest part, and the calibration factor (measured / the cut's estimate on the earlier release) that `rates` applies
  to this release's cut estimates. Every other release step: the earlier record times the protein ratio (the
  clustering steps), or the gene ratio (tokens, postings, dim, manifest, catalog);
- memory (the largest part's peak RSS of each call, scaled by the protein ratio, times the parts at once) and the
  clustering scratch (the earlier release's peak, given as --scratch-peak-gb, scaled the same way).

`rates`: every run record of the release that exists, against its estimate: a part `<call>.part-<p>` (or a call cut
into one part) at its cut record's estimate (`l2/cluster_plan/release=R/<call>.cut.json`) times the call's calibration
factor; any other step at `before`'s projection. Prints one JSON line (the judged rows outside 20%, the worst ratio, the
sums); exit 3 when a step estimated at 24 CPU-minutes or more is outside 20% of its estimate, else 0.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

TOL = 0.20
# DEC-43 (PI, 2026-11-15): a step estimated under this many CPU-minutes is reported but not judged.
MIN_JUDGED_CPU_MIN = 24.0
SPLIT_STEPS = ("tokens", "postings", "dim", "manifest", "catalog")  # scaled by genes; the rest by proteins


def gene_call_outputs(root: Path, stage: str, table: str, gcv: str) -> list[Path]:
    rec = json.loads((root / "runs" / f"stage={stage}" / "genecall.json").read_text())
    prefix = f"{table}/gene_call_version={gcv}/stage={stage}/"
    return [root / o["path"] for o in rec["outputs"] if o["path"].startswith(prefix) and o["path"].endswith(".parquet")]


def before(a) -> int:
    import duckdb
    import pyarrow.parquet as pq

    from biobox import profile as bprofile
    from biobox import stage as bs
    from biobox.l1.genecall import GENE_CALL_VERSION as GCV
    from biobox.l2 import cluster as l2

    root, stages = Path(a.root), a.stages.split(",")
    prof = bprofile.load(a.profile)
    calib = json.loads(Path(a.calib).read_text())
    t = pq.read_table(root / "release_def" / f"release={a.earlier}" / "part-0.parquet", partitioning=None).to_pylist()
    estages = [x["stage"] for x in sorted(t, key=lambda x: x["position"])]
    if estages != a.earlier_stages.split(","):
        raise SystemExit(f"release {a.earlier} is over {estages}, not {a.earlier_stages}")
    per_stage, files_of = {}, {}
    for s in stages:
        ps = gene_call_outputs(root, s, "l0seq/protein_seq", GCV)
        fp = gene_call_outputs(root, s, "l1/feature_protein", GCV)
        files_of[s] = ps
        per_stage[s] = {"protein_seq_rows": sum(pq.read_metadata(p).num_rows for p in ps),
                        "genes": sum(pq.read_metadata(p).num_rows for p in fp)}
    rows, missing = bs.release_proteins(root, stages, GCV)
    if missing is not None or rows != sum(v["protein_seq_rows"] for v in per_stage.values()):
        raise SystemExit(f"release_proteins: {rows} rows, stage without a record {missing}")
    Path(a.tmp).mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(config={"threads": a.duck_threads, "memory_limit": "6GB", "temp_directory": str(a.tmp),
                                 "preserve_insertion_order": False})

    def distinct(ss: list[str]) -> int:
        lit = "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for s in ss for p in files_of[s]) + "]"
        return con.execute(f"SELECT count(DISTINCT protein_uid) FROM read_parquet({lit}, hive_partitioning = false)"
                           ).fetchone()[0]

    n_all = distinct(stages)
    n_old = distinct(estages)
    new_by_stage = {}
    for s in stages:
        if s in estages:
            continue
        lit_new = "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for p in files_of[s]) + "]"
        lit_old = "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for t in estages for p in files_of[t]) + "]"
        new_by_stage[s] = con.execute(
            f"SELECT count(*) FROM (SELECT protein_uid FROM read_parquet({lit_new}, hive_partitioning = false) "
            f"EXCEPT SELECT protein_uid FROM read_parquet({lit_old}, hive_partitioning = false))").fetchone()[0]
    con.close()
    r = n_all / n_old
    genes_all = sum(per_stage[s]["genes"] for s in stages)
    genes_old = sum(per_stage[s]["genes"] for s in estages)
    rg = genes_all / genes_old
    t_one = bs.cluster_cpu_min(rows)
    rows_old = sum(per_stage[s]["protein_seq_rows"] for s in estages)
    split = a.parts is not None or t_one > bs.SHARD_MARGIN * prof.cpu_min_per_process
    cut = l2.cut_of(a.parts, None if a.parts is not None else float(repr(bs.SHARD_MARGIN * prof.cpu_min_per_process)),
                    None if prof.cluster_cpu_min_max is None else float(prof.cluster_cpu_min_max))
    inputs = l2.plan_inputs(root, a.release, GCV, stages, cut)
    reuse = {s: (l2._source_problem(root, a.release, s, inputs) or "qualifies") for s in l2._plan_releases(root)}
    # the earlier release's calls scaled to this one
    calls, part_max, total_proj, factors = [], 0.0, 0.0, {}
    for c in calib["calls"]:
        name = c["call"]
        rows_ = calib["part_rows"][name]
        u_old = sum(x["units"] or 0 for x in rows_)
        o_old = c["overhead_estimate"] or 0.0
        k = len(rows_) if a.parts is None else a.parts
        m_old = l2.call_cpu_min(name, n_old, u_old) + len(rows_) * o_old
        m_new = l2.call_cpu_min(name, n_all, round(u_old * r)) + k * o_old * r
        g = m_new / m_old if m_old else r
        proj = c["measured_sum"] * g
        largest = (c["max_part_cpu_min"] or 0) * g * len(rows_) / k
        factors[name] = round(c["measured_sum"] / c["estimate_sum"], 6) if c["estimate_sum"] else None
        part_max = max(part_max, largest)
        total_proj += proj
        calls.append({"call": name, "parts": k, "measured_earlier": c["measured_sum"], "growth": round(g, 4),
                      "projected": round(proj, 1), "largest_part_projected": round(largest, 1),
                      "calibration_factor": factors[name],
                      "peak_rss_gb_projected": round((c["max_peak_rss_gb"] or 0) * r, 1)})
    steps = {}
    for o in calib["others"]:
        scale = rg if o["step"] in SPLIT_STEPS else r
        steps[o["step"]] = round(o["cpu_min"] * scale, 3)
    others_total = sum(v for k_, v in steps.items() if k_.startswith("cluster"))
    mem = sorted(((c["peak_rss_gb_projected"], c["call"]) for c in calls), reverse=True)
    doc = {
        "release": a.release, "stages": stages, "earlier": a.earlier, "earlier_stages": estages,
        "measured": {"per_stage": per_stage, "protein_seq_rows": rows, "protein_seq_rows_earlier": rows_old,
                     "distinct_proteins": n_all, "distinct_proteins_earlier": n_old,
                     "distinct_new_not_in_earlier": n_all - n_old, "new_by_stage_not_in_earlier": new_by_stage,
                     "genes": genes_all, "genes_earlier": genes_old,
                     "label": "measured: Parquet footers and DuckDB count(DISTINCT) over the gene-call records' files"},
        "code": {"cluster_cpu_min_one_process": round(t_one, 1),
                 "cluster_cpu_min_one_process_earlier_rows": round(bs.cluster_cpu_min(rows_old), 1),
                 "margin_cpu_min": bs.SHARD_MARGIN * prof.cpu_min_per_process, "clusters_in_parts": split,
                 "cut": cut, "reuse": reuse, "label": "estimate: the code's functions"},
        "projection": {
            "protein_ratio": round(r, 6), "gene_ratio": round(rg, 6),
            "earlier_cluster_cpu_min": calib["cluster_cpu_min"],
            "cluster_parts_cpu_min": round(total_proj, 1), "cluster_other_steps_cpu_min": round(others_total, 1),
            "cluster_cpu_min": round(total_proj + others_total, 1),
            "cluster_cpu_min_by_one_process_ratio": round(calib["cluster_cpu_min"] * t_one / bs.cluster_cpu_min(rows_old), 1),
            "largest_part_cpu_min": round(part_max, 1), "every_part_under_240": part_max < 240,
            "calls": calls, "steps": steps, "factors": factors,
            "memory": {"largest_parts_gb": mem[:3], "jobs": a.jobs, "cpus": a.cpus,
                       "worst_parts_at_once_gb": round(mem[0][0] * a.jobs, 1) if a.jobs and mem else None},
            "scratch_peak_gb": round(a.scratch_peak_gb * r, 1) if a.scratch_peak_gb else None,
            "label": "derived: the earlier release's run records times the model's growth (clustering calls), "
                     "the protein ratio (other clustering steps) or the gene ratio (tokens, postings, dim, manifest, "
                     "catalog)"},
    }
    print(json.dumps(doc, indent=1, sort_keys=True))
    return 0


def rates(a) -> int:
    root, rel = Path(a.root), a.release
    est = json.loads(Path(a.est).read_text())
    proj = est["projection"]
    factors, steps = proj["factors"], proj["steps"]
    runs = root / "runs" / f"stage={rel}"
    plans = root / "l2" / "cluster_plan" / f"release={rel}"
    rows, bad, unknown = [], [], []
    cuts: dict[str, dict] = {}
    for f in sorted(runs.glob("*.json")) if runs.is_dir() else []:
        if f.name.endswith(".publishing.json"):
            continue
        name = f.name[: -len(".json")]
        try:
            r = json.loads(f.read_text())
        except (OSError, ValueError):
            continue
        e, basis = None, None
        call, _, part = name.partition(".part-")
        if part:
            if call not in cuts:
                cp = plans / f"{call}.cut.json"
                cuts[call] = json.loads(cp.read_text()) if cp.exists() else {}
            parts = cuts[call].get("parts", [])
            if parts and factors.get(call):
                e = parts[int(part)]["cpu_min"] * factors[call]
                basis = f"cut estimate {parts[int(part)]['cpu_min']:.3f} x factor {factors[call]}"
        elif (plans / f"{name}.cut.json").exists() and name not in cuts:
            cut = json.loads((plans / f"{name}.cut.json").read_text())
            if len(cut.get("parts", [])) == 1 and factors.get(name):  # a call cut into one part: the call itself
                e = cut["parts"][0]["cpu_min"] * factors[name]
                basis = f"one part: cut estimate x factor {factors[name]}"
        if e is None and name in steps:
            e, basis = steps[name], "earlier record x ratio"
        if e is None and name.endswith(".cut"):
            e, basis = steps.get(name), "earlier record x ratio"
        row = {"step": name, "cpu_min": round(r["cpu_min"], 3), "wall_s": r.get("wall_s"),
               "peak_rss_gb": round((r.get("peak_rss_bytes") or 0) / 1e9, 2)}
        if e is None:
            unknown.append(name)
        else:
            ratio = r["cpu_min"] / e if e > 0 else float("inf")
            judged = e >= MIN_JUDGED_CPU_MIN
            row.update({"estimate": round(e, 3), "ratio": round(ratio, 3), "judged": judged, "basis": basis})
            if judged and abs(ratio - 1) > TOL:
                bad.append(row)
        rows.append(row)
    judged = [x for x in rows if x.get("judged")]
    out = {"release": rel, "records": len(rows), "judged": len(judged),
           "cpu_min": round(sum(x["cpu_min"] for x in rows), 2),
           "estimate_of_recorded": round(sum(x.get("estimate", 0) for x in rows), 2),
           "worst": max(judged, key=lambda x: abs(x["ratio"] - 1), default=None),
           "largest": max(rows, key=lambda x: x["cpu_min"], default=None),
           "outside_20pct": bad, "no_estimate": unknown[-8:], "no_estimate_n": len(unknown)}
    if a.out:
        Path(a.out).write_text(json.dumps({**out, "rows": rows}, indent=1, sort_keys=True))
    print(json.dumps(out, sort_keys=True))
    return 3 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("before", "rates"))
    ap.add_argument("--root", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--stages")
    ap.add_argument("--earlier", default="rM.2")
    ap.add_argument("--earlier-stages", default="S1,P1,P2,P3,P4,P5,P6")
    ap.add_argument("--calib")
    ap.add_argument("--parts", type=int, default=None)
    ap.add_argument("--profile")
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--cpus", type=int, default=None)
    ap.add_argument("--scratch-peak-gb", type=float, default=None)
    ap.add_argument("--tmp", default="/tmp/if2453/c35/tmp/est")
    ap.add_argument("--duck-threads", type=int, default=2)
    ap.add_argument("--est")
    ap.add_argument("--out")
    a = ap.parse_args()
    return {"before": before, "rates": rates}[a.mode](a)


if __name__ == "__main__":
    sys.exit(main())

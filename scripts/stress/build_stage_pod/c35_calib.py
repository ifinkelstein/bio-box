"""Card C-35: what an earlier release's clustering in parts and its other release steps cost, call by call, against the
estimates its cut records held. Host side, standard library only (JSON); reads only:

    python3 -I c35_calib.py --root ROOT --release rM.2 [--out FILE.json]

For each clustering call: the parts (from `<call>.cut.json` under `l2/cluster_plan/release=<R>/`), the cut's estimate of
each part (`parts[].cpu_min`: the model's work share plus the per-part overhead), the run record of each part
(`runs/stage=<R>/<call>.part-<p>.json`, or the call's own record when it is one part) with its CPU-minutes, wall
seconds, threads and peak RSS; the sum, the largest part and the ratio measured / estimate of the call. Then every
other run record of the release (the plan, the unsplit calls, the joins, `cluster`, tokens, postings, dim, manifest,
catalog) with its CPU-minutes and peak. Prints one JSON document (and writes it to --out).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def load(p: Path) -> dict | None:
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return None


def rec_row(r: dict) -> dict:
    return {"cpu_min": round(r["cpu_min"], 3), "wall_s": round(r["wall_s"], 1), "threads": r.get("threads"),
            "peak_rss_gb": round((r.get("peak_rss_bytes") or 0) / 1e9, 2), "started_at": r.get("started_at")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--release", required=True)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    runs = a.root / "runs" / f"stage={a.release}"
    plans = a.root / "l2" / "cluster_plan" / f"release={a.release}"
    seen: set[str] = set()
    calls = []
    for cutf in sorted(plans.glob("*.cut.json")):
        call = cutf.name[: -len(".cut.json")]
        cut = load(cutf)
        parts = cut["parts"]
        rows = []
        for p, part in enumerate(parts):
            name = f"{call}.part-{p:05d}" if len(parts) > 1 else call
            r = load(runs / f"{name}.json")
            seen.add(name)
            row = {"part": p, "estimate": round(part["cpu_min"], 3), "units": part.get("units")}
            if r is not None:
                row.update(rec_row(r))
                row["ratio"] = round(r["cpu_min"] / part["cpu_min"], 3) if part["cpu_min"] else None
            rows.append(row)
        meas = [x["cpu_min"] for x in rows if "cpu_min" in x]
        est = sum(x["estimate"] for x in rows)
        join = load(runs / f"{call}.json") if len(parts) > 1 else None  # also listed under `others`
        calls.append({"call": call, "parts": len(parts), "work_estimate": cut.get("cpu_min"),
                      "overhead_estimate": cut.get("overhead_cpu_min"), "one_process": "one_process" in cut,
                      "estimate_sum": round(est, 2), "measured_sum": round(sum(meas), 2),
                      "ratio": round(sum(meas) / est, 3) if est and meas else None,
                      "max_part_cpu_min": max(meas) if meas else None,
                      "max_part_estimate": max(x["estimate"] for x in rows),
                      "min_part_ratio": min((x["ratio"] for x in rows if x.get("ratio") is not None), default=None),
                      "max_part_ratio": max((x["ratio"] for x in rows if x.get("ratio") is not None), default=None),
                      "max_peak_rss_gb": max((x["peak_rss_gb"] for x in rows if "peak_rss_gb" in x), default=None),
                      "threads": sorted({x["threads"] for x in rows if x.get("threads")}),
                      "cut_record_cpu_min": round(load(runs / f"{call}.cut.json")["cpu_min"], 3)
                      if (runs / f"{call}.cut.json").exists() else None,
                      "join": rec_row(join) if join else None, "rows": rows})
    others = []
    for f in sorted(runs.glob("*.json")):
        name = f.name[: -len(".json")]
        if name.endswith(".publishing") or name in seen:
            continue
        r = load(f)
        if r is not None:
            others.append({"step": name, **rec_row(r)})
    total = sum(x["measured_sum"] for x in calls)
    allrec = []
    for f in runs.glob("*.json"):
        if not f.name.endswith(".publishing.json"):
            r = load(f)
            if r is not None:
                allrec.append((f.name, r["cpu_min"], r.get("peak_rss_bytes") or 0))
    clu = [x for x in allrec if x[0].startswith("cluster")]
    doc = {"release": a.release, "records": len(allrec),
           "cluster_records": len(clu), "cluster_cpu_min": round(sum(x[1] for x in clu), 2),
           "all_cpu_min": round(sum(x[1] for x in allrec), 2),
           "largest": sorted(((n, round(c, 2)) for n, c, _ in allrec), key=lambda t: -t[1])[:5],
           "largest_peak_gb": sorted(((n, round(p / 1e9, 2)) for n, _, p in allrec), key=lambda t: -t[1])[:5],
           "parts_measured_cpu_min": round(total, 2),
           "calls": [{k: v for k, v in c.items() if k != "rows"} for c in calls],
           "others": others, "part_rows": {c["call"]: c["rows"] for c in calls}}
    text = json.dumps(doc, indent=1, sort_keys=True)
    if a.out:
        a.out.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())

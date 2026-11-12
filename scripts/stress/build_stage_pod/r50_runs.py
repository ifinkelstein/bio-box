"""Run records of one stage or release, as a table and totals: r50_runs.py ROOT STAGE [--json]. Standard library only.

50 Gbp build (planning/measurements/2026-11-11_build_r50.md). Reads ROOT/runs/stage=STAGE/*.json (not the intent
files). Prints one line per record in start order (step, start UTC, wall s, CPU-minutes, peak RSS MiB, threads,
over_cpu_budget, cpu_budget_reason), then totals: records, CPU-minutes, the largest CPU-minutes and peak RSS with
their steps, the wall time from the first start to the last end, and, for clustering in parts (card B-41), the
number of module-call steps, of `.cut` steps and of `.part-` steps, and the parts of each call that was cut into
more than one.
"""
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

root, stage = Path(sys.argv[1]), sys.argv[2]
recs = []
for f in sorted((root / "runs" / f"stage={stage}").glob("*.json")):
    if f.name.endswith(".publishing.json"):
        continue
    r = json.loads(f.read_text())
    recs.append(r)
recs.sort(key=lambda r: r["started_at"])
rows = []
for r in recs:
    rows.append({"step": r["step"], "started_at": r["started_at"], "wall_s": r["wall_s"], "cpu_min": r["cpu_min"],
                 "peak_rss_mib": None if r.get("peak_rss_bytes") is None else round(r["peak_rss_bytes"] / 2**20),
                 "threads": r["threads"], "over": r.get("over_cpu_budget"), "reason": r.get("cpu_budget_reason")})
if "--json" in sys.argv:
    print(json.dumps(rows))
else:
    for x in rows:
        print(f"{x['step']:50s} {x['started_at'][:19]} {x['wall_s']:10.1f} {x['cpu_min']:9.2f} "
              f"{x['peak_rss_mib'] if x['peak_rss_mib'] is not None else '-':>8} {x['threads']:3d} "
              f"{x['over']} {x['reason'] or ''}")


def t(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


tot = {"records": len(rows), "cpu_min": round(sum(x["cpu_min"] for x in rows), 2)}
if rows:
    big = max(rows, key=lambda x: x["cpu_min"])
    tot["max_cpu_min"] = [big["step"], big["cpu_min"]]
    mem = max(rows, key=lambda x: x["peak_rss_mib"] or 0)
    tot["max_peak_rss_mib"] = [mem["step"], mem["peak_rss_mib"]]
    first = min(t(x["started_at"]) for x in rows)
    last = max(t(x["started_at"]) + timedelta(seconds=x["wall_s"]) for x in rows)
    tot["first_start"], tot["last_end"] = first.isoformat(), last.isoformat()
    tot["span_s"] = round((last - first).total_seconds(), 1)
    tot["over_300"] = [[x["step"], x["cpu_min"]] for x in rows if x["cpu_min"] > 300]
cl = [x["step"] for x in rows if x["step"].startswith("cluster.")]
if cl:
    parts: dict[str, int] = {}
    for s in cl:
        if ".part-" in s:
            call = s.split(".part-")[0]
            parts[call] = parts.get(call, 0) + 1
    tot["cluster_steps"] = len(cl)
    tot["cluster_cut_steps"] = sum(1 for s in cl if s.endswith(".cut"))
    tot["cluster_part_steps"] = sum(parts.values())
    tot["calls_in_parts"] = parts
print("totals", json.dumps(tot))

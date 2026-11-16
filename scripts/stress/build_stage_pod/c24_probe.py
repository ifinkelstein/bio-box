"""Card C-24: read-only figures of the built stages for the estimates made before stage H1 is built. Host side,
standard library only (python3 -I), reads only:

    python3 -I c24_probe.py ROOT [S1,P1,P2,P3,...]
    python3 -I c24_probe.py --guard GUARD.JSONL

For each stage: the run records of its stage steps (`runs/stage=<S>/<step>.json`: CPU-minutes, wall seconds, peak RSS,
threads), the Pfam 38.2 scan's shards against its plan (`l1/domain_scan/library=pfam-38.2/stage=<S>/plan.json`: each
shard's residues, so CPU-minutes per million residues per shard), and the bytes of the root's files under a
`stage=<S>` directory. Prints one JSON document.

With --guard: the 50 Gbp build's guard lines (`r50_watch.py`, once a minute), the largest scratch seen while only stage
steps ran (no `release ...` command live), per stage named in the live steps, and over all lines.
"""
import json
import os
import sys
from pathlib import Path

LIB = "pfam-38.2"


def rec(p: Path) -> dict | None:
    try:
        r = json.loads(p.read_text())
    except (OSError, ValueError):
        return None
    return {k: r.get(k) for k in ("cpu_min", "wall_s", "peak_rss_bytes", "threads", "started_at", "image_digest")}


def stage_bytes(root: Path) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
    for top, _, files in os.walk(root):
        stages = [part[6:] for part in Path(top).relative_to(root).parts if part.startswith("stage=")]
        if not stages:
            continue
        s = stages[0]
        n = out.setdefault(s, [0, 0])
        for f in files:
            try:
                n[0] += os.lstat(os.path.join(top, f)).st_size
                n[1] += 1
            except OSError:
                pass
    return out


def guard(path: Path) -> int:
    peak: dict[str, list] = {}
    top = [0.0, None]
    for line in path.read_text().splitlines():
        try:
            g = json.loads(line)
        except ValueError:
            continue
        sc = g.get("scratch_gb") or 0.0
        if sc > top[0]:
            top = [sc, g.get("time")]
        live = [x["step"] for x in g.get("live", [])]
        if not live or any(s.startswith("release") for s in live):
            continue
        for s in {part[6:] for x in live for part in x.split() if part.startswith("stage=")}:
            if sc > peak.get(s, [0.0])[0]:
                peak[s] = [sc, g.get("time"), live]
    print(json.dumps({"stage_steps_only_scratch_peak_gb": peak, "all_lines_scratch_peak_gb": top}, indent=1))
    return 0


def main() -> int:
    if sys.argv[1] == "--guard":
        return guard(Path(sys.argv[2]))
    root = Path(sys.argv[1])
    stages = sys.argv[2].split(",") if len(sys.argv) > 2 else ["S1", "P1", "P2", "P3", "P4", "P5", "P6"]
    sizes = stage_bytes(root)
    doc = {"root": str(root), "stages": {}}
    for s in stages:
        runs = root / "runs" / f"stage={s}"
        d: dict = {"bytes_under_stage": sizes.get(s, [0, 0])[0], "files_under_stage": sizes.get(s, [0, 0])[1],
                   "steps": {}}
        shards = []
        for p in sorted(runs.glob("*.json")) if runs.exists() else []:
            if p.name.endswith(".publishing.json"):
                continue
            name = p.name[:-5]
            r = rec(p)
            if name.startswith(f"domains-{LIB}.shard-"):
                shards.append((int(name.rsplit("-", 1)[1]), r))
            else:
                d["steps"][name] = r
        plan_p = root / "l1" / "domain_scan" / f"library={LIB}" / f"stage={s}" / "plan.json"
        if plan_p.exists():
            plan = json.loads(plan_p.read_text())
            per = []
            for k, r in sorted(shards):
                res = plan["shards"][k]["residues"]
                per.append({"shard": k, "residues": res, "proteins": plan["shards"][k]["proteins"],
                            "cpu_min": r["cpu_min"], "wall_s": r["wall_s"], "threads": r["threads"],
                            "peak_rss_bytes": r["peak_rss_bytes"],
                            "cpu_min_per_m_residues": round(r["cpu_min"] / (res / 1e6), 3)})
            rates = [x["cpu_min_per_m_residues"] for x in per]
            d["scan"] = {"proteins": plan["proteins"], "residues": plan["residues"], "shards": len(plan["shards"]),
                         "shard_residues": plan["shard_residues"],
                         "cpu_min_shards": round(sum(x["cpu_min"] for x in per), 2),
                         "cpu_min_per_m_residues_all": round(sum(x["cpu_min"] for x in per) / (plan["residues"] / 1e6), 3)
                         if per else None,
                         "rate_min": min(rates) if rates else None, "rate_max": max(rates) if rates else None,
                         "cpu_min_max_shard": max((x["cpu_min"] for x in per), default=None),
                         "peak_rss_max_shard": max((x["peak_rss_bytes"] for x in per), default=None)}
        doc["stages"][s] = d
    print(json.dumps(doc, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

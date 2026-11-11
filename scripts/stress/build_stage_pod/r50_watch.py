"""Guard of the 50 Gbp build (planning/measurements/2026-11-11_build_r50.md). Host side, standard library only.

    python3 -I r50_watch.py --root R --scratch D --since FILE --state DIR [--unit-prefix r50-]
                            [--tmp-min-gb 200] [--scratch-max-gb 300] [--cpu-max 300]

Prints one JSON line and exits 0 when the build may go on, 3 when it must stop (the reasons are in the line):

- `tmp_free_gb`: free space of the filesystem of /tmp; a reason when under --tmp-min-gb.
- `scratch_gb`: bytes under D (the profile's scratch: MMseqs2's work directories, scan and ingest runs); a reason
  when over --scratch-max-gb. The largest value seen is kept in DIR/scratch_peak (bytes and the time).
- `records`: run records under R/runs written after FILE's mtime (this driver's), each with its CPU-minutes; a
  reason for each over --cpu-max (AGENTS §7.2).
- `live`: each step command running in a systemd unit whose name starts with --unit-prefix (a process whose
  command line runs `-m biobox.cli stage|release ...`, `-m biobox.l1.domains` or `-m biobox.l0seq.store`, counted
  once at the topmost such process), with the CPU-minutes of its whole process tree so far (utime + stime of each
  live process, plus cutime + cstime for the children each has reaped); a reason for each over --cpu-max.

Reads /proc and the root; writes only DIR/scratch_peak.
"""
import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

TICK = os.sysconf("SC_CLK_TCK")
MODULES = ("biobox.l1.domains", "biobox.l0seq.store")


def tree_bytes(d: Path) -> int:
    total = 0
    for top, _, files in os.walk(d):
        for f in files:
            try:
                total += os.lstat(os.path.join(top, f)).st_size
            except OSError:
                pass  # removed while walking
    return total


def procs() -> dict[int, dict]:
    out = {}
    uid = os.getuid()
    for p in Path("/proc").iterdir():
        if not p.name.isdigit():
            continue
        try:
            if p.stat().st_uid != uid:
                continue
            stat = (p / "stat").read_text()
            fields = stat[stat.rindex(")") + 2:].split()
            cmd = (p / "cmdline").read_bytes().split(b"\0")
            cg = (p / "cgroup").read_text()
        except OSError:
            continue
        out[int(p.name)] = {"ppid": int(fields[1]),
                            "ticks": sum(int(x) for x in fields[11:15]),  # utime stime cutime cstime
                            "cmd": [c.decode(errors="replace") for c in cmd if c],
                            "cgroup": cg}
    return out


def step_of(cmd: list[str]) -> str | None:
    """The step a command line runs: `<sub> <action-or-step>` of biobox.cli, or the module, else None."""
    for i, a in enumerate(cmd[:-1]):
        if a == "-m" and cmd[i + 1] in MODULES:
            return " ".join(cmd[i + 1:i + 3])
        if a == "-m" and cmd[i + 1] == "biobox.cli" and i + 2 < len(cmd) and cmd[i + 2] in ("stage", "release"):
            rest = cmd[i + 2:]
            name = rest[1] if len(rest) > 1 else ""
            if name == "verify":
                return None
            tag = []
            for flag in ("--stage", "--release", "--name", "--shard", "--batch"):
                if flag in rest and rest.index(flag) + 1 < len(rest):
                    tag.append(f"{flag[2:]}={rest[rest.index(flag) + 1]}")
            action = rest[2] if len(rest) > 2 and not rest[2].startswith("--") else ""
            return " ".join(x for x in (rest[0], name, action, *tag) if x)
    return None


def live(prefix: str) -> list[dict]:
    ps = procs()
    kids: dict[int, list[int]] = {}
    for pid, p in ps.items():
        kids.setdefault(p["ppid"], []).append(pid)

    def ours(p: dict) -> bool:
        return any(f"/{prefix}" in line for line in p["cgroup"].splitlines())

    def matched_above(pid: int) -> bool:
        q = ps[pid]["ppid"]
        while q in ps:
            if step_of(ps[q]["cmd"]) is not None:
                return True
            q = ps[q]["ppid"]
        return False

    res = []
    for pid, p in ps.items():
        name = step_of(p["cmd"])
        if name is None or not ours(p) or matched_above(pid):
            continue
        ticks, stack = 0, [pid]
        while stack:
            q = stack.pop()
            ticks += ps[q]["ticks"]
            stack.extend(kids.get(q, []))
        res.append({"pid": pid, "step": name, "cpu_min": round(ticks / TICK / 60, 2)})
    return sorted(res, key=lambda r: -r["cpu_min"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--scratch", required=True)
    ap.add_argument("--since", required=True)
    ap.add_argument("--state", required=True)
    ap.add_argument("--unit-prefix", default="r50-")
    ap.add_argument("--tmp-min-gb", type=float, default=200)
    ap.add_argument("--scratch-max-gb", type=float, default=300)
    ap.add_argument("--cpu-max", type=float, default=300)
    a = ap.parse_args()
    reasons = []
    out = {"time": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    free = shutil.disk_usage("/tmp").free / 1e9
    out["tmp_free_gb"] = round(free, 1)
    if free < a.tmp_min_gb:
        reasons.append(f"/tmp has {free:.1f} GB free, under {a.tmp_min_gb:g}")
    sb = tree_bytes(Path(a.scratch)) if Path(a.scratch).exists() else 0
    out["scratch_gb"] = round(sb / 1e9, 2)
    peak_file = Path(a.state) / "scratch_peak"
    try:
        peak = json.loads(peak_file.read_text())
    except (OSError, ValueError):
        peak = {"bytes": 0}
    if sb > peak["bytes"]:
        peak = {"bytes": sb, "time": out["time"]}
        tmp = peak_file.with_name(peak_file.name + ".tmp")
        tmp.write_text(json.dumps(peak))
        os.replace(tmp, peak_file)
    out["scratch_peak_gb"] = round(peak["bytes"] / 1e9, 2)
    if sb > a.scratch_max_gb * 1e9:
        reasons.append(f"scratch {a.scratch} holds {sb / 1e9:.1f} GB, over {a.scratch_max_gb:g}")
    since = Path(a.since).stat().st_mtime
    recs = []
    for f in Path(a.root, "runs").glob("stage=*/*.json"):
        if f.name.endswith(".publishing.json"):
            continue
        try:
            if f.stat().st_mtime < since:
                continue
            r = json.loads(f.read_text())
        except (OSError, ValueError):
            continue
        cm = r.get("cpu_min")
        if cm is not None and cm > a.cpu_max:
            reasons.append(f"run record {f.parent.name}/{f.name}: {cm:.1f} CPU-minutes, over {a.cpu_max:g}")
        recs.append(cm or 0)
    out["records"] = {"n": len(recs), "max_cpu_min": round(max(recs), 2) if recs else None}
    lv = live(a.unit_prefix)
    out["live"] = lv[:4]
    for r in lv:
        if r["cpu_min"] > a.cpu_max:
            reasons.append(f"step {r['step']} (pid {r['pid']}) has used {r['cpu_min']:.1f} CPU-minutes, over {a.cpu_max:g}")
    out["reasons"] = reasons
    print(json.dumps(out), flush=True)
    return 3 if reasons else 0


if __name__ == "__main__":
    sys.exit(main())

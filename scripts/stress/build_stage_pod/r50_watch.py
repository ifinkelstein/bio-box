"""Guard of the 50 Gbp build (planning/measurements/2026-11-11_build_r50.md). Host side, standard library only.

    python3 -I r50_watch.py --root R --scratch D --since FILE --state DIR [--unit-prefix r50-]
                            [--tmp-min-gb 200] [--scratch-max-gb 300] [--cpu-max 300] [--cluster-cpu-max N]
                            [--stages S,REL,...]

--cluster-cpu-max is the limit for clustering steps (default --cpu-max); --stages limits the run records read to
those stages and releases.

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

Clustering's exemptions (DEC-29, DEC-77) raise a step's limit to the one-process cap only where build-stage grants
them. A live shared-index `<call>.cut` whose cut record is not yet written gets DEC-77's cap only when build-stage's
own test holds (card R-05, `stage._cluster_in_parts`): the shared `createindex_cpu_min` over the call's target in the
plan's work directory (under the process's `--tmp`, else D) is over the plan's part budget and within its
one-process cap.

Reads /proc, the root and the target index of such a cut; writes only DIR/scratch_peak and DIR/createindex_counts.json
(the counts of that index, kept while the index is unchanged).
"""
import argparse
import hashlib
import json
import math
import os
import shutil
import stat
import sys
import time
from pathlib import Path

TICK = os.sysconf("SC_CLK_TCK")
MODULES = ("biobox.l1.domains", "biobox.l0seq.store")
DEC29_EXEMPT = "exempt: DEC-29"
DEC77_EXEMPT = "exempt: DEC-77"

# Standard-library mirrors of biobox.l2.cluster (the watch runs on the host with `python3 -I`); card R-05's tests hold
# each equal to the library.
DEFAULT_CLUSTERING_MODE = "cascade"  # l2.cluster.DEFAULT_CLUSTERING_MODE: a plan without "method" (PLAN_FORMAT 1)
PLAN_DIR = "l2/cluster_plan"  # l2.cluster.PLAN_DIR
PLAN_STEP = "cluster.plan"
WORK_DIR = "biobox-cluster"  # l2.cluster.WORK_DIR, under the scratch directory
PLAN_INPUT_KEYS_V1 = ("format", "release", "stages", "gene_call_version", "mmseqs2", "files", "cut")
PLAN_INPUT_KEYS = (*PLAN_INPUT_KEYS_V1, "method")  # l2.cluster._INPUT_KEYS (PLAN_FORMAT 2)
# each method's shared-index prefilter cut steps: (the call's target database in the work directory, its -s), as
# l2.cluster.calls(mode) makes them (`c.split and c.shares_index`; `c.positional[1]`, `-s`)
_CASCADE_NR50_CUTS = {
    "cluster.nr50.15.prefilter.cut": ("tmp50/w/input_step1", "2"),
    "cluster.nr50.19.prefilter.cut": ("tmp50/w/input_step2", "3"),
}
SHARED_INDEX_CUTS: dict[str, dict[str, tuple[str, str]]] = {
    "cascade": dict(_CASCADE_NR50_CUTS),
    "hybrid-a": dict(_CASCADE_NR50_CUTS),
    "hybrid-b1": {"cluster.nr50.15.prefilter.cut": _CASCADE_NR50_CUTS["cluster.nr50.15.prefilter.cut"]},
    "linclust": {},
}
SHARED_INDEX_CUT_STEPS = frozenset(s for cuts in SHARED_INDEX_CUTS.values() for s in cuts)
# l2.cluster.createindex_cpu_min (card B-68): CPU-minutes per 10^9 residues at 10M targets, by -s
CREATEINDEX_T0 = 10_000_000
CREATEINDEX_GROWTH = 0.48
CREATEINDEX_RATE_BY_SENS: dict[str, float] = {"2": 11.80, "3": 15.75}
INDEX_COUNTS = "createindex_counts.json"  # under DIR


def createindex_cpu_min(residues: int, targets: int, sensitivity: str = "2") -> float:
    """l2.cluster.createindex_cpu_min: the estimate of one shared-index `mmseqs createindex` `.cut` step."""
    if residues <= 0 or targets <= 0:
        return 0.0
    s = str(sensitivity)
    rate = CREATEINDEX_RATE_BY_SENS.get(s)
    if rate is None:
        try:
            s_val = float(s)
            rate = CREATEINDEX_RATE_BY_SENS["2"] * (s_val / 2.0) if s_val > 0 else CREATEINDEX_RATE_BY_SENS["2"]
        except ValueError:
            rate = CREATEINDEX_RATE_BY_SENS["2"]
    return rate * 1e-9 * residues * (targets / CREATEINDEX_T0) ** CREATEINDEX_GROWTH


def _safe_segment(s: str) -> bool:
    return bool(s) and "/" not in s and "\\" not in s and ".." not in s


def _valid_cap(mx: object, default_cap: float) -> float:
    if isinstance(mx, bool) or not isinstance(mx, (int, float)):
        return default_cap
    mx_f = float(mx)
    return mx_f if math.isfinite(mx_f) and mx_f > default_cap else default_cap


def dec29_cap(root: Path, release: str, call: str, default_cap: float) -> float:
    """Return DEC-29's one-process cap for a single clustering call from its cut record, else default_cap."""
    if (
        not _safe_segment(release)
        or not _safe_segment(call)
        or not call.startswith("cluster.")
        or ".part-" in call
        or call.endswith(".cut")
    ):
        return default_cap
    cut_path = root / "l2" / "cluster_plan" / f"release={release}" / f"{call}.cut.json"
    try:
        cut = json.loads(cut_path.read_text())
    except (OSError, ValueError):
        return default_cap
    if not isinstance(cut, dict):
        return default_cap
    parts = cut.get("parts")
    op = cut.get("one_process")
    if not isinstance(parts, list) or len(parts) != 1 or not isinstance(op, dict):
        return default_cap
    if op.get("decision") != "DEC-29":
        return default_cap
    return _valid_cap(op.get("cpu_min_max"), default_cap)


def _number(x: object) -> float | None:
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return None
    x_f = float(x)
    return x_f if math.isfinite(x_f) else None


def _json_bytes(doc: dict) -> bytes:
    return json.dumps(doc, sort_keys=True, indent=1).encode() + b"\n"  # l2.cluster._json_bytes


def work_dir(tmp: Path, root: Path, plan: dict) -> Path | None:
    """The plan's MMseqs2 work directory under the scratch directory `tmp`, as l2.cluster._work_dir names it from the
    root and the plan's input fields (`_plan_input_dict`: without "method" for PLAN_FORMAT 1); None for a plan
    without its input fields or with an unsafe release."""
    keys = PLAN_INPUT_KEYS_V1 if plan.get("format") == 1 else PLAN_INPUT_KEYS
    if any(k not in plan for k in keys):
        return None
    inputs = {k: plan[k] for k in keys}
    rel = inputs["release"]
    if not isinstance(rel, str) or not _safe_segment(rel):
        return None
    try:
        key = hashlib.sha256(f"{Path(root).resolve()}\n".encode() + _json_bytes(inputs)).hexdigest()[:16]
    except (TypeError, ValueError, OSError):
        return None
    return Path(tmp) / WORK_DIR / f"{rel}-{key}"


def _index_ident(st: os.stat_result, index: Path) -> dict:
    return {"path": str(index), "dev": st.st_dev, "ino": st.st_ino, "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def index_counts(index: Path, cache: Path | None = None) -> tuple[int, int] | None:
    """(entries, residues) of the MMseqs2 sequence database whose index is `index`, counted as
    l2.cluster.cut_createindex_estimate counts them (an entry per non-blank line; its length less the newline and the
    NUL), or None when it is missing or not an index of three tab-separated integers. The counts are kept in `cache`
    for that file while its device, inode, size and mtime are unchanged; a cache that does not match is not used."""
    try:
        st = index.stat()
    except OSError:
        return None
    if not stat.S_ISREG(st.st_mode):
        return None
    ident = _index_ident(st, index)
    if cache is not None:
        try:
            memo = json.loads(cache.read_text())
        except (OSError, ValueError):
            memo = None
        if isinstance(memo, dict) and memo.get("index") == ident:
            t, u = memo.get("targets"), memo.get("residues")
            if all(isinstance(v, int) and not isinstance(v, bool) and v >= 0 for v in (t, u)):
                return t, u
    targets = residues = 0
    try:
        with open(index, "rb") as fh:
            for line in fh:
                if not line.strip():
                    continue
                f = line.split(b"\t")
                if len(f) != 3:
                    return None
                k, o, n = int(f[0]), int(f[1]), int(f[2])
                if k < 0 or o < 0 or n < 0:
                    return None
                targets += 1
                residues += max(0, n - 2)
        st2 = index.stat()
    except (OSError, ValueError):
        return None
    if cache is not None and _index_ident(st2, index) == ident:
        tmp = cache.with_name(cache.name + ".tmp")
        try:
            tmp.write_text(json.dumps({"index": ident, "targets": targets, "residues": residues}))
            os.replace(tmp, cache)
        except OSError:
            pass  # a cache that cannot be written is counted again next time
    return targets, residues


def _live_cut_cap(root: Path, release: str, step: str, default_cap: float, tmp: Path | None,
                  cache: Path | None) -> float:
    """DEC-77's cap for a live shared-index `.cut` whose cut record is not yet written, by build-stage's own test
    (`stage._cluster_in_parts` with `l2.cluster.cut_createindex_estimate`): the plan clusters by a part budget
    (`cut.part_cpu_min`) and holds a one-process cap (`cut.one_process_cpu_min_max`, the profile's
    cluster_cpu_min_max), `step` is a shared-index cut of the plan's method, and the createindex estimate over its
    target in the plan's work directory under `tmp` is over the budget (else build-stage runs it without the
    exemption) and within the cap (else build-stage refuses to launch it); never for a plan of PLAN_FORMAT 1.
    Else default_cap."""
    if tmp is None:
        return default_cap
    try:
        plan = json.loads((root / PLAN_DIR / f"release={release}" / f"{PLAN_STEP}.json").read_text())
    except (OSError, ValueError):
        return default_cap
    if not isinstance(plan, dict):
        return default_cap
    if plan.get("format") == 1:
        # card R-02 (l2.cluster._format_1): a PLAN_FORMAT 1 cut keeps the pre-B-71 record and refusal, so
        # cut_createindex_estimate is None for it and build-stage never grants it exempt: DEC-77
        return default_cap
    method = plan.get("method", DEFAULT_CLUSTERING_MODE)
    cuts = SHARED_INDEX_CUTS.get(method) if isinstance(method, str) else None
    plan_cut = plan.get("cut")
    if not cuts or step not in cuts or not isinstance(plan_cut, dict):
        return default_cap
    budget = _number(plan_cut.get("part_cpu_min"))
    cap = _valid_cap(plan_cut.get("one_process_cpu_min_max"), default_cap)
    if budget is None or budget <= 0 or cap <= default_cap:
        return default_cap
    w = work_dir(tmp, root, plan)
    if w is None:
        return default_cap
    target, sens = cuts[step]
    counts = index_counts(Path(f"{w / target}.index"), cache)
    if counts is None:
        return default_cap
    est = createindex_cpu_min(counts[1], counts[0], sens)
    return cap if budget < est <= cap else default_cap


def dec77_cut_cap(root: Path, release: str, step: str, default_cap: float, *, live_proc: bool = False,
                  tmp: Path | None = None, cache: Path | None = None) -> float:
    """Return DEC-77's one-process cap for a shared-index `.cut` (`createindex`) step, else default_cap.

    When `<step>.json` exists under `l2/cluster_plan/release=<release>/` (a published cut, or a live redo of one), it
    decides: `len(parts) > 1` and `one_process.decision == "DEC-77"`. A live `.cut` process (`live_proc=True`) whose
    `<step>.json` is not yet written gets the cap only by build-stage's own test (`_live_cut_cap`, card R-05) over the
    work directory under `tmp` (the process's scratch directory); `cache` keeps the target index's counts.
    """
    if not _safe_segment(release) or not _safe_segment(step) or step not in SHARED_INDEX_CUT_STEPS:
        return default_cap
    pdir = root / PLAN_DIR / f"release={release}"
    cut_path = pdir / f"{step}.json"
    if cut_path.exists():
        try:
            cut = json.loads(cut_path.read_text())
        except (OSError, ValueError):
            return default_cap
        if not isinstance(cut, dict):
            return default_cap
        parts = cut.get("parts")
        op = cut.get("one_process")
        if not isinstance(parts, list) or len(parts) <= 1 or not isinstance(op, dict):
            return default_cap
        if op.get("decision") != "DEC-77":
            return default_cap
        return _valid_cap(op.get("cpu_min_max"), default_cap)
    if not live_proc:
        return default_cap
    return _live_cut_cap(root, release, step, default_cap, tmp, cache)


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
        row = {"pid": pid, "step": name, "cpu_min": round(ticks / TICK / 60, 2)}
        for flag in ("--root", "--tmp"):  # where the step itself finds its plan and its work directory (card R-05)
            v = _flag_value(p["cmd"], flag)
            if v and os.path.isabs(v):  # build-stage passes both resolved; a relative one is the process's cwd's
                row[flag[2:]] = v
        res.append(row)
    return sorted(res, key=lambda r: -r["cpu_min"])


def _flag_value(cmd: list[str], flag: str) -> str | None:
    """The value after the first `flag` of a command line, else None."""
    return cmd[cmd.index(flag) + 1] if flag in cmd and cmd.index(flag) + 1 < len(cmd) else None


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
    ap.add_argument("--cluster-cpu-max", type=float, default=None,
                    help="the limit for clustering steps (`cluster.*` records, `release cluster` commands); "
                         "default --cpu-max")
    ap.add_argument("--stages", default=None,
                    help="comma list: read only run records of these stages and releases (runs/stage=<S>/)")
    a = ap.parse_args()
    only = None if a.stages is None else {f"stage={s}" for s in a.stages.split(",")}
    ccap = a.cpu_max if a.cluster_cpu_max is None else a.cluster_cpu_max
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
        if f.name.endswith(".publishing.json") or (only is not None and f.parent.name not in only):
            continue
        try:
            if f.stat().st_mtime < since:
                continue
            r = json.loads(f.read_text())
        except (OSError, ValueError):
            continue
        cm = r.get("cpu_min")
        cap = ccap if f.name.startswith("cluster") else a.cpu_max
        if f.parent.name.startswith("stage=") and f.stem.startswith("cluster."):
            rel = f.parent.name[len("stage="):]
            reason = r.get("cpu_budget_reason")
            if reason == DEC29_EXEMPT:
                cap = dec29_cap(Path(a.root), rel, f.stem, ccap)
            elif reason == DEC77_EXEMPT:
                cap = dec77_cut_cap(Path(a.root), rel, f.stem, ccap, live_proc=False)
        if cm is not None and cm > cap:
            reasons.append(f"run record {f.parent.name}/{f.name}: {cm:.1f} CPU-minutes, over {cap:g}")
        recs.append(cm or 0)
    out["records"] = {"n": len(recs), "max_cpu_min": round(max(recs), 2) if recs else None}
    lv = live(a.unit_prefix)
    out["live"] = [{k: r[k] for k in ("pid", "step", "cpu_min") if k in r} for r in lv[:4]]
    for r in lv:
        step = r["step"]
        cap = ccap if step.startswith("release cluster") else a.cpu_max
        # an exemption only raises the cap, so it is looked up only for a step over the default (card R-05)
        if r["cpu_min"] > cap and step.startswith(("release cluster step ", "release cluster redo ")):
            kv = dict(tok.split("=", 1) for tok in step.split()[3:] if "=" in tok)
            rel, call = kv.get("release", ""), kv.get("name", "")
            if rel and call:
                if call.endswith(".cut"):
                    cap = dec77_cut_cap(Path(r.get("root") or a.root), rel, call, ccap, live_proc=True,
                                        tmp=Path(r.get("tmp") or a.scratch), cache=Path(a.state) / INDEX_COUNTS)
                else:
                    cap = dec29_cap(Path(a.root), rel, call, ccap)
        if r["cpu_min"] > cap:
            reasons.append(f"step {step} (pid {r['pid']}) has used {r['cpu_min']:.1f} CPU-minutes, over {cap:g}")
    out["reasons"] = reasons
    print(json.dumps(out), flush=True)
    return 3 if reasons else 0


if __name__ == "__main__":
    sys.exit(main())

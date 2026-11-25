"""Card C-46: release-train planner (DEC-52, DEC-54, C-38a, B-49).

Computes the next rolling release name (`r<YY>.<MM>.<n>`), selects the ordered stage list
(preserving the latest defined release's stage order and appending newly complete stages at
depth `called` or `annotated` in `inputs.json` timestamp order, excluding `stored`-only and
incomplete stages with explicit reasons via `biobox.stage.stage_incomplete`), estimates
clustering or detects `B-49` reuse, and prints the thread plan by the hour (`biobox.runstep`).

Read-only on the data root; writes nothing.

    python train_plan.py --root ROOT --profile PROFILE [--date YYYY-MM-DD] [--release R]
                         [--parts K] [--same-stages] [--json] [--no-live-cpu]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import date, datetime, timezone
from pathlib import Path

from biobox import catalog, profile as bprofile, runstep, stage as bs
from biobox.l1 import ncrna
from biobox.l1.genecall import GENE_CALL_VERSION as GCV
from biobox.l2 import cluster as l2

RELEASE_STEP_FILES = {
    "cluster.json",
    "cluster.plan.json",
    "tokens.json",
    "postings.json",
    "dim.json",
    "manifest.json",
    "catalog.json",
}
_DOMAIN_STEP_RE = re.compile(r"^domains-([^-]+)-(.+?)(?:\.plan|\.shard-\d+|\.publishing)?\.json$")


class ReleaseTrainError(ValueError):
    """Refusal condition for the release train (exits with code 2)."""


def defined_releases(root: Path) -> list[str]:
    """Releases with a valid `release_def/release=<R>/part-0.parquet`, in `catalog.release_order`."""
    d = Path(root) / "release_def"
    if not d.is_dir():
        return []
    names = []
    for p in d.iterdir():
        if p.is_dir() and p.name.startswith("release="):
            rel = p.name.removeprefix("release=")
            if (p / "part-0.parquet").is_file():
                names.append(rel)
    return sorted(names, key=catalog.release_order)


def existing_releases(root: Path) -> set[str]:
    """Union of release names present in `release_def/` (any dir) or `manifests/`."""
    out: set[str] = set()
    rd = Path(root) / "release_def"
    if rd.is_dir():
        for p in rd.iterdir():
            if p.is_dir() and p.name.startswith("release="):
                out.add(p.name.removeprefix("release="))
    out.update(catalog.manifest_releases(root))
    return out


def next_release_name(root: Path, on_date: date | None = None) -> str:
    """DEC-52 (1): `r<YY>.<MM>.<n>`, counted from 1 within `(YY, MM)` across existing releases."""
    d = date.today() if on_date is None else on_date
    yy, mm = d.year % 100, d.month
    max_n = 0
    for name in existing_releases(root):
        m = catalog.ROLLING_NAME.fullmatch(name)
        if m is not None and int(m[1]) == yy and int(m[2]) == mm:
            max_n = max(max_n, int(m[3]))
    return f"r{yy:02d}.{mm:02d}.{max_n + 1}"


def resolve_release_name(root: Path, on_date: date | None = None, explicit: str | None = None) -> str:
    """Return the next release name, or validate an explicit `--release` (refusing `latest` or existing)."""
    if explicit is None:
        return next_release_name(root, on_date)
    if explicit == catalog.LATEST:
        raise ReleaseTrainError("refused: release name 'latest' is reserved (card W-37)")
    if explicit in existing_releases(root):
        raise ReleaseTrainError(f"refused: release {explicit} already exists under {root}")
    return explicit


def _arg_after(argv: list[str], flag: str) -> str | None:
    try:
        i = argv.index(flag)
    except ValueError:
        return None
    return argv[i + 1] if i + 1 < len(argv) else None


def _parse_ts(ts: str | None) -> float | None:
    if not ts:
        return None
    try:
        dt = datetime.fromisoformat(ts.replace("Z", "+00:00"))
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt.timestamp()
    except ValueError:
        return None


def _read_json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def _is_release_run_dir(root: Path, name: str, rel_names: set[str]) -> bool:
    if name in rel_names:
        return True
    if (Path(root) / "l2" / "cluster_plan" / f"release={name}").is_dir():
        return True
    rdir = Path(root) / "runs" / f"stage={name}"
    if not rdir.is_dir():
        return False
    files = {p.name for p in rdir.iterdir() if p.is_file()}
    if "inputs.json" not in files and (files & RELEASE_STEP_FILES):
        return True
    return False


def infer_stage_spec(root: Path, stage: str) -> dict:
    """Infer a stage's built depth and build-stage options from its run records under `runs/stage=<stage>/`."""
    rdir = Path(root) / "runs" / f"stage={stage}"
    files = sorted(p for p in rdir.iterdir() if p.is_file() and p.name.endswith(".json")) if rdir.is_dir() else []
    names = [p.name for p in files]

    inputs_rec = _read_json(rdir / "inputs.json") or {}
    ingest_rec = _read_json(rdir / "ingest.json") or {}
    inputs_argv = inputs_rec.get("argv") or []
    ingest_argv = ingest_rec.get("argv") or []
    gcv = _arg_after(inputs_argv, "--gene-call-version") or _arg_after(ingest_argv, "--gene-call-version") or GCV
    genomes = _arg_after(inputs_argv, "--genomes")
    objects = _arg_after(ingest_argv, "--objects")

    ts_str = inputs_rec.get("started_at") or inputs_rec.get("finished_at")
    epoch = _parse_ts(ts_str)
    if epoch is None:
        inputs_p = rdir / "inputs.json"
        epoch = inputs_p.stat().st_mtime if inputs_p.is_file() else float("inf")

    lib_map: dict[tuple[str, str], dict] = {}
    for p in files:
        m = _DOMAIN_STEP_RE.match(p.name)
        if not m:
            continue
        rec = _read_json(p) if not p.name.endswith(".publishing.json") else None
        argv = (rec or {}).get("argv") or []
        lname = _arg_after(argv, "--library") or m.group(1)
        lver = _arg_after(argv, "--library-version") or m.group(2)
        hmm = _arg_after(argv, "--hmm")
        cutoff = _arg_after(argv, "--cutoff") or "gathering"
        key = (lname, lver)
        prev = lib_map.get(key)
        if prev is None or (prev.get("hmm") is None and hmm is not None):
            lib_map[key] = {"name": lname, "version": lver, "hmm": hmm, "cutoff": cutoff}

    crispr_enabled = any(n.startswith("crispr-") for n in names)

    ncrna_cm: str | None = None
    ncrna_set: str = ncrna.DEFAULT_SET
    for p in files:
        if p.name.startswith("ncrna-") and not p.name.endswith(".publishing.json"):
            argv = (_read_json(p) or {}).get("argv") or []
            ncrna_cm = _arg_after(argv, "--cm") or ncrna_cm
            ncrna_set = _arg_after(argv, "--ncrna-set") or ncrna_set

    non_stored_markers = {
        n for n in names
        if n not in ("inputs.json", "ingest.json", "labels.json")
        and not n.startswith("inputs.")
        and not n.startswith("ingest.")
        and not n.startswith("labels.")
    }
    if lib_map or ncrna_cm is not None:
        depth = "annotated"
    elif non_stored_markers:
        depth = "called"
    elif (rdir / "inputs.json").is_file() and (rdir / "ingest.json").is_file():
        depth = "stored"
    else:
        depth = "called"

    libraries = sorted(lib_map)
    return {
        "stage": stage,
        "depth": depth,
        "libraries": libraries,
        "library_specs": [lib_map[k] for k in libraries],
        "crispr": crispr_enabled,
        "ncrna_cm": ncrna_cm,
        "ncrna_set": ncrna_set,
        "gene_call_version": gcv,
        "genomes": genomes,
        "objects": objects,
        "inputs_started_at": ts_str,
        "inputs_epoch": epoch,
    }


def select_stages(root: Path) -> dict:
    """DEC-52 (2) + DEC-54 (2): preserve latest release's stage order, append newly complete stages."""
    root = Path(root)
    defs = defined_releases(root)
    latest = defs[-1] if defs else None
    prior_stages = catalog.release_stages(root, latest) if latest is not None else []
    prior_set = set(prior_stages)
    rel_names = existing_releases(root)

    runs_root = root / "runs"
    stage_dirs = sorted(
        p.name.removeprefix("stage=")
        for p in runs_root.iterdir()
        if p.is_dir() and p.name.startswith("stage=")
    ) if runs_root.is_dir() else []

    candidates: list[tuple[float, str, dict]] = []
    excluded: list[dict] = []

    for s in stage_dirs:
        if s in prior_set or _is_release_run_dir(root, s, rel_names):
            continue
        spec = infer_stage_spec(root, s)
        if spec["depth"] == "stored":
            problems = bs.stage_incomplete(root, s, "stored", gene_call_version=spec["gene_call_version"])
            reason = "stored only (depth=stored; DEC-54 requires depth >= called)"
            if problems:
                reason += "; " + "; ".join(problems)
            excluded.append({"stage": s, "depth": "stored", "problems": problems, "reason": reason})
            continue
        problems = bs.stage_incomplete(
            root,
            s,
            spec["depth"],
            libraries=spec["libraries"],
            crispr=spec["crispr"],
            ncrna_cm=Path(spec["ncrna_cm"]) if spec["ncrna_cm"] else None,
            ncrna_set=spec["ncrna_set"],
            gene_call_version=spec["gene_call_version"],
        )
        if problems:
            excluded.append({
                "stage": s,
                "depth": spec["depth"],
                "problems": problems,
                "reason": "; ".join(problems),
            })
        else:
            candidates.append((spec["inputs_epoch"], s, spec))

    candidates.sort(key=lambda item: (item[0], item[1]))
    new_stages = [c[1] for c in candidates]
    stages = prior_stages + new_stages
    last_stage = infer_stage_spec(root, stages[-1]) if stages else None
    return {
        "latest_release": latest,
        "prior_stages": prior_stages,
        "new_stages": new_stages,
        "stages": stages,
        "excluded": excluded,
        "nothing_new": len(new_stages) == 0,
        "last_stage": last_stage,
        "earlier_manifest_releases": catalog.manifest_releases(root),
    }


def estimate_and_reuse(root: Path, release: str, stages: list[str], prof: bprofile.Profile,
                       parts: int | None = None, *,
                       cluster_method: str = l2.DEFAULT_CLUSTERING_MODE) -> dict:
    """Clustering estimate from protein count and B-49 reuse detection (read-only)."""
    l2._check_mode(cluster_method)
    if not stages:
        return {
            "protein_seq_rows": 0,
            "cluster_cpu_min_one_process": 0.0,
            "clusters_in_parts": False,
            "cut": None,
            "cluster_method": cluster_method,
            "reuse_source": None,
            "reuse_reasons": {},
            "label": "estimate: no stages",
        }
    rows, missing = bs.release_proteins(root, stages, GCV)
    cut = l2.cut_of(
        parts,
        None if parts is not None else float(repr(bs.SHARD_MARGIN * prof.cpu_min_per_process)),
        None if prof.cluster_cpu_min_max is None else float(prof.cluster_cpu_min_max),
    )
    reuse_source: str | None = None
    reuse_reasons: dict[str, str] = {}
    plan_rels = l2._plan_releases(root)
    if plan_rels:
        try:
            inputs = l2.plan_inputs(root, release, GCV, stages, cut, mode=cluster_method)
            for src in plan_rels:
                why = l2._source_problem(root, release, src, inputs)
                reuse_reasons[src] = "qualifies" if why is None else why
                if why is None and reuse_source is None:
                    reuse_source = src
        except (FileNotFoundError, ValueError) as exc:
            reuse_reasons["_check_error"] = str(exc)

    if missing is not None or rows is None:
        return {
            "protein_seq_rows": None,
            "missing_protein_stage": missing,
            "cluster_cpu_min_one_process": None,
            "clusters_in_parts": parts is not None,
            "cut": cut,
            "cluster_method": cluster_method,
            "reuse_source": reuse_source,
            "reuse_reasons": reuse_reasons,
            "label": "estimate: protein_seq count unavailable (no genecall protein_seq outputs)",
        }
    t_one = round(bs.cluster_cpu_min(rows), 2)
    split = parts is not None or t_one > bs.SHARD_MARGIN * prof.cpu_min_per_process
    return {
        "protein_seq_rows": rows,
        "cluster_cpu_min_one_process": t_one,
        "margin_cpu_min": round(bs.SHARD_MARGIN * prof.cpu_min_per_process, 2),
        "clusters_in_parts": split,
        "cut": cut,
        "cluster_method": cluster_method,
        "reuse_source": reuse_source,
        "reuse_reasons": reuse_reasons,
        "label": "estimate: biobox.stage.cluster_cpu_min(protein_seq_rows)",
    }


def thread_plan(prof: bprofile.Profile, *, check_live: bool = True) -> dict:
    """Thread plan by the hour from the profile and `biobox.runstep`."""
    plan = {
        "day_window": f"{prof.day_start_hour:02d}:00-{prof.day_end_hour:02d}:00 weekdays",
        "max_threads_day": prof.max_threads_day,
        "night_window": (
            f"{prof.day_end_hour:02d}:00-{prof.day_start_hour:02d}:00 weekdays & weekends "
            f"(when other users < {runstep.OTHERS_BUSY_CPUS:g} cores)"
        ),
        "max_threads_night": prof.max_threads_night,
        "mem_scope": prof.mem_scope,
        "cpu_min_per_process": prof.cpu_min_per_process,
    }
    if check_live:
        try:
            now = datetime.now().astimezone()
            use = runstep.measure_cpu(window_s=min(runstep.CPU_WINDOW_S, 2.0))
            allow = runstep.allowance(prof, now, use)
            plan["now_local"] = now.strftime("%Y-%m-%d %H:%M %Z")
            plan["now_cap"] = allow.cap
            plan["now_allowance"] = allow.threads
            plan["now_ours_cpus"] = round(use.own_other + use.own_here, 2)
            plan["now_others_cpus"] = round(use.others, 2)
        except Exception as exc:  # noqa: BLE001
            plan["now_error"] = str(exc)
    return plan


def build_plan(root: Path, prof: bprofile.Profile, *, on_date: date | None = None,
               explicit_release: str | None = None, parts: int | None = None,
               same_stages: bool = False, check_live: bool = True,
               cluster_method: str = l2.DEFAULT_CLUSTERING_MODE) -> dict:
    """Full dry-run plan for the release train."""
    rel = resolve_release_name(root, on_date, explicit=explicit_release)
    sel = select_stages(root)
    est = estimate_and_reuse(root, rel, sel["stages"], prof, parts=parts, cluster_method=cluster_method)
    tplan = thread_plan(prof, check_live=check_live)
    return {
        "release": rel,
        "date": (date.today() if on_date is None else on_date).isoformat(),
        "same_stages": same_stages,
        "cluster_method": cluster_method,
        "action": (
            "noop (0 new complete stages since latest release; pass --same-stages to cut anyway)"
            if (sel["nothing_new"] and not same_stages)
            else ("cut" if sel["stages"] else "noop (no complete stages in root)")
        ),
        **sel,
        "clustering": est,
        "threads": tplan,
    }


def format_plan(doc: dict) -> str:
    """Human-readable summary for `--dry-run`."""
    lines = [
        f"Next release: {doc['release']} (date {doc['date']})",
        f"Latest defined release: {doc['latest_release'] or 'none'}",
        f"Stages ({len(doc['stages'])}): {','.join(doc['stages']) or '(none)'}",
        f"New stages ({len(doc['new_stages'])}): {','.join(doc['new_stages']) or '(none)'}",
        f"Action: {doc['action']}",
    ]
    if doc["excluded"]:
        lines.append(f"Excluded stages ({len(doc['excluded'])}; left out, never waited for):")
        for ex in doc["excluded"]:
            lines.append(f"  - {ex['stage']} ({ex['depth']}): {ex['reason']}")
    else:
        lines.append("Excluded stages: none")

    c = doc["clustering"]
    if c.get("reuse_source"):
        lines.append(
            f"Clustering: REUSE from {c['reuse_source']} (card B-49; no MMseqs2 run; "
            f"protein_seq rows: {c['protein_seq_rows']})"
        )
    elif c.get("cluster_cpu_min_one_process") is not None:
        mode = "in parts" if c["clusters_in_parts"] else "one process"
        lines.append(
            f"Clustering: {c['cluster_cpu_min_one_process']} CPU-min one-process (*estimate*) "
            f"from {c['protein_seq_rows']:,} protein_seq rows; mode={mode}; cut={json.dumps(c['cut'], sort_keys=True)}"
        )
    else:
        lines.append(f"Clustering: {c.get('label', 'unavailable')}")

    t = doc["threads"]
    lines.append(
        f"Thread plan by the hour: day ({t['day_window']}) <= {t['max_threads_day']}; "
        f"night/weekend ({t['night_window']}) <= {t['max_threads_night']}; "
        f"mem_scope={t['mem_scope']}, cpu_min_per_process={t['cpu_min_per_process']}"
    )
    if "now_allowance" in t:
        lines.append(
            f"Now ({t['now_local']}): cap={t['now_cap']}, allowance={t['now_allowance']} "
            f"(ours={t['now_ours_cpus']} cores, others={t['now_others_cpus']} cores)"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-46: release-train planner (dry-run, read-only)")
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--date", default=None, help="YYYY-MM-DD (default: today on host)")
    ap.add_argument("--release", default=None, help="explicit release name (refused if existing or 'latest')")
    ap.add_argument("--parts", type=int, default=None, help="override --cluster-parts")
    ap.add_argument("--cluster-method", choices=l2.CLUSTERING_MODES, default=l2.DEFAULT_CLUSTERING_MODE,
                    help="clustering workflow mode (card B-71, DEC-77)")
    ap.add_argument("--same-stages", action="store_true", help="allow cutting a release with 0 new stages (DEC-52)")
    ap.add_argument("--json", action="store_true", dest="as_json", help="emit JSON instead of text")
    ap.add_argument("--no-live-cpu", action="store_true", help="skip live /proc CPU measurement")
    a = ap.parse_args(argv)

    on_date = date.fromisoformat(a.date) if a.date else None
    prof = bprofile.load(a.profile)
    try:
        doc = build_plan(
            a.root,
            prof,
            on_date=on_date,
            explicit_release=a.release,
            parts=a.parts,
            same_stages=a.same_stages,
            check_live=not a.no_live_cpu,
            cluster_method=a.cluster_method,
        )
    except ReleaseTrainError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    if a.as_json:
        print(json.dumps(doc, indent=1, sort_keys=True))
    else:
        print(format_plan(doc))
    return 0


if __name__ == "__main__":
    sys.exit(main())

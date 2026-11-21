"""Card C-46: release-train estimates before launch and rate guard (DEC-43, DEC-53, DEC-64, AGENTS §6).

    python train_est.py before --root ROOT --release R --stages S0,S1 --profile P
                               [--earlier PREV] [--parts K] [--jobs J] [--cpus C] [--tmp D]
    python3 -I train_est.py rates --root ROOT --release R --est BEFORE.json [--out OUT.json]

`before` runs in the runtime image (read-only on the root):
- Counts `protein_seq` rows and `feature_protein` genes per stage from each stage's `genecall.json` outputs.
- Computes `stage.cluster_cpu_min(rows)` (*estimate*), the split/cut policy (`l2.cut_of`), and checks `B-49` reuse.
- When `--earlier` names an earlier release with run records (and optional `l2/cluster_plan/release=<PREV>/*.cut.json`),
  calibrates call-by-call factors from that earlier release and scales them to this release (*derived*).
  When no earlier calibrated release exists, uses cut estimates at factor 1.0 (*estimate*).

`rates` runs on the host with the standard library only:
- Judges every published step of `release` estimated at `>= 24` CPU-minutes (DEC-43).
- `ratio > 1.20` (`RECORD_ABOVE`): outside 20%, listed in `outside_20pct` (DEC-53).
- Card A-16 (DEC-64): a judged step whose reported three-decimal ratio is over 1.50 (`STOP_ABOVE`) is listed in
  `over_20pct` and exits 3 (stop unit at once), as a step over 1.20 did before; one over 1.20 and at or under 1.50 is
  listed under `recorded_over_20pct` (its row: step, ratio, estimate, basis), stays in `outside_20pct`, and exits 0.
- `ratio < 0.80` (`under_20pct`): recorded and reported, exit 0 (does not stop; DEC-53).
- Card A-13: a judged step over 1.20 that an entry of `ACCEPTED` matches (its release exactly, its step by
  `fnmatchcase`, its reported three-decimal ratio at or under the entry's bound) is listed under `accepted_over_20pct`
  with the entry's decision row and stays in `outside_20pct`; `over_20pct` and `worst_over` hold only what is not
  accepted. Since DEC-64 an entry decides only above 1.50: a step it matches at or under 1.50 is in
  `recorded_over_20pct` too, and goes on with or without the entry. Where two entries match a step, the tighter bound
  decides. A step over 1.50 and over the bound of the entry that matches it stops, its row naming the entry under
  `accepted_bound_exceeded`. The key `accepted_over_20pct` appears only for a release an entry names.
  A malformed table is refused with exit 3, the one status on which `release_train.sh` stopped the unit before card
  A-14 (it now stops on any status but 0).
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from fnmatch import fnmatchcase
from pathlib import Path
from typing import NamedTuple

TOL = 0.20
MIN_JUDGED_CPU_MIN = 24.0  # DEC-43
RECORD_ABOVE = 1.20  # DEC-64 (amending DEC-53): a judged step over this ratio is outside 20% and recorded
STOP_ABOVE = 1.50  # DEC-64: a judged step whose reported (three-decimal) ratio is over this stops the train
GENE_SCALED_STEPS = ("tokens", "postings", "dim", "manifest", "catalog")


class Accept(NamedTuple):
    """One record class the PI accepted over 20%: a release, a step pattern, the greatest ratio, the decision row."""

    release: str
    step: str
    max_ratio: float
    decision: str


# Card A-13: records accepted over 20% by a decision row of planning/DECISIONS.md. An entry is added only in the
# commit that adds its row, and only here: no flag, no variable of the shell and no file adds one. Empty: every
# judged step more than 1.50 times its estimate stops the train (DEC-43, DEC-64). Since DEC-64 a new entry's bound
# must be over STOP_ABOVE: at or under it, the rule already lets the step go on.
ACCEPTED: tuple[Accept, ...] = (
    # DEC-61 (PI, 2026-11-21: "accept up to 1.5"): the parts of the last NR50 alignment of r26.11.2; part-00001 ran at
    # 1.215 of an estimate whose factor card B-54 fitted on r26.11.1 alone.
    Accept("r26.11.2", "cluster.nr50.20.align.part-*", 1.50, "DEC-61"),
)

# The greatest bound an entry may carry: DEC-55's 1.60, the largest ratio the PI has accepted (NR50 parts of about
# 56 CPU-minutes against 45, far under the 280 CPU-minute process cap). A step further off than that needs its
# estimate refitted, not an entry in this table.
ACCEPT_CEILING = 1.60

# Card A-16: the entries of ACCEPTED made before DEC-64, each exactly as its decision row gave it, kept as the record
# of that decision. Only these may carry a bound over RECORD_ABOVE and at or under STOP_ABOVE; every other entry's
# bound must be over STOP_ABOVE. DEC-61's bound is STOP_ABOVE itself, so it accepts no step that the rule does not
# already let go on (a step it matches at or under 1.50 is in the band; one over 1.50 is over its bound and stops).
KEPT_BEFORE_DEC64: tuple[Accept, ...] = (
    Accept("r26.11.2", "cluster.nr50.20.align.part-*", 1.50, "DEC-61"),
)
_DECISION = re.compile(r"DEC-[1-9][0-9]*")


class AcceptedTableError(ValueError):
    """`ACCEPTED` holds an entry the guard refuses to read."""


def accepted_table_problems(table) -> list[str]:
    """Every reason the table is refused, or an empty list."""
    problems: list[str] = []
    seen: set[tuple[str, str]] = set()
    for k, e in enumerate(table):
        if not isinstance(e, tuple) or len(e) != 4:
            problems.append(f"entry {k}: {e!r} is not (release, step pattern, greatest ratio, decision row)")
            continue
        rel, pat, bound, dec = e
        if not isinstance(rel, str) or not rel or any(c in rel for c in "*?[]"):
            problems.append(f"entry {k}: release {rel!r} must be one release's exact name")
        if not isinstance(pat, str) or not pat.strip("*?"):
            problems.append(f"entry {k}: step pattern {pat!r} must name steps, not be empty or wildcards alone")
        if isinstance(bound, bool) or not isinstance(bound, (int, float)) or not math.isfinite(bound):
            problems.append(f"entry {k}: greatest ratio {bound!r} is not a finite number")
        elif not (RECORD_ABOVE < bound <= ACCEPT_CEILING):
            problems.append(f"entry {k}: greatest ratio {bound!r} must be over {STOP_ABOVE:.2f} (DEC-64: at or under "
                            f"it the rule already lets a step go on) and at most {ACCEPT_CEILING:.2f} (DEC-55)")
        elif bound <= STOP_ABOVE and not any(e == kept for kept in KEPT_BEFORE_DEC64):
            problems.append(f"entry {k}: greatest ratio {bound!r} must be over {STOP_ABOVE:.2f} (DEC-64: at or under "
                            f"it the rule already lets a step go on; only an entry made before DEC-64, as its row "
                            f"gave it, keeps a lower bound)")
        if not isinstance(dec, str) or not _DECISION.fullmatch(dec):
            problems.append(f"entry {k}: decision {dec!r} must be a row's name, DEC-<number>")
        if isinstance(rel, str) and isinstance(pat, str):
            if (rel, pat) in seen:
                problems.append(f"entry {k}: release {rel} and step pattern {pat} are named twice")
            seen.add((rel, pat))
    return problems


def _check_accepted() -> None:
    problems = accepted_table_problems(ACCEPTED)
    if problems:
        raise AcceptedTableError("; ".join(problems))


def _matching(release: str, step: str) -> list[Accept]:
    """The entries that match `step` of `release`, the tightest bound first (table order among equal bounds)."""
    return sorted((Accept(*e) for e in ACCEPTED if e[0] == release and fnmatchcase(step, e[1])),
                  key=lambda e: e.max_ratio)


def _load_json(p: Path) -> dict | None:
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return None


def _gene_call_outputs(root: Path, stage: str, table: str, gcv: str) -> list[Path]:
    rec = _load_json(root / "runs" / f"stage={stage}" / "genecall.json")
    if not rec:
        return []
    prefix = f"{table}/gene_call_version={gcv}/stage={stage}/"
    return [root / o["path"] for o in rec.get("outputs", ())
            if o["path"].startswith(prefix) and o["path"].endswith(".parquet")]


def calibrate_earlier(root: Path, release: str) -> dict:
    """Summarize an earlier release's clustering calls and other release steps from its run records."""
    runs = root / "runs" / f"stage={release}"
    plans = root / "l2" / "cluster_plan" / f"release={release}"
    seen: set[str] = set()
    calls = []
    part_rows: dict[str, list[dict]] = {}
    if plans.is_dir():
        for cutf in sorted(plans.glob("*.cut.json")):
            call = cutf.name[:-len(".cut.json")]
            cut = _load_json(cutf) or {}
            parts = cut.get("parts") or []
            rows = []
            for p, part in enumerate(parts):
                name = f"{call}.part-{p:05d}" if len(parts) > 1 else call
                r = _load_json(runs / f"{name}.json")
                seen.add(name)
                row = {"part": p, "estimate": round(part["cpu_min"], 3), "units": part.get("units")}
                if r is not None:
                    row["cpu_min"] = round(r["cpu_min"], 3)
                    row["peak_rss_gb"] = round((r.get("peak_rss_bytes") or 0) / 1e9, 2)
                rows.append(row)
            meas = [x["cpu_min"] for x in rows if "cpu_min" in x]
            est = sum(x["estimate"] for x in rows)
            part_rows[call] = rows
            calls.append({
                "call": call,
                "parts": len(parts),
                "overhead_estimate": cut.get("overhead_cpu_min"),
                "estimate_sum": round(est, 2),
                "measured_sum": round(sum(meas), 2),
                "max_part_cpu_min": max(meas) if meas else None,
                "max_peak_rss_gb": max((x["peak_rss_gb"] for x in rows if "peak_rss_gb" in x), default=None),
            })
    others = []
    if runs.is_dir():
        for f in sorted(runs.glob("*.json")):
            name = f.name[:-len(".json")]
            if name.endswith(".publishing") or name in seen:
                continue
            r = _load_json(f)
            if r is not None and "cpu_min" in r:
                others.append({"step": name, "cpu_min": round(r["cpu_min"], 3)})
    return {"release": release, "calls": calls, "others": others, "part_rows": part_rows}


def before(a) -> int:
    import duckdb
    import pyarrow.parquet as pq

    from biobox import catalog, profile as bprofile, stage as bs
    from biobox.l1.genecall import GENE_CALL_VERSION as GCV
    from biobox.l2 import cluster as l2

    root = Path(a.root)
    stages = [s for s in a.stages.split(",") if s]
    prof = bprofile.load(a.profile)

    per_stage: dict[str, dict] = {}
    files_of: dict[str, list[Path]] = {}
    for s in stages:
        ps = _gene_call_outputs(root, s, "l0seq/protein_seq", GCV)
        fp = _gene_call_outputs(root, s, "l1/feature_protein", GCV)
        files_of[s] = ps
        per_stage[s] = {
            "protein_seq_rows": sum(pq.read_metadata(p).num_rows for p in ps if p.is_file()),
            "genes": sum(pq.read_metadata(p).num_rows for p in fp if p.is_file()),
        }

    rows, missing = bs.release_proteins(root, stages, GCV)
    if missing is not None or rows is None:
        raise SystemExit(f"release_proteins: missing gene-call record for stage {missing}")

    t_one = bs.cluster_cpu_min(rows)
    split = a.parts is not None or t_one > bs.SHARD_MARGIN * prof.cpu_min_per_process
    cut = l2.cut_of(
        a.parts,
        None if a.parts is not None else float(repr(bs.SHARD_MARGIN * prof.cpu_min_per_process)),
        None if prof.cluster_cpu_min_max is None else float(prof.cluster_cpu_min_max),
    )

    reuse_source = None
    reuse_reasons: dict[str, str] = {}
    plan_rels = l2._plan_releases(root)
    if plan_rels:
        inputs = l2.plan_inputs(root, a.release, GCV, stages, cut)
        for src in plan_rels:
            why = l2._source_problem(root, a.release, src, inputs)
            reuse_reasons[src] = "qualifies" if why is None else why
            if why is None and reuse_source is None:
                reuse_source = src

    genes_all = sum(per_stage[s]["genes"] for s in stages)
    earlier = a.earlier
    estages: list[str] = []
    if earlier:
        estages = catalog.read_release_def(root, earlier) or []

    factors: dict[str, float | None] = {}
    steps: dict[str, float] = {}
    calls_out: list[dict] = []
    part_max = 0.0
    total_proj = 0.0
    n_all: int | None = None
    n_old: int | None = None

    if earlier and estages and all(s in per_stage for s in estages):
        calib = calibrate_earlier(root, earlier)
        if any(files_of[s] for s in stages):
            Path(a.tmp).mkdir(parents=True, exist_ok=True)
            con = duckdb.connect(config={
                "threads": a.duck_threads, "memory_limit": "4GB",
                "temp_directory": str(a.tmp), "preserve_insertion_order": False,
            })

            def distinct(ss: list[str]) -> int:
                lit = "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for s in ss for p in files_of[s]) + "]"
                return con.execute(
                    f"SELECT count(DISTINCT protein_uid) FROM read_parquet({lit}, hive_partitioning = false)"
                ).fetchone()[0]

            n_all = distinct(stages)
            n_old = distinct(estages)
            con.close()

        r = (n_all / n_old) if (n_all and n_old) else (rows / max(1, sum(per_stage[s]["protein_seq_rows"] for s in estages)))
        genes_old = max(1, sum(per_stage[s]["genes"] for s in estages))
        rg = genes_all / genes_old

        for c in calib["calls"]:
            name = c["call"]
            rows_ = calib["part_rows"][name]
            u_old = sum(x["units"] or 0 for x in rows_)
            o_old = c["overhead_estimate"] or 0.0
            k = len(rows_) if a.parts is None else a.parts
            if n_all and n_old:
                m_old = l2.call_cpu_min(name, n_old, u_old) + len(rows_) * o_old
                m_new = l2.call_cpu_min(name, n_all, round(u_old * r)) + k * o_old * r
                g = m_new / m_old if m_old else r
            else:
                g = r
            proj = c["measured_sum"] * g
            largest = (c["max_part_cpu_min"] or 0) * g * len(rows_) / max(1, k)
            factors[name] = round(c["measured_sum"] / c["estimate_sum"], 6) if c["estimate_sum"] else 1.0
            part_max = max(part_max, largest)
            total_proj += proj
            calls_out.append({
                "call": name, "parts": k, "measured_earlier": c["measured_sum"],
                "growth": round(g, 4), "projected": round(proj, 1),
                "largest_part_projected": round(largest, 1),
                "calibration_factor": factors[name],
                "peak_rss_gb_projected": round((c["max_peak_rss_gb"] or 0) * r, 1),
            })
        for o in calib["others"]:
            scale = rg if o["step"] in GENE_SCALED_STEPS else r
            steps[o["step"]] = round(o["cpu_min"] * scale, 3)

    others_total = sum(v for k_, v in steps.items() if k_.startswith("cluster"))
    doc = {
        "release": a.release,
        "stages": stages,
        "earlier": earlier,
        "earlier_stages": estages,
        "measured": {
            "per_stage": per_stage,
            "protein_seq_rows": rows,
            "distinct_proteins": n_all,
            "distinct_proteins_earlier": n_old,
            "genes": genes_all,
            "label": "measured: Parquet footers and DuckDB count(DISTINCT) over genecall outputs",
        },
        "code": {
            "cluster_cpu_min_one_process": round(t_one, 2),
            "margin_cpu_min": round(bs.SHARD_MARGIN * prof.cpu_min_per_process, 2),
            "clusters_in_parts": split,
            "cut": cut,
            "reuse_source": reuse_source,
            "reuse_reasons": reuse_reasons,
            "label": "estimate: biobox.stage and biobox.l2.cluster",
        },
        "projection": {
            "cluster_cpu_min": 0.0 if reuse_source else round(total_proj + others_total if calls_out else t_one, 1),
            "largest_part_cpu_min": 0.0 if reuse_source else round(part_max, 1),
            "every_part_under_240": True if reuse_source else (part_max < 240),
            "calls": calls_out,
            "steps": steps,
            "factors": factors,
            "label": (
                f"derived: B-49 clustering reuse from {reuse_source}" if reuse_source
                else ("derived: scaled from earlier release" if calls_out else "estimate: uncalibrated cut estimates")
            ),
        },
    }
    print(json.dumps(doc, indent=1, sort_keys=True))
    return 0


def check_rates(root: Path, release: str, est_path: Path, *, out_path: Path | None = None) -> dict:
    """Judge published run records of `release` under DEC-43 (>= 24 CPU-min), DEC-53 (outside 20%: over 1.20 or
    under 0.80) and DEC-64 (over 1.50 stops; over 1.20 and at or under 1.50 is recorded), less the records `ACCEPTED`
    names (card A-13)."""
    _check_accepted()
    named = any(e[0] == release for e in ACCEPTED)
    est = json.loads(Path(est_path).read_text())
    proj = est.get("projection", {})
    factors = proj.get("factors", {})
    steps = proj.get("steps", {})
    runs = Path(root) / "runs" / f"stage={release}"
    plans = Path(root) / "l2" / "cluster_plan" / f"release={release}"

    rows, over_bad, under_bad, bad, unknown, accepted, recorded = [], [], [], [], [], [], []
    cuts: dict[str, dict] = {}
    for f in sorted(runs.glob("*.json")) if runs.is_dir() else []:
        if f.name.endswith(".publishing.json"):
            continue
        name = f.name[:-len(".json")]
        r = _load_json(f)
        if r is None or "cpu_min" not in r:
            continue
        e, basis = None, None
        call, _, part = name.partition(".part-")
        if part:
            if call not in cuts:
                cp = plans / f"{call}.cut.json"
                cuts[call] = _load_json(cp) or {}
            parts = cuts[call].get("parts", [])
            if parts and int(part) < len(parts):
                fac = factors.get(call) or 1.0
                e = parts[int(part)]["cpu_min"] * fac
                basis = f"cut estimate {parts[int(part)]['cpu_min']:.3f} x factor {fac}"
        elif (plans / f"{name}.cut.json").is_file() and name not in cuts:
            cut = _load_json(plans / f"{name}.cut.json") or {}
            if len(cut.get("parts", [])) == 1:
                fac = factors.get(name) or 1.0
                e = cut["parts"][0]["cpu_min"] * fac
                basis = f"one part: cut estimate x factor {fac}"
        if e is None and name in steps:
            e, basis = steps[name], "earlier record x ratio"

        row = {
            "step": name,
            "cpu_min": round(r["cpu_min"], 3),
            "wall_s": r.get("wall_s"),
            "peak_rss_gb": round((r.get("peak_rss_bytes") or 0) / 1e9, 2),
        }
        if e is None:
            unknown.append(name)
        else:
            ratio = r["cpu_min"] / e if e > 0 else float("inf")
            judged = e >= MIN_JUDGED_CPU_MIN
            row.update({"estimate": round(e, 3), "ratio": round(ratio, 3), "judged": judged, "basis": basis})
            if judged and ratio > RECORD_ABOVE:
                match = _matching(release, name) if named else []
                if match and row["ratio"] <= match[0].max_ratio:  # the reported ratio, as the decision row reads it
                    row["accepted"] = match[0].decision
                    accepted.append(row)
                if row["ratio"] <= STOP_ABOVE:  # DEC-64, on the reported ratio: recorded, the train goes on
                    recorded.append(row)
                elif "accepted" not in row:
                    if match:
                        row["accepted_bound_exceeded"] = {"decision": match[0].decision, "step": match[0].step,
                                                          "max_ratio": match[0].max_ratio}
                    over_bad.append(row)
                bad.append(row)
            elif judged and ratio < 1 - TOL:
                under_bad.append(row)
                bad.append(row)
        rows.append(row)

    judged = [x for x in rows if x.get("judged")]
    out = {
        "release": release,
        "records": len(rows),
        "judged": len(judged),
        "cpu_min": round(sum(x["cpu_min"] for x in rows), 2),
        "estimate_of_recorded": round(sum(x.get("estimate", 0) for x in rows), 2),
        "worst_over": max((x for x in judged if "accepted" not in x), key=lambda x: x["ratio"], default=None),
        "largest": max(rows, key=lambda x: x["cpu_min"], default=None),
        "over_20pct": over_bad,
        "recorded_over_20pct": recorded,
        "under_20pct": under_bad,
        "outside_20pct": bad,
        "no_estimate": unknown[-8:],
        "no_estimate_n": len(unknown),
        "exit_code": 3 if over_bad else 0,
    }
    if named:
        out["accepted_over_20pct"] = accepted
    if out_path is not None:
        Path(out_path).write_text(json.dumps({**out, "rows": rows}, indent=1, sort_keys=True))
    return out


def rates(a) -> int:
    try:
        _check_accepted()
    except AcceptedTableError as exc:
        print(f"train_est.py rates: the ACCEPTED table is refused, the train stops: {exc}", file=sys.stderr)
        print(json.dumps({"release": a.release, "accepted_table_error": str(exc)}, sort_keys=True))
        return 3
    out = check_rates(Path(a.root), a.release, Path(a.est), out_path=Path(a.out) if a.out else None)
    print(json.dumps({k: v for k, v in out.items() if k != "exit_code"}, sort_keys=True))
    return out["exit_code"]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-46: release-train estimates and rate guard")
    ap.add_argument("mode", choices=("before", "rates"))
    ap.add_argument("--root", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--stages", default="")
    ap.add_argument("--earlier", default=None)
    ap.add_argument("--parts", type=int, default=None)
    ap.add_argument("--profile")
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--cpus", type=int, default=None)
    ap.add_argument("--tmp", default="/tmp/if2453/train/tmp/est")
    ap.add_argument("--duck-threads", type=int, default=2)
    ap.add_argument("--est")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    return {"before": before, "rates": rates}[a.mode](a)


if __name__ == "__main__":
    sys.exit(main())

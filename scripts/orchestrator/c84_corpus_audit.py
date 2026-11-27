#!/usr/bin/env python3
"""Card C-84 validator for the three-tier corpus inventory and acquisition roadmap.

Verifies exhaustive disposition mapping across pod stages/inputs/staging, `tower`
`sources_incoming/` and `sources/` directories, stage `.conf` and template files,
`scripts/orchestrator/tower_jobs/*.tsv` job lists, and `sources.yaml` entries, plus
`AGENTS.md` §6 number labelling (`measured`, `derived`, `estimate` with non-empty `source`).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
DEFAULT_INVENTORY = REPO / "planning" / "measurements" / "2026-11-27_corpus_inventory_and_roadmap.json"
DEFAULT_SOURCES = REPO / "sources.yaml"
DEFAULT_JOBS_DIR = REPO / "scripts" / "orchestrator" / "tower_jobs"
DEFAULT_STAGES_DIR = REPO / "scripts" / "stress" / "build_stage_pod" / "stages"

DISPOSITIONS: tuple[str, ...] = (
    "built",
    "building_or_queued_r26_11_4",
    "staged_r26_11_5",
    "staged_later_tranche",
    "landed_on_tower_unstaged",
    "actively_fetching_on_tower",
    "failed_or_corrupt_archive",
    "wrong_payload_non_dna",
    "non_sequence_aux",
    "superseded_duplicate",
)
DISPOSITION_SET = frozenset(DISPOSITIONS)

LABELS: tuple[str, ...] = ("measured", "derived", "estimate")
LABEL_SET = frozenset(LABELS)

_SOURCE_NAME_RE = re.compile(r"^  - name:\s*(\S+)\s*(?:#.*)?$")


def _is_number(val: Any) -> bool:
    return not isinstance(val, bool) and isinstance(val, (int, float))


def _is_number_wrapper(obj: dict[str, Any]) -> bool:
    return "value" in obj and ("label" in obj or "source" in obj)


def check_numbers(obj: Any, path: str = "root") -> list[str]:
    """Return `AGENTS.md` §6 labelling problems across all numeric leaves in `obj`."""
    problems: list[str] = []
    if _is_number(obj):
        problems.append(f"{path}: bare unlabelled number {obj!r}")
        return problems
    if isinstance(obj, list):
        for idx, item in enumerate(obj):
            problems.extend(check_numbers(item, f"{path}[{idx}]"))
        return problems
    if isinstance(obj, dict):
        if _is_number_wrapper(obj):
            val = obj.get("value")
            if not _is_number(val):
                problems.append(f"{path}.value: expected int or float, got {val!r}")
            label = obj.get("label")
            if not isinstance(label, str) or label not in LABEL_SET:
                problems.append(
                    f"{path}.label: invalid label {label!r} (must be one of {LABELS})"
                )
            source = obj.get("source")
            if not isinstance(source, str) or not source.strip():
                problems.append(f"{path}.source: missing or empty source string")
            low = obj.get("low")
            if low is not None and not _is_number(low):
                problems.append(f"{path}.low: expected number, got {low!r}")
            high = obj.get("high")
            if high is not None and not _is_number(high):
                problems.append(f"{path}.high: expected number, got {high!r}")
            if _is_number(val) and _is_number(low) and low > val:
                problems.append(f"{path}: low ({low}) > value ({val})")
            if _is_number(val) and _is_number(high) and val > high:
                problems.append(f"{path}: value ({val}) > high ({high})")
            for k, v in obj.items():
                if k not in ("value", "low", "high", "label", "source"):
                    problems.extend(check_numbers(v, f"{path}.{k}"))
            return problems
        for k, v in obj.items():
            problems.extend(check_numbers(v, f"{path}.{k}"))
    return problems


def parse_sources_yaml_names(sources_path: Path) -> list[str]:
    """Extract ordered `sources[*].name` entries from `sources.yaml`."""
    text = sources_path.read_text(encoding="utf-8")
    try:
        import yaml  # type: ignore[import-not-found]

        parsed = yaml.safe_load(text)
        if isinstance(parsed, dict) and isinstance(parsed.get("sources"), list):
            return [str(entry["name"]) for entry in parsed["sources"] if isinstance(entry, dict) and "name" in entry]
    except ModuleNotFoundError:
        pass

    names: list[str] = []
    in_sources = False
    for raw_line in text.splitlines():
        if raw_line.startswith("sources:"):
            in_sources = True
            continue
        if in_sources and raw_line and not raw_line.startswith((" ", "#")):
            in_sources = False
        if in_sources:
            m = _SOURCE_NAME_RE.match(raw_line)
            if m:
                names.append(m.group(1))
    return names


def count_job_lines(tsv_path: Path) -> int:
    """Count non-empty, non-comment lines in a `tower_jobs/*.tsv` file."""
    count = 0
    for raw in tsv_path.read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        count += 1
    return count


def _check_disposition(disp: Any, context: str, problems: list[str]) -> None:
    if not isinstance(disp, str) or disp not in DISPOSITION_SET:
        problems.append(
            f"{context}: invalid disposition {disp!r} (allowed: {', '.join(DISPOSITIONS)})"
        )


def validate_inventory(
    data: dict[str, Any],
    *,
    stages_dir: Path,
    jobs_dir: Path,
    sources_path: Path | None,
) -> tuple[list[str], list[str]]:
    """Validate classification completeness and numeric labels; returns `(unclassified, number_problems)`."""
    unclassified: list[str] = []
    number_problems = check_numbers(data)

    observed = data.get("observed")
    if not isinstance(observed, dict):
        unclassified.append("observed: missing or non-object 'observed' block")
        return unclassified, number_problems

    # 1. tower sources_incoming directories (including nested batch*/* directories)
    obs_incoming = observed.get("tower_sources_incoming_dirs")
    incoming_rows = data.get("tower_sources_incoming")
    if not isinstance(obs_incoming, list) or not isinstance(incoming_rows, list):
        unclassified.append("tower_sources_incoming: missing observed list or classification rows")
    else:
        seen_inc: set[str] = set()
        for idx, row in enumerate(incoming_rows):
            if not isinstance(row, dict) or not isinstance(row.get("dir"), str):
                unclassified.append(f"tower_sources_incoming[{idx}]: missing string 'dir'")
                continue
            d = row["dir"]
            if d in seen_inc:
                unclassified.append(f"tower_sources_incoming: duplicate dir {d!r}")
            seen_inc.add(d)
            _check_disposition(row.get("disposition"), f"tower_sources_incoming[{d}]", unclassified)
            if not isinstance(row.get("summary"), str) or not row["summary"].strip():
                unclassified.append(f"tower_sources_incoming[{d}]: missing summary")
        obs_inc_set = {str(x) for x in obs_incoming}
        for missing in sorted(obs_inc_set - seen_inc):
            unclassified.append(f"tower_sources_incoming: unclassified observed dir {missing!r}")
        for extra in sorted(seen_inc - obs_inc_set):
            unclassified.append(f"tower_sources_incoming: entry {extra!r} not in observed.tower_sources_incoming_dirs")

    # 2. tower sources/ mirror directories
    obs_mirrors = observed.get("tower_sources_dirs")
    mirror_rows = data.get("tower_sources_mirrors")
    if not isinstance(obs_mirrors, list) or not isinstance(mirror_rows, list):
        unclassified.append("tower_sources_mirrors: missing observed list or classification rows")
    else:
        seen_mir: set[str] = set()
        for idx, row in enumerate(mirror_rows):
            if not isinstance(row, dict) or not isinstance(row.get("dir"), str):
                unclassified.append(f"tower_sources_mirrors[{idx}]: missing string 'dir'")
                continue
            d = row["dir"]
            if d in seen_mir:
                unclassified.append(f"tower_sources_mirrors: duplicate dir {d!r}")
            seen_mir.add(d)
            _check_disposition(row.get("disposition"), f"tower_sources_mirrors[{d}]", unclassified)
        obs_mir_set = {str(x) for x in obs_mirrors}
        for missing in sorted(obs_mir_set - seen_mir):
            unclassified.append(f"tower_sources_mirrors: unclassified observed dir {missing!r}")
        for extra in sorted(seen_mir - obs_mir_set):
            unclassified.append(f"tower_sources_mirrors: entry {extra!r} not in observed.tower_sources_dirs")

    # 3. pod stages, s2/inputs, bio-box-staging, stages/*.conf, stages/templates/*.conf
    pod_rows = data.get("pod_stages_and_assets")
    if not isinstance(pod_rows, list):
        unclassified.append("pod_stages_and_assets: missing list")
    else:
        seen_pod: set[str] = set()
        for idx, row in enumerate(pod_rows):
            if not isinstance(row, dict) or not isinstance(row.get("name"), str):
                unclassified.append(f"pod_stages_and_assets[{idx}]: missing string 'name'")
                continue
            nm = row["name"]
            if nm in seen_pod:
                unclassified.append(f"pod_stages_and_assets: duplicate name {nm!r}")
            seen_pod.add(nm)
            _check_disposition(row.get("disposition"), f"pod_stages_and_assets[{nm}]", unclassified)

        for key in ("pod_root_stages", "pod_s2_inputs_dirs", "pod_biobox_staging_dirs"):
            obs_list = observed.get(key)
            if not isinstance(obs_list, list):
                unclassified.append(f"observed.{key}: missing list")
                continue
            for item in sorted({str(x) for x in obs_list} - seen_pod):
                unclassified.append(f"observed.{key}: unclassified pod stage/dir {item!r}")

        if not stages_dir.is_dir():
            unclassified.append(f"stages_dir not found: {stages_dir}")
        else:
            conf_stems = {p.stem for p in stages_dir.glob("*.conf")}
            for missing in sorted(conf_stems - seen_pod):
                unclassified.append(f"stages/*.conf: unclassified stage conf {missing!r}")
            tpl_dir = stages_dir / "templates"
            if tpl_dir.is_dir():
                tpl_stems = {p.stem for p in tpl_dir.glob("*.conf")}
                for missing in sorted(tpl_stems - seen_pod):
                    unclassified.append(f"stages/templates/*.conf: unclassified template conf {missing!r}")

    # 4. tower_jobs/*.tsv job lists and exact job counts
    job_rows = data.get("tower_job_lists")
    if not isinstance(job_rows, list):
        unclassified.append("tower_job_lists: missing list")
    elif not jobs_dir.is_dir():
        unclassified.append(f"jobs_dir not found: {jobs_dir}")
    else:
        disk_jobs = {p.name: count_job_lines(p) for p in jobs_dir.glob("*.tsv")}
        seen_jobs: dict[str, dict[str, Any]] = {}
        for idx, row in enumerate(job_rows):
            if not isinstance(row, dict) or not isinstance(row.get("file"), str):
                unclassified.append(f"tower_job_lists[{idx}]: missing string 'file'")
                continue
            fn = row["file"]
            if fn in seen_jobs:
                unclassified.append(f"tower_job_lists: duplicate file {fn!r}")
            seen_jobs[fn] = row
            _check_disposition(row.get("disposition"), f"tower_job_lists[{fn}]", unclassified)

        for missing in sorted(set(disk_jobs) - set(seen_jobs)):
            unclassified.append(f"tower_job_lists: unclassified job list {missing!r}")
        for extra in sorted(set(seen_jobs) - set(disk_jobs)):
            unclassified.append(f"tower_job_lists: unknown job list {extra!r} not on disk")
        for fn, expected_cnt in sorted(disk_jobs.items()):
            if fn in seen_jobs:
                jobs_obj = seen_jobs[fn].get("jobs")
                rec_cnt = jobs_obj.get("value") if isinstance(jobs_obj, dict) else None
                if rec_cnt != expected_cnt:
                    unclassified.append(
                        f"tower_job_lists[{fn}]: recorded jobs={rec_cnt!r} != disk jobs={expected_cnt}"
                    )

    # 5. sources.yaml reconciliation
    src_rows = data.get("sources_yaml_reconciliation")
    if not isinstance(src_rows, list):
        unclassified.append("sources_yaml_reconciliation: missing list")
    else:
        seen_src: set[str] = set()
        for idx, row in enumerate(src_rows):
            if not isinstance(row, dict) or not isinstance(row.get("name"), str):
                unclassified.append(f"sources_yaml_reconciliation[{idx}]: missing string 'name'")
                continue
            nm = row["name"]
            if nm in seen_src:
                unclassified.append(f"sources_yaml_reconciliation: duplicate source {nm!r}")
            seen_src.add(nm)
            _check_disposition(row.get("disposition"), f"sources_yaml_reconciliation[{nm}]", unclassified)
            if not isinstance(row.get("reason_code"), str) or not row["reason_code"].strip():
                unclassified.append(f"sources_yaml_reconciliation[{nm}]: missing reason_code")
            if not isinstance(row.get("detail"), str) or not row["detail"].strip():
                unclassified.append(f"sources_yaml_reconciliation[{nm}]: missing detail")

        if sources_path is not None:
            if not sources_path.is_file():
                unclassified.append(f"sources_path not found: {sources_path}")
            else:
                yaml_names = set(parse_sources_yaml_names(sources_path))
                for missing in sorted(yaml_names - seen_src):
                    unclassified.append(f"sources.yaml: unclassified source entry {missing!r}")
                for extra in sorted(seen_src - yaml_names):
                    unclassified.append(f"sources.yaml: reconciliation entry {extra!r} not in sources.yaml")

    # 6. Summary tiers & roadmap priorities presence
    if not isinstance(data.get("summary_tiers"), list) or not data["summary_tiers"]:
        unclassified.append("summary_tiers: must be a non-empty list")
    if not isinstance(data.get("roadmap_priorities"), list) or not data["roadmap_priorities"]:
        unclassified.append("roadmap_priorities: must be a non-empty list")

    return unclassified, number_problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Run exhaustive inventory and label checks")
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    parser.add_argument("--no-require-sources", action="store_true", help="Skip on-disk sources.yaml comparison")
    parser.add_argument("--jobs-dir", type=Path, default=DEFAULT_JOBS_DIR)
    parser.add_argument("--stages-dir", type=Path, default=DEFAULT_STAGES_DIR)
    args = parser.parse_args(argv)

    if not args.inventory.is_file():
        print(f"c84_corpus_audit: inventory file not found: {args.inventory}", file=sys.stderr)
        return 1

    try:
        data = json.loads(args.inventory.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        print(f"c84_corpus_audit: invalid JSON in {args.inventory}: {exc}", file=sys.stderr)
        return 1

    if not isinstance(data, dict):
        print("c84_corpus_audit: top-level JSON value must be an object", file=sys.stderr)
        return 1

    sources_path = None if args.no_require_sources else args.sources
    unclassified, number_problems = validate_inventory(
        data,
        stages_dir=args.stages_dir,
        jobs_dir=args.jobs_dir,
        sources_path=sources_path,
    )

    if unclassified or number_problems:
        for msg in unclassified:
            print(f"UNCLASSIFIED: {msg}", file=sys.stderr)
        for msg in number_problems:
            print(f"NUMBER_LABEL: {msg}", file=sys.stderr)
        print(
            f"c84_corpus_audit: FAILED unclassified={len(unclassified)} number_problems={len(number_problems)}",
            file=sys.stderr,
        )
        return 1

    print(
        f"c84_corpus_audit: OK unclassified=0 number_problems=0 "
        f"incoming={len(data.get('tower_sources_incoming', []))} "
        f"mirrors={len(data.get('tower_sources_mirrors', []))} "
        f"pod_assets={len(data.get('pod_stages_and_assets', []))} "
        f"job_lists={len(data.get('tower_job_lists', []))} "
        f"sources={len(data.get('sources_yaml_reconciliation', []))}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Card C-42 step (2): SPIRE sample metadata merge, per-habitat summaries, and sqrt-stratified ~1,000-sample selection.

    python -I c42_select.py run --tables-dir /tmp/if2453/s3/spire/tables --out-dir /tmp/if2453/s3/spire/select
        [--target 1000] [--seed 20261118]

Reads only SPIRE's metadata tables (no sequence):
1. `study_tsv/*.tsv` (`spire_study_api`): every SPIRE study's `sample_id, mags, lat, lon, microntology`
2. `api_samples/page_*.tsv` (`spire_pilot_assemblies` `url_samples_api`): per-sample `assembly_bp`,
   `assembly_n50`, `binned_fraction`, `bases`, `reads`, `n_mags`, `n_hq_mags`, `n_mq_mags`
3. `spire_v1_genome_metadata.tsv.gz` (`spire_tables`): all SPIRE v1 MQ/HQ MAGs, summed per sample for
   exact `mag_bp` and `mag_count`.

Assigns each sample to its SPIRE habitat (`habitat` level 1 and `habitat2` level 2) by the precedence
rule of card C-43 (`host-associated`, `anthropogenic`, `terrestrial`, `aquatic`, `air`, else `unlabelled`),
summarizes samples and assembled bases per habitat, and selects `target` samples (default 1,000) across
`habitat2` strata in proportion to `sqrt(N_h)` (Hamilton largest-remainder apportionment with floor >= 1
per non-empty stratum and cap <= N_h), using a deterministic per-stratum PRNG seeded with `f"{seed}:{h}"`.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import json
import math
import random
import statistics
import sys
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

HABITAT_ORDER = ("host-associated", "anthropogenic", "terrestrial", "aquatic", "air")
UNLABELLED = "unlabelled"


def parse_micron(cell: str) -> list[str]:
    """Parse SPIRE's Python-literal microntology list safely as data."""
    cell = (cell or "").strip()
    if not cell:
        return []
    val = ast.literal_eval(cell)
    if not isinstance(val, list) or not all(isinstance(x, str) for x in val):
        raise ValueError(f"invalid microntology cell: {cell[:80]!r}")
    return val


def spire_habitat(terms: list[str]) -> tuple[str, str]:
    """Return `(level1, level2)` habitat by card C-43's stated precedence."""
    tops = {t.split(":", 1)[0] for t in terms}
    for top in HABITAT_ORDER:
        if top in tops:
            seconds = sorted(":".join(t.split(":")[:2]) for t in terms if t.startswith(top + ":"))
            return top, (seconds[0] if seconds else top)
    return UNLABELLED, UNLABELLED


def allocate_sqrt(counts: dict[str, int], target: int) -> dict[str, int]:
    """Allocate `target` samples across non-empty strata in proportion to `sqrt(N_h)`.

    - Order-independent (sorts keys first).
    - If `target >= sum(counts.values())`, returns every available item `{h: N_h}`.
    - If `target >= len(non_empty)`, guarantees `1 <= k_h <= N_h` for every non-empty stratum.
    - Uses iterative water-filling + Hamilton largest-remainder rounding, breaking remainder ties by
      larger `N_h` and then lexicographically by stratum key.
    """
    if target < 0:
        raise ValueError(f"target must be >= 0, got {target}")
    clean = {k: int(v) for k, v in sorted(counts.items()) if int(v) > 0}
    if not clean or target == 0:
        return {k: 0 for k in clean}
    total_avail = sum(clean.values())
    if target >= total_avail:
        return dict(clean)

    keys = sorted(clean)
    n_strata = len(keys)
    if target < n_strata:
        # Fewer slots than strata: give 1 each to largest strata (ties broken by name).
        ranked = sorted(keys, key=lambda k: (-clean[k], k))
        return {k: (1 if i < target else 0) for i, k in enumerate(ranked)}

    # Start with floor 1 per stratum if N_h == 1 or if sqrt quota would fall below 1; water-fill caps/floors.
    weights = {k: math.sqrt(clean[k]) for k in keys}
    fixed: dict[str, int] = {}
    free = set(keys)
    while free:
        rem = target - sum(fixed.values())
        w_sum = sum(weights[k] for k in free)
        raw = {k: rem * weights[k] / w_sum for k in free}
        changed = False
        # Clamp strata whose raw share is <= 1.0 to 1, or >= N_h to N_h
        for k in sorted(free):
            if raw[k] >= clean[k]:
                fixed[k] = clean[k]
                free.remove(k)
                changed = True
            elif raw[k] < 1.0:
                fixed[k] = 1
                free.remove(k)
                changed = True
        if not changed:
            floors = {k: max(1, min(clean[k], int(math.floor(raw[k])))) for k in free}
            leftover = rem - sum(floors.values())
            order = sorted(
                free,
                key=lambda k: (-(raw[k] - floors[k]), -clean[k], k),
            )
            for k in order:
                if leftover <= 0:
                    break
                if floors[k] < clean[k]:
                    floors[k] += 1
                    leftover -= 1
            fixed.update(floors)
            break

    return {k: fixed[k] for k in keys}


def _opt_int(v) -> int | None:
    if v is None or v == "" or v == "NA" or v == "None":
        return None
    return int(float(v))


def _opt_float(v) -> float | None:
    if v is None or v == "" or v == "NA" or v == "None":
        return None
    return float(v)


def merge_sample_tables(
    study_rows: Iterable[dict],
    api_rows: Iterable[dict],
    mag_bp_by_sample: dict[str, tuple[int, int]],
) -> list[dict]:
    """Deterministically merge study TSV rows, API sample rows, and MAG genome sizes by `sample_id`."""
    by_study: dict[str, dict] = {}
    for r in study_rows:
        sid = r["sample_id"].strip()
        study = r.get("study", "").strip()
        mic = r.get("microntology", "")
        prev = by_study.get(sid)
        if prev is None or (study, mic) < (prev["study"], prev["microntology"]):
            by_study[sid] = {
                "sample_id": sid,
                "study": study,
                "microntology": mic,
                "study_mags": _opt_int(r.get("mags")) or 0,
                "lat": r.get("lat", "") or "",
                "lon": r.get("lon", "") or "",
            }

    by_api: dict[str, dict] = {}
    for r in api_rows:
        sid = (r.get("sample_name") or r.get("sample") or r.get("sample_id") or "").strip()
        if not sid:
            continue
        prev = by_api.get(sid)
        study = (r.get("study") or "").strip()
        if prev is None or study < prev["api_study"]:
            by_api[sid] = {
                "api_study": study,
                "assembly_bp": _opt_int(r.get("assembly_bp")),
                "assembly_n50": _opt_int(r.get("assembly_n50")),
                "binned_fraction": _opt_float(r.get("binned_fraction")),
                "read_bases_qc": _opt_int(r.get("bases")),
                "reads_qc": _opt_int(r.get("reads")),
                "n_mags_api": _opt_int(r.get("n_mags")),
                "n_hq_mags": _opt_int(r.get("n_hq_mags")),
                "n_mq_mags": _opt_int(r.get("n_mq_mags")),
                "microntology_api": (r.get("microntology") or r.get("environment") or "").strip(),
            }

    all_ids = sorted(set(by_study) | set(by_api))
    out: list[dict] = []
    for sid in all_ids:
        st = by_study.get(sid)
        ap = by_api.get(sid)
        if st is not None:
            terms = parse_micron(st["microntology"])
            study = st["study"]
        else:
            raw_m = ap["microntology_api"] if ap else ""
            if raw_m.startswith("["):
                terms = parse_micron(raw_m)
            elif raw_m:
                terms = [t.strip() for t in raw_m.replace(";", ",").split(",") if t.strip()]
            else:
                terms = []
            study = ap["api_study"] if ap else ""
        hab1, hab2 = spire_habitat(terms)
        n_mags_tbl, mag_bp = mag_bp_by_sample.get(sid, (0, 0))
        out.append({
            "sample_id": sid,
            "study": study,
            "in_study_tsv": st is not None,
            "in_samples_api": ap is not None,
            "habitat": hab1,
            "habitat2": hab2,
            "assembly_bp": ap["assembly_bp"] if ap else None,
            "assembly_n50": ap["assembly_n50"] if ap else None,
            "binned_fraction": ap["binned_fraction"] if ap else None,
            "read_bases_qc": ap["read_bases_qc"] if ap else None,
            "mag_count": n_mags_tbl,
            "mag_bp": mag_bp,
        })
    return out


def summarize_by_habitat(samples: Iterable[dict], key: str = "habitat") -> list[dict]:
    """Aggregate sample counts, assembled bases, and MAG/binned bases per habitat key."""
    groups: dict[str, list[dict]] = defaultdict(list)
    for r in samples:
        groups[r[key]].append(r)

    out: list[dict] = []
    for hab in sorted(groups):
        rows = groups[hab]
        asm_vals = [r["assembly_bp"] for r in rows if r["assembly_bp"] is not None]
        asm_sum = sum(asm_vals)
        mag_sum = sum(r["mag_bp"] for r in rows)
        bf_num = sum(
            r["assembly_bp"] * r["binned_fraction"]
            for r in rows
            if r["assembly_bp"] is not None and r["binned_fraction"] is not None
        )
        bf_den = sum(
            r["assembly_bp"]
            for r in rows
            if r["assembly_bp"] is not None and r["binned_fraction"] is not None
        )
        # Estimate total assembled bp when a minority of samples lack assembly_bp in metalog/API
        mean_asm = (asm_sum / len(asm_vals)) if asm_vals else 0.0
        est_total_bp = int(round(asm_sum + mean_asm * (len(rows) - len(asm_vals)))) if asm_vals else 0
        out.append({
            key: hab,
            "samples": len(rows),
            "samples_with_assembly_bp": len(asm_vals),
            "assembly_bp": asm_sum,
            "assembly_bp_estimated_all": est_total_bp,
            "mean_assembly_bp": round(mean_asm, 1) if asm_vals else None,
            "median_assembly_bp": float(statistics.median(asm_vals)) if asm_vals else None,
            "mag_count": sum(r["mag_count"] for r in rows),
            "mag_bp": mag_sum,
            "mag_bp_share_of_known_asm": round(mag_sum / asm_sum, 6) if asm_sum > 0 else None,
            "weighted_binned_fraction": round(bf_num / bf_den, 6) if bf_den > 0 else None,
        })
    return out


def select_pilot(
    samples: Iterable[dict],
    target: int = 1000,
    seed: int = 20261118,
    stratum_key: str = "habitat2",
) -> list[dict]:
    """Deterministic sqrt-proportional stratified sample selection."""
    by_stratum: dict[str, list[dict]] = defaultdict(list)
    for r in samples:
        by_stratum[r[stratum_key]].append(r)

    counts = {k: len(v) for k, v in by_stratum.items()}
    alloc = allocate_sqrt(counts, target=target)

    chosen: list[dict] = []
    for stratum in sorted(by_stratum):
        k_h = alloc.get(stratum, 0)
        if k_h <= 0:
            continue
        ordered = sorted(by_stratum[stratum], key=lambda r: r["sample_id"])
        rng = random.Random(f"{seed}:{stratum}")
        shuffled = list(ordered)
        rng.shuffle(shuffled)
        for rank, row in enumerate(shuffled[:k_h]):
            chosen.append({
                **row,
                "stratum": stratum,
                "stratum_n": len(ordered),
                "stratum_k": k_h,
                "rank_in_stratum": rank,
                "seed": seed,
            })
    chosen.sort(key=lambda r: (r["habitat"], r["habitat2"], r["rank_in_stratum"], r["sample_id"]))
    return chosen


def _load_study_tsvs(study_dir: Path, studies_html: Path) -> list[dict]:
    """Read `studies.html` study names -> `study_tsv/*.tsv` files mapped by URL-safe filename."""
    import re
    import urllib.parse

    html = studies_html.read_text()
    start = html.index("var studiesData")
    end = html.index("];", start)
    studies = [s for s, _ in re.findall(r'\{\s*study:\s*"((?:[^"\\]|\\.)*)",\s*samples:\s*(\d+)\s*\}', html[start:end])]
    rows: list[dict] = []
    for study in sorted(set(studies)):
        safe = urllib.parse.quote(study, safe="") + ".tsv"
        p = study_dir / safe
        if not p.is_file():
            raise FileNotFoundError(f"missing study TSV {p} for study {study!r}")
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows.append({"study": study, **r})
    return rows


def _load_api_samples(api_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for p in sorted(api_dir.glob("page_*.tsv")):
        with open(p, newline="") as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    return rows


def _load_mag_bp(genome_meta_gz: Path) -> dict[str, tuple[int, int]]:
    totals: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    with gzip.open(genome_meta_gz, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            sid = (r.get("derived_from_sample") or r.get("sample_id") or "").strip()
            gsize = _opt_int(r.get("genome_size")) or 0
            acc = totals[sid]
            acc[0] += 1
            acc[1] += gsize
    return {k: (v[0], v[1]) for k, v in totals.items()}


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    tmp = path.with_name(f".{path.name}.tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in fieldnames})
    tmp.replace(path)


def cmd_run(args: argparse.Namespace) -> int:
    tdir = Path(args.tables_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    study_rows = _load_study_tsvs(tdir / "study_tsv", tdir / "studies.html")
    api_rows = _load_api_samples(tdir / "api_samples")
    mag_bp = _load_mag_bp(tdir / "spire_v1_genome_metadata.tsv.gz")

    # Restrict universe to SPIRE's listed study samples (`in_study_tsv == True`) so every chosen sample
    # is a canonical SPIRE v1 study sample; also record any API-only rows in the summary.
    merged_all = merge_sample_tables(study_rows, api_rows, mag_bp)
    spire_samples = [r for r in merged_all if r["in_study_tsv"]]
    eligible_samples = [r for r in spire_samples if r["assembly_bp"] is None or r["assembly_bp"] > 0]

    hab1_all = summarize_by_habitat(spire_samples, key="habitat")
    hab2_all = summarize_by_habitat(spire_samples, key="habitat2")
    pilot = select_pilot(eligible_samples, target=args.target, seed=args.seed, stratum_key="habitat2")
    hab1_pilot = summarize_by_habitat(pilot, key="habitat")
    hab2_pilot = summarize_by_habitat(pilot, key="habitat2")

    h1_cols = [
        "habitat", "samples", "samples_with_assembly_bp", "assembly_bp", "assembly_bp_estimated_all",
        "mean_assembly_bp", "median_assembly_bp", "mag_count", "mag_bp",
        "mag_bp_share_of_known_asm", "weighted_binned_fraction",
    ]
    h2_cols = ["habitat2" if c == "habitat" else c for c in h1_cols]
    _write_tsv(out / "habitat_level1_all.tsv", hab1_all, h1_cols)
    _write_tsv(out / "habitat_level2_all.tsv", hab2_all, h2_cols)
    _write_tsv(out / "pilot_by_habitat1.tsv", hab1_pilot, h1_cols)
    _write_tsv(out / "pilot_by_habitat2.tsv", hab2_pilot, h2_cols)

    sel_cols = [
        "sample_id", "study", "habitat", "habitat2", "stratum_n", "stratum_k", "rank_in_stratum",
        "seed", "assembly_bp", "assembly_n50", "binned_fraction", "read_bases_qc", "mag_count", "mag_bp",
    ]
    _write_tsv(out / "pilot_selection.tsv", pilot, sel_cols)

    summary = {
        "seed": args.seed,
        "target": args.target,
        "stratum_key": "habitat2",
        "rule": "Hamilton largest-remainder proportional to sqrt(N_h) over SPIRE level-2 habitats, 1 <= k_h <= N_h, per-stratum PRNG random.Random(f'{seed}:{habitat2}') over sorted sample_ids",
        "study_tsv_samples": len(spire_samples),
        "eligible_nonempty_samples": len(eligible_samples),
        "zero_assembly_bp_samples": len(spire_samples) - len(eligible_samples),
        "api_samples_rows": len(api_rows),
        "merged_union_samples": len(merged_all),
        "api_only_samples": sum(1 for r in merged_all if not r["in_study_tsv"]),
        "study_samples_with_api_assembly_bp": sum(1 for r in spire_samples if r["assembly_bp"] is not None),
        "total_known_assembly_bp": sum(r["assembly_bp"] or 0 for r in spire_samples),
        "total_estimated_assembly_bp": sum(r["assembly_bp_estimated_all"] for r in hab2_all),
        "total_mag_count": sum(r["mag_count"] for r in spire_samples),
        "total_mag_bp": sum(r["mag_bp"] for r in spire_samples),
        "pilot_selected": len(pilot),
        "pilot_strata": len({r["habitat2"] for r in pilot}),
        "pilot_studies": len({r["study"] for r in pilot}),
        "pilot_known_assembly_bp": sum(r["assembly_bp"] or 0 for r in pilot),
        "pilot_estimated_assembly_bp": sum(r["assembly_bp_estimated_all"] for r in hab2_pilot),
        "pilot_mag_bp": sum(r["mag_bp"] for r in pilot),
    }
    (out / "select_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    run_p = sub.add_parser("run")
    run_p.add_argument("--tables-dir", required=True)
    run_p.add_argument("--out-dir", required=True)
    run_p.add_argument("--target", type=int, default=1000)
    run_p.add_argument("--seed", type=int, default=20261118)
    args = ap.parse_args(argv)
    if args.cmd == "run":
        return cmd_run(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())

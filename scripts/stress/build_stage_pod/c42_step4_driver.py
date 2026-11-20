"""Card C-42 step (4) resumable depth-1 (`--depth called --crispr`) driver.

    python -I c42_step4_driver.py run --spire-dir /tmp/if2453/s3/spire
        [--threads 35] [--scan-jobs 5] [--scan-cpus 7] [--max-scratch-gb 170]

Runs inside the pinned runtime image (`bio-box_faec1f4.sif`) under a `MemoryMax=128G` systemd scope:
1. Part A (`SP0`, unfiltered `--min-contig-len 0` across all 32 level-2 habitats, ~5-7 Gbp):
   measures exact unfiltered depth-1 CPU-min/Gbp, bytes/Gbp, genes/Gbp, and fragment-gene share
   (`dist_to_contig_end_bp == 0`) stratified across 6 contig-length bands (`200..499` through `>= 10,000 bp`),
   saves summaries to `depth1_unfiltered/`, then removes `root_unfiltered/` to reclaim scratch.
2. Part B (`SPA1..SPT2`, all 990 non-empty pilot samples at `--min-contig-len 1500`, `114.42 Gbp`):
   stages and builds one habitat stage at a time into `root/`, unlinking temporary `filtered/*.fa` and
   consumed raw `dl/assemblies/*.fa.gz` after each stage's `tapes.json` run record is published so total
   scratch stays flat (~85-125 GB, strictly < 180 GB).
3. Part C (`SPA1` search-pack probe, `DEC-54 (3)`):
   builds a scratch release + search pack on `SPA1` (`air`, `60.5 Mbp`) in `pack_probe/`, records exact
   `unique`, `igr`, and `dark_orf` bytes per Gbp, and removes `pack_probe/`.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PROFILE_TEMPLATE = """\
max_threads_day: 56
max_threads_night: 56
mem_scope: 180G
tmp: {tmp_dir}
durable: /stor/work/Wilke/if2453/all-the-data
bulk: /stor/scratch/Wilke/if2453/all-the-data
nice_monolithic: 19
nice_fanout: 19
day_start_hour: 8
day_end_hour: 17
cpu_min_per_process: 45
cluster_cpu_min_max: 20000
small_file_floor_bytes: 1048576
require_mem_cap: true
"""


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _log(msg: str) -> None:
    print(f"{_now()} {msg}", flush=True)


def _check_stop(spire_dir: Path) -> None:
    for p in (Path("/tmp/if2453/r10/root/STOP"), spire_dir / "STOP"):
        if p.exists():
            _log(f"STOP file found at {p}; exiting cleanly (75)")
            raise SystemExit(75)


def _dir_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        return path.stat().st_size
    total = 0
    for r, _, files in os.walk(path):
        for fn in files:
            try:
                total += os.path.getsize(os.path.join(r, fn))
            except OSError:
                pass
    return total


def _check_scratch(spire_dir: Path, cap_gb: float, label: str) -> float:
    used_gb = _dir_bytes(spire_dir) / 1e9
    _log(f"scratch [{label}]: {used_gb:.2f} GB (cap {cap_gb:.1f} GB)")
    if used_gb >= cap_gb:
        raise RuntimeError(f"scratch {used_gb:.2f} GB >= cap {cap_gb:.1f} GB at {label}")
    return used_gb


def _run_py(argv: list[str], label: str) -> None:
    _log(f"start {label}: {' '.join(argv)}")
    rc = subprocess.run([sys.executable, *argv], check=False).returncode
    if rc != 0:
        raise RuntimeError(f"{label} exited with code {rc}")
    _log(f"done {label}")


def _select_sp0_samples(sel_rows: list[dict], stats_by_s: dict[str, dict]) -> list[dict]:
    """Pick 1 representative non-empty sample per level-2 habitat (all 32 strata), preferring <= 400k contigs."""
    by_h2: dict[str, list[dict]] = defaultdict(list)
    for r in sel_rows:
        st = stats_by_s.get(r["sample_id"], {})
        if int(st.get("records", 0)) > 0:
            by_h2[r["habitat2"]].append(r)

    chosen: list[dict] = []
    for h2 in sorted(by_h2):
        cands = sorted(by_h2[h2], key=lambda x: (int(x["rank_in_stratum"]), x["sample_id"]))
        moderate = [c for c in cands if int(stats_by_s[c["sample_id"]]["records"]) <= 400_000]
        if moderate:
            chosen.append(moderate[0])
        else:
            smallest = min(cands, key=lambda x: (int(stats_by_s[x["sample_id"]]["records"]), x["sample_id"]))
            chosen.append(smallest)
    return chosen


def _write_single_stage_inputs_subproc(
    bin_dir: Path,
    samples: list[dict],
    fasta_dir: Path,
    fetch_log: Path,
    stage_name: str,
    habitat: str,
    out_dir: Path,
    min_contig_len: int,
    exclude_uids_bin: Path | None = None,
) -> dict:
    """Invoke `c42_inputs.py` in a child process so parent RSS/CPU stays near zero."""
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    helper_json = out_dir.parent / f".{stage_name}_req.json"
    helper_json.write_text(json.dumps({
        "samples": samples,
        "fasta_dir": str(fasta_dir),
        "fetch_log": str(fetch_log),
        "out_dir": str(out_dir),
        "min_contig_len": min_contig_len,
        "exclude_uids_bin": str(exclude_uids_bin) if exclude_uids_bin is not None else None,
    }))
    code = (
        "import json, sys; from pathlib import Path; "
        f"sys.path.insert(0, {str(bin_dir)!r}); import c42_inputs as m; "
        f"req = json.loads(Path({str(helper_json)!r}).read_text()); "
        "ret = {}; "
        "flog = Path(req['fetch_log']); "
        "[ret.__setitem__(r['sample_id'], r['fetched_utc']) for r in (json.loads(l) for l in flog.read_text().splitlines() if l.strip()) if r.get('sample_id') and r.get('fetched_utc')] if flog.is_file() else None; "
        "samps = [{**s, 'retrieved_at': ret.get(s['sample_id'], m.DEFAULT_RETRIEVED)} for s in req['samples']]; "
        "fastas = {s['sample_id']: Path(req['fasta_dir']) / f\"{s['sample_id']}.fa.gz\" for s in samps}; "
        "ex = Path(req['exclude_uids_bin']) if req.get('exclude_uids_bin') else None; "
        "m.build_stage_inputs(samps, fastas, Path(req['out_dir']), min_contig_len=req['min_contig_len'], exclude_uids_bin=ex)"
    )
    try:
        _run_py(["-c", code], f"inputs-{stage_name}")
    finally:
        helper_json.unlink(missing_ok=True)

    rep = json.loads((out_dir / "input_report.json").read_text())
    (out_dir / "samples.tsv").write_text(
        "sample_id\thabitat\thabitat2\tstudy\n"
        + "".join(f"{s['sample_id']}\t{s['habitat']}\t{s['habitat2']}\t{s['study']}\n" for s in samples)
    )
    return {"stage": stage_name, "habitat": habitat, **rep}


def _run_pack_probe_spa1(spire_dir: Path, profile_path: Path, root: Path, threads: int) -> None:
    out_json = spire_dir / "depth1" / "pack_probe_SPA1.json"
    if out_json.is_file():
        _log("pack_probe_SPA1.json already exists; skipping Part C")
        return
    probe = spire_dir / "pack_probe"
    if probe.exists():
        shutil.rmtree(probe)
    p_root = probe / "root"
    p_packs = probe / "packs"
    p_packs.mkdir(parents=True, exist_ok=True)
    # Copy SPA1 stage files (~45 MB total) into isolated probe root
    shutil.copytree(root / "l0" / "l0_manifest" / "stage=SPA1", p_root / "l0" / "l0_manifest" / "stage=SPA1")
    if (root / "l0seq" / "store_params.json").is_file():
        (p_root / "l0seq").mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / "l0seq" / "store_params.json", p_root / "l0seq" / "store_params.json")
    for sub in ("contig_index", "contig_label", "contig_occurrence", "genome_input", "seq_locator", "seq_shards"):
        src = root / "l0seq" / sub / "stage=SPA1"
        if src.exists():
            shutil.copytree(src, p_root / "l0seq" / sub / "stage=SPA1")
    for sub_glob in (
        ("l0seq/igr", "igr_version=*/stage=SPA1"),
        ("l0seq/protein_seq", "gene_call_version=*/stage=SPA1"),
        ("l1/feature", "producer=*/stage=SPA1"),
        ("l1/feature_protein", "gene_call_version=*/stage=SPA1"),
        ("l1/protein_occurrence", "gene_call_version=*/stage=SPA1"),
        ("l1/crispr", "producer=*/stage=SPA1"),
        ("l3/tapes", "tape_version=*/stage=SPA1"),
        ("l3/tape_version_def", "tape_version=*"),
    ):
        base = root / sub_glob[0]
        for match in base.glob(sub_glob[1]):
            rel = match.relative_to(root)
            shutil.copytree(match, p_root / rel)
    shutil.copytree(root / "ledger" / "stage=SPA1", p_root / "ledger" / "stage=SPA1")
    shutil.copytree(root / "runs" / "stage=SPA1", p_root / "runs" / "stage=SPA1")

    in_spa1 = spire_dir / "inputs" / "SPA1"
    try:
        _run_py([
            "-m", "biobox.cli", "build-stage",
            "--profile", str(profile_path),
            "--root", str(p_root),
            "--stage", "SPA1",
            "--depth", "called",
            "--release", "rSPA1",
            "--stages", "SPA1",
            "--genomes", str(in_spa1 / "genomes.json"),
            "--objects", str(in_spa1 / "objects.json"),
            "--crispr",
            "--threads", str(min(threads, 16)),
            "--scan-jobs", "2",
            "--scan-cpus", "4",
        ], "pack-probe-release-rSPA1")

        _run_py([
            "-m", "biobox.cli", "pack", "build",
            "--profile", str(profile_path),
            "--root", str(p_root),
            "--release", "rSPA1",
            "--pack-dir", str(p_packs),
        ], "pack-probe-build-rSPA1")

        # Measure pack target sizes vs SPA1 stored_bp
        d1_sum = json.loads((spire_dir / "depth1" / "depth1_summary.json").read_text())
        spa1_info = next(s for s in d1_sum["per_stage"] if s["stage"] == "SPA1")
        gbp = spa1_info["stored_bp"] / 1e9
        pack_dirs = list((p_packs / "release=rSPA1").iterdir())
        pdir = pack_dirs[0]
        by_file = {f.name: _dir_bytes(f) for f in sorted(pdir.iterdir())}
        tot_pack = _dir_bytes(pdir)
        man = json.loads((pdir / "manifest.json").read_text()) if (pdir / "manifest.json").is_file() else {}
        probe_summary = {
            "stage": "SPA1",
            "stored_bp": spa1_info["stored_bp"],
            "pack_id": pdir.name,
            "pack_total_bytes": tot_pack,
            "pack_gb_per_gbp": round((tot_pack / 1e9) / gbp, 4) if gbp else 0.0,
            "entries_bytes": by_file,
            "entries_gb_per_gbp": {k: round((v / 1e9) / gbp, 4) for k, v in by_file.items()} if gbp else {},
            "manifest_targets": man.get("targets", {}),
        }
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(probe_summary, indent=2) + "\n")
        _log(f"pack probe SPA1: {json.dumps(probe_summary)}")
    finally:
        shutil.rmtree(probe, ignore_errors=True)


def cmd_run(args: argparse.Namespace) -> int:
    spire_dir = Path(args.spire_dir).resolve()
    bin_dir = spire_dir / "bin"
    sys.path.insert(0, str(bin_dir))
    import c42_inputs as m_in

    tmp_dir = spire_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    profile_path = spire_dir / "profile.yaml"
    profile_path.write_text(PROFILE_TEMPLATE.format(tmp_dir=tmp_dir))

    sel_rows = list(csv.DictReader(open(spire_dir / "select" / "pilot_selection.tsv", newline=""), delimiter="\t"))
    stats_by_s = {
        r["sample_id"]: r
        for r in csv.DictReader(open(spire_dir / "contigs" / "sample_contig_stats.tsv", newline=""), delimiter="\t")
    }
    fasta_dir = spire_dir / "dl" / "assemblies"
    fetch_log = fasta_dir / "fetch_assemblies.jsonl"

    # --- Part A: Unfiltered (`--min-contig-len 0`) calibration stage SP0 across all 32 strata ---
    unf_summary = spire_dir / "depth1_unfiltered" / "depth1_summary.json"
    if not unf_summary.is_file():
        _check_stop(spire_dir)
        _check_scratch(spire_dir, args.max_scratch_gb, "before-SP0")
        sp0_samples = _select_sp0_samples(sel_rows, stats_by_s)
        in_unf = spire_dir / "inputs_unfiltered" / "SP0"
        root_unf = spire_dir / "root_unfiltered"
        if not (in_unf / "genomes.json").is_file():
            if in_unf.exists():
                shutil.rmtree(in_unf)
            sp0_rep = _write_single_stage_inputs_subproc(
                bin_dir, sp0_samples, fasta_dir, fetch_log, "SP0", "all_32_strata", in_unf, min_contig_len=0
            )
            (spire_dir / "inputs_unfiltered" / "stages_manifest.json").write_text(json.dumps([sp0_rep], indent=2) + "\n")

        _check_stop(spire_dir)
        _run_py([
            "-m", "biobox.cli", "build-stage",
            "--profile", str(profile_path),
            "--root", str(root_unf),
            "--stage", "SP0",
            "--depth", "called",
            "--genomes", str(in_unf / "genomes.json"),
            "--objects", str(in_unf / "objects.json"),
            "--threads", str(args.threads),
            "--scan-jobs", str(args.scan_jobs),
            "--scan-cpus", str(args.scan_cpus),
        ], "build-stage-SP0-unfiltered")

        _run_py([
            str(bin_dir / "c42_inputs.py"), "summarize-depth1",
            "--scratch-root", str(root_unf),
            "--stages-manifest", str(spire_dir / "inputs_unfiltered" / "stages_manifest.json"),
            "--out-dir", str(spire_dir / "depth1_unfiltered"),
        ], "summarize-SP0-unfiltered")

        shutil.rmtree(root_unf, ignore_errors=True)
        shutil.rmtree(tmp_dir / "scratch", ignore_errors=True)
        _check_scratch(spire_dir, args.max_scratch_gb, "after-SP0-cleanup")

    # --- Part B: Full 990-sample pilot at `--min-contig-len 1500` across habitat stages ---
    stages = m_in.partition_samples_into_stages(
        sel_rows,
        stats_by_s,
        min_contig_len=1500,
        max_stage_bp=19_000_000_000,
        max_stage_contigs=5_500_000,
    )
    _log(f"planned {len(stages)} stages at min_contig_len=1500: {[s[0] for s in stages]}")

    roots_dir = spire_dir / "roots"
    roots_dir.mkdir(parents=True, exist_ok=True)
    inputs_base = spire_dir / "inputs"
    inputs_base.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict] = []
    completed_stages: list[str] = []

    for st_name, hab, st_samples in stages:
        _check_stop(spire_dir)
        _check_scratch(spire_dir, args.max_scratch_gb, f"before-{st_name}")
        st_in = inputs_base / st_name
        st_root = roots_dir / st_name
        tapes_rec = st_root / "runs" / f"stage={st_name}" / "tapes.json"

        if not tapes_rec.is_file():
            owned_bin = inputs_base / ".earlier_owned_uids.bin"
            try:
                if completed_stages:
                    n_owned = m_in.write_owned_uids_bin(roots_dir, completed_stages, owned_bin)
                    _log(f"wrote {n_owned} earlier-owned contig_uids from {completed_stages} to {owned_bin.name}")
                else:
                    owned_bin.unlink(missing_ok=True)

                if (
                    not (st_in / "genomes.json").is_file()
                    or not (st_in / "filtered").is_dir()
                    or not (st_in / "cross_stage_occurrences.tsv").is_file()
                ):
                    if st_in.exists():
                        shutil.rmtree(st_in)
                    rep = _write_single_stage_inputs_subproc(
                        bin_dir,
                        st_samples,
                        fasta_dir,
                        fetch_log,
                        st_name,
                        hab,
                        st_in,
                        min_contig_len=1500,
                        exclude_uids_bin=owned_bin if owned_bin.is_file() else None,
                    )
                else:
                    rep = {"stage": st_name, "habitat": hab, **json.loads((st_in / "input_report.json").read_text())}
            finally:
                owned_bin.unlink(missing_ok=True)

            _check_stop(spire_dir)
            _check_scratch(spire_dir, args.max_scratch_gb, f"inputs-ready-{st_name}")
            # Build each stage in its own isolated scratch root `roots/<STAGE>` (`<= 45` row groups `<< 64`)
            # so KeyLookup's 64-row-group LRU cache has 100% hit rate on random 1-contig genome_uid lookups.
            # Run --crispr on stages with <= 350k contigs (SPA1, SPU1); larger stages (> 350k contigs)
            # avoid crispr.py line 566's 1-row-per-row-group Parquet footer bug (measured on SP0: 3.93M row groups,
            # 2.19 GB contigs.parquet, 29.9 GB RSS).
            crispr_flag = ["--crispr"] if int(rep["distinct_contigs"]) <= 350_000 else []
            _run_py([
                "-m", "biobox.cli", "build-stage",
                "--profile", str(profile_path),
                "--root", str(st_root),
                "--stage", st_name,
                "--depth", "called",
                "--genomes", str(st_in / "genomes.json"),
                "--objects", str(st_in / "objects.json"),
                *crispr_flag,
                "--threads", str(args.threads),
                "--scan-jobs", str(args.scan_jobs),
                "--scan-cpus", str(args.scan_cpus),
            ], f"build-stage-{st_name}")
            if not tapes_rec.is_file():
                raise RuntimeError(f"missing {tapes_rec} after build-stage {st_name}")
        else:
            rep = {"stage": st_name, "habitat": hab, **json.loads((st_in / "input_report.json").read_text())}
            _log(f"stage {st_name} already completed ({tapes_rec} present)")

        completed_stages.append(st_name)
        manifest_rows.append(rep)
        (inputs_base / "stages_manifest.json").write_text(json.dumps(manifest_rows, indent=2) + "\n")

        # Reclaim temporary filtered FASTAs and consumed raw `.fa.gz` files (keep tiny SPA1 intact for pack probe)
        shutil.rmtree(tmp_dir / "scratch", ignore_errors=True)
        if st_name != "SPA1":
            shutil.rmtree(st_in / "filtered", ignore_errors=True)
            for s in st_samples:
                (fasta_dir / f"{s['sample_id']}.fa.gz").unlink(missing_ok=True)
        _check_scratch(spire_dir, args.max_scratch_gb, f"after-{st_name}-cleanup")

    _run_py([
        str(bin_dir / "c42_inputs.py"), "summarize-depth1",
        "--scratch-root", str(roots_dir),
        "--stages-manifest", str(inputs_base / "stages_manifest.json"),
        "--out-dir", str(spire_dir / "depth1"),
    ], "summarize-depth1-all-stages")

    # --- Part C: Search-pack bytes/Gbp probe on SPA1 (DEC-54 (3)) ---
    _check_stop(spire_dir)
    _run_pack_probe_spa1(spire_dir, profile_path, roots_dir / "SPA1", args.threads)

    (spire_dir / "depth1" / "STEP4_DONE").write_text(_now() + "\n")
    _log("STEP 4 ALL PARTS COMPLETE")
    return 0


def cmd_migrate_shared_root(args: argparse.Namespace) -> int:
    """Move already-completed stages out of legacy shared `root/` into per-stage `roots/<STAGE>/`."""
    from biobox import stage as bstage

    base = Path(args.spire_dir).resolve()
    old = base / "root"
    roots = base / "roots"
    roots.mkdir(parents=True, exist_ok=True)
    for st, has_cr in (("SPA1", True), ("SPN1", False)):
        dst = roots / st
        if old.is_dir() and (old / "runs" / f"stage={st}" / "tapes.json").is_file():
            dst.mkdir(parents=True, exist_ok=True)
            if (old / "l0seq" / "store_params.json").is_file():
                (dst / "l0seq").mkdir(parents=True, exist_ok=True)
                shutil.copy2(old / "l0seq" / "store_params.json", dst / "l0seq" / "store_params.json")
            if (old / "l3" / "tape_version_def").is_dir():
                shutil.copytree(old / "l3" / "tape_version_def", dst / "l3" / "tape_version_def", dirs_exist_ok=True)
            for p in sorted(old.rglob(f"stage={st}")):
                if p.is_dir():
                    rel = p.relative_to(old)
                    target = dst / rel
                    target.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(p, target)
        if (dst / "runs" / f"stage={st}" / "tapes.json").is_file():
            problems = bstage.stage_incomplete(dst, st, "called", crispr=has_cr)
            if problems:
                raise RuntimeError(f"{st} incomplete after migration: {problems}")
            _log(f"verified {st} complete at {dst}")
    if old.is_dir():
        shutil.rmtree(old)
    shutil.rmtree(base / "tmp", ignore_errors=True)
    _log("removed legacy shared root/ and tmp/")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--spire-dir", required=True)
    p.add_argument("--threads", type=int, default=28)
    p.add_argument("--scan-jobs", type=int, default=4)
    p.add_argument("--scan-cpus", type=int, default=7)
    p.add_argument("--max-scratch-gb", type=float, default=170.0)
    pm = sub.add_parser("migrate-shared-root")
    pm.add_argument("--spire-dir", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "run":
        return cmd_run(args)
    if args.cmd == "migrate-shared-root":
        return cmd_migrate_shared_root(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())

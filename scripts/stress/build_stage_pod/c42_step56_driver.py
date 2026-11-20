"""Card C-42 steps (5) and (6) driver at 56 threads under MemoryMax=180G.

    python c42_step56_driver.py run --spire-dir /tmp/if2453/s3/spire
        --ref-root /tmp/if2453/r10/root --ref-release r26.11.1
        --pfam-hmm /stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz
        --threads 56 --align-parts 20 --max-scratch-gb 178

Checks `/tmp/if2453/r10/root/STOP` and `<spire-dir>/STOP` before every sub-step, keeps every
single child process strictly `< 280` CPU-min (`--align-parts 20` on MMseqs2 `align`, `<= 8M`
residues per Pfam round shard), and cleans bulky intermediate MMseqs2/Parquet files as soon
as summary tables are verified.
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
from pathlib import Path


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


def cmd_run(args: argparse.Namespace) -> int:
    spire_dir = Path(args.spire_dir).resolve()
    bin_dir = spire_dir / "bin"
    roots_dir = spire_dir / "roots"
    inputs_dir = spire_dir / "inputs"
    nov_dir = spire_dir / "novelty"
    tniq_dir = spire_dir / "tniq"

    # Reclaim SPA1 temporary FASTAs once pack_probe_SPA1.json is written
    if (spire_dir / "depth1" / "pack_probe_SPA1.json").is_file():
        shutil.rmtree(inputs_dir / "SPA1" / "filtered", ignore_errors=True)
        spa1_tsv = inputs_dir / "SPA1" / "samples.tsv"
        if spa1_tsv.is_file():
            with open(spa1_tsv, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    (spire_dir / "dl" / "assemblies" / f"{r['sample_id']}.fa.gz").unlink(missing_ok=True)

    stages_man = json.loads((inputs_dir / "stages_manifest.json").read_text())
    pilot_stages = ",".join(s["stage"] for s in stages_man)

    # --- Step 5: Novelty vs r26.11.1 via 50%-id MMseqs2 linclust ---
    _check_stop(spire_dir)
    _check_scratch(spire_dir, args.max_scratch_gb, "before-step5")
    nov_dir.mkdir(parents=True, exist_ok=True)

    if not (nov_dir / "export_meta.json").is_file():
        _run_py([
            str(bin_dir / "c42_novelty.py"), "export-faa",
            "--ref-root", args.ref_root,
            "--ref-release", args.ref_release,
            "--pilot-root", str(roots_dir),
            "--pilot-stages", pilot_stages,
            "--work-dir", str(nov_dir),
            "--direct-createdb",
        ], "step5-export-faa")
        _check_scratch(spire_dir, args.max_scratch_gb, "after-step5-export")

    _check_stop(spire_dir)
    if not (nov_dir / "linclust_steps.json").is_file():
        _run_py([
            str(bin_dir / "c42_novelty.py"), "linclust-50",
            "--work-dir", str(nov_dir),
            "--threads", str(args.threads),
            "--align-parts", str(args.align_parts),
            "--split-mem", args.split_mem,
            "--cleanup-db",
        ], "step5-linclust-50")
        _check_scratch(spire_dir, args.max_scratch_gb, "after-step5-linclust")

    _check_stop(spire_dir)
    if not (nov_dir / "novelty_summary.json").is_file():
        _run_py([
            str(bin_dir / "c42_novelty.py"), "analyze",
            "--work-dir", str(nov_dir),
            "--pilot-root", str(roots_dir),
            "--inputs-dir", str(inputs_dir),
            "--selection", str(spire_dir / "select" / "pilot_selection.tsv"),
            "--threads", str(min(args.threads, 28)),
            "--duck-mem", "48GB",
            "--cleanup-tsv",
        ], "step5-analyze")

    # Verify Step 5 outputs and remove bulky intermediate Parquet files
    req_nov = [
        "export_meta.json",
        "linclust_steps.json",
        "novelty_summary.json",
        "novelty_by_habitat1.tsv",
        "novelty_by_habitat2.tsv",
        "accumulation_global.tsv",
        "accumulation_by_habitat1.tsv",
    ]
    for fn in req_nov:
        if not (nov_dir / fn).is_file():
            raise RuntimeError(f"missing Step 5 output: {nov_dir / fn}")
    for bulky in (
        "clu50.parquet",
        "pilot_novel_map.parquet",
        "sample_novel_cids.parquet",
        "ref_uids.parquet",
        "pilot_uids.parquet",
        "cross_stage_occ.parquet",
    ):
        (nov_dir / bulky).unlink(missing_ok=True)
    _check_scratch(spire_dir, args.max_scratch_gb, "after-step5-cleanup")

    # --- Step 6: End-to-end TniQ search + +-10-gene Pfam 38.2 round + SPIRE MAG check ---
    _check_stop(spire_dir)
    tniq_dir.mkdir(parents=True, exist_ok=True)

    if not (tniq_dir / "tniq_summary.json").is_file():
        _run_py([
            str(bin_dir / "c42_tniq.py"), "run",
            "--pilot-root", str(roots_dir),
            "--inputs-dir", str(inputs_dir),
            "--selection", str(spire_dir / "select" / "pilot_selection.tsv"),
            "--pfam-hmm", args.pfam_hmm,
            "--out-dir", str(tniq_dir),
            "--threads", str(args.threads),
        ], "step6-tniq-initial")

    _check_stop(spire_dir)
    mag_json = tniq_dir / "hit_sample_mag_uids.json"
    if not mag_json.is_file() and (tniq_dir / "hit_samples.txt").is_file():
        _run_py([
            str(bin_dir / "c42_fetch.py"), "sample-mags",
            "--samples-file", str(tniq_dir / "hit_samples.txt"),
            "--genome-meta", str(spire_dir / "tables" / "spire_v1_genome_metadata.tsv.gz"),
            "--spire-dir", str(spire_dir),
            "--out-json", str(mag_json),
            "--interval", "0.5",
            "--unlink-after",
        ], "step6-fetch-hit-sample-mags")

    _check_stop(spire_dir)
    if mag_json.is_file():
        cur_sum = json.loads((tniq_dir / "tniq_summary.json").read_text()) if (tniq_dir / "tniq_summary.json").is_file() else {}
        if cur_sum.get("tniq_hit_loci_in_spire_mags") is None:
            _run_py([
                str(bin_dir / "c42_tniq.py"), "run",
                "--pilot-root", str(roots_dir),
                "--inputs-dir", str(inputs_dir),
                "--selection", str(spire_dir / "select" / "pilot_selection.tsv"),
                "--pfam-hmm", args.pfam_hmm,
                "--out-dir", str(tniq_dir),
                "--mag-uids-json", str(mag_json),
                "--threads", str(args.threads),
            ], "step6-tniq-with-mag-uids")

    (spire_dir / "STEP56_DONE").write_text(_now() + "\n")
    _log("STEPS 5 AND 6 COMPLETE")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--spire-dir", required=True)
    p.add_argument("--ref-root", default="/tmp/if2453/r10/root")
    p.add_argument("--ref-release", default="r26.11.1")
    p.add_argument("--pfam-hmm", default="/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz")
    p.add_argument("--threads", type=int, default=56)
    p.add_argument("--align-parts", type=int, default=20)
    p.add_argument("--split-mem", default="110G")
    p.add_argument("--max-scratch-gb", type=float, default=178.0)
    args = ap.parse_args(argv)
    if args.cmd == "run":
        return cmd_run(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())

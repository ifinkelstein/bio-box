"""Card B-63: meta-mode gene-calling cost and batch sizing before PL3.

Controlled, interleaved scratch measurements on samples of IMG-PR, mMGE, gut and RefSeq
from /tmp/if2453/s2/inputs/PL2 (bound read-only). Writes only under --work-dir
(default /tmp/if2453/b63). Never writes to the data root or /stor.

Subcommands:
  arm-a   Measure GeneFinder(meta=True) / ViralGeneFinder(meta=True) per-genome
          instantiation overhead vs per-thread/process reuse, stratified by source
          and contig length (1-5 kb, 10-30 kb, >=50 kb), plus unstratified source
          samples and exact call-for-call equality check.
  arm-b   Ingest a representative 1-in-N IMG-PR + mMGE slice (~80 Mbp at N=200) into
          a scratch root, compute cut_plan's estimate, and run `stage genecall`
          (4 workers, real subprocess) to measure total step CPU-minutes (self +
          children rusage), wall time, peak RSS, and `_batches` fragmentation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_PL2 = Path("/tmp/if2453/s2/inputs/PL2")
DEFAULT_WORK = Path("/tmp/if2453/b63")

SOURCE_FASTAS = {
    "IMG-PR": ("plasmidscope_IMG-PR", "plasmidscope_IMG-PR.fasta"),
    "mMGE": ("plasmidscope_mMGE", "plasmidscope_mMGE.fasta"),
    "gut": ("gut_plasmids_yu2024", "gut_plasmids_yu2024.fasta"),
    "RefSeq": ("plasmidscope_RefSeq", "plasmidscope_RefSeq.fasta"),
}

LENGTH_BINS = (
    ("short_1_5kb", 1_000, 5_000),
    ("med_10_30kb", 10_000, 30_000),
    ("long_50kb+", 50_000, 10**12),
)


def _load_avg() -> list[float]:
    try:
        return [round(x, 2) for x in os.getloadavg()]
    except OSError:
        return []


def _call_seq_list(seqs: list[tuple[str, str, bool]], mode: str, reuse: bool) -> tuple[float, int, str]:
    """Call `seqs` = [(name, seq, circular), ...] in `mode` ('meta' or 'gv-meta').

    Returns (cpu_s, total_genes, blake3_digest_of_all_calls).
    """
    from blake3 import blake3
    import pyrodigal
    import pyrodigal_gv
    from biobox.l1.genecall import _genes

    t0 = time.process_time()
    h = blake3()
    n_genes = 0
    cached = None
    if reuse:
        cached = pyrodigal_gv.ViralGeneFinder(meta=True) if mode == "gv-meta" else pyrodigal.GeneFinder(meta=True)
    for name, seq, circular in seqs:
        if reuse:
            finder = cached
        else:
            finder = pyrodigal_gv.ViralGeneFinder(meta=True) if mode == "gv-meta" else pyrodigal.GeneFinder(meta=True)
        genes = sorted(_genes(finder, seq, circular), key=lambda t: (t[0], t[1], t[2]))
        n_genes += len(genes)
        for start, end, strand, aa, code, partial in genes:
            h.update(f"{name}:{start}:{end}:{strand}:{code}:{int(partial)}:{aa}\n".encode())
    cpu_s = time.process_time() - t0
    return cpu_s, n_genes, h.hexdigest(length=16)


def _sample_source_bins(
    fasta_path: Path,
    per_bin_target_bp: int,
    per_bin_max_seqs: int,
    unstrat_every: int,
    unstrat_max_bp: int,
) -> tuple[dict[str, list[tuple[str, str, bool]]], list[tuple[str, str, bool]]]:
    """Stream `fasta_path` once to collect length-stratified bins and an unstratified sample."""
    from biobox import fasta

    bins: dict[str, list[tuple[str, str, bool]]] = {name: [] for name, _, _ in LENGTH_BINS}
    bin_bp: dict[str, int] = {name: 0 for name, _, _ in LENGTH_BINS}
    unstrat: list[tuple[str, str, bool]] = []
    unstrat_bp = 0

    for idx, rec in enumerate(fasta.iter_records(fasta_path)):
        seq_str = rec.seq.upper().decode("ascii")
        n = len(seq_str)
        if idx % unstrat_every == 0 and unstrat_bp < unstrat_max_bp:
            unstrat.append((rec.header.split()[0], seq_str, False))
            unstrat_bp += n
        for bname, lo, hi in LENGTH_BINS:
            if lo <= n < hi and bin_bp[bname] < per_bin_target_bp and len(bins[bname]) < per_bin_max_seqs:
                bins[bname].append((rec.header.split()[0], seq_str, False))
                bin_bp[bname] += n
                break
        if unstrat_bp >= unstrat_max_bp and all(
            bin_bp[b] >= per_bin_target_bp or len(bins[b]) >= per_bin_max_seqs for b, _, _ in LENGTH_BINS
        ):
            break
    return bins, unstrat


def run_arm_a(args: argparse.Namespace) -> None:
    import pyrodigal
    import pyrodigal_gv

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    pl2 = Path(args.pl2_dir)

    # 1. Measure standalone instantiation cost
    t0 = time.process_time()
    for _ in range(args.init_iters):
        pyrodigal.GeneFinder(meta=True)
    meta_init_us = (time.process_time() - t0) * 1e6 / args.init_iters

    t0 = time.process_time()
    for _ in range(args.init_iters):
        pyrodigal_gv.ViralGeneFinder(meta=True)
    gv_init_us = (time.process_time() - t0) * 1e6 / args.init_iters

    print(
        f"Standalone init ({args.init_iters} iters): "
        f"GeneFinder(meta=True) = {meta_init_us:.1f} us, "
        f"ViralGeneFinder(meta=True) = {gv_init_us:.1f} us",
        flush=True,
    )

    # 2. Sample sources
    sources = args.sources.split(",")
    samples: dict[str, dict[str, list[tuple[str, str, bool]]]] = {}
    for src in sources:
        _, fname = SOURCE_FASTAS[src]
        fpath = pl2 / fname
        bins, unstrat = _sample_source_bins(
            fpath,
            per_bin_target_bp=args.bin_bp,
            per_bin_max_seqs=args.bin_max_seqs,
            unstrat_every=args.unstrat_every,
            unstrat_max_bp=args.unstrat_bp,
        )
        samples[src] = {**bins, "unstrat": unstrat}
        summary = ", ".join(
            f"{k}: {len(v)} seqs / {sum(len(s) for _, s, _ in v)/1e6:.2f} Mbp"
            for k, v in samples[src].items()
        )
        print(f"Sampled {src}: {summary}", flush=True)

    # 3. Interleaved reps: (mode, reuse)
    modes_to_test = args.modes.split(",")
    records: list[dict] = []
    digests: dict[tuple[str, str, str], str] = {}

    for rep in range(1, args.reps + 1):
        for src in sources:
            for bname, seqs in samples[src].items():
                if not seqs:
                    continue
                total_bp = sum(len(s) for _, s, _ in seqs)
                for mode in modes_to_test:
                    order = [False, True] if rep % 2 == 1 else [True, False]
                    for reuse in order:
                        load_before = _load_avg()
                        cpu_s, n_genes, dig = _call_seq_list(seqs, mode, reuse)
                        load_after = _load_avg()
                        key = (src, bname, mode)
                        if key in digests:
                            if digests[key] != dig:
                                raise RuntimeError(
                                    f"Output mismatch for {key} reuse={reuse}: {dig} != {digests[key]}"
                                )
                        else:
                            digests[key] = dig
                        rec = {
                            "rep": rep,
                            "source": src,
                            "bin": bname,
                            "mode": mode,
                            "reuse": reuse,
                            "n_seqs": len(seqs),
                            "total_bp": total_bp,
                            "mean_bp": round(total_bp / len(seqs), 1),
                            "n_genes": n_genes,
                            "cpu_s": round(cpu_s, 6),
                            "cpu_s_per_mbp": round(cpu_s / (total_bp / 1e6), 4),
                            "cpu_ms_per_genome": round(cpu_s * 1e3 / len(seqs), 4),
                            "digest": dig,
                            "load_before": load_before,
                            "load_after": load_after,
                        }
                        records.append(rec)

    # 4. Also verify circular equality on a small subset (first 10 seqs of mMGE short & med)
    circ_check = []
    for bname in ("short_1_5kb", "med_10_30kb"):
        sub = [(n, s, True) for n, s, _ in samples.get("mMGE", {}).get(bname, [])[:10]]
        if not sub:
            continue
        for mode in modes_to_test:
            _, ng_f, dig_f = _call_seq_list(sub, mode, False)
            _, ng_r, dig_r = _call_seq_list(sub, mode, True)
            if dig_f != dig_r or ng_f != ng_r:
                raise RuntimeError(f"Circular mismatch on {bname} {mode}: {dig_f} vs {dig_r}")
            circ_check.append({"bin": bname, "mode": mode, "n_seqs": len(sub), "n_genes": ng_f, "digest": dig_f})

    # 5. Aggregate summary
    agg: list[dict] = []
    for src in sources:
        for bname in (*[b for b, _, _ in LENGTH_BINS], "unstrat"):
            for mode in modes_to_test:
                fresh = [r for r in records if r["source"] == src and r["bin"] == bname and r["mode"] == mode and not r["reuse"]]
                reused = [r for r in records if r["source"] == src and r["bin"] == bname and r["mode"] == mode and r["reuse"]]
                if not fresh or not reused:
                    continue
                f_mbp = statistics.median(r["cpu_s_per_mbp"] for r in fresh)
                r_mbp = statistics.median(r["cpu_s_per_mbp"] for r in reused)
                f_ms = statistics.median(r["cpu_ms_per_genome"] for r in fresh)
                r_ms = statistics.median(r["cpu_ms_per_genome"] for r in reused)
                agg.append({
                    "source": src,
                    "bin": bname,
                    "mode": mode,
                    "n_seqs": fresh[0]["n_seqs"],
                    "mean_bp": fresh[0]["mean_bp"],
                    "total_mbp": round(fresh[0]["total_bp"] / 1e6, 3),
                    "fresh_s_per_mbp": round(f_mbp, 4),
                    "reused_s_per_mbp": round(r_mbp, 4),
                    "delta_pct": round((r_mbp - f_mbp) / f_mbp * 100.0, 2),
                    "fresh_ms_per_genome": round(f_ms, 4),
                    "reused_ms_per_genome": round(r_ms, 4),
                    "saved_ms_per_genome": round(f_ms - r_ms, 4),
                })

    out = {
        "init_iters": args.init_iters,
        "meta_init_us": round(meta_init_us, 2),
        "gv_init_us": round(gv_init_us, 2),
        "reps": args.reps,
        "exact_equality_verified": True,
        "circular_equality_check": circ_check,
        "summary": agg,
        "records": records,
    }
    out_path = work / "arm_a.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nWrote {out_path}", flush=True)
    print(f"{'source':<8} {'bin':<12} {'mode':<8} {'n':>5} {'mean_bp':>8} {'fresh_s/M':>10} {'reuse_s/M':>10} {'delta%':>8} {'save_ms/g':>10}")
    for row in agg:
        print(
            f"{row['source']:<8} {row['bin']:<12} {row['mode']:<8} {row['n_seqs']:>5} "
            f"{row['mean_bp']:>8.0f} {row['fresh_s_per_mbp']:>10.4f} {row['reused_s_per_mbp']:>10.4f} "
            f"{row['delta_pct']:>+7.2f}% {row['saved_ms_per_genome']:>10.4f}"
        )


def _tree_sha256(root: Path) -> dict[str, str]:
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def _build_arm_b_slice(pl2: Path, slice_dir: Path, every: int, dtr_circular: bool) -> dict:
    """Create a 1-in-`every` IMG-PR + mMGE slice (FASTA + genomes.json + objects.json) under `slice_dir`."""
    from biobox import fasta, ids

    slice_dir.mkdir(parents=True, exist_ok=True)
    wanted_sources = {"plasmidscope_IMG-PR", "plasmidscope_mMGE"}

    # Read records.tsv to pick every Nth new (unowned, first-occurrence) record of IMG-PR and mMGE
    selected: dict[str, dict[str, tuple[bool, str]]] = {s: {} for s in wanted_sources}
    counts = {s: 0 for s in wanted_sources}
    with open(pl2 / "records.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            src = parts[idx["source"]]
            if src not in wanted_sources:
                continue
            if parts[idx["first_source"]] or parts[idx["owner_stage"]]:
                continue  # skip repeats and earlier-stage owned contigs
            c = counts[src]
            counts[src] = c + 1
            if c % every != 0:
                continue
            acc = parts[idx["accession"]]
            circ = parts[idx["circular"]] == "true"
            if dtr_circular and src == "plasmidscope_IMG-PR" and parts[idx["topology_value"]] == "direct terminal repeat":
                circ = True
            qual = parts[idx["assembly_quality"]]
            selected[src][acc] = (circ, qual)

    # Stream each FASTA and write the selected slice FASTA + genomes list
    genomes = []
    objs = []
    total_bp = 0
    total_circ = 0
    by_src = {}
    for src in ("plasmidscope_IMG-PR", "plasmidscope_mMGE"):
        want = selected[src]
        src_fa = pl2 / f"{src}.fasta"
        out_fa = slice_dir / f"{src}.slice.fasta"
        n_kept = 0
        bp_kept = 0
        circ_kept = 0
        with open(out_fa, "wb") as out_fh:
            for rec in fasta.iter_records(src_fa):
                acc = rec.header.split()[0]
                meta = want.get(acc)
                if meta is None:
                    continue
                circ, qual = meta
                cu = ids.contig_uid(rec.seq)
                out_fh.write(b"\n".join(fasta.format_record(rec.header, rec.seq, 80)) + b"\n")
                genomes.append({
                    "contig_uids": [cu.hex()],
                    "genetic_code": 11,
                    "kind": "plasmid",
                    "circular": [cu.hex()] if circ else [],
                    "assembly_quality": qual,
                })
                n_kept += 1
                bp_kept += len(rec.seq)
                circ_kept += int(circ)
        objs.append({
            "path": str(out_fa),
            "source_name": src,
            "source_release": "b63_slice",
            "license": "local-only",
            "redistributable": False,
        })
        by_src[src] = {"genomes": n_kept, "bp": bp_kept, "circular": circ_kept}
        total_bp += bp_kept
        total_circ += circ_kept

    g_path = slice_dir / "genomes.json"
    o_path = slice_dir / "objects.json"
    g_path.write_text(json.dumps(genomes))
    o_path.write_text(json.dumps(objs, indent=2))
    info = {
        "every": every,
        "dtr_circular": dtr_circular,
        "genomes": len(genomes),
        "total_bp": total_bp,
        "circular": total_circ,
        "by_source": by_src,
        "genomes_json": str(g_path),
        "objects_json": str(o_path),
    }
    (slice_dir / "slice_info.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


def _inspect_batches(root: Path, stage: str) -> dict:
    """Count how `_batches` partitions the ingested stage (long vs short, short batch sizes)."""
    from biobox.l0seq.inputs import stage_genomes
    from biobox.l0seq.lookup import KeyLookup
    from biobox.l0seq.store import SeqStore
    from biobox.l1 import genecall

    store = SeqStore(root, [stage])
    try:
        genomes, _ = stage_genomes(root, stage, None)
        lengths = KeyLookup(store.index_files, "contig_uid", ["contig_len", "is_default_variant"])
        n_long = 0
        short_batch_lens: list[int] = []
        short_batch_bases: list[int] = []
        for is_long, batch in genecall._batches(
            genomes,
            lambda g: genecall._genome_bases(store, lengths, g),
            genecall.BATCH_BASES,
            genecall.THREAD_MIN_BASES,
        ):
            if is_long:
                n_long += 1
            else:
                short_batch_lens.append(len(batch))
                short_batch_bases.append(sum(genecall._genome_bases(store, lengths, g) for g in batch))
        return {
            "n_long_tasks": n_long,
            "n_short_batches": len(short_batch_lens),
            "short_genomes_total": sum(short_batch_lens),
            "short_batch_median_genomes": statistics.median(short_batch_lens) if short_batch_lens else 0,
            "short_batch_mean_genomes": round(statistics.mean(short_batch_lens), 2) if short_batch_lens else 0,
            "short_batch_mean_bp": round(statistics.mean(short_batch_bases), 1) if short_batch_bases else 0,
            "short_batch_singletons": sum(1 for x in short_batch_lens if x == 1),
        }
    finally:
        store.close()


def run_arm_b(args: argparse.Namespace) -> None:
    from biobox.l0seq import SourceObject, ingest
    from biobox.l0seq.inputs import GenomeInput, write_genome_input
    from biobox.l1 import genecall

    work = Path(args.work_dir)
    tag = args.tag
    base = work / f"arm_b_{tag}"
    slice_dir = work / f"slice_every{args.every}{'_dtr' if args.dtr_circular else ''}"
    root = slice_dir / "root"
    stage = "PL3S"

    if not (slice_dir / "slice_info.json").is_file() or not root.exists():
        print(f"Building 1-in-{args.every} PL3 slice in {slice_dir}...", flush=True)
        info = _build_arm_b_slice(Path(args.pl2_dir), slice_dir, args.every, args.dtr_circular)
        print(
            f"Slice built: {info['genomes']:,} genomes, {info['total_bp']/1e6:.2f} Mbp, "
            f"{info['circular']:,} circular",
            flush=True,
        )
        raw_g = json.loads(Path(info["genomes_json"]).read_text())
        g_objs = [
            GenomeInput(
                tuple(bytes.fromhex(u) for u in g["contig_uids"]),
                g["genetic_code"],
                kind=g["kind"],
                circular=tuple(bytes.fromhex(u) for u in g["circular"]),
                assembly_quality=g["assembly_quality"],
            )
            for g in raw_g
        ]
        raw_o = json.loads(Path(info["objects_json"]).read_text())
        s_objs = [
            SourceObject(Path(o["path"]), o["source_name"], o["source_release"], o["license"], o["redistributable"])
            for o in raw_o
        ]
        if root.exists():
            shutil.rmtree(root)
        write_genome_input(root, stage, g_objs)
        ingest(s_objs, root, stage)
        print(f"Ingested slice into {root}", flush=True)
    else:
        info = json.loads((slice_dir / "slice_info.json").read_text())
        print(f"Reusing existing ingested slice at {root} ({info['genomes']:,} genomes, {info['total_bp']/1e6:.2f} Mbp)", flush=True)

    # Plan estimate from genecall.plan
    plan_tmp = base / "plan_tmp"
    plan_tmp.mkdir(parents=True, exist_ok=True)
    p = genecall.plan(root, stage, batch_cpu_min=240.0, tmp=plan_tmp)
    batch_stats = _inspect_batches(root, stage)
    print(
        f"Plan estimate: {p['estimate_cpu_min']:.3f} CPU-min "
        f"(bases={p['bases_called']:,}, genomes={p['genomes']:,}); "
        f"_batches: {batch_stats}",
        flush=True,
    )

    runs = []
    first_hashes = None
    for rep in range(1, args.reps + 1):
        out_dir = base / f"out_rep{rep}"
        tmp_dir = base / f"tmp_rep{rep}"
        for d in (out_dir, tmp_dir):
            if d.exists():
                shutil.rmtree(d)
            d.mkdir(parents=True)

        # Run `python -m biobox.cli stage genecall` in a child process and measure its + its children's CPU via rusage
        # Wrap with a tiny python snippet that runs cli.main and prints RUSAGE_SELF + RUSAGE_CHILDREN so we capture
        # exact CPU-minutes of this subprocess AND its worker processes independently of earlier reps.
        wrapper = (
            "import json, resource, sys, time; "
            "from biobox import cli; "
            "t0 = time.monotonic(); "
            "rc = cli.main(sys.argv[1:]); "
            "me = resource.getrusage(resource.RUSAGE_SELF); "
            "kids = resource.getrusage(resource.RUSAGE_CHILDREN); "
            "u = {'rc': rc, 'wall_s': round(time.monotonic() - t0, 3), "
            "'self_cpu_s': round(me.ru_utime + me.ru_stime, 4), "
            "'kids_cpu_s': round(kids.ru_utime + kids.ru_stime, 4), "
            "'total_cpu_min': round((me.ru_utime + me.ru_stime + kids.ru_utime + kids.ru_stime) / 60.0, 4), "
            "'peak_rss_mib': round(max(me.ru_maxrss, kids.ru_maxrss) / 1024.0, 1)}; "
            "print('B63_USAGE=' + json.dumps(u), file=sys.stderr, flush=True); "
            "sys.exit(rc)"
        )
        cmd = [
            sys.executable, "-c", wrapper,
            "stage", "genecall",
            "--root", str(root),
            "--stage", stage,
            "--out", str(out_dir),
            "--tmp", str(tmp_dir),
            "--workers", str(args.workers),
        ]
        load_before = _load_avg()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        load_after = _load_avg()
        if proc.returncode != 0:
            raise RuntimeError(f"stage genecall failed (rc={proc.returncode}):\n{proc.stderr[-2000:]}")
        usage = None
        for line in proc.stderr.splitlines():
            if line.startswith("B63_USAGE="):
                usage = json.loads(line[len("B63_USAGE="):])
        if usage is None:
            raise RuntimeError(f"Missing B63_USAGE in stderr:\n{proc.stderr[-1000:]}")
        stdout_s = proc.stdout.strip()
        try:
            cli_res = json.loads(stdout_s)["result"]
        except json.JSONDecodeError:
            # Handle any leading non-JSON lines before the JSON object
            idx = stdout_s.find("{")
            if idx < 0:
                raise RuntimeError(f"No JSON object in stdout: {stdout_s[:400]!r}")
            cli_res = json.loads(stdout_s[idx:])["result"]
        hashes = _tree_sha256(out_dir)
        if first_hashes is None:
            first_hashes = hashes
        elif hashes != first_hashes:
            raise RuntimeError(f"Rep {rep} output files differ from rep 1!")

        ratio = usage["total_cpu_min"] / p["estimate_cpu_min"]
        step_s_per_mbp = usage["total_cpu_min"] * 60.0 / (info["total_bp"] / 1e6)
        run_rec = {
            "rep": rep,
            "workers": args.workers,
            "estimate_cpu_min": p["estimate_cpu_min"],
            "measured_cpu_min": usage["total_cpu_min"],
            "ratio": round(ratio, 4),
            "step_s_per_mbp": round(step_s_per_mbp, 4),
            "self_cpu_s": usage["self_cpu_s"],
            "kids_cpu_s": usage["kids_cpu_s"],
            "wall_s": usage["wall_s"],
            "peak_rss_mib": usage["peak_rss_mib"],
            "genes": cli_res["genes"],
            "distinct_proteins": cli_res["distinct_proteins"],
            "wrapping_genes": cli_res["wrapping_genes"],
            "load_before": load_before,
            "load_after": load_after,
        }
        runs.append(run_rec)
        print(
            f"Rep {rep}: measured={usage['total_cpu_min']:.3f} CPU-min vs est={p['estimate_cpu_min']:.3f} "
            f"(ratio={ratio:.4f}, {step_s_per_mbp:.3f} s/Mbp, wall={usage['wall_s']:.1f}s, "
            f"self={usage['self_cpu_s']:.1f}s, kids={usage['kids_cpu_s']:.1f}s, rss={usage['peak_rss_mib']:.1f} MiB)",
            flush=True,
        )
        shutil.rmtree(tmp_dir, ignore_errors=True)
        if rep > 1:
            shutil.rmtree(out_dir, ignore_errors=True)

    result = {
        "tag": tag,
        "slice": info,
        "plan_estimate_cpu_min": p["estimate_cpu_min"],
        "plan_rates": p["rates"],
        "batch_stats": batch_stats,
        "workers": args.workers,
        "runs": runs,
        "median_cpu_min": round(statistics.median(r["measured_cpu_min"] for r in runs), 4),
        "median_ratio": round(statistics.median(r["ratio"] for r in runs), 4),
        "output_hashes": first_hashes,
    }
    out_json = work / f"arm_b_{tag}.json"
    out_json.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {out_json}", flush=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card B-63: meta-mode gene-calling cost and batch sizing arms.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("arm-a", help="Measure GeneFinder/ViralGeneFinder(meta=True) instantiation vs reuse overhead.")
    a.add_argument("--pl2-dir", default=str(DEFAULT_PL2), help="Read-only PL2 inputs directory.")
    a.add_argument("--work-dir", default=str(DEFAULT_WORK), help="Scratch output directory.")
    a.add_argument("--sources", default="IMG-PR,mMGE,gut,RefSeq", help="Comma-separated source names.")
    a.add_argument("--modes", default="meta,gv-meta", help="Comma-separated modes to test (meta,gv-meta).")
    a.add_argument("--bin-bp", type=int, default=3_000_000, help="Target bp per length bin per source.")
    a.add_argument("--bin-max-seqs", type=int, default=400, help="Max sequences per length bin per source.")
    a.add_argument("--unstrat-every", type=int, default=50, help="Stride for unstratified sample.")
    a.add_argument("--unstrat-bp", type=int, default=8_000_000, help="Target bp for unstratified sample per source.")
    a.add_argument("--reps", type=int, default=3, help="Number of interleaved repetitions.")
    a.add_argument("--init-iters", type=int, default=500, help="Iterations for standalone init timing.")

    b = sub.add_parser("arm-b", help="Run full `stage genecall` on a representative IMG-PR + mMGE scratch slice.")
    b.add_argument("--pl2-dir", default=str(DEFAULT_PL2), help="Read-only PL2 inputs directory.")
    b.add_argument("--work-dir", default=str(DEFAULT_WORK), help="Scratch output directory.")
    b.add_argument("--every", type=int, default=200, help="Take 1 in N new IMG-PR/mMGE records.")
    b.add_argument("--dtr-circular", action="store_true", help="Mark IMG-PR DTR records circular (DEC-63).")
    b.add_argument("--workers", type=int, default=4, help="Worker count for stage genecall.")
    b.add_argument("--reps", type=int, default=2, help="Number of repetitions.")
    b.add_argument("--tag", default="baseline", help="Label for output JSON and directories.")

    args = ap.parse_args(argv)
    if args.cmd == "arm-a":
        run_arm_a(args)
    elif args.cmd == "arm-b":
        run_arm_b(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

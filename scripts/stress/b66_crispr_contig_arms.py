"""Card B-66: pod scratch measurement arms for O(1)-stage `SeqStore.get` and per-stored-contig CRISPR CPU-minutes.

Read-only on `/tmp/if2453/r10/root` and `/tmp/if2453/s2/inputs/SPX1`; all scratch under `--scratch` (e.g.
`/tmp/if2453/b66/`).

Modes:
  summary   Inspect all stages in `--root` (stored contigs, stored bases, circular flank bases, CRISPR run records).
  seqstore  Benchmark before (linear-stage + double locate) vs after (find_first + single locate) on 20,000 G1 contigs
            with all 12 stages opened.
  g1        Build a deterministic 1-in-20 contig scratch view of stage G1 (~61,360 contigs, ~502 Mbp) and measure
            end-to-end `crispr.run` and `crispr.run_near_cas` (`cpus=4`, all 12 stages opened).
  spx1      Ingest a deterministic 1-in-40 assembly sample of SPX1 (`--min-contig-len 5000`) into scratch and measure
            end-to-end `crispr.run` (`cpus=4`, with the 12 existing stages also opened).
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox import fasta, ids, io
from biobox import stage as bs
from biobox.l0seq import SeqStore, SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l0seq.store import INDEX_DIR, LOCATOR_DIR, stages_on_disk
from biobox.l1 import crispr
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.registry.schemas import table


def _stage_contigs_local(root: Path, stage: str) -> int:
    rec = json.loads((root / "runs" / f"stage={stage}" / "ingest.json").read_text())
    prefix = f"{INDEX_DIR}/stage={stage}/"
    total = 0
    for o in rec["outputs"]:
        if not o["path"].startswith(prefix):
            continue
        pf = pq.ParquetFile(root / o["path"])
        for i in range(pf.num_row_groups):
            t = pf.read_row_group(i, columns=["is_default_variant"])
            total += pc.sum(pc.cast(t.column("is_default_variant"), pa.int64())).as_py() or 0
    return int(total)


def cmd_summary(args) -> int:
    root = Path(args.root)
    stages = stages_on_disk(root)
    def_step = f"crispr-{crispr.producer()}"
    nc_step = f"crispr-{crispr.near_cas_producer('38.2')}"
    rows = []
    for st in stages:
        ing = root / "runs" / f"stage={st}" / "ingest.json"
        if not ing.is_file():
            continue
        bases = bs.stage_bases(root, st)
        flank = bs.crispr_flank_bases(root, st)
        contigs = _stage_contigs_local(root, st)
        drec = root / "runs" / f"stage={st}" / f"{def_step}.json"
        nrec = root / "runs" / f"stage={st}" / f"{nc_step}.json"
        ddoc = json.loads(drec.read_text()) if drec.is_file() else None
        ndoc = json.loads(nrec.read_text()) if nrec.is_file() else None
        rows.append({
            "stage": st,
            "contigs_stored": contigs,
            "bases_stored": bases,
            "flank_bases": flank,
            "scanned_bases": bases + flank,
            "default_cpu_min": None if ddoc is None else ddoc.get("cpu_min"),
            "default_wall_s": None if ddoc is None else ddoc.get("wall_s"),
            "near_cas_cpu_min": None if ndoc is None else ndoc.get("cpu_min"),
            "near_cas_wall_s": None if ndoc is None else ndoc.get("wall_s"),
        })
    print(json.dumps({"stages": stages, "rows": rows}, indent=2))
    return 0


def cmd_seqstore(args) -> int:
    root = Path(args.root)
    stages = ["G1", *(s for s in stages_on_disk(root) if s != "G1")]
    # Collect first N default contigs from G1's contig_index
    n_want = args.n
    cus: list[bytes] = []
    for p in io.part_files(root / INDEX_DIR / "stage=G1"):
        pf = pq.ParquetFile(p)
        for i in range(pf.num_row_groups):
            t = pf.read_row_group(i, columns=["contig_uid", "is_default_variant"])
            for cu, d in zip(t.column("contig_uid").to_pylist(), t.column("is_default_variant").to_pylist()):
                if d:
                    cus.append(cu)
                    if len(cus) >= n_want:
                        break
            if len(cus) >= n_want:
                break
        if len(cus) >= n_want:
            break

    # 1. Pre-B-66 behavior (linear scan via .get() and double locator lookup: length_of + raw)
    s_old = SeqStore(root, stages)
    t0 = time.perf_counter()
    c0 = time.process_time()
    old_bases = 0
    for cu in cus:
        rd = next(r["raw_digest"] for r in s_old.index.get(cu) if r["is_default_variant"])
        row1 = s_old.locator.get(rd)[0]
        _len0 = row1["contig_len"]
        row2 = s_old.locator.get(rd)[0]
        rel, _len1 = row2["shard"], row2["contig_len"]
        old_bases += len(s_old._shard(rel).fetch(rd.hex(), 0, _len1))
    old_wall = time.perf_counter() - t0
    old_cpu = time.process_time() - c0
    old_idx_rg, old_loc_rg = s_old.index.row_groups_read, s_old.locator.row_groups_read
    s_old.close()

    # 2. Post-B-66 behavior (SeqStore.get with find_first + single _where)
    s_new = SeqStore(root, stages)
    t0 = time.perf_counter()
    c0 = time.process_time()
    new_bases = 0
    for cu in cus:
        new_bases += len(s_new.get(cu))
    new_wall = time.perf_counter() - t0
    new_cpu = time.process_time() - c0
    new_idx_rg, new_loc_rg = s_new.index.row_groups_read, s_new.locator.row_groups_read
    s_new.close()

    assert old_bases == new_bases
    print(json.dumps({
        "stages": stages,
        "n_contigs": len(cus),
        "bases": new_bases,
        "before": {
            "wall_s": round(old_wall, 4),
            "cpu_s": round(old_cpu, 4),
            "ms_per_contig": round(old_cpu * 1000.0 / len(cus), 4),
            "index_rg_read": old_idx_rg,
            "locator_rg_read": old_loc_rg,
        },
        "after": {
            "wall_s": round(new_wall, 4),
            "cpu_s": round(new_cpu, 4),
            "ms_per_contig": round(new_cpu * 1000.0 / len(cus), 4),
            "index_rg_read": new_idx_rg,
            "locator_rg_read": new_loc_rg,
        },
    }, indent=2))
    return 0


def _timed_subproc(code: str, argv: list[str]) -> dict:
    """Run `code` in a fresh Python child process and return wall_s, total CPU-s (self+children), and parsed JSON."""
    wrapper = (
        "import json, resource, sys, time\n"
        "t0 = time.perf_counter()\n"
        f"{code}\n"
        "wall = time.perf_counter() - t0\n"
        "rs = resource.getrusage(resource.RUSAGE_SELF)\n"
        "rc = resource.getrusage(resource.RUSAGE_CHILDREN)\n"
        "cpu_s = rs.ru_utime + rs.ru_stime + rc.ru_utime + rc.ru_stime\n"
        "print('__B66_RES__' + json.dumps({'wall_s': wall, 'cpu_s': cpu_s, 'cpu_min': cpu_s / 60.0, 'res': res}))\n"
    )
    proc = subprocess.run([sys.executable, "-c", wrapper, *argv], capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"child failed ({proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    for line in reversed(proc.stdout.splitlines()):
        if line.startswith("__B66_RES__"):
            return json.loads(line[len("__B66_RES__"):])
    raise RuntimeError(f"missing __B66_RES__ in stdout:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


def _symlink_l0seq(src_root: Path, dst_root: Path, stages: list[str]) -> None:
    (dst_root / "l0seq").mkdir(parents=True, exist_ok=True)
    sp = src_root / "l0seq" / "store_params.json"
    if sp.exists() and not (dst_root / "l0seq" / "store_params.json").exists():
        os.symlink(sp, dst_root / "l0seq" / "store_params.json")
    for sub in ("seq_shards", "seq_locator", "contig_index"):
        (dst_root / "l0seq" / sub).mkdir(parents=True, exist_ok=True)
        for st in stages:
            src = src_root / "l0seq" / sub / f"stage={st}"
            dst = dst_root / "l0seq" / sub / f"stage={st}"
            if src.is_dir() and not dst.exists():
                os.symlink(src, dst)


def cmd_g1(args) -> int:
    root = Path(args.root)
    scratch = Path(args.scratch) / "g1_arm"
    if scratch.exists():
        shutil.rmtree(scratch)
    view = scratch / "view"
    work = scratch / "work"
    out_def = scratch / "out_def"
    out_nc = scratch / "out_nc"
    work.mkdir(parents=True)

    all_stages = stages_on_disk(root)
    other_stages = [s for s in all_stages if s != "G1"]
    _symlink_l0seq(root, view, all_stages)

    # Deterministic 1-in-stride sample of G1's default-variant contigs
    stride = args.stride
    sampled_cus: list[bytes] = []
    sampled_bases = 0
    idx = 0
    for p in io.part_files(root / INDEX_DIR / "stage=G1"):
        pf = pq.ParquetFile(p)
        for rg in range(pf.num_row_groups):
            t = pf.read_row_group(rg, columns=["contig_uid", "contig_len", "is_default_variant"])
            for cu, clen, default in zip(
                t.column("contig_uid").to_pylist(),
                t.column("contig_len").to_pylist(),
                t.column("is_default_variant").to_pylist(),
            ):
                if default:
                    if idx % stride == 0:
                        sampled_cus.append(cu)
                        sampled_bases += int(clen)
                    idx += 1

    # Write sampled genome_input for G1 (each contig as a linear prokaryote genome, matching G1's linear contigs)
    genomes = [GenomeInput((cu,), 11, kind="prokaryote", circular=()) for cu in sampled_cus]
    write_genome_input(view, "G1", genomes)

    # Arm 1: default crispr.run
    def_code = (
        "from pathlib import Path\n"
        "from biobox.l1 import crispr\n"
        "res = crispr.run(Path(sys.argv[1]), 'G1', int(sys.argv[2]), tmp=Path(sys.argv[3]), "
        "out=Path(sys.argv[4]), seq_stages=sys.argv[5].split(','))\n"
    )
    def_m = _timed_subproc(def_code, [str(view), str(args.cpus), str(work), str(out_def), ",".join(other_stages)])

    # Wire up inputs for crispr.run_near_cas: default feature, G1 feature_protein, and empty/synthetic l1_domains
    prod = crispr.producer()
    feat_dst = view / "l1" / "feature" / f"producer={prod}"
    feat_dst.mkdir(parents=True, exist_ok=True)
    os.symlink(out_def / "l1" / "feature" / f"producer={prod}" / "stage=G1", feat_dst / "stage=G1")

    fp_src = root / "l1" / "feature_protein" / f"gene_call_version={GENE_CALL_VERSION}" / "stage=G1"
    fp_dst = view / "l1" / "feature_protein" / f"gene_call_version={GENE_CALL_VERSION}"
    fp_dst.mkdir(parents=True, exist_ok=True)
    os.symlink(fp_src, fp_dst / "stage=G1")

    dom_src = root / "l1" / "domains" / "library=pfam-38.2" / "stage=G1"
    dom_dir = view / "l1" / "domains" / "library=pfam-38.2" / "stage=G1"
    if dom_src.is_dir():
        dom_dir.parent.mkdir(parents=True, exist_ok=True)
        os.symlink(dom_src, dom_dir)
    else:
        # Write a valid empty/minimal l1_domains table with the schema of l1_domains
        dom_spec = table("l1_domains")
        empty_tbl = pa.table({f.name: pa.array([], type=f.type) for f in dom_spec.schema}, schema=dom_spec.schema)
        io.write_parts(empty_tbl, dom_spec, dom_dir, 64, 64 * 1024 * 1024)

    nc_code = (
        "from pathlib import Path\n"
        "from biobox.l1 import crispr\n"
        "res = crispr.run_near_cas(Path(sys.argv[1]), 'G1', int(sys.argv[2]), library_version='38.2', "
        "tmp=Path(sys.argv[3]), out=Path(sys.argv[4]), seq_stages=sys.argv[5].split(','))\n"
    )
    nc_m = _timed_subproc(nc_code, [str(view), str(args.cpus), str(work), str(out_nc), ",".join(other_stages)])

    print(json.dumps({
        "arm": "G1",
        "stride": stride,
        "cpus": args.cpus,
        "total_g1_contigs": idx,
        "sampled_contigs": len(sampled_cus),
        "sampled_bases": sampled_bases,
        "stages_opened": len(all_stages),
        "default": {
            "wall_s": round(def_m["wall_s"], 3),
            "cpu_s": round(def_m["cpu_s"], 3),
            "cpu_min": round(def_m["cpu_min"], 4),
            "arrays": def_m["res"]["arrays"],
            "contigs": def_m["res"]["contigs"],
            "bases": def_m["res"]["bases"],
        },
        "near_cas": {
            "wall_s": round(nc_m["wall_s"], 3),
            "cpu_s": round(nc_m["cpu_s"], 3),
            "cpu_min": round(nc_m["cpu_min"], 4),
            "kept": nc_m["res"]["kept"],
            "candidates": nc_m["res"].get("candidates"),
            "contigs": nc_m["res"]["contigs"],
            "bases": nc_m["res"]["bases"],
        },
    }, indent=2))
    return 0


def cmd_spx1(args) -> int:
    root = Path(args.root)
    spx_in = Path(args.spx_inputs)
    scratch = Path(args.scratch) / "spx1_arm"
    if scratch.exists():
        shutil.rmtree(scratch)
    sroot = scratch / "root"
    work = scratch / "work"
    out_def = scratch / "out_def"
    work.mkdir(parents=True)

    all_existing = stages_on_disk(root)
    _symlink_l0seq(root, sroot, all_existing)

    raw_objs = json.loads((spx_in / "objects.json").read_text())
    stride = args.stride
    picked = raw_objs[::stride]
    if args.max_objects:
        picked = picked[:args.max_objects]

    source_objs = []
    kept_cus: list[bytes] = []
    kept_bases = 0
    seen_cus: set[bytes] = set()
    for o in picked:
        p = Path(o["path"])
        if not p.is_absolute():
            p = spx_in / p
        source_objs.append(SourceObject(
            path=p,
            source_name=o["source_name"],
            source_release=str(o["source_release"]),
            license=o.get("license", "unknown"),
            redistributable=bool(o.get("redistributable", False)),
        ))
        for rec in fasta.iter_records(p):
            if len(rec.seq) >= 5000:
                cu = ids.contig_uid(rec.seq)
                if cu not in seen_cus:
                    seen_cus.add(cu)
                    kept_cus.append(cu)
                    kept_bases += len(rec.seq)

    genomes = [GenomeInput((cu,), 11, kind="prokaryote", circular=()) for cu in sorted(kept_cus)]
    write_genome_input(sroot, "SPX1S", genomes)
    t_ing0 = time.perf_counter()
    ing_rep = ingest(source_objs, sroot, "SPX1S", tmp=work, min_contig_len=5000)
    t_ing = time.perf_counter() - t_ing0

    # Count actually stored default contigs in SPX1S
    n_stored = 0
    b_stored = 0
    for p in io.part_files(sroot / INDEX_DIR / "stage=SPX1S"):
        pf = pq.ParquetFile(p)
        for rg in range(pf.num_row_groups):
            t = pf.read_row_group(rg, columns=["contig_len", "is_default_variant"])
            for clen, d in zip(t.column("contig_len").to_pylist(), t.column("is_default_variant").to_pylist()):
                if d:
                    n_stored += 1
                    b_stored += int(clen)

    def_code = (
        "from pathlib import Path\n"
        "from biobox.l1 import crispr\n"
        "res = crispr.run(Path(sys.argv[1]), 'SPX1S', int(sys.argv[2]), tmp=Path(sys.argv[3]), "
        "out=Path(sys.argv[4]), seq_stages=sys.argv[5].split(','))\n"
    )
    def_m = _timed_subproc(def_code, [str(sroot), str(args.cpus), str(work), str(out_def), ",".join(all_existing)])

    print(json.dumps({
        "arm": "SPX1",
        "stride": stride,
        "n_objects_total": len(raw_objs),
        "n_objects_picked": len(picked),
        "cpus": args.cpus,
        "stages_opened": len(all_existing) + 1,
        "ingest_wall_s": round(t_ing, 2),
        "contigs_stored": n_stored,
        "bases_stored": b_stored,
        "default": {
            "wall_s": round(def_m["wall_s"], 3),
            "cpu_s": round(def_m["cpu_s"], 3),
            "cpu_min": round(def_m["cpu_min"], 4),
            "arrays": def_m["res"]["arrays"],
            "contigs": def_m["res"]["contigs"],
            "bases": def_m["res"]["bases"],
        },
    }, indent=2))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("summary")
    p.add_argument("--root", default="/tmp/if2453/r10/root")
    p = sub.add_parser("seqstore")
    p.add_argument("--root", default="/tmp/if2453/r10/root")
    p.add_argument("--n", type=int, default=20_000)
    p = sub.add_parser("g1")
    p.add_argument("--root", default="/tmp/if2453/r10/root")
    p.add_argument("--scratch", default="/tmp/if2453/b66")
    p.add_argument("--stride", type=int, default=20)
    p.add_argument("--cpus", type=int, default=4)
    p = sub.add_parser("spx1")
    p.add_argument("--root", default="/tmp/if2453/r10/root")
    p.add_argument("--spx-inputs", default="/tmp/if2453/s2/inputs/SPX1")
    p.add_argument("--scratch", default="/tmp/if2453/b66")
    p.add_argument("--stride", type=int, default=40)
    p.add_argument("--max-objects", type=int, default=8)
    p.add_argument("--cpus", type=int, default=4)
    a = ap.parse_args()
    return {"summary": cmd_summary, "seqstore": cmd_seqstore, "g1": cmd_g1, "spx1": cmd_spx1}[a.mode](a)


if __name__ == "__main__":
    sys.exit(main())

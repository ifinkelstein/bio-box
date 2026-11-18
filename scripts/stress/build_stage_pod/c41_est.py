"""Card C-41: overlap, step estimates, sample measurements (genecall, NR90 linclust, CRISPR) and full-set cost.

    python c41_est.py --root /tmp/if2453/r10/root --inputs /tmp/if2453/s2/inputs/G1 \
        --lock /tmp/if2453/s3/globdb/genomes/G1/lock.json \
        --selection /tmp/if2453/s3/globdb/genomes/G1/selection.tsv \
        --table-summary /tmp/if2453/s3/globdb/tables/table_summary.json \
        --profile /tmp/if2453/ci/C-41/profiles/wilkcomp02.yaml \
        --pl2-inputs /tmp/if2453/s2/inputs/PL2 --v1-inputs /tmp/if2453/s2/inputs/V1 \
        --every 20 --crispr-every 100 --out /tmp/if2453/s3/globdb/est/G1
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import resource
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import fasta, ids, io, stage
from biobox.l0seq import store
from biobox.l0seq.inputs import GenomeInput
from biobox.l1 import crispr, genecall
from biobox.profile import load

csv.field_size_limit(1 << 30)


def cpu_self() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def cpu_children() -> float:
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def stage_disk_bytes(root: Path, st: str) -> dict[str, int]:
    """Bytes on disk in `root` under `stage=<st>` partitions, grouped by top-level directory."""
    by_top: dict[str, int] = defaultdict(int)
    for p in root.rglob(f"stage={st}"):
        if not p.is_dir():
            continue
        rel = p.relative_to(root).parts[0]
        for f in p.rglob("*"):
            if f.is_file():
                by_top[rel] += f.stat().st_size
    return dict(sorted(by_top.items()))


def p2_measured(root: Path) -> dict:
    runs = root / "runs" / "stage=P2"
    by = defaultdict(float)
    for f in sorted(runs.glob("*.json")):
        if f.name.endswith(".publishing.json"):
            continue
        d = json.loads(f.read_text())
        step = d.get("step", f.stem)
        key = (
            "genecall" if step.startswith("genecall")
            else "domains-pfam-38.2" if step.startswith("domains")
            else "crispr" if step.startswith("crispr")
            else step
        )
        by[key] += float(d.get("cpu_min") or 0.0)
    gcv = genecall.GENE_CALL_VERSION
    fp = io.part_files(root / "l1/feature_protein" / f"gene_call_version={gcv}" / "stage=P2")
    ps = io.part_files(root / "l0seq/protein_seq" / f"gene_call_version={gcv}" / "stage=P2",
                       "protein_prefix=*/part-*.parquet")
    ci = io.part_files(root / store.INDEX_DIR / "stage=P2")
    con = duckdb.connect(config={"threads": 2, "memory_limit": "4GB"})
    genes = sum(pq.ParquetFile(f).metadata.num_rows for f in fp)
    (nd, res), = con.execute(
        f"SELECT count(*), sum(length) FROM read_parquet({sql_list(ps)}, hive_partitioning = false)"
    ).fetchall()
    (nc, nb), = con.execute(
        f"SELECT count(*), sum(contig_len) FROM read_parquet({sql_list(ci)}, hive_partitioning = false) "
        "WHERE is_default_variant"
    ).fetchall()
    con.close()
    disk = stage_disk_bytes(root, "P2")
    return {
        "cpu_min_by_step": {k: round(v, 3) for k, v in sorted(by.items())},
        "genes": genes,
        "distinct_proteins": int(nd),
        "residues": int(res),
        "stored_contigs": int(nc),
        "stored_bases": int(nb),
        "genes_per_mbp": round(genes / nb * 1e6, 2),
        "distinct_per_gene": round(nd / genes, 4),
        "residues_per_distinct": round(res / nd, 2),
        "disk_bytes_by_dir": disk,
        "disk_bytes_total": sum(disk.values()),
        "disk_gb_per_gbp": round((sum(disk.values()) / 1e9) / (nb / 1e9), 4),
    }


def known_proteins(root: Path, uids: set[bytes]) -> dict:
    """How many of `uids` are in any stage's protein_seq on disk, overall and by stage."""
    gcv = genecall.GENE_CALL_VERSION
    stages = store.stages_on_disk(root)
    con = duckdb.connect(config={"threads": 2, "memory_limit": "6GB"})
    con.execute("CREATE TABLE q (protein_uid BLOB)")
    con.executemany("INSERT INTO q VALUES (?)", [(u,) for u in sorted(uids)])
    all_files = []
    by_stage: dict[str, int] = {}
    for s in stages:
        fs = io.part_files(
            root / "l0seq/protein_seq" / f"gene_call_version={gcv}" / f"stage={s}",
            "protein_prefix=*/part-*.parquet",
        )
        if not fs:
            continue
        all_files.extend(fs)
        (cnt,), = con.execute(
            f"SELECT count(DISTINCT protein_uid) FROM read_parquet({sql_list(fs)}, hive_partitioning = false) "
            "WHERE protein_uid IN (SELECT protein_uid FROM q)"
        ).fetchall()
        if cnt:
            by_stage[s] = int(cnt)
    (tot,), = con.execute(
        f"SELECT count(DISTINCT protein_uid) FROM read_parquet({sql_list(all_files)}, hive_partitioning = false) "
        "WHERE protein_uid IN (SELECT protein_uid FROM q)"
    ).fetchall()
    con.close()
    return {"total_known": int(tot), "by_stage": by_stage}


def lookup_shared_labels(root: Path, shared_rows: list[dict], sel_by_acc: dict[str, dict]) -> list[dict]:
    """Look up accession/original_header in `contig_label` (or `contig_occurrence`) for contigs shared with stored stages."""
    by_stage: dict[str, list[bytes]] = defaultdict(list)
    for r in shared_rows:
        by_stage[r["owner_stage"]].append(bytes.fromhex(r["contig_uid"]))
    label_map: dict[tuple[str, str], tuple[str, str]] = {}
    con = duckdb.connect(config={"threads": 1, "memory_limit": "2GB"})
    for st, uids in sorted(by_stage.items()):
        lfiles = io.part_files(root / "l0seq/contig_label" / f"stage={st}")
        if not lfiles or not uids:
            continue
        con.execute("CREATE OR REPLACE TABLE q (contig_uid BLOB)")
        con.executemany("INSERT INTO q VALUES (?)", [(u,) for u in sorted(set(uids))])
        rows = con.execute(
            f"SELECT hex(contig_uid), accession, original_header "
            f"FROM read_parquet({sql_list(lfiles)}, hive_partitioning = false) "
            "WHERE contig_uid IN (SELECT contig_uid FROM q)"
        ).fetchall()
        for h, acc, hdr in rows:
            label_map[(st, h.lower())] = (acc, hdr[:140])
    con.close()
    details = []
    for r in shared_rows:
        st = r["owner_stage"]
        cu_hex = r["contig_uid"]
        g1_acc = r["accession"]
        s_acc, s_hdr = label_map.get((st, cu_hex), ("", ""))
        details.append({
            "g1_accession": g1_acc,
            "g1_record": r["record_accession"],
            "g1_resource": sel_by_acc[g1_acc]["resource"],
            "g1_taxonomy": sel_by_acc[g1_acc]["taxonomy"],
            "contig_uid": cu_hex,
            "length": int(r["length"]),
            "owner_stage": st,
            "owner_accession": s_acc,
            "owner_header": s_hdr,
        })
    return details


def unbuilt_stage_overlap(stage_uids: set[str], inputs_dir: Path | None) -> int | None:
    if inputs_dir is None or not (inputs_dir / "genomes.json").is_file():
        return None
    gs = json.loads((inputs_dir / "genomes.json").read_text())
    other = {u for g in gs for u in g["contig_uids"]}
    return len(stage_uids & other)


def run_sample_linclust(prot: dict[bytes, str], threads: int = 2) -> dict:
    """Run MMseqs2 linclust at NR90 parameters (--min-seq-id 0.9 -c 0.8 --cov-mode 0 --cluster-mode 2) on `prot`."""
    work = Path(tempfile.mkdtemp(prefix="c41_nr90_"))
    try:
        fa = work / "sample_proteins.faa"
        with open(fa, "w") as fh:
            for u in sorted(prot):
                fh.write(f">{u.hex()}\n{prot[u]}\n")
        out_prefix = work / "clu"
        tmp_dir = work / "tmp"
        tmp_dir.mkdir()
        t0 = cpu_children()
        subprocess.run(
            [
                "mmseqs", "easy-linclust", str(fa), str(out_prefix), str(tmp_dir),
                "--min-seq-id", "0.9", "-c", "0.8", "--cov-mode", "0",
                "--cluster-mode", "2", "--kmer-per-seq", "80",
                "--threads", str(threads), "-v", "0",
            ],
            check=True,
            stdin=subprocess.DEVNULL,
        )
        dt = cpu_children() - t0
        tsv_path = work / "clu_cluster.tsv"
        reps: set[bytes] = set()
        len_mismatch_gt_20pct = 0
        members_non_rep = 0
        with open(tsv_path) as fh:
            for line in fh:
                r_hex, m_hex = line.rstrip("\n").split("\t")
                ru = bytes.fromhex(r_hex)
                mu = bytes.fromhex(m_hex)
                reps.add(ru)
                if mu != ru:
                    members_non_rep += 1
                    lr = len(prot[ru])
                    lm = len(prot[mu])
                    if abs(lm - lr) > 0.2 * lr:
                        len_mismatch_gt_20pct += 1
        rep_res = sum(len(prot[u]) for u in reps)
        tot_res = sum(len(s) for s in prot.values())
        return {
            "input_distinct_proteins": len(prot),
            "input_residues": tot_res,
            "nr90_clusters": len(reps),
            "non_rep_members": members_non_rep,
            "members_length_diff_gt_20pct": len_mismatch_gt_20pct,
            "rep_protein_fraction": round(len(reps) / len(prot), 4) if prot else 0.0,
            "rep_residues": rep_res,
            "rep_residue_fraction": round(rep_res / tot_res, 4) if tot_res else 0.0,
            "linclust_cpu_s": round(dt, 3),
        }
    finally:
        shutil.rmtree(work, ignore_errors=True)


def run_sample_crispr(sample_fastas: list[Path]) -> dict:
    """Run MinCED (-minNR 3) in <=50 Mbp batches as `biobox.l1.crispr` does, timing child JVM CPU."""
    java, jar = crispr._tool()
    work = Path(tempfile.mkdtemp(prefix="c41_crispr_"))
    try:
        batches: list[tuple[Path, int, int]] = []
        cur_idx = 0
        cur_bases = 0
        cur_recs = 0
        cur_fh = None
        cur_path = work / f"batch_{cur_idx}.fna"
        for fp in sample_fastas:
            for rec in fasta.iter_records(fp):
                if cur_fh is None:
                    cur_path = work / f"batch_{cur_idx}.fna"
                    cur_fh = open(cur_path, "wb")
                cur_fh.write(b">" + rec.accession.encode() + b"\n" + rec.seq.upper() + b"\n")
                cur_bases += len(rec.seq)
                cur_recs += 1
                if cur_bases >= crispr.BATCH_BASES:
                    cur_fh.close()
                    batches.append((cur_path, cur_recs, cur_bases))
                    cur_fh = None
                    cur_idx += 1
                    cur_bases = 0
                    cur_recs = 0
        if cur_fh is not None:
            cur_fh.close()
            batches.append((cur_path, cur_recs, cur_bases))

        tot_cpu_s = 0.0
        tot_arrays = 0
        tot_bases = sum(b[2] for b in batches)
        tot_recs = sum(b[1] for b in batches)
        for i, (bpath, _, _) in enumerate(batches):
            txt = work / f"b{i}.txt"
            gff = work / f"b{i}.gff"
            sc = work / f"tmp_{i}"
            sc.mkdir()
            argv = crispr._argv(java, jar, bpath, txt, gff, sc, min_repeats=crispr.MIN_REPEATS)
            t0 = cpu_children()
            subprocess.run(argv, check=True, stdin=subprocess.DEVNULL, capture_output=True)
            tot_cpu_s += cpu_children() - t0
            tot_arrays += sum(1 for _ in crispr.parse_gff(gff.read_text().splitlines()))
        return {
            "genomes": len(sample_fastas),
            "contigs": tot_recs,
            "bases": tot_bases,
            "batches_50mbp": len(batches),
            "arrays_found": tot_arrays,
            "arrays_per_gbp": round(tot_arrays / (tot_bases / 1e9), 2) if tot_bases else 0.0,
            "cpu_s": round(tot_cpu_s, 3),
            "cpu_s_per_gbp": round(tot_cpu_s / (tot_bases / 1e9), 2) if tot_bases else 0.0,
            "code_constant_cpu_s_per_gbp": stage.CRISPR_CPU_S_PER_GBP,
        }
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--inputs", required=True, type=Path)
    ap.add_argument("--lock", required=True, type=Path)
    ap.add_argument("--selection", required=True, type=Path)
    ap.add_argument("--table-summary", required=True, type=Path)
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--pl2-inputs", type=Path, default=None)
    ap.add_argument("--v1-inputs", type=Path, default=None)
    ap.add_argument("--every", type=int, default=20)
    ap.add_argument("--crispr-every", type=int, default=100)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    prof = load(a.profile)
    lock = json.loads(a.lock.read_text())["genomes"]
    sel_by_acc: dict[str, dict] = {}
    with open(a.selection, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sel_by_acc[r["accession"]] = r
    genomes_json = json.loads((a.inputs / "genomes.json").read_text())
    objects_json = json.loads((a.inputs / "objects.json").read_text())
    tsum = json.loads(a.table_summary.read_text())

    # Map contig_uid -> length and intra-stage genome memberships by scanning the stage FASTAs once
    contig_len: dict[bytes, int] = {}
    contig_accs: dict[bytes, list[str]] = defaultdict(list)
    acc_uids: dict[str, list[bytes]] = {}
    for acc in sorted(lock):
        fp = Path(lock[acc]["file"])
        uids = []
        seen_in_g: set[bytes] = set()
        for rec in fasta.iter_records(fp):
            cu = ids.contig_uid(rec.seq)
            contig_len.setdefault(cu, len(rec.seq))
            if cu not in seen_in_g:
                seen_in_g.add(cu)
                uids.append(cu)
                contig_accs[cu].append(acc)
        acc_uids[acc] = uids

    # Read shared_records.tsv (written by make_inputs.py prok --root)
    shared_rows = []
    shared_uids_by_stage: dict[str, list[bytes]] = defaultdict(list)
    stored_owned_uids: dict[bytes, str] = {}
    with open(a.inputs / "shared_records.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            shared_rows.append(r)
            cu = bytes.fromhex(r["contig_uid"])
            st = r["owner_stage"]
            shared_uids_by_stage[st].append(cu)
            stored_owned_uids[cu] = st

    distinct_stage_uids = {bytes.fromhex(u) for g in genomes_json for u in g["contig_uids"]}
    distinct_stage_hex = {u.hex() for u in distinct_stage_uids}
    held_bases = sum(contig_len[cu] for cu in distinct_stage_uids if cu in stored_owned_uids)
    new_bases = sum(contig_len[cu] for cu in distinct_stage_uids if cu not in stored_owned_uids)
    total_stage_bases = held_bases + new_bases

    # Intra-stage shared contigs (contigs in >1 genome of G1)
    intra_shared = [
        {"contig_uid": cu.hex(), "length": contig_len[cu], "genomes": accs,
         "resources": [sel_by_acc[x]["resource"] for x in accs]}
        for cu, accs in sorted(contig_accs.items(), key=lambda kv: (-contig_len[kv[0]], kv[0]))
        if len(accs) > 1
    ]

    overlap_out = {
        "stored_stages_on_disk": store.stages_on_disk(a.root),
        "shared_occurrences_with_stored_stages": len(shared_rows),
        "shared_distinct_contigs_with_stored_stages": len(stored_owned_uids),
        "shared_bases_with_stored_stages": held_bases,
        "shared_contigs_by_owner_stage": {
            st: len(set(us)) for st, us in sorted(shared_uids_by_stage.items())
        },
        "shared_bases_by_owner_stage": {
            st: sum(contig_len[cu] for cu in set(us)) for st, us in sorted(shared_uids_by_stage.items())
        },
        "shared_contig_labels_sample": lookup_shared_labels(a.root, shared_rows, sel_by_acc),
        "intra_stage_shared_contigs": len(intra_shared),
        "intra_stage_shared_bases": sum(x["length"] for x in intra_shared),
        "intra_stage_shared_examples": intra_shared[:25],
        "shared_with_PL2_inputs": unbuilt_stage_overlap(distinct_stage_hex, a.pl2_inputs),
        "shared_with_V1_inputs": unbuilt_stage_overlap(distinct_stage_hex, a.v1_inputs),
        "new_contigs": len(distinct_stage_uids) - len(stored_owned_uids),
        "new_bases": new_bases,
        "total_distinct_contigs": len(distinct_stage_uids),
        "total_distinct_bases": total_stage_bases,
    }

    # Spot-check stage.fasta_bytes on first 10 objects against lock.json
    for o in objects_json[:10]:
        p = Path(o["path"])
        acc = p.name[:-6]
        got_fb = stage.fasta_bytes(p)
        assert got_fb == lock[acc]["bytes"], f"{acc}: fasta_bytes {got_fb} != lock {lock[acc]['bytes']}"

    tot_fasta_bytes = sum(lock[acc]["bytes"] for acc in lock)
    tot_gz_bytes = sum(lock[acc]["gz_bytes"] for acc in lock)
    max_obj_acc = max(lock, key=lambda k: lock[k]["bytes"])
    max_obj_bytes = lock[max_obj_acc]["bytes"]

    ingest_out = {
        "objects": len(objects_json),
        "total_gz_bytes": tot_gz_bytes,
        "total_fasta_bytes": tot_fasta_bytes,
        "held_bases": held_bases,
        "cpu_min_one_process_b52": round(stage.ingest_cpu_min(tot_fasta_bytes, held_bases), 3),
        "cpu_min_one_process_unadjusted": round(stage.ingest_cpu_min(tot_fasta_bytes, 0), 3),
        "max_object_accession": max_obj_acc,
        "max_object_fasta_bytes": max_obj_bytes,
        "max_object_cpu_min": round(stage.ingest_cpu_min(max_obj_bytes, 0), 4),
        "shard_margin_cpu_min": stage.SHARD_MARGIN * prof.cpu_min_per_process,
        "runs_in_one_process": stage.ingest_cpu_min(tot_fasta_bytes, held_bases) <= stage.SHARD_MARGIN * prof.cpu_min_per_process,
    }

    # Gene-call plan (exact DEC-26 + DEC-27 ownership rules)
    # Least genome_uid in G1 owns each contig not owned by an earlier stage
    contig_least_guid: dict[bytes, bytes] = {}
    parsed_genomes: list[tuple[bytes, dict]] = []
    for g in genomes_json:
        cuids = [bytes.fromhex(u) for u in g["contig_uids"]]
        guid = ids.genome_uid(cuids)
        parsed_genomes.append((guid, g))
        for cu in cuids:
            if cu not in stored_owned_uids:
                if cu not in contig_least_guid or guid < contig_least_guid[cu]:
                    contig_least_guid[cu] = guid
    parsed_genomes.sort(key=lambda t: t[0])

    mode_counts = Counter()
    items = []
    for guid, g in parsed_genomes:
        cuids = [bytes.fromhex(u) for u in g["contig_uids"]]
        contigs = [(contig_len[cu], False) for cu in cuids]
        bases = sum(n for n, _ in contigs)
        mode_counts[genecall.mode_of(g["kind"], bases)] += 1
        writes_any = any((cu not in stored_owned_uids and contig_least_guid[cu] == guid) for cu in cuids)
        if not writes_any:
            items.append((guid, 0, bases, 0.0, False))
        else:
            items.append((
                guid,
                genecall.genome_cost_us("prokaryote", contigs),
                bases,
                genecall.genome_genes("prokaryote", bases),
                True,
            ))

    plan = genecall.cut_plan(lambda: iter(items), stage.SHARD_MARGIN * prof.cpu_min_per_process)
    genecall_out = {
        "modes": dict(mode_counts),
        "genomes": plan["genomes"],
        "genomes_not_called": plan["genomes_not_called"],
        "bases": plan["bases"],
        "bases_called": plan["bases_called"],
        "estimate_cpu_min": plan["estimate_cpu_min"],
        "batches": len(plan["batches"]),
        "largest_batch_cpu_min": max(b["estimate_cpu_min"] for b in plan["batches"]),
        "largest_genome": plan["largest_genome"],
        "genes_expected": plan["genes_expected"],
        "merge_estimate_cpu_min": plan["merge_estimate_cpu_min"],
        "rates": plan["rates"],
    }

    # Sample every Nth callable genome (ordered by accession for reproducibility)
    seen_guids: set[bytes] = set()
    callable_accs: list[str] = []
    for acc in sorted(lock):
        cuids = acc_uids[acc]
        guid = ids.genome_uid(cuids)
        if guid in seen_guids:
            continue
        seen_guids.add(guid)
        if any((cu not in stored_owned_uids and contig_least_guid.get(cu) == guid) for cu in cuids):
            callable_accs.append(acc)

    tot = defaultdict(lambda: {"genomes": 0, "contigs": 0, "bases": 0, "genes": 0, "partial_genes": 0, "cpu_s": 0.0})
    prot_seq: dict[bytes, str] = {}
    sample_genes_written = 0
    for idx_k, acc in enumerate(callable_accs):
        if idx_k % a.every != 0:
            continue
        cuids = acc_uids[acc]
        guid = ids.genome_uid(cuids)
        fp = Path(lock[acc]["file"])
        seqs: dict[bytes, str] = {}
        for rec in fasta.iter_records(fp):
            cu = ids.contig_uid(rec.seq)
            if cu not in seqs:
                seqs[cu] = rec.seq.upper().decode("ascii")
        g_in = GenomeInput(tuple(cuids), 11, kind="prokaryote", assembly_quality=sel_by_acc[acc]["assembly_quality"])
        t0 = cpu_self()
        calls = genecall._call(g_in, seqs)
        dt = cpu_self() - t0
        # Keep only contigs owned by this genome (matching genecall step)
        owned_calls = [
            c for c in calls
            if c["contig_uid"] not in stored_owned_uids and contig_least_guid.get(c["contig_uid"]) == guid
        ]
        partials = sum(1 for c in owned_calls if c["dist_to_contig_end_bp"] == 0)
        sample_genes_written += len(owned_calls)
        for c in owned_calls:
            aa = c["aa"]
            pu = ids.protein_uid(aa)
            prot_seq.setdefault(pu, aa)
        gbases = sum(len(s) for s in seqs.values())
        size_bin = "under_2mbp" if gbases < 2_000_000 else "ge_2mbp"
        res_name = sel_by_acc[acc]["resource"]
        for key in ("all", size_bin, f"res:{res_name}"):
            b = tot[key]
            b["genomes"] += 1
            b["contigs"] += len(cuids)
            b["bases"] += gbases
            b["genes"] += len(owned_calls)
            b["partial_genes"] += partials
            b["cpu_s"] += dt

    for b in tot.values():
        if b["bases"]:
            b["s_per_mbp"] = round(b["cpu_s"] / b["bases"] * 1e6, 4)
            b["genes_per_mbp"] = round(b["genes"] / b["bases"] * 1e6, 2)
        if b["genes"]:
            b["partial_gene_share"] = round(b["partial_genes"] / b["genes"], 4)
        b["cpu_s"] = round(b["cpu_s"], 3)

    known_info = known_proteins(a.root, set(prot_seq))
    nr90_info = run_sample_linclust(prot_seq, threads=2)

    crispr_fastas = [Path(lock[acc]["file"]) for idx_k, acc in enumerate(callable_accs) if idx_k % a.crispr_every == 0]
    crispr_sample = run_sample_crispr(crispr_fastas)

    tot_sample_res = sum(len(s) for s in prot_seq.values())
    sample_out = {
        "every": a.every,
        "callable_genomes_in_stage": len(callable_accs),
        "by": dict(tot),
        "genes_written": sample_genes_written,
        "distinct_proteins": len(prot_seq),
        "residues_distinct": tot_sample_res,
        "distinct_per_gene": round(len(prot_seq) / sample_genes_written, 4) if sample_genes_written else 0.0,
        "residues_per_distinct": round(tot_sample_res / len(prot_seq), 2) if prot_seq else 0.0,
        "known_in_stored_stages": known_info,
        "known_protein_share": round(known_info["total_known"] / len(prot_seq), 4) if prot_seq else 0.0,
        "nr90_linclust": nr90_info,
        "crispr_sample": crispr_sample,
    }

    # Derived estimates for stage G1 and for the full 41-stage (410.33 Gbp, 162,653-genome) non-GTDB set
    s_all = sample_out["by"]["all"]
    g_per_mbp = s_all["genes_per_mbp"]
    dpg = sample_out["distinct_per_gene"]
    rpd = sample_out["residues_per_distinct"]
    known_share = sample_out["known_protein_share"]
    rep_res_frac_sample = nr90_info["rep_residue_fraction"]
    # Also report rN.3's measured one-per-species NR90 residue scan share (0.831 on P1+P2, B-47_inferred_domains.md §2c)
    rep_res_frac_rn3 = 0.831
    sr = stage.shard_residues(prof)

    g1_genes = g_per_mbp * (new_bases / 1e6)
    g1_distinct_hi = g1_genes * dpg
    g1_res_hi = g1_distinct_hi * rpd
    g1_res_after_copy = g1_res_hi * (1.0 - known_share)
    # Stage-local NR90 `represent` step CPU-min scaled from rN.3's 992.2 CPU-min at 10,620,684 proteins (exponent 1.244)
    g1_represent_cpu_min = 992.2 * ((g1_distinct_hi / 10_620_684) ** stage.CLUSTER_EXPONENT)
    # IGR + tapes CPU-min scaled per Gbp from P2's measured run records
    p2 = p2_measured(a.root)
    p2_igr_tapes_per_gbp = (p2["cpu_min_by_step"].get("igr", 0.0) + p2["cpu_min_by_step"].get("tapes", 0.0)) / (
        p2["stored_bases"] / 1e9
    )
    g1_igr_tapes_cpu_min = p2_igr_tapes_per_gbp * (new_bases / 1e9)

    g1_pfam_full_cpu_min = stage.scan_cpu_min(int(round(g1_res_hi)))
    g1_pfam_after_b47b_cpu_min = stage.scan_cpu_min(int(round(g1_res_after_copy)))
    g1_pfam_nr90_sample_cpu_min = stage.scan_cpu_min(int(round(g1_res_after_copy * rep_res_frac_sample)))
    g1_pfam_nr90_rn3_cpu_min = stage.scan_cpu_min(int(round(g1_res_after_copy * rep_res_frac_rn3)))
    g1_crispr_def_code = stage.crispr_cpu_min(new_bases)
    g1_crispr_cas_code = g1_crispr_def_code * stage.CRISPR_NEAR_CAS_COST
    g1_crispr_def_measured_rate = (crispr_sample["cpu_s_per_gbp"] * (new_bases / 1e9)) / 60.0
    g1_crispr_cas_measured_rate = g1_crispr_def_measured_rate * stage.CRISPR_NEAR_CAS_COST

    # Root storage for G1: base P2 rate (0.78 GB/Gbp) + extra contig index/occurrence/label rows for 1.23M contigs
    # P2 had ~2,300 contigs in 9.50 Gbp; G1 has 1,227,287 contigs (~224 B/contig across occurrence 120 + index 56 + label 48)
    extra_contig_bytes_g1 = max(0, len(distinct_stage_uids) - 2300) * (120 + 56 + 48)
    g1_root_gb_annotated = p2["disk_gb_per_gbp"] * (new_bases / 1e9) + extra_contig_bytes_g1 / 1e9
    # Without l1_domains (~0.08 GB/Gbp on P2: 758 MB / 9.50 Gbp)
    p2_domains_gb_per_gbp = (p2["disk_bytes_by_dir"].get("l1", 0) / 1e9) / (p2["stored_bases"] / 1e9)
    g1_domain_only_gb = 0.0798 * (new_bases / 1e9)  # P2 l1_domains 757.9 MB / 9.50 Gbp

    derived_g1 = {
        "new_bases": new_bases,
        "genes_est": round(g1_genes),
        "distinct_proteins_upper": round(g1_distinct_hi),
        "residues_upper": round(g1_res_hi),
        "residues_after_b47b_copy": round(g1_res_after_copy),
        "shard_residues": sr,
        "ingest_cpu_min": ingest_out["cpu_min_one_process_b52"],
        "genecall_plan_cpu_min": genecall_out["estimate_cpu_min"],
        "genecall_at_sample_rate_cpu_min": round(
            s_all["s_per_mbp"] * genecall.STEP_OVERHEAD * (new_bases / 1e6) / 60.0, 2
        ),
        "igr_and_tapes_cpu_min": round(g1_igr_tapes_cpu_min, 2),
        "crispr_default_code_cpu_min": round(g1_crispr_def_code, 2),
        "crispr_near_cas_code_cpu_min": round(g1_crispr_cas_code, 2),
        "crispr_default_at_sample_rate_cpu_min": round(g1_crispr_def_measured_rate, 2),
        "crispr_near_cas_at_sample_rate_cpu_min": round(g1_crispr_cas_measured_rate, 2),
        "represent_nr90_cpu_min": round(g1_represent_cpu_min, 2),
        "pfam_full_cpu_min": round(g1_pfam_full_cpu_min, 1),
        "pfam_full_shards": math.ceil(g1_res_hi / sr),
        "pfam_after_b47b_cpu_min": round(g1_pfam_after_b47b_cpu_min, 1),
        "pfam_after_b47b_shards": math.ceil(g1_res_after_copy / sr),
        "pfam_nr90_at_sample_bound_cpu_min": round(g1_pfam_nr90_sample_cpu_min, 1),
        "pfam_nr90_at_sample_bound_shards": math.ceil((g1_res_after_copy * rep_res_frac_sample) / sr),
        "pfam_nr90_at_rn3_rate_cpu_min": round(g1_pfam_nr90_rn3_cpu_min, 1),
        "pfam_nr90_at_rn3_rate_shards": math.ceil((g1_res_after_copy * rep_res_frac_rn3) / sr),
        "root_gb_annotated_est": round(g1_root_gb_annotated, 2),
        "root_gb_called_est": round(g1_root_gb_annotated - g1_domain_only_gb, 2),
        "extra_contig_metadata_gb": round(extra_contig_bytes_g1 / 1e9, 3),
    }

    # Full non-GTDB set (162,653 genomes, 410,331,702,192 bp, 48,999,006 contigs, 41 stages G1..G41)
    full_bp = tsum["resources"]["ALL_NON_GTDB"]["bp"]
    full_genomes = tsum["resources"]["ALL_NON_GTDB"]["genomes"]
    full_contigs = tsum["resources"]["ALL_NON_GTDB"]["contigs"]
    n_stages = tsum["stage_proposal_genus_coherent"]["n_actual"]
    # FASTA bytes ratio from G1 (headers `><ID>_000000000001\n` + line-wrapped bases)
    fasta_byte_factor = tot_fasta_bytes / total_stage_bases
    full_fasta_bytes_est = round(full_bp * fasta_byte_factor)

    # Note: G1 is Archaea-heavy (mean genome 1.41 Mbp, higher coding density/shorter proteins than 2.52 Mbp overall).
    # Compute full-set numbers BOTH from P2/code prokaryote constants (917.4 genes/Mbp or P2's 899.6 genes/Mbp,
    # 0.9384 distinct/gene, 323.6 aa/protein — typical 2.5–3.2 Mbp Bacteria+Archaea) AND from G1's sample!
    p2_genes_full = p2["genes_per_mbp"] * (full_bp / 1e6)
    p2_distinct_per_stage = (p2_genes_full / n_stages) * p2["distinct_per_gene"]
    p2_tot_distinct_stage_sum = p2_distinct_per_stage * n_stages
    p2_tot_res = p2_tot_distinct_stage_sum * p2["residues_per_distinct"]

    g1rate_genes_full = g_per_mbp * (full_bp / 1e6)
    g1rate_tot_distinct_stage_sum = g1rate_genes_full * dpg
    g1rate_tot_res = g1rate_tot_distinct_stage_sum * rpd

    full_ingest_cpu_h = (n_stages * stage.INGEST_FIXED_CPU_MIN + full_fasta_bytes_est * stage.INGEST_CPU_MIN / stage.INGEST_BP) / 60.0
    full_genecall_code_cpu_h = (
        n_stages * genecall.PROCESS_S
        + (full_bp / 1e6) * genecall.CALL_S_PER_MBP[("single", False)] * genecall.STEP_OVERHEAD
    ) / 3600.0
    full_genecall_sample_cpu_h = (
        n_stages * genecall.PROCESS_S
        + (full_bp / 1e6) * s_all["s_per_mbp"] * genecall.STEP_OVERHEAD
    ) / 3600.0
    full_igr_tapes_cpu_h = (p2_igr_tapes_per_gbp * (full_bp / 1e9)) / 60.0
    full_crispr_def_code_cpu_h = stage.crispr_cpu_min(full_bp) / 60.0
    full_crispr_cas_code_cpu_h = full_crispr_def_code_cpu_h * stage.CRISPR_NEAR_CAS_COST
    full_crispr_def_sample_cpu_h = (crispr_sample["cpu_s_per_gbp"] * (full_bp / 1e9)) / 3600.0
    full_crispr_cas_sample_cpu_h = full_crispr_def_sample_cpu_h * stage.CRISPR_NEAR_CAS_COST

    def pfam_scenarios(tot_res: float, distinct_per_stage: float) -> dict:
        rep_step_cpu_h = (
            n_stages * 992.2 * ((distinct_per_stage / 10_620_684) ** stage.CLUSTER_EXPONENT)
        ) / 60.0
        full_scan_h = stage.scan_cpu_min(int(round(tot_res))) / 60.0
        after_copy_res = tot_res * (1.0 - known_share)
        after_copy_h = stage.scan_cpu_min(int(round(after_copy_res))) / 60.0
        nr90_rn3_h = stage.scan_cpu_min(int(round(after_copy_res * rep_res_frac_rn3))) / 60.0
        nr90_sample_h = stage.scan_cpu_min(int(round(after_copy_res * rep_res_frac_sample))) / 60.0
        return {
            "total_stage_distinct_proteins": round(distinct_per_stage * n_stages),
            "total_stage_residues": round(tot_res),
            "represent_41_stages_cpu_h": round(rep_step_cpu_h, 1),
            "pfam_full_unadjusted_cpu_h": round(full_scan_h, 1),
            "pfam_after_b47b_copy_cpu_h": round(after_copy_h, 1),
            "pfam_nr90_rn3_rate_cpu_h": round(nr90_rn3_h, 1),
            "pfam_nr90_sample_bound_cpu_h": round(nr90_sample_h, 1),
        }

    sc_p2 = pfam_scenarios(p2_tot_res, p2_distinct_per_stage)
    sc_g1 = pfam_scenarios(g1rate_tot_res, g1rate_tot_distinct_stage_sum / n_stages)

    # Depth 1 (`called`, DEC-54): ingest + genecall + default CRISPR + IGR + tapes
    depth_called_code_cpu_h = full_ingest_cpu_h + full_genecall_code_cpu_h + full_crispr_def_code_cpu_h + full_igr_tapes_cpu_h
    depth_called_sample_cpu_h = full_ingest_cpu_h + full_genecall_sample_cpu_h + full_crispr_def_sample_cpu_h + full_igr_tapes_cpu_h

    depth_called_lo, depth_called_hi = sorted([depth_called_sample_cpu_h, depth_called_code_cpu_h])

    # Depth 3 (`annotated` with NR90 representatives, B-47a..e): depth_called + represent + Pfam(NR90) + near-Cas CRISPR
    depth3_nr90_lo_cpu_h = (
        depth_called_lo + min(full_crispr_cas_sample_cpu_h, full_crispr_cas_code_cpu_h)
        + min(sc_p2["represent_41_stages_cpu_h"] + sc_p2["pfam_nr90_rn3_rate_cpu_h"],
              sc_g1["represent_41_stages_cpu_h"] + sc_g1["pfam_nr90_rn3_rate_cpu_h"])
    )
    depth3_nr90_hi_cpu_h = (
        depth_called_hi + max(full_crispr_cas_sample_cpu_h, full_crispr_cas_code_cpu_h)
        + max(sc_p2["represent_41_stages_cpu_h"] + sc_p2["pfam_nr90_sample_bound_cpu_h"],
              sc_g1["represent_41_stages_cpu_h"] + sc_g1["pfam_nr90_sample_bound_cpu_h"])
    )
    depth3_full_lo_cpu_h = (
        depth_called_lo + min(full_crispr_cas_sample_cpu_h, full_crispr_cas_code_cpu_h)
        + min(sc_p2["pfam_after_b47b_copy_cpu_h"], sc_g1["pfam_after_b47b_copy_cpu_h"])
    )
    depth3_full_hi_cpu_h = (
        depth_called_hi + max(full_crispr_cas_sample_cpu_h, full_crispr_cas_code_cpu_h)
        + max(sc_p2["pfam_full_unadjusted_cpu_h"], sc_g1["pfam_full_unadjusted_cpu_h"])
    )

    # Pod thread-hours/day: 5*(12*40 + 12*56) + 2*(24*56) = 8,448 thread-h/week = 1,206.86 thread-h/day (or 960 at 40 threads flat)
    pod_th_per_day_sched = (5 * (12 * 40 + 12 * 56) + 2 * (24 * 56)) / 7.0
    pod_th_per_day_40 = 40.0 * 24.0
    tacc_cores_per_node = 128.0  # TACC Lonestar6: 2x AMD EPYC 7763, 128 cores/node, 1 SU/node-h (card C-37)

    extra_contig_gb_full = max(0, full_contigs - n_stages * 5018) * (120 + 56 + 48) / 1e9
    full_root_annotated_gb = p2["disk_gb_per_gbp"] * (full_bp / 1e9) + extra_contig_gb_full
    # Stage-rep table (~48 B/protein for ~367M proteins ≈ 17.6 GB) + inferred_from column (~7 GB)
    full_stage_rep_plus_inferred_gb = 24.6
    full_domains_gb = 0.0790 * (full_bp / 1e9)  # P2 l1_domains 757.9 MB / 9.596 Gbp = 0.0790 GB/Gbp

    full_set_out = {
        "genomes": full_genomes,
        "bp": full_bp,
        "gbp": round(full_bp / 1e9, 3),
        "contigs": full_contigs,
        "stages": n_stages,
        "est_fasta_bytes": full_fasta_bytes_est,
        "step_cpu_hours": {
            "ingest_b52": round(full_ingest_cpu_h, 2),
            "genecall_code_rate": round(full_genecall_code_cpu_h, 2),
            "genecall_sample_rate": round(full_genecall_sample_cpu_h, 2),
            "igr_and_tapes": round(full_igr_tapes_cpu_h, 2),
            "crispr_default_code_rate": round(full_crispr_def_code_cpu_h, 2),
            "crispr_default_sample_rate": round(full_crispr_def_sample_cpu_h, 2),
            "crispr_near_cas_code_rate": round(full_crispr_cas_code_cpu_h, 2),
            "crispr_near_cas_sample_rate": round(full_crispr_cas_sample_cpu_h, 2),
            "pfam_model_from_P2_composition": sc_p2,
            "pfam_model_from_G1_archaea_sample": sc_g1,
        },
        "totals": {
            "depth_called_cpu_h_range": [round(depth_called_lo, 1), round(depth_called_hi, 1)],
            "depth_called_pod_days_40_56": [
                round(depth_called_lo / pod_th_per_day_sched, 2),
                round(depth_called_hi / pod_th_per_day_sched, 2),
            ],
            "depth_called_pod_days_flat40": [
                round(depth_called_lo / pod_th_per_day_40, 2),
                round(depth_called_hi / pod_th_per_day_40, 2),
            ],
            "depth_called_lonestar6_node_h_SUs": [
                round(depth_called_lo / tacc_cores_per_node, 2),
                round(depth_called_hi / tacc_cores_per_node, 2),
            ],
            "depth3_nr90_cpu_h_range": [round(depth3_nr90_lo_cpu_h, 1), round(depth3_nr90_hi_cpu_h, 1)],
            "depth3_nr90_pod_days_40_56": [
                round(depth3_nr90_lo_cpu_h / pod_th_per_day_sched, 2),
                round(depth3_nr90_hi_cpu_h / pod_th_per_day_sched, 2),
            ],
            "depth3_nr90_pod_days_flat40": [
                round(depth3_nr90_lo_cpu_h / pod_th_per_day_40, 2),
                round(depth3_nr90_hi_cpu_h / pod_th_per_day_40, 2),
            ],
            "depth3_nr90_lonestar6_node_h_SUs": [
                round(depth3_nr90_lo_cpu_h / tacc_cores_per_node, 1),
                round(depth3_nr90_hi_cpu_h / tacc_cores_per_node, 1),
            ],
            "depth3_full_pfam_cpu_h_range": [round(depth3_full_lo_cpu_h, 1), round(depth3_full_hi_cpu_h, 1)],
            "depth3_full_pfam_pod_days_40_56": [
                round(depth3_full_lo_cpu_h / pod_th_per_day_sched, 2),
                round(depth3_full_hi_cpu_h / pod_th_per_day_sched, 2),
            ],
            "depth3_full_pfam_lonestar6_node_h_SUs": [
                round(depth3_full_lo_cpu_h / tacc_cores_per_node, 1),
                round(depth3_full_hi_cpu_h / tacc_cores_per_node, 1),
            ],
        },
        "storage_gb": {
            "base_410gbp_at_p2_rate": round(p2["disk_gb_per_gbp"] * (full_bp / 1e9), 1),
            "extra_49M_contig_metadata": round(extra_contig_gb_full, 1),
            "depth_called_root_gb": round(full_root_annotated_gb - full_domains_gb, 1),
            "depth3_full_pfam_root_gb": round(full_root_annotated_gb, 1),
            "depth3_nr90_root_gb": round(full_root_annotated_gb + full_stage_rep_plus_inferred_gb, 1),
        },
    }

    out = {
        "profile_cpu_min_per_process": prof.cpu_min_per_process,
        "shard_margin": stage.SHARD_MARGIN,
        "overlap": overlap_out,
        "ingest": ingest_out,
        "genecall": genecall_out,
        "sample": sample_out,
        "p2_measured": p2,
        "derived_stage_G1": derived_g1,
        "full_set_410gbp": full_set_out,
    }
    (a.out / "est.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({
        "shared_with_stored_stages": overlap_out["shared_distinct_contigs_with_stored_stages"],
        "intra_stage_shared_contigs": overlap_out["intra_stage_shared_contigs"],
        "ingest_cpu_min": ingest_out["cpu_min_one_process_b52"],
        "genecall_cpu_min": genecall_out["estimate_cpu_min"],
        "genecall_batches": genecall_out["batches"],
        "sample_s_per_mbp": s_all["s_per_mbp"],
        "sample_nr90_rep_res_frac": rep_res_frac_sample,
        "crispr_sample_s_per_gbp": crispr_sample["cpu_s_per_gbp"],
    }, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Card C-48d stage 3: skani 95-96% ANI species-representative selection,
~10 Gbp stage proposals (prokaryote vs eukaryote), and locked pilot manifest writer.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from biobox import fasta, ids


def _quality_rank(rec: dict) -> tuple:
    comp = float(rec.get("completeness", 0.0))
    cont = float(rec.get("contamination", 0.0))
    qscore = comp - 5.0 * cont
    n50 = int(rec.get("n50", 0))
    bp = int(rec.get("bp", 0))
    return (-qscore, -n50, -bp, str(rec["genome_id"]))


def select_species_reps(
    candidates: list[dict],
    skani_edges: list[tuple[str, str, float, float, float]],
    ani_threshold: float = 95.0,
    min_af: float = 50.0,
    ref_ids: set[str] | None = None,
) -> dict:
    """Dereplicate candidate genomes against `ref_ids` and against each other at ANI/AF.

    skani_edges: list of (id_a, id_b, ani, af_a, af_b).
    """
    ref_set = set(ref_ids or ())
    cand_by_id = {c["genome_id"]: c for c in candidates}
    dropped_by_ref: dict[str, tuple[str, float]] = {}
    adj: dict[str, set[str]] = {cid: set() for cid in cand_by_id}

    for a, b, ani, af_a, af_b in skani_edges:
        if ani < ani_threshold or max(af_a, af_b) < min_af:
            continue
        if a in ref_set and b in cand_by_id:
            prev = dropped_by_ref.get(b)
            if prev is None or ani > prev[1]:
                dropped_by_ref[b] = (a, ani)
        elif b in ref_set and a in cand_by_id:
            prev = dropped_by_ref.get(a)
            if prev is None or ani > prev[1]:
                dropped_by_ref[a] = (b, ani)
        elif a in cand_by_id and b in cand_by_id and a != b:
            adj[a].add(b)
            adj[b].add(a)

    surviving = [c for cid, c in cand_by_id.items() if cid not in dropped_by_ref]
    surviving.sort(key=_quality_rank)

    novel_reps: list[dict] = []
    clustered_into: dict[str, str] = {}
    assigned: set[str] = set()

    for rec in surviving:
        gid = rec["genome_id"]
        if gid in assigned:
            continue
        novel_reps.append(rec)
        assigned.add(gid)
        for nbr in sorted(adj.get(gid, ())):
            if nbr not in assigned:
                assigned.add(nbr)
                clustered_into[nbr] = gid

    novel_reps.sort(key=lambda r: (r.get("taxonomy", ""), r["genome_id"]))
    return {
        "input_candidates": len(candidates),
        "dropped_by_ref": sorted(dropped_by_ref.keys()),
        "dropped_by_ref_detail": {k: {"ref_id": v[0], "ani": v[1]} for k, v in sorted(dropped_by_ref.items())},
        "clustered_into_novel_rep": dict(sorted(clustered_into.items())),
        "novel_reps": novel_reps,
        "novel_rep_count": len(novel_reps),
        "novel_rep_bases": sum(int(r.get("bp", 0)) for r in novel_reps),
    }


def _stem_from_path(p: str) -> str:
    name = Path(p).name
    for suf in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if name.endswith(suf):
            return name[: -len(suf)]
    return Path(name).stem


def run_skani_dedup(
    candidate_fastas: dict[str, Path],
    ref_fastas: dict[str, Path],
    work_dir: Path,
    threads: int = 4,
    ani_threshold: float = 95.0,
    min_af: float = 50.0,
    skani_bin: str = "skani",
) -> dict:
    """Run `skani sketch --separate-sketches`, `skani dist` (cand vs ref), and `skani triangle -E` (cand vs cand)."""
    work_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    edges: list[tuple[str, str, float, float, float]] = []

    cand_list = work_dir / "cands.list"
    cand_list.write_text("\n".join(str(candidate_fastas[k]) for k in sorted(candidate_fastas)) + "\n")

    cand_sk_dir = work_dir / "cand_sketches"
    if cand_sk_dir.exists():
        shutil.rmtree(cand_sk_dir)
    subprocess.run(
        [skani_bin, "sketch", "--separate-sketches", "-t", str(threads), "-l", str(cand_list), "-o", str(cand_sk_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    cand_sk_list = work_dir / "cand_sketches.list"
    cand_sk_files = sorted(cand_sk_dir.glob("*.sketch"))
    cand_sk_list.write_text("\n".join(str(p) for p in cand_sk_files) + "\n")

    # 1) Candidate vs Reference
    if ref_fastas:
        ref_list = work_dir / "refs.list"
        ref_list.write_text("\n".join(str(ref_fastas[k]) for k in sorted(ref_fastas)) + "\n")
        ref_sk_dir = work_dir / "ref_sketches"
        if ref_sk_dir.exists():
            shutil.rmtree(ref_sk_dir)
        subprocess.run(
            [skani_bin, "sketch", "--separate-sketches", "-t", str(threads), "-l", str(ref_list), "-o", str(ref_sk_dir)],
            check=True,
            capture_output=True,
            text=True,
        )
        ref_sk_list = work_dir / "ref_sketches.list"
        ref_sk_files = sorted(ref_sk_dir.glob("*.sketch"))
        ref_sk_list.write_text("\n".join(str(p) for p in ref_sk_files) + "\n")

        dist_out = work_dir / "cand_vs_ref.tsv"
        subprocess.run(
            [
                skani_bin, "dist", "-t", str(threads), "--min-af", str(min_af),
                "--ql", str(cand_sk_list), "--rl", str(ref_sk_list), "-o", str(dist_out),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        if dist_out.is_file():
            with open(dist_out, encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    r_id = _stem_from_path(row["Ref_file"])
                    q_id = _stem_from_path(row["Query_file"])
                    edges.append((
                        r_id, q_id,
                        float(row["ANI"]),
                        float(row["Align_fraction_ref"]),
                        float(row["Align_fraction_query"]),
                    ))

    # 2) Candidate vs Candidate
    if len(cand_sk_files) >= 2:
        tri_out = work_dir / "cand_triangle.tsv"
        subprocess.run(
            [
                skani_bin, "triangle", "-E", "-t", str(threads), "--min-af", str(min_af),
                "-l", str(cand_sk_list), "-o", str(tri_out),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        if tri_out.is_file():
            with open(tri_out, encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    r_id = _stem_from_path(row["Ref_file"])
                    q_id = _stem_from_path(row["Query_file"])
                    if r_id != q_id:
                        edges.append((
                            r_id, q_id,
                            float(row["ANI"]),
                            float(row["Align_fraction_ref"]),
                            float(row["Align_fraction_query"]),
                        ))

    return {
        "candidates": len(candidate_fastas),
        "references": len(ref_fastas),
        "edges_count": len(edges),
        "edges": edges,
        "wall_sec": round(time.monotonic() - t0, 2),
    }


def _genus_key(tax: str) -> str:
    parts = [p for p in (tax or "").split(";") if p.strip()]
    for p in parts:
        if p.startswith("g__"):
            idx = parts.index(p)
            return ";".join(parts[: idx + 1])
    return ";".join(parts) if parts else "Unclassified"


def propose_stages(
    reps: list[dict],
    prefix: str,
    target_bp: int = 10_000_000_000,
    redistributable: bool = False,
    kind: str = "prokaryote",
    depth_recommendation: str = "called",
    gene_caller_note: str = "pyrodigal single mode (code 11; screen Patescibacteriota/Clostridia for code 25/4)",
) -> list[dict]:
    """Partition representatives into ~target_bp stages in deterministic (taxonomy, genome_id) order."""
    if not reps:
        return []
    ordered = sorted(reps, key=lambda r: (r.get("taxonomy", ""), r["genome_id"]))
    stages: list[dict] = []
    cur: list[dict] = []
    cur_bp = 0
    idx = 1

    # Group contiguous genera so genera stay whole unless a single genus exceeds target_bp
    groups: list[list[dict]] = []
    for r in ordered:
        gk = _genus_key(r.get("taxonomy", ""))
        if not groups or _genus_key(groups[-1][0].get("taxonomy", "")) != gk:
            groups.append([r])
        else:
            groups[-1].append(r)

    def _emit(items: list[dict], k: int) -> dict:
        return {
            "stage": f"{prefix}{k}",
            "kind": kind,
            "redistributable": redistributable,
            "depth_recommendation": depth_recommendation,
            "gene_caller_note": gene_caller_note,
            "genomes": len(items),
            "bp": sum(int(x.get("bp", 0)) for x in items),
            "gbp": round(sum(int(x.get("bp", 0)) for x in items) / 1e9, 4),
            "first_id": items[0]["genome_id"],
            "last_id": items[-1]["genome_id"],
            "genome_ids": [x["genome_id"] for x in items],
        }

    for grp in groups:
        grp_bp = sum(int(x.get("bp", 0)) for x in grp)
        if cur and cur_bp + grp_bp > target_bp:
            stages.append(_emit(cur, idx))
            idx += 1
            cur = []
            cur_bp = 0
        cur.extend(grp)
        cur_bp += grp_bp
    if cur:
        stages.append(_emit(cur, idx))
    return stages


def _mag_quality(rec: dict) -> str:
    comp = float(rec.get("completeness", 0.0))
    cont = float(rec.get("contamination", 0.0))
    if comp > 90.0 and cont < 5.0:
        return "mag_hq"
    if comp >= 50.0 and cont <= 10.0:
        return "mag_mq"
    return "contig"


def write_pilot_manifest(
    stage_rows: list[dict],
    extracted_files: dict[str, Path],
    out_dir: Path,
    source_name: str,
    source_release: str,
    license_str: str,
    redistributable: bool = False,
    retrieved_at: str = "2026-11-19T00:00:00Z",
) -> dict:
    """Write deterministic locked pilot stage manifest (lock.json, selection.tsv,
    genomes.json, objects.json, inputs_SUMS) with `redistributable=False` for JGI sets.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    ordered = sorted(stage_rows, key=lambda r: r["genome_id"])
    lock_genomes: dict[str, dict] = {}
    genomes_json: list[dict] = []
    objects_json: list[dict] = []
    seen_gids: set[bytes] = set()
    total_contigs = 0
    total_bp = 0

    sel_path = out_dir / "selection.tsv"
    with open(sel_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "source", "translation_table", "assembly_quality", "bp", "completeness", "contamination", "taxonomy", "ecosystem"])
        for rec in ordered:
            gid = rec["genome_id"]
            fa_path = extracted_files[gid]
            raw_gz = fa_path.read_bytes()
            gz_sha = hashlib.sha256(raw_gz).hexdigest()
            uids: list[bytes] = []
            bp = 0
            for r in fasta.iter_records(fa_path):
                uids.append(ids.contig_uid(r.seq))
                bp += len(r.seq)
            kept = list(dict.fromkeys(uids))
            g_uid = ids.genome_uid(kept)
            qual = _mag_quality(rec)
            w.writerow([
                gid, rec.get("source", source_name), "11", qual, bp,
                rec.get("completeness", ""), rec.get("contamination", ""),
                rec.get("taxonomy", ""), rec.get("ecosystem", ""),
            ])
            lock_genomes[gid] = {
                "file": str(fa_path),
                "gz_bytes": len(raw_gz),
                "gz_sha256": gz_sha,
                "records": len(uids),
                "bp": bp,
                "retrieved_utc": retrieved_at,
            }
            if g_uid not in seen_gids:
                seen_gids.add(g_uid)
                genomes_json.append({
                    "contig_uids": [u.hex() for u in kept],
                    "genetic_code": 11,
                    "kind": "prokaryote",
                    "circular": [],
                    "assembly_quality": qual,
                })
                total_contigs += len(kept)
                total_bp += bp
            objects_json.append({
                "path": str(fa_path),
                "source_name": source_name,
                "source_release": source_release,
                "license": license_str,
                "redistributable": bool(redistributable),
                "retrieved_at": retrieved_at,
            })

    lock_path = out_dir / "lock.json"
    genomes_path = out_dir / "genomes.json"
    objects_path = out_dir / "objects.json"
    lock_path.write_text(json.dumps({"genomes": lock_genomes}, indent=2, sort_keys=True) + "\n")
    genomes_path.write_text(json.dumps(genomes_json) + "\n")
    objects_path.write_text(json.dumps(objects_json, indent=1) + "\n")

    sums_lines = []
    for fname in ("genomes.json", "lock.json", "objects.json", "selection.tsv"):
        digest = hashlib.sha256((out_dir / fname).read_bytes()).hexdigest()
        sums_lines.append(f"{digest}  {fname}")
    (out_dir / "inputs_SUMS").write_text("\n".join(sums_lines) + "\n")

    return {
        "genomes": len(genomes_json),
        "contigs": total_contigs,
        "bp": total_bp,
        "objects": len(objects_json),
        "redistributable": bool(redistributable),
        "sums": sums_lines,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-48d skani species-rep selection and stage proposal")
    ap.add_argument("--candidates-tsv", type=Path, required=True)
    ap.add_argument("--cand-dir", type=Path, required=True)
    ap.add_argument("--ref-dir", type=Path, default=None)
    ap.add_argument("--max-refs", type=int, default=500)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--pilot-out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--ani", type=float, default=95.0)
    ap.add_argument("--min-af", type=float, default=50.0)
    args = ap.parse_args(argv)

    cand_rows: list[dict] = []
    cand_fastas: dict[str, Path] = {}
    with open(args.candidates_tsv, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gid = row.get("genome_id") or row.get("bin_id") or ""
            fa = args.cand_dir / f"{gid}.fna.gz"
            if gid and fa.is_file():
                row["genome_id"] = gid
                row["bp"] = int(row.get("bp") or 0)
                row["n50"] = int(row.get("n50") or 0)
                row["completeness"] = float(row.get("completeness") or 0.0)
                row["contamination"] = float(row.get("contamination") or 0.0)
                cand_rows.append(row)
                cand_fastas[gid] = fa

    ref_fastas: dict[str, Path] = {}
    if args.ref_dir is not None and args.ref_dir.is_dir():
        for p in sorted(args.ref_dir.glob("*.fa.gz"))[: args.max_refs]:
            ref_fastas[_stem_from_path(p.name)] = p
        for p in sorted(args.ref_dir.glob("*.fna.gz"))[: max(0, args.max_refs - len(ref_fastas))]:
            ref_fastas[_stem_from_path(p.name)] = p

    sk_res = run_skani_dedup(
        cand_fastas, ref_fastas, args.work_dir / "skani", threads=args.threads, ani_threshold=args.ani, min_af=args.min_af
    )
    sel = select_species_reps(
        cand_rows, sk_res["edges"], ani_threshold=args.ani, min_af=args.min_af, ref_ids=set(ref_fastas.keys())
    )
    stages = propose_stages(sel["novel_reps"], prefix="J", target_bp=10_000_000_000, redistributable=False)
    pilot_meta = write_pilot_manifest(
        sel["novel_reps"],
        cand_fastas,
        args.pilot_out,
        source_name="jgi_gems_umag_iuig",
        source_release="2026-11-18",
        license_str="JGI Data Utilization Policy (redistributable=false)",
        redistributable=False,
    )
    out = {
        "skani": {k: v for k, v in sk_res.items() if k != "edges"},
        "selection": {k: v for k, v in sel.items() if k != "novel_reps"},
        "stages": [{k: v for k, v in s.items() if k != "genome_ids"} for s in stages],
        "pilot_manifest": pilot_meta,
    }
    (args.work_dir / "skani_stage_summary.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

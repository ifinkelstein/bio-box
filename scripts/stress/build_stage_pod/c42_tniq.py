"""Card C-42 step (6): end-to-end TniQ HMM search over the pilot's depth-1 proteins + +-10-gene Pfam round.

    python -I c42_tniq.py run --pilot-root /tmp/if2453/s3/spire/root --inputs-dir /tmp/if2453/s3/spire/inputs
        --selection /tmp/if2453/s3/spire/select/pilot_selection.tsv
        --pfam-hmm /stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz
        --out-dir /tmp/if2453/s3/spire/tniq [--mag-uids-json /tmp/if2453/s3/spire/tniq/hit_sample_mag_uids.json]
        --threads 16

Demonstrates card C-38 item (2) ("search first, annotate on the hit", DEC-49 & DEC-54):
1. Extracts `TniQ` (`PF06527`, gathering cutoff) from Pfam 38.2 and searches every distinct protein of
   the pilot's depth-1 (`called`) stages at `Z = domZ = 1e8` (DEC-04).
2. Locates all `TniQ` anchor genes on the pilot's `contig_tapes`, extracts their strand-oriented
   `[-10, +10]` gene windows and nearby `l1_crispr` arrays (`<= 10 kb`), and collects the distinct
   `protein_uid`s in those hit windows.
3. Runs a Pfam 38.2 domain scan (`biobox.l1.domains.scan`, gathering cutoffs, `Z = domZ = 1e8`) ONLY
   on the hit-window proteins, classifying every TniQ locus (`I-F3_CAST`, `I-B_CAST`, `I-C_CAST`,
   `V-K_CAST`, `other_Cas_with_adaptation_or_cas3`, `Tn7_non_Cas`) and reporting CPU savings vs full Pfam.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import io as pyio
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pyhmmer

from biobox import ids, io
from biobox.l0seq.inputs import GenomeInput
from biobox.l1 import crispr, domains
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l2 import cluster as l2cluster

# Pfam 38.2 accession stems from `examples/crispr_cas_neighbourhoods/README.md` §1 (plus type V-K Cas12k)
ROLES_BY_ACC = {
    "PF06527": "TniQ",
    "PF09611": "Cas8f",
    "PF09614": "Cas5f",
    "PF09615": "Cas7f",
    "PF09618": "Cas6f",
    "PF24121": "Cas8b",
    "PF24122": "Cas8b",
    "PF09709": "Cas8c",
    "PF05107": "Cas7",
    "PF09704": "Cas5d",
    "PF09559": "Cas6",
    "PF01881": "Cas6",
    "PF10040": "Cas6",
    "PF01867": "Cas1",
    "PF09827": "Cas2",
    "PF09707": "Cas2",
    "PF00271": "Cas3_helicase",
    "PF18019": "Cas3",
    "PF18395": "Cas3",
    "PF21384": "Cas2_Cas3",
    "PF21802": "Cas3",
    "PF22590": "Cas3",
    "PF30455": "Cas3",
    "PF01930": "Cas4",
    "PF08722": "TnsA",
    "PF08721": "TnsA",
    "PF27165": "TnsB",
    "PF29687": "TnsB",
    "PF29688": "TnsB",
    "PF27135": "TnsB",
    "PF00665": "TnsB_rve",
    "PF11426": "TnsC",
    "PF05621": "TniB",
    "PF15978": "TnsD",
    "PF18623": "TnsE",
    "PF00478": "IMPDH",
    "PF03061": "YciA_4HBT",
    "PF01380": "GlmS_SIS",
    "PF13522": "GlmS_GATase6",
    "PF01078": "YifB_MgChel",
}
ADAPT_INTERF = {"Cas1", "Cas2", "Cas3", "Cas2_Cas3"}
ANY_CAS_ROLES = {"Cas8f", "Cas5f", "Cas7f", "Cas6f", "Cas8b", "Cas8c", "Cas7", "Cas5d", "Cas6",
                 "Cas1", "Cas2", "Cas3", "Cas2_Cas3", "Cas4", "Cas12k"}


def _acc_stem(model_acc: str) -> str:
    return (model_acc or "").split(".", 1)[0]


def _role_of(model_acc: str, model_name: str) -> str | None:
    stem = _acc_stem(model_acc)
    if stem in ROLES_BY_ACC:
        return ROLES_BY_ACC[stem]
    low = (model_name or "").lower()
    if "cas12k" in low or "cas12" in low or "cpfl" in low or "c2c5" in low:
        return "Cas12k"
    if low.startswith("cas_") or low.startswith("csy") or low.startswith("csd") or low.startswith("cse"):
        return f"Cas_other:{model_name}"
    return None


def extract_hit_windows(
    tapes_rows: Iterable[dict],
    hit_protein_uids: set[bytes],
    *,
    crispr_rows: Iterable[dict] = (),
    window: int = 10,
    crispr_bp: int = 10_000,
) -> tuple[list[dict], list[bytes]]:
    """Extract strand-oriented `[-window, +window]` gene windows around every anchor hit on `contig_tapes`."""
    crispr_by_cu: dict[bytes, list[dict]] = defaultdict(list)
    for cr in crispr_rows:
        crispr_by_cu[cr["contig_uid"]].append(cr)
    for cu in crispr_by_cu:
        crispr_by_cu[cu].sort(key=lambda r: (r["start"], r["end"]))

    windows: list[dict] = []
    req_proteins: set[bytes] = set()

    for row in tapes_rows:
        cu = row["contig_uid"]
        p_uids = row["protein_uids"]
        if not any(pu in hit_protein_uids for pu in p_uids if pu is not None):
            continue
        f_uids = row["feature_uids"]
        starts = row["starts"]
        ends = row["ends"]
        strands = row["strands"]
        kinds = row["kinds"]
        c_len = int(row.get("contig_len") or (max(ends) if ends else 0))
        n_slots = len(p_uids)

        for idx, pu in enumerate(p_uids):
            if pu is None or pu not in hit_protein_uids:
                continue
            a_strand = int(strands[idx])
            lo = max(0, idx - window)
            hi = min(n_slots - 1, idx + window)

            # Truncation along anchor's 5'->3' orientation:
            # If anchor strand == +1: 5' is lower index (idx < window), 3' is higher index (idx + window >= n_slots)
            # If anchor strand == -1: 5' is higher index (idx + window >= n_slots), 3' is lower index (idx < window)
            trunc_low = idx < window
            trunc_high = (idx + window) >= n_slots
            trunc_5p = trunc_low if a_strand == 1 else trunc_high
            trunc_3p = trunc_high if a_strand == 1 else trunc_low

            slots: list[dict] = []
            w_min_bp = starts[lo]
            w_max_bp = ends[hi]
            for j in range(lo, hi + 1):
                spu = p_uids[j]
                if spu is not None:
                    req_proteins.add(spu)
                slots.append({
                    "offset_rel": (j - idx) * a_strand,
                    "feature_uid": f_uids[j],
                    "protein_uid": spu,
                    "start": int(starts[j]),
                    "end": int(ends[j]),
                    "strand": int(strands[j]),
                    "same_strand": int(strands[j]) == a_strand,
                    "kind": kinds[j],
                })
            slots.sort(key=lambda s: s["offset_rel"])

            nearby_arrays = [
                cr for cr in crispr_by_cu.get(cu, ())
                if cr["end"] >= (w_min_bp - crispr_bp) and cr["start"] <= (w_max_bp + crispr_bp)
            ]

            windows.append({
                "contig_uid": cu,
                "contig_len": c_len,
                "anchor_feature_uid": f_uids[idx],
                "anchor_protein_uid": pu,
                "anchor_start": int(starts[idx]),
                "anchor_end": int(ends[idx]),
                "anchor_strand": a_strand,
                "n_contig_slots": n_slots,
                "window_slots": len(slots),
                "truncated_5p": trunc_5p,
                "truncated_3p": trunc_3p,
                "slots": slots,
                "crispr_arrays": nearby_arrays,
            })

    windows.sort(key=lambda w: (w["contig_uid"], w["anchor_start"], w["anchor_end"], w["anchor_strand"]))
    return windows, sorted(req_proteins)


def classify_tniq_locus(
    window_slots: list[dict],
    pfam_by_protein: dict[bytes, list[tuple[str, str, float]]],
    crispr_arrays: list[dict] | tuple = (),
) -> dict:
    """Classify a TniQ +-10-gene neighborhood using Pfam 38.2 roles (Peters et al. 2017)."""
    roles_present: set[str] = set()
    role_by_offset: dict[int, list[str]] = defaultdict(list)
    arch_by_offset: dict[int, str] = {}

    for s in window_slots:
        off = s["offset_rel"]
        pu = s["protein_uid"]
        hits = pfam_by_protein.get(pu, []) if pu is not None else []
        if not hits:
            arch_by_offset[off] = "-"
            continue
        arch_by_offset[off] = "+".join(dict.fromkeys(h[1] for h in hits))
        for acc, name, _ in hits:
            role = _role_of(acc, name)
            if role:
                roles_present.add(role)
                if role not in role_by_offset[off]:
                    role_by_offset[off].append(role)

    has_adapt = bool(roles_present & ADAPT_INTERF)
    has_any_cas = bool(roles_present & ANY_CAS_ROLES) or any(r.startswith("Cas_other:") for r in roles_present)
    has_tns = bool(roles_present & {"TnsA", "TnsB", "TnsC", "TniB"})

    if not has_adapt and ("Cas7f" in roles_present or ("Cas5f" in roles_present and "Cas6f" in roles_present)):
        locus_class = "I-F3_CAST"
    elif not has_adapt and "Cas8b" in roles_present:
        locus_class = "I-B_CAST"
    elif not has_adapt and "Cas8c" in roles_present:
        locus_class = "I-C_CAST"
    elif not has_adapt and "Cas12k" in roles_present:
        locus_class = "V-K_CAST"
    elif not has_adapt and has_any_cas:
        locus_class = "minimal_Cas_other"
    elif has_any_cas and has_adapt:
        locus_class = "Cas_with_Cas1_Cas2_or_Cas3"
    elif has_tns:
        locus_class = "Tn7_non_Cas"
    else:
        locus_class = "TniQ_only_or_truncated"

    return {
        "locus_class": locus_class,
        "has_cas1_cas2_cas3": has_adapt,
        "has_any_cas": has_any_cas,
        "has_tns_abc": has_tns,
        "has_crispr_array": len(crispr_arrays) > 0,
        "cas5f_at_plus1": "Cas5f" in role_by_offset.get(1, ()),
        "roles": sorted(roles_present),
        "order_summary": "; ".join(
            f"{off:+d}:{arch_by_offset[off]}" for off in sorted(arch_by_offset) if arch_by_offset[off] != "-"
        ),
    }


def _cpu_s() -> float:
    r0 = resource.getrusage(resource.RUSAGE_SELF)
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r0.ru_utime + r0.ru_stime + r1.ru_utime + r1.ru_stime


def _extract_tniq_hmm(pfam_gz: Path, dst: Path) -> dict:
    with gzip.open(pfam_gz, "rb") as fh, pyhmmer.plan7.HMMFile(fh) as hf:
        for hmm in hf:
            acc = (hmm.accession.decode() if isinstance(hmm.accession, bytes) else (hmm.accession or ""))
            name = (hmm.name.decode() if isinstance(hmm.name, bytes) else (hmm.name or ""))
            if acc.startswith("PF06527") or name == "TniQ":
                with open(dst, "wb") as out:
                    hmm.write(out)
                return {
                    "name": name,
                    "accession": acc,
                    "length": hmm.M,
                    "ga_seq": hmm.cutoffs.gathering1 if hmm.cutoffs.gathering_available() else None,
                    "ga_dom": hmm.cutoffs.gathering2 if hmm.cutoffs.gathering_available() else None,
                }
    raise RuntimeError(f"TniQ (PF06527) not found in {pfam_gz}")


def _stage_root(base: Path, stage: str) -> Path:
    cand = base / stage
    return cand if (cand / "l0seq").exists() or (cand / "runs").exists() else base


def _search_tniq_streaming(
    pilot_root: Path,
    pilot_stages: list[str],
    tniq_hmm_path: Path,
    cpus: int,
    block_proteins: int = 500_000,
) -> tuple[dict[bytes, float], int, int]:
    """Stream distinct pilot proteins in blocks through pyhmmer against TniQ.hmm at GA and Z=1e8."""
    with open(tniq_hmm_path, "rb") as fh, pyhmmer.plan7.HMMFile(fh) as hf:
        tniq_hmm = next(iter(hf))

    alphabet = pyhmmer.easel.Alphabet.amino()
    files = [l2cluster._protein_files(_stage_root(pilot_root, s), GENE_CALL_VERSION, [s])[0] for s in pilot_stages]
    hits_by_uid: dict[bytes, float] = {}
    n_distinct = total_residues = 0
    last = None
    batch_seqs: list[pyhmmer.easel.TextSequence] = []

    def run_block() -> None:
        if not batch_seqs:
            return
        dig = pyhmmer.easel.TextSequenceBlock(batch_seqs).digitize(alphabet)
        for top_hits in pyhmmer.hmmsearch([tniq_hmm], dig, cpus=cpus, Z=domains.Z_PROT, domZ=domains.DOM_Z,
                                          bit_cutoffs="gathering"):
            for hit in top_hits:
                if hit.included:
                    uid = bytes.fromhex(hit.name.decode() if isinstance(hit.name, bytes) else hit.name)
                    hits_by_uid[uid] = float(hit.score)
        batch_seqs.clear()

    for u, seq in heapq.merge(*(l2cluster._stage_proteins(f, seq=True) for f in files), key=lambda r: r[0]):
        if u == last:
            continue
        last = u
        n_distinct += 1
        clean = seq.replace("*", "")
        total_residues += len(clean)
        batch_seqs.append(pyhmmer.easel.TextSequence(name=u.hex().encode(), sequence=clean))
        if len(batch_seqs) >= block_proteins:
            run_block()
    run_block()
    return hits_by_uid, n_distinct, total_residues


def cmd_run(args: argparse.Namespace) -> int:
    pilot_root = Path(args.pilot_root).resolve()
    inputs_dir = Path(args.inputs_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    spill = out_dir / "spill"
    spill.mkdir(parents=True, exist_ok=True)

    stages_man = json.loads((inputs_dir / "stages_manifest.json").read_text())
    pilot_stages = [s["stage"] for s in stages_man]

    # 1. Extract TniQ.hmm
    tniq_hmm_path = out_dir / "TniQ.hmm"
    t0_w, t0_c = time.perf_counter(), _cpu_s()
    hmm_info = _extract_tniq_hmm(Path(args.pfam_hmm), tniq_hmm_path)
    t_extract_w, t_extract_c = time.perf_counter() - t0_w, _cpu_s() - t0_c

    # 2. Search all distinct pilot proteins for TniQ at GA (or reuse cached hits on MAG-annotation rerun)
    import pyarrow as pa
    hit_pq = out_dir / "tniq_hit_proteins.parquet"
    search_cache = out_dir / "tniq_search_cache.json"
    if hit_pq.is_file() and search_cache.is_file():
        sc = json.loads(search_cache.read_text())
        n_distinct_prots = int(sc["n_distinct_prots"])
        total_residues = int(sc["total_residues"])
        t_search_w = float(sc["t_search_w"])
        t_search_c = float(sc["t_search_c"])
        ht = pq.read_table(hit_pq)
        tniq_hits = dict(zip(ht.column("protein_uid").to_pylist(), (float(x) for x in ht.column("tniq_bits").to_pylist())))
    else:
        t1_w, t1_c = time.perf_counter(), _cpu_s()
        tniq_hits, n_distinct_prots, total_residues = _search_tniq_streaming(
            pilot_root, pilot_stages, tniq_hmm_path, cpus=args.threads
        )
        t_search_w, t_search_c = time.perf_counter() - t1_w, _cpu_s() - t1_c
        pq.write_table(
            pa.table({
                "protein_uid": pa.array(sorted(tniq_hits), type=pa.binary(16)),
                "tniq_bits": pa.array([tniq_hits[u] for u in sorted(tniq_hits)], type=pa.float32()),
            }),
            hit_pq,
            compression="zstd",
        )
        search_cache.write_text(json.dumps({
            "n_distinct_prots": n_distinct_prots,
            "total_residues": total_residues,
            "t_search_w": t_search_w,
            "t_search_c": t_search_c,
        }, indent=2) + "\n")

    # 3. Pull matching contig_tapes, contig_len, crispr arrays, and sample/habitat metadata via DuckDB
    t2_w, t2_c = time.perf_counter(), _cpu_s()

    p_occ_files, tape_files, ci_files, crispr_files, c_occ_files, l0_man_files = [], [], [], [], [], []
    cross_rows: list[tuple[bytes, str, str]] = []
    for st in pilot_stages:
        sr = _stage_root(pilot_root, st)
        p_occ_files.extend(sorted((sr / "l1" / "protein_occurrence").glob(f"gene_call_version=*/stage={st}/part-*.parquet")))
        tape_files.extend(sorted((sr / "l3" / "tapes").glob(f"tape_version=*/stage={st}/part-*.parquet")))
        tape_files.extend(sorted((sr / "l3" / "contig_tapes").glob(f"tape_version=*/stage={st}/part-*.parquet")))
        ci_files.extend(io.part_files(sr / "l0seq" / "contig_index" / f"stage={st}"))
        crispr_files.extend(sorted((sr / "l1" / "crispr").glob(f"producer=*/stage={st}/part-*.parquet")))
        c_occ_files.extend(io.part_files(sr / "l0seq" / "contig_occurrence" / f"stage={st}", "source=*/part-*.parquet"))
        l0_man_files.extend(io.part_files(sr / "l0" / "l0_manifest" / f"stage={st}"))
        ctsv = inputs_dir / st / "cross_stage_occurrences.tsv"
        if ctsv.is_file():
            with open(ctsv, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    cross_rows.append((bytes.fromhex(r["contig_uid"]), r["sample_id"], r["accession"]))

    cross_pq = out_dir / "cross_stage_occ.parquet"
    pq.write_table(
        pa.table({
            "contig_uid": pa.array([x[0] for x in cross_rows], type=pa.binary(16)),
            "sample_id": pa.array([x[1] for x in cross_rows], type=pa.string()),
            "accession": pa.array([x[2] for x in cross_rows], type=pa.string()),
        }),
        cross_pq,
        compression="zstd",
    )

    con = duckdb.connect(config={"threads": min(args.threads, 28), "memory_limit": "32GB", "temp_directory": str(spill)})
    try:
        po_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in p_occ_files)
        tp_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in tape_files)
        ci_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in ci_files)
        co_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in c_occ_files)
        lm_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in l0_man_files)

        tapes_raw = con.execute(
            f"""
            WITH hit_contigs AS (
                SELECT DISTINCT po.contig_uid
                FROM read_parquet([{po_lits}], hive_partitioning=false) po
                JOIN read_parquet('{hit_pq}') hp USING (protein_uid)
            )
            SELECT
                t.contig_uid,
                t.contig_len::BIGINT AS contig_len,
                t.feature_uids,
                t.protein_uids,
                t.start_bps,
                t.end_bps,
                t.strands,
                t.feature_kinds
            FROM read_parquet([{tp_lits}], hive_partitioning=false) t
            JOIN hit_contigs hc USING (contig_uid)
            """
        ).fetchall()

        # Also map each hit contig to its owning stage so we can run MinCED on any hit contig whose stage skipped crispr
        hit_stage_rows = con.execute(
            f"""
            WITH hit_contigs AS (
                SELECT DISTINCT po.contig_uid
                FROM read_parquet([{po_lits}], hive_partitioning=false) po
                JOIN read_parquet('{hit_pq}') hp USING (protein_uid)
            )
            SELECT ci.contig_uid, regexp_extract(ci.filename, 'stage=([^/]+)', 1) AS stage
            FROM read_parquet([{ci_lits}], hive_partitioning=false, filename=true) ci
            JOIN hit_contigs hc USING (contig_uid)
            WHERE ci.is_default_variant
            """
        ).fetchall()

        by_owner_stage: dict[str, list[bytes]] = defaultdict(list)
        for cu, st_own in hit_stage_rows:
            by_owner_stage[st_own].append(cu)

        hit_crispr_out = out_dir / "hit_crispr_out"
        hit_crispr_tmp = out_dir / "hit_crispr_tmp"
        if hit_crispr_out.exists():
            shutil.rmtree(hit_crispr_out)
        hit_crispr_tmp.mkdir(parents=True, exist_ok=True)
        for st_own, cus in sorted(by_owner_stage.items()):
            sr_own = _stage_root(pilot_root, st_own)
            existing_st_cr = list((sr_own / "l1" / "crispr").glob(f"producer={crispr.producer()}/stage={st_own}/part-*.parquet"))
            if not existing_st_cr and cus:
                crispr.run(
                    sr_own,
                    st_own,
                    cpus=min(args.threads, 4),
                    genomes=[
                        GenomeInput(contig_uids=(cu,), genetic_code=11, kind="prokaryote", circular=(), assembly_quality="contig")
                        for cu in sorted(set(cus))
                    ],
                    tmp=hit_crispr_tmp,
                    out=hit_crispr_out,
                    seq_stages=[st_own] if sr_own != pilot_root else pilot_stages,
                )
        shutil.rmtree(hit_crispr_tmp, ignore_errors=True)
        crispr_files = sorted(set(crispr_files) | set((hit_crispr_out / "l1" / "crispr").glob("producer=*/stage=*/part-*.parquet")))

        crispr_raw = []
        if crispr_files:
            cr_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in crispr_files)
            crispr_raw = con.execute(
                f"""
                WITH hit_contigs AS (
                    SELECT DISTINCT po.contig_uid
                    FROM read_parquet([{po_lits}], hive_partitioning=false) po
                    JOIN read_parquet('{hit_pq}') hp USING (protein_uid)
                )
                SELECT
                    c.contig_uid,
                    min(c.start)::BIGINT AS start,
                    max(c."end")::BIGINT AS "end",
                    sum(CASE WHEN c.element = 'repeat' THEN 1 ELSE 0 END)::INTEGER AS n_repeats
                FROM read_parquet([{cr_lits}], hive_partitioning=false) c
                JOIN hit_contigs hc USING (contig_uid)
                GROUP BY c.contig_uid, c.feature_uid
                """
            ).fetchall()

        # Sample membership of hit contigs
        contig_sample_rows = con.execute(
            f"""
            WITH hit_contigs AS (
                SELECT DISTINCT po.contig_uid
                FROM read_parquet([{po_lits}], hive_partitioning=false) po
                JOIN read_parquet('{hit_pq}') hp USING (protein_uid)
            ),
            obj_sample AS (
                SELECT DISTINCT
                    object_sha256 AS l0_object_sha256,
                    regexp_replace(regexp_extract(object_path, '([^/]+)\\.fa(\\.gz)?$', 1), '\\.fa(\\.gz)?$', '') AS sample_id
                FROM read_parquet([{lm_lits}], hive_partitioning=false)
            ),
            all_occ AS (
                SELECT co.contig_uid, os.sample_id, co.accession
                FROM read_parquet([{co_lits}], hive_partitioning=false) co
                JOIN obj_sample os USING (l0_object_sha256)
                UNION
                SELECT contig_uid, sample_id, accession FROM read_parquet('{cross_pq}')
            )
            SELECT DISTINCT ao.contig_uid, ao.sample_id, ao.accession
            FROM all_occ ao
            JOIN hit_contigs hc USING (contig_uid)
            ORDER BY ao.contig_uid, ao.sample_id
            """
        ).fetchall()

        tapes_dicts = [
            {
                "contig_uid": r[0], "contig_len": int(r[1]), "feature_uids": list(r[2]),
                "protein_uids": list(r[3]), "starts": list(r[4]), "ends": list(r[5]),
                "strands": list(r[6]), "kinds": list(r[7]),
            }
            for r in tapes_raw
        ]
        crispr_dicts = [
            {"contig_uid": r[0], "start": int(r[1]), "end": int(r[2]), "n_repeats": int(r[3])}
            for r in crispr_raw
        ]

        windows, req_uids = extract_hit_windows(tapes_dicts, set(tniq_hits), crispr_rows=crispr_dicts, window=10)

        # Fetch sequences of ONLY the requested window proteins
        req_pq = out_dir / "window_protein_uids.parquet"
        pq.write_table(pa.table({"protein_uid": pa.array(req_uids, type=pa.binary(16))}), req_pq, compression="zstd")
        ps_files = []
        for st in pilot_stages:
            sr = _stage_root(pilot_root, st)
            ps_files.extend(sorted((sr / "l0seq" / "protein_seq").glob(f"gene_call_version=*/stage={st}/protein_prefix=*/part-*.parquet")))
        ps_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in ps_files)
        seq_rows = con.execute(
            f"""
            SELECT DISTINCT ps.protein_uid, ps.seq
            FROM read_parquet([{ps_lits}], hive_partitioning=false) ps
            JOIN read_parquet('{req_pq}') rq USING (protein_uid)
            """
        ).fetchall()
    finally:
        con.close()
        shutil.rmtree(spill, ignore_errors=True)

    window_proteins = {u: s for u, s in seq_rows}
    window_residues = sum(len(s) for s in window_proteins.values())
    t_win_w, t_win_c = time.perf_counter() - t2_w, _cpu_s() - t2_c

    # 4. Run Pfam 38.2 round scan (`biobox.l1.domains.scan`) ONLY on the +-10-gene window proteins!
    #    Shard into <= 8,000,000 residues per child process (~57 CPU-min << 280 CPU-min cap).
    round_pq = out_dir / "round_domains.parquet"
    pfam_cache = out_dir / "pfam_round_cache.json"
    dom_schema = domains.table("l1_domains").schema
    if round_pq.is_file() and pfam_cache.is_file():
        pcache = json.loads(pfam_cache.read_text())
        t_pfam_w = float(pcache["t_pfam_w"])
        t_pfam_c = float(pcache["t_pfam_c"])
        domain_rows = pq.read_table(round_pq).to_pylist()
    elif not window_proteins:
        t_pfam_w, t_pfam_c = 0.0, 0.0
        domain_rows = []
        pq.write_table(pa.Table.from_pylist([], schema=dom_schema), round_pq, compression="zstd")
        pfam_cache.write_text(json.dumps({"t_pfam_w": 0.0, "t_pfam_c": 0.0, "shards": 0}, indent=2) + "\n")
    elif window_residues <= 8_000_000:
        t3_w, t3_c = time.perf_counter(), _cpu_s()
        domain_rows = domains.scan(window_proteins, Path(args.pfam_hmm), "gathering", cpus=args.threads)
        t_pfam_w, t_pfam_c = time.perf_counter() - t3_w, _cpu_s() - t3_c
        pq.write_table(pa.Table.from_pylist(domain_rows, schema=dom_schema), round_pq, compression="zstd")
        pfam_cache.write_text(json.dumps({"t_pfam_w": t_pfam_w, "t_pfam_c": t_pfam_c, "shards": 1}, indent=2) + "\n")
    else:
        shards_dir = out_dir / "pfam_shards"
        shards_dir.mkdir(parents=True, exist_ok=True)
        uids_sorted = sorted(window_proteins)
        shard_inputs: list[Path] = []
        cur_u, cur_s, cur_res = [], [], 0
        for u in uids_sorted:
            s = window_proteins[u]
            if cur_u and cur_res + len(s) > 8_000_000:
                sp = shards_dir / f"in_{len(shard_inputs):04d}.parquet"
                pq.write_table(pa.table({"protein_uid": pa.array(cur_u, type=pa.binary(16)), "seq": pa.array(cur_s, type=pa.string())}), sp, compression="zstd")
                shard_inputs.append(sp)
                cur_u, cur_s, cur_res = [], [], 0
            cur_u.append(u)
            cur_s.append(s)
            cur_res += len(s)
        if cur_u:
            sp = shards_dir / f"in_{len(shard_inputs):04d}.parquet"
            pq.write_table(pa.table({"protein_uid": pa.array(cur_u, type=pa.binary(16)), "seq": pa.array(cur_s, type=pa.string())}), sp, compression="zstd")
            shard_inputs.append(sp)

        t3_w = time.perf_counter()
        t_pfam_c = 0.0
        shard_tables: list[pa.Table] = []
        for idx_s, sp in enumerate(shard_inputs):
            out_sp = shards_dir / f"out_{idx_s:04d}.parquet"
            proc = subprocess.Popen([
                sys.executable, str(Path(__file__).resolve()), "scan-shard",
                "--in-parquet", str(sp),
                "--out-parquet", str(out_sp),
                "--pfam-hmm", str(args.pfam_hmm),
                "--threads", str(args.threads),
            ])
            _, st_code, ru = os.wait4(proc.pid, 0)
            rc = os.waitstatus_to_exitcode(st_code)
            if rc != 0:
                raise RuntimeError(f"scan-shard {idx_s} exited with code {rc}")
            t_pfam_c += ru.ru_utime + ru.ru_stime
            shard_tables.append(pq.read_table(out_sp))
        t_pfam_w = time.perf_counter() - t3_w
        combined_tbl = pa.concat_tables(shard_tables) if shard_tables else pa.Table.from_pylist([], schema=dom_schema)
        pq.write_table(combined_tbl, round_pq, compression="zstd")
        domain_rows = combined_tbl.to_pylist()
        shutil.rmtree(shards_dir, ignore_errors=True)
        pfam_cache.write_text(json.dumps({"t_pfam_w": t_pfam_w, "t_pfam_c": t_pfam_c, "shards": len(shard_inputs)}, indent=2) + "\n")

    pfam_by_prot: dict[bytes, list[tuple[str, str, float]]] = defaultdict(list)
    for r in domain_rows:
        pfam_by_prot[r["protein_uid"]].append((r["model_acc"], r["model_name"], float(r["dom_bitscore"])))

    # Map sample -> habitat
    hab_by_sample = {}
    mag_count_by_sample = {}
    with open(args.selection, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            hab_by_sample[r["sample_id"]] = (r["habitat"], r["habitat2"], r["study"])
            mag_count_by_sample[r["sample_id"]] = int(r.get("mag_count") or 0)

    samples_of_cu: dict[bytes, list[tuple[str, str]]] = defaultdict(list)
    for cu, sid, acc in contig_sample_rows:
        samples_of_cu[cu].append((sid, acc))

    # Write distinct sample_ids of hit contigs so MAG membership can be checked cleanly if desired
    hit_sids = sorted({sid for _, sid, _ in contig_sample_rows if sid})
    (out_dir / "hit_samples.txt").write_text("".join(f"{s}\n" for s in hit_sids))

    mag_by_cu_hex: dict[str, str] = {}
    if args.mag_uids_json and Path(args.mag_uids_json).is_file():
        mag_by_cu_hex = json.loads(Path(args.mag_uids_json).read_text()).get("contig_uid_to_mag", {})

    locus_rows: list[dict] = []
    class_counts = Counter()
    hab_class_counts: dict[str, Counter] = defaultdict(Counter)
    len_band_stats: dict[str, list[bool]] = defaultdict(list)
    binned_loci = 0

    for w in windows:
        cu = w["contig_uid"]
        cls = classify_tniq_locus(w["slots"], pfam_by_prot, w["crispr_arrays"])
        s_list = samples_of_cu.get(cu, [("", "")])
        first_sid, first_acc = s_list[0]
        hab1, hab2, study = hab_by_sample.get(first_sid, ("unlabelled", "unlabelled", ""))
        clen = w["contig_len"]
        band = (
            "lt_1.5kb" if clen < 1500
            else "1.5_to_5kb" if clen < 5000
            else "5_to_10kb" if clen < 10000
            else "ge_10kb"
        )
        both_intact = not (w["truncated_5p"] or w["truncated_3p"])
        len_band_stats[band].append(both_intact)
        class_counts[cls["locus_class"]] += 1
        hab_class_counts[hab1][cls["locus_class"]] += 1

        mag_id = mag_by_cu_hex.get(cu.hex(), "")
        if mag_id:
            binned_loci += 1

        locus_rows.append({
            "contig_uid": cu.hex(),
            "contig_len": clen,
            "length_band": band,
            "sample_id": first_sid,
            "n_pilot_samples": len({s for s, _ in s_list if s}),
            "contig_accession": first_acc,
            "spire_mag_id": mag_id,
            "in_spire_mag": bool(mag_id),
            "study": study,
            "habitat": hab1,
            "habitat2": hab2,
            "anchor_feature_uid": w["anchor_feature_uid"].hex(),
            "anchor_protein_uid": w["anchor_protein_uid"].hex(),
            "tniq_bits": round(tniq_hits[w["anchor_protein_uid"]], 2),
            "anchor_start": w["anchor_start"],
            "anchor_end": w["anchor_end"],
            "anchor_strand": w["anchor_strand"],
            "window_slots": w["window_slots"],
            "truncated_5p": w["truncated_5p"],
            "truncated_3p": w["truncated_3p"],
            "n_crispr_arrays": len(w["crispr_arrays"]),
            "locus_class": cls["locus_class"],
            "has_cas1_cas2_cas3": cls["has_cas1_cas2_cas3"],
            "cas5f_at_plus1": cls["cas5f_at_plus1"],
            "roles": ",".join(cls["roles"]),
            "order_summary": cls["order_summary"],
        })

    if locus_rows:
        cols = list(locus_rows[0].keys())
        with open(out_dir / "tniq_loci.tsv", "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            wr.writeheader()
            wr.writerows(locus_rows)

    est_full_pfam_cpu_min = total_residues * 7.18 / 1e6
    actual_round_cpu_min = (t_search_c + t_win_c + t_pfam_c) / 60.0

    summary = {
        "hmm": hmm_info,
        "pilot_distinct_proteins_searched": n_distinct_prots,
        "pilot_total_residues": total_residues,
        "tniq_hit_distinct_proteins": len(tniq_hits),
        "tniq_hit_loci": len(windows),
        "tniq_hit_loci_in_spire_mags": binned_loci if args.mag_uids_json else None,
        "tniq_hit_loci_unbinned": (len(windows) - binned_loci) if args.mag_uids_json else None,
        "tniq_hit_contigs": len({w["contig_uid"] for w in windows}),
        "window_distinct_proteins_scanned_with_pfam": len(window_proteins),
        "window_residues_scanned_with_pfam": window_residues,
        "window_proteins_with_pfam_domain": len(pfam_by_prot),
        "round_domain_rows": len(domain_rows),
        "protein_reduction_factor": round(n_distinct_prots / max(1, len(window_proteins)), 1),
        "residue_reduction_factor": round(total_residues / max(1, window_residues), 1),
        "timing": {
            "hmm_extract_wall_s": round(t_extract_w, 2),
            "hmm_extract_cpu_min": round(t_extract_c / 60.0, 3),
            "tniq_search_wall_s": round(t_search_w, 2),
            "tniq_search_cpu_min": round(t_search_c / 60.0, 3),
            "tniq_search_cpu_min_per_10m_proteins": round((t_search_c / 60.0) / (n_distinct_prots / 1e7), 3) if n_distinct_prots else 0.0,
            "window_extract_wall_s": round(t_win_w, 2),
            "window_extract_cpu_min": round(t_win_c / 60.0, 3),
            "pfam_round_wall_s": round(t_pfam_w, 2),
            "pfam_round_cpu_min": round(t_pfam_c / 60.0, 3),
            "total_search_plus_round_cpu_min": round(actual_round_cpu_min, 3),
            "estimated_full_pfam_cpu_min_at_7_18": round(est_full_pfam_cpu_min, 1),
            "cpu_speedup_vs_full_pfam": round(est_full_pfam_cpu_min / max(0.001, actual_round_cpu_min), 1),
        },
        "locus_classes": dict(sorted(class_counts.items())),
        "by_habitat": {h: dict(sorted(c.items())) for h, c in sorted(hab_class_counts.items())},
        "intact_21_slot_window_by_contig_length": {
            b: {
                "loci": len(vals),
                "intact_both_sides": sum(vals),
                "truncated_share": round(1.0 - (sum(vals) / len(vals)), 4) if vals else 0.0,
            }
            for b, vals in sorted(len_band_stats.items())
        },
    }
    (out_dir / "tniq_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))
    return 0


def cmd_scan_shard(args: argparse.Namespace) -> int:
    import pyarrow as pa
    tbl = pq.read_table(args.in_parquet)
    prots = dict(zip(tbl.column("protein_uid").to_pylist(), tbl.column("seq").to_pylist()))
    rows = domains.scan(prots, Path(args.pfam_hmm), "gathering", cpus=args.threads) if prots else []
    tmp_out = Path(args.out_parquet).with_suffix(".parquet.tmp")
    pq.write_table(pa.Table.from_pylist(rows, schema=domains.table("l1_domains").schema), tmp_out, compression="zstd")
    tmp_out.replace(Path(args.out_parquet))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--pilot-root", required=True)
    p.add_argument("--inputs-dir", required=True)
    p.add_argument("--selection", required=True)
    p.add_argument("--pfam-hmm", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--mag-uids-json", default=None)
    p.add_argument("--threads", type=int, default=16)

    ps = sub.add_parser("scan-shard")
    ps.add_argument("--in-parquet", required=True)
    ps.add_argument("--out-parquet", required=True)
    ps.add_argument("--pfam-hmm", required=True)
    ps.add_argument("--threads", type=int, default=16)

    args = ap.parse_args(argv)
    if args.cmd == "run":
        return cmd_run(args)
    if args.cmd == "scan-shard":
        return cmd_scan_shard(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())

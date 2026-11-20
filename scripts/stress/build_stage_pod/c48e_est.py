"""Card C-48e part 1: the staged tranche measured against its estimate, its overlap, the cost of building it, the
search pack's size, and what the tranche leaves for later. Reads only; writes into `--out`.

    python c48e_est.py --root /tmp/if2453/r10/root --inputs /tmp/if2453/s2/inputs/SPX1 \
        --uids /tmp/if2453/c48e1/work/contigs_ge5000.parquet --selection SEL.tsv --select-dir SELDIR \
        --fetch-log DL/fetch_log.jsonl --pilot-depth1 T --pilot-depth1-unfiltered T --pilot-pack J \
        --profile profiles/wilkcomp02.yaml --other-inputs /tmp/if2453/s2/inputs --tmp TMP --out OUT \
        [--max-contigs 3000000 --max-bp 40000000000] [--sample-bits 8]

1. Measured contigs and bases at >= the length per sample (the inputs' `samples.tsv`) against the selection's
   estimate, by stratum, habitat and in all; the caps checked on the measured values; if over, the samples to drop
   from the end of the order (by order rank) until both hold, on the per-sample sums (an upper bound of the distinct).
2. Overlap by `contig_uid`: the tranche's distinct long contigs against each stage's `contig_index` (default
   variants) in the root, and against the `genomes.json` of each other staged inputs directory.
3. The cost of `build-stage --depth called --min-contig-len L` from the code's own estimates (`stage.ingest_cpu_min`,
   `stage.pass1_bytes`'s rule, `genecall.cut_plan` over the stage's one-contig genomes in genome_uid order,
   `stage.crispr_cpu_min`) and from the pilot's measured rates (card C-42: ingest, gene calling, IGR, tapes,
   CRISPR before card B-55), and the GB added to the root by a least-squares fit of the pilot's measured stage bytes
   on stored bases and contigs.
4. A measured sample: the long contigs whose contig_uid's first byte is under 2**(8 - sample_bits) (1/256 at 8),
   gene-called one contig at a time as the step calls them (`genecall._call`), their distinct proteins, IGRs
   (gaps of at least `pack.IGR_MIN_NT` between genes) and dark ORFs (`orfs.dark_orfs`); with the pilot's measured
   pack bytes per residue and per IGR (stage `SPA1`), the search pack with and without `dark_orf`.
5. The remainder (DEC-59 (3)): every sample of the habitats not staged, in the order a next tranche takes them
   (with the reason), the bases under the length in the staged assemblies by stratum, and the other habitats.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io as _io
import json
import resource
import shutil
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import fasta, ids, io, stage
from biobox.l0seq import store
from biobox.l0seq.inputs import GenomeInput
from biobox.l1 import genecall
from biobox.profile import load
from biobox.search import orfs, pack


def _sql(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def _tsv(path: Path) -> list[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict], cols: list[str] | None = None, gz: bool = False) -> None:
    cols = cols or (list(rows[0].keys()) if rows else [])
    buf = _io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})
    data = buf.getvalue().encode()
    if gz:
        with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as fh:
            fh.write(data)
    else:
        path.write_bytes(data)


def _cpu() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def _ratio(a, b):
    return round(a / b, 4) if b else None


# --- 1. measured against the estimate -----------------------------------------------------------------------------

def measured_vs_estimate(sel: list[dict], samples: list[dict], max_n: int, max_bp: int) -> dict:
    meas = {r["sample_id"]: r for r in samples}
    per = []
    for r in sel:
        m = meas.get(r["sample_id"])
        per.append({"order_rank": int(r["order_rank"]), "sample_id": r["sample_id"], "habitat": r["habitat"],
                    "habitat2": r["habitat2"], "est_source": r["est_source"], "assembly_bp": int(r["assembly_bp"]),
                    "est_n_ge_5000": int(r["est_n_ge_5000"]), "est_bp_ge_5000": int(r["est_bp_ge_5000"]),
                    "staged": m is not None,
                    "meas_bp": int(m["bp"]) if m else None,
                    "meas_n_ge_5000": int(m["records_at_least"]) if m else None,
                    "meas_bp_ge_5000": int(m["bp_at_least"]) if m else None})
    per.sort(key=lambda r: r["order_rank"])

    def group(key):
        g = defaultdict(lambda: defaultdict(int))
        for r in per:
            if not r["staged"]:
                continue
            for k in ([r[key]] if key else ["all"]):
                d = g[k]
                d["samples"] += 1
                d["assembly_bp"] += r["assembly_bp"]
                d["meas_bp"] += r["meas_bp"]
                d["est_n_ge_5000"] += r["est_n_ge_5000"]
                d["meas_n_ge_5000"] += r["meas_n_ge_5000"]
                d["est_bp_ge_5000"] += r["est_bp_ge_5000"]
                d["meas_bp_ge_5000"] += r["meas_bp_ge_5000"]
        out = []
        for k in sorted(g):
            d = g[k]
            out.append({"level": key or "all", "group": k, **d,
                        "ratio_n_meas_over_est": _ratio(d["meas_n_ge_5000"], d["est_n_ge_5000"]),
                        "ratio_bp_meas_over_est": _ratio(d["meas_bp_ge_5000"], d["est_bp_ge_5000"]),
                        "ratio_assembly_bp": _ratio(d["meas_bp"], d["assembly_bp"])})
        return out

    rows = group("habitat2") + group("habitat") + group(None)
    # caps on the per-sample sums, in order; samples to drop from the end until both hold
    staged = [r for r in per if r["staged"]]
    n = sum(r["meas_n_ge_5000"] for r in staged)
    bp = sum(r["meas_bp_ge_5000"] for r in staged)
    drop = []
    while staged and (n > max_n or bp > max_bp):
        r = staged.pop()
        drop.append(r["sample_id"])
        n -= r["meas_n_ge_5000"]
        bp -= r["meas_bp_ge_5000"]
    return {"per_sample": per, "groups": rows, "drop_from_end": drop,
            "kept_sum_n_ge_5000": n, "kept_sum_bp_ge_5000": bp}


# --- 2. overlap -----------------------------------------------------------------------------------------------------

def overlap(con, root: Path, uids: Path, other_inputs: Path, own_stage: str) -> dict:
    con.execute(f"CREATE OR REPLACE TABLE ours AS SELECT contig_uid, contig_len FROM read_parquet({_sql([uids])})")
    out = {"root_stages": {}, "inputs": {}}
    held_q = []
    for st in store.stages_on_disk(root):
        files = io.part_files(root / store.INDEX_DIR / f"stage={st}")
        if not files:
            continue
        q = (f"SELECT DISTINCT contig_uid FROM read_parquet({_sql(files)}, hive_partitioning = false) "
             "WHERE is_default_variant")
        held_q.append(q)
        (c, b), = con.execute(f"SELECT count(*), coalesce(sum(contig_len), 0) FROM ours "
                              f"SEMI JOIN ({q}) h USING (contig_uid)").fetchall()
        out["root_stages"][st] = {"contigs": int(c), "bp": int(b), "index_files": len(files)}
    if held_q:
        con.execute("CREATE OR REPLACE TABLE held AS SELECT DISTINCT contig_uid FROM (" + " UNION ALL ".join(held_q) + ")"
                    " SEMI JOIN ours USING (contig_uid)")
    else:
        con.execute("CREATE OR REPLACE TABLE held (contig_uid BLOB)")
    (c, b), = con.execute("SELECT count(*), coalesce(sum(contig_len), 0) FROM ours SEMI JOIN held USING (contig_uid)").fetchall()
    out["root_any"] = {"contigs": int(c), "bp": int(b)}
    for d in sorted(p for p in Path(other_inputs).iterdir() if p.is_dir()):
        if d.name == own_stage or not (d / "genomes.json").is_file():
            continue
        (c, b), = con.execute(
            f"WITH g AS (SELECT DISTINCT unnest(contig_uids) AS h FROM read_json({_sql([d / 'genomes.json'])}, "
            "format = 'array', columns = {'contig_uids': 'VARCHAR[]'})) "
            "SELECT count(*), coalesce(sum(o.contig_len), 0) FROM ours o JOIN g ON g.h = lower(hex(o.contig_uid))"
        ).fetchall()
        out["inputs"][d.name] = {"contigs": int(c), "bp": int(b)}
    (c, b), = con.execute("SELECT count(*), coalesce(sum(contig_len), 0) FROM ours").fetchall()
    out["ours"] = {"contigs": int(c), "bp": int(b)}
    return out


# --- 3. costs -------------------------------------------------------------------------------------------------------

def genecall_plan(con, work: Path, batch_cpu_min: float) -> dict:
    con.execute(f"COPY (SELECT o.contig_uid, o.contig_len, (h.contig_uid IS NOT NULL) AS held FROM ours o "
                f"LEFT JOIN held h USING (contig_uid) ORDER BY o.contig_uid) TO '{work / 'flags.parquet'}' (FORMAT PARQUET)")
    schema = pa.schema([("genome_uid", pa.binary(16)), ("cost_us", pa.int64()), ("bases", pa.int64()),
                        ("genes", pa.float64()), ("called", pa.bool_())])
    w = pq.ParquetWriter(work / "items.parquet", schema)
    modes = defaultdict(lambda: [0, 0])
    for b in pq.ParquetFile(work / "flags.parquet").iter_batches(batch_size=122_880):
        cu, ln, held = (b.column(i).to_pylist() for i in range(3))
        g, c, n, ge, ca = [], [], [], [], []
        for u, L, h in zip(cu, ln, held):
            g.append(ids.genome_uid([u]))
            n.append(L)
            ca.append(not h)
            if h:
                c.append(0)
                ge.append(0.0)
            else:
                c.append(genecall.genome_cost_us("prokaryote", [(L, False)]))
                ge.append(genecall.genome_genes("prokaryote", L))
                m = modes[genecall.mode_of("prokaryote", L)]
                m[0] += 1
                m[1] += L
        w.write_batch(pa.RecordBatch.from_arrays([pa.array(g, pa.binary(16)), pa.array(c, pa.int64()),
                                                  pa.array(n, pa.int64()), pa.array(ge, pa.float64()),
                                                  pa.array(ca, pa.bool_())], schema=schema))
    w.close()
    con.execute(f"COPY (SELECT * FROM read_parquet('{work / 'items.parquet'}') ORDER BY genome_uid) "
                f"TO '{work / 'items_sorted.parquet'}' (FORMAT PARQUET)")

    def items():
        for b in pq.ParquetFile(work / "items_sorted.parquet").iter_batches(batch_size=122_880):
            yield from zip(*(b.column(i).to_pylist() for i in range(5)))

    plan = genecall.cut_plan(items, batch_cpu_min)
    bt = plan.pop("batches")
    plan["batches_n"] = len(bt)
    plan["largest_batch_cpu_min"] = max(x["estimate_cpu_min"] for x in bt)
    plan["batches_cpu_min_sum"] = round(sum(x["estimate_cpu_min"] for x in bt), 3)
    plan["modes"] = {k: {"genomes": v[0], "bases": v[1]} for k, v in sorted(modes.items())}
    return plan


def fit_bytes(depth1: list[dict]) -> dict:
    """Least squares of the pilot's measured stage bytes on (stored bases, stored contigs), no intercept."""
    X = np.array([[float(r["stored_bp"]), float(r["owned_contigs"])] for r in depth1])
    out = {}
    for col in ("bytes_depth1_total", "bytes_depth0_stored", "bytes_without_nt_and_igr"):
        y = np.array([float(r[col]) for r in depth1])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coef
        out[col] = {"bytes_per_base": float(coef[0]), "bytes_per_contig": float(coef[1]),
                    "max_abs_rel_err": float(np.max(np.abs(pred - y) / y))}
    return out


def pilot_rates(depth1: list[dict]) -> dict:
    bp = sum(float(r["stored_bp"]) for r in depth1)
    inp = sum(float(r["input_bp"]) for r in depth1)
    contigs = sum(float(r["owned_contigs"]) for r in depth1)
    cr = [r for r in depth1 if r["crispr_ran"] == "True"]
    crbp = sum(float(r["stored_bp"]) for r in cr)
    return {
        "stages": [r["stage"] for r in depth1],
        "stored_gbp": bp / 1e9,
        "ingest_cpu_min_per_gbp": sum(float(r["cpu_min_ingest"]) for r in depth1) / (inp / 1e9),
        "genecall_cpu_min_per_gbp": sum(float(r["cpu_min_genecall"]) for r in depth1) / (bp / 1e9),
        "igr_cpu_min_per_gbp": sum(float(r["cpu_min_igr"]) for r in depth1) / (bp / 1e9),
        "tapes_cpu_min_per_gbp": sum(float(r["cpu_min_tapes"]) for r in depth1) / (bp / 1e9),
        "inputs_cpu_min_per_m_genomes": sum(float(r["cpu_min_inputs"]) for r in depth1) / (contigs / 1e6),
        "crispr_cpu_min_per_gbp_before_b55": (sum(float(r["cpu_min_crispr"]) for r in cr) / (crbp / 1e9)) if crbp else None,
        "gb_per_gbp_called_total": sum(float(r["bytes_depth1_total"]) for r in depth1) / bp,
        "genes_per_gbp": sum(float(r["genes"]) for r in depth1) / (bp / 1e9),
    }


# --- 4. the measured sample -----------------------------------------------------------------------------------------

def measured_sample(objects: list[dict], L: int, bits: int, held: set[bytes]) -> dict:
    lim = 1 << (8 - bits)
    seqs: dict[bytes, bytes] = {}
    for o in objects:
        for fr in fasta.iter_records(Path(o["path"])):
            if len(fr.seq) < L:
                continue
            cu = ids.contig_uid(fr.seq)
            if cu[0] < lim and cu not in held:
                seqs.setdefault(cu, fr.seq.upper())
    s = defaultdict(float)
    prots: dict[bytes, int] = {}
    t_call = t_dark = 0.0
    for cu in sorted(seqs):
        seq = seqs[cu]
        g = GenomeInput((cu,), 11, kind="prokaryote", assembly_quality="contig")
        t0 = _cpu()
        calls = genecall._call(g, {cu: seq.decode("ascii")})
        t_call += _cpu() - t0
        s["contigs"] += 1
        s["bp"] += len(seq)
        s["genes"] += len(calls)
        s["fragments"] += sum(1 for c in calls if c["dist_to_contig_end_bp"] == 0)
        for c in calls:
            prots.setdefault(ids.protein_uid(c["aa"]), len(c["aa"]))
        spans = sorted((c["start"], c["end"]) for c in calls)
        end = None
        for a, b in spans:
            if end is not None and a - end >= pack.IGR_MIN_NT:
                s["igrs"] += 1
                s["igr_nt"] += a - end
            end = b if end is None else max(end, b)
        t0 = _cpu()
        dark = orfs.dark_orfs(seq, False, [(c["start"], c["end"], c["strand"]) for c in calls], 11)
        t_dark += _cpu() - t0
        s["dark_orfs"] += len(dark)
        s["dark_aa"] += sum(len(o.aa) for o, _, _ in dark)
        s["dark_at_end"] += sum(1 for _, _, e in dark if e)
    s["distinct_proteins"] = len(prots)
    s["distinct_aa"] = sum(prots.values())
    s["call_cpu_s"] = round(t_call, 3)
    s["dark_cpu_s"] = round(t_dark, 3)
    s["call_s_per_mbp"] = round(t_call / (s["bp"] / 1e6), 4) if s["bp"] else None
    s["sample_fraction"] = f"contig_uid first byte < {lim} (1/{256 // lim})"
    return dict(s)


# --- main -----------------------------------------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    for k in ("root", "inputs", "uids", "selection", "select-dir", "fetch-log", "pilot-depth1", "pilot-depth1-unfiltered",
              "pilot-pack", "profile", "other-inputs", "tmp", "out"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--max-contigs", type=int, default=3_000_000)
    ap.add_argument("--max-bp", type=int, default=40_000_000_000)
    ap.add_argument("--sample-bits", type=int, default=8)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    inp = Path(a.inputs)
    rep = json.loads((inp / "input_report.json").read_text())
    L = int(rep["min_contig_len"])
    stage_name = rep["stage"]
    samples = _tsv(inp / "samples.tsv")
    sel = _tsv(Path(a.selection))
    prof = load(Path(a.profile))
    res: dict = {"stage": stage_name, "min_contig_len": L, "caps": {"max_contigs": a.max_contigs, "max_bp": a.max_bp}}

    # 1
    mv = measured_vs_estimate(sel, [r for r in samples], a.max_contigs, a.max_bp)
    _write_tsv(out / "per_sample_measured.tsv", mv["per_sample"])
    _write_tsv(out / "measured_vs_estimate.tsv", mv["groups"])
    res["measured"] = {k: rep[k] for k in ("samples_staged", "records", "bp", "records_below", "bp_below",
                                            "records_at_least", "bp_at_least", "distinct_contigs", "distinct_bp",
                                            "contigs_in_2plus_samples", "distinct_contigs_ge_100000",
                                            "distinct_bp_ge_100000", "gz_bytes", "fasta_bytes")}
    res["measured_vs_estimate_all"] = [g for g in mv["groups"] if g["level"] == "all"][0]
    res["caps_hold_on_distinct"] = rep["distinct_contigs"] <= a.max_contigs and rep["distinct_bp"] <= a.max_bp
    res["drop_from_end"] = mv["drop_from_end"]
    print(json.dumps({"measured": res["measured"], "vs_est": res["measured_vs_estimate_all"],
                      "drop": res["drop_from_end"]}), flush=True)

    # 2
    tmp = Path(tempfile.mkdtemp(prefix="c48e-est-", dir=a.tmp))
    con = duckdb.connect(config={"threads": 2, "memory_limit": "6GB", "temp_directory": str(tmp / "duck"),
                                 "preserve_insertion_order": False})
    try:
        t0 = time.time()
        ov = overlap(con, Path(a.root), Path(a.uids), Path(a.other_inputs), stage_name)
        ov["read_utc_seconds"] = round(time.time() - t0, 1)
        res["overlap"] = ov
        print(json.dumps({"overlap": ov}), flush=True)
        held = {bytes(r[0]) for r in con.execute("SELECT contig_uid FROM held").fetchall()}

        # 3
        new_bp = ov["ours"]["bp"] - ov["root_any"]["bp"]
        new_contigs = ov["ours"]["contigs"] - ov["root_any"]["contigs"]
        batch_cpu = stage.SHARD_MARGIN * prof.cpu_min_per_process
        gplan = genecall_plan(con, tmp, batch_cpu)
    finally:
        con.close()
    d1 = _tsv(Path(a.pilot_depth1))
    d0 = _tsv(Path(a.pilot_depth1_unfiltered))
    rates = pilot_rates(d1)
    unf = pilot_rates(d0)
    fit = fit_bytes(d1 + d0)
    fasta_bytes, gz_bytes = int(rep["fasta_bytes"]), int(rep["gz_bytes"])
    ratio = max(1.0, fasta_bytes / gz_bytes)
    per_shard_disk = int((stage.SHARD_MARGIN * prof.cpu_min_per_process - stage.INGEST_FIXED_CPU_MIN)
                         * stage.INGEST_BP / stage.INGEST_CPU_MIN / ratio)
    # pass-1 shards: maximal runs of objects (in pass-1 order, by sha256) of at most per_shard_disk bytes on disk
    order = sorted(((r["sha256"], int(r["file_bytes"]), int(r["bp"])) for r in samples))
    shards, cur, cur_b, cur_bp = [], 0, 0, 0
    for _, fb, bpp in order:
        if cur and cur_b + fb > per_shard_disk:
            shards.append((cur, cur_b, cur_bp))
            cur, cur_b, cur_bp = 0, 0, 0
        cur, cur_b, cur_bp = cur + 1, cur_b + fb, cur_bp + bpp
    if cur:
        shards.append((cur, cur_b, cur_bp))
    largest_shard_bp = max(s[2] for s in shards)
    ingest = {
        "fasta_bytes": fasta_bytes, "gz_bytes": gz_bytes, "fasta_per_gz": round(ratio, 4),
        "held_bp": ov["root_any"]["bp"],
        "code_cpu_min_one_process": round(stage.ingest_cpu_min(fasta_bytes, ov["root_any"]["bp"]), 1),
        "code_with_c47_gzip_factor_1_2035": round(1.2035 * stage.ingest_cpu_min(fasta_bytes, ov["root_any"]["bp"]), 1),
        "pilot_sp0_unfiltered_gz_rate_cpu_min_per_gbp": round(unf["ingest_cpu_min_per_gbp"], 3),
        "pilot_sp0_rate_cpu_min": round(unf["ingest_cpu_min_per_gbp"] * int(rep["bp"]) / 1e9, 1),
        "pass1_bytes_on_disk_per_shard": per_shard_disk,
        "pass1_shards": len(shards),
        "largest_pass1_shard_bp": largest_shard_bp,
        "largest_pass1_shard_cpu_min_at_pilot_sp0_rate": round(unf["ingest_cpu_min_per_gbp"] * largest_shard_bp / 1e9, 1),
        "budget_cpu_min_per_process": prof.cpu_min_per_process,
    }
    gbp_new = new_bp / 1e9
    costs = {
        "inputs_cpu_min_pilot_rate": round(rates["inputs_cpu_min_per_m_genomes"] * rep["genomes"] / 1e6, 2),
        "ingest": ingest,
        "genecall_code_plan": gplan,
        "genecall_pilot_rate_cpu_min": round(rates["genecall_cpu_min_per_gbp"] * gbp_new, 1),
        "crispr_code_cpu_min": round(stage.crispr_cpu_min(new_bp), 1),
        "crispr_pilot_before_b55_cpu_min": (round(rates["crispr_cpu_min_per_gbp_before_b55"] * gbp_new, 1)
                                            if rates["crispr_cpu_min_per_gbp_before_b55"] else None),
        "igr_pilot_rate_cpu_min": round(rates["igr_cpu_min_per_gbp"] * gbp_new, 1),
        "tapes_pilot_rate_cpu_min": round(rates["tapes_cpu_min_per_gbp"] * gbp_new, 1),
        "new_bp": new_bp, "new_contigs": new_contigs,
        "pilot_rates": rates, "pilot_unfiltered_rates": unf,
    }
    costs["total_cpu_min_code"] = round(costs["inputs_cpu_min_pilot_rate"] + ingest["code_cpu_min_one_process"]
                                        + gplan["batches_cpu_min_sum"] + gplan["merge_estimate_cpu_min"]
                                        + costs["crispr_code_cpu_min"] + costs["igr_pilot_rate_cpu_min"]
                                        + costs["tapes_pilot_rate_cpu_min"], 1)
    costs["total_cpu_min_pilot_rates"] = round(costs["inputs_cpu_min_pilot_rate"] + ingest["pilot_sp0_rate_cpu_min"]
                                               + costs["genecall_pilot_rate_cpu_min"]
                                               + (costs["crispr_pilot_before_b55_cpu_min"] or 0)
                                               + costs["igr_pilot_rate_cpu_min"] + costs["tapes_pilot_rate_cpu_min"], 1)
    root_gb = {k: round((v["bytes_per_base"] * new_bp + v["bytes_per_contig"] * new_contigs) / 1e9, 2)
               for k, v in fit.items()}
    costs["root_gb_fit"] = root_gb
    costs["root_gb_fit_coefficients"] = fit
    costs["root_gb_at_pilot_1500_rate"] = round(rates["gb_per_gbp_called_total"] * gbp_new, 2)
    res["costs"] = costs
    print(json.dumps({"costs": {k: v for k, v in costs.items() if k not in ("pilot_rates", "pilot_unfiltered_rates")}}),
          flush=True)

    # 4
    if a.sample_bits > 0:
        objs = json.loads((inp / "objects.json").read_text())
        t0 = time.time()
        ms = measured_sample(objs, L, a.sample_bits, held)
        ms["wall_s"] = round(time.time() - t0, 1)
        pk = json.loads(Path(a.pilot_pack).read_text())
        spa1 = {"unique_bytes": pk["entries_bytes"]["unique"], "unique_aa": 9_996_937,
                "dark_bytes": pk["entries_bytes"]["dark_orf"], "dark_aa": 39_399_560,
                "igr_bytes": pk["entries_bytes"]["igr"], "igrs": 34_205,
                "source": "card C-42 depth1/pack_probe_SPA1.json; residues and IGR count from its note §4.3"}
        gs = ms["bp"] / 1e9
        per_gbp = {k: ms[k] / gs for k in ("genes", "distinct_proteins", "distinct_aa", "igrs", "dark_orfs", "dark_aa")}
        packs = {
            "unique_gb": round(per_gbp["distinct_aa"] * gbp_new * spa1["unique_bytes"] / spa1["unique_aa"] / 1e9, 2),
            "igr_gb": round(per_gbp["igrs"] * gbp_new * spa1["igr_bytes"] / spa1["igrs"] / 1e9, 2),
            "dark_orf_gb": round(per_gbp["dark_aa"] * gbp_new * spa1["dark_bytes"] / spa1["dark_aa"] / 1e9, 2),
        }
        packs["with_dark_orf_gb"] = round(packs["unique_gb"] + packs["igr_gb"] + packs["dark_orf_gb"], 2)
        packs["without_dark_orf_gb"] = round(packs["unique_gb"] + packs["igr_gb"], 2)
        packs["pilot_rate_with_dark_gb"] = round(pk["pack_gb_per_gbp"] * gbp_new, 2)
        res["sample"] = ms
        res["sample_per_gbp"] = {k: round(v, 1) for k, v in per_gbp.items()}
        res["sample_fragment_share"] = round(ms["fragments"] / ms["genes"], 4) if ms["genes"] else None
        res["pack"] = packs
        res["pack_bytes_rates_spa1"] = spa1
        print(json.dumps({"sample": ms, "pack": packs}), flush=True)

    # 5
    not_staged = {d["sample_id"]: d["status"] for d in rep["not_staged"]}
    rem = []
    for r in sel:
        sid = r["sample_id"]
        if sid in not_staged:
            rem.append({**r, "status": f"taken_{not_staged[sid]}"})
        elif sid in mv["drop_from_end"]:
            rem.append({**r, "status": "dropped_over_cap"})
    for r in _tsv(Path(a.select_dir) / "remainder.tsv.gz"):
        rem.append(r)
    rem.sort(key=lambda r: (r["status"] == "no_assembly_size", int(r["order_rank"]) if r.get("order_rank") else 0,
                            r["sample_id"]))
    cols = ["order_rank", "sample_id", "study", "habitat", "habitat2", "stratum_n", "stratum_rank", "assembly_bp",
            "assembly_n50", "est_n_ge_5000", "est_bp_ge_5000", "est_source", "status"]
    _write_tsv(out / "remainder_aquatic_terrestrial.tsv.gz", rem, cols, gz=True)
    below = []
    for key in ("habitat2", "habitat"):
        g = defaultdict(lambda: defaultdict(int))
        for r in samples:
            d = g[r[key]]
            d["samples"] += 1
            for k in ("records", "bp", "records_below", "bp_below", "records_at_least", "bp_at_least"):
                d[k] += int(r[k])
        for k in sorted(g):
            d = g[k]
            below.append({"level": key, "group": k, **d, "share_records_below": _ratio(d["records_below"], d["records"]),
                          "share_bp_below": _ratio(d["bp_below"], d["bp"])})
    _write_tsv(out / "below_length_by_stratum.tsv", below)
    shutil.copyfile(Path(a.select_dir) / "other_habitats.tsv", out / "other_habitats.tsv")
    res["remainder"] = {
        "rows": len(rem),
        "by_status": {s: sum(1 for r in rem if r["status"] == s) for s in sorted({r["status"] for r in rem})},
        "not_taken_assembly_bp": sum(int(r["assembly_bp"]) for r in rem if r["assembly_bp"]),
        "not_taken_est_n_ge_5000": sum(int(r["est_n_ge_5000"]) for r in rem if r["est_n_ge_5000"]),
        "not_taken_est_bp_ge_5000": sum(int(r["est_bp_ge_5000"]) for r in rem if r["est_bp_ge_5000"]),
        "staged_below_length": {"records": rep["records_below"], "bp": rep["bp_below"],
                                "share_records": _ratio(rep["records_below"], rep["records"]),
                                "share_bp": _ratio(rep["bp_below"], rep["bp"])},
    }
    shutil.rmtree(tmp, ignore_errors=True)
    (out / "est.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps({"remainder": res["remainder"]}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

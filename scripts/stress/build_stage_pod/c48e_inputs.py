"""Card C-48e part 1: a SPIRE tranche's `build-stage` inputs, staged whole, as the pilot's depth-1 stages had them.

    python c48e_inputs.py build --selection SEL.tsv --fetch-log DL/fetch_log.jsonl --dl-dir DL --stage SPX1 \
        --final-dir /tmp/if2453/s2/inputs/SPX1 --out-dir OUT --min-contig-len 5000 \
        --uids-out /tmp/if2453/c48e1/work/contigs_ge5000.parquet --tmp /tmp/if2453/c48e1/tmp
    python c48e_inputs.py check-summary CHECK_FULL.json CHECK.json

What `build` writes into OUT (first as `.<name>.tmp` beside it, renamed when complete; refused if OUT exists):
- `assemblies/<sample_id>.fa.gz`: a hard link to each fetched assembly, whole (every contig the provider published).
  The length filter is ingest's (`build-stage --min-contig-len`, card C-38d), so that the ingest ledger records the
  records and bases it drops (`below:<L>`).
- `objects.json`: one source object per staged assembly with at least one record, its path under `--final-dir`
  (where OUT is renamed to; a second run into scratch with the same `--final-dir` gives the same bytes), source
  `spire` `v1`, CC BY-SA 4.0, redistributable, `retrieved_at` the fetch log's time; in sample_id order.
- `genomes.json`: as the pilot's `c42_inputs.py` (card C-42): one genome per distinct contig, kind `prokaryote`,
  genetic code 11, no circular contig, assembly quality `contig`; here only the contigs of at least
  `--min-contig-len` bases, which is what ingest with that length requires of `genome_input` (no shorter contig, every
  longer one); in contig_uid order. One genome per contig and not per sample: a sample is a mixture of organisms, and
  a sample-sized genome would be called in Prodigal's single mode, trained on the mixture; a contig under 100 kb is
  called in meta mode, one of 100 kb or more in single mode trained on itself (`genecall.mode_of`).
- `samples.tsv` (per sample: the file's bytes and sha256, records and bases in total, under and at least the length,
  at least 10 kb and 100 kb), `assemblies.sha256` (sha256sum form), `input_report.json` (totals, by habitat and
  stratum, contigs shared between samples). No time, host or absolute path of the run is written.
Bounded memory (AGENTS §4.1): one record at a time; the long contigs' uids go to a scratch Parquet file in row groups
and are aggregated and sorted by one DuckDB connection (one thread, 1 GB, spilling to `--tmp`). Each assembly's size
and sha256 are checked against the fetch log before it is read; a mismatch refuses the build.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import catalog, fasta, ids

ROW_GROUP = 122_880
SINGLE_MODE_MIN_BP = 100_000
_LONG = pa.schema([("contig_uid", pa.binary(16)), ("contig_len", pa.int64()), ("sample", pa.int32())])
_UIDS = pa.schema([("contig_uid", pa.binary(16)), ("contig_len", pa.int64()), ("n_samples", pa.int32()),
                   ("n_records", pa.int32())])


class Refused(Exception):
    """Inputs that cannot be staged (nothing is left in OUT)."""


def _sibling(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def _read_selection(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    ids_ = [r["sample_id"] for r in rows]
    if len(ids_) != len(set(ids_)):
        raise Refused(f"{path}: a sample is listed twice")
    return sorted(rows, key=lambda r: r["sample_id"])


SAMPLE_COLS = ["sample_id", "study", "habitat", "habitat2", "order_rank", "retrieved_at", "file_bytes", "sha256",
               "records", "bp", "records_below", "bp_below", "records_at_least", "bp_at_least", "records_ge_10000",
               "bp_ge_10000", "records_ge_100000", "bp_ge_100000"]


def build(selection: Path, fetch_log: Path, dl_dir: Path, stage: str, final_dir: Path, out_dir: Path, *,
          min_contig_len: int, uids_out: Path, tmp: Path) -> dict:
    """Write the stage's inputs into `out_dir` (module docstring); returns the report."""
    fetch_mod = _sibling("c48e_fetch")
    source = _sibling("c42_inputs").SPIRE_SOURCE
    catalog.check_stages([stage])
    out_dir, final_dir, dl_dir, uids_out = Path(out_dir), Path(final_dir), Path(dl_dir), Path(uids_out)
    if out_dir.exists():
        raise Refused(f"{out_dir} already exists")
    if min_contig_len < 1:
        raise Refused(f"min_contig_len must be >= 1, got {min_contig_len}")
    log = fetch_mod.read_log(Path(fetch_log))
    rows = _read_selection(Path(selection))
    staged = [r for r in rows if log.get(r["sample_id"], {}).get("status") == "ok"]
    not_staged = [{"sample_id": r["sample_id"], "status": log.get(r["sample_id"], {}).get("status", "not_fetched")}
                  for r in rows if log.get(r["sample_id"], {}).get("status") != "ok"]

    work = out_dir.with_name(f".{out_dir.name}.tmp")
    if work.exists():
        shutil.rmtree(work)
    (work / "assemblies").mkdir(parents=True)
    Path(tmp).mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix="c48e-inputs-", dir=tmp))
    try:
        rep = _build(staged, not_staged, log, dl_dir, stage, final_dir, work, scratch, min_contig_len, uids_out,
                     source)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    work.rename(out_dir)
    return rep


def _build(staged, not_staged, log, dl_dir, stage, final_dir, work, scratch, L, uids_out, source) -> dict:
    sample_rows, objects, sums = [], [], []
    writer = pq.ParquetWriter(scratch / "long.parquet", _LONG, compression="zstd", write_statistics=False)
    buf: list[list] = [[], [], []]

    def flush():
        if buf[0]:
            writer.write_batch(pa.RecordBatch.from_arrays([pa.array(buf[0], pa.binary(16)), pa.array(buf[1], pa.int64()),
                                                           pa.array(buf[2], pa.int32())], schema=_LONG))
            for b in buf:
                b.clear()

    tot = defaultdict(int)
    for i, r in enumerate(staged):
        sid = r["sample_id"]
        rec = log[sid]
        src = dl_dir / f"{sid}.fa.gz"
        if not src.is_file() or src.stat().st_size != rec["disk_bytes"] or _sha256(src) != rec["disk_sha256"]:
            raise Refused(f"{sid}: {src} is missing or differs from the fetch log (size or sha256)")
        s = defaultdict(int)
        for fr in fasta.iter_records(src):
            n = len(fr.seq)
            s["records"] += 1
            s["bp"] += n
            if n < L:
                s["records_below"] += 1
                s["bp_below"] += n
                continue
            s["records_at_least"] += 1
            s["bp_at_least"] += n
            if n >= 10_000:
                s["records_ge_10000"] += 1
                s["bp_ge_10000"] += n
            if n >= SINGLE_MODE_MIN_BP:
                s["records_ge_100000"] += 1
                s["bp_ge_100000"] += n
            buf[0].append(ids.contig_uid(fr.seq))
            buf[1].append(n)
            buf[2].append(i)
            if len(buf[0]) >= ROW_GROUP:
                flush()
        os.link(src, work / "assemblies" / f"{sid}.fa.gz")
        sums.append(f"{rec['disk_sha256']}  assemblies/{sid}.fa.gz\n")
        if s["records"]:
            objects.append({"path": str(final_dir / "assemblies" / f"{sid}.fa.gz"), **source,
                            "retrieved_at": rec["fetched_utc"]})
        sample_rows.append({**{k: r.get(k, "") for k in ("sample_id", "study", "habitat", "habitat2", "order_rank")},
                            "retrieved_at": rec["fetched_utc"], "file_bytes": rec["disk_bytes"],
                            "sha256": rec["disk_sha256"], **{k: s[k] for k in SAMPLE_COLS[8:]}})
        for k in SAMPLE_COLS[8:]:
            tot[k] += s[k]
        tot["fasta_bytes"] += int(rec.get("uncompressed_bytes") or 0)
        tot["gz_bytes"] += int(rec["disk_bytes"])
    flush()
    writer.close()

    # distinct long contigs: one DuckDB pass, sorted by contig_uid
    con = duckdb.connect(config={"threads": 1, "memory_limit": "1GB", "temp_directory": str(scratch / "duck"),
                                 "preserve_insertion_order": False})
    agg = dict(n=0, bp=0, multi=0, within=0, single=0, single_bp=0, lmismatch=0)
    uids_tmp = uids_out.with_name(f".{uids_out.name}.tmp")
    uids_out.parent.mkdir(parents=True, exist_ok=True)
    try:
        lit = "'" + str(scratch / "long.parquet").replace("'", "''") + "'"
        reader = con.sql(
            f"SELECT contig_uid, max(contig_len)::BIGINT AS contig_len, count(DISTINCT sample)::INTEGER AS n_samples, "
            f"count(*)::INTEGER AS n_records, min(contig_len)::BIGINT AS lmin FROM read_parquet({lit}) "
            "GROUP BY contig_uid ORDER BY contig_uid").to_arrow_reader(ROW_GROUP)
        uw = pq.ParquetWriter(uids_tmp, _UIDS, compression="zstd")
        with open(work / "genomes.json", "w") as gj:
            gj.write("[")
            first = True
            for b in reader:
                cu, ln, ns, nr, lmin = (b.column(j).to_pylist() for j in range(5))
                uw.write_batch(pa.RecordBatch.from_arrays([pa.array(cu, pa.binary(16)), pa.array(ln, pa.int64()),
                                                           pa.array(ns, pa.int32()), pa.array(nr, pa.int32())],
                                                          schema=_UIDS), row_group_size=ROW_GROUP)
                parts = []
                for u, n, s_, r_, lo in zip(cu, ln, ns, nr, lmin):
                    agg["n"] += 1
                    agg["bp"] += n
                    agg["multi"] += s_ > 1
                    agg["within"] += r_ - s_
                    agg["lmismatch"] += lo != n
                    if n >= SINGLE_MODE_MIN_BP:
                        agg["single"] += 1
                        agg["single_bp"] += n
                    parts.append(json.dumps({"contig_uids": [bytes(u).hex()], "genetic_code": 11, "kind": "prokaryote",
                                             "circular": [], "assembly_quality": "contig"}))
                if parts:
                    gj.write(("" if first else ", ") + ", ".join(parts))
                    first = False
            gj.write("]\n")
        uw.close()
    finally:
        con.close()
    if agg["lmismatch"]:
        raise Refused(f"{agg['lmismatch']} contig uids with two lengths")
    uids_tmp.replace(uids_out)

    (work / "objects.json").write_text(json.dumps(objects, indent=1) + "\n")
    (work / "assemblies.sha256").write_text("".join(sorted(sums)))
    with open(work / "samples.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SAMPLE_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sample_rows)

    def by(key):
        g: dict[str, dict] = {}
        for r in sample_rows:
            d = g.setdefault(r[key], defaultdict(int))
            d["samples"] += 1
            for k in ("records", "bp", "records_below", "bp_below", "records_at_least", "bp_at_least"):
                d[k] += r[k]
        return {k: dict(v) for k, v in sorted(g.items())}

    rep = {
        "stage": stage,
        "min_contig_len": L,
        "samples_selected": len(staged) + len(not_staged),
        "samples_staged": len(staged),
        "not_staged": not_staged,
        "objects": len(objects),
        "gz_bytes": tot["gz_bytes"],
        "fasta_bytes": tot["fasta_bytes"],
        **{k: tot[k] for k in SAMPLE_COLS[8:]},
        "distinct_contigs": agg["n"],
        "distinct_bp": agg["bp"],
        "genomes": agg["n"],
        "contigs_in_2plus_samples": agg["multi"],
        "repeated_records_between_samples": tot["records_at_least"] - agg["n"] - agg["within"],
        "within_sample_duplicate_records": agg["within"],
        "distinct_contigs_ge_100000": agg["single"],
        "distinct_bp_ge_100000": agg["single_bp"],
        "by_habitat": by("habitat"),
        "by_habitat2": by("habitat2"),
    }
    (work / "input_report.json").write_text(json.dumps(rep, indent=1) + "\n")
    return rep


def cmd_check_summary(full: Path, out: Path) -> int:
    d = json.loads(Path(full).read_text())
    objs = d.pop("objects")
    d["objects_summary"] = {
        "count": len(objs), "total_bytes": sum(o["bytes"] for o in objs),
        "source_names": sorted({o["source_name"] for o in objs}),
        "source_releases": sorted({o["source_release"] for o in objs}),
        "licenses": sorted({o["license"] for o in objs}),
        "redistributable": sorted({o["redistributable"] for o in objs}),
        "first": objs[0] if objs else None, "last": objs[-1] if objs else None,
    }
    Path(out).write_text(json.dumps(d, indent=1) + "\n")
    print(json.dumps({k: v for k, v in d.items() if k != "objects_summary"}))
    return 0 if not d["problems"] else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--selection", required=True)
    b.add_argument("--fetch-log", required=True)
    b.add_argument("--dl-dir", required=True)
    b.add_argument("--stage", required=True)
    b.add_argument("--final-dir", required=True)
    b.add_argument("--out-dir", required=True)
    b.add_argument("--min-contig-len", type=int, required=True)
    b.add_argument("--uids-out", required=True)
    b.add_argument("--tmp", required=True)
    c = sub.add_parser("check-summary")
    c.add_argument("full")
    c.add_argument("out")
    a = ap.parse_args(argv)
    if a.cmd == "check-summary":
        return cmd_check_summary(Path(a.full), Path(a.out))
    try:
        rep = build(Path(a.selection), Path(a.fetch_log), Path(a.dl_dir), a.stage, Path(a.final_dir), Path(a.out_dir),
                    min_contig_len=a.min_contig_len, uids_out=Path(a.uids_out), tmp=Path(a.tmp))
    except Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return 1
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

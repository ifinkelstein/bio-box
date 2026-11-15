"""Card C-33: the checks after the two CRISPR producers ran on the built stages of the data root. In the image
(pyarrow; `biobox.stage.stage_bases` from the tree's src), read-only on the root:

    python c33_check.py --root ROOT --before DIR --after DIR --units UNITS.TSV --logs DIR --out DIR \
        [--view NAME=DIR ...] [--oracle-stages S1,P1,P2]

1. Before against after (`listing.tsv`: path, bytes, mtime of every file of the root; `p1s1.sha256`: sha256 of
   every file of stages P1 and S1): no file of the root that existed before is missing, has another size or mtime,
   or (P1, S1) another sha256. The files that are new must be exactly the expected ones: for each stage and producer
   its `feature` and `l1_crispr` parts, its ledger part and its run record. Anything else new is listed as unexpected.
2. The oracle: for each of `--oracle-stages` and each producer, the root's `feature`, `l1_crispr` and ledger files are
   byte-identical (sha256, same file names) to those of each `--view` that holds that producer.
3. Per stage and producer: arrays kept (rows of the `feature` files), contigs searched and skipped, calls in and
   rejected by reason (the ledger), CPU-minutes, wall seconds, peak RSS, image digest and code id (the run record),
   the unit's MemoryPeak and CPUUsageNSec (`pod_run: <unit> usage` in its log), and the estimate: the stage's genome
   bases (`stage_bases`) at 182.0 CPU-seconds per Gbp, times 1.19 for the filtered producer (`CRISPR_CPU_S_PER_GBP`,
   `CRISPR_NEAR_CAS_COST` in src/biobox/stage.py).
4. `root_added.tsv`: every new file with its bytes and sha256, for the durable copy.

Writes into --out: report.txt, summary.json, root_added.tsv, per_stage.tsv, oracle.tsv. Exit 0 when every check
passes, 1 otherwise (the report says which failed).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import pyarrow.parquet as pq

from biobox.stage import CRISPR_CPU_S_PER_GBP, CRISPR_NEAR_CAS_COST, stage_bases

PRODUCERS = {"nr3": "minced-0.4.2-w8-nr3-r23x47-s26x50",
             "near": "minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2"}
STAGES = ["S1", "P1", "P2", "P3", "P4", "P5", "P6"]
RATE_TOLERANCE = 0.20  # AGENTS §6: a measured rate more than 20% from its estimate stops the card


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def listing(d: Path) -> dict[str, tuple[int, str]]:
    out = {}
    for line in (d / "listing.tsv").read_text().splitlines():
        path, size, mtime = line.split("\t")
        out[path] = (int(size), mtime)
    return out


def sums(d: Path) -> dict[str, str]:
    out = {}
    for line in (d / "p1s1.sha256").read_text().splitlines():
        h, path = line.split("  ", 1)
        out[path] = h
    return out


def dirs(stage: str, prod: str) -> list[str]:
    return [f"l1/feature/producer={prod}/stage={stage}", f"l1/crispr/producer={prod}/stage={stage}",
            f"ledger/stage={stage}/step=crispr-{prod}"]


def files_of(root: Path, stage: str, prod: str) -> list[str]:
    return sorted(p.relative_to(root).as_posix() for d in dirs(stage, prod) for p in (root / d).glob("*.parquet"))


def ledger_counts(path: Path) -> dict:
    t = pq.read_table(path, columns=["n_in", "n_kept", "n_rejected", "reason_code"]).to_pydict()
    c = {"contigs_searched": 0, "contigs_skipped_other_stage": 0, "contigs_shared_in_stage": 0, "calls_in": 0,
         "calls_kept": 0, "rejected": {}}
    for n_in, n_kept, n_rej, why in zip(t["n_in"], t["n_kept"], t["n_rejected"], t["reason_code"]):
        if why == "ok":
            c["contigs_searched"] += 1
            c["calls_in"] += n_in
            c["calls_kept"] += n_kept
        elif why.startswith("owned_by:"):
            c["contigs_shared_in_stage" if "/" in why else "contigs_skipped_other_stage"] += 1
        else:
            c["rejected"][why] = c["rejected"].get(why, 0) + n_rej
    return c


def unit_usage(log: Path, unit: str) -> dict:
    out = {}
    if not log.exists():
        return out
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith(f"pod_run: {unit} usage"):
            out = {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", line)}
        m = re.match(rf"pod_run: {re.escape(unit)} exited (\d+)$", line)
        if m:
            out["exit"] = int(m[1])
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--before", required=True, type=Path)
    ap.add_argument("--after", required=True, type=Path)
    ap.add_argument("--units", required=True, type=Path, help="TSV: unit, stage, kind (nr3|near), step")
    ap.add_argument("--logs", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--view", action="append", default=[], help="NAME=DIR: a scratch view for the oracle")
    ap.add_argument("--oracle-stages", default="S1,P1,P2")
    a = ap.parse_args()
    root, out = a.root, a.out
    out.mkdir(parents=True, exist_ok=True)
    report: list[str] = []
    failed: list[str] = []

    def say(s: str) -> None:
        report.append(s)
        print(s, flush=True)

    # 1. before against after
    b, f = listing(a.before), listing(a.after)
    gone = sorted(p for p in b if p not in f)
    changed = sorted(p for p in b if p in f and f[p] != b[p])
    new = sorted(p for p in f if p not in b)
    expected = set()
    for s in STAGES:
        for prod in PRODUCERS.values():
            expected.add(f"runs/stage={s}/crispr-{prod}.json")
    expected_dirs = tuple(d + "/" for s in STAGES for prod in PRODUCERS.values() for d in dirs(s, prod))
    unexpected = sorted(p for p in new if p not in expected and not (p.startswith(expected_dirs)
                                                                       and p.endswith(".parquet")))
    missing_new = sorted(p for p in expected if p not in f)
    sb, sa = sums(a.before), sums(a.after)
    sha_changed = sorted(p for p in sb if sa.get(p) != sb[p])
    say(f"before: {len(b)} files; after: {len(f)} files; new {len(new)}; gone {len(gone)}; size or mtime changed "
        f"{len(changed)}; P1 and S1 files hashed before {len(sb)}, after {len(sa)}, sha256 changed {len(sha_changed)}")
    for what, lst in (("gone", gone), ("changed", changed), ("sha256 changed", sha_changed),
                      ("unexpected new", unexpected), ("expected but absent", missing_new)):
        if lst:
            failed.append(f"{what}: {len(lst)}")
            for p in lst[:50]:
                say(f"  {what}: {p}")
    new_in_p1s1 = sorted(p for p in sa if p not in sb)
    say(f"new files among P1 and S1 in the after hashes: {len(new_in_p1s1)} (all must be the producers')")

    # 4. root_added.tsv
    added = []
    for p in new:
        added.append((p, (root / p).stat().st_size, sha256(root / p)))
    with open(out / "root_added.tsv", "w") as fh:
        fh.write("path\tbytes\tsha256\n")
        for p, n, h in added:
            fh.write(f"{p}\t{n}\t{h}\n")
    say(f"root_added.tsv: {len(added)} files, {sum(n for _, n, _ in added)} bytes")
    for p, n, h in added:
        if p in sa and sa[p] != h:
            failed.append(f"after-hash of {p} differs from its hash now")

    # 2. oracle
    views = dict(v.split("=", 1) for v in a.view)
    orows = []
    for s in a.oracle_stages.split(","):
        for kind, prod in PRODUCERS.items():
            mine = files_of(root, s, prod)
            for name, vd in views.items():
                vd = Path(vd)
                theirs = files_of(vd, s, prod)
                if not theirs:
                    continue  # this view does not hold the producer
                verdict = "equal"
                if mine != theirs:
                    verdict = f"file names differ: root {mine} view {theirs}"
                for rel in sorted(set(mine) | set(theirs)):
                    hr = sha256(root / rel) if (root / rel).exists() else "absent"
                    hv = sha256(vd / rel) if (vd / rel).exists() else "absent"
                    same = hr == hv
                    orows.append((s, kind, name, rel, hr, hv, "same" if same else "DIFFERENT"))
                    if not same:
                        verdict = "DIFFERENT"
                say(f"oracle {s} {kind} against {name}: {len(mine)} files, {verdict}")
                if verdict != "equal":
                    failed.append(f"oracle {s} {kind} {name}: {verdict}")
    with open(out / "oracle.tsv", "w") as fh:
        fh.write("stage\tproducer\tview\tpath\tsha256_root\tsha256_view\tverdict\n")
        for r in orows:
            fh.write("\t".join(r) + "\n")

    # 3. per stage and producer
    units = {}
    for line in a.units.read_text().splitlines():
        unit, s, kind, step = line.split("\t")
        units[(s, kind)] = unit  # the last launch of each wins
    rows = []
    for kind, prod in PRODUCERS.items():
        for s in STAGES:
            rec_p = root / "runs" / f"stage={s}" / f"crispr-{prod}.json"
            if not rec_p.exists():
                failed.append(f"no run record {rec_p}")
                continue
            rec = json.loads(rec_p.read_text())
            feats = [root / p for p in files_of(root, s, prod) if p.startswith("l1/feature/")]
            elems = [root / p for p in files_of(root, s, prod) if p.startswith("l1/crispr/")]
            led = [root / p for p in files_of(root, s, prod) if p.startswith("ledger/")]
            arrays = sum(pq.ParquetFile(p).metadata.num_rows for p in feats)
            elements = sum(pq.ParquetFile(p).metadata.num_rows for p in elems)
            lc = ledger_counts(led[0]) if len(led) == 1 else {}
            if len(led) != 1:
                failed.append(f"{s} {kind}: {len(led)} ledger parts")
            if lc and lc["calls_kept"] != arrays:
                failed.append(f"{s} {kind}: ledger kept {lc['calls_kept']}, feature rows {arrays}")
            bases = stage_bases(root, s)
            est = bases * CRISPR_CPU_S_PER_GBP / 1e9 / 60 * (CRISPR_NEAR_CAS_COST if kind == "near" else 1.0)
            ratio = rec["cpu_min"] / est
            unit = units.get((s, kind), "")
            use = unit_usage(a.logs / f"{unit}.log", unit) if unit else {}
            row = {"stage": s, "producer": kind, "bases": bases, "arrays_kept": arrays, "elements": elements,
                   **{k: v for k, v in lc.items() if k != "rejected"},
                   **{f"rejected_{k}": v for k, v in sorted(lc.get("rejected", {}).items())},
                   "cpu_min": rec["cpu_min"], "wall_s": rec["wall_s"], "peak_rss_bytes": rec["peak_rss_bytes"],
                   "estimate_cpu_min": round(est, 3), "measured_over_estimate": round(ratio, 3),
                   "cpu_s_per_gbp": round(rec["cpu_min"] * 60 / (bases / 1e9), 1), "threads": rec["threads"],
                   "image_digest": rec["image_digest"], "code_id": rec["code_id"], "code_commit": rec["code_commit"],
                   "started_at": rec["started_at"], "unit": unit, "unit_memory_peak": use.get("MemoryPeak"),
                   "unit_cpu_ns": use.get("CPUUsageNSec"), "unit_exit": use.get("exit"),
                   "over_cpu_budget": rec["over_cpu_budget"]}
            rows.append(row)
            flag = "" if abs(ratio - 1) <= RATE_TOLERANCE else "  OUTSIDE 20%"
            if flag:
                failed.append(f"{s} {kind}: CPU-minutes {rec['cpu_min']:.3f} against the estimate {est:.3f}")
            say(f"{s} {kind}: arrays {arrays}, calls in {lc.get('calls_in')}, rejected {lc.get('rejected')}, contigs "
                f"{lc.get('contigs_searched')} searched, {lc.get('contigs_skipped_other_stage')} another stage's; "
                f"{rec['cpu_min']:.3f} CPU-min (estimate {est:.3f}, ratio {ratio:.3f}){flag}, {rec['wall_s']:.1f} s, "
                f"peak {rec['peak_rss_bytes'] / 1e6:.1f} MB, unit peak {use.get('MemoryPeak')}")
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out / "per_stage.tsv", "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
    tot = {kind: round(sum(r["cpu_min"] for r in rows if r["producer"] == kind), 3) for kind in PRODUCERS}
    say(f"CPU-minutes in all: {tot}; images {sorted({r['image_digest'] for r in rows})}; code ids "
        f"{sorted({r['code_id'] for r in rows})}")
    (out / "summary.json").write_text(json.dumps(
        {"failed": failed, "new_files": len(new), "new_bytes": sum(n for _, n, _ in added), "gone": gone,
         "changed": changed, "sha256_changed": sha_changed, "unexpected_new": unexpected, "per_stage": rows,
         "cpu_min_total": tot}, indent=1, sort_keys=True))
    say("RESULT: " + ("all checks pass" if not failed else "FAILED: " + "; ".join(failed)))
    (out / "report.txt").write_text("\n".join(report) + "\n")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())

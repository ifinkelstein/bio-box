"""Export a stage's distinct proteins and its producer Pfam domain rows (card 0-09a). Reads only.

    python export_stage.py <root> <stage> <out_dir> [--library pfam-38.2] [--shard K] [--head N]

Files read are the producer's own lists, never globs: the proteins are the `protein_seq` files named in
the stage's shard plan (`l1/domain_scan/library=<L>/stage=<S>/plan.json`, which the shards scanned); the
truth is the `l1_domains` files named in the published step's run record
(`runs/stage=<S>/domains-<L>.json`, its `outputs`). `--shard K` keeps the plan's shard K (its uid range);
`--head N` then keeps its first N proteins by uid. The truth is cut to the proteins kept.

Writes, under <out_dir>:
  proteins.fa   ">" uid hex, the sequence without '*' (as `domains.scan` strips it)
  proteins.tsv  uid, length (residues after stripping)
  truth.tsv     uid, acc, env_start, env_end, ali_start, ali_end, seq_bits, dom_bits (0-based half-open)
  export.json   counts, the files read and their sha256
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pyarrow.parquet as pq


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("stage")
    ap.add_argument("out", type=Path)
    ap.add_argument("--library", default="pfam-38.2")
    ap.add_argument("--shard", type=int)
    ap.add_argument("--head", type=int)
    a = ap.parse_args()
    root, out = a.root, a.out
    out.mkdir(parents=True, exist_ok=True)
    plan = json.loads((root / f"l1/domain_scan/library={a.library}/stage={a.stage}/plan.json").read_text())
    record = json.loads((root / f"runs/stage={a.stage}/domains-{a.library}.json").read_text())
    if record.get("exit_code") != 0:
        raise SystemExit(f"run record of domains-{a.library} for stage {a.stage} has exit {record.get('exit_code')}")
    lo = hi = None
    if a.shard is not None:
        s = plan["shards"][a.shard]
        lo, hi = bytes.fromhex(s["first"]), bytes.fromhex(s["last"])
    kept: set[bytes] = set()
    n = residues = 0
    read = []
    with open(out / "proteins.fa", "w") as fa, open(out / "proteins.tsv", "w") as tsv:
        tsv.write("uid\tlength\n")
        for f in plan["files"]:
            path = root / f["path"]
            read.append({"path": f["path"], "bytes": path.stat().st_size, "sha256": sha256(path)})
            if path.stat().st_size != f["bytes"]:
                raise SystemExit(f"{path}: size differs from the plan")
            for b in pq.ParquetFile(path).iter_batches(batch_size=16384, columns=["protein_uid", "seq"]):
                for u, seq in zip(b.column(0).to_pylist(), b.column(1).to_pylist()):
                    if lo is not None and not lo <= u <= hi:
                        continue
                    if a.head is not None and n >= a.head:
                        break
                    seq = seq.replace("*", "")
                    h = u.hex()
                    fa.write(f">{h}\n{seq}\n")
                    tsv.write(f"{h}\t{len(seq)}\n")
                    kept.add(u)
                    n += 1
                    residues += len(seq)
    rows = 0
    truth_files = [o for o in record["outputs"] if o["path"].startswith("l1/domains/")]
    cols = ["protein_uid", "model_acc", "env_start", "env_end", "ali_start", "ali_end", "seq_bitscore", "dom_bitscore"]
    with open(out / "truth.tsv", "w") as fh:
        fh.write("uid\tacc\tenv_start\tenv_end\tali_start\tali_end\tseq_bits\tdom_bits\n")
        for o in truth_files:
            path = root / o["path"]
            got = sha256(path)
            if got != o["sha256"]:
                raise SystemExit(f"{path}: sha256 {got} is not the run record's {o['sha256']}")
            read.append({"path": o["path"], "bytes": path.stat().st_size, "sha256": got})
            for b in pq.ParquetFile(path).iter_batches(batch_size=65536, columns=cols):
                d = b.to_pydict()
                for i in range(b.num_rows):
                    u = d["protein_uid"][i]
                    if u not in kept:
                        continue
                    fh.write(f"{u.hex()}\t{d['model_acc'][i]}\t{d['env_start'][i]}\t{d['env_end'][i]}\t"
                             f"{d['ali_start'][i]}\t{d['ali_end'][i]}\t{d['seq_bitscore'][i]:.2f}\t"
                             f"{d['dom_bitscore'][i]:.2f}\n")
                    rows += 1
    doc = {"root": str(root), "stage": a.stage, "library": a.library, "shard": a.shard, "head": a.head,
           "plan_proteins": plan["proteins"], "plan_residues": plan["residues"],
           "proteins": n, "residues": residues, "truth_rows": rows, "files": read}
    (out / "export.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: doc[k] for k in ("stage", "shard", "head", "proteins", "residues", "truth_rows")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Name the stage records behind each overlap `s2_seqstats.py` found, read-only (S2 inputs, 2026-11-11).

    python s2_overlap_labels.py ROOT STAGE_UIDS_TSV OVERLAP_TSV OUT_TSV

For a stored stage, the accessions are those of `ROOT/l0seq/contig_occurrence/stage=<S>/source=*` (every source
record of that sequence in the stage, `contig_uid` and `accession` columns, files from
`biobox.io.part_files`); for a fetched stage (how `fasta`), the accession column of STAGE_UIDS_TSV.
Nothing in ROOT is written. Writes OUT_TSV: the overlap rows (match `same` only) with the stage's
accessions comma-joined, and prints counts per stage.
"""

import json
import sys
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

from biobox import io


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    root, stage_tsv, ov_tsv, out = Path(argv[0]), Path(argv[1]), Path(argv[2]), Path(argv[3])
    rows = []
    with open(ov_tsv) as fh:
        fh.readline()
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    rows = [r for r in rows if r[4] == "same"]
    want = {(r[3], bytes.fromhex(r[2])) for r in rows}
    names: dict[tuple[str, bytes], list[str]] = {}
    for stage in sorted({s.split("/")[0] for s, _ in want if s.endswith("/store")}):
        uids = {u for s, u in want if s == f"{stage}/store"}
        for f in io.part_files(root / "l0seq" / "contig_occurrence" / f"stage={stage}", "source=*/part-*.parquet"):
            t = pq.read_table(f, columns=["contig_uid", "accession"])
            for u, a in zip(t.column("contig_uid").to_pylist(), t.column("accession").to_pylist()):
                if u in uids:
                    names.setdefault((f"{stage}/store", u), []).append(a)
    with open(stage_tsv) as fh:
        fh.readline()
        for ln in fh:
            s, how, u, _, a = ln.rstrip("\n").split("\t")
            key = (f"{s}/{how}", bytes.fromhex(u))
            if how == "fasta" and key in want:
                names.setdefault(key, []).append(a)
    tmp = out.with_name(f".{out.name}.tmp")
    with open(tmp, "w") as w:
        w.write("ordinal\taccession\tcontig_uid\tstage\tstage_accessions\n")
        for r in rows:
            w.write("\t".join([r[0], r[1], r[2], r[3], ",".join(sorted(names.get((r[3], bytes.fromhex(r[2])), [])))])
                    + "\n")
    tmp.replace(out)
    print(json.dumps({"rows": len(rows), "per_stage": dict(sorted(Counter(r[3] for r in rows).items())),
                      "stage_records_named": {s: sum(len(v) for k, v in names.items() if k[0] == s)
                                              for s in sorted({k[0] for k in names})}}))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

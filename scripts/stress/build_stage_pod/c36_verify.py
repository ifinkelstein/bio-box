"""Card C-36 (2): every genome of a stage checked against its fetch's lock, and each record's contig uid and length.

    python c36_verify.py STAGE_DIR OUT_TSV

STAGE_DIR holds the stage's lock.json and selection.tsv (`c36_stages.py assign`). For each genome, in accession order:
the .fna.gz named by `file` must have the lock's `gz_bytes` and `gz_sha256` (the bytes archived and extracted), and
its decompressed bytes the lock's `bytes`, `sha256` (the FASTA NCBI's md5 was checked on at the fetch) and `md5`;
the records read by `biobox.fasta.iter_records` (the reader ingest and make_inputs.py use) must number the lock's
`records`, and their bases sum to the lock's `bp` and to GTDB's `genome_size` (selection.tsv). Each record's contig
uid is `biobox.ids.contig_uid` of its sequence, the store's key. Writes OUT_TSV (accession, record ordinal, contig_uid
hex, length, record accession), and prints a JSON report; exit 1 when any genome fails a check (each failure listed),
else 0.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

from biobox import fasta, ids

READ = 1 << 20


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    stage_dir, out = Path(argv[0]), Path(argv[1])
    lock = json.loads((stage_dir / "lock.json").read_text())["genomes"]
    with open(stage_dir / "selection.tsv", newline="") as fh:
        gsize = {r["accession"]: int(r["genome_size"]) for r in csv.DictReader(fh, delimiter="\t")}
    if set(gsize) != set(lock):
        print(f"c36_verify.py: selection and lock of {stage_dir} name different genomes", file=sys.stderr)
        return 2
    problems: list[str] = []
    rows: list[tuple] = []
    tot = {"genomes": 0, "records": 0, "bp": 0, "fasta_bytes": 0, "gz_bytes": 0}
    for acc in sorted(lock):
        e = lock[acc]
        p = Path(e["file"])
        gz = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
        if gz != e["gz_sha256"] or p.stat().st_size != e["gz_bytes"]:
            problems.append(f"{acc}: {p} gz sha256 {gz} or size is not the lock's")
            continue
        h, m, n = hashlib.sha256(), hashlib.md5(), 0
        with gzip.open(p, "rb") as fh:
            while chunk := fh.read(READ):
                h.update(chunk)
                m.update(chunk)
                n += len(chunk)
        if (h.hexdigest(), m.hexdigest(), n) != (e["sha256"], e["md5"], e["bytes"]):
            problems.append(f"{acc}: FASTA sha256/md5/bytes {h.hexdigest()} {m.hexdigest()} {n} are not the lock's")
            continue
        k = bp = 0
        for k, r in enumerate(fasta.iter_records(p), 1):
            rows.append((acc, k - 1, ids.contig_uid(r.seq).hex(), len(r.seq), r.accession))
            bp += len(r.seq)
        if (k, bp) != (e["records"], e["bp"]) or bp != gsize[acc]:
            problems.append(f"{acc}: {k} records and {bp} bp; lock {e['records']} and {e['bp']}, GTDB {gsize[acc]}")
            continue
        tot["genomes"] += 1
        tot["records"] += k
        tot["bp"] += bp
        tot["fasta_bytes"] += n
        tot["gz_bytes"] += e["gz_bytes"]
    tmp = out.with_name(f".{out.name}.tmp")
    with open(tmp, "w") as w:
        w.write("accession\tordinal\tcontig_uid\tlength\trecord\n")
        for r in rows:
            w.write("\t".join(map(str, r)) + "\n")
    tmp.replace(out)
    rep = {"stage_dir": str(stage_dir), "genomes_in_lock": len(lock), "verified": tot,
           "distinct_contigs": len({r[2] for r in rows}), "problems": problems,
           "contigs_tsv_sha256": hashlib.sha256(out.read_bytes()).hexdigest()}
    print(json.dumps(rep, indent=1, sort_keys=True))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

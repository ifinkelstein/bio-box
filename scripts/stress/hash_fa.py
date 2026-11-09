"""Pod side: sha256 of each record's sequence in a FASTA written by samtools faidx (stdlib only).

    python3 hash_fa.py in.fa > out.tsv
"""
import hashlib
import sys

name, parts = None, []


def emit():
    if name is not None:
        print(name + "\t" + hashlib.sha256(b"".join(parts)).hexdigest())


with open(sys.argv[1], "rb") as fh:
    for ln in fh:
        ln = ln.rstrip(b"\r\n")
        if ln.startswith(b">"):
            emit()
            name, parts = ln[1:].decode().split("/")[0], []
        else:
            parts.append(ln)
emit()

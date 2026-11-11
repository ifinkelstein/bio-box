"""Read Pfam-A.hmm(.gz) (HMMER3/f) for card 0-09a: the model table, an HH-suite HHM database, and
accession lines for `hmmemit -a` output.

    python pfam_models.py models <Pfam-A.hmm[.gz]> <models.tsv>
        acc, name, leng, ga_seq, ga_dom (the GA line's two bit-score cutoffs), in file order.
    python pfam_models.py hhm <Pfam-A.hmm[.gz]> <out_db>
        One HH-suite HHM text per model, written as an MMseqs2 database (<out_db> data, entries separated
        by NUL; <out_db>.index "key<TAB>offset<TAB>length", keys 0..n-1 in file order), the input
        `mmseqs convertprofiledb` takes. Emissions and transitions are HMMER's, converted from -ln(p) to
        HH-suite's -1000*log2(p) (rounded; '*' stays '*'); the master and consensus sequence is HMMER's
        CONS column, upper-cased. NAME is "<acc> <name>", so MMseqs2 reports the accession.
    python pfam_models.py sto <Pfam-A.hmm[.gz]> < hmmemit.sto > named.sto
        Inserts "#=GF ID <name>" and "#=GF AC <acc>" into each alignment `hmmemit -a` writes (it writes
        none), so `mmseqs convertmsa` names each MSA by accession. The model is the first sequence's
        name without its "-sample<N>" suffix.

Deterministic: output depends only on the input bytes.
"""

from __future__ import annotations

import gzip
import math
import sys
from collections.abc import Iterator

AA = "ACDEFGHIKLMNPQRSTVWY"
LOG2E_1000 = 1000.0 / math.log(2.0)


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def records(path: str) -> Iterator[list[str]]:
    """The lines of each model, '//' excluded."""
    cur: list[str] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("//"):
                yield cur
                cur = []
            else:
                cur.append(line.rstrip("\n"))


def header(lines: list[str]) -> dict:
    h = {}
    for line in lines:
        if line.startswith("HMM "):
            break
        key = line[:6].strip()
        h.setdefault(key, line[6:].strip())
    return h


def _hh(v: str) -> str:
    return "*" if v == "*" else str(int(round(float(v) * LOG2E_1000)))


def body(lines: list[str]) -> tuple[list[str], list[str], list[tuple[str, list[str], list[str]]]]:
    """(null emissions, begin transitions, [(consensus residue, 20 match emissions, 7 transitions)])."""
    i = next(k for k, line in enumerate(lines) if line.startswith("HMM "))
    rest = lines[i + 2:]  # skip the "HMM A C ..." and "m->m ..." lines
    if rest and rest[0].split()[0] == "COMPO":
        rest = rest[1:]
    null = rest[0].split()
    begin = rest[1].split()
    states = []
    for k in range(2, len(rest), 3):
        m = rest[k].split()
        if len(m) < 21:
            break
        emis = m[1:21]
        cons = m[22] if len(m) > 22 else "X"
        trans = rest[k + 2].split()
        states.append((cons.upper() if cons.isalpha() else "X", emis, trans))
    return null, begin, states


def models(path: str, out: str) -> None:
    with open(out, "w") as fh:
        fh.write("acc\tname\tleng\tga_seq\tga_dom\n")
        for lines in records(path):
            h = header(lines)
            ga = h.get("GA", "nan nan;").rstrip(";").split()
            fh.write(f"{h['ACC']}\t{h['NAME']}\t{h['LENG']}\t{ga[0]}\t{ga[1]}\n")


def hhm_text(lines: list[str]) -> str:
    h = header(lines)
    null, begin, states = body(lines)
    leng = int(h["LENG"])
    if len(states) != leng:
        raise ValueError(f"{h['ACC']}: {len(states)} match states read, LENG {leng}")
    cons = "".join(c for c, _, _ in states)
    out = ["HHsearch 1.5", f"NAME  {h['ACC']} {h['NAME']}", "FAM   ",
           f"LENG  {leng} match states, {leng} columns in multiple alignment", "NEFF  1.0", "SEQ",
           ">Consensus", cons, f">{h['ACC']}", cons, "#",
           "NULL   " + "\t".join(_hh(v) for v in null) + "\t",
           "HMM    " + "\t".join(AA) + "\t",
           "       M->M\tM->I\tM->D\tI->M\tI->I\tD->M\tD->D\tNeff\tNeff_I\tNeff_D",
           "       " + "\t".join(_hh(v) for v in begin) + "\t*\t*\t*"]
    for i, (c, emis, trans) in enumerate(states, start=1):
        out.append(f"{c} {i}    " + "\t".join(_hh(v) for v in emis) + f"\t{i}")
        out.append("       " + "\t".join(_hh(v) for v in trans) + "\t1000\t0\t0")
        out.append("")
    out.append("//")
    return "\n".join(out) + "\n"


def hhm(path: str, out_db: str) -> None:
    off = 0
    with open(out_db, "wb") as data, open(out_db + ".index", "w") as idx:
        for key, lines in enumerate(records(path)):
            b = hhm_text(lines).encode() + b"\0"
            data.write(b)
            idx.write(f"{key}\t{off}\t{len(b)}\n")
            off += len(b)


def sto(path: str) -> None:
    acc = {header(lines)["NAME"]: header(lines)["ACC"] for lines in records(path)}
    pending = None
    for line in sys.stdin:
        if line.startswith("# STOCKHOLM"):
            sys.stdout.write(line)
            pending = True
            continue
        if line.strip() and not line.startswith("#") and not line.startswith("//"):
            name, seq = line.split()
            if pending:
                sys.stdout.write(f"#=GF ID {name.rsplit('-sample', 1)[0]}\n#=GF AC {acc[name.rsplit('-sample', 1)[0]]}\n")
                pending = None
            # hmmemit writes insert-state residues in lower case and '.' for the other rows' gaps in insert
            # columns; MMseqs2 reads the alignment as FASTA, where '.' is not a gap: write '-'.
            line = f"{name} {seq.replace('.', '-').upper()}\n"
        sys.stdout.write(line)


def main(argv: list[str]) -> int:
    if len(argv) >= 2 and argv[0] == "models" and len(argv) == 3:
        models(argv[1], argv[2])
    elif len(argv) == 3 and argv[0] == "hhm":
        hhm(argv[1], argv[2])
    elif len(argv) == 2 and argv[0] == "sto":
        sto(argv[1])
    else:
        print(__doc__, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

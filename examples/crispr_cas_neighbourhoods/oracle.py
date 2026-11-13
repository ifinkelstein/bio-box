"""The hmmsearch oracle of the CRISPR-Cas example (standard library only; run in the image with `python -I`).

    python -I oracle.py ANCHORS_HMM TBLOUT NAME=HITS_TSV [NAME=HITS_TSV ...] > observed_oracle.tsv

ANCHORS_HMM is the HMMER3 file the runs searched with (one model per NAME, matched by the model's NAME line);
TBLOUT is `hmmsearch --tblout` of that file against the proteins of `sql/oracle_uids.sql` (every hit of the runs
and every protein the release's Pfam scan gave either model), run with `-Z 1e8 --domZ 1e8 -T 0`. For each NAME the
proteins hmmsearch scores at or above the model's sequence gathering cutoff (GA, first number) are compared with the
proteins of the run's `hits.tsv`. Writes TSV rows `check_id<TAB>value`:

    oracle_<name>_only_search      proteins in hits.tsv that hmmsearch scores below GA
    oracle_<name>_only_hmmsearch   proteins hmmsearch scores at or above GA that are not in hits.tsv
    oracle_<name>_max_bits_diff    the largest |hits.tsv seq_bitscore - hmmsearch score| over common proteins
"""

from __future__ import annotations

import csv
import sys


def gathering(path: str) -> dict[str, float]:
    """{model NAME: sequence GA} of an HMMER3 file."""
    out: dict[str, float] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split()[1]
            elif line.startswith("GA ") and name is not None:
                out[name] = float(line.split()[1].rstrip(";"))
    return out


def tblout(path: str) -> dict[str, dict[str, float]]:
    """{model NAME: {target: full-sequence score}} of `hmmsearch --tblout`."""
    out: dict[str, dict[str, float]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            out.setdefault(f[2], {})[f[0]] = float(f[5])
    return out


def run_hits(path: str) -> dict[str, float]:
    """{protein_uid: seq_bitscore} of a run's hits.tsv (one row per gene; a protein's rows share the score)."""
    out: dict[str, float] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            out[r["protein_uid"]] = float(r["seq_bitscore"])
    return out


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    ga = gathering(argv[0])
    scores = tblout(argv[1])
    rows = []
    for spec in argv[2:]:
        name, _, path = spec.partition("=")
        model = {"tniq": "TniQ", "cas7f": "Cas_Csy3"}.get(name, name)
        if model not in ga:
            print(f"oracle: no model {model} with a GA line in {argv[0]}", file=sys.stderr)
            return 2
        ours = run_hits(path)
        theirs = {u: s for u, s in scores.get(model, {}).items() if s >= ga[model]}
        both = ours.keys() & theirs.keys()
        diff = max((abs(ours[u] - theirs[u]) for u in both), default=0.0)
        rows += [(f"oracle_{name}_only_search", str(len(ours.keys() - theirs.keys()))),
                 (f"oracle_{name}_only_hmmsearch", str(len(theirs.keys() - ours.keys()))),
                 (f"oracle_{name}_max_bits_diff", f"{diff:.2f}")]
    print("check_id\tvalue")
    for k, v in rows:
        print(f"{k}\t{v}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

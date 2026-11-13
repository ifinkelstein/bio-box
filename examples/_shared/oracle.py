"""The hmmsearch oracle of a worked example: a run's hits against HMMER run by hand (standard library only; run in
the image with `python -I`; card E-11).

    python -I oracle.py --hmm MODELS_HMM --tblout TBLOUT --prefix NAME --model MODEL --bits BITS|ga HITS_TSV
        > observed_oracle.tsv

MODELS_HMM is the HMMER3 file the run searched with (MODEL is matched by its NAME line); TBLOUT is `hmmsearch
--tblout` of that file against the proteins the example chose (every hit of the run, and every protein the
release's Pfam scan gave the model), run with `-Z 1e8 --domZ 1e8 -T 0`, so the scores are those of the search at
the fixed reference Z. The proteins hmmsearch scores at or above the cutoff (`--bits`: a number of bits, or `ga`
for the model's sequence gathering cutoff) are compared with the proteins of the run's `hits.tsv`. Writes TSV rows
`check_id<TAB>value`:

    oracle_<NAME>_only_search      proteins in hits.tsv that hmmsearch scores below the cutoff
    oracle_<NAME>_only_hmmsearch   proteins hmmsearch scores at or above the cutoff that are not in hits.tsv
    oracle_<NAME>_max_bits_diff    the largest |hits.tsv seq_bitscore - hmmsearch score| over common proteins

The CRISPR-Cas example keeps its own `oracle.py` (two models at their gathering cutoffs, written before this one);
the two compare the same way.
"""

from __future__ import annotations

import argparse
import csv
import sys


def gathering(path: str) -> dict[str, float | None]:
    """{model NAME: sequence GA} of an HMMER3 file (None for a model without a GA line)."""
    out: dict[str, float | None] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split()[1]
                out[name] = None
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


def cutoff(text: str) -> str | float:
    if text == "ga":
        return text
    try:
        return float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{text!r} is neither a number of bits nor 'ga'") from None


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--tblout", required=True)
    ap.add_argument("--prefix", required=True, help="the name in the check ids: oracle_<prefix>_...")
    ap.add_argument("--model", required=True)
    ap.add_argument("--bits", required=True, type=cutoff, help="the inclusion cutoff in bits, or 'ga'")
    ap.add_argument("hits")
    a = ap.parse_args(argv)
    models = gathering(a.hmm)
    if a.model not in models:
        print(f"oracle: no model {a.model} in {a.hmm}", file=sys.stderr)
        return 2
    bits = models[a.model] if a.bits == "ga" else a.bits
    if bits is None:
        print(f"oracle: model {a.model} has no GA line in {a.hmm}", file=sys.stderr)
        return 2
    ours = run_hits(a.hits)
    theirs = {u: s for u, s in tblout(a.tblout).get(a.model, {}).items() if s >= bits}
    both = ours.keys() & theirs.keys()
    diff = max((abs(ours[u] - theirs[u]) for u in both), default=0.0)
    print("check_id\tvalue")
    print(f"oracle_{a.prefix}_only_search\t{len(ours.keys() - theirs.keys())}")
    print(f"oracle_{a.prefix}_only_hmmsearch\t{len(theirs.keys() - ours.keys())}")
    print(f"oracle_{a.prefix}_max_bits_diff\t{diff:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

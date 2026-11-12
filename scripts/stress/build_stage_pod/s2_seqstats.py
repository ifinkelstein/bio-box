"""Characterise a FASTA source without changing it, and count its sequences already in other stages (S2 inputs).

    python s2_seqstats.py FASTA STAGE_UIDS_TSV OUT_DIR

FASTA is read with `biobox.fasta.iter_records` (gzip by the `.gz` suffix), as ingest reads a source
object. Per record: the store's content key `ids.contig_uid` (BLAKE3-128 of the uppercase sequence), the
same of the reverse complement, `ids.raw_digest` (sha256 of the exact bytes), sha256 of the uppercase
sequence and of its reverse complement, and base counts. STAGE_UIDS_TSV is `s2_stage_uids.py`'s output.

Writes into OUT_DIR (each under a temporary name, then renamed):
- `records.tsv`: one row per record, in file order;
- `repeats.tsv`: every record whose uppercase sequence repeats an earlier record's (in file order),
  the records `make_inputs.py` does not list again, with the first record of its sequence;
- `overlap.tsv`: every (record, stage, how) where the record's uid, or its reverse complement's, is a
  contig uid of that stage;
- `summary.json`, also printed.
"""

import hashlib
import json
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

from biobox import fasta, ids

COMP = bytes.maketrans(b"ACGTRYKMBVDHNSW", b"TGCAYRMKVBHDNSW")
IUPAC_OTHER = b"RYKMBVDHSW"
LOWER = bytes(range(ord("a"), ord("z") + 1))
KEYWORDS = ("circular", "linear", "complete", "partial", "topology", "host", "sample", "plasmid", "phage",
            "chromosome", "=", "|", ";", ",")


def revcomp(up: bytes) -> bytes:
    return up.translate(COMP)[::-1]


def write(path: Path, text: str) -> None:
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    src, stage_tsv, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    stages: dict[str, set[bytes]] = {}
    with open(stage_tsv) as fh:
        fh.readline()
        for ln in fh:
            s, how, u, _, _ = ln.rstrip("\n").split("\t")
            stages.setdefault(f"{s}/{how}", set()).add(bytes.fromhex(u))
    recs, first_of = [], {}
    lines = ["ordinal\taccession\theader\tlength\tcontig_uid\trc_contig_uid\traw_sha256\tupper_sha256\t"
             "rc_upper_sha256\tA\tC\tG\tT\tN\tiupac_other\tnon_iupac\tlowercase"]
    repeats = ["ordinal\taccession\tcontig_uid\tfirst_ordinal\tfirst_accession"]
    overlap = ["ordinal\taccession\tcontig_uid\tstage\tmatch"]
    kw, ntok, shapes = Counter(), Counter(), Counter()
    tot = Counter()
    for i, r in enumerate(fasta.iter_records(src)):
        up = r.seq.upper()
        rc = revcomp(up)
        uid, rcu = ids.contig_uid(r.seq), ids.contig_uid(rc)
        cnt = {b: up.count(b.encode()) for b in "ACGTN"}
        other = sum(up.count(bytes([c])) for c in IUPAC_OTHER)
        non = len(up) - sum(cnt.values()) - other
        lower = len(r.seq) - len(r.seq.translate(None, LOWER))
        tot.update({**cnt, "iupac_other": other, "non_iupac": non, "lowercase": lower, "bp": len(r.seq),
                    "records_with_non_acgt": int(cnt["N"] + other + non > 0), "records_with_n": int(cnt["N"] > 0)})
        lines.append("\t".join(map(str, [i, r.accession, r.header, len(r.seq), uid.hex(), rcu.hex(),
                                         ids.raw_digest(r.seq).hex(), hashlib.sha256(up).hexdigest(),
                                         hashlib.sha256(rc).hexdigest(), cnt["A"], cnt["C"], cnt["G"], cnt["T"],
                                         cnt["N"], other, non, lower])))
        if uid in first_of:
            j = first_of[uid]
            repeats.append(f"{i}\t{r.accession}\t{uid.hex()}\t{j}\t{recs[j][0]}")
        else:
            first_of[uid] = i
        canon = min(uid, rcu)
        for name, uids in stages.items():
            if uid in uids:
                overlap.append(f"{i}\t{r.accession}\t{uid.hex()}\t{name}\tsame")
            elif rcu in uids:
                overlap.append(f"{i}\t{r.accession}\t{uid.hex()}\t{name}\treverse_complement")
        recs.append((r.accession, len(r.seq), uid, canon))
        ntok[len(r.header.split())] += 1
        shapes[re.sub(r"\d+", "#", r.header.split()[0] if r.header else "")[:60]] += 1
        low = r.header.lower()
        kw.update(k for k in KEYWORDS if k in low)
    write(out / "records.tsv", "\n".join(lines) + "\n")
    write(out / "repeats.tsv", "\n".join(repeats) + "\n")
    write(out / "overlap.tsv", "\n".join(overlap) + "\n")
    lens = sorted(n for _, n, _, _ in recs)
    by_uid, by_canon = Counter(u for _, _, u, _ in recs), Counter(c for _, _, _, c in recs)
    ov: dict[str, dict[str, int]] = {}
    for row in overlap[1:]:
        _, _, u, name, how = row.split("\t")
        d = ov.setdefault(name, {"same": 0, "reverse_complement": 0, "distinct_uids_same": set()})
        d[how] += 1
        if how == "same":
            d["distinct_uids_same"].add(u)
    for d in ov.values():
        d["distinct_uids_same"] = len(d["distinct_uids_same"])
    accs = Counter(a for a, _, _, _ in recs)
    rep = {
        "source": str(src), "records": len(recs), "bp": tot["bp"],
        "length": {"min": lens[0], "median": statistics.median(lens), "max": lens[-1],
                   "mean": round(tot["bp"] / len(recs), 1),
                   "under_200000": sum(1 for n in lens if n < 200_000),
                   "deciles": [lens[int(len(lens) * k / 10)] for k in range(1, 10)]},
        "repeats": {"distinct_contig_uids": len(by_uid),
                    "uids_seen_more_than_once": sum(1 for v in by_uid.values() if v > 1),
                    "records_repeating_an_earlier_record": len(repeats) - 1,
                    "distinct_by_either_strand": len(by_canon),
                    "records_repeating_an_earlier_record_on_either_strand": len(recs) - len(by_canon),
                    "distinct_sequences_whose_reverse_complement_is_another": len(by_uid) - len(by_canon)},
        "bases": {k: tot[k] for k in ("A", "C", "G", "T", "N", "iupac_other", "non_iupac", "lowercase")},
        "records_with_n": tot["records_with_n"], "records_with_non_acgt": tot["records_with_non_acgt"],
        "headers": {"tokens_per_header": dict(sorted(ntok.items())), "first_token_shapes": dict(shapes.most_common(10)),
                    "distinct_shapes": len(shapes), "keyword_counts": dict(sorted(kw.items())),
                    "accessions_distinct": len(accs), "accessions_repeated": sum(1 for v in accs.values() if v > 1)},
        "overlap": {name: ov.get(name, {"same": 0, "reverse_complement": 0, "distinct_uids_same": 0})
                    for name in sorted(stages)},
        "stage_sets": {name: len(v) for name, v in sorted(stages.items())},
    }
    text = json.dumps(rep, indent=1)
    write(out / "summary.json", text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

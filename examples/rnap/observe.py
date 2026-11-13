"""The RNA polymerase example's values read from the files its run wrote, not from SQL (standard library only; run
in the image with `python -I`; card E-11).

    python -I observe.py --hmm HMM --hits HITS_TSV --neighbors NEIGHBORS_TSV --flanks FLANKS_DIR \
        --control FEATURE_UID > observed_files.tsv

HMM is the model searched with; HITS_TSV and NEIGHBORS_TSV are `search hmm`'s `hits.tsv` and `neighbors.tsv`;
FLANKS_DIR is an `nbr flanks --out` directory of the same run (`groups.tsv`, `context_types.tsv`); FEATURE_UID is
the control gene. Writes TSV rows `check_id<TAB>value`:

    rnap_hmm_sha256               sha256 of the model file
    rnap_summary_loci, rnap_summary_contigs, rnap_summary_censored_loci, rnap_summary_neighbours
                                  the header of neighbors.tsv (n_anchor_loci, n_contigs, n_censored_loci, n_neighbors)
    rnap_top_pfam_neighbour, rnap_top_pfam_neighbour_loci
                                  the Pfam family next to most loci (the `unannotated` row is not a family; a tie goes
                                  to the least key) and its loci
    rnap_flank_groups, rnap_flank_leaves, rnap_flank_context_types
                                  groups.tsv: groups, leaves and context types in all
    rnap_flank_types_mixed_kind   context types (group, type) holding loci of a virus and of a non-virus genome
    rnap_flank_loci_missing       hit genes of hits.tsv with no row in context_types.tsv
    t7_flank_group_loci, t7_flank_group_types, t7_flank_type_kinds
                                  the control's group (its loci and context types) and the genome kinds of its context
                                  type; left out when context_types.tsv has no row for the control
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


def table(path: Path) -> list[dict[str, str]]:
    """The rows of a TSV whose first lines may be `# key<TAB>value` comments."""
    with open(path, newline="") as fh:
        lines = [line for line in fh.read().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE))


def summary(path: Path) -> tuple[dict[str, object], list[dict[str, str]]]:
    """neighbors.tsv: its `# key<TAB>json` header, and the rows of its pfam table."""
    head: dict[str, object] = {}
    blocks: list[list[str]] = [[]]
    with open(path, newline="") as fh:
        for line in fh.read().splitlines():
            if line.startswith("# "):
                k, _, v = line[2:].partition("\t")
                head[k] = json.loads(v)
            elif not line.strip():
                blocks.append([])
            else:
                blocks[-1].append(line)
    pfam: list[dict[str, str]] = []
    for b in blocks:
        rows = list(csv.DictReader(b, delimiter="\t", quoting=csv.QUOTE_NONE))
        pfam += [r for r in rows if r.get("level") == "pfam"]
    return head, pfam


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--neighbors", required=True)
    ap.add_argument("--flanks", required=True)
    ap.add_argument("--control", required=True)
    a = ap.parse_args(argv)
    rows: list[tuple[str, str]] = [("rnap_hmm_sha256", hashlib.sha256(Path(a.hmm).read_bytes()).hexdigest())]

    head, pfam = summary(Path(a.neighbors))
    rows += [("rnap_summary_loci", str(head["n_anchor_loci"])), ("rnap_summary_contigs", str(head["n_contigs"])),
             ("rnap_summary_censored_loci", str(head["n_censored_loci"])),
             ("rnap_summary_neighbours", str(head["n_neighbors"]))]
    named = [r for r in pfam if r["unannotated"] != "true"]
    if named:
        top = min(named, key=lambda r: (-int(r["n_loci"]), r["key"]))
        rows += [("rnap_top_pfam_neighbour", top["model_name"]), ("rnap_top_pfam_neighbour_loci", top["n_loci"])]

    flanks = Path(a.flanks)
    groups = table(flanks / "groups.tsv")
    ctx = table(flanks / "context_types.tsv")
    kinds: dict[tuple[str, str], set[str]] = {}
    for r in ctx:
        if r["context_type"]:
            kinds.setdefault((r["group"], r["context_type"]), set()).add(r["kind"])
    mixed = sum(1 for k in kinds.values() if "virus" in k and len(k) > 1)
    hit_genes = {r["feature_uid"] for r in table(Path(a.hits))}
    rows += [("rnap_flank_groups", str(len(groups))),
             ("rnap_flank_leaves", str(sum(int(g["n_leaves"]) for g in groups))),
             ("rnap_flank_context_types", str(sum(int(g["n_context_types"]) for g in groups))),
             ("rnap_flank_types_mixed_kind", str(mixed)),
             ("rnap_flank_loci_missing", str(len(hit_genes - {r["feature_uid"] for r in ctx})))]

    mine = [r for r in ctx if r["feature_uid"] == a.control]
    if mine:
        g = next(x for x in groups if x["group"] == mine[0]["group"])
        rows += [("t7_flank_group_loci", g["n_loci"]), ("t7_flank_group_types", g["n_context_types"]),
                 ("t7_flank_type_kinds", ",".join(sorted(kinds.get((g["group"], mine[0]["context_type"]), set()))))]

    print("check_id\tvalue")
    for k, v in rows:
        print(f"{k}\t{v}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

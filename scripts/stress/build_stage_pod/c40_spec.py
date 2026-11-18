"""Card C-40: the spec of `make_inputs.py plasmid-sources` for stage PL2.

    python c40_spec.py DL FASTA_DIR OUT_JSON

DL is the fetch directory (`gut/`, `plasmidscope/`, `SUMS` as c40_fetch.sh wrote them); FASTA_DIR holds
`<SRC>/<SRC>.fasta` from c40_concat.py for each PlasmidScope repository and `gut/predicted_plasmids_summary.txt`, the
gut table decompressed. The sources, in the order a record's genome is
first taken (the first source wins a sequence that several carry):

1. gut_plasmids_yu2024: predicted_plasmids.fa.gz (226,194 PlasX predictions) with predicted_plasmids_summary.txt.gz;
   only the 68,350 non-redundant representatives are taken (`in_dereplicated_plasmids` True). Topology from `circular`
   (True: detection in the sample's own assembly >= 0.95 and a reverse-forward read pair across the contig's ends):
   True circular, False linear. Completeness from the same column: True complete, False unknown (the table has no
   completeness; a representative is never a MobMess fragment cluster, `cluster_type` fragment has no representative).
2. PlasmidScope's metagenomic repositories, IMG-PR then mMGE, then its isolate repositories RefSeq, GenBank, ENA, DDBJ,
   COMPASS, TPA, Kraken2. Its PLSDB file is not taken: stage PL1 holds PLSDB 2024_05_31_v2. Topology from `Topology`:
   circular is circular; linear, inverted terminal repeat, concatemer and direct terminal repeat are linear (a direct
   terminal repeat is a circular molecule whose contig carries the repeat at both ends: wrapping it would join two
   copies of the repeat, so it is called linear; c40_analyze.py measures how many keep both copies). Completeness from
   `Completeness`: complete is complete, incomplete is fragment, '-' is unknown. Records of the FASTA without a row in
   the repository's table are skipped (the table is PlasmidScope's curated set; its archives hold a few more files).

Writes OUT_JSON and prints it.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

SCOPE_ORDER = ("IMG-PR", "mMGE", "RefSeq", "GenBank", "ENA", "DDBJ", "COMPASS", "TPA", "Kraken2")
SCOPE_TOPOLOGY = {"circular": "circular", "linear": "linear", "inverted terminal repeat": "linear",
                  "concatemer": "linear", "direct terminal repeat": "linear"}
SCOPE_COMPLETENESS = {"complete": "complete", "incomplete": "fragment", "-": "unknown"}
SCOPE_LICENSE = ("PlasmidScope states no licence; its article (CC BY 4.0) says all data are freely available without "
                 "registration (doi:10.1093/nar/gkae930)")
SCOPE_TERMS = {"IMG-PR": "IMG/PR records (JGI data use policy) through PlasmidScope's open download; " + SCOPE_LICENSE,
               "mMGE": "mMGE records, no terms stated; " + SCOPE_LICENSE}
NOT_REDISTRIBUTABLE = {"IMG-PR", "mMGE"}


def fetched_at(sums: Path) -> dict[str, str]:
    """file (relative to DL) -> the UTC time its fetch ended, from c40_fetch.sh's SUMS lines."""
    out = {}
    for ln in sums.read_text().splitlines():
        parts = ln.split()
        out[parts[0]] = parts[-1]
    return out


def spec(dl: Path, fasta_dir: Path) -> dict:
    t = fetched_at(dl / "SUMS")
    srcs = [{"name": "gut_plasmids_yu2024", "source_name": "gut_plasmids_yu2024", "source_release": "zenodo-5843600-v1.2",
             "license": "CC BY 4.0 (Zenodo record 5843600)", "redistributable": True,
             "retrieved_at": t["gut/predicted_plasmids.fa.gz"], "fasta": [str(dl / "gut/predicted_plasmids.fa.gz")],
             "table": str(fasta_dir / "gut/predicted_plasmids_summary.txt"), "id_column": "contig",
             "topology_column": "circular", "topology": {"True": "circular", "False": "linear"},
             "completeness_column": "circular", "completeness": {"True": "complete", "False": "unknown"},
             "select": {"column": "in_dereplicated_plasmids", "values": ["True"]}}]
    for s in SCOPE_ORDER:
        srcs.append({"name": f"plasmidscope_{s}", "source_name": "plasmidscope",
                     "source_release": f"{s}.download-{t[f'plasmidscope/{s}.fasta.tar.gz'][:10]}",
                     "license": SCOPE_TERMS.get(s, f"{s} records (INSDC or the repository's own); " + SCOPE_LICENSE),
                     "redistributable": s not in NOT_REDISTRIBUTABLE,
                     "retrieved_at": t[f"plasmidscope/{s}.fasta.tar.gz"],
                     "fasta": [str(fasta_dir / s / f"{s}.fasta")],
                     "table": str(dl / f"plasmidscope/{s}.plasmid_list.download.tsv"), "id_column": "Plasmid_ID",
                     "topology_column": "Topology", "topology": SCOPE_TOPOLOGY,
                     "completeness_column": "Completeness", "completeness": SCOPE_COMPLETENESS,
                     "without_row": "skip"})
    return {"card": "C-40", "stage": "PL2", "sources": srcs}


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    s = spec(Path(argv[0]), Path(argv[1]))
    Path(argv[2]).write_text(json.dumps(s, indent=1))
    print(json.dumps(s))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

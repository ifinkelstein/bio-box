"""Card C-39 step 6: what VIRE would add over UHGV and the soil atlas, from VIRE's genome table alone (stdlib only).

    python c39_vire.py VIRE_METADATA_TXT OUT

VIRE_METADATA_TXT is the response of VIRE's all-samples metadata URL (sources.yaml `vire_tables`), a TSV with one
row per genome. Read as data only. Writes OUT/summary.json (and prints it):

- genomes, species clusters (`species_cluster`, 95% ANI in the paper), genus clusters, bases; CheckV tiers
  (`checkv_checkv_quality`), `genomad_topology`, `genomad_genetic_code`, `source`;
- habitat from `microntology` (semicolon-separated labels): every label's genome count, and each genome put in one
  of: human gut (the label `...:mammalian host:human host` and a label `...:digestive tract:intestine`), soil (a
  label naming soil), other; a cluster is human gut when any member is, soil when any member is, else other;
- `alignment_to_imgvr` (VIRE's own flag, True when the genome aligns to IMG/VR): genomes and clusters with no
  member aligned, by habitat;
- for clusters, the bases of their longest member (an *estimate* of a representative's length; VIRE's table does
  not name representatives) summed by habitat, to size a stage of VIRE's species representatives.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

csv.field_size_limit(1 << 30)


def habitat(labels: list[str]) -> str:
    low = [x.lower() for x in labels]
    human = any(x.endswith(":human host") for x in low)  # the labels name host and body site separately
    if human and any("digestive tract:intestine" in x for x in low):
        return "human_gut"
    if any("soil" in x for x in low):
        return "soil"
    return "other"


def main(argv: list[str]) -> int:
    src, out = Path(argv[0]), Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    n = bases = 0
    quality, topo, code, source, label_n, hab_n = Counter(), Counter(), Counter(), Counter(), Counter(), Counter()
    imgvr, imgvr_by_hab = Counter(), Counter()
    clusters: dict[str, list] = {}  # species cluster -> [habitat set, any aligned to IMG/VR, longest, members, best q]
    genera: set[str] = set()
    with open(src, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n += 1
            ln = int(r["contig_length"])
            bases += ln
            quality[r["checkv_checkv_quality"]] += 1
            topo[r["genomad_topology"]] += 1
            code[r["genomad_genetic_code"]] += 1
            source[r["source"]] += 1
            labels = [x for x in (r["microntology"] or "").split(";") if x and x != "None"]
            label_n.update(set(labels))
            h = habitat(labels)
            hab_n[h] += 1
            al = r["alignment_to_imgvr"] == "True"
            imgvr[r["alignment_to_imgvr"]] += 1
            imgvr_by_hab[(h, al)] += 1
            genera.add(r["genus_cluster"])
            c = clusters.get(r["species_cluster"])
            if c is None:
                clusters[r["species_cluster"]] = [{h}, al, ln, 1, r["checkv_checkv_quality"]]
            else:
                c[0].add(h)
                c[1] = c[1] or al
                c[2] = max(c[2], ln)
                c[3] += 1
    chab = Counter()
    cnew = Counter()
    cbases = Counter()
    cbases_new = Counter()
    single = 0
    for hs, al, longest, members, _ in clusters.values():
        h = "human_gut" if "human_gut" in hs else "soil" if "soil" in hs else "other"
        chab[h] += 1
        cbases[h] += longest
        single += members == 1
        if not al:
            cnew[h] += 1
            cbases_new[h] += longest
    res = {
        "genomes": n, "bases": bases, "species_clusters": len(clusters), "genus_clusters": len(genera),
        "species_clusters_of_one_genome": single,
        "checkv_quality": dict(quality.most_common()), "genomad_topology": dict(topo.most_common()),
        "genomad_genetic_code": dict(code.most_common()), "source": dict(source.most_common()),
        "habitat_genomes": dict(hab_n), "habitat_species_clusters": dict(chab),
        "alignment_to_imgvr_genomes": dict(imgvr),
        "alignment_to_imgvr_by_habitat": {f"{h}:{a}": k for (h, a), k in sorted(imgvr_by_hab.items())},
        "species_clusters_no_member_aligned_to_imgvr": dict(cnew),
        "longest_member_bases_by_habitat": dict(cbases),
        "longest_member_bases_no_member_aligned_by_habitat": dict(cbases_new),
        "labels_top60": dict(label_n.most_common(60)),
        "labels_with_soil": {k: v for k, v in label_n.items() if "soil" in k.lower()},
        "labels_with_human": {k: v for k, v in label_n.items() if "human" in k.lower()},
    }
    (out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

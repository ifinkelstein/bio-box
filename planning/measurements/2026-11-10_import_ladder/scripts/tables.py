"""Per-table new rows/bytes per release and the query timings, from results_pod/<KEY>/result.json."""

import json
import sys
from pathlib import Path

d = Path(sys.argv[1])
for key in ("A", "B", "Bp"):
    R = json.loads((d / key / "result.json").read_text())
    s = R["steps"]
    c = s["catalog"]["result"]
    print(f"\n## {key} {R['release']}: files {c['files']} artifacts {c['artifacts']} total_bytes {c['total_bytes']}")
    for t, v in c["tables"].items():
        if v["new_files"]:
            print(f"  {t:28s} new_files {v['new_files']:3d} new_rows {v['new_rows']:>11,d} new_bytes {v['new_bytes']:>14,d}"
                  f"  (all: rows {v['rows']:,d} bytes {v['bytes']:,d})")
    q = s["queries"]["result"]
    print(f"  open: SeqStore {q['seqstore_open_s']} s +{q['seqstore_rss_growth_mb']} MB; open_catalog {q['open_catalog_s']} s "
          f"+{q['open_catalog_rss_growth_mb']} MB; row counts of {len(q['row_counts']['views'])} views {q['row_counts']['wall_s']} s")
    g = q["seq_get_1000"]
    print(f"  seq_get 1000: {g['wall_s']} s, {g['result']}, contigs {g['distinct_contigs']}, blocks {g['bgzf']['blocks_read']}, "
          f"index rg {g['index']}, warm {q['seq_get_1000_warm']['wall_s']} s")
    t = q["tapes_of_one"]
    print(f"  tapes_of: {t['wall_s']} s scans {[(x['row_groups_scanned'], x['row_groups_total']) for x in t['scans']]}")
    p = q["postings_lookup_nr50"]
    print(f"  postings nr50 {p['cluster_uid']}: {p['wall_s']} s rows {p['result']} scans "
          f"{[(x['row_groups_scanned'], x['row_groups_total']) for x in p['scans']]} rep_pfam {p.get('rep_pfam')}")
    for n in q["neighbours"]:
        print(f"  nbr {n['label']}: anchors {n['anchor_loci']} contigs {n['anchor_contigs']} with_partner {n['loci_with_partner']} "
              f"pairs {n['pair_rows']} {n['wall_s']} s scans {[(x['row_groups_scanned'], x['row_groups_total']) for x in n['scans']]} "
              f"top3 {n['top10_partners'][:3]}")
    if "prok_top_nr50" in q:
        print(f"  prok_top_nr50 {q['prok_top_nr50']}")

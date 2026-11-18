"""Card C-41: analyse GlobDB release 226's metadata tables, dictionaries and GTDB R226 crosswalk, and propose stages.

    python c41_tables.py analyze --dl DL_DIR --gtdb-tar GTDB_TAR --root ROOT \
        --p16-sel ALL_TSV --c16-stages STAGES_TSV --pl2-inputs PL2_DIR --v1-inputs V1_DIR --out OUT_DIR
    python c41_tables.py probe-tar --url URL --bytes N
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io as pyio
import json
import math
import re
import sys
import tarfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

GBP = 10_000_000_000
GC_RE = re.compile(r"(GC[AF]_(\d{9})(?:\.\d+)?)")

# Order of sequential dereplication in GlobDB r226 (website_pages/globdb_r226_methods and globdb_r226_citations).
RESOURCE_ORDER = (
    "GTDB",
    "mOTU",
    "RBG",
    "SPIRE",
    "GEM",
    "MGnify",
    "GOMC",
    "SMAG",
    "TPMC",
    "MRGM/HRGM2",
    "cFMD",
    "SHGO",
    "AMXMAG",
    "GFS",
)


def classify_id(gid: str) -> tuple[str, str]:
    """Return (resource_14, sub_source) from a GlobDB r226 genome ID."""
    if gid.startswith(("GB_GCA_", "RS_GCF_", "GCA_", "GCF_")):
        return "GTDB", "GTDB"
    if gid.startswith("MOTU40_") or gid.startswith("MOTU"):
        return "mOTU", "mOTU"
    if gid.startswith("BCRBG_"):
        return "RBG", "RBG"
    if gid.startswith("SPIREOTU_"):
        return "SPIRE", "SPIREOTU"
    if gid.startswith("SPECIV4_"):
        return "SPIRE", "SPECIV4"
    if gid.startswith("GEMOTU_"):
        return "GEM", "GEM"
    if gid.startswith("MGYG"):
        return "MGnify", "MGnify"
    if gid.startswith("GOMCOTU_"):
        return "GOMC", "GOMC"
    if gid.startswith("SMAGOTU_"):
        return "SMAG", "SMAG"
    if gid.startswith("TPMCOTU_"):
        return "TPMC", "TPMC"
    if gid.startswith("HRGM2_") or gid.startswith("HRGM"):
        return "MRGM/HRGM2", "HRGM2"
    if gid.startswith("MRGM_") or gid.startswith("MRGM"):
        return "MRGM/HRGM2", "MRGM"
    if gid.startswith("FMDMAG_"):
        return "cFMD", "cFMD"
    if gid.startswith("SHGOMAG_"):
        return "SHGO", "SHGO"
    if gid.startswith("AMXMAG_"):
        return "AMXMAG", "AMXMAG"
    if gid.startswith("GFSMAG_"):
        return "GFS", "GFS"
    return f"OTHER:{gid.split('_')[0]}", gid.split("_")[0]


def quantiles(vals: list[float]) -> dict:
    if not vals:
        return {}
    s = sorted(vals)
    n = len(s)

    def q(p: float) -> float:
        idx = p * (n - 1)
        lo = int(math.floor(idx))
        hi = int(math.ceil(idx))
        return s[lo] + (s[hi] - s[lo]) * (idx - lo)

    return {
        "count": n,
        "min": round(s[0], 4),
        "p01": round(q(0.01), 4),
        "p05": round(q(0.05), 4),
        "p10": round(q(0.10), 4),
        "p25": round(q(0.25), 4),
        "p50": round(q(0.50), 4),
        "p75": round(q(0.75), 4),
        "p90": round(q(0.90), 4),
        "p95": round(q(0.95), 4),
        "p99": round(q(0.99), 4),
        "max": round(s[-1], 4),
        "mean": round(sum(s) / n, 4),
    }


def load_dictionaries(dl: Path) -> tuple[dict[str, str], dict[str, dict]]:
    """Return (globdb_id -> original_id, per_file_summary)."""
    mapping: dict[str, str] = {}
    summary: dict[str, dict] = {}
    with tarfile.open(dl / "globdb_r226_dictionaries.tar.gz", "r:gz") as tf:
        for m in sorted(tf.getmembers(), key=lambda x: x.name):
            if not m.isfile() or m.name.endswith("/README"):
                continue
            fh = tf.extractfile(m)
            assert fh is not None
            text = pyio.TextIOWrapper(fh, encoding="utf-8")
            rd = csv.reader(text, delimiter="\t")
            rows = list(rd)
            name = Path(m.name).name
            has_header = bool(rows and rows[0][0] == "GlobDB_ID")
            start = 1 if has_header else 0
            cnt = 0
            gca_hits = 0
            bad_rows = 0
            two_col_rows = 0
            two_col_samples: list[list[str]] = []
            for r in rows[start:]:
                if not r or all(x == "" for x in r):
                    continue
                if has_header:
                    if len(r) >= 3 and r[2] != "":
                        gid, orig = r[0], r[2]
                    elif len(r) >= 2:
                        gid, orig = r[0], r[1]
                        two_col_rows += 1
                        if len(two_col_samples) < 5:
                            two_col_samples.append(r)
                    else:
                        bad_rows += 1
                        continue
                else:
                    if len(r) < 2:
                        bad_rows += 1
                        continue
                    orig, gid = r[0], r[1]
                cnt += 1
                mapping[gid] = orig
                if GC_RE.search(orig):
                    gca_hits += 1
            summary[name] = {
                "rows": cnt,
                "has_header": has_header,
                "orig_with_gca_gcf": gca_hits,
                "two_col_rows": two_col_rows,
                "two_col_samples": two_col_samples,
                "malformed_rows": bad_rows,
            }
    return mapping, summary


def load_gtdb(gtdb_tar: Path) -> tuple[dict[str, dict], dict[str, str], dict]:
    """Read GTDB R226 bac120 + ar53 metadata from the tar on /stor (streaming, bounded memory).

    Returns:
      reps: gtdb_accession (e.g. GB_GCA_... or RS_GCF_...) -> selected fields for all gtdb_representative == 't'
      all_stems: 9-digit numeric GCA/GCF stem -> gtdb_accession (for all 732,475 GTDB R226 genomes)
      counts: summary counts over all GTDB R226 rows
    """
    reps: dict[str, dict] = {}
    all_acc_to_rep: dict[str, tuple[str, str]] = {}  # stripped GCA_/GCF_ -> (gtdb_acc, gtdb_genome_representative)
    stem_to_gtdb: dict[str, list[str]] = defaultdict(list)
    total_rows = 0
    by_dom_total = Counter()
    by_dom_reps = Counter()
    reps_under_2mbp = 0
    reps_comp_lt_40 = 0
    reps_cont_gt_10 = 0
    keep_cols = (
        "accession",
        "gtdb_representative",
        "gtdb_genome_representative",
        "gtdb_taxonomy",
        "genome_size",
        "contig_count",
        "longest_contig",
        "n50_contigs",
        "checkm2_completeness",
        "checkm2_contamination",
        "ncbi_assembly_level",
        "ncbi_genome_category",
        "ncbi_translation_table",
        "protein_count",
    )
    with tarfile.open(gtdb_tar, "r:*") as tf:
        for m in sorted(tf.getmembers(), key=lambda x: x.name):
            if not m.isfile() or not m.name.endswith(".tsv.gz"):
                continue
            dom = "Archaea" if m.name.startswith("ar53") else "Bacteria"
            fh = tf.extractfile(m)
            assert fh is not None
            with gzip.open(fh, "rt", newline="") as gz:
                rd = csv.reader(gz, delimiter="\t")
                head = next(rd)
                idx = {c: i for i, c in enumerate(head)}
                c_acc = idx["accession"]
                c_rep = idx["gtdb_representative"]
                c_grep = idx["gtdb_genome_representative"]
                c_gs = idx["genome_size"]
                for row in rd:
                    total_rows += 1
                    by_dom_total[dom] += 1
                    gacc = row[c_acc]
                    stripped = gacc[3:] if gacc.startswith(("GB_", "RS_")) else gacc
                    all_acc_to_rep[stripped] = (gacc, row[c_grep])
                    mt = GC_RE.search(stripped)
                    if mt:
                        stem_to_gtdb[mt.group(2)].append(gacc)
                    if row[c_rep] == "t":
                        by_dom_reps[dom] += 1
                        gs = int(row[c_gs])
                        if gs < 2_000_000:
                            reps_under_2mbp += 1
                        if float(row[idx["checkm2_completeness"]]) < 40.0:
                            reps_comp_lt_40 += 1
                        if float(row[idx["checkm2_contamination"]]) > 10.0:
                            reps_cont_gt_10 += 1
                        reps[gacc] = {c: row[idx[c]] for c in keep_cols}
    counts = {
        "total_genomes": total_rows,
        "total_by_domain": dict(by_dom_total),
        "species_reps": len(reps),
        "species_reps_by_domain": dict(by_dom_reps),
        "species_reps_under_2mbp": reps_under_2mbp,
        "species_reps_checkm2_comp_lt_40": reps_comp_lt_40,
        "species_reps_checkm2_cont_gt_10": reps_cont_gt_10,
    }
    return reps, all_acc_to_rep, {"counts": counts, "stem_to_gtdb": stem_to_gtdb}


def load_our_stages(root: Path, p16_sel: Path, c16_stages: Path, pl2_dir: Path, v1_dir: Path) -> dict:
    """Load accessions and stems from our built stages (root) and selection/input directories."""
    p16_selected: dict[str, dict] = {}
    with open(p16_sel, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            p16_selected[r["accession"]] = r

    # Built P1..P6 objects in l0_manifest
    p16_built: dict[str, str] = {}  # accession -> stage
    for st in ("P1", "P2", "P3", "P4", "P5", "P6"):
        mdir = root / "l0" / "l0_manifest" / f"stage={st}"
        for pf in sorted(mdir.glob("*.parquet")):
            t = pq.read_table(pf, columns=["object_path"])
            for op in t.column("object_path").to_pylist():
                mt = GC_RE.search(Path(op).name)
                if mt:
                    p16_built[mt.group(1)] = st

    c16: dict[str, str] = {}
    with open(c16_stages, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c16[r["accession"]] = r["stage"]

    # Check if V1 or PL2 has any matching identifiers
    v1_ids: set[str] = set()
    if (v1_dir / "records.tsv").is_file():
        with open(v1_dir / "records.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                v1_ids.add(r["genome"])
    pl2_ids: set[str] = set()
    if (pl2_dir / "records.tsv").is_file():
        with open(pl2_dir / "records.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pl2_ids.add(r["accession"])

    return {
        "p16_selected": p16_selected,
        "p16_built": p16_built,
        "c16": c16,
        "v1_count": len(v1_ids),
        "pl2_count": len(pl2_ids),
        "v1_ids": v1_ids,
        "pl2_ids": pl2_ids,
    }


def propose_stages(non_gtdb_rows: list[dict], group_key: str, prefix: str = "G") -> dict:
    """Propose ~10 Gbp stages in taxonomy order keeping each `group_key` ('genus_tax' or 'family_tax') whole."""
    rows = sorted(non_gtdb_rows, key=lambda r: (r["taxonomy"], r["ID"]))
    total_bp = sum(r["genome_size"] for r in rows)
    n = max(1, round(total_bp / GBP))

    # Check contiguity of group_key along sorted(taxonomy, ID)
    group_stage: dict[str, int] = {}
    before = 0
    members: dict[int, list[dict]] = defaultdict(list)
    for r in rows:
        gk = r[group_key]
        if gk not in group_stage:
            group_stage[gk] = 1 + (n * before) // total_bp
        members[group_stage[gk]].append(r)
        before += r["genome_size"]

    seq = [group_stage[r[group_key]] for r in rows]
    contiguous = seq == sorted(seq)

    stage_rows = []
    for k in sorted(members):
        ms = members[k]
        bp = sum(r["genome_size"] for r in ms)
        contigs = sum(r["contig_number"] for r in ms)
        phyla = list(dict.fromkeys(f"{r['domain']}:{r['phylum']}" for r in ms))
        genera = len({r["genus_tax"] for r in ms})
        families = len({r["family_tax"] for r in ms})
        u2m = sum(1 for r in ms if r["genome_size"] < 2_000_000)
        u100k = sum(1 for r in ms if r["genome_size"] < 100_000)
        hq = sum(1 for r in ms if r["checkm2_completeness"] > 90.0 and r["checkm2_contamination"] < 5.0)
        by_src = Counter(r["resource"] for r in ms)
        stage_rows.append({
            "stage": f"{prefix}{k}",
            "index": k,
            "genomes": len(ms),
            "bp": bp,
            "gbp": round(bp / 1e9, 4),
            "contigs": contigs,
            "domains": sorted({r["domain"] for r in ms}),
            "phyla_count": len(phyla),
            "families_count": families,
            "genera_count": genera,
            "under_2mbp_genomes": u2m,
            "under_100kbp_genomes": u100k,
            "mimag_hq_genomes": hq,
            "first_id": ms[0]["ID"],
            "last_id": ms[-1]["ID"],
            "first_taxon": ms[0]["taxonomy"],
            "last_taxon": ms[-1]["taxonomy"],
            "by_resource": {res: by_src[res] for res in RESOURCE_ORDER if by_src[res]},
        })
    bps = [s["bp"] for s in stage_rows]
    return {
        "group_key": group_key,
        "n_target": n,
        "n_actual": len(stage_rows),
        "total_genomes": len(rows),
        "total_bp": total_bp,
        "contiguous_in_tax_order": contiguous,
        "min_stage_bp": min(bps),
        "max_stage_bp": max(bps),
        "mean_stage_bp": round(total_bp / len(stage_rows)),
        "stages": stage_rows,
    }


def analyze(dl: Path, gtdb_tar: Path, root: Path, p16_sel: Path, c16_stages: Path,
            pl2_dir: Path, v1_dir: Path, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    dict_map, dict_summary = load_dictionaries(dl)
    gtdb_reps, gtdb_all_acc, gtdb_extra = load_gtdb(gtdb_tar)
    stem_to_gtdb = gtdb_extra["stem_to_gtdb"]
    ours = load_our_stages(root, p16_sel, c16_stages, pl2_dir, v1_dir)

    tax_map: dict[str, str] = {}
    with gzip.open(dl / "globdb_r226_taxonomy.tsv.gz", "rt", newline="") as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if r:
                tax_map[r[0]] = r[1]

    p16_sel_stems = {GC_RE.search(a).group(2): a for a in ours["p16_selected"] if GC_RE.search(a)}
    p16_built_stems = {GC_RE.search(a).group(2): a for a in ours["p16_built"] if GC_RE.search(a)}
    c16_stems = {GC_RE.search(a).group(2): a for a in ours["c16"] if GC_RE.search(a)}
    p16_sel_species = {r["gtdb_taxonomy"] for r in ours["p16_selected"].values()}

    total_rows = 0
    by_res: dict[str, dict] = defaultdict(lambda: {
        "genomes": 0, "bp": 0, "contigs": 0, "under_2mbp": 0, "under_2mbp_bp": 0,
        "under_1mbp": 0, "under_500kbp": 0, "under_100kbp": 0, "under_100kbp_bp": 0,
        "hq_strict": 0, "hq_ge90_le5": 0, "mq": 0, "lq": 0, "qscore_lt_50": 0,
        "cultured": 0, "archaea": 0, "bacteria": 0, "unclassified": 0,
        "sizes": [], "comp": [], "cont": [], "n50": [], "contigs_list": [],
    })
    by_sub: Counter = Counter()
    by_phylum: dict[tuple[str, str], dict] = defaultdict(lambda: {
        "gtdb_genomes": 0, "gtdb_bp": 0, "non_gtdb_genomes": 0, "non_gtdb_bp": 0,
        "non_gtdb_under_2mbp": 0, "non_gtdb_under_2mbp_bp": 0, "non_gtdb_contigs": 0,
    })

    # Index GTDB reps by full gtdb_acc (RS_GCF_xxx.1), stripped (GCF_xxx.1), versionless (GCF_xxx), and 9-digit stem
    rep_by_key: dict[str, str] = {}
    for gacc in gtdb_reps:
        stripped = gacc[3:] if gacc.startswith(("GB_", "RS_")) else gacc
        nover = stripped.split(".")[0]
        mt = GC_RE.search(stripped)
        rep_by_key[gacc] = gacc
        rep_by_key[stripped] = gacc
        rep_by_key[nover] = gacc
        if mt:
            rep_by_key[mt.group(2)] = gacc

    matched_gtdb_reps: set[str] = set()
    sample_gtdb_ids: list[str] = []
    globdb_gtdb_ids: set[str] = set()
    globdb_gtdb_acc_mismatch = 0
    globdb_gtdb_size_mismatch = 0
    gtdb_overlap = Counter()
    non_gtdb_gca_rows = []
    non_gtdb_v1_pl2_hits = Counter()

    non_gtdb_rows: list[dict] = []
    comp_bands = Counter()
    cont_bands = Counter()
    size_bands = Counter()

    with gzip.open(dl / "globdb_r226_tax_plus_stats.tsv.gz", "rt", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            total_rows += 1
            gid = r["ID"]
            res, sub = classify_id(gid)
            by_sub[sub] += 1
            gs = int(r["genome_size"])
            nc = int(r["contig_number"])
            lc = int(r["largest_contig"])
            n50 = int(r["contig_N50"])
            comp = float(r["checkm2_completeness"])
            cont = float(r["checkm2_contamination"])
            qscore = comp - 5.0 * cont
            cult = r["culture_collection"] not in ("0", "", "None", "NA")
            raw_dom = r["domain"]
            phy = r["phylum"]
            tax_str = tax_map[gid]
            if tax_str.startswith("d__Archaea;"):
                dom = "Archaea"
            elif tax_str.startswith("d__Bacteria;"):
                dom = "Bacteria"
            else:
                dom = "Unclassified"

            # Record stats per resource and aggregate groups
            groups = [res, "ALL_GLOBDB", "ALL_GTDB" if res == "GTDB" else "ALL_NON_GTDB"]
            if res in ("MRGM/HRGM2", "SPIRE"):
                groups.append(f"SUB:{sub}")
            hq_s = comp > 90.0 and cont < 5.0
            hq_i = comp >= 90.0 and cont <= 5.0
            mq_s = (not hq_s) and (comp >= 50.0 and cont < 10.0)
            for g in groups:
                st = by_res[g]
                st["genomes"] += 1
                st["bp"] += gs
                st["contigs"] += nc
                if gs < 2_000_000:
                    st["under_2mbp"] += 1
                    st["under_2mbp_bp"] += gs
                if gs < 1_000_000:
                    st["under_1mbp"] += 1
                if gs < 500_000:
                    st["under_500kbp"] += 1
                if gs < 100_000:
                    st["under_100kbp"] += 1
                    st["under_100kbp_bp"] += gs
                st["hq_strict"] += hq_s
                st["hq_ge90_le5"] += hq_i
                st["mq"] += mq_s
                st["lq"] += (not hq_s and not mq_s)
                st["qscore_lt_50"] += qscore < 50.0
                st["cultured"] += cult
                if dom == "Archaea":
                    st["archaea"] += 1
                elif dom == "Bacteria":
                    st["bacteria"] += 1
                else:
                    st["unclassified"] += 1
                if g in ("ALL_GTDB", "ALL_NON_GTDB") or g in RESOURCE_ORDER or g.startswith("SUB:"):
                    st["sizes"].append(gs)
                    st["comp"].append(comp)
                    st["cont"].append(cont)
                    st["n50"].append(n50)
                    st["contigs_list"].append(nc)

            pst = by_phylum[(dom, phy if dom != "Unclassified" else f"Unclassified (raw domain={raw_dom})")]
            if res == "GTDB":
                pst["gtdb_genomes"] += 1
                pst["gtdb_bp"] += gs
            else:
                pst["non_gtdb_genomes"] += 1
                pst["non_gtdb_bp"] += gs
                pst["non_gtdb_contigs"] += nc
                if gs < 2_000_000:
                    pst["non_gtdb_under_2mbp"] += 1
                    pst["non_gtdb_under_2mbp_bp"] += gs

            if res == "GTDB":
                globdb_gtdb_ids.add(gid)
                if len(sample_gtdb_ids) < 5:
                    sample_gtdb_ids.append(gid)
                stripped = gid[3:] if gid.startswith(("GB_", "RS_")) else gid
                mt = GC_RE.search(stripped)
                stem = mt.group(2) if mt else ""
                matched_gacc = rep_by_key.get(gid) or rep_by_key.get(stripped) or rep_by_key.get(stem)
                if matched_gacc is None:
                    globdb_gtdb_acc_mismatch += 1
                else:
                    matched_gtdb_reps.add(matched_gacc)
                    if int(gtdb_reps[matched_gacc]["genome_size"]) != gs:
                        globdb_gtdb_size_mismatch += 1
                    # Use GTDB's full versioned accession for exact P1..P6 / C1..C6 comparison too
                    gtdb_ver_acc = matched_gacc[3:] if matched_gacc.startswith(("GB_", "RS_")) else matched_gacc
                    in_p16_b = gtdb_ver_acc in ours["p16_built"]
                    in_p16_s = gtdb_ver_acc in ours["p16_selected"]
                    in_c16 = gtdb_ver_acc in ours["c16"]
                    in_p16_b_stem = stem in p16_built_stems
                    in_p16_s_stem = stem in p16_sel_stems
                    in_c16_stem = stem in c16_stems
                    gtdb_overlap["in_p16_built_exact"] += in_p16_b
                    gtdb_overlap["in_p16_selected_exact"] += in_p16_s
                    gtdb_overlap["in_c16_exact"] += in_c16
                    gtdb_overlap["in_p16_or_c16_exact"] += (in_p16_b or in_c16)
                    gtdb_overlap["in_p16_built_stem"] += in_p16_b_stem
                    gtdb_overlap["in_p16_selected_stem"] += in_p16_s_stem
                    gtdb_overlap["in_c16_stem"] += in_c16_stem
                    gtdb_overlap["in_p16_or_c16_stem"] += (in_p16_b_stem or in_c16_stem)
                    gtdb_overlap["species_has_p16_complete_genome"] += (
                        gtdb_reps[matched_gacc]["gtdb_taxonomy"] in p16_sel_species
                    )
            else:
                # Non-GTDB distributions
                if comp >= 95:
                    comp_bands[">=95"] += 1
                elif comp >= 90:
                    comp_bands["90-<95"] += 1
                elif comp >= 70:
                    comp_bands["70-<90"] += 1
                elif comp >= 50:
                    comp_bands["50-<70"] += 1
                elif comp >= 40:
                    comp_bands["40-<50"] += 1
                else:
                    comp_bands["<40"] += 1

                if cont == 0:
                    cont_bands["0"] += 1
                elif cont <= 1:
                    cont_bands["(0,1]"] += 1
                elif cont <= 2:
                    cont_bands["(1,2]"] += 1
                elif cont <= 5:
                    cont_bands["(2,5]"] += 1
                elif cont <= 10:
                    cont_bands["(5,10]"] += 1
                else:
                    cont_bands[">10"] += 1

                if gs < 100_000:
                    size_bands["<100kb"] += 1
                elif gs < 500_000:
                    size_bands["100kb-<500kb"] += 1
                elif gs < 1_000_000:
                    size_bands["500kb-<1Mbp"] += 1
                elif gs < 2_000_000:
                    size_bands["1Mbp-<2Mbp"] += 1
                elif gs < 3_000_000:
                    size_bands["2Mbp-<3Mbp"] += 1
                elif gs < 5_000_000:
                    size_bands["3Mbp-<5Mbp"] += 1
                else:
                    size_bands[">=5Mbp"] += 1

                orig = dict_map.get(gid, gid)
                if gid in ours["v1_ids"] or orig in ours["v1_ids"]:
                    non_gtdb_v1_pl2_hits["v1"] += 1
                if gid in ours["pl2_ids"] or orig in ours["pl2_ids"]:
                    non_gtdb_v1_pl2_hits["pl2"] += 1
                mt = GC_RE.search(orig) or GC_RE.search(gid)
                if mt:
                    acc = mt.group(1)
                    stem = mt.group(2)
                    non_gtdb_gca_rows.append({
                        "globdb_id": gid,
                        "resource": res,
                        "original_id": orig,
                        "gca_acc": acc,
                        "stem": stem,
                        "in_p16_built_exact": acc in ours["p16_built"],
                        "in_p16_built_stem": stem in p16_built_stems,
                        "in_p16_selected_stem": stem in p16_sel_stems,
                        "in_c16_exact": acc in ours["c16"],
                        "in_c16_stem": stem in c16_stems,
                        "in_gtdb_r226_exact": acc in gtdb_all_acc,
                        "in_gtdb_r226_stem": stem in stem_to_gtdb,
                        "gtdb_r226_matches": stem_to_gtdb.get(stem, []),
                        "taxonomy": tax_str,
                        "genome_size": gs,
                    })

                parts = tax_str.split(";")
                non_gtdb_rows.append({
                    "ID": gid,
                    "resource": res,
                    "domain": dom,
                    "phylum": phy,
                    "taxonomy": tax_str,
                    "phylum_tax": ";".join(parts[:2]),
                    "family_tax": ";".join(parts[:5]),
                    "genus_tax": ";".join(parts[:6]),
                    "genome_size": gs,
                    "contig_number": nc,
                    "largest_contig": lc,
                    "contig_N50": n50,
                    "checkm2_completeness": comp,
                    "checkm2_contamination": cont,
                })

    # Find the exact GTDB R226 species reps absent from GlobDB r226
    missing_gtdb_ids = sorted(set(gtdb_reps) - matched_gtdb_reps)
    missing_gtdb = []
    for gacc in missing_gtdb_ids:
        info = dict(gtdb_reps[gacc])
        mt = GC_RE.search(gacc)
        stem = mt.group(2) if mt else ""
        info["in_p16_built"] = (gacc[3:] in ours["p16_built"])
        info["in_c16"] = (gacc[3:] in ours["c16"])
        info["stem_in_non_gtdb_gca"] = any(x["stem"] == stem for x in non_gtdb_gca_rows)
        missing_gtdb.append(info)

    # Summarise per-resource table cleanly
    res_out = {}
    for k in [*RESOURCE_ORDER, "SUB:SPIREOTU", "SUB:SPECIV4", "SUB:HRGM2", "SUB:MRGM",
              "ALL_GTDB", "ALL_NON_GTDB", "ALL_GLOBDB"]:
        if k not in by_res:
            continue
        st = by_res[k]
        res_out[k] = {
            "genomes": st["genomes"],
            "bp": st["bp"],
            "gbp": round(st["bp"] / 1e9, 4),
            "contigs": st["contigs"],
            "archaea": st["archaea"],
            "bacteria": st["bacteria"],
            "unclassified": st["unclassified"],
            "under_2mbp_genomes": st["under_2mbp"],
            "under_2mbp_bp": st["under_2mbp_bp"],
            "under_1mbp_genomes": st["under_1mbp"],
            "under_500kbp_genomes": st["under_500kbp"],
            "under_100kbp_genomes": st["under_100kbp"],
            "under_100kbp_bp": st["under_100kbp_bp"],
            "mimag_hq_strict": st["hq_strict"],
            "mimag_hq_ge90_le5": st["hq_ge90_le5"],
            "mimag_mq": st["mq"],
            "mimag_lq": st["lq"],
            "qscore_lt_50": st["qscore_lt_50"],
            "cultured": st["cultured"],
            "size_quantiles": quantiles(st["sizes"]),
            "completeness_quantiles": quantiles(st["comp"]),
            "contamination_quantiles": quantiles(st["cont"]),
            "contig_n50_quantiles": quantiles(st["n50"]),
            "contig_count_quantiles": quantiles(st["contigs_list"]),
        }

    # Phylum table
    phylum_rows = []
    for (dom, phy), pst in sorted(by_phylum.items(), key=lambda kv: (-kv[1]["non_gtdb_genomes"], kv[0][0], kv[0][1])):
        tot_g = pst["gtdb_genomes"] + pst["non_gtdb_genomes"]
        phylum_rows.append({
            "domain": dom,
            "phylum": phy,
            "non_gtdb_genomes": pst["non_gtdb_genomes"],
            "non_gtdb_bp": pst["non_gtdb_bp"],
            "non_gtdb_gbp": round(pst["non_gtdb_bp"] / 1e9, 4),
            "non_gtdb_contigs": pst["non_gtdb_contigs"],
            "non_gtdb_under_2mbp": pst["non_gtdb_under_2mbp"],
            "non_gtdb_under_2mbp_bp": pst["non_gtdb_under_2mbp_bp"],
            "gtdb_genomes": pst["gtdb_genomes"],
            "gtdb_bp": pst["gtdb_bp"],
            "total_globdb_genomes": tot_g,
            "non_gtdb_share": round(pst["non_gtdb_genomes"] / tot_g, 4) if tot_g else 0.0,
        })

    # Write full phylum TSV
    with open(out / "phyla.tsv", "w", newline="") as fh:
        cols = list(phylum_rows[0].keys())
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(phylum_rows)

    #Largest genus and family sizes in non-GTDB
    genus_bp = Counter()
    genus_n = Counter()
    fam_bp = Counter()
    fam_n = Counter()
    for r in non_gtdb_rows:
        genus_bp[r["genus_tax"]] += r["genome_size"]
        genus_n[r["genus_tax"]] += 1
        fam_bp[r["family_tax"]] += r["genome_size"]
        fam_n[r["family_tax"]] += 1

    prop_genus = propose_stages(non_gtdb_rows, "genus_tax", prefix="G")
    prop_family = propose_stages(non_gtdb_rows, "family_tax", prefix="G")

    # Write proposed genus-coherent stage table TSV and full per-genome assignment TSV
    with open(out / "stage_proposal.tsv", "w", newline="") as fh:
        cols = [
            "stage", "genomes", "bp", "contigs", "phyla_count", "families_count", "genera_count",
            "under_2mbp_genomes", "under_100kbp_genomes", "mimag_hq_genomes",
            "first_id", "last_id", "first_taxon", "last_taxon",
        ]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(prop_genus["stages"])

    # Also save genome -> stage assignment for the genus-coherent proposal
    sorted_non_gtdb = sorted(non_gtdb_rows, key=lambda r: (r["taxonomy"], r["ID"]))
    total_non_gtdb_bp = sum(r["genome_size"] for r in sorted_non_gtdb)
    n_stages = prop_genus["n_target"]
    g_stage_map: dict[str, int] = {}
    acc_bp = 0
    with open(out / "non_gtdb_stages.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ID", "stage", "position_in_stage", "resource", "genome_size", "contig_number",
                    "checkm2_completeness", "checkm2_contamination", "taxonomy"])
        pos_by_stage: Counter = Counter()
        for r in sorted_non_gtdb:
            gk = r["genus_tax"]
            if gk not in g_stage_map:
                g_stage_map[gk] = 1 + (n_stages * acc_bp) // total_non_gtdb_bp
            sk = f"G{g_stage_map[gk]}"
            pos = pos_by_stage[sk]
            pos_by_stage[sk] += 1
            w.writerow([r["ID"], sk, pos, r["resource"], r["genome_size"], r["contig_number"],
                        r["checkm2_completeness"], r["checkm2_contamination"], r["taxonomy"]])
            acc_bp += r["genome_size"]

    non_gtdb_gca_summary = {
        "total_non_gtdb_with_gca_gcf": len(non_gtdb_gca_rows),
        "by_resource": dict(Counter(x["resource"] for x in non_gtdb_gca_rows)),
        "in_p16_built_exact": sum(1 for x in non_gtdb_gca_rows if x["in_p16_built_exact"]),
        "in_p16_built_stem": sum(1 for x in non_gtdb_gca_rows if x["in_p16_built_stem"]),
        "in_p16_selected_stem": sum(1 for x in non_gtdb_gca_rows if x["in_p16_selected_stem"]),
        "in_c16_exact": sum(1 for x in non_gtdb_gca_rows if x["in_c16_exact"]),
        "in_c16_stem": sum(1 for x in non_gtdb_gca_rows if x["in_c16_stem"]),
        "in_gtdb_r226_exact": sum(1 for x in non_gtdb_gca_rows if x["in_gtdb_r226_exact"]),
        "in_gtdb_r226_stem": sum(1 for x in non_gtdb_gca_rows if x["in_gtdb_r226_stem"]),
        "examples_in_gtdb_or_ours": [
            x for x in non_gtdb_gca_rows
            if x["in_gtdb_r226_stem"] or x["in_p16_selected_stem"] or x["in_c16_stem"]
        ][:50],
        "first_10_gca": non_gtdb_gca_rows[:10],
    }

    qscores_non_gtdb = [
        r["checkm2_completeness"] - 5.0 * r["checkm2_contamination"] for r in non_gtdb_rows
    ]
    largest_contigs_non_gtdb = [r["largest_contig"] for r in non_gtdb_rows]

    result = {
        "total_globdb_rows": total_rows,
        "taxonomy_rows": len(tax_map),
        "gtdb_r226_metadata_counts": gtdb_extra["counts"],
        "globdb_gtdb_count": len(globdb_gtdb_ids),
        "sample_globdb_gtdb_ids": sample_gtdb_ids,
        "globdb_non_gtdb_count": len(non_gtdb_rows),
        "globdb_gtdb_not_in_gtdb_reps": globdb_gtdb_acc_mismatch,
        "globdb_gtdb_size_mismatch_vs_gtdb": globdb_gtdb_size_mismatch,
        "missing_7_gtdb_reps": missing_gtdb,
        "dictionaries": dict_summary,
        "by_sub_source": dict(by_sub),
        "resources": res_out,
        "non_gtdb_distributions": {
            "completeness_bands": dict(comp_bands),
            "contamination_bands": dict(cont_bands),
            "size_bands": dict(size_bands),
            "qscore_quantiles": quantiles(qscores_non_gtdb),
            "largest_contig_quantiles": quantiles(largest_contigs_non_gtdb),
        },
        "phyla_summary": {
            "total_phyla_in_globdb": len(phylum_rows),
            "archaea_phyla": sum(1 for p in phylum_rows if p["domain"] == "Archaea"),
            "bacteria_phyla": sum(1 for p in phylum_rows if p["domain"] == "Bacteria"),
            "unclassified_phyla": sum(1 for p in phylum_rows if p["domain"] == "Unclassified"),
            "novel_phyla_zero_gtdb": [p for p in phylum_rows if p["gtdb_genomes"] == 0],
            "top_25_phyla_by_non_gtdb_genomes": phylum_rows[:25],
            "archaea_all_phyla": [p for p in phylum_rows if p["domain"] == "Archaea"],
        },
        "accession_overlap": {
            "our_stage_counts": {
                "p16_selected": len(ours["p16_selected"]),
                "p16_built": len(ours["p16_built"]),
                "c16_fetched": len(ours["c16"]),
                "v1_records": ours["v1_count"],
                "pl2_records": ours["pl2_count"],
            },
            "globdb_gtdb_in_our_stages": dict(gtdb_overlap),
            "non_gtdb_gca_overlap": non_gtdb_gca_summary,
            "non_gtdb_id_in_v1_or_pl2": dict(non_gtdb_v1_pl2_hits),
        },
        "taxon_group_maxima_non_gtdb": {
            "distinct_genera": len(genus_bp),
            "distinct_families": len(fam_bp),
            "top_5_genera_by_bp": [
                {"genus": k, "genomes": genus_n[k], "bp": v}
                for k, v in genus_bp.most_common(5)
            ],
            "top_5_families_by_bp": [
                {"family": k, "genomes": fam_n[k], "bp": v}
                for k, v in fam_bp.most_common(5)
            ],
        },
        "stage_proposal_genus_coherent": prop_genus,
        "stage_proposal_family_coherent": {
            k: v for k, v in prop_family.items() if k != "stages"
        },
    }
    (out / "table_summary.json").write_text(json.dumps(result, indent=1) + "\n")
    return result


def probe_tar(url: str, nbytes: int) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "curl", "Range": f"bytes=0-{nbytes - 1}"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        raw = resp.read()
    members = []
    bio = pyio.BytesIO(raw)
    try:
        with tarfile.open(fileobj=bio, mode="r|gz") as tf:
            for m in tf:
                rec = {"name": m.name, "size": m.size, "isfile": m.isfile(), "isdir": m.isdir()}
                if m.isfile() and len(members) < 5:
                    fh = tf.extractfile(m)
                    if fh is not None:
                        data = fh.read()
                        with gzip.open(pyio.BytesIO(data), "rt") as gzf:
                            head_lines = [gzf.readline().rstrip("\r\n")[:120] for _ in range(4)]
                        rec["first_lines"] = head_lines
                members.append(rec)
                if len(members) >= 25:
                    break
    except (EOFError, tarfile.ReadError, gzip.BadGzipFile):
        pass
    return {"range_bytes": len(raw), "members_seen": len(members), "members": members}


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_an = sub.add_parser("analyze")
    p_an.add_argument("--dl", required=True, type=Path)
    p_an.add_argument("--gtdb-tar", required=True, type=Path)
    p_an.add_argument("--root", required=True, type=Path)
    p_an.add_argument("--p16-sel", required=True, type=Path)
    p_an.add_argument("--c16-stages", required=True, type=Path)
    p_an.add_argument("--pl2-inputs", required=True, type=Path)
    p_an.add_argument("--v1-inputs", required=True, type=Path)
    p_an.add_argument("--out", required=True, type=Path)

    p_tr = sub.add_parser("probe-tar")
    p_tr.add_argument("--url", required=True)
    p_tr.add_argument("--bytes", type=int, default=1_048_576)

    a = ap.parse_args(argv)
    if a.cmd == "analyze":
        res = analyze(a.dl, a.gtdb_tar, a.root, a.p16_sel, a.c16_stages, a.pl2_inputs, a.v1_inputs, a.out)
        print(json.dumps({
            "total_globdb_rows": res["total_globdb_rows"],
            "globdb_gtdb_count": res["globdb_gtdb_count"],
            "globdb_non_gtdb_count": res["globdb_non_gtdb_count"],
            "missing_7_count": len(res["missing_7_gtdb_reps"]),
            "stage_proposal_n": res["stage_proposal_genus_coherent"]["n_actual"],
            "total_non_gtdb_bp": res["stage_proposal_genus_coherent"]["total_bp"],
        }, indent=1))
    elif a.cmd == "probe-tar":
        print(json.dumps(probe_tar(a.url, a.bytes), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

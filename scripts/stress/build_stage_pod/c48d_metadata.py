#!/usr/bin/env python3
"""Card C-48d stage 1: metadata-level deduplication of GEMs, IUIG, and IMG-UMAG
against GTDB R226, GlobDB r226, and IMG/M + GOLD tables.

Streams directly from read-only archives on `/stor/work` and `/stor/scratch`
without unpacking large archives onto `/tmp` (AGENTS.md §4.1, §7, §8).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import statistics
import sys
import tarfile
import zipfile
from collections import Counter
from pathlib import Path

GCA_RE = re.compile(r"(GC[AF]_\d{9})(?:\.\d+)?")


def _stem(acc: str) -> str:
    m = GCA_RE.search(acc or "")
    return m.group(1) if m else ""


def _bare_acc(acc: str) -> str:
    s = (acc or "").strip()
    if s.startswith(("GB_", "RS_")):
        return s[3:]
    return s


def load_reference_sets(
    gtdb_tar: Path | None = None,
    globdb_dl_dir: Path | None = None,
    img_m_tar: Path | None = None,
) -> dict:
    """Load GTDB R226, GlobDB r226, and IMG/M reference lookup sets in bounded memory."""
    refs: dict = {
        "gtdb_accessions": set(),
        "gtdb_stems": set(),
        "gtdb_rep_stems": set(),
        "gtdb_ncbi_taxids": set(),
        "gtdb_species_names": set(),
        "globdb_all_ids": set(),
        "globdb_gem_ids": set(),
        "globdb_gem_dict_all": set(),
        "globdb_gca_stems": set(),
        "img_oid_to_gca": {},
        "img_oid_to_eco": {},
        "counts": {},
    }

    if gtdb_tar is not None and gtdb_tar.is_file():
        n_gtdb = 0
        n_gtdb_reps = 0
        with tarfile.open(gtdb_tar, "r") as tf:
            for member in tf.getmembers():
                if not member.name.endswith(".tsv.gz"):
                    continue
                fobj = tf.extractfile(member)
                if fobj is None:
                    continue
                with gzip.open(fobj, "rt", encoding="utf-8", errors="replace") as gz:
                    reader = csv.DictReader(gz, delimiter="\t")
                    for row in reader:
                        n_gtdb += 1
                        acc = _bare_acc(row.get("accession", ""))
                        if acc:
                            refs["gtdb_accessions"].add(acc)
                            st = _stem(acc)
                            if st:
                                refs["gtdb_stems"].add(st)
                                if row.get("gtdb_representative") == "t":
                                    n_gtdb_reps += 1
                                    refs["gtdb_rep_stems"].add(st)
                        taxid = (row.get("ncbi_taxid") or "").strip()
                        if taxid and taxid != "none":
                            refs["gtdb_ncbi_taxids"].add(taxid)
                        gtax = row.get("gtdb_taxonomy") or ""
                        for part in gtax.split(";"):
                            if part.startswith("s__") and len(part) > 3:
                                refs["gtdb_species_names"].add(part[3:].strip())
                        org = (row.get("ncbi_organism_name") or "").strip()
                        if org and org != "none":
                            toks = org.split()
                            if len(toks) >= 2:
                                refs["gtdb_species_names"].add(f"{toks[0]} {toks[1]}")
        refs["counts"]["gtdb_rows"] = n_gtdb
        refs["counts"]["gtdb_reps"] = n_gtdb_reps

    if globdb_dl_dir is not None and globdb_dl_dir.is_dir():
        stats_gz = globdb_dl_dir / "globdb_r226_tax_plus_stats.tsv.gz"
        if stats_gz.is_file():
            with gzip.open(stats_gz, "rt", encoding="utf-8", errors="replace") as gz:
                reader = csv.DictReader(gz, delimiter="\t")
                for row in reader:
                    gid = (row.get("ID") or "").strip()
                    if gid:
                        refs["globdb_all_ids"].add(gid)
                        st = _stem(gid)
                        if st:
                            refs["globdb_gca_stems"].add(st)
        dict_tgz = globdb_dl_dir / "globdb_r226_dictionaries.tar.gz"
        if dict_tgz.is_file():
            with tarfile.open(dict_tgz, "r:gz") as tf:
                for member in tf.getmembers():
                    if not member.isfile() or not member.name.endswith(".tsv"):
                        continue
                    fobj = tf.extractfile(member)
                    if fobj is None:
                        continue
                    is_gem = "GEM" in member.name.upper()
                    text = io.TextIOWrapper(fobj, encoding="utf-8", errors="replace")
                    for line in text:
                        parts = line.rstrip("\r\n").split("\t")
                        if len(parts) < 2:
                            continue
                        glob_id = parts[0].strip()
                        if glob_id == "GlobDB_ID":
                            continue
                        orig_ids = [p.strip() for p in parts[1:] if p.strip()]
                        for orig_id in orig_ids:
                            st = _stem(orig_id)
                            if st and glob_id in refs["globdb_all_ids"]:
                                refs["globdb_gca_stems"].add(st)
                            if is_gem:
                                refs["globdb_gem_dict_all"].add(orig_id)
                                if glob_id in refs["globdb_all_ids"]:
                                    refs["globdb_gem_ids"].add(orig_id)
                                    refs["globdb_gem_ids"].add(glob_id)
        refs["counts"]["globdb_all_ids"] = len(refs["globdb_all_ids"])
        refs["counts"]["globdb_gem_kept_orig_ids"] = len(refs["globdb_gem_ids"])
        refs["counts"]["globdb_gem_dict_all"] = len(refs["globdb_gem_dict_all"])

    if img_m_tar is not None and img_m_tar.is_file():
        n_img_rows = 0
        with tarfile.open(img_m_tar, "r") as tf:
            for member in tf.getmembers():
                if not member.name.endswith(".tsv.gz"):
                    continue
                fobj = tf.extractfile(member)
                if fobj is None:
                    continue
                with gzip.open(fobj, "rt", encoding="utf-8", errors="replace") as gz:
                    reader = csv.DictReader(gz, delimiter="\t")
                    fnames = reader.fieldnames or []
                    oid_col = next((c for c in fnames if "IMG Genome ID" in c or c == "taxon_oid"), None)
                    acc_col = next((c for c in fnames if "NCBI Assembly Accession" in c or "GenBank Assembly" in c), None)
                    eco_cols = [c for c in ("Ecosystem", "Ecosystem Category", "Ecosystem Type") if c in fnames]
                    if not oid_col:
                        continue
                    for row in reader:
                        n_img_rows += 1
                        oid = (row.get(oid_col) or "").strip()
                        if not oid:
                            continue
                        if acc_col:
                            val = (row.get(acc_col) or "").strip()
                            m = GCA_RE.search(val)
                            if m:
                                refs["img_oid_to_gca"][oid] = m.group(0)
                        if eco_cols:
                            eco = " / ".join((row.get(c) or "").strip() for c in eco_cols if (row.get(c) or "").strip())
                            if eco:
                                refs["img_oid_to_eco"][oid] = eco
        refs["counts"]["img_m_rows_indexed"] = n_img_rows
        refs["counts"]["img_oid_to_gca"] = len(refs["img_oid_to_gca"])
        refs["counts"]["img_oid_to_eco"] = len(refs["img_oid_to_eco"])

    return refs


def analyze_gems_metadata(gems_meta_tar: Path, refs: dict) -> dict:
    """Analyze GEMs 52,515 MAGs, 45,599 OTUs, and 24,313 VirSorter prophages against GTDB + GlobDB."""
    otu_rep_map: dict[str, str] = {}
    otu_tax_map: dict[str, dict] = {}
    mags: list[dict] = []
    virsorter_rows: list[dict] = []

    with tarfile.open(gems_meta_tar, "r") as tf:
        members = {m.name: m for m in tf.getmembers() if m.isfile()}

        if "gems/otus/otu_representative_genomes.tsv" in members:
            f = tf.extractfile(members["gems/otus/otu_representative_genomes.tsv"])
            assert f is not None
            for row in csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t"):
                otu_rep_map[row["otu_id"]] = row["genome_id"]

        if "gems/otus/otu_taxonomy.tsv" in members:
            f = tf.extractfile(members["gems/otus/otu_taxonomy.tsv"])
            assert f is not None
            for row in csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t"):
                otu_tax_map[row["otu_id"]] = {
                    "taxonomy": row.get("taxonomy", ""),
                    "domain": row.get("t_domain", ""),
                    "phylum": row.get("t_phylum", ""),
                }

        if "gems/genomes/genome_metadata.tsv" in members:
            f = tf.extractfile(members["gems/genomes/genome_metadata.tsv"])
            assert f is not None
            for row in csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t"):
                mags.append(row)

        if "gems/prophages/virsorter.tsv" in members:
            f = tf.extractfile(members["gems/prophages/virsorter.tsv"])
            assert f is not None
            for row in csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t"):
                virsorter_rows.append(row)

    globdb_gem_ids: set[str] = set(refs.get("globdb_gem_ids", set()))
    globdb_all_ids: set[str] = set(refs.get("globdb_all_ids", set()))
    gtdb_accs: set[str] = set(refs.get("gtdb_accessions", set()))
    gtdb_stems: set[str] = set(refs.get("gtdb_stems", set()))
    globdb_gem_dict_all: set[str] = set(refs.get("globdb_gem_dict_all", set()))

    total_mags = len(mags)
    hq_mq_mags = 0
    lq_mags = 0
    high_contam_mags = 0
    already_in_globdb_or_gtdb = 0
    dropped_non_rep_same_otu = 0
    globdb_evaluated_dropped = 0
    by_phylum: Counter[str] = Counter()
    by_ecosystem: Counter[str] = Counter()
    all_metagenome_ids: set[str] = set()

    # Keep best HQ/MQ MAG per GEM OTU not already represented in GlobDB/GTDB
    best_by_otu: dict[str, dict] = {}
    covered_bins: set[str] = set()
    novel_candidate_bins: set[str] = set()

    for r in mags:
        gid = r["genome_id"].strip()
        mg_id = (r.get("metagenome_id") or "").strip()
        if mg_id:
            all_metagenome_ids.add(mg_id)
        otu_id = (r.get("otu_id") or "").strip()
        comp = float(r["completeness"])
        cont = float(r["contamination"])
        qscore = float(r.get("quality_score") or (comp - 5.0 * cont))
        bp = int(r["genome_length"])
        n50 = int(r["n50"])
        n_contigs = int(r["num_contigs"])
        eco_raw = (r.get("ecosystem") or "").strip()
        # Real NERSC GEM genome_metadata.tsv omits a "taxonomy" header column between otu_id and ecosystem
        # (20 headers vs 21 fields), shifting ecosystem -> ecosystem_category in DictReader.
        if eco_raw.startswith(("d__", "Unclassified")):
            eco = (r.get("ecosystem_category") or "").strip() or "Unclassified"
            eco_cat = (r.get("ecosystem_type") or "").strip()
        else:
            eco = eco_raw or "Unclassified"
            eco_cat = (r.get("ecosystem_category") or "").strip()
        eco_path = f"{eco} / {eco_cat}" if eco_cat else eco

        if comp < 50.0:
            lq_mags += 1
            continue
        if cont > 10.0:
            high_contam_mags += 1
            continue
        hq_mq_mags += 1

        # Check if this MAG or its OTU representative is already in GlobDB r226 or GTDB R226
        rep_id = otu_rep_map.get(otu_id, "")
        rep_stem = _stem(rep_id)
        in_ref = (
            gid in globdb_gem_ids
            or gid in globdb_all_ids
            or otu_id in globdb_gem_ids
            or rep_id in globdb_gem_ids
            or rep_id in globdb_all_ids
            or rep_id in gtdb_accs
            or (rep_stem and rep_stem in gtdb_stems)
        )
        if in_ref:
            already_in_globdb_or_gtdb += 1
            covered_bins.add(gid)
            continue

        # Also track whether GlobDB r226 evaluated this OTU rep in 05_GEMOTU_dictionary.tsv
        # (meaning GlobDB clustered it into GTDB/mOTUs/SPIRE/RBG at >=96% ANI or dropped it <40% CheckM2)
        in_globdb_dict = (gid in globdb_gem_dict_all) or (otu_id in globdb_gem_dict_all) or (rep_id in globdb_gem_dict_all)
        if in_globdb_dict:
            globdb_evaluated_dropped += 1

        tax_info = otu_tax_map.get(otu_id, {})
        rec = {
            "genome_id": gid,
            "metagenome_id": mg_id,
            "otu_id": otu_id,
            "otu_rep_id": rep_id,
            "bp": bp,
            "num_contigs": n_contigs,
            "n50": n50,
            "completeness": comp,
            "contamination": cont,
            "quality_score": qscore,
            "mimag_quality": (r.get("mimag_quality") or "").strip(),
            "taxonomy": tax_info.get("taxonomy", ""),
            "domain": tax_info.get("domain", ""),
            "phylum": tax_info.get("phylum", "") or "Unclassified",
            "ecosystem": eco_path,
            "in_globdb_dict_dropped": in_globdb_dict,
            "source": "gems",
        }
        prev = best_by_otu.get(otu_id)
        if prev is None or (rec["quality_score"], rec["n50"], -len(gid), gid) > (
            prev["quality_score"], prev["n50"], -len(prev["genome_id"]), prev["genome_id"]
        ):
            if prev is not None:
                dropped_non_rep_same_otu += 1
            best_by_otu[otu_id] = rec
        else:
            dropped_non_rep_same_otu += 1

    candidates = sorted(best_by_otu.values(), key=lambda x: (x["taxonomy"], x["genome_id"]))
    for c in candidates:
        novel_candidate_bins.add(c["genome_id"])
        by_phylum[c["phylum"]] += 1
        by_ecosystem[c["ecosystem"]] += 1

    # VirSorter prophage sidecar summary
    vp_bins: set[str] = set()
    vp_contigs: set[str] = set()
    vp_lens: list[int] = []
    vp_on_globdb = 0
    vp_on_novel = 0
    for vr in virsorter_rows:
        bid = vr["bin_id"].strip()
        cid = vr["contig_id"].strip()
        start = int(vr["contig_start"])
        stop = int(vr["contig_stop"])
        vp_bins.add(bid)
        vp_contigs.add(cid)
        vp_lens.append(max(0, stop - start))
        if bid in covered_bins:
            vp_on_globdb += 1
        elif bid in novel_candidate_bins:
            vp_on_novel += 1

    vp_summary = {
        "total_prophages": len(virsorter_rows),
        "unique_bins": len(vp_bins),
        "unique_contigs": len(vp_contigs),
        "on_globdb_bins": vp_on_globdb,
        "on_novel_candidate_bins": vp_on_novel,
        "median_len_bp": int(statistics.median(vp_lens)) if vp_lens else 0,
        "total_len_bp": sum(vp_lens),
    }

    return {
        "total_mags": total_mags,
        "total_otus_in_rep_table": len(otu_rep_map),
        "hq_mq_mags": hq_mq_mags,
        "lq_mags_below_50_comp": lq_mags,
        "high_contam_above_10": high_contam_mags,
        "already_in_globdb_or_gtdb": already_in_globdb_or_gtdb,
        "globdb_evaluated_and_collapsed": globdb_evaluated_dropped,
        "dropped_non_rep_same_otu": dropped_non_rep_same_otu,
        "novel_otu_candidates": len(candidates),
        "novel_otu_not_in_globdb_dict": sum(1 for c in candidates if not c["in_globdb_dict_dropped"]),
        "novel_otu_bases": sum(c["bp"] for c in candidates),
        "by_phylum": dict(by_phylum.most_common(30)),
        "by_ecosystem": dict(by_ecosystem.most_common(20)),
        "metagenome_ids_count": len(all_metagenome_ids),
        "metagenome_ids": sorted(all_metagenome_ids),
        "virsorter": vp_summary,
        "candidates": candidates,
    }


def analyze_iuig_umag_metadata(meta_zip: Path, refs: dict) -> dict:
    """Analyze IUIG (143,797 rows) and IMG-UMAG (18,912 rows) metadata against GTDB + GlobDB."""
    gtdb_taxids: set[str] = set(refs.get("gtdb_ncbi_taxids", set()))
    gtdb_species: set[str] = set(refs.get("gtdb_species_names", set()))
    gtdb_stems: set[str] = set(refs.get("gtdb_stems", set()))
    globdb_stems: set[str] = set(refs.get("globdb_gca_stems", set()))
    img_oid_to_gca: dict[str, str] = dict(refs.get("img_oid_to_gca", {}))
    img_oid_to_eco: dict[str, str] = dict(refs.get("img_oid_to_eco", {}))
    gems_mgs: set[str] = set(refs.get("gems_metagenome_ids", set()))

    iuig_by_domain: Counter[str] = Counter()
    iuig_by_status: Counter[str] = Counter()
    iuig_prok_phylum: Counter[str] = Counter()
    iuig_matched = 0
    iuig_matched_by_gca = 0
    iuig_matched_by_taxid = 0
    iuig_matched_by_species = 0
    iuig_prok_novel_oids: list[str] = []
    iuig_prok_novel_rows: list[dict] = []
    iuig_total = 0

    umag_total = 0
    umag_by_eco: Counter[str] = Counter()
    umag_shared_gems = 0
    umag_oids: set[str] = set()

    with zipfile.ZipFile(meta_zip, "r") as zf:
        iuig_name = next((n for n in zf.namelist() if "img_genomes_tarballs_metadata" in n), None)
        umag_name = next((n for n in zf.namelist() if "unrestricted_metagenomes_for_bin_tarballs" in n), None)

        if iuig_name:
            with zf.open(iuig_name) as f:
                reader = csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t")
                for row in reader:
                    iuig_total += 1
                    oid = (row.get("taxon_oid") or "").strip()
                    dom = (row.get("domain") or "").strip() or "Unclassified"
                    status = (row.get("seq_status") or "").strip() or "Unknown"
                    phylum = (row.get("phylum") or "").strip() or "Unclassified"
                    species = (row.get("species") or "").strip()
                    taxid = (row.get("ncbi_taxon_id") or "").strip()
                    iuig_by_domain[dom] += 1
                    iuig_by_status[status] += 1
                    if dom not in ("Bacteria", "Archaea"):
                        continue

                    gca = img_oid_to_gca.get(oid, "")
                    st = _stem(gca)
                    if st and (st in gtdb_stems or st in globdb_stems):
                        iuig_matched += 1
                        iuig_matched_by_gca += 1
                        continue
                    if taxid and taxid in gtdb_taxids:
                        iuig_matched += 1
                        iuig_matched_by_taxid += 1
                        continue
                    sp_toks = species.split()
                    sp_bin = f"{sp_toks[0]} {sp_toks[1]}" if len(sp_toks) >= 2 else species
                    if sp_bin and sp_bin in gtdb_species:
                        iuig_matched += 1
                        iuig_matched_by_species += 1
                        continue

                    iuig_prok_novel_oids.append(oid)
                    iuig_prok_phylum[phylum] += 1
                    iuig_prok_novel_rows.append({
                        "taxon_oid": oid,
                        "domain": dom,
                        "phylum": phylum,
                        "genus": (row.get("genus") or "").strip(),
                        "species": species,
                        "seq_status": status,
                        "analysis_project_id": (row.get("analysis_project_id") or "").strip(),
                        "ecosystem": img_oid_to_eco.get(oid, "Unclassified"),
                    })

        if umag_name:
            with zf.open(umag_name) as f:
                reader = csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t")
                for row in reader:
                    umag_total += 1
                    oid = (row.get("taxon_oid") or "").strip()
                    umag_oids.add(oid)
                    if oid in gems_mgs:
                        umag_shared_gems += 1
                    eco = (row.get("phylum") or "").strip() or "Unclassified"
                    eco_cat = (row.get("ir_class") or "").strip()
                    umag_by_eco[f"{eco} / {eco_cat}" if eco_cat else eco] += 1

    iuig_prok_novel_oids.sort()
    iuig_prok_novel_rows.sort(key=lambda x: x["taxon_oid"])

    return {
        "iuig": {
            "total_rows": iuig_total,
            "by_domain": dict(iuig_by_domain.most_common()),
            "by_seq_status": dict(iuig_by_status.most_common()),
            "prok_total": iuig_by_domain.get("Bacteria", 0) + iuig_by_domain.get("Archaea", 0),
            "prok_matched_gtdb_or_globdb": iuig_matched,
            "prok_matched_by_gca": iuig_matched_by_gca,
            "prok_matched_by_taxid": iuig_matched_by_taxid,
            "prok_matched_by_species": iuig_matched_by_species,
            "prok_novel_candidates_count": len(iuig_prok_novel_oids),
            "prok_novel_by_phylum": dict(iuig_prok_phylum.most_common(30)),
            "prok_novel_candidate_oids": iuig_prok_novel_oids,
            "prok_novel_rows": iuig_prok_novel_rows,
        },
        "umag": {
            "total_metagenomes": umag_total,
            "shared_with_gems_metagenomes": umag_shared_gems,
            "by_ecosystem": dict(umag_by_eco.most_common(25)),
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-48d metadata deduplication pass")
    ap.add_argument("--gems-meta-tar", type=Path, required=True)
    ap.add_argument("--umag-iuig-meta-zip", type=Path, required=True)
    ap.add_argument("--gtdb-tar", type=Path, required=True)
    ap.add_argument("--globdb-dl-dir", type=Path, required=True)
    ap.add_argument("--img-m-tar", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    refs = load_reference_sets(args.gtdb_tar, args.globdb_dl_dir, args.img_m_tar)
    gems_res = analyze_gems_metadata(args.gems_meta_tar, refs)
    refs["gems_metagenome_ids"] = set(gems_res["metagenome_ids"])
    iuig_umag_res = analyze_iuig_umag_metadata(args.umag_iuig_meta_zip, refs)

    # Write candidate ID lists deterministically
    gems_cand_tsv = args.out_dir / "gems_candidates.tsv"
    with open(gems_cand_tsv, "w", encoding="utf-8", newline="") as fh:
        cols = [
            "genome_id", "metagenome_id", "otu_id", "bp", "num_contigs", "n50",
            "completeness", "contamination", "quality_score", "mimag_quality",
            "domain", "phylum", "taxonomy", "ecosystem", "in_globdb_dict_dropped",
        ]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in gems_res["candidates"]:
            w.writerow(row)

    iuig_cand_tsv = args.out_dir / "iuig_prok_candidates.tsv"
    with open(iuig_cand_tsv, "w", encoding="utf-8", newline="") as fh:
        cols = ["taxon_oid", "domain", "phylum", "genus", "species", "seq_status", "analysis_project_id", "ecosystem"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in iuig_umag_res["iuig"]["prok_novel_rows"]:
            w.writerow(row)

    summary = {
        "reference_counts": refs["counts"],
        "gems": {k: v for k, v in gems_res.items() if k not in ("candidates", "metagenome_ids")},
        "iuig": {k: v for k, v in iuig_umag_res["iuig"].items() if k not in ("prok_novel_candidate_oids", "prok_novel_rows")},
        "umag": iuig_umag_res["umag"],
    }
    (args.out_dir / "metadata_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

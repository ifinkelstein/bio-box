#!/usr/bin/env python3
"""Card C-85 (DEC-88 (1), (2), DEC-80, AGENTS §8): prepare gentle next-wave `tower` fetch queues and
stream-extract target non-GTDB species subsets from `globdb_r232_genome_fasta.tar.gz`.

Subcommands:
    extract-globdb-subset --tar-gz FILE_OR_DASH --ids IDS_FILE --out-dir DIR
    inspect-pod           --out-json FILE [paths...]
    build-queues          --config-json FILE --out-dir DIR
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import sys
import tarfile
import zipfile
from pathlib import Path
from typing import Iterable

ACC_VER_RE = re.compile(r"^(GC[AF])_(\d{3})(\d{3})(\d{3})\.(\d+)$")
ACC_NOVER_RE = re.compile(r"(GC[AF]_\d{9})")
NCBI_FTP_DIR_RE = re.compile(
    r"^(?:https|ftp)://ftp\.ncbi\.nlm\.nih\.gov/genomes/all/"
    r"(GC[AF])/(\d{3})/(\d{3})/(\d{3})/([A-Za-z0-9_.-]+?)/?$"
)
WGS_PREFIX_RE = re.compile(r"^([A-Z]{4,6}\d{2})")

# Exact validation patterns mirrored from scripts/orchestrator/fetch_to_tower.sh lines 157-163
NAME_RE = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*(/[A-Za-z0-9_][A-Za-z0-9_.-]*)*$")
FILE_RE = re.compile(r"^[A-Za-z0-9_., ()=+-]+$")
URL_RE = re.compile(r"^(https|ftp)://[A-Za-z0-9._~:/?#@%&=+,;!*()-]+$")
MD5_RE = re.compile(r"^[0-9a-f]{32}$")
LIC_RE = re.compile(r"^[A-Za-z0-9._ {}()=+,/-]+$")

STREAM_CHUNK = 65_536


class Refused(Exception):
    """Raised when an input accession, URL, or job row violates repository rules."""


def nover(acc: str) -> str:
    """Normalize `RS_GCF_000005845.2`, `GB_GCA_053992195.1`, or `GCA_053992195` to `GC[AF]_ddddddddd`."""
    m = ACC_NOVER_RE.search(acc.strip())
    return m.group(1) if m else acc.strip()


def wgs_prefix(val: str) -> str | None:
    """Extract 4-6 uppercase letters + 2 digits WGS master prefix (e.g. `CAJVWK010000000.1` -> `CAJVWK01`)."""
    s = val.strip()
    if not s or s == "na" or s.startswith(("GCA_", "GCF_", "RS_", "GB_")):
        return None
    m = WGS_PREFIX_RE.match(s)
    return m.group(1) if m else None


def build_ncbi_fna_url(accession_or_ftp_path: str, asm_name: str | None = None) -> tuple[str, str]:
    """Build `(local_filename, https_url)` for an NCBI GenBank/RefSeq genomic FASTA (`.fna.gz`).

    Accepts either an NCBI `assembly_summary.txt` `ftp_path` (`ftp://` or `https://` under
    `ftp.ncbi.nlm.nih.gov/genomes/all/...`) or `(versioned_accession, asm_name)`.
    """
    raw = accession_or_ftp_path.strip()
    if not raw or raw.lower() == "na":
        raise Refused(f"empty or 'na' NCBI path/accession: {accession_or_ftp_path!r}")

    if "://" in raw:
        m = NCBI_FTP_DIR_RE.match(raw)
        if not m:
            raise Refused(f"malformed or non-NCBI assembly directory URL: {raw!r}")
        prefix, d1, d2, d3, asm_dir = m.groups()
        acc_prefix = f"{prefix}_{d1}{d2}{d3}."
        if not asm_dir.startswith(acc_prefix):
            raise Refused(f"assembly_dir {asm_dir!r} does not match prefix {acc_prefix!r}")
        ver_acc = asm_dir.split("_", 2)
        if len(ver_acc) < 2:
            raise Refused(f"assembly_dir missing versioned accession: {asm_dir!r}")
        acc_ver = f"{ver_acc[0]}_{ver_acc[1]}"
        if not ACC_VER_RE.match(acc_ver):
            raise Refused(f"invalid versioned accession {acc_ver!r} in {asm_dir!r}")
    else:
        m_acc = ACC_VER_RE.match(raw)
        if not m_acc:
            raise Refused(f"expected versioned GCA_/GCF_ accession, got {raw!r}")
        if not asm_name or not asm_name.strip() or asm_name.strip().lower() == "na":
            raise Refused(f"asm_name required when building URL from bare accession {raw!r}")
        prefix, d1, d2, d3, _ver = m_acc.groups()
        acc_ver = raw
        clean_asm = re.sub(r"[^A-Za-z0-9_.-]+", "_", asm_name.strip())
        asm_dir = f"{acc_ver}_{clean_asm}"

    local_file = f"{acc_ver}_genomic.fna.gz"
    https_url = (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{d1}/{d2}/{d3}/"
        f"{asm_dir}/{asm_dir}_genomic.fna.gz"
    )
    validate_job_row(("check", local_file, https_url, "-", "INSDC"))
    return local_file, https_url


def validate_job_row(row: tuple[str, str, str, str, str]) -> tuple[str, str, str, str, str]:
    """Validate a five-column `fetch_to_tower.sh` row `(name, file, url, md5, licence)`."""
    if len(row) != 5:
        raise Refused(f"expected 5 fields, got {len(row)}: {row!r}")
    name, file, url, md5, lic = (str(x) for x in row)
    if not NAME_RE.match(name) or any(part in (".", "..") for part in name.split("/")):
        raise Refused(f"refused destination name: {name!r}")
    if (
        not FILE_RE.match(file)
        or file.startswith("-")
        or file in (".", "..", "SOURCE_URLS.tsv", "SHA256SUMS")
        or file.endswith(".partial")
    ):
        raise Refused(f"refused file name: {file!r}")
    if not URL_RE.match(url):
        raise Refused(f"refused URL: {url!r}")
    if not (md5 == "-" or MD5_RE.match(md5)):
        raise Refused(f"refused md5 (expected '-' or 32 lowercase hex digits): {md5!r}")
    if not lic or not LIC_RE.match(lic):
        raise Refused(f"refused licence text: {lic!r}")
    return (name, file, url, md5, lic)


def write_jobs_tsv(
    path: Path,
    header_lines: Iterable[str],
    rows: Iterable[tuple[str, str, str, str, str]],
    *,
    sort_rows: bool = True,
) -> int:
    """Validate all rows, reject duplicate `(name, file)` targets, and write atomically."""
    validated = [validate_job_row(r) for r in rows]
    if sort_rows:
        validated.sort(key=lambda r: (r[0], r[1]))
    seen: set[tuple[str, str]] = set()
    for r in validated:
        key = (r[0], r[1])
        if key in seen:
            raise Refused(f"duplicate destination {r[0]}/{r[1]} in {path}")
        seen.add(key)
    if not validated:
        raise Refused(f"refusing to write empty job list: {path}")

    out_lines: list[str] = []
    for h in header_lines:
        h_clean = h.rstrip("\n")
        if h_clean and not h_clean.startswith("#"):
            raise Refused(f"header comment must start with '#': {h_clean!r}")
        out_lines.append(h_clean)
    for r in validated:
        out_lines.append("\t".join(r))

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.partial")
    tmp.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    tmp.replace(path)
    return len(validated)


def is_excluded_mgnify_biome(name_or_biome: str) -> bool:
    """Exclude human body-site catalogs and eukaryote catalogs (gated on card C-60)."""
    slug = Path(name_or_biome).name.lower()
    if slug.startswith("human-") or slug.startswith("human_"):
        return True
    if "eukaryote" in slug:
        return True
    return False


def _iter_mgnify_lines(source: str | Path) -> Iterable[str]:
    if isinstance(source, Path):
        with source.open("r", encoding="utf-8") as fh:
            yield from fh
    else:
        yield from io.StringIO(source)


def filter_mgnify_catalogs(
    catalog_texts: dict[str, str | Path],
    *,
    known_mgyg_ids: set[str],
    known_insdc_accessions: set[str],
    known_wgs_prefixes: set[str],
    known_biosamples: set[str] | None = None,
    url_suffix_override: str | None = None,
) -> tuple[list[dict], dict]:
    """Stream-filter MGnify `genomes-all_metadata.tsv` tables to novel non-human prokaryotic species reps."""
    biosamples = known_biosamples or set()
    kept: list[dict] = []
    seen_reps: set[str] = set()
    stats = {
        "total_catalogs": len(catalog_texts),
        "excluded_biomes": 0,
        "included_biomes": 0,
        "total_rows_scanned": 0,
        "species_reps_scanned": 0,
        "excluded_in_globdb_mgyg": 0,
        "excluded_in_insdc_or_wgs": 0,
        "excluded_in_biosample": 0,
        "excluded_duplicate_rep": 0,
        "kept_reps": 0,
        "kept_bp": 0,
        "by_biome": {},
    }

    for label in sorted(catalog_texts):
        if is_excluded_mgnify_biome(label):
            stats["excluded_biomes"] += 1
            continue
        stats["included_biomes"] += 1
        biome_key = Path(label).name.replace("_genomes-all_metadata.tsv", "")
        b_scanned = 0
        b_reps = 0
        b_kept = 0
        b_bp = 0

        it = iter(_iter_mgnify_lines(catalog_texts[label]))
        header_line = next(it, None)
        if header_line is None:
            continue
        cols = header_line.rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(cols)}
        for req in ("Genome", "Length", "Genome_accession", "Species_rep", "Lineage", "Sample_accession", "FTP_download"):
            if req not in idx:
                raise Refused(f"{label}: missing required column {req!r}")

        i_gen = idx["Genome"]
        i_len = idx["Length"]
        i_acc = idx["Genome_accession"]
        i_rep = idx["Species_rep"]
        i_lin = idx["Lineage"]
        i_sam = idx["Sample_accession"]
        i_ftp = idx["FTP_download"]

        for raw_line in it:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            b_scanned += 1
            stats["total_rows_scanned"] += 1
            genome = parts[i_gen].strip()
            sp_rep = parts[i_rep].strip()
            if genome != sp_rep:
                continue
            b_reps += 1
            stats["species_reps_scanned"] += 1

            if genome in seen_reps:
                stats["excluded_duplicate_rep"] += 1
                continue
            if genome in known_mgyg_ids:
                stats["excluded_in_globdb_mgyg"] += 1
                continue

            g_acc = parts[i_acc].strip()
            g_nover = nover(g_acc)
            g_wgs = wgs_prefix(g_acc)
            if (g_nover and g_nover in known_insdc_accessions) or (g_wgs and g_wgs in known_wgs_prefixes):
                stats["excluded_in_insdc_or_wgs"] += 1
                continue

            sam = parts[i_sam].strip()
            if sam and sam in biosamples:
                stats["excluded_in_biosample"] += 1
                continue

            ftp_url = parts[i_ftp].strip()
            if ftp_url.startswith("ftp://ftp.ebi.ac.uk/"):
                https_url = "https://ftp.ebi.ac.uk/" + ftp_url[len("ftp://ftp.ebi.ac.uk/") :]
            elif ftp_url.startswith("https://ftp.ebi.ac.uk/"):
                https_url = ftp_url
            else:
                raise Refused(f"{label}: unexpected FTP_download host for {genome}: {ftp_url!r}")

            if url_suffix_override is not None and https_url.endswith(".gff.gz"):
                https_url = https_url[: -len(".gff.gz")] + url_suffix_override

            length_bp = int(parts[i_len])
            seen_reps.add(genome)
            b_kept += 1
            b_bp += length_bp
            stats["kept_reps"] += 1
            stats["kept_bp"] += length_bp
            kept.append(
                {
                    "biome": biome_key,
                    "genome": genome,
                    "length": length_bp,
                    "genome_accession": g_acc,
                    "sample_accession": sam,
                    "lineage": parts[i_lin].strip(),
                    "https_url": https_url,
                    "filename": https_url.rsplit("/", 1)[-1],
                }
            )

        stats["by_biome"][biome_key] = {
            "rows_scanned": b_scanned,
            "species_reps": b_reps,
            "kept_reps": b_kept,
            "kept_bp": b_bp,
        }

    kept.sort(key=lambda r: (r["biome"], r["genome"]))
    return kept, stats


def select_gtdb_r232_delta(
    *,
    gtdb_r232_reps: Iterable[dict],
    all_gtdb_r232_accessions: set[str],
    r226_species: set[str],
    r226_rep_accessions: set[str],
    built_accessions: set[str],
    globdb_covered_accessions: set[str],
    post_r226_candidates: Iterable[dict] = (),
) -> dict:
    """Compute unique new GTDB R232 species representatives and post-R232 RefSeq/GenBank species delta."""
    new_species_total = 0
    excluded_built = 0
    excluded_globdb = 0
    unique_reps: list[dict] = []

    for r in gtdb_r232_reps:
        if not r.get("is_representative", False):
            continue
        acc = nover(r["accession"])
        gb_acc = nover(r.get("genbank_accession") or acc)
        sp = r.get("species")
        if (sp and sp in r226_species) or acc in r226_rep_accessions or gb_acc in r226_rep_accessions:
            continue
        new_species_total += 1
        if acc in built_accessions or gb_acc in built_accessions:
            excluded_built += 1
            continue
        if acc in globdb_covered_accessions or gb_acc in globdb_covered_accessions:
            excluded_globdb += 1
            continue
        unique_reps.append(
            {
                "accession": acc,
                "accession_versioned": r.get("accession_versioned", acc),
                "genbank_accession": gb_acc,
                "species": sp,
                "domain": r.get("domain"),
            }
        )

    unique_reps.sort(key=lambda x: x["accession"])
    chosen_accs = {x["accession"] for x in unique_reps} | {x["genbank_accession"] for x in unique_reps}

    post_delta: list[dict] = []
    post_scanned = 0
    post_in_gtdb232 = 0
    post_in_globdb = 0
    post_in_built = 0
    for cand in post_r226_candidates:
        post_scanned += 1
        raw_acc = cand["accession"].strip()
        acc = nover(raw_acc)
        if acc in all_gtdb_r232_accessions or acc in chosen_accs:
            post_in_gtdb232 += 1
            continue
        if acc in globdb_covered_accessions:
            post_in_globdb += 1
            continue
        if acc in built_accessions:
            post_in_built += 1
            continue
        item = dict(cand)
        item["accession"] = acc
        item["accession_versioned"] = raw_acc
        post_delta.append(item)

    post_delta.sort(key=lambda x: x["accession"])
    return {
        "gtdb_unique_new_reps": unique_reps,
        "post_r232_ncbi_delta": post_delta,
        "counts": {
            "new_r232_species_total": new_species_total,
            "excluded_already_built": excluded_built,
            "excluded_globdb_covered": excluded_globdb,
            "gtdb_unique_new_reps": len(unique_reps),
            "post_r226_scanned": post_scanned,
            "post_r226_in_gtdb_r232": post_in_gtdb232,
            "post_r226_in_globdb": post_in_globdb,
            "post_r226_in_built": post_in_built,
            "post_r232_ncbi_delta": len(post_delta),
        },
    }


def _member_genome_id(member_name: str) -> tuple[str, str]:
    base = member_name.rsplit("/", 1)[-1]
    stem = base
    for suf in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    return stem, base


def extract_globdb_subset(tar_gz_path: str | Path, target_ids: set[str], out_dir: Path) -> dict:
    """Stream a `.tar.gz` archive (`r|gz`) in bounded memory, extracting only members in `target_ids`."""
    if not target_ids:
        raise Refused("target_ids must be non-empty")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records: list[tuple[str, str, int, str]] = []
    found_ids: set[str] = set()
    total_bytes = 0

    if str(tar_gz_path) == "-":
        tf_ctx = tarfile.open(fileobj=sys.stdin.buffer, mode="r|gz")
    else:
        p = Path(tar_gz_path)
        if not p.is_file():
            raise Refused(f"missing archive: {p}")
        tf_ctx = tarfile.open(p, mode="r|gz")

    with tf_ctx as tf:
        for member in tf:
            if not member.isfile():
                continue
            gid, base = _member_genome_id(member.name)
            if gid not in target_ids or gid in found_ids:
                continue
            if not FILE_RE.match(base) or base.startswith("-"):
                raise Refused(f"unsafe member filename in archive: {member.name!r}")
            src = tf.extractfile(member)
            if src is None:
                continue
            dst = out_dir / base
            partial = out_dir / f".{base}.partial"
            hasher = hashlib.sha256()
            nbytes = 0
            with partial.open("wb") as out_fh:
                while True:
                    chunk = src.read(STREAM_CHUNK)
                    if not chunk:
                        break
                    out_fh.write(chunk)
                    hasher.update(chunk)
                    nbytes += len(chunk)
            partial.replace(dst)
            digest = hasher.hexdigest()
            found_ids.add(gid)
            total_bytes += nbytes
            records.append((gid, base, nbytes, digest))

    records.sort(key=lambda r: r[0])
    manifest_tmp = out_dir / ".extracted_manifest.tsv.partial"
    manifest_lines = ["genome_id\tfilename\tbytes\tsha256"]
    manifest_lines.extend(f"{gid}\t{base}\t{nb}\t{dig}" for gid, base, nb, dig in records)
    manifest_tmp.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    manifest_tmp.replace(out_dir / "extracted_manifest.tsv")

    missing = sorted(target_ids - found_ids)
    return {
        "total_targets": len(target_ids),
        "matched": len(found_ids),
        "missing_count": len(missing),
        "missing": missing,
        "total_bytes": total_bytes,
    }


def norm_bin(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def strip_clean_suffix(s: str) -> str:
    return re.sub(r"(?i)[_-]clean$", "", s.strip())


def map_asg1_study_bins(study_bins: Iterable[str], ncbi_genbank_rows: Iterable[dict]) -> tuple[list[dict], list[str]]:
    """Map `00_ThisStudy.tsv` Asgardarchaeota bin names to canonical `latest` GenBank `GCA_` assemblies."""
    bins = [b.strip() for b in study_bins if b.strip()]
    exact_map: dict[str, str] = {}
    norm_map: dict[str, str] = {}
    for orig in bins:
        st = orig.split("__", 1)[-1] if "__" in orig else orig
        exact_map[orig.lower()] = orig
        exact_map[st.lower()] = orig
        norm_map[norm_bin(orig)] = orig
        norm_map[norm_bin(st)] = orig

    hits: dict[str, list[dict]] = {}
    for row in ncbi_genbank_rows:
        if row.get("version_status", "latest") != "latest":
            continue
        ftp_path = row.get("ftp_path", "")
        if not ftp_path or ftp_path == "na":
            continue
        matched_orig = None
        for key in ("isolate", "infraspecific_name", "asm_name"):
            val = row.get(key, "")
            if not val or val == "na":
                continue
            for cand in (val.strip(), strip_clean_suffix(val)):
                cl = cand.lower()
                cn = norm_bin(cand)
                matched_orig = exact_map.get(cl) or (norm_map.get(cn) if len(cn) >= 5 else None)
                if matched_orig:
                    break
            if matched_orig:
                break
        if matched_orig:
            hits.setdefault(matched_orig, []).append(row)

    chosen: list[dict] = []
    seen_accs: set[str] = set()
    for bin_name in sorted(hits):
        cands = sorted(
            hits[bin_name],
            key=lambda x: (
                "texas" in str(x.get("submitter", "")).lower(),
                str(x.get("seq_rel_date", "")),
                str(x.get("acc_ver", "")),
            ),
            reverse=True,
        )
        best = dict(cands[0])
        acc = nover(best["acc_ver"])
        if acc not in seen_accs:
            seen_accs.add(acc)
            best["study_bin"] = bin_name
            best["accession"] = acc
            chosen.append(best)

    unmapped = sorted(set(bins) - set(hits))
    return chosen, unmapped


def _cmd_extract_globdb_subset(args: argparse.Namespace) -> int:
    ids_path = Path(args.ids)
    if not ids_path.is_file():
        raise Refused(f"missing --ids file: {ids_path}")
    targets = {
        line.split("\t", 1)[0].strip()
        for line in ids_path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }
    summary = extract_globdb_subset(args.tar_gz, targets, Path(args.out_dir))
    sys.stdout.write(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_ext = sub.add_parser("extract-globdb-subset", help="Stream-extract target genome IDs from GlobDB r232 tar.gz")
    p_ext.add_argument("--tar-gz", required=True, help="Path to globdb_r232_genome_fasta.tar.gz or '-' for stdin")
    p_ext.add_argument("--ids", required=True, help="File of target GlobDB genome IDs (first column / one per line)")
    p_ext.add_argument("--out-dir", required=True, help="Destination directory for extracted FASTA files")
    p_ext.set_defaults(func=_cmd_extract_globdb_subset)

    return ap


def main(argv: list[str] | None = None) -> int:
    ap = build_parser()
    try:
        args = ap.parse_args(argv)
        return int(args.func(args))
    except Refused as exc:
        sys.stderr.write(f"c85_prepare_tower_queues: {exc}\n")
        return 2


if __name__ == "__main__":
    sys.exit(main())

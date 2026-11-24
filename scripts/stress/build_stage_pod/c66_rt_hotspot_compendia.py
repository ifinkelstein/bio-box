"""Card C-66 (DEC-69, DEC-71): Stage post-GlobDB r226 long-read peatland/sludge/rumen (LR2)
and hadal/mat/oral/holobiont (NX1) RT-hotspot compendia.

Subcommands:
    fetch-lr2   STOR_DIR META_DIR [--gtdb-tar TAR] [--globdb-tax GZ] [--globdb-dict-tar TAR] [--include-mfd]
    stage-lr2   MANIFEST_JSON OUT_META_DIR OUT_FASTA_DIR [--final-fasta-dir DIR]
                [--root ROOT] [--root-stages ...] [--owner-inputs NAME=DIR ...]
                [--skani-dedup] [--skani-threads 8] [--workers 2]
    fetch-nx1   STOR_DIR META_DIR --gtdb-tar TAR --globdb-tax GZ [--globdb-dict-tar TAR]
    stage-nx1   MANIFEST_JSON OUT_META_DIR OUT_FASTA_DIR [--final-fasta-dir DIR]
                [--root ROOT] [--root-stages ...] [--owner-inputs NAME=DIR ...]
                [--skani-dedup] [--skani-threads 8] [--workers 2]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path

import duckdb
import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402
from biobox.l0seq import store  # noqa: E402
from biobox.l1 import genecall  # noqa: E402

Refused = make_inputs.Refused

UA = "curl/8.5.0"
MIN_STOR_FILE_BYTES = 64 * 1024 * 1024  # AGENTS.md §7.3

_CIRCULAR_RE = re.compile(
    r"(?:(?:^|[\s|;,])(?:circular\s*=\s*(?:true|yes|y|1)|circ\s*=\s*(?:true|yes|y|1)|"
    r"is_circular\s*=\s*(?:true|yes|y|1)|topology\s*=\s*circular|circular_contig|dtr\b))",
    re.IGNORECASE,
)
_CODE25_TAX_RE = re.compile(
    r"(gracilibacteria|absconditabacteria|absconditabacterales|bd1-5|sr1|jaedam01|patescibacteri)",
    re.IGNORECASE,
)
_CODE4_TAX_RE = re.compile(
    r"(hodgkinia|nasuia|zinderia|mycoplasma|ureaplasma|spiroplasma|mesoplasma|entomoplasma|"
    r"mycoplasmoides|mycoplasmopsis|metamycoplasma|malacoplasma|acholeplasma|phytoplasma)",
    re.IGNORECASE,
)

FASTA_EXTS = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def _sql_lit(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _sql_list(files: Iterable[Path]) -> str:
    return "[" + ", ".join(_sql_lit(f) for f in files) + "]"


def strip_acc_version(acc: str) -> str:
    acc = acc.strip()
    if acc.startswith(("RS_", "GB_")):
        acc = acc[3:]
    if "." in acc:
        return acc.split(".", 1)[0]
    return acc


def _is_fasta_name(name: str) -> bool:
    low = name.lower()
    base = os.path.basename(low)
    if base.startswith(".") or base.startswith("._"):
        return False
    return any(low.endswith(ext) for ext in FASTA_EXTS)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def _download_file(
    url: str,
    dest: Path,
    *,
    expected_md5: str | None = None,
    expected_bytes: int | None = None,
    min_interval: float = 0.55,
    last_req: list[float] | None = None,
    max_refusals: int = 5,
    timeout: int = 120,
) -> dict:
    """Download URL atomically with rate limit (<= 2 req/s), User-Agent curl, and MD5/SHA256 verification."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file():
        sz = dest.stat().st_size
        if expected_bytes is None or sz == expected_bytes:
            md5 = _md5_file(dest)
            if expected_md5 is None or md5 == expected_md5:
                sha = _sha256_file(dest)
                return {"url": url, "path": str(dest), "bytes": sz, "md5": md5, "sha256": sha, "cached": True}

    if last_req is None:
        last_req = [0.0]

    refusals = 0
    while True:
        wait = min_interval - (time.monotonic() - last_req[0])
        if wait > 0:
            time.sleep(wait)
        last_req[0] = time.monotonic()

        tmp = dest.with_name(dest.name + ".tmp")
        h_md5 = hashlib.md5()
        h_sha = hashlib.sha256()
        nbytes = 0
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as out:
                while chunk := resp.read(1 << 20):
                    out.write(chunk)
                    h_md5.update(chunk)
                    h_sha.update(chunk)
                    nbytes += len(chunk)
            md5 = h_md5.hexdigest()
            sha = h_sha.hexdigest()
            if expected_bytes is not None and nbytes != expected_bytes:
                tmp.unlink(missing_ok=True)
                raise Refused(f"{url}: size mismatch got {nbytes}, expected {expected_bytes}")
            if expected_md5 is not None and md5 != expected_md5:
                tmp.unlink(missing_ok=True)
                raise Refused(f"{url}: md5 mismatch got {md5}, expected {expected_md5}")
            tmp.replace(dest)
            return {"url": url, "path": str(dest), "bytes": nbytes, "md5": md5, "sha256": sha, "cached": False}
        except urllib.error.HTTPError as e:
            tmp.unlink(missing_ok=True)
            if e.code in (403, 429, 500, 502, 503, 504):
                refusals += 1
                if refusals >= max_refusals:
                    raise Refused(f"{url}: stopped after {refusals} consecutive refusals (last HTTP {e.code})") from e
                time.sleep(min(2.0 * refusals, 10.0))
                continue
            raise Refused(f"{url}: HTTP {e.code}") from e
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            tmp.unlink(missing_ok=True)
            refusals += 1
            if refusals >= max_refusals:
                raise Refused(f"{url}: stopped after {refusals} consecutive errors: {e}") from e
            time.sleep(min(2.0 * refusals, 10.0))


def load_gtdb_and_globdb_refs(
    gtdb_tar: Path | None,
    globdb_tax_gz: Path | None,
    globdb_dict_tar: Path | None = None,
) -> tuple[set[str], set[str], set[str]]:
    """Return (gtdb_accs, globdb_accs, held_species_tokens) from GTDB R226 and GlobDB r226."""
    gtdb_accs: set[str] = set()
    globdb_accs: set[str] = set()
    held_species: set[str] = set()

    if gtdb_tar is not None and gtdb_tar.is_file():
        with tarfile.open(gtdb_tar, "r:*") as tf:
            for m in tf:
                if not m.isfile() or not m.name.endswith(".tsv.gz"):
                    continue
                ef = tf.extractfile(m)
                if ef is None:
                    continue
                with gzip.open(ef, "rt", encoding="utf-8", errors="replace") as fh:
                    hdr = fh.readline().rstrip("\r\n").split("\t")
                    i_acc = hdr.index("accession")
                    i_gb = hdr.index("ncbi_genbank_assembly_accession") if "ncbi_genbank_assembly_accession" in hdr else -1
                    i_tax = hdr.index("gtdb_taxonomy") if "gtdb_taxonomy" in hdr else -1
                    for line in fh:
                        parts = line.rstrip("\r\n").split("\t")
                        if i_acc < len(parts) and parts[i_acc]:
                            gtdb_accs.add(strip_acc_version(parts[i_acc]))
                        if i_gb >= 0 and i_gb < len(parts) and parts[i_gb] and parts[i_gb] != "none":
                            gtdb_accs.add(strip_acc_version(parts[i_gb]))
                        if i_tax >= 0 and i_tax < len(parts):
                            for tok in parts[i_tax].split(";"):
                                tok = tok.strip()
                                if tok.startswith("s__") and len(tok) > 3:
                                    held_species.add(tok)

    if globdb_tax_gz is not None and globdb_tax_gz.is_file():
        with gzip.open(globdb_tax_gz, "rt", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline().rstrip("\r\n").split("\t")
            i_tax = next((i for i, c in enumerate(hdr) if "taxonomy" in c.lower()), -1)
            for line in fh:
                parts = line.rstrip("\r\n").split("\t")
                if not parts:
                    continue
                gid = parts[0].strip()
                if gid.startswith(("GCA_", "GCF_", "RS_GCF_", "GB_GCA_")):
                    globdb_accs.add(strip_acc_version(gid))
                if i_tax >= 0 and i_tax < len(parts):
                    for tok in parts[i_tax].split(";"):
                        tok = tok.strip()
                        if tok.startswith("s__") and len(tok) > 3:
                            held_species.add(tok)

    if globdb_dict_tar is not None and globdb_dict_tar.is_file():
        with tarfile.open(globdb_dict_tar, "r:*") as tf:
            for m in tf:
                if not m.isfile():
                    continue
                ef = tf.extractfile(m)
                if ef is None:
                    continue
                for raw in ef:
                    for tok in raw.decode("utf-8", "replace").split():
                        if tok.startswith(("GCA_", "GCF_", "RS_GCF_", "GB_GCA_")):
                            globdb_accs.add(strip_acc_version(tok))

    return gtdb_accs, globdb_accs, held_species


def _extract_species_token(taxonomy: str) -> str:
    for tok in taxonomy.split(";"):
        tok = tok.strip()
        if tok.startswith("s__") and len(tok) > 3:
            return tok
    return ""


def _safe_float(val: str | None, default: float = 0.0) -> float:
    if val is None or val == "" or val.lower() in ("na", "nan", "none"):
        return default
    try:
        return float(val)
    except ValueError:
        return default


def _is_truthy(val: str | None) -> bool:
    return (val or "").strip().lower() in ("true", "t", "1", "yes", "y", "rep", "representative")


def filter_emerge_v2_metadata(
    tsv_path: Path,
    *,
    held_species: set[str],
    min_completeness: float = 50.0,
    max_contamination: float = 10.0,
) -> tuple[list[dict], list[dict]]:
    """Filter EMERGE v2.0.1 metadata TSV to species-dereplicated novel MIMAG MQ/HQ representatives.

    Handles flexible column naming across EMERGE metadata releases.
    """
    selected: list[dict] = []
    deferred: list[dict] = []

    with open(tsv_path, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        col_low = {c.lower(): c for c in cols}

        c_genome = next((col_low[k] for k in ("genome", "mag", "mag_id", "genome_id", "user_genome", "name") if k in col_low), cols[0])
        c_comp = next((col_low[k] for k in ("completeness", "checkm2_completeness", "checkm_completeness") if k in col_low), None)
        c_cont = next((col_low[k] for k in ("contamination", "checkm2_contamination", "checkm_contamination") if k in col_low), None)
        c_tax = next((col_low[k] for k in ("gtdb_taxonomy", "taxonomy", "classification", "gtdb2_taxonomy", "gtdb_r220_taxonomy", "gtdb_r214_taxonomy") if k in col_low), None)
        c_rep = next((col_low[k] for k in ("is_representative", "representative", "is_rep", "species_rep", "cluster_rep") if k in col_low), None)
        c_clust = next((col_low[k] for k in ("species_cluster", "cluster", "secondary_cluster", "rep_genome") if k in col_low), None)

        rows = list(reader)

    # If no explicit is_representative flag exists, pick best-quality genome per species_cluster
    best_in_cluster: dict[str, str] = {}
    if c_rep is None and c_clust is not None:
        by_cl: dict[str, list[dict]] = defaultdict(list)
        for r in rows:
            cl = (r.get(c_clust) or "").strip() or (r.get(c_genome) or "").strip()
            by_cl[cl].append(r)
        for cl, grp in by_cl.items():
            grp.sort(
                key=lambda x: (
                    -(_safe_float(x.get(c_comp)) - 5.0 * _safe_float(x.get(c_cont))),
                    (x.get(c_genome) or "").strip(),
                )
            )
            best_in_cluster[cl] = (grp[0].get(c_genome) or "").strip()

    for r in rows:
        gid = (r.get(c_genome) or "").strip()
        if not gid:
            continue
        comp = _safe_float(r.get(c_comp)) if c_comp else 100.0
        cont = _safe_float(r.get(c_cont)) if c_cont else 0.0
        tax = (r.get(c_tax) or "").strip() if c_tax else ""
        cl = (r.get(c_clust) or "").strip() if c_clust else gid

        if c_rep is not None:
            is_rep = _is_truthy(r.get(c_rep))
        elif c_clust is not None:
            is_rep = (best_in_cluster.get(cl) == gid)
        else:
            is_rep = True

        if not is_rep:
            deferred.append({"genome": gid, "reason": "not_representative", "taxonomy": tax, "completeness": comp, "contamination": cont})
            continue
        if comp < min_completeness or cont >= max_contamination:
            deferred.append({"genome": gid, "reason": "below_mimag_mq", "taxonomy": tax, "completeness": comp, "contamination": cont})
            continue
        sp = _extract_species_token(tax)
        if sp and sp in held_species:
            deferred.append({"genome": gid, "reason": "species_in_gtdb_or_globdb", "taxonomy": tax, "completeness": comp, "contamination": cont})
            continue

        qual = "mag_hq" if (comp > 90.0 and cont < 5.0) else "mag_mq"
        selected.append({
            "genome": gid,
            "completeness": comp,
            "contamination": cont,
            "assembly_quality": qual,
            "taxonomy": tax,
            "species_cluster": cl,
        })

    selected.sort(key=lambda x: x["genome"])
    return selected, deferred


def parse_skani_matches(
    tsv_path: Path,
    *,
    min_ani: float = 96.0,
    min_af: float = 50.0,
) -> dict[str, dict]:
    """Parse skani search/dist TSV and return query_name -> best hit with ANI >= min_ani and max(AF_ref, AF_query) >= min_af."""
    best: dict[str, dict] = {}
    with open(tsv_path, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            qname = (r.get("Query_name") or Path(r.get("Query_file", "")).stem).strip()
            rname = (r.get("Ref_name") or Path(r.get("Ref_file", "")).stem).strip()
            ani = _safe_float(r.get("ANI"))
            af_ref = _safe_float(r.get("Align_fraction_ref"))
            af_q = _safe_float(r.get("Align_fraction_query"))
            af = max(af_ref, af_q)
            if ani >= min_ani and af >= min_af:
                prev = best.get(qname)
                if prev is None or (ani, af) > (prev["ani"], prev["af"]):
                    best[qname] = {"ref": rname, "ani": ani, "af": af, "af_ref": af_ref, "af_query": af_q}
    return best


def filter_ehomd_novel(
    tsv_path: Path,
    *,
    gtdb_accs: set[str],
    globdb_accs: set[str],
) -> tuple[list[dict], list[dict]]:
    """Filter eHOMD GCA_ID_info.txt to assemblies absent from GTDB R226 and GlobDB r226."""
    novel: list[dict] = []
    skipped: list[dict] = []
    seen: set[str] = set()
    with open(tsv_path, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            gca = (r.get("GCA-ID") or "").strip()
            if not gca or gca in seen:
                continue
            seen.add(gca)
            base = strip_acc_version(gca)
            status = (r.get("Status") or "").strip()
            tax = f"g__{r.get('Genus', '').strip()};s__{r.get('Species', '').strip()}"
            item = {
                "gca_id": gca,
                "base_acc": base,
                "hmt_id": (r.get("HMT-ID") or "").strip(),
                "taxonomy": tax,
                "strain": (r.get("Strain") or "").strip(),
                "status": status,
                "combined_size": int(_safe_float(r.get("Combined_Size"))),
            }
            if base in gtdb_accs or base in globdb_accs:
                skipped.append({**item, "reason": "in_gtdb_or_globdb"})
            else:
                novel.append(item)
    novel.sort(key=lambda x: x["gca_id"])
    return novel, skipped


def detect_exact_dtr(seq: bytes, min_k: int = 20, max_k: int = 1000) -> int:
    """Return length k of exact terminal direct repeat (20..1,000 bp) per DEC-65, or 0 if none/low-complexity."""
    n = len(seq)
    if n < 2 * min_k + 1:
        return 0
    upper = seq.upper()
    hi = min(max_k, n // 2 - 1)
    # Search from longest to shortest exact prefix==suffix match
    for k in range(hi, min_k - 1, -1):
        if upper[:k] == upper[n - k :]:
            rep = upper[:k]
            # Reject homopolymer / dinucleotide low-complexity repeats
            if len(set(rep)) <= 2:
                return 0
            return k
    return 0


def _call_stats(seq_strs: list[str], total_bp: int, table: int) -> tuple[float, float, int]:
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*seq_strs, translation_table=table)
    coding_bp = 0
    n_genes = 0
    for s in seq_strs:
        for g in finder.find_genes(s):
            coding_bp += g.end - g.begin + 1
            n_genes += 1
    cd = coding_bp / total_bp if total_bp > 0 else 0.0
    mean_orf = coding_bp / n_genes if n_genes > 0 else 0.0
    return cd, mean_orf, n_genes


def screen_genetic_code(
    seqs: list[bytes],
    taxonomy: str = "",
    *,
    cd_gain_min: float = 0.10,
    orf_ratio_min: float = 1.25,
    cd11_trigger: float = 0.82,
    orf11_trigger: float = 600.0,
    screen_max_bp: int = 500_000,
) -> dict:
    """Screen a prokaryotic genome (>= 100 kb) for UGA recoding (code 11 vs 25 vs 4) using pyrodigal."""
    total_bp = sum(len(s) for s in seqs)
    if total_bp < genecall.SINGLE_MODE_MIN_BP:
        return {"genetic_code": 11, "reason": "below_100kb", "tested_alt": False, "cd_11": 0.0, "mean_orf_11": 0.0}

    ordered = sorted(seqs, key=len, reverse=True)
    pre_strs: list[str] = []
    pre_bp = 0
    for s in ordered:
        pre_strs.append(s.upper().decode("ascii", "replace"))
        pre_bp += len(s)
        if screen_max_bp > 0 and pre_bp >= screen_max_bp and pre_bp >= genecall.SINGLE_MODE_MIN_BP:
            break

    cd_11, mean_orf_11, _ = _call_stats(pre_strs, pre_bp, 11)
    force_tax = bool(_CODE25_TAX_RE.search(taxonomy) or _CODE4_TAX_RE.search(taxonomy))
    if not (force_tax or cd_11 < cd11_trigger or mean_orf_11 < orf11_trigger):
        return {
            "genetic_code": 11,
            "reason": "code11_pass",
            "tested_alt": False,
            "cd_11": round(cd_11, 5),
            "mean_orf_11": round(mean_orf_11, 2),
        }

    cd_25, mean_orf_25, _ = _call_stats(pre_strs, pre_bp, 25)
    cd_4, mean_orf_4, _ = _call_stats(pre_strs, pre_bp, 4)

    gain_25 = cd_25 - cd_11
    ratio_25 = (mean_orf_25 / mean_orf_11) if mean_orf_11 > 0 else 0.0
    gain_4 = cd_4 - cd_11
    ratio_4 = (mean_orf_4 / mean_orf_11) if mean_orf_11 > 0 else 0.0

    recoded = (max(gain_25, gain_4) >= cd_gain_min) and (max(ratio_25, ratio_4) >= orf_ratio_min)
    if not recoded:
        return {
            "genetic_code": 11,
            "reason": "code11_alt_checked",
            "tested_alt": True,
            "cd_11": round(cd_11, 5),
            "mean_orf_11": round(mean_orf_11, 2),
            "cd_25": round(cd_25, 5),
            "cd_4": round(cd_4, 5),
        }

    if _CODE4_TAX_RE.search(taxonomy) and not _CODE25_TAX_RE.search(taxonomy):
        code, reason = 4, "code4_tax_and_gain"
    elif _CODE25_TAX_RE.search(taxonomy) and not _CODE4_TAX_RE.search(taxonomy):
        code, reason = 25, "code25_tax_and_gain"
    elif (cd_25, mean_orf_25) >= (cd_4, mean_orf_4):
        code, reason = 25, "code25_gain"
    else:
        code, reason = 4, "code4_gain"

    return {
        "genetic_code": code,
        "reason": reason,
        "tested_alt": True,
        "cd_11": round(cd_11, 5),
        "mean_orf_11": round(mean_orf_11, 2),
        "cd_25": round(cd_25, 5),
        "cd_4": round(cd_4, 5),
    }


def _iter_fasta_file(path: Path) -> Iterable[tuple[str, bytes]]:
    for rec in fasta.iter_records(path):
        yield rec.header, rec.seq


def _iter_fasta_stream(fh: io.BufferedIOBase | io.RawIOBase) -> Iterable[tuple[str, bytes]]:
    header: str | None = None
    seq_lines: list[bytes] = []
    for raw in fh:
        if not raw.strip():
            continue
        ln = raw.rstrip(b"\r\n")
        if ln.startswith(b">"):
            if header is not None:
                yield header, b"".join(seq_lines)
                seq_lines = []
            header = ln[1:].decode("utf-8", "replace")
        else:
            if header is not None:
                seq_lines.append(ln)
    if header is not None:
        yield header, b"".join(seq_lines)


def _load_earlier_contig_owners(
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path] | None,
    *,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict[bytes, str]:
    """Map contig_uid (16 bytes) -> earliest owning stage across built root stages and staged owner_inputs."""
    owners: dict[bytes, str] = {}
    if root is not None:
        stages = list(root_stages)
        if not stages:
            ci_dir = root / "l0seq" / "contig_index"
            if ci_dir.is_dir():
                stages = sorted(p.name.split("=", 1)[1] for p in ci_dir.glob("stage=*") if p.is_dir())
        if stages:
            con = duckdb.connect()
            con.execute(f"SET threads={int(threads)}")
            con.execute(f"SET memory_limit='{memory_limit}'")
            for st in stages:
                files = sorted((root / "l0seq" / "contig_index" / f"stage={st}").glob("*.parquet"))
                if not files:
                    continue
                for (uid,) in con.execute(f"SELECT contig_uid FROM read_parquet({_sql_list(files)})").fetchall():
                    b = bytes(uid)
                    if b not in owners:
                        owners[b] = st
            con.close()

    for st, d in (owner_inputs or {}).items():
        gjson = Path(d) / "genomes.json"
        if not gjson.is_file():
            continue
        with open(gjson, encoding="utf-8") as fh:
            for g in json.load(fh):
                for h in g.get("contig_uids", ()):
                    b = bytes.fromhex(h)
                    if b not in owners:
                        owners[b] = st
    return owners


def stage_compendium(
    *,
    stage: str,
    mag_entries: list[dict],
    unbinned_entries: list[dict],
    out_meta_dir: Path,
    out_fasta_dir: Path,
    final_fasta_dir: Path | None,
    fasta_name: str,
    source_meta: dict,
    min_unbinned_len: int = 5_000,
    root: Path | None = None,
    root_stages: Iterable[str] = (),
    owner_inputs: dict[str, Path] | None = None,
    retrieved_at: str = "2026-11-24T00:00:00Z",
    screen_codes: bool = True,
    screen_max_bp: int = 500_000,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict:
    """Consolidate MAG bins and optional unbinned >= 5 kb contigs into one deterministic .fna.gz stage."""
    out_meta_dir.mkdir(parents=True, exist_ok=True)
    out_fasta_dir.mkdir(parents=True, exist_ok=True)

    earlier_owners = _load_earlier_contig_owners(
        root, root_stages, owner_inputs, threads=threads, memory_limit=memory_limit
    )

    fasta_path = out_fasta_dir / fasta_name
    tmp_fasta = fasta_path.with_name(fasta_path.name + ".tmp")

    written_uids: set[bytes] = set()
    binned_uids: set[bytes] = set()
    contig_props: dict[bytes, tuple[str, int, bool]] = {}  # uid -> (kind, code, circular)
    seen_genome_keys: set[tuple[bytes, ...]] = set()

    genomes_out: list[dict] = []
    records_rows: list[list] = []
    repeated_rows: list[list] = []
    trims_rows: list[list] = []
    gc_rows: list[list] = []

    counts = Counter()
    bases_by_source = Counter()
    genomes_by_source = Counter()
    code_counts = Counter()

    with gzip.open(tmp_fasta, "wb", compresslevel=6) as gz_out:
        # 1. Process MAG bins (sorted deterministically by genome_id)
        for entry in sorted(mag_entries, key=lambda x: x["genome_id"]):
            gid = entry["genome_id"]
            src = entry["source_name"]
            fpath = Path(entry["path"])
            qual = entry.get("assembly_quality", "mag_mq")
            tax = entry.get("taxonomy", "")
            kind = entry.get("kind", "prokaryote")

            raw_recs: list[tuple[str, bytes]] = []
            for hdr, seq in _iter_fasta_file(fpath):
                s = seq.strip()
                if not s:
                    continue
                raw_recs.append((hdr, s))
            if not raw_recs:
                counts["empty_mags_skipped"] += 1
                continue

            if "genetic_code" in entry and not screen_codes:
                code = int(entry["genetic_code"])
                gc_info = {"genetic_code": code, "reason": "preset", "tested_alt": False, "cd_11": 0.0, "mean_orf_11": 0.0}
            elif screen_codes and kind == "prokaryote":
                gc_info = screen_genetic_code([s for _, s in raw_recs], tax, screen_max_bp=screen_max_bp)
                code = int(gc_info["genetic_code"])
            else:
                code = int(entry.get("genetic_code", 11))
                gc_info = {"genetic_code": code, "reason": "default", "tested_alt": False, "cd_11": 0.0, "mean_orf_11": 0.0}

            g_uids: list[bytes] = []
            g_circ: list[bytes] = []
            g_seen: set[bytes] = set()
            g_bp = 0
            conflict = False

            # First pass: determine uids and check for within-stage DEC-27 property conflict
            processed_contigs: list[tuple[str, bytes, bytes, bool]] = []
            for hdr, seq in raw_recs:
                is_circ_flag = bool(_CIRCULAR_RE.search(hdr)) and len(raw_recs) == 1
                trimmed_seq = seq
                is_circ = False
                if is_circ_flag:
                    k = detect_exact_dtr(seq)
                    if k > 0:
                        trimmed_seq = seq[:-k]
                        is_circ = True
                        trims_rows.append([gid, hdr.split()[0], len(seq), k, len(trimmed_seq)])
                        counts["dtr_trimmed"] += 1
                    else:
                        is_circ = True
                uid = ids.contig_uid(trimmed_seq)
                if uid in g_seen:
                    counts["within_mag_dup_contigs"] += 1
                    continue
                g_seen.add(uid)
                if uid in contig_props:
                    prev_kind, prev_code, prev_circ = contig_props[uid]
                    if (prev_kind, prev_code, prev_circ) != (kind, code, is_circ):
                        # Align to first-seen property if only code differs on shared contig, else skip genome
                        conflict = True
                        break
                processed_contigs.append((hdr, trimmed_seq, uid, is_circ))

            if conflict:
                counts["dec27_property_conflicts_skipped"] += 1
                repeated_rows.append([gid, src, "dec27_property_conflict", ""])
                continue

            g_key = tuple(sorted(uid for _, _, uid, _ in processed_contigs))
            if g_key in seen_genome_keys:
                counts["repeated_genomes"] += 1
                repeated_rows.append([gid, src, "duplicate_genome_contig_set", g_key[0].hex()])
                continue
            seen_genome_keys.add(g_key)

            for rec_ord, (hdr, seq, uid, is_circ) in enumerate(processed_contigs):
                g_uids.append(uid)
                binned_uids.add(uid)
                if is_circ:
                    g_circ.append(uid)
                    counts["circular_contigs"] += 1
                if uid not in written_uids:
                    written_uids.add(uid)
                    contig_props[uid] = (kind, code, is_circ)
                    gz_out.write(b">" + uid.hex().encode("ascii") + b"\n")
                    for i in range(0, len(seq), 80):
                        gz_out.write(seq[i : i + 80] + b"\n")
                    bases_by_source[src] += len(seq)
                g_bp += len(seq)
                records_rows.append([gid, src, rec_ord, uid.hex(), len(seq), int(is_circ), kind, qual, code, tax])

            genomes_out.append({
                "contig_uids": [u.hex() for u in g_uids],
                "genetic_code": code,
                "kind": kind,
                "circular": [u.hex() for u in g_circ],
                "assembly_quality": qual,
            })
            counts["mag_genomes"] += 1
            genomes_by_source[src] += 1
            code_counts[code] += 1
            gc_rows.append([
                gid, src, g_bp, len(g_uids), qual, code,
                gc_info.get("reason", ""), gc_info.get("tested_alt", False),
                gc_info.get("cd_11", ""), gc_info.get("mean_orf_11", ""),
                gc_info.get("cd_25", ""), gc_info.get("cd_4", ""),
            ])

        # 2. Process unbinned >= min_unbinned_len contigs (e.g. Bickhart HiFi unbinned assembly)
        for u_entry in sorted(unbinned_entries, key=lambda x: (x["source_name"], x["path"])):
            src = u_entry["source_name"]
            fpath = Path(u_entry["path"])
            for rec_ord, (hdr, seq) in enumerate(_iter_fasta_file(fpath)):
                s = seq.strip()
                if len(s) < min_unbinned_len:
                    counts["unbinned_short_dropped"] += 1
                    counts["unbinned_short_bases"] += len(s)
                    continue
                is_circ_flag = bool(_CIRCULAR_RE.search(hdr))
                trimmed_seq = s
                is_circ = False
                if is_circ_flag:
                    k = detect_exact_dtr(s)
                    if k > 0:
                        trimmed_seq = s[:-k]
                        is_circ = True
                        trims_rows.append([f"{src}:unbinned:{rec_ord}", hdr.split()[0], len(s), k, len(trimmed_seq)])
                        counts["dtr_trimmed"] += 1
                    else:
                        is_circ = True
                if len(trimmed_seq) < min_unbinned_len:
                    counts["unbinned_short_dropped"] += 1
                    continue
                uid = ids.contig_uid(trimmed_seq)
                if uid in binned_uids:
                    counts["unbinned_already_binned_skipped"] += 1
                    continue
                if uid in written_uids:
                    counts["unbinned_dup_skipped"] += 1
                    continue
                g_key = (uid,)
                if g_key in seen_genome_keys:
                    counts["unbinned_dup_skipped"] += 1
                    continue
                seen_genome_keys.add(g_key)
                written_uids.add(uid)
                contig_props[uid] = ("prokaryote", 11, is_circ)
                gz_out.write(b">" + uid.hex().encode("ascii") + b"\n")
                for i in range(0, len(trimmed_seq), 80):
                    gz_out.write(trimmed_seq[i : i + 80] + b"\n")
                if is_circ:
                    counts["circular_contigs"] += 1
                bases_by_source[src] += len(trimmed_seq)
                genomes_out.append({
                    "contig_uids": [uid.hex()],
                    "genetic_code": 11,
                    "kind": "prokaryote",
                    "circular": [uid.hex()] if is_circ else [],
                    "assembly_quality": "contig",
                })
                counts["unbinned_contigs_kept"] += 1
                genomes_by_source[src] += 1
                code_counts[11] += 1
                records_rows.append([
                    f"{src}:unbinned:{rec_ord}", src, 0, uid.hex(), len(trimmed_seq),
                    int(is_circ), "prokaryote", "contig", 11, "",
                ])

    tmp_fasta.replace(fasta_path)
    fasta_sha256 = _sha256_file(fasta_path)
    fasta_bytes = fasta_path.stat().st_size

    # Compute contig_uid overlap against earlier stages
    shared_rows: list[list] = []
    shared_by_stage = Counter()
    shared_bases_by_stage = Counter()
    total_bp = sum(bases_by_source.values())
    uid_len: dict[bytes, int] = {}
    for row in records_rows:
        u = bytes.fromhex(row[3])
        if u not in uid_len:
            uid_len[u] = int(row[4])

    for uid in sorted(written_uids):
        owner = earlier_owners.get(uid)
        if owner is not None:
            ln = uid_len.get(uid, 0)
            shared_rows.append([uid.hex(), owner, ln])
            shared_by_stage[owner] += 1
            shared_bases_by_stage[owner] += ln

    # Write genomes.json, objects_staging.json, objects.json
    genomes_path = out_meta_dir / "genomes.json"
    genomes_path.write_text(json.dumps(genomes_out, separators=(",", ":")) + "\n", encoding="utf-8")

    obj_staging = [{
        "path": str(fasta_path),
        "source_name": source_meta["source_name"],
        "source_release": source_meta["source_release"],
        "license": source_meta["license"],
        "redistributable": bool(source_meta["redistributable"]),
        "retrieved_at": retrieved_at,
    }]
    obj_staging_path = out_meta_dir / "objects_staging.json"
    obj_staging_path.write_text(json.dumps(obj_staging, indent=2) + "\n", encoding="utf-8")

    final_path = (final_fasta_dir / fasta_name) if final_fasta_dir is not None else fasta_path
    obj_final = [{**obj_staging[0], "path": str(final_path)}]
    obj_final_path = out_meta_dir / "objects.json"
    obj_final_path.write_text(json.dumps(obj_final, indent=2) + "\n", encoding="utf-8")

    def _write_tsv(path: Path, header: list[str], rows: list[list]) -> str:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)
        return _sha256_file(path)

    rec_sha = _write_tsv(
        out_meta_dir / "records.tsv",
        ["genome_id", "source_name", "record_ordinal", "contig_uid", "length", "circular", "kind", "assembly_quality", "genetic_code", "taxonomy"],
        records_rows,
    )
    rep_sha = _write_tsv(
        out_meta_dir / "repeated_records.tsv",
        ["genome_id", "source_name", "reason", "first_contig_uid"],
        repeated_rows,
    )
    shr_sha = _write_tsv(
        out_meta_dir / "shared_records.tsv",
        ["contig_uid", "owner_stage", "length"],
        shared_rows,
    )
    trm_sha = _write_tsv(
        out_meta_dir / "trims.tsv",
        ["genome_id", "contig_name", "raw_len", "trimmed_k", "final_len"],
        trims_rows,
    )
    gc_sha = _write_tsv(
        out_meta_dir / "genetic_codes.tsv",
        ["genome_id", "source_name", "genome_bp", "n_contigs", "assembly_quality", "genetic_code", "reason", "tested_alt", "cd_11", "mean_orf_11", "cd_25", "cd_4"],
        gc_rows,
    )

    shared_contigs = len(shared_rows)
    shared_bp = sum(shared_bases_by_stage.values())

    report = {
        "stage": stage,
        "fasta_path": str(fasta_path),
        "fasta_bytes": fasta_bytes,
        "fasta_sha256": fasta_sha256,
        "genomes_json_sha256": _sha256_file(genomes_path),
        "objects_json_sha256": _sha256_file(obj_final_path),
        "objects_staging_json_sha256": _sha256_file(obj_staging_path),
        "records_tsv_sha256": rec_sha,
        "repeated_records_tsv_sha256": rep_sha,
        "shared_records_tsv_sha256": shr_sha,
        "trims_tsv_sha256": trm_sha,
        "genetic_codes_tsv_sha256": gc_sha,
        "total_genomes": len(genomes_out),
        "mag_genomes": counts["mag_genomes"],
        "unbinned_contigs_kept": counts["unbinned_contigs_kept"],
        "unbinned_short_dropped": counts["unbinned_short_dropped"],
        "unbinned_short_bases": counts["unbinned_short_bases"],
        "unbinned_already_binned_skipped": counts["unbinned_already_binned_skipped"],
        "unbinned_dup_skipped": counts["unbinned_dup_skipped"],
        "repeated_genomes": counts["repeated_genomes"],
        "dec27_property_conflicts_skipped": counts["dec27_property_conflicts_skipped"],
        "dtr_trimmed": counts["dtr_trimmed"],
        "circular_contigs": counts["circular_contigs"],
        "distinct_contigs": len(written_uids),
        "total_bp": total_bp,
        "shared_contigs": shared_contigs,
        "shared_bp": shared_bp,
        "net_new_contigs": len(written_uids) - shared_contigs,
        "net_new_bp": total_bp - shared_bp,
        "shared_by_stage": dict(shared_by_stage),
        "shared_bases_by_stage": dict(shared_bases_by_stage),
        "genomes_by_source": dict(genomes_by_source),
        "bases_by_source": dict(bases_by_source),
        "genetic_codes": {str(k): v for k, v in sorted(code_counts.items())},
    }
    (out_meta_dir / "stage_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def run_skani_self_dedup(
    mag_entries: list[dict],
    work_dir: Path,
    *,
    min_ani: float = 96.0,
    min_af: float = 50.0,
    threads: int = 8,
) -> tuple[list[dict], list[dict]]:
    """Run skani triangle on MAG bins at <= 8 threads (nice 19) to dereplicate within-stage >= 96% ANI / 50% AF."""
    if len(mag_entries) <= 1:
        return mag_entries, []
    work_dir.mkdir(parents=True, exist_ok=True)
    flist = work_dir / "skani_files.txt"
    path_to_entry: dict[str, dict] = {}
    ordered_entries = sorted(
        mag_entries,
        key=lambda e: (
            0 if e.get("assembly_quality") == "mag_hq" else 1,
            -int(e.get("size_bytes", 0)),
            e["genome_id"],
        ),
    )
    with open(flist, "w", encoding="utf-8") as fh:
        for e in ordered_entries:
            p = str(Path(e["path"]).resolve())
            path_to_entry[p] = e
            fh.write(p + "\n")

    out_tsv = work_dir / "skani_triangle.tsv"
    cmd = [
        "skani", "triangle",
        "-l", str(flist),
        "-t", str(threads),
        "--sparse",
        "--min-af", str(min_af),
        "-o", str(out_tsv),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise Refused(f"skani triangle failed (rc={proc.returncode}): {proc.stderr[:500]}")

    # Greedy clustering in priority order (HQ first, largest first)
    rank = {str(Path(e["path"]).resolve()): idx for idx, e in enumerate(ordered_entries)}
    adj: dict[str, set[str]] = defaultdict(set)
    with open(out_tsv, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rf = str(Path(r.get("Ref_file", "")).resolve())
            qf = str(Path(r.get("Query_file", "")).resolve())
            if rf == qf or rf not in rank or qf not in rank:
                continue
            ani = _safe_float(r.get("ANI"))
            af = max(_safe_float(r.get("Align_fraction_ref")), _safe_float(r.get("Align_fraction_query")))
            if ani >= min_ani and af >= min_af:
                adj[rf].add(qf)
                adj[qf].add(rf)

    kept: list[dict] = []
    dropped: list[dict] = []
    eliminated: set[str] = set()
    for e in ordered_entries:
        p = str(Path(e["path"]).resolve())
        if p in eliminated:
            continue
        kept.append(e)
        for nb in adj.get(p, ()):
            if nb not in eliminated and rank[nb] > rank[p]:
                eliminated.add(nb)
                dropped.append({
                    "genome_id": path_to_entry[nb]["genome_id"],
                    "source_name": path_to_entry[nb]["source_name"],
                    "reason": f"skani_within_stage_ge_{min_ani}ani_rep:{e['genome_id']}",
                })
    return kept, dropped


def _stream_extract_tar_gz(
    url: str,
    wanted_predicate,
    extract_dir: Path,
    *,
    expected_md5: str | None = None,
    expected_bytes: int | None = None,
    timeout: int = 300,
) -> dict:
    """Stream a remote .tar.gz in-flight without writing the full archive to disk, extracting matching members."""
    extract_dir.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    h_md5 = hashlib.md5()
    h_sha = hashlib.sha256()
    nbytes = 0
    extracted: list[str] = []
    skipped = 0

    class _HashingReader(io.RawIOBase):
        def __init__(self, raw):
            self._raw = raw
        def readable(self):
            return True
        def readinto(self, b):
            nonlocal nbytes
            data = self._raw.read(len(b))
            if not data:
                return 0
            n = len(data)
            b[:n] = data
            h_md5.update(data)
            h_sha.update(data)
            nbytes += n
            return n

    with urllib.request.urlopen(req, timeout=timeout) as resp:
        hr = io.BufferedReader(_HashingReader(resp), buffer_size=1 << 20)
        with tarfile.open(fileobj=hr, mode="r|gz") as tf:
            for m in tf:
                if not m.isfile():
                    continue
                if wanted_predicate(m.name):
                    base = os.path.basename(m.name)
                    out_p = extract_dir / base
                    ef = tf.extractfile(m)
                    if ef is not None:
                        with open(out_p, "wb") as out:
                            shutil.copyfileobj(ef, out, length=1 << 20)
                        extracted.append(str(out_p))
                else:
                    skipped += 1
        # Drain any trailing gzip padding bytes so full-file md5/sha256 matches
        while hr.read(1 << 20):
            pass

    md5 = h_md5.hexdigest()
    sha = h_sha.hexdigest()
    if expected_bytes is not None and nbytes != expected_bytes:
        raise Refused(f"{url}: streamed size {nbytes} != expected {expected_bytes}")
    if expected_md5 is not None and md5 != expected_md5:
        raise Refused(f"{url}: streamed md5 {md5} != expected {expected_md5}")
    return {
        "url": url,
        "bytes": nbytes,
        "md5": md5,
        "sha256": sha,
        "extracted_count": len(extracted),
        "skipped_count": skipped,
        "extracted_files": extracted,
    }


def cmd_stage(args: argparse.Namespace) -> None:
    manifest = json.loads(Path(args.manifest_json).read_text(encoding="utf-8"))
    mag_entries: list[dict] = manifest.get("mags", [])
    unbinned_entries: list[dict] = manifest.get("unbinned", [])
    source_meta: dict = manifest["source_meta"]

    dropped_skani: list[dict] = []
    if args.skani_dedup and mag_entries:
        mag_entries, dropped_skani = run_skani_self_dedup(
            mag_entries,
            Path(args.out_meta_dir) / "_skani_work",
            min_ani=args.skani_ani,
            min_af=args.skani_af,
            threads=args.skani_threads,
        )

    owner_inputs = {}
    for item in args.owner_inputs or []:
        k, v = item.split("=", 1)
        owner_inputs[k] = Path(v)

    rep = stage_compendium(
        stage=args.stage,
        mag_entries=mag_entries,
        unbinned_entries=unbinned_entries,
        out_meta_dir=Path(args.out_meta_dir),
        out_fasta_dir=Path(args.out_fasta_dir),
        final_fasta_dir=Path(args.final_fasta_dir) if args.final_fasta_dir else None,
        fasta_name=args.fasta_name,
        source_meta=source_meta,
        min_unbinned_len=args.min_unbinned_len,
        root=Path(args.root) if args.root else None,
        root_stages=args.root_stages.split(",") if args.root_stages else (),
        owner_inputs=owner_inputs,
        retrieved_at=args.retrieved_at,
        screen_codes=not args.no_screen_codes,
        screen_max_bp=args.screen_max_bp,
        threads=args.threads,
        memory_limit=args.memory_limit,
    )
    if dropped_skani:
        skani_drop_path = Path(args.out_meta_dir) / "skani_within_dropped.tsv"
        with open(skani_drop_path, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["genome_id", "source_name", "reason"], delimiter="\t")
            w.writeheader()
            w.writerows(dropped_skani)
        rep["skani_within_dropped"] = len(dropped_skani)
        (Path(args.out_meta_dir) / "stage_report.json").write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(rep, indent=2))


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Card C-66: LR2 and NX1 RT-hotspot compendia staging")
    sub = p.add_subparsers(dest="cmd", required=True)

    ps = sub.add_parser("stage")
    ps.add_argument("stage", choices=["LR2", "NX1"])
    ps.add_argument("manifest_json")
    ps.add_argument("out_meta_dir")
    ps.add_argument("out_fasta_dir")
    ps.add_argument("--final-fasta-dir", default=None)
    ps.add_argument("--fasta-name", required=True)
    ps.add_argument("--min-unbinned-len", type=int, default=5_000)
    ps.add_argument("--root", default=None)
    ps.add_argument("--root-stages", default="")
    ps.add_argument("--owner-inputs", action="append", default=[])
    ps.add_argument("--retrieved-at", default="2026-11-24T00:00:00Z")
    ps.add_argument("--no-screen-codes", action="store_true")
    ps.add_argument("--screen-max-bp", type=int, default=500_000)
    ps.add_argument("--skani-dedup", action="store_true")
    ps.add_argument("--skani-ani", type=float, default=96.0)
    ps.add_argument("--skani-af", type=float, default=50.0)
    ps.add_argument("--skani-threads", type=int, default=8)
    ps.add_argument("--threads", type=int, default=2)
    ps.add_argument("--memory-limit", default="4GB")
    ps.set_defaults(func=cmd_stage)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

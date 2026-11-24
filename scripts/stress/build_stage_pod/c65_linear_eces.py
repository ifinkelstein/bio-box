"""Card C-65 (DEC-69, DEC-71, DEC-65, DEC-50, DEC-51, DEC-26, DEC-27):
fetch, deduplicate, and stage giant linear extrachromosomal elements (`ECE1`):
  - archaeal `Borgs` & *Ca. Methanoperedens* linear/circular ECEs (`PRJNA866293`, `PRJNA850006`, NCBI nuccore);
  - *Streptomyces* / *Rhodococcus* / Actinomycetia linear plasmids & megaplasmids (`PRJNA747871`, NCBI linear
    Actinomycetota plasmids, Zenodo `7816421` BGC-encoding prokaryotic plasmids);
  - curated phage-plasmids & non-gut megaplasmids (`EpfeiferNutri/Phage-plasmids` tyPPing tables, Zenodo
    `16616313` complete phage-plasmid assemblies, and NCBI complete bacterial megaplasmids).

Hard rules enforced:
- Every source in `ECE1` must be `redistributable: true` (`DEC-51`).
- Linear replicons and public-database records are preserved **byte-identically as deposited** (`DEC-65`):
  legitimate biological terminal inverted repeats (TIRs) are detected and recorded in `ece_metadata.tsv`,
  NEVER trimmed (`trims.tsv` always has 0 rows).
- All `ECE1` replicons are staged with `kind="plasmid"`. `circular` is set ONLY where explicit public
  metadata (`LOCUS` / `esummary` `topology == "circular"`) marks `circular`, and NEVER for Borgs or
  linear plasmids.
- Within-stage `contig_uid` deduplication (`DEC-27`) and cross-stage build-order ownership (`DEC-26`).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # c48c_virus_inputs and make_inputs beside this script

import c48c_virus_inputs as v2  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

ECE_CLASSES = ("borg", "linear_megaplasmid", "phage_plasmid")
VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
_COMP_TABLE = bytes.maketrans(
    b"ACGTacgtRYKMSWBDHVNrykmswbdhvn",
    b"TGCAtgcaYRMKWSVHDBNyrmkwsvhdbn",
)

RECORD_HEAD = [
    "ordinal", "source", "ece_class", "accession", "contig_uid", "length",
    "circular", "tir_length", "host_taxon", "owner_stage",
]
SIDECAR_HEAD = [
    "stage", "source", "ece_class", "genome_id", "contig_uid", "accession",
    "length", "circular", "tir_length", "host_taxon", "owner_stage",
]
TRIM_HEAD = [
    "source", "genome", "accession", "source_contig_uid", "contig_uid",
    "source_length", "repeat_length", "bases_removed", "length",
]
SKIPPED_HEAD = ["source", "accession", "reason", "length"]

PFEIFER_SHA = "fb10c7fe5ba82da7faed0a2c0387f4e689f043f7"
PRJNA866293_ACCESSIONS = [
    "GCA_027460385.1", "GCA_027460405.1", "GCA_027460425.1", "GCA_027460435.1",
    "GCA_027460465.1", "GCA_027460485.1", "GCA_027460495.1", "GCA_027460525.1",
    "GCA_027460535.1", "GCA_027460565.1", "GCA_027460585.1", "GCA_027460605.1",
    "GCA_027460615.1",
]


def _base_acc(acc: str) -> str:
    acc = acc.strip()
    return acc.split(".", 1)[0] if "." in acc else acc


def detect_tir(seq: bytes, min_len: int = 20, max_len: int = 100_000) -> int:
    """Detect exact terminal inverted repeat (TIR) length at the ends of a linear replicon.

    Compares `seq[i]` against the Watson-Crick complement of `seq[n - 1 - i]` from `i = 0`
    inward up to `min(max_len, len(seq) // 2)`. Returns the exact match length if `>= min_len`,
    otherwise 0. Read-only: NEVER trims or modifies `seq` (`DEC-65`).
    """
    n = len(seq)
    if n < 2 * min_len:
        return 0
    limit = min(max_len, n // 2)
    upper = seq[:limit].upper()
    tail_rc = seq[n - limit :].upper().translate(_COMP_TABLE)[::-1]
    k = 0
    while k < limit and upper[k] == tail_rc[k]:
        k += 1
    return k if k >= min_len else 0


def _extract_host_taxon(header: str) -> str:
    """Extract a clean host/organism taxon string from a FASTA header description."""
    parts = header.split(None, 1)
    if len(parts) < 2:
        return ""
    desc = parts[1].strip()
    for prefix in ("MAG: ", "TPA_asm: ", "Uncultured ", "UNVERIFIED: "):
        if desc.startswith(prefix):
            desc = desc[len(prefix) :].strip()
    # Cut at common plasmid/strain/isolate/complete markers
    m = re.split(
        r"\b(?:strain|isolate|plasmid|chromosome|clone|megaplasmid|linear plasmid|extrachromosomal|whole genome|complete|contig|scaffold)\b|,",
        desc,
        maxsplit=1,
        flags=re.IGNORECASE,
    )
    taxon = m[0].strip().rstrip(":;")
    return taxon[:120]


def classify_ece_record(header: str, default_class: str, topology: str = "linear") -> dict:
    """Classify an ECE1 FASTA record into `ece_class`, `circular`, `accession`, and `host_taxon`.

    Hard rules:
    - `ece_class` is one of `("borg", "linear_megaplasmid", "phage_plasmid")`.
    - `circular` is True ONLY when `topology.lower() == "circular"` AND the record is not a
      Borg or explicitly linear replicon.
    """
    tok = header.split(None, 1)[0] if header else ""
    ece_class = default_class if default_class in ECE_CLASSES else "linear_megaplasmid"
    low_h = header.lower()
    topo = (topology or "linear").strip().lower()

    if ece_class == "borg":
        # Borgs are giant linear ECEs; only explicitly circular Methanoperedens plasmids (e.g. Mp_ECE_* circular)
        # if NOT labeled Borg and explicitly marked circular in LOCUS/esummary could be circular,
        # but per card C-65 & brief: Borgs are always staged linear.
        circular = False
    elif "linear plasmid" in low_h or "linear replicon" in low_h:
        circular = False
    else:
        circular = topo == "circular"

    return {
        "accession": tok,
        "ece_class": ece_class,
        "circular": circular,
        "host_taxon": _extract_host_taxon(header),
    }


def _open_fasta_in(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return open(path, "rb")


class _DetGzipWriter:
    def __init__(self, path: Path):
        self._path = path
        self._raw = None
        self._gz = None

    def __enter__(self):
        self._raw = open(self._path, "wb")
        if self._path.suffix == ".gz":
            self._gz = gzip.GzipFile(filename="", mode="wb", fileobj=self._raw, compresslevel=6, mtime=0)
            return self._gz
        return self._raw

    def __exit__(self, exc_type, exc, tb):
        if self._gz is not None:
            self._gz.close()
        if self._raw is not None:
            self._raw.close()


def stage_ece1(
    sources: list[dict],
    out_dir: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_dir: Path | None = None,
    stage_name: str = "ECE1",
) -> dict:
    """Stage `ECE1` inputs byte-identically (0 trims, `DEC-65`), refusing non-redistributable sources."""
    if not sources:
        raise Refused(f"{stage_name}: no sources to stage")
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir} exists")
    for s in sources:
        if not s.get("redistributable", False):
            raise Refused(
                f"{stage_name}: refusing non-redistributable source {s.get('source_name')!r} (DEC-51)"
            )
        fa = Path(s["fasta"])
        if not fa.is_file():
            raise Refused(f"{stage_name}: FASTA file not found: {fa}")

    work = out_dir.with_name(f".{out_dir.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out_dir

    try:
        rep = _stage_ece1_work(
            sources, work, obj_dir, stage_name, root, list(root_stages), owner_inputs, retrieved
        )
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise

    work.rename(out_dir)
    sum_files = sorted(
        p.name for p in out_dir.iterdir()
        if p.is_file() and p.name not in ("inputs_SUMS", "stage_report.json")
    )
    sha_map = {n: ids.sha256_file(out_dir / n).hex() for n in sum_files}
    bytes_map = {n: (out_dir / n).stat().st_size for n in sum_files}
    with open(out_dir / "inputs_SUMS", "w") as w:
        for n in sum_files:
            w.write(f"{sha_map[n]}  {n}\n")
    rep["sha256"] = sha_map
    rep["bytes"] = bytes_map
    (out_dir / "stage_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True))
    return rep


def _stage_ece1_work(
    sources: list[dict],
    work: Path,
    obj_dir: Path,
    stage_name: str,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
) -> dict:
    # (source_name, ece_class, acc, cu, length, circ, tir_len, host_taxon)
    staged_records: list[tuple[str, str, str, bytes, int, bool, int, str]] = []
    skipped: list[tuple[str, str, str, int]] = []
    trims: list[tuple] = []  # Always empty under DEC-65
    objects_json: list[dict] = []
    by_source: dict[str, dict] = {}

    for src in sources:
        sname = src["source_name"]
        out_name = src["output_name"]
        fa_path = Path(src["fasta"])
        default_class = src.get("ece_class", "linear_megaplasmid")
        topo_map: dict[str, str] = src.get("topology_map") or {}
        base_topo_map = {_base_acc(k): v for k, v in topo_map.items()}
        dst_fa = work / out_name
        rec_count = 0
        base_count = 0
        skip_count = 0

        with _open_fasta_in(fa_path) as fh, _DetGzipWriter(dst_fa) as w:
            for r in v2.iter_fasta_stream(fh):
                if len(r.seq) == 0:
                    skipped.append((sname, r.accession, "empty_sequence", 0))
                    skip_count += 1
                    continue
                if not VALID_DNA_RE.match(r.seq):
                    skipped.append((sname, r.accession, "non_iupac_dna", len(r.seq)))
                    skip_count += 1
                    continue
                topo = topo_map.get(r.accession) or base_topo_map.get(_base_acc(r.accession)) or "linear"
                meta = classify_ece_record(r.header, default_class, topo)
                seq = r.seq  # Byte-identical, NEVER trimmed (DEC-65)
                circ = bool(meta["circular"])
                tir_len = 0 if circ else detect_tir(seq)
                w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
                cu = ids.contig_uid(seq)
                staged_records.append((
                    sname, meta["ece_class"], meta["accession"], cu, len(seq),
                    circ, tir_len, meta["host_taxon"],
                ))
                rec_count += 1
                base_count += len(seq)

        if rec_count == 0:
            dst_fa.unlink(missing_ok=True)
            continue

        objects_json.append({
            "path": str(obj_dir / out_name),
            "source_name": sname,
            "source_release": src["source_release"],
            "license": src["license"],
            "redistributable": True,
            "retrieved_at": retrieved,
        })
        by_source[out_name] = {
            "source_name": sname,
            "ece_class": default_class,
            "records": rec_count,
            "bases": base_count,
            "skipped": skip_count,
        }

    if not staged_records:
        raise Refused(f"{stage_name}: no valid records staged")

    all_uids = [x[3] for x in staged_records]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    sidecar_rows: list[tuple] = []
    repeated_rows: list[tuple] = []
    shared_rows: list[tuple] = []
    first_of: dict[bytes, str] = {}
    class_counts: Counter[str] = Counter()
    genome_class_counts: Counter[str] = Counter()
    tir_count = 0
    circular_genomes = 0

    for ordinal, (sname, ece_cls, acc, cu, length, circ, tir_len, host_tax) in enumerate(staged_records):
        owner = owners.get(cu, "")
        rec_rows.append((
            ordinal, sname, ece_cls, acc, cu.hex(), length,
            str(circ).lower(), tir_len, host_tax, owner,
        ))
        if owner:
            shared_rows.append((ordinal, acc, cu.hex(), owner))
        class_counts[ece_cls] += 1

        g = listing.add([cu], 11, "plasmid", "contig")
        guid = ids.genome_uid([cu])
        if g is None:
            repeated_rows.append((ordinal, acc, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first_of[guid] = acc
        genome_class_counts[ece_cls] += 1
        if circ:
            circular_genomes += 1
        if tir_len > 0:
            tir_count += 1
        sidecar_rows.append((
            stage_name, sname, ece_cls, guid.hex(), cu.hex(), acc,
            length, str(circ).lower(), tir_len, host_tax, owner,
        ))

    sidecar_rows.sort(key=lambda t: (t[2], t[1], t[5], t[4]))
    make_inputs._write_tsv(work / "records.tsv", RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work / "ece_metadata.tsv", SIDECAR_HEAD, sidecar_rows)
    make_inputs._write_tsv(
        work / "repeated_records.tsv",
        ["ordinal", "accession", "genome_uid", "first_accession"],
        repeated_rows,
    )
    make_inputs._write_tsv(
        work / "shared_records.tsv",
        ["ordinal", "accession", "contig_uid", "owner_stage"],
        shared_rows,
    )
    make_inputs._write_tsv(work / "trims.tsv", TRIM_HEAD, trims)
    make_inputs._write_tsv(work / "skipped_records.tsv", SKIPPED_HEAD, skipped)

    (work / "objects.json").write_text(json.dumps(objects_json, indent=1))
    (work / "genomes.json").write_text(json.dumps(genomes))

    contig_len: dict[bytes, int] = {}
    contig_cls: dict[bytes, str] = {}
    for _, ece_cls, _, cu, length, *_ in staged_records:
        contig_len.setdefault(cu, length)
        contig_cls.setdefault(cu, ece_cls)

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    bases_by_class: Counter[str] = Counter()
    for cu, length in contig_len.items():
        bases_by_class[contig_cls[cu]] += length
        if cu in owners:
            st = owners[cu]
            owner_contigs[st] += 1
            owner_bases[st] += length

    total_contigs = len(contig_len)
    total_bases = sum(contig_len.values())
    held_contigs = sum(owner_contigs.values())
    held_bases = sum(owner_bases.values())

    return {
        "stage": stage_name,
        "redistributable": True,
        "records": len(staged_records),
        "skipped_records": len(skipped),
        "genomes": len(genomes),
        "contigs": total_contigs,
        "bases": total_bases,
        "circular_genomes": circular_genomes,
        "linear_genomes": len(genomes) - circular_genomes,
        "genomes_with_tir": tir_count,
        "repeated_records": len(repeated_rows),
        "trims": len(trims),
        "owned_by_stage": dict(sorted(owner_contigs.items())),
        "owned_bases_by_stage": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": total_contigs - held_contigs,
        "new_bases": total_bases - held_bases,
        "records_by_class": dict(sorted(class_counts.items())),
        "genomes_by_class": dict(sorted(genome_class_counts.items())),
        "bases_by_class": dict(sorted(bases_by_class.items())),
        "by_source": by_source,
        "root_stages": root_stages,
        "owner_inputs": list(owner_inputs),
    }


# --------------------------------------------------------------------------------------------------------------
# Network fetch helpers (DEC-50: <= 2 req/s, User-Agent: curl, stop after 5 consecutive refusals)


class _RateLimiter:
    def __init__(self, min_interval_s: float = 0.55):
        self.min_interval = min_interval_s
        self.last = 0.0
        self.refusals = 0

    def get(self, url: str, timeout: int = 120) -> bytes:
        wait = self.min_interval - (time.monotonic() - self.last)
        if wait > 0:
            time.sleep(wait)
        self.last = time.monotonic()
        req = urllib.request.Request(url, headers={"User-Agent": "curl"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
                self.refusals = 0
                return data
        except urllib.error.HTTPError as e:
            if e.code in (403, 429, 500, 502, 503, 504):
                self.refusals += 1
                if self.refusals >= 5:
                    raise Refused(f"5 consecutive refusals; last HTTP {e.code} on {url}") from e
            raise

    def download_file(self, url: str, dst: Path, expected_md5: str | None = None, timeout: int = 600) -> tuple[int, str, str]:
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp = dst.with_name(dst.name + ".part")
        wait = self.min_interval - (time.monotonic() - self.last)
        if wait > 0:
            time.sleep(wait)
        self.last = time.monotonic()
        req = urllib.request.Request(url, headers={"User-Agent": "curl"})
        h_sha = hashlib.sha256()
        h_md5 = hashlib.md5()
        nbytes = 0
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as w:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    w.write(chunk)
                    h_sha.update(chunk)
                    h_md5.update(chunk)
                    nbytes += len(chunk)
            self.refusals = 0
        except urllib.error.HTTPError as e:
            tmp.unlink(missing_ok=True)
            if e.code in (403, 429, 500, 502, 503, 504):
                self.refusals += 1
                if self.refusals >= 5:
                    raise Refused(f"5 consecutive refusals; last HTTP {e.code} on {url}") from e
            raise
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
        got_md5 = h_md5.hexdigest()
        if expected_md5 and got_md5 != expected_md5:
            tmp.unlink(missing_ok=True)
            raise Refused(f"{dst.name}: md5 mismatch (got {got_md5}, want {expected_md5})")
        tmp.rename(dst)
        return nbytes, h_sha.hexdigest(), got_md5


def _esearch_all(rl: _RateLimiter, term: str) -> list[str]:
    ids_out: list[str] = []
    retstart = 0
    retmax = 10000
    while True:
        u = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=nuccore&term={urllib.parse.quote(term)}&retmax={retmax}&retstart={retstart}&retmode=json"
        )
        d = json.loads(rl.get(u, timeout=60))["esearchresult"]
        batch = d.get("idlist") or []
        ids_out.extend(batch)
        total = int(d.get("count") or 0)
        retstart += len(batch)
        if not batch or retstart >= total:
            break
    return ids_out


def _esummary_map(rl: _RateLimiter, id_list: list[str], batch_size: int = 200) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for i in range(0, len(id_list), batch_size):
        chunk = id_list[i : i + batch_size]
        u = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=nuccore&id={','.join(chunk)}&retmode=json"
        )
        for attempt in range(3):
            try:
                res = json.loads(rl.get(u, timeout=60)).get("result") or {}
                for k in chunk:
                    if k in res and isinstance(res[k], dict):
                        out[k] = res[k]
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2.0)
    return out


def _efetch_fasta(rl: _RateLimiter, id_list: list[str], out_path: Path, batch_size: int = 80) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".part")
    total_bytes = 0
    with open(tmp, "wb") as w:
        for i in range(0, len(id_list), batch_size):
            chunk = id_list[i : i + batch_size]
            u = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                f"?db=nuccore&id={','.join(chunk)}&rettype=fasta&retmode=text"
            )
            for attempt in range(3):
                try:
                    data = rl.get(u, timeout=180)
                    w.write(data)
                    if data and not data.endswith(b"\n"):
                        w.write(b"\n")
                    total_bytes += len(data)
                    break
                except Exception:
                    if attempt == 2:
                        tmp.unlink(missing_ok=True)
                        raise
                    time.sleep(3.0)
    tmp.rename(out_path)
    return total_bytes


def fetch_all(dl_dir: Path) -> dict:
    """Fetch all `sources.yaml`-listed ECE1 datasets into `dl_dir`."""
    dl_dir = Path(dl_dir)
    dl_dir.mkdir(parents=True, exist_ok=True)
    rl = _RateLimiter(0.55)
    fetched: dict[str, dict] = {}

    # 1. NCBI Datasets PRJNA866293 (13 Ca. Methanoperedens / Borg assemblies)
    ds_zip = dl_dir / "prjna866293_datasets.zip"
    if not ds_zip.is_file():
        url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            + ",".join(PRJNA866293_ACCESSIONS)
            + "/download?include_annotation_type=GENOME_FASTA"
        )
        nb, sha, md5 = rl.download_file(url, ds_zip, timeout=300)
        fetched[ds_zip.name] = {"bytes": nb, "sha256": sha, "md5": md5}
    else:
        fetched[ds_zip.name] = {"bytes": ds_zip.stat().st_size, "sha256": ids.sha256_file(ds_zip).hex()}

    # 2. NCBI nuccore archaeal Borgs & Methanoperedens ECEs
    arch_fa = dl_dir / "ncbi_archaeal_borgs_eces.fna"
    arch_topo_p = dl_dir / "ncbi_archaeal_borgs_eces.topology.json"
    if not arch_fa.is_file():
        id_set: list[str] = []
        seen: set[str] = set()
        for term in (
            "PRJNA850006[BioProject] AND plasmid[Title]",
            "Methanoperedens[Organism] AND plasmid[Title]",
            "extrachromosomal[Title] AND archaea[Organism]",
        ):
            for x in _esearch_all(rl, term):
                if x not in seen:
                    seen.add(x)
                    id_set.append(x)
        smap = _esummary_map(rl, id_set)
        topo = {v.get("accessionversion", k): v.get("topology", "linear") for k, v in smap.items()}
        arch_topo_p.write_text(json.dumps(topo, indent=1, sort_keys=True))
        _efetch_fasta(rl, id_set, arch_fa, batch_size=50)
    fetched[arch_fa.name] = {"bytes": arch_fa.stat().st_size, "sha256": ids.sha256_file(arch_fa).hex()}

    # 3. NCBI PRJNA747871 Actinomycetia plasmids (699 plasmids)
    p747_fa = dl_dir / "prjna747871_actinomycetia_plasmids.fna"
    p747_topo_p = dl_dir / "prjna747871_actinomycetia_plasmids.topology.json"
    p747_accs: set[str] = set()
    if not p747_fa.is_file():
        p747_ids = _esearch_all(rl, "PRJNA747871[BioProject] AND plasmid[Title]")
        smap = _esummary_map(rl, p747_ids)
        topo = {v.get("accessionversion", k): v.get("topology", "linear") for k, v in smap.items()}
        p747_topo_p.write_text(json.dumps(topo, indent=1, sort_keys=True))
        p747_accs = {_base_acc(a) for a in topo}
        _efetch_fasta(rl, p747_ids, p747_fa, batch_size=60)
    else:
        p747_accs = {_base_acc(a) for a in json.loads(p747_topo_p.read_text())}
    fetched[p747_fa.name] = {"bytes": p747_fa.stat().st_size, "sha256": ids.sha256_file(p747_fa).hex()}

    # 4. NCBI Actinomycetota linear plasmids (>= 10 kb) + bacterial linear plasmids (>= 5 kb)
    lin_fa = dl_dir / "ncbi_actinomycetota_linear_plasmids.fna"
    lin_topo_p = dl_dir / "ncbi_actinomycetota_linear_plasmids.topology.json"
    if not lin_fa.is_file():
        act_ids = _esearch_all(
            rl,
            'Actinomycetota[Organism] AND plasmid[Title] AND "complete sequence"[Title] AND 10000:5000000[SLEN]',
        )
        smap = _esummary_map(rl, act_ids)
        keep_ids: list[str] = []
        topo: dict[str, str] = {}
        seen_base = set(p747_accs)
        for uid, rec in smap.items():
            if (rec.get("topology") or "").lower() != "linear":
                continue
            acc = rec.get("accessionversion") or uid
            bacc = _base_acc(acc)
            if bacc in seen_base:
                continue
            seen_base.add(bacc)
            keep_ids.append(uid)
            topo[acc] = "linear"
        # Also include explicit bacterial "linear plasmid"[Title] >= 5 kb
        bact_lin_ids = _esearch_all(
            rl,
            'Bacteria[Organism] AND "linear plasmid"[Title] AND 5000:5000000[SLEN]',
        )
        smap2 = _esummary_map(rl, bact_lin_ids)
        for uid, rec in smap2.items():
            acc = rec.get("accessionversion") or uid
            bacc = _base_acc(acc)
            if bacc in seen_base:
                continue
            seen_base.add(bacc)
            keep_ids.append(uid)
            topo[acc] = "linear"
        lin_topo_p.write_text(json.dumps(topo, indent=1, sort_keys=True))
        _efetch_fasta(rl, keep_ids, lin_fa, batch_size=60)
    fetched[lin_fa.name] = {"bytes": lin_fa.stat().st_size, "sha256": ids.sha256_file(lin_fa).hex()}

    # 5. Zenodo 7816421 prokaryotic BGC plasmids (283 MB zip, md5 e129115952c84b8a368857f954441791)
    z781 = dl_dir / "plasmid_fasta_final_filtered.zip"
    if not z781.is_file():
        nb, sha, md5 = rl.download_file(
            "https://zenodo.org/api/records/7816421/files/plasmid_fasta_final_filtered.zip/content",
            z781,
            expected_md5="e129115952c84b8a368857f954441791",
            timeout=600,
        )
        fetched[z781.name] = {"bytes": nb, "sha256": sha, "md5": md5}
    else:
        fetched[z781.name] = {"bytes": z781.stat().st_size, "sha256": ids.sha256_file(z781).hex()}

    # 6. Pfeifer/Rocha GitHub tyPPing tables + NCBI nuccore fetch for phage-plasmids
    pp_accs: list[str] = []
    seen_pp: set[str] = set()
    for fname in ("Final_prediction_table_0523_by_tyPPing.tsv", "Final_prediction_table_0321_by_tyPPing.tsv"):
        dst = dl_dir / fname
        if not dst.is_file():
            url = (
                f"https://raw.githubusercontent.com/EpfeiferNutri/Phage-plasmids/{PFEIFER_SHA}"
                f"/Publication_related_data/Supplementary_data/{fname}"
            )
            nb, sha, md5 = rl.download_file(url, dst, timeout=60)
            fetched[dst.name] = {"bytes": nb, "sha256": sha, "md5": md5}
        else:
            fetched[dst.name] = {"bytes": dst.stat().st_size, "sha256": ids.sha256_file(dst).hex()}
        with open(dst, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gid = (row.get("Genome ID") or "").strip()
                if gid and gid not in seen_pp:
                    seen_pp.add(gid)
                    pp_accs.append(gid)

    pp_fa = dl_dir / "pfeifer_phage_plasmids.fna"
    pp_topo_p = dl_dir / "pfeifer_phage_plasmids.topology.json"
    if not pp_fa.is_file():
        # Also include explicit NCBI phage-plasmid title matches
        extra_ids = _esearch_all(
            rl,
            '"phage-plasmid"[Title] OR "phage plasmid"[Title] OR "prophage-plasmid"[Title]',
        )
        smap_extra = _esummary_map(rl, extra_ids)
        topo: dict[str, str] = {}
        for rec in smap_extra.values():
            acc = rec.get("accessionversion") or ""
            if acc:
                topo[acc] = rec.get("topology", "linear")
                if _base_acc(acc) not in seen_pp and acc not in seen_pp:
                    seen_pp.add(acc)
                    pp_accs.append(acc)
        # Resolve esummary topology for pp_accs in batches of 40 (safe GET URL length)
        for i in range(0, len(pp_accs), 40):
            chunk = pp_accs[i : i + 40]
            q = " OR ".join(f"{a}[Accession]" for a in chunk)
            uids = _esearch_all(rl, q)
            sm = _esummary_map(rl, uids)
            for rec in sm.values():
                acc = rec.get("accessionversion") or ""
                if acc:
                    topo[acc] = rec.get("topology", "linear")
        pp_topo_p.write_text(json.dumps(topo, indent=1, sort_keys=True))
        _efetch_fasta(rl, pp_accs, pp_fa, batch_size=80)
    fetched[pp_fa.name] = {"bytes": pp_fa.stat().st_size, "sha256": ids.sha256_file(pp_fa).hex()}

    # 7. Zenodo 16616313 PPs_complete_assemblies.zip (294 KB)
    z166 = dl_dir / "PPs_complete_assemblies.zip"
    if not z166.is_file():
        nb, sha, md5 = rl.download_file(
            "https://zenodo.org/api/records/16616313/files/PPs_complete_assemblies.zip/content",
            z166,
            expected_md5="0bba2fa3fcd214c36b9e7ba02c2cc0e5",
            timeout=120,
        )
        fetched[z166.name] = {"bytes": nb, "sha256": sha, "md5": md5}
    else:
        fetched[z166.name] = {"bytes": z166.stat().st_size, "sha256": ids.sha256_file(z166).hex()}

    # 8. NCBI complete bacterial megaplasmids (>= 100 kb)
    mp_fa = dl_dir / "ncbi_bacterial_megaplasmids.fna"
    mp_topo_p = dl_dir / "ncbi_bacterial_megaplasmids.topology.json"
    if not mp_fa.is_file():
        mp_ids = _esearch_all(
            rl,
            'Bacteria[Organism] AND megaplasmid[Title] AND "complete sequence"[Title] AND 100000:5000000[SLEN]',
        )
        smap = _esummary_map(rl, mp_ids)
        topo = {v.get("accessionversion", k): v.get("topology", "linear") for k, v in smap.items()}
        mp_topo_p.write_text(json.dumps(topo, indent=1, sort_keys=True))
        _efetch_fasta(rl, mp_ids, mp_fa, batch_size=40)
    fetched[mp_fa.name] = {"bytes": mp_fa.stat().st_size, "sha256": ids.sha256_file(mp_fa).hex()}

    with open(dl_dir / "fetch_SUMS", "w") as w:
        for name in sorted(fetched):
            w.write(f"{fetched[name]['sha256']}  {name}\n")
    (dl_dir / "fetch_report.json").write_text(json.dumps(fetched, indent=1, sort_keys=True))
    return fetched


def _prepare_sources(dl_dir: Path, scratch: Path) -> list[dict]:
    scratch.mkdir(parents=True, exist_ok=True)
    sources: list[dict] = []

    # 1. PRJNA866293 datasets zip -> prjna866293_borgs.fna
    ds_zip = dl_dir / "prjna866293_datasets.zip"
    if ds_zip.is_file():
        out_fa = scratch / "prjna866293_borgs.fna"
        with zipfile.ZipFile(ds_zip) as zf, open(out_fa, "wb") as w:
            for m in sorted(zf.namelist()):
                if m.endswith(".fna"):
                    with zf.open(m) as fh:
                        shutil.copyfileobj(fh, w)
        sources.append({
            "source_name": "ncbi_datasets_c65_borg_hosts",
            "source_release": "PRJNA866293 (Al-Shayeb et al. 2022, Nature 610:310-316)",
            "license": "INSDC (public)",
            "redistributable": True,
            "fasta": out_fa,
            "output_name": "prjna866293_borgs.fna.gz",
            "ece_class": "borg",
        })

    # 2. NCBI archaeal Borgs & Methanoperedens ECEs
    arch_fa = dl_dir / "ncbi_archaeal_borgs_eces.fna"
    if arch_fa.is_file():
        topo = json.loads((dl_dir / "ncbi_archaeal_borgs_eces.topology.json").read_text())
        sources.append({
            "source_name": "ncbi_eutils_c65_ece1",
            "source_release": "2026-11-24 (PRJNA850006 + Methanoperedens/archaeal ECEs)",
            "license": "INSDC (public)",
            "redistributable": True,
            "fasta": arch_fa,
            "output_name": "ncbi_archaeal_borgs_eces.fna.gz",
            "ece_class": "borg",
            "topology_map": topo,
        })

    # 3. NCBI PRJNA747871 Actinomycetia plasmids
    p747_fa = dl_dir / "prjna747871_actinomycetia_plasmids.fna"
    if p747_fa.is_file():
        topo = json.loads((dl_dir / "prjna747871_actinomycetia_plasmids.topology.json").read_text())
        sources.append({
            "source_name": "ncbi_eutils_c65_ece1",
            "source_release": "2026-11-24 (PRJNA747871 Actinomycetia non-chromosomal plasmids)",
            "license": "INSDC (public)",
            "redistributable": True,
            "fasta": p747_fa,
            "output_name": "prjna747871_actinomycetia_plasmids.fna.gz",
            "ece_class": "linear_megaplasmid",
            "topology_map": topo,
        })

    # 4. NCBI Actinomycetota & bacterial linear plasmids
    lin_fa = dl_dir / "ncbi_actinomycetota_linear_plasmids.fna"
    if lin_fa.is_file():
        topo = json.loads((dl_dir / "ncbi_actinomycetota_linear_plasmids.topology.json").read_text())
        sources.append({
            "source_name": "ncbi_eutils_c65_ece1",
            "source_release": "2026-11-24 (Actinomycetota & bacterial complete linear plasmids)",
            "license": "INSDC (public)",
            "redistributable": True,
            "fasta": lin_fa,
            "output_name": "ncbi_actinomycetota_linear_plasmids.fna.gz",
            "ece_class": "linear_megaplasmid",
            "topology_map": topo,
        })

    # 5. Zenodo 7816421 BGC plasmids zip -> zenodo_7816421_bgc_plasmids.fna
    z781 = dl_dir / "plasmid_fasta_final_filtered.zip"
    if z781.is_file():
        out_fa = scratch / "zenodo_7816421_bgc_plasmids.fna"
        with zipfile.ZipFile(z781) as zf, open(out_fa, "wb") as w:
            for m in sorted(zf.namelist()):
                if m.endswith((".fa", ".fna", ".fasta")) and not m.startswith("__MACOSX"):
                    with zf.open(m) as fh:
                        data = fh.read()
                        if data:
                            w.write(data)
                            if not data.endswith(b"\n"):
                                w.write(b"\n")
        if out_fa.stat().st_size > 0:
            sources.append({
                "source_name": "zenodo_7816421_bgc_plasmids",
                "source_release": "10.5281/zenodo.7816421 (plasmid_fasta_final_filtered.zip)",
                "license": "CC BY 4.0",
                "redistributable": True,
                "fasta": out_fa,
                "output_name": "zenodo_7816421_bgc_plasmids.fna.gz",
                "ece_class": "linear_megaplasmid",
            })

    # 6. Pfeifer/Rocha phage-plasmids (NCBI fetched from tyPPing tables)
    pp_fa = dl_dir / "pfeifer_phage_plasmids.fna"
    if pp_fa.is_file():
        topo = json.loads((dl_dir / "pfeifer_phage_plasmids.topology.json").read_text())
        sources.append({
            "source_name": "pfeifer_phage_plasmids_github",
            "source_release": f"commit {PFEIFER_SHA} (tyPPing 0523+0321 + NCBI phage-plasmids)",
            "license": "MIT / INSDC (public)",
            "redistributable": True,
            "fasta": pp_fa,
            "output_name": "pfeifer_phage_plasmids.fna.gz",
            "ece_class": "phage_plasmid",
            "topology_map": topo,
        })

    # 7. Zenodo 16616313 complete phage-plasmid assemblies
    z166 = dl_dir / "PPs_complete_assemblies.zip"
    if z166.is_file():
        out_fa = scratch / "pfeifer_complete_assemblies.fna"
        with zipfile.ZipFile(z166) as zf, open(out_fa, "wb") as w:
            for m in sorted(zf.namelist()):
                if m.endswith((".fa", ".fna", ".fasta")) and not m.startswith("__MACOSX"):
                    with zf.open(m) as fh:
                        data = fh.read()
                        if data:
                            w.write(data)
                            if not data.endswith(b"\n"):
                                w.write(b"\n")
        if out_fa.stat().st_size > 0:
            sources.append({
                "source_name": "pfeifer_phage_plasmids_zenodo",
                "source_release": "10.5281/zenodo.16616313 (PPs_complete_assemblies.zip)",
                "license": "CC BY 4.0",
                "redistributable": True,
                "fasta": out_fa,
                "output_name": "pfeifer_complete_assemblies.fna.gz",
                "ece_class": "phage_plasmid",
            })

    # 8. NCBI complete bacterial megaplasmids (>= 100 kb)
    mp_fa = dl_dir / "ncbi_bacterial_megaplasmids.fna"
    if mp_fa.is_file():
        topo = json.loads((dl_dir / "ncbi_bacterial_megaplasmids.topology.json").read_text())
        sources.append({
            "source_name": "ncbi_eutils_c65_ece1",
            "source_release": "2026-11-24 (complete bacterial megaplasmids >= 100 kb)",
            "license": "INSDC (public)",
            "redistributable": True,
            "fasta": mp_fa,
            "output_name": "ncbi_bacterial_megaplasmids.fna.gz",
            "ece_class": "linear_megaplasmid",
            "topology_map": topo,
        })

    return sources


def _parse_kv(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        k, v = item.split("=", 1)
        out[k] = Path(v)
    return out


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Card C-65: fetch and stage ECE1 giant linear ECEs.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    fp = sub.add_parser("fetch")
    fp.add_argument("--dl", type=Path, required=True)

    sp = sub.add_parser("stage")
    sp.add_argument("--dl", type=Path, required=True)
    sp.add_argument("--out", type=Path, required=True)
    sp.add_argument("--retrieved", required=True)
    sp.add_argument("--root", type=Path, default=None)
    sp.add_argument("--root-stages", default="S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5,V1,G1,SPX1")
    sp.add_argument("--owner-inputs", action="append", default=[])
    sp.add_argument("--final-dir", type=Path, default=None)
    sp.add_argument("--out-report", type=Path, default=None)

    args = ap.parse_args(argv)
    if args.cmd == "fetch":
        rep = fetch_all(args.dl)
        print(json.dumps(rep, indent=1, sort_keys=True))
        return 0

    root_stages = [s.strip() for s in args.root_stages.split(",") if s.strip()]
    owner_inputs = _parse_kv(args.owner_inputs)
    scratch = args.dl / ".extract_scratch"
    shutil.rmtree(scratch, ignore_errors=True)
    try:
        sources = _prepare_sources(args.dl, scratch)
        rep = stage_ece1(
            sources,
            args.out,
            root=args.root,
            root_stages=root_stages,
            owner_inputs=owner_inputs,
            retrieved=args.retrieved,
            final_dir=args.final_dir,
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    if args.out_report is not None:
        args.out_report.write_text(json.dumps(rep, indent=1, sort_keys=True))
    print(json.dumps(rep, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

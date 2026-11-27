"""Card C-81 (DEC-79): stage the Groundwater Virome Catalogue (`GWV1`), the OMZ single-cell compendium (`SAG3`) and
the GWM wetland MAGs (`WET1`).

    python c81_gwv1_sag3_wet1.py wet1-stage ZIP XLSX INDEX_DIR GTDB_META_DIR META_OUT FASTA_OUT [as sag3-stage]
    python c81_gwv1_sag3_wet1.py wet1-ani META_DIR FASTA OUT [--whole NAME=META:FASTA[,FASTA]] [--by-genus ...]
        [--taxonomy non_gtdb_stages.tsv] [--archaea FASTA:LIST] [--threads 8]

`WET1` (zenodo 20263536): `Dereplicated genomes.zip` (on the pod as `Dereplicated_genomes.zip`, card R-18's relay
mapping) holds `dereplicated_genomes/<MAG id>.fa`; `MAG_meta.xlsx` describes every MAG (bin name, sample, GTDB-Tk
classification) and gives no accession, so the accession test takes GTDB R232's rows of BioProject `PRJNA1430164` by
bin name (`ncbi_isolate`, `ncbi_strain_identifiers`); every MAG is staged as card C-77 staged single cells
(`assembly_quality="mag"`, every contig name given its MAG id and `|`); `wet1-ani` records ANI >= 96 / AF >= 50 against
staged MAG sets on the pod and card C-75's r226 archaea of the same genera, one skani process per group; nothing is cut.

    python c81_gwv1_sag3_wet1.py gwv1-stage FNA_GZ INFO_TXT META_OUT FASTA_OUT [--root ROOT] [--root-stages ...]
        [--chain-home /tmp/if2453/s2/inputs] [--later NAME=DIR ...] [--inputs-home ...] [--work DIR] [--report-out F]
    python c81_gwv1_sag3_wet1.py gwv1-ani META_DIR FASTA OUT --ref NAME=FASTA ... [--every 25] [--threads 8]
    python c81_gwv1_sag3_wet1.py sag3-census TAR
    python c81_gwv1_sag3_wet1.py sag3-stage TAR XLSX INDEX_DIR GTDB_META_DIR META_OUT FASTA_OUT [--root ROOT]
        [--chain-home ...] [--later NAME=DIR ...] [--inputs-home ...] [--work DIR] [--workers N] [--report-out F]
    python c81_gwv1_sag3_wet1.py sag3-ani META_DIR FASTA OUT --ref NAME=META_DIR:FASTA ... [--threads 8]

`GWV1` (zenodo 11230969): `GWVC.vOTUs.fna.gz`, one record per vOTU representative, described row by row by
`GWVC.vOTUs.information.txt` (CheckV's quality and completeness method, geNomad's taxonomy, the well and zone). Every
record is staged, in the file's order, `kind="virus"` (called by pyrodigal-gv), `assembly_quality="vOTU"`. DEC-65 (1)
(a metagenome-assembled set): a record whose CheckV completeness method is a direct terminal repeat (`DTR (...)`) is
classified by `c62_metavr_inputs.apply_dec65`: an exact 20 to 1,000 bp repeat has its second copy trimmed and is
staged circular, `none` circular untrimmed, anything else linear untrimmed; every other record (ITR, provirus, AAI,
HMM) linear and untrimmed. A record the table does not describe, a table row with no record, or a length that differs
from the table's is refused. DEC-27: a repeated contig set is listed once (`repeated_records.tsv`).
`shared_records.tsv` holds the contigs that the root's stored stages or the `r26.11.4` stages hold (the template's
`EXPECT=held_by` lines), `later_records.tsv` those that only a stage of the `r26.11.5` chain before this one holds
(`--later`, in chain order). `gwv1-ani` is card C-48c's vOTU threshold on a sample (every EVERY-th single-record,
non-repeated genome in FASTA order) against whole held stages, one skani run per held stage, recorded and not cut.

`SAG3` (figshare 20459526, 20481603): `OMZ_Compendium_SAGs_Full.tar.gz` holds one FASTA per single cell
(`Final_Catalog/<name>[_contigs][_output_clean...].fna`); Table S1 (an xlsx, read with the standard library) names
the cells in its `Outputs-Genome_Sequencing` sheet (with their BioSample and GenBank (WGS) accessions) and in its
`Outputs_GTDB-TK` sheet (GTDB-Tk's classification, kept only for a name that sheet gives once: it cuts some pairs of
names, `<cell>_AB-905` and `<cell>_AB-908`, to `<cell>` twice). A member is matched to the longest Table S1 name equal to its stem or followed in it by
`_`; a name two members would take matches neither (labelled so). The accession test: GTDB R232's own metadata
(`ncbi_biosample`, `ncbi_wgs_master`) maps a cell's WGS project or BioSample to a GTDB R232 genome, whose accession
card C-74's `lookup` then places in GTDB R232 and GlobDB r232. No cell is removed for being listed there (the
orchestrating session, 2026-11-26): every cell's status is written to `r232_status.tsv`. GlobDB r232 is seen only
through GTDB R232's accessions here (Table S1 gives no GCA_ accession). Staging is card C-77's: whole contigs, the
genetic code screened per cell (`c58.screen_genetic_code`, 11, 25 or 4; GTDB-Tk's classification as the taxonomic
hint), a contig name that does not begin with the cell's name and `_` given the cell's name and `|`, `kind="prokaryote"`,
`assembly_quality="contig"`. One BGZF FASTA in the archive's member order; metadata in name order. `sag3-ani` records
skani ANI >= 96 over an aligned fraction >= 50 (the larger of the two, as card C-66's dedup) against `SAG1` and `SAG2`;
nothing is cut by it.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
import zipfile
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c39_ani  # noqa: E402
import c58_sags_symbionts_niche as c58  # noqa: E402
import c62_metavr_inputs as c62  # noqa: E402
import c77_sag2_inputs as c77  # noqa: E402
import make_inputs  # noqa: E402
import r232_index  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

ROOT_STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
# the r26.11.5 chain of the plan's §5, in order; a stage is read as a later holder only when its genomes.json is given
R26_11_5_CHAIN = ("GRARE1", "GRARE2", "GRARE3", "LR2", "SAG1", "NX1", "GV1", "EGV1", "ECE2", "MIR1", "ASG1", "SAG2",
                  "DPN1", "SYM2", "LRP1", "MIR2", "MGE4", "MV2", "MVEUKV", "KAVHIFI", "INPH2604", "GWV1", "SAG3",
                  "WET1")


def _write_tsv(path: Path, head: list[str], rows) -> None:
    make_inputs._write_tsv(path, head, rows)


def _tally(contig_len: dict[bytes, int], holders: dict[bytes, str], skip: dict[bytes, str] | None = None):
    c, b = Counter(), Counter()
    for cu, n in contig_len.items():
        if cu in holders and not (skip and cu in skip):
            c[holders[cu]] += 1
            b[holders[cu]] += n
    return dict(sorted(c.items())), dict(sorted(b.items()))


def _check_later(later_inputs: dict[str, Path]) -> None:
    for name, d in later_inputs.items():
        if not (Path(d) / "genomes.json").is_file():
            raise Refused(f"--later {name}: no genomes.json under {d}")


def _owners(records: Path, root, root_stages, owner_inputs, later_inputs, work, threads, memory_limit):
    owners = c58._owners_from_records(records, root, list(root_stages), owner_inputs, work=work, threads=threads,
                                      memory_limit=memory_limit)
    later = c58._owners_from_records(records, None, [], later_inputs, work=work, threads=threads,
                                     memory_limit=memory_limit) if later_inputs else {}
    return owners, later


# === GWV1 ==========================================================================================================

GWV1 = "GWV1"
GWV1_FASTA = "gwv1_gwvc.fna.gz"
GWV1_SOURCE = {
    "source_name": "tower_gwv1",
    "source_release": "Groundwater Virome Catalogue (GWVC), zenodo 11230969 (published 2024-05-21), "
                      "GWVC.vOTUs.fna.gz; DEC-65 terminal repeats by c81_gwv1_sag3_wet1.py",
    "license": "cc-by-4.0 as the Zenodo record states it (read 2026-11-26, card C-72 part 1)",
    "redistributable": True,  # DEC-80 (3): a recorded licence (CC BY 4.0) permits it
}
GWV1_RETRIEVED = "2026-11-26T18:37:00Z"  # tower's SOURCE_URLS.tsv, GWVC.vOTUs.fna.gz
INFO_NEED = ("contig_id", "contig_length", "checkv_quality", "miuvig_quality", "completeness", "completeness_method",
             "taxonomy", "lifestyle", "sample", "zone", "biosample_accession")
GWV1_RECORDS_HEAD = ["ordinal", "genome", "contig_uid", "length", "checkv_quality", "miuvig_quality", "completeness",
                     "completeness_method", "taxonomy", "lifestyle", "sample", "zone", "biosample", "circular",
                     "dtr_class", "source_contig_uid", "source_length", "trimmed_bases"]
GWV1_META = ("genomes.json", "objects.json", "objects_staging.json", "records.tsv", "repeated_records.tsv",
             "shared_records.tsv", "later_records.tsv", "trims.tsv")


def read_gwvc_info(path: Path) -> dict[str, tuple]:
    """contig_id -> (length, the other INFO_NEED fields in order); refuses a missing column or a repeated id."""
    out: dict[str, tuple] = {}
    with open(path, encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        head = next(rd, None) or []
        miss = [c for c in INFO_NEED if c not in head]
        if miss:
            raise Refused(f"{path}: no column {miss[0]}")
        ix = [head.index(c) for c in INFO_NEED]
        for row in rd:
            if not row:
                continue
            f = [row[i] if i < len(row) else "" for i in ix]
            if f[0] in out:
                raise Refused(f"{path}: contig {f[0]} has two rows")
            try:
                n = int(f[1])
            except ValueError:
                raise Refused(f"{path}: contig {f[0]}: contig_length {f[1]!r} is not a number") from None
            out[f[0]] = (n, *f[2:])
    return out


def is_dtr(method: str) -> bool:
    """CheckV's completeness method names a direct terminal repeat (`DTR (high-confidence)`, `DTR (medium-...)`)."""
    return method.startswith("DTR")


def stage_gwv1(fna: Path, info_txt: Path, meta_out: Path, fasta_out: Path, *, root: Path | None, root_stages,
               owner_inputs: dict[str, Path], later_inputs: dict[str, Path], inputs_home: Path,
               retrieved_at: str = GWV1_RETRIEVED, work: Path | None = None, threads: int = 2,
               memory_limit: str = "4GB") -> dict:
    """Stage GWV1; returns the report (module docstring). Nothing is left behind when it refuses."""
    meta_dir, fa_dir = Path(meta_out) / GWV1, Path(fasta_out) / GWV1
    final_fa = fa_dir / GWV1_FASTA
    for p in (meta_dir, final_fa):
        if p.exists():
            raise Refused(f"{p} exists")
    _check_later(later_inputs)
    info = read_gwvc_info(Path(info_txt))
    tmp = Path(meta_out) / f"{GWV1}.tmp"
    partial = fa_dir / f".{GWV1_FASTA}.partial"
    work_dir = Path(work) if work is not None else Path(meta_out) / f"{GWV1}.work"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    fa_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        rep = _stage_gwv1(Path(fna), info, tmp, partial, work_dir, root, list(root_stages), owner_inputs,
                          later_inputs, Path(inputs_home), final_fa, retrieved_at, threads, memory_limit)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        partial.unlink(missing_ok=True)
        if fa_dir.exists() and not any(fa_dir.iterdir()):
            fa_dir.rmdir()
        raise
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    partial.rename(final_fa)
    tmp.rename(meta_dir)
    files = {n: meta_dir / n for n in GWV1_META}
    files[GWV1_FASTA] = final_fa
    rep["sha256"] = {n: ids.sha256_file(p).hex() for n, p in files.items()}
    rep["bytes"] = {n: p.stat().st_size for n, p in files.items()}
    return rep


def _stage_gwv1(fna, info, tmp, partial, work_dir, root, root_stages, owner_inputs, later_inputs, inputs_home,
                final_fa, retrieved_at, threads, memory_limit) -> dict:
    recs: list[tuple] = []  # (name, cu, n, src_n, src_cu, circular, dtr_class, k)
    seen: set[str] = set()
    trims: list[tuple] = []
    dtr = Counter()
    w = c77.BgzfWriter(partial)
    try:
        for r in fasta.iter_records(fna):
            name = r.accession
            if name not in info:
                raise Refused(f"{fna}: record {name} is not in the information table")
            if name in seen:
                raise Refused(f"{fna}: record {name} occurs twice")
            seen.add(name)
            want = info[name][0]
            if len(r.seq) != want:
                raise Refused(f"{name}: {len(r.seq)} bases in the FASTA, {want} in the information table")
            src_cu = ids.contig_uid(r.seq)
            if is_dtr(info[name][4]):  # (length, checkv_quality, miuvig_quality, completeness, method, ...)
                seq, circ, cls, k = c62.apply_dec65("Direct terminal repeat", r.seq)
                dtr[cls] += 1
            else:
                seq, circ, cls, k = r.seq, False, "", 0
            cu = ids.contig_uid(seq) if k else src_cu
            if k:
                trims.append(("tower_gwv1", name, name, src_cu.hex(), cu.hex(), len(r.seq), k, k, len(seq)))
            w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
            recs.append((name, cu, len(seq), len(r.seq), src_cu, circ, cls, k))
    finally:
        w.close()
    missing = sorted(set(info) - seen)
    if missing:
        raise Refused(f"{len(missing)} contigs of the information table have no FASTA record; first {missing[0]}")

    records = tmp / "records.tsv"
    with open(records, "w", encoding="utf-8", buffering=1 << 20) as fh:
        fh.write("\t".join(GWV1_RECORDS_HEAD) + "\n")
        for i, (name, cu, n, src_n, src_cu, circ, cls, k) in enumerate(recs):
            (_, cq, mq, comp, meth, tax, life, sample, zone, bios) = info[name]
            fh.write("\t".join([str(i), name, cu.hex(), str(n), cq, mq, comp, meth, tax, life, sample, zone, bios,
                                str(circ).lower(), cls, src_cu.hex(), str(src_n), str(k)]) + "\n")
    owners, later = _owners(records, root, root_stages, owner_inputs, later_inputs, work_dir, threads, memory_limit)

    listing = make_inputs.Listing()
    genomes, repeated, shared, held_later = [], [], [], []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    by_quality: dict[str, list[int]] = {}
    for i, (name, cu, n, _src_n, _src_cu, circ, _cls, _k) in enumerate(recs):
        contig_len.setdefault(cu, n)
        q = by_quality.setdefault(info[name][1], [0, 0])
        q[0] += 1
        q[1] += n
        if cu in owners:
            shared.append((i, name, cu.hex(), n, owners[cu]))
        elif cu in later:
            held_later.append((i, name, cu.hex(), n, later[cu]))
        g = listing.add([cu], 11, "virus", "vOTU")
        if g is None:
            repeated.append((i, name, ids.genome_uid([cu]).hex(), first_of[cu]))
            continue
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first_of[cu] = name
    _write_tsv(tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    _write_tsv(tmp / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "length", "owner_stage"], shared)
    _write_tsv(tmp / "later_records.tsv", ["ordinal", "genome", "contig_uid", "length", "holder"], held_later)
    _write_tsv(tmp / "trims.tsv", c62.TRIM_HEAD, trims)
    obj = {**GWV1_SOURCE, "retrieved_at": retrieved_at}
    (tmp / "objects.json").write_text(json.dumps([{"path": str(inputs_home / GWV1 / GWV1_FASTA), **obj}], indent=1))
    (tmp / "objects_staging.json").write_text(json.dumps([{"path": str(final_fa), **obj}], indent=1))
    (tmp / "genomes.json").write_text(json.dumps(genomes))

    own_c, own_b = _tally(contig_len, owners)
    lat_c, lat_b = _tally(contig_len, later, owners)
    return {
        "records": len(recs), "source_bases": sum(r[3] for r in recs), "staged_bases": sum(r[2] for r in recs),
        "by_checkv_quality": {k: {"records": v[0], "bp": v[1]} for k, v in sorted(by_quality.items())},
        "dtr_records": sum(dtr.values()), "dtr_classes": dict(sorted(dtr.items())), "trims": len(trims),
        "trimmed_bases": sum(t[7] for t in trims), "circular_genomes": sum(1 for g in genomes if g["circular"]),
        "genomes": len(genomes), "genomes_repeated": len(repeated), "contigs_distinct": len(contig_len),
        "bases_distinct": sum(contig_len.values()),
        "held_by": own_c, "held_bases_by": own_b, "later_held_by": lat_c, "later_held_bases_by": lat_b,
        "contigs_new": sum(1 for cu in contig_len if cu not in owners),
        "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
        "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
        "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
        **{"listing_" + k: v for k, v in listing.report().items()},
    }


# --- the vOTU threshold on a sample (card C-48c's method) ----------------------------------------------------------

def pair_meets(ani: float, af_ref: float, af_query: float) -> dict[str, bool]:
    """Card C-39's three rules for one skani pair. `shorter` (the vOTU threshold, ANI >= 95 over >= 85% of the
    shorter genome) is read here as the larger of the two aligned fractions: both count the same aligned bases, so the
    shorter genome's fraction is the larger one, and no length of the held contig is needed."""
    return {"shorter": ani >= c39_ani.ANI_MIN and max(af_ref, af_query) >= c39_ani.AF_MIN,
            "both": ani >= c39_ani.ANI_MIN and min(af_ref, af_query) >= c39_ani.AF_MIN,
            "near_identical": ani >= 99.0 and min(af_ref, af_query) >= 90.0}


def sample_hits(table: Path, sample: set[str]) -> dict[str, set[str]]:
    """rule -> the sample records that meet it with a record outside the sample, in a skani table whichever side
    (query or reference) the sample was on."""
    out: dict[str, set[str]] = {"shorter": set(), "both": set(), "near_identical": set()}
    with open(table, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, ref = r["Query_name"], r["Ref_name"]
            s = q if q in sample and ref not in sample else ref if ref in sample and q not in sample else None
            if s is None:
                continue
            for rule, ok in pair_meets(float(r["ANI"]), float(r["Align_fraction_ref"]),
                                       float(r["Align_fraction_query"])).items():
                if ok:
                    out[rule].add(s)
    return out


def gwv1_sample(meta: Path, fasta_path: Path, out_fa: Path, every: int) -> tuple[dict[str, int], dict[str, str]]:
    """Every EVERY-th non-repeated genome of GWV1 in FASTA order: name -> length, and name -> CheckV quality."""
    with open(meta / "repeated_records.tsv", newline="") as fh:
        repeated = {r["genome"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(meta / "records.tsv", newline="") as fh:
        quality = {r["genome"]: r["checkv_quality"] for r in csv.DictReader(fh, delimiter="\t")}
    sample, k = {}, 0
    with open(out_fa, "w") as w:
        for r in fasta.iter_records(Path(fasta_path)):
            if r.accession in repeated:
                continue
            k += 1
            if k % every:
                continue
            w.write(f">{r.accession}\n{r.seq.decode()}\n")
            sample[r.accession] = len(r.seq)
    return sample, {s: quality.get(s, "") for s in sample}


def ani_summary(hits: dict[str, set[str]], sample: dict[str, int], quality: dict[str, str]) -> dict:
    out: dict = {"queries": len(sample)}
    for rule, met in hits.items():
        by_q, by_l = {}, {}
        for s, n in sample.items():
            for d, key in ((by_q, quality.get(s, "")), (by_l, c39_ani.length_class(n))):
                e = d.setdefault(key, {"queries": 0, "meet": 0})
                e["queries"] += 1
                e["meet"] += s in met
        out[rule] = {"meet": len(met), "by_checkv_quality": dict(sorted(by_q.items())),
                     "by_length": dict(sorted(by_l.items()))}
    return out


def fasta_chunks(paths: list[Path], chunk_bytes: int, work: Path):
    """The records of `paths` in order, in chunk files of about `chunk_bytes` of sequence each (a record is never
    split), yielded one at a time; each chunk file is removed when the caller asks for the next."""
    k, cur, n, path = 0, None, 0, None
    try:
        for p in paths:
            for r in fasta.iter_records(Path(p)):
                if cur is None:
                    path = work / f"chunk_{k:03d}.fna"
                    cur, n = open(path, "w"), 0
                cur.write(f">{r.accession}\n{r.seq.decode()}\n")
                n += len(r.seq)
                if n >= chunk_bytes:
                    cur.close()
                    cur = None
                    yield path
                    path.unlink(missing_ok=True)
                    k += 1
        if cur is not None:
            cur.close()
            cur = None
            yield path
            path.unlink(missing_ok=True)
    finally:
        if cur is not None:
            cur.close()
        if path is not None:
            path.unlink(missing_ok=True)


def gwv1_ani(meta: Path, fasta_path: Path, out: Path, refs: dict[str, list[Path]], *, every: int = 25,
             threads: int = 8, chunk_bytes: int = 4_000_000_000) -> dict:
    """The GWV1 sample against each held stage's whole FASTA (the sample as skani's reference set, the held stage's
    records as queries, so that skani holds the small side's sketches; the held stage is read in chunks of about
    CHUNK_BYTES, one skani process each, so that no process holds or spends much), and the union over held stages."""
    out = Path(out)
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)
    sample_fa = work / "gwv1_sample.fna"
    sample, quality = gwv1_sample(Path(meta), Path(fasta_path), sample_fa, every)
    res: dict = {"every": every, "threads": threads, "samples": len(sample), "sample_bases": sum(sample.values()),
                 "chunk_bytes": chunk_bytes,
                 "threshold": {"ani_min": c39_ani.ANI_MIN, "af_min": c39_ani.AF_MIN}, "by_stage": {}}
    union: dict[str, set[str]] = {"shorter": set(), "both": set(), "near_identical": set()}
    for name, paths in refs.items():
        h: dict[str, set[str]] = {"shorter": set(), "both": set(), "near_identical": set()}
        cpus = []
        for i, chunk in enumerate(fasta_chunks([Path(p) for p in paths], chunk_bytes, work)):
            table = out / f"{name}.{i:03d}.tsv"
            cpus.append(round(c39_ani.skani(threads, [chunk], [sample_fa], table), 1))
            for rule, met in sample_hits(table, set(sample)).items():
                h[rule] |= met
        for rule in union:
            union[rule] |= h[rule]
        res["by_stage"][name] = {"fasta": [str(p) for p in paths], "skani_runs": len(cpus), "cpu_s": cpus,
                                 **ani_summary(h, sample, quality)}
    res["any_held_stage"] = ani_summary(union, sample, quality)
    (out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    return res


# === SAG3 ==========================================================================================================

SAG3 = "SAG3"
SAG3_FASTA = "sag3_omz.fna.gz"
SAG3_SOURCE = {
    "source_name": "tower_sag3",
    "source_release": "OMZ single-cell compendium, figshare 20459526 (2023-04-06) OMZ_Compendium_SAGs_Full.tar.gz and "
                      "figshare 20481603 (2023-05-23) Table S1",
    "license": "CC BY 4.0 as both figshare articles state it (read 2026-11-26, card C-72 part 1)",
    "redistributable": True,  # DEC-80 (3): a recorded licence (CC BY 4.0) permits it
}
SAG3_RETRIEVED = "2026-11-26T18:27:24Z"  # tower's SOURCE_URLS.tsv, OMZ_Compendium_SAGs_Full.tar.gz
MEMBER_DIR = "Final_Catalog/"
SHEET_GTDBTK = "Outputs_GTDB-TK"
SHEET_SEQ = "Outputs-Genome_Sequencing"
SHEET_QUAL = "Outputs-Genome_Quality_Classifi"
# Outputs-Genome_Sequencing repeats its header names, so its columns are taken by position and their names checked
SEQ_COLS = {"name": (0, "Sample Name"), "biosample_repo": (15, "Data Repository"),
            "biosample": (16, "Data Record Accession"), "genbank_repo": (19, "Data Repository"),
            "genbank": (20, "Data Record Accession"), "img_repo": (23, "Data Repository"),
            "img": (24, "Data Record Accession")}
BIOSAMPLE_RE = re.compile(r"SAM[NED][A-Z]?\d+")
WGS_RE = re.compile(r"([A-Z]{4,6})\d{8,10}(\.\d+)?")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "p": "http://schemas.openxmlformats.org/package/2006/relationships"}
SAG3_META = ("genomes.json", "objects.json", "objects_staging.json", "records.tsv", "repeated_records.tsv",
             "shared_records.tsv", "later_records.tsv", "genetic_codes.tsv", "sag_match.tsv", "r232_status.tsv")


def _col_index(ref: str) -> int:
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_xlsx(path: Path) -> dict[str, list[list[str]]]:
    """Every sheet of an xlsx as rows of strings (shared and inline strings resolved, missing cells empty), by sheet
    name through the workbook's relationships; the standard library only (zipfile, xml.etree)."""
    try:
        z = zipfile.ZipFile(path)
    except (zipfile.BadZipFile, OSError) as exc:
        raise Refused(f"{path}: not an xlsx: {exc}") from exc
    with z:
        names = set(z.namelist())
        ss: list[str] = []
        if "xl/sharedStrings.xml" in names:
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
                ss.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
        rels = {r.get("Id"): r.get("Target") for r in
                ET.fromstring(z.read("xl/_rels/workbook.xml.rels")).findall("p:Relationship", NS)}
        out: dict[str, list[list[str]]] = {}
        for s in ET.fromstring(z.read("xl/workbook.xml")).find("m:sheets", NS):
            target = rels[s.get(f"{{{NS['r']}}}id")].lstrip("/")
            target = target if target.startswith("xl/") else "xl/" + target
            rows: list[list[str]] = []
            for r in ET.fromstring(z.read(target)).find("m:sheetData", NS).findall("m:row", NS):
                vals: dict[int, str] = {}
                for c in r.findall("m:c", NS):
                    v = c.find("m:v", NS)
                    t = c.get("t")
                    x = v.text if v is not None and v.text is not None else ""
                    if t == "s" and x:
                        x = ss[int(x)]
                    elif t == "inlineStr":
                        x = "".join(tt.text or "" for tt in c.iter(f"{{{NS['m']}}}t"))
                    vals[_col_index(c.get("r"))] = x
                rows.append([vals.get(k, "") for k in range(max(vals) + 1)] if vals else [])
            out[s.get("name")] = rows
    return out


def _cell(row: list[str], i: int) -> str:
    return row[i].strip() if i < len(row) else ""


def table_s1(xlsx: Path) -> dict[str, dict]:
    """Table S1, per cell named in its GTDB-Tk sheet: the classification, BioSample, WGS accession and IMG id
    (Outputs-Genome_Sequencing) and the quality class; refuses a sheet or a column that is not where it should be."""
    sheets = read_xlsx(Path(xlsx))
    for s in (SHEET_GTDBTK, SHEET_SEQ, SHEET_QUAL):
        if s not in sheets:
            raise Refused(f"{xlsx}: no sheet {s}")
    out: dict[str, dict] = {}

    def entry(name: str) -> dict:
        return out.setdefault(name, {"classification": "", "fastani_reference": "", "biosample": "", "wgs": "",
                                     "img": "", "quality": "", "seq_rows": 0, "gtdbtk_rows": 0})

    gt = sheets[SHEET_GTDBTK]
    h = gt[0] if gt else []
    if "user_genome" not in h or "classification" not in h:
        raise Refused(f"{xlsx}: sheet {SHEET_GTDBTK} has no user_genome or classification column")
    iu, ic = h.index("user_genome"), h.index("classification")
    ifr = h.index("fastani_reference") if "fastani_reference" in h else None
    for r in gt[1:]:
        name = _cell(r, iu)
        if not name:
            continue
        d = entry(name)
        d["gtdbtk_rows"] += 1
        fr = _cell(r, ifr) if ifr is not None else ""
        if d["gtdbtk_rows"] == 1:
            d["classification"], d["fastani_reference"] = _cell(r, ic), ("" if fr in ("N/A", "NA") else fr)
        else:  # the sheet gives this name to two cells (it cuts `<cell>_AB-905` and `<cell>_AB-908` short): neither
            d["classification"] = d["fastani_reference"] = ""
    seq = sheets[SHEET_SEQ]
    sh = seq[0] if seq else []
    for key, (i, want) in SEQ_COLS.items():
        if _cell(sh, i) != want:
            raise Refused(f"{xlsx}: {SHEET_SEQ} column {i} is {_cell(sh, i)!r}, not {want!r}")
    for r in seq[1:]:
        name = _cell(r, 0)
        if not name:
            continue
        d = entry(name)
        d["seq_rows"] += 1
        bs, gb, img = _cell(r, SEQ_COLS["biosample"][0]), _cell(r, SEQ_COLS["genbank"][0]), _cell(r, SEQ_COLS["img"][0])
        if BIOSAMPLE_RE.fullmatch(bs) and not d["biosample"]:
            d["biosample"] = bs
        if WGS_RE.fullmatch(gb) and not d["wgs"]:
            d["wgs"] = gb
        if img.isdigit() and not d["img"]:
            d["img"] = img
    qs = sheets[SHEET_QUAL]
    qh = qs[0] if qs else []
    if "Sample ID" in qh and "Quality Classification" in qh:
        i_n, i_q = qh.index("Sample ID"), qh.index("Quality Classification")
        rows = [(_cell(r, i_n), _cell(r, i_q)) for r in qs[1:] if _cell(r, i_n)]
        once = Counter(n for n, _ in rows)
        for n, q in rows:
            if n in out and once[n] == 1:
                out[n]["quality"] = q
    return out


def wgs_project(acc: str) -> str:
    """`AFHZ00000000`, `AQSI01000000`, `CABWYC000000000.1` -> the WGS project's letters (`AFHZ`, `AQSI`, `CABWYC`)."""
    m = WGS_RE.fullmatch(acc.strip())
    return m.group(1) if m else ""


def sag3_members(tar_path: Path) -> tuple[list[tuple[str, str, int]], dict]:
    """(member name, stem, size) of each single-cell FASTA in archive order, and the census of every entry."""
    mem, other, dirs = [], [], 0
    try:
        with tarfile.open(tar_path, "r:gz") as tf:
            for m in tf:
                if m.isdir():
                    dirs += 1
                elif m.isfile() and m.name.startswith(MEMBER_DIR) and m.name.endswith(".fna") and "/" not in \
                        m.name[len(MEMBER_DIR):]:
                    mem.append((m.name, m.name[len(MEMBER_DIR):-len(".fna")], m.size))
                else:
                    other.append(m.name)
    except (tarfile.TarError, OSError, EOFError, zlib.error) as exc:
        raise Refused(f"{tar_path}: not a readable tar.gz: {exc}") from exc
    stems = Counter(s for _, s, _ in mem)
    dup = sorted(s for s, n in stems.items() if n > 1)
    if dup:
        raise Refused(f"{tar_path}: two members for {dup[0]}")
    return mem, {"entries": len(mem) + len(other) + dirs, "directories": dirs, "cells": len(mem),
                 "other": sorted(other), "bytes": sum(s for _, _, s in mem)}


def match_members(stems: list[str], names: set[str]) -> dict[str, tuple[str | None, str]]:
    """stem -> (Table S1 name, basis): the longest name equal to the stem (`exact`) or followed in it by `_`
    (`prefix`); a name that two stems would take is given to neither (`tie`); `none` when no name fits."""
    pick: dict[str, str | None] = {}
    for s in stems:
        best = None
        for n in names:
            if (s == n or s.startswith(n + "_")) and (best is None or len(n) > len(best)):
                best = n
        pick[s] = best
    taken = Counter(n for n in pick.values() if n is not None)
    out = {}
    for s, n in pick.items():
        if n is None:
            out[s] = (None, "none")
        elif taken[n] > 1:
            out[s] = (None, "tie")
        else:
            out[s] = (n, "exact" if s == n else "prefix")
    return out


def gtdb_crosswalk(gtdb_meta_dir: Path, biosamples: set[str], projects: set[str]) -> tuple[dict, dict]:
    """From GTDB R232's metadata (ar53, bac120): BioSample -> sorted versionless accessions, WGS project ->
    sorted versionless accessions, for the keys asked only (one pass over each table, two columns kept)."""
    import pyarrow as pa
    import pyarrow.csv as pacsv

    by_bs: dict[str, set[str]] = defaultdict(set)
    by_wgs: dict[str, set[str]] = defaultdict(set)
    cols = ["accession", "ncbi_biosample", "ncbi_wgs_master"]
    for key in ("ar53", "bac120"):
        p = Path(gtdb_meta_dir) / r232_index.FILES[key]
        if not p.is_file():
            raise Refused(f"missing GTDB R232 metadata {p}")
        rd = pacsv.open_csv(str(p), read_options=pacsv.ReadOptions(block_size=1 << 24),
                            parse_options=pacsv.ParseOptions(delimiter="\t"),
                            convert_options=pacsv.ConvertOptions(include_columns=cols,
                                                                 column_types={c: pa.string() for c in cols}))
        for b in rd:
            for acc, bs, wgs in zip(*(b.column(i).to_pylist() for i in range(3))):
                a = r232_index.nover(acc or "")
                if bs and bs in biosamples:
                    by_bs[bs].add(a)
                pr = wgs_project(wgs or "")
                if pr and pr in projects:
                    by_wgs[pr].add(a)
    return {k: sorted(v) for k, v in by_bs.items()}, {k: sorted(v) for k, v in by_wgs.items()}


STATUS_HEAD = ["genome", "biosample", "wgs", "gtdb_r232_accession", "via", "r232_status", "in_gtdb_r232",
               "is_gtdb_representative", "gtdb_representative", "gtdb_species", "in_globdb_r232", "globdb_r232_id",
               "gtdbtk_classification", "gtdbtk_fastani_reference", "contigs", "bp"]
MATCH_HEAD = ["member", "stem", "genome", "basis", "contigs", "bp", "biosample", "wgs", "img", "quality"]


def r232_status(x: dict | None, has_accession: bool) -> str:
    """`not_visible`: Table S1 gives the cell no BioSample or WGS accession (or the member matched no row);
    `not_in_gtdb_r232`: it does, and no GTDB R232 genome carries it (GlobDB r232 cannot be asked without a GCA_);
    else C-74's `lookup` of the GTDB R232 accession."""
    if not has_accession:
        return "not_visible"
    if x is None:
        return "not_in_gtdb_r232"
    if x["in_gtdb_r232"] and x["in_globdb_r232"]:
        return "in_gtdb_and_globdb_r232"
    if x["in_gtdb_r232"]:
        return "in_gtdb_r232"
    return "in_globdb_r232" if x["in_globdb_r232"] else "not_in_r232"


def _read_records(fh, name: str) -> list[tuple[str, bytes, int]]:
    out: list[tuple[str, bytes, int]] = []
    header: str | None = None
    lines: list[bytes] = []
    for ln in fh:
        ln = ln.rstrip(b"\r\n")
        if not ln:
            continue
        if ln.startswith(b">"):
            if header is not None and lines:
                out.append((header, b"".join(lines), fasta._width(lines)))
            header, lines = ln[1:].decode("utf-8", "replace"), []
        elif header is None:
            raise Refused(f"{name}: sequence before the first header")
        else:
            lines.append(ln)
    if header is not None and lines:
        out.append((header, b"".join(lines), fasta._width(lines)))
    return out


def stage_sag3(tar_path: Path, xlsx: Path, index_dir: Path, gtdb_meta_dir: Path, meta_out: Path, fasta_out: Path, *,
               root: Path | None, root_stages, owner_inputs: dict[str, Path], later_inputs: dict[str, Path],
               inputs_home: Path, retrieved_at: str = SAG3_RETRIEVED, work: Path | None = None, threads: int = 2,
               memory_limit: str = "4GB", workers: int = 1) -> dict:
    """Stage SAG3; returns the report (module docstring). Nothing is left behind when it refuses."""
    tar_path, index_dir = Path(tar_path), Path(index_dir)
    meta_dir, fa_dir = Path(meta_out) / SAG3, Path(fasta_out) / SAG3
    final_fa = fa_dir / SAG3_FASTA
    for p in (meta_dir, final_fa):
        if p.exists():
            raise Refused(f"{p} exists")
    if not (index_dir / "gtdb_r232.parquet").is_file() or not (index_dir / "globdb_r232.parquet").is_file():
        raise Refused(f"{index_dir} is not card C-74's index (no gtdb_r232.parquet or globdb_r232.parquet)")
    _check_later(later_inputs)
    s1 = table_s1(Path(xlsx))
    mem, census = sag3_members(tar_path)
    if not mem:
        raise Refused(f"{tar_path}: no single-cell member ({MEMBER_DIR}*.fna)")
    matched = match_members([s for _, s, _ in mem], set(s1))
    genome_of = {s: (matched[s][0] or s) for _, s, _ in mem}
    dup = [g for g, n in Counter(genome_of.values()).items() if n > 1]
    if dup:
        raise Refused(f"two members would stage as genome {sorted(dup)[0]}")
    tmp = Path(meta_out) / f"{SAG3}.tmp"
    partial = fa_dir / f".{SAG3_FASTA}.partial"
    work_dir = Path(work) if work is not None else Path(meta_out) / f"{SAG3}.work"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    fa_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        rep = _stage_sag3(tar_path, s1, mem, census, matched, genome_of, index_dir, Path(gtdb_meta_dir), tmp,
                          partial, work_dir, root, list(root_stages), owner_inputs, later_inputs, Path(inputs_home),
                          final_fa, retrieved_at, threads, memory_limit, workers)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        partial.unlink(missing_ok=True)
        if fa_dir.exists() and not any(fa_dir.iterdir()):
            fa_dir.rmdir()
        raise
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    partial.rename(final_fa)
    tmp.rename(meta_dir)
    files = {n: meta_dir / n for n in SAG3_META}
    files[SAG3_FASTA] = final_fa
    rep["sha256"] = {n: ids.sha256_file(p).hex() for n, p in files.items()}
    rep["bytes"] = {n: p.stat().st_size for n, p in files.items()}
    return rep


def _stage_sag3(tar_path, s1, mem, census, matched, genome_of, index_dir, gtdb_meta_dir, tmp, partial, work_dir,
                root, root_stages, owner_inputs, later_inputs, inputs_home, final_fa, retrieved_at, threads,
                memory_limit, workers) -> dict:
    # Table S1's row of each genome, through the member's match only (a tied stem may equal a Table S1 name)
    row = {genome_of[s]: s1[n] for s, (n, _) in matched.items() if n is not None}
    hint = {g: (row[g]["classification"] if g in row else "") for g in genome_of.values()}
    gc_results: dict[str, dict] = {}
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    pending: list[tuple[str, object]] = []

    def submit(g: str, seqs: list[bytes]) -> None:
        if pool is None:
            gc_results[g] = c58.screen_genetic_code(seqs, hint[g], screen_max_bp=500_000)
            return
        pending.append((g, pool.submit(c58._screen_worker, (g, seqs, hint[g], 500_000))))
        while len(pending) >= 4 * workers:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]

    contigs: dict[str, tuple[bytes, list[int]]] = {}
    member_of: dict[str, str] = {}
    names: set[str] = set()
    n_renamed = 0
    want = {m: s for m, s, _ in mem}
    w = c77.BgzfWriter(partial)
    try:
        with tarfile.open(tar_path, "r:gz") as tf:
            for m in tf:
                if m.name not in want:
                    continue
                g = genome_of[want[m.name]]
                fh = tf.extractfile(m)
                if fh is None:
                    raise Refused(f"{m.name}: not a file")
                recs = _read_records(fh, m.name)
                for hdr, seq, width in recs:
                    first = hdr.split(None, 1)[0] if hdr.strip() else ""
                    if not first.startswith(g + "_"):
                        hdr = f"{g}|{hdr}"
                        n_renamed += 1
                    first = hdr.split(None, 1)[0]
                    if first in names:
                        raise Refused(f"contig name {first} occurs twice in the stage")
                    names.add(first)
                    w.write(b"\n".join(fasta.format_record(hdr, seq, width)) + b"\n")
                contigs[g] = (b"".join(ids.contig_uid(s) for _, s, _ in recs), [len(s) for _, s, _ in recs])
                member_of[g] = m.name
                submit(g, [s for _, s, _ in recs])
        while pending:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]
    except (tarfile.TarError, OSError, EOFError, zlib.error) as exc:
        raise Refused(f"{tar_path}: not a readable tar.gz: {exc}") from exc
    finally:
        w.close()
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    lost = sorted(set(genome_of.values()) - set(contigs))
    if lost:
        raise Refused(f"{len(lost)} members were listed but not read; first {lost[0]}")
    empty = sorted(g for g, (_, lens) in contigs.items() if not lens)
    if empty:
        raise Refused(f"{len(empty)} members hold no sequence; first {empty[0]}")
    kept = sorted(contigs)

    # the accession test
    bios = {row[g]["biosample"] for g in kept if g in row and row[g]["biosample"]}
    projs = {wgs_project(row[g]["wgs"]) for g in kept if g in row and row[g]["wgs"]}
    by_bs, by_wgs = gtdb_crosswalk(gtdb_meta_dir, bios, projs)
    via: dict[str, tuple[str, str]] = {}
    for g in kept:
        d = row.get(g)
        if d is None:
            continue
        pr = wgs_project(d["wgs"]) if d["wgs"] else ""
        if pr and by_wgs.get(pr):
            via[g] = (by_wgs[pr][0], "wgs")
        elif d["biosample"] and by_bs.get(d["biosample"]):
            via[g] = (by_bs[d["biosample"]][0], "biosample")
    look = {x["accession"]: x for x in r232_index.lookup(index_dir, sorted({a for a, _ in via.values()}))}
    has_acc = {g: bool(g in row and (row[g]["biosample"] or row[g]["wgs"])) for g in kept}
    status = {g: r232_status(look.get(via[g][0]) if g in via else None, has_acc[g]) for g in kept}

    resolved, genome_uids, conflicts = c58._resolve_within_stage_code_conflicts(kept, contigs, gc_results)
    records_path = tmp / "records.tsv"
    ordinal = 0
    with open(records_path, "w", encoding="utf-8", buffering=1 << 20) as rf, open(
            tmp / "genetic_codes.tsv", "w", encoding="utf-8", buffering=1 << 20) as gf:
        rf.write("ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\tbiosample\n")
        gf.write("\t".join(c58.GC_HEAD) + "\n")
        for g in kept:
            raw, lens = contigs[g]
            gc = gc_results[g]
            gf.write("\t".join([
                g, SAG3, hint[g], str(sum(lens)), str(len(lens)), "contig", str(resolved[g]), gc["reason"],
                str(gc["tested_alt"]).lower(),
                c58._fmt_opt(gc["cd_11"]), c58._fmt_opt(gc["mean_orf_11"]), c58._fmt_opt(gc["n_genes_11"]),
                c58._fmt_opt(gc["cd_25"]), c58._fmt_opt(gc["mean_orf_25"]), c58._fmt_opt(gc["n_genes_25"]),
                c58._fmt_opt(gc["cd_4"]), c58._fmt_opt(gc["mean_orf_4"]), c58._fmt_opt(gc["n_genes_4"]),
                c58._fmt_opt(gc["tga_per_kb"])]) + "\n")
            bs = row[g]["biosample"] if g in row else ""
            for k, n in enumerate(lens):
                rf.write(f"{ordinal}\t{g}\t{k}\t{raw[k * 16:(k + 1) * 16].hex()}\t{n}\tcontig\t{resolved[g]}\t{bs}\n")
                ordinal += 1
    owners, later = _owners(records_path, root, root_stages, owner_inputs, later_inputs, work_dir, threads,
                            memory_limit)

    listing = make_inputs.Listing()
    genomes, repeated, shared, held_later = [], [], [], []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ordinal = 0
    for g in kept:
        raw, lens = contigs[g]
        uids = [raw[i: i + 16] for i in range(0, len(raw), 16)]
        for k, (cu, n) in enumerate(zip(uids, lens)):
            contig_len.setdefault(cu, n)
            if cu in owners:
                shared.append((ordinal, g, k, cu.hex(), n, owners[cu]))
            elif cu in later:
                held_later.append((ordinal, g, k, cu.hex(), n, later[cu]))
            ordinal += 1
        e = listing.add(uids, resolved[g], "prokaryote", "contig")
        if e is None:
            repeated.append((ordinal - len(uids), g, genome_uids[g].hex(), first_of[genome_uids[g]]))
            continue
        genomes.append(e)
        first_of[genome_uids[g]] = g
    _write_tsv(tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    _write_tsv(tmp / "shared_records.tsv",
               ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"], shared)
    _write_tsv(tmp / "later_records.tsv", ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "holder"],
               held_later)
    stem_of = {genome_of[s]: s for _, s, _ in mem}
    _write_tsv(tmp / "sag_match.tsv", MATCH_HEAD, [
        (member_of[g], stem_of[g], g, matched[stem_of[g]][1], len(contigs[g][1]), sum(contigs[g][1]),
         *(row.get(g, {}).get(k, "") for k in ("biosample", "wgs", "img", "quality"))) for g in kept])

    def yn(v: bool) -> str:
        return "yes" if v else "no"

    rows = []
    for g in kept:
        x = look.get(via[g][0]) if g in via else None
        d = row.get(g, {})
        rows.append((g, d.get("biosample", ""), d.get("wgs", ""), via[g][0] if g in via else "",
                     via[g][1] if g in via else "", status[g], yn(bool(x and x["in_gtdb_r232"])),
                     str(x["is_representative"]) if x and x["in_gtdb_r232"] else "",
                     (x["representative"] or "") if x else "", (x["species"] or "") if x else "",
                     yn(bool(x and x["in_globdb_r232"])), (x["globdb_r232_id"] or "") if x else "",
                     d.get("classification", ""), d.get("fastani_reference", ""), len(contigs[g][1]),
                     sum(contigs[g][1])))
    _write_tsv(tmp / "r232_status.tsv", STATUS_HEAD, rows)

    obj = {**SAG3_SOURCE, "retrieved_at": retrieved_at}
    (tmp / "objects.json").write_text(json.dumps([{"path": str(inputs_home / SAG3 / SAG3_FASTA), **obj}], indent=1))
    (tmp / "objects_staging.json").write_text(json.dumps([{"path": str(final_fa), **obj}], indent=1))
    (tmp / "genomes.json").write_text(json.dumps(genomes))

    own_c, own_b = _tally(contig_len, owners)
    lat_c, lat_b = _tally(contig_len, later, owners)
    bp = {g: sum(contigs[g][1]) for g in kept}
    return {
        "census": census, "cells": len(kept), "contigs": sum(len(contigs[g][1]) for g in kept),
        "bases": sum(bp.values()), "contigs_renamed": n_renamed,
        "table_s1_cells": len(s1), "match_basis": dict(sorted(Counter(b for _, b in matched.values()).items())),
        "table_s1_unmatched": len(set(s1) - {n for n, _ in matched.values() if n is not None}),
        "with_biosample": sum(1 for g in kept if g in row and row[g]["biosample"]),
        "with_wgs": sum(1 for g in kept if g in row and row[g]["wgs"]),
        "without_classification": sum(1 for g in kept if not hint[g]),
        "via": dict(sorted(Counter(v for _, v in via.values()).items())),
        "r232_status": dict(sorted(Counter(status.values()).items())),
        "r232_status_bases": {k: sum(bp[g] for g in kept if status[g] == k) for k in sorted(set(status.values()))},
        "gtdb_representatives": sum(1 for g in kept if g in via and look[via[g][0]]["is_representative"]),
        "genetic_codes": dict(sorted(Counter(str(resolved[g]) for g in kept).items())),
        "code_conflicts": len(conflicts),
        "genomes": len(genomes), "genomes_repeated": len(repeated), "contigs_distinct": len(contig_len),
        "bases_distinct": sum(contig_len.values()),
        "held_by": own_c, "held_bases_by": own_b, "later_held_by": lat_c, "later_held_bases_by": lat_b,
        "contigs_new": sum(1 for cu in contig_len if cu not in owners),
        "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
        "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
        "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
        **{"listing_" + k: v for k, v in listing.report().items()},
    }


# --- skani against SAG1 and SAG2 (recorded, not cut) ---------------------------------------------------------------

ANI_SAG, AF_SAG = 96.0, 50.0


def _rows_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _records_tsv(meta: Path) -> Path:
    """A staged set's records.tsv: in META, or (a file over 64 MB, AGENTS §7.3) beside its FASTA on /stor."""
    p = Path(meta) / "records.tsv"
    if p.is_file():
        return p
    for o in json.loads((Path(meta) / "objects_staging.json").read_text()):
        q = Path(o["path"]).parent / "records.tsv"
        if q.is_file():
            return q
    raise Refused(f"{meta}: no records.tsv")


def split_genomes(meta: Path, fasta_path, out_dir: Path, keep: set[str] | None = None,
                  unlisted: list[str] | None = None) -> dict[str, Path]:
    """One FASTA per genome of a staged set (contig -> genome by records.tsv's contig_uid, first genome listed),
    written under out_dir, only the genomes in KEEP when it is given; FASTA_PATH is one file or a list; returns
    genome -> file."""
    paths = [Path(fasta_path)] if isinstance(fasta_path, (str, Path)) else [Path(p) for p in fasta_path]
    first: dict[str, str] = {}
    with open(_records_tsv(Path(meta)), newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        col = "genome" if "genome" in (rd.fieldnames or []) else "genome_id"  # cards C-58 and C-66 name it genome_id
        for r in rd:
            first.setdefault(r["contig_uid"], r[col])
    out_dir.mkdir(parents=True, exist_ok=True)
    files: dict[str, Path] = {}
    used: set[str] = set()
    handles: dict[str, object] = {}
    try:
        for r in (rec for p in paths for rec in fasta.iter_records(p)):
            g = first.get(ids.contig_uid(r.seq).hex())
            if g is None:
                if unlisted is None:
                    raise Refused(f"{paths[0]}: record {r.accession} is in no genome of {meta}/records.tsv")
                unlisted.append(r.accession)  # the caller asked to count, not refuse (an ANI reference set)
                continue
            if keep is not None and g not in keep:
                continue
            if g not in files:
                fname = re.sub(r"[^A-Za-z0-9_.-]", "_", g) + ".fna"
                if fname in used:
                    raise Refused(f"two genomes would share the file {fname}")
                used.add(fname)
                files[g] = out_dir / fname
            fh = handles.get(g)
            if fh is None:
                fh = handles[g] = open(files[g], "ab")
                if len(handles) > 64:
                    k0 = next(iter(handles))
                    handles.pop(k0).close()
            fh.write(b">" + r.accession.encode() + b"\n" + r.seq + b"\n")  # bytes as staged, decoded never
    finally:
        for fh in handles.values():
            fh.close()
    return files


def best_hits(table: Path, query_of: dict[str, str], ref_of: dict[str, tuple[str, str]]) -> dict[str, dict]:
    """query genome -> its best reference (largest ANI, then larger aligned fraction) among the pairs at ANI >= 96
    and an aligned fraction (the larger of the two) >= 50, from a skani table of genome files."""
    out: dict[str, dict] = {}
    with open(table, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, ref = query_of.get(r["Query_file"]), ref_of.get(r["Ref_file"])
            if q is None or ref is None:
                continue
            ani, af = float(r["ANI"]), max(float(r["Align_fraction_ref"]), float(r["Align_fraction_query"]))
            if ani < ANI_SAG or af < AF_SAG:
                continue
            cur = out.get(q)
            if cur is None or (-ani, -af, ref[0], ref[1]) < (-cur["ani"], -cur["af"], cur["stage"], cur["genome"]):
                out[q] = {"ani": ani, "af": af, "stage": ref[0], "genome": ref[1]}
    return out


def genus_of(lineage: str) -> str:
    """`...;g__Nitrosopelagicus;s__...` -> `Nitrosopelagicus`; '' when the lineage names no genus."""
    for tok in lineage.split(";"):
        tok = tok.strip()
        if tok.startswith("g__"):
            return tok[3:]
    return ""


def archaea_refs(fna: Path, list_tsv: Path, genera: set[str], out_dir: Path) -> dict[str, Path]:
    """Card C-75's r226 archaea (one FASTA in list order, `first_record` and `records` per genome): the genomes of
    the given GTDB genera written one file each under out_dir; returns genome -> file."""
    want: dict[int, str] = {}
    with open(list_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if genus_of(r["taxonomy"]) in genera:
                for i in range(int(r["first_record"]), int(r["first_record"]) + int(r["records"])):
                    want[i] = r["genome"]
    out_dir.mkdir(parents=True, exist_ok=True)
    files: dict[str, Path] = {}
    last = max(want) if want else -1
    for i, r in enumerate(fasta.iter_records(Path(fna))):
        if i > last:
            break
        g = want.get(i)
        if g is None:
            continue
        p = files.setdefault(g, out_dir / (re.sub(r"[^A-Za-z0-9_.-]", "_", g) + ".fna"))
        with open(p, "ab") as fh:
            fh.write(b">" + r.accession.encode() + b"\n" + r.seq + b"\n")
    return files


def sag3_ani(meta: Path, fasta_path: Path, out: Path, refs: dict[str, tuple[Path, Path]], *, threads: int = 8,
             archaea: tuple[Path, Path] | None = None) -> dict:
    """SAG3's cells against whole staged sets (REFS: name -> (meta dir, FASTA)) and, with ARCHAEA (card C-75's
    r226 archaeal FASTA and its list), against the r226 archaea of the genera GTDB-Tk gives SAG3's archaeal cells in
    Table S1 (the card's same-genus references, here for archaea only: the bacterial ones need the r226 archive)."""
    out = Path(out)
    work = out / "work"
    shutil.rmtree(work, ignore_errors=True)
    q_files = split_genomes(Path(meta), Path(fasta_path), work / "SAG3")
    ref_of: dict[str, tuple[str, str]] = {}
    for name, (m, fa) in refs.items():
        for g, p in split_genomes(m, fa, work / name).items():
            ref_of[str(p)] = (name, g)
    genera: set[str] = set()
    if archaea is not None:
        for r in _rows_tsv(Path(meta) / "r232_status.tsv"):
            if r["gtdbtk_classification"].startswith("d__Archaea") and genus_of(r["gtdbtk_classification"]):
                genera.add(genus_of(r["gtdbtk_classification"]))
        for g, p in archaea_refs(archaea[0], archaea[1], genera, work / "r226_archaea").items():
            ref_of[str(p)] = ("r226_archaea", g)
    (work / "q.txt").write_text("".join(f"{p}\n" for p in sorted(map(str, q_files.values()))))
    (work / "r.txt").write_text("".join(f"{p}\n" for p in sorted(ref_of)))
    table = out / "sag3_vs_refs.tsv"
    t0 = c39_ani.child_cpu()
    subprocess.run(["skani", "dist", "-t", str(threads), "--min-af", str(AF_SAG), "--ql", str(work / "q.txt"),
                    "--rl", str(work / "r.txt"), "-o", str(table)], check=True)
    cpu = c39_ani.child_cpu() - t0
    hits = best_hits(table, {str(p): g for g, p in q_files.items()}, ref_of)
    _write_tsv(out / "sag3_ani_hits.tsv", ["genome", "stage", "ref_genome", "ani", "af"],
               [(g, h["stage"], h["genome"], h["ani"], h["af"]) for g, h in sorted(hits.items())])
    res = {"queries": len(q_files), "references": dict(sorted(Counter(s for s, _ in ref_of.values()).items())),
           "archaeal_genera": sorted(genera),
           "threshold": {"ani_min": ANI_SAG, "af_min_larger_of_two": AF_SAG}, "threads": threads,
           "cpu_s": round(cpu, 1), "meet": len(hits),
           "meet_by_stage": dict(sorted(Counter(h["stage"] for h in hits.values()).items()))}
    (out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    shutil.rmtree(work, ignore_errors=True)
    return res


# === WET1 ==========================================================================================================

WET1 = "WET1"
WET1_FASTA = "wet1_gwm.fna.gz"
WET1_SOURCE = {
    "source_name": "tower_wet1",
    "source_release": "GWM global wetland catalogue, zenodo 20263536 (2026-05): 'Dereplicated genomes.zip' (on the pod "
                      "as Dereplicated_genomes.zip, relay_to_pod.sh --space-as-underscore) and MAG_meta.xlsx",
    "license": "cc-by-4.0 as the Zenodo record states it (read 2026-11-26, card C-72 part 1)",
    "redistributable": True,  # DEC-80 (3): a recorded licence (CC BY 4.0) permits it
}
WET1_RETRIEVED = "2026-11-26T19:04:00Z"  # tower's SOURCE_URLS.tsv, Dereplicated genomes.zip
WET1_DIR = "dereplicated_genomes/"
WET1_BIOPROJECT = "PRJNA1430164"
WET1_META = ("genomes.json", "objects.json", "objects_staging.json", "records.tsv", "repeated_records.tsv",
             "shared_records.tsv", "later_records.tsv", "genetic_codes.tsv", "r232_status.tsv")
MAG_NEED = {"id": "ID", "name": "Name", "sample": "Sample ID", "classification": "GTDB_classification",
            "completeness": "Completeness (%)", "contamination": "Contamination (%)", "quality": "Quality_tier (MIMAG)"}


def mag_meta(xlsx: Path) -> dict[str, dict]:
    """MAG_meta.xlsx (its one sheet: a title row, then the header row that starts with `ID`): MAG id -> name, sample,
    GTDB classification, completeness, contamination, quality tier. Refuses a missing column or a repeated id."""
    sheets = read_xlsx(Path(xlsx))
    rows = next(iter(sheets.values()), [])
    hi = next((i for i, r in enumerate(rows) if r and r[0].strip() == "ID"), None)
    if hi is None:
        raise Refused(f"{xlsx}: no header row starting with ID")
    head = [c.strip() for c in rows[hi]]
    miss = [c for c in MAG_NEED.values() if c not in head]
    if miss:
        raise Refused(f"{xlsx}: no column {miss[0]!r}")
    ix = {k: head.index(c) for k, c in MAG_NEED.items()}
    out: dict[str, dict] = {}
    for r in rows[hi + 1:]:
        mid = _cell(r, ix["id"])
        if not mid:
            continue
        if mid in out:
            raise Refused(f"{xlsx}: MAG {mid} has two rows")
        out[mid] = {k: _cell(r, i) for k, i in ix.items() if k != "id"}
    return out


def mimag_quality(tier: str) -> str:
    """MAG_meta's MIMAG tier as the stager's `assembly_quality` (`mag_hq`, `mag_mq`); anything else is refused (the
    publication keeps only MAGs of quality score > 50, high or medium)."""
    q = {"High-quality": "mag_hq", "Medium-quality": "mag_mq"}.get(tier.strip())
    if q is None:
        raise Refused(f"MIMAG tier {tier!r} is neither High-quality nor Medium-quality")
    return q


def wet1_members(zip_path: Path) -> tuple[list[tuple[str, str]], dict]:
    """(member name, MAG id) of each `dereplicated_genomes/<id>.fa` in archive order, and the census of every entry."""
    try:
        with zipfile.ZipFile(zip_path) as z:
            infos = z.infolist()
    except (zipfile.BadZipFile, OSError) as exc:
        raise Refused(f"{zip_path}: not a readable zip: {exc}") from exc
    mem, other, dirs = [], [], 0
    for i in infos:
        n = i.filename
        if i.is_dir():
            dirs += 1
        elif n.startswith(WET1_DIR) and n.endswith(".fa") and "/" not in n[len(WET1_DIR):]:
            mem.append((n, n[len(WET1_DIR):-3]))
        else:
            other.append(n)
    dup = [g for g, k in Counter(g for _, g in mem).items() if k > 1]
    if dup:
        raise Refused(f"{zip_path}: two members for {sorted(dup)[0]}")
    return mem, {"entries": len(infos), "directories": dirs, "mags": len(mem), "other": sorted(other)}


def gtdb_project_crosswalk(gtdb_meta_dir: Path, bioproject: str, names: set[str]) -> tuple[dict, int]:
    """GTDB R232's genomes of one BioProject (`ncbi_bioproject`), matched to the given bin names through
    `ncbi_isolate` or `ncbi_strain_identifiers`: name -> sorted versionless accessions, and the number of R232 rows of
    that BioProject."""
    import pyarrow as pa
    import pyarrow.csv as pacsv

    cols = ["accession", "ncbi_bioproject", "ncbi_isolate", "ncbi_strain_identifiers"]
    by: dict[str, set[str]] = defaultdict(set)
    n = 0
    for key in ("ar53", "bac120"):
        p = Path(gtdb_meta_dir) / r232_index.FILES[key]
        if not p.is_file():
            raise Refused(f"missing GTDB R232 metadata {p}")
        rd = pacsv.open_csv(str(p), read_options=pacsv.ReadOptions(block_size=1 << 24),
                            parse_options=pacsv.ParseOptions(delimiter="\t"),
                            convert_options=pacsv.ConvertOptions(include_columns=cols,
                                                                 column_types={c: pa.string() for c in cols}))
        for b in rd:
            for acc, bp, iso, strain in zip(*(b.column(i).to_pylist() for i in range(4))):
                if (bp or "").strip() != bioproject:
                    continue
                n += 1
                for v in [iso or "", *(strain or "").split(",")]:
                    v = v.strip()
                    if v in names:
                        by[v].add(r232_index.nover(acc or ""))
    return {k: sorted(v) for k, v in by.items()}, n


WET1_STATUS_HEAD = ["genome", "name", "sample", "gtdb_r232_accession", "r232_status", "in_gtdb_r232",
                    "is_gtdb_representative", "gtdb_representative", "gtdb_species", "in_globdb_r232",
                    "globdb_r232_id", "gtdb_classification", "completeness", "contamination", "quality_tier",
                    "contigs", "bp"]


def stage_wet1(zip_path: Path, xlsx: Path, index_dir: Path, gtdb_meta_dir: Path, meta_out: Path, fasta_out: Path, *,
               root: Path | None, root_stages, owner_inputs: dict[str, Path], later_inputs: dict[str, Path],
               inputs_home: Path, retrieved_at: str = WET1_RETRIEVED, work: Path | None = None, threads: int = 2,
               memory_limit: str = "4GB", workers: int = 1) -> dict:
    """Stage WET1: every dereplicated MAG of the archive, as card C-77 staged single cells; returns the report."""
    zip_path, index_dir = Path(zip_path), Path(index_dir)
    meta_dir, fa_dir = Path(meta_out) / WET1, Path(fasta_out) / WET1
    final_fa = fa_dir / WET1_FASTA
    for p in (meta_dir, final_fa):
        if p.exists():
            raise Refused(f"{p} exists")
    if not (index_dir / "gtdb_r232.parquet").is_file() or not (index_dir / "globdb_r232.parquet").is_file():
        raise Refused(f"{index_dir} is not card C-74's index (no gtdb_r232.parquet or globdb_r232.parquet)")
    _check_later(later_inputs)
    meta = mag_meta(Path(xlsx))
    mem, census = wet1_members(zip_path)
    if not mem:
        raise Refused(f"{zip_path}: no MAG member ({WET1_DIR}*.fa)")
    lacking = sorted(g for _, g in mem if g not in meta)
    if lacking:
        raise Refused(f"{len(lacking)} members are not in {xlsx}; first {lacking[0]}")
    for _, g in mem:
        meta[g]["assembly_quality"] = mimag_quality(meta[g]["quality"])
    tmp = Path(meta_out) / f"{WET1}.tmp"
    partial = fa_dir / f".{WET1_FASTA}.partial"
    work_dir = Path(work) if work is not None else Path(meta_out) / f"{WET1}.work"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    fa_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        rep = _stage_wet1(zip_path, meta, mem, census, index_dir, Path(gtdb_meta_dir), tmp, partial, work_dir, root,
                          list(root_stages), owner_inputs, later_inputs, Path(inputs_home), final_fa, retrieved_at,
                          threads, memory_limit, workers)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        partial.unlink(missing_ok=True)
        if fa_dir.exists() and not any(fa_dir.iterdir()):
            fa_dir.rmdir()
        raise
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    partial.rename(final_fa)
    tmp.rename(meta_dir)
    files = {n: meta_dir / n for n in WET1_META}
    files[WET1_FASTA] = final_fa
    rep["sha256"] = {n: ids.sha256_file(p).hex() for n, p in files.items()}
    rep["bytes"] = {n: p.stat().st_size for n, p in files.items()}
    return rep


def _stage_wet1(zip_path, meta, mem, census, index_dir, gtdb_meta_dir, tmp, partial, work_dir, root, root_stages,
                owner_inputs, later_inputs, inputs_home, final_fa, retrieved_at, threads, memory_limit,
                workers) -> dict:
    hint = {g: meta[g]["classification"] for _, g in mem}
    gc_results: dict[str, dict] = {}
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    pending: list[tuple[str, object]] = []

    def submit(g: str, seqs: list[bytes]) -> None:
        if pool is None:
            gc_results[g] = c58.screen_genetic_code(seqs, hint[g], screen_max_bp=500_000)
            return
        pending.append((g, pool.submit(c58._screen_worker, (g, seqs, hint[g], 500_000))))
        while len(pending) >= 4 * workers:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]

    contigs: dict[str, tuple[bytes, list[int]]] = {}
    names: set[str] = set()
    n_renamed = 0
    w = c77.BgzfWriter(partial)
    try:
        with zipfile.ZipFile(zip_path) as z:
            for member, g in mem:
                with z.open(member) as fh:
                    recs = _read_records(fh, member)
                for hdr, seq, width in recs:
                    first = hdr.split(None, 1)[0] if hdr.strip() else ""
                    if not first.startswith(g + "_"):
                        hdr = f"{g}|{hdr}"
                        n_renamed += 1
                    first = hdr.split(None, 1)[0]
                    if first in names:
                        raise Refused(f"contig name {first} occurs twice in the stage")
                    names.add(first)
                    w.write(b"\n".join(fasta.format_record(hdr, seq, width)) + b"\n")
                if not recs:
                    raise Refused(f"{member} holds no sequence")
                contigs[g] = (b"".join(ids.contig_uid(s) for _, s, _ in recs), [len(s) for _, s, _ in recs])
                submit(g, [s for _, s, _ in recs])
        while pending:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]
    except (zipfile.BadZipFile, OSError, EOFError, zlib.error) as exc:
        raise Refused(f"{zip_path}: not a readable zip: {exc}") from exc
    finally:
        w.close()
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    kept = sorted(contigs)

    # the accession test: the table gives no accession, so GTDB R232's rows of the BioProject, by bin name
    by_name, n_project = gtdb_project_crosswalk(gtdb_meta_dir, WET1_BIOPROJECT, {meta[g]["name"] for g in kept})
    via = {g: by_name[meta[g]["name"]][0] for g in kept if by_name.get(meta[g]["name"])}
    look = {x["accession"]: x for x in r232_index.lookup(index_dir, sorted(set(via.values())))}
    status = {g: r232_status(look.get(via[g]) if g in via else None, True) for g in kept}

    resolved, genome_uids, conflicts = c58._resolve_within_stage_code_conflicts(kept, contigs, gc_results)
    records_path = tmp / "records.tsv"
    ordinal = 0
    with open(records_path, "w", encoding="utf-8", buffering=1 << 20) as rf, open(
            tmp / "genetic_codes.tsv", "w", encoding="utf-8", buffering=1 << 20) as gf:
        rf.write("ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\tname\n")
        gf.write("\t".join(c58.GC_HEAD) + "\n")
        for g in kept:
            raw, lens = contigs[g]
            gc = gc_results[g]
            gf.write("\t".join([
                g, WET1, hint[g], str(sum(lens)), str(len(lens)), meta[g]["assembly_quality"], str(resolved[g]),
                gc["reason"],
                str(gc["tested_alt"]).lower(),
                c58._fmt_opt(gc["cd_11"]), c58._fmt_opt(gc["mean_orf_11"]), c58._fmt_opt(gc["n_genes_11"]),
                c58._fmt_opt(gc["cd_25"]), c58._fmt_opt(gc["mean_orf_25"]), c58._fmt_opt(gc["n_genes_25"]),
                c58._fmt_opt(gc["cd_4"]), c58._fmt_opt(gc["mean_orf_4"]), c58._fmt_opt(gc["n_genes_4"]),
                c58._fmt_opt(gc["tga_per_kb"])]) + "\n")
            for k, n in enumerate(lens):
                rf.write(f"{ordinal}\t{g}\t{k}\t{raw[k * 16:(k + 1) * 16].hex()}\t{n}\t"
                         f"{meta[g]['assembly_quality']}\t{resolved[g]}\t{meta[g]['name']}\n")
                ordinal += 1
    owners, later = _owners(records_path, root, root_stages, owner_inputs, later_inputs, work_dir, threads,
                            memory_limit)

    listing = make_inputs.Listing()
    genomes, repeated, shared, held_later = [], [], [], []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ordinal = 0
    for g in kept:
        raw, lens = contigs[g]
        uids = [raw[i: i + 16] for i in range(0, len(raw), 16)]
        for k, (cu, n) in enumerate(zip(uids, lens)):
            contig_len.setdefault(cu, n)
            if cu in owners:
                shared.append((ordinal, g, k, cu.hex(), n, owners[cu]))
            elif cu in later:
                held_later.append((ordinal, g, k, cu.hex(), n, later[cu]))
            ordinal += 1
        e = listing.add(uids, resolved[g], "prokaryote", meta[g]["assembly_quality"])
        if e is None:
            repeated.append((ordinal - len(uids), g, genome_uids[g].hex(), first_of[genome_uids[g]]))
            continue
        genomes.append(e)
        first_of[genome_uids[g]] = g
    _write_tsv(tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    _write_tsv(tmp / "shared_records.tsv",
               ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"], shared)
    _write_tsv(tmp / "later_records.tsv", ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "holder"],
               held_later)

    def yn(v: bool) -> str:
        return "yes" if v else "no"

    rows = []
    for g in kept:
        x = look.get(via[g]) if g in via else None
        d = meta[g]
        rows.append((g, d["name"], d["sample"], via.get(g, ""), status[g], yn(bool(x and x["in_gtdb_r232"])),
                     str(x["is_representative"]) if x and x["in_gtdb_r232"] else "",
                     (x["representative"] or "") if x else "", (x["species"] or "") if x else "",
                     yn(bool(x and x["in_globdb_r232"])), (x["globdb_r232_id"] or "") if x else "",
                     d["classification"], d["completeness"], d["contamination"], d["quality"],
                     len(contigs[g][1]), sum(contigs[g][1])))
    _write_tsv(tmp / "r232_status.tsv", WET1_STATUS_HEAD, rows)

    obj = {**WET1_SOURCE, "retrieved_at": retrieved_at}
    (tmp / "objects.json").write_text(json.dumps([{"path": str(inputs_home / WET1 / WET1_FASTA), **obj}], indent=1))
    (tmp / "objects_staging.json").write_text(json.dumps([{"path": str(final_fa), **obj}], indent=1))
    (tmp / "genomes.json").write_text(json.dumps(genomes))

    own_c, own_b = _tally(contig_len, owners)
    lat_c, lat_b = _tally(contig_len, later, owners)
    bp = {g: sum(contigs[g][1]) for g in kept}
    return {
        "census": census, "mags": len(kept), "contigs": sum(len(contigs[g][1]) for g in kept),
        "bases": sum(bp.values()), "contigs_renamed": n_renamed, "table_rows": len(meta),
        "table_not_in_archive": len(set(meta) - set(kept)), "gtdb_rows_of_bioproject": n_project,
        "matched_in_gtdb_r232": len(via),
        "r232_status": dict(sorted(Counter(status.values()).items())),
        "r232_status_bases": {k: sum(bp[g] for g in kept if status[g] == k) for k in sorted(set(status.values()))},
        "genetic_codes": dict(sorted(Counter(str(resolved[g]) for g in kept).items())),
        "code_conflicts": len(conflicts),
        "genomes": len(genomes), "genomes_repeated": len(repeated), "contigs_distinct": len(contig_len),
        "bases_distinct": sum(contig_len.values()),
        "held_by": own_c, "held_bases_by": own_b, "later_held_by": lat_c, "later_held_bases_by": lat_b,
        "contigs_new": sum(1 for cu in contig_len if cu not in owners),
        "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
        "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
        "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
        **{"listing_" + k: v for k, v in listing.report().items()},
    }


# --- skani against MAG stages on the pod and C-75's r226 archaea (recorded, not cut) -------------------------------

def genus_genomes(taxonomy_tsv: Path, genera: set[str]) -> set[str]:
    """The ids of a table with `ID` and `taxonomy` columns (card C-63's non_gtdb_stages.tsv) whose GTDB genus is
    one of GENERA."""
    out: set[str] = set()
    with open(taxonomy_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if genus_of(r["taxonomy"]) in genera:
                out.add(r["ID"])
    return out


def merge_best(groups: list[dict[str, dict]]) -> dict[str, dict]:
    """The best hit of each genome over several `best_hits` results (largest ANI, then aligned fraction, then names)."""
    out: dict[str, dict] = {}
    for h in groups:
        for g, x in h.items():
            cur = out.get(g)
            if cur is None or (-x["ani"], -x["af"], x["stage"], x["genome"]) < (
                    -cur["ani"], -cur["af"], cur["stage"], cur["genome"]):
                out[g] = x
    return out


def wet1_ani(meta: Path, fasta_path: Path, out: Path, *, whole: dict[str, tuple[Path, list[Path]]],
             by_genus: dict[str, tuple[Path, list[Path]]], taxonomy_tsv: Path | None,
             archaea: tuple[Path, Path] | None, threads: int = 8) -> dict:
    """WET1's MAGs against reference groups, one skani process per group (so that no process holds every sketch):
    WHOLE staged sets (name -> (meta dir, FASTA files)), BY_GENUS staged sets cut to the genomes whose GlobDB r226
    genus (TAXONOMY_TSV) is a genus of some WET1 MAG, and C-75's r226 archaea of those genera. Recorded, not cut."""
    out = Path(out)
    work = out / "work"
    shutil.rmtree(work, ignore_errors=True)
    q_files = split_genomes(Path(meta), Path(fasta_path), work / "WET1")
    (work / "q.txt").write_text("".join(f"{p}\n" for p in sorted(map(str, q_files.values()))))
    query_of = {str(p): g for g, p in q_files.items()}
    genera = {genus_of(r["gtdb_classification"]) for r in _rows_tsv(Path(meta) / "r232_status.tsv")} - {""}
    groups: list[tuple[str, object]] = [(n, ("whole", v)) for n, v in whole.items()]
    groups += [(n, ("genus", v)) for n, v in by_genus.items()]
    if archaea is not None:
        groups.append(("r226_archaea", ("archaea", archaea)))
    keep = genus_genomes(taxonomy_tsv, genera) if by_genus and taxonomy_tsv is not None else set()
    res: dict = {"queries": len(q_files), "genera": len(genera), "threads": threads,
                 "threshold": {"ani_min": ANI_SAG, "af_min_larger_of_two": AF_SAG}, "groups": {}}
    hits_all = []
    for name, (kind, v) in groups:
        d = work / name
        unlisted: list[str] = []
        if kind == "archaea":
            files = archaea_refs(v[0], v[1], genera, d)
        else:
            files = split_genomes(v[0], v[1], d, keep=keep if kind == "genus" else None, unlisted=unlisted)
        ref_of = {str(p): (name, g) for g, p in files.items()}
        entry = {"kind": kind, "references": len(files), "records_in_no_genome": len(unlisted),
                 "records_in_no_genome_first": unlisted[:5]}
        if files:
            (work / f"r_{name}.txt").write_text("".join(f"{p}\n" for p in sorted(ref_of)))
            table = out / f"wet1_vs_{name}.tsv"
            t0 = c39_ani.child_cpu()
            subprocess.run(["skani", "dist", "-t", str(threads), "--min-af", str(AF_SAG), "--ql", str(work / "q.txt"),
                            "--rl", str(work / f"r_{name}.txt"), "-o", str(table)], check=True)
            entry["cpu_s"] = round(c39_ani.child_cpu() - t0, 1)
            h = best_hits(table, query_of, ref_of)
            entry["meet"] = len(h)
            hits_all.append(h)
        shutil.rmtree(d, ignore_errors=True)
        res["groups"][name] = entry
    best = merge_best(hits_all)
    _write_tsv(out / "wet1_ani_hits.tsv", ["genome", "stage", "ref_genome", "ani", "af"],
               [(g, h["stage"], h["genome"], h["ani"], h["af"]) for g, h in sorted(best.items())])
    res["meet"] = len(best)
    res["meet_by_stage"] = dict(sorted(Counter(h["stage"] for h in best.values()).items()))
    (out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    shutil.rmtree(work, ignore_errors=True)
    return res


# === command line ==================================================================================================

def _check(stage_name: str, meta: Path) -> dict:
    p = subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "s2_check_inputs.py"), stage_name,
                        str(meta / "genomes.json"), str(meta / "objects_staging.json")],
                       capture_output=True, text=True, check=False)
    doc = json.loads(p.stdout) if p.stdout.strip().startswith("{") else {"stdout": p.stdout}
    return {"exit": p.returncode, "problems": doc.get("problems"), "gene_call_mode": doc.get("gene_call_mode"),
            "genomes": doc.get("genomes"), "genome_input_rows": doc.get("genome_input_rows"),
            "codes": doc.get("codes"), "kinds": doc.get("kinds"), "stderr": p.stderr[-2000:]}


def _pairs(items: list[str], what: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in items:
        k, _, v = item.partition("=")
        if not k or not v:
            raise Refused(f"{what} {item!r}: expected NAME=VALUE")
        if k in out:
            raise Refused(f"{what} {k} given twice")
        out[k] = v
    return out


def _stage_args(p: argparse.ArgumentParser, retrieved: str) -> None:
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--root-stages", default=",".join(ROOT_STAGES))
    p.add_argument("--chain-home", type=Path, default=None,
                   help="the r26.11.4 stages' inputs, read in a26's chain order as owners")
    p.add_argument("--later", action="append", default=[], metavar="NAME=DIR",
                   help="an r26.11.5-chain stage before this one whose genomes.json is under DIR (in chain order)")
    p.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    p.add_argument("--retrieved-at", default=retrieved)
    p.add_argument("--work", type=Path, default=None)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--memory-limit", default="4GB")
    p.add_argument("--report-out", type=Path, default=None)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage GWV1 and SAG3 inputs (card C-81, DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gwv1-stage")
    for n in ("fna", "info", "meta_out", "fasta_out"):
        g.add_argument(n, type=Path)
    _stage_args(g, GWV1_RETRIEVED)
    ga = sub.add_parser("gwv1-ani")
    for n in ("meta", "fasta", "out"):
        ga.add_argument(n, type=Path)
    ga.add_argument("--ref", action="append", default=[], metavar="NAME=FASTA[,FASTA...]")
    ga.add_argument("--every", type=int, default=25)
    ga.add_argument("--threads", type=int, default=8)
    ga.add_argument("--chunk-bytes", type=int, default=4_000_000_000)
    sc = sub.add_parser("sag3-census")
    sc.add_argument("tar", type=Path)
    s = sub.add_parser("sag3-stage")
    for n in ("tar", "xlsx", "index_dir", "gtdb_meta_dir", "meta_out", "fasta_out"):
        s.add_argument(n, type=Path)
    _stage_args(s, SAG3_RETRIEVED)
    s.add_argument("--workers", type=int, default=1)
    sa = sub.add_parser("sag3-ani")
    for n in ("meta", "fasta", "out"):
        sa.add_argument(n, type=Path)
    sa.add_argument("--ref", action="append", default=[], metavar="NAME=META_DIR:FASTA")
    sa.add_argument("--threads", type=int, default=8)
    sa.add_argument("--archaea", default="", metavar="FASTA:LIST_TSV", help="card C-75's r226 archaea and list")
    w = sub.add_parser("wet1-stage")
    for n in ("zip", "xlsx", "index_dir", "gtdb_meta_dir", "meta_out", "fasta_out"):
        w.add_argument(n, type=Path)
    _stage_args(w, WET1_RETRIEVED)
    w.add_argument("--workers", type=int, default=1)
    wa = sub.add_parser("wet1-ani")
    for n in ("meta", "fasta", "out"):
        wa.add_argument(n, type=Path)
    wa.add_argument("--whole", action="append", default=[], metavar="NAME=META_DIR:FASTA[,FASTA]")
    wa.add_argument("--by-genus", action="append", default=[], metavar="NAME=META_DIR:FASTA[,FASTA]")
    wa.add_argument("--taxonomy", type=Path, default=None, help="ID and taxonomy per genome (non_gtdb_stages.tsv)")
    wa.add_argument("--archaea", default="", metavar="FASTA:LIST_TSV")
    wa.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "sag3-census":
            mem, census = sag3_members(a.tar)
            print(json.dumps(census, indent=1, sort_keys=True))
            return 0
        if a.cmd == "gwv1-ani":
            refs = {k: [Path(x) for x in v.split(",") if x] for k, v in _pairs(a.ref, "--ref").items()}
            print(json.dumps(gwv1_ani(a.meta, a.fasta, a.out, refs, every=a.every, threads=a.threads,
                                      chunk_bytes=a.chunk_bytes), sort_keys=True))
            return 0
        if a.cmd == "sag3-ani":
            refs = {}
            for k, v in _pairs(a.ref, "--ref").items():
                m, _, fa = v.partition(":")
                if not m or not fa:
                    raise Refused(f"--ref {k}={v}: expected NAME=META_DIR:FASTA")
                refs[k] = (Path(m), Path(fa))
            arc = None
            if a.archaea:
                fa, _, lst = a.archaea.partition(":")
                if not fa or not lst:
                    raise Refused(f"--archaea {a.archaea}: expected FASTA:LIST_TSV")
                arc = (Path(fa), Path(lst))
            print(json.dumps(sag3_ani(a.meta, a.fasta, a.out, refs, threads=a.threads, archaea=arc), sort_keys=True))
            return 0
        if a.cmd == "wet1-ani":
            def sets(items, what):
                o = {}
                for k, v in _pairs(items, what).items():
                    m, _, fas = v.partition(":")
                    if not m or not fas:
                        raise Refused(f"{what} {k}={v}: expected NAME=META_DIR:FASTA[,FASTA]")
                    o[k] = (Path(m), [Path(x) for x in fas.split(",") if x])
                return o
            arc = None
            if a.archaea:
                fa, _, lst = a.archaea.partition(":")
                if not fa or not lst:
                    raise Refused(f"--archaea {a.archaea}: expected FASTA:LIST_TSV")
                arc = (Path(fa), Path(lst))
            print(json.dumps(wet1_ani(a.meta, a.fasta, a.out, whole=sets(a.whole, "--whole"),
                                      by_genus=sets(a.by_genus, "--by-genus"), taxonomy_tsv=a.taxonomy, archaea=arc,
                                      threads=a.threads), sort_keys=True))
            return 0
        owners = {st: a.chain_home / st for st, _ in a26.R26_11_4_CHAIN} if a.chain_home is not None else {}
        later = {k: Path(v) for k, v in _pairs(a.later, "--later").items()}
        roots = [x for x in a.root_stages.split(",") if x]
        if a.cmd == "gwv1-stage":
            name = GWV1
            rep = stage_gwv1(a.fna, a.info, a.meta_out, a.fasta_out, root=a.root, root_stages=roots,
                             owner_inputs=owners, later_inputs=later, inputs_home=a.inputs_home,
                             retrieved_at=a.retrieved_at, work=a.work, threads=a.threads, memory_limit=a.memory_limit)
        elif a.cmd == "wet1-stage":
            name = WET1
            rep = stage_wet1(a.zip, a.xlsx, a.index_dir, a.gtdb_meta_dir, a.meta_out, a.fasta_out, root=a.root,
                             root_stages=roots, owner_inputs=owners, later_inputs=later, inputs_home=a.inputs_home,
                             retrieved_at=a.retrieved_at, work=a.work, threads=a.threads,
                             memory_limit=a.memory_limit, workers=a.workers)
        else:
            name = SAG3
            rep = stage_sag3(a.tar, a.xlsx, a.index_dir, a.gtdb_meta_dir, a.meta_out, a.fasta_out, root=a.root,
                             root_stages=roots, owner_inputs=owners, later_inputs=later, inputs_home=a.inputs_home,
                             retrieved_at=a.retrieved_at, work=a.work, threads=a.threads,
                             memory_limit=a.memory_limit, workers=a.workers)
        rep["owner_stages"] = [*(roots if a.root else []), *owners]
        rep["later_stages"] = [*later]
        rep["s2_check"] = _check(name, a.meta_out / name)
        txt = json.dumps(rep, indent=1, sort_keys=True)
        if a.report_out is not None:
            a.report_out.parent.mkdir(parents=True, exist_ok=True)
            a.report_out.write_text(txt + "\n")
        print(txt)
        return 0 if rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [] else 1
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())

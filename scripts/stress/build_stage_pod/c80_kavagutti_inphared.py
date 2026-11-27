"""Card C-80 (DEC-79, DEC-80): stage Kavagutti's HiFi MAGs and NCLDVs as `KAVHIFI` and INPHARED's 7Apr2026 delta as
`INPH2604`.

    python c80_kavagutti_inphared.py kav  SUPPORT_TGZ META_OUT FASTA_OUT [--root ROOT] [--root-stages ...]
        [--chain-home DIR] [--later NAME=DIR ...] [--index INDEX_DIR] [--inputs-home DIR] [--workers 2]
        [--work DIR] [--report-out F]
    python c80_kavagutti_inphared.py inph GENOMES_FA META_OUT FASTA_OUT [--root ROOT] [--root-stages ...]
        [--chain-home DIR] [--later NAME=DIR ...] [--known-root-stages S1] [--known-input NAME=DIR ...]
        [--inputs-home DIR] [--work DIR] [--report-out F]

`kav` reads `Support-data.gz` (Figshare 20431332, `kavagutti_freshwater_hifi_phages`) once, as a stream:
- `bac-MAGs-fnas.tar.gz`: one genome per FASTA file, `kind="prokaryote"`, quality `contig` (the archive has no
  completeness table, as GV1), linear (no topology is given); the genetic code of each (11, 25 or 4) is screened with
  card C-66's `screen_genetic_code` (pyrodigal on the longest contigs); a MAG with a GCA_/GCF_ accession in its file
  name or headers is looked up in the C-74 index (`r232_index.lookup`) and not staged when GTDB R232 or GlobDB r232
  holds it. No ANI test is made.
- `NCLDV-data.tar.gz` (`NCLDV/NCLDV-fna/*.fna`): one genome per file, `kind="virus"` (called by gv-meta in any case),
  DEC-65 (1) exact 20..1,000 bp terminal repeat trimmed and staged circular, as card C-67 did for the phages.
- `phages-fna.tar.gz` was taken by `LRV1` (card C-67). `euk-MAGs-fna.tar.gz` (card C-60) and
  `mitochondrias-fna.tar.gz` are copied byte for byte to `deferred/`, and listed in `euk_deferred.tsv` and
  `mito_deferred.tsv`: the archive carries no per-record genetic code for the mitochondria (only file names), so they
  are not staged (card C-59a's `O2` routes by the RefSeq record's `transl_table`) and wait for a later organelle stage.
Ownership: `shared_records.tsv` the contigs the root's stored stages or the stages of OWNER_INPUTS hold (the
`r26.11.4` chain, in its order); `later_records.tsv` those that only a stage of `--later` (the stages of `r26.11.5`
before this one) holds. The stage keeps them (DEC-26: the first holder owns the contig's calls).

`inph` reads INPHARED's `7Apr2026_genomes.fa.gz`. A record is staged when its accession (first header token, version
dropped) is in neither `S1` (the root's `contig_label`) nor `SAT1` (its `records.tsv`), and its `contig_uid` is held by
neither `S1` nor `SAT1`. Records are byte-identical to the archive (DEC-65 (2)), `kind="virus"`, linear, quality
`complete`. The four cells of (accession known, sequence held) are counted; `skipped_records.tsv` lists the two
crossed ones.
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c48c_virus_inputs as v2  # noqa: E402
import c62_metavr_inputs as c62  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

ROOT_STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
ACC_RE = re.compile(r"GC[AF]_\d{9}")
FASTA_EXTS = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")

KAV_BAC = "figshare/bac-MAGs-fnas.tar.gz"
KAV_NCLDV = "figshare/NCLDV-data.tar.gz"
KAV_MITO = "figshare/mitochondrias-fna.tar.gz"
KAV_EUK = "figshare/euk-MAGs-fna.tar.gz"
KAV_FASTA = {"bac_mag": "kavagutti_hifi_bacterial_mags.fna", "ncldv": "kavagutti_hifi_ncldv.fna"}
KAV_SOURCE = {
    "source_name": "kavagutti_freshwater_hifi_phages",
    "source_release": "Kavagutti et al. 2023, Figshare 20431332 (doi:10.6084/m9.figshare.20431332.v1) "
                      "Support-data.gz; KAVHIFI selection by c80_kavagutti_inphared.py (DEC-65, DEC-79)",
    "license": "CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/, Figshare API v2, as sources.yaml records it)",
    "redistributable": True,
}
INPH_FASTA = "inphared_7apr2026_delta.fna"
INPH_SOURCE = {
    "source_name": "inphared_7apr2026",
    "source_release": "INPHARED 7Apr2026 (7Apr2026_genomes.fa.gz) delta against S1 (14Apr2025) and SAT1 "
                      "(9Aug2025 delta); INPH2604 selection by c80_kavagutti_inphared.py (DEC-65 (2), DEC-79)",
    "license": "INSDC (public); INPHARED scripts AGPL-3.0, as the inphared_9aug2025 entry records",
    "redistributable": True,
}

RECORDS_HEAD = ["ordinal", "source", "category", "genome", "accession", "contig_uid", "length", "genetic_code",
                "circular", "dtr_class", "source_contig_uid", "source_length", "trimmed_bases", "owner_stage",
                "later_holder"]
META_FILES = ("genomes.json", "objects.json", "objects_staging.json", "records.tsv", "repeated_records.tsv",
              "shared_records.tsv", "later_records.tsv", "trims.tsv")


def base_acc(acc: str) -> str:
    acc = acc.strip()
    return acc.split(".", 1)[0] if "." in acc else acc


def _is_fasta_name(name: str) -> bool:
    base = os.path.basename(name.lower())
    return not base.startswith(".") and base.endswith(FASTA_EXTS)


def _stem(name: str) -> str:
    base = os.path.basename(name)
    for ext in (".gz", ".fasta", ".fna", ".fa"):
        if base.lower().endswith(ext):
            base = base[: -len(ext)]
    return base


def _meta_parts(tmp: Path, final_dir: Path) -> dict[Path, Path]:
    """Every file under TMP (the scratch directory `_work` aside) at its place under FINAL_DIR."""
    return {final_dir / p.relative_to(tmp): p for p in sorted(tmp.rglob("*"))
            if p.is_file() and "_work" not in p.relative_to(tmp).parts}


def _hashes(finals: list[Path]) -> tuple[dict, dict]:
    ps = sorted(finals, key=lambda p: p.name)
    return {p.name: ids.sha256_file(p).hex() for p in ps}, {p.name: p.stat().st_size for p in ps}


def _commit(parts: dict[Path, Path]) -> None:
    """Rename each partial file onto its final path (all refused first if one exists)."""
    for final in parts:
        if final.exists():
            raise Refused(f"{final} exists")
    for final, part in parts.items():
        final.parent.mkdir(parents=True, exist_ok=True)
        os.replace(part, final)


# --- genetic code screen ------------------------------------------------------------------------------------------

def screen_codes(seqs: list[bytes]) -> dict:
    """The genetic code (11, 25 or 4) of a prokaryotic genome: card C-66's screen, no taxonomy hint."""
    import c66_rt_hotspot_compendia as c66

    return c66.screen_genetic_code(seqs, "")


def _screen_job(args: tuple[str, list[bytes]]) -> tuple[str, dict]:
    name, seqs = args
    return name, screen_codes(seqs)


# --- the accession test (card C-74's index) -----------------------------------------------------------------------

def accessions_in(*texts: str) -> list[str]:
    """The distinct versionless GCA_/GCF_ accessions in the texts, in order of appearance."""
    out: dict[str, None] = {}
    for t in texts:
        for m in ACC_RE.finditer(t):
            out.setdefault(m.group(0), None)
    return list(out)


def accession_lookup(index_dir: Path | str):
    """A function (accessions) -> rows of `r232_index.lookup` for the index at INDEX_DIR."""
    import r232_index

    def look(accs: list[str]) -> list[dict]:
        return r232_index.lookup(index_dir, accs)

    return look


# --- ownership ----------------------------------------------------------------------------------------------------

def _owners(all_uids: set[bytes], root, root_stages, owner_inputs, later_inputs, work: Path, threads: int,
            memory_limit: str):
    work.mkdir(parents=True, exist_ok=True)
    f = work / "all_uids.tsv"
    with open(f, "w") as fh:
        fh.write("contig_uid\n")
        for cu in sorted(all_uids):
            fh.write(cu.hex() + "\n")
    owners = c62._owners_from_records(f, root, list(root_stages), owner_inputs, work=work, threads=threads,
                                      memory_limit=memory_limit)
    later = c62._owners_from_records(f, None, [], later_inputs, work=work, threads=threads,
                                     memory_limit=memory_limit) if later_inputs else {}
    return owners, later


def known_accessions(root, stages, record_dirs: dict[str, Path], *, work: Path, threads: int = 2,
                     memory_limit: str = "2GB") -> set[str]:
    """Versionless accessions of the root stages' `contig_label` (`accession`) and of the `records.tsv` of each
    staged input (`accession`)."""
    import duckdb

    from biobox import io as bbio
    from biobox.l0seq import store

    work.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(config={"threads": int(threads), "memory_limit": memory_limit,
                                 "temp_directory": str(work / "duck_acc")})
    out: set[str] = set()
    try:
        for s in stages:
            files = bbio.part_files(Path(root) / store.LABEL_DIR / f"stage={s}") if root is not None else []
            if root is not None and not files:
                raise Refused(f"stage {s}: no contig_label under {root}")
            if files:
                lit = "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for p in files) + "]"
                out |= {base_acc(a) for (a,) in con.execute(
                    f"SELECT DISTINCT accession FROM read_parquet({lit}, hive_partitioning = false) "
                    "WHERE accession IS NOT NULL").fetchall()}
        for name, d in record_dirs.items():
            rec = Path(d) / "records.tsv"
            if not rec.is_file():
                raise Refused(f"{name}: no records.tsv in {d}")
            lit = "'" + str(rec).replace("'", "''") + "'"
            out |= {base_acc(a) for (a,) in con.execute(
                f"SELECT DISTINCT accession FROM read_csv({lit}, delim='\t', header=true, quote='', "
                "all_varchar=true) WHERE accession IS NOT NULL AND accession <> ''").fetchall()}
    finally:
        con.close()
    return out


# --- KAVHIFI ------------------------------------------------------------------------------------------------------

def _dec65(seq: bytes) -> tuple[bytes, int, str, bool]:
    """DEC-65 (1) on a metagenome-assembled record: an exact 20..1,000 bp terminal repeat is trimmed once and the
    record staged circular; every other class stays whole and linear."""
    e = v2.ends(seq)
    cls = v2.classify(e)
    if cls == "exact":
        return v2.trim(seq, e.exact), e.exact, cls, True
    return seq, 0, cls, False


def _read_member_fasta(fobj) -> list[tuple[str, bytes, int]]:
    return [(r.header, r.seq, r.line_width or 70) for r in v2.iter_fasta_stream(fobj)]


def _copy_member(outer: tarfile.TarFile, member: tarfile.TarInfo, dst: Path) -> tuple[Path, str, int]:
    if dst.exists():
        raise Refused(f"{dst} exists")
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_name(f".{dst.name}.partial")
    src = outer.extractfile(member)
    assert src is not None
    with open(part, "wb") as w:
        shutil.copyfileobj(src, w, 1 << 22)
    os.replace(part, dst)
    return dst, ids.sha256_file(dst).hex(), dst.stat().st_size


def _list_archive(path: Path) -> list[tuple[str, int, int]]:
    """(file name, records, bases) of each FASTA file of a tar.gz."""
    out = []
    with tarfile.open(path, "r|gz") as tf:
        for m in tf:
            if not m.isfile() or not _is_fasta_name(m.name):
                continue
            f = tf.extractfile(m)
            assert f is not None
            n = bp = 0
            for r in v2.iter_fasta_stream(f):
                n += 1
                bp += len(r.seq)
            out.append((m.name, n, bp))
    return sorted(out)


def stage_kavagutti(support: Path, meta_out: Path, fasta_out: Path, *, root: Path | None, root_stages,
                    owner_inputs: dict[str, Path], later_inputs: dict[str, Path], inputs_home: Path,
                    lookup=None, retrieved_at: str = "2026-11-24T07:08:00Z", work: Path | None = None,
                    threads: int = 2, workers: int = 2, memory_limit: str = "4GB", screen=None) -> dict:
    """Stream SUPPORT once: stage the bacterial MAGs and the NCLDV files, copy the eukaryotic MAGs and the
    mitochondria aside. `lookup(accessions)` is the accession test; `screen(seqs)` the genetic-code screen."""
    stage = "KAVHIFI"
    meta_dir, fasta_dir = Path(meta_out) / stage, Path(fasta_out) / stage
    if meta_dir.exists() and any((meta_dir / n).exists() for n in META_FILES):
        raise Refused(f"{meta_dir} holds a staged KAVHIFI")
    for n in KAV_FASTA.values():
        if (fasta_dir / n).exists():
            raise Refused(f"{fasta_dir / n} exists")
    meta_tmp = Path(meta_out) / f"{stage}.tmp"
    shutil.rmtree(meta_tmp, ignore_errors=True)
    meta_tmp.mkdir(parents=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(work) if work is not None else meta_tmp / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    partial = {k: fasta_dir / f".{n}.partial" for k, n in KAV_FASTA.items()}
    created: list[Path] = []
    try:
        rep = _stage_kav(Path(support), meta_tmp, fasta_dir, partial, work_dir, root, list(root_stages), owner_inputs,
                         later_inputs, Path(inputs_home), lookup, retrieved_at, threads, workers, memory_limit,
                         screen or screen_codes, created)
    except BaseException:
        shutil.rmtree(meta_tmp, ignore_errors=True)
        for p in [*partial.values(), *created]:
            p.unlink(missing_ok=True)
        raise
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    parts = _meta_parts(meta_tmp, meta_dir)
    parts.update({fasta_dir / KAV_FASTA[k]: p for k, p in partial.items()})
    _commit(parts)
    shutil.rmtree(meta_tmp, ignore_errors=True)
    if "mito" in rep["deferred"]:
        rep["deferred"]["mito"]["path"] = str(meta_dir / "deferred" / os.path.basename(KAV_MITO))
    rep["sha256"], rep["bytes"] = _hashes([*parts])
    return rep


def _stage_kav(support, tmp, fasta_dir, partial, work_dir, root, root_stages, owner_inputs, later_inputs,
               inputs_home, lookup, retrieved_at, threads, workers, memory_limit, screen, created) -> dict:
    # genome: dict(name, category, kind, code, reason, circular, contigs=[(cu, n, src_n, acc, dtr, src_cu, k)])
    genomes: list[dict] = []
    skipped: list[tuple] = []
    trims: list[tuple] = []
    dtr_classes: Counter = Counter()
    code_rows: list[tuple] = []
    acc_rows: list[tuple] = []
    seen_members = Counter()
    deferred: dict[str, dict] = {}
    writers = {k: open(p, "wb", buffering=1 << 22) for k, p in partial.items()}
    pool = concurrent.futures.ProcessPoolExecutor(max_workers=workers) if workers > 0 and screen is screen_codes \
        else None
    pending: collections.deque = collections.deque()

    def settle(g: dict, fut) -> None:
        name, res = fut.result() if hasattr(fut, "result") else fut
        g["code"], g["reason"] = int(res["genetic_code"]), res.get("reason", "")
        code_rows.append((name, g["code"], g["reason"], res.get("cd_11", ""), res.get("mean_orf_11", ""),
                          res.get("cd_25", ""), res.get("cd_4", "")))

    def bounded(limit: int) -> None:
        while len(pending) > limit:
            settle(*pending.popleft())

    try:
        with tarfile.open(support, "r|gz") as outer:
            for m in outer:
                if not m.isfile():
                    continue
                seen_members[m.name] += 1
                if m.name == KAV_BAC:
                    with tarfile.open(fileobj=outer.extractfile(m), mode="r|gz") as inner:
                        for im in inner:
                            if not im.isfile() or not _is_fasta_name(im.name):
                                continue
                            name = _stem(im.name)
                            f = inner.extractfile(im)
                            recs = [x for x in _read_member_fasta(f) if x[1] and VALID_DNA_RE.match(x[1])]
                            if not recs:
                                skipped.append((name, "bac_mag", "no_valid_dna_record", ""))
                                continue
                            accs = accessions_in(name, *(h for h, _, _ in recs))
                            hit = None
                            if accs and lookup is not None:
                                for row in lookup(accs):
                                    acc_rows.append((name, row["accession"], row["in_gtdb_r232"],
                                                     row["in_globdb_r232"]))
                                    if row["in_gtdb_r232"] or row["in_globdb_r232"]:
                                        hit = hit or row["accession"]
                            elif not accs:
                                acc_rows.append((name, "", "", ""))
                            if hit:
                                skipped.append((name, "bac_mag", "accession_in_gtdb_r232_or_globdb_r232", hit))
                                continue
                            contigs = []
                            for h, seq, width in recs:
                                writers["bac_mag"].write(b"\n".join(fasta.format_record(h, seq, width)) + b"\n")
                                cu = ids.contig_uid(seq)
                                contigs.append((cu, len(seq), len(seq), h.split(None, 1)[0] if h else "", "", cu, 0))
                            g = {"name": name, "category": "bac_mag", "kind": "prokaryote", "quality": "contig",
                                 "circular": [], "contigs": contigs, "code": 11, "reason": ""}
                            genomes.append(g)
                            seqs = [s for _, s, _ in recs]
                            if pool is not None:
                                pending.append((g, pool.submit(_screen_job, (name, seqs))))
                                bounded(2 * workers)
                            else:
                                pending.append((g, (name, screen(seqs))))
                                bounded(0)
                elif m.name == KAV_NCLDV:
                    with tarfile.open(fileobj=outer.extractfile(m), mode="r|gz") as inner:
                        for im in inner:
                            if not im.isfile() or not _is_fasta_name(im.name) or "/NCLDV-fna/" not in im.name:
                                continue
                            name = _stem(im.name)
                            recs = [x for x in _read_member_fasta(inner.extractfile(im))
                                    if x[1] and VALID_DNA_RE.match(x[1])]
                            if not recs:
                                skipped.append((name, "ncldv", "no_valid_dna_record", ""))
                                continue
                            contigs, circ = [], []
                            for h, seq, width in recs:
                                src_cu = ids.contig_uid(seq)
                                out, k, cls, is_c = _dec65(seq)
                                dtr_classes[cls] += 1
                                cu = ids.contig_uid(out) if k else src_cu
                                writers["ncldv"].write(b"\n".join(fasta.format_record(h, out, width)) + b"\n")
                                acc = h.split(None, 1)[0] if h else ""
                                if k:
                                    trims.append(("kavagutti_ncldv", name, acc, src_cu.hex(), cu.hex(), len(seq), k,
                                                  k, len(out)))
                                contigs.append((cu, len(out), len(seq), acc, cls, src_cu, k))
                                if is_c:
                                    circ.append(cu)
                            genomes.append({"name": name, "category": "ncldv", "kind": "virus", "quality": "contig",
                                            "circular": circ, "contigs": contigs, "code": 11,
                                            "reason": "gv-meta decides"})
                elif m.name in (KAV_EUK, KAV_MITO):
                    short = "euk" if m.name == KAV_EUK else "mito"
                    dst = (fasta_dir / "deferred" if m.name == KAV_EUK else tmp / "deferred") / os.path.basename(m.name)
                    created.append(dst)
                    deferred[short] = dict(zip(("path", "sha256", "bytes"), _copy_member(outer, m, dst)))
            bounded(0)
    finally:
        for w in writers.values():
            w.close()
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    for need in (KAV_BAC, KAV_NCLDV):
        if not seen_members[need]:
            raise Refused(f"{need} not found in {support}")
    for short, d in deferred.items():
        rows = _list_archive(Path(d["path"]))
        make_inputs._write_tsv(tmp / f"{short}_deferred.tsv", ["file", "records", "bases"], rows)
        d["files"], d["records"], d["bases"] = len(rows), sum(r[1] for r in rows), sum(r[2] for r in rows)

    genomes.sort(key=lambda g: (g["category"], g["name"]))
    all_uids = {c[0] for g in genomes for c in g["contigs"]}
    owners, later = _owners(all_uids, root, root_stages, owner_inputs, later_inputs, work_dir, threads, memory_limit)
    make_inputs._write_tsv(tmp / "genetic_codes.tsv",
                           ["genome", "genetic_code", "reason", "cd_11", "mean_orf_11", "cd_25", "cd_4"],
                           sorted(code_rows))
    make_inputs._write_tsv(tmp / "mag_accessions.tsv", ["genome", "accession", "in_gtdb_r232", "in_globdb_r232"],
                           sorted(acc_rows))
    make_inputs._write_tsv(tmp / "skipped_members.tsv", ["genome", "category", "reason", "detail"], sorted(skipped))
    rep = _write_stage_meta("KAVHIFI", genomes, trims, owners, later, tmp, inputs_home, fasta_dir, retrieved_at,
                            KAV_SOURCE, KAV_FASTA, "kavagutti_freshwater_hifi_phages")
    rep["dtr_classes"] = dict(sorted(dtr_classes.items()))
    rep["members_skipped"] = Counter(s[2] for s in skipped)
    rep["genetic_codes"] = dict(sorted(Counter(g["code"] for g in genomes if g["category"] == "bac_mag").items()))
    rep["accession_test"] = {
        "mags_staged_or_tested": sum(1 for g in genomes if g["category"] == "bac_mag") + sum(
            1 for s in skipped if s[1] == "bac_mag" and s[2] != "no_valid_dna_record"),
        "mags_with_an_insdc_accession": len({r[0] for r in acc_rows if r[1]}),
        "mags_without_one_not_visible_to_the_index": len({r[0] for r in acc_rows if not r[1]}),
        "found_in_gtdb_r232": len({r[0] for r in acc_rows if r[2] is True}),
        "found_in_globdb_r232": len({r[0] for r in acc_rows if r[3] is True}),
    }
    rep["deferred"] = deferred
    return rep


def _write_stage_meta(stage, genomes, trims, owners, later, out: Path, inputs_home: Path, fasta_dir: Path,
                      retrieved_at: str, source: dict, fasta_names: dict[str, str], source_label: str) -> dict:
    """records.tsv, genomes.json, repeated/shared/later lists, trims, objects (one per category's FASTA)."""
    listing = make_inputs.Listing()
    out_genomes, repeated, shared, held_later, rows = [], [], [], [], []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    contig_cat: dict[bytes, str] = {}
    ordinal = 0
    for g in genomes:
        uids = [c[0] for c in g["contigs"]]
        circ = set(g["circular"])
        for cu, n, src_n, acc, cls, src_cu, k in g["contigs"]:
            own, lat = owners.get(cu, ""), ("" if cu in owners else later.get(cu, ""))
            rows.append((ordinal, source_label, g["category"], g["name"], acc, cu.hex(), n, g["code"],
                         str(cu in circ).lower(), cls, src_cu.hex(), src_n, k, own, lat))
            if own:
                shared.append((ordinal, g["name"], cu.hex(), own))
            elif lat:
                held_later.append((ordinal, g["name"], cu.hex(), lat))
            contig_len.setdefault(cu, n)
            contig_cat.setdefault(cu, g["category"])
            ordinal += 1
        entry = listing.add(uids, g["code"], g["kind"], g["quality"])
        guid = ids.genome_uid(list(dict.fromkeys(uids)))
        if entry is None:
            repeated.append((ordinal - len(uids), g["name"], guid.hex(), first_of[guid]))
            continue
        entry["circular"] = [u.hex() for u in g["circular"]]
        out_genomes.append(entry)
        first_of[guid] = g["name"]
    make_inputs._write_tsv(out / "records.tsv", RECORDS_HEAD, rows)
    make_inputs._write_tsv(out / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(out / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"], shared)
    make_inputs._write_tsv(out / "later_records.tsv", ["ordinal", "genome", "contig_uid", "holder"], held_later)
    make_inputs._write_tsv(out / "trims.tsv", c62.TRIM_HEAD, sorted(trims, key=lambda t: (t[0], t[1], t[2])))
    obj = {**source, "retrieved_at": retrieved_at}
    used = [k for k in fasta_names if any(g["category"] == k for g in genomes)]
    (out / "objects.json").write_text(json.dumps(
        [{"path": str(inputs_home / stage / fasta_names[k]), **obj} for k in used], indent=1))
    (out / "objects_staging.json").write_text(json.dumps(
        [{"path": str(fasta_dir / fasta_names[k]), **obj} for k in used], indent=1))
    (out / "genomes.json").write_text(json.dumps(out_genomes))

    own_c, own_b, lat_c, lat_b = Counter(), Counter(), Counter(), Counter()
    for cu, n in contig_len.items():
        if cu in owners:
            own_c[owners[cu]] += 1
            own_b[owners[cu]] += n
        elif cu in later:
            lat_c[later[cu]] += 1
            lat_b[later[cu]] += n
    by_cat: dict[str, dict] = {}
    for g in genomes:
        e = by_cat.setdefault(g["category"], {"genomes": 0, "records": 0, "bases": 0})
        e["genomes"] += 1
        e["records"] += len(g["contigs"])
        e["bases"] += sum(c[1] for c in g["contigs"])
    return {
        "stage": stage, "by_category": dict(sorted(by_cat.items())),
        "records": len(rows), "genomes": len(out_genomes), "repeated_genomes": len(repeated),
        "contigs": len(contig_len), "contig_bases": sum(contig_len.values()),
        "distinct_by_category": {c: {"contigs": sum(1 for cu in contig_len if contig_cat[cu] == c),
                                     "bases": sum(n for cu, n in contig_len.items() if contig_cat[cu] == c)}
                                 for c in sorted(set(contig_cat.values()))},
        "circular_genomes": sum(1 for g in out_genomes if g["circular"]),
        "trims": len(trims), "trimmed_bases": sum(t[7] for t in trims),
        "shared_records": len(shared), "later_records": len(held_later),
        "owner_contigs": dict(sorted(own_c.items())), "owner_bases": dict(sorted(own_b.items())),
        "later_contigs": dict(sorted(lat_c.items())), "later_bases": dict(sorted(lat_b.items())),
        "new_contigs": sum(1 for cu in contig_len if cu not in owners),
        "new_bases": sum(n for cu, n in contig_len.items() if cu not in owners),
        "new_after_later_contigs": sum(1 for cu in contig_len if cu not in owners and cu not in later),
        "new_after_later_bases": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
        **listing.report(),
    }


# --- INPH2604 -----------------------------------------------------------------------------------------------------

JUMBO_BP = 200_000


def stage_inphared(genomes_fa: Path, meta_out: Path, fasta_out: Path, *, root: Path | None, root_stages,
                   owner_inputs: dict[str, Path], later_inputs: dict[str, Path], known_root_stages,
                   known_inputs: dict[str, Path], inputs_home: Path, retrieved_at: str = "2026-04-07T00:00:00Z",
                   work: Path | None = None, threads: int = 2, memory_limit: str = "4GB") -> dict:
    stage = "INPH2604"
    meta_dir, fasta_dir = Path(meta_out) / stage, Path(fasta_out) / stage
    if (fasta_dir / INPH_FASTA).exists() or (meta_dir / "genomes.json").exists():
        raise Refused(f"{stage} is already staged under {meta_dir} / {fasta_dir}")
    meta_tmp = Path(meta_out) / f"{stage}.tmp"
    shutil.rmtree(meta_tmp, ignore_errors=True)
    meta_tmp.mkdir(parents=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(work) if work is not None else meta_tmp / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    partial = fasta_dir / f".{INPH_FASTA}.partial"
    try:
        rep = _stage_inph(Path(genomes_fa), meta_tmp, fasta_dir, partial, work_dir, root, list(root_stages),
                          owner_inputs, later_inputs, list(known_root_stages), known_inputs, Path(inputs_home),
                          retrieved_at, threads, memory_limit)
    except BaseException:
        shutil.rmtree(meta_tmp, ignore_errors=True)
        partial.unlink(missing_ok=True)
        raise
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    parts = _meta_parts(meta_tmp, meta_dir)
    parts[fasta_dir / INPH_FASTA] = partial
    _commit(parts)
    shutil.rmtree(meta_tmp, ignore_errors=True)
    rep["sha256"], rep["bytes"] = _hashes([*parts])
    return rep


def _open_fa(path: Path):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")


def _stage_inph(genomes_fa, tmp, fasta_dir, partial, work_dir, root, root_stages, owner_inputs, later_inputs,
                known_root_stages, known_inputs, inputs_home, retrieved_at, threads, memory_limit) -> dict:
    known = known_accessions(root, known_root_stages, known_inputs, work=work_dir, threads=threads,
                             memory_limit="2GB")
    recs: list[tuple[str, bytes, int, str, fasta.Record | None]] = []  # (acc, uid, length, header, record if unknown)
    with _open_fa(genomes_fa) as fh:
        for r in v2.iter_fasta_stream(fh):
            if not r.seq or not VALID_DNA_RE.match(r.seq):
                recs.append((base_acc(r.accession), b"", 0, r.header, None))
                continue
            acc = base_acc(r.accession)
            recs.append((acc, ids.contig_uid(r.seq), len(r.seq), r.header, None if acc in known else r))
    valid = [x for x in recs if x[1]]
    all_uids = {x[1] for x in valid}
    work_dir.mkdir(parents=True, exist_ok=True)
    f = work_dir / "all_uids.tsv"
    with open(f, "w") as fh:
        fh.write("contig_uid\n" + "".join(u.hex() + "\n" for u in sorted(all_uids)))
    held = c62._owners_from_records(f, root, list(known_root_stages), known_inputs, work=work_dir, threads=threads,
                                    memory_limit=memory_limit)
    cells: Counter = Counter()
    staged, skipped = [], []
    for acc, uid, n, header, r in valid:
        acc_known, seq_held = r is None, uid in held
        cells[("accession_known" if acc_known else "accession_new", "sequence_held" if seq_held else "sequence_new")] += 1
        if not acc_known and not seq_held:
            staged.append((acc, uid, n, header, r))
        elif acc_known != seq_held:
            skipped.append((acc, uid.hex(), n, "accession_new_sequence_held_by_S1_or_SAT1" if seq_held
                            else "accession_known_sequence_not_held"))
    if not staged:
        raise Refused("no INPHARED record is new against the known accessions and sequences")
    genomes = []
    with open(partial, "wb", buffering=1 << 22) as w:
        for acc, uid, n, header, r in staged:
            w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width)) + b"\n")
            genomes.append({"name": acc, "category": "inphared_delta", "kind": "virus", "quality": "complete",
                            "circular": [], "code": 11, "contigs": [(uid, n, n, acc, "", uid, 0)]})
    all_staged = {g["contigs"][0][0] for g in genomes}
    owners, later = _owners(all_staged, root, root_stages, owner_inputs, later_inputs, work_dir, threads, memory_limit)
    make_inputs._write_tsv(tmp / "skipped_records.tsv", ["accession", "contig_uid", "length", "reason"],
                           sorted(skipped))
    rep = _write_stage_meta("INPH2604", genomes, [], owners, later, tmp, inputs_home, fasta_dir, retrieved_at,
                            INPH_SOURCE, {"inphared_delta": INPH_FASTA}, "inphared_7apr2026")
    desc = {a: h for a, _, _, h, _ in staged}
    rep.update({
        "archive_records": len(recs), "archive_invalid_records": len(recs) - len(valid),
        "known_accessions": len(known), "known_sequences_in_archive": len(held),
        "cells": {f"{a}+{b}": c for (a, b), c in sorted(cells.items())},
        "jumbo_ge_200kb": sum(1 for g in genomes if g["contigs"][0][1] >= JUMBO_BP),
        "jumbo_bases": sum(g["contigs"][0][1] for g in genomes if g["contigs"][0][1] >= JUMBO_BP),
        "wolbachia": sum(1 for h in desc.values() if "wolbachia" in h.lower()),
        "synthetic_construct": sum(1 for h in desc.values() if "synthetic construct" in h.lower()),
        "skipped_records": len(skipped),
    })
    return rep


# --- command line -------------------------------------------------------------------------------------------------

def _check(stage_name: str, meta: Path) -> dict:
    p = subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "s2_check_inputs.py"), stage_name,
                        str(meta / "genomes.json"), str(meta / "objects_staging.json")],
                       capture_output=True, text=True, check=False)
    doc = json.loads(p.stdout) if p.stdout.strip().startswith("{") else {"stdout": p.stdout}
    return {"exit": p.returncode, "problems": doc.get("problems"), "gene_call_mode": doc.get("gene_call_mode"),
            "genomes": doc.get("genomes"), "genome_input_rows": doc.get("genome_input_rows"),
            "kinds": doc.get("kinds"), "codes": doc.get("codes"), "stderr": p.stderr[-2000:]}


def _kv(items: list[str]) -> dict[str, Path]:
    out = {}
    for it in items:
        name, sep, path = it.partition("=")
        if not sep or not name or not path:
            raise Refused(f"expected NAME=DIR, got {it!r}")
        out[name] = Path(path)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage KAVHIFI and INPH2604 inputs (card C-80, DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for cmd in ("kav", "inph"):
        p = sub.add_parser(cmd)
        p.add_argument("source", type=Path)
        p.add_argument("meta_out", type=Path)
        p.add_argument("fasta_out", type=Path)
        p.add_argument("--root", type=Path, default=None)
        p.add_argument("--root-stages", default=",".join(ROOT_STAGES))
        p.add_argument("--chain-home", type=Path, default=None,
                       help="the r26.11.4 stages' inputs, read in a26's chain order as owners")
        p.add_argument("--later", action="append", default=[], help="NAME=DIR of an r26.11.5 stage before this one")
        p.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
        p.add_argument("--work", type=Path, default=None)
        p.add_argument("--threads", type=int, default=2)
        p.add_argument("--memory-limit", default="4GB")
        p.add_argument("--report-out", type=Path, default=None)
        if cmd == "kav":
            p.add_argument("--index", type=Path, default=None, help="card C-74's index directory")
            p.add_argument("--workers", type=int, default=2)
        else:
            p.add_argument("--known-root-stages", default="S1")
            p.add_argument("--known-input", action="append", default=[], help="NAME=DIR with a records.tsv")
    a = ap.parse_args(argv)
    try:
        owners = {s: a.chain_home / s for s, _ in a26.R26_11_4_CHAIN} if a.chain_home is not None else {}
        later = _kv(a.later)
        root_stages = [s for s in a.root_stages.split(",") if s]
        if a.cmd == "kav":
            lookup = accession_lookup(a.index) if a.index is not None else None
            rep = stage_kavagutti(a.source, a.meta_out, a.fasta_out, root=a.root, root_stages=root_stages,
                                  owner_inputs=owners, later_inputs=later, inputs_home=a.inputs_home, lookup=lookup,
                                  work=a.work, threads=a.threads, workers=a.workers, memory_limit=a.memory_limit)
            stage_name = "KAVHIFI"
        else:
            rep = stage_inphared(a.source, a.meta_out, a.fasta_out, root=a.root, root_stages=root_stages,
                                 owner_inputs=owners, later_inputs=later,
                                 known_root_stages=[s for s in a.known_root_stages.split(",") if s],
                                 known_inputs=_kv(a.known_input), inputs_home=a.inputs_home, work=a.work,
                                 threads=a.threads, memory_limit=a.memory_limit)
            stage_name = "INPH2604"
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    rep["owner_stages"] = [*(root_stages if a.root else []), *owners]
    rep["later_stages"] = [*later]
    rep["s2_check"] = _check(stage_name, a.meta_out / stage_name)
    txt = json.dumps(rep, indent=1, sort_keys=True, default=str)
    if a.report_out is not None:
        a.report_out.parent.mkdir(parents=True, exist_ok=True)
        a.report_out.write_text(txt + "\n")
    print(txt)
    return 0 if rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [] else 1


if __name__ == "__main__":
    sys.exit(main())

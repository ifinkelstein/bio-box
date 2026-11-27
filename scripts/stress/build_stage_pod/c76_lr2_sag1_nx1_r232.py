"""Card C-76 (DEC-79 (3)): the R232 / r232 status of the staged `LR2`, `SAG1` and `NX1`.

    c76_lr2_sag1_nx1_r232.py run STAGE META_DIR FASTA INDEX_DIR OUT_META_DIR
                             [--names ENA_TABLE] [--root ROOT --chain-home DIR] [--later STAGE=DIR ...]
                             [--no-check] [--work DIR] [--report-out FILE]

Cards C-66 and C-58 staged the three stages against GTDB R226 and GlobDB r226. DEC-79 (3) makes GTDB R232 and
GlobDB r232 the bases of the "absent from" tests of `r26.11.5`; the orchestrating session then decided (2026-11-26)
that a genome r232 also lists is NOT removed: it is in hand, staged and verified, and the later r232 cut is
deduplicated against what `r26.11.5` holds. So this script removes nothing and writes nothing but reports. For every
genome of the stage it records its r232 status: the INSDC accession is read from the genome id
(`ehomd:GCA_000406805.1`) or, for single cells, from an ENA assembly table by the assembly's name (`--names`); card
C-74's `lookup` (`r232_index.py`) says whether GTDB R232 or GlobDB r232 holds it. A genome with no accession (a MAG
named by its submitter, an unbinned contig) is `not_visible`: the accession test cannot see it (no r232 genome is on
the pod for an ANI test).

Written under OUT_META_DIR/STAGE: `r232_status.tsv` (one row per genome, sorted by genome id, deterministic) and
`r232_report.json` (counts, the stage's file table with sha256, exact `contig_uid` ownership when `--root` and
`--chain-home` are given, the `s2_check_inputs.py` result). The stage's own files are only read.

Two re-stages under new file names (the staged files are only read): `restage-lr2` (card R-19: `LR2` from its four
archives) and `restage-nx1` (card R-20: `NX1` without the three 16S alignment files of the Lost City archive):

    c76_lr2_sag1_nx1_r232.py restage-nx1 OLD_META OLD_FASTA OUT_META OUT_FASTA FINAL_DIR [--drop GENOME_ID ...]
                                         [--compressor bgzip|gzip] [--bgzip-threads N] [--source-release DATE]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import re
import subprocess
import sys
import tarfile
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import r232_index  # noqa: E402

ACC = re.compile(r"GC[AF]_\d{9}")
GENOME_COLS = ("genome_id", "genome")
ROOT_STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
STATUS_HEAD = ["genome_id", "source_name", "accession", "how", "status", "in_gtdb_r232", "is_representative",
               "representative", "species", "in_globdb_r232", "globdb_r232_id"]
STATUSES = ("in_gtdb_r232+in_globdb_r232", "in_gtdb_r232", "in_globdb_r232", "not_in_r232", "not_visible")


class Refused(Exception):
    pass


def accession_of(genome_id: str) -> str | None:
    """The versionless GCA_/GCF_ accession in a genome id, or None."""
    m = ACC.search(genome_id)
    return m.group(0) if m else None


def read_names(path) -> dict[str, str]:
    """ENA assembly table (`accession`, `assembly_name`, ...): name to versionless accession. A name that two
    accessions share is ambiguous and is left out."""
    names: dict[str, set[str]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            names.setdefault(row["assembly_name"], set()).add(r232_index.nover(row["accession"]))
    return {n: next(iter(a)) for n, a in sorted(names.items()) if len(a) == 1}


def index_lookup(index_dir):
    """`lookup(accessions)` of card C-74 over its index."""
    return lambda accs: r232_index.lookup(index_dir, accs)


def _read_records(records_tsv) -> tuple[list[str], list[list[str]]]:
    with open(records_tsv, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        head = next(rd)
        return head, [r for r in rd]


def _genome_col(head: list[str]) -> int:
    for c in GENOME_COLS:
        if c in head:
            return head.index(c)
    raise Refused(f"records.tsv has no genome column ({' or '.join(GENOME_COLS)})")


def genome_accessions(records_tsv, names: dict[str, str] | None = None) -> list[dict]:
    """One row per genome of records.tsv in first-seen order: genome_id, source_name, accession, how."""
    head, rows = _read_records(records_tsv)
    gi = _genome_col(head)
    si = head.index("source_name") if "source_name" in head else None
    out, seen = [], set()
    for r in rows:
        g = r[gi]
        if g in seen:
            continue
        seen.add(g)
        acc, how = accession_of(g), "genome_id"
        if acc is None and names:
            acc = names.get(g) or names.get(g.split(":", 1)[-1])
            how = "ena_name"
        out.append({"genome_id": g, "source_name": r[si] if si is not None else "-", "accession": acc,
                    "how": how if acc else None})
    return out


def status(lookup, records_tsv, names: dict[str, str] | None = None) -> dict:
    """The r232 status of every genome (`rows`, sorted by genome id) and the counts, by status and by source."""
    gens = genome_accessions(records_tsv, names)
    hits = {}
    accs = sorted({g["accession"] for g in gens if g["accession"]})
    if accs:
        hits = {h["accession"]: h for h in lookup(accs)}
    rows, by_source = [], {}
    for g in gens:
        h = hits.get(g["accession"]) if g["accession"] else None
        if h is None:
            st = "not_visible"
        else:
            st = "+".join(k for k, v in (("in_gtdb_r232", h["in_gtdb_r232"]),
                                         ("in_globdb_r232", h["in_globdb_r232"])) if v) or "not_in_r232"
        rows.append({**g, "status": st,
                     "in_gtdb_r232": bool(h and h["in_gtdb_r232"]), "is_representative": h["is_representative"] if h else None,
                     "representative": h["representative"] if h else None, "species": h["species"] if h else None,
                     "in_globdb_r232": bool(h and h["in_globdb_r232"]), "globdb_r232_id": h["globdb_r232_id"] if h else None})
        by_source.setdefault(g["source_name"], Counter())[st] += 1
    rows.sort(key=lambda r: r["genome_id"])
    return {"rows": rows, "genomes": len(rows),
            "by_status": {s: sum(1 for r in rows if r["status"] == s) for s in STATUSES},
            "with_accession": sum(1 for r in rows if r["accession"]),
            "without_accession": sum(1 for r in rows if not r["accession"]),
            "in_gtdb_r232": sum(1 for r in rows if r["in_gtdb_r232"]),
            "in_globdb_r232": sum(1 for r in rows if r["in_globdb_r232"]),
            "by_source": {s: {k: c[k] for k in sorted(c)} for s, c in sorted(by_source.items())}}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _atomic_text(path: Path, text: str) -> None:
    tmp = path.with_name("." + path.name + ".partial")
    tmp.write_text(text)
    tmp.replace(path)


def _status_text(rows: list[dict]) -> str:
    def cell(v):
        return "" if v is None else ("1" if v is True else "0" if v is False else str(v))
    return "\t".join(STATUS_HEAD) + "\n" + "".join(
        "\t".join(cell(r[k]) for k in STATUS_HEAD) + "\n" for r in rows)


def _stats(head: list[str], rows: list[list[str]]) -> dict:
    gi, ci, li = _genome_col(head), head.index("contig_uid"), head.index("length")
    contigs: dict[str, int] = {}
    for r in rows:
        contigs.setdefault(r[ci], int(r[li]))
    return {"genomes": len({r[gi] for r in rows}), "contigs": len(contigs), "bases": sum(contigs.values())}


# ---- exact ownership (card C-79's method) -----------------------------------------------------------------------

def owners(records_tsv: Path, root: Path | None, chain_home: Path | None, later: dict[str, Path], work: Path, *,
           threads: int = 2, memory_limit: str = "4GB") -> dict:
    import a26_r26_11_4_chain as a26
    import c62_metavr_inputs as c62

    work.mkdir(parents=True, exist_ok=True)
    chain = {s: chain_home / s for s, _ in a26.R26_11_4_CHAIN} if chain_home is not None else {}
    own = c62._owners_from_records(records_tsv, root, list(ROOT_STAGES) if root is not None else [], chain,
                                   work=work, threads=threads, memory_limit=memory_limit)
    lat = c62._owners_from_records(records_tsv, None, [], later, work=work, threads=threads,
                                   memory_limit=memory_limit) if later else {}
    head, rows = _read_records(records_tsv)
    ci, li = head.index("contig_uid"), head.index("length")
    contigs: dict[bytes, int] = {}
    for r in rows:
        contigs.setdefault(bytes.fromhex(r[ci]), int(r[li]))
    c, b, lc, lb = Counter(), Counter(), Counter(), Counter()
    for cu, n in contigs.items():
        if cu in own:
            c[own[cu]] += 1
            b[own[cu]] += n
        elif cu in lat:
            lc[lat[cu]] += 1
            lb[lat[cu]] += n
    new = [n for cu, n in contigs.items() if cu not in own]
    return {"owner_stages": [*(ROOT_STAGES if root is not None else ()), *chain], "later_stages": [*later],
            "held_contigs": dict(sorted(c.items())), "held_bases": dict(sorted(b.items())),
            "held_total_contigs": sum(c.values()), "held_total_bases": sum(b.values()),
            "later_contigs": dict(sorted(lc.items())), "later_bases": dict(sorted(lb.items())),
            "distinct_contigs": len(contigs), "distinct_bases": sum(contigs.values()),
            "new_contigs": len(new), "new_bases": sum(new)}


def _check(stage: str, meta: Path) -> dict:
    p = subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "s2_check_inputs.py"), stage,
                        str(meta / "genomes.json"), str(meta / "objects_staging.json")],
                       capture_output=True, text=True, check=False)
    doc = json.loads(p.stdout) if p.stdout.strip().startswith("{") else {"stdout": p.stdout}
    return {"exit": p.returncode, "problems": doc.get("problems"), "gene_call_mode": doc.get("gene_call_mode"),
            "genomes": doc.get("genomes"), "genome_input_rows": doc.get("genome_input_rows"),
            "stderr": p.stderr[-2000:]}


# ---- card R-19: LR2 re-staged from its four archives ---------------------------------------------------------------
#
# Card C-66's stager (`/tmp/if2453/c66/run_stage_lr2.py` on the pod, not in the tree) gunzipped a member of the
# Singleton archive when its name ended in `.gz` and read the members of the other archives as plain text. The members
# of the Bickhart bins archive are gzip files, so their compressed bytes went into `lr2_compendia.fna.gz` as sequence
# (card R-19's note). This re-stage follows C-66's logic step for step, with three changes: a member is decompressed
# when its first two bytes are the gzip magic, whatever its name; every record is refused unless its header is UTF-8
# and its sequence holds IUPAC letters only; and a genome of the two sources C-66 read correctly keeps the genetic code
# C-66's screen gave it (its `genetic_codes.tsv` row), so only the bins are screened again.

LR2_ARCHIVES = {
    "emerge": "EMERGE_MAGs_2019_long-short-hybrid.tar.gz",
    "singleton": "MQ_HQ_MAGs_20200502.tar.gz",
    "bins": "hifi_das.bin3c.bins.tar.gz",
    "flye": "metaflye_hifi_asm.tar.gz",
}
LR2_PRESET_SOURCES = ("emerge_stordalen_2019_hifi", "singleton_2021_wwtp_ont")
LR2_BINS_SOURCE = "bickhart_2022_rumen_hifi"  # its bins are bins of the metaFlye assembly the unbinned set is cut from
LR2_RECORDS_HEAD = ["genome_id", "source_name", "record_ordinal", "contig_uid", "length", "circular", "kind",
                    "assembly_quality", "genetic_code", "taxonomy"]
LR2_GC_HEAD = ["genome_id", "source_name", "genome_bp", "n_contigs", "assembly_quality", "genetic_code", "reason",
               "tested_alt", "cd_11", "mean_orf_11", "cd_25", "cd_4"]
LR2_GC_KEYS = ("reason", "tested_alt", "cd_11", "mean_orf_11", "cd_25", "cd_4")
LR2_MIN_UNBINNED = 5_000
_IUPAC = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn"
_FASTA_EXTS = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def _is_fasta_member(name: str) -> bool:
    low = name.lower()
    base = low.rsplit("/", 1)[-1]
    return not base.startswith(".") and any(low.endswith(e) for e in _FASTA_EXTS)


def member_stream(tf, member):
    """A tar member's bytes, gunzipped when its first two bytes are the gzip magic (not by its name)."""
    import gzip as _gzip

    ef = tf.extractfile(member)
    if ef is None:
        raise Refused(f"{member.name}: not a regular file")
    if ef.peek(2)[:2] == b"\x1f\x8b":
        return _gzip.GzipFile(fileobj=ef, mode="rb")
    return ef


def strict_records(fh, where: str):
    """(header, sequence) of a FASTA stream, split as card C-66's reader split it (CR/LF stripped, blank lines
    skipped); refused at a line before the first header, a header that is not UTF-8 or a byte that is not an IUPAC
    nucleotide letter."""
    header, seq, n = None, [], 0
    for raw in fh:
        if not raw.strip():
            continue
        ln = raw.rstrip(b"\r\n")
        if ln.startswith(b">"):
            if header is not None:
                yield header, b"".join(seq)
            n += 1
            try:
                header = ln[1:].decode("utf-8")
            except UnicodeDecodeError:
                raise Refused(f"{where}: record {n}: header is not UTF-8") from None
            seq = []
        else:
            if header is None:
                raise Refused(f"{where}: a line before the first header")
            if ln.translate(None, _IUPAC):
                raise Refused(f"{where}: record {n} ({header.split()[0] if header.split() else ''}): "
                              "a byte that is not an IUPAC nucleotide letter")
            seq.append(ln)
    if header is not None:
        yield header, b"".join(seq)


def read_preset_codes(genetic_codes_tsv) -> dict[str, list[str]]:
    """C-66's `genetic_codes.tsv`: genome id to its row (strings, as written)."""
    with open(genetic_codes_tsv, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        if head != LR2_GC_HEAD:
            raise Refused(f"{genetic_codes_tsv}: unexpected header {head}")
        return {r[0]: r for r in rd}


def _screen(args):
    import c66_rt_hotspot_compendia as c66

    seqs, tax = args
    return c66.screen_genetic_code(seqs, tax, screen_max_bp=250_000)


class _LR2Writer:
    """The stage's FASTA and tables, written as each genome is decided (C-66's `flush_batch`, kept step for step)."""

    def __init__(self, fasta_out, records_w, gc_w, genomes_fh):
        from biobox import ids

        self.ids = ids
        self.fa, self.rec_w, self.gc_w, self.gj = fasta_out, records_w, gc_w, genomes_fh
        self.written: set[bytes] = set()
        self.binned: set[bytes] = set()
        self.binned_names: set[str] = set()
        self.props: dict[bytes, tuple] = {}
        self.keys: set[tuple] = set()
        self.uid_len: dict[bytes, int] = {}
        self.counts = Counter()
        self.bases = Counter()
        self.genomes = Counter()
        self.codes = Counter()
        self.trims: list[list] = []
        self.repeated: list[list] = []
        self.n_genomes = 0

    def _contig(self, uid: bytes, seq: bytes) -> None:
        self.fa.write(b">" + uid.hex().encode("ascii") + b"\n")
        for i in range(0, len(seq), 80):
            self.fa.write(seq[i:i + 80] + b"\n")

    def _genome(self, g: dict) -> None:
        self.gj.write(("[" if self.n_genomes == 0 else ",") + json.dumps(g, separators=(",", ":")))
        self.n_genomes += 1

    def mag(self, gid, src, recs, qual, tax, gc_info: dict, preset_row: list[str] | None) -> None:
        import c66_rt_hotspot_compendia as c66

        code = int(preset_row[5]) if preset_row is not None else int(gc_info["genetic_code"])
        g_seen, proc = set(), []
        for hdr, seq in recs:
            is_circ_flag = bool(c66._CIRCULAR_RE.search(hdr)) and len(recs) == 1
            tseq, is_circ = seq, False
            if is_circ_flag:
                k = c66.detect_exact_dtr(seq)
                if k > 0:
                    tseq = seq[:-k]
                    self.trims.append([gid, hdr.split()[0], len(seq), k, len(tseq)])
                    self.counts["dtr_trimmed"] += 1
                is_circ = True
            uid = self.ids.contig_uid(tseq)
            if uid in g_seen:
                self.counts["within_mag_dup_contigs"] += 1
                continue
            g_seen.add(uid)
            if uid in self.props and self.props[uid] != ("prokaryote", code, is_circ):
                self.counts["dec27_property_conflicts_skipped"] += 1
                self.repeated.append([gid, src, "dec27_property_conflict", uid.hex()])
                return
            proc.append((hdr, tseq, uid, is_circ))
        if not proc:
            return
        key = tuple(sorted(u for _, _, u, _ in proc))
        if key in self.keys:
            self.counts["repeated_genomes"] += 1
            self.repeated.append([gid, src, "duplicate_genome_contig_set", key[0].hex()])
            return
        self.keys.add(key)
        g_uids, g_circ, g_bp = [], [], 0
        for rec_ord, (hdr, seq, uid, is_circ) in enumerate(proc):
            g_uids.append(uid)
            self.binned.add(uid)
            if src == LR2_BINS_SOURCE and hdr.split():
                self.binned_names.add(hdr.split()[0])
            if is_circ:
                g_circ.append(uid)
                self.counts["circular_contigs"] += 1
            if uid not in self.written:
                self.written.add(uid)
                self.props[uid] = ("prokaryote", code, is_circ)
                self.uid_len[uid] = len(seq)
                self._contig(uid, seq)
                self.bases[src] += len(seq)
            g_bp += len(seq)
            self.rec_w.writerow([gid, src, rec_ord, uid.hex(), len(seq), int(is_circ), "prokaryote", qual, code, tax])
        self._genome({"contig_uids": [u.hex() for u in g_uids], "genetic_code": code, "kind": "prokaryote",
                      "circular": [u.hex() for u in g_circ], "assembly_quality": qual})
        self.counts["mag_genomes"] += 1
        self.genomes[src] += 1
        self.codes[code] += 1
        if preset_row is not None:
            self.gc_w.writerow([gid, src, g_bp, len(g_uids), qual, code, *preset_row[6:]])
        else:
            self.gc_w.writerow([gid, src, g_bp, len(g_uids), qual, code, *(gc_info.get(k, "") for k in LR2_GC_KEYS)])

    def unbinned(self, src, cname, seq, is_circ_flag) -> None:
        import c66_rt_hotspot_compendia as c66

        s = seq.strip()
        if len(s) < LR2_MIN_UNBINNED:
            self.counts["unbinned_short_dropped"] += 1
            self.counts["unbinned_short_bases"] += len(s)
            return
        tseq, is_circ, trim = s, False, None
        if is_circ_flag:
            k = c66.detect_exact_dtr(s)
            if k > 0:
                tseq = s[:-k]
                trim = [f"bickhart_hifi:unbinned:{cname}", cname, len(s), k, len(tseq)]
            is_circ = True
        if len(tseq) < LR2_MIN_UNBINNED:
            self.counts["unbinned_short_dropped"] += 1
            return
        uid = self.ids.contig_uid(tseq)
        if uid in self.binned:
            self.counts["unbinned_already_binned_skipped"] += 1
            return
        if cname in self.binned_names:
            # a binned contig of the same assembly whose terminal repeat was trimmed here (circular in
            # assembly_info.txt) and not in its bin: the bin holds it, so it is not an unbinned contig
            self.counts["unbinned_binned_by_name_skipped"] += 1
            return
        if uid in self.written or (uid,) in self.keys:
            self.counts["unbinned_dup_skipped"] += 1
            return
        if trim is not None:  # a trim is recorded for a contig that is staged (C-66 recorded it before its skips)
            self.trims.append(trim)
            self.counts["dtr_trimmed"] += 1
        self.keys.add((uid,))
        self.written.add(uid)
        self.props[uid] = ("prokaryote", 11, is_circ)
        self.uid_len[uid] = len(tseq)
        self._contig(uid, tseq)
        if is_circ:
            self.counts["circular_contigs"] += 1
        self.bases[src] += len(tseq)
        self._genome({"contig_uids": [uid.hex()], "genetic_code": 11, "kind": "prokaryote",
                      "circular": [uid.hex()] if is_circ else [], "assembly_quality": "contig"})
        self.counts["unbinned_contigs_kept"] += 1
        self.genomes[src] += 1
        self.codes[11] += 1
        self.rec_w.writerow([f"bickhart_hifi:unbinned:{cname}", src, 0, uid.hex(), len(tseq), int(is_circ),
                             "prokaryote", "contig", 11, ""])


def _open_fasta_out(path: Path, compressor: str, threads: int):
    """(writer, closer) for the stage's FASTA: `bgzip` as card C-66 wrote it, or Python's gzip (the tests)."""
    import gzip as _gzip

    raw = open(path, "wb")
    if compressor == "bgzip":
        proc = subprocess.Popen(["bgzip", "-@", str(threads), "-c"], stdin=subprocess.PIPE, stdout=raw)

        def close():
            proc.stdin.close()
            rc = proc.wait()
            raw.close()
            if rc != 0:
                raise Refused(f"bgzip failed rc={rc}")
        return proc.stdin, close
    gz = _gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0)

    def close():
        gz.close()
        raw.close()
    return gz, close


def _rows_of(records_tsv: Path, sources) -> Iterable:
    with open(records_tsv, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for r in rd:
            if r[1] in sources:
                yield r


def restage_lr2(dl_dir, emerge_meta_tsv, old_meta, out_meta, out_fasta, final_dir, *, workers: int = 1,
                compressor: str = "gzip", bgzip_threads: int = 2, root: Path | None = None,
                chain_home: Path | None = None, work: Path | None = None,
                source_release: str = "2026-11-27") -> dict:
    """Re-stage LR2 from C-66's four archives into OUT_FASTA (a new name) and OUT_META; C-66's files are only read."""
    from concurrent.futures import ProcessPoolExecutor

    dl, old, out_meta, out_fasta = Path(dl_dir), Path(old_meta), Path(out_meta), Path(out_fasta)
    if out_fasta.exists():
        raise Refused(f"{out_fasta} exists")
    if out_meta.exists() and any(out_meta.iterdir()):
        raise Refused(f"{out_meta} exists and is not empty")
    arch = {k: dl / v for k, v in LR2_ARCHIVES.items()}
    for p in [*arch.values(), Path(emerge_meta_tsv), old / "genetic_codes.tsv", old / "objects_staging.json"]:
        if not Path(p).is_file():
            raise Refused(f"missing {p}")
    preset = read_preset_codes(old / "genetic_codes.tsv")
    e19 = {}
    with open(emerge_meta_tsv, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            e19[(r.get("filename") or "").strip()] = r
    import c66_rt_hotspot_compendia as c66

    out_meta.mkdir(parents=True, exist_ok=True)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    tmp_fasta = out_fasta.with_name(out_fasta.name + ".tmp")
    fa, close_fa = _open_fasta_out(tmp_fasta, compressor, bgzip_threads)
    rec_fh = open(out_meta / ".records.tsv.partial", "w", encoding="utf-8", newline="")
    gc_fh = open(out_meta / ".genetic_codes.tsv.partial", "w", encoding="utf-8", newline="")
    gj_fh = open(out_meta / ".genomes.json.partial", "w", encoding="utf-8")
    rec_w, gc_w = csv.writer(rec_fh, delimiter="\t"), csv.writer(gc_fh, delimiter="\t")
    rec_w.writerow(LR2_RECORDS_HEAD)
    gc_w.writerow(LR2_GC_HEAD)
    w = _LR2Writer(fa, rec_w, gc_w, gj_fh)
    members = Counter()
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None

    def flush(batch):
        need = [i for i, b in enumerate(batch) if b[5] is None]
        args = [([s for _, s in batch[i][2]], batch[i][4]) for i in need]
        res = list(pool.map(_screen, args)) if pool is not None else [_screen(a) for a in args]
        got = dict(zip(need, res))
        for i, (gid, src, recs, qual, tax, prow) in enumerate(batch):
            w.mag(gid, src, recs, qual, tax, got.get(i, {}), prow)
        batch.clear()

    def mags(key, src, ident):
        batch = []
        with tarfile.open(arch[key], "r:*") as tf:
            for m in tf:
                if not m.isfile() or not _is_fasta_member(m.name):
                    continue
                base = m.name.rsplit("/", 1)[-1]
                fh = member_stream(tf, m)
                members[(key, "gzip" if isinstance(fh, gzip_type()) else "plain")] += 1
                recs = [(h, s.strip()) for h, s in strict_records(fh, f"{arch[key].name}:{m.name}") if s.strip()]
                gid, qual, tax = ident(base)
                prow = preset.get(gid) if src in LR2_PRESET_SOURCES else None
                if recs:
                    if src in LR2_PRESET_SOURCES and prow is None:
                        raise Refused(f"{gid}: no genetic_codes.tsv row of card C-66")
                    batch.append((gid, src, recs, qual, tax, prow))
                if len(batch) >= 48:
                    flush(batch)
        flush(batch)

    def e19_ident(base):
        meta = e19.get(base, {})
        comp = c66._safe_float(meta.get("Completeness"), 80.0)
        cont = c66._safe_float(meta.get("Contamination"), 0.0)
        qual = "mag_hq" if (comp > 90.0 and cont < 5.0) else "mag_mq"
        return f"emerge2019:{base.rsplit('.', 1)[0]}", qual, (meta.get("classification") or "").strip()

    def stem(base):
        return base.replace(".gz", "").rsplit(".", 1)[0]

    try:
        mags("emerge", "emerge_stordalen_2019_hifi", e19_ident)
        mags("singleton", "singleton_2021_wwtp_ont", lambda b: (f"singleton2021:{stem(b)}", "mag_hq", ""))
        mags("bins", LR2_BINS_SOURCE, lambda b: (f"bickhart_bin:{stem(b)}", "mag_hq", ""))
        flye_circ = set()
        with tarfile.open(arch["flye"], "r:*") as tf:
            for m in tf:
                if m.isfile() and m.name.endswith("assembly_info.txt"):
                    for line in io.TextIOWrapper(member_stream(tf, m), encoding="utf-8", errors="replace"):
                        if line.startswith("#") or not line.strip():
                            continue
                        p = line.rstrip("\r\n").split("\t")
                        if len(p) >= 4 and p[3].strip().upper() == "Y":
                            flye_circ.add(p[0].strip())
                    break
        with tarfile.open(arch["flye"], "r:*") as tf:
            for m in tf:
                if m.isfile() and m.name.endswith("assembly.fasta"):
                    where = f"{arch['flye'].name}:{m.name}"
                    for hdr, seq in strict_records(member_stream(tf, m), where):
                        cname = hdr.split()[0]
                        w.unbinned("bickhart_2022_rumen_hifi_unbinned", cname, seq,
                                   (cname in flye_circ) or bool(c66._CIRCULAR_RE.search(hdr)))
                    break
    finally:
        if pool is not None:
            pool.shutdown()
        close_fa()
        gj_fh.write(("[" if w.n_genomes == 0 else "") + "]\n")
        for f in (rec_fh, gc_fh, gj_fh):
            f.close()
    for n in ("records.tsv", "genetic_codes.tsv", "genomes.json"):
        (out_meta / f".{n}.partial").replace(out_meta / n)
    tmp_fasta.replace(out_fasta)

    def wtsv(name, head, rows):
        with open(out_meta / name, "w", encoding="utf-8", newline="") as fh:
            wr = csv.writer(fh, delimiter="\t")
            wr.writerow(head)
            wr.writerows(rows)

    wtsv("repeated_records.tsv", ["genome_id", "source_name", "reason", "first_contig_uid"], w.repeated)
    wtsv("trims.tsv", ["genome_id", "contig_name", "raw_len", "trimmed_k", "final_len"], w.trims)
    shared: list[list] = []
    if root is not None or chain_home is not None:
        import a26_r26_11_4_chain as a26
        import c62_metavr_inputs as c62

        chain = {s: chain_home / s for s, _ in a26.R26_11_4_CHAIN} if chain_home is not None else {}
        own = c62._owners_from_records(out_meta / "records.tsv", root, list(ROOT_STAGES) if root is not None else [],
                                       chain, work=Path(work) if work else out_meta.parent / "_work_r19")
        shared = [[u.hex(), own[u], w.uid_len[u]] for u in sorted(own) if u in w.uid_len]
    wtsv("shared_records.tsv", ["contig_uid", "owner_stage", "length"], shared)
    obj = json.loads((old / "objects_staging.json").read_text(encoding="utf-8"))
    if len(obj) != 1:
        raise Refused(f"{old / 'objects_staging.json'}: expected one object")
    o = {**obj[0], "path": str(out_fasta), "source_release": source_release}
    (out_meta / "objects_staging.json").write_text(json.dumps([o], indent=2) + "\n", encoding="utf-8")
    (out_meta / "objects.json").write_text(
        json.dumps([{**o, "path": str(Path(final_dir) / out_fasta.name)}], indent=2) + "\n", encoding="utf-8")
    # the two sources C-66 read correctly come out row for row as C-66 wrote them
    old_rows = _rows_of(old / "records.tsv", LR2_PRESET_SOURCES)
    new_rows = _rows_of(out_meta / "records.tsv", LR2_PRESET_SOURCES)
    compared, first_diff = 0, None
    for a, b in itertools.zip_longest(old_rows, new_rows):
        if a != b:
            first_diff = {"row": compared + 1, "old": a, "new": b}
            break
        compared += 1
    shared_by = Counter(r[1] for r in shared)
    total_bp = sum(w.bases.values())
    rep = {
        "card": "R-19", "stage": "LR2", "fasta_path": str(out_fasta), "fasta_bytes": out_fasta.stat().st_size,
        "fasta_sha256": _sha256(out_fasta),
        "files_sha256": {p.name: _sha256(p) for p in sorted(out_meta.iterdir()) if p.is_file()
                         and not p.name.startswith(".")},
        "members": {f"{k}:{c}": n for (k, c), n in sorted(members.items())},
        "total_genomes": w.n_genomes, "mag_genomes": w.counts["mag_genomes"],
        "unbinned_contigs_kept": w.counts["unbinned_contigs_kept"],
        "unbinned_short_dropped": w.counts["unbinned_short_dropped"],
        "unbinned_short_bases": w.counts["unbinned_short_bases"],
        "unbinned_already_binned_skipped": w.counts["unbinned_already_binned_skipped"],
        "unbinned_binned_by_name_skipped": w.counts["unbinned_binned_by_name_skipped"],
        "unbinned_dup_skipped": w.counts["unbinned_dup_skipped"],
        "repeated_genomes": w.counts["repeated_genomes"],
        "dec27_property_conflicts_skipped": w.counts["dec27_property_conflicts_skipped"],
        "within_mag_dup_contigs": w.counts["within_mag_dup_contigs"],
        "dtr_trimmed": w.counts["dtr_trimmed"], "circular_contigs": w.counts["circular_contigs"],
        "distinct_contigs": len(w.written), "total_bp": total_bp,
        "shared_contigs": len(shared), "shared_bp": sum(r[2] for r in shared),
        "net_new_contigs": len(w.written) - len(shared), "net_new_bp": total_bp - sum(r[2] for r in shared),
        "shared_by_stage": dict(sorted(shared_by.items())),
        "genomes_by_source": dict(sorted(w.genomes.items())), "bases_by_source": dict(sorted(w.bases.items())),
        "genetic_codes": {str(k): v for k, v in sorted(w.codes.items())},
        "preset_rows_compared": compared, "preset_rows_first_difference": first_diff,
    }
    (out_meta / "stage_report.json").write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    return rep


def gzip_type():
    import gzip as _gzip

    return _gzip.GzipFile


# ---- card R-20: NX1 re-staged without the three 16S alignment files of the Lost City archive ----------------------
#
# Card C-66 took every FASTA of the Lost City zip as a MAG, so three "genomes" of NX1 are 16S rRNA alignments (their
# records hold gap bytes `-`; card R-19's scan). The re-stage reads C-66's staged FASTA and tables and writes them again
# under a new name without those genomes: a record is left out only when no kept genome holds its contig, every kept
# record and table row is copied byte for byte, and the kept records are checked (header UTF-8, IUPAC letters only,
# header the contig_uid of the bytes) as they are read back from the new file.

NX1_ALIGNMENTS = ("lost_city:Bipolaricaulota-16S-phylogeny", "lost_city:Methanosarcinales-16S-phylogeny",
                  "lost_city:Thermodesulfovibrionales-16S-phylogeny")
NX1_REASON = "16S rRNA alignment file of the source archive, not genome sequence (card R-20)"
NX1_GENOME_TABLES = ("genetic_codes.tsv", "trims.tsv", "repeated_records.tsv")  # first column genome_id
NX1_OPTIONAL_TABLES = ("skani_within_dropped.tsv",)
_SKANI_REP = b"_rep:"


def _filter_lines(src: Path, dst: Path, keep) -> tuple[int, int]:
    """Copy the header line and every line `keep(fields)` accepts, bytes and line ends as they are; (kept, left out)."""
    k = d = 0
    with open(src, "rb") as fi, open(dst, "wb") as fo:
        fo.write(fi.readline())
        for ln in fi:
            if keep(ln.rstrip(b"\r\n").split(b"\t")):
                fo.write(ln)
                k += 1
            else:
                d += 1
    return k, d


def _record_lines(fh):
    """(header line, [sequence lines]) of a FASTA stream, every line with its own line end, blank lines kept."""
    head, body = None, []
    for ln in fh:
        if ln.startswith(b">"):
            if head is not None:
                yield head, body
            head, body = ln, []
        elif head is None:
            if ln.strip():
                raise Refused("a line before the first header")
        else:
            body.append(ln)
    if head is not None:
        yield head, body


def restage_nx1(old_meta, old_fasta, out_meta, out_fasta, final_dir, *, drop=NX1_ALIGNMENTS,
                compressor: str = "gzip", bgzip_threads: int = 1, source_release: str = "2026-11-27") -> dict:
    """NX1 without the genomes `drop`, into OUT_FASTA (a new name) and OUT_META; C-66's files are only read."""
    import gzip as _gzip

    from biobox import ids

    old, old_fasta, out_meta, out_fasta = Path(old_meta), Path(old_fasta), Path(out_meta), Path(out_fasta)
    if out_fasta.exists():
        raise Refused(f"{out_fasta} exists")
    if out_meta.exists() and any(out_meta.iterdir()):
        raise Refused(f"{out_meta} exists and is not empty")
    need = ["records.tsv", "genomes.json", "shared_records.tsv", "objects_staging.json", *NX1_GENOME_TABLES]
    for n in need:
        if not (old / n).is_file():
            raise Refused(f"missing {old / n}")
    if not old_fasta.is_file():
        raise Refused(f"missing {old_fasta}")
    drop = tuple(sorted(set(drop)))
    drop_b = {g.encode() for g in drop}

    # records.tsv: each genome's contigs in C-66's order, the contigs to leave out, the counts of what stays
    order: list[bytes] = []
    contigs_of: dict[bytes, list[bytes]] = {}
    src_of: dict[bytes, bytes] = {}
    kept_uids: set[bytes] = set()
    drop_rows: Counter = Counter()
    drop_row_bases: Counter = Counter()
    drop_uids: dict[bytes, set[bytes]] = {}
    with open(old / "records.tsv", "rb") as fh:
        head = fh.readline().rstrip(b"\r\n").split(b"\t")
        gi, si, ci, li = (head.index(c) for c in (b"genome_id", b"source_name", b"contig_uid", b"length"))
        for ln in fh:
            f = ln.rstrip(b"\r\n").split(b"\t")
            g, u = f[gi], f[ci]
            if g not in contigs_of:
                order.append(g)
                contigs_of[g] = []
                src_of[g] = f[si]
            contigs_of[g].append(u)
            if g in drop_b:
                drop_rows[g] += 1
                drop_row_bases[g] += int(f[li])
                drop_uids.setdefault(g, set()).add(u)
            else:
                kept_uids.add(u)
    for g in drop:
        if g.encode() not in contigs_of:
            raise Refused(f"{g}: not in {old / 'records.tsv'}")
    gone = set().union(*drop_uids.values()) - kept_uids if drop_uids else set()
    shared_with_kept = set().union(*drop_uids.values()) & kept_uids if drop_uids else set()

    # genomes.json follows records.tsv genome for genome, and is written back as C-66 wrote it
    raw = (old / "genomes.json").read_bytes()
    gj = json.loads(raw)
    if len(gj) != len(order) or any([u.encode() for u in e["contig_uids"]] != contigs_of[g]
                                    for e, g in zip(gj, order)):
        raise Refused(f"{old / 'genomes.json'} does not follow records.tsv genome for genome")
    body = json.dumps(gj, separators=(",", ":")).encode()
    gj_end = raw[len(body):]  # C-66 ended the file with a newline
    if not raw.startswith(body) or gj_end.strip():
        raise Refused(f"{old / 'genomes.json'} is not in the compact form this re-stage writes")
    # a genome C-66's skani step dropped for a removed representative would be lost with it
    for n in NX1_OPTIONAL_TABLES:
        if (old / n).is_file():
            with open(old / n, "rb") as fh:
                next(fh)
                for ln in fh:
                    f = ln.rstrip(b"\r\n").split(b"\t")
                    rep = f[-1].split(_SKANI_REP, 1)[-1] if _SKANI_REP in f[-1] else b""
                    if rep in drop_b:
                        raise Refused(f"{old / n}: {f[0].decode()} was dropped by skani for {rep.decode()}, "
                                      "a genome this re-stage removes")

    out_meta.mkdir(parents=True, exist_ok=True)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    tmp_fasta = out_fasta.with_name(out_fasta.name + ".tmp")

    # the FASTA: every record whose contig stays, its lines as they are
    gap_of: dict[bytes, int] = {}
    kept_sha = hashlib.sha256()
    n_old = old_bp = n_kept = kept_bp = removed_bp = removed_gaps = 0
    seen: set[bytes] = set()
    def check(hl: bytes, body: list[bytes], where: str, k: int) -> str:
        """The record's header (refused unless UTF-8) and its sequence (refused unless IUPAC letters only)."""
        try:
            h = hl[1:].rstrip(b"\r\n").decode("utf-8")
        except UnicodeDecodeError:
            raise Refused(f"{where}: record {k}: header is not UTF-8") from None
        name = h.split()[0] if h.split() else ""
        if any(x.rstrip(b"\r\n").translate(None, _IUPAC) for x in body):
            raise Refused(f"{where}: record {k} ({name}): a byte that is not an IUPAC nucleotide letter")
        return name

    fa, close_fa = _open_fasta_out(tmp_fasta, compressor, bgzip_threads)
    try:
        try:
            with _gzip.open(old_fasta, "rb") as fh:
                for hl, body in _record_lines(fh):
                    n_old += 1
                    u = hl[1:].rstrip(b"\r\n").split(None, 1)[0] if hl[1:].strip() else b""
                    n = sum(len(x.rstrip(b"\r\n")) for x in body)
                    old_bp += n
                    if u in gone:
                        g = sum(x.count(b"-") for x in body)
                        gap_of[u] = g
                        removed_bp += n
                        removed_gaps += g
                        continue
                    check(hl, body, str(old_fasta), n_old)
                    if u not in kept_uids:
                        raise Refused(f"{old_fasta}: record {n_old} ({u.decode()}) is in no genome")
                    if u in seen:
                        raise Refused(f"{old_fasta}: record {n_old} ({u.decode()}) appears twice")
                    seen.add(u)
                    gap_of[u] = 0
                    for x in (hl, *body):
                        fa.write(x)
                        kept_sha.update(x)
                    n_kept += 1
                    kept_bp += n
        finally:
            close_fa()
        if seen != kept_uids:
            raise Refused(f"{len(kept_uids - seen)} contigs of the kept genomes are not in {old_fasta}")
        if gone - set(gap_of):
            raise Refused(f"{len(gone - set(gap_of))} contigs of the removed genomes are not in {old_fasta}")

        # read the new file back: the same bytes as the kept records of the old one, every record clean and keyed
        # by the contig_uid of its bytes
        back, mismatches, n_back = hashlib.sha256(), 0, 0
        with _gzip.open(tmp_fasta, "rb") as fh:
            for hl, body in _record_lines(fh):
                for x in (hl, *body):
                    back.update(x)
                n_back += 1
                name = check(hl, body, str(out_fasta), n_back)
                if ids.contig_uid(b"".join(x.rstrip(b"\r\n") for x in body)).hex() != name:
                    mismatches += 1
    except BaseException:
        tmp_fasta.unlink(missing_ok=True)
        raise
    identical = back.hexdigest() == kept_sha.hexdigest() and n_back == n_kept
    tmp_fasta.replace(out_fasta)

    # the tables, less the removed genomes' rows (and the removed contigs' rows of shared_records.tsv)
    part = {}
    for n in ("records.tsv", *NX1_GENOME_TABLES):
        p = out_meta / f".{n}.partial"
        part[n] = _filter_lines(old / n, p, lambda f: f[0] not in drop_b)
    part["shared_records.tsv"] = _filter_lines(old / "shared_records.tsv", out_meta / ".shared_records.tsv.partial",
                                               lambda f: f[0] not in gone)
    for n in NX1_OPTIONAL_TABLES:
        if (old / n).is_file():
            part[n] = _filter_lines(old / n, out_meta / f".{n}.partial", lambda f: f[0] not in drop_b)
    kept_gj = [e for e, g in zip(gj, order) if g not in drop_b]
    (out_meta / ".genomes.json.partial").write_bytes(json.dumps(kept_gj, separators=(",", ":")).encode() + gj_end)
    for n in [*part, "genomes.json"]:
        (out_meta / f".{n}.partial").replace(out_meta / n)
    removed_rows = [[g, src_of[g.encode()].decode(), str(drop_rows[g.encode()]), str(drop_row_bases[g.encode()]),
                     str(sum(gap_of.get(u, 0) for u in drop_uids[g.encode()])), NX1_REASON] for g in drop]
    with open(out_meta / "removed_genomes.tsv", "w", encoding="utf-8", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["genome_id", "source_name", "records", "bases", "gap_bytes", "reason"])
        wr.writerows(removed_rows)
    obj = json.loads((old / "objects_staging.json").read_text(encoding="utf-8"))
    if len(obj) != 1:
        raise Refused(f"{old / 'objects_staging.json'}: expected one object")
    o = {**obj[0], "path": str(out_fasta), "source_release": source_release}
    (out_meta / "objects_staging.json").write_text(json.dumps([o], indent=2) + "\n", encoding="utf-8")
    (out_meta / "objects.json").write_text(
        json.dumps([{**o, "path": str(Path(final_dir) / out_fasta.name)}], indent=2) + "\n", encoding="utf-8")

    with open(out_meta / "shared_records.tsv", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        shared = [(r[0], r[1], int(r[2])) for r in rd]
    kept_genomes = [g for g in order if g not in drop_b]
    rep = {
        "card": "R-20", "stage": "NX1",
        "old_fasta_path": str(old_fasta), "old_fasta_bytes": old_fasta.stat().st_size,
        "old_fasta_sha256": _sha256(old_fasta),
        "fasta_path": str(out_fasta), "fasta_bytes": out_fasta.stat().st_size, "fasta_sha256": _sha256(out_fasta),
        "files_sha256": {p.name: _sha256(p) for p in sorted(out_meta.iterdir()) if p.is_file()
                         and not p.name.startswith(".")},
        "removed_genome_ids": list(drop), "removed_genomes": len(drop),
        "removed_records_tsv_rows": sum(drop_rows.values()),
        "removed_contigs": len(gone), "removed_bases": removed_bp, "removed_gap_bytes": removed_gaps,
        "removed_contigs_kept_for_another_genome": len(shared_with_kept),
        "old_genomes": len(order), "old_distinct_contigs": n_old, "old_total_bp": old_bp,
        "total_genomes": len(kept_genomes), "distinct_contigs": n_kept, "total_bp": kept_bp,
        "records_tsv_rows": part["records.tsv"][0],
        "kept_stream_sha256": kept_sha.hexdigest(), "read_back_sha256": back.hexdigest(),
        "kept_records_byte_identical": identical, "uid_mismatches": mismatches,
        "shared_contigs": len(shared), "shared_bp": sum(r[2] for r in shared),
        "shared_by_stage": dict(sorted(Counter(r[1] for r in shared).items())),
        "net_new_contigs": n_kept - len(shared), "net_new_bp": kept_bp - sum(r[2] for r in shared),
        "genomes_by_source": dict(sorted(Counter(src_of[g].decode() for g in kept_genomes).items())),
        "genetic_codes": {str(k): v for k, v in sorted(Counter(e["genetic_code"] for e in kept_gj).items())},
        "table_rows_kept_and_removed": {n: list(v) for n, v in sorted(part.items())},
    }
    (out_meta / "stage_report.json").write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    return rep


# ---- the command -------------------------------------------------------------------------------------------------

def run(stage: str, meta_dir, fasta, index_dir, out_meta_root, *, names: dict[str, str] | None = None,
        check: bool = True, lookup=None, root: Path | None = None, chain_home: Path | None = None,
        later: dict[str, Path] | None = None, work: Path | None = None, hash_files: bool = True) -> dict:
    """Record the r232 status of every genome of the stage; the stage's files are read, never written."""
    meta, fasta, out_dir = Path(meta_dir), Path(fasta), Path(out_meta_root) / stage
    targets = [out_dir / "r232_status.tsv", out_dir / "r232_report.json"]
    if any(t.exists() for t in targets):
        raise Refused(f"{out_dir} exists")
    for n in ("records.tsv", "genomes.json", "objects.json", "objects_staging.json"):
        if not (meta / n).is_file():
            raise Refused(f"missing {meta / n}")
    if not fasta.is_file():
        raise Refused(f"missing {fasta}")
    st = status(lookup or index_lookup(index_dir), meta / "records.tsv", names)
    head, rows = _read_records(meta / "records.tsv")
    rep: dict = {"stage": stage, "removed": 0, "genomes": st["genomes"], "with_accession": st["with_accession"],
                 "without_accession": st["without_accession"], "in_gtdb_r232": st["in_gtdb_r232"],
                 "in_globdb_r232": st["in_globdb_r232"], "by_status": st["by_status"],
                 "by_source": st["by_source"], "stage_totals": _stats(head, rows)}
    out_dir.mkdir(parents=True, exist_ok=True)
    _atomic_text(out_dir / "r232_status.tsv", _status_text(st["rows"]))
    if root is not None or chain_home is not None or later:
        rep["ownership"] = owners(meta / "records.tsv", root, chain_home, later or {},
                                  Path(work) if work else out_dir / "_work")
    if check:
        rep["s2_check"] = _check(stage, meta)
    if hash_files:
        files = {p.name: {"bytes": p.stat().st_size, "sha256": _sha256(p)}
                 for p in sorted(meta.iterdir()) if p.is_file() and not p.name.startswith(".")}
        files[fasta.name] = {"bytes": fasta.stat().st_size, "sha256": _sha256(fasta)}
        rep["files"] = files
    rep["status_tsv_sha256"] = _sha256(out_dir / "r232_status.tsv")
    _atomic_text(out_dir / "r232_report.json", json.dumps(rep, indent=1, sort_keys=True) + "\n")
    return rep


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="R232 / r232 status of a staged LR2, SAG1 or NX1 (card C-76, DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("stage")
    for n in ("meta", "fasta", "index", "out_meta"):
        r.add_argument(n, type=Path)
    r.add_argument("--names", type=Path, default=None)
    r.add_argument("--root", type=Path, default=None)
    r.add_argument("--chain-home", type=Path, default=None)
    r.add_argument("--later", action="append", default=[], help="STAGE=DIR of an earlier r26.11.5 stage")
    r.add_argument("--no-check", action="store_true")
    r.add_argument("--work", type=Path, default=None)
    r.add_argument("--report-out", type=Path, default=None)
    s = sub.add_parser("restage-lr2", help="card R-19: LR2 from C-66's four archives, under a new FASTA name")
    for n in ("dl_dir", "emerge_meta", "old_meta", "out_meta", "out_fasta", "final_dir"):
        s.add_argument(n, type=Path)
    s.add_argument("--workers", type=int, default=1)
    s.add_argument("--compressor", choices=("bgzip", "gzip"), default="bgzip")
    s.add_argument("--bgzip-threads", type=int, default=2)
    s.add_argument("--root", type=Path, default=None)
    s.add_argument("--chain-home", type=Path, default=None)
    s.add_argument("--work", type=Path, default=None)
    s.add_argument("--source-release", default="2026-11-27")
    x = sub.add_parser("restage-nx1", help="card R-20: NX1 without its three 16S alignment files, under a new name")
    for n in ("old_meta", "old_fasta", "out_meta", "out_fasta", "final_dir"):
        x.add_argument(n, type=Path)
    x.add_argument("--drop", action="append", default=None,
                   help="a genome id to leave out (default: the three Lost City 16S alignment files)")
    x.add_argument("--compressor", choices=("bgzip", "gzip"), default="bgzip")
    x.add_argument("--bgzip-threads", type=int, default=1)
    x.add_argument("--source-release", default="2026-11-27")
    a = ap.parse_args(argv)
    if a.cmd == "restage-nx1":
        try:
            rep = restage_nx1(a.old_meta, a.old_fasta, a.out_meta, a.out_fasta, a.final_dir,
                              drop=tuple(a.drop) if a.drop else NX1_ALIGNMENTS, compressor=a.compressor,
                              bgzip_threads=a.bgzip_threads, source_release=a.source_release)
        except Refused as exc:
            print(f"refused: {exc}", file=sys.stderr)
            return 2
        print(json.dumps(rep, indent=1, sort_keys=True))
        return 0 if rep["kept_records_byte_identical"] and rep["uid_mismatches"] == 0 else 1
    if a.cmd == "restage-lr2":
        try:
            rep = restage_lr2(a.dl_dir, a.emerge_meta, a.old_meta, a.out_meta, a.out_fasta, a.final_dir,
                              workers=a.workers, compressor=a.compressor, bgzip_threads=a.bgzip_threads,
                              root=a.root, chain_home=a.chain_home, work=a.work, source_release=a.source_release)
        except Refused as exc:
            print(f"refused: {exc}", file=sys.stderr)
            return 2
        print(json.dumps(rep, indent=1, sort_keys=True))
        return 0 if rep["preset_rows_first_difference"] is None else 1
    try:
        later = {}
        for item in a.later:
            k, _, v = item.partition("=")
            later[k] = Path(v)
        rep = run(a.stage, a.meta, a.fasta, a.index, a.out_meta, names=read_names(a.names) if a.names else None,
                  check=not a.no_check, root=a.root, chain_home=a.chain_home, later=later, work=a.work)
    except (Refused, r232_index.Refused) as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    txt = json.dumps(rep, indent=1, sort_keys=True)
    if a.report_out is not None:
        a.report_out.parent.mkdir(parents=True, exist_ok=True)
        a.report_out.write_text(txt + "\n")
    print(txt)
    s = rep.get("s2_check")
    return 0 if s is None or (s["exit"] == 0 and s["problems"] == []) else 1


if __name__ == "__main__":
    sys.exit(main())

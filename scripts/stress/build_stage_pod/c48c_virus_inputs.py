"""Card C-48c (part 2): the `build-stage` inputs of stage V2, IMG/VR v4.1's high-confidence, high-quality vOTU
representatives plus the Inovirus genomes, streamed from the archives on /stor (read only), in the form of card C-39.

    python c48c_virus_inputs.py v2 HC_ZIP META_ZIP OUT_DIR [--root ROOT] [--root-stages S1,P1,...]
        [--owner-inputs NAME=DIR ...]
    python c48c_virus_inputs.py gv1 GVMAGS_TAR OUT_DIR [--root ROOT] [--root-stages S1,P1,...] [--owner-inputs NAME=DIR ...]

gv1 (card C-48c part 3; stage GV1): every genome of GVMAGs V2, the `.fna` members of all_18kgenomes_fna.tar.gz inside
GVMAGS_TAR (gvmags_v2_20261118.tar), streamed one member at a time; the member's md5 is checked against sources.yaml. A
genome is one `.fna` (its records its contigs, in order; kind virus, gv-meta; code 11 recorded; assembly_quality
`contig`, the schema's value for an assembly of unknown completeness, since the archive has no quality table; every
contig linear, the archive records no topology). A member that is not `.fna` is not read as a genome and is written
to skipped_members.tsv. A genome whose contig set an earlier genome has is that genome (DEC-27, repeated_records.tsv);
a contig in several genomes of different sets stays in each (DEC-27; the genome with the least genome_uid writes its
calls); a record whose contig an owner holds is written to shared_records.tsv (owners as for v2, V2's inputs among
them). Writes gvmags_v2_all_18kgenomes.fna (every record, header and line width as in the archive, archive order),
genomes.json, objects.json (CC BY-SA 4.0, redistributable), records.tsv, repeated_records.tsv, shared_records.tsv,
skipped_members.tsv.

HC_ZIP is jgi_img_vr_v4.1_hc_20261118.zip, META_ZIP jgi_inovirus_umag_iuig_meta_20261118.zip. Nothing is unpacked:
the table and the FASTA are read as streams of their zip members (AGENTS §4.1: one record at a time; what is held is
one row per high-quality vOTU and one line per record taken). The md5 of each member read is compared with the
`sources.yaml` record (`jgi_img_vr` files); a difference is refused (exit 2) and nothing is written.

The selection (stated before it is applied; `imgvr_representatives`, a pure function of the table):
- the table is IMGVR_all_Sequence_information-high_confidence.tsv; every row must be `High-confidence` (it is the
  high-confidence file; another value is refused);
- a vOTU is taken when at least one of its rows has `MIUViG quality` = `High-quality` (estimated >= 90% complete,
  the MIUViG tier); `Reference` (RefSeq), `Genome fragment` and `Unsure` rows are never representatives;
- IMG/VR names no representative, so one is chosen among the vOTU's high-quality rows: the highest `Estimated
  completeness` capped at 100 (CheckV values above 100 say the sequence is longer than the genomes it was compared
  with, not more complete), then the lowest `Estimated contamination`, then the longest `Length`, then the least
  `UVIG`. A high-quality row whose completeness or contamination is not a number is refused. The result depends on
  the rows, not their order.
- Inovirus (Roux et al. 2019, Custom_MPI-Inovirus): every record of every GenBank file in Gb_files_inoviruses.zip, the
  files in name order.

The stage's genomes, kind virus (gene calling in gv-meta mode, `genecall.mode_of`), genetic code 11 (neither source
gives one; in gv-meta mode the code of each gene comes from the caller's models, so the value is only recorded),
assembly_quality `vOTU` (the schema's value for a vOTU representative; the MIUViG tier is kept in records.tsv):
- an IMG/VR representative is one genome of every FASTA record whose header's first `|` field is its UViG (a GVMAG is a
  bin of several records); the records' bases must sum to the table's `Length` (refused otherwise);
- an Inovirus record is one genome;
- circular (`is_circular`) only for a one-record genome whose topology is IMG/VR's `Direct terminal repeat` or
  GenBank's `circular` and whose sequence no longer carries the repeat (its first 20 bases are not found again in its
  last 2,000; card C-39's test, `c39_inputs.has_terminal_repeat`): the gene caller wraps a circular contig end to
  start, which is right only when the second copy is gone. Every other genome is linear.
No record is dropped (DEC-27): a genome whose contig set an earlier genome has is that genome, listed once, and written
to repeated_records.tsv. Owners (DEC-26): with --root, the default variants of the `contig_index` of each stage of
--root-stages (explicit part files, read only; never a stage being built), then each --owner-inputs NAME=DIR, the
contig uids of DIR/records.tsv (staged inputs of a stage built before this one, V1): a record whose contig one of them
holds stays in its genome and is written to shared_records.tsv with the first owner in that order.

Writes, under OUT_DIR.tmp renamed to OUT_DIR when complete (an existing OUT_DIR is refused): imgvr_v4.1_hc_hq_votu_reps.fna
(the records taken, header and line width as in the source, FASTA order), inovirus_roux2019.fna (each GenBank record
as FASTA: header its LOCUS name and family/subfamily, 60 bases a line, case as in the file), genomes.json, objects.json
(redistributable false, JGI's terms), records.tsv, repeated_records.tsv, shared_records.tsv. Prints a JSON report.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import sys
import zipfile
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # make_inputs and c39_inputs beside this script

import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

HC = "Cus_VR/IMG_VR_2022-12-19_7.1/"
SEQ_TSV = HC + "IMGVR_all_Sequence_information-high_confidence.tsv"
FNA = HC + "IMGVR_all_nucleotides-high_confidence.fna.gz"
GB_ZIP = "Cusrus/data/Gb_files_inoviruses.zip"
# sources.yaml, entries jgi_img_vr and jgi_inovirus (JGI API md5sum values)
MD5 = {"seq_tsv": "3c516db128082fa29dc2c2f60520da1b", "fna": "83301c9c6dfefea3305a53ee2a41bac3",
       "gb_zip": "a075958a353ec9a36c01264513f329e1"}
IMGVR_FASTA = "imgvr_v4.1_hc_hq_votu_reps.fna"
INOVIRUS_FASTA = "inovirus_roux2019.fna"
JGI_TERMS = ("JGI Data Utilization Policy (dataset marked Unrestricted on the JGI Data Portal; accepted by the PI "
             "2026-11-18); held local, not redistributed (DEC-50, DEC-51)")
IMGVR_SOURCE = {"source_name": "jgi_img_vr",
                "source_release": "IMG_VR_2022-12-19_7.1 (IMG/VR v4.1) high-confidence; high-quality vOTU "
                                  "representatives chosen by c48c_virus_inputs.py",
                "license": JGI_TERMS, "redistributable": False}
INOVIRUS_SOURCE = {"source_name": "jgi_inovirus",
                   "source_release": "Custom_MPI-Inovirus Gb_files_inoviruses.zip (Roux et al. 2019, "
                                     "doi:10.1038/s41564-019-0510-x)",
                   "license": JGI_TERMS, "redistributable": False}
NEEDED = ("UVIG", "vOTU", "Length", "Topology", "Confidence", "Estimated completeness", "Estimated contamination",
          "MIUViG quality", "Ecosystem classification", "Host taxonomy prediction", "Host prediction method",
          "Sequence origin (doi)")
KEPT = ("UVIG", "vOTU", "Length", "Topology", "MIUViG quality", "Estimated completeness", "Estimated contamination",
        "Ecosystem classification", "Host prediction method", "Host taxonomy prediction", "Sequence origin (doi)")
LINE = 60  # bases per line of the Inovirus FASTA


# --- parsing --------------------------------------------------------------------------------------------------------

def iter_fasta_stream(fh) -> Iterator[fasta.Record]:
    """`fasta.iter_records` over an open binary stream (a zip or tar member): the same records, line widths and
    errors, one record held at a time."""
    header: str | None = None
    seq_lines: list[bytes] = []
    for raw in fh:
        if not raw.strip():
            continue
        ln = raw.rstrip(b"\r\n")
        if ln.startswith(b">"):
            if header is not None:
                yield fasta.Record(header, b"".join(seq_lines), fasta._width(seq_lines))
            header = ln[1:].decode()
            seq_lines = []
        else:
            if header is None:
                raise ValueError("sequence before the first header")
            seq_lines.append(ln)
    if header is not None:
        yield fasta.Record(header, b"".join(seq_lines), fasta._width(seq_lines))


def iter_fasta_taken(fh, take) -> Iterator[tuple[str, fasta.Record | None]]:
    """(header, record) for each record of a FASTA stream, the record only when `take(header)` is true (the others'
    lines are skipped, not joined): the scan of a 82 Gbp file for the 4% it keeps."""
    header: str | None = None
    keep = False
    seq_lines: list[bytes] = []
    for raw in fh:
        if raw[:1] == b">":
            if header is not None:
                yield header, (fasta.Record(header, b"".join(seq_lines), fasta._width(seq_lines)) if keep else None)
            header = raw[1:].rstrip(b"\r\n").decode()
            keep = take(header)
            seq_lines = []
        elif keep:
            ln = raw.rstrip(b"\r\n")
            if ln:
                seq_lines.append(ln)
        elif header is None and raw.strip():
            raise ValueError("sequence before the first header")
    if header is not None:
        yield header, (fasta.Record(header, b"".join(seq_lines), fasta._width(seq_lines)) if keep else None)


_SEQ_JUNK = re.compile(rb"[\s0-9]")


def iter_genbank(fh) -> Iterator[dict]:
    """The records of a GenBank flat file (binary stream): LOCUS name, its topology word (linear, circular or ''),
    the LOCUS length, ACCESSION's first token, DEFINITION, and the ORIGIN sequence (digits and spaces removed, case
    kept). A record without ORIGIN or `//` is refused (ValueError)."""
    rec: dict | None = None
    in_seq = False
    key = ""
    for raw in fh:
        ln = raw.rstrip(b"\r\n")
        if ln.startswith(b"LOCUS"):
            if rec is not None:
                raise ValueError(f"record {rec['locus']}: no // before the next LOCUS")
            f = ln.split()
            topo = next((t.decode() for t in f if t in (b"linear", b"circular")), "")
            length = next((int(f[i - 1]) for i, t in enumerate(f) if t in (b"bp", b"aa") and i > 0
                           and f[i - 1].isdigit()), -1)
            rec = {"locus": f[1].decode() if len(f) > 1 else "", "topology": topo, "locus_length": length,
                   "accession": "", "definition": "", "seq": []}
            in_seq = False
            key = "LOCUS"
            continue
        if rec is None:
            if ln.strip():
                raise ValueError(f"text before the first LOCUS: {ln[:60]!r}")
            continue
        if ln.startswith(b"//"):
            if not in_seq:
                raise ValueError(f"record {rec['locus']}: no ORIGIN")
            rec["seq"] = b"".join(rec["seq"])
            yield rec
            rec = None
            continue
        if in_seq:
            rec["seq"].append(_SEQ_JUNK.sub(b"", ln))
            continue
        if ln.startswith(b"ORIGIN"):
            in_seq = True
            continue
        if ln[:1] not in (b" ", b""):
            key = ln[:12].strip().decode()
            val = ln[12:].strip().decode("utf-8", "replace")
            if key == "ACCESSION":
                rec["accession"] = val.split()[0] if val.split() else ""
            elif key == "DEFINITION":
                rec["definition"] = val
        elif key == "DEFINITION":
            rec["definition"] += " " + ln.strip().decode("utf-8", "replace")
    if rec is not None:
        raise ValueError(f"record {rec['locus']}: the file ends before //")


class HashReader(io.RawIOBase):
    """A read-through stream that hashes (md5) every byte read from `raw`."""

    def __init__(self, raw):
        self.raw, self.md5, self.n = raw, hashlib.md5(), 0

    def readable(self) -> bool:
        return True

    def readinto(self, b) -> int:
        data = self.raw.read(len(b))
        n = len(data)
        b[:n] = data
        self.md5.update(data)
        self.n += n
        return n

    def drain(self) -> str:
        """Read (and hash) whatever the reader above left unread; the md5 of every byte of `raw`."""
        while True:
            data = self.raw.read(1 << 22)
            if not data:
                break
            self.md5.update(data)
            self.n += len(data)
        return self.md5.hexdigest()


# --- the rules --------------------------------------------------------------------------------------------------------

def _number(r: dict, col: str) -> float:
    try:
        return float(r[col])
    except (TypeError, ValueError):
        raise Refused(f"{r.get('UVIG')}: {col} {r.get(col)!r} is not a number") from None


def rep_key(r: dict) -> tuple:
    """The order of a vOTU's high-quality rows; the least is its representative (module docstring)."""
    return (-min(_number(r, "Estimated completeness"), 100.0), _number(r, "Estimated contamination"),
            -int(r["Length"]), r["UVIG"])


def imgvr_representatives(rows: Iterable[dict]) -> tuple[dict[str, dict], dict]:
    """vOTU -> its representative's row (the KEPT columns) for every vOTU with a high-quality row, and counts."""
    best: dict[str, tuple[tuple, dict]] = {}
    stats = Counter()
    for r in rows:
        stats["rows"] += 1
        if r["Confidence"] != "High-confidence":
            raise Refused(f"{r['UVIG']}: Confidence {r['Confidence']!r} in the high-confidence table")
        stats[f"miuvig:{r['MIUViG quality']}"] += 1
        if r["MIUViG quality"] != "High-quality":
            continue
        stats["hq_rows"] += 1
        k = rep_key(r)
        b = best.get(r["vOTU"])
        if b is None or k < b[0]:
            best[r["vOTU"]] = (k, {c: r[c] for c in KEPT})
    reps = {v: best[v][1] for v in sorted(best)}
    return reps, {"rows": stats.pop("rows", 0), "hq_rows": stats.pop("hq_rows", 0), "hq_votus": len(reps),
                  "miuvig": {k.split(":", 1)[1]: n for k, n in sorted(stats.items())}}


def has_terminal_repeat(seq: bytes, k: int = 20, window: int = 2000) -> bool:
    """Card C-39's test (`c39_inputs.has_terminal_repeat`): the first `k` bases found again in the last `window`."""
    s = seq.upper()
    return len(s) >= 2 * k and s[:k] in s[max(k, len(s) - window):]


def is_circular(topology: str, seq: bytes) -> bool:
    """Listed circular: a direct terminal repeat (IMG/VR) or a circular GenBank record whose repeat was trimmed."""
    return topology in ("Direct terminal repeat", "circular") and not has_terminal_repeat(seq)


# --- owners (DEC-26) ------------------------------------------------------------------------------------------------

def _sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def owners_of(uids: list[bytes], root: Path | None, root_stages: Iterable[str],
              owner_inputs: dict[str, Path]) -> dict[bytes, str]:
    """The first owner, in the order given, of each of `uids` that a stage holds: the default variants of each root
    stage's contig_index (explicit part files), then each staged input's records.tsv. One DuckDB join per owner; only
    the matches come back."""
    import duckdb

    from biobox import io as bbio
    from biobox.l0seq import store

    con = duckdb.connect(config={"threads": 1, "memory_limit": "2GB"})
    con.execute("CREATE TABLE v (contig_uid BLOB)")
    con.executemany("INSERT INTO v VALUES (?)", [(u,) for u in sorted(set(uids))])
    out: dict[bytes, str] = {}
    for s in root_stages:
        files = bbio.part_files(Path(root) / store.INDEX_DIR / f"stage={s}") if root is not None else []
        if root is not None and not files:
            raise Refused(f"stage {s}: no contig_index under {root}")
        if not files:
            continue
        for (cu,) in con.execute(f"SELECT DISTINCT contig_uid FROM read_parquet({_sql_list(files)}, "
                                 "hive_partitioning = false) WHERE is_default_variant AND contig_uid IN "
                                 "(SELECT contig_uid FROM v) ORDER BY 1").fetchall():
            out.setdefault(bytes(cu), s)
    for name, d in owner_inputs.items():
        rec = Path(d) / "records.tsv"
        for (cu,) in con.execute(f"SELECT DISTINCT unhex(contig_uid) FROM read_csv({_sql_list([rec])}, delim='\t', "
                                 "header=true, quote='', all_varchar=true) WHERE unhex(contig_uid) IN "
                                 "(SELECT contig_uid FROM v) ORDER BY 1").fetchall():
            out.setdefault(bytes(cu), name)
    con.close()
    return out


# --- staging ------------------------------------------------------------------------------------------------------------

def _md5_check(rep: dict, name: str, got: str, expected: dict[str, str] | None) -> None:
    want = (expected or {}).get(name)
    rep[name] = {"got": got, "expected": want, "equal": want is None or got == want}
    if want is not None and got != want:
        raise Refused(f"md5 of {name} is {got}, sources.yaml says {want}")


def _imgvr(hc_zip: Path, work: Path, expected_md5: dict | None, rep: dict):
    """Pass 1 (the table) and pass 2 (the FASTA): the representatives' records written to work/IMGVR_FASTA."""
    with zipfile.ZipFile(hc_zip) as z:
        with z.open(SEQ_TSV) as raw:
            h = HashReader(raw)
            text = io.TextIOWrapper(io.BufferedReader(h, 1 << 20), encoding="utf-8", newline="")
            reader = csv.DictReader(text, delimiter="\t", quoting=csv.QUOTE_NONE)
            missing = [c for c in NEEDED if c not in (reader.fieldnames or ())]
            if missing:
                raise Refused(f"{SEQ_TSV}: no column {missing[0]!r}")
            reps, stats = imgvr_representatives(reader)
            _md5_check(rep["md5"], "seq_tsv", h.drain(), expected_md5)
        rep.update(stats)
        chosen = {r["UVIG"]: r for r in reps.values()}
        parts: dict[str, list[tuple[bytes, int, int]]] = {}  # uvig -> (contig uid, length, record ordinal)
        circ: dict[str, bool] = {}  # judged on a genome's first record; a second record (a bin) makes it linear
        n_recs = n_fasta = 0
        with z.open(FNA) as raw, open(work / IMGVR_FASTA, "wb") as w:
            h = HashReader(raw)
            gz = gzip.GzipFile(fileobj=io.BufferedReader(h, 1 << 22))
            for header, r in iter_fasta_taken(gz, lambda hd: _uvig(hd) in chosen):
                n_fasta += 1
                if r is None:
                    continue
                uvig = _uvig(header)
                w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width)) + b"\n")
                cu = ids.contig_uid(r.seq)
                if uvig in parts:
                    circ[uvig] = False
                else:
                    circ[uvig] = is_circular(chosen[uvig]["Topology"], r.seq)
                parts.setdefault(uvig, []).append((cu, len(r.seq), n_recs))
                n_recs += 1
            gz.close()
            _md5_check(rep["md5"], "fna", h.drain(), expected_md5)
        rep["imgvr_fasta_records"] = n_fasta
    missing = sorted(set(chosen) - set(parts))
    if missing:
        raise Refused(f"{len(missing)} representatives have no FASTA record; first {missing[0]}")
    bad = [u for u, p in parts.items() if sum(n for _, n, _ in p) != int(chosen[u]["Length"])]
    if bad:
        raise Refused(f"{len(bad)} representatives' records do not sum to the table's Length; first {bad[0]} "
                      f"({sum(n for _, n, _ in parts[bad[0]])} against {chosen[bad[0]]['Length']})")
    return chosen, parts, circ


def _uvig(header: str) -> str:
    """The UViG of an IMG/VR FASTA header: its first `|` field (`IMGVR_UViG_<taxon>_<n>|taxon|scaffold[|coords]`)."""
    return header.split("|", 1)[0].split(None, 1)[0] if header else ""


def _inovirus(meta_zip: Path, work: Path, expected_md5: dict | None, rep: dict) -> list[dict]:
    with zipfile.ZipFile(meta_zip) as z0:
        inner = z0.read(GB_ZIP)  # 40 MB held: seeking in a deflated member would re-inflate it from its start
    _md5_check(rep["md5"], "gb_zip", hashlib.md5(inner).hexdigest(), expected_md5)
    out, seen = [], set()
    with zipfile.ZipFile(io.BytesIO(inner)) as z, open(work / INOVIRUS_FASTA, "wb") as w:
        for info in sorted((i for i in z.infolist() if not i.is_dir()), key=lambda i: i.filename):
            parts = info.filename.split("/")
            group = "/".join(parts[:-1])
            with z.open(info) as fh:
                for r in iter_genbank(fh):
                    if r["locus"] in seen:
                        raise Refused(f"Inovirus LOCUS {r['locus']} twice")
                    if r["locus_length"] != len(r["seq"]):
                        raise Refused(f"Inovirus {r['locus']}: LOCUS says {r['locus_length']} bp, ORIGIN has "
                                      f"{len(r['seq'])}")
                    seen.add(r["locus"])
                    header = f"{r['locus']} {group}"
                    w.write(b"\n".join(fasta.format_record(header, r["seq"], LINE)) + b"\n")
                    out.append({"locus": r["locus"], "group": group, "topology": r["topology"],
                                "seq_len": len(r["seq"]), "cu": ids.contig_uid(r["seq"]),
                                "circular": is_circular(r["topology"], r["seq"])})
    return out


def stage_v2(hc_zip: Path, meta_zip: Path, out: Path, *, root: Path | None, root_stages: Iterable[str],
             owner_inputs: dict[str, Path], expected_md5: dict[str, str] | None = MD5,
             retrieved_imgvr: str, retrieved_inovirus: str) -> dict:
    out = Path(out)
    if out.exists():
        raise Refused(f"{out} exists")
    work = out.with_name(f"{out.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        rep = _stage_v2(hc_zip, meta_zip, out, work, root, list(root_stages), owner_inputs, expected_md5,
                        retrieved_imgvr, retrieved_inovirus)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    work.rename(out)
    names = ["genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
             IMGVR_FASTA, INOVIRUS_FASTA]
    rep["sha256"] = {n: ids.sha256_file(out / n).hex() for n in names}
    rep["bytes"] = {n: (out / n).stat().st_size for n in names}
    return rep


def _stage_v2(hc_zip, meta_zip, out, work, root, root_stages, owner_inputs, expected_md5, retrieved_imgvr,
              retrieved_inovirus) -> dict:
    rep: dict = {"md5": {}}
    chosen, parts, circ = _imgvr(Path(hc_zip), work, expected_md5, rep)
    ino = _inovirus(Path(meta_zip), work, expected_md5, rep)
    owners = owners_of([cu for p in parts.values() for cu, _, _ in p] + [x["cu"] for x in ino], root, root_stages,
                       owner_inputs)

    listing = make_inputs.Listing()
    genomes, repeated, shared, rows = [], [], [], []
    first_of: dict[bytes, str] = {}  # genome uid -> the genome that listed it
    ordinal = 0
    order = sorted(parts, key=lambda u: parts[u][0][2])  # FASTA order of each genome's first record
    for uvig in order:
        row = chosen[uvig]
        uids = [cu for cu, _, _ in parts[uvig]]
        for cu, n, _ in parts[uvig]:
            rows.append((ordinal, "imgvr", uvig, row["vOTU"], cu.hex(), n, row["MIUViG quality"],
                         row["Estimated completeness"], row["Estimated contamination"], row["Topology"], 11,
                         str(circ[uvig]).lower(), row["Ecosystem classification"], row["Host prediction method"],
                         row["Host taxonomy prediction"], row["Sequence origin (doi)"]))
            if cu in owners:
                shared.append((ordinal, uvig, cu.hex(), owners[cu]))
            ordinal += 1
        _add(listing, genomes, repeated, first_of, uids, circ[uvig], uvig, ordinal - len(uids))
    for x in ino:
        rows.append((ordinal, "inovirus", x["locus"], x["group"], x["cu"].hex(), x["seq_len"], "", "", "",
                     x["topology"], 11, str(x["circular"]).lower(), "", "", "", "Inoviruses (10.1038/s41564-019-0510-x)"))
        if x["cu"] in owners:
            shared.append((ordinal, x["locus"], x["cu"].hex(), owners[x["cu"]]))
        _add(listing, genomes, repeated, first_of, [x["cu"]], x["circular"], x["locus"], ordinal)
        ordinal += 1

    make_inputs._write_tsv(work / "records.tsv", [
        "ordinal", "source", "genome", "votu", "contig_uid", "length", "miuvig_quality", "completeness",
        "contamination", "topology", "genetic_code", "circular", "ecosystem", "host_method", "host_taxonomy",
        "origin"], rows)
    make_inputs._write_tsv(work / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"],
                           repeated)
    make_inputs._write_tsv(work / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"], shared)
    objects = [{"path": str(out / IMGVR_FASTA), **IMGVR_SOURCE, "retrieved_at": retrieved_imgvr},
               {"path": str(out / INOVIRUS_FASTA), **INOVIRUS_SOURCE, "retrieved_at": retrieved_inovirus}]
    (work / "objects.json").write_text(json.dumps(objects, indent=1))
    (work / "genomes.json").write_text(json.dumps(genomes))

    contig_len: dict[bytes, int] = {}
    for r in rows:
        contig_len.setdefault(bytes.fromhex(r[4]), r[5])
    by = lambda rs: dict(sorted(Counter(r[3] for r in rs).items()))  # noqa: E731
    rep.update({
        "imgvr_records": sum(len(p) for p in parts.values()), "imgvr_genomes": len(parts),
        "imgvr_multi_record_genomes": sum(1 for p in parts.values() if len(p) > 1),
        "imgvr_bases": sum(n for p in parts.values() for _, n, _ in p),
        "imgvr_topology": dict(Counter(chosen[u]["Topology"] for u in parts).most_common()),
        "imgvr_dtr_genomes": sum(1 for u in parts if chosen[u]["Topology"] == "Direct terminal repeat"),
        "imgvr_circular_genomes": sum(1 for u in parts if circ[u]),
        "imgvr_origin": dict(Counter(chosen[u]["Sequence origin (doi)"] for u in parts).most_common()),
        "imgvr_ecosystem_l1": dict(Counter(chosen[u]["Ecosystem classification"].split(";")[0]
                                           for u in parts).most_common()),
        "imgvr_host_method": dict(Counter(chosen[u]["Host prediction method"] for u in parts).most_common()),
        "inovirus_records": len(ino), "inovirus_bases": sum(x["seq_len"] for x in ino),
        "inovirus_circular_word": sum(1 for x in ino if x["topology"] == "circular"),
        "inovirus_circular": sum(1 for x in ino if x["circular"]),
        "genomes": len(genomes), "contigs": len(contig_len), "contig_bases": sum(contig_len.values()),
        "circular_genomes": sum(1 for g in genomes if g["circular"]),
        "repeated_genomes": len(repeated), "shared_records": len(shared), "owner_records": by(shared),
        "owner_contigs": by(list({r[2]: r for r in shared}.values())),
        "new_contigs": sum(1 for cu in contig_len if cu not in owners),
        "new_bases": sum(n for cu, n in contig_len.items() if cu not in owners),
        **listing.report(), "objects": len(objects)})
    return rep


def _add(listing, genomes, repeated, first_of, uids, circular, name, ordinal) -> None:
    g = listing.add(uids, 11, "virus", "vOTU")
    guid = ids.genome_uid(list(dict.fromkeys(uids)))
    if g is None:
        repeated.append((ordinal, name, guid.hex(), first_of[guid]))
        return
    g["circular"] = [uids[0].hex()] if circular else []
    genomes.append(g)
    first_of[guid] = name


# --- stage GV1 (card C-48c part 3) ------------------------------------------------------------------------------------

GV_MEMBER = "gvmags_v2/all_18kgenomes_fna.tar.gz"
GV_MD5 = "077f125d844c85ca6169c331d52640b9"  # sources.yaml, gvmags_v2 (Zenodo's md5)
GV_FASTA = "gvmags_v2_all_18kgenomes.fna"
GV_QUALITY = "contig"  # unknown completeness: the archive has no quality table (part 3's decision 3; card C-40's value)
GV_SOURCE = {"source_name": "gvmags_v2",
             "source_release": "zenodo:10.5281/zenodo.20126600 GVMAGs V2, all_18kgenomes_fna.tar.gz (every genome)",
             "license": "CC BY-SA 4.0", "redistributable": True}


def stage_gv1(gv_tar: Path, out: Path, *, root: Path | None, root_stages: Iterable[str],
              owner_inputs: dict[str, Path], expected_md5: str | None = GV_MD5, retrieved: str) -> dict:
    """GV1's inputs (`gv1` mode, module docstring): every `.fna` member of GV_MEMBER one genome, under OUT.tmp, renamed
    to OUT when complete; an existing OUT is refused."""
    out = Path(out)
    if out.exists():
        raise Refused(f"{out} exists")
    work = out.with_name(f"{out.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        rep = _stage_gv1(Path(gv_tar), out, work, root, list(root_stages), owner_inputs, expected_md5, retrieved)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    work.rename(out)
    names = ["genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
             "skipped_members.tsv", GV_FASTA]
    rep["sha256"] = {n: ids.sha256_file(out / n).hex() for n in names}
    rep["bytes"] = {n: (out / n).stat().st_size for n in names}
    return rep


def _stage_gv1(gv_tar, out, work, root, root_stages, owner_inputs, expected_md5, retrieved) -> dict:
    import tarfile

    genomes_in: list[tuple[str, str, list[tuple[bytes, int, str]]]] = []  # (genome, member, [(uid, len, header)])
    skipped: list[tuple] = []
    rep: dict = {"md5": {}}
    with tarfile.open(gv_tar, mode="r:") as outer:
        raw = outer.extractfile(GV_MEMBER)
        if raw is None:
            raise Refused(f"{GV_MEMBER}: not a file in {gv_tar}")
        h = HashReader(raw)
        with tarfile.open(fileobj=io.BufferedReader(h, 1 << 22), mode="r|gz") as tf, \
                open(work / GV_FASTA, "wb") as w:
            for m in tf:
                if not m.isfile():
                    continue
                if not m.name.endswith(".fna"):
                    skipped.append((m.name, m.size, "not_fasta"))  # read as data, never run
                    continue
                g = m.name.rsplit("/", 1)[-1][:-4]
                recs = []
                for r in iter_fasta_stream(io.BytesIO(tf.extractfile(m).read())):  # one genome (a few MB) at a time
                    w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width)) + b"\n")
                    recs.append((ids.contig_uid(r.seq), len(r.seq), r.accession))
                if not recs:
                    raise Refused(f"{m.name}: no FASTA record")
                genomes_in.append((g, m.name, recs))
        _md5_check(rep["md5"], "gv_member", h.drain(), {"gv_member": expected_md5} if expected_md5 else None)
    names = [g for g, _, _ in genomes_in]
    if len(set(names)) != len(names):
        raise Refused("a genome name is in two members")
    owners = owners_of([u for _, _, rs in genomes_in for u, _, _ in rs], root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes, repeated, shared, rows = [], [], [], []
    first_of: dict[bytes, str] = {}
    held_whole = 0
    ordinal = 0
    for g, member, recs in genomes_in:
        uids = [u for u, _, _ in recs]
        for k, (u, n, acc) in enumerate(recs):
            rows.append((ordinal, g, k, acc, u.hex(), n, member))
            if u in owners:
                shared.append((ordinal, g, k, u.hex(), owners[u]))
            ordinal += 1
        guid = ids.genome_uid(list(dict.fromkeys(uids)))
        gi = listing.add(uids, 11, "virus", GV_QUALITY)
        if gi is None:
            repeated.append((ordinal - len(recs), g, guid.hex(), first_of[guid]))
            continue
        genomes.append(gi)  # every contig linear: the archive records no topology (part 3's decision 4)
        first_of[guid] = g
        held_whole += all(u in owners for u in uids)

    make_inputs._write_tsv(work / "records.tsv", ["ordinal", "genome", "record", "accession", "contig_uid", "length",
                                                  "member"], rows)
    make_inputs._write_tsv(work / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"],
                           repeated)
    make_inputs._write_tsv(work / "shared_records.tsv", ["ordinal", "genome", "record", "contig_uid", "owner_stage"],
                           shared)
    make_inputs._write_tsv(work / "skipped_members.tsv", ["member", "bytes", "reason"], skipped)
    objects = [{"path": str(out / GV_FASTA), **GV_SOURCE, "retrieved_at": retrieved}]
    (work / "objects.json").write_text(json.dumps(objects, indent=1))
    (work / "genomes.json").write_text(json.dumps(genomes))
    contig_len: dict[bytes, int] = {}
    for r in rows:
        contig_len.setdefault(bytes.fromhex(r[4]), r[5])
    by = lambda rs: dict(sorted(Counter(r[4] for r in rs).items()))  # noqa: E731
    rep.update({
        "fna_members": len(genomes_in), "skipped_members": len(skipped), "records": len(rows),
        "record_bases": sum(r[5] for r in rows), "genomes": len(genomes), "repeated_genomes": len(repeated),
        "contigs": len(contig_len), "contig_bases": sum(contig_len.values()),
        "contigs_per_genome_median": sorted(len(g["contig_uids"]) for g in genomes)[len(genomes) // 2] if genomes else 0,
        "genomes_all_contigs_held_earlier": held_whole,
        "shared_records": len(shared), "owner_records": by(shared),
        "owner_contigs": by(list({r[3]: r for r in shared}.values())),
        "new_contigs": sum(1 for u in contig_len if u not in owners),
        "new_bases": sum(n for u, n in contig_len.items() if u not in owners),
        **listing.report(), "objects": len(objects)})
    return rep


def main(argv: list[str]) -> int:
    try:
        argv, root = make_inputs._option(argv, "--root")
        argv, stages = make_inputs._option(argv, "--root-stages")
        owner_inputs = {}
        while "--owner-inputs" in argv:
            argv, spec = make_inputs._option(argv, "--owner-inputs")
            name, _, d = spec.partition("=")
            if not name or not d:
                raise Refused(f"--owner-inputs {spec!r}: expected NAME=DIR")
            owner_inputs[name] = Path(d)
    except Refused as e:
        print(f"c48c_virus_inputs.py: refused: {e}\n{__doc__.split(chr(10) * 2)[1]}", file=sys.stderr)
        return 2
    ok = (len(argv) == 4 and argv[:1] == ["v2"]) or (len(argv) == 3 and argv[:1] == ["gv1"])
    if not ok or (stages and not root):
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    common = {"root": None if root is None else Path(root), "root_stages": stages.split(",") if stages else (),
              "owner_inputs": owner_inputs}
    try:
        if argv[0] == "gv1":  # retrieved: the end of the Zenodo fetch, sources.yaml gvmags_v2
            rep = stage_gv1(Path(argv[1]), Path(argv[2]), retrieved="2026-11-19T10:42:14Z", **common)
        else:
            rep = stage_v2(Path(argv[1]), Path(argv[2]), Path(argv[3]), retrieved_imgvr="2026-11-18T23:15:37Z",
                           retrieved_inovirus="2026-11-18T23:16:08Z", **common)
    except Refused as e:
        print(f"c48c_virus_inputs.py: refused: {e}", file=sys.stderr)
        return 2
    print(json.dumps(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

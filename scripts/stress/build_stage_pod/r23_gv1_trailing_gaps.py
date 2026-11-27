"""Card R-23: `GV1` re-staged with the trailing gap bytes of its three contigs `Hoko_1`, `Hoko_3` and `Hoko_4` removed.

    r23_gv1_trailing_gaps.py restage OLD_DIR OUT_META OUT_FASTA FINAL_DIR [--expect ACC ...]
                                     [--root ROOT] [--chain-home DIR] [--work DIR]

Card R-19's scan found the last two bytes of three `GV1` contigs to be `--` (gap bytes of the source's FASTA), and card
R-20 read the ingest path: ingest stores a `-` as a base, counts it in the length and hashes it into the `contig_uid`.
This re-stage reads card C-48c's staged `GV1` in OLD_DIR (the FASTA its `objects.json` names, `records.tsv`,
`genomes.json`, `repeated_records.tsv`, `shared_records.tsv`, `skipped_members.tsv`) and writes it again under a new
FASTA name, OUT_FASTA, and OUT_META; OLD_DIR is only read.

- The FASTA is read record by record, each with its own lines and line ends. A record whose sequence holds a byte that
  is not an IUPAC nucleotide letter (either case) is taken only when those bytes are a run of `-` at its end and the
  rest is IUPAC letters: the run is removed from the end of its last lines (a line left empty goes with its line end).
  Any other such byte refuses the run, as does a set of trimmed records other than `--expect` (default the three).
  Every other record is written byte for byte.
- `records.tsv` holds one row per FASTA record in file order (C-48c's layout); each row is checked against its record
  (ordinal, accession, `contig_uid` of the bytes, length), so a table out of step with the file refuses the run. A
  trimmed record's row gets the uid and length of the trimmed bytes; every other row is copied as it is.
- `genomes.json` (C-48c's `json.dumps` form, no final newline, checked) has each old uid of a trimmed record replaced
  by the new one. Refused: a trimmed contig that becomes equal to another contig of its own genome, a trimmed record's
  genome named in `repeated_records.tsv`, an old uid of a trimmed record in `shared_records.tsv`. A new uid equal to a
  contig of another genome of the stage is reported (`new_uids_already_in_stage`; DEC-27 lists it in both). A trimmed
  genome whose contig set becomes an earlier genome's is, by C-48c's own DEC-27 rule, that genome: it leaves
  `genomes.json` and gets a row in `repeated_records.tsv` (`genomes_now_repeats`); its records stay in the FASTA and in
  `records.tsv`, as C-48c kept every repeat's records.
- `shared_records.tsv` and `skipped_members.tsv` are copied byte for byte, and `repeated_records.tsv` too unless a
  genome became a repeat (its row is then inserted in ordinal order, every old row as it was); `gap_trims.tsv` lists
  each trimmed record (old and new uid and length, the bytes removed); `objects.json` is C-48c's object with its path
  FINAL_DIR/<OUT_FASTA's name> and `objects_staging.json` the same with OUT_FASTA (both in C-48c's `indent=1` form).
- The new file is read back beside the old one: every header and every untrimmed record byte-identical, each trimmed
  record its old bytes less the gaps, every record's `contig_uid` and length equal to the new `records.tsv`.
  `removed_ranges` gives the removed bytes as (offset in the old file, length), so `cmp` can check the rest.
- With `--root` and/or `--chain-home`, the exact `contig_uid` ownership of the old and the new `records.tsv` over the
  root's stored stages and the `r26.11.4` chain (card C-76's `owners`, DEC-26), and which new uids an earlier stage
  holds (`trimmed_uids_held`).

Writes `stage_report.json` in OUT_META. Exit 0, or 2 when refused (nothing is left at OUT_FASTA).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

IUPAC = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn"
GAP = b"-"
EXPECT = ("Hoko_1", "Hoko_3", "Hoko_4")
RECORDS_HEAD = [b"ordinal", b"genome", b"record", b"accession", b"contig_uid", b"length", b"member"]
COPIED = ("repeated_records.tsv", "shared_records.tsv", "skipped_members.tsv")
GAP_TRIMS_HEAD = ["ordinal", "genome", "accession", "old_contig_uid", "old_length", "removed", "new_contig_uid",
                  "new_length"]


class Refused(Exception):
    pass


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _record_lines(fh):
    """(header line, [body lines]) of a FASTA stream, every line with its own line end, blank lines kept."""
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


def _content(ln: bytes) -> bytes:
    return ln.rstrip(b"\r\n")


def seq_of(body: list[bytes]) -> bytes:
    """The sequence as `biobox.fasta.iter_records` reads it: blank lines skipped, CR/LF stripped."""
    return b"".join(_content(x) for x in body if x.strip())


def trim_trailing_gaps(body: list[bytes], start: int) -> tuple[list[bytes], int, list[list[int]]]:
    """BODY less the run of `-` at the end of its sequence: (new lines, bytes of sequence removed, removed ranges as
    [offset, length] with START the offset of the body's first byte). A line left empty goes with its line end."""
    offs, at = [], start
    for x in body:
        offs.append(at)
        at += len(x)
    out = list(body)
    ranges: list[list[int]] = []
    removed = 0
    j = len(out) - 1
    while j >= 0:
        ln = out[j]
        if not ln.strip():
            j -= 1
            continue
        c = _content(ln)
        n = len(c) - len(c.rstrip(GAP))
        if n == 0:
            break
        removed += n
        if n == len(c):
            ranges.append([offs[j], len(ln)])
            del out[j]
            j -= 1
            continue
        ranges.append([offs[j] + len(c) - n, n])
        out[j] = c[: len(c) - n] + ln[len(c):]
        break
    return out, removed, sorted(ranges)


def _copy_bytes(src: Path, dst: Path) -> None:
    with open(src, "rb") as fi, open(dst, "wb") as fo:
        for b in iter(lambda: fi.read(1 << 22), b""):
            fo.write(b)


def _tsv_rows(path: Path) -> tuple[list[bytes], list[list[bytes]]]:
    with open(path, "rb") as fh:
        head = fh.readline().rstrip(b"\r\n").split(b"\t")
        return head, [ln.rstrip(b"\r\n").split(b"\t") for ln in fh if ln.strip()]


def _stats(records_tsv: Path, genomes: list[dict]) -> dict:
    contigs: dict[bytes, int] = {}
    rows = bases = 0
    with open(records_tsv, "rb") as fh:
        fh.readline()
        for ln in fh:
            f = ln.rstrip(b"\r\n").split(b"\t")
            rows += 1
            bases += int(f[5])
            contigs.setdefault(f[4], int(f[5]))
    return {"records": rows, "record_bases": bases, "contigs": len(contigs), "contig_bases": sum(contigs.values()),
            "genomes": len(genomes), "memberships": sum(len(g["contig_uids"]) for g in genomes)}


def restage(old_dir, out_meta, out_fasta, final_dir, *, expect=EXPECT, root: Path | None = None,
            chain_home: Path | None = None, work: Path | None = None) -> dict:
    """GV1 with the trailing gaps of `expect` removed, into OUT_FASTA (a new name) and OUT_META (module docstring)."""
    from biobox import ids

    old, out_meta, out_fasta = Path(old_dir), Path(out_meta), Path(out_fasta)
    if out_fasta.exists():
        raise Refused(f"{out_fasta} exists")
    if out_meta.exists() and any(out_meta.iterdir()):
        raise Refused(f"{out_meta} exists and is not empty")
    for n in ("records.tsv", "genomes.json", "objects.json", *COPIED):
        if not (old / n).is_file():
            raise Refused(f"missing {old / n}")
    objs = json.loads((old / "objects.json").read_text(encoding="utf-8"))
    if len(objs) != 1:
        raise Refused(f"{old / 'objects.json'}: expected one object")
    old_fasta = Path(objs[0]["path"])
    if not old_fasta.is_file():
        raise Refused(f"missing {old_fasta}")
    if old_fasta.resolve() == out_fasta.resolve():
        raise Refused("the new FASTA must have a new name")
    gj_raw = (old / "genomes.json").read_bytes()
    gj = json.loads(gj_raw)
    if json.dumps(gj).encode() != gj_raw:
        raise Refused(f"{old / 'genomes.json'} is not in C-48c's json.dumps form")
    old_stats_genomes = [dict(e) for e in gj]

    out_meta.mkdir(parents=True, exist_ok=True)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    tmp_fasta = out_fasta.with_name(out_fasta.name + ".tmp")
    part = {n: out_meta / f".{n}.partial" for n in ("records.tsv", "genomes.json", "gap_trims.tsv", *COPIED,
                                                     "objects.json", "objects_staging.json")}

    trims: list[dict] = []  # one per trimmed record
    ranges: list[list[int]] = []
    all_uids: set[str] = set()
    order: list[str] = []  # genome names in records.tsv order (C-48c's archive order)
    first_ordinal: dict[str, int] = {}
    n_rec = 0
    try:
        with open(old_fasta, "rb") as fi, open(tmp_fasta, "wb") as fo, open(old / "records.tsv", "rb") as ri, \
                open(part["records.tsv"], "wb") as ro:
            head_ln = ri.readline()
            if head_ln.rstrip(b"\r\n").split(b"\t") != RECORDS_HEAD:
                raise Refused(f"{old / 'records.tsv'}: unexpected header {head_ln!r}")
            ro.write(head_ln)
            at = 0  # offset in the old file
            for hl, body in _record_lines(fi):
                n_rec += 1
                body_at = at + len(hl)
                at = body_at + sum(len(x) for x in body)
                try:
                    acc = hl[1:].rstrip(b"\r\n").decode("utf-8").split(None, 1)[0] if hl[1:].strip() else ""
                except UnicodeDecodeError:
                    raise Refused(f"{old_fasta}: record {n_rec}: header is not UTF-8") from None
                seq = seq_of(body)
                row_ln = ri.readline()
                if not row_ln:
                    raise Refused(f"{old / 'records.tsv'} has fewer rows than {old_fasta} has records")
                row = row_ln.rstrip(b"\r\n").split(b"\t")
                u_old = ids.contig_uid(seq).hex()
                if (row[0] != str(n_rec - 1).encode() or row[3].decode() != acc or row[4].decode() != u_old
                        or int(row[5]) != len(seq)):
                    raise Refused(f"{old_fasta}: record {n_rec} ({acc}) is not records.tsv row {n_rec} "
                                  f"({b' '.join(row[:6]).decode(errors='replace')})")
                all_uids.add(u_old)
                gname = row[1].decode()
                if gname not in first_ordinal:
                    first_ordinal[gname] = n_rec - 1
                    order.append(gname)
                if seq.translate(None, IUPAC):
                    core = seq.rstrip(GAP)
                    if not core or core.translate(None, IUPAC):
                        raise Refused(f"{old_fasta}: record {n_rec} ({acc}): a byte that is not a nucleotide letter "
                                      "and not only trailing gap bytes")
                    new_body, k, rr = trim_trailing_gaps(body, body_at)
                    if seq_of(new_body) != core or k != len(seq) - len(core):
                        raise Refused(f"{old_fasta}: record {n_rec} ({acc}): the gaps could not be cut from its lines")
                    u_new = ids.contig_uid(core).hex()
                    trims.append({"ordinal": n_rec - 1, "genome": row[1].decode(), "accession": acc,
                                  "old_contig_uid": u_old, "old_length": len(seq), "removed": seq[len(core):].decode(),
                                  "new_contig_uid": u_new, "new_length": len(core), "ranges": rr})
                    ranges.extend(rr)
                    fo.write(hl)
                    for x in new_body:
                        fo.write(x)
                    row[4], row[5] = u_new.encode(), str(len(core)).encode()
                    ro.write(b"\t".join(row) + row_ln[len(row_ln.rstrip(b"\r\n")):])
                else:
                    fo.write(hl)
                    for x in body:
                        fo.write(x)
                    ro.write(row_ln)
            if ri.readline().strip():
                raise Refused(f"{old / 'records.tsv'} has more rows than {old_fasta} has records")

        got = tuple(sorted(t["accession"] for t in trims))
        if got != tuple(sorted(expect)):
            raise Refused(f"trimmed records {list(got)}, expected {sorted(expect)}")
        old_to_new = {t["old_contig_uid"]: t["new_contig_uid"] for t in trims}
        trimmed_genomes = {t["genome"] for t in trims}

        # the tables the re-stage copies must not name what it changes
        _, rep_rows = _tsv_rows(old / "repeated_records.tsv")
        for r in rep_rows:
            if {r[1].decode(), r[3].decode()} & trimmed_genomes:
                raise Refused(f"genome {sorted({r[1].decode(), r[3].decode()} & trimmed_genomes)[0]} of a trimmed "
                              "record is named in repeated_records.tsv")
        sh_head, sh_rows = _tsv_rows(old / "shared_records.tsv")
        ci = sh_head.index(b"contig_uid")
        for r in sh_rows:
            if r[ci].decode() in old_to_new:
                raise Refused(f"an old uid of a trimmed record ({r[ci].decode()}) is in shared_records.tsv")

        # genomes.json: the old uids replaced; DEC-27's genome uids stay distinct
        def guid(e):
            return ids.genome_uid([bytes.fromhex(u) for u in e["contig_uids"]])

        before = Counter(guid(e) for e in gj)
        if max(before.values(), default=1) > 1:
            raise Refused(f"{old / 'genomes.json'} lists one contig set twice")
        # genomes.json lists the genomes of records.tsv in its order, less the repeats (C-48c's Listing)
        repeated_names = {r[1].decode() for r in rep_rows}
        listed = [g for g in order if g not in repeated_names]
        if len(listed) != len(gj):
            raise Refused(f"{old / 'genomes.json'} has {len(gj)} genomes, records.tsv less its repeats {len(listed)}")
        found: set[str] = set()
        trim_by_old = {t["old_contig_uid"]: t for t in trims}
        for name, e in zip(listed, gj):
            if not any(u in old_to_new for u in e["contig_uids"]):
                continue
            if any(trim_by_old[u]["genome"] != name for u in e["contig_uids"] if u in trim_by_old):
                raise Refused(f"{old / 'genomes.json'} does not follow records.tsv at genome {name}")
            new_list = [old_to_new.get(u, u) for u in e["contig_uids"]]
            if len(set(new_list)) != len(new_list):
                dup = next(t for u, t in trim_by_old.items() if u in e["contig_uids"]
                           and new_list.count(t["new_contig_uid"]) > 1)
                raise Refused(f"{dup['accession']} without its gaps is listed twice in genome {dup['genome']}")
            found |= {u for u in e["contig_uids"] if u in old_to_new}
            e["contig_uids"] = new_list
        if found != set(old_to_new):
            raise Refused(f"{len(set(old_to_new) - found)} trimmed contigs are in no genome of genomes.json")
        # DEC-27 as C-48c's stager applies it: a genome whose contig set an earlier genome has is that genome; it leaves
        # genomes.json for repeated_records.tsv (its records stay in the FASTA and records.tsv)
        first_of: dict[bytes, str] = {}
        kept, new_repeats = [], []
        for name, e in zip(listed, gj):
            g = guid(e)
            if g in first_of:
                if name not in trimmed_genomes and first_of[g] not in trimmed_genomes:
                    raise Refused(f"genomes {first_of[g]} and {name} share a contig set in {old / 'genomes.json'}")
                new_repeats.append((first_ordinal[name], name, g.hex(), first_of[g]))
                continue
            first_of[g] = name
            kept.append(e)
        gj = kept
        part["genomes.json"].write_bytes(json.dumps(gj).encode())
        with open(old / "repeated_records.tsv", "rb") as fh:
            rep_lines = fh.readlines()
        rows_out = [(int(ln.split(b"\t", 1)[0]), ln) for ln in rep_lines[1:] if ln.strip()]
        rows_out += [(o, f"{o}\t{n}\t{g}\t{f}\n".encode()) for o, n, g, f in new_repeats]
        rows_out.sort(key=lambda x: x[0])
        new_in_stage = sorted({t["new_contig_uid"] for t in trims} & all_uids)

        with open(part["gap_trims.tsv"], "w", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(GAP_TRIMS_HEAD) + "\n")
            for t in trims:
                fh.write("\t".join(str(t[k]) for k in GAP_TRIMS_HEAD) + "\n")
        for n in COPIED:
            if n == "repeated_records.tsv" and new_repeats:
                part[n].write_bytes(rep_lines[0] + b"".join(ln for _, ln in rows_out))
            else:
                _copy_bytes(old / n, part[n])
        o = dict(objs[0])
        part["objects.json"].write_text(json.dumps([{**o, "path": str(Path(final_dir) / out_fasta.name)}], indent=1),
                                        encoding="utf-8")
        part["objects_staging.json"].write_text(json.dumps([{**o, "path": str(out_fasta)}], indent=1),
                                                encoding="utf-8")

        # read the new file back beside the old one, and against the new records.tsv
        trim_at = {t["ordinal"] for t in trims}
        mismatches = identical_other = n_back = 0
        same = True
        with open(old_fasta, "rb") as fa, open(tmp_fasta, "rb") as fb, open(part["records.tsv"], "rb") as rr_:
            rr_.readline()
            for (ha, ba), (hb, bb) in zip(_record_lines(fa), _record_lines(fb), strict=True):
                k = n_back
                n_back += 1
                sb = seq_of(bb)
                row = rr_.readline().rstrip(b"\r\n").split(b"\t")
                if ids.contig_uid(sb).hex() != row[4].decode() or len(sb) != int(row[5]):
                    mismatches += 1
                if ha != hb:
                    same = False
                elif k in trim_at:
                    sa = seq_of(ba)
                    if sa.rstrip(GAP) != sb or sb.translate(None, IUPAC):
                        same = False
                elif ba == bb:
                    identical_other += 1
                else:
                    same = False
    except BaseException:
        tmp_fasta.unlink(missing_ok=True)
        for p in part.values():
            p.unlink(missing_ok=True)
        raise
    os.replace(tmp_fasta, out_fasta)
    for n, p in part.items():
        p.replace(out_meta / n)

    rep: dict = {
        "card": "R-23", "stage": "GV1",
        "old_dir": str(old), "old_fasta_path": str(old_fasta), "old_fasta_bytes": old_fasta.stat().st_size,
        "old_fasta_sha256": _sha256(old_fasta),
        "fasta_path": str(out_fasta), "fasta_bytes": out_fasta.stat().st_size, "fasta_sha256": _sha256(out_fasta),
        "fasta_bytes_removed": old_fasta.stat().st_size - out_fasta.stat().st_size,
        "trimmed_records": len(trims), "removed_gap_bytes": sum(t["old_length"] - t["new_length"] for t in trims),
        "trims": [{k: v for k, v in t.items()} for t in trims],
        "removed_ranges": sorted(ranges),
        "records_read_back": n_back, "other_records_byte_identical": same and identical_other == n_back - len(trims),
        "uid_mismatches": mismatches, "new_uids_already_in_stage": new_in_stage,
        "genomes_now_repeats": [{"ordinal": o, "genome": n, "genome_uid": g, "first_genome": f}
                                for o, n, g, f in new_repeats],
        "old": _stats(old / "records.tsv", old_stats_genomes), "new": _stats(out_meta / "records.tsv", gj),
        "files_sha256": {p.name: _sha256(p) for p in sorted(out_meta.iterdir()) if p.is_file()
                         and not p.name.startswith(".") and p.name != "stage_report.json"},
    }
    if root is not None or chain_home is not None:
        import c76_lr2_sag1_nx1_r232 as c76

        w = Path(work) if work else out_meta.parent / "_work_r23"
        rep["ownership_old"] = c76.owners(old / "records.tsv", root, chain_home, {}, w / "old")
        rep["ownership_new"] = c76.owners(out_meta / "records.tsv", root, chain_home, {}, w / "new")
        rep["ownership_stages"] = rep["ownership_new"]["owner_stages"]
        import a26_r26_11_4_chain as a26
        import c62_metavr_inputs as c62

        chain = {s: chain_home / s for s, _ in a26.R26_11_4_CHAIN} if chain_home is not None else {}
        own = c62._owners_from_records(_uids_tsv(trims, w), root, list(c76.ROOT_STAGES) if root is not None else [],
                                       chain, work=w / "t")
        rep["trimmed_uids_held"] = {u.hex(): s for u, s in sorted(own.items())}
    (out_meta / "stage_report.json").write_text(json.dumps(rep, indent=1) + "\n", encoding="utf-8")
    return rep


def _uids_tsv(trims: list[dict], work: Path) -> Path:
    """A one-column `contig_uid` table of the trimmed records' new uids (for the owners query)."""
    work.mkdir(parents=True, exist_ok=True)
    p = work / "trimmed_new_uids.tsv"
    p.write_text("contig_uid\n" + "".join(t["new_contig_uid"] + "\n" for t in trims), encoding="utf-8")
    return p


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card R-23: GV1 re-staged without the trailing gaps of three contigs")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("restage")
    for n in ("old_dir", "out_meta", "out_fasta", "final_dir"):
        s.add_argument(n, type=Path)
    s.add_argument("--expect", action="append", default=None, help="an accession expected to be trimmed (repeat)")
    s.add_argument("--root", type=Path, default=None)
    s.add_argument("--chain-home", type=Path, default=None)
    s.add_argument("--work", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        rep = restage(a.old_dir, a.out_meta, a.out_fasta, a.final_dir, expect=tuple(a.expect) if a.expect else EXPECT,
                      root=a.root, chain_home=a.chain_home, work=a.work)
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(rep, indent=1))
    return 0 if rep["other_records_byte_identical"] and rep["uid_mismatches"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

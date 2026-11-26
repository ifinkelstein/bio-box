"""Card C-76 (DEC-79 (3)): the R232 / r232 re-diff of the staged `LR2`, `SAG1` and `NX1`.

    c76_lr2_sag1_nx1_r232.py run STAGE META_DIR FASTA INDEX_DIR OUT_META_DIR OUT_FASTA_DIR
                             [--names ENA_TABLE] [--inputs-home DIR] [--root ROOT --chain-home DIR]
                             [--later STAGE=DIR ...] [--no-check] [--report-out FILE]

Cards C-66 and C-58 staged the three stages against GTDB R226 and GlobDB r226. DEC-79 (3) moves every "absent from"
test of `r26.11.5` to GTDB R232 and GlobDB r232. The test is by accession (no r232 genome is on the pod): each genome's
INSDC accession is read from its genome id (`ehomd:GCA_000406805.1`) or, for single cells, from an ENA assembly table
by the assembly's name (`--names`), and card C-74's `lookup` (`r232_index.py`) says whether GTDB R232 or GlobDB r232
holds it. A genome it finds is removed and listed with its reason. A genome with no accession (a MAG named by its
submitter, an unbinned contig) cannot be seen by this test and stays; the report counts both.

Nothing is overwritten. An unchanged stage keeps its files byte for byte: only `r232_rediff.json`, `r232_rediff.tsv`
and `r232_removed.tsv` are written. A changed stage is written under new names beside the old ones: its metadata
under OUT_META_DIR/STAGE (genomes.json, records.tsv, shared_records.tsv, repeated_records.tsv, trims.tsv,
genetic_codes.tsv, objects.json, objects_staging.json) and its FASTA as `<stem>_r232.fna.gz` in OUT_FASTA_DIR
(gzip, mtime 0, level 6, so a second run gives the same bytes: AGENTS 4.3). Records, contigs and FASTA records keep
their order; a FASTA record is identified by the contig_uid of its sequence, not by its header (SAG1's headers are
the assemblers' names); a contig that a surviving genome still lists stays (DEC-27); `ordinal` columns are
renumbered.

With `--root` and `--chain-home` the final stage's exact `contig_uid` ownership is computed, as card C-79 did, over
the root's stored stages and the 19 `r26.11.4` stages, and over the `r26.11.5` stages named by `--later`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from blake3 import blake3  # noqa: E402

import r232_index  # noqa: E402

ACC = re.compile(r"GC[AF]_\d{9}")
GENOME_COLS = ("genome_id", "genome")
ROOT_STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
TSV_FILES = ("shared_records.tsv", "repeated_records.tsv", "trims.tsv", "genetic_codes.tsv")


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
    """`lookup(accessions)` of card C-74 over its index, loading it once."""
    return lambda accs: r232_index.lookup(index_dir, accs)


# ---- the records -----------------------------------------------------------------------------------------------

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
            bare = g.split(":", 1)[-1]
            acc = names.get(g) or names.get(bare)
            how = "ena_name"
        out.append({"genome_id": g, "source_name": r[si] if si is not None else "-", "accession": acc,
                    "how": how if acc else None})
    return out


def rediff(records_tsv, lookup, names: dict[str, str] | None = None) -> dict:
    """The genomes an accession test finds in GTDB R232 or GlobDB r232 (`removed`, with the reason) and a count of
    what the test could and could not see."""
    gens = genome_accessions(records_tsv, names)
    with_acc = [g for g in gens if g["accession"]]
    hits = {}
    if with_acc:
        got = lookup(sorted({g["accession"] for g in with_acc}))
        hits = {h["accession"]: h for h in got}
    removed, by_source = [], {}
    for g in gens:
        s = by_source.setdefault(g["source_name"], {"genomes": 0, "with_accession": 0, "removed": 0})
        s["genomes"] += 1
        if not g["accession"]:
            continue
        s["with_accession"] += 1
        h = hits[g["accession"]]
        if h["in_gtdb_r232"] or h["in_globdb_r232"]:
            s["removed"] += 1
            reason = "+".join(k for k, v in (("in_gtdb_r232", h["in_gtdb_r232"]),
                                             ("in_globdb_r232", h["in_globdb_r232"])) if v)
            removed.append({**g, "reason": reason, "in_gtdb_r232": h["in_gtdb_r232"],
                            "in_globdb_r232": h["in_globdb_r232"], "representative": h["representative"],
                            "species": h["species"], "globdb_r232_id": h["globdb_r232_id"]})
    return {"genomes": len(gens), "with_accession": len(with_acc), "without_accession": len(gens) - len(with_acc),
            "by_pattern": sum(1 for g in with_acc if g["how"] == "genome_id"),
            "by_name": sum(1 for g in with_acc if g["how"] == "ena_name"),
            "removed": removed, "by_source": dict(sorted(by_source.items())), "all": gens}


# ---- the rewrite -----------------------------------------------------------------------------------------------

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


def _tsv_text(head: list[str], rows) -> str:
    return "\t".join(head) + "\n" + "".join("\t".join(r) + "\n" for r in rows)


def _fasta_new_name(fasta: Path) -> str:
    n = fasta.name
    for suf in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if n.endswith(suf):
            return n[: -len(suf)] + "_r232" + suf
    return n + "_r232"


def _stats(head: list[str], rows: list[list[str]]) -> dict:
    gi, ci, li = _genome_col(head), head.index("contig_uid"), head.index("length")
    contigs: dict[str, int] = {}
    for r in rows:
        contigs.setdefault(r[ci], int(r[li]))
    return {"genomes": len({r[gi] for r in rows}), "contigs": len(contigs), "bases": sum(contigs.values())}


def _filter_side_table(path: Path, removed: set[str], kept_uids: set[str], ordinal_map: dict[str, str]):
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        head = next(rd)
        gcol = next((head.index(c) for c in GENOME_COLS if c in head), None)
        ucol = head.index("contig_uid") if "contig_uid" in head else None
        first_col = head.index("first_contig_uid") if "first_contig_uid" in head else None
        ocol = head.index("ordinal") if "ordinal" in head else None
        out = []
        for r in rd:
            if gcol is not None:
                if r[gcol] in removed:
                    continue
            elif ucol is not None and r[ucol] not in kept_uids:
                continue
            elif first_col is not None and r[first_col] not in kept_uids:
                continue
            if ocol is not None:
                if r[ocol] not in ordinal_map:
                    continue
                r = list(r)
                r[ocol] = ordinal_map[r[ocol]]
            out.append(r)
    return head, out


def _json_style(text: str) -> tuple[str, str]:
    return (",", ":") if '": ' not in text[:4000] and '", "' not in text[:4000] else (", ", ": ")


def _rewrite(stage: str, meta: Path, fasta: Path, removed: set[str], out_meta: Path, out_fasta_dir: Path,
             objects_path: str, objects_staging_path: str) -> dict:
    head, rows = _read_records(meta / "records.tsv")
    gi, ci = _genome_col(head), head.index("contig_uid")
    before = _stats(head, rows)
    order: list[str] = []
    uids_of: dict[str, list[str]] = {}
    for r in rows:
        if r[gi] not in uids_of:
            order.append(r[gi])
            uids_of[r[gi]] = []
        uids_of[r[gi]].append(r[ci])
    gtext = (meta / "genomes.json").read_text()
    genomes = json.loads(gtext)
    if len(genomes) != len(order):
        raise Refused(f"genomes.json has {len(genomes)} genomes, records.tsv {len(order)}")
    for g, name in zip(genomes, order):
        if g["contig_uids"] != uids_of[name]:
            raise Refused(f"genomes.json and records.tsv disagree on the contigs of genome {name}")
    unknown = removed - set(order)
    if unknown:
        raise Refused(f"removed genomes not in records.tsv: {sorted(unknown)[:3]}")
    kept_rows = [r for r in rows if r[gi] not in removed]
    kept_uids = {r[ci] for r in kept_rows}
    ordinal_map: dict[str, str] = {}
    oi = head.index("ordinal") if "ordinal" in head else None
    if oi is not None:
        for n, r in enumerate(kept_rows):
            ordinal_map[r[oi]] = str(n)
    else:
        ordinal_map = {}
    new_rows = []
    for n, r in enumerate(kept_rows):
        if oi is not None:
            r = list(r)
            r[oi] = str(n)
        new_rows.append(r)
    # FASTA first: it is the only long step, and it can refuse (a kept contig with no record)
    new_name = _fasta_new_name(fasta)
    out_fasta = Path(out_fasta_dir) / new_name
    if out_fasta.exists():
        raise Refused(f"{out_fasta} exists")
    out_fasta_dir.mkdir(parents=True, exist_ok=True)
    partial = out_fasta.with_name("." + new_name + ".partial")
    all_uids = {r[ci] for r in rows}
    written: set[str] = set()
    n_in = n_out = 0
    try:
        with open(partial, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0,
                                                       compresslevel=6) as gz, \
                (gzip.open(fasta, "rb") if str(fasta).endswith(".gz") else open(fasta, "rb")) as src:
            # A record is identified by the contig_uid of its sequence (the stage's own definition), not by its
            # header: SAG1's headers are the assemblers' names, NX1's are the uids.
            def flush(header: bytes | None, body: list[bytes]) -> None:
                nonlocal n_in, n_out
                if header is None:
                    return
                n_in += 1
                seq = b"".join(ln.strip() for ln in body)
                u = blake3(seq.upper()).digest(length=16).hex()
                if u not in all_uids:
                    raise Refused(f"FASTA record {header[1:60].decode(errors='replace').strip()} (contig_uid {u}) "
                                  f"is in no genome of records.tsv")
                if u in kept_uids and u not in written:
                    written.add(u)
                    n_out += 1
                    gz.write(header)
                    gz.writelines(body)

            header, body = None, []
            for line in src:
                if line.startswith(b">"):
                    flush(header, body)
                    header, body = line, []
                else:
                    body.append(line)
            flush(header, body)
        lost = kept_uids - written
        if lost:
            raise Refused(f"{len(lost)} kept contigs are missing from the FASTA; first {sorted(lost)[0]}")
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    partial.replace(out_fasta)
    # metadata
    out_meta.mkdir(parents=True, exist_ok=True)
    kept_genomes = [g for g, name in zip(genomes, order) if name not in removed]
    sep = _json_style(gtext)
    _atomic_text(out_meta / "genomes.json", json.dumps(kept_genomes, separators=sep))
    _atomic_text(out_meta / "records.tsv", _tsv_text(head, new_rows))
    for name in TSV_FILES:
        p = meta / name
        if p.is_file():
            h, r = _filter_side_table(p, removed, kept_uids, ordinal_map)
            _atomic_text(out_meta / name, _tsv_text(h, r))
    for name, newpath in (("objects.json", objects_path), ("objects_staging.json", objects_staging_path)):
        ot = (meta / name).read_text()
        obj = json.loads(ot)
        if len(obj) != 1:
            raise Refused(f"{name} lists {len(obj)} objects, expected 1")
        obj[0]["path"] = newpath
        _atomic_text(out_meta / name, json.dumps(obj, indent=2 if '\n  "' in ot else 1))
    after = _stats(head, new_rows)
    return {"before": before, "after": after, "fasta": str(out_fasta), "fasta_records_in": n_in,
            "fasta_records_out": n_out}


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


# ---- the command -------------------------------------------------------------------------------------------------

def run(stage: str, meta_dir, fasta, index_dir, out_meta_root, out_fasta_dir, *, names: dict[str, str] | None = None,
        inputs_home=Path("/tmp/if2453/s2/inputs"), check: bool = True, lookup=None, root: Path | None = None,
        chain_home: Path | None = None, later: dict[str, Path] | None = None, work: Path | None = None) -> dict:
    meta, fasta, out_dir = Path(meta_dir), Path(fasta), Path(out_meta_root) / stage
    targets = [out_dir / n for n in ("r232_rediff.json", "r232_rediff.tsv", "r232_removed.tsv")]
    if any(t.exists() for t in targets):
        raise Refused(f"{out_dir} exists")
    for n in ("records.tsv", "genomes.json", "objects.json", "objects_staging.json"):
        if not (meta / n).is_file():
            raise Refused(f"missing {meta / n}")
    if not fasta.is_file():
        raise Refused(f"missing {fasta}")
    d = rediff(meta / "records.tsv", lookup or index_lookup(index_dir), names)
    removed = d["removed"]
    changed = bool(removed)
    head, rows = _read_records(meta / "records.tsv")
    rep: dict = {"stage": stage, "changed": changed, "genomes": d["genomes"], "with_accession": d["with_accession"],
                 "without_accession": d["without_accession"], "by_pattern": d["by_pattern"],
                 "by_name": d["by_name"], "removed": len(removed), "by_source": d["by_source"],
                 "removed_by_reason": dict(sorted(Counter(r["reason"] for r in removed).items())),
                 "before": _stats(head, rows)}
    final_meta, final_fasta = meta, fasta
    if changed:
        if out_dir.exists() and any(out_dir.iterdir()):
            raise Refused(f"{out_dir} exists")
        new_name = _fasta_new_name(fasta)
        rw = _rewrite(stage, meta, fasta, {r["genome_id"] for r in removed}, out_dir, Path(out_fasta_dir),
                      str(Path(inputs_home) / stage / new_name), str(Path(out_fasta_dir) / new_name))
        rep.update(after=rw["after"], fasta_records_in=rw["fasta_records_in"],
                   fasta_records_out=rw["fasta_records_out"])
        final_meta, final_fasta = out_dir, Path(rw["fasta"])
    else:
        rep["after"] = rep["before"]
    out_dir.mkdir(parents=True, exist_ok=True)
    _atomic_text(out_dir / "r232_removed.tsv", _tsv_text(
        ["genome_id", "accession", "reason", "how", "source_name", "representative", "species", "globdb_r232_id"],
        [[str(r[k] or "") for k in ("genome_id", "accession", "reason", "how", "source_name", "representative",
                                    "species", "globdb_r232_id")] for r in removed]))
    _atomic_text(out_dir / "r232_rediff.tsv", _tsv_text(
        ["genome_id", "source_name", "accession", "how", "removed"],
        [[g["genome_id"], g["source_name"], g["accession"] or "", g["how"] or "",
          "1" if g["genome_id"] in {r["genome_id"] for r in removed} else "0"] for g in d["all"]]))
    rep["final_meta_dir"], rep["final_fasta"] = str(final_meta), str(final_fasta)
    if root is not None or chain_home is not None or later:
        rep["ownership"] = owners(final_meta / "records.tsv", root, chain_home, later or {},
                                  Path(work) if work else out_dir / "_work")
    if check:
        rep["s2_check"] = _check(stage, final_meta)
    files = {}
    for p in sorted(final_meta.iterdir()):
        if p.is_file() and not p.name.startswith("."):
            files[p.name] = {"bytes": p.stat().st_size, "sha256": _sha256(p)}
    files[final_fasta.name] = {"bytes": final_fasta.stat().st_size, "sha256": _sha256(final_fasta)}
    rep["files"] = files
    # the file holds no path, so that a second run to another directory gives the same bytes
    _atomic_text(out_dir / "r232_rediff.json", json.dumps(
        {k: v for k, v in rep.items() if k not in ("final_meta_dir", "final_fasta")}, indent=1, sort_keys=True) + "\n")
    return rep


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="R232 / r232 re-diff of a staged LR2, SAG1 or NX1 (card C-76, DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("stage")
    for n in ("meta", "fasta", "index", "out_meta", "out_fasta"):
        r.add_argument(n, type=Path)
    r.add_argument("--names", type=Path, default=None)
    r.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    r.add_argument("--root", type=Path, default=None)
    r.add_argument("--chain-home", type=Path, default=None)
    r.add_argument("--later", action="append", default=[], help="STAGE=DIR of an earlier r26.11.5 stage")
    r.add_argument("--no-check", action="store_true")
    r.add_argument("--work", type=Path, default=None)
    r.add_argument("--report-out", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        later = {}
        for item in a.later:
            k, _, v = item.partition("=")
            later[k] = Path(v)
        rep = run(a.stage, a.meta, a.fasta, a.index, a.out_meta, a.out_fasta,
                  names=read_names(a.names) if a.names else None, inputs_home=a.inputs_home,
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

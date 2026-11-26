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
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
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
    a = ap.parse_args(argv)
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

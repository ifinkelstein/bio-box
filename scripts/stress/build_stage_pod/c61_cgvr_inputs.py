"""Card C-61 (DEC-71, DEC-65, DEC-27): stage CGVR (120,568 gut viral genomes) from on-pod archive as CG1.

    python c61_cgvr_inputs.py stage SRC_FA OUT_DIR --retrieved ISO [--expected-sha256 HEX] \
        [--root ROOT] [--root-stages S1,P1,...] [--owner-inputs NAME=DIR ...] \
        [--final-dir DIR] [--compresslevel 1] [--out-report JSON]

Streams `cgvr_120568_viral_bins.fa.gz` in bounded memory (two streaming passes: pass 1 checks sha256 and counts
records per bin_id; pass 2 classifies and trims exact 20..1,000 bp direct terminal repeats on single-contig viral
records per DEC-65, writing deterministic gzip output `cgvr_120568_viral_bins.fna.gz`).

Because CGVR ships no topology metadata table:
- Single-contig viral records with an exact 20..1,000 bp terminal repeat (`classify(ends(seq)) == 'exact'`) have
  their trailing repeat copy trimmed once (recorded in `trims.tsv`) and are staged `circular`.
- All other single-contig records (`none`, `partial`, `long`, `low_complexity`, `unresolved`) and any multi-contig
  bin (`bin_counts[bin_id] > 1`) are staged untrimmed and `linear`.
- No record is dropped (DEC-27): genomes whose trimmed contig set repeats an earlier genome are listed once in
  `genomes.json` and logged in `repeated_records.tsv`.
- Ownership (`DEC-26`) is computed in build order across `--root-stages` (read-only `contig_index` default variants)
  followed by `--owner-inputs NAME=DIR` (`genomes.json`).
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # c48c_virus_inputs and make_inputs beside this script

import c48c_virus_inputs as v2  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

CGVR_FASTA = "cgvr_120568_viral_bins.fna.gz"
CGVR_SOURCE = {
    "source_name": "cgvr",
    "source_release": "OMIX007646-01 (2024-10-24; Wang et al. 2025, Microbiome 13:188, doi:10.1186/s40168-025-02185-9)",
    "license": "open access (NGDC CNCB OMIX007646, BioProject PRJCA031232; Wang et al. 2025 Microbiome 13:188 CC BY 4.0)",
    "redistributable": True,
}
TRIM_HEAD = [
    "source", "genome", "accession", "source_contig_uid", "contig_uid",
    "source_length", "repeat_length", "bases_removed", "length",
]
RECORD_HEAD = [
    "ordinal", "genome", "accession", "contig_uid", "length", "dtr_class",
    "circular", "source_contig_uid", "source_length", "trimmed_bases",
]


def parse_bin_id(header: str) -> str:
    """Extract the viral bin identifier from a FASTA header (first whitespace token, before any `|`)."""
    tok = header.split(None, 1)[0] if header else ""
    return tok.split("|", 1)[0] if "|" in tok else tok


def _open_fasta_in(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return open(path, "rb")


def _count_bins(src_fa: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with _open_fasta_in(src_fa) as fh:
        for raw in fh:
            if raw[:1] == b">":
                hd = raw[1:].rstrip(b"\r\n").decode("utf-8", "replace")
                bid = parse_bin_id(hd)
                if not bid:
                    raise Refused(f"{src_fa}: empty FASTA header")
                counts[bid] += 1
    if not counts:
        raise Refused(f"{src_fa}: no FASTA records found")
    return counts


def stage_cg1(
    src_fa: Path,
    out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    expected_sha256: str | None,
    retrieved: str,
    final_dir: Path | None = None,
    compresslevel: int = 1,
) -> dict:
    src_fa = Path(src_fa)
    out = Path(out)
    if out.exists():
        raise Refused(f"{out} exists")
    if not src_fa.is_file():
        raise Refused(f"{src_fa} is not a file")
    got_sha = ids.sha256_file(src_fa).hex()
    if expected_sha256 is not None and got_sha != expected_sha256:
        raise Refused(f"sha256 of {src_fa} is {got_sha}, expected {expected_sha256}")

    work = out.with_name(f"{out.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out
    try:
        rep = _stage_cg1(
            src_fa, obj_dir, work, root, list(root_stages), owner_inputs,
            got_sha, retrieved, compresslevel,
        )
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    work.rename(out)
    names = [
        "genomes.json", "objects.json", "records.tsv", "repeated_records.tsv",
        "shared_records.tsv", "trims.tsv", CGVR_FASTA,
    ]
    rep["sha256"] = {n: ids.sha256_file(out / n).hex() for n in names}
    rep["bytes"] = {n: (out / n).stat().st_size for n in names}
    return rep


def _stage_cg1(
    src_fa: Path,
    obj_dir: Path,
    work: Path,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    src_sha256: str,
    retrieved: str,
    compresslevel: int,
) -> dict:
    bin_counts = _count_bins(src_fa)

    # bin_id -> [(cu, staged_len, src_len, rec_ord, accession, dtr_cls, src_cu, trimmed_k, circ)]
    parts: dict[str, list[tuple[bytes, int, int, int, str, str, bytes, int, bool]]] = {}
    dtr_classes: Counter[str] = Counter()
    trims: list[tuple] = []
    n_recs = 0
    src_bases = 0
    staged_bases = 0

    with _open_fasta_in(src_fa) as fh, open(work / CGVR_FASTA, "wb") as raw_w:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_w, compresslevel=compresslevel, mtime=0) as w:
            for r in v2.iter_fasta_stream(fh):
                bid = parse_bin_id(r.header)
                src_len = len(r.seq)
                src_cu = ids.contig_uid(r.seq)
                src_bases += src_len
                if bin_counts[bid] == 1:
                    e = v2.ends(r.seq)
                    dtr_cls = v2.classify(e)
                    dtr_classes[dtr_cls] += 1
                    if dtr_cls == "exact":
                        trimmed_k = e.exact
                        seq = v2.trim(r.seq, trimmed_k)
                        circ = True
                    else:
                        trimmed_k = 0
                        seq = r.seq
                        circ = False
                else:
                    dtr_cls = "multi_contig"
                    dtr_classes[dtr_cls] += 1
                    trimmed_k = 0
                    seq = r.seq
                    circ = False

                w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
                staged_len = len(seq)
                staged_bases += staged_len
                cu = ids.contig_uid(seq) if trimmed_k else src_cu
                if trimmed_k:
                    trims.append((
                        "cgvr", bid, r.accession, src_cu.hex(), cu.hex(),
                        src_len, trimmed_k, trimmed_k, staged_len,
                    ))
                parts.setdefault(bid, []).append((
                    cu, staged_len, src_len, n_recs, r.accession, dtr_cls, src_cu, trimmed_k, circ,
                ))
                n_recs += 1

    trims.sort(key=lambda t: (t[0], t[1], t[2]))
    all_uids = [cu for p in parts.values() for cu, *_ in p]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    repeated: list[tuple] = []
    shared: list[tuple] = []
    rows: list[tuple] = []
    first_of: dict[bytes, str] = {}
    ordinal = 0

    order = sorted(parts, key=lambda b: parts[b][0][3])  # FASTA order of first record per bin
    for bid in order:
        recs = parts[bid]
        uids = [x[0] for x in recs]
        circ = len(recs) == 1 and recs[0][8]
        for cu, n, src_n, _, acc, dtr_cls, src_cu, trimmed_k, cflag in recs:
            rows.append((
                ordinal, bid, acc, cu.hex(), n, dtr_cls, str(cflag).lower(),
                src_cu.hex(), src_n, trimmed_k,
            ))
            if cu in owners:
                shared.append((ordinal, bid, cu.hex(), owners[cu]))
            ordinal += 1
        g = listing.add(uids, 11, "virus", "vOTU")
        guid = ids.genome_uid(list(dict.fromkeys(uids)))
        if g is None:
            repeated.append((ordinal - len(recs), bid, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [uids[0].hex()] if circ else []
        genomes.append(g)
        first_of[guid] = bid

    make_inputs._write_tsv(work / "records.tsv", RECORD_HEAD, rows)
    make_inputs._write_tsv(work / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(work / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"], shared)
    make_inputs._write_tsv(work / "trims.tsv", TRIM_HEAD, trims)

    objects = [{"path": str(obj_dir / CGVR_FASTA), **CGVR_SOURCE, "retrieved_at": retrieved}]
    (work / "objects.json").write_text(json.dumps(objects, indent=1))
    (work / "genomes.json").write_text(json.dumps(genomes))

    contig_len: dict[bytes, int] = {}
    for r in rows:
        contig_len.setdefault(bytes.fromhex(r[3]), r[4])

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    for cu, n in contig_len.items():
        if cu in owners:
            st = owners[cu]
            owner_contigs[st] += 1
            owner_bases[st] += n

    by_rec = dict(sorted(Counter(r[3] for r in shared).items()))
    held_contigs = sum(owner_contigs.values())
    held_bases = sum(owner_bases.values())

    return {
        "source_fasta": str(src_fa),
        "source_sha256": src_sha256,
        "source_bytes": src_fa.stat().st_size,
        "records": n_recs,
        "bins": len(parts),
        "multi_contig_bins": sum(1 for p in parts.values() if len(p) > 1),
        "source_bases": src_bases,
        "staged_bases": staged_bases,
        "dtr_classes": dict(sorted(dtr_classes.items())),
        "trims": len(trims),
        "trimmed_bases": sum(t[7] for t in trims),
        "genomes": len(genomes),
        "contigs": len(contig_len),
        "contig_bases": sum(contig_len.values()),
        "circular_genomes": sum(1 for g in genomes if g["circular"]),
        "repeated_genomes": len(repeated),
        "shared_records": len(shared),
        "owner_records": by_rec,
        "owner_contigs": dict(sorted(owner_contigs.items())),
        "owner_bases": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": len(contig_len) - held_contigs,
        "new_bases": sum(contig_len.values()) - held_bases,
        "root_stages": root_stages,
        "owner_inputs": list(owner_inputs),
        "build_order": [*root_stages, *owner_inputs, "CG1"],
        **listing.report(),
        "objects": len(objects),
    }


def _parse_owner_inputs(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        name, path = item.split("=", 1)
        if not name or not path:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        out[name] = Path(path)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage CGVR viral genomes as CG1 (card C-61).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sp = sub.add_parser("stage")
    sp.add_argument("src_fa", type=Path)
    sp.add_argument("out_dir", type=Path)
    sp.add_argument("--retrieved", required=True)
    sp.add_argument("--expected-sha256", default=None)
    sp.add_argument("--root", type=Path, default=None)
    sp.add_argument("--root-stages", default="")
    sp.add_argument("--owner-inputs", action="append", default=[])
    sp.add_argument("--final-dir", type=Path, default=None)
    sp.add_argument("--compresslevel", type=int, default=1)
    sp.add_argument("--out-report", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        root_stages = [s for s in a.root_stages.split(",") if s]
        owners = _parse_owner_inputs(a.owner_inputs)
        rep = stage_cg1(
            a.src_fa,
            a.out_dir,
            root=a.root,
            root_stages=root_stages,
            owner_inputs=owners,
            expected_sha256=a.expected_sha256,
            retrieved=a.retrieved,
            final_dir=a.final_dir,
            compresslevel=a.compresslevel,
        )
    except Refused as e:
        print(f"c61_cgvr_inputs.py: refused: {e}", file=sys.stderr)
        return 2
    txt = json.dumps(rep, indent=1, sort_keys=True) + "\n"
    if a.out_report is not None:
        a.out_report.parent.mkdir(parents=True, exist_ok=True)
        a.out_report.write_text(txt)
    print(json.dumps(rep, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

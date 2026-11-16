"""Card C-39: the `build-stage` inputs of stage V1 (UHGV v1.0 vOTU representatives), in the form of card C-26 part 1.

    python c39_inputs.py uhgv FNA_ZST VOTU_META_ZST OUT_DIR RETRIEVED_AT [--root ROOT]

FNA_ZST is UHGV's votus_full.fna.zst, decompressed byte for byte (Python's compression.zstd) to
OUT_DIR/uhgv_v1.0_votus_full.fna, the one source object (written under a temporary name, then moved). VOTU_META_ZST
is votus_metadata_extended.tsv.zst, one row per vOTU keyed by `uhgv_genome` (the representative genome). Every
record is joined to its row by the header's first token; a record without a row, a row without a record, or a
value outside the maps below is refused (exit 2) before genomes.json is written.

One genome per record, kind virus (gene calling in gv-meta mode, `genecall.mode_of`):
- genetic_code: the row's `genetic_code` (prodigal-gv's pick); refused unless it is an NCBI table number.
- assembly_quality: `complete` when `checkv_quality` is Complete, else `vOTU` (the schema's value for a vOTU
  representative whose completeness is short of complete or not determined); the CheckV tiers are counted.
- circular: when `checkv_completeness_method` is DTR (a direct terminal repeat: the molecule is circular or
  terminally redundant) and the sequence no longer carries the repeat (its first 20 bases are not found again in its
  last 2,000). The gene caller wraps a circular contig end to start (`genecall._genes`), which is right for a DTR
  genome whose second copy of the repeat was removed and would read the repeat twice where it was not; a DTR record
  that still carries a terminal repeat is listed linear, as is every other record (ITR, provirus, fragments).
No record is dropped (DEC-27): a record whose bytes repeat an earlier record's is that record's genome, listed once,
written to OUT_DIR/repeated_records.tsv. With --root, every stage on disk is read (contig_index, default variants,
`make_inputs.stage_owners`): a record whose contig a stage owns stays a genome of its own (that stage owns its calls,
DEC-26) and is written to OUT_DIR/shared_records.tsv with the owner. Also OUT_DIR/records.tsv: per record its
ordinal, genome id, vOTU, contig_uid, length, CheckV quality and method, code (the overlap and estimate scripts read
it). genomes.json, objects.json and the tables are written under OUT_DIR.tmp and renamed when complete; an existing
OUT_DIR is refused. Prints a JSON report.
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
from collections import Counter
from compression import zstd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # make_inputs beside this script

import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

ZOPT = {zstd.DecompressionParameter.window_log_max: 31}  # UHGV's files use a long window (zstd --long)

NCBI_CODES = {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33}
CHECKV = {"Complete": "complete", "High-quality": "vOTU", "Medium-quality": "vOTU", "Low-quality": "vOTU",
          "Not-determined": "vOTU"}
SOURCE = {"source_name": "uhgv", "source_release": "zenodo:10.5281/zenodo.17402089 v1.0",
          "license": "CC BY 4.0", "redistributable": True}
FASTA_NAME = "uhgv_v1.0_votus_full.fna"


def zst_rows(path: Path) -> list[dict]:
    with zstd.open(path, "rt", newline="", options=ZOPT) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def decompress(src: Path, dst: Path) -> None:
    tmp = dst.with_name(f".{dst.name}.tmp")
    with zstd.open(src, "rb", options=ZOPT) as r, open(tmp, "wb") as w:
        shutil.copyfileobj(r, w, 1 << 22)
    tmp.replace(dst)


def has_terminal_repeat(seq: bytes, k: int = 20, window: int = 2000) -> bool:
    s = seq.upper()
    return len(s) >= 2 * k and s[:k] in s[max(k, len(s) - window):]


def uhgv(fna_zst: Path, meta_zst: Path, out: Path, retrieved: str, root: Path | None) -> dict:
    if out.exists():
        raise make_inputs.Refused(f"{out} exists")
    meta = {}
    for r in zst_rows(meta_zst):
        if r["uhgv_genome"] in meta:
            raise make_inputs.Refused(f"{meta_zst}: genome {r['uhgv_genome']} has two rows")
        meta[r["uhgv_genome"]] = r
    work = out.with_name(f"{out.name}.tmp")
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    fa = out.parent / f".{out.name}.{FASTA_NAME}"
    final_fa = out / FASTA_NAME
    decompress(fna_zst, fa)
    owners = make_inputs.stage_owners(root) if root is not None else None
    listing = make_inputs.Listing()
    genomes, repeated, shared, recs = [], [], [], []
    first: dict[bytes, tuple] = {}
    seen: set[str] = set()
    quality, method, codes = Counter(), Counter(), Counter()
    dtr = dtr_repeat = record_bases = 0
    for k, r in enumerate(fasta.iter_records(fa)):
        gid = r.accession
        m = meta.get(gid)
        if m is None:
            raise make_inputs.Refused(f"record {k} {gid!r}: no row in {meta_zst}")
        if gid in seen:
            raise make_inputs.Refused(f"record {k}: genome {gid} is a second record")
        seen.add(gid)
        if int(m["genome_length"]) != len(r.seq):
            raise make_inputs.Refused(f"{gid}: genome_length {m['genome_length']} but {len(r.seq)} bases")
        try:
            code = int(m["genetic_code"])
        except ValueError:
            raise make_inputs.Refused(f"{gid}: genetic_code {m['genetic_code']!r}") from None
        if code not in NCBI_CODES:
            raise make_inputs.Refused(f"{gid}: genetic_code {code} is not an NCBI table")
        if m["checkv_quality"] not in CHECKV:
            raise make_inputs.Refused(f"{gid}: checkv_quality {m['checkv_quality']!r} not in {sorted(CHECKV)}")
        q = CHECKV[m["checkv_quality"]]
        quality[m["checkv_quality"]] += 1
        method[m["checkv_completeness_method"]] += 1
        codes[code] += 1
        circ = False
        if m["checkv_completeness_method"] == "DTR":
            dtr += 1
            kept_repeat = has_terminal_repeat(r.seq)
            dtr_repeat += kept_repeat
            circ = not kept_repeat
        record_bases += len(r.seq)
        cu = ids.contig_uid(r.seq)
        recs.append((k, gid, m["uhgv_votu"], cu.hex(), len(r.seq), m["checkv_quality"],
                     m["checkv_completeness_method"], code, str(circ).lower()))
        if owners is not None and cu in owners:
            shared.append((k, gid, cu.hex(), owners[cu]))
        if cu in first:
            f = first[cu]
            repeated.append((k, gid, cu.hex(), f[0], f[1], str(code == f[2]).lower(), str(q == f[3]).lower(),
                             str(circ == f[5]).lower()))
            listing.add([cu], code, "virus", q)
            continue
        g = listing.add([cu], code, "virus", q)
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first[cu] = (k, gid, code, q, len(r.seq), circ)
    missing = sorted(set(meta) - seen)
    if missing:
        raise make_inputs.Refused(f"{len(missing)} rows of {meta_zst} have no record; first {missing[0]}")
    make_inputs._write_tsv(work / "repeated_records.tsv", ["ordinal", "genome", "contig_uid", "first_ordinal",
                                                           "first_genome", "same_code", "same_quality",
                                                           "same_topology"], repeated)
    make_inputs._write_tsv(work / "records.tsv", ["ordinal", "genome", "votu", "contig_uid", "length",
                                                  "checkv_quality", "checkv_method", "genetic_code", "circular"], recs)
    if owners is not None:
        make_inputs._write_tsv(work / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"],
                               shared)
    objects = [{"path": str(final_fa), **SOURCE, "retrieved_at": retrieved}]
    (work / "objects.json").write_text(json.dumps(objects, indent=1))
    (work / "genomes.json").write_text(json.dumps(genomes))
    fa.replace(work / FASTA_NAME)
    work.rename(out)
    rep = {"records": len(seen), "genomes": len(genomes), "contigs": len(first), "record_bases": record_bases,
           "genome_bases": sum(f[4] for f in first.values()), "checkv_quality": dict(sorted(quality.items())),
           "assembly_quality": dict(Counter(g["assembly_quality"] for g in genomes)),
           "checkv_completeness_method": dict(sorted(method.items())),
           "codes": {str(c): n for c, n in sorted(codes.items())},
           "genome_codes": {str(c): n for c, n in sorted(Counter(g["genetic_code"] for g in genomes).items())},
           "dtr_records": dtr, "dtr_records_with_terminal_repeat": dtr_repeat,
           "circular_genomes": sum(1 for g in genomes if g["circular"]),
           "circular_genome_bases": sum(f[4] for f in first.values() if f[5]),
           "repeated_records": len(repeated), **listing.report(), "objects": len(objects),
           "fasta": str(final_fa), "fasta_bytes": final_fa.stat().st_size,
           "fasta_sha256": ids.sha256_file(final_fa).hex()}
    if owners is not None:
        by_stage = lambda rows: {s: sum(1 for r in rows if r[3] == s) for s in sorted({r[3] for r in rows})}  # noqa: E731
        rep.update({"shared_records": len(shared), "owner_records": by_stage(shared),
                    "owner_contigs": by_stage(list({r[2]: r for r in shared}.values())),
                    "new_contigs": sum(1 for cu in first if cu not in owners),
                    "new_bases": sum(f[4] for cu, f in first.items() if cu not in owners)})
    names = ["genomes.json", "objects.json", "repeated_records.tsv", "records.tsv"] + \
        (["shared_records.tsv"] if owners is not None else [])
    rep["sha256"] = {n: ids.sha256_file(out / n).hex() for n in names}
    return rep


def main(argv: list[str]) -> int:
    try:
        argv, root = make_inputs._option(argv, "--root")
    except make_inputs.Refused:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    if len(argv) != 5 or argv[0] != "uhgv":
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    try:
        rep = uhgv(Path(argv[1]), Path(argv[2]), Path(argv[3]), argv[4], None if root is None else Path(root))
    except make_inputs.Refused as e:
        print(f"c39_inputs.py: refused: {e}", file=sys.stderr)
        return 2
    print(json.dumps(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

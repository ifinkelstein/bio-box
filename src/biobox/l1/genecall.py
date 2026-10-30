"""Gene calling (operation O-CALL) and the protein tables it feeds.

Policy p2 (approved 2026-10-29, PLAN.md §3.1; field practice in
planning/DECISIONS.md DEC-13):

- viral contigs: pyrodigal-gv in meta mode (prodigal-gv models for giant
  viruses and alternative codes, including code 15), code taken per gene;
- standalone plasmids: meta mode;
- prokaryotic genomes >= 100 kb: single mode trained on all contigs of the
  genome together with the genome's code; below 100 kb: meta mode;
- circular contigs: genes may wrap across the origin. Genes are called on the
  contig with `WRAP_FLANK` bp of the other end added on both sides and kept
  when their start lies on the contig, so each gene appears once with full
  context; a wrapping gene has `end > contig_len`.

Writes, per gene-calling version and stage (paths under `<root>/`):

    l0seq/protein_seq/gene_call_version=<v>/stage=<S>/protein_prefix=<xx>/part-0.parquet
    l1/slot_protein/gene_call_version=<v>/stage=<S>/part-0.parquet
    l1/protein_occurrence/gene_call_version=<v>/stage=<S>/part-0.parquet

Coordinates are 0-based half-open (pyrodigal reports 1-based inclusive).
`dist_to_contig_end_bp` is 0 for genes running off a linear contig end and
equals the contig length on circular contigs, which have no edge.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pyrodigal
import pyrodigal_gv

from biobox import ids, ledger
from biobox.l0seq.store import SeqStore
from biobox.registry.schemas import table

SINGLE_MODE_MIN_BP = 100_000  # Prodigal: below this, use anonymous/meta mode
WRAP_FLANK = 50_000
GENE_CALL_VERSION = f"pyrodigal-{pyrodigal.__version__}-gv-{pyrodigal_gv.__version__}-p2"
KINDS = ("prokaryote", "virus", "plasmid")


@dataclass(frozen=True)
class GenomeInput:
    contig_uids: tuple[bytes, ...]
    genetic_code: int = 11  # used in single mode; meta modes report the code per gene
    kind: str = "prokaryote"
    circular: tuple[bytes, ...] = ()  # contig_uids that are circular
    assembly_quality: str = "complete"
    tier: int = 1
    is_species_rep: bool = False


def _write(tbl: pa.Table, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(tbl, path, compression="zstd", row_group_size=122_880)


def _finder(genome: GenomeInput, seqs: dict[bytes, str]):
    if genome.kind not in KINDS:
        raise ValueError(f"unknown genome kind {genome.kind!r}")
    if genome.kind == "virus":
        return pyrodigal_gv.ViralGeneFinder(meta=True), "gv-meta"
    if genome.kind == "plasmid" or sum(len(s) for s in seqs.values()) < SINGLE_MODE_MIN_BP:
        return pyrodigal.GeneFinder(meta=True), "meta"
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*seqs.values(), translation_table=genome.genetic_code)
    return finder, "single"


def _genes(finder, seq: str, circular: bool):
    """Yield (start, end, strand, aa, code, partial) with 0-based half-open coordinates."""
    n = len(seq)
    if not circular:
        for g in finder.find_genes(seq):
            yield g.begin - 1, g.end, g.strand, g.translate(include_stop=False), g.translation_table, bool(
                g.partial_begin or g.partial_end
            )
        return
    flank = min(n, WRAP_FLANK)
    ext = seq[n - flank :] + seq + seq[:flank]
    for g in finder.find_genes(ext):
        start = g.begin - 1 - flank
        if not 0 <= start < n or g.partial_begin or g.partial_end:
            continue
        end = g.end - flank
        yield start, end, g.strand, g.translate(include_stop=False), g.translation_table, False


def call_genome(store: SeqStore, genome: GenomeInput) -> list[dict]:
    seqs = {cu: store.get(cu).upper().decode("ascii") for cu in genome.contig_uids}
    finder, mode = _finder(genome, seqs)
    guid = ids.genome_uid(list(genome.contig_uids))
    calls: list[dict] = []
    for cu in sorted(seqs):
        seq = seqs[cu]
        circular = cu in genome.circular
        genes = sorted(_genes(finder, seq, circular), key=lambda t: (t[0], t[1], t[2]))
        for ordinal, (start, end, strand, aa, code, partial) in enumerate(genes):
            if circular:
                dist = len(seq)
            else:
                dist = 0 if partial else min(start, len(seq) - end)
            calls.append({
                "contig_uid": cu,
                "slot_ordinal": ordinal,
                "genome_uid": guid,
                "start": start,
                "end": end,
                "strand": strand,
                "aa": aa,
                "genetic_code": code,
                "mode": mode,
                "wraps_origin": end > len(seq),
                "dist_to_contig_end_bp": dist,
                "assembly_quality": genome.assembly_quality,
                "tier": genome.tier,
                "is_species_rep": genome.is_species_rep,
            })
    return calls


def run(root: Path, stage: str, genomes: list[GenomeInput], version: str = GENE_CALL_VERSION) -> dict:
    root = Path(root)
    store = SeqStore(root)
    calls = [c for g in genomes for c in call_genome(store, g)]

    slot_spec = table("slot_protein")
    slot_rows = sorted(
        (
            {
                "contig_uid": c["contig_uid"],
                "slot_ordinal": c["slot_ordinal"],
                "protein_uid": ids.protein_uid(c["aa"]),
                "genetic_code": c["genetic_code"],
                "edit_type": "none",
            }
            for c in calls
        ),
        key=lambda r: (r["contig_uid"], r["slot_ordinal"]),
    )
    _write(
        pa.Table.from_pylist(slot_rows, schema=slot_spec.schema),
        root / f"l1/slot_protein/gene_call_version={version}/stage={stage}/part-0.parquet",
    )

    occ_spec = table("protein_occurrence")
    occ_rows = sorted(
        (
            {
                "protein_uid": ids.protein_uid(c["aa"]),
                "contig_uid": c["contig_uid"],
                "slot_ordinal": c["slot_ordinal"],
                "genome_uid": c["genome_uid"],
                "start": c["start"],
                "end": c["end"],
                "strand": c["strand"],
                "tax_order": None,
                "is_species_rep": c["is_species_rep"],
                "tier": c["tier"],
                "assembly_quality": c["assembly_quality"],
                "dist_to_contig_end_bp": c["dist_to_contig_end_bp"],
            }
            for c in calls
        ),
        key=lambda r: (r["protein_uid"], r["contig_uid"], r["slot_ordinal"]),
    )
    _write(
        pa.Table.from_pylist(occ_rows, schema=occ_spec.schema),
        root / f"l1/protein_occurrence/gene_call_version={version}/stage={stage}/part-0.parquet",
    )

    prot_spec = table("protein_seq")
    distinct = {ids.protein_uid(c["aa"]): c["aa"] for c in calls}
    by_prefix: dict[str, list[bytes]] = defaultdict(list)
    for pu in distinct:
        by_prefix[pu[:1].hex()].append(pu)
    for prefix in sorted(by_prefix):
        rows = [
            {
                "protein_uid": pu,
                "protein_digest": ids.raw_digest(distinct[pu].encode("ascii")),
                "length": len(distinct[pu]),
                "seq": distinct[pu],
            }
            for pu in sorted(by_prefix[prefix])
        ]
        _write(
            pa.Table.from_pylist(rows, schema=prot_spec.schema),
            root / f"l0seq/protein_seq/gene_call_version={version}/stage={stage}/protein_prefix={prefix}/part-0.parquet",
        )
    # ledger: per contig, genes in; kept = first gene of its protein sequence in (contig, slot) order
    seen: set[bytes] = set()
    per: dict[bytes, list[int]] = {cu: [0, 0] for g in genomes for cu in g.contig_uids}
    for r in slot_rows:
        c = per[r["contig_uid"]]
        c[0] += 1
        if r["protein_uid"] not in seen:
            seen.add(r["protein_uid"])
            c[1] += 1
    ledger.write(root, stage, "genecall",
                 [ledger.row(cu, "contig", n, k, n - k) for cu, (n, k) in per.items()])
    return {
        "genes": len(calls),
        "distinct_proteins": len(distinct),
        "modes": sorted({c["mode"] for c in calls}),
        "wrapping_genes": sum(c["wraps_origin"] for c in calls),
    }

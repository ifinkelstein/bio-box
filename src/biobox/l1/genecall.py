"""Gene calling with pyrodigal (operation O-CALL) and the protein tables it feeds.

Writes, per gene-calling version and stage (paths under `<root>/`):

    l0seq/protein_seq/gene_call_version=<v>/stage=<S>/protein_prefix=<xx>/part-0.parquet
    l1/slot_protein/gene_call_version=<v>/stage=<S>/part-0.parquet
    l1/protein_occurrence/gene_call_version=<v>/stage=<S>/part-0.parquet

Coordinates are 0-based half-open on the stored contig (pyrodigal reports
1-based inclusive). Slots are numbered per contig in coordinate order. Genes
running off a contig end keep `dist_to_contig_end_bp = 0`.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pyrodigal

from biobox import ids
from biobox.l0seq.store import SeqStore
from biobox.registry.schemas import table

# Policy p1: single mode trained per genome when the genome is >= 100 kb
# (Prodigal's training minimum); meta mode below that. Viral contigs move to
# pyrodigal-gv once geNomad regimes exist (PLAN.md §3.1).
SINGLE_MODE_MIN_BP = 100_000
GENE_CALL_VERSION = f"pyrodigal-{pyrodigal.__version__}-p1"


@dataclass(frozen=True)
class GenomeInput:
    contig_uids: tuple[bytes, ...]
    genetic_code: int
    assembly_quality: str = "complete"
    tier: int = 1
    is_species_rep: bool = False


def _write(tbl: pa.Table, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(tbl, path, compression="zstd", row_group_size=122_880)


def call_genome(store: SeqStore, genome: GenomeInput) -> list[dict]:
    """Train on all contigs of one genome, then call genes per contig."""
    seqs = {cu: store.get(cu).upper().decode("ascii") for cu in genome.contig_uids}
    if sum(len(s) for s in seqs.values()) >= SINGLE_MODE_MIN_BP:
        finder = pyrodigal.GeneFinder(meta=False)
        finder.train(*seqs.values(), translation_table=genome.genetic_code)
    else:
        finder = pyrodigal.GeneFinder(meta=True)
    guid = ids.genome_uid(list(genome.contig_uids))
    calls: list[dict] = []
    for cu in sorted(seqs):
        seq = seqs[cu]
        genes = sorted(finder.find_genes(seq), key=lambda g: (g.begin, g.end, g.strand))
        for ordinal, g in enumerate(genes):
            aa = g.translate(include_stop=False)
            start, end = g.begin - 1, g.end
            calls.append({
                "contig_uid": cu,
                "slot_ordinal": ordinal,
                "genome_uid": guid,
                "start": start,
                "end": end,
                "strand": g.strand,
                "aa": aa,
                "genetic_code": g.translation_table,
                "partial": bool(g.partial_begin or g.partial_end),
                "dist_to_contig_end_bp": 0 if (g.partial_begin or g.partial_end) else min(start, len(seq) - end),
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

    # protein_seq: one row per distinct protein, partitioned by uid prefix
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
    return {"genes": len(calls), "distinct_proteins": len(distinct)}

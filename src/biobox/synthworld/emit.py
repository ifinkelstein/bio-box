"""synthworld emit: per-source FASTA, the ingest manifest and the build-stage inputs, the truth tables in
emitted-contig coordinates, a convenience GFF and the world manifest (planning/tests/synthworld_spec.md §2.1 S11,
§3, §6.2; card F-02 session 8). Reads a world whose S1-S10 are built and writes ``world/emit/`` and
``world/manifest.json``. Nothing of an earlier stage is rewritten.

**Sources.** Every source of ``sources/sources.tsv`` is one FASTA file, ``fasta/<accession>.fna``, its contigs
in ``source_contigs.tsv`` order under their contig accessions (``>SWC00000001.1``), uncompressed at 80 columns
(§6.2), the bytes those of the S9 contig, the S10 contaminant or the S10 duplicate (their ``contig_uid`` is
checked). The 12 phages of S6 (``mobile/phages.fna``) are emitted too, each a source of its own (``src_a``, kind
virus), with accessions and contig accessions numbered after the host sources' and one BioProject for the set.

**Ingest manifest** (``ingest.tsv``, one row per contig occurrence): ``source_class`` isolate (complete or draft
assembly), MAG or ``MGE database`` (a phage); ``replicon_class`` chromosome, plasmid or virus (PLAN.md §2's
vocabulary; a contaminant takes its donor contig's); ``molecule_type`` DNA; ``topology`` circular only for a
complete circular contig; the genome's genetic code; ``kind`` prokaryote or virus and ``assembly_quality``
complete, contig (draft), mag_hq (MAG completeness >= 0.9) or mag_mq, as the builders' ``genome_input`` takes them
(phages complete); ``tax_lineage`` the synthetic lineage, mislabelled for a source with probability
``emit.mislabel_rate`` (stream ``S11/mislabel/<source>/<genome>``: the lineage of a genome of another genus, drawn
uniformly); ``biome`` one of ``emit.biomes`` per taxon at ``emit.biome_rank`` (stream ``S11/biome/<taxon>``),
empty for phages. ``genomes.json`` and ``objects.json`` are ``bio-box build-stage``'s ``--genomes`` and
``--objects`` for all sources in one stage: one genome per source, its distinct contigs in file order, a genome
whose contig set an earlier source already listed (an exact duplicate) listed once (DEC-27, as
``scripts/stress/build_stage_pod/make_inputs.py``'s ``Listing``); object paths relative to ``emit/``.

**Truth** (``truth/<table>.tsv``, §3; every row starts with ``world_id`` and ``release_tag``). Coordinates are
0-based half-open on the emitted contig, strand +1/-1, a feature across a circular contig's origin has end >
length. A feature of an earlier stage is lifted to the S8 frame (S5 -> S6 by ``mobile/insertions.tsv``, S6 -> S7
by ``fusions/edits.tsv``, S7 -> S8 by ``assembly/edits.tsv``), placed on the S9 contigs (``assembly/contigs.tsv``)
and from there on every emitted occurrence of the contig: the contig itself, a near duplicate's copy (same
coordinates, another ``contig_uid``) and a contaminant window cut from it (clipped and flipped as the window
was). An exact duplicate shares the ``contig_uid``, so its rows are the original's. One row per placed piece;
``status`` says whole (one piece, nothing clipped), cut (a piece of a feature split by a contig edge or clipped
by a contaminant window; ``partial_5``/``partial_3`` on genes), or, with no coordinates, why the feature is not on
an emitted contig: dropped (on contigs not emitted), s6_insertion, s7_event or s8_edit (a later stage's edit
touched it), partial (a call error's gene not whole on any emitted contig). ``occurrence`` is assembly, near_dup,
contaminant or phage. Stable keys: ``contig_uid`` (BLAKE3-128 of the contig), ``feature_uid`` =
``ids.feature_uid(contig_uid, start, end, strand, kind)`` (cds for genes, call errors and fused ORFs; ncrna for
RNA motifs; crispr_array for arrays), ``region_uid`` = ``ids.region_uid(contig_uid, start, end, +1)`` for IGRs,
``genome_uid`` over a source's distinct contigs; no ordinal is a key. ``gene_id`` is the S8 gene's copy label, so
a gene's pieces and copies share it (the row key is ``feature_uid``). ``stop_bp`` is end on + and start on -
(empty when the 3' end is cut off). An IGR of zero or negative length has start = end = the upper end of the
overlap on the contig (``overlap_bp``), and ``left_gene_id`` is the gene at the lower contig coordinate.
``igr_uid`` is BLAKE3-128 of the IGR's bases on the contig's + strand. List columns are ``;``-separated
(``domain_arch`` ``acc:start-end``, ``elements`` ``kind:start-end:strand``, ``expected_calls``
``start-end:strand``, ``gene_ids``, ``partner_family_ids``).

Not computed (empty): ``ref_bitscore`` and ``ref_bitscore_median`` (no scan at build), ``gene_tree_extant``,
``origin_node``, ``strength_class``, ``realized_independent_transfers``, ``realized_vehicle_independent_loci``,
``powered_expected``; ``truth_igr`` holds the S5 junctions only (junctions made by S6-S8 insertions are not
listed); domain architectures are given for S4 and S7 genes and phage genes, not for Cas genes or smORFs.

Shape on disk (``build(world)``)::

    emit/fasta/<accession>.fna, .fai   one file per source
    emit/ingest.tsv                    the ingest manifest (above)
    emit/genomes.json, objects.json    build-stage inputs
    emit/truth/truth_<name>.tsv        genome contig edges gene igr rna spacer_target fusion callerr dup contam
                                       family module traps
    emit/truth.gff                     genes and RNA motifs, GFF3, one seqid per contig accession
    emit/emit.json                     format, seed, the emit preset and its digest, world_id, release_tag,
                                       upstream digests, versions, stream paths, counts, sha256 of every file,
                                       digest
    manifest.json                      format, world_id, release_tag, preset name and digest, seed, versions,
                                       panel and motif sources, every stage's digest and format, world_digest

``world_id`` is ``<preset>-<seed>``; ``release_tag`` ``sw-`` and BLAKE3-64 over the sorted (stage, digest) lines
of S1-S10; ``world_digest`` BLAKE3-128 over the sorted (stage, digest) lines of every stage, emit included. Each
stage's json lists the sha256 of its files. One genome's contigs are in memory at a time (the placement tables
of every genome are small and held whole).
"""

from __future__ import annotations

import collections
import csv
import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
from blake3 import blake3

from biobox import ids
from biobox.synthworld import assembly, backbone, fusions, intergenic, mobile, seeds, sequences, sources

FORMAT = "synthworld-emit/1"
WORLD_FORMAT = "synthworld-world/1"
UPSTREAM = (("backbone", backbone), ("sequences", sequences), ("replicons", intergenic), ("mobile", mobile),
            ("fusions", fusions), ("assembly", assembly), ("sources", sources))
KIND = {"cds": ids.KIND["cds"], "ncrna": ids.KIND["ncrna"], "crispr_array": ids.KIND["crispr_array"]}
ASSEMBLY_QUALITY = {"complete": "complete", "draft": "contig"}
SOURCE_CLASS = {"complete": "isolate", "draft": "isolate", "mag": "MAG"}
q6 = sequences.q6
revcomp = intergenic.revcomp
_tsv = intergenic._tsv
W = ("world_id", "release_tag")


# --- tables ------------------------------------------------------------------------------------------------------

class TruthGenome(NamedTuple):
    world_id: str
    release_tag: str
    genome_id: str
    species_node: str
    source_label: str
    accession: str
    bioproject: str
    genome_uid: str
    kind: str
    source_class: str
    tax_lineage: str
    tax_lineage_true: str
    mislabeled: bool
    biome: str
    genetic_code: int
    assembly_class: str
    assembly_quality: str
    completeness_true: float
    gc_true: float
    is_phage: bool
    n_contigs: int
    length: int


class TruthContig(NamedTuple):
    world_id: str
    release_tag: str
    contig_uid: str
    accession: str
    source_label: str
    genome_id: str
    contig_id: str
    origin: str
    length: int
    replicon_id: str
    replicon_class: str
    topology: str
    is_circular: bool
    is_contaminant: bool
    is_dropped: bool
    prophage_bearing: bool
    gc: float
    file: str
    record: str


class TruthEdges(NamedTuple):
    world_id: str
    release_tag: str
    contig_uid: str
    genome_id: str
    contig_id: str
    origin: str
    replicon_id: str
    replicon_start: int
    replicon_end: int
    orientation: int
    rotation_offset: int
    left_edge: str
    right_edge: str
    repeat_family_left: str
    repeat_family_right: str
    n_partial_genes: int


class TruthGene(NamedTuple):
    world_id: str
    release_tag: str
    feature_uid: str
    gene_id: str
    contig_uid: str
    start: int
    end: int
    strand: int
    stop_bp: object
    partial_5: bool
    partial_3: bool
    origin_spanning: bool
    feature_kind: str
    family_id: str
    copy_label: str
    protein_uid: str
    domain_arch: str
    module_id: str
    vehicle_id: str
    operon_id: str
    genetic_code: int
    hgt_on_lineage: object
    genome_id: str
    source: str
    callerr_id: str
    status: str
    occurrence: str


class TruthIgr(NamedTuple):
    world_id: str
    release_tag: str
    region_uid: str
    contig_uid: str
    start: object
    end: object
    overlap_bp: int
    igr_uid: str
    left_gene_id: str
    right_gene_id: str
    junction: str
    rule: str
    elements: str
    genome_id: str
    replicon: int
    s5_start: int
    s5_end: int
    status: str
    occurrence: str


class TruthRna(NamedTuple):
    world_id: str
    release_tag: str
    feature_uid: str
    rna_id: str
    rfam_acc: str
    cm_name: str
    clade_emission_id: str
    exp: float
    local_mode: bool
    attached_gene_id: str
    attached_family_id: str
    contig_uid: str
    start: object
    end: object
    strand: object
    dist_to_start: int
    n_subst_from_emission: int
    n_compensatory: int
    ref_bitscore: str
    genome_id: str
    status: str
    occurrence: str


class TruthSpacerTarget(NamedTuple):
    world_id: str
    release_tag: str
    array_feature_uid: str
    array_id: str
    genome_id: str
    contig_uid: str
    array_start: object
    array_end: object
    array_strand: object
    repeat_consensus: str
    cas_module_id: str
    spacer_id: str
    spacer_index: int
    spacer_start: object
    spacer_end: object
    spacer_seq: str
    acquisition_node: str
    target_contig_uid: str
    target_kind: str
    target_id: str
    proto_start: object
    proto_end: object
    proto_strand: object
    n_mismatches: int
    pam_seq: str
    pam_present: bool
    designed: bool
    status: str
    target_status: str


class TruthFusion(NamedTuple):
    world_id: str
    release_tag: str
    feature_uid: str
    fusion_id: str
    kind: str
    family_a: str
    family_b: str
    split_domain_index: int
    clade_node: str
    genome_id: str
    contig_uid: str
    start: object
    end: object
    strand: object
    junction_aa: int
    linker_aa: int
    gene_ids: str
    expected_arch: str
    in_module: str
    status: str
    occurrence: str


class TruthCallerr(NamedTuple):
    world_id: str
    release_tag: str
    feature_uid: str
    callerr_id: str
    kind: str
    genome_id: str
    contig_uid: str
    true_start: object
    true_end: object
    strand: object
    error_pos: object
    alt_start: object
    true_gene_id: str
    true_protein_uid: str
    expected_calls: str
    expected_arch: str
    genetic_code: int
    status: str
    occurrence: str


TruthDup = NamedTuple("TruthDup", [(f, str) for f in W] + list(sources.DupRow.__annotations__.items()))
TruthContam = NamedTuple("TruthContam", [(f, str) for f in W] + list(sources.ContamRow.__annotations__.items())
                         + [("contig_accession", str)])


class TruthFamily(NamedTuple):
    world_id: str
    release_tag: str
    family_id: str
    zombi_family: object
    kind: str
    backbone_kind: str
    domain_hmms: str
    in_scan_panel: bool
    root_source: str
    rate_multiplier: object
    target_identity: object
    realized_mean_identity: object
    realized_min_identity: object
    ref_bitscore_median: str
    n_extant_copies: int
    n_genomes_present: int
    n_transfers: object
    gene_tree_extant: str


class TruthModule(NamedTuple):
    world_id: str
    release_tag: str
    module_id: str
    kind: str
    family_id: str
    member_role: str
    nogd_clade_node: str
    order_index: int
    together: object
    fates: str
    adjacent: bool
    origin_node: str
    strength_class: str
    realized_loci: object
    realized_independent_transfers: str
    realized_vehicle_independent_loci: str
    powered_expected: str


class TruthTrap(NamedTuple):
    world_id: str
    release_tag: str
    trap_id: str
    kind: str
    family_id: str
    vehicle_id: str
    partner_family_ids: str
    expected_label: str
    rationale: str


class IngestRow(NamedTuple):
    source_label: str
    genome_id: str
    accession: str
    bioproject: str
    contig_accession: str
    contig_uid: str
    length: int
    file: str
    record: str
    source_class: str
    replicon_class: str
    molecule_type: str
    topology: str
    genetic_code: int
    kind: str
    assembly_quality: str
    tax_lineage: str
    biome: str


TABLES = {"truth_genome": TruthGenome, "truth_contig": TruthContig, "truth_edges": TruthEdges,
          "truth_gene": TruthGene, "truth_igr": TruthIgr, "truth_rna": TruthRna,
          "truth_spacer_target": TruthSpacerTarget, "truth_fusion": TruthFusion, "truth_callerr": TruthCallerr,
          "truth_dup": TruthDup, "truth_contam": TruthContam, "truth_family": TruthFamily,
          "truth_module": TruthModule, "truth_traps": TruthTrap}


def _cell(v) -> str:
    return "" if v is None else backbone._cell(v)


class _Out:
    """A TSV written row by row; None is an empty cell."""

    def __init__(self, path: Path, cls):
        self.fh = open(path, "w", encoding="utf-8", newline="")
        self.fh.write("\t".join(cls._fields) + "\n")
        self.n = 0

    def row(self, r) -> None:
        self.fh.write("\t".join(_cell(v) for v in r) + "\n")
        self.n += 1

    def close(self) -> None:
        self.fh.close()


def check_config(c: dict) -> None:
    if not 0.0 <= float(c["mislabel_rate"]) <= 1.0:
        raise ValueError("emit.mislabel_rate is a probability")
    if not c["biomes"] or len(set(c["biomes"])) != len(c["biomes"]):
        raise ValueError("emit.biomes is a non-empty list of distinct names")
    if c["biome_rank"] not in ("phylum", "class", "order", "family", "genus", "species"):
        raise ValueError(f"emit.biome_rank is a taxonomy rank, got {c['biome_rank']!r}")


def stop_bp(start: int, end: int, strand: int, partial_3: bool):
    """The slot join key (§3.1): just past the stop on +, the stop codon's lowest base on -; None when cut off."""
    if partial_3:
        return None
    return end if strand == 1 else start


def assembly_quality(klass: str, completeness: float) -> str:
    if klass == "mag":
        return "mag_hq" if completeness >= 0.9 else "mag_mq"
    return ASSEMBLY_QUALITY[klass]


def _rank(lineage: str, prefix: str) -> str:
    return next((x for x in lineage.split(";") if x.startswith(prefix)), "")


# --- placement ---------------------------------------------------------------------------------------------------

class Occ(NamedTuple):
    """An emitted contig derived from an S9 contig: the contig itself, a near duplicate's copy, or a contaminant
    window [x0, y0) of it (``o`` -1 when reverse-complemented)."""
    uid: str
    origin: str
    x0: int
    y0: int
    o: int
    host: str
    contig_id: str


class Place(NamedTuple):
    occ: Occ
    start: int
    end: int
    strand: int
    p5: bool
    p3: bool
    length: int


def window(cs: int, ce: int, cst: int, p5: bool, p3: bool, L: int, x0: int, y0: int, o: int):
    """The parts of a feature [cs, ce) (strand cst, end > L across a circular origin) on a contig of L bases inside
    the window [x0, y0), in window coordinates (reverse-complemented when o is -1), with the 5'/3' flags."""
    out = []
    W_ = y0 - x0
    for sh in (0, -L):
        fs, fe = cs + sh, ce + sh
        lo, hi = max(fs, x0), min(fe, y0)
        if lo >= hi:
            continue
        lc, hc = fs < x0, fe > y0
        q5, q3 = (p5 or lc, p3 or hc) if cst == 1 else (p5 or hc, p3 or lc)
        a, b = lo - x0, hi - x0
        out.append((a, b, cst, q5, q3) if o == 1 else (W_ - b, W_ - a, -cst, q5, q3))
    return out


class Placer:
    """Lifts features of the S5, S6, S7 and S8 frames to the S8 frame and places them on every emitted contig."""

    def __init__(self, world: Path, contam: list[dict], dup: list[dict]):
        self.world = world
        self.rows = collections.defaultdict(list)
        self.by_id = {}
        for r in _tsv(world / "assembly" / "contigs.tsv"):
            self.rows[r["record"]].append(r)
            self.by_id[r["contig_id"]] = r
        reps = {r["record"]: r for r in _tsv(world / "assembly" / "replicons.tsv")}
        self.recs = {}
        for rec, rows in self.rows.items():
            rows.sort(key=lambda r: int(r["replicon_start"]))
            ps = [assembly.Piece(0, int(r["replicon_start"]), int(r["replicon_end"]), r["is_circular"] == "true",
                                 "", "", "", "") for r in rows]
            self.recs[rec] = (ps, int(reps[rec]["length"]), reps[rec]["topology"] == "circular", rows)
        self.occs = collections.defaultdict(list)
        for r in self.by_id.values():
            if r["emitted"] == "true":
                L = int(r["length"])
                self.occs[r["contig_id"]].append(Occ(r["contig_uid"], "assembly", 0, L, 1, r["genome_id"],
                                                     r["contig_id"]))
        for d in dup:
            if d["kind"] == "near":
                L = int(self.by_id[d["rep_contig_id"]]["length"])
                self.occs[d["rep_contig_id"]].append(Occ(d["member_contig_uid"], "near_dup", 0, L, 1,
                                                         d["member_genome_id"], d["member_contig_id"]))
        for c in contam:
            self.occs[c["donor_contig_id"]].append(Occ(c["contig_uid"], "contaminant", int(c["donor_start"]),
                                                       int(c["donor_end"]), int(c["orientation"]),
                                                       c["host_genome_id"], c["contig_id"]))
        for v in self.occs.values():
            v.sort(key=lambda o: (o.origin, o.contig_id))
        # the lifts, per record
        self.len5 = {r["record"]: int(r["length"]) for r in _tsv(world / "replicons" / "replicons.tsv")}
        self.len6 = {r["record"]: int(r["length"]) for r in _tsv(world / "mobile" / "replicons.tsv")}
        self.len7 = {r["record"]: int(r["length"]) for r in _tsv(world / "fusions" / "replicons.tsv")}
        self.ins6 = collections.defaultdict(list)
        for r in _tsv(world / "mobile" / "insertions.tsv"):
            self.ins6[r["record"]].append((int(r["at"]), int(r["length"])))
        self.ed7 = collections.defaultdict(list)
        for r in _tsv(world / "fusions" / "edits.tsv"):
            self.ed7[r["record"]].append((int(r["start"]), int(r["end"]), int(r["length"])))
        self.ed8 = collections.defaultdict(list)
        for r in _tsv(world / "assembly" / "edits.tsv"):
            self.ed8[r["record"]].append((int(r["start"]), int(r["end"]), int(r["length"])))
        self._lifts: dict[str, tuple] = {}

    def _lift_objs(self, rec: str):
        if rec not in self._lifts:
            self._lifts[rec] = (mobile.Lift(self.len5[rec], self.ins6[rec]),
                                assembly.Lift(self.len6[rec], self.ed7[rec]),
                                assembly.Lift(self.len7[rec], self.ed8[rec]))
        return self._lifts[rec]

    def lift(self, rec: str, stage: int, s: int, e: int):
        """[s, e) of the stage-``stage`` frame (5 to 8) on the S8 frame, or the reason it is not there."""
        l56, l67, l78 = self._lift_objs(rec)
        if stage <= 5:
            try:
                s, e = l56.feature(s, e)
            except ValueError:
                return "s6_insertion"
        if stage <= 6:
            if l67.touched(s, e):
                return "s7_event"
            s, e = l67.feature(s, e)
        if stage <= 7:
            if l78.touched(s, e):
                return "s8_edit"
            s, e = l78.feature(s, e)
        return s, e

    def frame_parts(self, rec: str, s: int, e: int, strand: int):
        """The parts of S8-frame feature [s, e) on the S9 contigs: (contig row, start, end, strand, 5' cut, 3' cut)."""
        ps, N, circ, rows = self.recs[rec]
        out = []
        for j, x, y, cl, cr in assembly.place(ps, N, circ, s, e):
            row = rows[j]
            cs, ce, cst = assembly.to_contig(int(row["length"]), int(row["orientation"]), x, y, strand)
            p5, p3 = (cl, cr) if strand == 1 else (cr, cl)
            out.append((row, cs, ce, cst, p5, p3))
        return out

    def derive(self, row: dict, cs: int, ce: int, cst: int, p5: bool, p3: bool) -> list[Place]:
        L = int(row["length"])
        out = []
        for oc in self.occs.get(row["contig_id"], ()):
            if oc.origin != "contaminant":
                out.append(Place(oc, cs, ce, cst, p5, p3, L))
            else:
                for a, b, st, q5, q3 in window(cs, ce, cst, p5, p3, L, oc.x0, oc.y0, oc.o):
                    out.append(Place(oc, a, b, st, q5, q3, oc.y0 - oc.x0))
        return out

    def place(self, rec: str, s: int, e: int, strand: int) -> tuple[str, list[Place]]:
        """(status, emitted pieces) of S8-frame feature [s, e): whole, cut or dropped (module docstring)."""
        parts = self.frame_parts(rec, s, e, strand)
        em = [p for p in parts if p[0]["emitted"] == "true"]
        if not em:
            return "dropped", []
        status = "whole" if len(parts) == 1 and not (parts[0][4] or parts[0][5]) else "cut"
        return status, [pl for p in em for pl in self.derive(*p)]

    def lifted(self, rec: str, stage: int, s: int, e: int, strand: int) -> tuple[str, list[Place]]:
        x = self.lift(rec, stage, s, e)
        if isinstance(x, str):
            return x, []
        return self.place(rec, x[0], x[1], strand)


def piece_status(pl: Place) -> str:
    return "cut" if pl.p5 or pl.p3 else "whole"


def point(placer: Placer, rec: str, x: int, y: int, uid: str):
    """Lowest contig coordinate of S8-frame interval [x, y) on emitted contig ``uid`` when it lies there whole."""
    _, pls = placer.place(rec, x, y, 1)
    got = [pl for pl in pls if pl.occ.uid == uid and not (pl.p5 or pl.p3)]
    return got[0].start if len(got) == 1 else None


# --- sequences of one genome's emitted contigs --------------------------------------------------------------------

class _Seqs:
    def __init__(self, world: Path, contaminants: dict[str, str]):
        self.world = world
        self.gid = None
        self.seqs: dict[str, str] = {}
        self.contaminants = contaminants

    def genome(self, gid: str) -> dict[str, str]:
        if gid != self.gid:
            self.gid, self.seqs = gid, {}
            for name, s in intergenic._fasta_records(self.world / "assembly" / "contigs" / f"{gid}.fna"):
                self.seqs[name] = s
            for f in sorted((self.world / "sources" / "duplicates").glob(f"{gid}_*.fna")):
                self.seqs.update(dict(intergenic._fasta_records(f)))
            self.seqs.update(self.contaminants)
        return self.seqs


# --- build -------------------------------------------------------------------------------------------------------

@dataclass
class EmitWorld:
    seed: int
    world_id: str
    release_tag: str
    digest: str
    world_digest: str
    out: Path
    meta: dict = field(default_factory=dict, repr=False)


def _digest(d: dict[str, str], length: int) -> str:
    return blake3("".join(f"{k}\t{v}\n" for k, v in sorted(d.items())).encode()).hexdigest(length=length)


def build(world) -> EmitWorld:
    """Write ``world/emit/`` and ``world/manifest.json`` (module docstring) from a world whose S1-S10 are built."""
    world = Path(world)
    docs = {}
    for d, mod in UPSTREAM:
        docs[d] = json.loads((world / d / f"{d}.json").read_text(encoding="utf-8"))
        if docs[d]["format"] != mod.FORMAT:
            raise ValueError(f"{world}/{d}: format {docs[d]['format']!r}, expected {mod.FORMAT!r}")
    bb = docs["backbone"]
    preset = bb["preset"]
    cfg = preset["emit"]
    check_config(cfg)
    seed = int(bb["seed"])
    out = world / "emit"
    if out.exists() or (world / "manifest.json").exists():
        raise FileExistsError(f"{out} or {world}/manifest.json exists; the emit stage is written once")
    tmp = world / "emit.partial"
    if tmp.exists():
        shutil.rmtree(tmp)
    (tmp / "fasta").mkdir(parents=True)
    (tmp / "truth").mkdir()
    world_id = f"{preset['name']}-{seed}"
    release_tag = "sw-" + _digest({d: docs[d]["digest"] for d, _ in UPSTREAM}, 8)
    WV = (world_id, release_tag)
    counts = collections.Counter()

    # --- the earlier stages' small tables
    srcs = _tsv(world / "sources" / "sources.tsv")
    occ_rows = _tsv(world / "sources" / "source_contigs.tsv")
    contam = _tsv(world / "sources" / "contam.tsv")
    dup = _tsv(world / "sources" / "dup.tsv")
    genomes = {r["genome_id"]: r for r in _tsv(world / "assembly" / "genomes.tsv")}
    reps8 = {r["record"]: r for r in _tsv(world / "assembly" / "replicons.tsv")}
    host_genomes = {r["genome_id"]: r for r in _tsv(world / "backbone" / "genomes.tsv")}
    lineage = {r["genome_id"]: r["tax_lineage"] for r in _tsv(world / "backbone" / "taxonomy.tsv")}
    phages = _tsv(world / "mobile" / "phages.tsv")
    ph_world = world / "mobile" / "phages"
    ph_lineage = {r["genome_id"]: r["tax_lineage"] for r in _tsv(ph_world / "backbone" / "taxonomy.tsv")}
    placer = Placer(world, contam, dup)
    contig_of = placer.by_id

    # --- accessions of the phage sources, after the host sources'
    n_src = len(srcs)
    n_cacc = len(occ_rows)
    n_bp = max(int(r["bioproject"][5:]) for r in srcs)
    phage_src = []
    for k, p in enumerate(phages):
        phage_src.append({"source_label": "src_a", "genome_id": p["phage_id"], "kind": "phage", "rep_genome_id": "",
                          "accession": f"GCA_SW{n_src + k + 1:06d}.1", "bioproject": f"PRJSW{n_bp + 1:05d}",
                          "contig_accession": f"SWC{n_cacc + k + 1:08d}.1", "phage": p})

    # --- taxonomy labels and biomes (S11/mislabel, S11/biome)
    rank_prefix = {"phylum": "p__", "class": "c__", "order": "o__", "family": "f__", "genus": "g__",
                   "species": "s__"}[cfg["biome_rank"]]
    gids = sorted(genomes)
    biome = {}
    for g in gids:
        taxon = _rank(lineage[g], rank_prefix)
        r = seeds.generator(seed, f"S11/biome/{taxon}")
        biome[g] = str(cfg["biomes"][int(r.integers(len(cfg["biomes"])))])
    label = {}
    for s in srcs:
        g = s["genome_id"]
        r = seeds.generator(seed, f"S11/mislabel/{s['source_label']}/{g}")
        lab = lineage[g]
        if r.random() < float(cfg["mislabel_rate"]):
            other = [h for h in gids if _rank(lineage[h], "g__") != _rank(lineage[g], "g__")]
            if other:
                lab = lineage[other[int(r.integers(len(other)))]]
        label[(s["source_label"], g)] = lab

    # --- the occurrences: replicon class, topology, the source of each contig
    rep_of_member = {d["member_contig_id"]: d["rep_contig_id"] for d in dup}
    donor_of = {c["contig_id"]: c["donor_contig_id"] for c in contam}

    def base_contig(occ: dict) -> dict:
        cid = occ["contig_id"]
        if occ["origin"] in ("exact_dup", "near_dup"):
            cid = rep_of_member[cid]
        elif occ["origin"] == "contaminant":
            cid = donor_of[cid]
        return contig_of[cid]

    def replicon_id(record: str) -> str:
        k = int(record.rsplit("_r", 1)[1])
        return "chr" if k == 0 else f"p{k}"

    # --- FASTA, ingest manifest, build-stage inputs
    accs_of_uid = collections.defaultdict(list)       # contig_uid -> contig accessions (exact dups: several)
    ingest_rows: list[IngestRow] = []
    src_contigs = collections.defaultdict(list)
    for o in occ_rows:
        src_contigs[(o["source_label"], o["genome_id"])].append(o)
    contaminant_seqs = dict(intergenic._fasta_records(world / "sources" / "contaminants.fna"))
    cont_by_uid = {}
    for c in contam:
        cont_by_uid[c["contig_uid"]] = contaminant_seqs[c["contig_id"]]
    genome_entries, listed = [], set()
    objects = []
    truth_genome_rows = []
    gc_of = collections.defaultdict(lambda: [0, 0])
    for r in reps8.values():
        gc_of[r["genome_id"]][0] += float(r["gc"]) * int(r["length"])
        gc_of[r["genome_id"]][1] += int(r["length"])
    phage_seqs = dict(intergenic._fasta_records(world / "mobile" / "phages.fna"))
    gc_acc = {}
    for s in srcs:
        g, lab = s["genome_id"], s["source_label"]
        gr = genomes[g]
        code = int(gr["genetic_code"])
        klass = gr["assembly_class"]
        aq = assembly_quality(klass, float(gr["completeness"]))
        rel = f"fasta/{s['accession']}.fna"
        cache: dict[str, dict[str, str]] = {}
        fx = sequences._Fasta(tmp / rel)
        uids, circ = [], []
        for o in src_contigs[(lab, g)]:
            if o["file"] not in cache:
                cache[o["file"]] = dict(intergenic._fasta_records(world / o["file"]))
            seq = cache[o["file"]][o["record"]]
            uid = ids.contig_uid(seq.encode("ascii")).hex()
            if uid != o["contig_uid"] or len(seq) != int(o["length"]):
                raise AssertionError(f"{o['contig_accession']}: bytes disagree with source_contigs.tsv")
            fx.write(o["contig_accession"], seq)
            gc_acc[o["contig_accession"]] = q6((seq.count("G") + seq.count("C")) / len(seq))
            base = base_contig(o)
            topo = "circular" if base["is_circular"] == "true" and o["origin"] != "contaminant" else "linear"
            rclass = reps8[base["record"]]["kind"]
            accs_of_uid[uid].append(o["contig_accession"])
            ingest_rows.append(IngestRow(lab, g, s["accession"], s["bioproject"], o["contig_accession"], uid,
                                         len(seq), rel, o["contig_accession"], SOURCE_CLASS[klass], rclass, "DNA",
                                         topo, code, "prokaryote", aq, label[(lab, g)], biome[g]))
            if uid not in uids:
                uids.append(uid)
                if topo == "circular":
                    circ.append(uid)
        fx.close()
        guid = ids.genome_uid([bytes.fromhex(u) for u in uids]).hex()
        if guid not in listed:
            listed.add(guid)
            genome_entries.append({"contig_uids": uids, "genetic_code": code, "kind": "prokaryote",
                                   "circular": circ, "assembly_quality": aq})
        else:
            counts["genomes_listed_once"] += 1
        objects.append({"path": rel, "source_name": str(cfg["source_name"]), "source_release": world_id,
                        "license": str(cfg["license"]), "redistributable": True})
        gc = gc_of[g]
        truth_genome_rows.append(TruthGenome(*WV, g, host_genomes[g]["species_node"], lab, s["accession"],
                                             s["bioproject"], guid, "prokaryote", SOURCE_CLASS[klass],
                                             label[(lab, g)], lineage[g], label[(lab, g)] != lineage[g], biome[g],
                                             code, klass, aq, float(gr["completeness"]), q6(gc[0] / gc[1]), False,
                                             len(src_contigs[(lab, g)]),
                                             sum(int(o["length"]) for o in src_contigs[(lab, g)])))
    for ps in phage_src:
        p = ps["phage"]
        rel = f"fasta/{ps['accession']}.fna"
        seq = phage_seqs[p["record"]]
        uid = ids.contig_uid(seq.encode("ascii")).hex()
        if uid != p["contig_uid"]:
            raise AssertionError(f"{p['phage_id']}: bytes disagree with phages.tsv")
        fx = sequences._Fasta(tmp / rel)
        fx.write(ps["contig_accession"], seq)
        fx.close()
        accs_of_uid[uid].append(ps["contig_accession"])
        tax = ";".join(["d__Viruses"] + ph_lineage[p["genome_id"]].split(";")[1:])
        ingest_rows.append(IngestRow("src_a", p["phage_id"], ps["accession"], ps["bioproject"],
                                     ps["contig_accession"], uid, len(seq), rel, ps["contig_accession"],
                                     "MGE database", "virus", "DNA", "linear", 11, "virus", "complete", tax, ""))
        guid = ids.genome_uid([bytes.fromhex(uid)]).hex()
        if guid not in listed:
            listed.add(guid)
            genome_entries.append({"contig_uids": [uid], "genetic_code": 11, "kind": "virus", "circular": [],
                                   "assembly_quality": "complete"})
        objects.append({"path": rel, "source_name": str(cfg["source_name"]), "source_release": world_id,
                        "license": str(cfg["license"]), "redistributable": True})
        truth_genome_rows.append(TruthGenome(*WV, p["phage_id"], p["species_node"], "src_a", ps["accession"],
                                             ps["bioproject"], guid, "virus", "MGE database", tax, tax, False, "",
                                             11, "complete", "complete", 1.0, float(p["gc"]), True, 1, len(seq)))
    t = _Out(tmp / "ingest.tsv", IngestRow)
    for r in ingest_rows:
        t.row(r)
    t.close()
    (tmp / "genomes.json").write_text(json.dumps(genome_entries, indent=1) + "\n", encoding="utf-8")
    (tmp / "objects.json").write_text(json.dumps(objects, indent=1) + "\n", encoding="utf-8")
    t = _Out(tmp / "truth" / "truth_genome.tsv", TruthGenome)
    for r in truth_genome_rows:
        t.row(r)
    t.close()
    shared = collections.Counter()
    for e in genome_entries:
        for u in e["contig_uids"]:
            shared[u] += 1
    counts["contigs_in_several_genomes"] = sum(1 for v in shared.values() if v > 1)
    counts["sources"] = len(srcs) + len(phage_src)
    counts["genomes_listed"] = len(genome_entries)

    gff = open(tmp / "truth.gff", "w", encoding="utf-8", newline="")
    gff.write("##gff-version 3\n")
    for r in ingest_rows:
        gff.write(f"##sequence-region {r.contig_accession} 1 {r.length}\n")
    for r in ingest_rows:
        if r.topology == "circular":     # GFF3: a feature may end past the end of a circular region
            gff.write(f"{r.contig_accession}\tsynthworld\tregion\t1\t{r.length}\t.\t+\t.\t"
                      f"ID={r.contig_accession};Is_circular=true\n")

    def gff_rows(uid, start, end, strand, ftype, ident, attrs):
        for acc in accs_of_uid[uid]:
            gff.write(f"{acc}\tsynthworld\t{ftype}\t{start + 1}\t{end}\t.\t{'+' if strand == 1 else '-'}\t"
                      f"{'0' if ftype == 'CDS' else '.'}\tID={ident}@{acc};{attrs}\n")

    # --- genes (truth_gene)
    seq_genes = {}
    with open(world / "sequences" / "genes.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq_genes[r["copy_label"]] = (r["hgt_on_lineage"], r["domain_arch"])
    ph_genes = {}
    with open(ph_world / "sequences" / "genes.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ph_genes[r["copy_label"]] = (r["hgt_on_lineage"], r["domain_arch"])
    ev_genes = {r["copy_label"]: r for r in _tsv(world / "fusions" / "genes.tsv") if r["source"] in
                ("fused", "split", "nogd")}
    group_kind = {r["name"]: r["kind"] for r in _tsv(world / "backbone" / "groups.tsv")}
    unit_group = {}
    with open(world / "backbone" / "genes.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["group"]:
                unit_group[r["copy_label"]] = r["group"]
    operon = {}
    with open(world / "replicons" / "genes.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            operon[r["copy_label"]] = r["operon_id"]

    def gene_attrs(label: str, source: str):
        """(hgt_on_lineage, domain_arch, module_id, vehicle_id, operon_id) of an S8 gene."""
        if source == "prophage":
            base = label.split("@")[0][3:]
            hgt, arch = ph_genes.get(base, (None, ""))
            return hgt, arch, "", "", ""
        if label in ev_genes:
            arch = ev_genes[label]["domain_arch"]
            parents = ev_genes[label]["parents"].split(";") if ev_genes[label]["parents"] else []
            grp = unit_group.get(parents[0], "") if parents else ""
            hgt = None
        else:
            hgt, arch = seq_genes.get(label, (None, ""))
            grp = unit_group.get(label, "")
        kind = group_kind.get(grp, "")
        return (hgt, arch, grp if kind in ("planted", "profile_only") else "", grp if kind == "vehicle" else "",
                operon.get(label, ""))

    tg = _Out(tmp / "truth" / "truth_gene.tsv", TruthGene)
    n_partial = collections.Counter()
    fam_copies = collections.Counter()
    fam_genomes = collections.defaultdict(set)
    gene_stream = assembly._Stream(world / "assembly" / "genes.tsv")
    for gid, pieces in mobile._by_genome(world / "assembly" / "gene_contigs.tsv"):
        seen_uid_genes = set()
        g8 = {r["copy_label"]: r for r in gene_stream.rows(gid)}
        for r in g8.values():
            fam_copies[r["family_id"]] += 1
            fam_genomes[r["family_id"]].add(gid)
        rows = []
        for p in pieces:
            if p["emitted"] != "true":
                continue
            g = g8[p["copy_label"]]
            row = contig_of[p["contig_id"]]
            for pl in placer.derive(row, int(p["start"]), int(p["end"]), int(p["strand"]), p["partial_5"] == "true",
                                    p["partial_3"] == "true"):
                rows.append((pl, g))
        out_rows = []
        for pl, g in rows:
            key = (pl.occ.uid, pl.start, pl.end, pl.strand, g["copy_label"])
            if key in seen_uid_genes:
                continue
            seen_uid_genes.add(key)
            hgt, arch, mod, veh, op = gene_attrs(g["copy_label"], g["source"])
            fuid = ids.feature_uid(bytes.fromhex(pl.occ.uid), pl.start, pl.end, pl.strand, KIND["cds"]).hex()
            out_rows.append(TruthGene(*WV, fuid, g["copy_label"], pl.occ.uid, pl.start, pl.end, pl.strand,
                                      stop_bp(pl.start, pl.end, pl.strand, pl.p3), pl.p5, pl.p3, pl.end > pl.length,
                                      g["feature_kind"], g["family_id"], g["copy_label"], g["protein_uid"], arch,
                                      mod, veh, op, int(g["genetic_code"]), hgt, gid, g["source"], g["callerr_id"],
                                      piece_status(pl), pl.occ.origin))
            if pl.p5 or pl.p3:
                n_partial[pl.occ.uid] += 1
        out_rows.sort(key=lambda x: (x.contig_uid, x.start, x.end, x.strand, x.gene_id))
        for x in out_rows:
            tg.row(x)
            gff_rows(x.contig_uid, x.start, x.end, x.strand, "CDS" if x.feature_kind != "pseudogene" else "pseudogene",
                     x.gene_id, f"gene_id={x.gene_id};family_id={x.family_id};feature_kind={x.feature_kind};"
                                f"partial_5={_cell(x.partial_5)};partial_3={_cell(x.partial_3)}")
            counts[f"genes_{x.occurrence}_{x.status}"] += 1
    ph_by_id = {ps["genome_id"]: ps for ps in phage_src}
    ph_rows = []
    for r in _tsv(world / "mobile" / "phage_genes.tsv"):
        p = ph_by_id[r["phage_id"]]["phage"]
        uid, s, e, st = p["contig_uid"], int(r["start"]), int(r["end"]), int(r["strand"])
        hgt, arch = ph_genes.get(r["copy_label"][3:], (None, ""))
        fam_copies[r["family_id"]] += 1
        fam_genomes[r["family_id"]].add(r["phage_id"])
        fuid = ids.feature_uid(bytes.fromhex(uid), s, e, st, KIND["cds"]).hex()
        ph_rows.append(TruthGene(*WV, fuid, r["copy_label"], uid, s, e, st, stop_bp(s, e, st, False), False, False,
                                 False, "CDS", r["family_id"], r["copy_label"], r["protein_uid"], arch, "", "", "", 11,
                                 hgt, r["phage_id"], "phage", "", "whole", "phage"))
    ph_rows.sort(key=lambda x: (x.contig_uid, x.start, x.end, x.strand, x.gene_id))
    for x in ph_rows:
        tg.row(x)
        gff_rows(x.contig_uid, x.start, x.end, x.strand, "CDS", x.gene_id,
                 f"gene_id={x.gene_id};family_id={x.family_id};feature_kind=CDS")
        counts["genes_phage_whole"] += 1
    tg.close()

    # --- call errors (truth_callerr)
    tc = _Out(tmp / "truth" / "truth_callerr.tsv", TruthCallerr)
    for r in _tsv(world / "assembly" / "callerr.tsv"):
        rec, s, e, st = r["record"], int(r["true_start"]), int(r["true_end"]), int(r["strand"])
        status, pls = placer.place(rec, s, e, st)
        whole = [pl for pl in pls if not (pl.p5 or pl.p3)]
        base = (r["callerr_id"], r["kind"], r["genome_id"])
        if not whole:
            tc.row(TruthCallerr(*WV, "", *base, "", None, None, None, None, None, r["true_gene_id"],
                                r["true_protein_uid"], "", r["expected_arch"], int(r["genetic_code"]),
                                "dropped" if status == "dropped" else "partial", ""))
            counts[f"callerr_{'dropped' if status == 'dropped' else 'partial'}"] += 1
            continue
        N8 = int(reps8[rec]["length"])
        for pl in whole:
            uid = pl.occ.uid
            ep = int(r["error_pos"])
            err = point(placer, rec, ep, ep + 1, uid)
            alt = None
            if int(r["alt_start"]) >= 0:
                a = int(r["alt_start"]) % N8
                alt = point(placer, rec, a, a + 3, uid)
            calls = []
            for item in filter(None, r["expected_calls"].split(";")):
                span, cst = item.split(":")
                x, y = (int(v) for v in span.split("-"))
                if x >= N8:
                    x, y = x - N8, y - N8
                _, cp = placer.place(rec, x, y, int(cst))
                got = [c for c in cp if c.occ.uid == uid and not (c.p5 or c.p3)]
                if len(got) == 1:
                    calls.append(f"{got[0].start}-{got[0].end}:{got[0].strand}")
            fuid = ids.feature_uid(bytes.fromhex(uid), pl.start, pl.end, pl.strand, KIND["cds"]).hex()
            tc.row(TruthCallerr(*WV, fuid, *base, uid, pl.start, pl.end, pl.strand, err, alt, r["true_gene_id"],
                                r["true_protein_uid"], ";".join(calls), r["expected_arch"], int(r["genetic_code"]),
                                "whole", pl.occ.origin))
            counts[f"callerr_whole_{pl.occ.origin}"] += 1
    tc.close()

    # --- junctions and their elements (truth_igr)
    seqs = _Seqs(world, contaminant_seqs)
    elements = assembly._Stream(world / "replicons" / "elements.tsv")
    ti = _Out(tmp / "truth" / "truth_igr.tsv", TruthIgr)
    for gid, juncs in mobile._by_genome(world / "replicons" / "junctions.tsv"):
        els = collections.defaultdict(list)
        for r in elements.rows(gid):
            els[(r["replicon"], r["left_position"])].append(r)
        gseq = seqs.genome(gid)
        uid_seq = {}
        for cid in list(gseq):
            uid_seq[ids.contig_uid(gseq[cid].encode("ascii")).hex()] = gseq[cid]
        rows = []
        for j in juncs:
            rec = f"{gid}_r{j['replicon']}"
            s, e, ov = int(j["start"]), int(j["end"]), int(j["overlap_bp"])
            N5 = placer.len5[rec]
            zero = e == s
            if zero:
                fs, fe = (s - ov, s) if ov else (s - 1, s)
                if fs < 0:
                    fs, fe = fs + N5, fe + N5
            else:
                fs, fe = s, e
            st_, pls = placer.lifted(rec, 5, fs, fe, 1)
            common = (j["junction"], j["rule"])
            if not pls:
                rows.append(TruthIgr(*WV, "", "", None, None, ov, "", j["left_copy"], j["right_copy"], *common, "",
                                     gid, int(j["replicon"]), s, e, st_, ""))
                counts[f"igr_{st_}"] += 1
                continue
            placed_els = collections.defaultdict(list)
            for el in els[(j["replicon"], j["left_position"])]:
                _, epl = placer.lifted(rec, 5, int(el["start"]), int(el["end"]), int(el["strand"]))
                for x in epl:
                    if not (x.p5 or x.p3):
                        placed_els[x.occ.uid].append((x.start, x.end, x.strand, el["kind"]))
            for pl in pls:
                flip = pl.strand == -1          # placed on the frame's + strand: -1 means the contig is reversed
                left, right = (j["right_copy"], j["left_copy"]) if flip else (j["left_copy"], j["right_copy"])
                if zero:
                    a = b = (pl.start if flip else pl.end) if not ov else pl.end
                    igr_uid, region = "", ""
                    el_txt = ""
                else:
                    a, b = pl.start, pl.end
                    sq = uid_seq[pl.occ.uid]
                    bases = assembly.cut(sq, a, b, 1)
                    igr_uid = ids.contig_uid(bases.encode("ascii")).hex()
                    region = ids.region_uid(bytes.fromhex(pl.occ.uid), a, b, 1).hex()
                    el_txt = ";".join(f"{k}:{x}-{y}:{z}" for x, y, z, k in sorted(placed_els[pl.occ.uid])
                                      if x < b and y > a)
                rows.append(TruthIgr(*WV, region, pl.occ.uid, a, b, ov, igr_uid, left, right, *common, el_txt, gid,
                                     int(j["replicon"]), s, e, piece_status(pl), pl.occ.origin))
                counts[f"igr_{pl.occ.origin}_{piece_status(pl)}"] += 1
        for x in rows:
            ti.row(x)
    ti.close()

    # --- RNA motifs (truth_rna)
    tr = _Out(tmp / "truth" / "truth_rna.tsv", TruthRna)
    for r in _tsv(world / "replicons" / "rna.tsv"):
        rec = f"{r['genome_id']}_r{r['replicon']}"
        st_, pls = placer.lifted(rec, 5, int(r["start"]), int(r["end"]), int(r["strand"]))
        base = (r["rna_id"], r["model_acc"], r["model_name"], r["emission_id"], float(r["exp"]),
                r["local_mode"] == "true", r["attached_copy"], r["family_id"])
        tail = (int(r["dist_to_start"]), int(r["n_subst_from_emission"]), int(r["n_compensatory"]), "",
                r["genome_id"])
        if not pls:
            tr.row(TruthRna(*WV, "", *base, "", None, None, None, *tail, st_, ""))
            counts[f"rna_{st_}"] += 1
            continue
        for pl in pls:
            fuid = ids.feature_uid(bytes.fromhex(pl.occ.uid), pl.start, pl.end, pl.strand, KIND["ncrna"]).hex()
            tr.row(TruthRna(*WV, fuid, *base, pl.occ.uid, pl.start, pl.end, pl.strand, *tail, piece_status(pl),
                            pl.occ.origin))
            counts[f"rna_{pl.occ.origin}_{piece_status(pl)}"] += 1
            if piece_status(pl) == "whole":
                gff_rows(pl.occ.uid, pl.start, pl.end, pl.strand, "ncRNA", r["rna_id"],
                         f"rna_id={r['rna_id']};model={r['model_acc']}")
    tr.close()
    gff.close()

    # --- spacers and targets (truth_spacer_target)
    phage_uid = {p["phage_id"]: p["contig_uid"] for p in phages}
    ts = _Out(tmp / "truth" / "truth_spacer_target.tsv", TruthSpacerTarget)
    kind_map = {"phage": "phage_contig", "prophage": "prophage", "none": "none"}
    for r in _tsv(world / "mobile" / "spacer_targets.tsv"):
        rec, ast = r["record"], int(r["array_strand"])
        a_st, a_pls = placer.lifted(rec, 6, int(r["array_start"]), int(r["array_end"]), ast)
        s_st, s_pls = placer.lifted(rec, 6, int(r["spacer_start"]), int(r["spacer_end"]), ast)
        # the target side
        tk = r["target_kind"]
        if tk == "phage":
            targets = [(phage_uid[r["target_id"]], int(r["proto_start"]), int(r["proto_end"]),
                        int(r["proto_strand"]), "whole")]
        elif tk == "prophage":
            t_st, t_pls = placer.lifted(r["target_record"], 6, int(r["proto_start"]), int(r["proto_end"]),
                                        int(r["proto_strand"]))
            targets = [(x.occ.uid, x.start, x.end, x.strand, piece_status(x)) for x in t_pls] or \
                [("", None, None, None, t_st)]
        else:
            targets = [("", None, None, None, "none")]
        hosts = []
        for x in a_pls:
            sp = [y for y in s_pls if y.occ.uid == x.occ.uid and not (y.p5 or y.p3)]
            if not (x.p5 or x.p3) and len(sp) == 1:
                hosts.append((x, sp[0]))
        if not hosts:
            why = a_st if a_st not in ("whole",) else "cut"
            hosts = [(None, None)]
        for x, sp in hosts:
            for tu, ps_, pe_, pst, tstat in targets:
                if x is None:
                    host = ("", "", None, None, None)
                    spc = (None, None)
                    status = why
                else:
                    fuid = ids.feature_uid(bytes.fromhex(x.occ.uid), x.start, x.end, x.strand,
                                           KIND["crispr_array"]).hex()
                    host = (fuid, x.occ.uid, x.start, x.end, x.strand)
                    spc = (sp.start, sp.end)
                    status = "whole"
                ts.row(TruthSpacerTarget(*WV, host[0], r["array_id"], r["genome_id"], host[1], host[2], host[3],
                                         host[4], r["repeat_consensus"], r["cas_module_id"], r["spacer_id"],
                                         int(r["spacer_index"]), spc[0], spc[1], r["spacer_seq"],
                                         r["acquisition_node"], tu, kind_map[tk], r["target_id"], ps_, pe_, pst,
                                         int(r["n_mismatches"]), r["pam_seq"], r["pam_present"] == "true",
                                         r["designed"] == "true", status, tstat))
                counts[f"spacer_{status}_{tstat}"] += 1
    ts.close()

    # --- fusions (truth_fusion)
    tf = _Out(tmp / "truth" / "truth_fusion.tsv", TruthFusion)
    for r in _tsv(world / "fusions" / "fusions.tsv"):
        st_, pls = placer.lifted(r["record"], 7, int(r["start"]), int(r["end"]), int(r["strand"]))
        base = (r["fusion_id"], r["kind"], r["family_a"], r["family_b"], int(r["split_domain_index"]),
                r["clade_node"], r["genome_id"])
        tail = (int(r["junction_aa"]), int(r["linker_aa"]), r["gene_ids"], r["expected_arch"], r["in_module"])
        if not pls:
            tf.row(TruthFusion(*WV, "", *base, "", None, None, None, *tail, st_, ""))
            counts[f"fusion_{st_}"] += 1
            continue
        for pl in pls:
            fuid = ids.feature_uid(bytes.fromhex(pl.occ.uid), pl.start, pl.end, pl.strand, KIND["cds"]).hex() \
                if r["kind"] == "fusion" else ""
            tf.row(TruthFusion(*WV, fuid, *base, pl.occ.uid, pl.start, pl.end, pl.strand, *tail, piece_status(pl),
                               pl.occ.origin))
            counts[f"fusion_{pl.occ.origin}_{piece_status(pl)}"] += 1
    tf.close()

    # --- contigs, edges (truth_contig, truth_edges)
    prophage_uids = set()
    for r in _tsv(world / "mobile" / "prophages.tsv"):
        lifted = placer.lift(r["record"], 6, int(r["start"]), int(r["end"]))
        if isinstance(lifted, str):
            continue
        for row, *_ in placer.frame_parts(r["record"], lifted[0], lifted[1], 1):
            prophage_uids.add(row["contig_uid"])
            for oc in placer.occs.get(row["contig_id"], ()):
                if oc.origin != "contaminant":
                    prophage_uids.add(oc.uid)
    tcn = _Out(tmp / "truth" / "truth_contig.tsv", TruthContig)
    te = _Out(tmp / "truth" / "truth_edges.tsv", TruthEdges)
    edges_done = set()
    contam_by_id = {c["contig_id"]: c for c in contam}
    for o in occ_rows:
        base = base_contig(o)
        rec = base["record"]
        uid = o["contig_uid"]
        L = int(o["length"])
        seq_gc = gc_acc[o["contig_accession"]]
        is_circ = base["is_circular"] == "true" and o["origin"] != "contaminant"
        tcn.row(TruthContig(*WV, uid, o["contig_accession"], o["source_label"], o["genome_id"], o["contig_id"],
                            o["origin"], L, replicon_id(rec), reps8[rec]["kind"], "circular" if is_circ else "linear",
                            is_circ, o["origin"] == "contaminant", False, uid in prophage_uids, seq_gc, o["file"],
                            o["record"]))
        if uid in edges_done:
            continue
        edges_done.add(uid)
        if o["origin"] != "contaminant":
            te.row(TruthEdges(*WV, uid, base["genome_id"], o["contig_id"], o["origin"], replicon_id(rec),
                              int(base["replicon_start"]), int(base["replicon_end"]), int(base["orientation"]),
                              int(base["rotation_offset"]), base["left_edge"], base["right_edge"],
                              base["repeat_family_left"], base["repeat_family_right"], n_partial[uid]))
            continue
        c = contam_by_id[o["contig_id"]]
        x0, y0, oc = int(c["donor_start"]), int(c["donor_end"]), int(c["orientation"])
        Lb, ob, a = int(base["length"]), int(base["orientation"]), int(base["replicon_start"])
        fa, fb = (a + x0, a + y0) if ob == 1 else (a + Lb - y0, a + Lb - x0)
        # edges at the window's ends, in donor-contig direction, then in the contaminant's
        lo = base["left_edge"] if x0 == 0 else "break_random"
        hi = base["right_edge"] if y0 == Lb else "break_random"
        flo = base["repeat_family_left"] if x0 == 0 else ""
        fhi = base["repeat_family_right"] if y0 == Lb else ""
        le, re_, fl, fr = (lo, hi, flo, fhi) if oc == 1 else (hi, lo, fhi, flo)
        te.row(TruthEdges(*WV, uid, base["genome_id"], o["contig_id"], "contaminant", replicon_id(rec), fa, fb,
                          ob * oc, 0, le, re_, fl, fr, n_partial[uid]))
    for r in placer.by_id.values():
        if r["emitted"] == "false":
            tcn.row(TruthContig(*WV, r["contig_uid"], "", "src_a", r["genome_id"], r["contig_id"], "assembly",
                                int(r["length"]), replicon_id(r["record"]), reps8[r["record"]]["kind"],
                                "circular" if r["is_circular"] == "true" else "linear", r["is_circular"] == "true",
                                False, True, False, float(r["gc"]), "", ""))
            counts["contigs_dropped"] += 1
    for ps in phage_src:
        p = ps["phage"]
        tcn.row(TruthContig(*WV, p["contig_uid"], ps["contig_accession"], "src_a", p["phage_id"], p["phage_id"],
                            "phage", int(p["length"]), p["phage_id"], "virus", "linear", False, False, False, False,
                            float(p["gc"]), f"fasta/{ps['accession']}.fna", ps["contig_accession"]))
        te.row(TruthEdges(*WV, p["contig_uid"], p["phage_id"], p["phage_id"], "phage", p["phage_id"], 0,
                          int(p["length"]), 1, 0, "replicon_end", "replicon_end", "", "", 0))
    counts["contigs_emitted"] = len(occ_rows) + len(phage_src)
    tcn.close()
    te.close()

    # --- duplicates and contaminants (truth_dup, truth_contam): S10's rows with the world columns
    acc_of = {(o["source_label"], o["contig_id"]): o["contig_accession"] for o in occ_rows}
    t = _Out(tmp / "truth" / "truth_dup.tsv", TruthDup)
    for d in dup:
        t.row(TruthDup(*WV, *(d[f] for f in sources.DupRow._fields)))
    t.close()
    t = _Out(tmp / "truth" / "truth_contam.tsv", TruthContam)
    for c in contam:
        t.row(TruthContam(*WV, *(c[f] for f in sources.ContamRow._fields), acc_of[("src_a", c["contig_id"])]))
    t.close()

    # --- families, modules, traps
    fam_bb = {r["family_id"]: r for r in _tsv(world / "backbone" / "families.tsv")}
    fam_id_of = {int(r["family"]): r["family_id"] for r in fam_bb.values()}
    fam_seq = {r["family_id"]: r for r in _tsv(world / "sequences" / "families.tsv")}
    transfers = collections.Counter()
    with open(world / "backbone" / "events.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"].startswith("transfer"):
                transfers[int(r["family"])] += 1
    panel_kind = "pfam" if docs["sequences"]["panel"]["source"] == "pfam" else "panel"
    tfm = _Out(tmp / "truth" / "truth_family.tsv", TruthFamily)
    special = ("vehicle_chassis", "cargo", "hub", "decoy")
    for fid in sorted(fam_seq):
        q, b = fam_seq[fid], fam_bb.get(fid)
        bkind = b["kind"] if b else ""
        dark = q["dark"] == "true"
        kind = bkind if bkind in special else "dark" if dark else panel_kind
        tfm.row(TruthFamily(*WV, fid, int(q["family"]), kind, bkind, q["domain_hmms"], not dark, q["root_source"],
                            float(q["rate_multiplier"]), float(q["target_identity"]),
                            float(q["realized_mean_identity"]), float(q["realized_min_identity"]), "",
                            fam_copies[fid], len(fam_genomes[fid]), transfers[int(q["family"])], ""))
    for r in _tsv(world / "fusions" / "nogd_families.tsv"):
        tfm.row(TruthFamily(*WV, r["family_id"], int(r["family"]), "nogd", "", r["domain_hmms"], r["dark"] != "true",
                            "hmm_sample", float(r["rate_multiplier"]), float(r["target_identity"]), None, None, "",
                            fam_copies[r["family_id"]], len(fam_genomes[r["family_id"]]), None, ""))
    cas_rows = _tsv(world / "mobile" / "cas.tsv")
    for r in cas_rows:
        tfm.row(TruthFamily(*WV, r["family_id"], int(r["family"]), "cas", "", r["domain_hmms"], r["dark"] != "true",
                            "hmm_sample", float(r["rate_multiplier"]), float(r["target_identity"]), None, None, "",
                            fam_copies[r["family_id"]], len(fam_genomes[r["family_id"]]), None, ""))
    for r in _tsv(ph_world / "sequences" / "families.tsv"):
        fid = "P" + r["family_id"][1:]
        tfm.row(TruthFamily(*WV, fid, int(r["family"]), "phage", "", r["domain_hmms"], r["dark"] != "true",
                            r["root_source"], float(r["rate_multiplier"]), float(r["target_identity"]),
                            float(r["realized_mean_identity"]), float(r["realized_min_identity"]), "",
                            fam_copies[fid], len(fam_genomes[fid]), None, ""))
    counts["families"] = tfm.n
    tfm.close()

    loci = collections.Counter()
    for u in _tsv(world / "backbone" / "truth_units.tsv"):
        if int(u["n_members"]) >= 2 and u["contiguous"] == "true":
            loci[u["group"]] += 1
    groups = _tsv(world / "backbone" / "groups.tsv")
    tm = _Out(tmp / "truth" / "truth_module.tsv", TruthModule)
    mkind = {"planted": "planted", "profile_only": "profile_only", "vehicle": "vehicle_chassis"}
    for gr in groups:
        if gr["kind"] not in mkind:
            continue
        for k, fam in enumerate(int(x) for x in gr["members"].split(";")):
            tm.row(TruthModule(*WV, gr["name"], mkind[gr["kind"]], fam_id_of[fam], "core", "", k,
                               float(gr["together"]), gr["fates"], gr["adjacent"] == "true", "", "",
                               loci[gr["name"]] if gr["adjacent"] == "true" else None, "", "", ""))
    for ev in _tsv(world / "fusions" / "events.tsv"):
        if ev["kind"] != "nogd":
            continue
        k = int(fam_bb[ev["family_a"]]["order_index"]) if ev["family_a"] in fam_bb else -1
        for fam, role in ((ev["family_a"], "nogd_replaced"), (ev["replacement_family"], "nogd_replacement")):
            tm.row(TruthModule(*WV, ev["module_id"], "planted", fam, role, ev["clade_node"], k, None, "", True, "",
                               "", None, "", "", ""))
    for r in cas_rows:
        tm.row(TruthModule(*WV, r["cas_module_id"], "cas_operon", r["family_id"], "core", "", int(r["order_index"]),
                           None, "", True, "", "", None, "", "", ""))
    tm.close()

    tt = _Out(tmp / "truth" / "truth_traps.tsv", TruthTrap)
    n_trap = 0

    def trap(kind, fam, veh, partners, label_, why):
        nonlocal n_trap
        tt.row(TruthTrap(*WV, f"TR{n_trap:04d}", kind, fam, veh, ";".join(partners), label_, why))
        n_trap += 1

    prom = collections.defaultdict(set)
    for r in _tsv(world / "sequences" / "domains.tsv"):
        if r["promiscuous"] == "true":
            prom[r["accession"]].add(r["family_id"])
    for acc in sorted(prom):
        fams = sorted(prom[acc])
        for f in fams:
            trap("promiscuous_domain", f, "", [x for x in fams if x != f], "must_not_call",
                 f"shares the promiscuous domain {acc} with unrelated families by construction")
    for gr in groups:
        if gr["kind"] != "vehicle":
            continue
        mem = [fam_id_of[int(x)] for x in gr["members"].split(";")]
        for f in mem[1:-1]:
            trap("mge_cargo", f, gr["name"], [x for x in mem if x != f], "must_not_call",
                 "cargo always shares chassis flanks; vehicle dedup must collapse")
    for fid, b in sorted(fam_bb.items()):
        if b["kind"] == "hub":
            trap("hub_family", fid, "", [], "must_not_call",
                 "transfers at many times the genome rate into many contexts; its neighbours are not a module")
        elif b["kind"] == "decoy":
            trap("decoy_profile", fid, "", [fam_id_of[int(b["partner"])]], "must_not_call",
                 "its transfers land only where the partner is present: profile without adjacency or co-transfer")
    node_gid = {h["species_node"]: g for g, h in host_genomes.items()}
    burst = sorted(node_gid[x] for x in bb["burst_tips"])
    trap("clonal_burst", "", "", [], "UNKNOWN",
         "genomes " + ",".join(burst) + " are one clonal burst under one BioProject: counts by genome must not "
         "treat them as independent")
    tt.close()

    # --- the stage document and the world manifest
    files = {str(f.relative_to(tmp)): sequences._sha(f) for f in sorted(tmp.rglob("*")) if f.is_file()}
    digest = _digest(files, 16)
    versions = {}
    for d, _ in UPSTREAM:
        versions.update(docs[d].get("versions", {}))
    doc = {"format": FORMAT, "seed": seed, "preset": cfg, "params_digest": backbone.params_digest(cfg),
           "world_id": world_id, "release_tag": release_tag,
           **{f"{d}_digest": docs[d]["digest"] for d, _ in UPSTREAM},
           "versions": dict(sorted(versions.items())),
           "stream_paths": ["S11/biome/<taxon>", "S11/mislabel/<source>/<genome>"],
           "counts": dict(sorted(counts.items())), "files": files, "digest": digest}
    (tmp / "emit.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    stage_digests = {d: docs[d]["digest"] for d, _ in UPSTREAM}
    stage_digests["emit"] = digest
    manifest = {"format": WORLD_FORMAT, "world_id": world_id, "release_tag": release_tag,
                "preset": preset["name"], "params_digest": backbone.params_digest(preset), "seed": seed,
                "versions": dict(sorted(versions.items())), "panel": docs["sequences"]["panel"],
                "motif_source": docs["replicons"]["motif_source"],
                "stages": {d: {"format": docs[d]["format"] if d in docs else FORMAT, "digest": v}
                           for d, v in sorted(stage_digests.items())},
                "world_digest": _digest(stage_digests, 16)}
    tmp.rename(out)
    mtmp = world / "manifest.json.partial"
    mtmp.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    mtmp.rename(world / "manifest.json")
    return read(world)


def read(world) -> EmitWorld:
    """Load ``world/emit/emit.json`` and ``world/manifest.json``."""
    world = Path(world)
    doc = json.loads((world / "emit" / "emit.json").read_text(encoding="utf-8"))
    if doc["format"] != FORMAT:
        raise ValueError(f"{world}/emit: format {doc['format']!r}, expected {FORMAT!r}")
    man = json.loads((world / "manifest.json").read_text(encoding="utf-8"))
    return EmitWorld(doc["seed"], doc["world_id"], doc["release_tag"], doc["digest"], man["world_digest"], world,
                     {**doc, "manifest": man})


__all__ = ["FORMAT", "WORLD_FORMAT", "TABLES", "check_config", "stop_bp", "assembly_quality", "window", "Placer",
           "Occ", "Place", "build", "read", "EmitWorld"]

"""synthworld gene fusions, splits and non-orthologous gene displacement (planning/tests/synthworld_spec.md §2.1
S7, §3.3 ``truth_fusion`` and the NOGD rows of ``truth_module``, §5.1; card F-02 session 6). Reads a world whose
S1-S6 are built and writes ``world/fusions/``: the host replicons of ``mobile/`` with the events applied, and the
lift from the S6 frame.

**Events.** ``fusions.fusions`` fusions ``FU00`` ..., ``fusions.splits`` splits ``SP00`` ... and ``fusions.nogd``
displacements ``NG00`` ..., placed in that order. Each happens once, on a host species-tree branch with
``clade_size`` extant genomes below it, at a time uniform along the branch, and is inherited by every extant
genome below (the clade), as S6 places its insertions (`mobile.clade`): stream ``S7/placement`` ranks
``TRIES`` branches per event and the first branch with a free candidate is used; stream ``S7/choose/<event>``
draws the candidate uniformly among the free ones, in sorted order. Candidates use host genes of ``mobile/``
only (no prophage or Cas gene), never a gene across a replicon's origin, and every family at most once over
the world's events, so no two events touch one gene.

- *Fusion* of families (A, B): in every clade genome A and B each have one copy (in the genome) and sit next to
  each other on one replicon, co-directional, A transcriptionally upstream of B (no gene between them; the
  junction across a replicon's origin is not used). The fused ORF is A's coding sequence without its stop
  codon, a linker of ``linker_aa`` codons (uniform on the range; residues uniform over ``linker_residues``,
  codons from the event lineage's usage, S4), and B's coding sequence without its initiation codon: on the +
  strand the bases [A.end - 3, B.start + 3) are replaced by the linker (on the - strand [B.end - 3, A.start +
  3) by its reverse complement), so the IGR between them and any overlap go, and the fused protein is
  protein(A) + linker + protein(B)[1:]. The linker is drawn once per event (stream ``S7/<event>``) and is the
  same in every clade genome.
- *Split* of a family C with two or more domains (the promiscuous trap counts): one copy per clade genome; a
  domain boundary k (1 <= k < domains, drawn once per event) and, per copy, the protein position p = floor of the
  midpoint between domain k - 1's end and domain k's start (``sequences/genes.tsv`` ``domain_arch``, the copy's own
  coordinates). The block STOP + spacer + ATG is inserted after codon p - 1 (on the copy's strand), so the first
  half encodes protein[:p] and the second M + protein[p:]. The stop codon is TAA; the spacer (``split_spacer``
  bases, uniform on the range) is the event lineage's intergenic composition (S5) with an RBS (AGGAGG) ending
  r bases before the ATG (r uniform on ``intergenic.cassettes.rbs_spacer``); both are drawn once per event.
- *NOGD* of a planted-module family M (``backbone/families.tsv`` kind planted; its copy in a unit of its
  group): one copy per clade genome, not overlapping its neighbours. A new family N = ``NOGD_FAMILY_BASE`` +
  10 e + j (the first j whose domains share no accession with M's) is designed and evolved by the S3-S4 code
  on the clade's tree (streams ``S3/family/<id>``, ``S4/family/<id>``), and each clade genome's copy of M has
  its coding sequence replaced by N's copy for that lineage, in place (strand, neighbours, IGRs kept).

**Lift.** ``edits.tsv`` is the lift from the S6 frame: each edit replaces the bases [start, end) of a replicon of
S6 (start = end for an insertion) by ``length`` new bases, at [new_start, new_end) of the new replicon. With the
edits of a replicon sorted, a feature [s, e) of S6 that no edit touches (an edit [a, b) touches it when a < e
and b > s, or a = b with s < a < e) moves to s + d and e + d, d the sum of length - (b - a) over the edits
with b <= s; one across the origin (e > N) ends at N' + (e - N), and no edit touches either of its parts.
Every S6 feature that an edit touches belongs to that event (its genes, or the IGR and elements between a fused
pair).

Shape on disk (``build(world, pfam_hmm=None)``)::

    fasta/<genome_id>.fna, .fai   host replicons after the events (records as in S5 and S6)
    replicons.tsv   genome_id replicon record kind topology length s6_length n_genes n_edits gc contig_uid
    edits.tsv       edit_id event_id kind genome_id replicon record start end length new_start new_end
    genes.tsv       genome_id replicon copy_label family_id source insertion_id event_id parents start end strand
                    origin_spanning cds_edit protein_uid domain_arch
    events.tsv      event_id kind clade_node insertion_time family_a family_b split_domain_index module_id
                    replacement_family linker_aa linker spacer n_genomes genomes
    fusions.tsv     fusion_id kind family_a family_b split_domain_index clade_node genome_id record contig_uid
                    start end strand junction_aa linker_aa gene_ids expected_arch in_module
    nogd.tsv        event_id module_id family_replaced family_replacement clade_node genome_id record replaced_copy
                    copy_label start end strand protein_uid
    nogd_families.tsv event_id module_id family family_id dark n_domains domain_hmms target_identity
                    rate_multiplier
    fusions.json    format, seed, the ``fusions`` preset and its digest, upstream digests, versions, stream
                    paths, counts, sha256 of every file, digest

Coordinates are 0-based half-open on the + strand of their record; strand +1/-1; a feature across the origin has
end > length. ``genes.tsv`` lists every gene of the new replicons, in (genome, replicon, start) order: S6's
genes lifted (``event_id`` empty), and the genes the events made (``source`` fused, split or nogd, ``parents``
the S6 copy labels they come from). A fused gene is ``<A copy>+<B copy>`` of family ``<A>+<B>``; split halves
are ``<C copy>/1`` and ``/2`` of family C; a NOGD gene is ``<lineage>_<event>`` of family N. ``protein_uid`` is
BLAKE3-128 of the protein (initiator M, no stop) of every gene; ``domain_arch`` (``accession:start-end;...`` on
the gene's protein) is given for the genes the events made (S6's are in ``sequences/genes.tsv``).
``fusions.tsv`` is §3.3 ``truth_fusion`` on the S7 records: for a fusion its fused ORF (``gene_ids`` one gene,
``junction_aa`` the protein position where A ends, ``expected_arch`` fused); for a split the span of the two
halves (``family_a`` = ``family_b`` = C, ``junction_aa`` p, ``expected_arch`` split_adjacent); ``in_module`` the
planted group of A and B when both copies are in one unit of it. ``nogd.tsv`` gives the NOGD instances and
``events.tsv`` the module rows (``member_role`` nogd_replaced for ``family_a``, nogd_replacement for
``replacement_family``, ``nogd_clade_node`` the clade node). Genomes without an edit are copied byte for byte.
``digest`` is BLAKE3-128 over the sorted (path, sha256) list of every file but fusions.json.
"""

from __future__ import annotations

import collections
import json
import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
from blake3 import blake3

from biobox import ids
from biobox.synthworld import backbone, codons, intergenic, mobile, profiles, seeds, sequences

FORMAT = "synthworld-fusions/1"
NOGD_FAMILY_BASE = 910_000
TRIES = 8
SPLIT_STOP = "TAA"
RBS = intergenic.RBS
q6 = sequences.q6
revcomp = intergenic.revcomp
_ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)


# --- the lift ----------------------------------------------------------------------------------------------

class Edits:
    """Replacements ([a, b) -> l new bases, a = b an insertion) on a record of length N, non-overlapping."""

    def __init__(self, N: int, edits):
        self.N = int(N)
        ed = sorted((int(a), int(b), int(l)) for a, b, l in edits)
        for a, b, l in ed:
            if not 0 <= a <= b <= self.N or l < 0:
                raise ValueError(f"edit [{a}, {b}) -> {l} does not fit a record of {self.N}")
            if a == b == 0:
                raise ValueError("an insertion at 0 is ambiguous on a circular record")
            if a == b and l == 0:
                raise ValueError("an empty insertion is not an edit")
        for (a0, b0, _), (a1, b1, _) in zip(ed, ed[1:]):
            if b0 > a1 or (a0 == b0 == a1 == b1):
                raise ValueError(f"edits [{a0}, {b0}) and [{a1}, {b1}) overlap")
        self.a = np.array([x[0] for x in ed], dtype=np.int64)
        self.b = np.array([x[1] for x in ed], dtype=np.int64)
        self.l = np.array([x[2] for x in ed], dtype=np.int64)
        self.cum = np.concatenate([[0], np.cumsum(self.l - (self.b - self.a))])
        self.new_length = self.N + int(self.cum[-1])

    def new_starts(self) -> list[int]:
        """Where each edit's new bases start, in sorted edit order."""
        return [int(a) + int(self.cum[k]) for k, a in enumerate(self.a)]

    def shift(self, x: int) -> int:
        """The length change of the edits with b <= x (the edits wholly before position x)."""
        return int(self.cum[np.searchsorted(self.b, int(x), side="right")])

    def touched(self, s: int, e: int) -> bool:
        s, e = int(s), int(e)
        parts = [(s, e)] if e <= self.N else [(s, self.N + 1), (-1, e - self.N)]
        for x, y in parts:
            for a, b in zip(self.a.tolist(), self.b.tolist()):
                if (a < b and a < y and b > x) or (a == b and x < a < y):
                    return True
        return False

    def feature(self, s: int, e: int) -> tuple[int, int]:
        """[s, e) of the old record in the new one; refuses a feature an edit touches."""
        s, e = int(s), int(e)
        if self.touched(s, e):
            raise ValueError(f"feature [{s}, {e}) is touched by an edit")
        d = int(self.cum[np.searchsorted(self.b, s, side="right")])
        if e <= self.N:
            return s + d, e + d
        return s + d, self.new_length + (e - self.N)


def apply_edits(seq: str, edits) -> tuple[str, list[int]]:
    """``seq`` with the edits (a, b, block) applied; returns it and each block's new start, in the order given."""
    order = sorted(range(len(edits)), key=lambda k: (edits[k][0], edits[k][1]))
    parts, starts, prev, shift = [], [0] * len(edits), 0, 0
    for k in order:
        a, b, block = int(edits[k][0]), int(edits[k][1]), edits[k][2]
        if a < prev:
            raise ValueError("edits overlap")
        parts.append(seq[prev:a])
        starts[k] = a + shift
        parts.append(block)
        shift += len(block) - (b - a)
        prev = b
    parts.append(seq[prev:])
    return "".join(parts), starts


# --- the edits of one event in one genome ------------------------------------------------------------------

def fusion_edit(a: tuple[int, int, int], b: tuple[int, int, int], linker: str) -> tuple[int, int, str]:
    """The replacement that fuses gene a = (start, end, strand) with the transcriptionally downstream gene b."""
    (sa, ea, st), (sb, eb, stb) = a, b
    if st != stb:
        raise ValueError("a fusion joins co-directional genes")
    if st == 1:
        if not (sa < sb and ea - 3 <= sb + 3):
            raise ValueError("on +, gene a must end where gene b starts (an overlap of at most 4)")
        return ea - 3, sb + 3, linker
    if not (sb < sa and eb - 3 <= sa + 3):
        raise ValueError("on -, gene b must end where gene a starts (an overlap of at most 4)")
    return eb - 3, sa + 3, revcomp(linker)


def split_edit(g: tuple[int, int, int], p: int, block: str) -> tuple[int, int, str]:
    """The insertion that splits gene g after codon p - 1 (``block`` read on the gene's strand)."""
    s, e, st = g
    if not 1 <= p < (e - s) // 3 - 1:
        raise ValueError(f"split position {p} is not inside the gene")
    at = s + 3 * p if st == 1 else e - 3 * p
    return at, at, block if st == 1 else revcomp(block)


def split_point(arch: list[tuple[str, int, int]], k: int) -> int:
    """Protein position between domain k - 1 and domain k (floor of the midpoint of the gap)."""
    return (arch[k - 1][2] + arch[k][1]) // 2


def parse_arch(text: str) -> list[tuple[str, int, int]]:
    out = []
    for item in text.split(";") if text else []:
        acc, span = item.rsplit(":", 1)
        a, b = span.split("-")
        out.append((acc, int(a), int(b)))
    return out


def format_arch(arch) -> str:
    return ";".join(f"{acc}:{a}-{b}" for acc, a, b in arch)


# --- tables ------------------------------------------------------------------------------------------------

class RepliconRow(NamedTuple):
    genome_id: str
    replicon: int
    record: str
    kind: str
    topology: str
    length: int
    s6_length: int
    n_genes: int
    n_edits: int
    gc: float
    contig_uid: str


class EditRow(NamedTuple):
    edit_id: str
    event_id: str
    kind: str
    genome_id: str
    replicon: int
    record: str
    start: int
    end: int
    length: int
    new_start: int
    new_end: int


class GeneRow(NamedTuple):
    genome_id: str
    replicon: int
    copy_label: str
    family_id: str
    source: str
    insertion_id: str
    event_id: str
    parents: tuple[str, ...]
    start: int
    end: int
    strand: int
    origin_spanning: bool
    cds_edit: str
    protein_uid: str
    domain_arch: str


class EventRow(NamedTuple):
    event_id: str
    kind: str
    clade_node: str
    insertion_time: float
    family_a: str
    family_b: str
    split_domain_index: int
    module_id: str
    replacement_family: str
    linker_aa: int
    linker: str
    spacer: str
    n_genomes: int
    genomes: tuple[str, ...]


class FusionRow(NamedTuple):
    fusion_id: str
    kind: str
    family_a: str
    family_b: str
    split_domain_index: int
    clade_node: str
    genome_id: str
    record: str
    contig_uid: str
    start: int
    end: int
    strand: int
    junction_aa: int
    linker_aa: int
    gene_ids: tuple[str, ...]
    expected_arch: str
    in_module: str


class NogdRow(NamedTuple):
    event_id: str
    module_id: str
    family_replaced: str
    family_replacement: str
    clade_node: str
    genome_id: str
    record: str
    replaced_copy: str
    copy_label: str
    start: int
    end: int
    strand: int
    protein_uid: str


class NogdFamilyRow(NamedTuple):
    event_id: str
    module_id: str
    family: int
    family_id: str
    dark: bool
    n_domains: int
    domain_hmms: tuple[str, ...]
    target_identity: float
    rate_multiplier: float


_TABLES = {"replicons.tsv": RepliconRow, "edits.tsv": EditRow, "genes.tsv": GeneRow, "events.tsv": EventRow,
           "fusions.tsv": FusionRow, "nogd.tsv": NogdRow, "nogd_families.tsv": NogdFamilyRow}


def check_config(cfg: dict) -> None:
    for k in ("fusions", "splits", "nogd"):
        if int(cfg[k]) < 0:
            raise ValueError(f"fusions.{k} must be >= 0")
    if int(cfg["nogd"]) > 99:
        raise ValueError("at most 99 NOGD events (family numbers 910000 + 10 e + j)")
    for k in ("clade_size", "linker_aa", "split_spacer"):
        lo, hi = (int(x) for x in cfg[k])
        if not 0 <= lo <= hi:
            raise ValueError(f"fusions.{k} must be 0 <= lo <= hi, got {cfg[k]}")
    if int(cfg["clade_size"][0]) < 1:
        raise ValueError("fusions.clade_size starts at 1")
    if int(cfg["split_spacer"][0]) < len(RBS) + 10:
        raise ValueError(f"fusions.split_spacer must hold the RBS and its spacer (>= {len(RBS) + 10})")
    bad = set(str(cfg["linker_residues"])) - set(codons.AA)
    if bad or not cfg["linker_residues"]:
        raise ValueError(f"fusions.linker_residues are amino acids of {codons.AA}, got {cfg['linker_residues']!r}")


@dataclass
class Event:
    event_id: str
    kind: str                  # fusion | split | nogd
    tries: list[tuple[str, float]]
    node: str = ""
    t_insert: float = 0.0
    clade: mobile.Clade | None = None
    key: tuple = ()
    data: dict = field(default_factory=dict)


@dataclass
class FusionWorld:
    seed: int
    events: list[EventRow]
    replicons: list[RepliconRow]
    digest: str
    out: Path
    meta: dict = field(default_factory=dict, repr=False)


def _candidates(rows: list[dict], fam_kind: dict[int, str], n_arch: dict[int, int], group_of: dict[str, tuple[str, int]],
                fam_group: dict[int, str]) -> dict[str, set]:
    """One genome's candidate keys per event kind (module docstring)."""
    n = collections.Counter(r["family_id"] for r in rows)
    by_rep = collections.defaultdict(list)
    for r in rows:
        by_rep[int(r["replicon"])].append(r)
    fusion, split, nogd = set(), set(), set()

    def ok(r) -> bool:
        return r["source"] == "host" and r["origin_spanning"] == "false" and n[r["family_id"]] == 1

    for rep, rr in by_rep.items():
        for i, r in enumerate(rr):
            if not ok(r):
                continue
            f = int(r["family_id"][1:])
            if n_arch.get(f, 0) >= 2:
                split.add(f)
            if fam_kind.get(f) == "planted" and 0 < i < len(rr) - 1:
                g = group_of.get(r["copy_label"])
                if g and g[0] == fam_group[f] and int(rr[i - 1]["end"]) <= int(r["start"]) \
                        and int(r["end"]) <= int(rr[i + 1]["start"]):
                    nogd.add(f)
            if i + 1 < len(rr):
                q = rr[i + 1]
                if ok(q) and r["strand"] == q["strand"]:
                    a, b = (r, q) if r["strand"] == "1" else (q, r)
                    fusion.add((int(a["family_id"][1:]), int(b["family_id"][1:])))
    return {"fusion": fusion, "split": split, "nogd": nogd}


def _draw_codons(rng, residues: str, n: int, cum_row: np.ndarray) -> tuple[str, str]:
    """n codons for residues drawn uniformly from ``residues``, from one lineage's usage (``cum_row``)."""
    if n == 0:
        return "", ""
    aa = "".join(residues[int(i)] for i in rng.integers(len(residues), size=n))
    rows = np.array([codons.AA.index(x) for x in aa], dtype=np.int64)
    idx = codons.draw(cum_row, rows, rng.random(n))
    return codons.CODON_BYTES[idx].tobytes().decode("ascii"), aa


def _background(rng, n: int, comp_cum: np.ndarray) -> str:
    return _ACGT[np.minimum(np.count_nonzero(comp_cum[None, :] <= rng.random(n)[:, None], axis=1), 3)] \
        .tobytes().decode("ascii")


def build(world, *, pfam_hmm=None) -> FusionWorld:
    """Write ``world/fusions/`` (module docstring) from a world whose S1-S6 are built."""
    world = Path(world)
    bb = json.loads((world / "backbone" / "backbone.json").read_text(encoding="utf-8"))
    sq = json.loads((world / "sequences" / "sequences.json").read_text(encoding="utf-8"))
    mb = json.loads((world / "mobile" / "mobile.json").read_text(encoding="utf-8"))
    for doc, want in ((bb, backbone.FORMAT), (sq, sequences.FORMAT), (mb, mobile.FORMAT)):
        if doc["format"] != want:
            raise ValueError(f"{world}: format {doc['format']!r}, expected {want!r}")
    preset = bb["preset"]
    cfg = preset["fusions"]
    check_config(cfg)
    seed = int(bb["seed"])
    out = world / "fusions"
    if out.exists():
        raise FileExistsError(f"{out} exists; the fusions stage is written once")
    tmp = world / "fusions.partial"
    if tmp.exists():
        shutil.rmtree(tmp)
    (tmp / "fasta").mkdir(parents=True)

    branches = backbone.species_branches(world)
    host_genomes = mobile._tsv(world / "backbone" / "genomes.tsv")
    gid_of = {r["species_node"]: r["genome_id"] for r in host_genomes}
    node_of = {g: n for n, g in gid_of.items()}
    tree = mobile.species_tree(branches, gid_of)
    seq_cfg = sq["sequences_preset"]
    lin_seq = sequences.lineage_model(branches, seq_cfg, seed)
    lin_igr = intergenic.lineages(world, preset["intergenic"]["composition"])
    fams = mobile._tsv(world / "backbone" / "families.tsv")
    fam_kind = {int(r["family"]): r["kind"] for r in fams}
    fam_group = {int(r["family"]): r["group"] for r in fams}
    groups = {r["name"]: r for r in mobile._tsv(world / "backbone" / "groups.tsv")}
    planted = {g for g, r in groups.items() if r["kind"] == "planted"}
    group_of = {}
    for r in mobile._tsv(world / "backbone" / "genes.tsv"):
        if r["group"] in planted:
            group_of[r["copy_label"]] = (r["group"], int(r["unit"]))
    seq_genes = {}
    n_arch: dict[int, int] = {}
    for r in mobile._tsv(world / "sequences" / "genes.tsv"):
        seq_genes[r["copy_label"]] = (r["protein_uid"], r["domain_arch"])
        n_arch[int(r["family"])] = len(parse_arch(r["domain_arch"]))

    # --- placement (stream S7/placement)
    rng = seeds.generator(seed, "S7/placement")
    lo, hi = (int(x) for x in cfg["clade_size"])
    cand_nodes = [x for x in tree.labels[1:] if lo <= len(tree.below[x]) <= hi]
    if not cand_nodes:
        raise ValueError(f"no species branch has {lo} to {hi} extant genomes below it")
    events: list[Event] = []
    for kind, prefix, n_ev in (("fusion", "FU", int(cfg["fusions"])), ("split", "SP", int(cfg["splits"])),
                               ("nogd", "NG", int(cfg["nogd"]))):
        for i in range(n_ev):
            perm = rng.permutation(len(cand_nodes))[:TRIES]
            events.append(Event(f"{prefix}{i:02d}", kind, [(cand_nodes[int(k)], float(rng.random())) for k in perm]))

    # --- pass 1: per (event, try), the candidate keys in every clade genome
    gen_tries: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for ei, ev in enumerate(events):
        for ti, (node, _) in enumerate(ev.tries):
            for tip in tree.below[node]:
                gen_tries[gid_of[tip]].append((ei, ti))
    stats: dict[tuple[int, int], collections.Counter] = collections.defaultdict(collections.Counter)
    for gid, rows in mobile._by_genome(world / "mobile" / "genes.tsv"):
        if gid not in gen_tries:
            continue
        cands = _candidates(rows, fam_kind, n_arch, group_of, fam_group)
        for ei, ti in gen_tries[gid]:
            stats[(ei, ti)].update(cands[events[ei].kind])
    used: set[int] = set()
    for ei, ev in enumerate(events):
        for ti, (node, frac) in enumerate(ev.tries):
            k = len(tree.below[node])
            keys = sorted(x for x, c in stats[(ei, ti)].items() if c == k
                          and not (set(x if isinstance(x, tuple) else (x,)) & used))
            if keys:
                ev.key = keys[int(seeds.generator(seed, f"S7/choose/{ev.event_id}").integers(len(keys)))]
                ev.node = node
                t = q6(tree.start[node] + frac * (tree.end[node] - tree.start[node]))
                ev.t_insert = min(max(t, tree.start[node]), tree.end[node])
                ev.clade = mobile.clade(tree, node, ev.t_insert)
                used |= set(ev.key if isinstance(ev.key, tuple) else (ev.key,))
                break
        if ev.clade is None:
            raise ValueError(f"{ev.event_id}: none of {len(ev.tries)} branches has a free candidate")

    # --- per event content
    panel = profiles.panel(seq_cfg["profiles"], seed, pfam_hmm)
    rbs_lo, rbs_hi = (int(x) for x in preset["intergenic"]["cassettes"]["rbs_spacer"])
    nogd_rows: list[NogdFamilyRow] = []
    for ev in events:
        r = seeds.generator(seed, f"S7/{ev.event_id}")
        li_seq, li_igr = lin_seq.index[ev.clade.node], lin_igr.index[ev.clade.node]
        if ev.kind == "fusion":
            n = int(r.integers(int(cfg["linker_aa"][0]), int(cfg["linker_aa"][1]) + 1))
            ev.data["linker"], ev.data["linker_prot"] = _draw_codons(r, str(cfg["linker_residues"]), n,
                                                                     lin_seq.cum[li_seq])
        elif ev.kind == "split":
            ev.data["k"] = 1 + int(r.integers(n_arch[ev.key] - 1))
            L = int(r.integers(int(cfg["split_spacer"][0]), int(cfg["split_spacer"][1]) + 1))
            rr = int(r.integers(rbs_lo, rbs_hi + 1))
            comp = lin_igr.comp_cum[li_igr]
            sp = _background(r, L - len(RBS) - rr, comp) + RBS + _background(r, rr, comp)
            ev.data["spacer"] = sp
            ev.data["block"] = SPLIT_STOP + sp + "ATG"
        else:
            e = int(ev.event_id[2:])
            m_accs = _family_accessions(world, ev.key)
            for j in range(10):
                fam = NOGD_FAMILY_BASE + 10 * e + j
                d = sequences.design_family(seq_cfg, seed, fam, panel)
                if not ({x.accession for x in d.domains} & m_accs):
                    break
            else:
                raise ValueError(f"{ev.event_id}: no replacement family without the domains of F{ev.key:06d}")
            nodes = mobile.gene_tree(ev.clade, fam, ev.event_id)
            res = sequences.evolve_family(d, nodes, lin_seq, seq_cfg, seeds.generator(seed, f"S4/family/{d.family_id}"))
            ev.data["family"] = d.family_id
            ev.data["leaves"] = {lf.lineage: lf for lf in res.leaves}
            ev.data["accessions"] = [x.accession for x in d.domains]
            nogd_rows.append(NogdFamilyRow(ev.event_id, fam_group[ev.key], fam, d.family_id, d.dark, d.n_domains,
                                           tuple(x.accession for x in d.domains), d.target_identity, res.mu))

    # --- pass 2: one host genome at a time
    by_genome: dict[str, list[Event]] = collections.defaultdict(list)
    for ev in events:
        for tip in ev.clade.tips:
            by_genome[gid_of[tip]].append(ev)
    reps6 = collections.defaultdict(list)
    for r in mobile._tsv(world / "mobile" / "replicons.tsv"):
        reps6[r["genome_id"]].append(r)
    w = {n: mobile._Rows(tmp / n, c) for n, c in _TABLES.items() if n in ("replicons.tsv", "edits.tsv", "genes.tsv")}
    fusion_rows: list[FusionRow] = []
    nogd_inst: list[NogdRow] = []
    rep_rows: list[RepliconRow] = []
    counts = collections.Counter()
    for gid, rows in mobile._by_genome(world / "mobile" / "genes.tsv"):
        my = by_genome.get(gid, [])
        node = node_of[gid]
        src = world / "mobile" / "fasta" / f"{gid}.fna"
        by_rep = collections.defaultdict(list)
        for r in rows:
            by_rep[int(r["replicon"])].append(r)
        if not my:
            shutil.copyfile(src, tmp / "fasta" / f"{gid}.fna")
            shutil.copyfile(str(src) + ".fai", tmp / "fasta" / f"{gid}.fna.fai")
        recs = dict(intergenic._fasta_records(src)) if my else {}
        by_fam = {r["family_id"]: r for r in rows if r["source"] == "host"}
        fa = sequences._Fasta(tmp / "fasta" / f"{gid}.fna") if my else None
        for rr in sorted(reps6[gid], key=lambda r: int(r["replicon"])):
            rep, rec = int(rr["replicon"]), rr["record"]
            # the edits of this replicon: (a, b, block, event, payload)
            plan = []
            if my:
                seq = recs[rec]
                for ev in my:
                    _plan_event(ev, gid, node, rep, seq, by_fam, seq_genes, plan)
            if not plan:
                out_genes = [_lifted(gid, r, None, seq_genes) for r in by_rep[rep]]
                if fa is not None:
                    fa.write(rec, recs[rec])
                rrow = RepliconRow(gid, rep, rec, rr["kind"], rr["topology"], int(rr["length"]), int(rr["length"]),
                                   len(out_genes), 0, float(rr["gc"]), rr["contig_uid"])
            else:
                N = len(seq)
                new, starts = apply_edits(seq, [(a, b, blk) for a, b, blk, _, _ in plan])
                lift = Edits(N, [(a, b, len(blk)) for a, b, blk, _, _ in plan])
                if len(new) != lift.new_length:
                    raise AssertionError(f"{rec}: edited length disagrees with the lift")
                fa.write(rec, new)
                cuid = ids.contig_uid(new.encode("ascii")).hex()
                touched = set()
                event_genes = []
                for k, ((a, b, blk, ev, pay), b0) in enumerate(zip(plan, starts)):
                    eid = f"{ev.event_id}@{gid}"
                    w["edits.tsv"].row(EditRow(eid, ev.event_id, ev.kind, gid, rep, rec, a, b, len(blk), b0,
                                               b0 + len(blk)))
                    touched |= set(pay["copies"])
                    event_genes += _event_genes(ev, gid, rep, rec, cuid, new, a, b, blk, b0, pay, lift,
                                                fusion_rows, nogd_inst, group_of, seq_genes)
                    counts[ev.kind + "_instances"] += 1
                out_genes = [_lifted(gid, r, lift, seq_genes) for r in by_rep[rep] if r["copy_label"] not in touched]
                for r in by_rep[rep]:
                    if r["copy_label"] not in touched and lift.touched(int(r["start"]), int(r["end"])):
                        raise AssertionError(f"{r['copy_label']} is touched by an edit of another gene")
                out_genes += event_genes
                rrow = RepliconRow(gid, rep, rec, rr["kind"], rr["topology"], len(new), N, len(out_genes), len(plan),
                                   q6((new.count("G") + new.count("C")) / len(new)), cuid)
            out_genes.sort(key=lambda x: (x.start, x.end, x.copy_label))
            for x in out_genes:
                w["genes.tsv"].row(x)
            w["replicons.tsv"].row(rrow)
            rep_rows.append(rrow)
            counts["genes"] += len(out_genes)
            counts["replicons"] += 1
        if fa is not None:
            fa.close()
        counts["genomes"] += 1
        counts["genomes_edited"] += bool(my)
    for x in w.values():
        x.close()

    # --- event truth
    ev_rows = []
    for ev in events:
        fa_, fb_, k, mod, repl = "", "", -1, "", ""
        if ev.kind == "fusion":
            fa_, fb_ = f"F{ev.key[0]:06d}", f"F{ev.key[1]:06d}"
        elif ev.kind == "split":
            fa_ = fb_ = f"F{ev.key:06d}"
            k = ev.data["k"]
        else:
            fa_ = f"F{ev.key:06d}"
            mod, repl = fam_group[ev.key], ev.data["family"]
        ev_rows.append(EventRow(ev.event_id, ev.kind, ev.node, ev.t_insert, fa_, fb_, k, mod, repl,
                                len(ev.data.get("linker", "")) // 3, ev.data.get("linker", ""),
                                ev.data.get("spacer", ""), len(ev.clade.tips),
                                tuple(sorted(gid_of[x] for x in ev.clade.tips))))
    for name, rows_ in (("events.tsv", ev_rows), ("nogd_families.tsv", nogd_rows),
                        ("fusions.tsv", sorted(fusion_rows, key=lambda x: (x.fusion_id, x.genome_id))),
                        ("nogd.tsv", sorted(nogd_inst, key=lambda x: (x.event_id, x.genome_id)))):
        t = mobile._Rows(tmp / name, _TABLES[name])
        for x in rows_:
            t.row(x)
        t.close()
    counts.update(events=len(events))

    files = {str(f.relative_to(tmp)): sequences._sha(f) for f in sorted(tmp.rglob("*")) if f.is_file()}
    digest = blake3("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest(length=16)
    doc = {"format": FORMAT, "seed": seed, "preset": cfg, "params_digest": backbone.params_digest(cfg),
           "backbone_digest": bb["digest"], "sequences_digest": sq["digest"], "mobile_digest": mb["digest"],
           "versions": {"numpy": np.__version__},
           "stream_paths": ["S7/placement", "S7/choose/<event>", "S7/<event>", "S3/family/F91<ee><j>",
                            "S4/family/F91<ee><j>"],
           "counts": dict(sorted(counts.items())), "files": files, "digest": digest}
    (tmp / "fusions.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    tmp.rename(out)
    return read(world)


def _family_accessions(world: Path, family: int) -> set[str]:
    for r in mobile._tsv(world / "sequences" / "families.tsv"):
        if int(r["family"]) == family:
            return {x for x in r["domain_hmms"].split(";") if x}
    raise KeyError(family)


def _plan_event(ev: Event, gid: str, node: str, rep: int, seq: str, by_fam: dict, seq_genes: dict, plan: list) -> None:
    """Append the edit of ``ev`` in genome ``gid`` if it falls on replicon ``rep``."""
    def g(fid):
        r = by_fam[fid]
        return r, (int(r["start"]), int(r["end"]), int(r["strand"]))

    if ev.kind == "fusion":
        ra, a = g(f"F{ev.key[0]:06d}")
        rb, b = g(f"F{ev.key[1]:06d}")
        if int(ra["replicon"]) != rep:
            return
        if int(rb["replicon"]) != rep:
            raise AssertionError(f"{ev.event_id}@{gid}: the fused pair is on two replicons")
        x, y, blk = fusion_edit(a, b, ev.data["linker"])
        if a[2] == 1:
            stop, start = seq[a[1] - 3:a[1]], seq[b[0]:b[0] + 3]
        else:
            stop, start = revcomp(seq[a[0]:a[0] + 3]), revcomp(seq[b[1] - 3:b[1]])
        if stop not in intergenic.STOPS or start not in ("ATG", "GTG", "TTG"):
            raise AssertionError(f"{ev.event_id}@{gid}: the fused pair does not end and start as genes ({stop}, {start})")
        plan.append((x, y, blk, ev, {"copies": (ra["copy_label"], rb["copy_label"]), "a": a, "b": b, "ra": ra, "rb": rb}))
    elif ev.kind == "split":
        rc_, c = g(f"F{ev.key:06d}")
        if int(rc_["replicon"]) != rep:
            return
        arch = parse_arch(seq_genes[rc_["copy_label"]][1])
        p = split_point(arch, ev.data["k"])
        x, y, blk = split_edit(c, p, ev.data["block"])
        plan.append((x, y, blk, ev, {"copies": (rc_["copy_label"],), "c": c, "p": p, "arch": arch, "rc": rc_}))
    else:
        rm, m = g(f"F{ev.key:06d}")
        if int(rm["replicon"]) != rep:
            return
        lf = ev.data["leaves"][node]
        blk = lf.cds if m[2] == 1 else revcomp(lf.cds)
        plan.append((m[0], m[1], blk, ev, {"copies": (rm["copy_label"],), "m": m, "rm": rm, "leaf": lf}))


def _lifted(gid: str, r: dict, lift: Edits | None, seq_genes: dict) -> GeneRow:
    s, e = int(r["start"]), int(r["end"])
    if lift is not None:
        s, e = lift.feature(s, e)
    uid = r["protein_uid"] or seq_genes[r["copy_label"]][0]
    N = lift.new_length if lift is not None else None
    span = e > N if N is not None else r["origin_spanning"] == "true"
    return GeneRow(gid, int(r["replicon"]), r["copy_label"], r["family_id"], r["source"], r["insertion_id"], "", (),
                   s, e, int(r["strand"]), span, r["cds_edit"], uid, "")


def _protein(seq: str, s: int, e: int, strand: int) -> str:
    cds = seq[s:e] if strand == 1 else revcomp(seq[s:e])
    p = codons.translate(cds)
    if not p.endswith("*") or "*" in p[:-1]:
        raise AssertionError("an event gene does not read as one ORF")
    return p[:-1]


def _uid(prot: str) -> str:
    return blake3(prot.encode("ascii")).hexdigest(length=16)


def _event_genes(ev: Event, gid: str, rep: int, rec: str, cuid: str, new: str, a: int, b: int, blk: str, b0: int,
                 pay: dict, lift: Edits, fusion_rows: list, nogd_inst: list, group_of: dict,
                 seq_genes: dict) -> list[GeneRow]:
    d = len(blk) - (b - a)
    if ev.kind == "fusion":
        (sa, ea, st), (sb, eb, _) = pay["a"], pay["b"]
        ra, rb = pay["ra"], pay["rb"]
        lo, hi = (min(sa, sb), max(ea, eb))
        s2 = lo + lift.shift(lo)
        e2 = s2 + (hi - lo) + d
        prot = _protein(new, s2, e2, st)
        la = (ea - sa) // 3 - 1                                     # protein(A) length
        n_link = len(ev.data["linker"]) // 3
        label = f"{ra['copy_label']}+{rb['copy_label']}"
        fam = f"{ra['family_id']}+{rb['family_id']}"
        ga, gb = group_of.get(ra["copy_label"]), group_of.get(rb["copy_label"])
        in_mod = ga[0] if ga and ga == gb else ""
        off = la + n_link - 1                                       # protein(B)[1:] follows the linker
        arch = format_arch(parse_arch(seq_genes[ra["copy_label"]][1])
                           + [(acc, x + off, y + off) for acc, x, y in parse_arch(seq_genes[rb["copy_label"]][1])])
        edit = ";".join([x for x in ra["cds_edit"].split(";") if x.startswith("start:")]
                        + [x for x in rb["cds_edit"].split(";") if x.startswith("stop:")])
        fusion_rows.append(FusionRow(ev.event_id, "fusion", ra["family_id"], rb["family_id"], -1, ev.node, gid, rec,
                                     cuid, s2, e2, st, la, n_link, (label,), "fused", in_mod))
        return [GeneRow(gid, rep, label, fam, "fused", ra["insertion_id"], ev.event_id,
                        (ra["copy_label"], rb["copy_label"]), s2, e2, st, False, edit, _uid(prot), arch)]
    if ev.kind == "split":
        s, e, st = pay["c"]
        p, arch, rc_ = pay["p"], pay["arch"], pay["rc"]
        s2 = s + lift.shift(s)
        e2 = s2 + (e - s) + d
        k = ev.data["k"]
        L1 = 3 * p + 3                                             # first half with its stop codon
        L2 = 3 + (e - s) - 3 * p                                   # ATG + the rest
        if st == 1:
            h1, h2 = (s2, s2 + L1), (e2 - L2, e2)
        else:
            h1, h2 = (e2 - L1, e2), (s2, s2 + L2)
        p1, p2 = _protein(new, *h1, st), _protein(new, *h2, st)
        a1 = format_arch(arch[:k])
        a2 = format_arch([(acc, x - p + 1, y - p + 1) for acc, x, y in arch[k:]])
        lab = rc_["copy_label"]
        fusion_rows.append(FusionRow(ev.event_id, "split", rc_["family_id"], rc_["family_id"], k, ev.node, gid, rec,
                                     cuid, s2, e2, st, p, 0, (f"{lab}/1", f"{lab}/2"), "split_adjacent", ""))
        edit = rc_["cds_edit"]
        e1 = ";".join(x for x in edit.split(";") if x.startswith("start:"))
        e2_ = ";".join(x for x in edit.split(";") if x.startswith("stop:"))
        return [GeneRow(gid, rep, f"{lab}/1", rc_["family_id"], "split", rc_["insertion_id"], ev.event_id, (lab,),
                        *h1, st, False, e1, _uid(p1), a1),
                GeneRow(gid, rep, f"{lab}/2", rc_["family_id"], "split", rc_["insertion_id"], ev.event_id, (lab,),
                        *h2, st, False, e2_, _uid(p2), a2)]
    s, e, st = pay["m"]
    rm, lf = pay["rm"], pay["leaf"]
    s2, e2 = b0, b0 + len(blk)
    prot = _protein(new, s2, e2, st)
    if prot != lf.protein:
        raise AssertionError(f"{ev.event_id}@{gid}: the NOGD gene does not encode its protein")
    arch = format_arch([(acc, x, y) for acc, (x, y) in zip(ev.data["accessions"], lf.domain_spans)])
    label = f"{lf.lineage}_{ev.event_id}"
    nogd_inst.append(NogdRow(ev.event_id, group_of[rm["copy_label"]][0], rm["family_id"], ev.data["family"], ev.node,
                             gid, rec, rm["copy_label"], label, s2, e2, st, _uid(prot)))
    return [GeneRow(gid, rep, label, ev.data["family"], "nogd", "", ev.event_id, (rm["copy_label"],), s2, e2, st,
                    False, "", _uid(prot), arch)]


def read(world) -> FusionWorld:
    """Load ``world/fusions/``'s events and replicons."""
    world = Path(world)
    fdir = world / "fusions"
    doc = json.loads((fdir / "fusions.json").read_text(encoding="utf-8"))
    if doc["format"] != FORMAT:
        raise ValueError(f"{fdir}: format {doc['format']!r}, expected {FORMAT!r}")
    t = {}
    for name in ("events.tsv", "replicons.tsv"):
        cls = _TABLES[name]
        conv = backbone._parser(cls)
        rows = mobile._tsv(fdir / name)
        t[name] = [cls(*(c(r[f]) for c, f in zip(conv, cls._fields))) for r in rows]
    return FusionWorld(doc["seed"], t["events.tsv"], t["replicons.tsv"], doc["digest"], world, doc)


__all__ = ["FORMAT", "NOGD_FAMILY_BASE", "Edits", "apply_edits", "fusion_edit", "split_edit", "split_point",
           "parse_arch", "format_arch", "check_config", "build", "read", "FusionWorld"]

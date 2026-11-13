"""synthworld contamination, exact and near duplicates, and the mock FCS-GX sidecar (planning/tests/
synthworld_spec.md §2.1 S10, §3.3 ``truth_contam`` and ``truth_dup``, §5.1; card F-02 session 7). Reads a world
whose S1-S9 are built and writes ``world/sources/``.

**Sources.** Every genome's emitted contigs (``assembly/contigs/``) are its original source ``src_a``. Stream
``S10/roles`` permutes the genomes: the first ``contamination.genomes`` are contaminated hosts, the next
``duplicates.exact`` get an exact duplicate (source ``src_b``) and the next ``duplicates.near`` a near duplicate
(source ``src_c``), so no genome has two roles. Accessions ``GCA_SW<nnnnnn>.1`` are numbered over the sources in
(source_label, genome_id) order; BioProjects ``PRJSW<nnnnn>`` in the same order, one per original source except
that the clonal burst's genomes share one, and a new one for each duplicate source; contig occurrences get
``SWC<nnnnnnnn>.1`` in ``source_contigs.tsv`` order.

**Contamination** (stream ``S10/contam/<host>``). The donor is the genome farthest from the host on the species
tree (path length; the first genome id on a tie), so it is a distant clade of the main world (no outgroup world
is built: ``donor_world`` is always main and ``kind`` cross_clade). ``contigs`` windows (uniform on the range)
each of L bases (uniform on ``contig_length``) are cut from a donor emitted contig of at least L bases (uniform
among them, the start uniform; the longest donor contig, whole, when none is long enough) and
reverse-complemented with probability 1/2; they are added to the host's original source as contigs
``<host>_x<k>``. The mock FCS-GX sidecar (MR-19) flags every contaminant contig with action EXCLUDE.

**Duplicates.** An exact duplicate repeats its genome's emitted contigs base for base (equal ``contig_uid``,
MR-06). A near duplicate (stream ``S10/near/<genome>``) changes round((1 - ``near_ani``) x L) bases of each emitted
contig of L bases: round(``igr_share`` x n) on bases in no gene piece (to another base, uniform) and the rest on
synonymous sites, a base of a whole gene piece's internal codon covered by one gene whose change keeps the amino
acid under the genome's code (uniform over the alternatives); a pool too small passes the rest to the other.
Positions are drawn without replacement, uniformly in each pool. Proteins are therefore unchanged.

Shape on disk (``build(world)``)::

    contaminants.fna, .fai          contaminant contigs, records <host>_x<k>
    duplicates/<genome>_<source>.fna, .fai  a duplicate source's contigs, records <contig_id>_<source>
    sources.tsv        source_label genome_id kind rep_genome_id accession bioproject n_contigs length
    source_contigs.tsv source_label genome_id accession contig_accession contig_id origin file record length
                       contig_uid
    contam.tsv         contig_uid host_genome_id donor_genome_id donor_world tree_distance length gc kind
                       in_mock_mask contig_id donor_contig_id donor_start donor_end orientation
    dup.tsv            dup_group_id kind ani_target ani_realized rep_genome_id member_genome_id member_source_label
                       member_accession rep_contig_id member_contig_id rep_contig_uid member_contig_uid n_mutations
                       expected_collapse
    mutations.tsv      member_contig_id rep_contig_id position ref alt site
    sidecars/fcs_gx_mock.tsv  seq_id start_pos end_pos seq_len action div agg_cont_cov top_tax_name
    sources.json       format, seed, the contamination and duplicates presets and their digest, upstream digests,
                       versions, stream paths, roles, counts, sha256 of every file, digest

``origin`` is assembly, contaminant, exact_dup or near_dup; ``file`` and ``record`` say where the bases are,
relative to the world directory. ``donor_start``/``donor_end`` are on the donor contig (forward), ``orientation``
-1 when the contaminant is its reverse complement. ``mutations.position`` is on the contig (both contigs share
coordinates); ``site`` igr or syn. ``expected_collapse`` is true for exact and near duplicates alike (a 99.9% ANI
copy is the same genome for counting, MR-06 and MR-07). ``ani_realized`` = 1 - n_mutations / length. The sidecar's
``seq_id`` is the contaminant's contig accession, positions 1-based inclusive, ``div`` ``prok:<donor phylum>``,
``top_tax_name`` the donor's species. ``digest`` is BLAKE3-128 over the sorted (path, sha256) list of every file but
sources.json.
"""

from __future__ import annotations

import collections
import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
from blake3 import blake3

from biobox import ids
from biobox.synthworld import assembly, backbone, fusions, intergenic, mobile, seeds, sequences

FORMAT = "synthworld-sources/1"
SOURCE_LABELS = {"original": "src_a", "exact_dup": "src_b", "near_dup": "src_c"}
q6 = sequences.q6
revcomp = intergenic.revcomp


class SourceRow(NamedTuple):
    source_label: str
    genome_id: str
    kind: str
    rep_genome_id: str
    accession: str
    bioproject: str
    n_contigs: int
    length: int


class SourceContigRow(NamedTuple):
    source_label: str
    genome_id: str
    accession: str
    contig_accession: str
    contig_id: str
    origin: str
    file: str
    record: str
    length: int
    contig_uid: str


class ContamRow(NamedTuple):
    contig_uid: str
    host_genome_id: str
    donor_genome_id: str
    donor_world: str
    tree_distance: float
    length: int
    gc: float
    kind: str
    in_mock_mask: bool
    contig_id: str
    donor_contig_id: str
    donor_start: int
    donor_end: int
    orientation: int


class DupRow(NamedTuple):
    dup_group_id: str
    kind: str
    ani_target: float
    ani_realized: float
    rep_genome_id: str
    member_genome_id: str
    member_source_label: str
    member_accession: str
    rep_contig_id: str
    member_contig_id: str
    rep_contig_uid: str
    member_contig_uid: str
    n_mutations: int
    expected_collapse: bool


class MutationRow(NamedTuple):
    member_contig_id: str
    rep_contig_id: str
    position: int
    ref: str
    alt: str
    site: str


class FcsRow(NamedTuple):
    seq_id: str
    start_pos: int
    end_pos: int
    seq_len: int
    action: str
    div: str
    agg_cont_cov: int
    top_tax_name: str


_TABLES = {"sources.tsv": SourceRow, "source_contigs.tsv": SourceContigRow, "contam.tsv": ContamRow,
           "dup.tsv": DupRow, "mutations.tsv": MutationRow, "sidecars/fcs_gx_mock.tsv": FcsRow}


def check_config(c: dict, d: dict, n_genomes: int) -> None:
    n = int(c["genomes"]) + int(d["exact"]) + int(d["near"])
    if min(int(c["genomes"]), int(d["exact"]), int(d["near"])) < 0 or n > n_genomes:
        raise ValueError(f"{n} roles (contaminated, exact, near) need as many genomes, the world has {n_genomes}")
    for k in ("contigs", "contig_length"):
        lo, hi = (int(x) for x in c[k])
        if not 1 <= lo <= hi:
            raise ValueError(f"contamination.{k} must be 1 <= lo <= hi, got {c[k]}")
    if not 0.0 < float(d["near_ani"]) < 1.0 or not 0.0 <= float(d["igr_share"]) <= 1.0:
        raise ValueError("duplicates.near_ani lies in (0, 1) and igr_share in [0, 1]")


def tree_distance(tree: mobile.SpeciesTree, x: str, y: str) -> float:
    """Path length between two tips of the species tree."""
    anc = set()
    a = x
    while a:
        anc.add(a)
        a = tree.parent[a]
    m = y
    while m not in anc:
        m = tree.parent[m]
    return (tree.end[x] - tree.end[m]) + (tree.end[y] - tree.end[m])


def synonymous_sites(seq: str, pieces: list[tuple[int, int, int]], code: int,
                     partial: list[tuple[int, int, int]] = ()) -> tuple[dict[int, list[str]], np.ndarray]:
    """Per contig position, the bases (on the contig's + strand) a substitution may take there and keep the amino
    acid: internal codons of the whole gene pieces (start, end, strand; end > len across a circular origin) whose
    three bases lie in one gene piece only; and the number of gene pieces over each base, the ``partial`` pieces
    (genes cut by a contig edge) included."""
    L = len(seq)
    cover = np.zeros(L, dtype=np.int16)
    for s, e, _ in list(pieces) + list(partial):
        if e <= L:
            cover[s:e] += 1
        else:
            cover[s:] += 1
            cover[:e - L] += 1
    t = assembly.code_table(code)
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    out: dict[int, list[str]] = {}
    for s, e, st in pieces:
        n = (e - s) // 3
        for c in range(1, n - 1):
            pos = [((s + 3 * c + i) if st == 1 else (e - 1 - 3 * c - i)) % L for i in range(3)]
            if any(cover[p] != 1 for p in pos):
                continue
            cod = "".join(seq[p] if st == 1 else comp[seq[p]] for p in pos)
            aa = t[cod]
            for i, p in enumerate(pos):
                alts = [b for b in "ACGT" if b != cod[i] and t[cod[:i] + b + cod[i + 1:]] == aa]
                if alts:
                    out[p] = [b if st == 1 else comp[b] for b in alts]
    return out, cover


@dataclass
class SourcesWorld:
    seed: int
    sources: list[SourceRow]
    digest: str
    out: Path
    meta: dict = field(default_factory=dict, repr=False)


def build(world) -> SourcesWorld:
    """Write ``world/sources/`` (module docstring) from a world whose S1-S9 are built."""
    world = Path(world)
    docs = {}
    for d, mod in (("backbone", backbone), ("sequences", sequences), ("replicons", intergenic), ("mobile", mobile),
                   ("fusions", fusions), ("assembly", assembly)):
        docs[d] = json.loads((world / d / f"{d}.json").read_text(encoding="utf-8"))
        if docs[d]["format"] != mod.FORMAT:
            raise ValueError(f"{world}/{d}: format {docs[d]['format']!r}, expected {mod.FORMAT!r}")
    bb = docs["backbone"]
    preset = bb["preset"]
    cc, dc = preset["contamination"], preset["duplicates"]
    seed = int(bb["seed"])
    genomes = mobile._tsv(world / "assembly" / "genomes.tsv")
    gids = sorted(g["genome_id"] for g in genomes)
    code = {g["genome_id"]: int(g["genetic_code"]) for g in genomes}
    check_config(cc, dc, len(gids))
    out = world / "sources"
    if out.exists():
        raise FileExistsError(f"{out} exists; the sources stage is written once")
    tmp = world / "sources.partial"
    if tmp.exists():
        shutil.rmtree(tmp)
    (tmp / "duplicates").mkdir(parents=True)
    (tmp / "sidecars").mkdir()

    branches = backbone.species_branches(world)
    host_genomes = mobile._tsv(world / "backbone" / "genomes.tsv")
    gid_of = {r["species_node"]: r["genome_id"] for r in host_genomes}
    node_of = {g: n for n, g in gid_of.items()}
    tree = mobile.species_tree(branches, gid_of)
    lineage = {r["genome_id"]: r["tax_lineage"] for r in mobile._tsv(world / "backbone" / "taxonomy.tsv")}
    burst = {gid_of[x] for x in bb["burst_tips"]}

    rng = seeds.generator(seed, "S10/roles")
    order = [gids[int(i)] for i in rng.permutation(len(gids))]
    nh, ne, nn = int(cc["genomes"]), int(dc["exact"]), int(dc["near"])
    hosts, exact, near = sorted(order[:nh]), sorted(order[nh:nh + ne]), sorted(order[nh + ne:nh + ne + nn])

    contigs = collections.defaultdict(list)
    for r in mobile._tsv(world / "assembly" / "contigs.tsv"):
        if r["emitted"] == "true":
            contigs[r["genome_id"]].append(r)
    need_seq = set(hosts) | set(exact) | set(near)
    donors = {}
    for h in hosts:
        donors[h] = max((g for g in gids if g != h), key=lambda g: (tree_distance(tree, node_of[h], node_of[g]),
                                                                    [-ord(ch) for ch in g]))
        need_seq.add(donors[h])
    seqs = {}
    for g in sorted(need_seq):
        seqs.update(intergenic._fasta_records(world / "assembly" / "contigs" / f"{g}.fna"))

    # --- contamination
    contam: dict[str, list[tuple[str, str, ContamRow]]] = collections.defaultdict(list)
    fx = sequences._Fasta(tmp / "contaminants.fna")
    for h in hosts:
        r = seeds.generator(seed, f"S10/contam/{h}")
        d = donors[h]
        dist = q6(tree_distance(tree, node_of[h], node_of[d]))
        n = int(r.integers(int(cc["contigs"][0]), int(cc["contigs"][1]) + 1))
        dcs = sorted(contigs[d], key=lambda x: x["contig_id"])
        for k in range(n):
            L = int(r.integers(int(cc["contig_length"][0]), int(cc["contig_length"][1]) + 1))
            ok = [x for x in dcs if int(x["length"]) >= L]
            if ok:
                src = ok[int(r.integers(len(ok)))]
                a = int(r.integers(int(src["length"]) - L + 1))
            else:
                src = max(dcs, key=lambda x: (int(x["length"]), x["contig_id"]))
                L, a = int(src["length"]), 0
            o = -1 if r.random() < 0.5 else 1
            s = seqs[src["contig_id"]][a:a + L]
            s = s if o == 1 else revcomp(s)
            cid = f"{h}_x{k}"
            fx.write(cid, s)
            row = ContamRow(ids.contig_uid(s.encode()).hex(), h, d, "main", dist, L,
                            q6((s.count("G") + s.count("C")) / L), "cross_clade", True, cid, src["contig_id"], a,
                            a + L, o)
            contam[h].append((cid, s, row))
    fx.close()

    # --- duplicates
    dup_seqs: dict[str, list[tuple[str, str, str, str, int]]] = {}   # genome -> (member id, seq, rep id, rep seq, n)
    mut_rows: list[MutationRow] = []
    pieces = collections.defaultdict(list)
    partial = collections.defaultdict(list)
    if near:
        want = set(near)
        for p in mobile._tsv(world / "assembly" / "gene_contigs.tsv"):
            if p["genome_id"] in want and p["emitted"] == "true":
                whole = p["partial_5"] == p["partial_3"] == "false"
                (pieces if whole else partial)[p["contig_id"]].append((int(p["start"]), int(p["end"]),
                                                                        int(p["strand"])))
    gap = 1.0 - float(dc["near_ani"])
    for kind, group in (("exact_dup", exact), ("near_dup", near)):
        lab = SOURCE_LABELS[kind]
        for g in group:
            fd = sequences._Fasta(tmp / "duplicates" / f"{g}_{lab}.fna")
            rows = []
            r = seeds.generator(seed, f"S10/near/{g}") if kind == "near_dup" else None
            for c in sorted(contigs[g], key=lambda x: x["contig_id"]):
                rep = seqs[c["contig_id"]]
                mid = f"{c['contig_id']}_{lab}"
                if kind == "exact_dup":
                    new, n_mut = rep, 0
                else:
                    L = len(rep)
                    n_mut = int(round(gap * L))
                    syn, cover = synonymous_sites(rep, pieces[c["contig_id"]], code[g], partial[c["contig_id"]])
                    igr = np.flatnonzero(cover == 0)
                    syn_pos = np.array(sorted(syn), dtype=np.int64)
                    n_igr = min(int(round(float(dc["igr_share"]) * n_mut)), len(igr))
                    n_syn = min(n_mut - n_igr, len(syn_pos))
                    n_igr = min(n_mut - n_syn, len(igr))
                    pick_i = np.sort(r.choice(igr, size=n_igr, replace=False)) if n_igr else np.array([], np.int64)
                    pick_s = np.sort(r.choice(syn_pos, size=n_syn, replace=False)) if n_syn else np.array([], np.int64)
                    buf = bytearray(rep.encode())
                    muts = []
                    for p in pick_i.tolist():
                        alt = [b for b in "ACGT" if b != rep[p]][int(r.integers(3))]
                        muts.append((p, alt, "igr"))
                    for p in pick_s.tolist():
                        alts = syn[p]
                        muts.append((p, alts[int(r.integers(len(alts)))], "syn"))
                    for p, alt, site in sorted(muts):
                        buf[p] = ord(alt)
                        mut_rows.append(MutationRow(mid, c["contig_id"], p, rep[p], alt, site))
                    new = buf.decode()
                    n_mut = len(muts)
                fd.write(mid, new)
                rows.append((mid, new, c["contig_id"], rep, n_mut))
            fd.close()
            dup_seqs[g] = rows

    # --- sources, accessions, BioProjects
    srcs = [("src_a", g, "original", "") for g in gids]
    srcs += [(SOURCE_LABELS["exact_dup"], g, "exact_dup", g) for g in exact]
    srcs += [(SOURCE_LABELS["near_dup"], g, "near_dup", g) for g in near]
    srcs.sort(key=lambda x: (x[0], x[1]))
    bp, burst_bp, n_bp = {}, None, 0
    for lab, g, kind, _ in srcs:
        if kind == "original" and g in burst:
            if burst_bp is None:
                n_bp += 1
                burst_bp = f"PRJSW{n_bp:05d}"
            bp[(lab, g)] = burst_bp
        else:
            n_bp += 1
            bp[(lab, g)] = f"PRJSW{n_bp:05d}"
    w_src = mobile._Rows(tmp / "sources.tsv", SourceRow)
    w_sc = mobile._Rows(tmp / "source_contigs.tsv", SourceContigRow)
    dup_rows, fcs_rows, contam_rows = [], [], []
    n_acc = 0
    counts = collections.Counter()
    for i, (lab, g, kind, rep) in enumerate(srcs):
        acc = f"GCA_SW{i + 1:06d}.1"
        occ = []
        if kind == "original":
            for c in sorted(contigs[g], key=lambda x: x["contig_id"]):
                occ.append((c["contig_id"], "assembly", f"assembly/contigs/{g}.fna", c["contig_id"], int(c["length"]),
                            c["contig_uid"], None))
            for cid, s, row in contam[g]:
                occ.append((cid, "contaminant", "sources/contaminants.fna", cid, len(s), row.contig_uid, row))
        else:
            for mid, s, rid, rseq, n_mut in dup_seqs[g]:
                occ.append((mid, kind, f"sources/duplicates/{g}_{lab}.fna", mid, len(s),
                            ids.contig_uid(s.encode()).hex(), (rid, rseq, n_mut)))
        for cid, origin, f, rec, L, cuid, extra in occ:
            n_acc += 1
            cacc = f"SWC{n_acc:08d}.1"
            w_sc.row(SourceContigRow(lab, g, acc, cacc, cid, origin, f, rec, L, cuid))
            counts[origin] += 1
            if origin == "contaminant":
                contam_rows.append(extra)
                d = extra.donor_genome_id
                tax = lineage[d].split(";")
                phylum = next((x for x in tax if x.startswith("p__")), "p__")
                fcs_rows.append(FcsRow(cacc, 1, L, L, "EXCLUDE", f"prok:{phylum[3:]}", 100, tax[-1]))
            elif origin in ("exact_dup", "near_dup"):
                rid, rseq, n_mut = extra
                gi = (exact if origin == "exact_dup" else near).index(g)
                dup_rows.append(DupRow(f"{'EX' if origin == 'exact_dup' else 'NE'}{gi:02d}",
                                       "exact" if origin == "exact_dup" else "near",
                                       1.0 if origin == "exact_dup" else float(dc["near_ani"]),
                                       q6(1.0 - n_mut / L), g, g, lab, cacc, rid, cid,
                                       ids.contig_uid(rseq.encode()).hex(), cuid, n_mut, True))
        w_src.row(SourceRow(lab, g, kind, rep, acc, bp[(lab, g)], len(occ), sum(x[4] for x in occ)))
    w_src.close()
    w_sc.close()
    for name, rows in (("contam.tsv", contam_rows), ("dup.tsv", dup_rows), ("mutations.tsv", mut_rows),
                       ("sidecars/fcs_gx_mock.tsv", fcs_rows)):
        t = mobile._Rows(tmp / name, _TABLES[name])
        for x in rows:
            t.row(x)
        t.close()
    counts.update(sources=len(srcs), mutations=len(mut_rows))

    files = {str(f.relative_to(tmp)): sequences._sha(f) for f in sorted(tmp.rglob("*")) if f.is_file()}
    digest = blake3("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest(length=16)
    cfg = {"contamination": cc, "duplicates": dc}
    doc = {"format": FORMAT, "seed": seed, "preset": cfg, "params_digest": backbone.params_digest(cfg),
           **{f"{d}_digest": docs[d]["digest"] for d in docs},
           "versions": {"numpy": np.__version__},
           "roles": {"contaminated": hosts, "donors": donors, "exact_dup": exact, "near_dup": near},
           "stream_paths": ["S10/roles", "S10/contam/<genome>", "S10/near/<genome>"],
           "counts": dict(sorted(counts.items())), "files": files, "digest": digest}
    (tmp / "sources.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    tmp.rename(out)
    return read(world)


def read(world) -> SourcesWorld:
    """Load ``world/sources/``'s source table."""
    world = Path(world)
    sdir = world / "sources"
    doc = json.loads((sdir / "sources.json").read_text(encoding="utf-8"))
    if doc["format"] != FORMAT:
        raise ValueError(f"{sdir}: format {doc['format']!r}, expected {FORMAT!r}")
    conv = backbone._parser(SourceRow)
    rows = [SourceRow(*(c(r[f]) for c, f in zip(conv, SourceRow._fields))) for r in mobile._tsv(sdir / "sources.tsv")]
    return SourcesWorld(doc["seed"], rows, doc["digest"], world, doc)


__all__ = ["FORMAT", "SOURCE_LABELS", "check_config", "tree_distance", "synonymous_sites", "build", "read",
           "SourcesWorld"]

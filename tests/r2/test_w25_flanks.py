"""Ring R2: `bio-box nbr flanks` groups a search's loci by their flanks (card W-25; planning/design/W-25_flanks.md).

A constructed release with planted contexts, written table by table (tapes, protein occurrences, clusters, domains,
genome inputs, contig_dim, labels) over two stages, S1 (phages, kind virus) and P1 (prokaryotes). Every gene carries a
protein of its own; a family is a cluster at NR50 and NR90 and one Pfam-like domain of that name. Window 5 genes.

    T1      up A1+ A2+ A3- A4+ A5+, down B1+ B2- B3+ B4+ B5+ (strands relative to the anchor), on two phage and
            two prokaryote contigs, one of each with the anchor on the minus strand
    ins     T1 with X inserted between A1 and A2 (A5 leaves the window)
    un      T1 with A3 replaced by a gene with no cluster and no domain
    cens    T1 cut by the start of a linear contig after A2 (two genes upstream)
    rearr   T1 with the down side in the order B3 B1 B5 B2 B4
    T2      up C1..C5, down D1..D5 (D4 on the opposite strand), two prokaryote contigs
    U1, U2  ten genes with no cluster and no domain each; their anchors are another NR90 cluster (R2) of the same
            NR50 cluster as the others (R)

The oracle (`oracle_*` below) is a brute-force loop over the planted lists, written apart from `biobox.search.flanks`:
it builds the tokens from the planted family names, truncates as the censoring rule says and counts sets by loops. Its
average linkage recomputes every cluster distance from the loci at each step (no update formula, no collapsing).
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pytest

from biobox import catalog, cli, ids
from biobox import io as bio
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l3 import tapes
from biobox.registry.schemas import table
from biobox.search import flanks

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-25")]

REL = "rW25.1"
STAGES = ["S1", "P1"]
TV = tapes.tape_version([GENE_CALL_VERSION])
W = 5
GENE, GAP = 300, 60

A = [("A1", 1), ("A2", 1), ("A3", -1), ("A4", 1), ("A5", 1)]
B = [("B1", 1), ("B2", -1), ("B3", 1), ("B4", 1), ("B5", 1)]
C = [(f"C{k}", 1) for k in range(1, 6)]
D = [("D1", 1), ("D2", 1), ("D3", 1), ("D4", -1), ("D5", 1)]
NONE5 = [(None, 1), (None, -1), (None, 1), (None, 1), (None, -1)]
EXTRA = [("E", 1)]  # one gene beyond the window on a side that is not cut

# name: (stage, kind, anchor family, anchor strand, up outward, down outward, up cut by the contig start)
PLANTED = {
    "t1_v1": ("S1", "virus", "R", 1, A + EXTRA, B + EXTRA, False),
    "t1_v2": ("S1", "virus", "R", -1, A + EXTRA, B + EXTRA, False),
    "t1_p1": ("P1", "prokaryote", "R", 1, A + EXTRA, B + EXTRA, False),
    "t1_p2": ("P1", "prokaryote", "R", -1, A + EXTRA, B + EXTRA, False),
    "t1_ins": ("P1", "prokaryote", "R", 1, [A[0], ("X", 1)] + A[1:] + EXTRA, B + EXTRA, False),
    "t1_un": ("S1", "virus", "R", 1, A[:2] + [(None, -1)] + A[3:] + EXTRA, B + EXTRA, False),
    "t1_cens": ("P1", "prokaryote", "R", 1, A[:2], B + EXTRA, True),
    "t1_rearr": ("S1", "virus", "R", 1, A + EXTRA, [B[2], B[0], B[4], B[1], B[3]] + EXTRA, False),
    "t2_p1": ("P1", "prokaryote", "R", 1, C + EXTRA, D + EXTRA, False),
    "t2_p2": ("P1", "prokaryote", "R", -1, C + EXTRA, D + EXTRA, False),
    "u1": ("S1", "virus", "R2", 1, NONE5 + EXTRA, NONE5 + EXTRA, False),
    "u2": ("P1", "prokaryote", "R2", 1, NONE5 + EXTRA, NONE5 + EXTRA, False),
}
TRUTH = {"T1": ["t1_v1", "t1_v2", "t1_p1", "t1_p2", "t1_ins", "t1_un", "t1_cens"], "T1r": ["t1_rearr"],
         "T2": ["t2_p1", "t2_p2"], "U1": ["u1"], "U2": ["u2"]}
SHARED_GENOME = ("t1_v1", "t1_p1")  # one more genome (in P1) holds these two contigs: genome counts go through it


def _uid(text: str) -> bytes:
    return ids.protein_uid(text)


def _nr90(fam: str) -> bytes:
    return _uid(f"nr90:{fam}")


def _nr50(fam: str) -> bytes:
    return _uid("nr50:R" if fam in ("R", "R2") else f"nr50:{fam}")


def _write(name: str, rows: list[dict], path: Path, key) -> None:
    spec = table(name)
    bio.write_table([pa.Table.from_pylist(sorted(rows, key=key), schema=spec.schema)], spec, path)


def build(root: Path) -> dict:
    """Write the constructed release under `root`; return {locus name: facts} for the oracle."""
    loci, tapes_rows, occ, gin, labels, domains, members, dims = {}, {s: [] for s in STAGES}, {s: [] for s in STAGES}, \
        {s: [] for s in STAGES}, {s: [] for s in STAGES}, {s: [] for s in STAGES}, [], []
    families = sorted({f for v in PLANTED.values() for side in (v[4], v[5]) for f, _ in side if f} | {"R", "R2"})
    fam_members: dict[str, list[bytes]] = {f: [] for f in families}
    genomes: dict[bytes, list[bytes]] = {}
    for name, (stage, kind, afam, astrand, up, down, cut) in PLANTED.items():
        cu = ids.contig_uid(f"contig {name}".encode())
        # coordinate order: on a + anchor the up side lies at lower coordinates
        genes = ([(f, s * astrand, "u", k) for k, (f, s) in enumerate(up)][::-1] + [(afam, astrand, "anchor", 0)]
                 + [(f, s * astrand, "d", k) for k, (f, s) in enumerate(down)])
        if astrand == -1:
            genes = genes[::-1]
        pos = 0 if cut else 100
        slots, pus, fus = [], [], []
        anchor = None
        for i, (fam, strand, role, k) in enumerate(genes):
            start, end = pos, pos + GENE
            pos = end + GAP
            pu = _uid(f"M{name}G{i}")
            fu = ids.feature_uid(cu, start, end, strand, 1)
            slots.append((start, end, strand))
            pus.append(pu)
            fus.append(fu)
            if role == "anchor":
                anchor = (pu, fu, start, end)
            if fam:
                fam_members[fam].append(pu)
            occ[stage].append({"protein_uid": pu, "contig_uid": cu, "start": start, "end": end, "strand": strand,
                               "dist_to_contig_end_bp": 0, "feature_uid": fu, "_fam": fam})
        length = pos - GAP + 100
        tape = tapes.build_tape(slots, length, False)
        tapes_rows[stage].append({"contig_uid": cu, "contig_len": length, "is_circular": False, "n_slots": len(slots),
                                  **tape, "protein_uids": pus, "feature_uids": fus})
        labels[stage].append({"contig_uid": cu, "accession": name.upper(), "original_header": f"{name} constructed",
                              "n_records": 1, "l0_object_sha256": bytes(32), "record_ordinal": 0})
        genomes.setdefault(cu, []).append(ids.genome_uid([cu]))
        gin[stage].append({"contig_uid": cu, "genome_uid": ids.genome_uid([cu]), "is_circular": False, "kind": kind,
                           "genetic_code": 11, "assembly_quality": "complete"})
        loci[name] = {"contig_uid": cu, "anchor_protein_uid": anchor[0], "feature_uid": anchor[1], "stage": stage,
                      "kind": kind, "anchor_family": afam, "up": up[:W], "down": down[:W], "cut_up": cut,
                      "contig_len": length}
    shared = [loci[n]["contig_uid"] for n in SHARED_GENOME]
    gx = ids.genome_uid(shared)
    for cu in shared:
        genomes[cu].append(gx)
        gin["P1"].append({"contig_uid": cu, "genome_uid": gx, "is_circular": False, "kind": "prokaryote",
                          "genetic_code": 11, "assembly_quality": "complete"})
    for name, f in loci.items():
        f["genomes"] = sorted(genomes[f["contig_uid"]])
    for i, fam in enumerate(families):
        rep = min(fam_members[fam])
        for pu in fam_members[fam]:
            members.append({"member_uid": pu, "level": "nr90", "cluster_uid": _nr90(fam), "rep_uid": rep})
            members.append({"member_uid": pu, "level": "nr50", "cluster_uid": _nr50(fam), "rep_uid": rep})
    for stage in STAGES:
        for r in occ[stage]:
            fam = r.pop("_fam")
            if fam:
                domains[stage].append({"protein_uid": r["protein_uid"], "model_acc": f"PF{families.index(fam):05d}.1",
                                       "model_name": fam, "model_len": 100, "seq_bitscore": 50.0,
                                       "seq_evalue_zref": 1e-10, "dom_bitscore": 50.0, "dom_ievalue_zref": 1e-10,
                                       "env_start": 1, "env_end": 90, "ali_start": 1, "ali_end": 90, "hmm_from": 1,
                                       "hmm_to": 100, "cutoff": "gathering", "z": 1e8, "dom_z": 1e8})
    for cu, gs in genomes.items():
        dims.append({"contig_uid": cu, "genome_uid": min(gs),
                     "contig_len": next(f["contig_len"] for f in loci.values() if f["contig_uid"] == cu),
                     "is_circular": False, "tier": None, "assembly_quality": "complete", "tax_order": None,
                     "is_species_rep": False, "n_genomes": len(gs)})
    for stage in STAGES:
        _write("contig_tapes", tapes_rows[stage], root / f"l3/tapes/tape_version={TV}/stage={stage}/part-0.parquet",
               lambda r: r["contig_uid"])
        _write("protein_occurrence", occ[stage],
               root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage={stage}/part-0.parquet",
               lambda r: (r["protein_uid"], r["contig_uid"], r["start"]))
        _write("genome_input", gin[stage], root / f"l0seq/genome_input/stage={stage}/part-0.parquet",
               lambda r: (r["contig_uid"], r["genome_uid"]))
        _write("contig_label", labels[stage], root / f"l0seq/contig_label/stage={stage}/part-0.parquet",
               lambda r: r["contig_uid"])
        _write("l1_domains", domains[stage], root / f"l1/domains/library=pfam-38.2/stage={stage}/part-0.parquet",
               lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
    _write("member_cluster", members, root / f"l2/member_cluster/release={REL}/part-0.parquet",
           lambda r: (r["level"], r["member_uid"]))
    _write("contig_dim", dims, root / f"dim/contig_dim/release={REL}/part-0.parquet", lambda r: r["contig_uid"])
    catalog.define_release(root, REL, STAGES)
    return loci


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    root = tmp_path_factory.mktemp("w25") / "root"
    return root, build(root)


def anchors(loci: dict, names=None) -> str:
    return ",".join(sorted(loci[n]["anchor_protein_uid"].hex() for n in (names or loci)))


def run_flanks(world, out: Path, *extra: str, names=None) -> int:
    root, loci = world
    return cli.main(["nbr", "flanks", anchors(loci, names), "--root", str(root), "--release", REL,
                     "--out", str(out), *extra])


def read_tsv(path: Path) -> list[dict]:
    lines = [x for x in path.read_text().splitlines() if x and not x.startswith("#")]
    cols = lines[0].split("\t")
    return [dict(zip(cols, x.split("\t"))) for x in lines[1:]]


def read_matrix(path: Path) -> dict[tuple[str, str], float]:
    lines = path.read_text().splitlines()
    names = lines[0].split("\t")[1:]
    out = {}
    for line in lines[1:]:
        cells = line.split("\t")
        for n, v in zip(names, cells[1:]):
            out[(cells[0], n)] = float(v)
    return out


def by_name(loci: dict) -> dict[str, str]:
    return {f["feature_uid"].hex(): n for n, f in loci.items()}


# --- the oracle: brute force over the planted lists -----------------------------------------------------------------


def oracle_tokens(loci: dict, name: str, side: str) -> list[tuple]:
    out = []
    for k, (fam, strand) in enumerate(loci[name]["up" if side == "u" else "down"]):
        out.append(((fam if fam else ("unannotated", name, side, k)), side, strand))
    return out


def oracle_obs(loci: dict, name: str, side: str) -> int:
    if side == "u" and loci[name]["cut_up"]:
        return len(loci[name]["up"])
    return W


def _distinct(items: list) -> list:
    seen = []
    for x in items:
        if x not in seen:
            seen.append(x)
    return seen


def _jaccard_distance(xs: list, ys: list) -> float:
    xs, ys = _distinct(xs), _distinct(ys)
    inter = 0
    for x in xs:
        for y in ys:
            if x == y:
                inter += 1
    union = len(xs) + len(ys) - inter
    return 1.0 if union == 0 else 1 - inter / union


def oracle_pair(loci: dict, x: str, y: str) -> tuple[float, float, int, list, list]:
    """(J, A, positions compared, X's pairs, Y's pairs)."""
    toks = {x: [], y: []}
    pairs = {x: [], y: []}
    positions = 0
    for side in ("u", "d"):
        m = min(oracle_obs(loci, x, side), oracle_obs(loci, y, side))
        positions += m
        for n in (x, y):
            prev = "anchor"
            for t in oracle_tokens(loci, n, side)[:m]:
                toks[n].append(t)
                pairs[n].append((side, prev, t))
                prev = t
    return (_jaccard_distance(toks[x], toks[y]), _jaccard_distance(pairs[x], pairs[y]), positions, pairs[x],
            pairs[y])


def oracle_distance(loci: dict, x: str, y: str, kind: str, w: float = 0.0) -> float:
    j, a, *_ = oracle_pair(loci, x, y)
    return {"j": j, "a": a, "combined": 1 - (w * (1 - j) + (1 - w) * (1 - a))}[kind]


def oracle_types(loci: dict, names: list[str], kind: str, cut: float, w: float = 0.0) -> set[frozenset]:
    """Average linkage, every cluster distance recomputed from its loci, merged while the least is <= cut."""
    clusters = [[n] for n in sorted(names)]
    while len(clusters) > 1:
        best = None
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                total = 0.0
                for p in clusters[i]:
                    for q in clusters[j]:
                        total += oracle_distance(loci, p, q, kind, w)
                d = total / (len(clusters[i]) * len(clusters[j]))
                if best is None or d < best[0]:
                    best = (d, i, j)
        if best[0] > cut + 1e-9:
            break
        _, i, j = best
        clusters[i] = clusters[i] + clusters[j]
        del clusters[j]
    return {frozenset(c) for c in clusters}


def types_of(out: Path, loci: dict) -> dict[tuple[str, str], set[str]]:
    names = by_name(loci)
    got: dict[tuple[str, str], set[str]] = {}
    for r in read_tsv(out / "context_types.tsv"):
        if r["in_tree"] == "true":
            got.setdefault((r["group"], r["context_type"]), set()).add(names[r["feature_uid"]])
    return got


# --- tests ----------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("level", ["nr50", "pfam"])
def test_w25_i_two_flank_types_one_shared_by_phages_and_prokaryotes_are_two_context_types(world, tmp_path, level):
    """(i) Around the R anchors the T1 loci (with an insertion, an unannotated gene, a cut and a rearranged copy) and
    the T2 loci are two context types at the default cut and distance, equal to the oracle's average linkage; T1 holds
    phage and prokaryote loci. The U loci (anchors in NR90 cluster R2) are their own group, one type each."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out, "--level", level) == 0
    got = types_of(out, loci)
    group_r, group_r2 = _nr90("R").hex(), _nr90("R2").hex()
    want_r = oracle_types(loci, [n for n in loci if loci[n]["anchor_family"] == "R"], "j", flanks.DEFAULT_CUT)
    assert {frozenset(v) for (g, _), v in got.items() if g == group_r} == want_r
    assert want_r == {frozenset(TRUTH["T1"] + TRUTH["T1r"]), frozenset(TRUTH["T2"])}
    assert {frozenset(v) for (g, _), v in got.items() if g == group_r2} == {frozenset(["u1"]), frozenset(["u2"])}
    rows = {by_name(loci)[r["feature_uid"]]: r for r in read_tsv(out / "context_types.tsv")}
    t1 = {rows[n]["context_type"] for n in TRUTH["T1"]}
    assert len(t1) == 1 and {rows[n]["kind"] for n in TRUTH["T1"]} == {"virus", "prokaryote"}
    assert {rows[n]["stage"] for n in TRUTH["T1"]} == {"S1", "P1"}
    assert rows["t1_v1"]["n_genomes"] == "2" and rows["t1_v2"]["n_genomes"] == "1"
    assert all("ordinal" not in c for c in rows["t1_v1"])  # AGENTS §3.3: genes are named by feature_uid


def _pair_values(out: Path, loci: dict, group: str):
    rows = read_tsv(out / "context_types.tsv")
    leaf = {by_name(loci)[r["feature_uid"]]: r["leaf"] for r in rows}
    d = out / "groups" / group
    return leaf, read_matrix(d / "distance_j.tsv"), read_matrix(d / "distance_a.tsv"), read_matrix(d / "positions.tsv")


def test_w25_distances_equal_the_brute_force_oracle_for_every_pair(world, tmp_path):
    """Every pair of loci of the R group: J, A and the positions compared equal the oracle's (4 decimals), and two
    loci in one leaf are at distance 0 by the oracle."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out) == 0
    names = [n for n in loci if loci[n]["anchor_family"] == "R"]
    leaf, mj, ma, mp = _pair_values(out, loci, _nr90("R").hex())
    for x in names:
        for y in names:
            j, a, pos, *_ = oracle_pair(loci, x, y)
            if leaf[x] == leaf[y]:
                assert (j, a) == (0, 0), (x, y)
                continue
            assert mj[(leaf[x], leaf[y])] == pytest.approx(j, abs=5e-5), (x, y)
            assert ma[(leaf[x], leaf[y])] == pytest.approx(a, abs=5e-5), (x, y)
            assert mp[(leaf[x], leaf[y])] == pos, (x, y)


def test_w25_ii_a_rearranged_order_is_separated_by_a_and_not_by_j(world, tmp_path):
    """(ii) The rearranged copy has T1's families on the same sides and strands: J = 0, A = 1 - 5/15. With
    --distance a it is a context type of its own, as the oracle's average linkage on A says."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out) == 0
    leaf, mj, ma, _ = _pair_values(out, loci, _nr90("R").hex())
    assert mj[(leaf["t1_v1"], leaf["t1_rearr"])] == 0
    assert ma[(leaf["t1_v1"], leaf["t1_rearr"])] == pytest.approx(1 - 5 / 15, abs=5e-5)
    out_a = tmp_path / "out_a"
    assert run_flanks(world, out_a, "--distance", "a") == 0
    names = [n for n in loci if loci[n]["anchor_family"] == "R"]
    want = oracle_types(loci, names, "a", flanks.DEFAULT_CUT)
    assert frozenset(["t1_rearr"]) in want
    assert {frozenset(v) for (g, _), v in types_of(out_a, loci).items() if g == _nr90("R").hex()} == want


def test_w25_iii_one_inserted_gene_changes_a_by_at_most_the_pairs_it_touches(world, tmp_path):
    """(iii) X between A1 and A2: the pair sets differ by the broken pair (A1, A2), the two new pairs and the pair
    (A4, A5) that leaves the window at its far end, four pairs; A = 1 - 8/12 and J = 1 - 9/11, both as the oracle.
    A positional comparison would lose every position beyond the insertion."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out) == 0
    leaf, mj, ma, _ = _pair_values(out, loci, _nr90("R").hex())
    j, a, _, px, py = oracle_pair(loci, "t1_v1", "t1_ins")
    touched = [p for p in px if p not in py] + [p for p in py if p not in px]
    assert len(touched) <= 4
    union = len(_distinct(px + py))
    assert a <= len(touched) / union + 1e-12
    assert ma[(leaf["t1_v1"], leaf["t1_ins"])] == pytest.approx(1 - 8 / 12, abs=5e-5) == pytest.approx(a, abs=5e-5)
    assert mj[(leaf["t1_v1"], leaf["t1_ins"])] == pytest.approx(1 - 9 / 11, abs=5e-5)


def test_w25_iv_a_censored_locus_is_compared_only_over_the_positions_it_observes(world, tmp_path):
    """(iv) The cut locus observes two upstream genes: against T1 it is compared over 2 + 5 positions and is at
    distance 0 by J and A; T1 against the insertion copy is compared over 10. The locus row says the side is cut."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out) == 0
    leaf, mj, ma, mp = _pair_values(out, loci, _nr90("R").hex())
    assert mp[(leaf["t1_v1"], leaf["t1_cens"])] == 7 and mp[(leaf["t1_v1"], leaf["t1_ins"])] == 10
    assert mj[(leaf["t1_v1"], leaf["t1_cens"])] == 0 and ma[(leaf["t1_v1"], leaf["t1_cens"])] == 0
    assert leaf["t1_v1"] != leaf["t1_cens"]  # a cut locus never collapses with a whole one
    rows = {by_name(loci)[r["feature_uid"]]: r for r in read_tsv(out / "context_types.tsv")}
    assert (rows["t1_cens"]["cut_up"], rows["t1_cens"]["n_up"], rows["t1_cens"]["cut_down"]) == ("true", "2", "false")
    assert (rows["t1_v1"]["cut_up"], rows["t1_v1"]["n_up"]) == ("false", "5")


@pytest.mark.parametrize("level", ["nr50", "pfam"])
def test_w25_v_unannotated_neighbours_do_not_make_two_loci_similar(world, tmp_path, level):
    """(v) U1 and U2 hold ten genes each with no cluster and no domain: J = A = 1 between them (each such gene is a
    token of its own), they stay two leaves, and in one group with every locus neither joins anything."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out, "--level", level, "--group-by", "none") == 0
    leaf, mj, ma, _ = _pair_values(out, loci, "all")
    assert leaf["u1"] != leaf["u2"]
    assert mj[(leaf["u1"], leaf["u2"])] == 1 and ma[(leaf["u1"], leaf["u2"])] == 1
    got = {frozenset(v) for v in types_of(out, loci).values()}
    assert frozenset(["u1"]) in got and frozenset(["u2"]) in got


def test_w25_w_is_the_largest_weight_on_j_that_still_separates_the_planted_types(world):
    """The combined distance's w (design note): the largest w on 0, 0.05, ..., 1 at which every between-type distance
    of the planted truth (T1 with its variants, the rearranged copy, T2, U1, U2) exceeds every within-type distance,
    computed with the oracle."""
    _, loci = world
    names = [n for v in TRUTH.values() for n in v]
    truth = {n: t for t, v in TRUTH.items() for n in v}
    best = None
    for k in range(21):
        w = k / 20
        within, between = [], []
        for i, x in enumerate(names):
            for y in names[i + 1:]:
                (within if truth[x] == truth[y] else between).append(oracle_distance(loci, x, y, "combined", w))
        if min(between) - max(within) > 0:
            best = w
    assert best == flanks.W_COMBINED


def test_w25_conservation_table_equals_the_oracle_counts(world, tmp_path):
    """Per context type and family: loci holding it, share, modal offset along the anchor's strand, same-strand share
    and distinct genomes (the shared genome counted through genome_input), from a loop over the planted lists."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out, "--level", "pfam") == 0
    rows = read_tsv(out / "conservation.tsv")
    ct = {by_name(loci)[r["feature_uid"]]: r["context_type"] for r in read_tsv(out / "context_types.tsv")}
    t1 = TRUTH["T1"] + TRUTH["T1r"]
    t1_type = ct["t1_v1"]
    group = _nr90("R").hex()
    for fam in ("A1", "A3", "B2", "X"):
        holders, offsets, same, genomes = set(), {}, 0, set()
        for n in t1:
            for side, sign in (("up", -1), ("down", 1)):
                for k, (f, strand) in enumerate(loci[n][side]):
                    if f == fam:
                        holders.add(n)
                        offsets[sign * (k + 1)] = offsets.get(sign * (k + 1), 0) + 1
                        same += strand == 1
                        genomes.update(loci[n]["genomes"])
        modal = max(offsets, key=lambda o: (offsets[o], -abs(o), o))
        got = [r for r in rows if r["group"] == group and r["context_type"] == t1_type and r["key"] == fam]
        assert len(got) == 1, fam
        r = got[0]
        assert int(r["n_loci"]) == len(holders) and int(r["n_loci_type"]) == len(t1), fam
        assert float(r["share_of_loci"]) == pytest.approx(round(len(holders) / len(t1), 4)), fam
        assert int(r["modal_offset"]) == modal, fam
        assert float(r["same_strand_share"]) == pytest.approx(round(same / sum(offsets.values()), 4)), fam
        assert int(r["n_genomes"]) == len(genomes), fam
    unannotated = [r for r in rows if r["group"] == group and r["context_type"] == t1_type and r["unannotated"] == "true"]
    assert len(unannotated) == 1 and unannotated[0]["n_loci"] == "1"


def test_w25_above_k_leaves_the_sample_is_the_least_leaf_names_and_says_so(world, tmp_path, capsys):
    """--max-loci 3: of the R group's six leaves the three with the least names (each the least feature_uid of its
    loci) are in the tree, the others in_tree=false without a context type; groups.tsv and stderr say so."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out, "--max-loci", "3") == 0
    err = capsys.readouterr().err
    group = _nr90("R").hex()
    rows = [r for r in read_tsv(out / "context_types.tsv") if r["group"] == group]
    leaves = sorted({r["leaf"] for r in rows})
    assert len(leaves) == 6
    for r in rows:
        assert (r["in_tree"] == "true") == (r["leaf"] in leaves[:3])
        assert (r["context_type"] == "") == (r["in_tree"] == "false")
    leaf_of = {}
    for r in rows:
        leaf_of.setdefault(r["leaf"], []).append(r["feature_uid"])
    assert all(name == min(fus) for name, fus in leaf_of.items())
    g = {r["group"]: r for r in read_tsv(out / "groups.tsv")}[group]
    assert (g["n_leaves"], g["n_leaves_in_tree"], g["sampled"]) == ("6", "3", "true")
    assert group in err and "--max-loci 3" in err
    assert len(read_matrix(out / "groups" / group / "distance_j.tsv")) == 9


def test_w25_tree_is_newick_over_the_leaves_in_tree_with_ultrametric_heights(world, tmp_path):
    """The Newick tree names every leaf in the tree once, and every leaf's path to the root has the same length."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out) == 0
    group = _nr90("R").hex()
    text = (out / "groups" / group / "tree.nwk").read_text().strip()
    assert text.endswith(";") and text.count("(") == text.count(")")
    names = re.findall(r"[(,]([0-9a-f]{32}):", text)
    leaves = [r["leaf"] for r in read_tsv(out / "groups" / group / "leaves.tsv") if r["in_tree"] == "true"]
    assert sorted(names) == sorted(leaves)
    # the path length of each leaf: its own branch length plus those of the clades that hold it
    spans, lengths = [], {}
    for m in re.finditer(r"\):([0-9.]+)", text):
        close = m.start()
        level, i = 0, close
        while True:
            if text[i] == ")":
                level += 1
            elif text[i] == "(":
                level -= 1
                if level == 0:
                    break
            i -= 1
        spans.append((i, close, float(m.group(1))))
    for m in re.finditer(r"([0-9a-f]{32}):([0-9.]+)", text):
        lengths[m.group(1)] = float(m.group(2)) + sum(bl for a, b, bl in spans if a < m.start() < b)
    assert max(lengths.values()) - min(lengths.values()) < 1e-5


def test_w25_gff_export_has_one_file_per_locus_with_tokens_and_a_protein_group_table(world, tmp_path):
    """--export gff: one GFF3 per locus; the window's genes and the anchor as CDS rows, 1-based, with the gene's
    feature_uid, protein_uid and token; and protein_uid -> group for every protein in them."""
    root, loci = world
    out = tmp_path / "out"
    assert run_flanks(world, out, "--export", "gff", "--level", "pfam") == 0
    files = sorted((out / "gff").glob("locus_*.gff3"))
    assert len(files) == len(loci)
    f = loci["t1_v2"]
    path = out / "gff" / f"locus_T1_V2_{f['feature_uid'].hex()}.gff3"
    lines = path.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    cds = [x.split("\t") for x in lines if x and not x.startswith("#") and x.split("\t")[2] == "CDS"]
    assert len(cds) == 11  # five up, the anchor, five down
    anchor = [c for c in cds if f"ID={f['feature_uid'].hex()};" in c[8]]
    assert len(anchor) == 1 and anchor[0][6] == "-" and "anchor=true" in anchor[0][8]
    tokens = [re.search(r"token=([^;]+)", c[8]).group(1) for c in cds]
    assert "A3" in tokens and "B2" in tokens
    groups = read_tsv(out / "gff" / "protein_groups.tsv")
    assert {r["protein_uid"] for r in groups} >= {re.search(r"protein_uid=([0-9a-f]+)", c[8]).group(1) for c in cds}


def test_w25_output_is_byte_identical_across_threads_and_hash_seeds(world, tmp_path):
    """Two runs in separate processes, one with 1 DuckDB thread and PYTHONHASHSEED=0, the other with 3 threads and
    PYTHONHASHSEED=12345: every output file is byte-identical (AGENTS §4.3)."""
    root, loci = world
    outs = []
    for k, (threads, seed) in enumerate([("1", "0"), ("3", "12345")]):
        out = tmp_path / f"out{k}"
        env = {**os.environ, "PYTHONHASHSEED": seed}
        got = subprocess.run([sys.executable, "-m", "biobox.cli", "nbr", "flanks", anchors(loci), "--root", str(root),
                              "--release", REL, "--out", str(out), "--threads", threads, "--export", "gff",
                              "--level", "pfam"], env=env, capture_output=True, text=True)
        assert got.returncode == 0, got.stderr
        outs.append(out)
    a = sorted(p.relative_to(outs[0]) for p in outs[0].rglob("*") if p.is_file())
    b = sorted(p.relative_to(outs[1]) for p in outs[1].rglob("*") if p.is_file())
    assert a == b and len(a) > 10
    for p in a:
        assert (outs[0] / p).read_bytes() == (outs[1] / p).read_bytes(), p


@pytest.mark.parametrize("extra, needle", [
    (["--window", "5kb"], "gene window"),
    (["--cut", "1.5"], "--cut"),
    (["--level", "nr70"], "--level"),
])
def test_w25_refusals(world, tmp_path, capsys, extra, needle):
    """A bp window, a cut outside [0, 1] and an unknown level are refused (exit 2) and nothing is written."""
    out = tmp_path / "out"
    assert run_flanks(world, out, *extra) == 2
    assert needle in capsys.readouterr().err
    assert not out.exists()


def test_w25_a_non_empty_output_directory_is_refused(world, tmp_path, capsys):
    out = tmp_path / "out"
    out.mkdir()
    (out / "x").write_text("x")
    assert run_flanks(world, out) == 2
    assert "not empty" in capsys.readouterr().err

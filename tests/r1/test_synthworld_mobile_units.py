"""Ring R1: the pure parts of synthworld's phages, prophages and CRISPR arrays (synthworld_spec §2.1 S6, §3.3
truth_spacer_target, §7.3; card F-02 session 5). No ZOMBI2 here: the micro world is tested in
tests/r2/test_synthworld_mobile.py.

The oracles are written in this file: string slicing and reverse complements, the genetic code, a
brute-force Hamming scan over every window, the PAM rule read base by base, and the ancestry of a hand-made
tree.
"""

from __future__ import annotations

import collections
import copy

import numpy as np
import pytest

from biobox.synthworld import backbone, mobile

pytestmark = pytest.mark.covers("F-02")  # ring R1 (by directory)

CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
SENSE = sorted(c for c, a in CODE11.items() if a != "*")


def rc(s: str) -> str:
    return "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[c] for c in reversed(s))


def cut(seq: str, s: int, e: int) -> str:
    """Bases [s, e) of a circular record, e > len across the origin."""
    return (seq + seq)[s:e] if e > len(seq) else seq[s:e]


def protein(cds: str) -> str:
    return "".join(CODE11[cds[i:i + 3]] for i in range(0, len(cds), 3))


def rand_seq(rng, n: int) -> str:
    return "".join(rng.choice(list("ACGT"), size=n))


def rand_cds(rng, n_codons: int) -> str:
    return "ATG" + "".join(SENSE[int(i)] for i in rng.integers(len(SENSE), size=n_codons)) + "TAA"


# --- the lift from the S5 frame -------------------------------------------------------------------------------

def test_lift_moves_every_feature_with_its_bases():
    rng = np.random.default_rng(11)
    n_cross = 0
    for trial in range(30):
        N = int(rng.integers(50, 400))
        old = rand_seq(rng, N)
        pts = sorted({int(x) for x in rng.integers(0, N, size=int(rng.integers(1, 5)))})
        blocks = [(a, rand_seq(rng, int(rng.integers(1, 40)))) for a in pts]
        new, starts = mobile.insert_blocks(old, blocks)
        lift = mobile.Lift(N, [(a, len(b)) for a, b in blocks])
        assert len(new) == lift.new_length == N + sum(len(b) for _, b in blocks)
        for (a, b), s in zip(blocks, starts):
            assert new[s:s + len(b)] == b
        boundaries = {a if a else N for a in pts}
        n_checked = 0
        for _ in range(200):
            s = int(rng.integers(0, N))
            e = s + int(rng.integers(1, N))
            # straddling, written here: a point strictly inside the covered stretch, or N inside a wrap
            inside = any(s < a < e or s < a + N < e for a in boundaries) or (e > N and N in boundaries)
            if inside:
                with pytest.raises(ValueError):
                    lift.feature(s, e)
                continue
            s2, e2 = lift.feature(s, e)
            assert e2 - s2 == e - s and 0 <= s2 < len(new)
            assert cut(new, s2, e2) == cut(old, s, e), (trial, s, e, pts)
            n_checked += 1
            n_cross += e > N
        assert n_checked > 20
    assert n_cross > 20           # features across the origin were lifted too


def test_a_point_at_zero_counts_as_the_end():
    old = "ACGTACGTAC"
    new, starts = mobile.insert_blocks(old, [(0, "TTT")])
    assert new == old + "TTT" and starts == [10]
    lift = mobile.Lift(10, [(0, 3)])
    assert lift.feature(0, 4) == (0, 4) and lift.feature(6, 10) == (6, 10)
    with pytest.raises(ValueError):
        lift.feature(8, 12)            # across the origin, where the block sits


def test_flip_reverses_a_feature_with_its_block():
    block = "AAACCCGGGTTTACGT"
    s, e = 3, 9
    for strand in (1, -1):
        a, b, st = mobile.flip((s, e), strand, len(block), -1)
        assert st == -strand and rc(block)[a:b] == rc(block[s:e])
        assert mobile.flip((s, e), strand, len(block), 1) == (s, e, strand)


# --- the target scan ------------------------------------------------------------------------------------------

def _brute(seq: str, q: str, max_mm: int):
    n = len(q)
    out = set()
    for w in range(len(seq) - n + 1):
        win = seq[w:w + n]
        d = sum(x != y for x, y in zip(win, q))
        if d <= max_mm:
            out.add((w, 1, d))
        d = sum(x != y for x, y in zip(rc(win), q))
        if d <= max_mm:
            out.add((w, -1, d))
    return out


def test_find_targets_equals_a_brute_force_scan():
    rng = np.random.default_rng(5)
    seq = list(rand_seq(rng, 3000))
    queries = [rand_seq(rng, 30) for _ in range(24)]
    # plant copies with 0..7 mismatches on both strands
    for i, q in enumerate(queries[:16]):
        m = i % 8
        x = list(q)
        for p in rng.choice(30, size=m, replace=False):
            x[p] = "ACGT"[("ACGT".index(x[p]) + 1 + int(rng.integers(3))) % 4]
        x = "".join(x)
        if i % 2:
            x = rc(x)
        w = 100 + 170 * i
        seq[w:w + 30] = list(x)
    seq = "".join(seq)
    got = mobile.find_targets(seq, queries, 5, "NGG")
    by_q = collections.defaultdict(set)
    for qi, s, e, st, mm, pam in got:
        assert e - s == 30
        by_q[qi].add((s, st, mm))
        want_pam = seq[e:e + 3] if st == 1 else rc(seq[s - 3:s]) if s >= 3 else ""
        if st == 1 and e + 3 > len(seq):
            want_pam = ""
        assert pam == want_pam
    n_found = 0
    for i, q in enumerate(queries):
        want = _brute(seq, q, 5)
        assert by_q[i] == want, i
        n_found += len(want)
    assert n_found >= 12         # the planted copies with at most 5 mismatches


def test_scan_refuses_blocks_too_long_for_the_pigeonhole():
    ix = mobile.KmerIndex(mobile.to_codes("ACGT" * 20), 6)
    with pytest.raises(ValueError):
        mobile.scan(ix, mobile.to_codes("A" * 30)[None, :], 5)


# --- PAM and spacers ------------------------------------------------------------------------------------------

def test_phage_sites_follow_the_pam_rule_base_by_base():
    rng = np.random.default_rng(3)
    seqs = [rand_seq(rng, 900), rand_seq(rng, 400)]
    sites = mobile.PhageSites(seqs, 30, "NGG")
    for j, s in enumerate(seqs):
        L = len(s)
        d = sites.sites[j]
        assert sorted(set(d[(1, True)]) | set(d[(1, False)])) == list(range(0, L - 33 + 1))
        assert sorted(set(d[(-1, True)]) | set(d[(-1, False)])) == list(range(3, L - 30 + 1))
        plus_p = set(d[(1, True)].tolist())
        for w in range(0, L - 33 + 1):
            ok = s[w + 31] == "G" and s[w + 32] == "G"
            assert (w in plus_p) == ok
        plus_m = set(d[(-1, True)].tolist())
        for w in range(3, L - 30 + 1):
            pam = rc(s[w - 3:w])          # 3' of the protospacer on the minus strand
            assert (w in plus_m) == (pam[1:] == "GG")


def test_drawn_spacers_differ_from_their_protospacers_at_exactly_n_mismatches():
    rng = np.random.default_rng(8)
    seqs = [rand_seq(rng, 2000) for _ in range(3)]
    cfg = copy.deepcopy(backbone.load_preset("micro")["crispr"])
    sites = mobile.PhageSites(seqs, int(cfg["spacer_length"]), cfg["pam"])
    draw = np.random.default_rng(21)
    kinds = collections.Counter()
    mms = collections.Counter()
    for k in range(2000):
        sp = mobile.draw_spacer(cfg, sites, draw, f"x{k}", "X", "n1", 0.0)
        kinds[sp.source] += 1
        if sp.source == "decoy":
            assert sp.phage == -1 and len(sp.seq) == 30
            continue
        s = seqs[sp.phage]
        proto = s[sp.start:sp.end] if sp.strand == 1 else rc(s[sp.start:sp.end])
        pam = s[sp.end:sp.end + 3] if sp.strand == 1 else rc(s[sp.start - 3:sp.start])
        assert sp.end - sp.start == 30 and len(proto) == 30
        assert sum(a != b for a, b in zip(proto, sp.seq)) == sp.n_mismatches
        assert sp.pam_seq == pam and sp.pam_present == (pam[1:] == "GG")
        kinds["pam"] += sp.pam_present
        mms[sp.n_mismatches] += 1
    n = kinds["phage"]
    # 2,000 draws: decoys 0.05, PAM-less 0.1 of the rest, mismatches by the weights (4 sd bounds)
    assert abs(kinds["decoy"] / 2000 - cfg["decoy_fraction"]) <= 4 * (0.05 * 0.95 / 2000) ** 0.5
    assert abs(kinds["pam"] / n - (1 - cfg["no_pam_fraction"])) <= 4 * (0.09 / n) ** 0.5
    for m, w in enumerate(cfg["mismatch_weights"]):
        assert abs(mms[m] / n - w) <= 4 * (w * (1 - w) / n) ** 0.5, (m, mms)


# --- trees, arrays and prophages ------------------------------------------------------------------------------

def _tree():
    """A hand-made complete tree: stem n0; n1 (two tips and an extinct e5) and the tip n2."""
    br = [backbone.SpeciesBranch("n0", "", 1.0), backbone.SpeciesBranch("n1", "n0", 1.0),
          backbone.SpeciesBranch("n3", "n1", 1.0), backbone.SpeciesBranch("e5", "n1", 0.4),
          backbone.SpeciesBranch("n4", "n1", 1.0), backbone.SpeciesBranch("n2", "n0", 2.0)]
    return mobile.species_tree(br, {"n2", "n3", "n4"})


ANCESTORS = {"n3": {"n0", "n1", "n3"}, "n4": {"n0", "n1", "n4"}, "n2": {"n0", "n2"}}


def test_species_tree_clade_keeps_only_branches_with_extant_descendants():
    t = _tree()
    assert t.below["n1"] == ("n3", "n4") and t.below["e5"] == () and t.below["n0"] == ("n3", "n4", "n2")
    assert t.end["n3"] == t.end["n2"] == t.present == 3.0
    c = mobile.clade(t, "n1", 1.25)
    assert [n.lineage for n in c.nodes] == ["n1", "n3", "n4"] and c.tips == ["n3", "n4"]
    assert c.nodes[0].branch_length == pytest.approx(0.75) and c.nodes[0].start == 1.25
    assert c.nodes[1].parent == 0 and c.nodes[1].start == 2.0
    nodes = mobile.gene_tree(c, 900001, "CR00.cas1")
    assert [n.copy_label for n in nodes if n.extant] == ["n3_CR00.cas1", "n4_CR00.cas1"]
    with pytest.raises(ValueError):
        mobile.clade(t, "n1", 2.5)


def test_arrays_acquire_at_the_leader_end_and_delete_inside():
    t = _tree()
    cfg = copy.deepcopy(backbone.load_preset("micro")["crispr"])
    rng = np.random.default_rng(4)
    acquired_counts = []
    for rep in range(300):
        c = mobile.clade(t, "n0", 0.5)
        reg = []

        def acquire(node, time):
            reg.append((node, time))
            return len(reg) - 1

        tips = mobile.evolve_array(c, cfg, acquire, rng)
        assert set(tips) == {"n2", "n3", "n4"}
        order = {}
        for tip, ids in tips.items():
            assert cfg["min_spacers"] <= len(ids) <= cfg["max_spacers"]
            assert all(a > b for a, b in zip(ids, ids[1:]))        # index 0 is the newest
            for i in ids:
                assert reg[i][0] in ANCESTORS[tip]                  # acquired on the tip's own path
            times = [reg[i][1] for i in ids]
            assert all(a >= b for a, b in zip(times, times[1:]))
            order[tip] = ids
        for a in order:
            for b in order:
                shared = set(order[a]) & set(order[b])
                assert [x for x in order[a] if x in shared] == [x for x in order[b] if x in shared]
        # without deletions or a ceiling the acquisitions along a path are Poisson with the preset mean
        cfg2 = dict(cfg, deletions_per_spacer=0.0, max_spacers=10_000)
        reg = []
        tips2 = mobile.evolve_array(c, cfg2, acquire, rng)
        acquired_counts += [len(v) - cfg["initial_spacers"] for v in tips2.values()]
    mean = float(np.mean(acquired_counts))
    # 900 paths (three per run, correlated within a run: about 300 independent), sd sqrt(15 / 300) = 0.22
    assert abs(mean - cfg["acquisitions_mean"]) <= 4 * (cfg["acquisitions_mean"] / 300) ** 0.5, mean


def test_coding_sites_leave_out_start_stop_and_shared_codons():
    # gene A + [0, 30), gene B - [27, 51): they share bases 27-29 (A's stop, B's start on minus)
    igr, pos, strand = mobile.coding_sites(60, [(0, 30, 1), (27, 51, -1)])
    assert igr.tolist() == list(range(51, 60))
    a = sorted(p for p, s in zip(pos.tolist(), strand.tolist()) if s == 1)
    b = sorted(p for p, s in zip(pos.tolist(), strand.tolist()) if s == -1)
    assert a == [3, 6, 9, 12, 15, 18, 21, 24]            # codons 1..8 of A (its stop at 27 shared and left out)
    assert b == [30, 33, 36, 39, 42, 45]                  # internal codons of B: not its start (48) or stop (27)


def test_codons_read_and_written_on_either_strand():
    seq = np.frombuffer(b"ATGAAACCCGGGTTT", dtype=np.uint8).copy()
    pos = np.array([3, 6], dtype=np.int64)
    idx = mobile.read_codons(seq, pos, np.array([1, -1], dtype=np.int8))
    from biobox.synthworld import codons
    assert [codons.CODONS[i] for i in idx] == ["AAA", "GGG"]      # CCC read on minus is GGG
    mobile.write_codons(seq, pos, np.array([1, -1], dtype=np.int8), np.array([codons.CODONS.index("AAG"),
                                                                              codons.CODONS.index("GGA")]))
    assert seq.tobytes() == b"ATGAAGTCCGGGTTT"


def test_prophage_evolution_keeps_every_protein():
    rng = np.random.default_rng(9)
    a = rand_cds(rng, 60)                        # + [10, 10 + len)
    b = rand_cds(rng, 45)
    # A on plus, B on minus, 1,000 bases outside them
    seq = rand_seq(rng, 300) + a + rand_seq(rng, 400) + rc(b) + rand_seq(rng, 300)
    ga = (300, 300 + len(a), 1)
    sb = 300 + len(a) + 400
    gb = (sb, sb + len(b), -1)
    t = _tree()
    c = mobile.clade(t, "n0", 0.0)
    n = len(t.labels)
    index = {lab: i for i, lab in enumerate(t.labels)}
    from biobox.synthworld import codons
    usage = codons.cumulative(codons.usage_table(0.5, {"ATG": 1.0}, {"TAA": 1.0}))
    comp = np.cumsum(np.full(4, 0.25))
    tips = mobile.evolve_prophage(seq, [ga, gb], c, np.ones(n), np.stack([usage] * n), np.stack([comp] * n), index,
                                  1.0, np.random.default_rng(2))
    assert set(tips) == {"n2", "n3", "n4"}
    for tip, s in tips.items():
        assert len(s) == len(seq) and s != seq
        ca = s[ga[0]:ga[1]]
        cb = rc(s[gb[0]:gb[1]])
        assert protein(ca) == protein(a) and protein(cb) == protein(b)
        assert ca[:3] == a[:3] and ca[-3:] == a[-3:] and cb[:3] == b[:3] and cb[-3:] == b[-3:]
        assert ca != a or cb != b                 # synonymous changes happened at rate 1 over 3 time units
    # bases outside the genes change: expected fraction 3/4 (1 - exp(-r t)) with t = 3 on every path
    outside = [i for i in range(len(seq)) if not (ga[0] <= i < ga[1] or gb[0] <= i < gb[1])]
    frac = np.mean([np.mean([tips[x][i] != seq[i] for i in outside]) for x in tips])
    assert abs(frac - 0.75 * (1 - np.exp(-3.0))) < 0.06       # sd about 0.014 per tip (1,000 bases)


# --- presets --------------------------------------------------------------------------------------------------

def test_phage_preset_is_one_linear_chromosome_without_groups_or_motifs():
    p = backbone.load_preset("micro")
    q = mobile.phage_preset(p)
    assert "phages" not in q and "crispr" not in q
    assert q["genomes"]["topology"] == "linear" and q["species"]["n_extant"] == p["phages"]["count"] == 12
    assert q["species"]["burst_size"] == 1 and q["modules"]["count"] == 0 and q["hubs"]["count"] == 0
    assert q["plasmids"]["sizes"] == [] and q["vehicles"]["cargo"] == [] and q["decoys"]["partner_modules"] == []
    assert q["rna_motifs"]["n_models"] == 0 and q["rna_motifs"]["source"] == "synthetic"
    assert q["sequences"] == p["sequences"] and q["intergenic"] == p["intergenic"]
    assert p["genomes"]["topology"] == "circular"                # the host preset is untouched
    backbone.check_preset(q)
    assert mobile.phage_seed(p["default_seed"]) != p["default_seed"]
    assert mobile.phage_seed(p["default_seed"]) == mobile.phage_seed(p["default_seed"])


def test_micro_preset_has_the_section_5_1_counts_and_check_config_refuses_bad_ones():
    p = backbone.load_preset("micro")
    assert p["phages"]["count"] == 12 and p["prophages"]["count"] == 4 and p["crispr"]["arrays"] == 10
    assert (p["crispr"]["min_spacers"], p["crispr"]["max_spacers"]) == (5, 40)
    mobile.check_config(p["phages"], p["prophages"], p["crispr"])
    for section, key, value in (("prophages", "count", 13), ("crispr", "mismatch_weights", [0.2] * 7),
                                ("crispr", "initial_spacers", 50), ("crispr", "pam", "NGX"),
                                ("crispr", "clade_size", [3, 2])):
        q = copy.deepcopy(p)
        q[section][key] = value
        with pytest.raises(ValueError):
            mobile.check_config(q["phages"], q["prophages"], q["crispr"])

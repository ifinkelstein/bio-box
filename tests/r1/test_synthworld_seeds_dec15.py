"""Ring R1: synthworld seed derivation (synthworld_spec §6.1) and the DEC-15 preset guards (card F-02).

These need no ZOMBI2: the seeds are numpy and blake3, and the DEC-15 guards check a resolved preset
before anything is simulated. The worlds themselves are tested in tests/r2/test_synthworld_backbone.py.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest
from blake3 import blake3

from biobox.synthworld import backbone, seeds

pytestmark = pytest.mark.covers("F-02")  # ring R1 (by directory)


def _key(path: str) -> tuple[int, ...]:
    """§6.1 written out independently: blake3(path)[:16] read as four big-endian uint32 words."""
    d = blake3(path.encode("utf-8")).digest(length=16)
    return tuple(int.from_bytes(d[i:i + 4], "big") for i in (0, 4, 8, 12))


def test_spawn_key_is_blake3_of_the_component_path():
    assert seeds.spawn_key("S2/genomes") == _key("S2/genomes")
    assert seeds.spawn_key("S1/species_tree") != seeds.spawn_key("S2/genomes")


def test_generator_depends_only_on_seed_and_path():
    a = seeds.generator(7, "S1/clonal_burst").integers(0, 2**32, size=8)
    _ = seeds.generator(7, "S2/genomes").integers(0, 2**32, size=100)  # another stream drawn in between
    b = seeds.generator(7, "S1/clonal_burst").integers(0, 2**32, size=8)
    want = np.random.Generator(np.random.PCG64(np.random.SeedSequence(7, spawn_key=_key("S1/clonal_burst"))))
    assert a.tolist() == b.tolist() == want.integers(0, 2**32, size=8).tolist()
    assert a.tolist() != seeds.generator(8, "S1/clonal_burst").integers(0, 2**32, size=8).tolist()


def test_tool_seeds_are_nonzero_31_bit():
    got = {seeds.tool_seed(s, p) for s in range(200) for p in ("S1/species_tree", "S2/genomes")}
    assert len(got) == 400
    assert all(1 <= x <= 2**31 - 1 for x in got)


@pytest.mark.parametrize("bad", [-1, 2**64, 1.5, True])
def test_master_seed_is_u64(bad):
    with pytest.raises(ValueError):
        seeds.tool_seed(bad, "S2/genomes")


def test_micro_preset_meets_dec15():
    p = backbone.load_preset("micro")
    backbone.check_dec15(p)
    assert p["genomes"]["replacement"] == 0.0
    assert p["modules"]["fates"] == "group" and p["modules"]["scored"] is True
    assert p["root_scramble"]["transposition_passes"] > 0


def test_preset_overrides_merge_without_touching_the_rest():
    p = backbone.load_preset("micro", overrides={"modules": {"fates": "family", "scored": False}})
    q = backbone.load_preset("micro")
    assert p["modules"]["fates"] == "family" and p["modules"]["scored"] is False
    assert p["modules"]["count"] == q["modules"]["count"]
    assert p["genomes"] == q["genomes"]
    with pytest.raises(KeyError):
        backbone.load_preset("micro", overrides={"modules": {"no_such_key": 1}})


def test_dec15_a_refuses_a_scored_group_whose_fates_is_not_group():
    p = copy.deepcopy(backbone.load_preset("micro"))
    p["modules"]["fates"] = "family"           # still scored for adjacency
    with pytest.raises(backbone.Dec15Error, match=r"\(a\)"):
        backbone.check_dec15(p)
    p["modules"]["scored"] = False             # truth from the realized order only: allowed
    backbone.check_dec15(p)


def test_dec15_b_refuses_replacement_in_a_world_with_groups():
    p = copy.deepcopy(backbone.load_preset("micro"))
    p["genomes"]["replacement"] = 0.5
    with pytest.raises(backbone.Dec15Error, match=r"\(b\)"):
        backbone.check_dec15(p)
    p["modules"]["count"] = 0                  # a world with no groups at all may replace
    p["profile_only"]["count"] = 0
    backbone.check_dec15(p)


@pytest.mark.parametrize("key", ["transposition_passes", "stem_length"])
def test_dec15_c_refuses_an_unscrambled_root(key):
    p = copy.deepcopy(backbone.load_preset("micro"))
    p["root_scramble"][key] = 0
    with pytest.raises(backbone.Dec15Error, match=r"\(c\)"):
        backbone.check_dec15(p)


def test_dec15_c_scramble_check_reads_an_unscrambled_layout_as_unscrambled():
    # ZOMBI2's default root: every gene on +, declared families in one block at the end
    initial = [(f, 1) for f in range(100)] + [(f, 1) for f in range(100, 110)]
    stats = backbone.scramble_stats(initial, initial, declared=set(range(100, 110)), grouped_pairs=set())
    assert stats["plus_fraction"] == 1.0
    assert stats["initial_adjacency_kept"] == 1.0
    with pytest.raises(backbone.Dec15Error, match=r"\(c\)"):
        backbone.check_scrambled(stats)

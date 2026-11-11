"""Ring R1: which prefilter calls of clustering in parts search a shared index (card B-43), without MMseqs2.

A prefilter part may search one `createindex` index of its call's target instead of building the target's k-mer
table itself, but only at the call settings where that was proven to give the rebuilt result entry for entry
(MMseqs2 18.8cc5c). The choice is the hard-coded table `cluster.INDEX_SETTINGS`, checked here against the module
list: it names exactly the NR50 cascade's sensitive steps (`-s 2` and `-s 3`, with `--comp-bias-corr 1
--diag-score 1`), never the `-s 1 --comp-bias-corr 0 --diag-score 0` calls (the NR90 step and NR50 step 0, where
MMseqs2 aborts), and holds no setting that no call has. The index is built with the call's own settings, and a
part's per-part overhead in the cut is the index's load where it is shared, else the table's rebuild.
"""

from __future__ import annotations

import pytest

from biobox.l2 import cluster

pytestmark = pytest.mark.covers("MR-16")  # ring R1 (by directory)

SENSITIVE = ("cluster.nr50.15.prefilter", "cluster.nr50.19.prefilter")  # NR50 cascade steps 1 and 2
FAST = ("cluster.nr90.11.prefilter", "cluster.nr50.11.prefilter")  # the NR90 step and NR50 step 0


def prefilters() -> list[cluster.Call]:
    return [c for c in cluster.calls() if c.module == "prefilter"]


def test_b43_the_index_table_names_exactly_the_proven_calls_of_the_module_list():
    pf = prefilters()
    assert [c.name for c in pf] == [FAST[0], FAST[1], *SENSITIVE]
    assert sorted(c.name for c in cluster.calls() if c.shares_index) == sorted(SENSITIVE)
    by = {c.name: c for c in pf}
    assert {cluster.prefilter_settings(by[n]) for n in SENSITIVE} == {("2", "1", "1", "15"), ("3", "1", "1", "15")}
    assert {cluster.prefilter_settings(by[n]) for n in FAST} == {("1", "0", "0", "0")}
    # no stale row: every setting in the table is a call's
    assert cluster.INDEX_SETTINGS <= {cluster.prefilter_settings(c) for c in pf}
    assert not any(cluster.prefilter_settings(by[n]) in cluster.INDEX_SETTINGS for n in FAST)
    # only a prefilter can share an index
    assert not any(c.shares_index for c in cluster.calls() if c.module != "prefilter")


def test_b43_the_index_is_built_with_the_calls_own_settings():
    for c in prefilters():
        if not c.shares_index:
            continue
        argv = cluster.index_argv(c)
        assert argv[0] == "createindex"
        assert argv[1] == f"{cluster.part_dir(c)}/target" and argv[2] == f"{cluster.part_dir(c)}/index-tmp"
        flags = dict(zip(argv[3::2], argv[4::2]))
        for f in ("-s", "-k", "--comp-bias-corr", "--comp-bias-corr-scale", "--mask", "--mask-prob",
                  "--mask-lower-case", "--mask-n-repeat", "--spaced-kmer-mode", "--alph-size", "--seed-sub-mat",
                  "--max-seq-len", "--k-score", "--split", "--split-memory-limit"):
            assert flags[f] == c.argv[c.argv.index(f) + 1], f
        assert flags["--index-subset"] == "0"
        assert cluster.index_path(c) == f"{cluster.part_dir(c)}/target.idx"
    with pytest.raises(ValueError, match="shared index"):
        cluster.index_argv(next(c for c in prefilters() if not c.shares_index))


def test_b43_a_part_pays_the_index_load_where_it_is_shared_and_the_rebuild_elsewhere():
    targets = 1_000_000
    for c in prefilters():
        got = cluster.part_overhead_cpu_min(c, targets)
        rate = (cluster.PREFILTER_INDEX_LOAD_CPU_S_PER_TARGET if c.shares_index
                else cluster.PREFILTER_INDEX_CPU_S_PER_TARGET)
        assert got == pytest.approx(rate * targets / 60)
    assert cluster.PREFILTER_INDEX_LOAD_CPU_S_PER_TARGET < cluster.PREFILTER_INDEX_CPU_S_PER_TARGET / 4
    # the index itself is one process: its estimate grows with the target
    assert cluster.index_cpu_min(targets) == pytest.approx(cluster.CREATEINDEX_CPU_S_PER_TARGET * targets / 60)

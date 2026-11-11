"""Ring R1: the pieces of clustering in parts (card B-41) that need no MMseqs2.

- A step's record digests an MMseqs2 database over its type and its (key, entry) pairs in key order, not over
  its files: MMseqs2 writes one data file with one thread and `.0`, `.1`, ... with more, and its index need not
  be in key order. The same entries in either layout give one digest; one changed byte, or one key fewer, gives
  another. The databases here are written by hand in MMseqs2's format (data: NUL-terminated entries; index:
  `key<TAB>offset<TAB>length` lines, offsets into the data files taken as one run; dbtype: a 4-byte integer).
- The join of part databases holds every part's entries, in key order (the hand-written entries are the oracle).
- The cut of a per-query step: consecutive key ranges, a pure function of the keys, their sizes and the budget;
  every part's estimate fits the part budget (80% of the CPU-minute budget in build-stage); a smaller budget
  cuts more parts; a forced count cuts that many; a budget the per-part overhead alone fills is refused.
- A prefilter log without `Process prefiltering step 1 of 1` is refused: MMseqs2's own target split changes
  results (B-37 design §3, measured).
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from biobox.l2 import cluster

pytestmark = pytest.mark.covers("MR-16")  # ring R1 (by directory)

RESULT = (7).to_bytes(4, "little")  # an MMseqs2 result database's dbtype


def write_db(base: Path, entries: dict[int, bytes], files: int = 1, shuffle: int | None = None,
             dbtype: bytes = RESULT) -> None:
    """`entries` (without their NUL) as an MMseqs2 database at `base`: one data file, or `files` numbered ones
    each holding a consecutive share of the entries; the index in key order, or shuffled with seed `shuffle`."""
    keys = sorted(entries)
    chunks = [keys[i * len(keys) // files:(i + 1) * len(keys) // files] for i in range(files)]
    lines, off = [], 0
    for i, chunk in enumerate(chunks):
        data = b""
        for k in chunk:
            e = entries[k] + b"\0"
            lines.append(f"{k}\t{off + len(data)}\t{len(e)}\n")
            data += e
        Path(f"{base}.{i}" if files > 1 else base).write_bytes(data)
        off += len(data)
    if shuffle is not None:
        random.Random(shuffle).shuffle(lines)
    Path(f"{base}.index").write_text("".join(lines))
    Path(f"{base}.dbtype").write_bytes(dbtype)


def entries_of(n: int, seed: int) -> dict[int, bytes]:
    rng = random.Random(seed)
    return {k: bytes(rng.choice(b"ACDEFGHIKLMNPQRSTVWY\t\n0123456789") for _ in range(rng.randint(0, 40)))
            for k in rng.sample(range(10 * n), n)}


def test_b41_a_databases_digest_is_over_its_entries_not_its_files(tmp_path):
    e = entries_of(200, 1)
    write_db(tmp_path / "one", e)
    write_db(tmp_path / "two", e, files=2)
    write_db(tmp_path / "four", e, files=4, shuffle=7)
    d = cluster.db_digest(tmp_path / "one")
    assert cluster.db_digest(tmp_path / "two") == d == cluster.db_digest(tmp_path / "four")
    assert (tmp_path / "one").read_bytes() != (tmp_path / "two.0").read_bytes()  # the files differ
    assert list(cluster.db_entries(tmp_path / "four")) == [(k, e[k] + b"\0") for k in sorted(e)]
    k0 = sorted(e)[17]
    changed = e | {k0: e[k0] + b"x"}
    write_db(tmp_path / "changed", changed, files=2)
    fewer = {k: v for k, v in e.items() if k != k0}
    write_db(tmp_path / "fewer", fewer)
    write_db(tmp_path / "typed", e, dbtype=(0).to_bytes(4, "little"))
    others = {cluster.db_digest(tmp_path / n) for n in ("changed", "fewer", "typed")}
    assert d not in others and len(others) == 3


def test_b41_an_entry_that_is_not_nul_terminated_is_refused(tmp_path):
    write_db(tmp_path / "db", {1: b"ab", 2: b"cd"})
    (tmp_path / "db.index").write_text("1\t0\t3\n2\t3\t2\n")  # the second entry stops before its NUL
    with pytest.raises(ValueError, match="NUL"):
        cluster.db_digest(tmp_path / "db")


def test_b41_the_join_of_parts_holds_their_entries_in_key_order(tmp_path):
    e = entries_of(300, 2)
    keys = sorted(e)
    cuts = [keys[:50], keys[50:51], keys[51:]]
    parts = []
    for i, ks in enumerate(cuts):
        parts.append(tmp_path / f"part{i}")
        write_db(parts[-1], {k: e[k] for k in ks}, files=1 + i, shuffle=i)
    cluster.join_dbs(parts, tmp_path / "joined")
    assert list(cluster.db_entries(tmp_path / "joined")) == [(k, e[k] + b"\0") for k in keys]
    assert (tmp_path / "joined.dbtype").read_bytes() == RESULT
    write_db(tmp_path / "whole", e)
    assert cluster.db_digest(tmp_path / "joined") == cluster.db_digest(tmp_path / "whole")
    # parts whose keys are not in order, part after part, are refused
    with pytest.raises(ValueError, match="key"):
        cluster.join_dbs([parts[2], parts[0]], tmp_path / "bad")


def rows_of(n: int, seed: int) -> list[tuple[int, int]]:
    rng = random.Random(seed)
    return [(k, rng.randint(50, 600)) for k in sorted(rng.sample(range(5 * n), n))]


def check_cut(rows, parts) -> None:
    """Consecutive, non-empty key ranges that hold every key once, in order."""
    keys = [k for k, _ in rows]
    i = 0
    for p in parts:
        assert p["keys"] >= 1 and keys[i] == p["first"]
        i += p["keys"]
        assert keys[i - 1] == p["last"]
        assert p["units"] == sum(u for _, u in rows[i - p["keys"]:i])
    assert i == len(keys)


def test_b41_a_parts_estimate_over_80_percent_of_the_budget_is_cut_smaller():
    rows = rows_of(5_000, 3)
    budget = 300.0
    allowance = 0.8 * budget  # SHARD_MARGIN of build-stage
    one = cluster.cut(lambda: iter(rows), 100.0, part_cpu_min=allowance, overhead_cpu_min=5.0)
    assert len(one) == 1 and one[0]["cpu_min"] == pytest.approx(105.0)
    est = 1_000.0
    got = cluster.cut(lambda: iter(rows), est, part_cpu_min=allowance, overhead_cpu_min=5.0)
    check_cut(rows, got)
    assert len(got) == 5  # ceil(1000 / (240 - 5))
    assert all(p["cpu_min"] <= allowance for p in got)
    assert sum(p["cpu_min"] for p in got) == pytest.approx(est + 5.0 * len(got))
    smaller = cluster.cut(lambda: iter(rows), est, part_cpu_min=0.8 * 100.0, overhead_cpu_min=5.0)
    check_cut(rows, smaller)
    assert len(smaller) > len(got) and all(p["cpu_min"] <= 80.0 for p in smaller)
    # a pure function of the keys, their sizes and the budget
    assert cluster.cut(lambda: iter(list(rows)), est, part_cpu_min=allowance, overhead_cpu_min=5.0) == got
    # a budget that the per-part overhead alone fills cannot be met by more parts
    with pytest.raises(ValueError, match="overhead"):
        cluster.cut(lambda: iter(rows), est, part_cpu_min=5.0, overhead_cpu_min=5.0)


def test_b41_a_forced_count_cuts_that_many_parts_and_never_more_than_the_keys():
    rows = rows_of(1_000, 4)
    for n in (1, 2, 4, 7):
        got = cluster.cut(lambda: iter(rows), 1.0, parts=n)
        check_cut(rows, got)
        assert len(got) == n
        units = [p["units"] for p in got]
        assert max(units) - min(units) <= 2 * max(u for _, u in rows)  # equal work, to one key's size
    assert len(cluster.cut(lambda: iter(rows[:3]), 1.0, parts=8)) == 3
    zeros = [(k, 0) for k, _ in rows[:10]]  # empty entries: cut by count
    got = cluster.cut(lambda: iter(zeros), 1.0, parts=4)
    check_cut(zeros, got)
    assert [p["keys"] for p in got] == [3, 2, 3, 2]
    with pytest.raises(ValueError, match="order"):
        cluster.cut(lambda: iter([(2, 1), (1, 1)]), 1.0, parts=2)


@pytest.mark.parametrize("text, ok", [
    ("Query db start 1 to 10\nProcess prefiltering step 1 of 1\n\nk-mer similarity threshold: 154\n", True),
    ("Process prefiltering step 1 of 3\nProcess prefiltering step 2 of 3\nProcess prefiltering step 3 of 3\n", False),
    ("Index table: counting k-mers\n", False),
    ("", False),
])
def test_b41_a_prefilter_log_without_step_1_of_1_is_refused(text, ok):
    if ok:
        cluster.check_prefilter_log(text, "cluster.nr90.11.prefilter.part-00000")
    else:
        with pytest.raises(ValueError, match="step 1 of 1"):
            cluster.check_prefilter_log(text, "cluster.nr90.11.prefilter.part-00000")

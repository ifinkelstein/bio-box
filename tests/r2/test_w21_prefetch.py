"""Ring R2: a search reads the next block of its target while it searches the current one (card W-21).

The reader and the consumer are fakes with events, not clocks: the consumer of block 0 waits (a condition polled
under a timeout) until the reader has started on block 1, which can only happen if the two overlap. Nothing here
depends on wall time. The identical-bytes claim is tested against the unchanged W-19/W-07 tests, which run the
same searches through the same `Target`.

Part 2 of the card: `worker._families` reads `l1_domains` one row group at a time and keeps only the proteins asked
for; the oracle is a brute-force read of the whole table written here.
"""

from __future__ import annotations

import threading
import time
import weakref
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox.search import engine, worker

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-21")]

WAIT = 20.0  # seconds a condition is polled for before the test fails


def _wait(cond, what: str, timeout: float = WAIT) -> None:
    end = time.monotonic() + timeout
    while not cond():
        assert time.monotonic() < end, f"timed out waiting for: {what}"
        time.sleep(0.005)


class _Block:
    """A stand-in for a digitized block; counted while alive."""
    live = 0
    peak = 0
    lock = threading.Lock()

    def __init__(self, i: int):
        self.i = i
        with _Block.lock:
            _Block.live += 1
            _Block.peak = max(_Block.peak, _Block.live)
        weakref.finalize(self, _Block._gone)

    @staticmethod
    def _gone() -> None:
        with _Block.lock:
            _Block.live -= 1


@pytest.fixture(autouse=True)
def _reset_blocks():
    _Block.live = _Block.peak = 0
    yield


class _Reader:
    """A fake `Target._blocks`: records when each block's read starts and ends, and when the reader closed."""

    def __init__(self, n: int | None, fail_at: int | None = None):
        self.n, self.fail_at = n, fail_at
        self.started: list[int] = []
        self.ended: list[int] = []
        self.closed = threading.Event()

    def __call__(self):
        try:
            i = 0
            while self.n is None or i < self.n:
                self.started.append(i)
                if i == self.fail_at:
                    raise OSError(f"disk error reading block {i}")
                b = _Block(i)
                self.ended.append(i)
                yield b
                del b
                i += 1
        finally:
            self.closed.set()


def _target(reader) -> engine.Target:
    return engine.Target("protein", 0, reader)


def test_w21_the_next_block_is_read_while_the_current_one_is_searched():
    r = _Reader(4)
    seen = []
    t = _target(r)
    for block in t.blocks():
        if block.i == 0:
            # today the reader is a generator that runs only when the consumer asks: block 1 is never started
            _wait(lambda: 1 in r.started, "the reader to start block 1 while block 0 is being searched")
        seen.append(block.i)
        del block
    assert seen == [0, 1, 2, 3] and t.blocks_read == 4


def test_w21_no_more_than_two_blocks_are_alive_and_the_reader_stays_one_ahead():
    r = _Reader(6)
    t = _target(r)
    for block in t.blocks():
        i = block.i
        _wait(lambda: len(r.ended) >= min(i + 2, 6), f"block {i + 1} to be read ahead")
        end = time.monotonic() + 0.1  # a negative cannot be waited for: watch a reader that is not bounded run on
        while time.monotonic() < end:
            assert len(r.started) <= i + 2, f"the reader ran {len(r.started) - i - 1} blocks ahead of block {i}"
            time.sleep(0.005)
        del block
    assert _Block.peak <= 2, f"{_Block.peak} blocks were alive at once; the bound is two"


def test_w21_a_reader_error_surfaces_in_the_caller_and_does_not_hang_it():
    r = _Reader(5, fail_at=2)
    got = []
    with pytest.raises(OSError, match="disk error reading block 2"):
        for block in _target(r).blocks():
            got.append(block.i)
            del block
    assert got == [0, 1]
    _wait(r.closed.is_set, "the reader to close")


def test_w21_a_consumer_that_stops_stops_the_reader():
    r = _Reader(None)  # a target that never ends
    gen = _target(r).blocks()
    first = next(gen)
    assert first.i == 0
    _wait(lambda: 1 in r.started, "the reader to run ahead")
    gen.close()  # a cancelled search stops iterating
    _wait(r.closed.is_set, "the reader to stop after the consumer closed")
    _wait(lambda: not [t for t in threading.enumerate() if t.name == "biobox-prefetch"], "the reader thread to exit")


def test_w21_a_failed_search_stops_the_reader():
    """`engine.search` on a target whose blocks the search cannot use: the error reaches the caller and the reader,
    which would otherwise read a target that never ends, is closed."""
    r = _Reader(None)
    hmm = engine.build_from_sequence("w21", "MKVLAAGIVGLLLAQPAFAEEKTHRPSLVDSY")
    with pytest.raises(Exception):
        engine.search(hmm, _target(r), 10.0, 5.0, threads=1)  # a fake block is not a DigitalSequenceBlock
    _wait(r.closed.is_set, "the reader to close after the search failed")
    _wait(lambda: not [t for t in threading.enumerate() if t.name == "biobox-prefetch"], "the reader thread to exit")


def test_w21_an_empty_target_and_a_target_read_twice_give_each_pass_its_own_reader():
    t = _target(_Reader(0))
    assert list(t.blocks()) == [] and t.blocks_read == 0
    r = _Reader(3)
    t = _target(r)
    assert [b.i for b in t.blocks()] == [0, 1, 2]
    assert [b.i for b in t.blocks()] == [0, 1, 2] and t.blocks_read == 3


# ---- part 2: worker._families streams l1_domains ------------------------------------------------------------------


class _Files:
    """The one method of `ReleaseFiles` that `_families` uses."""

    def __init__(self, paths: list[Path]):
        self.paths = paths

    def table(self, name: str, **_):
        assert name == "l1_domains"
        return self.paths


def _domains(tmp_path: Path, n_files=2, groups=5, per_group=40):
    """An `l1_domains`-shaped table in several files of several row groups; proteins sorted, 3 models cycling."""
    rows = []
    paths = []
    k = 0
    for f in range(n_files):
        tabs = []
        for g in range(groups):
            uid = [(k + j).to_bytes(16, "big") for j in range(per_group)]
            acc = [f"PF{(k + j) % 3:05d}" for j in range(per_group)]
            rows += list(zip(uid, acc))
            tabs.append(pa.table({"protein_uid": pa.array(uid, pa.binary(16)), "model_acc": acc}))
            k += per_group
        p = tmp_path / f"l1_domains_{f}.parquet"
        with pq.ParquetWriter(p, tabs[0].schema) as w:
            for t in tabs:
                w.write_table(t)  # one row group per table
        paths.append(p)
    return _Files(paths), rows, groups * per_group


def test_w21_families_reads_one_row_group_at_a_time_and_keeps_only_the_proteins_asked_for(tmp_path, monkeypatch):
    files, rows, group_rows = _domains(tmp_path)
    want_uids = [rows[i][0] for i in (3, 77, 150, 399)] + [b"\xff" * 16]  # the last one has no domain
    oracle: dict[bytes, set[str]] = {}
    for u, a in rows:  # brute force, from the rows written, not from the function under test
        if u in want_uids:
            oracle.setdefault(u, set()).add(a)

    def boom(*a, **k):
        raise AssertionError("pq.read_table: the whole of l1_domains is read at once")

    monkeypatch.setattr(pq, "read_table", boom)
    reads: list[int] = []
    orig = pq.ParquetFile.read_row_group

    def spy(self, i, *a, **k):
        t = orig(self, i, *a, **k)
        reads.append(t.num_rows)
        return t

    monkeypatch.setattr(pq.ParquetFile, "read_row_group", spy)
    fam = worker._Families(files)
    got = fam.of(want_uids)
    assert got == oracle
    assert reads and max(reads) == 40 < len(rows), reads  # a row group at a time, never the 400 rows
    assert fam.rows_materialized <= 8, fam.rows_materialized  # only matching rows became Python objects
    assert fam.of([]) == {} and fam.of([b"\xff" * 16]) == {}

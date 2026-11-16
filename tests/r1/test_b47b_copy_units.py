"""Ring R1, card B-47b (DEC-35 (2)): the parts of a domain scan that copies measured rows, without a scan.

- `domains._assign`: which of several sorted uid streams first holds each uid of a sorted stream, by a merge that
  holds one uid per stream. The oracle is a loop over Python sets written here; streams of random 2-byte uids
  (many collisions) with a fixed seed, empty streams, and a stream that holds every uid.
- `domains._Rows.take`: the rows of a sorted domain file for uids asked in increasing order, across row groups (a
  protein's rows that straddle a row-group boundary, row groups with no asked uid), against a filter written here.
- `domains._Rows.until` and `rest`: the rows before a uid and the rest, which together give back every row once.
- The CLI takes `--seq-stages` for `stage domains plan` only: the plan records what is copied, so the other actions
  refuse it (exit 2) before reading anything.
"""

from __future__ import annotations

import random

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import cli
from biobox import io as bio
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("MR-12", "MR-08")  # ring R1 by its directory (tests/conftest.py)


def _domains():
    from biobox.l1 import domains  # imported here so that a missing name fails the tests, not the collection

    return domains


def _brute(mine: list[bytes], sources: list[list[bytes]]) -> list[int]:
    sets = [set(s) for s in sources]
    out = []
    for u in mine:
        hit = [k for k, s in enumerate(sets) if u in s]
        out.append(hit[0] if hit else -1)
    return out


def _uids(rng: random.Random, n: int) -> list[bytes]:
    return sorted({rng.randrange(0, 1500).to_bytes(2, "big") for _ in range(n)})


@pytest.mark.parametrize("seed", range(6))
def test_assign_gives_the_first_stream_that_holds_each_uid(seed):
    rng = random.Random(seed)
    mine = _uids(rng, 800)
    sources = [_uids(rng, rng.choice([0, 50, 400, 1200])) for _ in range(rng.randrange(1, 5))]
    got = [(u, n, k) for u, n, k in _domains()._assign(((u, len(u) + i) for i, u in enumerate(mine)),
                                                         [iter(s) for s in sources])]
    assert [u for u, _, _ in got] == mine
    assert [n for _, n, _ in got] == [len(u) + i for i, u in enumerate(mine)]
    assert [k for _, _, k in got] == _brute(mine, sources)


def test_assign_with_no_stream_or_a_stream_holding_everything():
    d = _domains()
    mine = [bytes([i]) for i in range(0, 200, 3)]
    assert [k for _, _, k in d._assign(((u, 1) for u in mine), [])] == [-1] * len(mine)
    every = [bytes([i]) for i in range(256)]
    got = [k for _, _, k in d._assign(((u, 1) for u in mine), [iter([]), iter(every), iter(every)])]
    assert got == [1] * len(mine)


def _domain_file(path, uids_rows: list[tuple[bytes, int]], row_group_rows: int) -> list[dict]:
    """A sorted l1_domains file: `n` rows for each (uid, n), row groups of `row_group_rows` rows. Returns the rows."""
    rows = []
    for u, n in uids_rows:
        for j in range(n):
            rows.append({"protein_uid": u, "model_acc": f"M{j:03d}", "model_name": f"m{j}", "model_len": 100,
                         "seq_bitscore": 30.0, "seq_evalue_zref": 1e-5, "dom_bitscore": 30.0, "dom_ievalue_zref": 1e-5,
                         "env_start": j, "env_end": j + 50, "ali_start": j, "ali_end": j + 40, "hmm_from": 0,
                         "hmm_to": 90, "cutoff": "gathering", "z": 1e8, "dom_z": 1e8})
    spec = table("l1_domains")
    bio.write_table(pa.Table.from_pylist(rows, schema=spec.schema), spec, path, row_group_rows=row_group_rows)
    return rows


def _uid16(i: int) -> bytes:
    return i.to_bytes(16, "big")


def test_rows_take_finds_every_row_of_each_asked_uid_across_row_groups(tmp_path):
    rng = random.Random(7)
    present = sorted(rng.sample(range(10_000), 300))
    spec_rows = [(_uid16(i), rng.choice([1, 1, 2, 3, 7])) for i in present]
    rows = _domain_file(tmp_path / "part-0.parquet", spec_rows[:150], row_group_rows=5)
    rows += _domain_file(tmp_path / "part-1.parquet", spec_rows[150:], row_group_rows=4)
    assert pq.ParquetFile(tmp_path / "part-0.parquet").metadata.num_row_groups > 20
    asked = sorted(set(rng.sample(present, 90)) | set(rng.sample(range(10_000), 60)))
    cur = _domains()._Rows([tmp_path / "part-0.parquet", tmp_path / "part-1.parquet"], None)
    for i in asked:
        got = [r for t in cur.take(_uid16(i)) for r in t.to_pylist()]
        want = [r for r in rows if r["protein_uid"] == _uid16(i)]
        assert got == want, i


def test_rows_until_and_rest_give_back_every_row_once_in_order(tmp_path):
    rng = random.Random(11)
    spec_rows = [(_uid16(i), rng.choice([1, 2, 5])) for i in sorted(rng.sample(range(5_000), 200))]
    rows = _domain_file(tmp_path / "part-0.parquet", spec_rows, row_group_rows=6)
    cur = _domains()._Rows([tmp_path / "part-0.parquet"], None)
    cuts = sorted(rng.sample(range(5_000), 40))
    got = []
    for c in cuts:
        before = [r for t in cur.until(_uid16(c)) for r in t.to_pylist()]
        assert all(r["protein_uid"] < _uid16(c) for r in before)
        got += before
        assert len(got) == sum(1 for r in rows if r["protein_uid"] < _uid16(c))
    got += [r for t in cur.rest() for r in t.to_pylist()]
    assert got == rows


@pytest.mark.parametrize("action", [None, "shard", "merge"])
def test_seq_stages_is_refused_by_every_domains_action_but_plan(tmp_path, capsys, action):
    argv = ["stage", "domains", *([] if action is None else [action]), "--root", str(tmp_path / "root"),
            "--stage", "B", "--library", "t", "--library-version", "1", "--hmm", str(tmp_path / "x.hmm"),
            "--cutoff", "gathering", "--shard", "0", "--tmp", str(tmp_path / "scratch"), "--seq-stages", "B,A"]
    assert cli.main(argv) == 2
    err = capsys.readouterr().err
    assert "--seq-stages" in err and "plan" in err
    assert not (tmp_path / "root").exists()

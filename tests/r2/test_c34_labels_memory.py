"""Ring R2: `bio-box stage labels` in bounded memory (card C-34; AGENTS §4.1, §4.2), in the form of card E-12's test.

Two stages, each labelled by `store.write_labels` in a process of its own, launched through `runstep.measured`, its
peak memory read with `runstep.step_peak` (card C-15):

- the 100x synthetic store of `scripts/stress/synth_scale.py`, run unchanged: its stage S0 was ingested with label
  files, so a view root holding only S0's `contig_occurrence` and ingest ledger (symlinks) stands for the stage before
  card B-27, and the step's label file must equal the one ingest wrote there, by sha256;
- a constructed stage of 1,000,000 records over 600,000 contigs in 64 source objects (30 times the 32,815 contigs of
  stages S1, P1 and P2 together, *derived* from card W-36's count), written through `io.write_table`, with an ingest
  ledger that counts them; its label rows are checked against a separate loop over the same rows.

The bound. The step holds no table in Python: DuckDB sorts the stage's records under `store.DUCK_MEMORY_LIMIT` (1 GB,
spilling to scratch beyond it) and the label rows come one batch of `io.ROW_GROUP_ROWS` at a time into the writer. So
the peak is DuckDB's limit plus the interpreter, pyarrow and one batch: `PEAK_BOUND` = 1 GiB + 256 MiB = 1.25 GiB
(derived), for both stages.
"""

from __future__ import annotations

import contextlib
import hashlib
import io as stdio
import os
import runpy
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import io, ledger, runstep
from biobox.registry.schemas import table

pytestmark = [pytest.mark.r2, pytest.mark.covers("C-34")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
COPIES = 100
PEAK_BOUND = (1 << 30) + (256 << 20)
N_RECORDS, N_CONTIGS, N_OBJECTS = 1_000_000, 600_000, 64

CHILD = r"""
import json, sys
from pathlib import Path
from biobox.l0seq import store
res = store.write_labels(Path(sys.argv[1]), "S0", out=Path(sys.argv[2]), tmp=Path(sys.argv[3]))
print(json.dumps({k: res[k] for k in ("records", "contigs", "objects")}))
"""


def _synth(out: Path) -> Path:
    argv = sys.argv
    sys.argv = [str(SCRIPT), str(REPO), str(out), str(COPIES), "--shard-bases", str(2**26), "--seed", "1"]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()):
            g = runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
    g["con"].close()
    return Path(g["root"])


def _view(src: Path, dst: Path) -> Path:
    """A root holding only stage S0's `contig_occurrence` and ingest ledger of `src`, as symlinks."""
    for rel in ("l0seq/contig_occurrence/stage=S0", "ledger/stage=S0/step=ingest"):
        (dst / rel).parent.mkdir(parents=True, exist_ok=True)
        os.symlink(src / rel, dst / rel)
    return dst


def _constructed(root: Path) -> list[tuple]:
    """Stage S0 of N_RECORDS records over N_CONTIGS contigs in N_OBJECTS objects, written as `contig_occurrence` (one
    source) with its ingest ledger; returns the label oracle: per contig, (contig_uid, first accession, records)."""
    rng = np.random.default_rng(34)
    shas = sorted(rng.bytes(32) for _ in range(N_OBJECTS))
    contig = rng.integers(0, N_CONTIGS, N_RECORDS)
    contig[:N_CONTIGS] = np.arange(N_CONTIGS)  # every contig has a record
    rng.shuffle(contig)
    uids = [hashlib.blake2b(int(i).to_bytes(4, "big"), digest_size=16).digest() for i in range(N_CONTIGS)]
    per = np.array_split(np.arange(N_RECORDS), N_OBJECTS)
    spec = table("contig_occurrence")

    def batches():
        for k, idx in enumerate(per):
            yield pa.table({
                "source_uid": pa.array(np.full(len(idx), 7, np.uint64)),
                "source_release": pa.array(["1"] * len(idx)),
                "l0_object_sha256": pa.array([shas[k]] * len(idx), pa.binary(32)),
                "record_ordinal": pa.array(np.arange(len(idx), dtype=np.uint64)),
                "accession": pa.array([f"acc{int(i)}" for i in idx]),
                "original_header": pa.array([f"acc{int(i)} contig {int(contig[i])}" for i in idx]),
                "raw_digest": pa.array([hashlib.sha256(int(contig[i]).to_bytes(4, "big")).digest() for i in idx],
                                       pa.binary(32)),
                "contig_uid": pa.array([uids[contig[i]] for i in idx], pa.binary(16)),
                "line_width": pa.array(np.full(len(idx), 80, np.uint32)),
            }, schema=spec.schema)

    io.write_table(batches(), spec, root / "l0seq/contig_occurrence/stage=S0/source=syn/part-0.parquet")
    ledger.write(root, "S0", "ingest", (ledger.row(shas[k], "l0_object", len(idx), len(idx))
                                        for k, idx in enumerate(per)))
    first: dict[int, int] = {}
    n = np.bincount(contig, minlength=N_CONTIGS)
    for k, idx in enumerate(per):  # objects in sha order, records in ordinal order: the first is the first seen
        for i in idx:
            first.setdefault(int(contig[i]), int(i))
    return sorted((uids[c], f"acc{first[c]}", int(n[c])) for c in range(N_CONTIGS))


def _measured(root: Path, tmp: Path, name: str) -> tuple[int, str, Path]:
    peak, out, scratch = tmp / f"peak-{name}.json", tmp / f"out-{name}", tmp / f"scratch-{name}"
    scratch.mkdir()
    cmd = runstep.measured([sys.executable, "-c", CHILD, str(root), str(out), str(scratch)], peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    return rss, source, out


def test_c34_labels_of_the_100x_store_in_bounded_memory(tmp_path_factory, tmp_path, record_property):
    store_root = _synth(tmp_path_factory.mktemp("c34x100"))
    view = _view(store_root, tmp_path / "view")
    rss, source, out = _measured(view, tmp_path, "x100")
    record_property("c34_x100_peak_bytes", rss)
    print(f"C34 100x: peak {rss / 2**20:.1f} MiB ({source})")
    got = io.part_files(out / "l0seq/contig_label/stage=S0")
    want = io.part_files(store_root / "l0seq/contig_label/stage=S0")
    assert want and [p.name for p in got] == [p.name for p in want]
    assert [p.read_bytes() for p in got] == [p.read_bytes() for p in want]
    assert rss <= PEAK_BOUND, (rss, source)


def test_c34_labels_of_a_million_records_in_bounded_memory(tmp_path, record_property):
    root = tmp_path / "root"
    oracle = _constructed(root)
    rss, source, out = _measured(root, tmp_path, "m1")
    record_property("c34_m1_peak_bytes", rss)
    print(f"C34 {N_RECORDS} records: peak {rss / 2**20:.1f} MiB ({source})")
    files = io.part_files(out / "l0seq/contig_label/stage=S0")
    got = []
    for f in files:
        md = pq.read_metadata(f)
        assert all(md.row_group(i).num_rows <= io.ROW_GROUP_ROWS for i in range(md.num_row_groups))
        t = pq.read_table(f, columns=["contig_uid", "accession", "n_records"], partitioning=None)
        got += list(zip(*(t.column(c).to_pylist() for c in ("contig_uid", "accession", "n_records"))))
    assert got == oracle
    assert rss <= PEAK_BOUND, (rss, source)

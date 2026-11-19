"""Ring R2, card C-38f (§3.8 test (7)): round plan and merge stay in bounded memory on the 100x synthetic store
(AGENTS §4.1, §4.2), modeled on `tests/r2/test_b47b_copy_memory.py`.

`scripts/stress/synth_scale.py` is run at 10x and 100x copies of the micro-S0 fixture. Stage S0 is prepared as a
completed annotated source for `synfam-1`; stage T holds the same proteins at depth `called`, and release `rT.1`
lists only `T` (so every protein of `T` is `not scanned` in `rT.1`, while `S0` on disk already measured them).
`domains.plan_round` and `domains.merge_round` run in a child process measured by `runstep.measured`: every protein
of `T` is copied from `S0`, `round_domains` equals `S0`'s `l1_domains` byte for byte, `round_coverage` holds one
row per protein of `T` with exact `n_rows`, and peak RSS stays within `PEAK_BOUND` (512 MiB) and grows by at most
`GROWTH_BOUND` (64 MiB) from 10x to 100x.
"""

from __future__ import annotations

import contextlib
import hashlib
import io as stdio
import json
import runpy
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, ids, ledger, runstep
from biobox import io as bio
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.registry.schemas import table

pytestmark = [pytest.mark.r2, pytest.mark.covers("C-38f", "MR-08")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
SMALL, BIG = 10, 100
GROWTH_BOUND = 64 * 2**20
PEAK_BOUND = 512 * 2**20
LIB, VER = "synfam", "1"
STEP = f"domains-{LIB}-{VER}"

CHILD = r"""
import json, sys
from pathlib import Path
from biobox.l1 import domains
root, hmm, req, out, tmp = (Path(x) for x in sys.argv[1:6])
p = domains.plan_round(root, "rT.1", sys.argv[6], "synfam", "1", hmm, "gathering", req, 10**12,
                       seq_stages=["S0"], tmp=tmp, out=out)
m = domains.merge_round(root, p["round_id"], "synfam", "1", tmp=tmp, out=out)
print(json.dumps({"plan": p, "merge": m}, sort_keys=True))
"""


def _store(out: Path, copies: int) -> Path:
    argv = sys.argv
    sys.argv = [str(SCRIPT), str(REPO), str(out), str(copies), "--shard-bases", str(2**26), "--seed", "1"]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()):
            g = runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
    g["con"].close()
    return Path(g["root"])


def _protein_seq(root: Path, stage: str) -> list[bytes]:
    fp = bio.part_files(root / f"l1/feature_protein/gene_call_version={V}/stage=S0")
    uids = sorted({u for f in fp for u in pq.read_table(f, columns=["protein_uid"], partitioning=None)
                   .column(0).to_pylist()})
    seqs = ["M" * (60 + u[0] % 200) for u in uids]
    t = pa.table({"protein_uid": uids, "protein_digest": [hashlib.sha256(s.encode()).digest() for s in seqs],
                  "length": [len(s) for s in seqs], "seq": seqs})
    bio.write_table(t, table("protein_seq"),
                    root / f"l0seq/protein_seq/gene_call_version={V}/stage={stage}/protein_prefix=all/part-0.parquet")
    return uids


def _as_completed_source(root: Path, hmm: Path) -> None:
    from biobox.l1 import domains

    side = root.parent / "plan-side"
    seq = f"l0seq/protein_seq/gene_call_version={V}/stage=S0"
    (side / seq).parent.mkdir(parents=True)
    (side / seq).symlink_to(root / seq)
    domains.plan(side, "S0", V, LIB, VER, hmm, "gathering", 10**12)
    rel = domains.scan_dir("S0", LIB, VER) + "/plan.json"
    (root / rel).parent.mkdir(parents=True)
    (root / rel).write_bytes((side / rel).read_bytes())
    files = bio.part_files(root / f"l1/domains/library={LIB}-{VER}/stage=S0")
    n_in = pq.read_metadata(bio.part_files(root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0/protein_prefix=all")[0]).num_rows
    hit = len({u for f in files for u in pq.read_table(f, columns=["protein_uid"], partitioning=None).column(0).to_pylist()})
    tools = domains._tools()
    lib = f"{LIB}-{VER}".encode()
    ledger.write(root, "S0", STEP, [
        ledger.row(lib, "library", n_in, hit, n_rejected=n_in - hit, reason="no_hit", tools=tools),
        ledger.row(lib, "library", n_in, hit, tools=tools),
    ])
    led = bio.part_files(root / f"ledger/stage=S0/step={STEP}")
    outputs = [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": ids.sha256_file(p).hex()}
               for p in files + led]
    rec = root / f"runs/stage=S0/{STEP}.json"
    rec.parent.mkdir(parents=True, exist_ok=True)
    rec.write_text(json.dumps({"stage": "S0", "step": STEP, "outputs": sorted(outputs, key=lambda o: o["path"])},
                              sort_keys=True, indent=1) + "\n")


def _measured(root: Path, hmm: Path, req: Path, tmp: Path) -> tuple[int, str, dict]:
    name = root.parent.name
    peak, out, scratch = tmp / f"peak-{name}.json", tmp / f"out-{name}", tmp / f"scratch-{name}"
    scratch.mkdir()
    cmd = runstep.measured([sys.executable, "-c", CHILD, str(root), str(hmm), str(req), str(out), str(scratch), V],
                           peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    res = json.loads(p.stdout.strip().splitlines()[-1])
    rid = res["plan"]["round_id"]
    src = bio.part_files(root / f"l1/domains/library={LIB}-{VER}/stage=S0")
    got = bio.part_files(out / f"l1/round_domains/round={rid}/library={LIB}-{VER}")
    assert [f.read_bytes() for f in got] == [f.read_bytes() for f in src]
    cov = bio.part_files(out / f"l1/round_coverage/round={rid}/library={LIB}-{VER}")
    assert sum(pq.read_metadata(f).num_rows for f in cov) == res["plan"]["proteins"]
    return rss, source, res


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    out = {}
    hmm = tmp_path_factory.mktemp("c38f-hmm") / "stand-in.hmm"
    hmm.write_bytes(b"a stand-in for an HMM file: only its sha256 is read by a plan\n")
    for n in (SMALL, BIG):
        base = tmp_path_factory.mktemp(f"c38fx{n}")
        root = _store(base, n)
        _protein_seq(root, "S0")
        uids = _protein_seq(root, "T")
        _as_completed_source(root, hmm)
        catalog.define_release(root, "rT.1", ["T"])
        req = base / "req.txt"
        req.write_text("\n".join(u.hex() for u in uids) + "\n", encoding="utf-8")
        out[n] = (root, req)
    return out, hmm


def test_c38f_round_plan_and_merge_peak_memory_bounded_100x(stores, tmp_path, record_property):
    roots, hmm = stores
    got = {n: _measured(root, hmm, req, tmp_path) for n, (root, req) in roots.items()}
    record_property("c38f_peak_bytes", {n: g[0] for n, g in got.items()})
    record_property("c38f_runs", {n: g[2] for n, g in got.items()})
    print(f"C-38f round peaks {{copies: (bytes, source, result)}}: { {n: g for n, g in got.items()} }")
    for n, (_, _, res) in got.items():
        assert res["plan"]["shards"] == 0 and res["plan"]["copied_proteins"] == res["plan"]["proteins"] > 0, n
    assert got[BIG][2]["plan"]["proteins"] >= 9 * got[SMALL][2]["plan"]["proteins"]
    growth = got[BIG][0] - got[SMALL][0]
    assert growth <= GROWTH_BOUND, (f"peak RSS {got[SMALL][0] / 2**20:.1f} MiB at {SMALL}x, {got[BIG][0] / 2**20:.1f} "
                                    f"MiB at {BIG}x: grew {growth / 2**20:.1f} MiB, bound {GROWTH_BOUND / 2**20:.0f}")
    assert max(g[0] for g in got.values()) <= PEAK_BOUND, got

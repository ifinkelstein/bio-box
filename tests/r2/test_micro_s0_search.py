"""Ring R2: the search pack, runner and exports on micro-S0, with fault injection (C-12, S-04..S-07, S-13).

Micro-S0 is built as in test_micro_s0_l3 (plus IGRs), then the release's pack. Runs use the
real worker under the real limits; faults are injected through the runner's test hook
(`fault=`), and on Linux the OOM test allocates past a real RLIMIT_AS. S-07 runs on a
synthetic world with a planted trap: a seed family A, a library family B that A's seed does
not carry, and bridge sequences between them, so an unguarded jackhmmer drifts into B.
"""

from __future__ import annotations

import hashlib
import io
import json
import platform
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer
import pytest

from biobox import catalog, cli, contracts, ids
from biobox.l1 import domains, igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l3 import tapes
from biobox.registry.schemas import conformance_problems, table
from biobox.search import engine, export, fastz, nbr, orfs, pack, runner
from biobox.search import recipe as rcp
from biobox.search.recipe import RecipeError

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_l3 as micro  # noqa: E402

pytestmark = pytest.mark.r2
REL = micro.REL
REPO = Path(__file__).resolve().parents[2]
LIMITS = runner.Limits(cpus=2, wall_s=120)


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("search")
    root = d / "a"
    circ = micro.build(root)
    igr.run(root, "S0", circ)
    man = pack.build(root, REL, models={"test-1": d / "a.hmm"})
    return root, d, man


def _seed(root: Path) -> tuple[bytes, bytes]:
    """A library-hit protein and an unrelated protein (no domain hit, not similar)."""
    dom = pq.read_table(next((root / "l1/domains").rglob("*.parquet"))).to_pylist()
    hit = sorted({r["protein_uid"] for r in dom})
    others = sorted(set(domains.load_proteins(root, V, "S0")) - set(hit))
    return hit[0], others[len(others) // 2]


def _recipe(seed_uid: bytes, positive: bytes, tool="phmmer", mx=1, extra="", target="unique", bits=25) -> str:
    return (f"recipe: micro_{tool}\nwritten_against: {REL}\nseed: [{{protein_uid: {seed_uid.hex()}}}]\n"
            f"target: {target}\ntool: {tool}\niterations: {{max: {mx}, inclusion_bits: {bits}}}\n"
            f"controls: {{positive: [{{protein_uid: {positive.hex()}}}]}}\n{extra}")


def _run(world, text, **kw):
    root, d, _ = world
    return runner.run(root, text, REL, work_base=d / "work", limits=kw.pop("limits", LIMITS), **kw)


# --- C-12 ------------------------------------------------------------------------------------------


def test_c12_pack_verifies_and_counts_equal_sources(world):
    root, d, man = world
    assert pack.verify(root, REL, bytes.fromhex(man["pack_id"])) == []
    t = man["targets"]
    assert t["unique"]["n_seq"] == len(domains.load_proteins(root, V, "S0"))
    igrs = pq.read_table(next((root / "l0seq/igr").rglob("*.parquet"))).to_pylist()
    assert t["igr"]["n_seq"] == sum(1 for r in igrs if r["end"] - r["start"] >= 30) > 0
    assert t["models"]["n_seq"] == 6 and t["dark_orf"]["n_seq"] > 0
    assert "spacers" in man["missing_targets"]


def _copy(world, tmp_path) -> Path:
    root = world[0]
    dst = tmp_path / "a"
    shutil.copytree(root, dst)
    return dst


@pytest.mark.parametrize("fault", ["byte", "stray", "missing_meta", "count", "pack_id"])
def test_c12_pack_faults_are_caught(world, tmp_path, fault):
    root = _copy(world, tmp_path)
    pdir = pack.pack_dir(root, REL)
    pid = None
    if fault == "byte":
        p = pdir / "unique/part-0.faa.zst"
        b = bytearray(p.read_bytes())
        b[-5] ^= 0xFF
        p.write_bytes(bytes(b))
    elif fault == "stray":
        (pdir / "unique/part-9.faa.zst").write_bytes(b"x")
    elif fault == "missing_meta":
        (pdir / "igr/igr_meta.parquet").unlink()
    elif fault == "count":
        m = json.loads((pdir / "pack_manifest.json").read_text())
        m["targets"]["unique"]["n_seq"] += 1
        (pdir / "pack_manifest.json").write_text(json.dumps(m))
    elif fault == "pack_id":
        pid = b"\0" * 16
    assert pack.verify(root, REL, pid)


def test_c12_small_shards_fail_unless_last_and_pack_is_deterministic(world, tmp_path):
    root = _copy(world, tmp_path)
    man = pack.build(root, REL, models={"test-1": world[1] / "a.hmm"}, shard_text_bytes=20_000)
    assert sum(1 for s in man["shards"] if s["target"] == "unique") > 1
    assert any("not the last shard" in p for p in pack.verify(root, REL))
    assert pack.verify(root, REL, min_shard_bytes=0) == []
    again = pack.build(root, REL, models={"test-1": world[1] / "a.hmm"}, shard_text_bytes=20_000)
    assert again["pack_id"] == man["pack_id"]  # D-01: no timestamps in the pack


def test_pack_is_outside_the_release_manifest(world):
    man = catalog.build_manifest(world[0], REL)
    assert not any(f["path"].startswith(("search_pack/", "search/")) for f in man["files"])


def test_dark_orfs_translate_their_own_coordinates(world):
    """Each dark ORF's protein equals the translation of its slice from the sequence store."""
    from biobox.l0seq.store import SeqStore

    root = world[0]
    pdir = pack.pack_dir(root, REL)
    meta = {r["region_uid"].hex(): r for r in pq.read_table(pdir / "dark_orf/dark_orf_meta.parquet").to_pylist()}
    seqs = dict(r for s in sorted(pdir.glob("dark_orf/part-*.faa.zst")) for r in fastz.read_records(s))
    store = SeqStore(root)
    circ = {r["contig_uid"] for r in pq.read_table(next((root / "l3/tapes").rglob("*.parquet"))).to_pylist()
            if r["is_circular"]}
    for name in sorted(seqs)[::7]:
        m = meta[name]
        nt = store.get(m["contig_uid"], m["start"], m["end"], m["strand"], circular=m["contig_uid"] in circ)
        aa = orfs._translate(nt.decode().upper(), m["genetic_code"])
        assert aa.rstrip("*") == seqs[name] and "*" not in seqs[name]
        assert m["region_uid"] == ids.region_uid(m["contig_uid"], m["start"], m["end"], m["strand"])
    store.close()


def _declare_code(root: Path, contig_prefix: str, code: int) -> bytes:
    """Rewrite the genetic code of every slot of the contig whose uid starts with the hex prefix."""
    target = None
    for f in sorted((root / "l1/slot_protein").rglob("*.parquet")):
        t = pq.read_table(f)
        rows = t.to_pylist()
        for r in rows:
            if r["contig_uid"].hex().startswith(contig_prefix):
                r["genetic_code"] = code
                target = r["contig_uid"]
        pq.write_table(pa.Table.from_pylist(rows, schema=t.schema), f, compression="zstd")
    assert target is not None
    return target


def _oracle_dark(root: Path, cu: bytes, circular: bool, code: int) -> dict:
    """{(start, end, strand): aa} of the dark ORFs of one contig, from the r1 per-position ORF oracle."""
    from biobox.l0seq.store import SeqStore

    sys.path.insert(0, str(REPO / "tests/r1"))
    from test_search_units import _oracle_orfs

    genes = []
    for f in pack._occurrence_files(root, V):
        genes += [(r["start"], r["end"], r["strand"]) for r in pq.read_table(f).to_pylist() if r["contig_uid"] == cu]
    store = SeqStore(root)
    seq = store.get(cu).decode().upper()
    store.close()
    L = len(seq)
    covered = set()
    for g in genes:
        covered |= {x % L for x in range(g[0], g[1])} if circular else set(range(g[0], g[1]))
    own = {(1, e % L) for s_, e, st in genes if st == 1} | {(-1, s_ % L) for s_, e, st in genes if st == -1}
    out = {}
    for start, end, strand, aa, _ in _oracle_orfs(seq, circular, code, orfs.MIN_AA):
        if (strand, (end if strand == 1 else start) % L) in own:
            continue
        frac = len({x % L for x in range(start, end)} & covered) / (end - start)
        at_end = not circular and (start < orfs.END_BP or end > L - orfs.END_BP)
        if frac >= orfs.MAX_GENE_OVERLAP and not at_end:
            continue
        out[(start, end, strand)] = aa
    return out


# (contig 014272 carries IUPAC ambiguity codes; the plain-table oracle covers ACGT contigs only)
@pytest.mark.parametrize("prefix", ["be7014", "b4db86"])
@pytest.mark.covers("W-01")
def test_w01_pack_builds_with_a_contig_declared_code_15(world, tmp_path, prefix):
    root = _copy(world, tmp_path)
    cu = _declare_code(root, prefix, 15)
    pack.build(root, REL, models={"test-1": world[1] / "a.hmm"})
    pdir = pack.pack_dir(root, REL)
    meta = pq.read_table(pdir / "dark_orf/dark_orf_meta.parquet").to_pylist()
    seqs = dict(r for s in sorted(pdir.glob("dark_orf/part-*.faa.zst")) for r in fastz.read_records(s))
    circ = {r["contig_uid"] for r in pq.read_table(next((root / "l3/tapes").rglob("*.parquet"))).to_pylist()
            if r["is_circular"]}
    got = {(m["start"], m["end"], m["strand"]): seqs[m["region_uid"].hex()] for m in meta if m["contig_uid"] == cu}
    assert all(m["genetic_code"] == 15 for m in meta if m["contig_uid"] == cu)
    want = _oracle_dark(root, cu, cu in circ, 15)
    assert got == want and want
    assert want != _oracle_dark(root, cu, cu in circ, 11)  # TAG matters for this contig


@pytest.mark.parametrize("bad", [7, 26])
@pytest.mark.covers("W-01")
def test_w01_unsupported_code_fails_before_any_pack_file_is_written(world, tmp_path, bad):
    root = _copy(world, tmp_path)
    cu = _declare_code(root, "be7014", bad)
    pdir = pack.pack_dir(root, REL)
    shutil.rmtree(pdir)
    with pytest.raises(ValueError) as e:
        pack.build(root, REL, models={"test-1": world[1] / "a.hmm"})
    assert str(bad) in str(e.value) and cu.hex() in str(e.value)
    assert not pdir.exists()


# --- runs, controls (S-04) -----------------------------------------------------------------------------


def test_s04_positive_recovered_and_negative_within_expectation(world):
    hit, _ = _seed(world[0])
    row = _run(world, _recipe(hit, hit))
    assert row["status"] == "ok", row["error"]
    ctl = json.loads(row["controls_json"])
    assert ctl["positive"][0]["recovered"]
    assert ctl["negative"]["pass"] and ctl["negative"]["n_hits"] <= ctl["negative"]["allowed"]
    assert row["target_size"] == world[2]["targets"]["unique"]["n_seq"]
    assert row["z_prot"] == 1e8 and row["pack_id"] == bytes.fromhex(world[2]["pack_id"])


def test_s04_unrecovered_positive_fails_the_run(world):
    hit, other = _seed(world[0])
    row = _run(world, _recipe(hit, other))
    assert row["status"] == "failed" and "positive control" in row["error"]


@pytest.mark.parametrize("case", ["not_in_target", "unknown_model", "igr_target", "unknown_seed"])
def test_s04_submit_refusals(world, case):
    hit, _ = _seed(world[0])
    text = {
        "not_in_target": _recipe(hit, hit).replace(f"positive: [{{protein_uid: {hit.hex()}", "positive: [{protein_uid: "
                                                    + "ab" * 16),
        "unknown_model": _recipe(hit, hit, tool="hmmsearch").replace(f"seed: [{{protein_uid: {hit.hex()}}}]",
                                                                    "seed: [{hmm: PF99999}]"),
        "igr_target": _recipe(hit, hit, tool="nhmmer", target="igr"),
        "unknown_seed": _recipe(b"\xcd" * 16, hit),
    }[case]
    with pytest.raises(RecipeError):
        _run(world, text)


def test_hmm_seed_from_the_library_finds_its_stored_hits(world):
    """An hmmsearch with a library model finds every protein L1 stored for it (same Z, same cutoff)."""
    root = world[0]
    dom = pq.read_table(next((root / "l1/domains").rglob("*.parquet"))).to_pylist()
    stored = {r["protein_uid"] for r in dom if r["model_acc"] == "TEST0000"}
    pos = sorted(stored)[0]
    row = _run(world, _recipe(pos, pos, tool="hmmsearch").replace(f"seed: [{{protein_uid: {pos.hex()}}}]",
                                                                 "seed: [{hmm: TEST0000}]"))
    assert row["status"] == "ok", row["error"]
    assert row["rounds"] == 1
    assert {h["target_uid"] for h in runner.read_hits(root, row["run_id"].hex())} >= stored


# --- provenance and failure states (S-05), envelope (S-06) -----------------------------------------------


@pytest.mark.parametrize("fault,status", [("oom", "oom"), ("timeout", "timeout"), ("crash", "failed")])
def test_s05_injected_failures_are_states_never_zero_hits(world, fault, status):
    hit, _ = _seed(world[0])
    lim = runner.Limits(cpus=1, wall_s=3, mem_gb=4 if platform.system() == "Linux" else None)
    row = _run(world, _recipe(hit, hit), fault=fault, limits=lim)
    assert row["status"] == status, row["error"]
    d = runner.run_dir(world[0], row["run_id"])
    assert sorted(p.name for p in d.iterdir()) == ["run.parquet"]  # no hits file, so never "0 hits"
    assert runner.read_hits(world[0], row["run_id"].hex()) == []


def test_s05_every_hit_resolves_to_a_run_and_catalog_shows_runs(world):
    root = world[0]
    hit, _ = _seed(root)
    _run(world, _recipe(hit, hit, tool="jackhmmer", mx=3))
    man = catalog.build_manifest(root, REL)
    con = catalog.open_catalog(catalog.generate_catalog(man, world[1] / "c.duckdb"))
    try:
        orphans = con.execute("SELECT count(*) FROM search_hit h ANTI JOIN search_run r USING (run_id)").fetchone()[0]
        assert orphans == 0
        n_hits = con.execute("SELECT count(*) FROM search_hit").fetchone()[0]
        assert n_hits == sum(pq.read_metadata(p).num_rows for p in (root / "search").glob("run=*/hits.parquet"))
        # SD-02: every protein hit has locus context through protein_occurrence
        ctx = con.execute("SELECT count(*) FROM search_hit h JOIN protein_occurrence o ON o.protein_uid = h.target_uid "
                          "WHERE h.target_kind = 'protein'").fetchone()[0]
        assert ctx >= con.execute("SELECT count(DISTINCT target_uid) FROM search_hit WHERE target_kind = 'protein'"
                                  ).fetchone()[0] > 0
        assert con.execute("SELECT count(*) FROM search_region").fetchone()[0] == (
            world[2]["targets"]["igr"]["n_seq"] + world[2]["targets"]["dark_orf"]["n_seq"])
    finally:
        con.close()


def test_s06_cpu_cap_quota_and_envelope(world):
    hit, _ = _seed(world[0])
    row = _run(world, _recipe(hit, hit), fault="cpu", limits=runner.Limits(cpus=1, cpu_min=1 / 60, wall_s=60))
    assert row["status"] == "timeout" and "CPU" in row["error"]
    assert row["cpu_min_max_pid"] < 0.2
    row = _run(world, _recipe(hit, hit), limits=runner.Limits(tmp_quota_gb=1e-9))
    assert row["status"] == "failed" and "quota" in row["error"]
    row = _run(world, _recipe(hit, hit, tool="jackhmmer", mx=3))
    assert row["status"] == "ok"
    lim = json.loads(row["limits_json"])
    assert lim["cpu_min"] <= 300 and lim["nice"] == 15 and lim["tmp_quota_gb"] <= 150
    assert row["cpu_min_max_pid"] < 300 and row["wall_s"] < 60  # SEARCH_PLAN §5: S0 is seconds


def test_output_cap_keeps_partial_outputs_flagged(world):
    hit, _ = _seed(world[0])
    row = _run(world, _recipe(hit, hit, extra="max_output_gb: 0.000000001\n"))
    assert row["status"] == "output_cap"
    assert (runner.run_dir(world[0], row["run_id"]) / "hits.parquet").exists()


def test_dark_orf_target_runs(world):
    root = world[0]
    pdir = pack.pack_dir(root, REL)
    name, aa = next(fastz.read_records(next(pdir.glob("dark_orf/part-*.faa.zst"))))
    text = (f"recipe: dark\nwritten_against: {REL}\nseed: [{{sequence: {aa}}}]\ntarget: dark_orf\ntool: phmmer\n"
            f"iterations: {{max: 1, inclusion_bits: 20}}\ncontrols: {{positive: [{{region_uid: {name}}}]}}\n")
    row = _run(world, text)
    assert row["status"] == "ok", row["error"]
    hits = runner.read_hits(root, row["run_id"].hex())
    assert {h["target_kind"] for h in hits} == {"dark_orf"}


# --- exports and CLI (S-13) ---------------------------------------------------------------------------------


def test_s13_cli_json_records_match_the_locked_schemas(world, tmp_path):
    root = world[0]
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit, tool="jackhmmer", mx=2))
    cmd = [sys.executable, "-m", "biobox.cli"]
    out = subprocess.run(cmd + ["search", "run", str(rec), "--root", str(root), "--release", REL, "--wait",
                                "--json", "--cpus", "1", "--work", str(world[1] / "work")],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    run = json.loads(out[0])
    locked = json.loads(export.lock_path(REPO).read_text())
    assert run["schema"] == "search_run/1" and list(run) == locked["search_run"]["keys"]
    lines = subprocess.run(cmd + ["search", "hits", run["run_id"], "--root", str(root), "--json"],
                           capture_output=True, text=True, check=True).stdout.splitlines()
    hits = [json.loads(x) for x in lines]
    assert hits and all(list(h) == locked["search_hit"]["keys"] for h in hits)
    assert all(h["n_occurrences"] >= 1 for h in hits)
    lines = subprocess.run(cmd + ["nbr", "--anchors", run["run_id"], "--root", str(root), "--window", "3g", "--json"],
                           capture_output=True, text=True, check=True).stdout.splitlines()
    nb = [json.loads(x) for x in lines]
    assert nb and all(list(r) == locked["neighborhood"]["keys"] for r in nb)
    assert {r["anchor_protein_uid"] for r in nb} == {h["target_uid"] for h in hits}
    # card B-16 (AGENTS §3.3): every exported ordinal has its feature_uid beside it, in the run's tape version
    tv = tapes.tape_version([V])
    tps, feats = _tape_rows(root, tv), _features(root)
    assert all(h["tape_version"] == tv for h in hits)
    for h in hits:
        for o in h["occurrences"]:
            tape = tps[bytes.fromhex(o["contig_uid"])]
            assert tape["feature_uids"][o["slot_ordinal"]].hex() == o["feature_uid"]
            assert feats[bytes.fromhex(o["feature_uid"])]["start"] == o["start"]
    for r in nb:
        tape = tps[bytes.fromhex(r["contig_uid"])]
        assert r["tape_version"] == tv
        assert tape["feature_uids"][r["ordinal"]].hex() == r["feature_uid"]
        assert tape["feature_uids"][r["anchor_ordinal"]].hex() == r["anchor_feature_uid"]
        assert bytes.fromhex(r["feature_uid"]) in feats and bytes.fromhex(r["anchor_feature_uid"]) in feats


def test_async_submit_then_status(world, tmp_path):
    import time

    root = world[0]
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit))
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1), world[1] / "work")
    deadline = time.time() + 60
    while runner.status(root, rid)["status"] == "running" and time.time() < deadline:
        time.sleep(0.2)
    assert runner.status(root, rid)["status"] == "ok"


def test_neighborhoods_written_for_final_hits(world):
    root = world[0]
    hit, _ = _seed(root)
    row = _run(world, _recipe(hit, hit, extra="context: {window: 5g, edge: truncate}\n"))
    nb = pq.read_table(runner.run_dir(root, row["run_id"]) / "neighborhoods.parquet").to_pylist()
    final = {h["target_uid"] for h in runner.read_hits(root, row["run_id"].hex())}
    assert {r["anchor_protein_uid"] for r in nb} == final
    assert all(abs(r["offset"]) <= 5 for r in nb)


# --- card B-16: search outputs name slots by feature_uid and say which tape version -------------------


def _tape_rows(root: Path, tv: str) -> dict[bytes, dict]:
    return {r["contig_uid"]: r for f in sorted((root / "l3/tapes").glob(f"tape_version={tv}/stage=*/part-*.parquet"))
            for r in pq.read_table(f, partitioning=None).to_pylist()}


def _features(root: Path) -> dict[bytes, dict]:
    return {r["feature_uid"]: r for f in sorted((root / "l1/feature").glob("producer=*/stage=*/part-*.parquet"))
            for r in pq.read_table(f, partitioning=None).to_pylist()}


def _check_neighborhood_rows(root: Path, rows: list[dict], tv: str) -> None:
    """Each row's feature uids are the tape's at its ordinals (in tape version tv) and resolve in `feature`."""
    tps, feats = _tape_rows(root, tv), _features(root)
    assert rows
    for r in rows:
        tape = tps[r["contig_uid"]]
        assert r["tape_version"] == tv
        assert r["feature_uid"] == tape["feature_uids"][r["ordinal"]]
        assert r["anchor_feature_uid"] == tape["feature_uids"][r["anchor_ordinal"]]
        f = feats[r["feature_uid"]]
        assert (f["contig_uid"], f["start"], f["end"], f["strand"]) == (r["contig_uid"], r["start"], r["end"], r["strand"])
        assert feats[r["anchor_feature_uid"]]["contig_uid"] == r["contig_uid"]
    keys = [(r["contig_uid"], r["anchor_feature_uid"], r["feature_uid"], r["policy_hash"]) for r in rows]
    assert len(set(keys)) == len(keys)


@pytest.mark.covers("S-09")
def test_b16_run_neighborhoods_name_slots_by_feature_uid_in_the_packs_tape_version(world):
    root, _, man = world
    tv = tapes.tape_version([V])
    assert man["tape_version"] == tv
    assert table("neighborhoods").primary_key == ("contig_uid", "anchor_feature_uid", "feature_uid", "policy_hash")
    hit, _ = _seed(root)
    row = _run(world, _recipe(hit, hit, extra="context: {window: 5g, edge: truncate}\n"))
    t = pq.read_table(runner.run_dir(root, row["run_id"]) / "neighborhoods.parquet")
    assert conformance_problems(t, table("neighborhoods")) == []
    _check_neighborhood_rows(root, t.to_pylist(), tv)


@pytest.mark.covers("C-12")
def test_b16_search_region_flanks_are_the_igr_tables_feature_uids(world):
    root, _, man = world
    tv = tapes.tape_version([V])
    assert (man["tape_version"], man["igr_version"]) == (tv, igr.igr_version([V]))
    feats, tape = _features(root), _tape_rows(root, tv)
    src = {(r["contig_uid"], r["start"], r["end"]): r
           for f in sorted((root / "l0seq/igr").glob(f"igr_version={man['igr_version']}/stage=*/part-*.parquet"))
           for r in pq.read_table(f, partitioning=None).to_pylist()}
    pdir = pack.pack_dir(root, REL)
    n_flanks = 0
    for target in ("igr", "dark_orf"):
        t = pq.read_table(pdir / f"{target}/{target}_meta.parquet")
        assert conformance_problems(t, table("search_region")) == []
        for m in t.to_pylist():
            for side in ("left", "right"):
                fu, slot = m[f"{side}_feature_uid"], m[f"{side}_slot"]
                assert (fu is None) == (slot is None), (target, side)
                if fu is not None:
                    n_flanks += 1
                    assert feats[fu]["contig_uid"] == m["contig_uid"]
                    assert tape[m["contig_uid"]]["feature_uids"][slot] == fu  # the legacy slot names the same feature
            if target == "igr":
                s = src[(m["contig_uid"], m["start"], m["end"])]
                assert (m["left_feature_uid"], m["right_feature_uid"]) == (s["left_feature_uid"], s["right_feature_uid"])
    assert n_flanks > 0


NEW = "test-ncrna"
TWO = "rS0.two"


@pytest.fixture(scope="module")
def two(world, tmp_path_factory):
    """The micro-S0 store with a second producer: one kind-2 feature in an IGR of M. genitalium; IGRs, tapes
    and the pack of release TWO are built for both producers, beside the one-producer ones."""
    d = tmp_path_factory.mktemp("two")
    root = d / "a"
    shutil.copytree(world[0], root)
    circ = {cu for cu, t in _tape_rows(root, tapes.tape_version([V])).items() if t["is_circular"]}
    mg = next(iter(circ))  # the only circular contig is M. genitalium
    fu = micro._add_ncrna(root, mg)
    producers = [V, NEW]
    igr.run(root, "S0", circ, producers=producers)
    tapes.run(root, "S0", circ, producers=producers)
    man = pack.build(root, TWO, models={"test-1": world[1] / "a.hmm"}, producers=producers)
    return root, d, man, mg, fu


@pytest.mark.covers("MR-20", "S-09")
def test_b16_neighbor_feature_uids_stay_when_a_producer_adds_a_feature_and_ordinals_move(two):
    root, _, man, mg, fu = two
    tv1, tv2 = tapes.tape_version([V]), tapes.tape_version([V, NEW])
    assert man["tape_version"] == tv2 != tv1
    t2 = _tape_rows(root, tv2)[mg]
    i, n, k = t2["feature_uids"].index(fu), t2["n_slots"], 3
    # anchors on M. genitalium whose k-gene window in the two-producer tape does not reach the new feature
    far = [j for j in range(n) if min((j - i) % n, (i - j) % n) > k and t2["protein_uids"][j] is not None]
    far_fus = {t2["feature_uids"][j] for j in far}
    uids = {t2["protein_uids"][j] for j in far}
    policy = nbr.Policy.parse(window=f"{k}g")
    one, _ = nbr.fetch(root, uids, policy)  # default: the gene-call producer's tape version
    other, _ = nbr.fetch(root, uids, policy, tape_version=tv2)
    _check_neighborhood_rows(root, one, tv1)
    _check_neighborhood_rows(root, other, tv2)

    def kept(rows):  # coding neighbors of the far anchors (a protein may also occur next to the new feature)
        return [r for r in rows if r["protein_uid"] is not None and (r["contig_uid"] != mg or r["anchor_feature_uid"] in far_fus)]

    def triples(rows):
        return {(r["contig_uid"], r["anchor_feature_uid"], r["feature_uid"]) for r in kept(rows)}

    assert triples(one) == triples(other) and any(r["contig_uid"] == mg for r in kept(one))
    assert fu not in {r["feature_uid"] for r in other}
    ords = [{(r["anchor_feature_uid"], r["feature_uid"]): (r["anchor_ordinal"], r["ordinal"])
             for r in kept(rows) if r["contig_uid"] == mg} for rows in (one, other)]
    assert set(ords[0]) == set(ords[1])
    # the ordinals of the same features move by one downstream of the new feature (at index i of tape 2)
    for key, (a, o) in ords[0].items():
        assert ords[1][key] == (a + (a >= i), o + (o >= i))
    assert any(ords[0][key] != ords[1][key] for key in ords[0])


@pytest.mark.covers("S-05")
def test_b16_a_run_is_refused_on_a_tape_version_other_than_the_packs(two, capsys):
    root, d, man, mg, fu = two
    tv1, tv2 = tapes.tape_version([V]), tapes.tape_version([V, NEW])
    hit, _ = _seed(root)
    text = _recipe(hit, hit, extra="context: {window: 3g, edge: truncate}\n").replace(
        f"written_against: {REL}", f"written_against: {TWO}")
    before = sorted(runner.runs_dir(root).glob("run=*"))
    with pytest.raises(RecipeError, match=f"tape version {tv1}.*{tv2}"):
        runner.run(root, text, TWO, work_base=d / "work", limits=LIMITS, tape_version=tv1)
    assert sorted(runner.runs_dir(root).glob("run=*")) == before  # a refusal records nothing
    # the pack's own tape version is the run's: neighborhoods and hit occurrences use its ordinals
    row = runner.run(root, text, TWO, work_base=d / "work", limits=LIMITS)
    assert row["status"] == "ok", row["error"]
    nb = pq.read_table(runner.run_dir(root, row["run_id"]) / "neighborhoods.parquet").to_pylist()
    _check_neighborhood_rows(root, nb, tv2)
    assert cli.main(["search", "hits", row["run_id"].hex(), "--root", str(root), "--json"]) == 0
    hits = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    tps = _tape_rows(root, tv2)
    assert hits and all(h["tape_version"] == tv2 for h in hits)
    for o in (o for h in hits for o in h["occurrences"]):
        assert tps[bytes.fromhex(o["contig_uid"])]["feature_uids"][o["slot_ordinal"]].hex() == o["feature_uid"]
    assert cli.main(["nbr", "--anchors", row["run_id"].hex(), "--root", str(root), "--json",
                     "--tape-version", tv1]) == 2
    assert "tape version" in capsys.readouterr().err


# --- W-06: neighbor summary against a brute-force loop over tapes ------------------------------------------


def _r1():
    sys.path.insert(0, str(REPO / "tests/r1"))
    import test_search_units as r1

    return r1


def _w06_anchors(root: Path) -> set[bytes]:
    """Anchors near the origin of the circular contig (both sides), next to each linear contig end, beside each
    domain-carrying lambda gene, and the library-hit seed."""
    tps = _tape_rows(root, tapes.tape_version([V]))
    doms = {r["protein_uid"] for f in sorted((root / "l1/domains").rglob("*.parquet")) for r in pq.read_table(f).to_pylist()}
    out = {_seed(root)[0]}
    for t in tps.values():
        p = t["protein_uids"]
        out |= {p[1], p[-1]}  # circular: windows across the origin; linear: windows cut by the contig end
        if not t["is_circular"]:
            out |= {p[j + 2] for j, u in enumerate(p) if u in doms and j + 2 < len(p)}
    assert None not in out
    return out


def _w06_store(root: Path, release: str = REL) -> dict:
    """Domains, clusters, protein lengths and genomes straight from the Parquet shards."""
    dom: dict[bytes, dict[tuple[str, str], str]] = {}
    for f in sorted((root / "l1/domains").glob("library=*/stage=*/part-*.parquet")):
        lib = f.parent.parent.name.split("=", 1)[1]
        for r in pq.read_table(f, partitioning=None).to_pylist():
            dom.setdefault(r["protein_uid"], {})[(lib, r["model_acc"])] = r["model_name"]
    clus: dict[str, dict[bytes, tuple[bytes, bytes]]] = {"nr50": {}, "nr90": {}}
    for f in sorted((root / f"l2/cluster_member/release={release}").glob("part-*.parquet")):
        for r in pq.read_table(f, partitioning=None).to_pylist():
            if r["level"] in clus:
                clus[r["level"]][r["member_uid"]] = (r["cluster_uid"], r["rep_uid"])
    lens = {r["protein_uid"]: r["length"] for f in sorted((root / "l0seq/protein_seq").rglob("part-*.parquet"))
            for r in pq.read_table(f, columns=["protein_uid", "length"]).to_pylist()}
    gen = {r["contig_uid"]: r["genome_uid"] for f in sorted((root / f"dim/contig_dim/release={release}").glob("*.parquet"))
           for r in pq.read_table(f, partitioning=None).to_pylist()}
    return {"dom": dom, "clus": clus, "len": lens, "genome": gen}


def _oracle_summary(root: Path, anchors: set[bytes], window: str, levels=("pfam", "nr50"), release: str = REL) -> dict:
    """The summary by a loop over every tape slot: each anchor slot's window from the S-09 coordinate oracle
    (tests/r1 `_oracle_nbr`), then every count by hand."""
    r1 = _r1()
    policy = nbr.Policy.parse(window=window)
    st_ = _w06_store(root, release)
    tps = _tape_rows(root, tapes.tape_version([V]))
    loci, contigs, censored, placed, nbrs = set(), set(), 0, set(), []
    for cu in sorted(tps):
        t = tps[cu]
        slots = list(zip(t["start_bps"], t["end_bps"], t["strands"]))
        for i, pu in enumerate(t["protein_uids"]):
            if pu not in anchors:
                continue
            got = r1._oracle_nbr(slots, t["contig_len"], t["is_circular"], i, policy)
            if got is None:
                continue
            loci.add((cu, i))
            contigs.add(cu)
            placed.add(pu)
            censored += any(d == 0 and fl & 1 for _, d, _, _, fl in got)
            for j, d, rel, _, _ in got:
                if d != 0:
                    nbrs.append({"locus": (cu, i), "protein": t["protein_uids"][j], "rel": rel,
                                 "same": slots[j][2] == slots[i][2], "contig": cu})
    basis = "genome_uid" if st_["genome"] and all(c in st_["genome"] for c in contigs) else "contig_uid"
    genome = (lambda c: st_["genome"][c]) if basis == "genome_uid" else (lambda c: c)
    n_loci = len(loci)
    rows, n_keys = [], {}
    for level in levels:
        groups: dict[str, list[dict]] = {}
        meta: dict[str, dict] = {}
        for nb in nbrs:
            u = nb["protein"]
            if level == "pfam":
                keys = {f"{lib}:{acc}": {"library": lib, "model_acc": acc, "model_name": name}
                        for (lib, acc), name in st_["dom"].get(u, {}).items()}
            else:
                keys = {}
                if u in st_["clus"][level]:
                    c, rep = st_["clus"][level][u]
                    keys = {c.hex(): {"cluster_uid": c.hex(), "rep_uid": rep.hex(), "rep_length": st_["len"].get(rep)}}
            for k, m in keys.items():
                groups.setdefault(k, []).append(nb)
                meta[k] = m
            if not keys:
                groups.setdefault("", []).append(nb)
        if level != "pfam":
            for k, m in meta.items():
                names, reps = {}, []
                for member, (c, rep) in st_["clus"][level].items():
                    if c.hex() == k:
                        reps.append(rep)
                        for name in set(st_["dom"].get(member, {}).values()):
                            names[name] = names.get(name, 0) + 1
                m["model_name"] = min(names, key=lambda x: (-names[x], x)) if names else ""
                # the stand-in clusters name each member its own representative: take the least, as documented
                m["rep_uid"] = min(reps).hex()
                m["rep_length"] = st_["len"].get(min(reps))

        def row(k, g):
            hist: dict[int, int] = {}
            for nb in g:
                hist[nb["rel"]] = hist.get(nb["rel"], 0) + 1
            n = len(g)
            ls = {nb["locus"] for nb in g}
            return {"level": level, "key": k, "unannotated": k == "", **meta.get(k, {}),
                    "n_loci": len(ls), "share_of_loci": round(len(ls) / n_loci, 4) if n_loci else None,
                    "n_neighbors": n, "n_genomes": len({genome(nb["contig"]) for nb in g}), "genome_basis": basis,
                    "offsets": {str(o): hist[o] for o in sorted(hist)},
                    "modal_offset": min(hist, key=lambda o: (-hist[o], abs(o), -o)) if hist else None,
                    "same_strand_share": round(sum(nb["same"] for nb in g) / n, 4) if n else None}

        keyed = sorted((row(k, g) for k, g in groups.items() if k), key=lambda r: (-r["n_loci"], r["key"]))
        n_keys[level] = len(keyed)
        rows += keyed + [row("", groups.get("", []))]
    return {"header": {"n_anchor_proteins": len(placed), "n_anchor_loci": n_loci, "n_contigs": len(contigs),
                       "n_censored_loci": censored, "n_neighbors": len(nbrs), "genome_basis": basis, "n_keys": n_keys},
            "rows": rows}


def _compare(got: dict, want: dict) -> None:
    for k, v in want["header"].items():
        assert got["header"][k] == v, k
    assert [(r["level"], r["key"]) for r in got["rows"]] == [(r["level"], r["key"]) for r in want["rows"]]
    for g, w in zip(got["rows"], want["rows"]):
        assert {k: g[k] for k in w} == w, (w["level"], w["key"])


@pytest.mark.covers("W-06", "S-09")
@pytest.mark.parametrize("window", ["3g", "5g"])
def test_w06_summary_equals_a_brute_force_loop_over_tapes(world, window):
    root = world[0]
    anchors = _w06_anchors(root)
    want = _oracle_summary(root, anchors, window)
    got = nbr.summarize(root, anchors, nbr.Policy.parse(window=window), by=("pfam", "nr50"), release=REL, top_n=0)
    _compare(got, want)
    h = got["header"]
    assert (h["release"], h["tape_version"], h["policy_hash"]) == (REL, tapes.tape_version([V]),
                                                                  nbr.Policy.parse(window=window).hash())
    # the world is what the card asks for: an anchor window across the origin, windows cut at linear ends,
    # neighbors with and without a family and two neighbors of one family or cluster at one locus
    assert want["header"]["n_censored_loci"] >= 2 and h["genome_basis"] == "contig_uid"
    rows, _ = nbr.fetch(root, anchors, nbr.Policy.parse(window=window))
    assert any(r["flags"] & nbr.FLAG_CROSSES_ORIGIN for r in rows)
    for level in ("pfam", "nr50"):
        lv = [r for r in want["rows"] if r["level"] == level]
        assert lv[-1]["unannotated"] and len(lv) > 2
    assert [r for r in want["rows"] if r["level"] == "pfam"][-1]["n_neighbors"] > 0
    nr = [r for r in want["rows"] if r["level"] == "nr50"]
    assert sum(r["n_neighbors"] for r in nr) == want["header"]["n_neighbors"]
    # top-N keeps totals: the cap drops rows only
    capped = nbr.summarize(root, anchors, nbr.Policy.parse(window=window), release=REL, top_n=2)
    assert capped["header"]["n_keys"] == want["header"]["n_keys"]
    assert [(r["level"], r["key"]) for r in capped["rows"]] == [
        (r["level"], r["key"]) for lv in ("pfam", "nr50")
        for r in [x for x in want["rows"] if x["level"] == lv and not x["unannotated"]][:2]
        + [x for x in want["rows"] if x["level"] == lv and x["unannotated"]]]


@pytest.mark.covers("W-06")
def test_w06_summary_counts_genomes_from_contig_dim_when_present(world, tmp_path):
    root = _copy(world, tmp_path)
    tps = _tape_rows(root, tapes.tape_version([V]))
    lin = sorted(cu for cu, t in tps.items() if not t["is_circular"])
    one = ids.genome_uid(lin)  # both phages in one made-up genome, so genomes < contigs
    rows = [{"contig_uid": cu, "genome_uid": one if cu in lin else ids.genome_uid([cu]), "contig_len": t["contig_len"],
             "is_circular": t["is_circular"], "tier": None, "assembly_quality": None, "tax_order": None,
             "is_species_rep": False} for cu, t in sorted(tps.items())]
    d = root / f"dim/contig_dim/release={REL}"
    d.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=table("contig_dim").schema), d / "part-0.parquet")
    anchors = _w06_anchors(root)
    want = _oracle_summary(root, anchors, "5g", levels=("pfam", "nr50", "nr90"))
    got = nbr.summarize(root, anchors, nbr.Policy.parse(window="5g"), by=("pfam", "nr50", "nr90"), release=REL,
                        top_n=0)
    assert want["header"]["genome_basis"] == "genome_uid"
    _compare(got, want)
    # the made-up genome changes some counts (the world without contig_dim counts contigs)
    by_contig = _oracle_summary(world[0], anchors, "5g", levels=("pfam", "nr50", "nr90"))
    assert by_contig["header"]["genome_basis"] == "contig_uid"
    assert [r["n_genomes"] for r in want["rows"]] != [r["n_genomes"] for r in by_contig["rows"]]
    # file lists from the release manifest give the same summary as the store's layout
    man = catalog.build_manifest(root, REL)
    assert nbr.summarize(root, anchors, nbr.Policy.parse(window="5g"), by=("pfam", "nr50", "nr90"), release=REL,
                         top_n=0, manifest=man) == got


@pytest.mark.covers("W-06", "S-13")
def test_w06_cli_summary_json_matches_the_lock_and_tsv_reads_back(world, capsys):
    root = world[0]
    hit, _ = _seed(root)
    row = _run(world, _recipe(hit, hit, tool="jackhmmer", mx=2))
    assert row["status"] == "ok", row["error"]
    rid = row["run_id"].hex()
    args = ["nbr", "summarize", rid, "--root", str(root), "--window", "5g", "--by", "pfam,nr50"]
    assert cli.main(args + ["--json"]) == 0
    recs = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    locked = json.loads(export.lock_path(REPO).read_text())
    assert recs[0]["schema"] == "nbr_summary_header/1" and list(recs[0]) == locked["nbr_summary_header"]["keys"]
    assert len(recs) > 3 and all(list(r) == locked["nbr_summary"]["keys"] for r in recs[1:])
    assert {"nbr_summary", "nbr_summary_header"} <= set(locked)
    assert export.lock_problems(locked) == [] and locked == json.loads(json.dumps(export.SCHEMAS, sort_keys=True))
    # the run's hits are the anchors, its release and tape version the summary's
    anchors = {h["target_uid"] for h in runner.read_hits(root, rid) if h["target_kind"] == "protein"}
    s = nbr.summarize(root, anchors, nbr.Policy.parse(window="5g"), release=REL)
    assert recs == list(export.nbr_summary_records(s))
    assert recs[0]["release"] == REL and recs[0]["tape_version"] == tapes.tape_version([V])
    assert recs[0]["top_n"] == 50 and recs[0]["levels"] == ["pfam", "nr50"]
    assert cli.main(args + ["--format", "tsv"]) == 0
    r1 = _r1()
    back = r1.parse_summary_tsv(capsys.readouterr().out)
    assert r1.same_numbers(back["rows"], recs[1:])
    for k in ("n_anchor_proteins", "n_anchor_loci", "n_contigs", "n_censored_loci", "release", "tape_version",
              "policy", "n_keys"):
        assert back["header"][k] == recs[0][k], k
    # anchors as a list of protein uids take the same path; a uid list needs the release named
    uids = ",".join(sorted(u.hex() for u in anchors))
    assert cli.main(["nbr", "summarize", uids, "--root", str(root), "--window", "5g", "--release", REL, "--json"]) == 0
    assert [json.loads(x) for x in capsys.readouterr().out.splitlines()] == recs
    assert cli.main(["nbr", "summarize", uids, "--root", str(root), "--window", "5g"]) == 2
    assert "--release" in capsys.readouterr().err


# --- S-07: planted drift trap ---------------------------------------------------------------------------------

AA = "ACDEFGHIKLMNPQRSTVWY"


def _mutate(s, frac, rng):
    s = list(s)
    for i in rng.sample(range(len(s)), int(frac * len(s))):
        s[i] = rng.choice(AA)
    return "".join(s)


def _hmm_of(name: str, acc: str, members: list[str], ga: float, path: Path) -> None:
    alpha = pyhmmer.easel.Alphabet.amino()
    msa = pyhmmer.easel.TextMSA(name=name.encode(), sequences=[
        pyhmmer.easel.TextSequence(name=f"m{i}".encode(), sequence=s) for i, s in enumerate(members)]).digitize(alpha)
    hmm, _, _ = pyhmmer.plan7.Builder(alpha, seed=1).build_msa(msa, pyhmmer.plan7.Background(alpha))
    hmm.name, hmm.accession = name, acc
    hmm.cutoffs.gathering = (ga, ga)
    with open(path, "wb") as fh:
        hmm.write(fh)


def _place(root: Path, loci: list[tuple[str, str, str, bool]]) -> None:
    """One linear 1.7 kb contig per locus (left, target, right), all + strand; `partial` puts
    the target at the contig end (dist_to_contig_end_bp = 0, as for a partial gene call)."""
    from biobox.l3.tapes import build_tape, tape_version

    occ, tapes_ = [], []
    for k, (left, target, right, partial) in enumerate(loci):
        cu = ids.contig_uid(f"synthetic-contig-{k}".encode())
        gu = ids.genome_uid([cu])
        coords = [(100, 460, 1), (600, 1080, 1), (1200, 1560, 1)]
        prots = [ids.protein_uid(x) for x in (left, target, right)]
        for i, ((st_, en, sd), pu) in enumerate(zip(coords, prots)):
            occ.append({"protein_uid": pu, "contig_uid": cu, "slot_ordinal": i, "genome_uid": gu, "start": st_,
                        "end": en, "strand": sd, "tax_order": None, "is_species_rep": True, "tier": 1,
                        "assembly_quality": "complete",
                        "dist_to_contig_end_bp": 0 if partial and i == 1 else min(st_, 1700 - en),
                        "feature_uid": ids.feature_uid(cu, st_, en, sd, ids.KIND["cds"])})
        tapes_.append({"contig_uid": cu, "tax_order_bucket": None, "genome_uid": gu, "contig_len": 1700,
                       "is_circular": False, "n_slots": 3, "tier": 1, "assembly_quality": "complete",
                       **build_tape(coords, 1700, False), "protein_uids": prots,
                       "feature_uids": [ids.feature_uid(cu, a, b, sd, ids.KIND["cds"]) for a, b, sd in coords]})
    occ.sort(key=lambda r: (r["protein_uid"], r["contig_uid"], r["slot_ordinal"]))
    tapes_.sort(key=lambda r: r["contig_uid"])
    for rel, rows, name in ((f"l1/protein_occurrence/gene_call_version={V}/stage=S0", occ, "protein_occurrence"),
                            (f"l3/tapes/tape_version={tape_version([V])}/stage=S0", tapes_, "contig_tapes")):
        (root / rel).mkdir(parents=True)
        pq.write_table(pa.Table.from_pylist(rows, schema=table(name).schema), root / rel / "part-0.parquet")


def trap_world(root: Path, lib: Path, seed: int = 7, placed: bool = False) -> dict:
    """Seed family A (10), library family B (30, GA 60) and 5 bridges half way from A to B on
    the sites where they differ, plus 50 random decoys; A shares ~24% identity with B.

    `placed` adds locus context: every A member, bridge and B member on its own contig; A
    members and bridges flanked by members of a library family M (`ctx-1`, GA 60), B members by
    proteins with no library hit; three A fragments (the first 40%) flanked by M; and two A
    members placed as partial genes. The base sequences do not change."""
    rng = random.Random(seed)
    A = "".join(rng.choice(AA) for _ in range(160))
    sites = rng.sample(range(160), 128)
    B = "".join(rng.choice(AA) if i in set(sites) else c for i, c in enumerate(A))
    fam_a = [_mutate(A, 0.15, rng) for _ in range(10)]
    fam_b = [_mutate(B, 0.15, rng) for _ in range(30)]
    bridges = []
    for _ in range(5):
        s = list(A)
        for i in rng.sample(sites, 64):
            s[i] = B[i]
        bridges.append(_mutate("".join(s), 0.075, rng))
    decoys = ["".join(rng.choice(AA) for _ in range(160)) for _ in range(50)]
    world = {"A": A, "a": fam_a, "b": {ids.protein_uid(s) for s in fam_b}}
    extra: list[str] = []
    if placed:
        r2 = random.Random(seed + 1)
        M = "".join(r2.choice(AA) for _ in range(120))
        frags = [fam_a[i][:64] for i in (7, 8, 9)]
        partial = fam_a[5:7]
        m_members = [_mutate(M, 0.15, r2) for _ in range(2 * (len(fam_a) + len(bridges) + len(frags)))]
        plain = ["".join(r2.choice(AA) for _ in range(120)) for _ in range(2 * len(fam_b))]
        it, ip = iter(m_members), iter(plain)
        loci = [(next(it), t, next(it), t in partial) for t in fam_a + bridges + frags]
        loci += [(next(ip), t, next(ip), False) for t in fam_b]
        _place(root, loci)
        extra = frags + m_members + plain
        world.update(fragments={ids.protein_uid(x) for x in frags}, partial={ids.protein_uid(x) for x in partial})
    rows = sorted(({"protein_uid": ids.protein_uid(s), "protein_digest": ids.raw_digest(s.encode()), "length": len(s),
                    "seq": s} for s in set(fam_a + fam_b + bridges + decoys + extra)), key=lambda r: r["protein_uid"])
    out = root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0/protein_prefix=all/part-0.parquet"
    out.parent.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=table("protein_seq").schema), out)
    _hmm_of("trapB", "TRAP0001", fam_b[:8], 60.0, lib)
    domains.run(root, "S0", V, "trap", "1", lib, "gathering")
    world["models"] = {"trap-1": lib}
    if placed:
        ctx = lib.with_name("ctx.hmm")
        _hmm_of("ctxM", "CTX0001", m_members[:8], 60.0, ctx)
        domains.run(root, "S0", V, "ctx", "1", ctx, "gathering")
        world["models"]["ctx-1"] = ctx
    return world


def test_s07_drift_alarm_reverts_a_planted_corruption(tmp_path):
    root = tmp_path / "a"
    w = trap_world(root, tmp_path / "trap.hmm")
    pack.build(root, "rSyn", models={"trap-1": tmp_path / "trap.hmm"})
    rows = {}
    for name, stop in (("guarded", "[no_new_reps, drift_alarm]"), ("unguarded", "[no_new_reps]")):
        text = (f"recipe: trap\nwritten_against: rSyn\nseed: [{{sequence: {w['A']}}}]\ntarget: unique\n"
                f"tool: jackhmmer\niterations: {{max: 5, inclusion_bits: 30, stop: {stop}}}\n"
                f"controls: {{positive: [{{sequence: {w['a'][0]}}}]}}\n")
        rows[name] = runner.run(root, text, "rSyn", work_base=tmp_path / "work", limits=LIMITS)
    g, u = rows["guarded"], rows["unguarded"]
    leaked = {n: len({h["target_uid"] for h in runner.read_hits(root, r["run_id"].hex())} & w["b"])
              for n, r in rows.items()}
    # the trap is real: without the alarm the profile drifts into family B
    assert u["status"] == "ok" and leaked["unguarded"] >= 25
    # with it, the round that drifts is logged, dropped, and the run reverts to the round before
    assert g["status"] == "ok", g["error"]
    assert g["stop_reason"].startswith("drift_alarm: TRAP0001") and g["final_round"] == g["rounds"] - 1
    assert leaked["guarded"] <= 3
    dropped = [h for h in runner.read_hits(root, g["run_id"].hex(), round_=g["rounds"], final_only=False)]
    assert dropped and all(h["curation_label"] == "drop" and h["curation_reason"] == "drift_alarm" for h in dropped)


def _warnings(root, run):
    import json

    row = pq.read_table(next((root / "search").glob(f"run={run['run_id'].hex()}/run.parquet"))).to_pylist()[0]
    return json.loads(row["warnings_json"])


def test_w08_drift_alarm_warns_when_the_release_has_no_domain_table(tmp_path):
    import shutil

    texts = {}
    for name in ("with", "without"):
        root = tmp_path / name
        w = trap_world(root, tmp_path / f"{name}.hmm")
        pack.build(root, "rSyn", models={"trap-1": tmp_path / f"{name}.hmm"})
        if name == "without":
            shutil.rmtree(root / "l1/domains")
        text = (f"recipe: trap\nwritten_against: rSyn\nseed: [{{sequence: {w['A']}}}]\ntarget: unique\n"
                f"tool: jackhmmer\niterations: {{max: 3, inclusion_bits: 30, stop: [no_new_reps, drift_alarm]}}\n"
                f"controls: {{positive: [{{sequence: {w['a'][0]}}}]}}\n")
        texts[name] = runner.run(root, text, "rSyn", work_base=tmp_path / f"work-{name}", limits=LIMITS)
        texts[name] = (root, texts[name])
    (r1, with_), (r2, without) = texts["with"], texts["without"]
    assert not any("domain annotations" in x for x in _warnings(r1, with_))
    assert any("drift_alarm could not be evaluated" in x and "no domain annotations" in x for x in _warnings(r2, without))
    assert without["status"] == with_["status"] == "ok"


def test_search_tables_conform_to_the_registry(world):
    root = world[0]
    for p in sorted((root / "search").glob("run=*/*.parquet")):
        name = {"run.parquet": "search_run", "hits.parquet": "search_hit", "profiles.parquet": "search_profile",
                "neighborhoods.parquet": "neighborhoods"}[p.name]
        assert contracts is not None
        from biobox.registry.schemas import conformance_problems

        assert conformance_problems(pq.read_table(p), table(name)) == []
        if name == "search_profile":  # blobs decode on their own: HMMER3 text and Stockholm
            for r in pq.read_table(p).to_pylist():
                assert fastz.unzst(r["profile_zst"]).startswith(b"HMMER3")
                assert r["msa_zst"] is None or fastz.unzst(r["msa_zst"]).startswith(b"# STOCKHOLM")


def _trap_recipe(w, stop="[no_new_reps]", tool="jackhmmer", mx=5, extra=""):
    return (f"recipe: trap\nwritten_against: rSyn\nseed: [{{sequence: {w['A']}}}]\ntarget: unique\n"
            f"tool: {tool}\niterations: {{max: {mx}, inclusion_bits: 30, stop: {stop}}}\n"
            f"controls: {{positive: [{{sequence: {w['a'][0]}}}]}}\n{extra}")


@pytest.fixture(scope="module")
def placed(tmp_path_factory):
    d = tmp_path_factory.mktemp("placed")
    root = d / "a"
    w = trap_world(root, d / "trap.hmm", placed=True)
    pack.build(root, "rSyn", models=w["models"])
    return root, d, w


def test_s07_context_lost_reverts_when_new_hits_lose_round1_neighbors(placed):
    root, d, w = placed
    row = runner.run(root, _trap_recipe(w, "[no_new_reps, context_lost]"), "rSyn", work_base=d / "work", limits=LIMITS)
    assert row["status"] == "ok", row["error"]
    assert row["stop_reason"].startswith("context_lost") and row["final_round"] == row["rounds"] - 1
    log = json.loads(row["rounds_json"])
    assert log[-1]["context_retention"] < 0.5
    final = {h["target_uid"] for h in runner.read_hits(root, row["run_id"].hex())}
    assert len(final & w["b"]) <= 3
    # the same world without the rule drifts into B: the check, not the world, stops it
    free = runner.run(root, _trap_recipe(w), "rSyn", work_base=d / "work", limits=LIMITS)
    assert len({h["target_uid"] for h in runner.read_hits(root, free["run_id"].hex())} & w["b"]) >= 25


@pytest.mark.parametrize("policy", ["exclude", "flag", "include"])
def test_curation_labels_for_coverage_and_partial_genes(placed, policy):
    root, d, w = placed
    text = _trap_recipe(w, tool="phmmer", mx=1,
                        extra=f"partial: {policy}\ncuration: {{min_query_coverage: 0.9}}\n")
    row = runner.run(root, text, "rSyn", work_base=d / "work", limits=LIMITS)
    assert row["status"] == "ok", row["error"]
    rows = runner.read_hits(root, row["run_id"].hex(), round_=1, final_only=False)
    label = {r["target_uid"]: (r["curation_label"], r["curation_reason"], r["is_partial"]) for r in rows}
    assert w["fragments"] <= set(label)
    for u in w["fragments"]:
        assert label[u][0] == "drop" and label[u][1].startswith("query_coverage")
    for u in w["partial"]:
        lab, reason, is_partial = label[u]
        assert is_partial
        assert (lab, reason) == (("drop", "partial") if policy == "exclude" else ("keep", None))
    kept = {u for u, (lab, _, _) in label.items() if lab == "keep"}
    assert not kept & w["fragments"]  # dropped targets never enter the next profile or the final set
    final = {h["target_uid"] for h in runner.read_hits(root, row["run_id"].hex())}
    assert final == kept


HMMSEARCH = shutil.which("hmmsearch")


@pytest.mark.skipif(HMMSEARCH is None, reason="hmmsearch binary not on PATH")
def test_rounds_equal_the_hmmsearch_binary(placed, tmp_path):
    """O-SEARCH oracle: each round's stored profile, searched by HMMER's own hmmsearch at the
    same Z and threshold, gives the same included targets, sequence scores and kept domains."""
    root, d, w = placed
    row = runner.run(root, _trap_recipe(w), "rSyn", work_base=d / "work", limits=LIMITS)
    assert row["status"] == "ok" and row["rounds"] >= 2
    faa = tmp_path / "unique.faa"
    with open(faa, "w") as fh:
        for sh in sorted(pack.pack_dir(root, "rSyn").glob("unique/part-*.faa.zst")):
            for name, seq in fastz.read_records(sh):
                fh.write(f">{name}\n{seq}\n")
    profiles = pq.read_table(runner.run_dir(root, row["run_id"]) / "profiles.parquet").to_pylist()
    hits = pq.read_table(runner.run_dir(root, row["run_id"]) / "hits.parquet").to_pylist()
    bits = 30.0
    for prof in profiles:
        k = prof["round"]
        hmm = tmp_path / f"r{k}.hmm"
        hmm.write_bytes(fastz.unzst(prof["profile_zst"]))
        tbl, dom = tmp_path / f"r{k}.tbl", tmp_path / f"r{k}.dom"
        subprocess.run([HMMSEARCH, "--noali", "--cpu", "1", "-Z", "1e8", "--domZ", "1e8", "-T", str(bits),
                        "--incT", str(bits), "--tblout", str(tbl), "--domtblout", str(dom), "-o", "/dev/null",
                        str(hmm), str(faa)], check=True)
        seq_scores = {}
        for line in tbl.read_text().splitlines():
            if not line.startswith("#"):
                f = line.split()
                if float(f[5]) >= bits:
                    seq_scores[f[0]] = float(f[5])
        doms: dict[str, list[tuple[int, int, float]]] = {}
        for line in dom.read_text().splitlines():
            if not line.startswith("#"):
                f = line.split()
                if f[0] in seq_scores:
                    doms.setdefault(f[0], []).append((int(f[19]) - 1, int(f[20]), float(f[13])))
        mine = [h for h in hits if h["round"] == k]
        got_seq = {h["target_uid"].hex(): h["seq_bitscore"] for h in mine}
        assert set(got_seq) == set(seq_scores), f"round {k}: included targets differ"
        assert all(abs(got_seq[t] - seq_scores[t]) <= 0.051 for t in seq_scores)
        for t, ds in doms.items():
            best = max(ds, key=lambda x: x[2])
            want = sorted((a, b) for a, b, sc in ds if sc >= bits / 2 or (a, b, sc) == best)
            have = sorted((h["env_start"], h["env_end"]) for h in mine if h["target_uid"].hex() == t)
            assert have == want, f"round {k}, {t}: kept domains differ"


# --- W-02: a user's HMM file and seed sequences as seeds -------------------------------------------------


def _library_models(path: Path, accs: list[str]) -> list:
    hmms = list(domains.iter_hmms(path))
    return [next(h for h in hmms if engine.text(h.accession) == a) for a in accs]


def _write_hmms(path: Path, hmms: list) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        for h in hmms:
            h.write(fh)
    return path


def _write_faa(path: Path, records: list[tuple[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{n}\n{s}\n" for n, s in records))
    return path


def _file_recipe(kind: str, path: Path, positive: bytes, rel: str = REL, tool: str = "hmmsearch", mx: int = 1,
                 bits: float = 25, seed_extra: str = "") -> str:
    return (f"recipe: w02\nwritten_against: {rel}\nseed: [{{{kind}: '{path}'{seed_extra}}}]\ntarget: unique\n"
            f"tool: {tool}\niterations: {{max: {mx}, inclusion_bits: {bits}}}\n"
            f"controls: {{positive: [{{protein_uid: {positive.hex()}}}]}}\n")


def _stored(root: Path, acc: str) -> set[bytes]:
    return {r["protein_uid"] for f in sorted((root / "l1/domains").rglob("*.parquet"))
            for r in pq.read_table(f).to_pylist() if r["model_acc"] == acc}


def _rows(root: Path, row: dict) -> list[tuple]:
    return sorted((h["round"], h["target_uid"], h["bitscore"], h["seq_bitscore"], h["evalue_zref"], h["seq_evalue_zref"],
                   h["env_start"], h["env_end"], h["ali_start"], h["ali_end"], h["hmm_from"], h["hmm_to"],
                   h["model_len"], h["curation_label"])
                  for h in runner.read_hits(root, row["run_id"].hex(), final_only=False))


def _pinned(root: Path, text: str, rel: str = REL) -> tuple[str, str, dict]:
    r, _ = rcp.load(text)
    res = runner.resolve(root, r, rel)
    return rcp.recipe_hash(res["recipe"]), runner.query_hash(res["seed"]), res


@pytest.mark.covers("W-02")
def test_w02_hmm_file_returns_exactly_the_library_models_hits(world, tmp_path):
    root, d, _ = world
    pos = sorted(_stored(root, "TEST0000"))[0]
    f = _write_hmms(tmp_path / "own" / "lam0.hmm", _library_models(d / "a.hmm", ["TEST0000"]))
    lib = _run(world, _recipe(pos, pos, tool="hmmsearch").replace(f"seed: [{{protein_uid: {pos.hex()}}}]",
                                                                 "seed: [{hmm: TEST0000}]"))
    own = _run(world, _file_recipe("hmm_file", f, pos))
    assert lib["status"] == own["status"] == "ok", (lib["error"], own["error"])
    assert (own["z_prot"], own["rounds"], own["final_round"]) == (lib["z_prot"], 1, 1)
    want = _rows(root, lib)
    assert want and _rows(root, own) == want
    assert pos in {r[1] for r in want}


def _score_line_edit(text: str) -> str:
    """One digit changed in the first COMPO line (a score line); the file still parses."""
    lines = text.splitlines(keepends=True)
    i = next(k for k, line in enumerate(lines) if line.split()[:1] == ["COMPO"])
    j = next(k for k, c in enumerate(lines[i]) if c.isdigit())
    lines[i] = lines[i][:j] + ("2" if lines[i][j] != "2" else "3") + lines[i][j + 1:]
    return "".join(lines)


@pytest.mark.covers("W-02")
@pytest.mark.parametrize("kind", ["hmm_file", "sequences", "msa"])
def test_w02_seed_file_hashes_follow_its_content_not_its_path(world, tmp_path, kind):
    root, d, _ = world
    pos = sorted(_stored(root, "TEST0000"))[0]
    prots = domains.load_proteins(root, V, "S0")
    three = [(f"p{i}", prots[u].rstrip("*")) for i, u in enumerate(sorted(prots)[:3])]
    if kind == "hmm_file":
        f = _write_hmms(tmp_path / "a" / "m.hmm", _library_models(d / "a.hmm", ["TEST0000"]))
        other = _write_hmms(tmp_path / "o" / "m.hmm", _library_models(d / "a.hmm", ["TEST0001"]))
        edited = _score_line_edit(f.read_text())
    elif kind == "sequences":
        f = _write_faa(tmp_path / "a" / "m.faa", three)
        other = _write_faa(tmp_path / "o" / "m.faa", three[:2])
        edited = f.read_text().replace(three[0][1][:5], three[0][1][:4] + ("W" if three[0][1][4] != "W" else "Y"), 1)
    else:
        s = three[0][1]
        f = tmp_path / "a" / "m.sto"
        f.parent.mkdir(parents=True)
        f.write_text(f"# STOCKHOLM 1.0\n\nx1 {s[:60]}\nx2 {s[1:61]}\n//\n")
        other = tmp_path / "o" / "m.sto"
        other.parent.mkdir(parents=True)
        other.write_text(f"# STOCKHOLM 1.0\n\nx1 {s[:60]}\nx2 {s[2:62]}\n//\n")
        edited = f.read_text().replace(f"x2 {s[1]}", f"x2 {'W' if s[1] != 'W' else 'Y'}")
    assert edited != f.read_text()
    moved = tmp_path / "moved" / "elsewhere" / f"renamed{f.suffix}"
    moved.parent.mkdir(parents=True)
    shutil.copyfile(f, moved)
    changed = tmp_path / "changed" / f.name
    changed.parent.mkdir(parents=True)
    changed.write_text(edited)
    tool = "hmmsearch" if kind == "hmm_file" else "jackhmmer"
    a, b, c, o = (_pinned(root, _file_recipe(kind, p, pos, tool=tool)) for p in (f, moved, changed, other))
    assert a[:2] == b[:2]  # the same content at another path: same recipe_hash and query_hash
    assert c[0] != a[0] and c[1] != a[1] and o[0] != a[0] and o[1] != a[1]
    sha = hashlib.sha256(f.read_bytes()).hexdigest()
    form = json.dumps(rcp.identity(a[2]["recipe"]), sort_keys=True)
    assert sha in form and str(tmp_path) not in form and f.name not in form
    assert a[2]["seed"]["sha256"] == sha and a[2]["seed"]["path"] == str(f)


def _dna_hmm(path: Path) -> Path:
    alpha = pyhmmer.easel.Alphabet.dna()
    seq = pyhmmer.easel.TextSequence(name=b"d", sequence="ACGTTGCAACGTAGCTAGCTAGGATCCA").digitize(alpha)
    hmm, _, _ = pyhmmer.plan7.Builder(alpha, seed=42).build(seq, pyhmmer.plan7.Background(alpha))
    return _write_hmms(path, [hmm])


@pytest.mark.covers("W-02")
@pytest.mark.parametrize("case,match", [
    ("missing_hmm", "not found"), ("two_models", "2 models"), ("absent_model", "no model 'NOPE'"),
    ("not_hmm", "not an HMMER3 profile"), ("dna_hmm", "amino-acid profile"), ("missing_faa", "not found"),
    ("one_sequence", "at least 2"), ("dna_letters", "nucleotide"), ("duplicate_names", "duplicate"),
])
def test_w02_seed_files_are_refused_at_submit(world, tmp_path, case, match):
    root, d, _ = world
    pos = sorted(_stored(root, "TEST0000"))[0]
    prots = domains.load_proteins(root, V, "S0")
    p = [prots[u].rstrip("*") for u in sorted(prots)[:3]]
    kind, extra = ("sequences" if case in ("missing_faa", "one_sequence", "dna_letters", "duplicate_names")
                   else "hmm_file"), ""
    two = lambda: _write_hmms(tmp_path / "two.hmm", _library_models(d / "a.hmm", ["TEST0000", "TEST0001"]))  # noqa: E731
    f = {
        "missing_hmm": lambda: tmp_path / "absent.hmm",
        "two_models": two,
        "absent_model": two,
        "not_hmm": lambda: _write_faa(tmp_path / "x.hmm", [("a", p[0]), ("b", p[1])]),
        "dna_hmm": lambda: _dna_hmm(tmp_path / "dna.hmm"),
        "missing_faa": lambda: tmp_path / "absent.faa",
        "one_sequence": lambda: _write_faa(tmp_path / "one.faa", [("a", p[0])]),
        "dna_letters": lambda: _write_faa(tmp_path / "dna.faa", [("a", "ACGTTGCAACGT" * 5), ("b", "ACGTAGCTAGGA" * 5)]),
        "duplicate_names": lambda: _write_faa(tmp_path / "dup.faa", [("a", p[0]), ("b", p[1]), ("a", p[2])]),
    }[case]()
    if case == "absent_model":
        extra = ", model: NOPE"
    before = sorted(runner.runs_dir(root).glob("run=*"))
    with pytest.raises(RecipeError, match=match):
        _run(world, _file_recipe(kind, f, pos, seed_extra=extra))
    assert sorted(runner.runs_dir(root).glob("run=*")) == before  # a refusal records nothing


@pytest.mark.covers("W-02")
def test_w02_model_key_picks_one_model_and_the_run_records_its_provenance(world, tmp_path):
    root, d, _ = world
    lam0, lam1 = _library_models(d / "a.hmm", ["TEST0000", "TEST0001"])
    f = _write_hmms(tmp_path / "lib" / "two.hmm", [lam0, lam1])
    pos = sorted(_stored(root, "TEST0001"))[0]
    by_name = _pinned(root, _file_recipe("hmm_file", f, pos, seed_extra=", model: lam1"))[2]["seed"]
    assert (by_name["model_name"], by_name["model_accession"]) == ("lam1", "TEST0001")
    text = _file_recipe("hmm_file", f, pos, seed_extra=", model: TEST0001")
    row = _run(world, text)
    assert row["status"] == "ok", row["error"]
    seed = json.loads(row["seed_json"])
    assert seed["kind"] == "hmm_file" and seed["path"] == str(f) and seed["model"] == "TEST0001"
    assert seed["sha256"] == hashlib.sha256(f.read_bytes()).hexdigest()
    assert (seed["model_name"], seed["model_accession"], seed["model_length"]) == ("lam1", "TEST0001", lam1.M)
    rh, qh, _ = _pinned(root, text)
    assert (row["recipe_hash"], row["query_hash"]) == (rh, qh)
    hits = runner.read_hits(root, row["run_id"].hex())
    assert pos in {h["target_uid"] for h in hits} and {h["model_len"] for h in hits} == {lam1.M}


def _b_seeds(placed) -> list[tuple[bytes, str, str]]:
    """Three members of family B, as the library model TRAP0001 stored them: (uid, name, sequence)."""
    root = placed[0]
    prots = domains.load_proteins(root, V, "S0")
    return [(u, f"bseed{i}", prots[u]) for i, u in enumerate(sorted(_stored(root, "TRAP0001"))[:3])]


@pytest.mark.covers("W-02")
def test_w02_sequences_seed_recovers_its_family_and_stores_the_alignment(placed, tmp_path):
    root, d, w = placed
    seeds = _b_seeds(placed)
    assert {u for u, _, _ in seeds} <= w["b"]
    faa = _write_faa(tmp_path / "seeds.faa", [(n, s) for _, n, s in seeds])
    row = runner.run(root, _file_recipe("sequences", faa, seeds[0][0], rel="rSyn", bits=30), "rSyn",
                     work_base=d / "work", limits=LIMITS)
    assert row["status"] == "ok", row["error"]
    final = {h["target_uid"] for h in runner.read_hits(root, row["run_id"].hex())}
    assert {u for u, _, _ in seeds} <= final
    assert len(final & w["b"]) >= 25  # three seeds find the family, not only themselves
    prof = [p for p in pq.read_table(runner.run_dir(root, row["run_id"]) / "profiles.parquet").to_pylist()
            if p["round"] == 1]
    assert len(prof) == 1 and prof[0]["n_seq"] == 3
    sto = fastz.unzst(prof[0]["msa_zst"]).decode()
    assert sto.startswith("# STOCKHOLM") and all(n in sto for _, n, _ in seeds)
    msa = pyhmmer.easel.MSAFile(io.BytesIO(sto.encode()), format="stockholm").read()
    assert sorted(engine.text(n) for n in msa.names) == sorted(n for _, n, _ in seeds)
    seed = json.loads(row["seed_json"])
    assert seed["path"] == str(faa) and seed["sha256"] == hashlib.sha256(faa.read_bytes()).hexdigest()
    assert seed["n_sequences"] == 3 and seed["names"] == [n for _, n, _ in seeds]


@pytest.mark.covers("W-02")
@pytest.mark.parametrize("kind", ["hmm_file", "sequences", "msa"])
def test_w02_file_seeds_iterate(placed, tmp_path, kind):
    root, d, w = placed
    seeds = _b_seeds(placed)
    if kind == "hmm_file":
        f = _write_hmms(tmp_path / "trap.hmm", _library_models(w["models"]["trap-1"], ["TRAP0001"]))
    elif kind == "sequences":
        f = _write_faa(tmp_path / "seeds.faa", [(n, s) for _, n, s in seeds])
    else:
        (trap,) = _library_models(w["models"]["trap-1"], ["TRAP0001"])
        block = pyhmmer.easel.TextSequenceBlock(
            pyhmmer.easel.TextSequence(name=n.encode(), sequence=s) for _, n, s in seeds).digitize(trap.alphabet)
        f = tmp_path / "seeds.sto"
        with open(f, "wb") as fh:
            pyhmmer.hmmer.hmmalign(trap, block, trim=True).write(fh, "stockholm")
    row = runner.run(root, _file_recipe(kind, f, seeds[0][0], rel="rSyn", tool="jackhmmer", mx=3, bits=30), "rSyn",
                     work_base=d / "work", limits=LIMITS)
    assert row["status"] == "ok", row["error"]
    profiles = pq.read_table(runner.run_dir(root, row["run_id"]) / "profiles.parquet").to_pylist()
    assert row["rounds"] >= 2 and len(profiles) == row["rounds"]
    assert {u for u, _, _ in seeds} <= {h["target_uid"] for h in runner.read_hits(root, row["run_id"].hex())}


HMMER = {b: shutil.which(b) for b in ("hmmbuild", "hmmalign", "hmmsearch")}


def _tblout(path: Path, bits: float) -> dict[str, float]:
    out = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            f = line.split()
            if float(f[5]) >= bits:
                out[f[0]] = float(f[5])
    return out


@pytest.mark.covers("W-02")
@pytest.mark.skipif(None in HMMER.values(), reason="hmmbuild, hmmalign or hmmsearch binary not on PATH")
def test_w02_sequences_profile_hits_equal_the_hmmer_binaries(placed, tmp_path):
    """AGENTS §5.2 oracle: HMMER's binaries build the profile the walk-through's way (a single-sequence
    profile of the first seed, `hmmalign --trim` of all seeds to it, `hmmbuild`) and `hmmsearch` it at
    the same Z and threshold; the included targets equal the run's. The run's stored round-1 profile,
    searched by `hmmsearch`, also gives the run's sequence scores (the tolerance of the round oracle)."""
    root, d, w = placed
    seeds = _b_seeds(placed)
    bits = 30.0
    faa = _write_faa(tmp_path / "seeds.faa", [(n, s) for _, n, s in seeds])
    row = runner.run(root, _file_recipe("sequences", faa, seeds[0][0], rel="rSyn", bits=bits), "rSyn",
                     work_base=d / "work", limits=LIMITS)
    assert row["status"] == "ok", row["error"]
    got = {h["target_uid"].hex(): h["seq_bitscore"] for h in runner.read_hits(root, row["run_id"].hex())}
    unique = tmp_path / "unique.faa"
    with open(unique, "w") as fh:
        for sh in sorted(pack.pack_dir(root, "rSyn").glob("unique/part-*.faa.zst")):
            for name, seq in fastz.read_records(sh):
                fh.write(f">{name}\n{seq}\n")
    first = tmp_path / "first.sto"
    first.write_text(f"# STOCKHOLM 1.0\n\n{seeds[0][1]} {seeds[0][2]}\n//\n")
    run = lambda *a: subprocess.run([str(x) for x in a], check=True, capture_output=True)  # noqa: E731
    run(HMMER["hmmbuild"], "--amino", "--singlemx", "-n", "first", "-o", "/dev/null", tmp_path / "first.hmm", first)
    run(HMMER["hmmalign"], "--trim", "--amino", "-o", tmp_path / "seeds.sto", tmp_path / "first.hmm", faa)
    run(HMMER["hmmbuild"], "--amino", "-n", "w02", "-o", "/dev/null", tmp_path / "w02.hmm", tmp_path / "seeds.sto")
    stored = tmp_path / "stored.hmm"
    prof = pq.read_table(runner.run_dir(root, row["run_id"]) / "profiles.parquet").to_pylist()
    stored.write_bytes(fastz.unzst(next(p for p in prof if p["round"] == 1)["profile_zst"]))
    tables = {}
    for name in ("w02", "stored"):
        tbl = tmp_path / f"{name}.tbl"
        run(HMMER["hmmsearch"], "--noali", "--cpu", "1", "-Z", "1e8", "--domZ", "1e8", "-T", bits, "--incT", bits,
            "--tblout", tbl, "-o", "/dev/null", tmp_path / f"{name}.hmm" if name == "w02" else stored, unique)
        tables[name] = _tblout(tbl, bits)
    assert set(tables["w02"]) == set(got), "the binaries' profile includes other targets"
    assert set(tables["stored"]) == set(got)
    assert all(abs(got[t] - tables["stored"][t]) <= 0.051 for t in got)


# --- cancel and runner death ------------------------------------------------------------------------


def _wait_for(pred, timeout=60.0):
    import time

    deadline = time.time() + timeout
    while time.time() < deadline:
        if pred():
            return True
        time.sleep(0.1)
    return False


def _alive(pid: int) -> bool:
    import os

    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    try:  # a zombie still answers kill(0); it is gone for our purposes
        st = subprocess.run(["ps", "-o", "stat=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
        return bool(st) and not st.startswith("Z")
    except FileNotFoundError:
        return True


def test_cancel_records_stopped_with_no_hits(world, tmp_path):
    root, d, _ = world
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit))
    work = tmp_path / "work"
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout")
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    runner.cancel(root, rid)
    assert _wait_for(lambda: runner.status(root, rid)["status"] != "running")
    st = runner.status(root, rid)
    assert st["status"] == "stopped", st
    assert sorted(p.name for p in runner.run_dir(root, rid).iterdir()) == ["run.parquet"]


def test_runner_death_is_reported_and_the_worker_does_not_outlive_it(world, tmp_path):
    import os
    import signal

    root, d, _ = world
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit))
    work = tmp_path / "work"
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout")
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    worker_pid = int((work / rid / "worker.pid").read_text())
    runner_pid = json.loads((runner.runs_dir(root) / "pending" / f"{rid}.json").read_text())["pid"]
    os.kill(runner_pid, signal.SIGKILL)
    assert _wait_for(lambda: not _alive(runner_pid), 10)
    st = runner.status(root, rid)
    assert st["status"] == "failed" and "without closing" in st["error"]
    assert _wait_for(lambda: not _alive(worker_pid), 10), "worker outlived its runner"
    assert not runner.run_dir(root, rid).exists()


@pytest.mark.covers("A-01")
def test_status_does_not_trust_a_pid_a_live_unrelated_process_may_own(tmp_path):
    """Pid reuse: the recorded pid answers kill(0), but nobody holds the runner lock: the run is dead."""
    import os

    pend = runner.runs_dir(tmp_path) / "pending"
    pend.mkdir(parents=True)
    lock = tmp_path / "work" / "ab" / "runner.lock"
    lock.parent.mkdir(parents=True)
    lock.write_text("")
    (pend / "ab.json").write_text(json.dumps({"run_id": "ab", "pid": os.getpid(), "lock": str(lock),
                                              "status": "running"}))
    st = runner.status(tmp_path, "ab")
    assert st["status"] == "failed" and "without closing" in st["error"], st
    import fcntl

    with open(lock, "a") as f:  # status must have released the lock it probed with
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)


@pytest.mark.covers("A-01")
def test_cancel_does_not_signal_a_pid_whose_runner_is_gone(tmp_path):
    """The runner died and its pid was reused: cancel must leave the new owner of the pid alone."""
    import sys

    pend = runner.runs_dir(tmp_path) / "pending"
    pend.mkdir(parents=True)
    lock = tmp_path / "work" / "cd" / "runner.lock"
    lock.parent.mkdir(parents=True)
    lock.write_text("")
    bystander = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        (pend / "cd.json").write_text(json.dumps({"run_id": "cd", "pid": bystander.pid, "lock": str(lock),
                                                  "status": "running"}))
        runner.cancel(tmp_path, "cd")
        with pytest.raises(subprocess.TimeoutExpired):
            bystander.wait(timeout=1.0)
    finally:
        bystander.kill()
        bystander.wait()


@pytest.mark.covers("A-01")
def test_runner_lock_is_held_before_the_pending_record_is_visible_and_status_frees_it(world, tmp_path, monkeypatch):
    import fcntl
    import os
    import signal

    def held(path):
        with open(path, "a") as f:
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return True
            return False

    root, d, _ = world
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit))
    work = tmp_path / "work"
    seen = {}
    real_popen = subprocess.Popen

    def spy(cmd, *a, **kw):
        rid = cmd[cmd.index("--run-id") + 1]
        lock = work / rid / "runner.lock"
        seen["lock_exists"] = lock.exists()
        seen["held_at_spawn"] = lock.exists() and held(lock)
        seen["pending_at_spawn"] = (runner.runs_dir(root) / "pending" / f"{rid}.json").exists()
        return real_popen(cmd, *a, **kw)

    monkeypatch.setattr(runner.subprocess, "Popen", spy)
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout")
    monkeypatch.undo()
    assert seen == {"lock_exists": True, "held_at_spawn": True, "pending_at_spawn": False}, seen
    lock = work / rid / "runner.lock"
    info = json.loads((runner.runs_dir(root) / "pending" / f"{rid}.json").read_text())
    assert info["lock"] == str(lock)
    assert runner.status(root, rid)["status"] == "running"
    assert held(lock)
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    runner.cancel(root, rid)
    assert _wait_for(lambda: runner.status(root, rid)["status"] != "running")
    assert runner.status(root, rid)["status"] == "stopped"


# --- W-04: worker threads, exit and kill labelling ------------------------------------------------------


@pytest.mark.covers("W-04", "S-06")
def test_w04_worker_thread_pools_follow_the_runs_cpu_limit(world, record_property):
    import os

    root, d, _ = world
    hit, _ = _seed(root)
    row = runner.run(root, _recipe(hit, hit, tool="jackhmmer", mx=3), REL, work_base=d / "work",
                     limits=runner.Limits(cpus=2, wall_s=120), keep_work=True)
    assert row["status"] == "ok", row["error"]
    t = json.loads((d / "work" / row["run_id"].hex() / "out" / "result.json").read_text())["threads"]
    for k in ("python_peak", "os_peak", "arrow_cpu_count", "arrow_io_thread_count"):
        record_property(f"worker_{k}", t[k])
    record_property("host_cpu_count", os.cpu_count())
    assert t["cpus"] == 2
    assert (t["arrow_cpu_count"], t["arrow_io_thread_count"]) == (2, 2)
    assert t["env"] == {v: "2" for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                          "POLARS_MAX_THREADS", "ARROW_IO_THREADS")}
    assert t["python_peak"] >= 1 and t["samples"] >= row["rounds"]
    if platform.system() == "Linux":
        assert t["os_peak"] is not None and t["os_peak"] >= t["python_peak"]
    if (os.cpu_count() or 1) >= 8:  # a host-sized pool would show here
        assert t["python_peak"] < 16
        if platform.system() == "Linux":  # /proc/self/status counts native threads too; skipped elsewhere
            assert t["os_peak"] < 16, t


@pytest.mark.covers("W-04", "S-05")
@pytest.mark.parametrize("fault,last_line", [
    ("round2", "RuntimeError: injected in round 2"),
    ("round2_thread", "RuntimeError: injected in a worker thread in round 2"),
])
def test_w04_an_exception_in_round_2_ends_the_run_failed_within_seconds(placed, fault, last_line, record_property):
    """In the main thread, or in a thread whose caller then waits on it forever (a pool that lost a
    worker): the process ends at once, non-zero, with its traceback; no hits from round 1 survive."""
    root, d, w = placed
    row = runner.run(root, _trap_recipe(w), "rSyn", work_base=d / "work", limits=runner.Limits(cpus=2, wall_s=60),
                     fault=fault, keep_work=True)
    record_property("run_wall_s", row["wall_s"])
    assert row["status"] == "failed", row["error"]
    assert row["error"].endswith(last_line), row["error"]
    assert row["wall_s"] < 20
    assert sorted(p.name for p in runner.run_dir(root, row["run_id"]).iterdir()) == ["run.parquet"]
    assert runner.read_hits(root, row["run_id"].hex()) == []
    work = d / "work" / row["run_id"].hex()
    assert (work / "out" / "hits.parquet").exists()  # round 1 was written: the fault came after it
    err = (work / "stderr.txt").read_text()
    assert "Traceback" in err and err.strip().splitlines()[-1] == last_line
    assert not _alive(int((work / "worker.pid").read_text()))


@pytest.mark.covers("W-04", "S-05")
def test_w04_a_sigkill_from_outside_is_failed_not_oom(world, tmp_path):
    import os
    import signal

    root, d, _ = world
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit))
    work = tmp_path / "work"
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout")
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    os.kill(int((work / rid / "worker.pid").read_text()), signal.SIGKILL)
    assert _wait_for(lambda: runner.status(root, rid)["status"] != "running")
    st = runner.status(root, rid)
    assert st["status"] == "failed", st
    assert "signal 9" in st["error"] and "outside" in st["error"]
    assert sorted(p.name for p in runner.run_dir(root, rid).iterdir()) == ["run.parquet"]


# --- C-03 and C-08 for search outputs -----------------------------------------------------------------


def test_c03_search_links_resolve_and_a_foreign_run_id_is_caught(world, tmp_path):
    root = world[0]
    hit, _ = _seed(root)
    _run(world, _recipe(hit, hit, tool="jackhmmer", mx=2))
    man = catalog.build_manifest(root, REL)
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    try:
        assert contracts.check_search_links(con) == []
        assert con.execute("SELECT count(*) FROM search_hit").fetchone()[0] > 0
    finally:
        con.close()
    bad = _copy(world, tmp_path)
    hits_file = sorted((bad / "search").glob("run=*/hits.parquet"))[0]
    t = pq.read_table(hits_file).to_pylist()
    for r in t:
        r["run_id"] = b"\x11" * 16
    pq.write_table(pa.Table.from_pylist(t, schema=table("search_hit").schema), hits_file)
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(bad, REL), tmp_path / "d.duckdb"))
    try:
        problems = contracts.check_search_links(con)
    finally:
        con.close()
    assert any(p.startswith("search_hit -> search_run") for p in problems)
    assert any(p.startswith("run_id equals its directory") for p in problems)


def test_c08_search_runs_and_pack_record_the_locked_versions(world):
    from biobox import provenance

    root = world[0]
    assert sorted((root / "search").glob("run=*/run.parquet"))
    assert provenance.check_search(root) == []
    lock = {**provenance.locked_versions(), "pyhmmer": "0.0.0"}
    problems = provenance.check_search(root, lock)
    assert any("pack_manifest.json" in p for p in problems) and any("run.parquet" in p for p in problems)


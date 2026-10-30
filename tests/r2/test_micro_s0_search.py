"""Ring R2: the search pack, runner and exports on micro-S0, with fault injection (C-12, S-04..S-07, S-13).

Micro-S0 is built as in test_micro_s0_l3 (plus IGRs), then the release's pack. Runs use the
real worker under the real limits; faults are injected through the runner's test hook
(`fault=`), and on Linux the OOM test allocates past a real RLIMIT_AS. S-07 runs on a
synthetic world with a planted trap: a seed family A, a library family B that A's seed does
not carry, and bridge sequences between them, so an unguarded jackhmmer drifts into B.
"""

from __future__ import annotations

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

from biobox import catalog, contracts, ids
from biobox.l1 import domains, igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.registry.schemas import table
from biobox.search import export, fastz, orfs, pack, runner
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
    from biobox.l3.tapes import build_tape

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
                        "dist_to_contig_end_bp": 0 if partial and i == 1 else min(st_, 1700 - en)})
        tapes_.append({"contig_uid": cu, "tax_order_bucket": None, "genome_uid": gu, "contig_len": 1700,
                       "is_circular": False, "n_slots": 3, "tier": 1, "assembly_quality": "complete",
                       **build_tape(coords, 1700, False), "protein_uids": prots})
    occ.sort(key=lambda r: (r["protein_uid"], r["contig_uid"], r["slot_ordinal"]))
    tapes_.sort(key=lambda r: r["contig_uid"])
    for rel, rows, name in ((f"l1/protein_occurrence/gene_call_version={V}/stage=S0", occ, "protein_occurrence"),
                            (f"l3/tapes/gene_call_version={V}/stage=S0", tapes_, "contig_tapes")):
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


def test_search_tables_conform_to_the_registry(world):
    root = world[0]
    for p in sorted((root / "search").glob("run=*/*.parquet")):
        name = {"run.parquet": "search_run", "hits.parquet": "search_hit", "profiles.parquet": "search_profile",
                "neighborhoods.parquet": "neighborhoods"}[p.name]
        assert contracts is not None
        from biobox.registry.schemas import conformance_problems

        assert conformance_problems(pq.read_table(p), table(name)) == []


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
    codec = pa.Codec("zstd")
    profiles = pq.read_table(runner.run_dir(root, row["run_id"]) / "profiles.parquet").to_pylist()
    hits = pq.read_table(runner.run_dir(root, row["run_id"]) / "hits.parquet").to_pylist()
    bits = 30.0
    for prof in profiles:
        k = prof["round"]
        hmm = tmp_path / f"r{k}.hmm"
        hmm.write_bytes(codec.decompress(prof["profile_zst"], asbytes=True))
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


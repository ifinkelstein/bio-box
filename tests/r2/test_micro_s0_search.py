"""Ring R2: the search pack, runner and exports on micro-S0, with fault injection (C-12, S-04..S-07, S-13).

Micro-S0 is built as in test_micro_s0_l3 (plus IGRs), then the release's pack, in a pack directory
beside the data root, never under it (card W-03); runs go to a runs directory beside it too. Runs use the
real worker under the real limits; faults are injected through the runner's test hook
(`fault=`), and on Linux the OOM test allocates past a real RLIMIT_AS. S-07 runs on a
synthetic world with a planted trap: a seed family A, a library family B that A's seed does
not carry, and bridge sequences between them, so an unguarded jackhmmer drifts into B.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import platform
import random
import re
import shutil
import signal
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
from biobox.search import engine, export, fastz, nbr, orfs, pack, runner, worker
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
    man = pack.build(root, REL, models={"test-1": d / "a.hmm"}, packs=_packs_of(root))
    return root, d, man


# Card W-03: a data root's pack and runs live beside it, never under it. Each test root `<dir>/<name>` has its
# pack directory `<dir>/<name>.packs` and its runs directory `<dir>/<name>.runs`.


def _packs_of(root: Path) -> Path:
    return Path(root).parent / f"{Path(root).name}.packs"


def _runs_of(root: Path) -> Path:
    return Path(root).parent / f"{Path(root).name}.runs"


def _rdir(root: Path, run_id) -> Path:
    """The directory of a run of the root (`run_id` as bytes or hex)."""
    return runner.run_dir(_runs_of(root), run_id)


def _closed_runs(root: Path) -> list[Path]:
    """The closed runs in the root's runs directory."""
    return sorted(p.parent for p in _runs_of(root).glob("run=*/run.parquet"))


def _at(root: Path) -> dict:
    """The runner's keyword arguments that name the root's pack and runs directories."""
    return {"packs": _packs_of(root), "runs": _runs_of(root)}


def _dirs(root: Path) -> list[str]:
    """The command-line options that name the root's runs and pack directories."""
    return ["--runs", str(_runs_of(root)), "--pack-dir", str(_packs_of(root))]


def _packs(world) -> Path:
    return _packs_of(world[0])


def _runs(world) -> Path:
    return _runs_of(world[0])


def _rd(world, run_id) -> Path:
    return _rdir(world[0], run_id)


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
    return runner.run(root, text, REL, work_base=d / "work", limits=kw.pop("limits", LIMITS), packs=_packs(world),
                      runs=_runs(world), **kw)


# --- C-12 ------------------------------------------------------------------------------------------


def test_c12_pack_verifies_and_counts_equal_sources(world):
    root, d, man = world
    assert pack.verify(root, REL, bytes.fromhex(man["pack_id"]), packs=_packs(world)) == []
    t = man["targets"]
    assert t["unique"]["n_seq"] == len(domains.load_proteins(root, V, "S0"))
    igrs = pq.read_table(next((root / "l0seq/igr").rglob("*.parquet"))).to_pylist()
    assert t["igr"]["n_seq"] == sum(1 for r in igrs if r["end"] - r["start"] >= 30) > 0
    assert t["models"]["n_seq"] == 6 and t["dark_orf"]["n_seq"] > 0
    assert "spacers" in man["missing_targets"]


def _copy(world, tmp_path) -> Path:
    """A copy of the world's root, with copies of its pack and runs directories beside it (card W-03: they are
    not under the root, so copying the root alone no longer takes them)."""
    root = world[0]
    dst = tmp_path / "a"
    shutil.copytree(root, dst)
    for src, to in ((_packs_of(root), _packs_of(dst)), (_runs_of(root), _runs_of(dst))):
        if src.exists():
            shutil.copytree(src, to)
    return dst


@pytest.mark.parametrize("fault", ["byte", "stray", "missing_meta", "count", "pack_id"])
def test_c12_pack_faults_are_caught(world, tmp_path, fault):
    root = _copy(world, tmp_path)
    packs = _packs_of(root)
    pdir = pack.find(packs, REL)
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
    assert pack.verify(root, REL, pid, packs=packs)


def _pack_bytes(packs: Path) -> dict[str, bytes]:
    return {p.relative_to(packs).as_posix(): p.read_bytes() for p in sorted(packs.rglob("*")) if p.is_file()}


def test_c12_small_shards_fail_unless_last_and_pack_is_deterministic(world, tmp_path):
    # card W-03: a pack is never overwritten, so the second build goes into a pack directory of its own (a
    # build into the same one is a no-op: test_w03_pack_build_twice_is_a_no_op_and_a_pack_is_never_overwritten)
    root = _copy(world, tmp_path)
    man = pack.build(root, REL, models={"test-1": world[1] / "a.hmm"}, shard_text_bytes=20_000, packs=tmp_path / "p1")
    assert sum(1 for s in man["shards"] if s["target"] == "unique") > 1
    assert any("not the last shard" in p for p in pack.verify(root, REL, packs=tmp_path / "p1"))
    assert pack.verify(root, REL, min_shard_bytes=0, packs=tmp_path / "p1") == []
    again = pack.build(root, REL, models={"test-1": world[1] / "a.hmm"}, shard_text_bytes=20_000, packs=tmp_path / "p2")
    assert again["pack_id"] == man["pack_id"]  # D-01: no timestamps in the pack
    assert _pack_bytes(tmp_path / "p1") == _pack_bytes(tmp_path / "p2")  # byte for byte, every file


def test_pack_is_outside_the_release_manifest(world):
    man = catalog.build_manifest(world[0], REL, ["S0"])
    assert not any(f["path"].startswith(("search_pack/", "search/")) for f in man["files"])


# --- W-16: the streaming build gives the bytes of the whole-release build -----------------------------------


def _w16_oracle(root: Path, release: str, producers: list[str] | None = None) -> dict[str, tuple[list, list | None]]:
    """The pack's records and meta rows by the rules of the module docstring, computed with every table of the
    release held in memory (the shape of the build before card W-16, written here separately): per target the
    FASTA records (name, sequence) and, for the regions, the `search_region` rows, both in pack order."""
    from collections import Counter, defaultdict

    from biobox.l0seq.store import SeqStore
    from biobox.l3.geometry import junction

    rf = catalog.ReleaseFiles.of(root, release)
    producers = sorted(producers) if producers else [V]
    tv, iv = tapes.tape_version(producers), igr.igr_version(producers)
    lengths, circular, strands = {}, set(), {}
    for f in rf.table("contig_tapes", tape_version=tv):
        for r in pq.read_table(f).to_pylist():
            lengths[r["contig_uid"]] = r["contig_len"]
            if r["is_circular"]:
                circular.add(r["contig_uid"])
            strands.update(zip(r["feature_uids"], r["strands"]))
    genes = defaultdict(list)
    for f in rf.table("protein_occurrence", gene_call_version=V):
        for r in pq.read_table(f).to_pylist():
            genes[r["contig_uid"]].append((r["start"], r["end"], r["strand"]))
    counts = defaultdict(Counter)
    for f in rf.table("feature_protein", gene_call_version=V):
        for r in pq.read_table(f).to_pylist():
            counts[r["contig_uid"]][r["genetic_code"]] += 1
    prots = {}
    for f in rf.table("protein_seq", gene_call_version=V):
        for r in pq.read_table(f).to_pylist():
            prots[r["protein_uid"]] = r["seq"]

    def region(kind, cu, s, e, strand, n, **kw):
        row = {"region_uid": ids.region_uid(cu, s, e, strand), "target_kind": kind, "contig_uid": cu, "start": s,
               "end": e, "strand": strand, "n_residues": n, "orientation": None, "gc": None,
               "dist_to_contig_end_bp": None, "partial": False, "at_contig_end": False, "gene_overlap_frac": None,
               "genetic_code": None, "left_feature_uid": None, "right_feature_uid": None}
        row.update(kw)
        return row

    def gc(seq: bytes) -> float:
        s = seq.upper()
        acgt = sum(s.count(b) for b in (b"A", b"C", b"G", b"T"))
        return (s.count(b"G") + s.count(b"C")) / acgt if acgt else 0.0

    store = SeqStore(root, rf)
    igr_rows, dark_rows = [], []
    for f in rf.table("igr", igr_version=iv):
        for r in pq.read_table(f, partitioning=None).to_pylist():
            cu, s, e = r["contig_uid"], r["start"], r["end"]
            if e - s < pack.IGR_MIN_NT:
                continue
            circ, L = cu in circular, lengths[cu]
            seq = store.get(cu, s, e, 1, circular=circ)
            lf, rt = r["left_feature_uid"], r["right_feature_uid"]
            igr_rows.append((region("igr", cu, s, e, 1, e - s, left_feature_uid=lf, right_feature_uid=rt,
                                    orientation=junction((0, 1, strands[lf]), (2, 3, strands[rt])) if lf and rt else None,
                                    gc=gc(seq), dist_to_contig_end_bp=None if circ else min(s, L - e),
                                    partial=not circ and (s == 0 or e == L)), seq.decode()))
    for cu in sorted(lengths):
        code = counts[cu].most_common(1)[0][0] if cu in counts else 11
        for o, frac, at_end in orfs.dark_orfs(store.get(cu), cu in circular, sorted(genes.get(cu, [])), code):
            dark_rows.append((region("dark_orf", cu, o.start, o.end, o.strand, len(o.aa), partial=o.partial,
                                     at_contig_end=at_end, gene_overlap_frac=frac, genetic_code=code), o.aa))
    store.close()
    out = {"unique": (sorted((u.hex(), s.replace("*", "")) for u, s in prots.items()), None)}
    for t, rows in (("igr", igr_rows), ("dark_orf", dark_rows)):
        rows.sort(key=lambda x: (x[0]["region_uid"], x[1]))
        out[t] = ([(m["region_uid"].hex(), s) for m, s in rows], [m for m, _ in rows])
    return out


def _w16_expected_files(oracle: dict, d: Path, shard_bytes: int, rg_rows: int) -> dict[str, bytes]:
    """The oracle written with the pack's two file writers: FASTA shards by `fastz.write_records` and the meta
    rows by one `pq.write_table` of the whole table."""
    suffix = {"unique": ".faa.zst", "igr": ".fna.zst", "dark_orf": ".faa.zst"}
    for t, (records, meta) in oracle.items():
        fastz.write_records(d / t, suffix[t], records, shard_bytes)
        if meta is not None:
            pq.write_table(pa.Table.from_pylist(meta, schema=table("search_region").schema), d / f"{t}/{t}_meta.parquet",
                           compression="zstd", row_group_size=rg_rows)
    return {p.relative_to(d).as_posix(): p.read_bytes() for p in sorted(d.rglob("*")) if p.is_file()}


def _w16_pack_files(pdir: Path) -> dict[str, bytes]:
    return {p.relative_to(pdir).as_posix(): p.read_bytes() for p in sorted(pdir.rglob("*"))
            if p.is_file() and p.parts[len(pdir.parts)] in ("unique", "igr", "dark_orf")}


@pytest.mark.covers("W-16", "D-01")
def test_w16_fixture_pack_equals_the_in_memory_oracle(world, tmp_path):
    root, d, man = world
    got = _w16_pack_files(pack.find(_packs(world), REL))
    want = _w16_expected_files(_w16_oracle(root, REL), tmp_path / "want", pack.SHARD_TEXT_BYTES, 122_880)
    assert sorted(got) == sorted(want)
    assert {k: v == want[k] for k, v in got.items()} == {k: True for k in want}
    assert man["targets"]["igr"]["n_seq"] and man["targets"]["dark_orf"]["n_seq"] and man["targets"]["unique"]["n_seq"]


@pytest.mark.covers("W-16", "MR-20")
def test_w16_two_producer_pack_equals_the_in_memory_oracle(two, tmp_path):
    root, _, man, _, _ = two
    got = _w16_pack_files(pack.pack_dir(_packs_of(root), TWO, man["pack_id"]))
    want = _w16_expected_files(_w16_oracle(root, TWO, [V, NEW]), tmp_path / "want", pack.SHARD_TEXT_BYTES, 122_880)
    assert {k: v == want.get(k) for k, v in got.items()} == {k: True for k in want}


@pytest.mark.covers("W-16", "D-01")
def test_w16_many_runs_shards_and_row_groups_give_the_oracles_bytes(world, tmp_path, monkeypatch):
    """Caps far below the fixture's size, so every part goes through several sorted runs and merge passes, input
    batches, shards and meta row groups: the files still equal the whole-table oracle written with the same shard
    size and row group size, and the pack verifies."""
    root, d, _ = world
    monkeypatch.setattr(pack, "RUN_ROWS", 50)
    monkeypatch.setattr(pack, "FAN_IN", 3)  # several merge passes
    monkeypatch.setattr(pack, "MERGE_BATCH", 16)
    monkeypatch.setattr(pack, "META_ROW_GROUP_ROWS", 64)
    monkeypatch.setattr(pack, "BATCH_ROWS", 7)
    monkeypatch.setattr(pack, "TAPE_BATCH_ROWS", 1)
    man = pack.build(root, REL, models={"test-1": d / "a.hmm"}, shard_text_bytes=20_000, packs=tmp_path / "p")
    pdir = pack.pack_dir(tmp_path / "p", REL, man["pack_id"])
    assert pq.ParquetFile(pdir / "dark_orf/dark_orf_meta.parquet").metadata.num_row_groups > 2
    assert sum(1 for s in man["shards"] if s["target"] == "dark_orf") > 1
    got = _w16_pack_files(pdir)
    want = _w16_expected_files(_w16_oracle(root, REL), tmp_path / "want", 20_000, 64)
    assert sorted(got) == sorted(want)
    assert {k: v == want[k] for k, v in got.items()} == {k: True for k in want}
    assert pack.verify(root, REL, packs=tmp_path / "p", min_shard_bytes=0) == []
    assert [p.name for p in pack.release_dir(tmp_path / "p", REL).iterdir()] == [pdir.name]  # no scratch left


_W16_STAMP = re.compile(r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d[+-]\d\d:\d\d pack build: ")


@pytest.mark.covers("W-16", "C-13")
def test_w16_progress_lines_one_per_part_and_per_shard_on_stderr(world, tmp_path, capsys):
    root, d, _ = world
    man = pack.build(root, REL, models={"test-1": d / "a.hmm"}, shard_text_bytes=20_000, packs=tmp_path / "p")
    cap = capsys.readouterr()
    assert cap.out == ""
    lines = cap.err.splitlines()
    assert lines and all(_W16_STAMP.match(x) for x in lines), [x for x in lines if not _W16_STAMP.match(x)]
    msgs = [_W16_STAMP.sub("", x) for x in lines]
    for t in ("unique", "igr", "dark_orf", "models"):
        assert sum(m.startswith(f"{t}: start") for m in msgs) == 1, t
    for t in ("unique", "igr", "dark_orf"):
        assert sum(m.startswith(f"{t}: done") for m in msgs) == 1, t
    for s in man["shards"]:
        t, name = s["target"], Path(s["path"]).name
        want = f"{t}: shard {name} written: {s['n_seq']:,} sequences" if t != "models" else f"models: {name} written"
        assert sum(m.startswith(want) for m in msgs) == 1, want
    assert msgs[-1].startswith(f"release {REL}: pack {man['pack_id']} written")


_W16_HELD = r"""
import signal, sys
from pathlib import Path
from biobox.search import pack
root, packs, release, models, flag, point = sys.argv[1:7]
log = pack._log


def held(msg):  # the build stops here until a signal comes
    log(msg)
    if msg.startswith(point):
        Path(flag).write_text(msg)
        signal.pause()


pack._log = held
pack.build(Path(root), release, models={"test-1": Path(models)}, packs=Path(packs))
"""


@pytest.mark.parametrize("point", ["igr: sorting", "dark_orf: start"])
@pytest.mark.covers("W-16")
def test_w16_sigterm_removes_the_building_and_scratch_directories_and_the_next_build_starts_again(
        world, tmp_path, point):
    root, d, man = world
    packs, flag = tmp_path / "packs", tmp_path / "held"
    p = subprocess.Popen([sys.executable, "-c", _W16_HELD, str(root), str(packs), REL, str(d / "a.hmm"), str(flag),
                          point], stderr=subprocess.PIPE, text=True)
    try:
        assert _wait_for(flag.exists, timeout=120), "the build never reached " + point
        base = pack.release_dir(packs, REL)
        assert sorted(x.name.split("-")[0] for x in base.iterdir()) == [".building", ".scratch"]
        p.send_signal(signal.SIGTERM)
        _, err = p.communicate(timeout=60)
    finally:
        if p.poll() is None:
            p.kill()
    assert p.returncode == 143, err
    assert "stopped by a signal (exit 143)" in err and "the next call starts again" in err, err
    assert list(base.iterdir()) == []  # neither directory is left, and no pack
    again = pack.build(root, REL, models={"test-1": d / "a.hmm"}, packs=packs)
    assert again["pack_id"] == man["pack_id"]


@pytest.mark.covers("W-16", "W-03")
def test_w16_a_missing_pack_names_the_build_command_with_model_optional(world, tmp_path):
    with pytest.raises(pack.NoPack) as e:
        pack.find(tmp_path / "none", REL, root=world[0])
    msg = str(e.value)
    cmd = f"bio-box pack build --root {world[0]} --release {REL} --pack-dir {tmp_path / 'none'} (add --model"
    assert cmd in msg and "if it has any" in msg, msg


def test_dark_orfs_translate_their_own_coordinates(world):
    """Each dark ORF's protein equals the translation of its slice from the sequence store."""
    from biobox.l0seq.store import SeqStore

    root = world[0]
    pdir = pack.find(_packs(world), REL)
    meta = {r["region_uid"].hex(): r for r in pq.read_table(pdir / "dark_orf/dark_orf_meta.parquet").to_pylist()}
    seqs = dict(r for s in sorted(pdir.glob("dark_orf/part-*.faa.zst")) for r in fastz.read_records(s))
    store = SeqStore(root, ["S0"])
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
    """Rewrite the genetic code of every CDS feature of the contig whose uid starts with the hex prefix."""
    target = None
    for f in sorted((root / "l1/feature_protein").rglob("*.parquet")):
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
    for f in pack._occurrence_files(catalog.ReleaseFiles.of(root, REL), V):
        genes += [(r["start"], r["end"], r["strand"]) for r in pq.read_table(f).to_pylist() if r["contig_uid"] == cu]
    store = SeqStore(root, ["S0"])
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
    # card W-03: the copy's pack is built into a pack directory of its own (the world's pack is not the copy's)
    man = pack.build(root, REL, models={"test-1": world[1] / "a.hmm"}, packs=tmp_path / "packs")
    pdir = pack.pack_dir(tmp_path / "packs", REL, man["pack_id"])
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
    packs = tmp_path / "packs"  # card W-03: the copy has no pack of its own; none may appear
    with pytest.raises(ValueError) as e:
        pack.build(root, REL, models={"test-1": world[1] / "a.hmm"}, packs=packs)
    assert str(bad) in str(e.value) and cu.hex() in str(e.value)
    assert not packs.exists()


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
    assert {h["target_uid"] for h in runner.read_hits(_rdir(root, row["run_id"]))} >= stored


# --- provenance and failure states (S-05), envelope (S-06) -----------------------------------------------


_S05_MEM_GB = 4 if platform.system() == "Linux" else None
# Each fault gets the limits its state is about (card A-07). The wall clock is the subject of the
# `timeout` fault only: it sleeps forever, so the wall must fire, and 3 s keeps the test short. For
# `oom` and `crash` a wall the fault can reach turns the test into a race between the fault and the
# clock: under a systemd scope (GitHub's runner) the kernel kills only after the worker has imported
# everything and touched mem_gb of pages, 4.5-5.2 s for 4 GB in a scope on the pod (measured
# 2026-11-10); on GitHub it timed out in 4 of the 7 failed runs card A-07 lists. Their work is bounded by construction
# (one allocation of twice the limit, which ends at the limit; one raise), so they run under the
# file's wall (LIMITS.wall_s), which only a hung worker reaches. The assertions are unchanged.
S05_LIMITS = {"timeout": runner.Limits(cpus=1, wall_s=3, mem_gb=_S05_MEM_GB),
              "oom": runner.Limits(cpus=1, wall_s=LIMITS.wall_s, mem_gb=_S05_MEM_GB),
              "crash": runner.Limits(cpus=1, wall_s=LIMITS.wall_s, mem_gb=_S05_MEM_GB)}


@pytest.mark.parametrize("fault,status", [("oom", "oom"), ("timeout", "timeout"), ("crash", "failed")])
def test_s05_injected_failures_are_states_never_zero_hits(world, fault, status):
    hit, _ = _seed(world[0])
    row = _run(world, _recipe(hit, hit), fault=fault, limits=S05_LIMITS[fault])
    # which path limited the worker (scope, inherited cap, RLIMIT_AS) is in argv and limits_json
    assert row["status"] == status, (row["error"], row["argv"][:4], row["limits_json"])
    d = _rdir(world[0], row["run_id"])
    assert sorted(p.name for p in d.iterdir()) == ["run.parquet"]  # no hits file, so never "0 hits"
    assert runner.read_hits(_rdir(world[0], row["run_id"])) == []


@pytest.mark.covers("A-07", "S-05")
def test_a07_oom_is_oom_however_long_the_fault_takes_to_reach_its_limit(world, monkeypatch):
    """CI flake (1): under a systemd scope the kernel kills only once the worker has touched mem_gb of
    pages (4.5-5.2 s for 4 GB on the pod, measured), so an oom test whose wall limit the fault can
    reach is a race. Injected condition: the fault starts 4 s late, longer than the old 3 s wall."""
    monkeypatch.setenv(worker.OOM_DELAY_VAR, "4")
    hit, _ = _seed(world[0])
    row = _run(world, _recipe(hit, hit), fault="oom", limits=S05_LIMITS["oom"])
    assert row["status"] == "oom", (row["error"], row["argv"][:4], row["limits_json"])
    assert sorted(p.name for p in _rdir(world[0], row["run_id"]).iterdir()) == ["run.parquet"]


def test_s05_every_hit_resolves_to_a_run_and_catalog_shows_runs(world):
    root = world[0]
    hit, _ = _seed(root)
    _run(world, _recipe(hit, hit, tool="jackhmmer", mx=3))
    man = catalog.build_manifest(root, REL, ["S0"])
    # card W-03: the runs and the pack are given to the catalog as directories outside the root
    con = catalog.open_catalog(catalog.generate_catalog(man, world[1] / "c.duckdb"), runs=_closed_runs(root),
                               pack=pack.find(_packs(world), REL))
    try:
        orphans = con.execute("SELECT count(*) FROM search_hit h ANTI JOIN search_run r USING (run_id)").fetchone()[0]
        assert orphans == 0
        n_hits = con.execute("SELECT count(*) FROM search_hit").fetchone()[0]
        assert n_hits == sum(pq.read_metadata(p).num_rows for p in _runs_of(root).glob("run=*/hits.parquet"))
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
    assert (_rdir(world[0], row["run_id"]) / "hits.parquet").exists()


def test_dark_orf_target_runs(world):
    root = world[0]
    pdir = pack.find(_packs_of(root), REL)
    name, aa = next(fastz.read_records(next(pdir.glob("dark_orf/part-*.faa.zst"))))
    text = (f"recipe: dark\nwritten_against: {REL}\nseed: [{{sequence: {aa}}}]\ntarget: dark_orf\ntool: phmmer\n"
            f"iterations: {{max: 1, inclusion_bits: 20}}\ncontrols: {{positive: [{{region_uid: {name}}}]}}\n")
    row = _run(world, text)
    assert row["status"] == "ok", row["error"]
    hits = runner.read_hits(_rdir(root, row["run_id"]))
    assert {h["target_kind"] for h in hits} == {"dark_orf"}


# --- exports and CLI (S-13) ---------------------------------------------------------------------------------


def test_s13_cli_json_records_match_the_locked_schemas(world, tmp_path):
    root = world[0]
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit, tool="jackhmmer", mx=2))
    cmd = [sys.executable, "-m", "biobox.cli"]
    out = subprocess.run(cmd + ["search", "run", str(rec), "--root", str(root), *_dirs(root), "--release", REL,
                                "--wait", "--json", "--cpus", "1", "--work", str(world[1] / "work")],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    run = json.loads(out[0])
    locked = json.loads(export.lock_path(REPO).read_text())
    assert run["schema"] == "search_run/1" and list(run) == locked["search_run"]["keys"]
    lines = subprocess.run(cmd + ["search", "hits", run["run_id"], "--root", str(root), *_dirs(root), "--json"],
                           capture_output=True, text=True, check=True).stdout.splitlines()
    hits = [json.loads(x) for x in lines]
    assert hits and all(list(h) == locked["search_hit"]["keys"] for h in hits)
    assert all(h["n_occurrences"] >= 1 for h in hits)
    lines = subprocess.run(cmd + ["nbr", "--anchors", run["run_id"], "--root", str(root), *_dirs(root), "--window", "3g",
                                  "--json"], capture_output=True, text=True, check=True).stdout.splitlines()
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
            # card B-17 renamed slot_ordinal to ordinal; card W-05 (search_hit/2) joins genome_uid and tier back
            # from contig_dim at export time
            assert set(o) == {"contig_uid", "ordinal", "feature_uid", "start", "end", "strand", "genome_uid", "tier"}
            assert tape["feature_uids"][o["ordinal"]].hex() == o["feature_uid"]
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
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1), world[1] / "work", **_at(root))
    deadline = time.time() + 60
    while runner.status(_rdir(root, rid))["status"] == "running" and time.time() < deadline:
        time.sleep(0.2)
    assert runner.status(_rdir(root, rid))["status"] == "ok"


def test_neighborhoods_written_for_final_hits(world):
    root = world[0]
    hit, _ = _seed(root)
    row = _run(world, _recipe(hit, hit, extra="context: {window: 5g, edge: truncate}\n"))
    nb = pq.read_table(_rdir(root, row["run_id"]) / "neighborhoods.parquet").to_pylist()
    final = {h["target_uid"] for h in runner.read_hits(_rdir(root, row["run_id"]))}
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
    t = pq.read_table(_rdir(root, row["run_id"]) / "neighborhoods.parquet")
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
    pdir = pack.find(_packs_of(root), REL)
    n_flanks = 0
    for target in ("igr", "dark_orf"):
        t = pq.read_table(pdir / f"{target}/{target}_meta.parquet")
        assert conformance_problems(t, table("search_region")) == []
        assert not {"left_slot", "right_slot"} & set(t.column_names)  # card B-17: flanks by feature uid only
        for m in t.to_pylist():
            for side in ("left", "right"):
                fu = m[f"{side}_feature_uid"]
                if fu is not None:
                    n_flanks += 1
                    assert feats[fu]["contig_uid"] == m["contig_uid"]
                    assert fu in tape[m["contig_uid"]]["feature_uids"]  # a slot of the pack's tape version
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
    catalog.define_release(root, TWO, ["S0"])  # the pack reads the release's stages (card B-28)
    man = pack.build(root, TWO, models={"test-1": world[1] / "a.hmm"}, producers=producers, packs=_packs_of(root))
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
    files = catalog.ReleaseFiles.of(root, TWO)
    one, _ = nbr.fetch(root, uids, policy, files=files)  # default: the gene-call producer's tape version
    other, _ = nbr.fetch(root, uids, policy, tape_version=tv2, files=files)
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
    before = sorted(_runs_of(root).glob("run=*"))
    with pytest.raises(RecipeError, match=f"tape version {tv1}.*{tv2}"):
        runner.run(root, text, TWO, work_base=d / "work", limits=LIMITS, tape_version=tv1, **_at(root))
    assert sorted(_runs_of(root).glob("run=*")) == before  # a refusal records nothing
    # the pack's own tape version is the run's: neighborhoods and hit occurrences use its ordinals
    row = runner.run(root, text, TWO, work_base=d / "work", limits=LIMITS, **_at(root))
    assert row["status"] == "ok", row["error"]
    nb = pq.read_table(_rdir(root, row["run_id"]) / "neighborhoods.parquet").to_pylist()
    _check_neighborhood_rows(root, nb, tv2)
    assert cli.main(["search", "hits", row["run_id"].hex(), "--root", str(root), *_dirs(root), "--json"]) == 0
    hits = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    tps = _tape_rows(root, tv2)
    assert hits and all(h["tape_version"] == tv2 for h in hits)
    for o in (o for h in hits for o in h["occurrences"]):
        assert tps[bytes.fromhex(o["contig_uid"])]["feature_uids"][o["ordinal"]].hex() == o["feature_uid"]
    assert cli.main(["nbr", "--anchors", row["run_id"].hex(), "--root", str(root), *_dirs(root), "--json",
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
    rows, _ = nbr.fetch(root, anchors, nbr.Policy.parse(window=window), files=catalog.ReleaseFiles.of(root, REL))
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
    man = catalog.build_manifest(root, REL, ["S0"])
    assert nbr.summarize(root, anchors, nbr.Policy.parse(window="5g"), by=("pfam", "nr50", "nr90"), release=REL,
                         top_n=0, manifest=man) == got


@pytest.mark.covers("W-06", "S-13")
def test_w06_cli_summary_json_matches_the_lock_and_tsv_reads_back(world, capsys):
    root = world[0]
    hit, _ = _seed(root)
    row = _run(world, _recipe(hit, hit, tool="jackhmmer", mx=2))
    assert row["status"] == "ok", row["error"]
    rid = row["run_id"].hex()
    args = ["nbr", "summarize", rid, "--root", str(root), *_dirs(root), "--window", "5g", "--by", "pfam,nr50"]
    assert cli.main(args + ["--json"]) == 0
    recs = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    locked = json.loads(export.lock_path(REPO).read_text())
    assert recs[0]["schema"] == "nbr_summary_header/1" and list(recs[0]) == locked["nbr_summary_header"]["keys"]
    assert len(recs) > 3 and all(list(r) == locked["nbr_summary"]["keys"] for r in recs[1:])
    assert {"nbr_summary", "nbr_summary_header"} <= set(locked)
    assert export.lock_problems(locked) == [] and locked == json.loads(json.dumps(export.SCHEMAS, sort_keys=True))
    # the run's hits are the anchors, its release and tape version the summary's
    anchors = {h["target_uid"] for h in runner.read_hits(_rdir(root, rid)) if h["target_kind"] == "protein"}
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
    assert cli.main(["nbr", "summarize", uids, "--root", str(root), *_dirs(root), "--window", "5g", "--release", REL,
                     "--json"]) == 0
    assert [json.loads(x) for x in capsys.readouterr().out.splitlines()] == recs
    assert cli.main(["nbr", "summarize", uids, "--root", str(root), *_dirs(root), "--window", "5g"]) == 2
    assert "--release" in capsys.readouterr().err


# --- W-17: a release built before card B-27 has no member_cluster; readers must not answer from an empty view ---


def _pre_b27(world, tmp_path, *, manifest: bool, labels: bool = False) -> Path:
    """A copy of the world's root as a release built before card B-27 leaves it: no `member_cluster` files (and,
    with `labels`, no `contig_label` files), listed by a manifest of the remaining files or, without one, by the
    release_def stages."""
    root = _copy(world, tmp_path)
    shutil.rmtree(root / "l2/member_cluster")
    if labels:
        shutil.rmtree(root / "l0seq/contig_label")
    assert (root / f"l2/cluster_member/release={REL}").is_dir()
    catalog.manifest_path(root, REL).unlink(missing_ok=True)
    if manifest:
        catalog.write_manifest(catalog.build_manifest(root, REL, ["S0"]))
        assert not [f for f in catalog.read_manifest(catalog.manifest_path(root, REL))["files"]
                    if f["table"] == "member_cluster"]
    return root


@pytest.mark.covers("W-17", "W-06")
@pytest.mark.parametrize("manifest", [False, True], ids=["release_def", "manifest"])
def test_w17_summary_of_a_release_without_member_cluster_equals_the_brute_force_loop(world, tmp_path, capsys, manifest):
    """The answer is the loop over the tapes and `cluster_member` (the oracle of W-06), not all `unannotated`; the
    slow path says so once on stderr, and a release with the table says nothing."""
    root = _pre_b27(world, tmp_path, manifest=manifest)
    anchors = _w06_anchors(root)
    want = _oracle_summary(root, anchors, "5g")
    assert [r for r in want["rows"] if r["level"] == "nr50" and not r["unannotated"]], "the world has no cluster"
    capsys.readouterr()
    got = nbr.summarize(root, anchors, nbr.Policy.parse(window="5g"), by=("pfam", "nr50"), release=REL, top_n=0,
                        manifest=catalog.read_manifest(catalog.manifest_path(root, REL)) if manifest else None)
    err = capsys.readouterr().err.splitlines()
    _compare(got, want)
    assert got["header"]["n_keys"]["nr50"] > 0
    assert len(err) == 1 and "member_cluster" in err[0] and "slow" in err[0], err
    nbr.summarize(world[0], anchors, nbr.Policy.parse(window="5g"), by=("pfam", "nr50"), release=REL, top_n=0)
    assert capsys.readouterr().err == ""


@pytest.mark.covers("W-17")
def test_w17_cli_summary_without_member_cluster_prints_what_the_release_with_it_prints(world, tmp_path, capsys):
    (tmp_path / "x").mkdir()
    root = _pre_b27(world, tmp_path / "x", manifest=True)
    uids = ",".join(sorted(u.hex() for u in _w06_anchors(root)))
    for fmt in (["--json"], ["--format", "tsv"]):
        args = ["nbr", "summarize", uids, "--window", "5g", "--by", "pfam,nr50", "--release", REL, *fmt]
        assert cli.main([*args, "--root", str(world[0]), *_dirs(world[0])]) == 0
        want = capsys.readouterr().out
        assert cli.main([*args, "--root", str(root), *_dirs(root)]) == 0
        out = capsys.readouterr()
        assert out.out == want and "member_cluster" in out.err


@pytest.mark.covers("W-17")
def test_w17_export_annotations_without_member_cluster_read_cluster_member(world, tmp_path, capsys):
    """`search hits` and `nbr --format tsv` take `nr90`/`nr50` from `member_cluster`: without it they read
    `cluster_member`, never an empty relation."""
    root = _pre_b27(world, tmp_path, manifest=True, labels=True)
    sqls = ("SELECT member_uid, level, cluster_uid, rep_uid FROM clus ORDER BY level, member_uid",
            "SELECT contig_uid, arg_min(accession, (l0_object_sha256, record_ordinal)), "
            "arg_min(original_header, (l0_object_sha256, record_ordinal)), sum(n_records)::BIGINT FROM lab "
            "GROUP BY contig_uid ORDER BY contig_uid")
    rows = {}
    for r in (world[0], root):
        with export.open_store(r, release=REL, gene_call_version=V, relations=("clus", "lab")) as st:
            rows[r] = [st.con.execute(sql).fetchall() for sql in sqls]
    assert rows[root][0] and rows[root][1] and rows[root] == rows[world[0]]
    err = capsys.readouterr().err
    assert "member_cluster" in err and "contig_label" in err


@pytest.mark.covers("W-17")
def test_w17_release_verify_warns_for_each_lookup_table_the_release_lacks(world, tmp_path, capsys):
    root = _pre_b27(world, tmp_path, manifest=True, labels=True)
    assert cli.main(["release", "verify", "--root", str(root), "--release", REL]) == 0  # a warning, not a failure
    err = capsys.readouterr().err.splitlines()
    assert sorted(t for t in ("member_cluster", "contig_label") if sum(t in x and x.startswith("warning") for x in err) == 1) \
        == ["contig_label", "member_cluster"], err
    # a release that has both says nothing
    (tmp_path / "g").mkdir()
    good = _copy(world, tmp_path / "g")
    catalog.manifest_path(good, REL).unlink(missing_ok=True)
    catalog.write_manifest(catalog.build_manifest(good, REL, ["S0"]))
    assert cli.main(["release", "verify", "--root", str(good), "--release", REL]) == 0
    assert capsys.readouterr().err == ""


@pytest.mark.covers("W-18")
def test_w18_nbr_takes_anchors_from_a_file_one_uid_per_line(world, capsys, tmp_path):
    """Card W-18 (3): `--anchors-file` gives the same answer as the comma-separated argument (1,957 uids made a 64 KB
    argument on rN.1); blank lines are skipped; a file together with another anchor source, a missing file and a line
    that is not a uid are refused with exit 2."""
    root = world[0]
    hit, _ = _seed(root)
    row = _run(world, _recipe(hit, hit, tool="jackhmmer", mx=2))
    assert row["status"] == "ok", row["error"]
    anchors = sorted(h["target_uid"].hex() for h in runner.read_hits(_rdir(root, row["run_id"].hex()))
                     if h["target_kind"] == "protein")
    f = tmp_path / "anchors.txt"
    f.write_text("\n".join(anchors[:1] + [""] + anchors[1:]) + "\n")
    where = ["--root", str(root), *_dirs(root), "--window", "5g", "--release", REL, "--json"]
    assert cli.main(["nbr", "summarize", ",".join(anchors), *where]) == 0
    want = capsys.readouterr().out
    assert want
    assert cli.main(["nbr", "summarize", "--anchors-file", str(f), *where]) == 0
    assert capsys.readouterr().out == want
    assert cli.main(["nbr", "--anchors-file", str(f), *where]) == 0  # the plain `nbr`, counts and rows
    plain = capsys.readouterr().out
    assert cli.main(["nbr", "--anchors", ",".join(anchors), *where]) == 0
    assert capsys.readouterr().out == plain
    assert cli.main(["nbr", "summarize", anchors[0], "--anchors-file", str(f), *where]) == 2
    assert "--anchors-file" in capsys.readouterr().err
    assert cli.main(["nbr", "summarize", "--anchors-file", str(tmp_path / "none.txt"), *where]) == 2
    assert "none.txt" in capsys.readouterr().err
    (tmp_path / "bad.txt").write_text(anchors[0] + "\nnot-a-uid\n")
    assert cli.main(["nbr", "summarize", "--anchors-file", str(tmp_path / "bad.txt"), *where]) == 2
    assert "not-a-uid" in capsys.readouterr().err


# --- W-05: exports a person can open, against brute-force joins over the Parquet files ----------------------

W05_HIT_COLUMNS = ["accession", "header", "contig_len", "start", "end", "strand", "length", "seq_bitscore",
                   "seq_evalue_zref", "pfam_architecture", "protein_uid", "nr90", "nr50", "nr50_rep", "genome_uid",
                   "tier", "contig_uid", "feature_uid", "ordinal", "tape_version", "dist_to_contig_end_bp", "n_records",
                   "stage", "kind"]  # card W-23: the gene's stage and its genome's kind, appended
W05_LOCUS_COLUMNS = ["offset", "offset_rel", "ordinal", "feature_uid", "start", "end", "strand", "protein_uid",
                     "length_aa", "pfam_architecture", "nr90", "nr50", "hit_seq_bitscore", "delta_bp", "crosses_origin"]
W05_BITS, W05_EVALUE = {"seq_bitscore", "hit_seq_bitscore"}, {"seq_evalue_zref"}


@pytest.fixture(scope="module")
def w05(world, tmp_path_factory):
    """A copy of the world with a finished jackhmmer run and a contig_dim of REL written by hand: the two linear
    contigs in one made-up genome and tiers 1, 2, 3 in contig order, so genome and tier are constructed truth."""
    hit, _ = _seed(world[0])
    row = _run(world, _recipe(hit, hit, tool="jackhmmer", mx=2))
    assert row["status"] == "ok", row["error"]
    root = tmp_path_factory.mktemp("w05") / "a"
    shutil.copytree(world[0], root)
    # the run and the pack are beside the world's root, not under it (card W-03): the copy gets its own copies
    shutil.copytree(_rd(world, row["run_id"]), _rdir(root, row["run_id"]))
    shutil.copytree(_packs(world), _packs_of(root))
    tps = _tape_rows(root, tapes.tape_version([V]))
    lin = sorted(cu for cu, t in tps.items() if not t["is_circular"])
    one = ids.genome_uid(lin)
    rows = [{"contig_uid": cu, "genome_uid": one if cu in lin else ids.genome_uid([cu]), "contig_len": t["contig_len"],
             "is_circular": t["is_circular"], "tier": k + 1, "assembly_quality": None, "tax_order": None,
             "is_species_rep": False} for k, (cu, t) in enumerate(sorted(tps.items()))]
    d = root / f"dim/contig_dim/release={REL}"
    d.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=table("contig_dim").schema), d / "part-0.parquet")
    return root, row["run_id"].hex(), {r["contig_uid"]: r for r in rows}


def _w05_cell(col: str, v) -> str:
    if v is None:
        return ""
    if col in W05_BITS:
        return f"{v:.1f}"
    if col in W05_EVALUE:
        return f"{v:.1e}"
    if col == "strand":  # card W-20: +1 or -1 everywhere (migrated from the bare `str(v)`)
        return f"{v:+d}"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, bytes):
        return v.hex()
    return str(v).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def _pq(root: Path, base: str) -> list[dict]:
    return [r for f in sorted((root / base).rglob("*.parquet")) for r in pq.read_table(f, partitioning=None).to_pylist()]


def _w05_final_hits(root: Path, rid: str) -> list[dict]:
    """The run's kept hits of its final round, read from its files."""
    run = pq.read_table(_rdir(root, rid) / "run.parquet").to_pylist()[0]
    return [h for h in pq.read_table(_rdir(root, rid) / "hits.parquet").to_pylist()
            if h["round"] == run["final_round"] and h["curation_label"] == "keep"]


def _w05_store(root: Path) -> dict:
    """Everything the exports join, as dicts built by loops over the Parquet shards."""
    recs: dict[bytes, list[dict]] = {}
    for r in _pq(root, "l0seq/contig_occurrence"):
        recs.setdefault(r["contig_uid"], []).append(r)
    label = {}
    for cu, rs in recs.items():
        first = min(rs, key=lambda r: (r["l0_object_sha256"], r["record_ordinal"]))
        label[cu] = (first["accession"], first["original_header"], len(rs))
    assert sorted(p.name for p in (root / "l1/domains").glob("library=*")) == ["library=test-1"]
    doms: dict[bytes, set] = {}
    for r in _pq(root, "l1/domains"):
        doms.setdefault(r["protein_uid"], set()).add((r["env_start"], r["model_acc"], r["model_name"]))
    arch = {}
    for pu, ds in doms.items():
        names: list[str] = []
        for _, _, name in sorted(ds):
            if not names or names[-1] != name:
                names.append(name)
        arch[pu] = "+".join(names)
    return {"label": label, "arch": arch,
            "clen": {r["contig_uid"]: r["contig_len"] for r in _pq(root, "l0seq/contig_index") if r["is_default_variant"]},
            "prot": {r["protein_uid"]: r for r in _pq(root, f"l0seq/protein_seq/gene_call_version={V}")},
            "clus": {(r["member_uid"], r["level"]): r for r in _pq(root, f"l2/cluster_member/release={REL}")},
            "tapes": _tape_rows(root, tapes.tape_version([V]))}


def _w05_scores(hits: list[dict]) -> dict[bytes, tuple[float, float]]:
    score: dict[bytes, tuple[float, float]] = {}
    for h in hits:
        b, e = score.get(h["target_uid"], (h["seq_bitscore"], h["seq_evalue_zref"]))
        score[h["target_uid"]] = (max(b, h["seq_bitscore"]), min(e, h["seq_evalue_zref"]))
    return score


def _w05_oracle_genes(root: Path, hits: list[dict], dims: dict) -> list[dict]:
    """One dict per gene of a hit protein, in the export's order (accession, contig, coordinates, protein)."""
    s, score, tv = _w05_store(root), _w05_scores(hits), tapes.tape_version([V])
    # card W-23: each gene's stage is its occurrence file's, its kind the genome_input kind of the contig in that stage
    kinds = {(r["contig_uid"], f.parent.name.removeprefix("stage=")): r["kind"]
             for f in sorted((root / "l0seq/genome_input").rglob("*.parquet"))
             for r in pq.read_table(f, partitioning=None).to_pylist()}
    occs = [{**o, "_stage": f.parent.name.removeprefix("stage=")}
            for f in sorted((root / f"l1/protein_occurrence/gene_call_version={V}").rglob("*.parquet"))
            for o in pq.read_table(f, partitioning=None).to_pylist()]
    out = []
    for o in occs:
        pu, cu = o["protein_uid"], o["contig_uid"]
        if pu not in score:
            continue
        acc, hdr, n = s["label"].get(cu, (None, None, None))
        c90, c50, d = s["clus"].get((pu, "nr90")), s["clus"].get((pu, "nr50")), dims.get(cu)
        row = {"accession": acc, "header": hdr, "contig_len": s["clen"][cu], "start": o["start"], "end": o["end"],
               "strand": o["strand"], "length": s["prot"][pu]["length"], "seq_bitscore": score[pu][0],
               "seq_evalue_zref": score[pu][1], "pfam_architecture": s["arch"].get(pu) or "-", "protein_uid": pu,
               "nr90": c90 and c90["cluster_uid"], "nr50": c50 and c50["cluster_uid"], "nr50_rep": c50 and c50["rep_uid"],
               "genome_uid": d and d["genome_uid"], "tier": d and d["tier"], "contig_uid": cu,
               "feature_uid": o["feature_uid"], "ordinal": s["tapes"][cu]["feature_uids"].index(o["feature_uid"]),
               "tape_version": tv, "dist_to_contig_end_bp": o["dist_to_contig_end_bp"], "n_records": n,
               "stage": o["_stage"], "kind": kinds.get((cu, o["_stage"]))}
        out.append(((acc is None, acc or "", cu, o["start"], o["end"], o["strand"], pu), row))
    return [r for _, r in sorted(out, key=lambda x: x[0])]


def _w05_tsv(cols: list[str], rows: list[dict]) -> list[str]:
    return ["\t".join(cols)] + ["\t".join(_w05_cell(c, r[c]) for c in cols) for r in rows]


def _w05_hits(root: Path) -> list[dict]:
    """Hits made by hand, as search_hit rows: the W-06 anchors (genes next to the circular contig's origin, at each
    linear contig end, beside the domain-carrying lambda genes, and the library-hit seed), each with its own
    score and E-value, two domain rows for the first (the export keeps one row per gene)."""
    rows = [{"target_kind": "protein", "target_uid": u, "seq_bitscore": 20.0 + 3.37 * k,
             "seq_evalue_zref": 1.234 * 10.0 ** (-k)} for k, u in enumerate(sorted(_w06_anchors(root)))]
    return rows + [dict(rows[0])]


@pytest.mark.covers("W-05")
def test_w05_search_hits_tsv_equals_a_brute_force_join_over_the_parquet_files(w05, world, capsys, tmp_path):
    root, rid, dims = w05
    args = ["search", "hits", rid, "--root", str(root), *_dirs(root), "--format", "tsv"]
    assert cli.main(args) == 0
    out = capsys.readouterr().out
    genes = _w05_oracle_genes(root, _w05_final_hits(root, rid), dims)
    assert genes and out == "\n".join(_w05_tsv(W05_HIT_COLUMNS, genes)) + "\n"
    # deterministic bytes, and --out writes the same bytes
    assert cli.main(args) == 0
    assert capsys.readouterr().out == out
    dst = tmp_path / "hits.tsv"
    assert cli.main(args + ["--out", str(dst)]) == 0
    assert dst.read_bytes() == out.encode()
    # a broader hit set on every contig, through the same export
    made = _w05_hits(root)
    genes = _w05_oracle_genes(root, made, dims)
    got = list(export.hits_tsv(root, made, release=REL, tape_version=tapes.tape_version([V]), gene_call_version=V))
    assert got == _w05_tsv(W05_HIT_COLUMNS, genes)
    # the world is what the card asks for: genes on every contig, with and without a domain architecture, genome
    # and tier from contig_dim, cluster ids that are cluster uids (16 bytes as hex), never token ids
    assert len({g["contig_uid"] for g in genes}) == 3 and len(genes) >= len(made) - 1
    assert {g["pfam_architecture"] == "-" for g in genes} == {True, False}
    assert all(g["tier"] is not None and g["genome_uid"] is not None for g in genes)
    assert all(len(_w05_cell("nr50", g["nr50"])) == 32 for g in genes)
    # without contig_dim (the world itself) genome and tier are empty, not guessed
    bare = list(export.hits_tsv(world[0], made, release=REL, tape_version=tapes.tape_version([V]),
                                gene_call_version=V))
    # (card W-20: and no `tier` column at all, since no gene has one; migrated from an always-empty column)
    assert bare == _w05_tsv([c for c in W05_HIT_COLUMNS if c != "tier"], _w05_oracle_genes(world[0], made, {}))


@pytest.mark.covers("W-05")
def test_w05_search_hits_faa_sequences_equal_protein_seq(w05, capsys):
    root, rid, dims = w05
    prot = {r["protein_uid"]: r for r in _pq(root, f"l0seq/protein_seq/gene_call_version={V}")}

    def oracle(hits):
        genes, want = _w05_oracle_genes(root, hits, dims), []
        for pu in sorted({g["protein_uid"] for g in genes}):
            mine = [g for g in genes if g["protein_uid"] == pu]
            g = mine[0]
            place = f"{g['accession'] or g['contig_uid'].hex()}:{g['start']}-{g['end']} strand={g['strand']:+d}"  # W-20
            desc = (f"length={prot[pu]['length']} seq_bitscore={_w05_cell('seq_bitscore', g['seq_bitscore'])} "
                    f"seq_evalue_zref={_w05_cell('seq_evalue_zref', g['seq_evalue_zref'])} n_genes={len(mine)} "
                    f"example={place}" + (f" {_w05_cell('header', g['header'])}" if g["header"] else ""))
            want += [f">{pu.hex()} {desc}", prot[pu]["seq"]]
        return want

    args = ["search", "hits", rid, "--root", str(root), *_dirs(root), "--format", "faa"]
    assert cli.main(args) == 0
    out = capsys.readouterr().out
    assert out == "\n".join(oracle(_w05_final_hits(root, rid))) + "\n"
    assert cli.main(args) == 0
    assert capsys.readouterr().out == out
    made = _w05_hits(root)
    got = list(export.hits_faa(root, made, release=REL, tape_version=tapes.tape_version([V]), gene_call_version=V))
    assert got == oracle(made)
    names = [x[1:].split(" ")[0] for x in got if x.startswith(">")]
    assert names == sorted(names) and len(names) == len({h["target_uid"] for h in made}) > 3


@pytest.mark.covers("W-05", "S-13")
def test_w05_hit_json_is_version_2_with_genome_and_tier_joined_from_contig_dim(w05, world, capsys):
    root, rid, dims = w05
    locked = json.loads(export.lock_path(REPO).read_text())
    assert locked["search_hit"]["version"] == 2 and export.SCHEMAS["search_hit"]["version"] == 2
    assert export.lock_problems(locked) == [] and locked == json.loads(json.dumps(export.SCHEMAS, sort_keys=True))
    for r, want_dim in ((root, True), (world[0], False)):  # the world has no contig_dim: null, not absent
        assert cli.main(["search", "hits", rid, "--root", str(r), *_dirs(r), "--json", "--members"]) == 0
        recs = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
        assert recs and all(h["schema"] == "search_hit/2" and list(h) == locked["search_hit"]["keys"] for h in recs)
        occ = [o for h in recs for o in h["occurrences"]]
        assert occ and all(list(o) == ["contig_uid", "ordinal", "feature_uid", "start", "end", "strand", "genome_uid",
                                       "tier"] for o in occ)
        for o in occ:
            d = dims[bytes.fromhex(o["contig_uid"])]
            assert (o["genome_uid"], o["tier"]) == ((d["genome_uid"].hex(), d["tier"]) if want_dim else (None, None))


def _w05_oracle_loci(root: Path, hits: list[dict], window: str) -> list[tuple[str, list[str]]]:
    """(file name, lines) per anchor locus: windows from the S-09 coordinate oracle (tests/r1 `_oracle_nbr`) over
    the tapes, annotations from `_w05_store`; loci by accession, contig and anchor ordinal; rows by offset."""
    r1 = _r1()
    policy = nbr.Policy.parse(window=window)
    s, score, tv = _w05_store(root), _w05_scores(hits), tapes.tape_version([V])
    loci = []
    for cu, t in s["tapes"].items():
        slots = list(zip(t["start_bps"], t["end_bps"], t["strands"]))
        for i, apu in enumerate(t["protein_uids"]):
            if apu not in score:
                continue
            got = r1._oracle_nbr(slots, t["contig_len"], t["is_circular"], i, policy)
            if got is None:
                continue
            acc, hdr, _ = s["label"].get(cu, (None, None, None))
            rows = []
            for j, d, rel, gap, fl in sorted(got, key=lambda x: x[1]):
                pu = t["protein_uids"][j]
                c90, c50 = s["clus"].get((pu, "nr90")), s["clus"].get((pu, "nr50"))
                rows.append({"offset": d, "offset_rel": rel, "ordinal": j, "feature_uid": t["feature_uids"][j],
                             "start": slots[j][0], "end": slots[j][1], "strand": slots[j][2], "protein_uid": pu,
                             "length_aa": s["prot"][pu]["length"] if pu else None,
                             "pfam_architecture": (s["arch"].get(pu) or "-") if pu else None,
                             "nr90": c90 and c90["cluster_uid"], "nr50": c50 and c50["cluster_uid"],
                             "hit_seq_bitscore": score[pu][0] if pu in score else None,
                             "delta_bp": gap if d != 0 else None, "crosses_origin": bool(fl & 2)})
            cens = any(fl & 1 for _, d, _, _, fl in got if d == 0)
            facts = [("contig_uid", cu.hex()), ("anchor_feature_uid", t["feature_uids"][i].hex()), ("anchor_ordinal", i),
                     ("anchor_protein_uid", apu.hex()), ("tape_version", tv), ("policy_hash", policy.hash()),
                     ("censored", "true" if cens else "false")]
            lines = [f"# {_w05_cell('header', hdr) if hdr else cu.hex()}",
                     "# " + "\t".join(f"{k}={v}" for k, v in facts)] + _w05_tsv(W05_LOCUS_COLUMNS, rows)
            name = f"locus_{re.sub(r'[^A-Za-z0-9._-]', '_', acc or cu.hex())}_{t['feature_uids'][i].hex()}.tsv"
            loci.append(((acc is None, acc or "", cu, i), name, lines))
    return [(name, lines) for _, name, lines in sorted(loci, key=lambda x: x[0])]


@pytest.mark.covers("W-05", "S-09")
def test_w05_nbr_tsv_writes_one_table_per_locus_equal_to_the_coordinate_oracle(w05, capsys, tmp_path):
    root, rid, _ = w05
    args = ["nbr", rid, "--root", str(root), *_dirs(root), "--window", "3g", "--format", "tsv"]
    assert cli.main(args) == 0
    out = capsys.readouterr().out
    want = _w05_oracle_loci(root, _w05_final_hits(root, rid), "3g")
    assert want and out == "\n\n".join("\n".join(lines) for _, lines in want) + "\n"  # tables, a blank line between
    blocks = out.rstrip("\n").split("\n\n")
    assert cli.main(args) == 0
    assert capsys.readouterr().out == out
    # --out DIR: one file per locus with the same table
    d = tmp_path / "loci"
    assert cli.main(args + ["--out", str(d)]) == 0
    assert sorted(p.name for p in d.iterdir()) == sorted(name for name, _ in want)
    for name, lines in want:
        assert (d / name).read_text() == "\n".join(lines) + "\n"
    # anchors given as protein uids need the release named, as nbr summarize does
    uids = ",".join(sorted({h["target_uid"].hex() for h in _w05_final_hits(root, rid)}))
    assert cli.main(["nbr", uids, "--root", str(root), *_dirs(root), "--window", "3g", "--format", "tsv"]) == 2
    assert "--release" in capsys.readouterr().err
    assert cli.main(["nbr", uids, "--root", str(root), *_dirs(root), "--window", "3g", "--format", "tsv",
                     "--release", REL]) == 0
    no_scores = capsys.readouterr().out.rstrip("\n").split("\n\n")
    assert len(no_scores) == len(blocks)
    # a broader anchor set through the same export: windows across the circular origin and cut at linear ends
    made = _w05_hits(root)
    for window in ("3g", "5g"):
        rows, _ = nbr.fetch(root, {h["target_uid"] for h in made}, nbr.Policy.parse(window=window),
                            files=catalog.ReleaseFiles.of(root, REL))
        got = list(export.locus_tables(root, rows, release=REL, gene_call_version=V,
                                       hit_bits={h["target_uid"]: h["seq_bitscore"] for h in made}))
        want = _w05_oracle_loci(root, made, window)
        assert got == want
        assert len({lines[1].split("\t")[0] for _, lines in want}) == 3  # loci on every contig
        assert any("censored=true" in lines[1] for _, lines in want)
        assert any(line.endswith("\ttrue") for _, lines in want for line in lines[3:])  # crosses the origin
        assert any(line.split("\t")[12] for _, lines in want for line in lines[3:] if not line.startswith("0\t"))


@pytest.mark.covers("W-05")
def test_w05_export_proteins_writes_every_distinct_protein_named_by_uid(w05, capsys, tmp_path):
    root = tmp_path / "a"
    shutil.copytree(w05[0], root)
    prot = {r["protein_uid"]: r["seq"] for r in _pq(root, f"l0seq/protein_seq/gene_call_version={V}")}
    want = "".join(f">{u.hex()}\n{prot[u]}\n" for u in sorted(prot))
    args = ["export", "proteins", "--release", REL, "--root", str(root)]
    # card W-11: without the release manifest the export refuses (see test_w11_export_proteins_without_...)
    assert cli.main(args) == 2
    capsys.readouterr()
    # with the release manifest the files come from it; every distinct protein, and --out writes the same bytes
    catalog.write_manifest(catalog.build_manifest(root, REL, ["S0"]), root / f"manifests/release={REL}.json")
    dst = tmp_path / "p.faa"
    assert cli.main(args + ["--out", str(dst)]) == 0
    assert dst.read_text() == want
    assert cli.main(args) == 0
    assert capsys.readouterr().out == want


@pytest.mark.covers("W-11")
def test_w11_export_proteins_without_the_release_manifest_refuses(w05, capsys, tmp_path):
    root = tmp_path / "a"
    shutil.copytree(w05[0], root)
    assert not (root / f"manifests/release={REL}.json").exists()
    dst = tmp_path / "p.faa"
    capsys.readouterr()
    assert cli.main(["export", "proteins", "--release", REL, "--root", str(root), "--out", str(dst)]) == 2
    cap = capsys.readouterr()
    assert cap.out == "" and not dst.exists()
    assert f"manifests/release={REL}.json" in cap.err and "bio-box release manifest" in cap.err
    # the library call refuses too, before any iteration, instead of reading every file it finds
    with pytest.raises(ValueError, match="manifest"):
        export.protein_fasta(root, release=REL, gene_call_version=V)


@pytest.mark.covers("W-11")
def test_w11_no_stored_protein_contains_a_star(w05):
    # rL.2 (2,053,315 proteins, measured 2026-11-10 with scripts/stress/protein_stars.py) holds none, so the
    # FASTA writers keep the sequence as stored; this pins the fixture to the same property.
    prot = _pq(w05[0], f"l0seq/protein_seq/gene_call_version={V}")
    assert prot and all("*" not in r["seq"] for r in prot)


# --- W-12: hit and locus exports read the release's files only -----------------------------------------------

import test_release_readers as rr  # noqa: E402
from test_release_readers import store  # noqa: E402,F401  (the two-stage fixture: stages S0 and S1, release R1 over S0)


@pytest.fixture(scope="module")
def w12(store):
    """The two-stage store with a finished search of R1 seeded with protein P, which also sits on an S1 contig."""
    row = runner.run(store["root"], rr._recipe(rr.R1, store["p"]), rr.R1, work_base=store["d"] / "work12",
                     limits=rr.LIMITS, **rr.at(store))
    assert row["status"] == "ok", row.get("error")
    return store, row["run_id"].hex()


def _w12_commands(store, root: Path, rid: str) -> dict[str, list[str]]:
    """The commands against the data root `root` (a copy of the store's); the run and the pack stay the store's,
    in its runs and pack directories (card W-03)."""
    where = ["--runs", str(store["runs"]), "--pack-dir", str(store["pack_dir"])]
    return {"hits-tsv": ["search", "hits", rid, "--root", str(root), *where, "--format", "tsv"],
            "hits-faa": ["search", "hits", rid, "--root", str(root), *where, "--format", "faa"],
            "nbr-tsv": ["nbr", rid, "--root", str(root), *where, "--window", "3g", "--format", "tsv"]}


def _w12_outputs(store, root: Path, rid: str, capsys) -> dict[str, str]:
    capsys.readouterr()
    out = {}
    for name, args in _w12_commands(store, root, rid).items():
        assert cli.main(args) == 0, name
        out[name] = capsys.readouterr().out
    return out


def _w12_copy(store, tmp_path: Path, manifest: bool) -> Path:
    root = tmp_path / "a"
    shutil.copytree(store["root"], root)
    path = catalog.manifest_path(root, rr.R1)
    path.unlink(missing_ok=True)
    if manifest:
        catalog.write_manifest(catalog.build_manifest(root, rr.R1, ["S0"]), path)
    return root


def _w12_check_s0_only(store, out: dict[str, str]) -> None:
    s1 = [c.hex() for c in store["s1_contigs"]]
    assert all(out.values())
    for name, text in out.items():
        assert not any(c in text for c in s1), name  # no S1 contig, nowhere
    rows = [line.split("\t") for line in out["hits-tsv"].rstrip("\n").split("\n")]
    col = rows[0].index("contig_uid")
    genes = [r for r in rows[1:] if r[rows[0].index("protein_uid")] == store["p"].hex()]
    assert [r[col] for r in genes] == [store["lam"].hex()]  # P has one gene in R1: S1's copy is not R1's
    faa = [line for line in out["hits-faa"].split("\n") if line.startswith(">")]
    mine = [line for line in faa if line.startswith(">" + store["p"].hex())]
    assert len(mine) == 1 and " n_genes=1 " in mine[0]
    assert {line.split()[0][1:] for line in faa}.isdisjoint(u.hex() for u in store["s1_only"])


@pytest.mark.parametrize("manifest", [False, True], ids=["release_def", "manifest"])
@pytest.mark.covers("W-12", "S-09")
def test_w12_hit_and_locus_exports_of_r1_return_nothing_from_s1(w12, tmp_path, capsys, manifest):
    store, rid = w12
    root = _w12_copy(store, tmp_path, manifest)
    assert catalog.manifest_path(root, rr.R1).exists() == manifest
    assert any(p.is_dir() for p in root.rglob("stage=S1"))  # the stage the release does not list is on disk
    _w12_check_s0_only(store, _w12_outputs(store, root, rid, capsys))


@pytest.mark.covers("W-12")
def test_w12_release_def_and_manifest_give_the_same_bytes(w12, tmp_path, capsys):
    store, rid = w12
    bare = _w12_outputs(store, _w12_copy(store, tmp_path / "x", False), rid, capsys)
    with_man = _w12_outputs(store, _w12_copy(store, tmp_path / "y", True), rid, capsys)
    assert bare == with_man


@pytest.mark.covers("W-12")
def test_w12_with_a_manifest_an_extra_stage_on_disk_changes_no_output_byte(w12, tmp_path, capsys):
    store, rid = w12
    full = _w12_copy(store, tmp_path / "full", True)
    lean = tmp_path / "lean"
    shutil.copytree(full, lean)
    gone = sorted(lean.rglob("stage=S1"))
    assert gone
    for d in gone:
        shutil.rmtree(d)
    assert not list(lean.rglob("stage=S1"))
    assert _w12_outputs(store, lean, rid, capsys) == _w12_outputs(store, full, rid, capsys)


@pytest.mark.covers("W-12")
def test_w12_an_undefined_release_exits_2_naming_the_release_and_release_define(w12, tmp_path, capsys):
    store, rid = w12
    root = _w12_copy(store, tmp_path, False)
    catalog.release_def_path(root, rr.R1).unlink()
    capsys.readouterr()
    for name, args in _w12_commands(store, root, rid).items():
        assert cli.main(args) == 2, name
        cap = capsys.readouterr()
        assert cap.out == "", name
        assert rr.R1 in cap.err and "bio-box release define" in cap.err, (name, cap.err)
    # the library calls refuse too, before reading any file
    with pytest.raises(ValueError, match="bio-box release define"):
        list(export.hits_tsv(root, [], release=rr.R1, tape_version=tapes.tape_version([V]), gene_call_version=V))


@pytest.mark.covers("W-12")
def test_w12_export_py_lists_no_files_by_globbing():
    src = (Path(export.__file__)).read_text()
    assert "rglob" not in src and "glob(" not in src


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
        coords = [(100, 460, 1), (600, 1080, 1), (1200, 1560, 1)]
        prots = [ids.protein_uid(x) for x in (left, target, right)]
        for i, ((st_, en, sd), pu) in enumerate(zip(coords, prots)):
            occ.append({"protein_uid": pu, "contig_uid": cu, "start": st_, "end": en, "strand": sd,
                        "dist_to_contig_end_bp": 0 if partial and i == 1 else min(st_, 1700 - en),
                        "feature_uid": ids.feature_uid(cu, st_, en, sd, ids.KIND["cds"])})
        tapes_.append({"contig_uid": cu, "contig_len": 1700, "is_circular": False, "n_slots": 3,
                       **build_tape(coords, 1700, False), "protein_uids": prots,
                       "feature_uids": [ids.feature_uid(cu, a, b, sd, ids.KIND["cds"]) for a, b, sd in coords]})
    occ.sort(key=lambda r: (r["protein_uid"], r["contig_uid"], r["start"], r["end"], r["strand"]))
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
    catalog.define_release(root, "rSyn", ["S0"])  # its pack and runs read the release's stages (card B-28)
    return world


def test_s07_drift_alarm_reverts_a_planted_corruption(tmp_path):
    root = tmp_path / "a"
    w = trap_world(root, tmp_path / "trap.hmm")
    pack.build(root, "rSyn", models={"trap-1": tmp_path / "trap.hmm"}, packs=_packs_of(root))
    rows = {}
    for name, stop in (("guarded", "[no_new_reps, drift_alarm]"), ("unguarded", "[no_new_reps]")):
        text = (f"recipe: trap\nwritten_against: rSyn\nseed: [{{sequence: {w['A']}}}]\ntarget: unique\n"
                f"tool: jackhmmer\niterations: {{max: 5, inclusion_bits: 30, stop: {stop}}}\n"
                f"controls: {{positive: [{{sequence: {w['a'][0]}}}]}}\n")
        rows[name] = runner.run(root, text, "rSyn", work_base=tmp_path / "work", limits=LIMITS, **_at(root))
    g, u = rows["guarded"], rows["unguarded"]
    leaked = {n: len({h["target_uid"] for h in runner.read_hits(_rdir(root, r["run_id"]))} & w["b"])
              for n, r in rows.items()}
    # the trap is real: without the alarm the profile drifts into family B
    assert u["status"] == "ok" and leaked["unguarded"] >= 25
    # with it, the round that drifts is logged, dropped, and the run reverts to the round before
    assert g["status"] == "ok", g["error"]
    assert g["stop_reason"].startswith("drift_alarm: TRAP0001") and g["final_round"] == g["rounds"] - 1
    assert leaked["guarded"] <= 3
    dropped = [h for h in runner.read_hits(_rdir(root, g["run_id"]), round_=g["rounds"], final_only=False)]
    assert dropped and all(h["curation_label"] == "drop" and h["curation_reason"] == "drift_alarm" for h in dropped)


def _warnings(root, run):
    import json

    row = pq.read_table(next(_runs_of(root).glob(f"run={run['run_id'].hex()}/run.parquet"))).to_pylist()[0]
    return json.loads(row["warnings_json"])


def test_w08_drift_alarm_warns_when_the_release_has_no_domain_table(tmp_path):
    import shutil

    texts = {}
    for name in ("with", "without"):
        root = tmp_path / name
        w = trap_world(root, tmp_path / f"{name}.hmm")
        pack.build(root, "rSyn", models={"trap-1": tmp_path / f"{name}.hmm"}, packs=_packs_of(root))
        if name == "without":
            shutil.rmtree(root / "l1/domains")
        text = (f"recipe: trap\nwritten_against: rSyn\nseed: [{{sequence: {w['A']}}}]\ntarget: unique\n"
                f"tool: jackhmmer\niterations: {{max: 3, inclusion_bits: 30, stop: [no_new_reps, drift_alarm]}}\n"
                f"controls: {{positive: [{{sequence: {w['a'][0]}}}]}}\n")
        texts[name] = runner.run(root, text, "rSyn", work_base=tmp_path / f"work-{name}", limits=LIMITS,
                                 **_at(root))
        texts[name] = (root, texts[name])
    (r1, with_), (r2, without) = texts["with"], texts["without"]
    assert not any("domain annotations" in x for x in _warnings(r1, with_))
    assert any("drift_alarm could not be evaluated" in x and "no domain annotations" in x for x in _warnings(r2, without))
    assert without["status"] == with_["status"] == "ok"


def test_search_tables_conform_to_the_registry(world):
    root = world[0]
    for p in sorted(_runs_of(root).glob("run=*/*.parquet")):
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
    pack.build(root, "rSyn", models=w["models"], packs=_packs_of(root))
    return root, d, w


def test_s07_context_lost_reverts_when_new_hits_lose_round1_neighbors(placed):
    root, d, w = placed
    row = runner.run(root, _trap_recipe(w, "[no_new_reps, context_lost]"), "rSyn", work_base=d / "work", limits=LIMITS,
                     **_at(root))
    assert row["status"] == "ok", row["error"]
    assert row["stop_reason"].startswith("context_lost") and row["final_round"] == row["rounds"] - 1
    log = json.loads(row["rounds_json"])
    assert log[-1]["context_retention"] < 0.5
    final = {h["target_uid"] for h in runner.read_hits(_rdir(root, row["run_id"]))}
    assert len(final & w["b"]) <= 3
    # the same world without the rule drifts into B: the check, not the world, stops it
    free = runner.run(root, _trap_recipe(w), "rSyn", work_base=d / "work", limits=LIMITS, **_at(root))
    assert len({h["target_uid"] for h in runner.read_hits(_rdir(root, free["run_id"]))} & w["b"]) >= 25


@pytest.mark.parametrize("policy", ["exclude", "flag", "include"])
def test_curation_labels_for_coverage_and_partial_genes(placed, policy):
    root, d, w = placed
    text = _trap_recipe(w, tool="phmmer", mx=1,
                        extra=f"partial: {policy}\ncuration: {{min_query_coverage: 0.9}}\n")
    row = runner.run(root, text, "rSyn", work_base=d / "work", limits=LIMITS, **_at(root))
    assert row["status"] == "ok", row["error"]
    rows = runner.read_hits(_rdir(root, row["run_id"]), round_=1, final_only=False)
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
    final = {h["target_uid"] for h in runner.read_hits(_rdir(root, row["run_id"]))}
    assert final == kept


HMMSEARCH = shutil.which("hmmsearch")


@pytest.mark.skipif(HMMSEARCH is None, reason="hmmsearch binary not on PATH")
def test_rounds_equal_the_hmmsearch_binary(placed, tmp_path):
    """O-SEARCH oracle: each round's stored profile, searched by HMMER's own hmmsearch at the
    same Z and threshold, gives the same included targets, sequence scores and kept domains."""
    root, d, w = placed
    row = runner.run(root, _trap_recipe(w), "rSyn", work_base=d / "work", limits=LIMITS, **_at(root))
    assert row["status"] == "ok" and row["rounds"] >= 2
    faa = tmp_path / "unique.faa"
    with open(faa, "w") as fh:
        for sh in sorted(pack.find(_packs_of(root), "rSyn").glob("unique/part-*.faa.zst")):
            for name, seq in fastz.read_records(sh):
                fh.write(f">{name}\n{seq}\n")
    profiles = pq.read_table(_rdir(root, row["run_id"]) / "profiles.parquet").to_pylist()
    hits = pq.read_table(_rdir(root, row["run_id"]) / "hits.parquet").to_pylist()
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
                  for h in runner.read_hits(_rdir(root, row["run_id"]), final_only=False))


def _pinned(root: Path, text: str, rel: str = REL) -> tuple[str, str, dict]:
    r, _ = rcp.load(text)
    res = runner.resolve(root, r, rel, _packs_of(root))
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
    before = sorted(_runs_of(root).glob("run=*"))
    with pytest.raises(RecipeError, match=match):
        _run(world, _file_recipe(kind, f, pos, seed_extra=extra))
    assert sorted(_runs_of(root).glob("run=*")) == before  # a refusal records nothing


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
    hits = runner.read_hits(_rdir(root, row["run_id"]))
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
                     work_base=d / "work", limits=LIMITS, **_at(root))
    assert row["status"] == "ok", row["error"]
    final = {h["target_uid"] for h in runner.read_hits(_rdir(root, row["run_id"]))}
    assert {u for u, _, _ in seeds} <= final
    assert len(final & w["b"]) >= 25  # three seeds find the family, not only themselves
    prof = [p for p in pq.read_table(_rdir(root, row["run_id"]) / "profiles.parquet").to_pylist()
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
                     work_base=d / "work", limits=LIMITS, **_at(root))
    assert row["status"] == "ok", row["error"]
    profiles = pq.read_table(_rdir(root, row["run_id"]) / "profiles.parquet").to_pylist()
    assert row["rounds"] >= 2 and len(profiles) == row["rounds"]
    assert {u for u, _, _ in seeds} <= {h["target_uid"] for h in runner.read_hits(_rdir(root, row["run_id"]))}


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
                     work_base=d / "work", limits=LIMITS, **_at(root))
    assert row["status"] == "ok", row["error"]
    got = {h["target_uid"].hex(): h["seq_bitscore"] for h in runner.read_hits(_rdir(root, row["run_id"]))}
    unique = tmp_path / "unique.faa"
    with open(unique, "w") as fh:
        for sh in sorted(pack.find(_packs_of(root), "rSyn").glob("unique/part-*.faa.zst")):
            for name, seq in fastz.read_records(sh):
                fh.write(f">{name}\n{seq}\n")
    first = tmp_path / "first.sto"
    first.write_text(f"# STOCKHOLM 1.0\n\n{seeds[0][1]} {seeds[0][2]}\n//\n")
    run = lambda *a: subprocess.run([str(x) for x in a], check=True, capture_output=True)  # noqa: E731
    run(HMMER["hmmbuild"], "--amino", "--singlemx", "-n", "first", "-o", "/dev/null", tmp_path / "first.hmm", first)
    run(HMMER["hmmalign"], "--trim", "--amino", "-o", tmp_path / "seeds.sto", tmp_path / "first.hmm", faa)
    run(HMMER["hmmbuild"], "--amino", "-n", "w02", "-o", "/dev/null", tmp_path / "w02.hmm", tmp_path / "seeds.sto")
    stored = tmp_path / "stored.hmm"
    prof = pq.read_table(_rdir(root, row["run_id"]) / "profiles.parquet").to_pylist()
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
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout", **_at(root))
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    runner.cancel(_rdir(root, rid))
    assert _wait_for(lambda: runner.status(_rdir(root, rid))["status"] != "running")
    st = runner.status(_rdir(root, rid))
    assert st["status"] == "stopped", st
    assert sorted(p.name for p in _rdir(root, rid).iterdir()) == ["run.parquet"]


def test_runner_death_is_reported_and_the_worker_does_not_outlive_it(world, tmp_path):
    import os
    import signal

    root, d, _ = world
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit))
    work = tmp_path / "work"
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout", **_at(root))
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    worker_pid = int((work / rid / "worker.pid").read_text())
    runner_pid = json.loads(runner.pending_path(_rdir(root, rid)).read_text())["pid"]
    os.kill(runner_pid, signal.SIGKILL)
    assert _wait_for(lambda: not _alive(runner_pid), 10)
    st = runner.status(_rdir(root, rid))
    assert st["status"] == "failed" and "without closing" in st["error"]
    assert _wait_for(lambda: not _alive(worker_pid), 10), "worker outlived its runner"
    assert not _rdir(root, rid).exists()


@pytest.mark.covers("A-01")
def test_status_does_not_trust_a_pid_a_live_unrelated_process_may_own(tmp_path):
    """Pid reuse: the recorded pid answers kill(0), but nobody holds the runner lock: the run is dead."""
    import os

    d = runner.run_dir(tmp_path / "runs", "ab")
    pend = runner.pending_path(d)
    pend.parent.mkdir(parents=True)
    lock = tmp_path / "work" / "ab" / "runner.lock"
    lock.parent.mkdir(parents=True)
    lock.write_text("")
    pend.write_text(json.dumps({"run_id": "ab", "pid": os.getpid(), "lock": str(lock), "status": "running"}))
    st = runner.status(d)
    assert st["status"] == "failed" and "without closing" in st["error"], st
    import fcntl

    with open(lock, "a") as f:  # status must have released the lock it probed with
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)


@pytest.mark.covers("A-01")
def test_cancel_does_not_signal_a_pid_whose_runner_is_gone(tmp_path):
    """The runner died and its pid was reused: cancel must leave the new owner of the pid alone."""
    import sys

    d = runner.run_dir(tmp_path / "runs", "cd")
    pend = runner.pending_path(d)
    pend.parent.mkdir(parents=True)
    lock = tmp_path / "work" / "cd" / "runner.lock"
    lock.parent.mkdir(parents=True)
    lock.write_text("")
    bystander = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        pend.write_text(json.dumps({"run_id": "cd", "pid": bystander.pid, "lock": str(lock), "status": "running"}))
        runner.cancel(d)
        with pytest.raises(subprocess.TimeoutExpired):
            bystander.wait(timeout=1.0)
    finally:
        bystander.kill()
        bystander.wait()


# A runner whose main thread has exited while another thread of it still runs: the lock it inherited
# stays held (the thread group's file table is open), and /proc shows the leader as a zombie. The leader
# leaves through the raw exit syscall, which ends one thread (exit_group would end them all).
_ZOMBIE_LEADER = """
import ctypes, os, sys, threading
threading.Thread(target=os.read, args=(0, 1)).start()  # holds the process open until stdin closes
print("ready", flush=True)
ctypes.CDLL(None).syscall({nr}, 0)
"""
_SYS_EXIT = {"x86_64": 60, "aarch64": 93}


def _proc_state(pid: int) -> str | None:
    try:
        s = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return None
    return s[s.rindex(")") + 2:].split()[0]


@pytest.mark.covers("A-07", "A-01")
@pytest.mark.skipif(platform.system() != "Linux" or platform.machine() not in _SYS_EXIT,
                    reason="a zombie thread-group leader is a Linux state; the syscall number is per arch")
def test_a07_status_is_failed_once_the_runner_leader_is_a_zombie_though_its_lock_is_still_held(tmp_path):
    """CI flake (2): after SIGKILL a multi-threaded runner's leader can be a zombie (what `_alive` and
    `ps` report as gone) while its last thread still tears the process down and holds runner.lock
    (2 of 20 kills on the pod: Z with the lock held for 42-43 ms; measured 2026-11-10). Injected
    condition: a stand-in runner whose leader thread exits while a second thread keeps the lock."""
    import fcntl
    import os

    lock = tmp_path / "work" / "ef" / runner.LOCK_NAME
    lock.parent.mkdir(parents=True)
    fd = os.open(lock, os.O_RDWR | os.O_CREAT, 0o644)
    fcntl.flock(fd, fcntl.LOCK_EX)  # taken and handed down as `submit` does
    proc = subprocess.Popen([sys.executable, "-c", _ZOMBIE_LEADER.format(nr=_SYS_EXIT[platform.machine()])],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, pass_fds=(fd,), text=True)
    os.close(fd)
    try:
        assert proc.stdout.readline().strip() == "ready"
        assert _wait_for(lambda: _proc_state(proc.pid) == "Z", 10)
        assert not _alive(proc.pid)  # gone as the runner-death test and `ps` see it
        assert runner._runner_holds(str(lock))  # the condition: the lock is still held
        d = runner.run_dir(tmp_path / "runs", "ef")
        pend = runner.pending_path(d)
        pend.parent.mkdir(parents=True)
        pend.write_text(json.dumps({"run_id": "ef", "pid": proc.pid, "lock": str(lock), "status": "running"}))
        st = runner.status(d)
        assert st["status"] == "failed" and "without closing" in st["error"], st
    finally:
        proc.stdin.close()
        proc.wait(timeout=30)


@pytest.mark.covers("A-07")
@pytest.mark.skipif(platform.system() != "Linux", reason="reads /proc/<pid>/fd")
def test_a07_only_the_runner_holds_a_descriptor_of_its_lock(world, tmp_path):
    """The lock is the runner's alone: no fd of the worker (or of systemd-run or nice, which exec it)
    refers to runner.lock, so the worker cannot keep a dead runner looking alive."""
    import os

    root, d, _ = world
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit))
    work = tmp_path / "work"
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout", **_at(root))
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    info = json.loads(runner.pending_path(_rdir(root, rid)).read_text())
    worker_pid = int((work / rid / "worker.pid").read_text())

    def lock_fds(pid):
        out = []
        for n in os.listdir(f"/proc/{pid}/fd"):
            try:
                out.append(os.readlink(f"/proc/{pid}/fd/{n}"))
            except OSError:
                pass
        return [p for p in out if p == info["lock"]]

    try:
        assert lock_fds(info["pid"]) == [info["lock"]]
        assert lock_fds(worker_pid) == []
    finally:
        runner.cancel(_rdir(root, rid))
        assert _wait_for(lambda: runner.status(_rdir(root, rid))["status"] != "running")


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
        seen["pending_at_spawn"] = runner.pending_path(_rdir(root, rid)).exists()
        return real_popen(cmd, *a, **kw)

    monkeypatch.setattr(runner.subprocess, "Popen", spy)
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout", **_at(root))
    monkeypatch.undo()
    assert seen == {"lock_exists": True, "held_at_spawn": True, "pending_at_spawn": False}, seen
    lock = work / rid / "runner.lock"
    info = json.loads(runner.pending_path(_rdir(root, rid)).read_text())
    assert info["lock"] == str(lock)
    assert runner.status(_rdir(root, rid))["status"] == "running"
    assert held(lock)
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    runner.cancel(_rdir(root, rid))
    assert _wait_for(lambda: runner.status(_rdir(root, rid))["status"] != "running")
    assert runner.status(_rdir(root, rid))["status"] == "stopped"


# --- W-04: worker threads, exit and kill labelling ------------------------------------------------------


@pytest.mark.covers("W-04", "S-06")
def test_w04_worker_thread_pools_follow_the_runs_cpu_limit(world, record_property):
    import os

    root, d, _ = world
    hit, _ = _seed(root)
    row = runner.run(root, _recipe(hit, hit, tool="jackhmmer", mx=3), REL, work_base=d / "work",
                     limits=runner.Limits(cpus=2, wall_s=120), keep_work=True, **_at(root))
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
                     fault=fault, keep_work=True, **_at(root))
    record_property("run_wall_s", row["wall_s"])
    assert row["status"] == "failed", row["error"]
    assert row["error"].endswith(last_line), row["error"]
    assert row["wall_s"] < 20
    assert sorted(p.name for p in _rdir(root, row["run_id"]).iterdir()) == ["run.parquet"]
    assert runner.read_hits(_rdir(root, row["run_id"])) == []
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
    rid = runner.submit(root, rec, REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout", **_at(root))
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    os.kill(int((work / rid / "worker.pid").read_text()), signal.SIGKILL)
    assert _wait_for(lambda: runner.status(_rdir(root, rid))["status"] != "running")
    st = runner.status(_rdir(root, rid))
    assert st["status"] == "failed", st
    assert "signal 9" in st["error"] and "outside" in st["error"]
    assert sorted(p.name for p in _rdir(root, rid).iterdir()) == ["run.parquet"]


# --- C-03 and C-08 for search outputs -----------------------------------------------------------------


def test_c03_search_links_resolve_and_a_foreign_run_id_is_caught(world, tmp_path):
    root = world[0]
    hit, _ = _seed(root)
    _run(world, _recipe(hit, hit, tool="jackhmmer", mx=2))
    man = catalog.build_manifest(root, REL, ["S0"])
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"), runs=_closed_runs(root),
                               pack=pack.find(_packs(world), REL))
    try:
        assert contracts.check_search_links(con) == []
        assert con.execute("SELECT count(*) FROM search_hit").fetchone()[0] > 0
    finally:
        con.close()
    bad = _copy(world, tmp_path)
    hits_file = sorted(_runs_of(bad).glob("run=*/hits.parquet"))[0]
    t = pq.read_table(hits_file).to_pylist()
    for r in t:
        r["run_id"] = b"\x11" * 16
    pq.write_table(pa.Table.from_pylist(t, schema=table("search_hit").schema), hits_file)
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(bad, REL, ["S0"]), tmp_path / "d.duckdb"),
                               runs=_closed_runs(bad), pack=pack.find(_packs_of(bad), REL))
    try:
        problems = contracts.check_search_links(con)
    finally:
        con.close()
    assert any(p.startswith("search_hit -> search_run") for p in problems)
    assert any(p.startswith("run_id equals its directory") for p in problems)


def test_c08_search_runs_and_pack_record_the_locked_versions(world):
    from biobox import provenance

    root = world[0]
    assert sorted(_runs_of(root).glob("run=*/run.parquet"))
    # card W-03: the packs and the runs are read from their own directories, outside the data root
    assert provenance.check_search(_packs(world), _runs(world)) == []
    lock = {**provenance.locked_versions(), "pyhmmer": "0.0.0"}
    problems = provenance.check_search(_packs(world), _runs(world), lock)
    assert any("pack_manifest.json" in p for p in problems) and any("run.parquet" in p for p in problems)


# --- W-03: the pack and the runs live outside the data root ------------------------------------------------------


def _tree(root: Path) -> dict[str, str]:
    """sha256 of every file under `root`, and every directory, by root-relative path."""
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else "dir"
            for p in sorted(root.rglob("*"))}


def _read_only(root: Path, request) -> None:
    """Take write permission from every file and directory under `root` (and `root`); a finalizer gives it back."""
    modes = {p: p.stat().st_mode for p in [root, *root.rglob("*")]}

    def restore():
        for p, m in modes.items():
            p.chmod(m)

    request.addfinalizer(restore)
    for p, m in modes.items():
        p.chmod(m & ~0o222)


def _w03_recipe(root: Path, tmp_path: Path) -> Path:
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit, tool="jackhmmer", mx=2, extra="context: {window: 3g, edge: truncate}\n"))
    return rec


def _stat_tree(d: Path) -> dict[str, tuple]:
    """(inode, mtime, sha256) of every file under `d`: a rewrite of the same bytes changes the first two."""
    return {p.relative_to(d).as_posix(): (p.stat().st_ino, p.stat().st_mtime_ns,
                                          hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None)
            for p in sorted(d.rglob("*"))}


@pytest.mark.covers("W-03", "S-13")
@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores file permissions, so a read-only data root cannot be made")
def test_w03_a_full_search_against_a_read_only_data_root_writes_nothing_under_it(world, tmp_path, capsys, request):
    root = _copy(world, tmp_path)
    for d in (_packs_of(root), _runs_of(root)):  # the copy starts with no pack and no run of its own
        shutil.rmtree(d, ignore_errors=True)
    packs, runs, work = tmp_path / "packs", tmp_path / "runs", tmp_path / "work"
    before = _tree(root)
    build = ["pack", "build", "--root", str(root), "--release", REL, "--model", f"test-1={world[1] / 'a.hmm'}",
             "--pack-dir", str(packs)]
    assert cli.main(build) == 0
    assert _tree(root) == before  # building the pack wrote nothing under the root
    _read_only(root, request)
    with pytest.raises(PermissionError):
        (root / "probe").write_text("")  # the root really is read-only
    with pytest.raises(PermissionError):
        (root / "l0seq" / "probe").write_text("")
    capsys.readouterr()
    assert cli.main(["pack", "verify", "--root", str(root), "--release", REL, "--pack-dir", str(packs)]) == 0
    assert capsys.readouterr().out.strip() == f"pack {REL}: verified"
    where =["--root", str(root), "--runs", str(runs), "--pack-dir", str(packs)]
    assert cli.main(["search", "run", str(_w03_recipe(root, tmp_path)), "--release", REL, "--wait", "--json",
                     "--cpus", "1", "--work", str(work), *where]) == 0
    run = json.loads(capsys.readouterr().out.splitlines()[0])
    rid = run["run_id"]
    assert run["status"] == "ok"
    assert cli.main(["search", "status", rid, *where]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "ok"
    d = runner.run_dir(runs, rid)
    assert cli.main(["search", "status", str(d), "--root", str(root)]) == 0  # a run is also found by its directory
    assert json.loads(capsys.readouterr().out)["run_id"] == rid
    assert cli.main(["search", "hits", rid, "--json", *where]) == 0
    hits = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    assert hits and all(h["schema"] == "search_hit/2" for h in hits)
    for fmt in ("tsv", "faa"):
        assert cli.main(["search", "hits", rid, "--format", fmt, *where]) == 0
        assert capsys.readouterr().out.strip(), fmt
    assert cli.main(["nbr", rid, "--window", "3g", "--format", "tsv", *where]) == 0
    assert capsys.readouterr().out.strip()
    assert cli.main(["nbr", "summarize", rid, "--window", "5g", "--by", "pfam,nr50", *where]) == 0
    assert capsys.readouterr().out.strip()
    assert _tree(root) == before  # every file byte-identical, and no file or directory added
    assert sorted(p.name for p in d.iterdir()) == ["hits.parquet", "neighborhoods.parquet", "profiles.parquet",
                                                   "run.parquet"]


@pytest.mark.covers("W-03")
def test_w03_status_and_cancel_find_a_run_in_its_own_runs_directory_only(world, tmp_path, capsys):
    root = world[0]
    other, work = tmp_path / "other-runs", tmp_path / "work"
    rid = runner.submit(root, _w03_recipe(root, tmp_path), REL, runner.Limits(cpus=1, wall_s=300), work, fault="timeout",
                        packs=_packs(world), runs=other)
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    d = runner.run_dir(other, rid)
    mine = ["--root", str(root), "--runs", str(_runs(world))]
    assert cli.main(["search", "status", rid, *mine]) == 2  # not in the world's runs directory
    assert cli.main(["search", "cancel", rid, *mine]) == 2
    assert "no pending run" in capsys.readouterr().err
    assert runner.status(d)["status"] == "running"  # the refused cancel signalled nothing
    assert cli.main(["search", "cancel", rid, "--root", str(root), "--runs", str(other)]) == 0
    assert _wait_for(lambda: runner.status(d)["status"] != "running")
    assert runner.status(d)["status"] == "stopped"
    assert sorted(p.name for p in d.iterdir()) == ["run.parquet"]
    assert runner.pending_path(d).parent == other / "pending"  # its pending record is beside it, not in the root
    assert not (root / "search").exists()


def _w03_profile(tmp_path: Path) -> Path:
    data = {"max_threads_day": 4, "max_threads_night": 4, "mem_scope": "4G", "tmp": str(tmp_path / "ptmp"),
            "durable": str(tmp_path / "durable"), "bulk": str(tmp_path / "bulk"), "nice_monolithic": 19,
            "nice_fanout": 15, "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": 300,
            "small_file_floor_bytes": 1000, "require_mem_cap": False}
    path = tmp_path / "profile.yaml"
    path.write_text(json.dumps(data))  # JSON is YAML
    return path


@pytest.mark.covers("W-03")
def test_w03_two_users_runs_directories_do_not_see_each_others_runs(world, tmp_path, capsys, monkeypatch):
    root, prof = world[0], _w03_profile(tmp_path)
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    monkeypatch.setenv("LOGNAME", "alice")
    assert cli.main(["search", "run", str(_w03_recipe(root, tmp_path)), "--root", str(root), "--release", REL,
                     "--wait", "--cpus", "1", "--profile", str(prof), "--pack-dir", str(_packs(world))]) == 0
    rid = capsys.readouterr().out.split()[0]
    ptmp = tmp_path / "ptmp"
    assert (ptmp / "search/alice/runs" / f"run={rid}" / "run.parquet").is_file()  # the per-user default
    assert (ptmp / "search/alice/work").is_dir()  # the work directory is per user under the profile's tmp too
    assert cli.main(["search", "status", rid, "--root", str(root), "--profile", str(prof)]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "ok"
    monkeypatch.setenv("LOGNAME", "bob")
    for action in ("status", "hits", "cancel"):
        assert cli.main(["search", action, rid, "--root", str(root), "--profile", str(prof)]) == 2, action
        assert str(ptmp / "search/bob/runs") in capsys.readouterr().err, action
    assert not (root / "search").exists()


@pytest.mark.covers("W-03", "C-12", "D-01")
def test_w03_pack_build_twice_is_a_no_op_and_a_pack_is_never_overwritten(world, tmp_path, capsys, monkeypatch):
    root, d, man = world
    packs = tmp_path / "packs"
    build = ["pack", "build", "--root", str(root), "--release", REL, "--model", f"test-1={d / 'a.hmm'}",
             "--pack-dir", str(packs)]
    verify = ["pack", "verify", "--root", str(root), "--release", REL, "--pack-dir", str(packs)]
    before = _tree(root)
    assert cli.main(build) == 0
    first = capsys.readouterr().out
    pdir = pack.find(packs, REL)
    assert pdir == pack.pack_dir(packs, REL, man["pack_id"])  # the same inputs give the world's pack_id
    assert _pack_bytes(packs) == _pack_bytes(_packs(world))  # and its bytes
    snap = _stat_tree(packs)

    def no_build(*a, **k):
        raise AssertionError("the second build built")

    monkeypatch.setattr(pack, "_build_into", no_build)
    assert cli.main(build) == 0  # a pack of the same inputs is there and verifies: nothing is built or written
    assert capsys.readouterr().out == first
    assert _stat_tree(packs) == snap
    monkeypatch.undo()
    assert cli.main(verify) == 0
    assert capsys.readouterr().out.strip() == f"pack {REL}: verified"
    # one pack file altered: verify fails, and the build will not write over it
    p = pdir / "unique/part-0.faa.zst"
    b = bytearray(p.read_bytes())
    b[-5] ^= 0xFF
    p.write_bytes(bytes(b))
    assert cli.main(verify) == 1
    assert "sha256 differs" in capsys.readouterr().out
    assert cli.main(build) == 2
    assert "never overwritten" in capsys.readouterr().err
    assert p.read_bytes() == bytes(b)
    assert [x.name for x in pack.release_dir(packs, REL).iterdir()] == [pdir.name]  # no temporary pack left
    assert _tree(root) == before


@pytest.mark.covers("W-03")
def test_w03_a_pack_run_or_work_directory_under_the_data_root_is_refused(world, tmp_path, capsys):
    root = world[0]
    before = _tree(root)
    assert cli.main(["pack", "build", "--root", str(root), "--release", REL, "--pack-dir", str(root / "packs")]) == 2
    assert "data root" in capsys.readouterr().err
    with pytest.raises(ValueError, match="data root"):
        pack.build(root, REL, packs=root / "l0seq" / "packs")
    rec = _w03_recipe(root, tmp_path)
    base = ["search", "run", str(rec), "--root", str(root), "--release", REL, "--cpus", "1", "--pack-dir",
            str(_packs(world)), "--runs", str(tmp_path / "runs"), "--work", str(tmp_path / "work")]
    for wait in (["--wait"], []):
        for flag, inside in (("--runs", root / "runs"), ("--out", root / "run"), ("--work", root / "work"),
                             ("--pack-dir", root / "packs")):
            assert cli.main(base + wait + [flag, str(inside)]) == 2, (flag, wait)
            assert "data root" in capsys.readouterr().err, (flag, wait)
    assert _tree(root) == before
    assert not (tmp_path / "runs").exists() and not (tmp_path / "work").exists()  # refused before anything


@pytest.mark.covers("W-03")
def test_w03_a_search_with_no_pack_exits_2_naming_the_command_that_builds_it(world, tmp_path, capsys):
    root = world[0]
    empty = tmp_path / "no-packs"
    rec = _w03_recipe(root, tmp_path)
    base = ["search", "run", str(rec), "--root", str(root), "--release", REL, "--cpus", "1", "--pack-dir", str(empty),
            "--runs", str(tmp_path / "runs"), "--work", str(tmp_path / "work")]
    for wait in (["--wait"], []):
        assert cli.main(base + wait) == 2
        err = capsys.readouterr().err
        assert f"bio-box pack build --root {root.resolve()} --release {REL} --pack-dir {empty}" in err, err
    assert not empty.exists()  # never built implicitly
    assert not (tmp_path / "runs").exists()


# --- W-13: a run started with --wait is watched and cancelled as a submitted one is ------------------------------

W13_RID = "13" * 16  # given with --run-id, so the test knows where the run's records go before it starts


def _w13_wait_command(root: Path, rec: Path, runs: Path, work: Path) -> list[str]:
    """`bio-box search run --wait` as a separate process, the way `scripts/pod_run.sh` starts it in a unit: its
    worker sleeps until it is stopped (`fault=timeout`), under a wall far beyond the test."""
    return [sys.executable, "-m", "biobox.cli", "search", "run", str(rec), "--root", str(root), "--release", REL,
            "--wait", "--run-id", W13_RID, "--cpus", "1", "--wall-s", "300", "--work", str(work), "--fault", "timeout",
            "--runs", str(runs), "--pack-dir", str(_packs_of(root))]


@pytest.mark.covers("W-13")
def test_w13_a_run_started_with_wait_has_a_pending_record_and_is_cancelled_from_another_shell(world, tmp_path, capsys):
    """Card W-13 (2): measured on the pod, a `--wait` run had no pending record while it ran, so `search status`
    found nothing and `search cancel` could not be used. Now it writes the record `submit` writes (the same keys:
    the runner's pid and start time and the lock it holds), so another process sees `running` and cancels it."""
    root = world[0]
    hit, _ = _seed(root)
    rec = tmp_path / "r.yaml"
    rec.write_text(_recipe(hit, hit))
    runs, work = tmp_path / "runs", tmp_path / "work"
    d = runner.run_dir(runs, W13_RID)
    where = ["--root", str(root), "--runs", str(runs)]
    with (tmp_path / "stderr.txt").open("w") as err:
        proc = subprocess.Popen(_w13_wait_command(root, rec, runs, work), stdout=subprocess.PIPE, stderr=err, text=True)
        try:
            assert _wait_for(lambda: (work / W13_RID / "worker.pid").exists()), (tmp_path / "stderr.txt").read_text()
            info = json.loads(runner.pending_path(d).read_text())
            assert sorted(info) == ["lock", "pid", "run_id", "start", "status"]  # the keys `submit` writes
            assert (info["run_id"], info["pid"], info["status"]) == (W13_RID, proc.pid, "running")
            assert Path(info["lock"]).resolve() == (work / W13_RID / runner.LOCK_NAME).resolve()
            if platform.system() == "Linux":
                assert info["start"] == runner.proc_start_time(proc.pid)
            capsys.readouterr()
            assert cli.main(["search", "status", W13_RID, *where]) == 0
            assert json.loads(capsys.readouterr().out) == {"run_id": W13_RID, "status": "running"}
            # a second runner of the same run id is refused, and leaves the live run's records alone
            assert cli.main(_w13_wait_command(root, rec, runs, work)[3:]) == 2
            assert "running" in capsys.readouterr().err
            assert json.loads(runner.pending_path(d).read_text()) == info
            assert (work / W13_RID / "worker.pid").exists()
            assert cli.main(["search", "cancel", W13_RID, *where]) == 0
            out, _ = proc.communicate(timeout=60)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
    assert out.splitlines() == [f"{W13_RID}\tstopped"]
    assert proc.returncode == 1  # a run that did not end `ok`
    st = runner.status(d)
    assert (st["status"], st["error"]) == ("stopped", "cancelled"), st
    assert sorted(p.name for p in d.iterdir()) == ["run.parquet"]
    # the run id is in the log as soon as the run is pending, so a unit's log names it while it runs
    assert f"search run: {W13_RID} running" in (tmp_path / "stderr.txt").read_text()


@pytest.mark.covers("W-13")
def test_w13_a_closed_wait_run_reads_as_closed_and_cancel_signals_nothing(world, tmp_path, capsys):
    """The pending record goes when the run has closed (card W-20; migrated: it used to stay): status reads the
    closed run from its own record, and cancel finds no live runner, signals nothing and exits 0."""
    root = world[0]
    runs, work = tmp_path / "runs", tmp_path / "work"
    where = ["--root", str(root), "--runs", str(runs)]
    assert cli.main(["search", "run", str(_w03_recipe(root, tmp_path)), "--release", REL, "--wait", "--cpus", "1",
                     "--work", str(work), "--pack-dir", str(_packs(world)), *where]) == 0
    rid = capsys.readouterr().out.split()[0]
    d = runner.run_dir(runs, rid)
    assert not runner.pending_path(d).exists()
    assert not (work / rid).exists()
    assert cli.main(["search", "status", rid, *where]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "ok"
    assert cli.main(["search", "cancel", rid, *where]) == 0  # this process is the recorded pid: not signalled
    assert runner.status(d)["status"] == "ok"


@pytest.mark.covers("W-13")
def test_w13_the_runner_submit_starts_is_told_the_pending_record_is_submits(world, tmp_path, monkeypatch):
    """`submit` writes the pending record of the runner it starts (`search run --wait`); it passes `--submitted`, so
    that runner does not write one too. The record is the one `submit` wrote, and the run cancels as before."""
    root = world[0]
    runs, work, seen = tmp_path / "runs", tmp_path / "work", {}
    real_popen = subprocess.Popen

    def spy(cmd, *a, **kw):
        seen["cmd"] = list(cmd)
        return real_popen(cmd, *a, **kw)

    monkeypatch.setattr(runner.subprocess, "Popen", spy)
    rid = runner.submit(root, _w03_recipe(root, tmp_path), REL, runner.Limits(cpus=1, wall_s=300), work,
                        fault="timeout", packs=_packs(world), runs=runs)
    monkeypatch.undo()
    assert seen["cmd"][seen["cmd"].index("run"):].count("--wait") == 1 and "--submitted" in seen["cmd"]
    d = runner.run_dir(runs, rid)
    def stamp():
        st = runner.pending_path(d).stat()
        return st.st_ino, st.st_mtime_ns, runner.pending_path(d).read_text()

    written = stamp()
    assert _wait_for(lambda: (work / rid / "worker.pid").exists())
    assert stamp() == written  # the runner left it alone: a rewrite by rename would change the inode
    runner.cancel(d)
    assert _wait_for(lambda: runner.status(d)["status"] != "running")
    assert runner.status(d)["status"] == "stopped"


@pytest.mark.covers("W-13")
def test_w13_a_refused_wait_run_writes_no_pending_record(world, tmp_path, capsys):
    root = world[0]
    runs = tmp_path / "runs"
    d = runner.run_dir(runs, W13_RID)
    d.mkdir(parents=True)  # a run directory is written once
    assert cli.main(_w13_wait_command(root, _w03_recipe(root, tmp_path), runs, tmp_path / "work")[3:]) == 2
    assert "exists" in capsys.readouterr().err
    assert not runner.pending_path(d).exists()
    assert not (tmp_path / "work").exists()



@pytest.mark.covers("W-14")
def test_w14_peak_memory_is_the_workers_own_not_the_runners(world):
    """Linux hands a process's memory image at exec to the `ru_maxrss` of what it execs, so `wait4` on a worker
    forked from a big runner reports at least the runner's peak (card C-15). The runner here holds 768 MiB
    (written, so resident); the worker is a fresh interpreter that holds a small part of that."""
    held = b"\x01" * (768 << 20)
    hit, _ = _seed(world[0])
    row = _run(world, _recipe(hit, hit))
    assert len(held) == 768 << 20
    assert row["status"] == "ok", (row["error"], row["status"])
    print(f"recorded peak MiB: {row['peak_rss_bytes'] / 2**20:.0f}")
    assert 1 << 20 < row["peak_rss_bytes"] < 768 << 20, row["peak_rss_bytes"]


# --- W-07: one command ----------------------------------------------------------------------------------------


def _w07_labels(root: Path) -> dict[bytes, tuple[str, str]]:
    """(accession, header) of each contig's first source record by (l0_object_sha256, record_ordinal), by a loop
    over the contig_label files."""
    out: dict[bytes, tuple[str, str]] = {}
    for r in sorted(_pq(root, "l0seq/contig_label"), key=lambda r: (r["l0_object_sha256"], r["record_ordinal"])):
        out.setdefault(r["contig_uid"], (r["accession"], r["original_header"]))
    return out


def _w07_genes(root: Path) -> list[dict]:
    """Every gene of the store (protein_occurrence) with its contig's accession and header."""
    labels = _w07_labels(root)
    return [{**r, "accession": labels[r["contig_uid"]][0], "header": labels[r["contig_uid"]][1]}
            for r in _pq(root, "l1/protein_occurrence")]


def _w07_spec(root: Path, protein: bytes) -> str:
    """`ACCESSION:START-END` of the protein's first gene, as `search hits --format tsv` prints those columns."""
    g = min((g for g in _w07_genes(root) if g["protein_uid"] == protein), key=lambda g: (g["accession"], g["start"]))
    return f"{g['accession']}:{g['start']}-{g['end']}"


def _w07_hmm(world, tmp_path: Path) -> tuple[Path, bytes]:
    """The library's TEST0000 model as a user's HMM file, and a protein the library assigns to it."""
    root, d, _ = world
    return (_write_hmms(tmp_path / "in" / "lam0.hmm", _library_models(d / "a.hmm", ["TEST0000"])),
            sorted(_stored(root, "TEST0000"))[0])


def _w07_files(d: Path) -> dict[str, bytes]:
    return {p.relative_to(d).as_posix(): p.read_bytes() for p in sorted(d.rglob("*")) if p.is_file()}


def _w07_all_hits(d: Path) -> list[tuple]:
    """Every hit row of a run, every round and label, without the run id."""
    return sorted(tuple(sorted((k, v) for k, v in h.items() if k != "run_id"))
                  for h in runner.read_hits(d, final_only=False))


def _w07_separate(root: Path, packs: Path, runs: Path, rid: str, window: str, by: str, out: Path, capsys) -> None:
    """The separate commands of section 7 of the runbook, each writing what `search hmm` writes into `out`."""
    w = ["--root", str(root), "--pack-dir", str(packs), "--runs", str(runs)]
    assert cli.main(["search", "hits", rid, *w, "--format", "tsv", "--out", str(out / "hits.tsv")]) == 0
    assert cli.main(["search", "hits", rid, *w, "--format", "faa", "--out", str(out / "hits.faa")]) == 0
    assert cli.main(["nbr", rid, *w, "--window", window, "--format", "tsv", "--out", str(out / "loci")]) == 0
    capsys.readouterr()
    assert cli.main(["nbr", "summarize", rid, *w, "--window", window, "--by", by, "--format", "tsv"]) == 0
    (out / "neighbors.tsv").write_text(capsys.readouterr().out)  # as `> neighbors.tsv` in a shell


@pytest.mark.covers("W-07", "W-02", "W-05", "W-06")
def test_w07_search_hmm_writes_what_the_separate_commands_write_and_its_recipe_reruns(world, tmp_path, capsys):
    """The acceptance of card W-07 on the fixture: the one command's files are byte-identical to those of the
    separate commands (`search run`, `search hits` as TSV and FASTA, `nbr summarize`, `nbr --format tsv`) with the
    same recipe; and the recipe it writes into DIR, run again through `search run`, gives the same hits."""
    root = world[0]
    hmm, pos = _w07_hmm(world, tmp_path)
    spec = _w07_spec(root, pos)
    out, runs, work = tmp_path / "out", tmp_path / "runs", tmp_path / "work"
    code = cli.main(["search", "hmm", "--root", str(root), "--release", REL, "--hmm", str(hmm), "--bits", "25",
                     "--neighbors", "5", "--positive", spec, "--out", str(out), "--cpus", "1",
                     "--pack-dir", str(_packs(world)), "--runs", str(runs), "--work", str(work)])
    got = capsys.readouterr()
    assert code == 0, got.err
    rid, status = got.out.splitlines()[0].split("\t")
    assert status == "ok"
    assert sorted(p.name for p in out.iterdir()) == ["README.txt", "hits.faa", "hits.tsv", "loci", "neighbors.tsv",
                                                     "recipe.yaml"]  # card W-20: README.txt
    assert not list(tmp_path.glob(".out*"))  # no temporary directory left beside DIR
    # the recipe in DIR is the one the run ran: the file seed, the bits, the window, the resolved positive control
    run1 = runner.read_run(runner.run_dir(runs, rid))
    assert run1["recipe_yaml"] == (out / "recipe.yaml").read_text()
    r, _ = rcp.load(out / "recipe.yaml")
    assert r["seed"][0]["kind"] == "hmm_file" and Path(r["seed"][0]["value"]) == hmm.resolve()
    assert (r["tool"]["name"], r["iterations"]["max"], r["iterations"]["inclusion_bits"]) == ("hmmsearch", 1, 25.0)
    assert r["context"]["window"] == "5g"
    assert r["controls"]["positive"] == [{"kind": "protein_uid", "value": pos.hex(), "label": spec}]
    assert json.loads(run1["controls_json"])["positive"][0]["recovered"]
    # run again through `search run`: the same hits
    runs2 = tmp_path / "runs2"
    assert cli.main(["search", "run", str(out / "recipe.yaml"), "--root", str(root), "--release", REL, "--wait",
                     "--cpus", "1", "--pack-dir", str(_packs(world)), "--runs", str(runs2), "--work", str(work)]) == 0
    rid2, status2 = capsys.readouterr().out.split()
    assert status2 == "ok" and rid2 != rid
    hits1 = _w07_all_hits(runner.run_dir(runs, rid))
    assert hits1 and hits1 == _w07_all_hits(runner.run_dir(runs2, rid2))
    assert run1["recipe_hash"] == runner.read_run(runner.run_dir(runs2, rid2))["recipe_hash"]
    # the separate commands on that run: byte-identical files
    sep = tmp_path / "sep"
    _w07_separate(root, _packs(world), runs2, rid2, "5g", "pfam,nr50", sep, capsys)
    mine = _w07_files(out)
    assert mine.pop("recipe.yaml") and mine.pop("README.txt")  # card W-20: README.txt is `search hmm`'s own
    assert any(k.startswith("loci/") for k in mine) and len(mine["hits.tsv"].splitlines()) > 1
    assert pos.hex().encode() in mine["hits.faa"]
    assert mine == _w07_files(sep)


@pytest.mark.covers("W-07", "W-06")
def test_w07_search_hmm_on_a_release_without_a_domain_library_summarizes_at_the_cluster_levels(world, tmp_path, capsys):
    """A release may have no domain library (`rN.1` has none): the command still works, its neighbour summary is
    at the cluster levels only, and it says so."""
    root = _copy(world, tmp_path)
    shutil.rmtree(root / "l1/domains")
    assert not (root / f"manifests/release={REL}.json").exists()  # the release's files come from its release_def
    hmm, pos = _w07_hmm(world, tmp_path)
    out, runs = tmp_path / "out", tmp_path / "runs"
    code = cli.main(["search", "hmm", "--root", str(root), "--release", REL, "--hmm", str(hmm), "--bits", "25",
                     "--positive", _w07_spec(root, pos), "--out", str(out), "--cpus", "1",
                     "--pack-dir", str(_packs_of(root)), "--runs", str(runs), "--work", str(tmp_path / "work")])
    got = capsys.readouterr()
    assert code == 0, got.err
    rid = got.out.splitlines()[0].split("\t")[0]
    assert f"release {REL} has no domain library" in got.out and "nr50,nr90" in got.out, got.out
    summary = (out / "neighbors.tsv").read_text()
    assert '# levels\t["nr50","nr90"]' in summary.splitlines()
    tsv = [line.split("\t") for line in (out / "hits.tsv").read_text().splitlines()]
    col = tsv[0].index("pfam_architecture")
    assert len(tsv) > 1 and {row[col] for row in tsv[1:]} == {"-"}
    sep = tmp_path / "sep"
    _w07_separate(root, _packs_of(root), runs, rid, "5g", "nr50,nr90", sep, capsys)  # --neighbors defaults to 5
    mine = _w07_files(out)
    mine.pop("recipe.yaml")
    mine.pop("README.txt")  # card W-20
    assert mine == _w07_files(sep)


@pytest.mark.covers("W-07", "S-04")
def test_w07_a_positive_that_names_no_protein_or_several_is_refused_listing_what_matched(world, tmp_path, capsys):
    """The positive control is resolved to a protein uid through the release's own lookup (contig labels, then the
    contig's slots in the pack's tape version). Zero or several proteins: exit 2, nothing run or written, and the
    message lists what matched."""
    root = world[0]
    hmm, pos = _w07_hmm(world, tmp_path)
    out, runs = tmp_path / "out", tmp_path / "runs"
    base = ["search", "hmm", "--root", str(root), "--release", REL, "--hmm", str(hmm), "--bits", "25", "--out",
            str(out), "--cpus", "1", "--pack-dir", str(_packs(world)), "--runs", str(runs), "--work",
            str(tmp_path / "work")]
    genes = _w07_genes(root)
    lam = [g for g in genes if "phage lambda" in g["header"]]
    lam_proteins = {g["protein_uid"] for g in lam}
    assert len(lam_proteins) > 20

    def refused(spec: str) -> str:
        assert cli.main(base + ["--positive", spec]) == 2, spec
        assert not out.exists() and not runs.exists(), spec
        return capsys.readouterr().err

    err = refused("phage lambda")  # a header substring: every protein of the lambda contig
    assert f"matched {len(lam_proteins)} proteins" in err, err
    listed = set(re.findall(r"\b([0-9a-f]{32})\b", err))
    assert listed and listed <= {p.hex() for p in lam_proteins}
    g = min(lam, key=lambda g: g["start"])
    assert f"{g['accession']}:{g['start']}-{g['end']}" in err  # each line is the exact form to give instead
    err = refused(f"{g['accession']}:{g['start'] + 1}-{g['end']}")  # no gene there: the overlapping gene is listed
    assert "matched no protein" in err and f"{g['accession']}:{g['start']}-{g['end']}" in err, err
    assert "no contig" in refused("NC_999999.9:0-90")
    assert "matched no protein" in refused("no contig header holds this")
    # the exact gene resolves to its protein; at a threshold nothing passes the control is not recovered: the run
    # ends `failed` (exit 1) and DIR holds the recipe it ran, nothing else
    assert cli.main(base + ["--positive", _w07_spec(root, pos), "--bits", "100000"]) == 1
    assert "positive control" in capsys.readouterr().err
    assert sorted(p.name for p in out.iterdir()) == ["recipe.yaml"]
    r, _ = rcp.load(out / "recipe.yaml")
    assert r["controls"]["positive"][0]["value"] == pos.hex()


# --- W-19: a search holds its hits, not its target; it checks its pack by the build's record; it writes nothing
# into the pack directory; its pending record comes first ---------------------------------------------------------

_W19_RECORD = getattr(pack, "VERIFIED", "pack_verified.json")  # the build's verification record (card W-19)


@pytest.fixture(scope="module")
def w19_small(world, tmp_path_factory):
    """The world's release packed in small shards, so every searchable target has several."""
    root, d, _ = world
    packs = tmp_path_factory.mktemp("w19") / "packs"
    man = pack.build(root, REL, models={"test-1": d / "a.hmm"}, shard_text_bytes=20_000, packs=packs)
    return packs, man


def _w19_shards(packs: Path, man: dict, target: str) -> list[Path]:
    pdir = pack.pack_dir(packs, REL, man["pack_id"])
    return [pdir / s["path"] for s in man["shards"] if s["target"] == target]


def _w19_whole_block(hmm, records: list[tuple[str, str]], bits: float) -> dict:
    """The oracle: one pyhmmer search of every record held in one block (the worker's search before card W-19,
    written here separately), with every reported domain."""
    block = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=n.encode(), sequence=s) for n, s in records).digitize(engine.AMINO)
    out = {}
    for hits in pyhmmer.hmmsearch([hmm], block, cpus=1, Z=1e8, domZ=1e8, T=bits, incT=bits):
        for h in hits:
            if h.score >= bits:
                out[engine.text(h.name)] = (h.score, h.evalue, sorted(
                    (d.score, d.i_evalue, d.env_from, d.env_to, d.alignment.hmm_from, d.alignment.hmm_to)
                    for d in h.domains))
    return out


@pytest.mark.covers("W-19")
def test_w19_a_target_searched_block_by_block_scores_as_one_whole_block(world, w19_small):
    """Card W-19 (1): the target is read from the pack's shards a block at a time (`engine.stream_target`), not
    held. At the fixed Z every score and E-value is that of one search of the whole target: compared exactly, for
    a library model and a single-sequence profile, on both protein targets, with blocks of a few sequences."""
    root, d, _ = world
    packs, man = w19_small
    lib = engine.find_model([d / "a.hmm"], "TEST0000")
    prof = engine.build_from_sequence("w19", next(iter(domains.load_proteins(root, V, "S0").values())).rstrip("*"))
    for target, kind in (("unique", "protein"), ("dark_orf", "dark_orf")):
        paths = _w19_shards(packs, man, target)
        records = [r for p in paths for r in fastz.read_records(p)]
        assert len(paths) > 1 and len(records) == man["targets"][target]["n_seq"]
        for hmm, bits in ((lib, 10.0), (prof, 10.0)):
            want = _w19_whole_block(hmm, records, bits)
            t = engine.stream_target(kind, paths, man["targets"][target]["n_seq"], block_residues=500)
            got = engine.search(hmm, t, bits, -1e9, threads=2, keep_sequences=True)
            assert t.size == len(records) and t.blocks_read > 1, (target, t.blocks_read)
            assert {u.hex(): (h["seq_bitscore"], h["seq_evalue_zref"], sorted(
                (x["bitscore"], x["evalue_zref"], x["env_start"] + 1, x["env_end"], x["hmm_from"] + 1, x["hmm_to"])
                for x in h["domains"])) for u, h in got.items()} == want
            seqs = dict(records)
            assert {u.hex(): s for u, s in t.seqs.items()} == {n: seqs[n] for n in want}  # the hits' sequences only
        if target == "unique":
            assert want, "the library model finds no protein: the comparison would be empty"


@pytest.mark.covers("W-19")
def test_w19_a_run_in_blocks_of_a_few_sequences_gives_the_hits_of_one_block(world, tmp_path, monkeypatch):
    """Card W-19 (1), through the runner and the worker: an iterative run whose target is read in blocks of about
    300 residues has the hits, controls and rounds of the same run read in one block. Both runs have one run id,
    from which the negative control's shuffle is seeded, in runs directories of their own."""
    root, d, _ = world
    hit, _ = _seed(root)
    text = _recipe(hit, hit, tool="jackhmmer", mx=3)
    rid, rows = b"\x07" * 16, {}
    for block in (300, 1 << 30):
        monkeypatch.setattr(engine, "BLOCK_RESIDUES", block, raising=False)
        runs, work = tmp_path / f"runs{block}", tmp_path / f"work{block}"
        row = runner.run(root, text, REL, limits=LIMITS, work_base=work, packs=_packs(world), runs=runs, run_id=rid,
                         keep_work=True)
        assert row["status"] == "ok", row["error"]
        res = json.loads((work / rid.hex() / "out" / "result.json").read_text())
        rows[block] = (row, res, _w07_all_hits(runner.run_dir(runs, rid)))
    assert rows[300][1]["target_blocks"] > 1 and rows[1 << 30][1]["target_blocks"] == 1, \
        {b: r[1].get("target_blocks") for b, r in rows.items()}
    (a, _, ha), (b, _, hb) = rows[300], rows[1 << 30]
    assert ha == hb != []
    for k in ("rounds", "final_round", "stop_reason", "controls_json", "target_size", "rounds_json"):
        assert a[k] == b[k], k


def _w19_between(names: list[str]) -> str:
    """A name that sorts between two adjacent names of `names` (sorted), is none of them, and holds a letter (an
    all-digit name would be read as an integer by YAML)."""
    for x, y in zip(names, names[1:]):
        for k in range(1, min(int(y, 16) - int(x, 16), 64)):
            name = f"{int(x, 16) + k:032x}"
            if not name.isdigit():
                return name
    raise AssertionError("no gap between adjacent names")


@pytest.mark.covers("W-19", "S-04")
def test_w19_positive_controls_and_seeds_are_looked_up_without_a_set_of_every_name(world, w19_small, monkeypatch):
    """Card W-19 (2): the positive control, and a `protein_uid` seed, are looked up by key in the pack's shards
    (sorted by name: card W-16), not by reading every record into a set. The oracle reads every record whole;
    the runner then may not call `fastz.read_records` at all. Names at the start and end of each shard and in the
    middle are found, with their sequences; names before the first, after the last and between two records are
    refused, as before."""
    root = world[0]
    packs, man = w19_small
    oracle = {}
    for target in ("unique", "dark_orf"):
        recs = [list(fastz.read_records(p)) for p in _w19_shards(packs, man, target)]
        names = sorted(n for r in recs for n, _ in r)
        present = sorted({r[0][0] for r in recs} | {r[-1][0] for r in recs} | {recs[1][len(recs[1]) // 2][0]})
        absent = ["0" * 30 + "ab", "f" * 32, _w19_between(names)]  # before the first, after the last, between two
        assert not set(absent) & set(names)
        oracle[target] = (dict(x for r in recs for x in r), present, absent)
        for p in _w19_shards(packs, man, target):  # the shard-level lookup against the whole-shard oracle
            whole = dict(fastz.read_records(p))
            ask = set(present) | set(absent) | set(list(whole)[:3])
            assert fastz.lookup(p, ask) == {n: whole[n] for n in ask if n in whole}

    def whole_read(*a, **k):
        raise AssertionError("the runner read every record of a shard")

    monkeypatch.setattr(fastz, "read_records", whole_read)
    seqs, present, absent = oracle["unique"]
    other = present[0]
    for name in present:
        r, _ = rcp.load(_recipe(bytes.fromhex(name), bytes.fromhex(other)))
        res = runner.resolve(root, r, REL, packs)
        assert res["seed"]["sequence"] == seqs[name]
        r, _ = rcp.load(_recipe(bytes.fromhex(other), bytes.fromhex(name)))
        assert [p["uid"] for p in runner.resolve(root, r, REL, packs)["positives"]] == [name]
    for name in absent:
        r, _ = rcp.load(_recipe(bytes.fromhex(other), bytes.fromhex(name)))
        with pytest.raises(RecipeError, match="could never be recovered"):
            runner.resolve(root, r, REL, packs)
        r, _ = rcp.load(_recipe(bytes.fromhex(name), bytes.fromhex(other)))
        with pytest.raises(RecipeError, match="is not in release"):
            runner.resolve(root, r, REL, packs)
    seqs, present, absent = oracle["dark_orf"]
    for name, ok in [(n, True) for n in present] + [(n, False) for n in absent]:
        r, _ = rcp.load(f"recipe: dark\nwritten_against: {REL}\nseed: [{{sequence: {seqs[present[0]]}}}]\n"
                        f"target: dark_orf\ntool: phmmer\niterations: {{max: 1, inclusion_bits: 20}}\n"
                        f"controls: {{positive: [{{region_uid: {name}}}]}}\n")
        if ok:
            assert [p["uid"] for p in runner.resolve(root, r, REL, packs)["positives"]] == [name]
        else:
            with pytest.raises(RecipeError, match="could never be recovered"):
                runner.resolve(root, r, REL, packs)


def _w19_packs(world, tmp_path: Path, recorded: bool) -> Path:
    """A copy of the world's pack directory; without the build's verification record when not `recorded` (a pack
    built before card W-19)."""
    packs = tmp_path / "packs"
    shutil.copytree(_packs(world), packs)
    if not recorded:
        (pack.find(packs, REL) / _W19_RECORD).unlink(missing_ok=True)
    return packs


def _w19_snapshot(packs: Path) -> dict:
    """Every path under the pack directory with its mode, mtime and (for a file) bytes: a directory made and removed
    inside it changes its parent's mtime."""
    return {p.relative_to(packs).as_posix(): (p.stat().st_mode, p.stat().st_mtime_ns,
                                              p.read_bytes() if p.is_file() else None)
            for p in [packs, *sorted(packs.rglob("*"))]}


@pytest.mark.covers("W-19", "W-03", "C-12")
@pytest.mark.parametrize("recorded", [True, False], ids=["recorded", "built_before_w19"])
@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores file permissions, so a read-only pack cannot be made")
def test_w19_a_search_against_a_read_only_pack_directory_succeeds_and_writes_nothing_into_it(
        world, tmp_path, request, monkeypatch, capsys, recorded):
    """Card W-19 (3) and (4): a search checks a pack its build verified by the manifest and the file sizes only (no
    full verify: `pack.verify_dir` is not called); a pack built before W-19 has no record and is verified in full
    once for the run, with its scratch in the run's work directory, and stderr says so. Either way the pack
    directory is read-only and nothing in it changes."""
    root = world[0]
    packs = _w19_packs(world, tmp_path, recorded)
    _read_only(packs, request)
    with pytest.raises(PermissionError):
        (pack.release_dir(packs, REL) / "probe").write_text("")
    before = _w19_snapshot(packs)
    calls = []
    real = pack.verify_dir

    def counted(*a, **k):
        calls.append(1)
        return real(*a, **k)

    monkeypatch.setattr(pack, "verify_dir", counted)
    hit, _ = _seed(root)
    capsys.readouterr()
    row = runner.run(root, _recipe(hit, hit, tool="jackhmmer", mx=2), REL, limits=LIMITS, work_base=tmp_path / "work",
                     packs=packs, runs=tmp_path / "runs", pending=True)
    err = capsys.readouterr().err
    assert row["status"] == "ok", row["error"]
    assert _w19_snapshot(packs) == before
    assert len(calls) == (0 if recorded else 1)
    assert ("built before card W-19" in err) == (not recorded), err
    assert not (tmp_path / "work" / row["run_id"].hex()).exists()  # the full verify's scratch went with the run


@pytest.mark.covers("W-19", "C-12")
def test_w19_a_damaged_pack_is_still_refused_or_failed(world, tmp_path):
    """What the search-time check still catches. A file whose size differs from the manifest: refused, naming `pack
    verify`. A byte flipped in a target shard of a pack with a record (same size): refused naming `pack verify` when
    the key lookup cannot decompress the shard, else the run fails at staging, which hashes the staged shards
    against the manifest (zstd frames carry no checksum here, so which one depends on the byte). The same byte in a
    pack built before W-19: refused by the full verify."""
    root = world[0]
    hit, _ = _seed(root)
    text = _recipe(hit, hit)
    for recorded, how in ((True, "truncate"), (True, "flip"), (False, "flip")):
        sub = tmp_path / f"{recorded}-{how}"
        packs = _w19_packs(world, sub, recorded)
        p = pack.find(packs, REL) / "unique/part-0.faa.zst"
        b = bytearray(p.read_bytes())
        if how == "truncate":
            del b[-1]
        else:
            b[-5] ^= 0xFF
        p.write_bytes(bytes(b))
        kw = dict(limits=LIMITS, work_base=sub / "work", packs=packs, runs=sub / "runs")
        if recorded and how == "flip":
            try:
                row = runner.run(root, text, REL, **kw)
            except RecipeError as e:
                assert "pack verify" in str(e), e
            else:
                assert row["status"] == "failed" and "sha256 differs" in row["error"], row["error"]
        else:
            with pytest.raises(RecipeError, match="pack verify" if how == "truncate" else "does not verify"):
                runner.run(root, text, REL, **kw)


@pytest.mark.covers("W-19", "C-12", "D-01")
def test_w19_the_build_records_its_verification_beside_the_manifest_and_records_an_older_pack(
        world, tmp_path, monkeypatch):
    """Card W-19 (3): the build verifies the pack before it renames it into place and writes the record beside the
    manifest, not in it, so the manifest, and the pack id that is its hash, are those of a pack built before
    W-19. `pack build` on such a pack (the same inputs) builds nothing and adds the record; `pack verify` checks the
    record against the manifest."""
    root, d, man = world
    pdir = pack.find(_packs(world), REL)
    manifest, pid = pack.read_manifest(pdir)
    assert set(manifest) == {*pack.FINGERPRINT_KEYS, "targets", "shards", "meta"}  # no new key
    rec = json.loads((pdir / _W19_RECORD).read_text())
    assert rec["pack_id"] == pid.hex() == man["pack_id"] and rec["problems"] == []
    packs = _w19_packs(world, tmp_path, recorded=False)
    files = {k: v for k, v in _stat_tree(packs).items() if v[2] is not None}

    def no_build(*a, **k):
        raise AssertionError("built again")

    monkeypatch.setattr(pack, "_build_into", no_build)
    again = pack.build(root, REL, models={"test-1": d / "a.hmm"}, packs=packs)
    monkeypatch.undo()
    assert again["pack_id"] == man["pack_id"]
    after = {k: v for k, v in _stat_tree(packs).items() if v[2] is not None}
    new = f"release={REL}/{man['pack_id']}/{_W19_RECORD}"
    assert after.pop(new)[2] == hashlib.sha256((pdir / _W19_RECORD).read_bytes()).hexdigest()
    assert after == files  # every other file as it was: same inode, mtime and bytes
    assert pack.verify(root, REL, packs=packs) == []
    copy = pack.find(packs, REL) / _W19_RECORD
    copy.write_text(json.dumps({**rec, "pack_id": "00" * 16}))
    assert any(_W19_RECORD in p for p in pack.verify(root, REL, packs=packs))


@pytest.mark.covers("W-19", "W-13")
def test_w19_the_pending_record_exists_before_the_pack_check_and_the_lookups(world, tmp_path, monkeypatch):
    """Card W-19 (5): with a pending record (`search hmm`, `search run --wait`) the record, and so `search status`,
    exists before any long work: the full verify of a pack built before W-19 and the key lookups in the shards. A
    refusal found after the record is written removes it again."""
    root = world[0]
    packs = _w19_packs(world, tmp_path, recorded=False)
    runs, rid = tmp_path / "runs", b"\x19" * 16
    pend = runner.pending_path(runner.run_dir(runs, rid))
    seen = {"verify": [], "lookup": []}
    real_verify, real_lookup = pack.verify_dir, getattr(fastz, "lookup", None)

    def verify(*a, **k):
        seen["verify"].append(pend.exists())
        return real_verify(*a, **k)

    def lookup(*a, **k):
        seen["lookup"].append(pend.exists())
        return real_lookup(*a, **k)

    monkeypatch.setattr(pack, "verify_dir", verify)
    monkeypatch.setattr(fastz, "lookup", lookup, raising=False)
    hit, _ = _seed(root)
    row = runner.run(root, _recipe(hit, hit), REL, limits=LIMITS, work_base=tmp_path / "work", packs=packs, runs=runs,
                     run_id=rid, pending=True)
    assert row["status"] == "ok", row["error"]
    assert seen["verify"] == [True] and seen["lookup"] and all(seen["lookup"]), seen
    # a positive control not in the target, found once the record exists: refused, and the record is gone again
    rid2 = b"\x20" * 16
    with pytest.raises(RecipeError, match="could never be recovered"):
        runner.run(root, _recipe(hit, b"\xab" * 16), REL, limits=LIMITS, work_base=tmp_path / "work", packs=packs,
                   runs=runs, run_id=rid2, pending=True)
    assert not runner.pending_path(runner.run_dir(runs, rid2)).exists()
    assert not (tmp_path / "work" / rid2.hex()).exists() and not runner.run_dir(runs, rid2).exists()


# --- W-20: output wording and small defects from the second trial ---------------------------------------------


def _w20_run_hmm(world, tmp_path, root=None, packs=None):
    """`search hmm` on the world's (or `root`'s) release; returns (exit code, out dir, runs dir)."""
    root = world[0] if root is None else root
    hmm, pos = _w07_hmm(world, tmp_path)
    out, runs = tmp_path / "out", tmp_path / "runs"
    code = cli.main(["search", "hmm", "--root", str(root), "--release", REL, "--hmm", str(hmm), "--bits", "25",
                     "--positive", _w07_spec(world[0], pos), "--out", str(out), "--cpus", "1",
                     "--pack-dir", str(packs or _packs(world)), "--runs", str(runs), "--work", str(tmp_path / "work")])
    return code, out, runs


@pytest.mark.covers("W-20")
def test_w20_nbr_summarize_by_pfam_on_a_release_without_a_domain_library_says_so_and_drops_the_level(world, tmp_path,
                                                                                                    capsys):
    root = _copy(world, tmp_path)
    shutil.rmtree(root / "l1/domains")
    uids = ",".join(sorted(u.hex() for u in _w06_anchors(world[0])))
    args = ["nbr", "summarize", uids, "--window", "5g", "--release", REL, "--root", str(root), *_dirs(root),
            "--format", "tsv"]
    capsys.readouterr()
    assert cli.main([*args, "--by", "pfam,nr50"]) == 0
    got = capsys.readouterr()
    assert got.err.count("no domain library") == 1 and "pfam" in got.err, got.err
    assert '# levels\t["nr50"]' in got.out.splitlines() and "\tpfam\t" not in got.out
    assert cli.main([*args, "--by", "nr50"]) == 0  # the same answer without being told
    want = capsys.readouterr()
    assert want.out == got.out and "domain library" not in want.err
    assert cli.main([*args, "--by", "pfam"]) == 2  # nothing left to summarize
    assert "no domain library" in capsys.readouterr().err
    # a release with a library keeps the level and says nothing
    assert cli.main(["nbr", "summarize", uids, "--window", "5g", "--release", REL, "--root", str(world[0]),
                     *_dirs(world[0]), "--format", "tsv", "--by", "pfam,nr50"]) == 0
    got = capsys.readouterr()
    assert '# levels\t["pfam","nr50"]' in got.out.splitlines() and "domain library" not in got.err


@pytest.mark.covers("W-20")
def test_w20_the_slow_path_note_prints_once_per_command(world, tmp_path, capsys):
    root = _pre_b27(world, tmp_path, manifest=True, labels=True)
    capsys.readouterr()
    code, out, runs = _w20_run_hmm(world, tmp_path, root=root, packs=_packs_of(root))
    err = capsys.readouterr().err
    assert code == 0, err
    for name in ("member_cluster", "contig_label"):
        notes = [x for x in err.splitlines() if x.startswith("note:") and f"has no {name} files" in x]
        assert len(notes) == 1, (name, err)


@pytest.mark.covers("W-20")
def test_w20_search_hits_table_says_how_many_rows_of_how_many_it_shows(w05, capsys, monkeypatch):
    root, rid, _ = w05
    rows = _w05_final_hits(root, rid)
    n = 3 * len(rows)  # the fixture run has one kept hit: the table reads it three times over
    monkeypatch.setattr(runner, "read_hits", lambda d, round_=None, final_only=True: rows * 3)
    base = ["search", "hits", rid, "--root", str(root), *_dirs(root)]
    capsys.readouterr()
    assert cli.main([*base, "--top_n", "2"]) == 0
    got = capsys.readouterr()
    assert len(got.out.splitlines()) == 3  # the column line and two rows: stdout stays a clean table
    assert f"showing 2 of {n} rows" in got.err, got.err
    assert cli.main([*base, "--top_n", "0"]) == 0
    assert f"showing {n} of {n} rows" in capsys.readouterr().err


@pytest.mark.covers("W-20", "W-13")
def test_w20_the_pending_record_goes_when_the_run_closes_and_status_and_cancel_still_answer(world, tmp_path, capsys):
    root = world[0]
    runs, work = tmp_path / "runs", tmp_path / "work"
    where = ["--root", str(root), "--runs", str(runs)]
    assert cli.main(["search", "run", str(_w03_recipe(root, tmp_path)), "--release", REL, "--wait", "--cpus", "1",
                     "--work", str(work), "--pack-dir", str(_packs(world)), *where]) == 0
    rid = capsys.readouterr().out.split()[0]
    d = runner.run_dir(runs, rid)
    assert (d / "run.parquet").exists() and not runner.pending_path(d).exists()
    assert cli.main(["search", "status", rid, *where]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "ok"
    assert cli.main(["search", "cancel", rid, *where]) == 0  # a closed run: nothing to cancel, said so, not refused
    err = capsys.readouterr().err
    assert "closed" in err and "refused" not in err, err
    assert cli.main(["search", "cancel", "f" * 32, *where]) == 2  # a run that never existed is still refused
    assert "no pending run" in capsys.readouterr().err


@pytest.mark.covers("W-20", "W-07", "W-05")
def test_w20_search_hmm_writes_a_readme_and_one_strand_convention(world, tmp_path, capsys):
    capsys.readouterr()
    code, out, runs = _w20_run_hmm(world, tmp_path)
    err = capsys.readouterr().err
    assert code == 0, err
    readme = (out / "README.txt").read_text()
    # the facts a reader needs, and every column of every table named in it
    assert "0-based" in readme and "half-open" in readme and "+1" in readme and "-1" in readme
    assert "n_records" in readme and "anchor's strand" in readme and "dist_to_contig_end_bp" in readme
    words = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", readme))
    neighbors_cols = set(export.SUMMARY_TSV_COLUMNS["pfam"]) | set(export.SUMMARY_TSV_COLUMNS["nr"]) \
        | set(export.SUMMARY_TSV_STATS)
    for col in [*export.HITS_TSV_COLUMNS, *export.LOCUS_TSV_COLUMNS, *sorted(neighbors_cols)]:
        assert col in words, f"README.txt does not describe {col}"
    # strand is +1 or -1 in every exported file
    tsv = [x.split("\t") for x in (out / "hits.tsv").read_text().splitlines()]
    col = tsv[0].index("strand")
    assert len(tsv) > 1 and {r[col] for r in tsv[1:]} <= {"+1", "-1"}
    loci = sorted((out / "loci").iterdir())
    assert loci
    for f in loci:
        rows = [x.split("\t") for x in f.read_text().splitlines() if not x.startswith("#")]
        c = rows[0].index("strand")
        assert {r[c] for r in rows[1:]} <= {"+1", "-1"}, f
    faa = (out / "hits.faa").read_text()
    assert not re.search(r"\([+-]1\)", faa) and re.search(r" strand=[+-]1", faa)
    assert sorted(p.name for p in out.iterdir()) == ["README.txt", "hits.faa", "hits.tsv", "loci", "neighbors.tsv",
                                                     "recipe.yaml"]


@pytest.mark.covers("W-20", "W-05")
def test_w20_tier_is_a_column_only_when_the_release_has_tiers(w05, world):
    root, rid, dims = w05
    made = _w05_hits(root)
    tv = tapes.tape_version([V])
    with_tier = next(iter(export.hits_tsv(root, made, release=REL, tape_version=tv, gene_call_version=V))).split("\t")
    no_tier = next(iter(export.hits_tsv(world[0], made, release=REL, tape_version=tv, gene_call_version=V))).split("\t")
    assert "tier" in with_tier and "tier" not in no_tier
    assert [c for c in with_tier if c != "tier"] == no_tier
    # `dist_to_contig_end_bp` is the gene caller's: 0 for a gene that runs off a linear contig end, the contig's
    # length on a circular contig (README.txt says so), so it stays
    assert "dist_to_contig_end_bp" in no_tier


@pytest.mark.covers("W-20", "W-05")
def test_w20_the_fasta_example_names_the_positive_controls_contig(w05, tmp_path):
    """A protein on several contigs: its example is the first gene in sort order, unless a control names it. No
    protein of the fixture sits on two contigs, so a copy of it gets one: a gene of another contig is given the
    protein of a hit (protein_occurrence only; the exports read the tapes for ordinals and slots, not proteins)."""
    src, rid, dims = w05
    root = tmp_path / "a"
    shutil.copytree(src, root)
    made = _w05_hits(root)
    pu = made[0]["target_uid"]
    tables = {f: pq.read_table(f, partitioning=None)
              for f in sorted((root / "l1/protein_occurrence").rglob("*.parquet"))}
    home = {r["contig_uid"] for t in tables.values() for r in t.to_pylist() if r["protein_uid"] == pu}
    f, j = next((f, j) for f, t in tables.items() for j, r in enumerate(t.to_pylist())
                if r["contig_uid"] not in home and r["protein_uid"] != pu)
    col = tables[f].column("protein_uid").to_pylist()
    col[j] = pu
    tables[f] = tables[f].set_column(tables[f].schema.get_field_index("protein_uid"), "protein_uid",
                                     pa.array(col, tables[f].schema.field("protein_uid").type))
    pq.write_table(tables[f], f)
    gs = [g for g in _w05_oracle_genes(root, made, dims) if g["protein_uid"] == pu]
    assert len({g["accession"] for g in gs}) > 1, "the protein is on one contig"
    last = gs[-1]  # not the first in sort order
    assert last["accession"] != gs[0]["accession"]
    tv = tapes.tape_version([V])

    def example(prefer):
        lines = list(export.hits_faa(root, made, release=REL, tape_version=tv, gene_call_version=V, prefer=prefer))
        head = next(x for x in lines if x.startswith(">" + pu.hex()))
        return re.search(r"example=(\S+)", head).group(1)

    assert example(None).startswith(f"{gs[0]['accession']}:")
    assert example({pu: (last["accession"], None)}).startswith(f"{last['accession']}:")
    assert example({pu: (None, last["header"])}).startswith(f"{last['accession']}:")
    # what the CLI hands over: a control labelled with a span names its accession, any other label is header text
    controls = json.dumps({"positive": [{"kind": "protein_uid", "uid": pu.hex(), "label": "NC_001604:3170-5822"},
                                        {"kind": "protein_uid", "uid": (b"\x01" * 16).hex(), "label": "T7 gp1"},
                                        {"kind": "protein_uid", "uid": (b"\x02" * 16).hex(), "label": None}]})
    assert cli._positive_prefs(controls) == {pu: ("NC_001604", None), b"\x01" * 16: (None, "T7 gp1")}
    assert cli._positive_prefs(None) == {}

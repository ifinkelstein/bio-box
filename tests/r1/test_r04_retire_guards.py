"""Ring R1: `release retire` checks what the tar holds and what the copies say (card R-04, DEC-66 S-4, card C-54).

The oracle for "the tar holds the file" is `tarfile` in random-access mode and `hashlib` over the file the manifest
lists, written here independently of `catalog`'s streaming scan. Every refusal is checked to leave the root, and the
profile's pack directory, byte for byte as they were.
"""

from __future__ import annotations

import builtins
import hashlib
import io
import json
import os
import sys
import tarfile
import tracemalloc
from pathlib import Path

import pytest

from biobox import catalog, cli

sys.path.insert(0, str(Path(__file__).parent))
from test_c54_release_retire import _write_minimal_release, _write_minimal_stage  # noqa: E402

pytestmark = [pytest.mark.covers("R-04", "C-54", "LY-10")]

REL = "r26.11.1"
# the eight release= tables DEC-66 (S-4) names, written out here rather than read from catalog.RETIRED_TABLES
LARGE = {"cluster_member", "cluster_member_aln", "member_cluster", "token_dict", "tape_tokens", "token_postings",
         "contig_dim", "genome_dim"}
TOWER = "tower:/mnt/user/research-bulk/bio-box/root_backups"


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _file_sha(p: Path) -> str:
    return _sha(p.read_bytes())


def _world(tmp_path: Path, with_plan: bool = True) -> Path:
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")
    for rel in ("r26.11.1", "r26.11.2", "r26.11.3"):
        _write_minimal_release(root, rel, with_plan=with_plan)
    return root


def _profile(tmp_path: Path) -> tuple[Path, Path]:
    prof_tmp = tmp_path / "prof_tmp"
    pack = prof_tmp / "packs" / f"release={REL}" / ("a" * 32)
    pack.mkdir(parents=True)
    (pack / "pack_manifest.json").write_text("{}\n")
    prof = tmp_path / "profile.yaml"
    prof.write_text(
        f"max_threads_day: 2\nmax_threads_night: 2\nday_start_hour: 8\nday_end_hour: 20\n"
        f"mem_scope: 8G\nrequire_mem_cap: false\nnice_monolithic: 19\nnice_fanout: 15\n"
        f"cpu_min_per_process: 300\nsmall_file_floor_bytes: 1000\ntmp: {prof_tmp}\n"
        f"durable: {tmp_path / 'dur'}\nbulk: {tmp_path / 'bulk'}\n"
    )
    return prof, prof_tmp


def _snapshot(*dirs: Path) -> dict:
    """Every directory and every file's bytes (as sha256) under `dirs`."""
    out: dict[str, str] = {}
    for d in dirs:
        for cur, subdirs, files in os.walk(d):
            out[str(Path(cur))] = "<dir>"
            for name in files:
                p = Path(cur) / name
                out[str(p)] = _file_sha(p)
    return out


def _large_files(root: Path, release: str) -> list[dict]:
    man = catalog.read_manifest(catalog.manifest_path(root, release))
    return sorted(
        (f for f in man["files"]
         if f["table"] in LARGE and f["partitions"].get("release") == release and "stage" not in f["partitions"]),
        key=lambda f: f["path"],
    )


def _plan_files(root: Path, release: str) -> list[str]:
    d = root / "l2" / "cluster_plan" / f"release={release}"
    return sorted(p.relative_to(root).as_posix() for p in d.glob("*.json")) if d.is_dir() else []


def _release_paths(root: Path, release: str) -> list[str]:
    """Every file of the release's own partition plus its cluster plans, as a backup tar of the release holds them."""
    man = catalog.read_manifest(catalog.manifest_path(root, release))
    return sorted(f["path"] for f in man["files"] if f["partitions"].get("release") == release) + _plan_files(root, release)


def _tar(tar_path: Path, root: Path, paths, *, drop=(), replace=None, prefix: str = "", symlink=()) -> str:
    """Write a tar of `paths` under `root`; `drop` leaves members out, `replace` gives other bytes for a member,
    `symlink` stores a member as a symbolic link. Returns the tar's sha256."""
    replace = replace or {}
    with tarfile.open(tar_path, "w") as tf:
        for rel in paths:
            if rel in drop:
                continue
            name = prefix + rel
            if rel in symlink:
                info = tarfile.TarInfo(name)
                info.type = tarfile.SYMTYPE
                info.linkname = "elsewhere"
                tf.addfile(info)
            elif rel in replace:
                data = replace[rel]
                info = tarfile.TarInfo(name)
                info.size = len(data)
                tf.addfile(info, io.BytesIO(data))
            else:
                tf.add(root / rel, arcname=name)
    return _file_sha(tar_path)


def _retire(root: Path, tars, copies, *extra) -> int:
    argv = ["release", "retire", "--root", str(root), "--release", REL]
    for t in tars:
        argv += ["--tar", t]
    for c in copies:
        argv += ["--copy", c]
    return cli.main(argv + list(extra))


# --- the tar must hold every file to be removed, with the manifest's size and bytes -------------------------------


@pytest.mark.parametrize("fault", ["lacks_file", "other_size", "other_bytes", "other_release", "lacks_plan",
                                   "symlink_member", "partial_tar"])
def test_r04_a_tar_that_does_not_hold_the_release_is_refused_before_anything_is_removed(tmp_path, capsys, fault):
    root = _world(tmp_path)
    prof, prof_tmp = _profile(tmp_path)
    victim = next(f for f in _large_files(root, REL) if f["table"] == "token_dict")
    vbytes = (root / victim["path"]).read_bytes()
    tar_path = tmp_path / "backup.tar"
    paths = _release_paths(root, REL)
    if fault == "lacks_file":
        sha = _tar(tar_path, root, paths, drop={victim["path"]})
        want = "not in any --tar"
    elif fault == "other_size":
        sha = _tar(tar_path, root, paths, replace={victim["path"]: vbytes + b"\0"})
        want = f"holds {victim['bytes'] + 1} bytes, expected {victim['bytes']}"
    elif fault == "other_bytes":
        flipped = bytes([vbytes[0] ^ 0xFF]) + vbytes[1:]
        assert len(flipped) == victim["bytes"] and _sha(flipped) != victim["sha256"]
        sha = _tar(tar_path, root, paths, replace={victim["path"]: flipped})
        want = f"holds sha256 {_sha(flipped)}, expected {victim['sha256']}"
    elif fault == "other_release":
        sha = _tar(tar_path, root, _release_paths(root, "r26.11.2"))
        want = "not in any --tar"
    elif fault == "lacks_plan":
        plan = _plan_files(root, REL)
        assert plan, "the fixture release has a cluster plan"
        sha = _tar(tar_path, root, paths, drop={plan[0]})
        want = "not in any --tar"
    elif fault == "symlink_member":
        sha = _tar(tar_path, root, paths, symlink={victim["path"]})
        want = "not a regular file"
    else:  # a tar cut short in the middle of a member: a partial copy
        full = tmp_path / "full.tar"
        _tar(full, root, paths)
        data = full.read_bytes()
        # cut inside the data of one of the files to be removed
        with tarfile.open(full) as tf:
            m = tf.getmember(victim["path"])
            cut = m.offset_data + m.size // 2
        tar_path.write_bytes(data[:cut])
        sha = _file_sha(tar_path)
        want = "cannot be read as a tar"
    before = _snapshot(root, prof_tmp)
    rc = _retire(root, [f"{tar_path}={sha}"], [f"{TOWER}/backup.tar={sha}"], "--profile", str(prof))
    err = capsys.readouterr().err
    assert rc == 2, err
    assert want in err, err
    assert _snapshot(root, prof_tmp) == before
    assert not catalog.is_archived(root, REL)


def test_r04_a_file_held_twice_with_other_bytes_in_a_second_tar_is_refused(tmp_path, capsys):
    """Restoring extracts every recorded tar in turn, so a second tar holding other bytes would overwrite the first."""
    root = _world(tmp_path)
    victim = _large_files(root, REL)[0]
    good, bad = tmp_path / "a.tar", tmp_path / "b.tar"
    sha_good = _tar(good, root, _release_paths(root, REL))
    sha_bad = _tar(bad, root, [victim["path"]], replace={victim["path"]: b"x" * victim["bytes"]})
    before = _snapshot(root)
    rc = _retire(root, [f"{good}={sha_good}", f"{bad}={sha_bad}"], [f"{TOWER}/a.tar={sha_good}"])
    err = capsys.readouterr().err
    assert rc == 2 and f"{bad} holds sha256 {_sha(b'x' * victim['bytes'])}" in err, err
    assert _snapshot(root) == before


# --- --copy <location>=<sha256> ---------------------------------------------------------------------------------


@pytest.mark.parametrize("copies", [
    [],                                                  # no --copy at all
    [f"{TOWER}/backup.tar"],                             # no =<sha256>
    [f"{TOWER}/backup.tar=not-a-sha"],                   # not 64 hex
    [f"{TOWER}/backup.tar={'a' * 64}"],                  # a sha256 that is no --tar's
], ids=["none", "no_sha", "malformed_sha", "not_a_tars_sha"])
def test_r04_a_copy_without_a_tars_sha256_is_refused(tmp_path, capsys, copies):
    root = _world(tmp_path)
    prof, prof_tmp = _profile(tmp_path)
    tar_path = tmp_path / "backup.tar"
    sha = _tar(tar_path, root, _release_paths(root, REL))
    before = _snapshot(root, prof_tmp)
    rc = _retire(root, [f"{tar_path}={sha}"], copies, "--profile", str(prof))
    err = capsys.readouterr().err
    assert rc == 2, err
    assert "--copy" in err, err
    assert _snapshot(root, prof_tmp) == before


# --- a correct call ---------------------------------------------------------------------------------------------


def test_r04_a_correct_call_over_two_tars_retires_and_records_the_copies_as_the_operators_statement(tmp_path, capsys):
    root = _world(tmp_path)
    paths = _release_paths(root, REL)
    l2 = [p for p in paths if p.startswith("l2/") or p.startswith("dim/")]
    rest = [p for p in paths if p not in l2]
    ta, tb = tmp_path / "l2.tar", tmp_path / "l3.tar"
    sa = _tar(ta, root, l2)
    sb = _tar(tb, root, rest, prefix="./")  # as `tar -C <root> -cf l3.tar ./l3 ...` names them
    large = _large_files(root, REL)
    plan = _plan_files(root, REL)
    plan_bytes = {p: (root / p).read_bytes() for p in plan}
    h = catalog.manifest_hash(catalog.read_manifest(catalog.manifest_path(root, REL)))
    copies = [f"{TOWER}/l3.tar={sb}", f"{TOWER}/l2.tar={sa.upper()}"]
    assert _retire(root, [f"{ta}={sa}", f"{tb}={sb}"], copies) == 0
    capsys.readouterr()

    raw = catalog.archived_path(root, REL).read_text()
    rec = json.loads(raw)
    assert raw == json.dumps(rec, indent=2, sort_keys=True) + "\n"
    assert set(rec) == {"copies", "manifest_hash", "release", "removed_files", "removed_plan_files",
                        "schema_version", "tars"}
    assert rec["copies"] == [
        {"location": f"{TOWER}/l2.tar", "sha256": sa, "stated_by": "operator"},
        {"location": f"{TOWER}/l3.tar", "sha256": sb, "stated_by": "operator"},
    ]
    assert rec["tars"] == sorted([{"path": str(ta), "sha256": sa}, {"path": str(tb), "sha256": sb}],
                                 key=lambda t: t["path"])
    assert rec["manifest_hash"] == h and rec["schema_version"] == 1
    assert rec["removed_files"] == [{"path": f["path"], "sha256": f["sha256"], "size_bytes": f["bytes"]} for f in large]
    assert rec["removed_plan_files"] == [
        {"path": p, "sha256": _sha(b), "size_bytes": len(b)} for p, b in sorted(plan_bytes.items())
    ]
    for f in large:
        assert not (root / f["path"]).exists()
    assert catalog.verify_manifest(catalog.read_manifest(catalog.manifest_path(root, REL)), root) == []
    assert cli.main(["release", "verify", "--root", str(root), "--release", REL]) == 0
    assert f"manifest hash {h}" in capsys.readouterr().out

    # The same call again is the idempotent path: accepted, the record is not rewritten
    assert _retire(root, [f"{ta}={sa}", f"{tb}={sb}"], copies) == 0
    assert catalog.archived_path(root, REL).read_text() == raw
    capsys.readouterr()


def test_r04_each_tar_is_opened_once(tmp_path, capsys, monkeypatch):
    """One pass over a tar yields both its sha256 and every member's: the tar file is opened exactly once."""
    root = _world(tmp_path)
    tar_path = tmp_path / "backup.tar"
    sha = _tar(tar_path, root, _release_paths(root, REL))
    opened: list[str] = []
    real_open = builtins.open

    def counting_open(file, *a, **k):
        if isinstance(file, (str, os.PathLike)) and os.fspath(file) == str(tar_path):
            opened.append(os.fspath(file))
        return real_open(file, *a, **k)

    monkeypatch.setattr(builtins, "open", counting_open)
    monkeypatch.setattr(io, "open", counting_open)
    rc = _retire(root, [f"{tar_path}={sha}"], [f"{TOWER}/backup.tar={sha}"])
    monkeypatch.undo()
    assert rc == 0, capsys.readouterr().err
    assert len(opened) == 1, opened


def test_r04_the_tar_scan_reads_members_by_blocks(tmp_path):
    """Bounded memory: a 48 MiB member is hashed with a Python heap peak far below its size."""
    block = hashlib.sha256(b"r04").digest() * (1 << 15)  # 1 MiB, fixed
    n_blocks = 48
    size = len(block) * n_blocks

    class Gen(io.RawIOBase):
        def __init__(self):
            self.left = size

        def readable(self):
            return True

        def readinto(self, b):
            k = min(len(b), self.left, len(block))
            if k <= 0:
                return 0
            b[:k] = block[:k]
            self.left -= k
            return k

    tar_path = tmp_path / "big.tar"
    with tarfile.open(tar_path, "w") as tf:
        info = tarfile.TarInfo("l2/cluster_member/release=r26.11.1/part-0.parquet")
        info.size = size
        tf.addfile(info, io.BufferedReader(Gen(), buffer_size=1 << 16))
    want = hashlib.sha256()
    for _ in range(n_blocks):
        want.update(block)
    tracemalloc.start()
    try:
        tar_sha, found = catalog._scan_tar(tar_path, {"l2/cluster_member/release=r26.11.1/part-0.parquet": size})
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert tar_sha == _file_sha(tar_path)
    assert found == {"l2/cluster_member/release=r26.11.1/part-0.parquet": [("file", size, want.hexdigest())]}
    assert peak < 8 * (1 << 20), peak


# --- records written before this card ---------------------------------------------------------------------------


def _legacy_record(root: Path, release: str, tars: list[dict], copies: list[str]) -> str:
    """An archive record as card C-54's code wrote it (copies as plain strings), written by hand; removes the files."""
    removed = [{"path": f["path"], "sha256": f["sha256"], "size_bytes": f["bytes"]} for f in _large_files(root, release)]
    plans = [{"path": p, "sha256": _file_sha(root / p), "size_bytes": (root / p).stat().st_size}
             for p in _plan_files(root, release)]
    rec = {
        "copies": copies,
        "manifest_hash": catalog.manifest_hash(catalog.read_manifest(catalog.manifest_path(root, release))),
        "release": release,
        "removed_files": removed,
        "removed_plan_files": plans,
        "schema_version": 1,
        "tars": tars,
    }
    raw = json.dumps(rec, indent=2, sort_keys=True) + "\n"
    p = catalog.archived_path(root, release)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(raw)
    for item in removed + plans:
        (root / item["path"]).unlink()
    return raw


@pytest.mark.parametrize("shape", ["rN_three_tars_sha_in_copy", "r26_11_1_copy_without_sha"])
def test_r04_an_archive_record_written_before_this_card_still_verifies(tmp_path, capsys, shape):
    """The shapes of the records on the pod (2026-11-25): rN.1a to rN.3 and rM.2 (three tars, each copy a string
    `<location>=<sha256>`, no plan files) and r26.11.1 (one tar, the copy a bare location, 26 plan files)."""
    if shape.startswith("rN"):
        root = _world(tmp_path, with_plan=False)
        tars = [{"path": f"/stor/work/x/root_backups/{n}.tar", "sha256": c * 64}
                for n, c in (("l1_and_records", "6"), ("l2", "e"), ("l3", "5"))]
        copies = [f"{TOWER}/{n}.tar={c * 64}" for n, c in (("l1_and_records", "6"), ("l2", "e"), ("l3", "5"))]
    else:
        root = _world(tmp_path, with_plan=True)
        tars = [{"path": "/stor/work/x/root_backups/added_r26_11_1.tar", "sha256": "3" * 64}]
        copies = [f"{TOWER}/added_r26_11_1.tar"]
    raw = _legacy_record(root, REL, tars, copies)
    man = catalog.read_manifest(catalog.manifest_path(root, REL))
    assert catalog.verify_manifest(man, root) == []
    assert cli.main(["release", "verify", "--root", str(root), "--release", REL]) == 0
    out = capsys.readouterr().out
    assert f"release {REL}: verified (archived: " in out
    assert f"manifest hash {catalog.manifest_hash(man)}" in out
    msg = catalog.archived_error(root, REL)
    assert f"archived/release={REL}.json" in msg
    for t in tars:
        assert f"tar -C {root} -xf {t['path']}" in msg
    assert catalog.archived_path(root, REL).read_text() == raw


def test_r04_a_rerun_over_an_old_record_removes_a_straggler_only_when_the_tar_holds_it(tmp_path, capsys):
    """The idempotent path of a record written before this card: its tar was never checked, so a rerun that would
    remove a file left behind checks the tar first."""
    for name, drop, want_rc in (("partial", True, 2), ("full", False, 0)):
        d = tmp_path / name
        d.mkdir()
        root = _world(d)
        victim = _large_files(root, REL)[0]
        vbytes = (root / victim["path"]).read_bytes()
        tar_path = d / f"{name}.tar"
        sha = _tar(tar_path, root, _release_paths(root, REL), drop={victim["path"]} if drop else ())
        raw = _legacy_record(root, REL, [{"path": str(tar_path), "sha256": sha}], [f"{TOWER}/{name}.tar"])
        straggler = root / victim["path"]
        straggler.parent.mkdir(parents=True, exist_ok=True)
        straggler.write_bytes(vbytes)
        before = _snapshot(root)
        rc = _retire(root, [f"{tar_path}={sha}"], [f"{TOWER}/{name}.tar={sha}"])
        err = capsys.readouterr().err
        assert rc == want_rc, err
        if drop:  # the tar lacks the file: the straggler stays and nothing changes
            assert "not in any --tar" in err, err
            assert _snapshot(root) == before
        else:  # the tar holds it: the straggler is removed and the record is unchanged
            assert not straggler.exists()
            assert catalog.archived_path(root, REL).read_text() == raw
            assert catalog.verify_manifest(catalog.read_manifest(catalog.manifest_path(root, REL)), root) == []


def test_r04_verify_rejects_a_copy_whose_sha256_is_no_tars(tmp_path, capsys):
    root = _world(tmp_path)
    tar_path = tmp_path / "backup.tar"
    sha = _tar(tar_path, root, _release_paths(root, REL))
    assert _retire(root, [f"{tar_path}={sha}"], [f"{TOWER}/backup.tar={sha}"]) == 0
    capsys.readouterr()
    apath = catalog.archived_path(root, REL)
    raw = apath.read_text()
    rec = json.loads(raw)
    man = catalog.read_manifest(catalog.manifest_path(root, REL))
    assert catalog.verify_manifest(man, root) == []
    for bad in (
        [{**rec["copies"][0], "sha256": "b" * 64}],
        [{"location": rec["copies"][0]["location"], "sha256": sha}],  # no stated_by
        [{**rec["copies"][0], "location": ""}],
    ):
        apath.write_text(json.dumps({**rec, "copies": bad}, indent=2, sort_keys=True) + "\n")
        assert catalog.verify_manifest(man, root), bad
        assert cli.main(["release", "verify", "--root", str(root), "--release", REL]) == 1
    apath.write_text(raw)
    assert catalog.verify_manifest(man, root) == []
    capsys.readouterr()

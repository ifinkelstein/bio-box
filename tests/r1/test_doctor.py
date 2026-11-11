"""Ring R1: `biobox.doctor` produces a complete, JSON-serializable report."""

import json
import os

import pytest

from biobox import doctor


def test_report_has_every_section():
    r = doctor.report()
    assert set(r) >= {"python", "platform", "cpu", "memory", "ulimits", "filesystems", "tools", "packages"}
    assert r["packages"]["pyarrow"]
    json.dumps(r, default=str)


# ---- A-04: cgroup memory limit and path probes (fake /proc and /sys/fs/cgroup trees) ----

GIB = 2**30


def _fake_cgroup(tmp_path, cgroup_line, files):
    """Build fake `proc` and `cg` roots. `files` maps 'rel/dir/memory.max' to its text."""
    proc, cg = tmp_path / "proc", tmp_path / "cg"
    (proc / "self").mkdir(parents=True)
    cg.mkdir()
    (proc / "self" / "cgroup").write_text(cgroup_line)
    for rel, text in files.items():
        f = cg / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(text)
    return proc, cg


def test_cgroup_child_max_parent_number(tmp_path):
    proc, cg = _fake_cgroup(
        tmp_path,
        "0::/a/b\n",
        {"a/b/memory.max": "max\n", "a/memory.max": f"{4 * GIB}\n", "a/b/memory.current": "12345\n"},
    )
    r = doctor.cgroup_memory(proc_root=proc, cgroup_root=cg)
    assert r["available"] is True
    assert r["version"] == 2
    assert r["memory_max_bytes"] == 4 * GIB
    assert r["limited_at"] == "/a"
    assert r["memory_current_bytes"] == 12345


def test_cgroup_child_smaller_than_parent(tmp_path):
    proc, cg = _fake_cgroup(
        tmp_path,
        "0::/a/b\n",
        {"a/b/memory.max": f"{2 * GIB}\n", "a/memory.max": f"{8 * GIB}\n", "a/b/memory.high": f"{GIB}\n"},
    )
    r = doctor.cgroup_memory(proc_root=proc, cgroup_root=cg)
    assert r["memory_max_bytes"] == 2 * GIB
    assert r["limited_at"] == "/a/b"
    assert r["memory_high_bytes"] == GIB


def test_cgroup_all_max_is_unlimited(tmp_path):
    proc, cg = _fake_cgroup(
        tmp_path, "0::/a/b\n", {"a/b/memory.max": "max\n", "a/memory.max": "max\n", "a/b/memory.high": "max\n"}
    )
    r = doctor.cgroup_memory(proc_root=proc, cgroup_root=cg)
    assert r["available"] is True
    assert r["memory_max_bytes"] is None
    assert r["memory_high_bytes"] is None
    assert "unlimited" in r["memory_max"]


def test_cgroup_v1_is_unavailable_with_reason(tmp_path):
    line = "12:memory:/docker/abc\n1:name=systemd:/docker/abc\n"
    proc, cg = _fake_cgroup(tmp_path, line, {"memory/docker/abc/memory.limit_in_bytes": "1\n"})
    r = doctor.cgroup_memory(proc_root=proc, cgroup_root=cg)
    assert r["available"] is False
    assert "v1" in r["reason"]


def test_cgroup_missing_tree_is_unavailable(tmp_path):
    r = doctor.cgroup_memory(proc_root=tmp_path / "nope", cgroup_root=tmp_path / "nada")
    assert r["available"] is False
    assert r["reason"]


def test_cgroup_v2_path_without_files_is_unavailable(tmp_path):
    proc, cg = _fake_cgroup(tmp_path, "0::/a\n", {})
    r = doctor.cgroup_memory(proc_root=proc, cgroup_root=cg)
    assert r["available"] is False


def test_probe_writable_directory_leaves_nothing(tmp_path):
    d = tmp_path / "w"
    d.mkdir()
    r = doctor.probe_path(str(d), mounts_file=tmp_path / "no_mounts")
    assert r["exists"] and r["is_dir"] and r["readable"] and r["writable"]
    assert r["free_gib"] is not None
    assert list(d.iterdir()) == []


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores permission bits")
def test_probe_read_only_directory(tmp_path):
    d = tmp_path / "ro"
    d.mkdir()
    (d / "x").write_text("x")
    d.chmod(0o555)
    try:
        r = doctor.probe_path(str(d), mounts_file=tmp_path / "no_mounts")
        assert r["exists"] and r["is_dir"] and r["readable"]
        assert r["writable"] is False
    finally:
        d.chmod(0o755)
    assert [p.name for p in d.iterdir()] == ["x"]


def test_probe_missing_path_is_not_an_error_and_is_not_created(tmp_path):
    p = tmp_path / "stor" / "absent"
    r = doctor.probe_path(str(p), mounts_file=tmp_path / "no_mounts")
    assert r["exists"] is False
    assert r["status"] == "missing"
    assert r["writable"] is False
    assert not p.exists()
    assert not p.parent.exists()


def test_probe_file_is_not_a_directory(tmp_path):
    f = tmp_path / "f.txt"
    f.write_text("hi")
    r = doctor.probe_path(str(f), mounts_file=tmp_path / "no_mounts")
    assert r["exists"] is True
    assert r["is_dir"] is False
    assert r["readable"] is True
    assert sorted(p.name for p in tmp_path.iterdir()) == ["f.txt"]


def test_fs_type_longest_matching_mount_point(tmp_path):
    mounts = tmp_path / "mounts"
    mounts.write_text(
        "rootfs / ext4 rw 0 0\n"
        "nas:/vol /stor nfs4 rw 0 0\n"
        "nas:/vol2 /stor/scratch nfs rw 0 0\n"
        "tmpfs /stor/scratch/with\\040space tmpfs rw 0 0\n"
        "tmpfs /storage tmpfs rw 0 0\n"
    )
    assert doctor.fs_type("/etc/passwd", mounts) == "ext4"
    assert doctor.fs_type("/stor", mounts) == "nfs4"
    assert doctor.fs_type("/stor/work/x", mounts) == "nfs4"
    assert doctor.fs_type("/stor/scratch/Wilke", mounts) == "nfs"
    assert doctor.fs_type("/stor/scratch/with space/a", mounts) == "tmpfs"
    assert doctor.fs_type("/storage/a", mounts) == "tmpfs"
    assert doctor.fs_type("/x", tmp_path / "absent") is None


def test_probe_reports_fs_type_from_mounts(tmp_path):
    d = tmp_path / "w"
    d.mkdir()
    mounts = tmp_path / "mounts"
    mounts.write_text(f"rootfs / ext4 rw 0 0\nnas:/v {tmp_path} nfs4 rw 0 0\n")
    assert doctor.probe_path(str(d), mounts_file=mounts)["fs_type"] == "nfs4"


def _profile_text(tmp, durable, bulk):
    return (
        "max_threads_day: 40\nmax_threads_night: 56\nmem_scope: 200G\n"
        f"tmp: {tmp}\ndurable: {durable}\nbulk: {bulk}\n"
        "nice_monolithic: 19\nnice_fanout: 15\nday_start_hour: 8\nday_end_hour: 20\n"
        "cpu_min_per_process: 300\nsmall_file_floor_bytes: 67108864\n"
    )


def test_profile_supplies_the_three_paths(tmp_path):
    t, d = tmp_path / "t", tmp_path / "d"
    t.mkdir()
    d.mkdir()
    b = tmp_path / "missing_bulk"
    prof = tmp_path / "p.yaml"
    prof.write_text(_profile_text(t, d, b))
    assert doctor.profile_paths(prof) == [str(t), str(d), str(b)]
    r = doctor.report(profile=prof)
    assert set(r["paths"]) == {str(t), str(d), str(b)}
    assert r["paths"][str(t)]["writable"] is True
    assert r["paths"][str(b)]["status"] == "missing"
    assert not b.exists()
    json.dumps(r, default=str)


def test_report_has_cgroup_section():
    r = doctor.report()
    assert "cgroup" in r
    assert "available" in r["cgroup"]
    assert r["paths"] == {}


def test_main_parses_paths_and_profile(tmp_path, capsys):
    d = tmp_path / "d"
    d.mkdir()
    doctor.main(["--paths", str(d), str(tmp_path / "gone")])
    out = json.loads(capsys.readouterr().out)
    assert out["paths"][str(d)]["writable"] is True
    assert out["paths"][str(tmp_path / "gone")]["status"] == "missing"

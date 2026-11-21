"""Ring R1: a stage driver's mode, the release train and the release waiter report a failure (card A-11, 2026-11-21).

On 2026-11-20 `c47_driver.sh check` exited 0 after `c47_check.py` exited 1: the script ended `esac` and a `log` line,
whose status became the script's. The queue behind it ran `backup` and wrote "all passed", and a waiter started a
release on that line. Card C-47 fixed its own driver (tests/r1/test_c47_driver.py); this file does the same for every
other driver of the same form, for `release_train.sh` and for `rel_wait.sh`.

Each driver runs here from a copy whose `/tmp/if2453` (and the backup directory on /stor) is rewritten to the test's
temporary directory, with fake `systemd-run`, `systemctl`, `singularity`, `python3` (for the scripts only; inline code
runs on the real interpreter), `pod_run.sh`, `sleep`, `nice` and `ionice` first on PATH. A pattern in `FAKE_FAIL` makes
every stub command whose arguments match it exit 1: each mode's inner command fails in turn. Nothing of the pod's
`/tmp/if2453`, and nothing under /stor, is read or written.

The success paths print `A11CALLS <name> <sha256>`: the sha256 of every command line the stubs saw (the temporary
directory replaced by `<T>`). The card's acceptance compares those lines from the tree before the change with the tree
after it: a success path launches the same commands.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-11")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
ORCH = REPO / "scripts" / "orchestrator"
POD_BACKUPS = "/stor/work/Wilke/if2453/all-the-data/root_backups"

SYSTEMD_RUN = """#!/bin/sh
while [ $# -gt 0 ]; do
  case "$1" in --user|--scope|--quiet|--collect|--unit=*) shift ;; -p) shift 2 ;; *) break ;; esac
done
case "$*" in *xargs*) [ -n "${FAKE_HASH_RC:-}" ] && exit "$FAKE_HASH_RC" ;; esac
exec "$@"
"""
NICE = """#!/bin/sh
[ "$1" = -n ] && shift 2
exec "$@"
"""
IONICE = """#!/bin/sh
[ "$1" = -c ] && shift 2
exec "$@"
"""
NOOP = "#!/bin/sh\nexit 0\n"
INACTIVE = """#!/bin/sh
echo "systemctl $*" >> "$FAKE_LOG"
case "$*" in *is-active*) exit 3 ;; esac
exit 0
"""
SINGULARITY = """#!/bin/sh
echo "singularity $*" >> "$FAKE_LOG"
if [ -n "${FAKE_FAIL:-}" ] && echo "$*" | grep -Eq -- "$FAKE_FAIL"; then echo "fake failure: $*" >&2; exit 1; fi
case "$*" in
  *" release verify "*) echo "release R: verified, manifest hash abc123"; exit 0 ;;
  *"cli releases "*) echo "r26.11.2 listed"; exit 0 ;;
  *train_plan.py*--json*) cat "$FAKE_PLAN"; exit 0 ;;
  *train_plan.py*) echo "Action: cut"; exit 0 ;;
  *train_est.py*) echo '{}'; exit 0 ;;
  *" pack verify "*) echo "pack verified"; exit 0 ;;
esac
out=""; prev=""
for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
if [ -n "$out" ] && [ -z "${FAKE_NO_ADDED:-}" ]; then
  mkdir -p "$out"
  printf 'path\\tbytes\\tsha256\\nstage=X/f\\t1\\tab\\n' > "$out/root_added.tsv"
fi
echo "RESULT: fake"
exit 0
"""
PYTHON3 = """#!/bin/sh
case "${2:-}" in -c|-) exec "$REAL_PYTHON" "$@" ;; esac
echo "python3 $*" >> "$FAKE_LOG"
if [ -n "${FAKE_FAIL:-}" ] && echo "$*" | grep -Eq -- "$FAKE_FAIL"; then echo "fake failure: $*" >&2; exit 1; fi
case "$*" in
  *train_check.py*)
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
    mkdir -p "$out"
    if [ -z "${FAKE_NO_ADDED:-}" ]; then printf 'path\\tbytes\\tsha256\\nx\\t1\\tab\\n' > "$out/root_added.tsv"; fi
    echo "# announcement" > "$out/announcement.md"
    echo "check ok" ;;
  *) echo '{}' ;;
esac
exit 0
"""
POD_RUN = """#!/bin/sh
log=""; unit=""
while [ $# -gt 0 ]; do
  case "$1" in
    --log) log=$2; shift 2 ;;
    --*) shift 2 ;;
    *) unit=$1; shift; break ;;
  esac
done
echo "pod_run $unit $*" >> "$FAKE_LOG"
rel=""; prev=""
for a in "$@"; do [ "$prev" = "--release" ] && rel=$a; prev=$a; done
rc=0
if [ -n "${FAKE_POD_FAIL:-}" ] && echo "$unit" | grep -Eq -- "$FAKE_POD_FAIL"; then rc=1; fi
mkdir -p "$(dirname "$log")"
{
  [ "$rc" = 0 ] && [ -n "$rel" ] && echo "build-stage: stage S and release $rel: all 9 steps done"
  echo "pod_run: $unit exited $rc"
} >> "$log"
exit 0
"""
EX = """#!/bin/sh
echo "release_examples $*" >> "$FAKE_LOG"
[ -n "${FAKE_FAIL:-}" ] && echo "release_examples" | grep -Eq -- "$FAKE_FAIL" && exit 1
exit 0
"""

PLAN = (
    '{"release": "r26.11.2", "stages": ["S1", "P1"], "latest_release": "r26.11.1", "nothing_new": false,'
    ' "earlier_manifest_releases": ["r26.11.1"], "last_stage": {"stage": "P1", "depth": "called",'
    ' "genomes": "g.json", "objects": "o.json", "crispr": false}}'
)

# stage directories the drivers hash ("earlier" stages), so that every snapshot hashes at least one file
STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5")


class World:
    """One test's temporary pod: stubs, the rewritten driver copies, a root with files, the fake logs."""

    def __init__(self, tmp: Path):
        self.tmp = tmp
        self.pod = tmp / "if2453"
        self.fake = tmp / "fakebin"
        self.bin = tmp / "bin"
        self.tree = tmp / "tree"
        self.root = self.pod / "r10" / "root"
        self.log = tmp / "fake.log"
        self.stor = tmp / "stor"
        for d in (self.fake, self.bin, self.tree / "scripts", self.tree / "profiles", self.root):
            d.mkdir(parents=True, exist_ok=True)
        for name, text in (("systemd-run", SYSTEMD_RUN), ("nice", NICE), ("ionice", IONICE), ("sleep", NOOP),
                           ("top", NOOP), ("uptime", NOOP), ("systemctl", INACTIVE), ("singularity", SINGULARITY),
                           ("python3", PYTHON3)):
            (self.fake / name).write_text(text)
            (self.fake / name).chmod(0o755)
        (self.tree / "scripts" / "pod_run.sh").write_text(POD_RUN)
        (self.tree / "scripts" / "release_examples.sh").write_text(EX)
        (self.tree / "profiles" / "wilkcomp02.yaml").write_text("tmp: /x\nmax_threads_day: 40\n")
        (self.tmp / "plan.json").write_text(PLAN)
        for s in STAGES:
            d = self.root / "l0seq" / "x" / f"stage={s}"
            d.mkdir(parents=True, exist_ok=True)
            (d / "part-0.parquet").write_bytes(b"x")
            (self.root / "l0seq" / "seq_locator" / f"stage={s}").mkdir(parents=True, exist_ok=True)
        self.log.write_text("")

    def driver(self, name: str) -> Path:
        """A copy of scripts/.../<name> with the pod's paths rewritten to this test's directory."""
        src = (BSP / name).read_text()
        text = src.replace("/tmp/if2453", str(self.pod)).replace(POD_BACKUPS, str(self.stor))
        dst = self.bin / name
        dst.write_text(text)
        for other in BSP.iterdir():  # scripts the driver names beside itself ($BIN/x.py) stay reachable
            link = self.bin / other.name
            if other.name != name and not link.exists():
                link.symlink_to(other)
        return dst

    def env(self, **extra: str) -> dict:
        return {**os.environ, "PATH": f"{self.fake}:{os.environ['PATH']}", "FAKE_LOG": str(self.log),
                "FAKE_PLAN": str(self.tmp / "plan.json"), "REAL_PYTHON": sys.executable, "TREE": str(self.tree),
                "START_MIN_GB": "0", **extra}

    def run(self, script: Path, *args: str, **extra: str) -> subprocess.CompletedProcess:
        return subprocess.run(["bash", str(script), *args], env=self.env(**extra), capture_output=True, text=True,
                              timeout=300)

    def calls(self) -> str:
        return self.log.read_text().replace(str(self.tmp), "<T>")

    def record(self, name: str) -> None:
        print("A11CALLS", name, hashlib.sha256(self.calls().encode()).hexdigest())


def _ok(proc: subprocess.CompletedProcess, w: World | None = None, name: str = "") -> None:
    if w is not None and name:  # recorded first: the card compares these lines before and after its change
        w.record(name)
    assert proc.returncode == 0, proc.stdout + proc.stderr


# --------------------------------------------------------------------------------------------------------------
# the drivers of the form `case "$MODE" in ... esac; log "done"`: c24, c26, c33, c34, c35

OLD = {
    # name: (script, C directory name, the check script, modes that fail on a failed verify)
    "c24": ("c24_driver.sh", "c24", "c24_check.py"),
    "c26": ("c26_driver.sh", "c26", "c26_check.py"),
    "c33": ("c33_driver.sh", "c33", "c33_check.py"),
    "c34": ("c34_driver.sh", "c34", "c34_check.py"),
    "c35": ("c35_driver.sh", "c35", "c35_check.py"),
}
RELS = "rN.1a rN.1b rN.1 rN.2 rN.3 rM.2".split()


def _old(tmp_path: Path, key: str):
    script, cdir, chk = OLD[key]
    w = World(tmp_path)
    drv = w.driver(script)
    c = w.pod / cdir
    c.mkdir(parents=True, exist_ok=True)
    (c / "units").write_text("u-1\n")
    (c / "rehearsal").mkdir(exist_ok=True)
    (c / "estimate").mkdir(exist_ok=True)
    (c / "estimate" / "before.json").write_text("{}")
    (c / "before").mkdir(exist_ok=True)
    for r in RELS + ["r26.11.1"]:  # c35's after compares the manifest hash with the one in before/
        (c / "before" / f"verify-{r}.log").write_text("release R: verified, manifest hash abc123\n")
    return w, drv, c, chk


@pytest.mark.parametrize("key", sorted(OLD))
def test_old_driver_a_failed_check_is_the_exit_status(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    (c / "after").mkdir(exist_ok=True)
    (c / "after" / "DONE").touch()
    got = w.run(drv, "check", FAKE_FAIL=chk)
    assert got.returncode != 0, "check exited 0 after its check script exited 1"
    assert not (c / "check" / "DONE").exists()


@pytest.mark.parametrize("key", sorted(OLD))
def test_old_driver_a_passed_check_writes_its_marker(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    (c / "after").mkdir(exist_ok=True)
    (c / "after" / "DONE").touch()
    got = w.run(drv, "check")
    _ok(got, w, f"{key}-check")
    assert (c / "check" / "DONE").exists()


@pytest.mark.parametrize("key", sorted(OLD))
def test_old_driver_check_refuses_without_the_after_marker(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    got = w.run(drv, "check")
    assert got.returncode != 0
    assert not (c / "check" / "DONE").exists()


@pytest.mark.parametrize("key", sorted(OLD))
def test_old_driver_a_check_that_wrote_no_list_is_a_failure(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    (c / "after").mkdir(exist_ok=True)
    (c / "after" / "DONE").touch()
    got = w.run(drv, "check", FAKE_NO_ADDED="1")
    assert got.returncode != 0
    assert not (c / "check" / "DONE").exists()


@pytest.mark.parametrize("key", sorted(OLD))
def test_old_driver_a_failed_verify_after_the_build_is_the_exit_status(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    got = w.run(drv, "after", FAKE_FAIL=r" release verify ")
    assert got.returncode != 0
    assert not (c / "after" / "DONE").exists()


@pytest.mark.parametrize("key", sorted(OLD))
def test_old_driver_a_passed_after_writes_its_marker(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    got = w.run(drv, "after")
    _ok(got, w, f"{key}-after")
    assert (c / "after" / "DONE").exists()


@pytest.mark.parametrize("key", ["c24", "c26", "c33", "c34"])
def test_old_driver_a_failed_hash_in_the_snapshot_fails_the_mode(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    got = w.run(drv, "after", FAKE_HASH_RC="1")
    assert got.returncode != 0
    assert not (c / "after" / "DONE").exists()


@pytest.mark.parametrize("key", ["c24", "c26", "c34", "c35"])
def test_old_driver_a_failed_releases_listing_fails_after(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    got = w.run(drv, "after", FAKE_FAIL=r"cli releases ")
    assert got.returncode != 0
    assert not (c / "after" / "DONE").exists()


@pytest.mark.parametrize("key", sorted(OLD))
def test_old_driver_an_unknown_mode_exits_2(tmp_path, key):
    w, drv, c, chk = _old(tmp_path, key)
    assert w.run(drv, "nosuchmode").returncode == 2


def test_c35_a_failed_rates_read_fails_check(tmp_path):
    w, drv, c, chk = _old(tmp_path, "c35")
    (c / "after").mkdir(exist_ok=True)
    (c / "after" / "DONE").touch()
    got = w.run(drv, "check", FAKE_FAIL="c35_est.py")
    assert got.returncode != 0
    assert not (c / "check" / "DONE").exists()


# b48 and b49 are other forms: b48's `compare` ends on a log line, b49's `check` on `du`

def test_b48_a_failed_compare_is_the_exit_status(tmp_path):
    w = World(tmp_path)
    drv = w.driver("b48_driver.sh")
    a, b = w.pod / "ra", w.pod / "rb"
    a.mkdir()
    b.mkdir()
    got = w.run(drv, "compare", str(a), str(b), "n1", FAKE_FAIL="b48_compare.py")
    assert got.returncode != 0


def test_b48_a_passed_compare_exits_0(tmp_path):
    w = World(tmp_path)
    drv = w.driver("b48_driver.sh")
    a, b = w.pod / "ra", w.pod / "rb"
    a.mkdir()
    b.mkdir()
    _ok(w.run(drv, "compare", str(a), str(b), "n1"), w, "b48-compare")


def test_b49_a_failed_check_is_the_exit_status(tmp_path):
    w = World(tmp_path)
    drv = w.driver("b49_driver.sh")
    (w.pod / "b49").mkdir()
    got = w.run(drv, "check", FAKE_FAIL="b49_check.py")
    assert got.returncode != 0, "check exited 0 after b49_check.py exited 1 (the last command is `du`)"


def test_b49_a_passed_check_exits_0(tmp_path):
    w = World(tmp_path)
    drv = w.driver("b49_driver.sh")
    (w.pod / "b49").mkdir()
    _ok(w.run(drv, "check"), w, "b49-check")


# --------------------------------------------------------------------------------------------------------------
# c48c_v1_driver.sh (stage V1): every mode

V1_STEPS = ("inputs ingest genecall crispr-minced-0.4.2-w8-nr3-r23x47-s26x50 igr tapes domains-pfam-38.2 "
            "domains-pfam-38.2.plan crispr-minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2").split()


def _v1(tmp_path: Path):
    w = World(tmp_path)
    drv = w.driver("c48c_v1_driver.sh")
    c = w.pod / "c48c"
    c.mkdir(parents=True, exist_ok=True)
    (c / "units").write_text("c48c-v1-1\n")
    w.stor.mkdir(exist_ok=True)
    return w, drv, c


def _mark(c: Path, *names: str) -> None:
    for n in names:
        (c / n).parent.mkdir(parents=True, exist_ok=True)
        (c / n).touch()


def test_v1_before_a_failed_verify_is_the_exit_status(tmp_path):
    w, drv, c = _v1(tmp_path)
    got = w.run(drv, "before", FAKE_FAIL=r" release verify ")
    assert got.returncode != 0
    assert not (c / "before" / "DONE").exists()


def test_v1_before_a_failed_hash_in_the_snapshot_fails_the_mode(tmp_path):
    w, drv, c = _v1(tmp_path)
    got = w.run(drv, "before", FAKE_HASH_RC="1")
    assert got.returncode != 0, "the snapshot's sha256sum failed and `before` wrote DONE"
    assert not (c / "before" / "DONE").exists()


def test_v1_before_writes_its_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _ok(w.run(drv, "before"), w, "v1-before")
    assert (c / "before" / "DONE").exists()


def test_v1_estimate_with_inputs_that_differ_fails_and_writes_no_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    (w.pod / "s2" / "inputs" / "V1").mkdir(parents=True)
    got = w.run(drv, "estimate")
    assert got.returncode != 0
    assert not (c / "estimate" / "DONE").exists()


def test_v1_a_failed_estimate_withdraws_an_earlier_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _mark(c, "estimate/DONE")
    (w.pod / "s2" / "inputs" / "V1").mkdir(parents=True)
    assert w.run(drv, "estimate").returncode != 0
    assert not (c / "estimate" / "DONE").exists(), "run would start on the marker of an estimate that failed"


def test_v1_after_refuses_without_the_run_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    got = w.run(drv, "after")
    assert got.returncode != 0
    assert not (c / "after" / "DONE").exists()


@pytest.mark.parametrize("pattern,why", [
    (r" release verify ", "a release does not verify"),
    (r"cli releases ", "bio-box releases fails"),
])
def test_v1_after_failures_are_the_exit_status(tmp_path, pattern, why):
    w, drv, c = _v1(tmp_path)
    _mark(c, "run/DONE")
    got = w.run(drv, "after", FAKE_FAIL=pattern)
    assert got.returncode != 0, why
    assert not (c / "after" / "DONE").exists()


def test_v1_after_a_failed_hash_in_the_snapshot_fails_the_mode(tmp_path):
    w, drv, c = _v1(tmp_path)
    _mark(c, "run/DONE")
    got = w.run(drv, "after", FAKE_HASH_RC="1")
    assert got.returncode != 0
    assert not (c / "after" / "DONE").exists()


def test_v1_after_writes_its_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _mark(c, "run/DONE")
    _ok(w.run(drv, "after"), w, "v1-after")
    assert (c / "after" / "DONE").exists()


def test_v1_a_failed_after_withdraws_an_earlier_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _mark(c, "run/DONE", "after/DONE")
    assert w.run(drv, "after", FAKE_FAIL=r" release verify ").returncode != 0
    assert not (c / "after" / "DONE").exists()


def test_v1_check_refuses_without_the_after_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    got = w.run(drv, "check")
    assert got.returncode != 0
    assert not (c / "check" / "DONE").exists()


def test_v1_a_failed_check_is_the_exit_status(tmp_path):
    w, drv, c = _v1(tmp_path)
    _mark(c, "after/DONE")
    got = w.run(drv, "check", FAKE_FAIL="c48c_v1_check.py")
    assert got.returncode != 0
    assert not (c / "check" / "DONE").exists()


def test_v1_a_check_that_wrote_no_list_is_a_failure(tmp_path):
    w, drv, c = _v1(tmp_path)
    _mark(c, "after/DONE")
    got = w.run(drv, "check", FAKE_NO_ADDED="1")
    assert got.returncode != 0
    assert not (c / "check" / "DONE").exists()


def test_v1_a_failed_check_withdraws_an_earlier_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _mark(c, "after/DONE", "check/DONE")
    assert w.run(drv, "check", FAKE_FAIL="c48c_v1_check.py").returncode != 0
    assert not (c / "check" / "DONE").exists(), "backup would start on the marker of a check that failed"


def test_v1_a_passed_check_writes_its_marker_and_the_list(tmp_path):
    w, drv, c = _v1(tmp_path)
    _mark(c, "after/DONE")
    _ok(w.run(drv, "check"), w, "v1-check")
    assert (c / "check" / "DONE").exists()
    assert (c / "root_added.tsv").read_text() == (c / "check" / "root_added.tsv").read_text()


def _backup_ready(w: World, c: Path) -> None:
    _mark(c, "check/DONE")
    f = w.root / "added.bin"
    with open(f, "wb") as fh:
        fh.truncate(70 * 1024 * 1024)
    (c / "root_added.tsv").write_text("path\tbytes\tsha256\nadded.bin\t1\tab\n")


def test_v1_backup_refuses_without_the_check_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    (c / "root_added.tsv").write_text("path\tbytes\tsha256\n")
    got = w.run(drv, "backup")
    assert got.returncode != 0
    assert not (c / "backup" / "DONE").exists()
    assert not any(w.stor.iterdir())


def test_v1_backup_with_a_missing_file_fails_and_writes_no_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _backup_ready(w, c)
    (c / "root_added.tsv").write_text("path\tbytes\tsha256\nnot-there.bin\t1\tab\n")
    got = w.run(drv, "backup")
    assert got.returncode != 0
    assert not (c / "backup" / "DONE").exists()


def test_v1_backup_that_cannot_write_the_copy_fails_and_writes_no_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _backup_ready(w, c)
    w.stor.rmdir()
    got = w.run(drv, "backup")
    assert got.returncode != 0
    assert not (c / "backup" / "DONE").exists()


def test_v1_a_passed_backup_writes_its_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _backup_ready(w, c)
    got = w.run(drv, "backup")
    _ok(got, w, "v1-backup")
    assert (c / "backup" / "DONE").exists()
    assert (w.stor / "added_v1.tar").exists() and (c / "backup.tsv").exists()


def _run_ready(w: World, c: Path) -> None:
    _mark(c, "before/DONE", "rh/OK", "estimate/DONE")
    for step in V1_STEPS:
        d = w.root / "runs" / "stage=V1"
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{step}.json").write_text("{}")


@pytest.mark.parametrize("marker", ["before/DONE", "rh/OK", "estimate/DONE"])
def test_v1_run_refuses_without_the_markers_of_the_modes_before_it(tmp_path, marker):
    w, drv, c = _v1(tmp_path)
    _run_ready(w, c)
    (c / marker).unlink()
    got = w.run(drv, "run")
    assert got.returncode != 0
    assert not (c / "run" / "DONE").exists()
    assert "pod_run" not in w.calls()


def test_v1_run_whose_unit_failed_is_the_exit_status(tmp_path):
    w, drv, c = _v1(tmp_path)
    _run_ready(w, c)
    got = w.run(drv, "run", FAKE_POD_FAIL="c48c-v1-")
    assert got.returncode != 0
    assert not (c / "run" / "DONE").exists()


def test_v1_run_with_a_missing_run_record_fails(tmp_path):
    w, drv, c = _v1(tmp_path)
    _run_ready(w, c)
    (w.root / "runs" / "stage=V1" / "tapes.json").unlink()
    got = w.run(drv, "run")
    assert got.returncode != 0
    assert not (c / "run" / "DONE").exists()


def test_v1_run_whose_rates_check_crashed_at_the_end_fails(tmp_path):
    w, drv, c = _v1(tmp_path)
    _run_ready(w, c)
    got = w.run(drv, "run", FAKE_FAIL="c48c_v1_est.py rates")
    assert got.returncode != 0, "the rates helper crashed (exit 1, not the guard's 3) and `run` went on"
    assert not (c / "run" / "DONE").exists()


def test_v1_run_writes_its_marker(tmp_path):
    w, drv, c = _v1(tmp_path)
    _run_ready(w, c)
    _ok(w.run(drv, "run"), w, "v1-run")
    assert (c / "run" / "DONE").exists()


def test_v1_rehearse_whose_view_failed_is_the_exit_status(tmp_path):
    w, drv, c = _v1(tmp_path)
    got = w.run(drv, "rehearse", FAKE_FAIL="c26_view.py")
    assert got.returncode != 0
    assert not (c / "rh" / "OK").exists()


def test_v1_an_unknown_mode_exits_2(tmp_path):
    w, drv, c = _v1(tmp_path)
    assert w.run(drv, "nosuchmode").returncode == 2


# --------------------------------------------------------------------------------------------------------------
# release_train.sh run

TRAIN = BSP / "release_train.sh"


def _train(tmp_path: Path):
    w = World(tmp_path)
    drv = w.driver("release_train.sh")
    work = w.pod / "c46" / "train"
    # release_train.sh names its real root in REAL_ROOT; the copy's rewrite points it at the test's root
    return w, drv, work


def _go(w: World, drv: Path, **extra) -> subprocess.CompletedProcess:
    return w.run(drv, "run", "--no-live-cpu", ROOT=str(w.root), **{"TMP_MIN_GB": "0", **extra})


def test_train_a_good_run_writes_its_marker(tmp_path):
    w, drv, work = _train(tmp_path)
    got = _go(w, drv)
    _ok(got, w, "train-run")
    assert (work / "DONE").read_text().split()[0] == "r26.11.2"


def test_train_a_second_run_withdraws_the_marker_of_the_first_when_it_fails(tmp_path):
    w, drv, work = _train(tmp_path)
    _ok(_go(w, drv))
    assert (work / "DONE").exists()
    got = _go(w, drv, FAKE_FAIL=r"cli pack verify|pack verify")
    assert got.returncode != 0
    assert not (work / "DONE").exists(), "a waiter would read the earlier run's marker as this run's"


@pytest.mark.parametrize("pattern,why", [
    (r"train_plan.py.*--json", "the plan"),
    (r"train_est.py before", "the estimate before launch"),
    (r" release verify ", "a release does not verify"),
    (r"pack verify", "the pack does not verify"),
    (r"train_check.py", "the post-release check"),
])
def test_train_each_inner_command_failing_is_the_exit_status(tmp_path, pattern, why):
    w, drv, work = _train(tmp_path)
    got = _go(w, drv, FAKE_FAIL=pattern)
    assert got.returncode != 0, why
    assert not (work / "DONE").exists()


@pytest.mark.parametrize("unit", ["rel-1", "pack"])
def test_train_a_unit_that_exited_nonzero_is_the_exit_status(tmp_path, unit):
    w, drv, work = _train(tmp_path)
    got = _go(w, drv, FAKE_POD_FAIL=unit)
    assert got.returncode != 0
    assert not (work / "DONE").exists()


def test_train_a_failed_releases_listing_is_the_exit_status(tmp_path):
    w, drv, work = _train(tmp_path)
    got = _go(w, drv, FAKE_FAIL=r"cli releases ")
    assert got.returncode != 0
    assert not (work / "DONE").exists()


def test_train_a_failed_examples_check_is_the_exit_status(tmp_path):
    w, drv, work = _train(tmp_path)
    got = _go(w, drv, FAKE_FAIL="release_examples")
    assert got.returncode != 0
    assert not (work / "DONE").exists()


def test_train_a_check_that_wrote_no_list_is_a_failure(tmp_path):
    w, drv, work = _train(tmp_path)
    got = _go(w, drv, FAKE_NO_ADDED="1")
    assert got.returncode != 0
    assert not (work / "DONE").exists()


def test_train_a_noop_writes_no_marker_and_exits_0(tmp_path):
    w, drv, work = _train(tmp_path)
    (w.tmp / "plan.json").write_text(PLAN.replace('"nothing_new": false', '"nothing_new": true'))
    _ok(_go(w, drv))
    assert not (work / "DONE").exists()


# --------------------------------------------------------------------------------------------------------------
# rel_wait.sh: starts a release only on the markers of the driver's modes

REL_WAIT = ORCH / "rel_wait.sh"
MODES = ("before/DONE", "run/DONE", "after/DONE", "check/DONE", "backup/DONE")


def _wait_world(tmp_path: Path, markers=MODES, line=True, backup_tsv=True, suite="passed"):
    w = World(tmp_path)
    d = w.pod / "c47"
    d.mkdir(parents=True)
    for m in markers:
        _mark(d, m)
    if backup_tsv:
        (d / "backup.tsv").write_text("x\t1\tab\n")
    (d / "queue.out").write_text("c47_queue: done: rebefore, run, after, check and backup all passed\n" if line
                                 else "c47_queue: STOPPED at mode check (exit 1)\n")
    sd = w.pod / "ci-tmp" / "rel1120-detached"
    sd.mkdir(parents=True)
    for k in range(8):
        (sd / f"log.{k}").write_text(f"10 {suite} in 1s\n")
    tree = w.pod / "ci" / "rel1120"
    (tree / "scripts" / "stress" / "build_stage_pod").mkdir(parents=True)
    stub = tree / "scripts" / "stress" / "build_stage_pod" / "release_train.sh"
    stub.write_text('#!/bin/bash\necho "release_train $*" >> "$FAKE_LOG"\n[ "$1" = --dry-run ] && echo "Action: cut"\nexit 0\n')
    (w.pod / "c46" / "wait").mkdir(parents=True)
    return w, d


def _wait(w: World, **extra) -> subprocess.CompletedProcess:
    sh = w.bin / "rel_wait.sh"
    sh.write_text(REL_WAIT.read_text().replace("/tmp/if2453", str(w.pod)))
    return subprocess.run(["bash", str(sh)], env=w.env(**extra), capture_output=True, text=True, timeout=120)


def test_rel_wait_starts_the_release_when_every_marker_is_there(tmp_path):
    w, d = _wait_world(tmp_path)
    got = _wait(w)
    _ok(got, w, "rel-wait")
    assert "release_train run" in w.calls()


@pytest.mark.parametrize("missing", MODES + ("backup.tsv",))
def test_rel_wait_does_not_start_on_a_log_line_when_a_marker_is_missing(tmp_path, missing):
    w, d = _wait_world(tmp_path, backup_tsv=missing != "backup.tsv",
                       markers=[m for m in MODES if m != missing])
    got = _wait(w)
    assert got.returncode != 0, "the queue's log line was taken at its word"
    assert "release_train" not in w.calls()


def test_rel_wait_does_not_start_when_the_queue_failed(tmp_path):
    w, d = _wait_world(tmp_path, line=False, markers=MODES[:3])
    got = _wait(w)
    assert got.returncode != 0
    assert "release_train" not in w.calls()


def test_rel_wait_does_not_start_when_a_suite_shard_failed(tmp_path):
    w, d = _wait_world(tmp_path, suite="failed")
    got = _wait(w)
    assert got.returncode != 0
    assert "release_train" not in w.calls()


def test_rel_wait_takes_the_driver_directory_and_the_suite_as_arguments(tmp_path):
    w, d = _wait_world(tmp_path)
    other = w.pod / "c48c"
    d.rename(other)
    (w.pod / "ci-tmp" / "rel1120-detached").rename(w.pod / "ci-tmp" / "other9-detached")
    (w.pod / "ci" / "rel1120").rename(w.pod / "ci" / "other9")
    sh = w.bin / "rel_wait.sh"
    sh.write_text(REL_WAIT.read_text().replace("/tmp/if2453", str(w.pod)))
    got = subprocess.run(["bash", str(sh), str(other), "other9"], env=w.env(), capture_output=True, text=True,
                         timeout=120)
    _ok(got)
    assert "release_train run" in w.calls()

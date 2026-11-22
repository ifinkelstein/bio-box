"""Ring R1: card A-20 — post-manifest resume, verified-pack reuse, and host-profile pack promotion in release_train.sh."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-20", "A-15", "C-46")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"

PACK_ID = "02dec8312fff239cfed68b255e02f60e"

SYSTEMD_RUN = """#!/bin/sh
while [ $# -gt 0 ]; do
  case "$1" in --user|--scope|--quiet|--collect|--unit=*) shift ;; -p) shift 2 ;; *) break ;; esac
done
exec "$@"
"""
NICE = '#!/bin/sh\n[ "$1" = -n ] && shift 2\nexec "$@"\n'
NOOP = "#!/bin/sh\nexit 0\n"
INACTIVE = """#!/bin/sh
echo "systemctl $*" >> "$FAKE_LOG"
case "$*" in
  *is-active*) exit 3 ;;
  *show*-p*ActiveState*) echo inactive; exit 0 ;;
esac
exit 0
"""
SINGULARITY = """#!/bin/sh
echo "singularity $*" >> "$FAKE_LOG"
if [ -n "${FAKE_FAIL:-}" ] && echo "$*" | grep -Eq -- "$FAKE_FAIL"; then echo "fake failure: $*" >&2; exit 1; fi
case "$*" in
  *" release verify "*) echo "release R: verified, manifest hash abc123"; exit 0 ;;
  *"cli releases "*)
    prof="" ; prev=""
    for a in "$@"; do [ "$prev" = "--profile" ] && prof=$a; prev=$a; done
    echo "r26.11.2 listed profile=$prof"
    exit 0 ;;
  *train_plan.py*--json*) cat "$FAKE_PLAN"; exit 0 ;;
  *train_plan.py*) echo "Action: cut"; exit 0 ;;
  *train_est.py*) echo '{}'; exit 0 ;;
  *" pack verify "*) echo "pack verified"; exit 0 ;;
esac
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
    printf 'path\\tbytes\\tsha256\\nx\\t1\\tab\\n' > "$out/root_added.tsv"
    echo "# announcement" > "$out/announcement.md"
    echo "check ok" ;;
  *) echo '{}' ;;
esac
exit 0
"""
POD_RUN = f"""#!/bin/sh
log=""; unit=""; prof=""
while [ $# -gt 0 ]; do
  case "$1" in
    --log) log=$2; shift 2 ;;
    --profile) prof=$2; shift 2 ;;
    --*) shift 2 ;;
    *) unit=$1; shift; break ;;
  esac
done
echo "pod_run $unit $*" >> "$FAKE_LOG"
rel=""; prev=""
for a in "$@"; do [ "$prev" = "--release" ] && rel=$a; prev=$a; done
rc=0
if [ -n "${{FAKE_POD_FAIL:-}}" ] && echo "$unit" | grep -Eq -- "$FAKE_POD_FAIL"; then rc=1; fi
mkdir -p "$(dirname "$log")"
if [ "$rc" = 0 ] && [ "$1" = "pack" ] && [ "${{2:-}}" = "build" ]; then
  tmp_dir=$(awk -F': *' '/^tmp:/{{print $2; exit}}' "$prof")
  pdir="$tmp_dir/packs/release=$rel/{PACK_ID}"
  mkdir -p "$pdir"
  printf '{{"pack_id": "{PACK_ID}"}}\\n' > "$pdir/pack_manifest.json"
  printf '{{"ok": true}}\\n' > "$pdir/pack_verified.json"
  printf 'pack-payload\\n' > "$pdir/data.bin"
  echo "pack {PACK_ID} built" >> "$log"
fi
{{
  [ "$rc" = 0 ] && [ "$1" = "build-stage" ] && [ -n "$rel" ] && echo "build-stage: stage P1 and release $rel: all 9 steps done"
  echo "pod_run: $unit exited $rc"
}} >> "$log"
exit 0
"""
EX = """#!/bin/sh
echo "release_examples $*" >> "$FAKE_LOG"
[ -n "${FAKE_FAIL:-}" ] && echo "release_examples" | grep -Eq -- "$FAKE_FAIL" && exit 1
exit 0
"""

PLAN_R2 = {
    "release": "r26.11.2",
    "stages": ["S1", "P1"],
    "new_stages": ["P1"],
    "latest_release": "r26.11.1",
    "nothing_new": False,
    "earlier_manifest_releases": ["r26.11.1"],
    "last_stage": {
        "stage": "P1", "depth": "called", "genomes": "g.json", "objects": "o.json", "crispr": False,
    },
}

PLAN_AFTER_R2_NOOP = {
    "release": "r26.11.3",
    "stages": ["S1", "P1"],
    "new_stages": [],
    "latest_release": "r26.11.2",
    "nothing_new": True,
    "earlier_manifest_releases": ["r26.11.1", "r26.11.2"],
    "last_stage": {
        "stage": "P1", "depth": "called", "genomes": "g.json", "objects": "o.json", "crispr": False,
    },
}

PLAN_AFTER_R2_NEW_S2 = {
    "release": "r26.11.3",
    "stages": ["S1", "P1", "S2"],
    "new_stages": ["S2"],
    "latest_release": "r26.11.2",
    "nothing_new": False,
    "earlier_manifest_releases": ["r26.11.1", "r26.11.2"],
    "last_stage": {
        "stage": "S2", "depth": "called", "genomes": "g2.json", "objects": "o2.json", "crispr": False,
    },
}

FIRST_VERIFY = "release r26.11.1: verified (stages S1, 1 files, 1 artifacts, manifest hash abc123)\n"
FIRST_ESTIMATE = '{"first_launch": true}\n'


class TrainWorld:
    def __init__(self, tmp: Path):
        self.tmp = tmp
        self.pod = tmp / "if2453"
        self.fake = tmp / "fakebin"
        self.bin = tmp / "bin"
        self.tree = tmp / "tree"
        self.root = self.pod / "r10" / "root"
        self.work = self.pod / "c46" / "train"
        self.host_tmp = self.pod / "bio-box"
        self.log = tmp / "fake.log"
        for d in (self.fake, self.bin, self.tree / "scripts", self.tree / "profiles", self.root, self.host_tmp):
            d.mkdir(parents=True, exist_ok=True)
        for name, text in (("systemd-run", SYSTEMD_RUN), ("nice", NICE), ("sleep", NOOP),
                           ("systemctl", INACTIVE), ("singularity", SINGULARITY), ("python3", PYTHON3)):
            p = self.fake / name
            p.write_text(text)
            p.chmod(0o755)
        (self.tree / "scripts" / "pod_run.sh").write_text(POD_RUN)
        (self.tree / "scripts" / "release_examples.sh").write_text(EX)
        self.profile = self.tree / "profiles" / "wilkcomp02.yaml"
        self.profile.write_text(f"tmp: {self.host_tmp}\nmax_threads_day: 40\n")
        (self.root / "manifests").mkdir(parents=True, exist_ok=True)
        (self.root / "manifests" / "release=r26.11.1.json").write_text("{}\n")
        self.plan_file = self.tmp / "plan.json"
        self.set_plan(PLAN_R2)
        src = (BSP / "release_train.sh").read_text().replace("/tmp/if2453", str(self.pod))
        self.drv = self.bin / "release_train.sh"
        self.drv.write_text(src)
        for other in BSP.iterdir():
            link = self.bin / other.name
            if other.name != "release_train.sh" and not link.exists():
                link.symlink_to(other)
        self.log.write_text("")

    def set_plan(self, obj: dict) -> None:
        self.plan_file.write_text(json.dumps(obj))

    def publish_release(self, rel: str = "r26.11.2", n_records: int = 2) -> None:
        rdir = self.root / "runs" / f"stage={rel}"
        rdir.mkdir(parents=True, exist_ok=True)
        for k in range(n_records):
            (rdir / f"cluster.step{k}.json").write_text('{"cpu_min": 1.0}\n')
        (self.root / "release_def" / f"release={rel}").mkdir(parents=True, exist_ok=True)
        mdir = self.root / "manifests"
        mdir.mkdir(parents=True, exist_ok=True)
        (mdir / f"release={rel}.json").write_text("{}\n")

    def run(self, *args: str, **extra: str) -> subprocess.CompletedProcess:
        env = {
            **os.environ,
            "PATH": f"{self.fake}:{os.environ['PATH']}",
            "FAKE_LOG": str(self.log),
            "FAKE_PLAN": str(self.plan_file),
            "REAL_PYTHON": sys.executable,
            "TREE": str(self.tree),
            "ROOT": str(self.root),
            "WORK": str(self.work),
            "START_MIN_GB": "0",
            "TMP_MIN_GB": "0",
            **extra,
        }
        return subprocess.run(["bash", str(self.drv), *args], env=env, capture_output=True, text=True, timeout=120)

    def calls(self) -> str:
        return self.log.read_text()


def _files(top: Path) -> dict[Path, tuple[bytes, int]]:
    return {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in sorted(top.rglob("*")) if p.is_file()}


def _post_manifest_failed_examples(tmp_path: Path) -> TrainWorld:
    """First launch of r26.11.2 publishes manifest + verified pack, then fails in release_examples.sh before DONE."""
    w = TrainWorld(tmp_path)
    got = w.run("run", "--no-live-cpu", FAKE_FAIL="release_examples")
    assert got.returncode != 0, got.stdout + got.stderr
    assert not (w.work / "DONE").exists()
    assert (w.work / "tmp" / "packs" / "release=r26.11.2" / PACK_ID / "pack_verified.json").is_file()
    (w.work / "before" / "verify-r26.11.1.log").write_text(FIRST_VERIFY)
    (w.work / "estimate" / "before.json").write_text(FIRST_ESTIMATE)
    w.publish_release("r26.11.2", n_records=2)
    w.set_plan(PLAN_AFTER_R2_NOOP)
    w.log.write_text("")
    return w


# (a) Post-manifest resume when DONE is absent -----------------------------------------------------------------


def test_a20_dry_run_post_manifest_without_done_reports_resume_and_writes_nothing(tmp_path):
    w = _post_manifest_failed_examples(tmp_path)
    snap = _files(w.pod)
    got = w.run("--dry-run", "--no-live-cpu")
    assert got.returncode == 0, got.stdout + got.stderr
    assert "Action: resume r26.11.2 (2 run records published; unit 2)" in got.stdout, got.stdout
    assert "manifest published; finishing post-build steps" in got.stdout, got.stdout
    assert "the same as the first launch's" in got.stdout, got.stdout
    assert _files(w.pod) == snap


def test_a20_run_post_manifest_resume_reuses_verified_pack_and_writes_done(tmp_path):
    w = _post_manifest_failed_examples(tmp_path)
    listing = w.work / "before" / "listing.tsv"
    listing_snap = (listing.read_bytes(), listing.stat().st_mtime_ns)

    got = w.run("run", "--no-live-cpu")
    assert got.returncode == 0, got.stdout + got.stderr
    assert "Action: resume r26.11.2 (2 run records published; unit 2)" in got.stdout
    assert (listing.read_bytes(), listing.stat().st_mtime_ns) == listing_snap
    assert (w.work / "before" / "verify-r26.11.1.log").read_text() == FIRST_VERIFY
    assert (w.work / "estimate" / "before.json").read_text() == FIRST_ESTIMATE
    assert (w.work / "resume-2" / "verify-r26.11.1.log").is_file()

    calls = w.calls()
    assert "train_est.py before" not in calls
    assert "pod_run c46-train-pack" not in calls, "verified pack on disk must be reused instead of rebuilding"
    assert " pack verify " in calls
    assert "release_examples r26.11.2" in calls
    assert f"--pack-id {PACK_ID}" in calls
    assert (w.work / "DONE").read_text().split()[0] == "r26.11.2"


def test_a20_resume_rebuilds_pack_when_pack_verified_missing(tmp_path):
    w = _post_manifest_failed_examples(tmp_path)
    # Remove pack_verified.json from both work and host pack dirs so resume must run pack build.
    for base in (w.work / "tmp" / "packs", w.host_tmp / "packs"):
        pv = base / "release=r26.11.2" / PACK_ID / "pack_verified.json"
        if pv.exists():
            pv.unlink()
    got = w.run("run", "--no-live-cpu")
    assert got.returncode == 0, got.stdout + got.stderr
    assert "pod_run c46-train-pack-2 " in w.calls(), w.calls()
    assert (w.work / "DONE").read_text().split()[0] == "r26.11.2"


# (b) Post-manifest without DONE refused when a new complete stage appeared ------------------------------------


def test_a20_post_manifest_without_done_refuses_when_new_stage_appeared(tmp_path):
    w = _post_manifest_failed_examples(tmp_path)
    w.set_plan(PLAN_AFTER_R2_NEW_S2)
    snap = _files(w.work / "before") | _files(w.work / "estimate") | {
        w.work / "plan.json": (w.work / "plan.json").read_bytes(),
    }

    dry = w.run("--dry-run", "--no-live-cpu")
    assert dry.returncode != 0
    assert "Action: refused:" in dry.stdout, dry.stdout
    assert "finish r26.11.2 before cutting over new stages" in dry.stdout, dry.stdout

    got = w.run("run", "--no-live-cpu")
    assert got.returncode != 0
    assert "finish r26.11.2 before cutting over new stages" in got.stderr, got.stderr
    assert "pod_run " not in w.calls()
    assert not (w.work / "DONE").exists()
    assert (
        _files(w.work / "before") | _files(w.work / "estimate") | {
            w.work / "plan.json": (w.work / "plan.json").read_bytes(),
        }
        == snap
    )


# (c) When DONE already holds r26.11.2: noop keeps DONE; new stage cuts r26.11.3 as first ---------------------


def test_a20_done_present_zero_new_stages_noop_preserves_done(tmp_path):
    w = TrainWorld(tmp_path)
    got1 = w.run("run", "--no-live-cpu")
    assert got1.returncode == 0, got1.stdout + got1.stderr
    assert (w.work / "DONE").read_text().split()[0] == "r26.11.2"

    w.publish_release("r26.11.2", n_records=2)
    w.set_plan(PLAN_AFTER_R2_NOOP)
    w.log.write_text("")

    dry = w.run("--dry-run", "--no-live-cpu")
    assert dry.returncode == 0, dry.stdout + dry.stderr
    assert "Action: resume" not in dry.stdout

    got2 = w.run("run", "--no-live-cpu")
    assert got2.returncode == 0, got2.stdout + got2.stderr
    assert "0 new complete stages since r26.11.2; nothing to cut" in got2.stdout
    assert "pod_run " not in w.calls()
    assert (w.work / "DONE").read_text().split()[0] == "r26.11.2"


def test_a20_done_present_new_stage_cuts_next_release_as_first(tmp_path):
    w = TrainWorld(tmp_path)
    got1 = w.run("run", "--no-live-cpu")
    assert got1.returncode == 0, got1.stdout + got1.stderr
    assert (w.work / "DONE").read_text().split()[0] == "r26.11.2"

    w.publish_release("r26.11.2", n_records=2)
    w.set_plan(PLAN_AFTER_R2_NEW_S2)
    w.log.write_text("")

    got2 = w.run("run", "--no-live-cpu")
    assert got2.returncode == 0, got2.stdout + got2.stderr
    assert "Action: resume" not in got2.stdout
    assert "train_est.py before" in w.calls()
    assert "pod_run c46-train-rel-2 " in w.calls()
    assert (w.work / "DONE").read_text().split()[0] == "r26.11.3"


# (d) Host-profile shared pack directory promotion when ROOT == REAL_ROOT --------------------------------------


def test_a20_promotes_pack_to_host_profile_dir_when_root_is_real_root(tmp_path):
    w = TrainWorld(tmp_path)
    got = w.run("run", "--no-live-cpu", REAL_ROOT=str(w.root))
    assert got.returncode == 0, got.stdout + got.stderr

    dst_rel = w.host_tmp / "packs" / "release=r26.11.2"
    dst_pack = dst_rel / PACK_ID
    assert (dst_pack / "pack_verified.json").is_file(), "pack was not promoted to host profile pack dir"
    assert (dst_pack / "data.bin").read_text() == "pack-payload\n"
    assert not list(dst_rel.glob(".*.tmp.*")), "dot-temp directory left behind"
    assert (w.work / "pack" / "verify-realprof.log").is_file()
    assert f"pack verify --root {w.root} --release r26.11.2 --profile {w.profile}" in w.calls()
    assert f"profile={w.profile}" in (w.work / "after" / "releases.txt").read_text()


def test_a20_does_not_overwrite_existing_host_pack(tmp_path):
    w = TrainWorld(tmp_path)
    dst_pack = w.host_tmp / "packs" / "release=r26.11.2" / PACK_ID
    dst_pack.mkdir(parents=True)
    (dst_pack / "sentinel.txt").write_text("existing-pack\n")

    got = w.run("run", "--no-live-cpu", REAL_ROOT=str(w.root))
    assert got.returncode == 0, got.stdout + got.stderr
    assert (dst_pack / "sentinel.txt").read_text() == "existing-pack\n"
    assert not (dst_pack / "data.bin").exists(), "existing host pack directory was overwritten"
    assert (w.work / "pack" / "verify-realprof.log").is_file()
    assert f"profile={w.profile}" in (w.work / "after" / "releases.txt").read_text()


def test_a20_skips_host_pack_promotion_when_root_is_not_real_root(tmp_path):
    w = TrainWorld(tmp_path)
    got = w.run("run", "--no-live-cpu", REAL_ROOT=str(tmp_path / "different_real_root"))
    assert got.returncode == 0, got.stdout + got.stderr
    assert not (w.host_tmp / "packs" / "release=r26.11.2").exists()
    assert not (w.work / "pack" / "verify-realprof.log").exists()
    assert f"profile={w.work / 'profile.yaml'}" in (w.work / "after" / "releases.txt").read_text()

"""Ring R1: one stage driver for every stage (card A-12, 2026-11-21).

`scripts/stress/build_stage_pod/stage_driver.sh STAGE MODE` replaces the per-stage copies (`c47_driver.sh`,
`c48c_v1_driver.sh`, ...): it reads `stages/<STAGE>.conf` and runs the modes estimate, rehearse, before, run, after,
check and backup. `stage_queue.sh STAGE` runs them in order on their markers. Every mode's own status is the driver's
exit status, and a mode writes `<mode>/DONE` only on success (card A-11, STATUS "6m").

The driver runs here from a copy whose `/tmp/if2453` is rewritten to the test's temporary directory, with fake
`systemd-run`, `systemctl`, `singularity`, `python3` (the conf parser and inline code run on the real interpreter),
`pod_run.sh`, `sleep`, `nice` and `ionice` first on PATH, as tests/r1/test_driver_exit_status.py does. A pattern in
`FAKE_FAIL` makes every stub command whose arguments match it exit 1. Nothing of the pod's `/tmp/if2453` and nothing
under /stor is read or written. The success paths print `A12CALLS <name> <sha256>` of the command lines the stubs saw.

The pure parts of `stage_est.py` (the conf, the rate guard of DEC-43 and DEC-53, the accepted records, the estimate
before launch) and of `stage_check.py` (the ledger's within-stage rows, card C-47) run on the real code.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-12")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
ORCH = REPO / "scripts" / "orchestrator"
STAGES = BSP / "stages"
PFAM = "/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz"
MODES = ("estimate", "rehearse", "before", "run", "after", "check", "backup")


def _module(name: str):
    spec = importlib.util.spec_from_file_location(name, BSP / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # a dataclass looks its module up by name
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------------------------------------------
# stubs

SYSTEMD_RUN = """#!/bin/sh
while [ $# -gt 0 ]; do
  case "$1" in --user|--scope|--quiet|--collect|--unit=*) shift ;; -p) shift 2 ;; *) break ;; esac
done
case "$*" in *xargs*) [ -n "${FAKE_HASH_RC:-}" ] && exit "$FAKE_HASH_RC" ;; esac
exec "$@"
"""
NICE = '#!/bin/sh\n[ "$1" = -n ] && shift 2\nexec "$@"\n'
IONICE = '#!/bin/sh\n[ "$1" = -c ] && shift 2\nexec "$@"\n'
NOOP = "#!/bin/sh\nexit 0\n"
SLEEP = '#!/bin/sh\n[ -n "${FAKE_SLEEP:-}" ] && exec /bin/sleep "$FAKE_SLEEP"\nexit 0\n'
DF = """#!/bin/sh
case "$*" in *--output=avail*) echo "Avail"; echo "${FAKE_DF_GB:-500}"; exit 0 ;; esac
echo "Filesystem Size Used Avail Use% Mounted on"; echo "/dev/fake 1T 1G ${FAKE_DF_GB:-500}G 1% /"
"""
DATE = """#!/bin/sh
if [ "${1:-}" = "+%H" ] && [ -n "${FAKE_HOURS:-}" ] && [ -s "$FAKE_HOURS" ]; then
  head -n 1 "$FAKE_HOURS"; sed -i 1d "$FAKE_HOURS"; exit 0
fi
exec /bin/date "$@"
"""
SHA256SUM = """#!/bin/sh
if [ -n "${FAKE_SHA_EMPTY:-}" ]; then case "$*" in *.tar*) exit 0 ;; esac; fi
exec /usr/bin/sha256sum "$@"
"""
# A unit's state: FAKE_STATES (one answer a line: active, inactive or fail, which is a failed systemctl call), else
# FAKE_ACTIVE_FILE (active while its count is above 0), else inactive. A unit that turns inactive gets the exit line
# pod_run.sh's wrapper writes ("pod_run: <unit> exited <rc>"), which the pod_run stub left pending.
SYSTEMCTL = """#!/bin/sh
echo "systemctl $*" >> "$FAKE_LOG"
d=$(dirname "$FAKE_LOG")
finish() { f="$d/pod.$1"; [ -f "$f" ] || return 0; read -r log rc < "$f"; echo "pod_run: $1 exited $rc" >> "$log"; rm -f "$f"; }
state() {
  if [ -n "${FAKE_STOP_ON_ACTIVE:-}" ]; then touch "$FAKE_STOP_ON_ACTIVE"; fi
  if [ -n "${FAKE_STATES:-}" ] && [ -s "$FAKE_STATES" ]; then
    s=$(head -n 1 "$FAKE_STATES"); sed -i 1d "$FAKE_STATES"
  elif [ -n "${FAKE_ACTIVE_FILE:-}" ] && [ -s "$FAKE_ACTIVE_FILE" ] && [ "$(cat "$FAKE_ACTIVE_FILE")" -gt 0 ]; then
    echo $(( $(cat "$FAKE_ACTIVE_FILE") - 1 )) > "$FAKE_ACTIVE_FILE"; s=active
  else
    s=inactive
  fi
  echo "$s"
}
last=""
for a in "$@"; do last=$a; done
case "$*" in
  *is-active*)
    s=$(state)
    case "$s" in active) exit 0 ;; fail) exit 1 ;; esac
    finish "$last"; exit 3 ;;
  *"show -p ActiveState"*)
    s=$(state)
    [ "$s" = fail ] && exit 1
    [ "$s" = inactive ] && finish "$last"
    echo "$s"; exit 0 ;;
  *list-units*)
    units="${FAKE_UNITS:-}"
    [ -n "${FAKE_UNITS_FILE:-}" ] && [ -f "$FAKE_UNITS_FILE" ] && units="$units $(cat "$FAKE_UNITS_FILE")"
    for u in $units; do
      for a in "$@"; do
        case "$a" in --*|list-units) continue ;; esac
        case "$u" in $a) echo "$u loaded active running x"; break ;; esac
      done
    done
    exit 0 ;;
  *" stop "*)
    if [ -n "${FAKE_STOP_FAIL:-}" ] && [ -s "$FAKE_STOP_FAIL" ] && [ "$(cat "$FAKE_STOP_FAIL")" -gt 0 ]; then
      echo $(( $(cat "$FAKE_STOP_FAIL") - 1 )) > "$FAKE_STOP_FAIL"; exit 1
    fi
    exit 0 ;;
  *set-property*) exit "${FAKE_SETPROP_RC:-0}" ;;
esac
exit 0
"""
SINGULARITY = r"""#!/bin/sh
echo "singularity $*" >> "$FAKE_LOG"
if [ -n "${FAKE_FAIL:-}" ] && echo "$*" | grep -Eq -- "$FAKE_FAIL"; then echo "fake failure: $*" >&2; exit 1; fi
arg() { prev=""; for a in "$@"; do [ "$prev" = "$want" ] && { echo "$a"; return; }; prev=$a; done; }
case "$*" in
  *" release verify "*)
    want=--root; r=$(arg "$@"); want=--release; rel=$(arg "$@")
    h=$(sha256sum "$r/manifests/release=$rel.json" | cut -c1-32)
    echo "release $rel: verified, manifest hash $h"; exit 0 ;;
  *"cli releases "*) echo "releases listed"; exit 0 ;;
  *"stage_est.py before"*)
    want=--out-plan; o=$(arg "$@"); [ -n "$o" ] && echo '{}' > "$o"
    echo '{"estimate": "fake"}'; exit "${FAKE_EST_RC:-0}" ;;
  *"stage_est.py rates"*)
    [ -n "${FAKE_RATES_SLEEP:-}" ] && /bin/sleep "$FAKE_RATES_SLEEP"
    echo '{"over_20pct": []}'; exit "${FAKE_RATES_RC:-0}" ;;
  *"stage_est.py cap"*)
    if [ -n "${FAKE_QUOTAS:-}" ] && [ -s "$FAKE_QUOTAS" ]; then head -n 1 "$FAKE_QUOTAS"; sed -i 1d "$FAKE_QUOTAS"
    else echo "${FAKE_QUOTA:-54}"; fi
    exit 0 ;;
  *"stage_est.py complete"*) exit "${FAKE_COMPLETE_RC:-0}" ;;
  *"stage_est.py fixture"*)
    want=--out; o=$(arg "$@"); mkdir -p "$o/E" "$o/R"
    for s in E R; do echo '[]' > "$o/$s/genomes.json"; echo '[]' > "$o/$s/objects.json"; done
    printf '>t4\nACGT\n' > "$o/R/stage.fna"
    echo '{}'; exit 0 ;;
  *stage_check.py*)
    want=--out; o=$(arg "$@"); want=--root; r=$(arg "$@"); want=--before; b=$(arg "$@"); want=--after; af=$(arg "$@")
    mkdir -p "$o"
    if [ -z "${FAKE_NO_ADDED:-}" ]; then
      printf 'path\tbytes\tsha256\n' > "$o/root_added.tsv"
      cut -f1 "$b/listing.tsv" | LC_ALL=C sort > "$o/.b"
      cut -f1 "$af/listing.tsv" | LC_ALL=C sort > "$o/.a"
      LC_ALL=C comm -13 "$o/.b" "$o/.a" | while read -r p; do
        printf '%s\t%s\t%s\n' "$p" "$(stat -c %s "$r/$p")" "$(sha256sum "$r/$p" | cut -d' ' -f1)"
      done >> "$o/root_added.tsv"
    fi
    echo "RESULT: fake"; exit "${FAKE_CHECK_RC:-0}" ;;
esac
echo "RESULT: fake"
exit 0
"""
PYTHON3 = """#!/bin/sh
case "$*" in *"stage_est.py conf"*) exec "$REAL_PYTHON" "$@" ;; esac
case "${2:-}" in -c|-) exec "$REAL_PYTHON" "$@" ;; esac
echo "python3 $*" >> "$FAKE_LOG"
if [ -n "${FAKE_FAIL:-}" ] && echo "$*" | grep -Eq -- "$FAKE_FAIL"; then echo "fake failure: $*" >&2; exit 1; fi
case "$*" in *r50_watch.py*) echo '{"reasons": []}'; exit "${FAKE_GUARD_RC:-0}" ;; esac
echo '{}'
exit 0
"""
POD_RUN = """#!/bin/sh
log=""; unit=""
while [ $# -gt 0 ]; do
  case "$1" in
    --log) log=$2; shift 2 ;;
    --*) echo "pod_run-opt $1 $2" >> "$FAKE_LOG"; shift 2 ;;
    *) unit=$1; shift; break ;;
  esac
done
echo "pod_run $unit $*" >> "$FAKE_LOG"
root=""; st=""; rel=""; prev=""
for a in "$@"; do
  [ "$prev" = "--root" ] && root=$a; [ "$prev" = "--stage" ] && st=$a; [ "$prev" = "--release" ] && rel=$a; prev=$a
done
rc=0
if [ -n "${FAKE_POD_FAIL:-}" ] && echo "$unit" | grep -Eq -- "$FAKE_POD_FAIL"; then rc=1; fi
if [ "$rc" = 0 ] && [ -n "$root" ] && [ -n "$st" ]; then
  mkdir -p "$root/l0seq/x/stage=$st"
  # a relaunch publishes nothing: the stage's file is written once
  [ -e "$root/l0seq/x/stage=$st/part-0.parquet" ] || echo "$st" > "$root/l0seq/x/stage=$st/part-0.parquet"
  if [ -n "$rel" ]; then mkdir -p "$root/manifests"; echo "{}" > "$root/manifests/release=$rel.json"; fi
fi
mkdir -p "$(dirname "$log")"
echo "pod_run: $unit started" >> "$log"
# the exit line comes when the unit ends: the systemctl stub writes it when it first answers that the unit is inactive
echo "$log $rc" > "$(dirname "$FAKE_LOG")/pod.$unit"
exit 0
"""
PROFILE = """max_threads_day: 40
max_threads_night: 56
mem_scope: 200G
tmp: /x
durable: /x
bulk: /x
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: 300
cluster_cpu_min_max: 20000
small_file_floor_bytes: 67108864
require_mem_cap: true
"""
OTHER_STAGES = ("S1", "P1", "V1")


def conf_text(stage: str, inputs: Path, sums: dict[str, str], *, depth="annotated", library="pfam 38.2 gathering",
              crispr="yes", minlen=None, after=None, extra="") -> str:
    lines = [f"# test conf of stage {stage}", f"STAGE={stage}", "NOTE=planning/measurements/2026-11-20_stage_V1.md",
             f"INPUTS={inputs}"]
    lines += [f"SUM={h} {n}" for n, h in sums.items()]
    lines.append(f"DEPTH={depth}")
    if library:
        lines.append(f"LIBRARY={library}")
    lines.append(f"CRISPR={crispr}")
    if minlen is not None:
        lines.append(f"MIN_CONTIG_LEN={minlen}")
    if after:
        lines.append(f"AFTER={after}")
    return "\n".join(lines) + "\n" + extra


class World:
    """One test's temporary pod: stubs, the rewritten driver copies, a root with releases and stages, inputs, a conf."""

    def __init__(self, tmp: Path, stage: str = "ZZ1", **conf_kw):
        self.tmp = tmp
        self.pod = tmp / "if2453"
        self.fake = tmp / "fakebin"
        self.bin = tmp / "bin"
        self.tree = tmp / "tree"
        self.root = self.pod / "r10" / "root"
        self.log = tmp / "fake.log"
        self.stor = tmp / "stor"
        self.stage = stage
        self.c = self.pod / "stages" / stage.lower()
        for d in (self.fake, self.bin / "stages", self.tree / "scripts", self.tree / "profiles", self.root, self.stor):
            d.mkdir(parents=True, exist_ok=True)
        for name, text in (("systemd-run", SYSTEMD_RUN), ("nice", NICE), ("ionice", IONICE), ("sleep", SLEEP),
                           ("top", NOOP), ("uptime", NOOP), ("systemctl", SYSTEMCTL), ("singularity", SINGULARITY),
                           ("python3", PYTHON3), ("df", DF), ("sha256sum", SHA256SUM), ("date", DATE)):
            (self.fake / name).write_text(text)
            (self.fake / name).chmod(0o755)
        (self.tree / "scripts" / "pod_run.sh").write_text(POD_RUN)
        (self.tree / "profiles" / "wilkcomp02.yaml").write_text(PROFILE)
        # the notes and decision rows a conf names are checked against the tree (card A-12 (8))
        (self.tree / "planning" / "measurements").mkdir(parents=True)
        (self.tree / "planning" / "measurements" / "2026-11-20_stage_V1.md").write_text("a note\n")
        (self.tree / "planning" / "DECISIONS.md").write_text("| DEC-56 | 2026-11-20 | accept C5's ingest |\n")
        # the root: three published releases (manifests), one release_def without a manifest, three other stages
        (self.root / "manifests").mkdir()
        for r in ("rA", "rB", "r26.11.1"):
            (self.root / "manifests" / f"release={r}.json").write_text(json.dumps({"release": r}))
            (self.root / "release_def" / f"release={r}").mkdir(parents=True)
        (self.root / "release_def" / "release=rOld").mkdir(parents=True)
        for s in OTHER_STAGES:
            d = self.root / "l0seq" / "x" / f"stage={s}"
            d.mkdir(parents=True)
            (d / "part-0.parquet").write_bytes(s.encode())
            (self.root / "runs" / f"stage={s}").mkdir(parents=True)
            (self.root / "runs" / f"stage={s}" / "ingest.json").write_text("{}")
        (self.root / "l0seq" / "store_params.json").write_text("{}")
        # the stage's inputs
        self.inputs = self.pod / "s2" / "inputs" / stage
        self.inputs.mkdir(parents=True)
        (self.inputs / "genomes.json").write_text("[]")
        (self.inputs / "objects.json").write_text("[]")
        (self.inputs / "x.fna").write_text(">a\nACGT\n")
        sums = {n: hashlib.sha256((self.inputs / n).read_bytes()).hexdigest()
                for n in ("genomes.json", "objects.json", "x.fna")}
        self.write_conf(stage, conf_text(stage, self.inputs, sums, **conf_kw))
        self.log.write_text("")

    def write_conf(self, stage: str, text: str) -> None:
        (self.bin / "stages" / f"{stage}.conf").write_text(text)

    def script(self, name: str) -> Path:
        """A copy of scripts/.../<name> with the pod's paths rewritten to this test's directory."""
        text = (BSP / name).read_text().replace("/tmp/if2453", str(self.pod))
        dst = self.bin / name
        if dst.is_symlink():  # a link an earlier copy made: writing through it would rewrite the tree's own script
            dst.unlink()
        dst.write_text(text)
        for other in BSP.iterdir():  # scripts the driver names beside itself ($BIN/x.py) stay reachable
            link = self.bin / other.name
            if other.is_file() and other.name != name and not link.exists():
                link.symlink_to(other)
        return dst

    def env(self, **extra: str) -> dict:
        return {**os.environ, "PATH": f"{self.fake}:{os.environ['PATH']}", "FAKE_LOG": str(self.log),
                "REAL_PYTHON": sys.executable, "TREE": str(self.tree), "BACKUP_DIR": str(self.stor), **extra}

    def drive(self, mode: str, *more: str, stage: str | None = None, **extra: str) -> subprocess.CompletedProcess:
        drv = self.script("stage_driver.sh")
        return subprocess.run(["bash", str(drv), self.stage if stage is None else stage, mode, *more],
                              env=self.env(**extra),
                              capture_output=True, text=True, timeout=300)

    def conf_sha(self, stage: str | None = None) -> str:
        return hashlib.sha256((self.bin / "stages" / f"{stage or self.stage}.conf").read_bytes()).hexdigest()

    def mark(self, *modes: str, c: Path | None = None) -> None:
        """Markers as the driver writes them: the conf's sha256 (card A-12 (4))."""
        for m in modes:
            d = (c or self.c) / m
            d.mkdir(parents=True, exist_ok=True)
            (d / "DONE").write_text(self.conf_sha() + "\n")

    def done(self, mode: str) -> bool:
        return (self.c / mode / "DONE").exists()

    def calls(self) -> str:
        return self.log.read_text().replace(str(self.tmp), "<T>")

    def record(self, name: str) -> None:
        print("A12CALLS", name, hashlib.sha256(self.calls().encode()).hexdigest())


def _ok(proc: subprocess.CompletedProcess, w: World | None = None, name: str = "") -> None:
    if w is not None and name:
        w.record(name)
    assert proc.returncode == 0, proc.stdout + proc.stderr + (w.c / "driver.log").read_text() if w else ""


def _driver_log(w: World) -> str:
    p = w.c / "driver.log"
    return p.read_text() if p.exists() else ""


# --------------------------------------------------------------------------------------------------------------
# the conf (stage_est.read_conf)

def _conf_file(tmp_path: Path, text: str, stage: str = "ZZ1") -> Path:
    p = tmp_path / f"{stage}.conf"
    p.write_text(text)
    return p


GOOD = {"genomes.json": "a" * 64, "objects.json": "b" * 64}


def test_a_conf_is_read_with_every_value(tmp_path):
    est = _module("stage_est")
    text = conf_text("ZZ1", Path("/tmp/if2453/s2/inputs/ZZ1"), GOOD, minlen=5000, after="V2,G1",
                     extra="JUDGE=genecall 236 planning/measurements/2026-11-20_virus_gv1_inputs.md\n"
                           "ACCEPTED=ingest 45.803 38.057 DEC-56\n")
    c = est.read_conf(_conf_file(tmp_path, text), "ZZ1")
    assert c.stage == "ZZ1" and c.depth == "annotated" and c.min_contig_len == 5000 and c.crispr
    assert c.library == ("pfam", "38.2", "gathering") and c.after == ("V2", "G1")
    assert c.judge == {"genecall": (236.0, "planning/measurements/2026-11-20_virus_gv1_inputs.md")}
    assert c.accepted == {("ingest", 45.803, 38.057): "DEC-56"}
    assert dict((n, h) for h, n in c.sums) == GOOD


@pytest.mark.parametrize("bad,why", [
    ("ACCEPTED=ingest 45.803 38.057\n", "an accepted record without its decision row"),
    ("ACCEPTED=ingest 45.803 38.057 DEC56\n", "a decision row not named DEC-<n>"),
    ("ACCEPTED=ingest 45.8 38.057 DEC-56\n", "an accepted value not at the precision of the record (3 decimals)"),
    ("THREADS=56\n", "an unknown key (a thread count belongs to the profile)"),
    ("INPUTS=/tmp/x\n", "a key given twice"),
    ("JUDGE=genecall many planning/x.md\n", "a figure that is not a number"),
    ("JUDGE=tapes 10 planning/x.md\n", "a judged step the driver does not group"),
    ("AFTER=ZZ1\n", "a stage after itself"),
    ("MIN_CONTIG_LEN=5kb\n", "a length that is not an integer"),
    ("SUM=abc genomes.json\n", "a sum that is not 64 hex digits"),
    ("SUM=" + "c" * 64 + " ../x.fna\n", "a summed file outside the inputs directory"),
])
def test_a_conf_with_a_bad_line_is_refused(tmp_path, bad, why):
    est = _module("stage_est")
    p = _conf_file(tmp_path, conf_text("ZZ1", Path("/tmp/x"), GOOD) + bad)
    with pytest.raises(est.ConfError):
        est.read_conf(p, "ZZ1")


@pytest.mark.parametrize("kw,why", [
    ({"depth": "annotated", "library": None}, "annotated needs a library"),
    ({"depth": "called"}, "called takes no library"),
    ({"depth": "stored", "library": None}, "a stored stage is not one this driver builds"),
    ({"library": "pfam 37.0 gathering"}, "a library whose HMM file the driver does not know"),
    ({"library": "pfam 38.2 loose"}, "a cutoff that build-stage does not take"),
    ({"crispr": "maybe"}, "CRISPR is yes or no"),
])
def test_a_conf_with_a_bad_combination_is_refused(tmp_path, kw, why):
    est = _module("stage_est")
    p = _conf_file(tmp_path, conf_text("ZZ1", Path("/tmp/x"), GOOD, **kw))
    with pytest.raises(est.ConfError):
        est.read_conf(p, "ZZ1")


def test_a_conf_must_sum_genomes_and_objects(tmp_path):
    est = _module("stage_est")
    p = _conf_file(tmp_path, conf_text("ZZ1", Path("/tmp/x"), {"genomes.json": "a" * 64}))
    with pytest.raises(est.ConfError):
        est.read_conf(p, "ZZ1")


def test_a_conf_of_another_stage_is_refused(tmp_path):
    est = _module("stage_est")
    p = _conf_file(tmp_path, conf_text("ZZ2", Path("/tmp/x"), GOOD), stage="ZZ1")
    with pytest.raises(est.ConfError):
        est.read_conf(p, "ZZ1")


def test_a_value_with_a_shell_metacharacter_is_refused(tmp_path):
    est = _module("stage_est")
    p = _conf_file(tmp_path, conf_text("ZZ1", Path("/tmp/x;rm -rf /"), GOOD))
    with pytest.raises(est.ConfError):
        est.read_conf(p, "ZZ1")


# the queued stages' confs in the repository: each value from its note

def _notes_text(conf) -> str:
    out = []
    for n in conf.notes:
        p = REPO / n
        out.append(p.read_text())
        d = p.with_suffix("")
        if d.is_dir():
            out += [f.read_text(errors="replace") for f in sorted(d.rglob("*")) if f.is_file() and f.stat().st_size < 2e6]
    return "\n".join(out)


def _repo_confs():
    return sorted(STAGES.glob("*.conf"))


def test_the_repository_has_a_conf_for_each_stage_whose_note_gives_every_value():
    assert {p.stem for p in _repo_confs()} == {"V2", "GV1"}


@pytest.mark.parametrize("path", _repo_confs(), ids=lambda p: p.stem)
def test_each_repository_conf_is_valid_and_its_sums_and_figures_are_its_notes(path):
    est = _module("stage_est")
    c = est.read_conf(path, path.stem, tree=REPO)  # every note it names is a file of the tree
    text = _notes_text(c)
    for h, name in c.sums:
        assert h in text, f"{path.name}: sum of {name} is not in its note"
    for step, (fig, src) in c.judge.items():
        assert f"{fig:g}" in (REPO / src).read_text(), f"{path.name}: {step} {fig} is not in {src}"
    for st, (n, src) in c.expect.items():
        note = (REPO / src).read_text()
        assert f"{n:,}" in note or f" {n} " in note or f"| {n} |" in note, f"{path.name}: held_by {st} {n} not in {src}"
    assert not c.accepted, "a new stage carries no accepted record"
    assert c.inputs == Path(f"/tmp/if2453/s2/inputs/{c.stage}")


def test_gv1_comes_after_v2_judges_its_gene_calling_against_the_notes_figure_and_expects_v2s_contigs():
    est = _module("stage_est")
    c = est.read_conf(STAGES / "GV1.conf", "GV1", tree=REPO)
    assert c.after == ("V2",)
    assert c.judge["genecall"][0] == 236.0
    assert c.expect["V2"][0] == 17664


@pytest.mark.parametrize("extra,why", [
    ("ACCEPTED=ingest 45.803 38.057 DEC-999\n", "a decision row the tree's DECISIONS.md does not have"),
    ("NOTE=planning/measurements/no_such_note.md\n", "a note that is not a file of the tree"),
    ("JUDGE=genecall 10 planning/measurements/no_such_note.md\n", "a JUDGE figure's note that is not a file"),
    ("EXPECT=held_by V2 5 planning/measurements/no_such_note.md\n", "an EXPECT note that is not a file"),
])
def test_a_conf_naming_what_the_tree_does_not_hold_is_refused(tmp_path, extra, why):
    est = _module("stage_est")
    tree = tmp_path / "tree"
    (tree / "planning" / "measurements").mkdir(parents=True)
    (tree / "planning" / "measurements" / "2026-11-20_stage_V1.md").write_text("x")
    (tree / "planning" / "DECISIONS.md").write_text("| DEC-56 | ... |\n")
    p = _conf_file(tmp_path, conf_text("ZZ1", Path("/tmp/x"), GOOD) + extra)
    with pytest.raises(est.ConfError):
        est.read_conf(p, "ZZ1", tree=tree)
    good = _conf_file(tmp_path, conf_text("ZZ1", Path("/tmp/x"), GOOD) + "ACCEPTED=ingest 45.803 38.057 DEC-56\n")
    assert est.read_conf(good, "ZZ1", tree=tree).accepted


def test_a_record_stopped_against_a_judge_figure_can_be_accepted_by_its_group_name(tmp_path):
    est = _module("stage_est")
    text = conf_text("ZZ1", Path("/tmp/x"), GOOD, extra="JUDGE=genecall 25 planning/measurements/2026-11-20_stage_V1.md\n"
                     "ACCEPTED=genecall.note 40.100 25.000 DEC-56\n")
    c = est.read_conf(_conf_file(tmp_path, text), "ZZ1")
    doc = est.judge_items([_item("genecall.note", 40.1, 25.0)], c.accepted)
    assert doc["over_20pct"] == [] and doc["accepted_over_20pct"] == {"genecall.note": "DEC-56"}


@pytest.mark.parametrize("bad", ["EXPECT=held_by ZZ1 5 planning/measurements/a.md\n",
                                 "EXPECT=held_by V2 many planning/measurements/a.md\n",
                                 "EXPECT=owned V2 5 planning/measurements/a.md\n"])
def test_a_bad_expect_line_is_refused(tmp_path, bad):
    est = _module("stage_est")
    with pytest.raises(est.ConfError):
        est.read_conf(_conf_file(tmp_path, conf_text("ZZ1", Path("/tmp/x"), GOOD) + bad), "ZZ1")


# --------------------------------------------------------------------------------------------------------------
# the rate guard (stage_est.judge_items): DEC-43, DEC-53, accepted records only by their decision row

def _item(name, cpu, est, complete=True):
    return {"name": name, "cpu_min": cpu, "estimate": est, "basis": "test", "complete": complete, "records": 1}


def test_a_step_estimated_under_24_cpu_minutes_is_never_judged():
    est = _module("stage_est")
    doc = est.judge_items([_item("tapes", 50.0, 10.0)], {})
    assert doc["over_20pct"] == [] and not doc["rows"][0]["judged"]


def test_a_judged_step_more_than_20_percent_over_stops():
    est = _module("stage_est")
    doc = est.judge_items([_item("ingest", 30.25, 25.0)], {})
    assert doc["over_20pct"] == ["ingest"] and est.rates_exit(doc) == 3


def test_a_judged_step_more_than_20_percent_under_is_recorded_and_does_not_stop():
    est = _module("stage_est")
    doc = est.judge_items([_item("crispr-x", 18.9, 25.0)], {})
    assert doc["under_20pct"] == ["crispr-x"] and doc["over_20pct"] == [] and est.rates_exit(doc) == 0


def test_a_group_still_running_is_judged_over_on_its_partial_sum_and_never_under():
    est = _module("stage_est")
    over = est.judge_items([_item("genecall:note", 300.0, 236.0, complete=False)], {})
    under = est.judge_items([_item("genecall:note", 10.0, 236.0, complete=False)], {})
    assert over["over_20pct"] == ["genecall:note"]
    assert under["under_20pct"] == [] and under["over_20pct"] == []


def test_an_accepted_record_is_reported_and_not_stopped_only_at_its_exact_values():
    est = _module("stage_est")
    acc = {("ingest", 45.803, 38.057): "DEC-56"}
    doc = est.judge_items([_item("ingest", 45.803, 38.057)], acc)
    assert doc["over_20pct"] == [] and doc["accepted_over_20pct"] == {"ingest": "DEC-56"}
    other = est.judge_items([_item("ingest", 45.9, 38.057)], acc)
    assert other["over_20pct"] == ["ingest"], "an accepted record covers its one record, at its values"


def test_the_group_cpu_of_ingest_adds_every_piece():
    est = _module("stage_est")
    recs = {"ingest.pass1-00000": {"cpu_min": 100.0}, "ingest.plan": {"cpu_min": 1.0},
            "ingest.batch-0000": {"cpu_min": 50.0}, "ingest": {"cpu_min": 2.0}, "inputs": {"cpu_min": 0.5}}
    names = est.group_names(recs, "ingest")
    assert sorted(names) == ["ingest", "ingest.batch-0000", "ingest.pass1-00000", "ingest.plan"]
    assert sum(recs[n]["cpu_min"] for n in names) == 153.0


# --------------------------------------------------------------------------------------------------------------
# the estimate before launch, on real code (stage_est.py before): no root is written

def _seq(rng: random.Random, n: int) -> bytes:
    return bytes(rng.choice(b"ACGT") for _ in range(n))


def test_the_estimate_before_launch_counts_flanks_as_the_crispr_step_does_and_drops_short_records(tmp_path):
    from biobox import ids
    from biobox import stage as bs
    from biobox.l1 import crispr, genecall

    rng = random.Random(12)
    circ, lin, short = _seq(rng, 150_000), _seq(rng, 30_000), _seq(rng, 2_000)
    inp = tmp_path / "inputs"
    inp.mkdir()
    with gzip.open(inp / "x.fna.gz", "wb") as fh:
        for h, s in (("c", circ), ("l", lin), ("s", short)):
            fh.write(b">" + h.encode() + b"\n" + s + b"\n")
    cu = {k: ids.contig_uid(v) for k, v in (("c", circ), ("l", lin), ("s", short))}
    genomes = [{"contig_uids": [cu["c"].hex()], "genetic_code": 11, "kind": "plasmid", "circular": [cu["c"].hex()],
                "assembly_quality": "contig"},
               {"contig_uids": [cu["l"].hex()], "genetic_code": 11, "kind": "virus", "circular": [],
                "assembly_quality": "contig"}]
    (inp / "genomes.json").write_text(json.dumps(genomes))
    (inp / "objects.json").write_text(json.dumps([{"path": str(inp / "x.fna.gz"), "source_name": "t",
                                                   "source_release": "1", "license": "t", "redistributable": True}]))
    sums = {n: hashlib.sha256((inp / n).read_bytes()).hexdigest() for n in ("genomes.json", "objects.json")}
    conf = _conf_file(tmp_path, conf_text("ZZ1", inp, sums, depth="called", library=None, minlen=5000))
    root = tmp_path / "root"
    (root / "l0seq" / "contig_index").mkdir(parents=True)
    prof = tmp_path / "profile.yaml"
    prof.write_text(PROFILE.replace("tmp: /x", f"tmp: {tmp_path}"))
    got = subprocess.run([sys.executable, str(BSP / "stage_est.py"), "before", "--conf", str(conf), "--stage", "ZZ1",
                          "--root", str(root), "--profile", str(prof), "--inputs", str(inp), "--tmp",
                          str(tmp_path / "t")], capture_output=True, text=True, timeout=600,
                         env={**os.environ, "PYTHONPATH": str(REPO / "src")})
    assert got.returncode == 0, got.stderr
    doc = json.loads(got.stdout)
    # an independent count: the code's formula for the CRISPR step's flanks (stage.crispr_flank_bases, card B-53)
    assert doc["flank_bases"] == 2 * min(crispr.CIRC_FLANK, 150_000)
    assert doc["bases_stored"] == 180_000 and doc["short_records"] == 1 and doc["short_bases"] == 2_000
    assert doc["held_bases"] == 0
    want = (genecall.PROCESS_S * 1e6 + genecall.genome_cost_us("plasmid", [(150_000, True)])
            + genecall.genome_cost_us("virus", [(30_000, False)])) / 60e6
    assert doc["genecall"]["estimate_cpu_min"] == pytest.approx(want, abs=1e-3)
    assert doc["crispr_default_cpu_min"] == pytest.approx(bs.crispr_cpu_min(180_000, 200_000), abs=1e-3)
    assert "crispr_near_cpu_min" not in doc  # no library: the near-Cas producer does not run


# --------------------------------------------------------------------------------------------------------------
# stage_check: the gene-call ledger's rows, as card C-47 corrected them

def test_the_check_counts_a_within_stage_row_with_a_genome_suffix_as_within_and_compares_whole_stage_names():
    import duckdb

    chk = _module("stage_check")
    con = duckdb.connect()
    con.execute("CREATE TABLE led (reason_code VARCHAR)")
    reasons = ["ok", "owned_by:PL1", "owned_by:C5/0011aa", "owned_by:C5", "owned_by:C50", "duplicate_contig:x"]
    con.execute("INSERT INTO led VALUES " + ", ".join(f"('{r}')" for r in reasons))
    got = [r[0] for r in con.execute(f"SELECT {chk.skip_kind_sql('reason_code', 'C5')} FROM led").fetchall()]
    assert got == [None, "earlier", "within", "within", "earlier", "within"]


def test_the_check_reads_tables_through_duckdb_and_not_whole_into_python():
    src = (BSP / "stage_check.py").read_text()
    assert "to_pylist()" not in src and "to_pydict()" not in src and "pq.read_table" not in src


# --------------------------------------------------------------------------------------------------------------
# the driver: refusals before anything runs

def test_an_unknown_mode_exits_2(tmp_path):
    w = World(tmp_path)
    assert w.drive("nosuchmode").returncode == 2


@pytest.mark.parametrize("stage", ["zz1", "Z;Z", "../ZZ1", ""])
def test_a_bad_stage_name_exits_2(tmp_path, stage):
    w = World(tmp_path)
    assert w.drive("before", stage=stage).returncode == 2


def test_a_stage_without_a_conf_is_refused(tmp_path):
    w = World(tmp_path)
    got = w.drive("estimate", stage="ZZ9")
    assert got.returncode != 0
    assert "singularity" not in w.calls()


@pytest.mark.parametrize("mode", MODES)
def test_every_mode_stops_on_a_stop_file_and_writes_no_marker(tmp_path, mode):
    w = World(tmp_path)
    w.mark(*MODES)
    (w.c / mode / "DONE").unlink()
    (w.root / "STOP").write_text("test")
    got = w.drive(mode)
    assert got.returncode != 0
    assert not w.done(mode)
    assert "pod_run" not in w.calls()


# --------------------------------------------------------------------------------------------------------------
# estimate

def test_estimate_with_inputs_that_differ_from_the_conf_fails_and_writes_no_marker(tmp_path):
    w = World(tmp_path)
    (w.inputs / "x.fna").write_text(">a\nACGA\n")
    got = w.drive("estimate")
    assert got.returncode != 0 and not w.done("estimate")
    assert "stage_est.py before" not in w.calls()


@pytest.mark.parametrize("rc,why", [("1", "the estimate crashed"), ("5", "a contig is owned by an incomplete stage"),
                                    ("7", "build-stage would refuse the stage")])
def test_estimate_whose_helper_fails_writes_no_marker(tmp_path, rc, why):
    w = World(tmp_path)
    got = w.drive("estimate", FAKE_EST_RC=rc)
    assert got.returncode != 0, why
    assert not w.done("estimate")


def test_estimate_writes_its_marker_and_withdraws_the_markers_after_it(tmp_path):
    w = World(tmp_path)
    w.mark("before", "run", "after", "check")
    _ok(w.drive("estimate"), w, "estimate")
    assert w.done("estimate")
    assert not any(w.done(m) for m in ("before", "run", "after", "check"))
    assert (w.c / "estimate" / "before.json").read_text().strip()


# --------------------------------------------------------------------------------------------------------------
# before: the releases come from the root, never from a list in the script

def _verified(w: World) -> list[str]:
    return sorted(re.findall(r" release verify --root \S+ --release (\S+)", w.calls()))


def test_before_verifies_every_release_the_root_has_a_manifest_for(tmp_path):
    w = World(tmp_path)
    _ok(w.drive("before"), w, "before")
    assert _verified(w) == ["r26.11.1", "rA", "rB"]
    assert "rOld" in _driver_log(w), "a release_def without a manifest is named in the log"
    assert w.done("before")


def test_before_verifies_a_release_added_to_the_root_without_any_change_to_the_scripts(tmp_path):
    w = World(tmp_path)
    (w.root / "manifests" / "release=r26.11.2.json").write_text("{}")
    _ok(w.drive("before"))
    assert "r26.11.2" in _verified(w)


def test_before_leaves_the_stages_own_files_out_of_its_listing_and_hashes_every_other_stage(tmp_path):
    w = World(tmp_path)
    own = w.root / "l0seq" / "x" / f"stage={w.stage}"
    own.mkdir(parents=True)
    (own / "part-0.parquet").write_text("published by an earlier launch")
    _ok(w.drive("before"))
    listing = (w.c / "before" / "listing.tsv").read_text()
    assert f"stage={w.stage}/" not in listing and "stage=S1/" in listing
    hashed = (w.c / "before" / "earlier.sha256").read_text()
    for s in OTHER_STAGES:
        assert f"stage={s}/part-0.parquet" in hashed
    assert f"stage={w.stage}" not in hashed


def test_before_a_failed_verify_is_the_exit_status(tmp_path):
    w = World(tmp_path)
    got = w.drive("before", FAKE_FAIL=r" release verify ")
    assert got.returncode != 0 and not w.done("before")


def test_before_a_failed_hash_fails_the_mode(tmp_path):
    w = World(tmp_path)
    got = w.drive("before", FAKE_HASH_RC="1")
    assert got.returncode != 0 and not w.done("before")


def test_before_with_no_published_release_is_refused(tmp_path):
    w = World(tmp_path)
    for p in (w.root / "manifests").iterdir():
        p.unlink()
    got = w.drive("before")
    assert got.returncode != 0 and not w.done("before")


@pytest.mark.parametrize("unit", ["c46-train-rel-1", "c47-c5-3", "zz2-1", "v2-4"])
def test_before_refuses_while_a_build_unit_in_the_root_is_active(tmp_path, unit):
    w = World(tmp_path)
    w.write_conf("ZZ2", conf_text("ZZ2", w.inputs, {"genomes.json": "a" * 64, "objects.json": "b" * 64}))
    w.write_conf("V2", conf_text("V2", w.inputs, {"genomes.json": "a" * 64, "objects.json": "b" * 64}))
    got = w.drive("before", FAKE_UNITS=unit)
    assert got.returncode != 0 and not w.done("before")
    assert "release verify" not in w.calls()


def test_before_is_not_held_by_the_stages_own_queue(tmp_path):
    w = World(tmp_path)
    _ok(w.drive("before", FAKE_UNITS=f"{w.stage.lower()}-queue"))


# --------------------------------------------------------------------------------------------------------------
# run

def _run_ready(w: World) -> None:
    w.mark("estimate", "rehearse", "before")


def _build_line(w: World) -> str:
    lines = [ln for ln in w.calls().splitlines() if ln.startswith("pod_run ")]
    assert len(lines) == 1, lines
    return lines[0]


def test_run_launches_build_stage_with_the_confs_values_and_writes_its_marker(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    got = w.drive("run")
    _ok(got, w, "run")
    assert w.done("run")
    c = "<T>/if2453/stages/zz1"
    want = (f"pod_run zz1-1 build-stage --profile {c}/profile.yaml --root <T>/if2453/r10/root --stage ZZ1 "
            f"--depth annotated --genomes <T>/if2453/s2/inputs/ZZ1/genomes.json "
            f"--objects <T>/if2453/s2/inputs/ZZ1/objects.json --library pfam 38.2 {PFAM} gathering --crispr "
            "--threads 6 --scan-cpus 4")
    assert _build_line(w) == want
    assert "pod_run-opt --threads 54" in w.calls(), "the unit's quota: the cap of the hour less the helpers' 2"
    assert "pod_run-opt --mem 200G" in w.calls()
    assert "stage_est.py complete" in w.calls()
    # the rate guard reads the FASTA bytes the estimate measured, never decompressing every object again
    assert re.search(r"stage_est.py rates .*--before \S+/stages/zz1/estimate/before.json", w.calls())


def test_run_of_a_called_stage_with_a_least_length_and_no_library(tmp_path):
    w = World(tmp_path, depth="called", library=None, minlen=5000)
    _run_ready(w)
    _ok(w.drive("run"))
    line = _build_line(w)
    assert "--depth called --min-contig-len 5000 --genomes" in line
    assert "--library" not in line and line.endswith("--crispr --threads 6 --scan-cpus 4")


@pytest.mark.parametrize("marker", ["estimate", "rehearse", "before"])
def test_run_refuses_without_the_markers_of_the_modes_before_it(tmp_path, marker):
    w = World(tmp_path)
    _run_ready(w)
    (w.c / marker / "DONE").unlink()
    got = w.drive("run")
    assert got.returncode != 0 and not w.done("run")
    assert "pod_run" not in w.calls()


def test_run_refuses_after_an_earlier_failure_until_a_person_moves_it_aside(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    (w.c / "run").mkdir(parents=True)
    (w.c / "run" / "FAILED").write_text("guard: over 280 CPU-minutes")
    got = w.drive("run")
    assert got.returncode != 0 and "pod_run" not in w.calls()
    assert (w.c / "run" / "FAILED").read_text() == "guard: over 280 CPU-minutes", "a stop stays a stop"


def test_run_refuses_while_a_build_unit_is_active(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    got = w.drive("run", FAKE_UNITS="c46-train-rel-1")
    assert got.returncode != 0 and "pod_run" not in w.calls()


def test_run_refuses_under_its_free_space_floor(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    got = w.drive("run", START_MIN_GB="999999")
    assert got.returncode != 0 and "pod_run" not in w.calls()


@pytest.mark.parametrize("env", [{"T": "57"}, {"CPUS": "41"}, {"Q": "64"}, {"CPU_MAX": "300"}])
def test_run_refuses_a_thread_count_above_the_profiles_or_a_process_guard_above_280(tmp_path, env):
    w = World(tmp_path)
    _run_ready(w)
    got = w.drive("run", **env)
    assert got.returncode != 0 and "pod_run" not in w.calls()


def test_run_waits_for_the_stages_it_comes_after(tmp_path):
    w = World(tmp_path, after="V2")
    _run_ready(w)
    got = w.drive("run")
    assert got.returncode != 0 and "pod_run" not in w.calls()
    (w.pod / "stages" / "v2" / "backup").mkdir(parents=True)
    (w.pod / "stages" / "v2" / "backup" / "DONE").touch()
    _ok(w.drive("run"))


def test_run_whose_unit_failed_is_the_exit_status(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    got = w.drive("run", FAKE_POD_FAIL="zz1-")
    assert got.returncode != 0 and not w.done("run")
    assert (w.c / "run" / "FAILED").exists()


def test_run_of_an_incomplete_stage_fails(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    got = w.drive("run", FAKE_COMPLETE_RC="1")
    assert got.returncode != 0 and not w.done("run")


def test_run_whose_rates_helper_crashed_at_the_end_fails(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    got = w.drive("run", FAKE_RATES_RC="1")
    assert got.returncode != 0 and not w.done("run")


@pytest.mark.parametrize("env,why", [({"FAKE_GUARD_RC": "3"}, "the 280 CPU-minute or free-space guard"),
                                     ({"FAKE_GUARD_RC": "1"}, "the guard itself crashed: a blind guard is a stop"),
                                     ({"FAKE_RATES_RC": "3"}, "DEC-53: a judged step more than 20% over")])
def test_run_stops_its_unit_on_a_guard(tmp_path, env, why):
    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("3")
    got = w.drive("run", FAKE_ACTIVE_FILE=str(active), **env)
    assert got.returncode != 0, why
    assert not w.done("run")
    assert "systemctl --user stop zz1-1" in w.calls()
    assert (w.c / "run" / "FAILED").exists()


def test_run_ends_at_a_stop_that_appears_and_writes_no_marker(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("2")
    got = w.drive("run", FAKE_ACTIVE_FILE=str(active), FAKE_STOP_ON_ACTIVE=str(w.root / "STOP"))
    assert got.returncode != 0 and not w.done("run")
    assert "systemctl --user stop" not in w.calls(), "build-stage itself ends at the step boundary"


def test_run_numbers_its_units_in_lower_case(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    (w.c).mkdir(parents=True, exist_ok=True)
    (w.c / "units").write_text("zz1-1\n")
    _ok(w.drive("run"))
    assert _build_line(w).startswith("pod_run zz1-2 build-stage")


# --------------------------------------------------------------------------------------------------------------
# after

def test_after_refuses_without_the_run_marker(tmp_path):
    w = World(tmp_path)
    got = w.drive("after")
    assert got.returncode != 0 and not w.done("after")


def _after_ready(w: World) -> None:
    w.mark("estimate", "rehearse", "before")
    _ok(w.drive("before"))
    w.mark("run")
    own = w.root / "l0seq" / "x" / f"stage={w.stage}"
    own.mkdir(parents=True)
    (own / "part-0.parquet").write_bytes(b"new")


def test_after_writes_its_marker_with_every_release_verified_and_its_hash_unchanged(tmp_path):
    w = World(tmp_path)
    _after_ready(w)
    _ok(w.drive("after"), w, "after")
    assert w.done("after")
    assert _verified(w).count("rA") == 2


@pytest.mark.parametrize("change,why", [
    ("hash", "a manifest changed: its hash differs from before"),
    ("new", "a release appeared between before and after"),
])
def test_after_fails_when_the_releases_differ_from_before(tmp_path, change, why):
    w = World(tmp_path)
    _after_ready(w)
    if change == "hash":
        (w.root / "manifests" / "release=rA.json").write_text('{"release": "rA", "changed": 1}')
    else:
        (w.root / "manifests" / "release=rC.json").write_text("{}")
    got = w.drive("after")
    assert got.returncode != 0, why
    assert not w.done("after")


@pytest.mark.parametrize("env", [{"FAKE_FAIL": r" release verify "}, {"FAKE_FAIL": r"cli releases "},
                                 {"FAKE_HASH_RC": "1"}])
def test_after_inner_failures_are_the_exit_status(tmp_path, env):
    w = World(tmp_path)
    _after_ready(w)
    got = w.drive("after", **env)
    assert got.returncode != 0 and not w.done("after")


# --------------------------------------------------------------------------------------------------------------
# check

def _check_ready(w: World) -> None:
    _after_ready(w)
    _ok(w.drive("after"))


def test_check_refuses_without_the_after_marker(tmp_path):
    w = World(tmp_path)
    got = w.drive("check")
    assert got.returncode != 0 and not w.done("check")


def test_check_writes_its_marker_and_the_list_of_files_added(tmp_path):
    w = World(tmp_path)
    _check_ready(w)
    _ok(w.drive("check"), w, "check")
    assert w.done("check")
    added = (w.c / "root_added.tsv").read_text()
    assert f"l0seq/x/stage={w.stage}/part-0.parquet\t3\t" in added
    assert added == (w.c / "check" / "1" / "root_added.tsv").read_text()


def test_check_passes_the_conf_to_the_check_script(tmp_path):
    w = World(tmp_path)
    _check_ready(w)
    _ok(w.drive("check"))
    line = next(ln for ln in w.calls().splitlines() if "stage_check.py" in ln)
    assert f"--conf <T>/bin/stages/{w.stage}.conf" in line and f"--stage {w.stage}" in line
    assert "--before-est <T>/if2453/stages/zz1/estimate/before.json" in line


@pytest.mark.parametrize("env", [{"FAKE_CHECK_RC": "1"}, {"FAKE_NO_ADDED": "1"}])
def test_a_failed_check_is_the_exit_status(tmp_path, env):
    w = World(tmp_path)
    _check_ready(w)
    got = w.drive("check", **env)
    assert got.returncode != 0 and not w.done("check")


def test_a_check_after_a_backup_must_list_the_files_the_backup_holds(tmp_path):
    w = World(tmp_path)
    _check_ready(w)
    _ok(w.drive("check"))
    (w.c / "backup.tsv").write_text("x\t1\tab\t1\n")
    extra = w.root / "l0seq" / "x" / f"stage={w.stage}" / "part-1.parquet"
    extra.write_bytes(b"later")
    _ok(w.drive("after"))
    got = w.drive("check")
    assert got.returncode != 0 and not w.done("check")


# --------------------------------------------------------------------------------------------------------------
# backup: a temporary name, read back, mv -n, never over an existing file

def _backup_ready(w: World) -> None:
    _check_ready(w)
    _ok(w.drive("check"))


def test_backup_refuses_without_the_check_marker(tmp_path):
    w = World(tmp_path)
    (w.c).mkdir(parents=True, exist_ok=True)
    (w.c / "root_added.tsv").write_text("path\tbytes\tsha256\n")
    got = w.drive("backup", TAR_FLOOR_BYTES="0")
    assert got.returncode != 0 and not w.done("backup")
    assert not any(w.stor.iterdir())


def test_backup_writes_the_tar_reads_it_back_and_writes_its_marker(tmp_path):
    w = World(tmp_path)
    _backup_ready(w)
    _ok(w.drive("backup", TAR_FLOOR_BYTES="0"), w, "backup")
    assert w.done("backup")
    tar = w.stor / f"added_{w.stage.lower()}.tar"
    rec = (w.c / "backup.tsv").read_text().split("\t")
    assert rec[0] == str(tar) and rec[2] == hashlib.sha256(tar.read_bytes()).hexdigest()
    assert not list(w.stor.glob(".*")), "no temporary file is left"


def test_backup_under_the_small_file_floor_fails(tmp_path):
    w = World(tmp_path)
    _backup_ready(w)
    got = w.drive("backup")
    assert got.returncode != 0 and not w.done("backup")
    assert not any(w.stor.iterdir())


def test_backup_never_overwrites_a_tar_that_exists(tmp_path):
    w = World(tmp_path)
    _backup_ready(w)
    tar = w.stor / f"added_{w.stage.lower()}.tar"
    tar.write_bytes(b"an earlier copy")
    got = w.drive("backup", TAR_FLOOR_BYTES="0")
    assert got.returncode != 0 and not w.done("backup")
    assert tar.read_bytes() == b"an earlier copy"


def test_backup_with_a_missing_file_fails(tmp_path):
    w = World(tmp_path)
    _backup_ready(w)
    with open(w.c / "root_added.tsv", "a") as fh:
        fh.write("not/there.parquet\t1\tab\n")
    got = w.drive("backup", TAR_FLOOR_BYTES="0")
    assert got.returncode != 0 and not w.done("backup")


def test_backup_that_cannot_write_its_copy_fails(tmp_path):
    w = World(tmp_path)
    _backup_ready(w)
    w.stor.chmod(0o500)
    try:
        got = w.drive("backup", TAR_FLOOR_BYTES="0")
    finally:
        w.stor.chmod(0o700)
    assert got.returncode != 0 and not w.done("backup")


def test_a_floor_of_zero_is_refused_for_a_backup_under_stor(tmp_path):
    w = World(tmp_path)
    _backup_ready(w)
    got = w.drive("backup", TAR_FLOOR_BYTES="0", BACKUP_DIR="/stor/work/Wilke/if2453/all-the-data/root_backups")
    assert got.returncode != 0 and not w.done("backup")
    assert "dd " not in _driver_log(w)


def test_a_second_backup_after_success_writes_nothing(tmp_path):
    w = World(tmp_path)
    _backup_ready(w)
    _ok(w.drive("backup", TAR_FLOOR_BYTES="0"))
    before = (w.c / "backup.tsv").read_text()
    _ok(w.drive("backup", TAR_FLOOR_BYTES="0"))
    assert (w.c / "backup.tsv").read_text() == before and w.done("backup")


# --------------------------------------------------------------------------------------------------------------
# rehearse: the fixture, in scratch, through the driver's own modes

def test_rehearse_runs_every_mode_on_the_fixture_in_scratch_and_writes_its_marker(tmp_path):
    w = World(tmp_path)
    _ok(w.drive("rehearse"), w, "rehearse")
    assert w.done("rehearse")
    rh = w.c / "rh"
    for m in ("estimate", "before", "run", "after", "check", "backup"):
        assert (rh / "drv" / m / "DONE").exists(), m
    calls = w.calls()
    assert not re.search(r"--root <T>/if2453/r10/root(\s|$)", calls), "nothing of the rehearsal writes the data root"
    assert (rh / "backups" / "added_zz1r.tar").exists()
    assert not any(w.stor.iterdir()), "the rehearsal's backup goes to its scratch directory"
    for unit in re.findall(r"^pod_run (\S+)", calls, flags=re.M):
        assert unit.startswith("zz1-rh-"), unit
    assert "pod_run-opt --threads 4" in calls and "pod_run-opt --threads 56" not in calls


@pytest.mark.parametrize("env", [{"FAKE_FAIL": "stage_est.py fixture"}, {"FAKE_POD_FAIL": "zz1-rh-1"},
                                 {"FAKE_CHECK_RC": "1"}, {"FAKE_POD_FAIL": "zz1-rh-e"}])
def test_rehearse_whose_inner_step_fails_writes_no_marker(tmp_path, env):
    w = World(tmp_path)
    got = w.drive("rehearse", **env)
    assert got.returncode != 0 and not w.done("rehearse")


def test_the_rehearsal_drivers_root_cannot_be_the_data_root(tmp_path):
    w = World(tmp_path)
    w.mark("estimate", "before")
    got = w.drive("run", STAGE_REHEARSAL="1")
    assert got.returncode != 0 and "pod_run" not in w.calls()


# --------------------------------------------------------------------------------------------------------------
# info: what the queue reads

def test_info_names_the_stages_to_wait_for_and_the_active_build_units(tmp_path):
    w = World(tmp_path, after="V2")
    got = w.drive("info", FAKE_UNITS="c46-train-rel-1 zz1-queue")
    _ok(got)
    lines = dict(ln.split("=", 1) for ln in got.stdout.splitlines() if "=" in ln)
    assert lines["AFTER"] == "V2" and lines["BUSY"] == "c46-train-rel-1"


# --------------------------------------------------------------------------------------------------------------
# the queue

DRIVER_STUB = """#!/bin/bash
C=$STAGES_HOME/$(echo "$1" | tr 'A-Z' 'a-z')
case "$2" in
  info) echo "AFTER=${STUB_AFTER:-}"; echo "BUSY=${STUB_BUSY:-}"; exit 0 ;;
esac
mkdir -p "$C"
echo "$2" >> "$C/modes"
[ -n "${STUB_NO_MARKER:-}" ] && [ "$2" = "$STUB_NO_MARKER" ] && exit 0
[ -n "${STUB_FAIL:-}" ] && [ "$2" = "$STUB_FAIL" ] && exit 1
mkdir -p "$C/$2"
touch "$C/$2/DONE"
exit 0
"""


def _queue(tmp_path: Path, **extra: str) -> subprocess.CompletedProcess:
    home = tmp_path / "stages"
    home.mkdir(exist_ok=True)
    b = tmp_path / "qbin"
    b.mkdir(exist_ok=True)
    (b / "stage_driver.sh").write_text(DRIVER_STUB)
    q = b / "stage_queue.sh"
    q.write_text((BSP / "stage_queue.sh").read_text().replace("/tmp/if2453", str(tmp_path / "pod")))
    (tmp_path / "pod" / "r10" / "root").mkdir(parents=True, exist_ok=True)
    fake = tmp_path / "fakebin"
    fake.mkdir(exist_ok=True)
    for name, text in (("systemctl", SYSTEMCTL), ("sleep", NOOP)):
        (fake / name).write_text(text)
        (fake / name).chmod(0o755)
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(tmp_path / "fake.log"),
           "STAGES_HOME": str(home), "POLL": "0", "MAX_WAIT_H": "0", **extra}
    return subprocess.run(["bash", str(q), "ZZ1"], env=env, capture_output=True, text=True, timeout=120)


def _modes(tmp_path: Path) -> list[str]:
    p = tmp_path / "stages" / "zz1" / "modes"
    return p.read_text().split() if p.exists() else []


def test_the_queue_runs_every_mode_in_order_when_each_writes_its_marker(tmp_path):
    got = _queue(tmp_path)
    assert got.returncode == 0, got.stdout + got.stderr
    assert _modes(tmp_path) == list(MODES)
    assert "all passed" in got.stdout


def test_the_queue_stops_when_a_mode_exits_0_without_its_marker(tmp_path):
    got = _queue(tmp_path, STUB_NO_MARKER="before")
    assert got.returncode != 0
    assert _modes(tmp_path) == ["estimate", "rehearse", "before"]
    assert "all passed" not in got.stdout


def test_the_queue_stops_at_a_failed_mode(tmp_path):
    got = _queue(tmp_path, STUB_FAIL="check")
    assert got.returncode != 0
    assert _modes(tmp_path)[-1] == "check" and "backup" not in _modes(tmp_path)


def test_the_queue_starts_nothing_on_a_stop_file(tmp_path):
    (tmp_path / "pod" / "r10" / "root").mkdir(parents=True)
    (tmp_path / "pod" / "r10" / "root" / "STOP").touch()
    got = _queue(tmp_path)
    assert got.returncode != 0 and _modes(tmp_path) == []


def test_the_queue_waits_for_the_stage_it_comes_after_and_gives_up_at_its_limit(tmp_path):
    got = _queue(tmp_path, STUB_AFTER="V2")
    assert got.returncode != 0 and _modes(tmp_path) == []


def test_the_queue_starts_once_the_stage_it_comes_after_is_backed_up(tmp_path):
    (tmp_path / "stages" / "v2" / "backup").mkdir(parents=True)
    (tmp_path / "stages" / "v2" / "backup" / "DONE").touch()
    got = _queue(tmp_path, STUB_AFTER="V2")
    assert got.returncode == 0, got.stdout
    assert _modes(tmp_path) == list(MODES)


def test_the_queue_stops_when_the_stage_it_comes_after_failed(tmp_path):
    (tmp_path / "stages" / "v2" / "check").mkdir(parents=True)
    (tmp_path / "stages" / "v2" / "check" / "FAILED").write_text("x")
    got = _queue(tmp_path, STUB_AFTER="V2", MAX_WAIT_H="48")
    assert got.returncode != 0 and _modes(tmp_path) == []


def test_the_queue_waits_while_a_build_unit_is_active(tmp_path):
    got = _queue(tmp_path, STUB_BUSY="c46-train-rel-1")
    assert got.returncode != 0 and _modes(tmp_path) == []


# --------------------------------------------------------------------------------------------------------------
# rel_wait.sh takes the unit pattern from its argument

def test_rel_wait_holds_a_release_while_a_unit_of_the_stage_it_waits_for_is_active(tmp_path):
    pod = tmp_path / "if2453"
    d = pod / "stages" / "v2"
    for m in ("before", "run", "after", "check", "backup"):
        (d / m).mkdir(parents=True)
        (d / m / "DONE").touch()
    (d / "backup.tsv").write_text("x\t1\tab\n")
    (pod / "c46" / "wait").mkdir(parents=True)
    (pod / "r10" / "root").mkdir(parents=True)
    confs = pod / "ci" / "s9" / "scripts" / "stress" / "build_stage_pod" / "stages"
    confs.mkdir(parents=True)
    (confs / "V2.conf").write_text("STAGE=V2\n")
    fake = tmp_path / "fakebin"
    fake.mkdir()
    for name, text in (("systemctl", SYSTEMCTL), ("sleep", NOOP)):
        (fake / name).write_text(text)
        (fake / name).chmod(0o755)
    sh = tmp_path / "rel_wait.sh"
    sh.write_text((ORCH / "rel_wait.sh").read_text().replace("/tmp/if2453", str(pod)))
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(tmp_path / "fake.log"),
           "FAKE_UNITS": "v2-7"}
    got = subprocess.run(["bash", str(sh), str(d), "s9"], env=env, capture_output=True, text=True, timeout=120)
    assert got.returncode == 6, got.stdout
    assert "v2-7" in got.stdout


# --------------------------------------------------------------------------------------------------------------
# the acceptance comparison: c47_driver.sh and stage_driver.sh launch the same build-stage command for C5

def test_c47_and_the_stage_driver_launch_the_same_build_stage_command_for_c5(tmp_path):
    w = World(tmp_path, stage="C5")
    c = w.pod / "c47"
    inputs = w.pod / "s2" / "inputs"
    # C5's values (c47_driver.sh: depth annotated, Pfam 38.2, CRISPR, DEC-56's accepted record); the sums are those of
    # this test's input files, since `run` checks them again (card A-12 (4)) and they enter no command line
    sums = {n: hashlib.sha256((inputs / "C5" / n).read_bytes()).hexdigest() for n in ("genomes.json", "objects.json")}
    w.write_conf("C5", conf_text("C5", inputs / "C5", sums, extra="ACCEPTED=ingest 45.803 38.057 DEC-56\n"))
    # c47_driver.sh run, with its own overrides and its markers
    w.mark("before", c=c)
    (c / "rh").mkdir(parents=True, exist_ok=True)
    (c / "rh" / "OK").touch()
    (c / "estimate").mkdir(parents=True, exist_ok=True)
    (c / "estimate" / "before.json").write_text("{}")
    runs = w.root / "runs" / "stage=C5"
    runs.mkdir(parents=True, exist_ok=True)
    nr3 = "minced-0.4.2-w8-nr3-r23x47-s26x50"
    near = "minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2"
    for step in ("inputs", "ingest", "genecall", f"crispr-{nr3}", "igr", "tapes", "domains-pfam-38.2",
                 "domains-pfam-38.2.plan", f"crispr-{near}"):
        (runs / f"{step}.json").write_text("{}")
    drv = w.script("c47_driver.sh")
    old = subprocess.run(["bash", str(drv), "run"], capture_output=True, text=True, timeout=300,
                         env=w.env(C47_DIR=str(c), C47_ROOT=str(w.root), C47_INPUTS=str(inputs)))
    assert old.returncode == 0, old.stdout + old.stderr
    old_line = next(ln for ln in w.calls().splitlines() if ln.startswith("pod_run "))
    w.log.write_text("")
    # stage_driver.sh run, in the same directory, so the profile path is the same
    w.mark("estimate", "rehearse", c=c)
    new = w.drive("run", STAGE_DIR=str(c), STAGE_ROOT=str(w.root))
    assert new.returncode == 0, new.stdout + new.stderr + (c / "driver.log").read_text()
    new_line = next(ln for ln in w.calls().splitlines() if ln.startswith("pod_run "))
    print("A12CALLS c47-run", old_line)
    print("A12CALLS stage-run", new_line)
    old_unit, old_cmd = old_line.split(" ", 2)[1:]
    new_unit, new_cmd = new_line.split(" ", 2)[1:]
    assert old_unit.startswith("c47-c5-") and new_unit.startswith("c5-")
    assert new_cmd == old_cmd
    for name in ("stage_driver.sh", "c47_driver.sh"):  # the copies never write through to the tree's scripts
        assert "/tmp/if2453/r10/root" in (BSP / name).read_text()
        assert str(w.pod) not in (BSP / name).read_text()


# ==============================================================================================================
# Card A-12's amendment after the independent review (2026-11-21): one block per point

# (1) `run` leaves its loop only on the unit's exit line or a real inactive answer; a trap stops the unit

def _lines(path: Path, *rows: str) -> str:
    path.write_text("".join(r + "\n" for r in rows))
    return str(path)


def test_1_a_failed_systemctl_answer_does_not_end_the_watch(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    states = _lines(tmp_path / "states", "fail", "active", "fail", "active", "inactive")
    _ok(w.drive("run", FAKE_STATES=states), w, "run-a-failed-answer")
    assert w.done("run")
    assert w.calls().count("show -p ActiveState --value zz1-1") == 5
    assert "no answer from systemd about zz1-1" in _driver_log(w)


def test_1_a_stop_that_fails_is_asked_again_and_the_driver_waits_for_the_unit_to_go(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("20")
    stopfail = tmp_path / "stopfail"
    stopfail.write_text("1")
    got = w.drive("run", FAKE_ACTIVE_FILE=str(active), FAKE_STOP_FAIL=str(stopfail), FAKE_GUARD_RC="3")
    assert got.returncode != 0 and (w.c / "run" / "FAILED").exists()
    assert w.calls().count("systemctl --user stop zz1-1") == 2, "the failed stop is asked again"
    assert active.read_text().strip() == "0", "the driver said 'stopped' only once the unit was gone"


def test_1_the_driver_killed_while_its_unit_runs_stops_the_unit(tmp_path):
    import signal
    import time

    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("300")
    drv = w.script("stage_driver.sh")
    p = subprocess.Popen(["bash", str(drv), w.stage, "run"], env=w.env(FAKE_ACTIVE_FILE=str(active), FAKE_SLEEP="0.1"),
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    t0 = time.monotonic()
    while "show -p ActiveState --value zz1-1" not in w.calls() and time.monotonic() - t0 < 60:
        time.sleep(0.2)
    p.send_signal(signal.SIGTERM)
    p.communicate(timeout=120)
    assert p.returncode != 0
    assert "systemctl --user stop zz1-1" in w.calls(), "the trap stops the unit the driver leaves"
    assert "the driver exits (status 143) while zz1-1 runs" in _driver_log(w)


# (2) the release train and the release waiter know the stage driver's units

def test_2_the_release_train_does_not_start_while_a_stage_unit_builds(tmp_path):
    w = World(tmp_path)
    train = w.script("release_train.sh")
    got = subprocess.run(["bash", str(train), "run", "--no-live-cpu"], capture_output=True, text=True, timeout=120,
                         env=w.env(ROOT=str(w.root), FAKE_UNITS="zz1-3"))
    assert got.returncode != 0
    assert "active build unit(s): zz1-3" in got.stderr
    other = subprocess.run(["bash", str(train), "run", "--no-live-cpu"], capture_output=True, text=True, timeout=120,
                           env=w.env(ROOT=str(w.root), FAKE_UNITS="zz1-queue"))
    assert "active build unit(s)" not in other.stderr


def test_2_rel_wait_looks_again_just_before_the_run(tmp_path):
    pod = tmp_path / "if2453"
    d = pod / "stages" / "v2"
    for m in ("before", "run", "after", "check", "backup"):
        (d / m).mkdir(parents=True)
        (d / m / "DONE").touch()
    (d / "backup.tsv").write_text("x\t1\tab\n")
    (pod / "c46" / "wait").mkdir(parents=True)
    (pod / "r10" / "root").mkdir(parents=True)
    tree = pod / "ci" / "s9"
    (tree / "scripts" / "stress" / "build_stage_pod" / "stages").mkdir(parents=True)
    (tree / "scripts" / "stress" / "build_stage_pod" / "stages" / "GV1.conf").write_text("STAGE=GV1\n")
    units = tmp_path / "units"
    units.write_text("")
    (tree / "scripts" / "stress" / "build_stage_pod" / "release_train.sh").write_text(
        '#!/bin/bash\necho "release_train $*" >> "$FAKE_LOG"\n'
        f'[ "$1" = --dry-run ] && {{ echo gv1-2 > {units}; echo "Action: cut"; }}\nexit 0\n')
    suite = pod / "ci-tmp" / "s9-detached"
    suite.mkdir(parents=True)
    for k in range(8):
        (suite / f"log.{k}").write_text("10 passed in 1s\n")
    fake = tmp_path / "fakebin"
    fake.mkdir()
    for name, text in (("systemctl", SYSTEMCTL), ("sleep", NOOP)):
        (fake / name).write_text(text)
        (fake / name).chmod(0o755)
    sh = tmp_path / "rel_wait.sh"
    sh.write_text((ORCH / "rel_wait.sh").read_text().replace("/tmp/if2453", str(pod)))
    log = tmp_path / "fake.log"
    log.write_text("")
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log), "FAKE_UNITS_FILE": str(units)}
    got = subprocess.run(["bash", str(sh), str(d), "s9"], env=env, capture_output=True, text=True, timeout=120)
    assert got.returncode == 6, got.stdout
    assert "gv1-2" in got.stdout
    assert "release_train run" not in log.read_text(), "a unit of stage GV1 began after the dry run: no run"


# (3), (6), (7) the check, the rate guard and the completeness check on a constructed root (no tool needed)

def _u(i: int) -> bytes:
    return bytes([i]) * 16


def _write(root: Path, rel: str, rows: list[dict]) -> int:
    import pyarrow as pa
    import pyarrow.parquet as pq

    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), p)
    return p.stat().st_size


def _dom(p: bytes, acc: str, score: float) -> dict:
    return {"protein_uid": p, "model_acc": acc, "model_name": acc, "model_len": 100, "seq_bitscore": score,
            "seq_evalue_zref": 1e-10, "dom_bitscore": score, "dom_ievalue_zref": 1e-10, "env_start": 1, "env_end": 50,
            "ali_start": 1, "ali_end": 50, "hmm_from": 1, "hmm_to": 50, "cutoff": "gathering", "z": 1000,
            "dom_z": 1000}


def build_root(tmp: Path, mut=None, post=None, **conf_kw) -> dict:
    """A stage ZQ as build-stage leaves it, by hand: contigs A and B stored (B listed by two genomes: one row within
    the stage, DEC-27), H held by stage OA (DEC-26); genes on A and B; a Pfam plan whose copy from OA is the row of the
    protein on H (cards B-47b, C-38h); both CRISPR producers; a published gene-call plan; the two snapshots."""
    from biobox.l1 import crispr, genecall

    gcv = genecall.GENE_CALL_VERSION
    p1, p2 = crispr.producer(), crispr.near_cas_producer("38.2")
    A, B, H = _u(1), _u(2), _u(3)
    g1, g2, g3 = _u(11), _u(12), _u(13)
    pr1, pr2, pr3 = _u(21), _u(22), _u(23)

    def gi(c, g):
        return {"contig_uid": c, "genome_uid": g, "kind": "virus", "genetic_code": 11, "is_circular": False,
                "assembly_quality": "complete"}

    def ci(c, n, dv=True):
        return {"contig_uid": c, "contig_len": n, "is_default_variant": dv}

    def led(c, why):
        return {"object_key": c, "reason_code": why}

    def cl(n_in, kept, why="ok"):
        return {"n_in": n_in, "n_kept": kept, "n_rejected": 0, "reason_code": why}

    S = {
        "T": {
            "gi": ["l0seq/genome_input/stage=ZQ/part-0.parquet", [gi(A, g1), gi(B, g1), gi(B, g2), gi(H, g3)]],
            "ci": ["l0seq/contig_index/stage=ZQ/part-0.parquet", [ci(A, 2000), ci(B, 1500)]],
            "ingled": ["ledger/stage=ZQ/step=ingest/part-0.parquet",
                       [{"reason_code": "ok", "n_bases_rejected": 0}, {"reason_code": "below:1000", "n_bases_rejected": 500}]],
            "occ": ["l0seq/contig_occurrence/stage=ZQ/part-0.parquet", [{"contig_uid": c} for c in (A, B, H)]],
            "label": ["l0seq/contig_label/stage=ZQ/part-0.parquet", [{"contig_uid": A}]],
            "gcled": ["ledger/stage=ZQ/step=genecall/part-0.parquet",
                      [led(A, "ok"), led(B, "ok"), led(B, f"owned_by:ZQ/{g1.hex()}"), led(H, "owned_by:OA")]],
            "fp": [f"l1/feature_protein/gene_call_version={gcv}/stage=ZQ/part-0.parquet",
                   [{"contig_uid": A, "protein_uid": pr1}, {"contig_uid": B, "protein_uid": pr2}]],
            "ps": [f"l0seq/protein_seq/gene_call_version={gcv}/stage=ZQ/part-0.parquet",
                   [{"protein_uid": pr1, "length": 100}, {"protein_uid": pr2, "length": 200}]],
            "dom": ["l1/domains/library=pfam-38.2/stage=ZQ/part-0.parquet", [_dom(pr1, "PF1", 30.0), _dom(pr3, "PF3", 40.0)]],
            "feat1": [f"l1/feature/producer={p1}/stage=ZQ/part-0.parquet", [{"feature_uid": _u(31)}]],
            "elem1": [f"l1/crispr/producer={p1}/stage=ZQ/part-0.parquet", [{"feature_uid": _u(31)}]],
            "cled1": [f"ledger/stage=ZQ/step=crispr-{p1}/part-0.parquet", [cl(1, 1), cl(0, 0), cl(0, 0, "owned_by:OA")]],
            "feat2": [f"l1/feature/producer={p2}/stage=ZQ/part-0.parquet", [{"feature_uid": _u(32)}]],
            "elem2": [f"l1/crispr/producer={p2}/stage=ZQ/part-0.parquet", [{"feature_uid": _u(32)}]],
            "cled2": [f"ledger/stage=ZQ/step=crispr-{p2}/part-0.parquet", [cl(1, 1), cl(0, 0), cl(0, 0, "owned_by:OA")]],
            "oaci": ["l0seq/contig_index/stage=OA/part-0.parquet", [ci(H, 3000)]],
            "oafp": [f"l1/feature_protein/gene_call_version={gcv}/stage=OA/part-0.parquet",
                     [{"contig_uid": H, "protein_uid": pr3}]],
            "oaps": [f"l0seq/protein_seq/gene_call_version={gcv}/stage=OA/part-0.parquet",
                     [{"protein_uid": pr3, "length": 300}]],
            "oadom": ["l1/domains/library=pfam-38.2/stage=OA/part-0.parquet", [_dom(pr3, "PF3", 40.0)]],
        },
        "records": {
            "ZQ": {"inputs": ["gi"], "ingest": ["ci", "ingled", "occ", "label"], "genecall": ["gcled", "fp", "ps"],
                   f"crispr-{p1}": ["feat1", "elem1", "cled1"], "igr": [], "tapes": [], "domains-pfam-38.2": ["dom"],
                   "domains-pfam-38.2.plan": [], "domains-pfam-38.2.shard-00000": [],
                   f"crispr-{p2}": ["feat2", "elem2", "cled2"]},
            "OA": {"ingest": ["oaci"], "genecall": ["oafp", "oaps"], "domains-pfam-38.2": ["oadom"]},
        },
        "cpu": {},
        "plan": None,  # filled below once the paths are known
        "gcplan": {"batches": [{"batch": 0, "genomes": 3, "bases": 3500, "estimate_cpu_min": 0.3}],
                   "merge_estimate_cpu_min": 0.25, "genes_expected": 3},
        "conf": {"depth": "annotated", "library": "pfam 38.2 gathering", "crispr": "yes", "minlen": 1000, **conf_kw},
        "extra_conf": "", "p1": p1, "p2": p2, "u": {"A": A, "B": B, "H": H, "pr3": pr3, "pr1": pr1},
    }
    T = S["T"]
    S["plan"] = {"stage": "ZQ", "proteins": 3, "residues": 600, "shards": [{"proteins": 2, "residues": 300}],
                 "copy": {"proteins": 1, "residues": 300, "rows": 1, "proteins_with_hit": 1,
                          "sources": [{"stage": "OA", "protein_files": [{"path": T["oaps"][0]}],
                                       "files": [{"path": T["oadom"][0]}], "proteins": 1, "rows": 1}]},
                 "owned": {"stages": [{"stage": "OA", "feature_protein": [{"path": T["oafp"][0]}]}]}}
    if mut is not None:
        mut(S)
    root = tmp / "root"
    sizes = {k: _write(root, rel, rows) for k, (rel, rows) in T.items()}
    for st, recs in S["records"].items():
        d = root / "runs" / f"stage={st}"
        d.mkdir(parents=True, exist_ok=True)
        for step, keys in recs.items():
            rec = {"stage": st, "step": step, "exit_code": 0, "cpu_min": S["cpu"].get(step, 0.1), "wall_s": 1.0,
                   "peak_rss_bytes": 1 << 20, "threads": 1, "started_at": "2026-11-21T00:00:00", "over_cpu_budget": False,
                   "image_digest": "sha256:x", "code_id": "x",
                   "outputs": [{"path": T[k][0], "bytes": sizes[k]} for k in keys]}
            (d / f"{step}.json").write_text(json.dumps(rec))
    if S["plan"] is not None:
        p = root / "l1/domain_scan/library=pfam-38.2/stage=ZQ/plan.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(S["plan"]))
    from biobox.l1 import genecall as gc

    if S["gcplan"] is not None:
        gp = root / gc.plan_record_rel("ZQ")
        gp.parent.mkdir(parents=True, exist_ok=True)
        gp.write_text(json.dumps(S["gcplan"]))
    before, after = tmp / "before", tmp / "after"
    before.mkdir(parents=True)
    after.mkdir(parents=True)
    lines = sorted(f"{p.relative_to(root)}\t{p.stat().st_size}\t{p.stat().st_mtime}" for p in root.rglob("*")
                   if p.is_file())
    (after / "listing.tsv").write_text("".join(ln + "\n" for ln in lines))
    (before / "listing.tsv").write_text("".join(ln + "\n" for ln in lines if "/stage=ZQ/" not in "/" + ln.split("\t")[0]))
    cf = S["conf"]
    conf = tmp / "ZQ.conf"
    conf.write_text(conf_text("ZQ", tmp / "inputs", GOOD, depth=cf["depth"], library=cf["library"], crispr=cf["crispr"],
                              minlen=cf["minlen"]) + S["extra_conf"])
    prof = tmp / "profile.yaml"
    prof.write_text(PROFILE.replace("tmp: /x", f"tmp: {tmp}"))
    bj = tmp / "before.json"
    bj.write_text(json.dumps({"fasta_bytes": 10_000, "held_bases": 3000}))
    r = {"root": root, "conf": conf, "before": before, "after": after, "profile": prof, "before_json": bj,
         "out": tmp / "out", "S": S}
    if post is not None:  # a change to the snapshots or the root after they were made
        post(r)
    return r


@pytest.fixture
def no_minced(monkeypatch):
    """The producer names without the MinCED binary (as tests/r1/test_c47_est_accept.py does): a fixed version."""
    from biobox.l1 import crispr

    monkeypatch.setattr(crispr, "minced_version", lambda: "0.4.2")


def run_check(r: dict, monkeypatch, *extra: str) -> tuple[int, dict]:
    chk = _module("stage_check")
    argv = ["stage_check.py", "--conf", str(r["conf"]), "--stage", "ZQ", "--root", str(r["root"]), "--before",
            str(r["before"]), "--after", str(r["after"]), "--objects", str(r["root"] / "no-objects.json"), "--out",
            str(r["out"]), "--tmp", str(r["out"].parent / "chk-tmp"), "--profile", str(r["profile"]), "--before-est",
            str(r["before_json"]), *extra]
    monkeypatch.setattr(sys, "argv", argv)
    rc = chk.main()
    return rc, json.loads((r["out"] / "summary.json").read_text())


def test_3_the_check_passes_on_a_stage_built_as_it_should_be(tmp_path, monkeypatch, no_minced):
    rc, s = run_check(build_root(tmp_path), monkeypatch)
    assert rc == 0, s["failed"]
    assert s["contigs"]["ledger_skipped_other"] == 1 and s["contigs"]["ledger_skipped_within"] == 1
    assert s["pfam"]["copy_check"]["copied_proteins"] == 1


def _drop_ledger_row_of_h(S):
    S["T"]["gcled"][1] = [x for x in S["T"]["gcled"][1] if x["object_key"] != S["u"]["H"]]


def _earlier_reason_with_a_suffix(S):
    S["T"]["gcled"][1] = [dict(x, reason_code="owned_by:OA/abc") if x["reason_code"] == "owned_by:OA" else x
                          for x in S["T"]["gcled"][1]]


def _uncovered_contig(S):
    D, g4 = _u(4), _u(14)
    S["T"]["ci"][1].append({"contig_uid": D, "contig_len": 1200, "is_default_variant": True})
    S["T"]["gi"][1].append({"contig_uid": D, "genome_uid": g4, "kind": "virus", "genetic_code": 11,
                            "is_circular": False, "assembly_quality": "complete"})


def _genes_on_a_contig_not_stored(S):
    S["T"]["fp"][1].append({"contig_uid": S["u"]["H"], "protein_uid": _u(24)})


def _copied_row_differs(S):
    S["T"]["dom"][1][1] = _dom(S["u"]["pr3"], "PF3", 41.0)


def _crispr_count_differs(S):
    S["T"]["cled1"][1][0] = {"n_in": 2, "n_kept": 2, "n_rejected": 0, "reason_code": "ok"}


def _no_scan_plan(S):
    S["plan"] = None


def _crispr_record_without_files(S):
    S["records"]["ZQ"][f"crispr-{S['p1']}"] = []


def _judged_size_step_unnamed(S):
    S["cpu"]["tapes"] = 30.0


@pytest.mark.parametrize("mut,conf_kw,expect", [
    (_drop_ledger_row_of_h, {}, "ledger skips for another stage that differ"),
    (_earlier_reason_with_a_suffix, {}, "ledger skips for another stage that differ"),
    (_uncovered_contig, {}, "stored contigs have no gene row and no ledger row"),
    (_genes_on_a_contig_not_stored, {}, "contigs have gene rows in ZQ but are not stored by it"),
    (_copied_row_differs, {}, "copied rows from OA differ"),
    (_crispr_count_differs, {}, "default: ledger keeps 2, feature files hold 1"),
    (None, {"depth": "called", "library": None}, "domain records without a library"),
    (None, {"minlen": 1600}, "stored contigs under the least length 1600"),
    (_no_scan_plan, {}, "no scan plan"),
    (_crispr_record_without_files, {}, "lists no ledger file"),
    (_judged_size_step_unnamed, {}, "records of judged size"),
], ids=["ledger-count", "ledger-reason-suffix", "uncovered", "genes-not-stored", "copy-differs", "crispr-count",
        "step-of-another-depth", "under-least-length", "no-scan-plan", "crispr-no-files", "judged-size-unnamed"])
def test_3_each_assertion_of_the_check_fails_on_its_defect(tmp_path, monkeypatch, no_minced, mut, conf_kw, expect):
    rc, s = run_check(build_root(tmp_path, mut, **conf_kw), monkeypatch)
    assert rc == 1
    assert any(expect in f for f in s["failed"]), s["failed"]


def test_7_a_crispr_producer_that_found_no_array_passes_when_its_ledger_says_so(tmp_path, monkeypatch, no_minced):
    def none_found(S):
        S["records"]["ZQ"][f"crispr-{S['p1']}"] = ["cled1"]
        S["T"]["cled1"][1][0] = {"n_in": 0, "n_kept": 0, "n_rejected": 0, "reason_code": "ok"}

    rc, s = run_check(build_root(tmp_path, none_found), monkeypatch)
    assert rc == 0, s["failed"]


def test_3_the_check_wants_a_copy_in_the_rehearsal(tmp_path, monkeypatch, no_minced):
    def no_copy(S):
        S["plan"]["copy"] = {"proteins": 0, "rows": 0, "sources": []}
        S["T"]["dom"][1] = S["T"]["dom"][1][:1]

    rc, s = run_check(build_root(tmp_path, no_copy), monkeypatch, "--expect-copy")
    assert rc == 1 and any("copied no protein" in f for f in s["failed"])


def test_3_rates_on_a_published_plan_with_shards_a_judge_line_and_an_accepted_record(tmp_path, no_minced):
    from biobox import profile as bprofile

    from biobox.l1 import genecall

    est = _module("stage_est")
    batch = f"genecall.{genecall.batch_name(0)}"

    def records(S):
        S["records"]["ZQ"][batch] = []
        S["cpu"].update({batch: 40.0, "domains-pfam-38.2.shard-00000": 0.5})
        S["gcplan"]["batches"][0]["estimate_cpu_min"] = 30.0

    r = build_root(tmp_path, records)
    judge = "JUDGE=genecall 25 planning/measurements/2026-11-20_stage_V1.md\n"
    r["conf"].write_text(r["conf"].read_text() + judge)
    doc = est.rates_doc(est.read_conf(r["conf"], "ZQ"), r["root"], bprofile.load(r["profile"]), tmp_path / "none",
                        str(r["before_json"]), tmp_path / "t")
    by = {x["step"]: x for x in doc["rows"]}
    assert by[batch]["over_20pct"], "40 against the plan's 30"
    assert by["genecall.note"]["cpu_min"] == pytest.approx(40.1) and by["genecall.note"]["over_20pct"]
    assert by["domains-pfam-38.2.shard-00000"]["judged"] is False, "an estimate under 24 is never judged (DEC-43)"
    assert est.rates_exit(doc) == 3
    r["conf"].write_text(r["conf"].read_text() + "ACCEPTED=genecall.note 40.100 25.000 DEC-56\n"
                         f"ACCEPTED={batch} 40.000 30.000 DEC-56\n")
    doc = est.rates_doc(est.read_conf(r["conf"], "ZQ"), r["root"], bprofile.load(r["profile"]), tmp_path / "none",
                        str(r["before_json"]), tmp_path / "t")
    assert doc["over_20pct"] == [] and set(doc["accepted_over_20pct"]) == {"genecall.note", batch}


def test_9_rates_takes_what_does_not_change_from_the_estimate(tmp_path, no_minced):
    from biobox import profile as bprofile

    est = _module("stage_est")
    r = build_root(tmp_path)
    doc = est.rates_doc(est.read_conf(r["conf"], "ZQ"), r["root"], bprofile.load(r["profile"]),
                        tmp_path / "no-objects.json", str(r["before_json"]), tmp_path / "t")
    ing = next(x for x in doc["rows"] if x["step"] == "ingest")
    assert "stage.ingest_cpu_min(10000, held=3000)" in ing["basis"], "no FASTA read, no held bases counted again"


def test_rates_of_a_stage_called_in_one_process_with_a_scratch_directory_not_yet_made(tmp_path, no_minced):
    """Found by the rehearsal of 529b6b0 (2026-11-21, ZX1R): with fasta_bytes and held_bases read from before.json,
    nothing made --tmp, and genecall.plan's mkdtemp under it failed; rates exited 1 and the guard stopped the build."""
    from biobox import profile as bprofile

    est = _module("stage_est")

    def one_process(S):
        S["gcplan"] = None  # gene calling in one process: rates asks genecall.plan for the estimate

    def stages_on_disk(r):
        for s in ("ZQ", "OA"):
            (r["root"] / "l0seq" / "seq_locator" / f"stage={s}").mkdir(parents=True)

    r = build_root(tmp_path, one_process, stages_on_disk)
    doc = est.rates_doc(est.read_conf(r["conf"], "ZQ"), r["root"], bprofile.load(r["profile"]), tmp_path / "none",
                        str(r["before_json"]), tmp_path / "not" / "made" / "yet")
    assert "genecall" in {x["step"] for x in doc["rows"]}


def test_3_complete_names_a_missing_step(tmp_path, monkeypatch, no_minced):
    est = _module("stage_est")
    r = build_root(tmp_path)
    (r["root"] / "runs" / "stage=ZQ" / "tapes.json").unlink()
    monkeypatch.setattr(sys, "argv", ["stage_est.py", "complete", "--conf", str(r["conf"]), "--stage", "ZQ", "--root",
                                      str(r["root"])])
    assert est.main() == 1
    r2 = build_root(tmp_path / "b")
    monkeypatch.setattr(sys, "argv", ["stage_est.py", "complete", "--conf", str(r2["conf"]), "--stage", "ZQ", "--root",
                                      str(r2["root"])])
    assert est.main() == 0


def test_6_the_check_reads_reasons_no_looser_than_c47_check():
    import duckdb

    chk = _module("stage_check")
    con = duckdb.connect()
    con.execute("CREATE TABLE led (reason_code VARCHAR)")
    con.execute("INSERT INTO led VALUES ('duplicateXcontig'), ('owned_byXOA'), ('duplicate_contig:a')")
    got = [r[0] for r in con.execute(f"SELECT {chk.skip_kind_sql('reason_code', 'C5')} FROM led").fetchall()]
    assert got == [None, None, "within"], "LIKE's `_` would match any character"


def test_6_estimate_stops_when_an_expect_line_differs(tmp_path):
    w = World(tmp_path)
    got = w.drive("estimate", FAKE_EST_RC="8")
    assert got.returncode != 0 and not w.done("estimate")
    assert "differ from the note's" in _driver_log(w)


def test_6_the_real_estimate_compares_an_expect_line_with_the_root(tmp_path):
    est = _module("stage_est")
    inp = tmp_path / "inputs"
    inp.mkdir()
    (inp / "x.fna").write_text(">a\n" + "ACGT" * 2000 + "\n")
    from biobox import ids

    cu = ids.contig_uid(b"ACGT" * 2000)
    (inp / "genomes.json").write_text(json.dumps([{"contig_uids": [cu.hex()], "genetic_code": 11, "kind": "virus",
                                                  "circular": [], "assembly_quality": "contig"}]))
    (inp / "objects.json").write_text(json.dumps([{"path": str(inp / "x.fna"), "source_name": "t",
                                                   "source_release": "1", "license": "t", "redistributable": True}]))
    conf = _conf_file(tmp_path, conf_text("ZZ1", inp, GOOD, depth="called", library=None,
                                          extra="EXPECT=held_by OA 1 planning/measurements/2026-11-20_stage_V1.md\n"))
    root = tmp_path / "root"
    (root / "l0seq" / "contig_index").mkdir(parents=True)
    prof = tmp_path / "profile.yaml"
    prof.write_text(PROFILE.replace("tmp: /x", f"tmp: {tmp_path}"))
    got = subprocess.run([sys.executable, str(BSP / "stage_est.py"), "before", "--conf", str(conf), "--stage", "ZZ1",
                          "--root", str(root), "--profile", str(prof), "--inputs", str(inp), "--tmp",
                          str(tmp_path / "t")], capture_output=True, text=True, timeout=300,
                         env={**os.environ, "PYTHONPATH": str(REPO / "src")})
    assert got.returncode == 8, got.stderr
    assert json.loads(got.stdout)["expect"]["OA"] == {"note": "planning/measurements/2026-11-20_stage_V1.md",
                                                      "expected": 1, "held": 0}
    assert est.read_conf(conf, "ZZ1").expect == {"OA": (1, "planning/measurements/2026-11-20_stage_V1.md")}


def test_6_each_check_keeps_its_report_in_a_numbered_directory(tmp_path):
    w = World(tmp_path)
    _check_ready(w)
    _ok(w.drive("check"))
    _ok(w.drive("check"))
    assert (w.c / "check" / "1" / "check.out").exists() and (w.c / "check" / "2" / "check.out").exists()


# (4) markers are tied to the conf they were made for; run checks the inputs' sums again

def test_4_a_marker_made_for_another_conf_is_refused(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    w.write_conf(w.stage, (w.bin / "stages" / f"{w.stage}.conf").read_text() + "# changed\n")
    got = w.drive("run")
    assert got.returncode != 0 and "pod_run" not in w.calls()
    assert "made for another conf" in _driver_log(w)


def test_4_run_checks_the_inputs_sums_again(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    (w.inputs / "x.fna").write_text(">a\nACGA\n")
    got = w.drive("run")
    assert got.returncode != 0 and "pod_run" not in w.calls()


def test_4_a_rehearsal_marker_of_another_conf_is_redone_and_withdraws_what_rests_on_it(tmp_path):
    w = World(tmp_path)
    (w.c / "rehearse").mkdir(parents=True)
    (w.c / "rehearse" / "DONE").write_text("0" * 64 + "\n")
    w.mark("run", "after")
    _ok(w.drive("rehearse"))
    assert "rehearse: mode estimate passed" in _driver_log(w)
    assert not w.done("run") and not w.done("after")
    assert (w.c / "rehearse" / "DONE").read_text().strip() == w.conf_sha()


# (5) the unit's quota: the profile's cap of the hour less the helpers', changed when the hour's cap changes

def test_5_the_cap_of_the_hour_comes_from_the_profile_by_the_code_of_run_step():
    from datetime import datetime

    est = _module("stage_est")
    prof = REPO / "profiles" / "wilkcomp02.yaml"
    assert est.cap_at(prof, datetime(2026, 11, 21, 3, 0)) == 56, "Saturday night"
    assert est.cap_at(prof, datetime(2026, 11, 23, 7, 59)) == 56, "Monday before 08:00"
    assert est.cap_at(prof, datetime(2026, 11, 23, 10, 0)) == 40, "Monday by day"
    assert est.cap_at(prof, datetime(2026, 11, 23, 20, 0)) == 56, "Monday evening"


def test_5_the_unit_quota_follows_the_hours_cap(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("10")
    quotas = _lines(tmp_path / "quotas", "54", "54", "38")
    _ok(w.drive("run", FAKE_ACTIVE_FILE=str(active), FAKE_QUOTAS=quotas), w, "run-quota")
    assert "pod_run-opt --threads 54" in w.calls()
    assert "systemctl --user set-property --runtime zz1-1 CPUQuota=3800%" in w.calls()


def test_5_q_only_lowers_the_quota(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    _ok(w.drive("run", Q="10"))
    assert "pod_run-opt --threads 10" in w.calls()


# (8) the environment: the data root is compared after realpath, and keeps its conf, memory and floors

def test_8_a_trailing_slash_does_not_make_the_data_root_a_scratch_root(tmp_path):
    w = World(tmp_path)
    w.mark("estimate", "before")
    got = w.drive("run", STAGE_REHEARSAL="1", STAGE_ROOT=str(w.root) + "/")
    assert got.returncode != 0 and "pod_run" not in w.calls()


@pytest.mark.parametrize("env", [{"MEM": "300G"}, {"TMP_MIN_GB": "10"}, {"START_MIN_GB": "100"},
                                 {"SCRATCH_MAX_GB": "500"}])
def test_8_the_data_root_never_loosens_its_limits(tmp_path, env):
    w = World(tmp_path)
    _run_ready(w)
    got = w.drive("run", **env)
    assert got.returncode != 0 and "pod_run" not in w.calls()


def test_8_in_the_data_root_the_conf_is_the_trees_own(tmp_path):
    w = World(tmp_path)
    other = tmp_path / "ZZ1.conf"
    other.write_text((w.bin / "stages" / "ZZ1.conf").read_text())
    got = w.drive("estimate", STAGE_CONF=str(other))
    assert got.returncode != 0 and "stage_est.py before" not in w.calls()


def test_8_an_accepted_record_without_its_decision_row_is_refused_by_the_driver(tmp_path):
    w = World(tmp_path)
    w.write_conf(w.stage, (w.bin / "stages" / "ZZ1.conf").read_text() + "ACCEPTED=ingest 1.000 1.000 DEC-999\n")
    got = w.drive("estimate")
    assert got.returncode != 0 and "stage_est.py before" not in w.calls()


# (9) smaller

def test_9_backup_wants_two_hashes_of_64_hex_digits(tmp_path):
    w = World(tmp_path)
    _backup_ready(w)
    got = w.drive("backup", TAR_FLOOR_BYTES="0", FAKE_SHA_EMPTY="1")
    assert got.returncode != 0 and not w.done("backup")


def test_9_a_stop_during_run_still_judges_the_records_published(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("2")
    got = w.drive("run", FAKE_ACTIVE_FILE=str(active), FAKE_STOP_ON_ACTIVE=str(w.root / "STOP"),
                  FAKE_RATES_RC="3")
    assert got.returncode == 1 and (w.c / "run" / "FAILED").exists(), "a record over 20% fails, STOP or not"


def test_9_each_rehearsal_starts_with_its_own_logs(tmp_path):
    w = World(tmp_path)
    _ok(w.drive("rehearse"))
    (w.c / "rehearse" / "DONE").unlink()
    _ok(w.drive("rehearse"))
    old = w.c / "rh-logs.1"
    assert (old / "inner-driver.log").exists() and (old / "zz1-rh-e.log").exists()
    assert (w.c / "rh" / "drv" / "driver.log").read_text().count("mode estimate done") == 1, "the new log only"
    assert (old / "inner-driver.log").read_text().count("mode estimate done") == 1


def test_9_the_queue_reads_run_failed_before_it_runs_estimate_again(tmp_path):
    (tmp_path / "stages" / "zz1" / "run").mkdir(parents=True)
    (tmp_path / "stages" / "zz1" / "run" / "FAILED").write_text("guard: over 280 CPU-minutes")
    got = _queue(tmp_path)
    assert got.returncode == 3 and _modes(tmp_path) == []


# ==============================================================================================================
# Card A-12, third round (2026-11-21): points A to F of the second review

# A. the release side refuses a tree whose lock could not see the stage driver's units

def test_A_the_release_train_refuses_a_tree_without_stage_confs(tmp_path):
    w = World(tmp_path)
    (w.bin / "stages" / "ZZ1.conf").unlink()
    train = w.script("release_train.sh")
    got = subprocess.run(["bash", str(train), "run", "--no-live-cpu"], capture_output=True, text=True, timeout=120,
                         env=w.env(ROOT=str(w.root)))
    assert got.returncode != 0
    assert "no stage confs" in got.stderr
    assert "train_plan.py" not in w.calls()


def test_A_rel_wait_takes_the_live_driver_trees_stages_through_a_variable(tmp_path):
    pod = tmp_path / "if2453"
    d = pod / "stages" / "v2"
    for m in ("before", "run", "after", "check", "backup"):
        (d / m).mkdir(parents=True)
        (d / m / "DONE").touch()
    (d / "backup.tsv").write_text("x\t1\tab\n")
    (pod / "c46" / "wait").mkdir(parents=True)
    (pod / "r10" / "root").mkdir(parents=True)
    (pod / "ci" / "old1").mkdir(parents=True)  # a frozen tree from before card A-12: no stages/
    live = tmp_path / "live" / "stages"
    live.mkdir(parents=True)
    (live / "GV1.conf").write_text("STAGE=GV1\n")
    fake = tmp_path / "fakebin"
    fake.mkdir()
    for name, text in (("systemctl", SYSTEMCTL), ("sleep", NOOP)):
        (fake / name).write_text(text)
        (fake / name).chmod(0o755)
    sh = tmp_path / "rel_wait.sh"
    sh.write_text((ORCH / "rel_wait.sh").read_text().replace("/tmp/if2453", str(pod)))
    log = tmp_path / "fake.log"
    log.write_text("")
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log), "FAKE_UNITS": "gv1-1",
           "REL_STAGES_DIR": str(live)}
    got = subprocess.run(["bash", str(sh), str(d), "old1"], env=env, capture_output=True, text=True, timeout=120)
    assert got.returncode == 6 and "gv1-1" in got.stdout


# B. the watch has a time limit without an answer, and so has stop_unit; neither claims the unit stopped

def test_B_no_answer_for_the_watch_limit_fails_saying_the_state_is_unknown(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    states = _lines(tmp_path / "states", *(["fail"] * 10))
    got = w.drive("run", FAKE_STATES=states, WATCH_LIMIT_S="0")
    assert got.returncode != 0 and not w.done("run")
    failed = (w.c / "run" / "FAILED").read_text()
    assert "state is unknown" in failed and "may still be running" in failed
    assert "stopped" not in failed.replace("may still be running", "")


def test_B_stop_unit_gives_up_loudly_at_its_bound_and_never_claims_the_unit_stopped(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("1000")
    got = w.drive("run", FAKE_ACTIVE_FILE=str(active), FAKE_GUARD_RC="3", STOP_LIMIT_S="0")
    assert got.returncode != 0
    failed = (w.c / "run" / "FAILED").read_text()
    assert "still not seen gone" in failed and "(unit zz1-1 stopped)" not in failed


def test_B_stop_unit_without_an_answer_ends_at_the_watch_limit(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    states = _lines(tmp_path / "states", "active", *(["fail"] * 20))
    got = w.drive("run", FAKE_STATES=states, FAKE_GUARD_RC="3", WATCH_LIMIT_S="0")
    assert got.returncode != 0
    assert "state is unknown" in (w.c / "run" / "FAILED").read_text()


# C. a TERM to the driver during a long helper call runs the trap at once

def test_C_a_stop_during_a_long_helper_call_stops_the_unit_at_once(tmp_path):
    import signal
    import time

    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("300")
    drv = w.script("stage_driver.sh")
    p = subprocess.Popen(["bash", str(drv), w.stage, "run"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                         env=w.env(FAKE_ACTIVE_FILE=str(active), FAKE_RATES_SLEEP="60", FAKE_SLEEP="0.05"))
    t0 = time.monotonic()
    while "stage_est.py rates" not in w.calls() and time.monotonic() - t0 < 60:
        time.sleep(0.2)
    t1 = time.monotonic()
    p.send_signal(signal.SIGTERM)
    p.communicate(timeout=120)
    assert time.monotonic() - t1 < 30, "the trap waited for the helper (a 60 s call)"
    assert "systemctl --user stop zz1-1" in w.calls()


# D. a conf that is not the tree's own never reaches a mode that reads or writes the data root; rehearse may take one

@pytest.mark.parametrize("mode", ["estimate", "before", "run", "after", "check", "backup"])
def test_D_a_conf_not_the_trees_own_never_reaches_a_data_root_mode(tmp_path, mode):
    w = World(tmp_path)
    w.mark("estimate", "rehearse", "before", "run", "after", "check")
    other = tmp_path / "confs" / "ZZ1.conf"
    other.parent.mkdir()
    other.write_text((w.bin / "stages" / "ZZ1.conf").read_text())
    got = w.drive(mode, STAGE_CONF=str(other), TAR_FLOOR_BYTES="0")
    assert got.returncode != 0
    assert "the conf is the tree's own" in _driver_log(w)
    assert "singularity" not in w.calls() and "pod_run" not in w.calls() and "sha256sum" not in _driver_log(w)


def test_D_rehearse_may_take_a_conf_that_is_not_the_trees_own(tmp_path):
    w = World(tmp_path)
    other = tmp_path / "confs" / "ZZ1.conf"
    other.parent.mkdir()
    other.write_text((w.bin / "stages" / "ZZ1.conf").read_text())
    (w.bin / "stages" / "ZZ1.conf").unlink()
    _ok(w.drive("rehearse", STAGE_CONF=str(other)), w, "rehearse-other-conf")
    assert (w.c / "rh" / "drv" / "backup" / "DONE").exists()


# E. the profile in the data root is the tree's own

def test_E_a_profile_not_the_trees_own_is_refused_in_the_data_root(tmp_path):
    w = World(tmp_path)
    other = tmp_path / "bigger.yaml"
    other.write_text(PROFILE.replace("max_threads_day: 40", "max_threads_day: 112"))
    got = w.drive("estimate", PROFILE=str(other))
    assert got.returncode != 0 and "the profile is the tree's own" in _driver_log(w)
    assert "stage_est.py before" not in w.calls()


# F. the check's other assertions, each failing on its own message

def _listing_line(r: dict, rel: str) -> str:
    p = r["root"] / rel
    return f"{rel}\t{p.stat().st_size}\t{p.stat().st_mtime}\n"


def _gone(r):
    with open(r["before"] / "listing.tsv", "a") as fh:
        fh.write("l0seq/contig_index/stage=OA/part-9.parquet\t5\t1.0\n")


def _changed(r):
    rel = "l0seq/contig_index/stage=OA/part-0.parquet"
    t = (r["before"] / "listing.tsv").read_text().replace(_listing_line(r, rel), f"{rel}\t1\t1.0\n")
    (r["before"] / "listing.tsv").write_text(t)


def _sha_changed(r):
    rel = "l0seq/contig_index/stage=OA/part-0.parquet"
    (r["before"] / "earlier.sha256").write_text(f"{'a' * 64}  {rel}\n")
    (r["after"] / "earlier.sha256").write_text(f"{'b' * 64}  {rel}\n")


def _unexpected_new(r):
    p = r["root"] / "l0seq" / "stray.parquet"
    p.write_bytes(b"x")
    with open(r["after"] / "listing.tsv", "a") as fh:
        fh.write(_listing_line(r, "l0seq/stray.parquet"))


def _held_by_two(S):
    S["T"]["obci"] = ["l0seq/contig_index/stage=OB/part-0.parquet",
                      [{"contig_uid": S["u"]["H"], "contig_len": 3000, "is_default_variant": True}]]
    S["records"]["OB"] = {"ingest": ["obci"]}


def _stored_and_held(S):
    S["T"]["oaci"][1].append({"contig_uid": S["u"]["A"], "contig_len": 2000, "is_default_variant": True})


def _neither(S):
    S["T"]["gi"][1].append({"contig_uid": _u(5), "genome_uid": _u(15), "kind": "virus", "genetic_code": 11,
                            "is_circular": False, "assembly_quality": "complete"})


def _no_gene_rows(S):
    S["records"]["ZQ"]["genecall"] = ["gcled", "ps"]


def _plan_count(S):
    S["plan"]["proteins"] = 4


@pytest.mark.parametrize("mut,post,expect", [
    (None, _gone, "gone: 1"),
    (None, _changed, "changed: 1"),
    (None, _sha_changed, "other stage's sha256 changed: 1"),
    (None, _unexpected_new, "unexpected new: 1"),
    (_held_by_two, None, "contigs held by two other stages (DEC-26): 1"),
    (_stored_and_held, None, "contigs both stored by the stage and held by another: 1"),
    (_neither, None, "genome_input contigs neither stored nor held: 1"),
    (_no_gene_rows, None, "no gene rows"),
    (_plan_count, None, "scan plan holds 4 proteins, expected 3"),
], ids=["file-gone", "file-changed", "sha-changed", "unexpected-new", "held-by-two", "stored-and-held", "neither",
        "no-gene-rows", "plan-count"])
def test_F_each_further_assertion_of_the_check_fails_on_its_defect(tmp_path, monkeypatch, no_minced, mut, post, expect):
    rc, s = run_check(build_root(tmp_path, mut, post), monkeypatch)
    assert rc == 1
    assert any(f.startswith(expect) or expect in f for f in s["failed"]), s["failed"]


# smaller: the quota is looked at as soon as the wall-clock hour changes

def test_the_quota_is_looked_at_within_one_poll_of_the_hour_changing(tmp_path):
    w = World(tmp_path)
    _run_ready(w)
    active = tmp_path / "active"
    active.write_text("3")  # three polls: far from the sixth
    quotas = _lines(tmp_path / "quotas", "54", "54", "38")
    hours = _lines(tmp_path / "hours", "07", "07", "08")
    _ok(w.drive("run", FAKE_ACTIVE_FILE=str(active), FAKE_QUOTAS=quotas, FAKE_HOURS=hours))
    assert "systemctl --user set-property --runtime zz1-1 CPUQuota=3800%" in w.calls()

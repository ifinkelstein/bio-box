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
SYSTEMCTL = """#!/bin/sh
echo "systemctl $*" >> "$FAKE_LOG"
case "$*" in
  *is-active*)
    if [ -n "${FAKE_STOP_ON_ACTIVE:-}" ]; then touch "$FAKE_STOP_ON_ACTIVE"; fi
    if [ -n "${FAKE_ACTIVE_FILE:-}" ] && [ -s "$FAKE_ACTIVE_FILE" ]; then
      n=$(cat "$FAKE_ACTIVE_FILE")
      if [ "$n" -gt 0 ]; then echo $((n - 1)) > "$FAKE_ACTIVE_FILE"; exit 0; fi
    fi
    exit 3 ;;
  *list-units*)
    for u in ${FAKE_UNITS:-}; do
      for a in "$@"; do
        case "$a" in --*|list-units) continue ;; esac
        case "$u" in $a) echo "$u loaded active running x"; break ;; esac
      done
    done
    exit 0 ;;
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
  *"stage_est.py rates"*) echo '{"over_20pct": []}'; exit "${FAKE_RATES_RC:-0}" ;;
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
echo "pod_run: $unit exited $rc" >> "$log"
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
        for name, text in (("systemd-run", SYSTEMD_RUN), ("nice", NICE), ("ionice", IONICE), ("sleep", NOOP),
                           ("top", NOOP), ("uptime", NOOP), ("systemctl", SYSTEMCTL), ("singularity", SINGULARITY),
                           ("python3", PYTHON3)):
            (self.fake / name).write_text(text)
            (self.fake / name).chmod(0o755)
        (self.tree / "scripts" / "pod_run.sh").write_text(POD_RUN)
        (self.tree / "profiles" / "wilkcomp02.yaml").write_text(PROFILE)
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
                "REAL_PYTHON": sys.executable, "TREE": str(self.tree), "BACKUP_DIR": str(self.stor),
                "START_MIN_GB": "0", "TMP_MIN_GB": "0", **extra}

    def drive(self, mode: str, *more: str, stage: str | None = None, **extra: str) -> subprocess.CompletedProcess:
        drv = self.script("stage_driver.sh")
        return subprocess.run(["bash", str(drv), self.stage if stage is None else stage, mode, *more],
                              env=self.env(**extra),
                              capture_output=True, text=True, timeout=300)

    def mark(self, *modes: str) -> None:
        for m in modes:
            (self.c / m).mkdir(parents=True, exist_ok=True)
            (self.c / m / "DONE").touch()

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
    c = est.read_conf(path, path.stem)
    text = _notes_text(c)
    for h, name in c.sums:
        assert h in text, f"{path.name}: sum of {name} is not in its note"
    for step, (fig, src) in c.judge.items():
        assert (REPO / src).is_file()
        assert f"{fig:g}" in (REPO / src).read_text(), f"{path.name}: {step} {fig} is not in {src}"
    assert not c.accepted, "a new stage carries no accepted record"
    assert c.inputs == Path(f"/tmp/if2453/s2/inputs/{c.stage}")


def test_gv1_comes_after_v2_and_judges_its_gene_calling_against_the_notes_figure():
    est = _module("stage_est")
    c = est.read_conf(STAGES / "GV1.conf", "GV1")
    assert c.after == ("V2",)
    assert c.judge["genecall"][0] == 236.0


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
    assert "pod_run-opt --threads 56" in w.calls(), "the unit's quota is the profile's night cap"
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
    assert added == (w.c / "check" / "root_added.tsv").read_text()


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

C5_GENOMES_SHA = "1c65060c80ef120d4550d702a76c4a9c5c3eff17ee25f2694d5a172873468573"
C5_OBJECTS_SHA = "09d742d70a05528e0bd91b735605b63e0f67ad3758e9f1d02107be30bf5483ad"


def test_c47_and_the_stage_driver_launch_the_same_build_stage_command_for_c5(tmp_path):
    w = World(tmp_path, stage="C5")
    c = w.pod / "c47"
    inputs = w.pod / "s2" / "inputs"
    w.write_conf("C5", conf_text("C5", inputs / "C5", {"genomes.json": C5_GENOMES_SHA, "objects.json": C5_OBJECTS_SHA},
                                 extra="ACCEPTED=ingest 45.803 38.057 DEC-56\n"))
    # c47_driver.sh run, with its own overrides and its markers
    for m in ("before/DONE", "rh/OK"):
        (c / m).parent.mkdir(parents=True, exist_ok=True)
        (c / m).touch()
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
    for m in ("estimate", "rehearse"):
        (c / m).mkdir(parents=True, exist_ok=True)
        (c / m / "DONE").touch()
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
        assert "REAL_ROOT=/tmp/if2453/r10/root" in (BSP / name).read_text() or name == "c47_driver.sh"
        assert str(w.pod) not in (BSP / name).read_text()

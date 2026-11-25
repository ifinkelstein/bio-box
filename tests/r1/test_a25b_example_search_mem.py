"""Ring R1: card A-25b — derive `search hmm --mem-gb` from `MEM` in worked examples and verify clean resume of `$WORK/examples`."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from biobox.search import runner

pytestmark = pytest.mark.covers("A-25", "E-11", "A-20")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
LAUNCH = REPO / "examples" / "_shared" / "launch.sh"
CRISPR_RUN = REPO / "examples" / "crispr_cas_neighbourhoods" / "run.sh"
RNAP_RUN = REPO / "examples" / "rnap" / "run.sh"


@pytest.mark.parametrize("env_extra,expected", [
    ({}, "14"),
    ({"MEM": "16G"}, "14"),
    ({"MEM": "64G"}, "62"),
    ({"MEM": "64Gi"}, "62"),
    ({"MEM": "64GB"}, "62"),
    ({"MEM": "64GiB"}, "62"),
    ({"MEM": "64g"}, "62"),
    ({"MEM": "08G"}, "6"),
    ({"MEM": "2G"}, "1"),
    ({"MEM": "1G"}, "1"),
    ({"MEM": "16384M"}, "14"),
    ({"MEM": "65536M"}, "62"),
    ({"MEM": "512M"}, "1"),
    ({"MEM": "1T"}, "1022"),
    ({"MEM": "64G", "SEARCH_MEM_GB": "48"}, "48"),
    ({"MEM": "16G", "SEARCH_MEM_GB": "10"}, "10"),
])
def test_search_mem_gb_derives_from_mem_and_respects_override(env_extra, expected):
    env = {k: v for k, v in os.environ.items() if k not in ("MEM", "SEARCH_MEM_GB")}
    env.update(env_extra)
    got = subprocess.run(
        ["bash", "-c", f'. "{LAUNCH}" && search_mem_gb'],
        env=env, capture_output=True, text=True, timeout=30,
    )
    assert got.returncode == 0, got.stderr
    assert got.stdout.strip() == expected, f"got {got.stdout!r}, stderr={got.stderr!r}"


@pytest.mark.parametrize("script_path", [CRISPR_RUN, RNAP_RUN])
def test_example_run_scripts_do_not_hardcode_mem_gb_14(script_path):
    text = script_path.read_text()
    assert "--mem-gb 14" not in text, f"{script_path} still hardcodes --mem-gb 14"
    assert "$(search_mem_gb)" in text, f"{script_path} must derive --mem-gb via search_mem_gb"


def _make_example_stub_world(tmp_path: Path, example: str) -> tuple[Path, Path, Path, dict[str, str]]:
    """Minimal stubbed host environment for running an example's `run.sh` through its `search hmm` call(s)."""
    tree = tmp_path / "tree"
    fake = tmp_path / "fakebin"
    root = tmp_path / "root"
    out = tmp_path / "out"
    pfam_dir = tmp_path / "pfam"
    capture = tmp_path / "pod_run.capture"
    for d in (tree / "scripts", tree / "examples" / "_shared", tree / "examples" / example, fake, root, pfam_dir):
        d.mkdir(parents=True, exist_ok=True)

    (tree / "examples" / "_shared" / "launch.sh").write_text(LAUNCH.read_text())
    (tree / "examples" / "_shared" / "check.py").write_text("#!/usr/bin/env python3\n")
    src_run = REPO / "examples" / example / "run.sh"
    (tree / "examples" / example / "run.sh").write_text(src_run.read_text())
    (tree / "examples" / example / "expected.tsv").write_text("check_id\tkind\texpected\n")
    pfam = pfam_dir / "Pfam-A.hmm.gz"
    pfam.write_text("fake pfam\n")
    capture.write_text("")

    # Stub pod_run.sh: record full argv, write a fake exit-0 log for `find`, and fail cleanly right after
    # recording `search hmm` so the test stops once the search invocation(s) have been captured.
    (tree / "scripts" / "pod_run.sh").write_text(
        "#!/bin/bash\n"
        "log=\"\"; unit=\"\"; mem=\"\"; threads=\"\"\n"
        "while [ $# -gt 0 ]; do\n"
        "  case \"$1\" in\n"
        "    --log) log=$2; shift 2 ;;\n"
        "    --mem) mem=$2; shift 2 ;;\n"
        "    --threads) threads=$2; shift 2 ;;\n"
        "    --*) shift 2 ;;\n"
        "    *) unit=$1; shift; break ;;\n"
        "  esac\n"
        "done\n"
        "echo \"mem=$mem threads=$threads cmd=$*\" >> \"$CAPTURE_LOG\"\n"
        "mkdir -p \"$(dirname \"$log\")\"\n"
        "if [ \"$1\" = \"search\" ]; then\n"
        "  echo \"pod_run: $unit exited 99\" >> \"$log\"\n"
        "else\n"
        "  echo \"pod_run: $unit exited 0\" >> \"$log\"\n"
        "fi\n"
        "exit 0\n"
    )

    # Stub systemd-run (used by `tool` for hmmfetch) to emit minimal HMM files with GA lines.
    (fake / "systemd-run").write_text(
        "#!/bin/bash\n"
        "out_hmm=\"\"\n"
        "prev=\"\"\n"
        "for a in \"$@\"; do\n"
        "  if [ \"$prev\" = \"-o\" ]; then out_hmm=$a; fi\n"
        "  prev=$a\n"
        "done\n"
        "if [ -n \"$out_hmm\" ]; then\n"
        "  printf 'HMMER3/f\\nNAME  Test\\nGA    25.00 25.00;\\n//\\n' > \"$out_hmm\"\n"
        "fi\n"
        "exit 0\n"
    )
    (fake / "systemd-run").chmod(0o755)
    (fake / "ps").write_text("#!/bin/sh\necho 0.0\n")
    (fake / "ps").chmod(0o755)
    (fake / "sha256sum").write_text("#!/bin/sh\necho '0000  fake.hmm'\n")
    (fake / "sha256sum").chmod(0o755)

    env = {
        k: v for k, v in os.environ.items()
        if k not in ("MEM", "SEARCH_MEM_GB", "THREADS", "CAP", "OUT", "ROOT", "TREE", "RELEASE")
    }
    env.update({
        "PATH": f"{fake}:{os.environ['PATH']}",
        "CAPTURE_LOG": str(capture),
        "TREE": str(tree),
        "PROFILE": str(tmp_path / "profile.yaml"),
        "IMAGE": str(tmp_path / "bio-box.sif"),
        "DIGEST": "sha256:" + "0" * 64,
        "ROOT": str(root),
        "RELEASE": "r26.11.3",
        "OUT": str(out),
        "PFAM": str(pfam),
        "BIND": str(tmp_path),
        "CAP": "80",
    })
    return tree / "examples" / example / "run.sh", capture, out, env


@pytest.mark.parametrize("example", ["crispr_cas_neighbourhoods", "rnap"])
@pytest.mark.parametrize("mem,override,expected_mem_gb", [
    ("64G", None, "62"),
    ("16G", None, "14"),
    ("64G", "50", "50"),
])
def test_example_run_passes_derived_mem_gb_to_search_hmm(tmp_path, example, mem, override, expected_mem_gb):
    run_sh, capture, _, env = _make_example_stub_world(tmp_path, example)
    env["MEM"] = mem
    env["THREADS"] = "24"
    if override is not None:
        env["SEARCH_MEM_GB"] = override

    # Stops with rc=2 after the stubbed search hmm exits 99; we inspect the captured search invocation.
    subprocess.run(["bash", str(run_sh)], env=env, capture_output=True, text=True, timeout=60)
    search_lines = [ln for ln in capture.read_text().splitlines() if "cmd=search hmm " in ln]
    assert search_lines, f"no `search hmm` invocation captured for {example}: {capture.read_text()!r}"
    for line in search_lines:
        assert f"mem={mem} threads=24 " in line, line
        assert f"--mem-gb {expected_mem_gb} " in line, line


def test_rlimits_with_derived_62gb_vs_hardcoded_14gb_under_64g_cgroup():
    cap_64gib = 64 * 1024**3
    old_rl = runner.rlimits(runner.Limits(mem_gb=14, cpus=24), "inherited", "Linux", cap_bytes=cap_64gib)
    new_rl = runner.rlimits(runner.Limits(mem_gb=62, cpus=24), "inherited", "Linux", cap_bytes=cap_64gib)
    import resource
    assert old_rl[resource.RLIMIT_AS] == (14 * runner.GB, 14 * runner.GB)
    assert new_rl[resource.RLIMIT_AS] == (62 * runner.GB, 62 * runner.GB)


def test_release_train_resume_clears_incomplete_examples_dir_before_rerun(tmp_path):
    """When release_examples.sh fails leaving a non-empty $WORK/examples, resume must clear it before rerunning."""
    from tests.r1 import test_a20_train_post_manifest_resume as a20

    w = a20.TrainWorld(tmp_path)
    # Stub release_examples.sh to create output on first (failing) run and refuse non-empty OUT on rerun,
    # mirroring the real `scripts/release_examples.sh` (`die "OUT $OUT is not empty"`).
    (w.tree / "scripts" / "release_examples.sh").write_text(
        "#!/bin/bash\n"
        "echo \"release_examples $* OUT=$OUT\" >> \"$FAKE_LOG\"\n"
        "if [ -e \"$OUT\" ] && [ -n \"$(ls -A \"$OUT\")\" ]; then\n"
        "  echo \"release_examples.sh: OUT $OUT is not empty\" >&2\n"
        "  exit 2\n"
        "fi\n"
        "mkdir -p \"$OUT/crispr_cas_neighbourhoods/logs\"\n"
        "echo partial > \"$OUT/crispr_cas_neighbourhoods/logs/search-tniq.log\"\n"
        "if [ -n \"${FAKE_FAIL:-}\" ] && echo \"release_examples\" | grep -Eq -- \"$FAKE_FAIL\"; then\n"
        "  exit 1\n"
        "fi\n"
        "exit 0\n"
    )

    got1 = w.run("run", "--no-live-cpu", FAKE_FAIL="release_examples")
    assert got1.returncode != 0
    assert (w.work / "examples" / "crispr_cas_neighbourhoods" / "logs" / "search-tniq.log").is_file()
    assert not (w.work / "DONE").exists()

    (w.work / "before" / "verify-r26.11.1.log").write_text(a20.FIRST_VERIFY)
    (w.work / "estimate" / "before.json").write_text(a20.FIRST_ESTIMATE)
    w.publish_release("r26.11.2", n_records=2)
    w.set_plan(a20.PLAN_AFTER_R2_NOOP)
    w.log.write_text("")

    got2 = w.run("run", "--no-live-cpu")
    assert got2.returncode == 0, got2.stdout + got2.stderr
    assert (w.work / "DONE").read_text().split()[0] == "r26.11.2"

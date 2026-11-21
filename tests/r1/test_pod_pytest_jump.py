"""Card A-09: `BIOBOX_POD_JUMP` puts `-J <user@host>` on every ssh call of the pod scripts.

Stub `ssh` and `rsync` executables come first on PATH and record their argument lists; nothing touches a network.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-09")

ROOT = Path(__file__).resolve().parents[2]
JUMP = "if2453@wilkcomp03.ccbb.utexas.edu"

STUB = """#!/bin/sh
printf '%s %s\\n' "$(basename "$0")" "$*" >> "$STUB_LOG"
case "$(basename "$0")" in
  ssh) echo "ladder-B=inactive/success" ;;
esac
exit 0
"""


@pytest.fixture
def stubs(tmp_path: Path) -> dict[str, str]:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for name in ("ssh", "rsync"):
        p = bindir / name
        p.write_text(STUB)
        p.chmod(0o755)
    tree = tmp_path / "tree"
    tree.mkdir()
    env = dict(os.environ)
    env.pop("BIOBOX_POD_JUMP", None)
    env["PATH"] = f"{bindir}{os.pathsep}{env['PATH']}"
    env["STUB_LOG"] = str(tmp_path / "log")
    env["BIOBOX_TREE"] = str(tree)
    return env


def calls(env: dict[str, str]) -> list[str]:
    log = Path(env["STUB_LOG"])
    return log.read_text().splitlines() if log.exists() else []


def scripts() -> list[list[str]]:
    return [
        ["sh", str(ROOT / "scripts" / "pod_pytest.sh"), "T-00", "tests/x.py"],
        ["sh", str(ROOT / "scripts" / "orchestrator" / "watch_until.sh"), "true", "true", "1"],
        ["sh", str(ROOT / "scripts" / "orchestrator" / "watch_ladder.sh"), "1"],
    ]


@pytest.mark.parametrize("idx", [0, 1, 2])
def test_jump_on_every_call(stubs: dict[str, str], idx: int) -> None:
    stubs["BIOBOX_POD_JUMP"] = JUMP
    r = subprocess.run(scripts()[idx], env=stubs, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
    got = calls(stubs)
    assert got
    for line in got:
        if line.startswith("ssh "):
            assert f" -J {JUMP} " in line + " ", line
        else:
            assert f"ssh -o BatchMode=yes -J {JUMP}" in line, line


@pytest.mark.parametrize("idx", [0, 1, 2])
def test_no_jump_when_unset_or_empty(stubs: dict[str, str], idx: int) -> None:
    for value in (None, ""):
        Path(stubs["STUB_LOG"]).unlink(missing_ok=True)
        if value is not None:
            stubs["BIOBOX_POD_JUMP"] = value
        r = subprocess.run(scripts()[idx], env=stubs, capture_output=True, text=True, timeout=60)
        assert r.returncode == 0, r.stderr
        got = calls(stubs)
        assert got
        assert not any("-J" in line for line in got)


@pytest.mark.parametrize("bad", ["a b@host", "u@h;rm -rf x", "u@h|x", "u@h$(x)", "u@h&", "u@h`x`", "u@h'x", "u@h\ty"])
@pytest.mark.parametrize("idx", [0, 1, 2])
def test_bad_value_refused(stubs: dict[str, str], idx: int, bad: str) -> None:
    stubs["BIOBOX_POD_JUMP"] = bad
    r = subprocess.run(scripts()[idx], env=stubs, capture_output=True, text=True, timeout=60)
    assert r.returncode == 2
    assert calls(stubs) == []

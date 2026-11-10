"""Ring R0: RUNBOOK.md names only `bio-box` subcommands and options that the parsers have (card C-06).

Static only: nothing in the runbook is run. Every command in a fenced `bash` block is read:

- a `bio-box ...` command is parsed by the real parser of `biobox.cli` (the top-level parser is
  captured while `cli.main` builds it; the command is never dispatched);
- a `scripts/pod_run.sh [options] UNIT bio-box-args` command is checked option by option against
  the options the script parses, then its `bio-box` arguments the same way;
- a `scripts/pod_pytest.sh NAME [options] ...` command is checked against the script's options;
- `python -m biobox.<module> ...` (or `biopy -m ...`) must name a module that exists, and
  `biobox.doctor` is parsed by its own parser.

Shell variables (`$ROOT`, `"${PROFILE}"`) stand for a value; the parser sees a placeholder.
Lines inside a here-document are not commands, and neither is a function or array definition. A
fenced block of another kind (`text`, `python`, `yaml`) is not read. Section 10's limits table must
also carry every key and value of `profiles/wilkcomp02.yaml`.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import re
import shlex
from pathlib import Path

import pytest
import yaml

from biobox import cli, doctor

pytestmark = [pytest.mark.r0, pytest.mark.covers("C-06")]

ROOT = Path(__file__).resolve().parents[2]
RUNBOOK = ROOT / "RUNBOOK.md"
POD_RUN = ROOT / "scripts" / "pod_run.sh"
POD_PYTEST = ROOT / "scripts" / "pod_pytest.sh"
PROFILE = ROOT / "profiles" / "wilkcomp02.yaml"
FENCE = re.compile(r"^```(\w*)\s*$")
VAR = re.compile(r"\$\{?[A-Za-z_][A-Za-z0-9_]*\}?")
HEREDOC = re.compile(r"<<-?\s*['\"]?(\w+)['\"]?")
PLACEHOLDER = "X"


class _Parsed(Exception):
    """Raised by the captured parser once it has parsed, so that nothing is dispatched."""


def parse_with(main, argv: list[str]) -> str | None:
    """Run `main(argv)` up to its `parse_args` call; return None when the parser accepts `argv`,
    else the parser's error message."""
    real = argparse.ArgumentParser.parse_args
    errors: list[str] = []

    def capture(self, args=None, namespace=None):
        try:
            real(self, args, namespace)
        except SystemExit:
            pass
        raise _Parsed

    def fail(self, message):
        errors.append(message)
        raise SystemExit(2)

    argparse.ArgumentParser.parse_args = capture
    real_error = argparse.ArgumentParser.error
    argparse.ArgumentParser.error = fail
    try:
        main(argv)
    except (_Parsed, SystemExit):
        pass
    finally:
        argparse.ArgumentParser.parse_args = real
        argparse.ArgumentParser.error = real_error
    return errors[0] if errors else None


def commands(text: str) -> list[tuple[int, str]]:
    """(line number, command) of every command line in the fenced `bash` blocks: a trailing
    backslash joins the next line; comments, blank lines and here-document bodies are dropped."""
    out: list[tuple[int, str]] = []
    kind: str | None = None
    pending, start, doc_end = "", 0, None
    for n, raw in enumerate(text.splitlines(), 1):
        m = FENCE.match(raw)
        if m and kind is None:
            kind = m.group(1) or "text"
            continue
        if m and kind is not None:
            kind, pending, doc_end = None, "", None
            continue
        if kind != "bash":
            continue
        if doc_end is not None:
            if raw.strip() == doc_end:
                doc_end = None
            continue
        line = raw.strip()
        if not pending and (not line or line.startswith("#")):
            continue
        start = start if pending else n
        if line.endswith("\\"):
            pending += line[:-1] + " "
            continue
        cmd = (pending + line).strip()
        pending = ""
        out.append((start, cmd))
        h = HEREDOC.search(cmd)
        if h:
            doc_end = h.group(1)
    return out


def tokens(cmd: str) -> list[str]:
    cmd = re.sub(r"\s#\s.*$", "", cmd)  # a trailing comment
    cmd = re.sub(r"^\$ ", "", cmd)
    return shlex.split(VAR.sub(PLACEHOLDER, cmd))


def options_of(script: Path) -> dict[str, bool]:
    """The `--option` branches of a script's `case` loop: {option: takes a value}."""
    found = {}
    for m in re.finditer(r"^\s*(--[a-z-]+)\)\s*(.*)$", script.read_text(), re.M):
        found[m.group(1)] = "shift 2" in m.group(2) or "+=(" in m.group(2)
    return found


def check(cmd: str) -> list[str]:
    """What is wrong with one command line (empty: nothing found)."""
    if re.match(r"^[\w-]+\(\)\s*\{", cmd) or re.match(r"^\w+=\(", cmd):  # a function or array definition
        return []
    t = tokens(cmd)
    problems: list[str] = []
    for i, tok in enumerate(t):
        if tok.endswith("pod_run.sh"):
            return problems + check_pod_run(t[i + 1:])
        if tok.endswith("pod_pytest.sh"):
            return problems + check_script_options(POD_PYTEST, t[i + 1:], positional=1)
        if tok == "bio-box":
            return problems + check_bio_box(t[i + 1:])
        if tok == "-m" and i + 1 < len(t) and t[i + 1].startswith("biobox."):
            return problems + check_module(t[i + 1], t[i + 2:])
    return problems


def check_bio_box(args: list[str]) -> list[str]:
    if not args:
        return ["bio-box without a subcommand"]
    if args[0] == "run-step" and "--" not in args:
        return ["bio-box run-step without `-- CMD`"]
    msg = parse_with(cli.main, list(args))
    return [] if msg is None else [f"bio-box {' '.join(args)}: {msg}"]


def check_script_options(script: Path, args: list[str], positional: int) -> list[str]:
    known, problems, i, seen = options_of(script), [], 0, 0
    while i < len(args):
        a = args[i]
        if a.startswith("--") and seen >= positional:
            if a not in known:
                problems.append(f"{script.name}: unknown option {a}")
                break
            i += 2 if known[a] else 1
        else:
            seen += 1
            i += 1
    return problems


def check_pod_run(args: list[str]) -> list[str]:
    known, i = options_of(POD_RUN), 0
    while i < len(args) and args[i].startswith("--"):
        if args[i] == "--":
            i += 1
            break
        if args[i] not in known:
            return [f"pod_run.sh: unknown option {args[i]}"]
        i += 2 if known[args[i]] else 1
    if i >= len(args):
        return ["pod_run.sh: no unit name"]
    return check_bio_box(args[i + 1:])  # args[i] is the unit


def check_module(module: str, args: list[str]) -> list[str]:
    if importlib.util.find_spec(module) is None:
        return [f"python -m {module}: no such module"]
    if module == "biobox.doctor":
        msg = parse_with(doctor.main, list(args))
        return [] if msg is None else [f"python -m {module} {' '.join(args)}: {msg}"]
    if module == "biobox.cli":
        return check_bio_box(args)
    return []


def read_commands() -> list[tuple[int, str]]:
    assert RUNBOOK.is_file(), f"{RUNBOOK} does not exist"
    return commands(RUNBOOK.read_text())


# --- the checker itself: it must reject what it is meant to reject ---------------------------------


@pytest.mark.parametrize("cmd", [
    "bio-box no-such-command --root R",
    "bio-box build-stage --profile P --root R --stage S --release X --genomes g --objects o --no-such-option",
    "bio-box build-stage --profile P --root R --stage S --release X --genomes g",  # --objects is required
    "bio-box search status",  # the run id is missing
    "bio-box run-step --profile P --root R --stage S --step T",  # no `-- CMD`
    "bash scripts/pod_run.sh --no-such-option c06-x build-stage --profile P",
    "bash scripts/pod_run.sh c06-x release nosuchstep --root R --release X",
    "biopy -m biobox.doctor --no-such-option",
    "biopy -m biobox.no_such_module",
    "sh scripts/pod_pytest.sh C-06 --no-such-option",
])
def test_the_checker_rejects_a_wrong_command(cmd):
    assert check(cmd), f"accepted: {cmd}"


@pytest.mark.parametrize("cmd", [
    "bio-box search status $RUN --root \"$ROOT\"",
    "bio-box release manifest --root R --release X",
    "bash scripts/pod_run.sh --profile \"$PROFILE\" --image \"$IMAGE\" c06-x export proteins --root R --release X",
    "biopy -m biobox.doctor --profile \"$PROFILE\"",
    "sh scripts/pod_pytest.sh C-06 --shards 8",
])
def test_the_checker_accepts_a_right_command(cmd):
    assert check(cmd) == []


def test_a_heredoc_body_is_not_a_command():
    text = "```bash\ncat > f <<'EOF'\nbio-box nonsense\nEOF\nbio-box pack verify --root R --release X\n```\n"
    assert [c for _, c in commands(text)] == ["cat > f <<'EOF'", "bio-box pack verify --root R --release X"]


# --- the runbook ------------------------------------------------------------------------------------


def test_the_runbook_has_commands():
    assert len(read_commands()) >= 20


def test_every_bio_box_command_in_the_runbook_is_one_the_parser_accepts():
    bad = [f"RUNBOOK.md:{n}: {p}" for n, cmd in read_commands() for p in check(cmd)]
    assert not bad, "\n".join(bad)


def test_the_limits_table_carries_every_value_of_the_host_profile():
    """Section 10 is copied from the profile; a changed profile value must show up in the table."""
    text = RUNBOOK.read_text()
    start = text.index("## 10. Limits")
    section = text[start:text.index("\n## ", start + 1)]
    lines = section.splitlines()
    missing = []
    for key, value in yaml.safe_load(PROFILE.read_text()).items():
        shown = str(value).lower() if isinstance(value, bool) else str(value)
        rows = [ln for ln in lines if f"`{key}`" in ln]
        if not rows:
            missing.append(f"{key}: no row names it")
        elif not any(shown in ln for ln in rows):
            missing.append(f"{key}: no row that names it shows {shown!r}")
    assert not missing, "\n".join(missing)


def test_every_script_the_runbook_names_exists():
    named = set(re.findall(r"scripts/[A-Za-z0-9_./-]+\.(?:sh|py)", RUNBOOK.read_text()))
    assert named, "the runbook names no script"
    missing = sorted(s for s in named if not (ROOT / s).is_file())
    assert not missing, f"RUNBOOK.md names scripts that do not exist: {missing}"

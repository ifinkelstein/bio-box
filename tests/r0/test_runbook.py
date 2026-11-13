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

Card C-29: the pod's test scripts and the documents name the image built from `main` and no dependency
directory; `scripts/pod_pytest.sh` is run with `ssh` and `rsync` replaced by stubs, so nothing leaves the host.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import os
import re
import shlex
import subprocess
from pathlib import Path

import pytest
import yaml

from biobox import cli, doctor

pytestmark = [pytest.mark.r0, pytest.mark.covers("C-06")]

ROOT = Path(__file__).resolve().parents[2]
RUNBOOK = ROOT / "RUNBOOK.md"
AGENT_GUIDE = ROOT / "AGENT_GUIDE.md"  # for agents that use releases; its commands are checked the same way
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


def read_commands(path: Path = RUNBOOK) -> list[tuple[int, str]]:
    assert path.is_file(), f"{path} does not exist"
    return commands(path.read_text())


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


@pytest.mark.parametrize("doc", [RUNBOOK, AGENT_GUIDE], ids=lambda p: p.name)
def test_every_bio_box_command_in_the_runbook_is_one_the_parser_accepts(doc):
    bad = [f"{doc.name}:{n}: {p}" for n, cmd in read_commands(doc) for p in check(cmd)]
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


def section_of(text: str, heading: str) -> str:
    start = text.index(heading)
    end = text.find("\n## ", start + 1)
    return text[start:] if end < 0 else text[start:end]


@pytest.mark.covers("C-11")
def test_the_runbook_shows_the_scan_options_of_build_stage():
    """Card C-11: the build section shows a `build-stage` command with `--scan-jobs` and `--scan-cpus` (the
    parser check above reads it), and the limits and known-limits sections say what bounds them."""
    text = RUNBOOK.read_text()
    build = [c for _, c in commands(section_of(text, "## 2. ")) if "build-stage" in c]
    assert any("--scan-jobs" in c and "--scan-cpus" in c for c in build), "no build-stage command with both options"
    for heading in ("## 10. Limits", "## 11. Known limits"):
        part = section_of(text, heading)
        assert "--scan-jobs" in part and "--scan-cpus" in part, f"{heading} does not name both options"
    appendix = section_of(text, "## Appendix A")
    rows = [ln for ln in appendix.splitlines() if "--scan-jobs" in ln]
    assert rows and all("tested on the fixture" in ln for ln in rows)


@pytest.mark.covers("C-12")
def test_the_unit_quota_is_the_night_cap_and_the_runbook_says_build_stage_holds_the_hour():
    """Card C-12: `pod_run.sh` sets the unit's CPU quota from the profile's night cap, the most any hour allows, and
    `build-stage` holds the limit of the hour itself, checked again before each step and each shard launch. The
    runbook says so where it names the quota, and no longer lists the allowance decided once as a known limit."""
    script = POD_RUN.read_text()
    assert re.search(r"^THREADS=\$\{THREADS:-\$\(prof max_threads_night\)\}$", script, re.M)
    header = script.split("set -euo pipefail", 1)[0]
    assert "max_threads_night" in header and "build-stage" in header
    text = RUNBOOK.read_text()
    assert "decided once" not in section_of(text, "## 11. Known limits")
    for heading in ("## 2. ", "## 10. Limits"):
        part = section_of(text, heading)
        assert "max_threads_night" in part and "before each shard launch" in part, heading


@pytest.mark.covers("C-28")
def test_c28_the_runbook_states_the_night_rule_as_implemented():
    """Card C-28: the night cap holds while no other user computes, read as other users' CPU over a window, not as a
    load average above 8; the allowance of a launch subtracts ours elsewhere; section 1 shows how to read it."""
    from biobox import runstep

    text = RUNBOOK.read_text()
    for heading in ("## 1. ", "## 2. ", "## 10. Limits"):
        part = section_of(text, heading)
        assert "load average above 8" not in part and "load average is at most 8" not in part, heading
    limits = section_of(text, "## 10. Limits")
    row = next(ln for ln in limits.splitlines() if ln.startswith('| "Other compute"'))
    assert f"{runstep.OTHERS_BUSY_CPUS:g} of a core" in row and f"{runstep.CPU_WINDOW_S:g} s window" in row
    assert f"uid above {runstep.SYSTEM_UID_MAX}" in row and "OTHERS_BUSY_CPUS" in row and "BUSY_LOAD" not in limits
    assert any(ln.startswith("| Allowance of a launch") and "ours elsewhere" in ln for ln in limits.splitlines())
    build = section_of(text, "## 2. ")
    assert "ours elsewhere" in build and "06:45" in build and "allowance ..." in build
    first = [c for _, c in commands(section_of(text, "## 1. "))]
    assert any(c.startswith("biopy -m biobox.runstep") for c in first), first


@pytest.mark.covers("C-13")
def test_c13_the_verify_section_runs_release_verify_with_and_without_a_previous_release():
    """Card C-13 (4): section 6 verifies a release with `bio-box release verify` (the parser check above reads the
    commands), once alone and once against the previous release."""
    cmds = [c for _, c in commands(section_of(RUNBOOK.read_text(), "## 6. "))]
    verify = [c for c in cmds if "bio-box release verify" in c]
    assert any("--previous" not in c for c in verify) and any("--previous" in c for c in verify), cmds
    assert not any("verify_manifest" in c for c in cmds)  # no python one-liner in its place


@pytest.mark.covers("C-13")
def test_c13_the_runbook_states_the_cost_of_a_stop_inside_a_step_and_the_logs_are_a_timeline():
    """Card C-13 (2), (3), (7): a `STOP` inside a step discards that step's finished work, with C-10's measured
    cost; `run-step` and `build-stage` lines are timestamped as they happen; a stopped unit's log still ends with the
    `pod_run: ... exited` line."""
    text = RUNBOOK.read_text()
    stop = section_of(text, "## 3. ")
    assert "57.7 s" in stop and "4.98 CPU-minutes" in stop and "measured" in stop
    build = section_of(text, "## 2. ")
    assert "buffered until it exits" not in build and "timestamp" in build
    assert "leaves no `pod_run: ... exited` line" not in stop
    assert "no `pod_run: ... exited` line" not in section_of(text, "## 11. ")


def test_every_script_the_runbook_names_exists():
    named = set(re.findall(r"scripts/[A-Za-z0-9_./-]+\.(?:sh|py)", RUNBOOK.read_text()))
    assert named, "the runbook names no script"
    missing = sorted(s for s in named if not (ROOT / s).is_file())
    assert not missing, f"RUNBOOK.md names scripts that do not exist: {missing}"


# --- W-13 (amended after W-03): the search commands name the directories outside the data root ----------------


def bio_box_args(cmd: str) -> list[str] | None:
    """The `bio-box` arguments of a command line (after `bio-box`, or after the unit of `pod_run.sh`), or None."""
    if re.match(r"^[\w-]+\(\)\s*\{", cmd) or re.match(r"^\w+=\(", cmd):
        return None
    t = tokens(cmd)
    for i, tok in enumerate(t):
        if tok.endswith("pod_run.sh"):
            known, j = options_of(POD_RUN), i + 1
            while j < len(t) and t[j].startswith("--"):
                j += 2 if known.get(t[j]) else 1
            return t[j + 1:]
        if tok == "bio-box":
            return t[i + 1:]
    return None


# What each search command needs when it has no --profile (card W-03): one option of each group.
SEARCH_DIRS = {
    ("pack", "build"): [("--pack-dir",)],
    ("pack", "verify"): [("--pack-dir",)],
    ("search", "run"): [("--pack-dir",), ("--runs", "--out"), ("--work",)],
    ("search", "status"): [("--runs",)],
    ("search", "cancel"): [("--runs",)],
    ("search", "hits"): [("--runs",), ("--pack-dir",)],
    ("search", "hmm"): [("--pack-dir",), ("--runs",), ("--work",)],  # card W-07
    ("nbr", None): [("--runs",), ("--pack-dir",)],
}


def search_dir_problems(args: list[str]) -> list[str]:
    key = (args[0], args[1] if args[0] != "nbr" else None) if args else None
    if key not in SEARCH_DIRS or "--profile" in args:
        return []
    return [f"{' '.join(args[:2])}: no --profile and none of {group}"
            for group in SEARCH_DIRS[key] if not any(o in args for o in group)]


@pytest.mark.covers("W-13", "W-03")
def test_w13_the_checker_of_search_directories_rejects_a_command_with_neither():
    assert search_dir_problems(["search", "status", "X", "--root", "X"])
    assert search_dir_problems(["pack", "build", "--root", "X", "--release", "X"])
    assert search_dir_problems(["nbr", "summarize", "X", "--root", "X", "--runs", "X"])  # no --pack-dir
    assert not search_dir_problems(["search", "status", "X", "--root", "X", "--profile", "X"])
    assert not search_dir_problems(["search", "cancel", "X", "--root", "X", "--runs", "X"])


@pytest.mark.covers("W-13", "W-03")
def test_w13_every_pack_search_and_nbr_command_names_its_directories_or_the_profile():
    """Card W-13 (c): since W-03 the pack and the runs live outside the data root, and `pod_run.sh` does not export
    BIOBOX_PROFILE, so every pack, search and nbr command of the runbook passes --profile or the directories."""
    found, bad = 0, []
    for n, cmd in read_commands():
        args = bio_box_args(cmd)
        if args and args[0] in ("pack", "search", "nbr"):
            found += 1
            bad += [f"RUNBOOK.md:{n}: {p}" for p in search_dir_problems(args)]
    assert found >= 10, found
    assert not bad, "\n".join(bad)


@pytest.mark.covers("W-13")
def test_w13_the_search_section_watches_and_cancels_a_wait_run_with_bio_box():
    """Card W-13: a `--wait` run has a pending record, so section 7 watches it with `search status` and stops it with
    `search cancel`; the runbook no longer says status finds a run only once it has ended, nor that the default
    work directory is /tmp/atd/search, and it says what `search run` without --wait does in a unit."""
    text = RUNBOOK.read_text()
    search = section_of(text, "## 7. ")
    cmds = [bio_box_args(c) or [] for _, c in commands(search)]
    runs = [a for a in cmds if a[:2] == ["search", "run"]]
    assert runs and all("--wait" in a for a in runs)
    assert any(a[:2] == ["search", "status"] for a in cmds) and any(a[:2] == ["search", "cancel"] for a in cmds)
    assert "finds a run only when it has ended" not in text
    assert "finds nothing until the run has ended" not in text
    assert "default is\n`/tmp/atd/search`" not in text and "the default is `/tmp/atd/search`" not in text
    assert "outside the data root" in search
    assert "without `--wait`" in search and "BIOBOX_POD_UNIT" in search
    clean = section_of(text, "## 9. ")  # what is ours under /tmp, and what may be cleaned
    assert "$BBTMP/packs" in clean and "$BBTMP/search/$USER" in clean
    known = section_of(text, "## 11. ")
    assert "runner exited without closing the run" not in known  # W-13 fixed both search notes
    appendix = section_of(text, "## Appendix A")
    assert "search cancel" in appendix and "w13-" in appendix


@pytest.mark.covers("W-07")
def test_w07_the_search_section_shows_the_one_command_with_its_positive_control_and_output_directory():
    """Card W-07: section 7 shows `search hmm` in a `pod_run.sh` unit with a positive control and `--out`, and says
    what it writes, that it waits for its run, what a missing pack does and what a release without a domain library
    gives."""
    search = section_of(RUNBOOK.read_text(), "## 7. ")
    hmm = [a for a in (bio_box_args(c) or [] for _, c in commands(search)) if a[:2] == ["search", "hmm"]]
    assert hmm and all({"--positive", "--out", "--hmm", "--bits", "--release"} <= set(a) for a in hmm), hmm
    assert any("pod_run.sh" in c for _, c in commands(search) if "search hmm" in c)
    for words in ("recipe.yaml", "hits.tsv", "hits.faa", "neighbors.tsv", "loci/", "no domain library",
                  "pack build", "foreground"):
        assert words in search, words


# --- C-14: errors found by the 10 Gbp rebuild ---------------------------------------------------------------


@pytest.mark.covers("C-14")
def test_c14_the_done_line_counts_12_steps_plus_one_per_library():
    """Card C-14 (1): `build-stage` prints `all <len(plan)> steps done`; the plan is 12 steps and one scan per
    `--library`, so 13 with one library and 12 with none. Section 2 says `--library` is optional and shows both."""
    build = section_of(RUNBOOK.read_text(), "## 2. ")
    assert "`--library` is optional" in build
    assert "all 12 steps done" in build and "all 13 steps done" in build
    for ln in build.splitlines():
        if "all 13 steps done" in ln:
            assert "one `--library`" in ln or "one library" in ln, ln


@pytest.mark.covers("C-14")
def test_c14_the_scan_recommendation_leaves_room_under_a_cap_that_is_a_total():
    """Card C-14 (2): J x C = 40 by day is the whole cap; the cap is the total over everything of ours on the pod, so
    section 2 gives J as a formula of what else runs, with a worked example, and not `--scan-jobs 10 --scan-cpus 4`."""
    text = RUNBOOK.read_text()
    build = section_of(text, "## 2. ")
    assert "--scan-jobs 10 --scan-cpus 4" not in text
    assert "total" in build and "everything of ours" in build
    assert "worked example" in build and "other threads" in build
    limits = section_of(text, "## 10. Limits")
    assert "total" in limits


@pytest.mark.covers("C-14")
def test_c14_the_stop_rule_for_running_shards_is_stated_once_and_both_sections_agree():
    """Card C-14 (3): a shard (a `run-step`) whose command ends while `STOP` exists publishes nothing and exits 75
    (runstep.py, after the command), and `build-stage` launches no more shards (stage.py `_launch`). So section 2 does
    not say shards at a `STOP` `finish and publish`; section 3 states the rule and section 2 points to it."""
    text = RUNBOOK.read_text()
    build, stop = section_of(text, "## 2. "), section_of(text, "## 3. ")
    assert "finish and publish" not in build and "finish and publish" not in stop
    assert "finish and publish" not in section_of(text, "## 11. ")
    assert "shard" in stop and "publishes nothing" in stop and "launches no more" in stop
    assert "(section 3)" in build.split("A shard that fails", 1)[1]
    assert stop.count("`runstep.py`") + stop.count("runstep.py") >= 1


@pytest.mark.covers("C-14")
def test_c14_no_command_block_uses_a_nested_quote_python_one_liner():
    """Card C-14 (4): `python -c '...'` with quotes inside survives no `ssh '...' bash -c "..."`; the runbook uses
    script files, which take their arguments without quoting."""
    bad = [f"RUNBOOK.md:{n}: {c[:60]}" for n, c in read_commands() if re.search(r"(?:python|biopy)\s+-c\b", c)]
    assert not bad, "\n".join(bad)


# --- GUIDE.md, the user guide: its commands are checked as the runbook's are --------------------------------


@pytest.mark.parametrize("doc", [ROOT / "GUIDE.md"], ids=lambda p: p.name)
def test_every_bio_box_command_in_the_user_guide_is_one_the_parser_accepts(doc):
    """GUIDE.md prints commands for a scientist to copy: each is read and checked as the runbook's are, and every pack,
    search and nbr command names its directories or the profile (W-03)."""
    cmds = commands(doc.read_text())
    assert cmds, f"{doc.name} has no command in a bash block"
    bad = [f"{doc.name}:{n}: {p}" for n, cmd in cmds for p in check(cmd)]
    for n, cmd in cmds:
        args = bio_box_args(cmd)
        if args and args[0] in ("pack", "search", "nbr"):
            bad += [f"{doc.name}:{n}: {p}" for p in search_dir_problems(args)]
    assert not bad, "\n".join(bad)


# --- W-22: one place for a release's search pack ------------------------------------------------------------------


def cli_args(cmd: str) -> list[str] | None:
    """The `bio-box` arguments of a command line, also when it calls `python -m biobox.cli` (the agent guide's form)."""
    args = bio_box_args(cmd)
    if args is not None:
        return args
    t = tokens(cmd)
    for i, tok in enumerate(t):
        if tok == "-m" and i + 1 < len(t) and t[i + 1] == "biobox.cli":
            return t[i + 2:]
    return None


def profile_pack_dir() -> Path:
    """The pack directory `bio-box` takes from the shipped host profile when given only `--profile`."""
    return cli._packs(argparse.Namespace(pack_dir=None, profile=str(PROFILE)))


@pytest.mark.covers("W-22")
@pytest.mark.parametrize("doc", [ROOT / "GUIDE.md", AGENT_GUIDE, RUNBOOK], ids=lambda p: p.name)
def test_w22_the_documents_find_the_pack_through_the_profile(doc):
    """Card W-22: the host profile names the host's one pack directory, so no command of the guides or the runbook
    that passes `--profile` also names a pack directory, and each document says where the profile's is."""
    assert profile_pack_dir() == Path("/tmp/if2453/bio-box/packs")
    bad = []
    for n, cmd in commands(doc.read_text()):
        args = cli_args(cmd)
        if args and "--profile" in args and "--pack-dir" in args:
            bad.append(f"{doc.name}:{n}: {' '.join(args[:2])} names --pack-dir beside --profile")
    assert not bad, "\n".join(bad)
    assert str(profile_pack_dir()) in doc.read_text(), f"{doc.name} does not say where the profile's pack directory is"


@pytest.mark.covers("W-22")
def test_w22_the_runbook_puts_the_scratch_of_a_pack_verify_under_the_profiles_tmp():
    """Card W-22: `pack verify` and a no-op `pack build` write their scratch under the profile's tmp, so the cleaning
    section lists `.verify-*` there and no longer in the pack directory, and section 7 says a pack directory may be
    read-only for all but its builder. The agent guide no longer warns that both write scratch into it."""
    text = RUNBOOK.read_text()
    clean = section_of(text, "## 9. ")
    assert "$BBTMP/search/$USER/work/.verify-" in clean
    assert not re.search(r"\$BBTMP/packs/release=<R>/[^\n]*\n?[^\n]*\.verify-", clean)
    assert "read-only" in section_of(text, "## 7. ")
    assert "both write scratch there" not in AGENT_GUIDE.read_text()


# --- C-29: tests and builds run in the image built from `main`, with no dependency directory ----------------------

DETACHED = ROOT / "scripts" / "orchestrator" / "detached_suite.sh"
GUIDE = ROOT / "GUIDE.md"
POD_IMAGES = "/stor/home/if2453/bio-box/images"
TEST_SIF = "bio-box-test_64f1bf0.sif"  # ghcr.io/ifinkelstein/bio-box-test@sha256:cb7cd0f1...
RUNTIME_SIF = "bio-box_64f1bf0.sif"  # ghcr.io/ifinkelstein/bio-box@sha256:2faf4b8d...
# sha256 of the .sif files on wilkcomp02 (measured, `sha256sum`, 2026-11-13): what `pod_run.sh --digest` compares
TEST_SIF_SHA256 = "3d77df6b7f0b894b5371e42eaf077dfda09ebc97b1f78f2ed1d2cfc5fa9fbda4"
RUNTIME_SIF_SHA256 = "cb7353f85316981be5268c2212100322b5fbf1086c32b07e2f9b47802b6afba3"


def run_pod_pytest(tmp_path: Path, env_extra: dict[str, str]) -> str:
    """Run scripts/pod_pytest.sh with `ssh` and `rsync` replaced by stubs; return the remote command of its last
    `ssh` call (the pytest launch). Nothing leaves the host."""
    stubs = tmp_path / "bin"
    stubs.mkdir(parents=True)
    log = tmp_path / "ssh.log"
    (stubs / "ssh").write_text(f'#!/bin/sh\nfor a in "$@"; do last=$a; done\nprintf "%s\\n" "$last" >> "{log}"\n')
    (stubs / "rsync").write_text("#!/bin/sh\nexit 0\n")
    for stub in stubs.iterdir():
        stub.chmod(0o755)
    tree = tmp_path / "tree"
    tree.mkdir()
    env = {"PATH": f"{stubs}:{os.environ['PATH']}", "BIOBOX_TREE": str(tree), "BIOBOX_POD_CI": "/ci", **env_extra}
    subprocess.run(["sh", str(POD_PYTEST), "T", "tests/r0/test_x.py"], env=env, check=True, capture_output=True)
    return log.read_text().splitlines()[-1]


@pytest.mark.covers("C-29")
def test_c29_pod_pytest_runs_the_new_image_with_only_the_tree_on_the_path(tmp_path):
    """Card C-29 (3): by default the suite runs in the test image built from `main` at 64f1bf0, whose environment is
    the locked one, with the tree's source alone on PYTHONPATH (a dependency directory first on the path would shadow
    the locked packages). The binds and the base temp of before stay."""
    remote = run_pod_pytest(tmp_path, {})
    assert f"'{POD_IMAGES}/{TEST_SIF}'" in remote, remote
    assert re.search(r"--env PYTHONPATH=/ci/T/src\s", remote), remote
    assert "pydeps" not in remote and "74b3e05" not in remote
    assert "--bind /tmp/if2453 " in remote and "--bind /sys/fs/cgroup:/sys/fs/cgroup:ro" in remote
    assert "--basetemp=/tmp/if2453/ci-tmp/T-0" in remote and "MemoryMax=32G" in remote


@pytest.mark.covers("C-29")
def test_c29_pod_pytest_keeps_the_image_and_dependency_overrides(tmp_path):
    remote = run_pod_pytest(tmp_path, {"BIOBOX_POD_IMAGE": "/x/old.sif", "BIOBOX_POD_PYDEPS": "/a:/b"})
    assert "'/x/old.sif'" in remote and TEST_SIF not in remote
    assert re.search(r"--env PYTHONPATH=/ci/T/src:/a:/b\s", remote), remote


@pytest.mark.covers("C-29")
def test_c29_the_detached_suite_defaults_to_the_new_image_and_no_dependency_directory():
    """Card C-29 (3): detached_suite.sh, the suite as detached units, takes the same defaults and overrides as
    pod_pytest.sh, adds a dependency directory only when one is given, passes the cross-host peer file on, and can run
    the image's own suite (`--own`: the source and tests built into the image, no PYTHONPATH)."""
    text = DETACHED.read_text()
    assert re.search(rf"^IMAGE=\$\{{BIOBOX_POD_IMAGE:-{POD_IMAGES}/{TEST_SIF}\}}$", text, re.M), "default image"
    assert re.search(r"^PYDEPS=\$\{BIOBOX_POD_PYDEPS:-\}$", text, re.M), "no default dependency directory"
    assert "pilot_s1" not in text and "74b3e05" not in text
    assert "${PYDEPS:+:$PYDEPS}" in text
    assert "--bind /tmp/if2453 --bind /sys/fs/cgroup:/sys/fs/cgroup:ro" in text
    assert '"--basetemp=/tmp/if2453/ci-tmp/$NAME-$k"' in text
    assert '"${1:-}" = --own' in text and "DST=/opt/bio-box" in text and "BIOBOX_PEER_MANIFEST=$" in text


@pytest.mark.covers("C-29")
def test_c29_the_runbook_sets_up_the_new_image_without_a_dependency_directory():
    """Card C-29 (3): section 0 names the runtime image built from `main` and its sha256, and binds only the tree's
    source; section 8 has the image check (the image's own suite) in the test image; the known-limits table no longer
    says the image is old; the runbook says how the next build takes the new image."""
    text = RUNBOOK.read_text()
    setup = "\n".join(c for _, c in commands(section_of(text, "## 0. ")))
    assert f'export IMAGE="$HOME/bio-box/images/{RUNTIME_SIF}"' in setup
    assert f"export DIGEST=sha256:{RUNTIME_SIF_SHA256}" in setup
    assert 'export PYTHONPATH_POD="$TREE/src"' in setup
    assert "pydeps" not in setup.lower() and "74b3e05" not in setup
    tests = section_of(text, "## 8. ")
    own = [c for _, c in commands(tests) if "detached_suite.sh" in c and "--own" in c]
    assert own, "section 8 has no image check"
    assert TEST_SIF in tests and TEST_SIF_SHA256 in tests and "1,720 passed" in tests
    known = section_of(text, "## 11. ")
    assert "The image on the pod is old" not in known
    assert "r50_driver.sh" in text and "`IMAGE` and `DIGEST`" in text and "`PYPATH`" in text


@pytest.mark.covers("C-29")
@pytest.mark.parametrize("doc", [GUIDE, AGENT_GUIDE, RUNBOOK], ids=lambda p: p.name)
def test_c29_no_document_puts_a_dependency_directory_on_the_path(doc):
    """Card C-29 (3): the image holds the locked environment, so no PYTHONPATH a document sets names a dependency
    directory (`pydeps`)."""
    bad = [f"{doc.name}:{n}: {line.strip()[:80]}" for n, line in enumerate(doc.read_text().splitlines(), 1)
           if "PYTHONPATH" in line and "pydeps" in line.lower()]
    assert not bad, "\n".join(bad)

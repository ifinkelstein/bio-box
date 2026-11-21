"""Card C-49: `scripts/orchestrator/relay_to_tower.sh` mirrors files from the pod to the NAS `tower` through the laptop.

Stub `ssh`, `scp` and `sha256sum` come first on PATH and record every call; nothing touches a network. The stub `ssh`
runs the remote command locally with `sh -c`, and two local directories stand for the pod's source directory and
`tower`'s destination directory, so the script's whole logic (the plan, the stream, the hash taken as the bytes arrive,
the size check, `mv -n`, `SHA256SUMS`) runs here. `STUB_CORRUPT=1` flips one byte of every stream on its way to
`tower`, as a fault in transit would.
"""

from __future__ import annotations

import hashlib
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-49")

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "orchestrator" / "relay_to_tower.sh"
JUMP = "if2453@wilkcomp03.ccbb.utexas.edu"

SSH = """#!{python}
import json, os, shutil, subprocess, sys
args = sys.argv[1:]
cfg, no_stdin, i = None, False, 0
while i < len(args) and args[i].startswith("-"):
    if args[i] in ("-F", "-o", "-J", "-i", "-p", "-l", "-S", "-E"):
        if args[i] == "-F":
            cfg = args[i + 1]
        i += 2
        continue
    no_stdin = no_stdin or args[i] == "-n"
    i += 1
host, cmd = args[i], " ".join(args[i + 1:])
with open(os.environ["STUB_LOG"], "a") as f:
    f.write(json.dumps({{"tool": "ssh", "host": host, "cmd": cmd}}) + "\\n")
if cfg and os.environ.get("STUB_CFG_LOG"):
    shutil.copyfile(cfg, os.environ["STUB_CFG_LOG"])
if host not in ("pod", "tower"):
    sys.exit(255)
if no_stdin:
    r = subprocess.run(["sh", "-c", cmd], stdin=subprocess.DEVNULL)
elif host == "tower" and os.environ.get("STUB_CORRUPT") and "tee" in cmd:
    data = sys.stdin.buffer.read()
    if data:
        data = bytes([data[0] ^ 1]) + data[1:]
    r = subprocess.run(["sh", "-c", cmd], input=data)
else:
    r = subprocess.run(["sh", "-c", cmd])
sys.exit(r.returncode)
"""

SHA256SUM = """#!{python}
import hashlib, json, os, sys
args = [a for a in sys.argv[1:] if a != "--"]
with open(os.environ["STUB_LOG"], "a") as f:
    f.write(json.dumps({{"tool": "sha256sum", "args": args}}) + "\\n")
rc = 0
for name in args or ["-"]:
    h = hashlib.sha256()
    try:
        src = sys.stdin.buffer if name == "-" else open(name, "rb")
    except OSError:
        sys.stderr.write("sha256sum: " + name + ": No such file or directory\\n")
        rc = 1
        continue
    for block in iter(lambda: src.read(1 << 16), b""):
        h.update(block)
    print(h.hexdigest() + "  " + name)
sys.exit(rc)
"""

SCP = """#!{python}
import json, os, sys
with open(os.environ["STUB_LOG"], "a") as f:
    f.write(json.dumps({{"tool": "scp", "args": sys.argv[1:]}}) + "\\n")
sys.exit(1)
"""


@pytest.fixture
def world(tmp_path: Path) -> dict:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for name, text in (("ssh", SSH), ("sha256sum", SHA256SUM), ("scp", SCP)):
        p = bindir / name
        p.write_text(text.format(python=sys.executable))
        p.chmod(0o755)
    pod = tmp_path / "pod" / "root_backups"
    (pod / "root_rM2").mkdir(parents=True)
    rnd = random.Random(49)
    files = {
        "root_rM2/l0seq.tar": bytes(rnd.getrandbits(8) for _ in range(300_000)),
        "root_rM2/l2.tar": bytes(rnd.getrandbits(8) for _ in range(170_000)),
        "added_r26_11_1.tar": bytes(rnd.getrandbits(8) for _ in range(90_001)),
    }
    for rel, data in files.items():
        (pod / rel).write_bytes(data)
    tower = tmp_path / "tower"
    tower.mkdir()
    home = tmp_path / "home"
    home.mkdir()
    env = dict(os.environ)
    for k in ("BIOBOX_POD", "BIOBOX_POD_JUMP", "BIOBOX_TOWER", "BIOBOX_TOWER_IDENTITY", "STUB_CORRUPT"):
        env.pop(k, None)
    env.update({
        "PATH": f"{bindir}{os.pathsep}{env['PATH']}",
        "HOME": str(home),
        "STUB_LOG": str(tmp_path / "calls.jsonl"),
        "STUB_CFG_LOG": str(tmp_path / "ssh_config.seen"),
        "BIOBOX_POD_JUMP": JUMP,
        "BIOBOX_TOWER_IDENTITY": str(home / ".ssh" / "id_test"),
        "BIOBOX_RELAY_SSH_INCLUDE": str(home / "no-such-config"),
    })
    return {"tmp": tmp_path, "pod": pod, "dst": tower / "bio-box" / "root_backups", "tower": tower,
            "files": files, "env": env, "log": tmp_path / "relay.log"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_list(w: dict, entries: list[tuple[str, str]], form: str = "tsv") -> Path:
    p = w["tmp"] / f"list.{form}"
    if form == "tsv":
        p.write_text("".join(f"{rel}\t{h}\n" for rel, h in entries))
    else:
        p.write_text("".join(f"{h}  {rel}\n" for rel, h in entries))
    return p


def all_entries(w: dict) -> list[tuple[str, str]]:
    return [(rel, sha(data)) for rel, data in w["files"].items()]


def run(w: dict, *args: str, src: str | None = None, dst: str | None = None, lst: Path | None = None,
        env: dict | None = None) -> subprocess.CompletedProcess:
    lst = lst if lst is not None else write_list(w, all_entries(w))
    cmd = ["bash", str(SCRIPT), "--log", str(w["log"]), "--pause", "0", *args,
           src if src is not None else str(w["pod"]), dst if dst is not None else str(w["dst"]), str(lst)]
    return subprocess.run(cmd, env=env or w["env"], capture_output=True, text=True, timeout=120)


def calls(w: dict) -> list[dict]:
    p = Path(w["env"]["STUB_LOG"])
    return [json.loads(line) for line in p.read_text().splitlines()] if p.exists() else []


def tree(base: Path) -> dict[str, tuple[int, int, str]]:
    """Every file under `base`: (size, mtime_ns, sha256), by relative path."""
    if not base.exists():
        return {}
    return {p.relative_to(base).as_posix(): (p.stat().st_size, p.stat().st_mtime_ns, sha(p.read_bytes()))
            for p in sorted(base.rglob("*")) if p.is_file()}


def sums_lines(w: dict) -> list[str]:
    p = w["dst"] / "SHA256SUMS"
    return p.read_text().splitlines() if p.exists() else []


def sums_of(entries: list[tuple[str, str]]) -> list[str]:
    """`SHA256SUMS` lines for `entries`, in path order (bytewise, as `LC_ALL=C sort`)."""
    return [f"{h}  {rel}" for rel, h in sorted(entries)]


WRITES = re.compile(r"(\b(rm|mv|tee|mkdir|touch|cp|dd|chmod|chown|ln|truncate)\b|>)")


def assert_pod_read_only(w: dict, before: dict) -> None:
    assert tree(w["pod"]) == before, "a file on the pod side changed"
    pod_cmds = [c["cmd"] for c in calls(w) if c["tool"] == "ssh" and c["host"] == "pod"]
    assert pod_cmds, "the pod was never read"
    bad = [c for c in pod_cmds if WRITES.search(c)]
    assert not bad, f"a command sent to the pod can write: {bad}"


@pytest.mark.parametrize("streams", ["1", "2"])
def test_a_good_copy_lands_is_listed_and_the_pod_is_untouched(world, streams):
    w = world
    before = tree(w["pod"])
    r = run(w, "--streams", streams)
    assert r.returncode == 0, r.stdout + r.stderr
    for rel, data in w["files"].items():
        assert (w["dst"] / rel).read_bytes() == data
    assert not [p for p in w["tower"].rglob("*") if p.name.endswith(".partial")]
    assert sums_lines(w) == sums_of(all_entries(w))
    log = w["log"].read_text()
    assert log.count(": OK:") == 3, log
    # the hash is taken from the stream as it arrives on tower: no file is read back unless --reverify
    reads = [c["args"] for c in calls(w) if c["tool"] == "sha256sum"]
    assert len(reads) == 3 and all(a in ([], ["-"]) for a in reads), reads
    assert not [c for c in calls(w) if c["tool"] == "scp"]
    cfg = Path(w["env"]["STUB_CFG_LOG"]).read_text()
    assert f"ProxyJump {JUMP}" in cfg
    assert f"IdentityFile {w['env']['BIOBOX_TOWER_IDENTITY']}" in cfg
    assert "Include" not in cfg  # the include named does not exist
    assert_pod_read_only(w, before)


def test_reverify_reads_each_file_back_on_tower(world):
    w = world
    r = run(w, "--reverify")
    assert r.returncode == 0, r.stdout + r.stderr
    reads = [c["args"] for c in calls(w) if c["tool"] == "sha256sum"]
    reread = [a for a in reads if a not in ([], ["-"])]
    assert len(reread) == 3, reads
    assert all(a[0].endswith(".partial") for a in reread), reread
    assert len(sums_lines(w)) == 3


def test_a_sha256sums_list_is_read_and_an_existing_sums_file_is_extended(world):
    w = world
    w["dst"].mkdir(parents=True)
    other = f"{'a' * 64}  older/kept.tar"
    (w["dst"] / "SHA256SUMS").write_text(other + "\n")
    lst = write_list(w, all_entries(w), form="sums")
    r = run(w, lst=lst)
    assert r.returncode == 0, r.stdout + r.stderr
    assert sums_lines(w) == sums_of(all_entries(w) + [("older/kept.tar", "a" * 64)])
    assert other in sums_lines(w)


@pytest.mark.parametrize("fault", ["stream", "record"])
def test_a_corrupted_stream_is_refused_its_partial_removed_and_the_exit_non_zero(world, fault):
    w = world
    before = tree(w["pod"])
    entries = all_entries(w)
    env = dict(w["env"])
    if fault == "stream":
        env["STUB_CORRUPT"] = "1"
    else:  # the record says other bytes than the pod holds
        entries[0] = (entries[0][0], sha(b"something else"))
    r = run(w, "--tries", "2", lst=write_list(w, entries), env=env)
    assert r.returncode == 1, r.stdout + r.stderr
    bad = entries[0][0]
    assert not (w["dst"] / bad).exists()
    assert not [p for p in w["tower"].rglob("*") if p.name.endswith(".partial")]
    assert "differs" in w["log"].read_text()
    assert not [line for line in sums_lines(w) if line.endswith("  " + bad)]
    if fault == "record":  # the other files still land
        assert (w["dst"] / entries[1][0]).read_bytes() == w["files"][entries[1][0]]
    assert_pod_read_only(w, before)


@pytest.mark.parametrize("reverify", [False, True])
def test_a_file_in_place_is_skipped_and_never_overwritten(world, reverify):
    w = world
    rel = "root_rM2/l2.tar"
    (w["dst"] / "root_rM2").mkdir(parents=True)
    (w["dst"] / rel).write_bytes(b"not the pod's bytes")
    kept = tree(w["dst"])[rel]
    r = run(w, *(["--reverify"] if reverify else []))
    assert tree(w["dst"])[rel] == kept, "a file in place was changed"
    assert "already in place" in w["log"].read_text()
    others = [(x, h) for x, h in all_entries(w) if x != rel]
    for x, _ in others:
        assert (w["dst"] / x).read_bytes() == w["files"][x]
    # an unverified file in place is not listed as if it were verified
    assert sums_lines(w) == sums_of(others)
    if reverify:  # read back, it differs from the record: reported, still not touched
        assert r.returncode == 1, r.stdout + r.stderr
        assert "differs" in w["log"].read_text()
    else:
        assert r.returncode == 0, r.stdout + r.stderr


def test_a_file_missing_on_the_pod_fails_and_the_rest_still_land(world):
    w = world
    entries = all_entries(w) + [("root_rM2/absent.tar", "b" * 64)]
    r = run(w, lst=write_list(w, entries))
    assert r.returncode == 1, r.stdout + r.stderr
    assert "root_rM2/absent.tar" in w["log"].read_text()
    for rel, data in w["files"].items():
        assert (w["dst"] / rel).read_bytes() == data


def test_dry_run_copies_nothing(world):
    w = world
    before_pod, before_tower = tree(w["pod"]), tree(w["tower"])
    r = run(w, "--dry-run")
    assert r.returncode == 0, r.stdout + r.stderr
    assert tree(w["tower"]) == before_tower
    assert not w["dst"].exists()
    for rel, data in w["files"].items():
        assert f"copy {rel} {len(data)}" in r.stdout
    assert f"3 files to copy, {sum(len(d) for d in w['files'].values())} bytes" in r.stdout
    tower_cmds = [c["cmd"] for c in calls(w) if c["tool"] == "ssh" and c["host"] == "tower"]
    assert not [c for c in tower_cmds if re.search(r"\b(rm|mv|tee|mkdir|touch|cp)\b", c)], tower_cmds
    assert not w["log"].exists()
    assert_pod_read_only(w, before_pod)


BAD_LINES = [
    "a b.tar\t" + "c" * 64,
    "a;b.tar\t" + "c" * 64,
    "a$(id).tar\t" + "c" * 64,
    "a`id`.tar\t" + "c" * 64,
    "a'b.tar\t" + "c" * 64,
    "a|b.tar\t" + "c" * 64,
    "a&b.tar\t" + "c" * 64,
    "../escape.tar\t" + "c" * 64,
    "x/../../escape.tar\t" + "c" * 64,
    "/abs.tar\t" + "c" * 64,
    "-rf.tar\t" + "c" * 64,
    "ok.tar\t" + "C" * 64,
    "ok.tar\t" + "c" * 63,
    "ok.tar\t" + "c" * 63 + ";",
    "ok.tar\t" + "c" * 64 + "\textra",
    "c" * 64 + "  a;b.tar",
    "not a list line",
]


@pytest.mark.parametrize("line", BAD_LINES)
def test_a_path_or_checksum_with_a_metacharacter_is_refused_before_any_call(world, line):
    w = world
    p = w["tmp"] / "bad.tsv"
    p.write_text(f"added_r26_11_1.tar\t{sha(w['files']['added_r26_11_1.tar'])}\n{line}\n")
    r = run(w, lst=p)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == []
    assert not w["dst"].exists()


@pytest.mark.parametrize("where", ["src", "dst", "jump", "tower", "pod", "identity"])
def test_a_bad_directory_or_host_is_refused_before_any_call(world, where):
    w = world
    env = dict(w["env"])
    src, dst = str(w["pod"]), str(w["dst"])
    if where == "src":
        src = str(w["pod"]) + "/x;rm"
    elif where == "dst":
        dst = "relative/dir"
    elif where == "jump":
        env["BIOBOX_POD_JUMP"] = "u@h;rm -rf x"
    elif where == "tower":
        env["BIOBOX_TOWER"] = "root@tower x"
    elif where == "pod":
        env["BIOBOX_POD"] = "u@h$(id)"
    else:
        env["BIOBOX_TOWER_IDENTITY"] = "/tmp/id key"
    r = run(w, src=src, dst=dst, env=env)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == []

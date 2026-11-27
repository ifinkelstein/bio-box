"""Card C-73: `scripts/orchestrator/relay_to_pod.sh` copies files from the NAS `tower` to the pod through the laptop.

It is the mirror of `relay_to_tower.sh` (card C-49), and these tests use that card's fakes: stub `ssh`, `scp` and
`sha256sum` come first on PATH and record every call; nothing touches a network. The stub `ssh` runs the remote
command locally with `sh -c`. A local directory stands for `tower`'s source directory; for the pod, the stub rewrites
the two roots a destination may lie under (`/tmp/if2453/` and `/stor/scratch/Wilke/if2453/bio-box-staging/`) to two
local directories, so the script sends the pod's real paths and its whole logic (the list check, the plan, the stream,
the hash taken as the bytes arrive, the size check, `mv -n`, `STOP`, `SHA256SUMS`) runs here. `STUB_CORRUPT=1` flips
one byte of every stream on its way into the pod, as a fault in transit would; `STUB_STOP_AFTER_MV=<file>` makes that
file after the first file is put in place, as the PI dropping a `STOP` file during a run would.

Every test runs twice: with the normal PATH, and with the tools hidden (the stubs, then `/usr/bin:/bin` only), as
GitHub Actions has none of the pod's tools.
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

pytestmark = pytest.mark.covers("C-73")

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "orchestrator" / "relay_to_pod.sh"
JUMP = "if2453@wilkcomp03.ccbb.utexas.edu"
POD_TMP = "/tmp/if2453/"
POD_STOR = "/stor/scratch/Wilke/if2453/bio-box-staging/"
MIB64 = 64 * 1024 * 1024

SSH = """#!{python}
import json, os, re, shutil, subprocess, sys
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
if host == "pod":  # the pod's two roots, in one pass, longest first
    roots = json.loads(os.environ["STUB_POD_ROOTS"])
    pat = "|".join(re.escape(k) for k in sorted(roots, key=len, reverse=True))
    cmd = re.sub(pat, lambda m: roots[m.group(0)], cmd)
if no_stdin:
    r = subprocess.run(["sh", "-c", cmd], stdin=subprocess.DEVNULL)
elif host == "pod" and os.environ.get("STUB_CORRUPT") and "tee" in cmd:
    data = sys.stdin.buffer.read()
    if data:
        data = bytes([data[0] ^ 1]) + data[1:]
    r = subprocess.run(["sh", "-c", cmd], input=data)
else:
    r = subprocess.run(["sh", "-c", cmd])
stop = os.environ.get("STUB_STOP_AFTER_MV")
if host == "pod" and stop and "mv -n" in cmd and r.returncode == 0 and not os.path.exists(stop):
    open(stop, "w").close()
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


@pytest.fixture(params=["normal-path", "tools-hidden"])
def world(request, tmp_path: Path) -> dict:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for name, text in (("ssh", SSH), ("sha256sum", SHA256SUM), ("scp", SCP)):
        p = bindir / name
        p.write_text(text.format(python=sys.executable))
        p.chmod(0o755)
    src = tmp_path / "tower" / "sources_incoming" / "EGV1"
    (src / "tables").mkdir(parents=True)
    rnd = random.Random(73)
    files = {
        "egovirales.zip": bytes(rnd.getrandbits(8) for _ in range(300_000)),
        "tables/S1.tsv": bytes(rnd.getrandbits(8) for _ in range(170_000)),
        "README.txt": bytes(rnd.getrandbits(8) for _ in range(90_001)),
    }
    for rel, data in files.items():
        (src / rel).write_bytes(data)
    pod_tmp, pod_stor = tmp_path / "pod_tmp", tmp_path / "pod_stor"
    pod_tmp.mkdir()
    pod_stor.mkdir()
    home = tmp_path / "home"
    home.mkdir()
    env = dict(os.environ)
    for k in ("BIOBOX_POD", "BIOBOX_POD_JUMP", "BIOBOX_TOWER", "BIOBOX_TOWER_IDENTITY", "STUB_CORRUPT",
              "STUB_STOP_AFTER_MV"):
        env.pop(k, None)
    rest = env["PATH"] if request.param == "normal-path" else "/usr/bin:/bin"
    env.update({
        "PATH": f"{bindir}{os.pathsep}{rest}",
        "HOME": str(home),
        "STUB_LOG": str(tmp_path / "calls.jsonl"),
        "STUB_CFG_LOG": str(tmp_path / "ssh_config.seen"),
        "STUB_POD_ROOTS": json.dumps({POD_TMP: f"{pod_tmp}/", POD_STOR: f"{pod_stor}/"}),
        "BIOBOX_POD_JUMP": JUMP,
        "BIOBOX_TOWER_IDENTITY": str(home / ".ssh" / "id_test"),
        "BIOBOX_RELAY_SSH_INCLUDE": str(home / "no-such-config"),
    })
    return {"tmp": tmp_path, "src": src, "pod_tmp": pod_tmp, "pod_stor": pod_stor,
            "dst": POD_TMP + "c73/EGV1", "dst_local": pod_tmp / "c73" / "EGV1",
            "files": files, "env": env, "log": tmp_path / "relay.log"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_list(w: dict, entries: list[tuple[str, str]], form: str = "tsv") -> Path:
    p = w["tmp"] / f"list.{form}"
    if form == "tsv":
        p.write_text("".join(f"{rel}\t{h}\n" for rel, h in entries))
    elif form == "sums-binary":
        p.write_text("".join(f"{h} *{rel}\n" for rel, h in entries))
    else:
        p.write_text("".join(f"{h}  {rel}\n" for rel, h in entries))
    return p


def all_entries(w: dict) -> list[tuple[str, str]]:
    return [(rel, sha(data)) for rel, data in w["files"].items()]


def run(w: dict, *args: str, src: str | None = None, dst: str | None = None, lst: Path | None = None,
        env: dict | None = None) -> subprocess.CompletedProcess:
    lst = lst if lst is not None else write_list(w, all_entries(w))
    cmd = ["bash", str(SCRIPT), "--log", str(w["log"]), "--pause", "0", *args,
           src if src is not None else str(w["src"]), dst if dst is not None else w["dst"], str(lst)]
    return subprocess.run(cmd, env=env or w["env"], capture_output=True, text=True, timeout=300)


def calls(w: dict) -> list[dict]:
    p = Path(w["env"]["STUB_LOG"])
    return [json.loads(line) for line in p.read_text().splitlines()] if p.exists() else []


def host_cmds(w: dict, host: str) -> list[str]:
    return [c["cmd"] for c in calls(w) if c["tool"] == "ssh" and c["host"] == host]


def tree(base: Path) -> dict[str, tuple[int, int, str]]:
    """Every file under `base`: (size, mtime_ns, sha256), by relative path."""
    if not base.exists():
        return {}
    return {p.relative_to(base).as_posix(): (p.stat().st_size, p.stat().st_mtime_ns, sha(p.read_bytes()))
            for p in sorted(base.rglob("*")) if p.is_file()}


def sums_lines(w: dict) -> list[str]:
    p = w["dst_local"] / "SHA256SUMS"
    return p.read_text().splitlines() if p.exists() else []


def sums_of(entries: list[tuple[str, str]]) -> list[str]:
    """`SHA256SUMS` lines for `entries`, in path order (bytewise, as `LC_ALL=C sort`)."""
    return [f"{h}  {rel}" for rel, h in sorted(entries)]


def partials(w: dict) -> list[Path]:
    return [p for base in (w["pod_tmp"], w["pod_stor"]) for p in base.rglob("*") if p.name.endswith(".partial")]


WRITES = re.compile(r"(\b(rm|mv|tee|mkdir|touch|cp|dd|chmod|chown|ln|truncate|sha256sum)\b|>)")


def assert_tower_read_only(w: dict, before: dict) -> None:
    assert tree(w["src"]) == before, "a file on tower changed"
    cmds = host_cmds(w, "tower")
    assert cmds, "tower was never read"
    bad = [c for c in cmds if WRITES.search(c)]
    assert not bad, f"a command sent to tower can write: {bad}"
    words = {m for c in cmds for m in re.findall(r"\b(stat|cat)\b", c)}
    assert words <= {"stat", "cat"} and "stat" in words, cmds


@pytest.mark.parametrize("streams", ["1", "2"])
def test_a_good_copy_lands_is_listed_and_tower_is_untouched(world, streams):
    w = world
    before = tree(w["src"])
    r = run(w, "--streams", streams)
    assert r.returncode == 0, r.stdout + r.stderr
    for rel, data in w["files"].items():
        assert (w["dst_local"] / rel).read_bytes() == data
    assert not partials(w)
    assert sums_lines(w) == sums_of(all_entries(w))
    log = w["log"].read_text()
    assert log.count(": OK:") == 3, log
    # the pod hashes the stream as it arrives: no file is read back unless --reverify
    reads = [c["args"] for c in calls(w) if c["tool"] == "sha256sum"]
    assert len(reads) == 3 and all(a in ([], ["-"]) for a in reads), reads
    # each stream is written as <dir>/.<name>.partial and put in place by mv -n
    tees = [c for c in host_cmds(w, "pod") if "tee" in c]
    assert sorted(tees) == sorted(
        f"tee -- '{w['dst']}/{rel.rsplit('/', 1)[0] + '/' if '/' in rel else ''}.{rel.rsplit('/', 1)[-1]}.partial'"
        " | sha256sum" for rel in w["files"]), tees
    assert len([c for c in host_cmds(w, "pod") if "mv -n" in c]) == 3
    assert not [c for c in calls(w) if c["tool"] == "scp"]
    cfg = Path(w["env"]["STUB_CFG_LOG"]).read_text()
    assert f"ProxyJump {JUMP}" in cfg
    assert f"IdentityFile {w['env']['BIOBOX_TOWER_IDENTITY']}" in cfg
    assert "Include" not in cfg  # the include named does not exist
    assert_tower_read_only(w, before)


def test_reverify_reads_each_file_back_on_the_pod(world):
    w = world
    before = tree(w["src"])
    r = run(w, "--reverify")
    assert r.returncode == 0, r.stdout + r.stderr
    reads = [c["args"] for c in calls(w) if c["tool"] == "sha256sum"]
    reread = [a for a in reads if a not in ([], ["-"])]
    assert len(reread) == 3, reads
    assert all(a[0].endswith(".partial") and a[0].startswith(str(w["pod_tmp"])) for a in reread), reread
    assert len(sums_lines(w)) == 3
    assert_tower_read_only(w, before)


@pytest.mark.parametrize("form", ["sums", "sums-binary"])
def test_a_sha256sums_list_is_read_and_an_existing_sums_file_is_extended(world, form):
    w = world
    w["dst_local"].mkdir(parents=True)
    other = f"{'a' * 64}  older/kept.tar"
    (w["dst_local"] / "SHA256SUMS").write_text(other + "\n")
    r = run(w, lst=write_list(w, all_entries(w), form=form))
    assert r.returncode == 0, r.stdout + r.stderr
    assert sums_lines(w) == sums_of(all_entries(w) + [("older/kept.tar", "a" * 64)])


def test_a_conflict_in_sha256sums_leaves_it_as_it_was(world):
    w = world
    w["dst_local"].mkdir(parents=True)
    old = f"{'d' * 64}  README.txt\n"
    (w["dst_local"] / "SHA256SUMS").write_text(old)
    r = run(w)
    assert r.returncode == 1, r.stdout + r.stderr
    assert (w["dst_local"] / "SHA256SUMS").read_text() == old
    assert "CONFLICT" in w["log"].read_text()


@pytest.mark.parametrize("fault", ["stream", "record"])
def test_a_corrupted_stream_is_refused_its_partial_removed_and_the_exit_non_zero(world, fault):
    w = world
    before = tree(w["src"])
    entries = all_entries(w)
    env = dict(w["env"])
    if fault == "stream":
        env["STUB_CORRUPT"] = "1"
    else:  # the record says other bytes than tower holds
        entries[0] = (entries[0][0], sha(b"something else"))
    r = run(w, "--tries", "2", lst=write_list(w, entries), env=env)
    assert r.returncode == 1, r.stdout + r.stderr
    bad = entries[0][0]
    assert not (w["dst_local"] / bad).exists()
    assert not partials(w)
    assert "differs" in w["log"].read_text()
    assert not [line for line in sums_lines(w) if line.endswith("  " + bad)]
    if fault == "record":  # the other files still land
        assert (w["dst_local"] / entries[1][0]).read_bytes() == w["files"][entries[1][0]]
    assert_tower_read_only(w, before)


def test_a_partial_left_by_a_killed_run_is_replaced(world):
    w = world
    (w["dst_local"]).mkdir(parents=True)
    (w["dst_local"] / ".README.txt.partial").write_bytes(b"half a file")
    r = run(w)
    assert r.returncode == 0, r.stdout + r.stderr
    assert (w["dst_local"] / "README.txt").read_bytes() == w["files"]["README.txt"]
    assert not partials(w)


@pytest.mark.parametrize("reverify", [False, True])
def test_a_file_in_place_is_skipped_and_never_overwritten(world, reverify):
    w = world
    rel = "tables/S1.tsv"
    (w["dst_local"] / "tables").mkdir(parents=True)
    (w["dst_local"] / rel).write_bytes(b"not tower's bytes")
    kept = tree(w["dst_local"])[rel]
    r = run(w, *(["--reverify"] if reverify else []))
    assert tree(w["dst_local"])[rel] == kept, "a file in place was changed"
    assert "already in place" in w["log"].read_text()
    assert not [c for c in host_cmds(w, "tower") if rel in c and "cat" in c], "a file in place was streamed"
    others = [(x, h) for x, h in all_entries(w) if x != rel]
    for x, _ in others:
        assert (w["dst_local"] / x).read_bytes() == w["files"][x]
    # an unverified file in place is not listed as if it were verified
    assert sums_lines(w) == sums_of(others)
    if reverify:  # read back, it differs from the record: reported, still not touched
        assert r.returncode == 1, r.stdout + r.stderr
        assert "differs" in w["log"].read_text()
    else:
        assert r.returncode == 0, r.stdout + r.stderr


def test_a_file_missing_on_tower_fails_and_the_rest_still_land(world):
    w = world
    entries = all_entries(w) + [("tables/absent.tsv", "b" * 64)]
    r = run(w, lst=write_list(w, entries))
    assert r.returncode == 1, r.stdout + r.stderr
    assert "tables/absent.tsv" in w["log"].read_text()
    for rel, data in w["files"].items():
        assert (w["dst_local"] / rel).read_bytes() == data


def test_dry_run_copies_nothing_and_writes_nothing(world):
    w = world
    before_src, before_tmp, before_stor = tree(w["src"]), tree(w["pod_tmp"]), tree(w["pod_stor"])
    r = run(w, "--dry-run")
    assert r.returncode == 0, r.stdout + r.stderr
    assert tree(w["pod_tmp"]) == before_tmp and tree(w["pod_stor"]) == before_stor
    assert not w["dst_local"].exists()
    for rel, data in w["files"].items():
        assert f"copy {rel} {len(data)}" in r.stdout
    assert f"3 files to copy, {sum(len(d) for d in w['files'].values())} bytes" in r.stdout
    pod = host_cmds(w, "pod")
    assert pod and not [c for c in pod if re.search(r"\b(rm|mv|tee|mkdir|touch|cp|sha256sum)\b", c)], pod
    assert not w["log"].exists()
    assert_tower_read_only(w, before_src)


def test_a_stop_file_in_the_destination_ends_the_run_before_the_next_file(world):
    w = world
    stop = w["dst_local"] / "STOP"
    env = dict(w["env"], STUB_STOP_AFTER_MV=str(stop))
    w["dst_local"].mkdir(parents=True)
    r = run(w, env=env)
    assert r.returncode == 1, r.stdout + r.stderr
    # largest first: the zip lands, then the STOP file is seen and nothing more is streamed
    assert (w["dst_local"] / "egovirales.zip").read_bytes() == w["files"]["egovirales.zip"]
    assert not (w["dst_local"] / "README.txt").exists() and not (w["dst_local"] / "tables").exists()
    assert len([c for c in host_cmds(w, "tower") if "cat --" in c]) == 1
    assert not partials(w)
    assert "STOP" in w["log"].read_text()
    assert sums_lines(w) == sums_of([("egovirales.zip", sha(w["files"]["egovirales.zip"]))])
    # the STOP file removed, a second run copies the rest and skips what is in place
    stop.unlink()
    r = run(w)
    assert r.returncode == 0, r.stdout + r.stderr
    for rel, data in w["files"].items():
        assert (w["dst_local"] / rel).read_bytes() == data
    assert sums_lines(w) == sums_of(all_entries(w))


def test_a_stop_file_present_at_the_start_copies_nothing(world):
    w = world
    w["dst_local"].mkdir(parents=True)
    (w["dst_local"] / "STOP").write_text("")
    r = run(w)
    assert r.returncode == 1, r.stdout + r.stderr
    assert sorted(tree(w["dst_local"])) == ["STOP"]
    assert not [c for c in host_cmds(w, "tower") if "cat --" in c]
    assert not [c for c in host_cmds(w, "pod") if re.search(r"\b(rm|mv|tee|mkdir)\b", c)]
    assert "STOP" in w["log"].read_text()


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
    "SHA256SUMS\t" + "c" * 64,
    "x/.y.tar.partial\t" + "c" * 64,
    "STOP\t" + "c" * 64,
    "README.txt\t" + "c" * 64,  # the same path twice
]


@pytest.mark.parametrize("line", BAD_LINES)
def test_a_bad_list_line_is_refused_before_any_call(world, line):
    w = world
    p = w["tmp"] / "bad.tsv"
    p.write_text(f"README.txt\t{sha(w['files']['README.txt'])}\n{line}\n")
    r = run(w, lst=p)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == []
    assert not w["dst_local"].exists()


@pytest.mark.parametrize("where", ["src", "dst", "jump", "tower", "pod", "identity", "include"])
def test_a_bad_directory_or_host_is_refused_before_any_call(world, where):
    w = world
    env = dict(w["env"])
    src, dst = str(w["src"]), w["dst"]
    if where == "src":
        src = str(w["src"]) + "/x;rm"
    elif where == "dst":
        dst = "relative/dir"
    elif where == "jump":
        env["BIOBOX_POD_JUMP"] = "u@h;rm -rf x"
    elif where == "tower":
        env["BIOBOX_TOWER"] = "root@tower x"
    elif where == "pod":
        env["BIOBOX_POD"] = "u@h$(id)"
    elif where == "identity":
        env["BIOBOX_TOWER_IDENTITY"] = "/tmp/id key"
    else:
        env["BIOBOX_RELAY_SSH_INCLUDE"] = "relative/config"
    r = run(w, src=src, dst=dst, env=env)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == []


@pytest.mark.parametrize("dst", [
    "/tmp/other/EGV1",
    "/tmp/if2453",
    "/tmp/if2453x/EGV1",
    "/tmp/if2453/./EGV1",
    "/tmp/if2453/../etc",
    "/stor/scratch/Wilke/if2453/bio-box-staging",
    "/stor/scratch/Wilke/if2453/all-the-data/EGV1",
    "/stor/work/Wilke/if2453/EGV1",
    "/stor/home/if2453/EGV1",
    "/home/if2453/EGV1",
])
def test_a_destination_outside_the_two_roots_is_refused_before_any_call(world, dst):
    w = world
    r = run(w, dst=dst)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == []
    assert tree(w["pod_tmp"]) == {} and tree(w["pod_stor"]) == {}


@pytest.mark.parametrize("dry", [False, True])
def test_a_file_under_64_mib_bound_for_stor_is_refused_before_any_call_to_the_pod(world, dry):
    w = world
    before = tree(w["src"])
    big = w["src"] / "big.zip"
    with open(big, "wb") as f:  # exactly 64 MiB, sparse on the tower side
        f.truncate(MIB64)
    entries = [("big.zip", sha(b"\0" * MIB64))] + all_entries(w)
    r = run(w, *(["--dry-run"] if dry else []), dst=POD_STOR + "EGV1", lst=write_list(w, entries))
    assert r.returncode == 2, r.stdout + r.stderr
    for rel in w["files"]:
        assert rel in r.stderr
    assert "big.zip" not in r.stderr
    assert host_cmds(w, "pod") == [], "the pod was called"
    assert tree(w["pod_stor"]) == {} and tree(w["pod_tmp"]) == {}
    assert not [c for c in host_cmds(w, "tower") if "cat --" in c]
    before["big.zip"] = tree(w["src"])["big.zip"]
    assert_tower_read_only(w, before)


def test_a_file_of_64_mib_lands_on_stor_and_no_small_file_is_written_there(world):
    w = world
    big = w["src"] / "big.zip"
    with open(big, "wb") as f:
        f.truncate(MIB64)
    entries = [("big.zip", sha(b"\0" * MIB64))]
    r = run(w, dst=POD_STOR + "EGV1", lst=write_list(w, entries))
    assert r.returncode == 0, r.stdout + r.stderr
    landed = tree(w["pod_stor"])
    assert sorted(landed) == ["EGV1/big.zip"], landed
    assert landed["EGV1/big.zip"][0] == MIB64 and landed["EGV1/big.zip"][2] == entries[0][1]
    assert all(size >= MIB64 for size, _, _ in landed.values())  # no SHA256SUMS on /stor
    log = w["log"].read_text()
    assert entries[0][1] in log and "SHA256SUMS: not written" in log


# --- card R-18: --space-as-underscore ---------------------------------------------------------------------------------

OPT = "--space-as-underscore"


def add_file(w: dict, rel: str, data: bytes) -> tuple[str, str]:
    (w["src"] / rel).write_bytes(data)
    return rel, sha(data)


def test_r18_a_space_is_refused_without_the_option_and_the_message_names_it(world):
    w = world
    rel, h = add_file(w, "a b.tar", b"spaced bytes")
    r = run(w, lst=write_list(w, [(rel, h)]))
    assert r.returncode == 2, r.stdout + r.stderr
    assert OPT in r.stderr
    assert calls(w) == []


@pytest.mark.parametrize("streams", ["1", "2"])
def test_r18_with_the_option_a_spaced_name_lands_under_its_underscore_name(world, streams):
    w = world
    data = bytes(random.Random(18).getrandbits(8) for _ in range(50_000))
    rel, h = add_file(w, "a b.tar", data)
    nested = "dir one/c d e.tsv"
    (w["src"] / "dir one").mkdir()
    entries = [add_file(w, nested, b"nested spaced")] + [(rel, h)]
    before = tree(w["src"])
    r = run(w, OPT, "--streams", streams, lst=write_list(w, entries))
    assert r.returncode == 0, r.stdout + r.stderr
    assert (w["dst_local"] / "a_b.tar").read_bytes() == data
    assert (w["dst_local"] / "dir_one" / "c_d_e.tsv").read_bytes() == b"nested spaced"
    assert not (w["dst_local"] / "a b.tar").exists()
    assert not partials(w)
    assert sums_lines(w) == sums_of([("a_b.tar", h), ("dir_one/c_d_e.tsv", sha(b"nested spaced"))])
    log = w["log"].read_text()
    assert "a b.tar" in log and "a_b.tar" in log and f"sha256 {h}" in log
    assert (w["src"] / "a b.tar").read_bytes() == data and tree(w["src"]) == before


def test_r18_dry_run_lists_the_pod_name_and_writes_nothing(world):
    w = world
    rel, h = add_file(w, "Dereplicated genomes.zip", b"zip bytes")
    r = run(w, OPT, "--dry-run", lst=write_list(w, [(rel, h)]))
    assert r.returncode == 0, r.stdout + r.stderr
    assert "Dereplicated genomes.zip" in r.stdout and "Dereplicated_genomes.zip" in r.stdout
    assert not w["dst_local"].exists() and not w["log"].exists()


@pytest.mark.parametrize("name", ["a  b.tar", " a.tar", "a.tar ", "a b.tar ", "a ;b.tar", "d/ e.tar", "d /e.tar", "a\tb.tar"])
def test_r18_with_the_option_other_spacings_are_still_refused_before_any_call(world, name):
    w = world
    p = w["tmp"] / "bad.tsv"
    p.write_text(f"README.txt\t{sha(w['files']['README.txt'])}\n{name}\t{'c' * 64}\n")
    r = run(w, OPT, lst=p)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == []
    assert not w["dst_local"].exists()


@pytest.mark.parametrize("pair", [("a b.tar", "a_b.tar"), ("a_b.tar", "a b.tar")])
def test_r18_two_lines_with_one_pod_name_are_refused_before_any_call(world, pair):
    w = world
    p = w["tmp"] / "dup.tsv"
    p.write_text(f"{pair[0]}\t{'c' * 64}\n{pair[1]}\t{'d' * 64}\n")
    r = run(w, OPT, lst=p)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == []


def test_r18_a_pod_name_that_is_the_relays_own_is_refused(world):
    w = world
    p = w["tmp"] / "own.tsv"
    p.write_text(f"x y.partial\t{'c' * 64}\n")
    r = run(w, OPT, lst=p)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == []


@pytest.mark.parametrize("name", ["a b.tar", "a =sh", "a ,b"])
def test_r18_the_tower_path_is_one_single_quoted_word(world, name):
    w = world
    data = name.encode() + b" content"
    rel, h = add_file(w, name, data)
    r = run(w, OPT, lst=write_list(w, [(rel, h)]))
    assert r.returncode == 0, r.stdout + r.stderr
    assert (w["dst_local"] / name.replace(" ", "_")).read_bytes() == data
    tower_cmds = host_cmds(w, "tower")
    assert f"cat -- '{w['src']}/{name}'" in tower_cmds, tower_cmds
    for c in host_cmds(w, "tower") + host_cmds(w, "pod"):
        assert name not in re.sub(r"'[^']*'", "''", c), c
    assert not [c for c in host_cmds(w, "pod") if name in c], "the tower name reached the pod"

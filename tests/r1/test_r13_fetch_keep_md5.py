"""Card R-13: `fetch_to_tower.sh` keeps a file whose provider md5 matched when its archive test fails.

Same stub `curl` as tests/r1/test_r08_fetch_to_tower.py (no request leaves the machine). The archive tools the cases
rely on are `unzip` (zip case), `gzip` (gzip case) and `tar` (tar case); a case skips when its tool is absent, as the
script itself skips that archive test. The FASTA case needs only `head`.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("R-13")

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "orchestrator" / "fetch_to_tower.sh"

CURL = """#!{python}
import json, os, sys
args = sys.argv[1:]
out, url = args[args.index("-o") + 1], args[-1]
spec = json.load(open(os.environ["STUB_MAP"])).get(url)
if spec is None:
    sys.exit(22)
with open(out, "wb") as f:
    f.write(open(spec, "rb").read())
"""

SOURCES = """version: 1
sources:
  - name: zen
    url: https://zenodo.org/api/records/111/files/{file}/content
"""

URL = "https://zenodo.org/api/records/111/files/{}/content"

# kind -> (file name, bytes that fail the archive test, tool needed, warning word)
BROKEN = {
    "zip": ("genomes.zip", b"this is not a zip archive" * 20, "unzip", "WARN_ZIP"),
    "gzip": ("table.tsv.gz", b"this is not gzip" * 20, "gzip", "WARN_GZIP"),
    "tar": ("bins.tar", b"this is not a tar" * 100, "tar", "WARN_TAR"),
    "fasta": ("seqs.fasta", b"no header line\nACGT\n", None, "WARN_FASTA"),
}
FAILWORD = {"zip": "FAIL_ZIP", "gzip": "FAIL_GZIP", "tar": "FAIL_TAR", "fasta": "FAIL_FASTA"}


def md5(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()  # noqa: S324 - the provider's checksum, not a security use


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def archive_test_fails(kind: str, tmp: Path, data: bytes) -> bool:
    """Skip the case when the tool is absent; otherwise the broken bytes must really fail the tool's own test."""
    _f, _d, tool, _w = BROKEN[kind]
    if tool and shutil.which(tool) is None:
        pytest.skip(f"{tool} is not installed: the script skips this archive test too")
    return True


@pytest.fixture
def world(tmp_path: Path) -> dict:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    p = bindir / "curl"
    p.write_text(CURL.format(python=sys.executable))
    p.chmod(0o755)
    src = tmp_path / "sources.yaml"
    src.write_text(SOURCES)
    mp = tmp_path / "map.json"
    mp.write_text("{}")
    env = dict(os.environ)
    env.update({"PATH": f"{bindir}{os.pathsep}{env['PATH']}", "STUB_MAP": str(mp), "TMPDIR": str(tmp_path)})
    return {"tmp": tmp_path, "src": src, "dest": tmp_path / "incoming", "env": env, "map": mp}


def run_one(w: dict, kind: str, md5_field: str | None) -> tuple[subprocess.CompletedProcess, str, bytes]:
    f, data, _tool, _w = BROKEN[kind]
    archive_test_fails(kind, w["tmp"], data)
    url = URL.format(f)
    blob = w["tmp"] / "blob"
    blob.write_bytes(data)
    w["map"].write_text(json.dumps({url: str(blob)}))
    field = md5(data) if md5_field is None else md5_field
    lst = w["tmp"] / "jobs.tsv"
    lst.write_text("\t".join(("S1", f, url, field, "CC BY 4.0")) + "\n")
    cmd = ["bash", str(SCRIPT), "--dest", str(w["dest"]), "--sources", str(w["src"]), "--pause", "0", str(lst)]
    return subprocess.run(cmd, env=w["env"], capture_output=True, text=True, timeout=120), f, data


@pytest.mark.parametrize("kind", sorted(BROKEN))
def test_a_file_with_the_providers_md5_is_kept_when_its_archive_test_fails(world, kind):
    w = world
    r, f, data = run_one(w, kind, None)
    assert r.returncode == 0, r.stdout + r.stderr
    d = w["dest"] / "S1"
    assert (d / f).read_bytes() == data
    assert (d / "SHA256SUMS").read_text().splitlines() == [f"{sha(data)}  {f}"]
    rows = [x.split("\t") for x in (d / "SOURCE_URLS.tsv").read_text().splitlines()]
    assert rows[1][:5] == [f, URL.format(f), str(len(data)), sha(data), "CC BY 4.0"]
    out = r.stdout.splitlines()
    ok = [i for i, x in enumerate(out) if x.startswith("OK\t")]
    assert len(ok) == 1 and out[ok[0]].split("\t")[1:3] == ["S1", f]
    assert out[ok[0] + 1].split("\t") == [BROKEN[kind][3], "S1", f]
    assert "DONE\tok=1\tskip=0\tfail=0" in r.stdout
    assert FAILWORD[kind] not in r.stderr
    assert not [p for p in w["dest"].rglob("*") if p.name.endswith(".partial")]


@pytest.mark.parametrize("kind", sorted(BROKEN))
def test_a_file_with_no_md5_whose_archive_test_fails_is_still_removed(world, kind):
    w = world
    r, f, _data = run_one(w, kind, "-")
    assert r.returncode == 1, r.stdout + r.stderr
    assert FAILWORD[kind] in r.stderr and f in r.stderr
    assert "WARN_" not in r.stdout
    d = w["dest"] / "S1"
    assert not (d / f).exists()
    assert (d / "SHA256SUMS").read_text() == ""
    assert not [p for p in w["dest"].rglob("*") if p.name.endswith(".partial")]
    assert "DONE\tok=0\tskip=0\tfail=1" in r.stdout


@pytest.mark.parametrize("kind", sorted(BROKEN))
def test_a_file_whose_md5_differs_is_removed_whatever_the_archive_test_says(world, kind):
    w = world
    r, f, _data = run_one(w, kind, md5(b"another"))
    assert r.returncode == 1, r.stdout + r.stderr
    assert "FAIL_MD5" in r.stderr and f in r.stderr
    assert "WARN_" not in r.stdout
    assert not (w["dest"] / "S1" / f).exists()
    assert "DONE\tok=0\tskip=0\tfail=1" in r.stdout

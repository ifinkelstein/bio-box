"""D-05 (card D-04): the micro-S0 fixture store has the same manifest hash on another host.

The peer is the GitHub Actions `cross-host` job (`.github/workflows/cross_host.yml`), which builds the
fixture with `scripts/stress/build_fixture.py` and prints the report; the laptop no longer computes.
Get the peer's hash from the run for the commit under test and give the pod a file that holds it:

    gh run view <run id> --log | grep -o '"manifest_hash": *"[0-9a-f]*"' | sort -u > peer.json

then run the test with BIOBOX_PEER_MANIFEST=<path of peer.json>. The file may also be the whole report
JSON of another host, or a bare 32-hex hash. This test builds the fixture itself with the same script,
so the two hashes are for the same code and fixture.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
ENV = "BIOBOX_PEER_MANIFEST"
_HASH = re.compile(r'"?manifest_hash"?\s*[:=]\s*"?([0-9a-f]{32})\b')
_BARE = re.compile(r"[0-9a-f]{32}")


def peer_hash(path: Path) -> str:
    """The manifest hash in a peer file: a report JSON, a `"manifest_hash": "..."` line, or a bare hash."""
    text = path.read_text()
    found = set(_HASH.findall(text))
    if not found and _BARE.fullmatch(text.strip()):
        found = {text.strip()}
    if len(found) != 1:
        raise ValueError(f"{path}: expected exactly one manifest hash, found {sorted(found)}")
    return found.pop()


@pytest.fixture(scope="module")
def own_hash(tmp_path_factory) -> str:
    root = tmp_path_factory.mktemp("xhost") / "fixture_root"
    out = subprocess.run(
        [sys.executable, str(ROOT / "scripts/stress/build_fixture.py"), str(ROOT), str(root), "rX.fixture"],
        capture_output=True, text=True, check=True,
    ).stdout
    report = json.loads(out)
    assert report["genomes"] == 3
    return report["manifest_hash"]


@pytest.mark.covers("D-05")
def test_d05_cross_host_manifest_hash_equal(own_hash, record_property):
    path = os.environ.get(ENV)
    if not path:
        pytest.skip(f"{ENV} not set: no manifest hash from another host to compare with")
    theirs = peer_hash(Path(path))
    record_property("own_manifest_hash", own_hash)
    record_property("peer_manifest_hash", theirs)
    assert own_hash == theirs, f"manifest hash differs from the peer's ({path}): own {own_hash}, peer {theirs}"


@pytest.mark.covers("D-05")
def test_d05_wrong_peer_hash_is_detected(own_hash, tmp_path):
    wrong = "0" * 32 if own_hash != "0" * 32 else "1" * 32
    p = tmp_path / "peer.json"
    p.write_text(json.dumps({"manifest_hash": wrong}))
    assert peer_hash(p) != own_hash
    p.write_text(f'"manifest_hash": "{own_hash}"\n')
    assert peer_hash(p) == own_hash
    p.write_text(own_hash + "\n")
    assert peer_hash(p) == own_hash
    p.write_text("no hash here\n")
    with pytest.raises(ValueError):
        peer_hash(p)

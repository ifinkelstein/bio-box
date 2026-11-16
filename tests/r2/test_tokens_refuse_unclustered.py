"""Card B-51: the tokens step refuses a release that was never clustered.

A release that is defined over a stage with tapes but has no `cluster_member` file would get tokens without
NR90 and NR50 levels, silently. `release tokens` exits 2 with a message that names the release and the missing
clustering, and writes nothing. A release whose tapes hold no protein slot builds as before (the check is on
the tapes' protein slots, not on the release's files alone).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

from biobox import catalog, cli
from biobox.l3 import tokens

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_release_def as rd  # noqa: E402

pytestmark = pytest.mark.r2

REL = "rB51.unclustered"


@pytest.fixture(scope="module")
def staged(tmp_path_factory):
    d = tmp_path_factory.mktemp("b51")
    root = d / "root"
    rd.add_stage(root, "S0", d / "lib.hmm")
    return root


def _files(root: Path) -> set[str]:
    return {str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()}


@pytest.mark.covers("D-01")
def test_b51_cli_refuses_a_defined_release_without_clustering_and_writes_nothing(staged, tmp_path, capsys):
    c = tmp_path / "c"
    shutil.copytree(staged, c)
    base = ["--root", str(c), "--release", REL]
    assert cli.main(["release", "define", *base, "--stages", "S0"]) == 0
    capsys.readouterr()
    before = _files(c)
    assert cli.main(["release", "tokens", *base, "--tmp", str(tmp_path)]) == 2
    err = capsys.readouterr().err
    assert REL in err and "cluster" in err and "release cluster" in err
    assert _files(c) == before  # nothing written: no token_dict, no tape_tokens, no scratch in the root


@pytest.mark.covers("D-01")
def test_b51_library_call_raises_before_the_definition_is_written(staged, tmp_path):
    c = tmp_path / "c"
    shutil.copytree(staged, c)
    before = _files(c)
    with pytest.raises(FileNotFoundError, match=REL) as e:
        tokens.run(c, REL, ["S0"], tmp=tmp_path)
    assert "cluster" in str(e.value)
    assert _files(c) == before
    assert not catalog.release_def_path(c, REL).exists()


@pytest.mark.covers("D-01")
def test_b51_a_release_without_protein_slots_needs_no_clustering(staged, tmp_path):
    c = tmp_path / "c"
    shutil.copytree(staged, c)
    # an unrelated tape version whose features are none: tapes without slots
    from biobox import io
    from biobox.l3 import tapes
    from biobox.registry.schemas import table

    io.write_table([], table("feature"), c / "l1/feature/producer=test-none/stage=S0/part-0.parquet")
    circ = ()
    tv = tapes.run(c, "S0", circ, producers=["test-none"])["tape_version"]
    res = tokens.run(c, REL + "-empty", ["S0"], producers=["test-none"])
    assert res["tokens"] == {"family": 0, "nr50": 0, "nr90": 0}
    assert (c / f"l3/tape_tokens/tape_version={tv}/release={REL}-empty").is_dir()

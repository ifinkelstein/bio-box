"""Card E-07: the frozen ncRNA model set `prok-1` is the configuration DEC-28 decided.

The oracle is the published design note `planning/design/0-08_ncrna_models.md` read here, not its rule script
(which generated the file): the set is Appendix A's 258 models with a hit (rule a) less the four eukaryotic rRNA
models of rule (a-), plus Appendix B's models added by rule (b), 588 in all; the whole-genome models are the 23
in-CDS models of the note's section 4, transcribed below. Every accession and name must be a model of Rfam 15.1 as
the note's model table lists it (`0-08_rfam15.1_models.tsv`).
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pytest

from biobox.l1 import ncrna

pytestmark = pytest.mark.covers("E-07")  # ring R1 (by directory)

REPO = Path(__file__).resolve().parents[2]
NOTE = REPO / "planning/design/0-08_ncrna_models.md"
MODELS = REPO / "planning/design/0-08_rfam15.1_models.tsv"
# section 2, rule (a-): SSU_rRNA_eukarya, SSU_rRNA_microsporidia, LSU_rRNA_eukarya, 5_8S_rRNA
EUK_RRNA = {"RF01960", "RF02542", "RF02543", "RF00002"}
# section 4, the in-CDS class: 23 models scanned on whole genomes
IN_CDS = {"RF00382", "RF01779", "RF01781", "RF01782", "RF01843", "RF01988", "RF01989", "RF02194", "RF02376",
          "RF02465", "RF02514", "RF02592", "RF02654", "RF02662", "RF02674", "RF02736", "RF02737", "RF02738",
          "RF02741", "RF02742", "RF02746", "RF02757", "RF02850"}


def note_sets() -> tuple[dict[str, str], dict[str, str]]:
    """(Appendix A's models, Appendix B's models), accession -> name, read from the note's text."""
    text = NOTE.read_text()
    a_text = text.split("## Appendix A", 1)[1].split("## Appendix B", 1)[0]
    b_text = text.split("## Appendix B", 1)[1]
    a = {m[1]: m[2] for m in re.finditer(r"^\| (RF\d{5}) \| ([^|]+?) \|", a_text, re.M)}
    b = {}
    for line in b_text.splitlines():
        if line.startswith("- (b"):
            b.update({m[1]: m[2] for m in re.finditer(r"(RF\d{5}) ([^,\s]+)", line.split(":", 1)[1])})
    return a, b


def test_prok1_is_the_decided_configuration():
    a, b = note_sets()
    assert len(a) == 258  # the note's rule (a)
    assert EUK_RRNA <= set(a) and not set(a) & set(b)
    want = {**{k: v for k, v in a.items() if k not in EUK_RRNA}, **b}
    assert len(want) == 588
    got = ncrna.read_set("prok-1")
    assert len(got) == 588 and len({m.acc for m in got}) == 588
    assert {m.acc: m.name for m in got} == want
    whole = {m.acc for m in got if m.whole_genome}
    assert len(whole) == 23 and whole == IN_CDS
    # every accession is a model of Rfam 15.1 with that name (the note's model table)
    table = {}
    for line in MODELS.read_text().splitlines()[1:]:
        f = line.split("\t")
        table[f[0]] = f[1]
    assert len(table) == 4227
    assert all(table.get(m.acc) == m.name for m in got)


def test_prok1_is_frozen_and_names_its_source():
    path = ncrna.set_path("prok-1")
    assert path == ncrna.SET_DIR / "prok-1.tsv"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == ncrna.FROZEN_SETS["prok-1"]
    assert ncrna.set_version("prok-1") == "prok1"
    head = "".join(ln for ln in path.read_text().splitlines(keepends=True) if ln.startswith("#"))
    assert "planning/design/0-08_ncrna_models.md" in head and "0-08_model_rules.py" in head
    assert re.search(r"at commit [0-9a-f]{40}", head) and re.search(r"20\d\d-\d\d-\d\d", head)


def test_a_frozen_set_that_changed_is_refused_and_another_file_is_named_by_its_sha256(tmp_path, monkeypatch):
    other = tmp_path / "mine.tsv"
    other.write_text("acc\tname\twhole_genome\nRF00005\ttRNA\t0\n")
    assert ncrna.set_version(other) == "sha" + hashlib.sha256(other.read_bytes()).hexdigest()[:12]
    assert [m.acc for m in ncrna.read_set(other)] == ["RF00005"]
    monkeypatch.setitem(ncrna.FROZEN_SETS, "prok-1", "0" * 64)
    with pytest.raises(ValueError, match="frozen model set"):
        ncrna.read_set("prok-1")
    for bad in ("acc\tname\n", "acc\tname\twhole_genome\nRF1\tx\t2\n", "acc\tname\twhole_genome\nRF1\tx\t0\nRF1\ty\t1\n"):
        other.write_text(bad)
        with pytest.raises(ValueError):
            ncrna.read_set(other)

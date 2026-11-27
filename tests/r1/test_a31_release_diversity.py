"""Card A-31: the release-diversity figure's data are traceable, add up, and draw the planned releases.

The data file is `planning/measurements/release_diversity.json`; the collector, checks and renderer are
`scripts/plot_release_diversity.py`. These tests read the tree only (no pod tools, no data root); the rendering case
skips without matplotlib.
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-31")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "plot_release_diversity.py"
DATA = REPO / "planning" / "measurements" / "release_diversity.json"
NOTE = REPO / "planning" / "measurements" / "2026-11-27_release_diversity.md"


def real_sources() -> Path:
    """The tree's sources.yaml; the test image (built from main, container/Dockerfile) does not ship it, so skip there."""
    p = REPO.joinpath("sources.yaml")
    if not p.is_file():
        pytest.skip("sources.yaml is not in the image: this check runs from a checkout (ci and the pod's mirror)")
    return p


def _mod():
    spec = importlib.util.spec_from_file_location("plot_release_diversity", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["plot_release_diversity"] = mod
    spec.loader.exec_module(mod)
    return mod


def _data() -> dict:
    return json.loads(DATA.read_text(encoding="utf-8"))


def test_every_number_has_a_label_and_a_source():
    m = _mod()
    nums = m.walk_numbers(_data())
    assert len(nums) > 80
    for where, num in nums:
        assert num["label"] in ("measured", "derived", "estimate"), where
        assert str(num.get("source", "")).strip(), where


def test_quotes_formulas_classes_and_releases_hold():
    real_sources()  # one quote (GBRMGD's archive bytes) is in sources.yaml
    assert _mod().data_problems(_data()) == []


def test_sums_by_class_equal_the_totals_and_each_stage_is_in_one_class_and_one_release():
    m = _mod()
    data = _data()
    cols = m.table(data)
    for rid, col in cols.items():
        assert sum(col["by_class"].values()) == pytest.approx(col["total"], abs=1e-9), rid
        assert sum(col["by_state"].values()) == pytest.approx(col["total"], abs=1e-9), rid
        assert col["low"] <= col["total"] + 1e-9 <= col["high"] + 2e-9, rid
    for r in data["releases"]:
        if r["status"] == "released":
            assert cols[r["id"]]["total"] * 1e9 == pytest.approx(r["stored_bases"]["value"], abs=0.5), r["id"]
    names = [s["stage"] for s in data["stages"]]
    assert len(names) == len(set(names))
    classes = {c["id"] for c in data["classes"]}
    releases = {r["id"] for r in data["releases"]}
    for s in data["stages"]:
        assert s["class"] in classes and s["release"] in releases, s["stage"]
    # a number that does not add up is found
    broken = json.loads(json.dumps(data))
    broken["releases"][0]["stored_bases"]["value"] += 1
    assert any("stages sum to" in p for p in m.data_problems(broken))


def test_measured_values_of_the_released_releases_equal_the_audit():
    data = _data()
    audit = data["audit"]
    rels = [r["id"] for r in data["releases"]]
    for r in data["releases"]:
        if r["status"] != "released":
            continue
        a = audit["releases"][r["id"]]
        assert r["stored_bases"]["value"] == a["stored_bases"] == sum(
            audit["stages"][s]["default_bases"] for s in a["stages"]), r["id"]
        assert r["proteins"]["value"] == a["proteins"], r["id"]
        assert a["cluster_member_rows"] == 2 * a["proteins"], r["id"]
        for lvl in ("nr50", "nr90"):
            if a[lvl] is not None:
                assert r[lvl]["value"] == a[lvl], (r["id"], lvl)
        drawn = sorted(s["stage"] for s in data["stages"] if rels.index(s["release"]) <= rels.index(r["id"]))
        assert drawn == sorted(a["stages"]), r["id"]
    for s in data["stages"]:
        if s["state"] in ("released", "built"):
            assert s["bases"]["value"] == audit["stages"][s["stage"]]["default_bases"], s["stage"]
    # every audited value is in the note, as the note's table writes it
    note = NOTE.read_text(encoding="utf-8")
    for st, a in audit["stages"].items():
        assert f"{a['default_bases']:,}" in note, st
    m = _mod()
    broken = json.loads(json.dumps(data))
    broken["releases"][4]["nr50"]["value"] += 1
    assert any("nr50" in p for p in m.data_problems(broken))


def _a26_chain() -> list[str]:
    text = (REPO / "planning/measurements/2026-11-24_a26_r26_11_4_confs.md").read_text(encoding="utf-8")
    line = next(ln for ln in text.splitlines() if ln.startswith("`V2 -> PL2"))
    return [s.strip() for s in line.strip("`").split("->")]


def _dec85_members() -> list[str]:
    text = (REPO / "planning/DECISIONS.md").read_text(encoding="utf-8")
    row = next(ln for ln in text.splitlines() if ln.startswith("| DEC-85 "))
    part = row.split("as they now stand:", 1)[1].split("; `ASG1`", 1)[0]
    return " ".join(re.findall(r"`([^`]+)`", part)).split()


@pytest.mark.covers("A-31", "A-32")
def test_the_planned_releases_hold_their_members_and_states():
    data = _data()
    by_rel: dict[str, list[dict]] = {}
    for s in data["stages"]:
        by_rel.setdefault(s["release"], []).append(s)
    r4 = by_rel["r26.11.4"]
    assert [s["stage"] for s in r4] == _a26_chain()
    assert len(r4) == 19
    built_r4 = [s["stage"] for s in r4 if s["state"] == "built"]
    assert {"V2", "PL2", "PL3", "CG1", "MGE1"} <= set(built_r4)
    assert [s["stage"] for s in r4[:len(built_r4)]] == built_r4  # wave-6 builds in strict chain prefix order
    assert all(s["state"] == "staged" for s in r4[len(built_r4):])
    r5 = by_rel["r26.11.5"]
    assert [s["stage"] for s in r5] == _dec85_members()
    assert len(r5) == 24
    assert all(s["state"] == "staged" for s in r5)  # RMD1 and GBRMGD staged by C-83 and promoted by A-32
    assert all(s["bases"]["label"] in ("measured", "derived") for s in r5)
    assert all(s["state"] == "planned" and s["bases"]["label"] == "estimate" for s in by_rel["r26.11.6"])
    assert data["r26_11_5_ownership"] == "planning/measurements/2026-11-27_release_diversity/r26_11_5_ownership.json"
    assert data["r26_11_5_upper_bound"] is False
    # the old r26.11.4 projection of NR50 families has no stated formula and is not drawn
    assert all("nr50" not in r for r in data["releases"] if r["status"] != "released")


@pytest.mark.covers("A-31", "A-32")
def test_a32_uniref_benchmarks_and_24_stage_ownership_are_verified():
    m = _mod()
    data = _data()
    own_rel = data["r26_11_5_ownership"]
    assert own_rel == "planning/measurements/2026-11-27_release_diversity/r26_11_5_ownership.json"
    own_path = REPO / own_rel
    assert own_path.is_file(), f"missing committed ownership JSON: {own_path}"
    own = json.loads(own_path.read_text(encoding="utf-8"))
    assert own["problems"] == [] and own["discrepancies"] == []
    assert own["chain"] == _dec85_members() and len(own["chain"]) == 24

    class Args:
        data = DATA
        ownership = None

    overlay = m._data_with_overlay(Args)
    assert m.data_problems(overlay) == []
    assert overlay["r26_11_5_upper_bound"] is False
    by_stage = {s["stage"]: s for s in overlay["stages"]}
    for st in ("RMD1", "GBRMGD"):
        assert by_stage[st]["state"] == "staged"
        assert by_stage[st]["bases"]["label"] == "measured"
        assert by_stage[st]["bases"]["value"] == own["stages"][st]["new_bases"]

    prefs = {r["id"]: r for r in data["benchmarks"]["protein_refs"]}
    for bid in ("uniref50", "uniref100"):
        ref = prefs[bid]
        assert ref.get("re_verified") is True, f"{bid} must be re-verified from UniProt release notes (DEC-87)"
        assert "citation" not in ref
        assert re.fullmatch(r"[0-9a-f]{64}", str(ref.get("sha256", ""))), bid
        assert re.fullmatch(r"[0-9a-f]{64}", str(ref.get("relnotes_sha256", ""))), bid

    if REPO.joinpath("sources.yaml").is_file():
        import yaml

        sdoc = yaml.safe_load(real_sources().read_text(encoding="utf-8"))
        entries = {e["name"]: e for e in sdoc["sources"] if isinstance(e, dict) and "name" in e}
        assert "uniprot_relnotes" in entries
        files = {f["file"]: f for f in entries["uniprot_relnotes"].get("files", [])}
        assert set(files) == {"relnotes.txt", "uniref50.release_note", "uniref100.release_note"}
        assert prefs["uniref50"]["sha256"] == files["uniref50.release_note"]["sha256"]
        assert prefs["uniref100"]["sha256"] == files["uniref100.release_note"]["sha256"]
        assert prefs["uniref50"]["relnotes_sha256"] == files["relnotes.txt"]["sha256"]
        assert prefs["uniref100"]["relnotes_sha256"] == files["relnotes.txt"]["sha256"]


def test_card_a29s_ownership_replaces_the_r26_11_5_new_bases(tmp_path):
    m = _mod()
    data = _data()
    own = {"stages": {"GRARE1": {"new_bases": 11_000_000_000}, "RMD1": {"new_bases": 33_000_000_000},
                      "V2": {"new_bases": 1}}}
    path = tmp_path / "own.json"
    path.write_text(json.dumps(own, indent=1, sort_keys=True), encoding="utf-8")
    out = m.apply_ownership(data, json.loads(path.read_text(encoding="utf-8")), str(path))
    st = {s["stage"]: s for s in out["stages"]}
    assert st["GRARE1"]["bases"]["value"] == 11_000_000_000 and st["GRARE1"]["bases"]["label"] == "measured"
    assert st["RMD1"]["state"] == "staged" and st["RMD1"]["bases"]["value"] == 33_000_000_000
    assert st["V2"]["bases"]["value"] == 8_199_089_073  # only r26.11.5 stages are replaced
    assert st["GRARE2"]["bases"] == {s["stage"]: s for s in data["stages"]}["GRARE2"]["bases"]
    assert out["r26_11_5_upper_bound"] is False
    assert m.table(out)["r26.11.5"]["total"] != m.table(data)["r26.11.5"]["total"]
    # the command reads it from the path the data file names
    data2 = json.loads(json.dumps(data))
    data2["r26_11_5_ownership"] = str(path)
    f = tmp_path / "data.json"
    f.write_text(json.dumps(data2), encoding="utf-8")

    class A:
        data = f
        ownership = None

    assert {s["stage"]: s for s in m._data_with_overlay(A)["stages"]}["GRARE1"]["bases"]["value"] == 11_000_000_000


def test_the_tower_remainder_is_read_from_sources_yaml(tmp_path):
    m = _mod()
    src = tmp_path / "sources.yaml"
    src.write_text(
        "version: 1\nsources:\n"
        "  - name: tower_a\n    files:\n      - {file: a.fna.gz, bytes: 100}\n      - {file: t.xlsx, bytes: 7}\n"
        "  - name: tower_b\n    files:\n      - {file: b.tar.gz, bytes: 50}\n"
        "  - name: tower_c\n    fetched: nothing yet\n"
        "  - name: other\n    files:\n      - {file: x.fna, bytes: 1000}\n",
        encoding="utf-8")
    rule = {"tower_remainder": {"assigned": {"tower_b": "B"}}}
    t = m.tower_remainder(rule, src)
    assert (t["sequence_bytes"], t["other_bytes"]) == (100, 7)
    assert t["landed_unassigned_entries"] == 1 and t["not_landed_unassigned_entries"] == 1
    real = m.tower_remainder(_data(), real_sources())
    assert real["unknown_assigned"] == []
    assert real["sequence_bytes"] > 0


FIGURE_TEXT_KEYS = ("figure_title", "footnote", "figure_label", "figure_text", "label", "novelty_caption", "upper",
                    "exact")


def _texts(obj, out):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in FIGURE_TEXT_KEYS and isinstance(v, str):
                out.append(v)
            _texts(v, out)
    elif isinstance(obj, list):
        for v in obj:
            _texts(v, out)
    return out


def test_figure_text_has_no_dashes_of_prose():
    m = _mod()
    texts = _texts(_data(), []) + list(m.STATE_TEXT.values())
    assert len(texts) > 20
    for t in texts:
        assert "—" not in t and "–" not in t, t


def test_render_writes_a_png_and_an_svg_without_a_date(tmp_path):
    pytest.importorskip("matplotlib")
    pytest.importorskip("yaml")
    m = _mod()
    data = _data()
    tower = m.tower_remainder(data, real_sources()) if REPO.joinpath("sources.yaml").is_file() else None
    png, svg = m.render(data, tmp_path / "a" / "fig", tower)
    assert png.stat().st_size > 50_000
    text = svg.read_text(encoding="utf-8")
    assert "<dc:date>" not in text
    assert "r26.11.6" in text and "Planned (estimate)" in text
    _, svg2 = m.render(data, tmp_path / "b" / "fig", tower)
    assert svg2.read_text(encoding="utf-8") == text

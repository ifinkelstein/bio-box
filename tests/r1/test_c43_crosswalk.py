"""Card C-43: `scripts/stress/c43_crosswalk.py`, one assembly per sequencing run across resources.

Small constructed tables with the truth written by hand: which run each resource holds, how a
sample or a WGS set is linked to its runs through the tables' own cross-reference columns, what
cannot be linked, and which runs only one resource holds. The oracle is the hand-written truth
below, not the function's own output.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-43")  # ring R1 by its directory (tests/conftest.py)

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "c43_crosswalk.py"


def _mod():
    spec = importlib.util.spec_from_file_location("c43_crosswalk", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    sys.modules["c43_crosswalk"] = m  # dataclasses resolve string annotations through sys.modules
    spec.loader.exec_module(m)
    return m


STUDIES_HTML = """
<script>
var studiesData = [
    { study: "PRJEB9584_gut", samples: 211 },

    { study: "Howe_2014_soil", samples: 2 },
    { study: "a b/c", samples: 1 },
];
new gridjs.Grid({ data: studiesData.map(function(d) { return [d.study, d.samples]; }) });
</script>
"""


def test_parse_studies_reads_the_embedded_list():
    m = _mod()
    assert m.parse_studies(STUDIES_HTML) == [("PRJEB9584_gut", 211), ("Howe_2014_soil", 2), ("a b/c", 1)]


def test_study_url_encodes_the_name():
    m = _mod()
    assert m.study_url("a b/c") == "https://spire.embl.de/spire/api/study/a%20b%2Fc?format=tsv"


def test_microntology_cell_is_parsed_as_data():
    m = _mod()
    assert m.parse_micron("['aquatic', 'aquatic:marine']") == ["aquatic", "aquatic:marine"]
    assert m.parse_micron("[]") == []
    assert m.parse_micron("") == []
    with pytest.raises(ValueError):
        m.parse_micron("__import__('os')")


@pytest.mark.parametrize("terms, label", [
    # sewage gets 'aquatic' pulled in; the more specific anthropogenic tag wins
    (["anthropogenic", "anthropogenic:wastewater", "anthropogenic:wastewater:sewage", "aquatic"],
     ("anthropogenic", "anthropogenic:wastewater")),
    # a mangrove: terrestrial before aquatic
    (["aquatic", "aquatic:littoral", "terrestrial", "terrestrial:wetland", "terrestrial:wetland:mangrove"],
     ("terrestrial", "terrestrial:wetland")),
    (["host-associated", "host-associated:animal host", "aquatic", "aquatic:marine"],
     ("host-associated", "host-associated:animal host")),
    (["aquatic", "aquatic:marine", "aquatic:marine:pelagic", "salinity:high salinity"], ("aquatic", "aquatic:marine")),
    (["air"], ("air", "air")),
    # only descriptors, no habitat
    (["salinity", "salinity:low salinity", "age group:adult"], ("unlabelled", "unlabelled")),
    ([], ("unlabelled", "unlabelled")),
])
def test_spire_habitat_by_stated_precedence(terms, label):
    m = _mod()
    assert m.spire_habitat(terms) == label


def test_mapped_label_longest_prefix_and_exact_names():
    m = _mod()
    table = {"root:Host-associated": "host-associated", "root:Environmental:Aquatic": "aquatic",
             "root:Environmental:Aquatic:Marine:Wetlands": "terrestrial", "marine metagenome": "aquatic"}
    assert m.mapped_label("root:Environmental:Aquatic:Marine", table) == "aquatic"
    assert m.mapped_label("root:Environmental:Aquatic:Marine:Wetlands:Sediment", table) == "terrestrial"
    assert m.mapped_label("root:Host-associated:Human:Skin", table) == "host-associated"
    # a lineage only shares a prefix at a ':' boundary
    assert m.mapped_label("root:Host-associatedX", table) == "unlabelled"
    assert m.mapped_label("marine metagenome", table) == "aquatic"
    assert m.mapped_label("marine metagenomes", table) == "unlabelled"
    assert m.mapped_label(None, table) == "unlabelled"


# A read-run table: run, primary sample, secondary sample, read bases
RUNS = [
    ("ERR1", "SAMEA1", "ERS1", 100),
    ("ERR2", "SAMEA1", "ERS1", 50),     # sample SAMEA1 has two runs
    ("SRR3", "SAMN3", "SRS3", 70),
    ("SRR4", "SAMN4", "SRS4", 30),
    ("SRR5", "SAMN5", "SRS5", 10),
    ("DRR6", "SAMD6", "DRS6", 5),
]


def test_link_by_run_then_by_either_sample_accession():
    m = _mod()
    idx = m.RunIndex(RUNS)
    # direct runs win and are kept even when the run table lacks them (bases unknown)
    assert idx.link(["SAMEA1"], ["ERR1", "ERR9"]) == ({"ERR1", "ERR9"}, {"SAMEA1"}, "run")
    # a sample given by its primary or its secondary accession reaches all of its runs
    assert idx.link(["SAMEA1"], []) == ({"ERR1", "ERR2"}, {"SAMEA1"}, "sample")
    assert idx.link(["ERS1"], []) == ({"ERR1", "ERR2"}, {"SAMEA1"}, "sample")
    # nothing in the table: unlinked, the sample kept as given
    assert idx.link(["SAMEA77"], []) == (set(), {"SAMEA77"}, "none")
    assert idx.link([], []) == (set(), set(), "none")
    assert idx.bases("ERR2") == 50
    assert idx.bases("ERR9") is None


def test_a_spire_sample_named_by_a_run_accession_is_that_run():
    m = _mod()
    idx = m.RunIndex(RUNS)
    assert m.spire_link("SRR3", idx) == ({"SRR3"}, set(), "run")
    assert m.spire_link("ERR1137728", idx) == ({"ERR1137728"}, set(), "run")
    assert m.spire_link("SAMEA1", idx) == ({"ERR1", "ERR2"}, {"SAMEA1"}, "sample")
    # a study's own sample name is not guessed at
    assert m.spire_link("MH0002", idx) == (set(), {"MH0002"}, "none")
    assert m.spire_link("SRR3x", idx) == (set(), {"SRR3x"}, "none")


def test_insdc_samples_for_lookup():
    m = _mod()
    assert m.insdc_samples(["SAMEA104389540", "SAMN02954240", "SAMD00008643", "ERS1", "SRS927152", "DRS6",
                            "MH0002", "G80150", "SAMX", ""]) == [
        "DRS6", "ERS1", "SAMD00008643", "SAMEA104389540", "SAMN02954240", "SRS927152"]


def test_lookup_query_names_each_sample_by_its_field():
    m = _mod()
    assert m.lookup_query(["ERS1", "SAMEA1", "SRS9"]) == (
        'secondary_sample_accession="ERS1" OR sample_accession="SAMEA1" OR secondary_sample_accession="SRS9"')
    url = m.lookup_url(["SAMN5"])
    assert url.startswith("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=")
    assert "query=sample_accession%3D%22SAMN5%22&fields=run_accession," in url
    assert url.endswith("&limit=0&format=tsv")


@pytest.mark.parametrize("atype, include, keep", [
    ("primary metagenome", {"primary metagenome"}, True),
    ("Primary Metagenome", {"primary metagenome"}, True),
    ("binned metagenome", {"primary metagenome"}, False),
    ("Metagenome-Assembled Genome (MAG)", {"primary metagenome"}, False),
    ("", {"primary metagenome"}, False),
    ("", {"primary metagenome", ""}, True),
])
def test_wgs_type_filter(atype, include, keep):
    m = _mod()
    assert m.wgs_included(atype, include) is keep


def _u(m, resource, uid, runs=(), samples=(), bases=None, habitat="unlabelled"):
    return m.Unit(resource=resource, uid=uid, runs=frozenset(runs), samples=frozenset(samples),
                  linked_by="run" if runs else ("sample" if samples else "none"), asm_bases=bases, habitat=habitat)


def test_choose_one_assembly_per_run_in_order_of_preference():
    m = _mod()
    units = [
        _u(m, "spire", "SAMEA1", runs={"ERR1", "ERR2"}, samples={"SAMEA1"}),
        _u(m, "spire", "SAMD6", samples={"SAMD6"}),                       # SPIRE sample with no run in the table
        _u(m, "mgnify", "ERZ1", runs={"ERR1"}),                           # same run as SPIRE: covered
        _u(m, "mgnify", "ERZ3", runs={"SRR3"}),                           # new
        _u(m, "mgnify", "ERZ4", runs={"SRR3", "SRR4"}),                   # partial: SRR3 already taken
        _u(m, "mgnify", "ERZ6", samples={"SAMD6"}),                       # no runs; sample matches SPIRE: covered
        _u(m, "wgs", "AAAA01", runs={"SRR3"}),                            # covered by ERZ3
        _u(m, "wgs", "BBBB01", runs={"SRR5"}),                            # new
        _u(m, "wgs", "CCCC01", runs={"SRR5"}),                            # second assembly of SRR5: covered
        _u(m, "wgs", "DDDD01", samples={"SAMN99"}),                       # unlinked and unmatched: taken, flagged
        _u(m, "wgs", "EEEE01", runs={"DRR6"}, samples={"SAMD6"}),         # run new, but its sample SPIRE assembled
    ]
    got = m.choose(units, ["spire", "mgnify", "wgs"])
    assert got == {
        "SAMEA1": "take", "SAMD6": "take_unlinked",
        "ERZ1": "covered", "ERZ3": "take", "ERZ4": "partial", "ERZ6": "covered",
        "AAAA01": "covered", "BBBB01": "take", "CCCC01": "covered", "DDDD01": "take_unlinked",
        "EEEE01": "covered",
    }


def test_choose_is_independent_of_input_order():
    m = _mod()
    units = [_u(m, "wgs", "CCCC01", runs={"SRR5"}), _u(m, "wgs", "BBBB01", runs={"SRR5"}),
             _u(m, "mgnify", "ERZ5", runs={"SRR5", "SRR4"})]
    a = m.choose(units, ["spire", "mgnify", "wgs"])
    b = m.choose(list(reversed(units)), ["spire", "mgnify", "wgs"])
    assert a == b == {"ERZ5": "take", "BBBB01": "covered", "CCCC01": "covered"}


def test_run_resources_and_exclusive_counts():
    m = _mod()
    idx = m.RunIndex(RUNS)
    units = [
        _u(m, "spire", "SAMEA1", runs={"ERR1", "ERR2"}, samples={"SAMEA1"}, habitat="aquatic"),
        _u(m, "mgnify", "ERZ1", runs={"ERR1"}, habitat="host-associated"),
        _u(m, "mgnify", "ERZ3", runs={"SRR3"}, habitat="soil-like"),
        _u(m, "wgs", "AAAA01", runs={"SRR3", "SRR4"}, habitat="terrestrial"),
        _u(m, "wgs", "BBBB01", runs={"ERR9"}, habitat="aquatic"),         # run not in the run table
    ]
    per_run = m.run_resources(units, ["spire", "mgnify", "wgs"])
    assert per_run == {
        "ERR1": (("spire", "mgnify"), "aquatic"),          # SPIRE's label first
        "ERR2": (("spire",), "aquatic"),
        "SRR3": (("mgnify", "wgs"), "soil-like"),          # then the earlier resource's label
        "SRR4": (("wgs",), "terrestrial"),
        "ERR9": (("wgs",), "aquatic"),
    }
    table = m.combo_table(per_run, idx)
    assert table[("spire",)]["aquatic"] == {"runs": 1, "read_bases": 50, "runs_without_bases": 0}
    assert table[("spire", "mgnify")]["aquatic"] == {"runs": 1, "read_bases": 100, "runs_without_bases": 0}
    assert table[("mgnify", "wgs")]["soil-like"] == {"runs": 1, "read_bases": 70, "runs_without_bases": 0}
    assert table[("wgs",)]["terrestrial"] == {"runs": 1, "read_bases": 30, "runs_without_bases": 0}
    assert table[("wgs",)]["aquatic"] == {"runs": 1, "read_bases": 0, "runs_without_bases": 1}


def _w(path, text):
    path.write_text(text)
    return str(path)


def test_build_end_to_end_on_small_tables(tmp_path):
    """Every reader and the outputs, on tables shaped as the providers give them."""
    import csv
    import json
    m = _mod()
    spire = tmp_path / "spire"
    spire.mkdir()
    _w(spire / "S1_gut.tsv", "sample_id\tmags\tlat\tlon\tmicrontology\n"
       "SAMEA1\t3\t\t\t['host-associated', 'host-associated:animal host']\n"
       "SAMD6\t0\t\t\t[]\n")
    _w(spire / "probe_S1_gut.tsv", "sample_id\tmags\tlat\tlon\tmicrontology\nSAMX\t0\t\t\t[]\n")
    rr = _w(tmp_path / "rr.tsv", "run_accession\tsample_accession\tsecondary_sample_accession\tstudy_accession\t"
            "secondary_study_accession\tbase_count\tlibrary_strategy\tinstrument_platform\ttax_id\tscientific_name\t"
            "first_public\n"
            "ERR1\tSAMEA1\tERS1\tP\tE\t100\tWGS\tI\t1\tx\t2020\n"
            "ERR2\tSAMEA1\tERS1\tP\tE\t50\tWGS\tI\t1\tx\t2020\n"
            "SRR3\tSAMN3\tSRS3\tP\tE\t70\tWGS\tI\t1\tx\t2020\n"
            "SRR4\tSAMN4\tSRS4\tP\tE\t30\tWGS\tI\t1\tx\t2020\n"
            "SRR5\tSAMN5\tSRS5\tP\tE\t10\tWGS\tI\t1\tx\t2020\n")
    # ENA's assembly table: MGnify's own assembly has no WGS set
    wgs = _w(tmp_path / "asm_ena.tsv", "assembly_accession\tassembly_type\tbase_count\trun_accession\tsample_accession\t"
             "secondary_sample_accession\tstudy_accession\twgs_set\tcenter_name\ttax_id\tscientific_name\tfirst_created\n"
             "GCA_900000001\tprimary metagenome\t900\t\tSAMN3\t\tP\t\tEMG\t1\tsoil metagenome\t2020\n"
             "GCA_000000002\tprimary metagenome\t800\t\tSAMN5\t\tP\tBBBB01\tX\t1\tmarine metagenome\t2020\n"
             "GCA_000000004\tprimary metagenome\t700\tERR2\tSAMEA1\t\tP\tDDDD01\tY\t1\tmetagenome\t2020\n"
             "GCA_000000005\tbinned metagenome\t3\tSRR5\tSAMN5\t\tP\tEEEE01\tZ\t1\tmetagenome\t2020\n")
    ncbi = _w(tmp_path / "asm.txt", "##  See README\n#assembly_accession\tbioproject\tbiosample\torganism_name\t"
              "genome_size\n"
              "GCA_000000002.1\tP\tSAMN5\tmarine metagenome\t800\n"
              "GCA_000000009.1\tP\tSAMN99\tfeces metagenome\t600\n")
    mga = tmp_path / "mga"
    mga.mkdir()
    page = {"links": {"next": None}, "data": [
        {"id": "ERZ3", "attributes": {"is-private": False, "legacy-accession": "GCA_900000001", "wgs-accession": None},
         "relationships": {"runs": {"data": [{"id": "SRR3"}]}, "samples": {"data": [{"id": "SRS3"}]}}},
        {"id": "ERZ1", "attributes": {"is-private": False, "legacy-accession": None, "wgs-accession": None},
         "relationships": {"runs": {"data": []}, "samples": {"data": [{"id": "ERS1"}]}}},
    ]}
    _w(mga / "page_000001.json", json.dumps(page))
    mgs = tmp_path / "mgs"
    mgs.mkdir()
    _w(mgs / "page_000001.json", json.dumps({"links": {"next": None}, "data": [
        {"id": "SRS3", "attributes": {"biosample": "SAMN3"},
         "relationships": {"biome": {"data": {"id": "root:Environmental:Terrestrial:Soil"}}}}]}))
    # ENA's analyses by MGnify: file sizes of its contig files; ERZ9 is MGnify's but not in its API listing
    emg = _w(tmp_path / "emg.tsv", "analysis_accession\tassembly_type\tsubmitted_bytes\tgenerated_bytes\t"
             "sample_accession\trun_accession\tstudy_accession\tfirst_public\n"
             "ERZ3\tprimary metagenome\t40;2\t\tSAMN3\tSRR3\tP\t2020\n"
             "ERZ9\tprimary metagenome\t11\t\tSAMN4\tSRR4\tP\t2020\n"
             "ERZ8\tbinned metagenome\t5\t\tSAMN4\tSRR4\tP\t2020\n")
    hm = _w(tmp_path / "map.tsv", "kind\tkey\tlabel\nmgnify_lineage\troot:Environmental:Terrestrial\tterrestrial\n"
            "organism_name\tmarine metagenome\taquatic\norganism_name\tfeces metagenome\thost-associated\n")
    out = tmp_path / "out"
    assert m.main(["build", "--spire-dir", str(spire), "--read-run", rr, "--insdc-assemblies", wgs,
                   "--ncbi-summary", ncbi,
                   "--mgnify-assemblies", str(mga), "--mgnify-samples", str(mgs), "--mgnify-ena-analyses", emg,
                   "--habitat-map", hm, "--out-dir", str(out)]) == 0

    def rows(name):
        with open(out / name) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    units = {r["uid"]: r for r in rows("units.tsv")}
    # the probe file is not a study; GCA_900000001 is MGnify's ERZ3 and not a unit of its own
    assert set(units) == {"SAMEA1", "SAMD6", "ERZ3", "ERZ1", "ERZ9", "GCA_000000002", "GCA_000000004",
                          "GCA_000000009"}
    assert units["ERZ3"]["gz_bytes"] == "42" and units["ERZ1"]["gz_bytes"] == ""
    assert units["ERZ9"]["runs"] == "SRR4" and units["ERZ9"]["gz_bytes"] == "11"
    assert units["ERZ9"]["decision"] == "take" and units["ERZ9"]["note"] == "source=ena_analysis_only"
    assert units["SAMEA1"]["runs"] == "ERR1;ERR2" and units["SAMEA1"]["habitat"] == "host-associated"
    assert units["SAMD6"]["decision"] == "take_unlinked" and units["SAMD6"]["linked_by"] == "none"
    assert units["ERZ1"]["decision"] == "covered" and units["ERZ1"]["linked_by"] == "sample"
    assert units["ERZ3"]["decision"] == "take" and units["ERZ3"]["asm_bases"] == "900"
    assert units["ERZ3"]["habitat"] == "terrestrial"
    g2, g4 = units["GCA_000000002"], units["GCA_000000004"]
    assert g2["decision"] == "take" and g2["habitat"] == "aquatic"
    assert g2["note"] == "ncbi_assembly=yes;wgs_set=BBBB01;center=X"
    assert g4["decision"] == "covered" and g4["habitat"] == "unlabelled"
    assert units["GCA_000000009"]["decision"] == "take_unlinked" and units["GCA_000000009"]["asm_bases"] == "600"
    res = {r["resource"]: r for r in rows("resources.tsv")}
    assert (res["spire"]["units"], res["spire"]["unlinked"], res["spire"]["read_bases"]) == ("2", "1", "150")
    assert (res["insdc"]["units"], res["insdc"]["distinct_runs"], res["insdc"]["asm_bases"]) == ("3", "2", "2100")
    assert (res["mgnify"]["units"], res["mgnify"]["units_with_gz_bytes"], res["mgnify"]["gz_bytes"]) == ("3", "2", "53")
    combos = {(r["resources"], r["habitat"]): r["runs"] for r in rows("runs_by_combination_and_habitat.tsv")}
    assert combos == {("spire+mgnify", "host-associated"): "1", ("spire+mgnify+insdc", "host-associated"): "1",
                      ("mgnify", "terrestrial"): "1", ("mgnify", "unlabelled"): "1", ("insdc", "aquatic"): "1"}
    assert [r["uid"] for r in rows("take_insdc.tsv")] == ["GCA_000000002", "GCA_000000009"]
    # a bin is not an assembly of the run; MGnify's own INSDC record is the MGnify unit
    notu = {(r["why"], r["assembly_type"]): (r["assemblies"], r["bases"]) for r in rows("insdc_not_units.tsv")}
    assert notu == {("mgnify's own", "primary metagenome"): ("1", "900"),
                    ("type not taken", "binned metagenome"): ("1", "3")}

"""Card C-74 (DEC-59, DEC-79): the accession index of GTDB R232 and GlobDB r232, and the r226-to-r232 re-key.

Every test builds small fixture tables whose truth is written out by hand in this file (constructed truth, not a
snapshot of the script's output). Each test runs twice, with the normal PATH and with PATH reduced to `/usr/bin:/bin`
(as tests/r1/test_relay_to_pod.py does): the index needs no tool beyond Python and pyarrow, and GitHub Actions has none
of the pod's tools.
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import subprocess
import sys
import tarfile
from pathlib import Path

import pyarrow.parquet as pq
import pytest

pytestmark = pytest.mark.covers("C-74")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "stress" / "build_stage_pod" / "r232_index.py"


def _load():
    spec = importlib.util.spec_from_file_location("r232_index", SCRIPT)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["r232_index"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def _gz(text: str) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as g:
        g.write(text.encode())
    return buf.getvalue()


def _tar(members: dict[str, str]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as t:
        for name, text in members.items():
            data = text.encode()
            info = tarfile.TarInfo("dicts/" + name)
            info.size = len(data)
            t.addfile(info, io.BytesIO(data))
    return buf.getvalue()


META_HEAD = ["accession", "gtdb_genome_representative", "gtdb_representative", "gtdb_taxonomy",
             "ncbi_genbank_assembly_accession", "genome_size"]


def _meta(rows: list[tuple]) -> bytes:
    lines = ["\t".join(META_HEAD)]
    lines += ["\t".join(r) for r in rows]
    return _gz("\n".join(lines) + "\n")


def _tax(domain: str, species: str) -> str:
    return f"d__{domain};p__P;c__C;o__O;f__F;g__G;s__{species}"


@pytest.fixture()
def world(tmp_path: Path) -> dict:
    gt = tmp_path / "GTDB_R232_META"
    gl = tmp_path / "GLOBDB_R232_TABLES"
    r6 = tmp_path / "r226"
    for d in (gt, gl, r6):
        d.mkdir()
    # GTDB R232: 5 genomes, 3 species, 3 representatives
    (gt / "ar53_metadata_r232.tsv.gz").write_bytes(_meta([
        ("GB_GCA_000000011.1", "GB_GCA_000000011.1", "t", _tax("Archaea", "Archaeum alpha"), "GCA_000000011.1", "1000"),
        ("RS_GCF_000000022.1", "GB_GCA_000000011.1", "f", _tax("Archaea", "Archaeum alpha"), "GCA_000000022.1", "1100"),
    ]))
    (gt / "bac120_metadata_r232.tsv.gz").write_bytes(_meta([
        ("RS_GCF_000000033.2", "RS_GCF_000000033.2", "t", _tax("Bacteria", "Bacterium beta"), "GCA_000000033.2", "2000"),
        ("GB_GCA_000000044.1", "RS_GCF_000000033.2", "f", _tax("Bacteria", "Bacterium beta"), "GCA_000000044.1", "2100"),
        ("GB_GCA_000000055.1", "GB_GCA_000000055.1", "t", _tax("Bacteria", "Bacterium gamma"), "GCA_000000055.1", "3000"),
    ]))
    (gt / "MD5SUM.txt").write_text("")
    # GlobDB r232: 3 GTDB representatives and 6 others
    ids = {
        "GCA_000000011": _tax("Archaea", "Archaeum alpha"),
        "GCF_000000033": _tax("Bacteria", "Bacterium beta"),
        "GCA_000000055": _tax("Bacteria", "Bacterium gamma"),
        "MOTU40_000001": _tax("Bacteria", "Mot one"),
        "MOTU40_000002": _tax("Archaea", "Mot two"),
        "BCRBG_00001": _tax("Bacteria", "Bin one"),
        "SPIREOTU_00000010": _tax("Bacteria", "Spire ten"),
        "SPECIV4_00000": _tax("Bacteria", "Spec zero"),
        "MGYG000000001": _tax("Bacteria", "Mgy one"),
    }
    (gl / "globdb_r232_taxonomy.tsv.gz").write_bytes(_gz("".join(f"{k}\t{v}\n" for k, v in ids.items())))
    ck = ["Name\tCompleteness\tContamination\tGenome_Size\tTotal_Contigs"]
    ck += [f"{k}\t90.5\t1.25\t{1000 + i}\t{i + 1}" for i, k in enumerate(ids)]
    (gl / "globdb_r232_checkm2.tsv.gz").write_bytes(_gz("\n".join(ck) + "\n"))
    (gl / "globdb_r232_dataset_list.tsv").write_text(
        "Name\tIdentifier\tTitle\tPaper link\tSpecies reps in dataset\tNumber of genomes in GlobDB r232\n"
        "GTDB\tGCA & GCF\tt\tl\t3\t3\n"
        "mOTU\tMOTU40\tt\tl\t9\t2\n"
        "SPIRE\tSPIREOTU & SPECIV4\tt\tl\t9\t2\n"
        "RBG\tBCRBG\tt\tl\t1\t1\n"
        "MGnify\tMGYG\tt\tl\t1\t1\n"
        "AMXMAG\tAMXMAG\tt\tl\t1\t0\n")
    (gl / "globdb_r232_dictionaries.tar.gz").write_bytes(_tar({
        "MOTU4_dictionary.tsv": "orgA\tMOTU40_000001\norgB\tMOTU40_000002\norgC\tMOTU40_000003\n",
        "BCRBG_dictionary.tsv": "binchicken_1.fna\tBCRBG_00001.fa\n",
        "SPIREOTU_dictionary.tsv": "GlobDB_ID\tspire_cluster\tgenome_id\n"
                                   "SPIREOTU_00000010\tspire_c1\tspire_mag_10\n"
                                   "SPECIV4_00000\tspecI_v4_00000 GCA_000000044.1\n",
        "AMXMAG_dictionary.tsv": "GCA_000000099.1\tAMXMAG_0001\n",
        "README": "not a dictionary\n",
    }))
    # r226
    r226_ids = ["MOTU40_000009", "MOTU40_000005", "BCRBG_00001", "SPIREOTU_00000011", "AMXMAG_0007", "GFSMAG_0001",
                "MGYG000000001", "GCA_000000077", "GCA_000000044", "SPECIV4_00000"]
    (r6 / "globdb_r226_taxonomy.tsv.gz").write_bytes(_gz("".join(f"{k}\td__X;s__x\n" for k in r226_ids)))
    (r6 / "globdb_r226_dictionaries.tar.gz").write_bytes(_tar({
        "02_MOTU4_dictionary.tsv": "orgA\tMOTU40_000009\norgC\tMOTU40_000005\n",
        "04_BCRBG_dictionary.tsv": "binchicken_1\tBCRBG_00001\n",
        "03_SPIREOTU_dictionary.tsv": "GlobDB_ID\tspire_cluster\tgenome_id\n"
                                     "SPIREOTU_00000011\tspire_c1\tspire_mag_11\n"
                                     "SPECIV4_00000\tspecI_v4_00000 GCA_000000044.1\n",
        "13_AMXMAG_dictionary.tsv": "GCA_000000022.1\tAMXMAG_0007\n",
        "14_GFSMAG_dictionary.tsv": "gfs_1\tGFSMAG_0001\n",
    }))
    return {"tmp": tmp_path, "gt": gt, "gl": gl, "r6": r6, "out": tmp_path / "index", "r226_ids": r226_ids}


def test_gtdb_r232_maps_every_genome_to_accession_species_and_representative(world):
    m = _load()
    m.build(world["gl"], world["gt"], world["r6"], world["out"])
    t = pq.read_table(world["out"] / "gtdb_r232.parquet").to_pylist()
    got = {r["accession"]: (r["species"], r["is_representative"], r["representative"], r["domain"]) for r in t}
    assert got == {
        "GCA_000000011": ("Archaeum alpha", True, "GCA_000000011", "Archaea"),
        "GCF_000000022": ("Archaeum alpha", False, "GCA_000000011", "Archaea"),
        "GCF_000000033": ("Bacterium beta", True, "GCF_000000033", "Bacteria"),
        "GCA_000000044": ("Bacterium beta", False, "GCF_000000033", "Bacteria"),
        "GCA_000000055": ("Bacterium gamma", True, "GCA_000000055", "Bacteria"),
    }
    s = json.loads((world["out"] / "summary.json").read_text())
    assert (s["gtdb_genomes"], s["gtdb_species"], s["gtdb_representatives"]) == (5, 3, 3)


def test_globdb_r232_maps_genome_to_resource_taxonomy_and_insdc_accession(world):
    m = _load()
    m.build(world["gl"], world["gt"], world["r6"], world["out"])
    t = pq.read_table(world["out"] / "globdb_r232.parquet").to_pylist()
    got = {r["globdb_id"]: (r["resource"], r["is_gtdb"], r["insdc_accession"], r["original_id"]) for r in t}
    assert got == {
        "GCA_000000011": ("GTDB", True, "GCA_000000011", None),
        "GCF_000000033": ("GTDB", True, "GCF_000000033", None),
        "GCA_000000055": ("GTDB", True, "GCA_000000055", None),
        "MOTU40_000001": ("mOTU", False, None, "orgA"),
        "MOTU40_000002": ("mOTU", False, None, "orgB"),
        "BCRBG_00001": ("RBG", False, None, "binchicken_1"),
        "SPIREOTU_00000010": ("SPIRE", False, None, "spire_mag_10"),
        "SPECIV4_00000": ("SPIRE", False, "GCA_000000044", "GCA_000000044.1"),
        "MGYG000000001": ("MGnify", False, None, None),
    }
    row = next(r for r in t if r["globdb_id"] == "MOTU40_000002")
    assert row["taxonomy"].startswith("d__Archaea") and row["domain"] == "Archaea"
    assert (row["completeness"], row["contamination"]) == (90.5, 1.25)
    s = json.loads((world["out"] / "summary.json").read_text())
    assert s["globdb_genomes"] == 9
    assert s["globdb_by_resource"] == {"GTDB": 3, "MGnify": 1, "RBG": 1, "SPIRE": 2, "mOTU": 2}
    # the provider's list says 3, 2, 2, 1, 1 and 0 genomes: AMXMAG has 0 here and no row, so it differs by nothing
    assert s["dataset_list_differences"] == {}


def test_a_count_that_differs_from_the_dataset_list_is_reported(world):
    m = _load()
    p = world["gl"] / "globdb_r232_dataset_list.tsv"
    p.write_text(p.read_text().replace("mOTU\tMOTU40\tt\tl\t9\t2", "mOTU\tMOTU40\tt\tl\t9\t5"))
    m.build(world["gl"], world["gt"], world["r6"], world["out"])
    s = json.loads((world["out"] / "summary.json").read_text())
    assert s["dataset_list_differences"] == {"mOTU": {"listed": 5, "found": 2}}


def test_lookup_returns_membership_and_representative(world):
    m = _load()
    m.build(world["gl"], world["gt"], world["r6"], world["out"])
    # versioned, prefixed, GenBank-form of a RefSeq accession, and unknown inputs are all accepted
    res = m.lookup(world["out"], ["GCA_000000011.1", "RS_GCF_000000033.2", "GCA_000000033.2", "GCA_000000044.1",
                                  "GCA_000000022.1", "GCA_999999999.1", "GCA_000000099.1"])
    by = {r["accession"]: r for r in res}
    assert [r["accession"] for r in res] == ["GCA_000000011.1", "RS_GCF_000000033.2", "GCA_000000033.2",
                                              "GCA_000000044.1", "GCA_000000022.1", "GCA_999999999.1",
                                              "GCA_000000099.1"]
    a = by["GCA_000000011.1"]
    assert (a["in_gtdb_r232"], a["in_globdb_r232"], a["representative"], a["is_representative"]) == (
        True, True, "GCA_000000011", True)
    b = by["RS_GCF_000000033.2"]
    assert (b["in_gtdb_r232"], b["in_globdb_r232"], b["representative"]) == (True, True, "GCF_000000033")
    assert by["GCA_000000033.2"]["in_gtdb_r232"] is True and by["GCA_000000033.2"]["in_globdb_r232"] is True
    # a non-representative GTDB genome is in GTDB R232 and, by its INSDC accession, in GlobDB r232 (SPECIV4_00000)
    c = by["GCA_000000044.1"]
    assert (c["in_gtdb_r232"], c["in_globdb_r232"], c["representative"], c["is_representative"]) == (
        True, True, "GCF_000000033", False)
    assert c["globdb_r232_id"] == "SPECIV4_00000"
    d = by["GCA_000000022.1"]
    assert (d["in_gtdb_r232"], d["in_globdb_r232"], d["representative"]) == (True, False, "GCA_000000011")
    e = by["GCA_999999999.1"]
    assert (e["in_gtdb_r232"], e["in_globdb_r232"], e["representative"]) == (False, False, None)
    # an INSDC accession that only a dictionary names is not in GlobDB r232 unless the genome is in the taxonomy
    assert by["GCA_000000099.1"]["in_globdb_r232"] is False


def test_rekey_returns_the_r232_identifier_replaced_by_or_dropped(world):
    m = _load()
    m.build(world["gl"], world["gt"], world["r6"], world["out"])
    res = {r["r226_id"]: r for r in m.rekey(world["out"], world["r226_ids"] + ["NOT_IN_R226_1"])}
    # same original under a new identifier
    assert (res["MOTU40_000009"]["status"], res["MOTU40_000009"]["r232_id"]) == ("kept", "MOTU40_000001")
    # the r232 dictionary lists the genome but dereplication left it out of GlobDB r232
    assert (res["MOTU40_000005"]["status"], res["MOTU40_000005"]["r232_id"]) == ("dropped", None)
    assert res["MOTU40_000005"]["reason"] == "dictionary_lists_it_but_not_in_r232"
    # the .fna suffix of r232's dictionary does not hide the genome
    assert (res["BCRBG_00001"]["status"], res["BCRBG_00001"]["r232_id"]) == ("kept", "BCRBG_00001")
    # the genome is gone, its cluster id is now held by another genome
    assert (res["SPIREOTU_00000011"]["status"], res["SPIREOTU_00000011"]["r232_id"]) == (
        "replaced", "SPIREOTU_00000010")
    assert res["SPIREOTU_00000011"]["reason"] == "same_cluster_id"
    # an INSDC original that is now a non-representative GTDB R232 genome: replaced by the representative
    assert (res["AMXMAG_0007"]["status"], res["AMXMAG_0007"]["r232_id"]) == ("replaced", "GCA_000000011")
    assert res["AMXMAG_0007"]["reason"] == "insdc_accession_in_gtdb_r232"
    assert (res["GFSMAG_0001"]["status"], res["GFSMAG_0001"]["reason"]) == ("dropped", "resource_not_in_r232")
    # no dictionary: the identifier is the key
    assert (res["MGYG000000001"]["status"], res["MGYG000000001"]["r232_id"]) == ("kept", "MGYG000000001")
    assert (res["GCA_000000077"]["status"], res["GCA_000000077"]["r232_id"]) == ("dropped", None)
    assert (res["GCA_000000044"]["status"], res["GCA_000000044"]["r232_id"]) == ("replaced", "GCF_000000033")
    assert (res["SPECIV4_00000"]["status"], res["SPECIV4_00000"]["r232_id"]) == ("kept", "SPECIV4_00000")
    # an identifier r226 never had is refused as unknown, not guessed
    assert res["NOT_IN_R226_1"]["status"] == "unknown"


def test_an_identifier_reused_for_another_genome_is_not_taken_as_kept(world):
    m = _load()
    # r226 says MOTU40_000001 is orgZ; r232 says it is orgA
    (world["r6"] / "globdb_r226_taxonomy.tsv.gz").write_bytes(_gz("MOTU40_000001\td__X;s__x\n"))
    (world["r6"] / "globdb_r226_dictionaries.tar.gz").write_bytes(_tar({"02_MOTU4_dictionary.tsv": "orgZ\tMOTU40_000001\n"}))
    m.build(world["gl"], world["gt"], world["r6"], world["out"])
    r = m.rekey(world["out"], ["MOTU40_000001"])[0]
    assert r["status"] == "dropped" and r["reason"] == "original_not_in_r232"


def test_unbuilt_r226_selection_leaves_out_g1_and_the_gart_phyla(world):
    m = _load()
    rows = [
        ("ID", "stage", "resource", "taxonomy"),
        ("A1", "G1", "RBG", "d__Archaea;p__Aenigmatarchaeota;c__c"),
        ("B1", "G5", "RBG", "d__Bacteria;p__Patescibacteriota;c__c"),
        ("B2", "G6", "mOTU", "d__Bacteria;p__Planctomycetota;c__c"),
        ("B3", "G7", "mOTU", "d__Bacteria;p__Chloroflexota;c__c"),
        ("B4", "G8", "SPIRE", "d__Bacteria;p__Myxococcota_A;c__c"),
        ("B5", "G9", "SPIRE", "d__Archaea;p__Thermoproteota;c__c"),
    ]
    p = world["tmp"] / "non_gtdb_stages.tsv"
    p.write_text("\n".join("\t".join(r) for r in rows) + "\n")
    assert m.unbuilt_ids(p) == ["B3", "B5"]


def test_the_index_is_sorted_has_key_statistics_and_is_byte_identical_when_rebuilt(world):
    m = _load()
    m.build(world["gl"], world["gt"], world["r6"], world["out"])
    first = {f.name: f.read_bytes() for f in sorted(world["out"].iterdir())}
    assert {"gtdb_r232.parquet", "globdb_r232.parquet", "globdb_r226_keys.parquet", "summary.json"} <= set(first)
    for name, col in (("gtdb_r232.parquet", "accession"), ("globdb_r232.parquet", "globdb_id"),
                      ("globdb_r226_keys.parquet", "r226_id")):
        pf = pq.ParquetFile(world["out"] / name)
        vals = pf.read().column(col).to_pylist()
        assert vals == sorted(vals), name
        st = pf.metadata.row_group(0).column(pf.schema_arrow.get_field_index(col)).statistics
        assert st.min == vals[0] and st.max == vals[-1]
    out2 = world["tmp"] / "index2"
    # another input order of the dictionary members and a rebuilt directory must give the same bytes
    m.build(world["gl"], world["gt"], world["r6"], out2)
    second = {f.name: f.read_bytes() for f in sorted(out2.iterdir())}
    assert first == second


def test_the_command_line_builds_looks_up_and_rekeys(world):
    out = world["out"]
    r = subprocess.run([sys.executable, "-I", str(SCRIPT), "build", str(world["gl"]), str(world["gt"]),
                        str(world["r6"]), str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    ids = world["tmp"] / "ids.txt"
    ids.write_text("GCA_000000011.1\nGCA_999999999.1\n")
    r = subprocess.run([sys.executable, "-I", str(SCRIPT), "lookup", str(out), str(ids)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lines = r.stdout.splitlines()
    assert lines[0].split("\t")[0] == "accession" and lines[1].split("\t")[:3] == ["GCA_000000011.1", "True", "True"]
    ids.write_text("MOTU40_000009\n")
    r = subprocess.run([sys.executable, "-I", str(SCRIPT), "rekey", str(out), str(ids)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.splitlines()[1].split("\t")[:3] == ["MOTU40_000009", "kept", "MOTU40_000001"]


def test_a_missing_input_file_is_refused(world):
    (world["gt"] / "ar53_metadata_r232.tsv.gz").unlink()
    r = subprocess.run([sys.executable, "-I", str(SCRIPT), "build", str(world["gl"]), str(world["gt"]),
                        str(world["r6"]), str(world["out"])], capture_output=True, text=True)
    assert r.returncode == 2 and "ar53_metadata_r232.tsv.gz" in r.stderr

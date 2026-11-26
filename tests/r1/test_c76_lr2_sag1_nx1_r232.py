"""Card C-76 (DEC-79): the R232 / r232 status of the staged LR2, SAG1 and NX1, and their conf templates.

No genome is removed because r232 lists it (the orchestrating session, 2026-11-26): the status is recorded, the
stage's files stay byte for byte what cards C-66 and C-58 wrote. Truth is written by hand in this file: which fixture
genome sits in which release. The accession test uses card C-74's `lookup` over a tiny hand-built index (never a snapshot of this script's
output). Every test runs twice, with the normal PATH and with PATH reduced to `/usr/bin:/bin` (as
tests/r1/test_r232_index.py does): nothing here needs a tool beyond Python and pyarrow, and GitHub Actions has none of
the pod's tools. No test reads a top-level file of the repository (tests/r0/test_image_files.py).
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-76")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-26_lr2_sag1_nx1_r232.md"


def _load(name: str, file: str):
    sys.path.insert(0, str(BSP))
    spec = importlib.util.spec_from_file_location(name, BSP / file)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(BSP / file)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _mod():
    return _load("c76_lr2_sag1_nx1_r232", "c76_lr2_sag1_nx1_r232.py")


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def seq_for(tag: str) -> str:
    return (tag * 40)[:80] + "ACGT" * 5


def uid_of(tag: str) -> str:
    return blake3(seq_for(tag).upper().encode()).digest(length=16).hex()


# --- a tiny r232 index (C-74's file layout, only the columns `lookup` reads) ----------------------------------------

@pytest.fixture()
def index(tmp_path: Path) -> Path:
    d = tmp_path / "index"
    d.mkdir()
    gt = pa.table({
        "accession": ["GCA_000000011", "GCF_000000077"],
        "genbank_accession": ["GCA_000000011", "GCA_000000078"],
        "is_representative": [True, False],
        "representative": ["GCA_000000011", "GCA_000000011"],
        "species": ["Archaeum alpha", "Archaeum alpha"],
    })
    pq.write_table(gt, d / "gtdb_r232.parquet")
    gl = pa.table({
        "globdb_id": ["GCA_000000011", "MOTU40_000001", "GEM_00007"],
        "insdc_accession": ["GCA_000000011", None, "GCF_000000033"],
    })
    pq.write_table(gl, d / "globdb_r232.parquet")
    return d


# --- a fixture stage in the layout of cards C-66 and C-58 -----------------------------------------------------------

# genome, source, contigs; G1 is in GTDB R232 (and GlobDB), G4 only in GlobDB r232 (its INSDC accession), G2 in neither,
# G3 carries no accession. G1 and G3 share contig b (DEC-27).
GENOMES = [
    ("ehomd:GCA_000000011.1", "ehomd_v11_03_genomes", ["a", "b"]),
    ("lichen:GCA_000000099.1", "ena_c66_nx1_holobiont_subsurface", ["c"]),
    ("lost_city:BinX", "lost_city_serpentinizing_mags", ["b", "d"]),
    ("samail:GCF_000000033.2", "ena_c66_nx1_holobiont_subsurface", ["e", "f"]),
]
RECORDS_HEAD = ["genome_id", "source_name", "record_ordinal", "contig_uid", "length", "circular", "kind",
                "assembly_quality", "genetic_code", "taxonomy"]


def make_stage(d: Path, *, genomes=GENOMES, compact: bool = True) -> Path:
    d.mkdir(parents=True)
    rows, gj = [], []
    for gid, src, tags in genomes:
        for i, t in enumerate(tags):
            rows.append("\t".join([gid, src, str(i), uid_of(t), str(len(seq_for(t))), "0", "prokaryote", "mag_mq",
                                   "11", "t"]))
        gj.append({"contig_uids": [uid_of(t) for t in tags], "genetic_code": 11, "kind": "prokaryote",
                   "circular": [], "assembly_quality": "mag_mq"})
    (d / "records.tsv").write_text("\t".join(RECORDS_HEAD) + "\n" + "\n".join(rows) + "\n")
    seps = (",", ":") if compact else (", ", ": ")
    (d / "genomes.json").write_text(json.dumps(gj, separators=seps))
    (d / "shared_records.tsv").write_text(
        "contig_uid\towner_stage\tlength\n"
        f"{uid_of('a')}\tPL1\t{len(seq_for('a'))}\n{uid_of('d')}\tPL3\t{len(seq_for('d'))}\n")
    (d / "repeated_records.tsv").write_text("genome_id\tsource_name\treason\tfirst_contig_uid\n")
    (d / "trims.tsv").write_text("genome_id\tcontig_name\traw_len\ttrimmed_k\tfinal_len\n"
                                 "ehomd:GCA_000000011.1\ta\t100\t20\t80\nlost_city:BinX\td\t100\t20\t80\n")
    (d / "genetic_codes.tsv").write_text(
        "genome_id\tsource_name\tgenome_bp\tgenetic_code\n" +
        "".join(f"{g}\t{s}\t1\t11\n" for g, s, _ in genomes))
    obj = [{"path": "/tmp/if2453/s2/inputs/NX1/nx1_compendia.fna.gz", "source_name": "c66_nx1_compendia",
            "source_release": "2026-11-24", "license": "CC BY 4.0", "redistributable": True,
            "retrieved_at": "2026-11-24T00:00:00Z"}]
    (d / "objects.json").write_text(json.dumps(obj, indent=2))
    obj[0]["path"] = "/stor/scratch/Wilke/if2453/bio-box-staging/NX1/nx1_compendia.fna.gz"
    (d / "objects_staging.json").write_text(json.dumps(obj, indent=2))
    return d


def make_fasta(path: Path, tags: list[str], *, assembler_names: bool = False) -> Path:
    """A gzip FASTA, records wrapped at 80; the header is the contig_uid, or (SAG1's case) the assembler's name."""
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(path, "wb"), mtime=0) as g:
        for t in tags:
            s = seq_for(t)
            h = f"NODE_{t};1655516;5094" if assembler_names else uid_of(t)
            g.write(f">{h}\n{s[:80]}\n{s[80:]}\n".encode())
    return path


@pytest.fixture()
def stage(tmp_path: Path) -> dict:
    meta = make_stage(tmp_path / "meta")
    fa = make_fasta(tmp_path / "nx1_compendia.fna.gz", ["a", "b", "c", "d", "e", "f"])
    return {"meta": meta, "fasta": fa, "tmp": tmp_path}


def fasta_uids(path: Path) -> list[str]:
    with gzip.open(path, "rt") as f:
        return [ln[1:].strip() for ln in f if ln.startswith(">")]


# --- accessions -----------------------------------------------------------------------------------------------------

def test_accession_of_reads_every_form_and_none_for_names():
    m = _mod()
    assert m.accession_of("ehomd:GCA_000000011.1") == "GCA_000000011"
    assert m.accession_of("GCF_000000033.2") == "GCF_000000033"
    assert m.accession_of("samail:GCA_039666335.1") == "GCA_039666335"
    assert m.accession_of("lost_city:Chloroflexi_metaspades_Bin_1_4-contigs") is None
    assert m.accession_of("emerge2019:cmr5.MA.201907_P_1_20to24_24") is None
    assert m.accession_of("AG-313-A04") is None


def test_names_table_maps_unique_names_only(tmp_path: Path):
    m = _mod()
    t = tmp_path / "ena.tsv"
    t.write_text("accession\tassembly_name\tbase_count\n"
                 "GCA_000000011\tAG-1\t5\nGCA_000000022\tAG-2\t5\nGCA_000000033\tAG-2\t5\n")
    assert m.read_names(t) == {"AG-1": "GCA_000000011"}


# --- the r232 status (a genome the index finds is kept and carries its status) ---------------------------------------

def stage_hashes(stage) -> dict:
    return {p.name: sha(p) for p in stage["meta"].iterdir()} | {stage["fasta"].name: sha(stage["fasta"])}


def test_status_of_each_genome_from_the_index(stage, index):
    m = _mod()
    r = m.status(m.index_lookup(index), stage["meta"] / "records.tsv")
    by = {x["genome_id"]: x for x in r["rows"]}
    assert [x["genome_id"] for x in r["rows"]] == sorted(by)  # sorted
    assert by["ehomd:GCA_000000011.1"]["status"] == "in_gtdb_r232+in_globdb_r232"
    assert by["ehomd:GCA_000000011.1"]["is_representative"] is True
    assert by["ehomd:GCA_000000011.1"]["species"] == "Archaeum alpha"
    assert by["samail:GCF_000000033.2"]["status"] == "in_globdb_r232"
    assert by["lichen:GCA_000000099.1"]["status"] == "not_in_r232"
    assert by["lost_city:BinX"]["status"] == "not_visible" and by["lost_city:BinX"]["accession"] is None
    assert r["by_status"] == {"in_gtdb_r232+in_globdb_r232": 1, "in_gtdb_r232": 0, "in_globdb_r232": 1,
                              "not_in_r232": 1, "not_visible": 1}
    assert (r["genomes"], r["with_accession"], r["without_accession"]) == (4, 3, 1)
    assert (r["in_gtdb_r232"], r["in_globdb_r232"]) == (1, 2)
    assert r["by_source"]["ena_c66_nx1_holobiont_subsurface"] == {"in_globdb_r232": 1, "not_in_r232": 1}
    assert r["by_source"]["lost_city_serpentinizing_mags"] == {"not_visible": 1}


def test_a_sag_is_seen_through_the_ena_names_table(tmp_path: Path, index):
    m = _mod()
    d = tmp_path / "sag"
    d.mkdir()
    (d / "records.tsv").write_text(
        "ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\n"
        "0\tAG-313-A04\t0\tu1\t10\tsag\t11\n1\tAG-313-A04\t1\tu2\t10\tsag\t11\n2\tAG-313-A07\t0\tu3\t10\tsag\t11\n"
        "3\tAG-313-Z99\t0\tu4\t10\tsag\t11\n")
    names = {"AG-313-A04": "GCA_000000011", "AG-313-A07": "GCA_000000555"}
    r = m.status(m.index_lookup(index), d / "records.tsv", names)
    by = {x["genome_id"]: x for x in r["rows"]}
    assert by["AG-313-A04"]["status"] == "in_gtdb_r232+in_globdb_r232" and by["AG-313-A04"]["how"] == "ena_name"
    assert by["AG-313-A07"]["status"] == "not_in_r232"
    assert by["AG-313-Z99"]["status"] == "not_visible"


def test_lookup_is_the_index_lookup_of_card_c74(index):
    m = _mod()
    got = m.index_lookup(index)(["GCA_000000011", "GCA_000000999", "GCF_000000033"])
    assert [(g["in_gtdb_r232"], g["in_globdb_r232"]) for g in got] == [(True, True), (False, False), (False, True)]


def test_run_removes_nothing_and_leaves_the_stages_files_byte_identical(stage, index):
    m = _mod()
    before = stage_hashes(stage)
    rep = m.run("NX1", stage["meta"], stage["fasta"], index, stage["tmp"] / "out", check=False)
    assert before == stage_hashes(stage)  # the stage's metadata and FASTA are untouched
    assert fasta_uids(stage["fasta"]) == [uid_of(t) for t in "abcdef"]
    assert rep["removed"] == 0
    assert rep["stage_totals"] == {"genomes": 4, "contigs": 6, "bases": 6 * len(seq_for("a"))}
    assert rep["by_status"]["in_gtdb_r232+in_globdb_r232"] == 1 and rep["by_status"]["in_globdb_r232"] == 1
    # nothing but the two reports is written
    assert sorted(p.name for p in (stage["tmp"] / "out" / "NX1").iterdir()) == ["r232_report.json", "r232_status.tsv"]
    # a genome the index finds is still in the stage and carries its status in the table
    lines = (stage["tmp"] / "out" / "NX1" / "r232_status.tsv").read_text().splitlines()
    assert lines[0].split("\t")[:5] == ["genome_id", "source_name", "accession", "how", "status"]
    assert [ln.split("\t")[0] for ln in lines[1:]] == sorted(g for g, _, _ in GENOMES)
    row = dict(zip(lines[0].split("\t"), next(ln for ln in lines[1:] if ln.startswith("ehomd:")).split("\t")))
    assert (row["status"], row["in_gtdb_r232"], row["is_representative"], row["in_globdb_r232"],
            row["globdb_r232_id"]) == ("in_gtdb_r232+in_globdb_r232", "1", "1", "1", "GCA_000000011")
    assert "ehomd:GCA_000000011.1" in (stage["meta"] / "records.tsv").read_text()
    assert len(json.loads((stage["meta"] / "genomes.json").read_text())) == 4


def test_the_status_table_is_byte_identical_when_written_again(stage, index):
    m = _mod()
    m.run("NX1", stage["meta"], stage["fasta"], index, stage["tmp"] / "a", check=False)
    m.run("NX1", stage["meta"], stage["fasta"], index, stage["tmp"] / "b", check=False)
    for n in ("r232_status.tsv", "r232_report.json"):
        assert sha(stage["tmp"] / "a" / "NX1" / n) == sha(stage["tmp"] / "b" / "NX1" / n), n


def test_refusals_existing_target_and_missing_input(stage, index):
    m = _mod()
    m.run("NX1", stage["meta"], stage["fasta"], index, stage["tmp"] / "o", check=False)
    with pytest.raises(m.Refused, match="exists"):
        m.run("NX1", stage["meta"], stage["fasta"], index, stage["tmp"] / "o", check=False)
    with pytest.raises(m.Refused, match="missing"):
        m.run("NX1", stage["meta"], stage["tmp"] / "nofile.gz", index, stage["tmp"] / "o2", check=False)


def test_the_command_line_runs_a_stage_and_exits_zero(stage, index, capsys):
    m = _mod()
    rc = m.main(["run", "NX1", str(stage["meta"]), str(stage["fasta"]), str(index), str(stage["tmp"] / "o"),
                 "--no-check"])
    assert rc == 0
    rep = json.loads(capsys.readouterr().out)
    assert rep["removed"] == 0 and rep["genomes"] == 4 and rep["in_gtdb_r232"] == 1


# --- the conf templates ---------------------------------------------------------------------------------------------

def test_templates_pass_stage_est_and_carry_the_cards_values():
    est = _load("stage_est", "stage_est.py")
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    confs = {s: est.read_conf(TEMPLATES / f"{s}.conf", s, tree=REPO) for s in ("LR2", "SAG1", "NX1")}
    for s, c in confs.items():
        assert (c.depth, c.crispr, c.library) == ("called", True, None), s
        assert str(c.inputs) == f"/tmp/if2453/s2/inputs/{s}"
        assert c.min_contig_len is None, s  # DEC-60 (1): MAGs and SAGs whole; LR2's unbinned contigs were cut at 5 kb
        assert NOTE in c.notes
    assert confs["LR2"].after == ("GRARE3",)
    assert confs["SAG1"].after == ("LR2",)
    assert confs["NX1"].after == ("SAG1",)
    assert confs["LR2"].genecall_mode == "gv-meta-c25"
    assert confs["SAG1"].genecall_mode == "default" and confs["NX1"].genecall_mode == "default"
    fasta = {"LR2": "lr2_compendia.fna.gz", "SAG1": "sag1_gorg.fna.gz", "NX1": "nx1_compendia.fna.gz"}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    allowed = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}
    allowed |= {s for s, _ in a26.R26_11_4_CHAIN}
    for s, c in confs.items():
        assert {n for _, n in c.sums} == {"genomes.json", "objects.json", fasta[s]}, s
        for h, _ in c.sums:
            assert h in note, f"{s}: SUM {h} not in the note"
        assert set(c.expect) <= allowed
        for owner, (n, _) in c.expect.items():
            assert n > 0

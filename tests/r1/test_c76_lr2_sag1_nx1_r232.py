"""Card C-76 (DEC-79): the R232 / r232 re-diff of the staged LR2, SAG1 and NX1, and their conf templates.

Truth is written by hand in this file: which fixture genome sits in which release, which records and contigs must
survive. The accession test uses card C-74's `lookup` over a tiny hand-built index (never a snapshot of this script's
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


# --- the re-diff ----------------------------------------------------------------------------------------------------

def test_rediff_removes_a_genome_found_in_gtdb_or_globdb_and_lists_the_reason(stage, index):
    m = _mod()
    r = m.rediff(stage["meta"] / "records.tsv", m.index_lookup(index))
    removed = {x["genome_id"]: x for x in r["removed"]}
    assert set(removed) == {"ehomd:GCA_000000011.1", "samail:GCF_000000033.2"}
    assert removed["ehomd:GCA_000000011.1"]["reason"] == "in_gtdb_r232+in_globdb_r232"
    assert removed["samail:GCF_000000033.2"]["reason"] == "in_globdb_r232"
    assert removed["samail:GCF_000000033.2"]["accession"] == "GCF_000000033"
    # what the test could and could not see: three of four carry an accession, one does not
    assert (r["genomes"], r["with_accession"], r["without_accession"]) == (4, 3, 1)
    assert r["by_source"]["lost_city_serpentinizing_mags"] == {"genomes": 1, "with_accession": 0, "removed": 0}
    assert r["by_source"]["ena_c66_nx1_holobiont_subsurface"] == {"genomes": 2, "with_accession": 2, "removed": 1}


def test_rediff_sees_a_sag_through_the_ena_names_table(tmp_path: Path, index):
    m = _mod()
    d = tmp_path / "sag"
    d.mkdir()
    (d / "records.tsv").write_text(
        "ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\n"
        "0\tAG-313-A04\t0\tu1\t10\tsag\t11\n1\tAG-313-A04\t1\tu2\t10\tsag\t11\n2\tAG-313-A07\t0\tu3\t10\tsag\t11\n")
    names = {"AG-313-A04": "GCA_000000011", "AG-313-A07": "GCA_000000555"}
    r = m.rediff(d / "records.tsv", m.index_lookup(index), names)
    assert [x["genome_id"] for x in r["removed"]] == ["AG-313-A04"]
    assert r["removed"][0]["how"] == "ena_name"
    assert (r["genomes"], r["with_accession"], r["without_accession"]) == (2, 2, 0)


def test_lookup_is_the_index_lookup_of_card_c74(index):
    m = _mod()
    got = m.index_lookup(index)(["GCA_000000011", "GCA_000000999", "GCF_000000033"])
    assert [(g["in_gtdb_r232"], g["in_globdb_r232"]) for g in got] == [(True, True), (False, False), (False, True)]


# --- the rewrite ----------------------------------------------------------------------------------------------------

def run_changed(stage, index, tag="a"):
    m = _mod()
    out = stage["tmp"] / f"out_{tag}"
    fa_out = stage["tmp"] / f"fa_{tag}"
    rep = m.run("NX1", stage["meta"], stage["fasta"], index, out, fa_out, inputs_home=Path("/tmp/if2453/s2/inputs"),
                check=False)
    return rep, out / "NX1", fa_out


def test_a_changed_stage_is_rewritten_with_the_right_genomes_records_and_contigs(stage, index):
    rep, out, fa_out = run_changed(stage, index)
    assert rep["changed"] is True
    new_fa = fa_out / "nx1_compendia_r232.fna.gz"
    # genomes.json: G2 and G3 survive, byte-for-byte the style of the input (compact)
    old = json.loads((stage["meta"] / "genomes.json").read_text())
    assert (out / "genomes.json").read_text() == json.dumps([old[1], old[2]], separators=(",", ":"))
    # records: only G2 and G3, in order
    lines = (out / "records.tsv").read_text().splitlines()
    assert lines[0].split("\t") == RECORDS_HEAD
    assert [ln.split("\t")[0] for ln in lines[1:]] == ["lichen:GCA_000000099.1", "lost_city:BinX", "lost_city:BinX"]
    # FASTA: c, b (shared with the removed G1, still in G3), d; a, e, f gone; the old order kept
    assert fasta_uids(new_fa) == [uid_of("b"), uid_of("c"), uid_of("d")]
    # shared records: a (held by PL1) belonged to G1 only and goes; d stays
    assert (out / "shared_records.tsv").read_text().splitlines()[1:] == [
        f"{uid_of('d')}\tPL3\t{len(seq_for('d'))}"]
    assert (out / "trims.tsv").read_text().splitlines()[1:] == ["lost_city:BinX\td\t100\t20\t80"]
    assert [ln.split("\t")[0] for ln in (out / "genetic_codes.tsv").read_text().splitlines()[1:]] == [
        "lichen:GCA_000000099.1", "lost_city:BinX"]
    # objects: the new names, the old fields
    o = json.loads((out / "objects.json").read_text())[0]
    s = json.loads((out / "objects_staging.json").read_text())[0]
    assert o["path"] == "/tmp/if2453/s2/inputs/NX1/nx1_compendia_r232.fna.gz"
    assert s["path"] == str(new_fa)
    assert o["license"] == "CC BY 4.0" and o["redistributable"] is True
    # the old files are untouched and the removed genomes are listed
    assert fasta_uids(stage["fasta"]) == [uid_of(t) for t in "abcdef"]
    rem = (out / "r232_removed.tsv").read_text().splitlines()
    assert rem[0].split("\t")[:3] == ["genome_id", "accession", "reason"]
    assert {ln.split("\t")[0] for ln in rem[1:]} == {"ehomd:GCA_000000011.1", "samail:GCF_000000033.2"}
    assert rep["after"] == {"genomes": 2, "contigs": 3, "bases": 3 * len(seq_for("a"))}
    assert rep["before"] == {"genomes": 4, "contigs": 6, "bases": 6 * len(seq_for("a"))}


def test_a_changed_stage_is_byte_identical_when_rewritten_again(stage, index):
    run_changed(stage, index, "a")
    run_changed(stage, index, "b")
    t = stage["tmp"]
    names = sorted(p.name for p in (t / "out_a" / "NX1").iterdir())
    assert names == sorted(p.name for p in (t / "out_b" / "NX1").iterdir())
    for n in names:
        if n == "objects_staging.json":  # it names the FASTA's own directory, which differs between the two runs
            continue
        if n == "r232_rediff.json":  # its file table holds that file's checksum
            ra, rb = (json.loads((t / d / "NX1" / n).read_text()) for d in ("out_a", "out_b"))
            ra["files"].pop("objects_staging.json")
            rb["files"].pop("objects_staging.json")
            assert ra == rb
            continue
        assert sha(t / "out_a" / "NX1" / n) == sha(t / "out_b" / "NX1" / n), n
    a = json.loads((t / "out_a" / "NX1" / "objects_staging.json").read_text())
    b = json.loads((t / "out_b" / "NX1" / "objects_staging.json").read_text())
    assert a[0]["path"].endswith("fa_a/nx1_compendia_r232.fna.gz")
    assert {**a[0], "path": ""} == {**b[0], "path": ""}
    assert sha(t / "fa_a" / "nx1_compendia_r232.fna.gz") == sha(t / "fa_b" / "nx1_compendia_r232.fna.gz")


def test_an_unchanged_stage_keeps_its_files_and_writes_only_the_report(tmp_path: Path, index):
    m = _mod()
    genomes = [("lichen:GCA_000000099.1", "s1", ["c"]), ("lost_city:BinX", "s2", ["d"])]
    meta = make_stage(tmp_path / "meta", genomes=genomes)
    fa = make_fasta(tmp_path / "nx1_compendia.fna.gz", ["c", "d"])
    before = {p.name: sha(p) for p in meta.iterdir()} | {fa.name: sha(fa)}
    rep = m.run("NX1", meta, fa, index, tmp_path / "out", tmp_path / "fa", inputs_home=Path("/x"), check=False)
    assert rep["changed"] is False and rep["removed"] == 0
    assert before == {p.name: sha(p) for p in meta.iterdir()} | {fa.name: sha(fa)}
    assert sorted(p.name for p in (tmp_path / "out" / "NX1").iterdir()) == [
        "r232_rediff.json", "r232_rediff.tsv", "r232_removed.tsv"]
    assert not (tmp_path / "fa").exists() or not any((tmp_path / "fa").iterdir())


def test_ordinal_columns_are_renumbered_after_a_removal(tmp_path: Path, index):
    m = _mod()
    d = tmp_path / "sag"
    d.mkdir()
    tags = {"AG-1": ["a", "b"], "AG-2": ["c"], "AG-3": ["d"]}
    rows = []
    n = 0
    for g, ts in tags.items():
        for i, t in enumerate(ts):
            rows.append(f"{n}\t{g}\t{i}\t{uid_of(t)}\t{len(seq_for(t))}\tsag\t11")
            n += 1
    (d / "records.tsv").write_text(
        "ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\n" + "\n".join(rows) + "\n")
    (d / "genomes.json").write_text(json.dumps([
        {"contig_uids": [uid_of(t) for t in ts], "genetic_code": 11, "kind": "prokaryote", "circular": [],
         "assembly_quality": "contig"} for ts in tags.values()]))
    (d / "shared_records.tsv").write_text(
        "ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\towner_stage\n"
        f"1\tAG-1\t1\t{uid_of('b')}\t{len(seq_for('b'))}\tMV1\n3\tAG-3\t0\t{uid_of('d')}\t{len(seq_for('d'))}\tMV1\n")
    (d / "repeated_records.tsv").write_text("ordinal\tgenome\tgenome_uid\tfirst_genome\n")
    (d / "genetic_codes.tsv").write_text("genome\tstage\tgenetic_code\nAG-1\tSAG1\t11\nAG-2\tSAG1\t11\nAG-3\tSAG1\t11\n")
    (d / "objects.json").write_text(json.dumps([{"path": "/tmp/x/sag1_gorg.fna.gz", "license": "l"}]))
    (d / "objects_staging.json").write_text(json.dumps([{"path": "/stor/x/sag1_gorg.fna.gz", "license": "l"}]))
    fa = make_fasta(tmp_path / "sag1_gorg.fna.gz", ["a", "b", "c", "d"], assembler_names=True)
    names = tmp_path / "ena.tsv"
    names.write_text("accession\tassembly_name\nGCA_000000011\tAG-1\nGCA_000000555\tAG-2\n")
    rep = m.run("SAG1", d, fa, index, tmp_path / "out", tmp_path / "fa", inputs_home=Path("/tmp/if2453/s2/inputs"),
                names=m.read_names(names), check=False)
    out = tmp_path / "out" / "SAG1"
    assert rep["changed"] is True and rep["removed"] == 1
    recs = [ln.split("\t") for ln in (out / "records.tsv").read_text().splitlines()[1:]]
    assert [(r[0], r[1]) for r in recs] == [("0", "AG-2"), ("1", "AG-3")]
    shared = [ln.split("\t") for ln in (out / "shared_records.tsv").read_text().splitlines()[1:]]
    assert [(s[0], s[1]) for s in shared] == [("1", "AG-3")]  # AG-3's record was ordinal 3 and is now 1
    # the records are found by the contig_uid of their sequence; the headers (the assembler's names) are kept as they are
    assert fasta_uids(tmp_path / "fa" / "sag1_gorg_r232.fna.gz") == ["NODE_c;1655516;5094", "NODE_d;1655516;5094"]
    assert len(json.loads((out / "genomes.json").read_text())) == 2


def test_the_rewrite_refuses_inputs_that_do_not_agree(stage, index):
    m = _mod()
    # a genomes.json that lists another number of genomes than records.tsv
    g = json.loads((stage["meta"] / "genomes.json").read_text())
    (stage["meta"] / "genomes.json").write_text(json.dumps(g[:3]))
    with pytest.raises(m.Refused, match="genomes"):
        m.run("NX1", stage["meta"], stage["fasta"], index, stage["tmp"] / "o", stage["tmp"] / "f",
              inputs_home=Path("/x"), check=False)


def test_the_rewrite_refuses_a_fasta_that_lacks_a_kept_contig(stage, index):
    m = _mod()
    short = make_fasta(stage["tmp"] / "short.fna.gz", ["a", "b", "c", "e", "f"])  # no d
    with pytest.raises(m.Refused, match="missing"):
        m.run("NX1", stage["meta"], short, index, stage["tmp"] / "o", stage["tmp"] / "f", inputs_home=Path("/x"),
              check=False)
    assert not (stage["tmp"] / "f").exists() or not list((stage["tmp"] / "f").glob("*.fna.gz"))


def test_nothing_is_overwritten(stage, index):
    m = _mod()
    run_changed(stage, index, "a")
    with pytest.raises(m.Refused, match="exists"):
        m.run("NX1", stage["meta"], stage["fasta"], index, stage["tmp"] / "out_a", stage["tmp"] / "fa_a",
              inputs_home=Path("/tmp/if2453/s2/inputs"), check=False)


def test_the_command_line_runs_a_stage_and_exits_zero(stage, index, capsys):
    m = _mod()
    rc = m.main(["run", "NX1", str(stage["meta"]), str(stage["fasta"]), str(index), str(stage["tmp"] / "o"),
                 str(stage["tmp"] / "f"), "--no-check", "--inputs-home", "/tmp/if2453/s2/inputs"])
    assert rc == 0
    rep = json.loads(capsys.readouterr().out)
    assert rep["removed"] == 2 and rep["after"]["genomes"] == 2


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
    fasta = {"LR2": "lr2_compendia.fna.gz", "SAG1": "sag1_gorg_r232.fna.gz", "NX1": "nx1_compendia_r232.fna.gz"}
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

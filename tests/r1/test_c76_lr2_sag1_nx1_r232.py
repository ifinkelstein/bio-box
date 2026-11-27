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
R19_NOTE = "planning/measurements/2026-11-27_staged_fasta_scan.md"


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
    # card R-19: LR2 is re-staged under a new FASTA name; its SUMs are recorded in R-19's note, the others in C-76's
    # card R-20: NX1 is re-staged without its three 16S alignment files; its SUMs are in R-19's note (a dated section)
    fasta = {"LR2": "lr2_compendia_r19.fna.gz", "SAG1": "sag1_gorg.fna.gz", "NX1": "nx1_compendia_r20.fna.gz"}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    r19 = (REPO / R19_NOTE).read_text(encoding="utf-8")
    notes = {"LR2": r19, "SAG1": note, "NX1": r19}
    assert R19_NOTE in confs["LR2"].notes
    assert R19_NOTE in confs["NX1"].notes
    allowed = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}
    allowed |= {s for s, _ in a26.R26_11_4_CHAIN}
    for s, c in confs.items():
        assert {n for _, n in c.sums} == {"genomes.json", "objects.json", fasta[s]}, s
        for h, _ in c.sums:
            assert h in notes[s], f"{s}: SUM {h} not in the note"
        assert set(c.expect) <= allowed
        for owner, (n, _) in c.expect.items():
            assert n > 0


# --- card R-19: LR2 re-staged from its four archives ----------------------------------------------------------------
#
# The fixture archives have the layout of the four LR2 inputs on the pod. Truth by hand: which contig is in which
# member, which member is gzip-compressed (the Bickhart bins' `bin.1.fa.gz`, read as plain text by card C-66), which
# assembly contig is binned, short or circular, and the genetic codes C-66 recorded for the two sources it read right.

pytestmark_r19 = pytest.mark.covers("R-19")


def dna(tag: str, n: int) -> bytes:
    out, h = bytearray(), hashlib.sha256(tag.encode()).digest()
    while len(out) < n:
        out += bytes(b"ACGT"[x & 3] for x in h)
        h = hashlib.sha256(h).digest()
    return bytes(out[:n])


def wrap(seq: bytes, w: int = 60) -> bytes:
    return b"".join(seq[i:i + w] + b"\n" for i in range(0, len(seq), w))


def fasta_bytes(recs) -> bytes:
    return b"".join(b">" + h.encode() + b"\n" + wrap(s) for h, s in recs)


def _tar(path: Path, members: dict[str, bytes]) -> None:
    import io
    import tarfile

    with tarfile.open(path, "w:gz") as t:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            t.addfile(info, io.BytesIO(data))


SEQ = {"e1a": dna("e1a", 3000), "e1b": dna("e1b", 2000), "s1": dna("s1", 2500), "s2": dna("s2", 2600),
       "c1": dna("c1", 6000), "c2": dna("c2", 3000), "c3": dna("c3", 5500) + dna("c3", 50), "c4": dna("c4", 7000),
       "c5": dna("c5", 4000)}  # contig_3 ends with a copy of its first 50 bases (a terminal direct repeat)
CIRCULAR = {3, 4}  # assembly_info.txt: contig_3 (binned, its repeat trimmed only in the unbinned pass) and contig_4


def uid(seq: bytes) -> str:
    return blake3(seq.upper()).digest(length=16).hex()


def lr2_world(tmp: Path, *, bad_bin: bool = False) -> dict:
    dl = tmp / "dl"
    dl.mkdir(parents=True)
    _tar(dl / "EMERGE_MAGs_2019_long-short-hybrid.tar.gz",
         {"MAGs/e1.fa": fasta_bytes([("k141_1 len=3000", SEQ["e1a"]), ("k141_2", SEQ["e1b"])]),
          "MAGs/README.txt": b"not a genome\n"})
    _tar(dl / "MQ_HQ_MAGs_20200502.tar.gz",
         {"MAGs/s1.fa.gz": gzip.compress(fasta_bytes([("s1_c1", SEQ["s1"])]), mtime=0),
          "MAGs/s2.fa": fasta_bytes([("s2_c1", SEQ["s2"])])})
    c3 = SEQ["c3"][:100] + b"\xe5" + SEQ["c3"][101:] if bad_bin else SEQ["c3"]
    _tar(dl / "hifi_das.bin3c.bins.tar.gz",
         {"bins/bin.1.fa.gz": gzip.compress(fasta_bytes([("contig_1", SEQ["c1"]), ("contig_2", SEQ["c2"])]), mtime=0),
          "bins/bin.2.fa": fasta_bytes([("contig_3", c3)])})
    info = "#seq_name\tlength\tcov.\tcirc.\trepeat\n" + "".join(
        f"contig_{i}\t{len(SEQ[f'c{i}'])}\t10\t{'Y' if i in CIRCULAR else 'N'}\tN\n" for i in range(1, 6))
    _tar(dl / "metaflye_hifi_asm.tar.gz",
         {"asm/assembly_info.txt": info.encode(),
          "asm/assembly.fasta": fasta_bytes([(f"contig_{i}", SEQ[f"c{i}"]) for i in range(1, 6)])})
    meta = tmp / "metadata_MAGs_2019_EMERGE.tsv"
    meta.write_text("filename\tclassification\tCompleteness\tContamination\ne1.fa\td__Bacteria;s__x\t95\t1\n")
    old = tmp / "old"
    old.mkdir()
    gc_rows = [["emerge2019:e1", "emerge_stordalen_2019_hifi", "5000", "2", "mag_hq", "11", "below_100kb", "False",
                "0.0", "0.0", "", ""],
               ["singleton2021:s1", "singleton_2021_wwtp_ont", "2500", "1", "mag_hq", "25", "code25_gain", "True",
                "0.7", "300.0", "0.9", "0.7"],
               ["singleton2021:s2", "singleton_2021_wwtp_ont", "2600", "1", "mag_hq", "11", "below_100kb", "False",
                "0.0", "0.0", "", ""]]
    head = ["genome_id", "source_name", "genome_bp", "n_contigs", "assembly_quality", "genetic_code", "reason",
            "tested_alt", "cd_11", "mean_orf_11", "cd_25", "cd_4"]
    (old / "genetic_codes.tsv").write_bytes(("\r\n".join("\t".join(r) for r in [head, *gc_rows]) + "\r\n").encode())
    # C-66's records of the two sources it read right, as C-66 wrote them (csv, CRLF)
    rec = [["emerge2019:e1", "emerge_stordalen_2019_hifi", "0", uid(SEQ["e1a"]), "3000", "0", "prokaryote", "mag_hq",
            "11", "d__Bacteria;s__x"],
           ["emerge2019:e1", "emerge_stordalen_2019_hifi", "1", uid(SEQ["e1b"]), "2000", "0", "prokaryote", "mag_hq",
            "11", "d__Bacteria;s__x"],
           ["singleton2021:s1", "singleton_2021_wwtp_ont", "0", uid(SEQ["s1"]), "2500", "0", "prokaryote", "mag_hq",
            "25", ""],
           ["singleton2021:s2", "singleton_2021_wwtp_ont", "0", uid(SEQ["s2"]), "2600", "0", "prokaryote", "mag_hq",
            "11", ""],
           ["bickhart_bin:bin.1", "bickhart_2022_rumen_hifi", "0", "00" * 16, "9", "0", "prokaryote", "mag_hq", "11",
            ""]]
    (old / "records.tsv").write_bytes(("\r\n".join("\t".join(r) for r in [RECORDS_HEAD, *rec]) + "\r\n").encode())
    (old / "objects_staging.json").write_text(json.dumps([{
        "path": "/stor/x/LR2/lr2_compendia.fna.gz", "source_name": "c66_lr2_compendia", "source_release": "2026-11-24",
        "license": "CC BY 4.0", "redistributable": True, "retrieved_at": "2026-11-24T00:00:00Z"}], indent=2))
    return {"dl": dl, "meta": meta, "old": old, "tmp": tmp}


def _restage(m, w, name="out"):
    return m.restage_lr2(w["dl"], w["meta"], w["old"], w["tmp"] / name, w["tmp"] / name / "fa" / "lr2_r19.fna.gz",
                         "/tmp/if2453/s2/inputs/LR2", workers=1, compressor="gzip")


def _read_fasta(path: Path) -> list[tuple[str, bytes]]:
    recs, h, s = [], None, []
    with gzip.open(path, "rb") as f:
        for ln in f:
            ln = ln.rstrip(b"\n")
            if ln.startswith(b">"):
                if h is not None:
                    recs.append((h, b"".join(s)))
                h, s = ln[1:].decode(), []
            else:
                s.append(ln)
    if h is not None:
        recs.append((h, b"".join(s)))
    return recs


@pytestmark_r19
def test_r19_a_gzip_member_is_read_by_its_magic_whatever_its_name(tmp_path: Path):
    import io
    import tarfile

    m = _mod()
    p = tmp_path / "a.tar.gz"
    _tar(p, {"x.fa": gzip.compress(b">a\nACGT\n", mtime=0), "y.fa.gz": b">b\nGG\n"})
    with tarfile.open(p, "r:*") as tf:
        got = {mb.name: list(m.strict_records(m.member_stream(tf, mb), mb.name)) for mb in tf if mb.isfile()}
    assert got == {"x.fa": [("a", b"ACGT")], "y.fa.gz": [("b", b"GG")]}
    with pytest.raises(m.Refused, match="not UTF-8"):
        list(m.strict_records(io.BytesIO(b">\xff\nAC\n"), "t"))
    with pytest.raises(m.Refused, match="IUPAC"):
        list(m.strict_records(io.BytesIO(b">a\nAC\x1f\x8b\n"), "t"))
    with pytest.raises(m.Refused, match="before the first header"):
        list(m.strict_records(io.BytesIO(b"\x1f\x8b\x08\n>a\nAC\n"), "t"))


@pytestmark_r19
def test_r19_restage_decompresses_the_bins_and_writes_only_sequence(tmp_path: Path):
    m = _mod()
    w = lr2_world(tmp_path)
    rep = _restage(m, w)
    out = tmp_path / "out"
    recs = _read_fasta(Path(rep["fasta_path"]))
    want = [SEQ[k] for k in ("e1a", "e1b", "s1", "s2", "c1", "c2", "c3", "c4")]
    assert [s for _, s in recs] == want  # the bins' sequences, not their compressed bytes
    assert [h for h, _ in recs] == [uid(s) for s in want]
    scan = _load("r19_fasta_scan", "r19_fasta_scan.py").scan_file(Path(rep["fasta_path"]))
    assert scan["status"] == "clean" and scan["records"] == 8
    assert rep["members"] == {"bins:gzip": 1, "bins:plain": 1, "emerge:plain": 1, "singleton:gzip": 1,
                              "singleton:plain": 1}
    g = json.loads((out / "genomes.json").read_text())
    assert [x["contig_uids"] for x in g] == [[uid(SEQ["e1a"]), uid(SEQ["e1b"])], [uid(SEQ["s1"])], [uid(SEQ["s2"])],
                                             [uid(SEQ["c1"]), uid(SEQ["c2"])], [uid(SEQ["c3"])], [uid(SEQ["c4"])]]
    # contig_1 and contig_3 are binned and left out of the unbinned set (contig_3 by its name: its repeat is trimmed
    # in the unbinned pass, so its uid is not its bin's); contig_2 and contig_5 are under 5 kb; contig_4 is circular
    # (assembly_info) and has no terminal repeat
    assert (rep["unbinned_already_binned_skipped"], rep["unbinned_binned_by_name_skipped"],
            rep["unbinned_short_dropped"], rep["unbinned_contigs_kept"]) == (1, 1, 2, 1)
    assert rep["dtr_trimmed"] == 0  # contig_3's trim is not recorded: it is not staged as unbinned
    assert (out / "trims.tsv").read_text().splitlines() == ["genome_id\tcontig_name\traw_len\ttrimmed_k\tfinal_len"]
    assert g[-1]["circular"] == [uid(SEQ["c4"])] and g[-1]["assembly_quality"] == "contig"
    # the codes C-66 gave the sources it read right are kept (s1: 25); the bins are screened (under 100 kb: 11)
    assert [x["genetic_code"] for x in g] == [11, 25, 11, 11, 11, 11]
    assert rep["genomes_by_source"] == {"bickhart_2022_rumen_hifi": 2, "bickhart_2022_rumen_hifi_unbinned": 1,
                                        "emerge_stordalen_2019_hifi": 1, "singleton_2021_wwtp_ont": 2}
    assert rep["total_bp"] == sum(len(s) for s in want)
    # C-66's rows of the two sources come out identical; its bins rows are not compared
    assert rep["preset_rows_compared"] == 4 and rep["preset_rows_first_difference"] is None
    o = json.loads((out / "objects.json").read_text())
    assert o[0]["path"] == "/tmp/if2453/s2/inputs/LR2/lr2_r19.fna.gz" and o[0]["source_release"] == "2026-11-27"
    assert json.loads((out / "objects_staging.json").read_text())[0]["path"] == rep["fasta_path"]


@pytestmark_r19
def test_r19_restage_is_byte_identical_when_run_again(tmp_path: Path):
    m = _mod()
    w = lr2_world(tmp_path)
    a, b = _restage(m, w, "a"), _restage(m, w, "b")
    assert a["fasta_sha256"] == b["fasta_sha256"]
    for n in ("genomes.json", "records.tsv", "genetic_codes.tsv", "trims.tsv", "repeated_records.tsv",
              "shared_records.tsv"):
        assert sha(tmp_path / "a" / n) == sha(tmp_path / "b" / n), n


@pytestmark_r19
def test_r19_restage_refuses_a_byte_that_is_not_sequence_and_an_existing_target(tmp_path: Path):
    m = _mod()
    w = lr2_world(tmp_path, bad_bin=True)
    with pytest.raises(m.Refused, match="bin.2.fa.*IUPAC"):
        _restage(m, w)
    w2 = lr2_world(tmp_path / "w2")
    _restage(m, w2)
    with pytest.raises(m.Refused, match="exists"):
        _restage(m, w2)


@pytestmark_r19
def test_r19_restage_reports_a_row_of_c66_that_comes_out_different(tmp_path: Path):
    m = _mod()
    w = lr2_world(tmp_path)
    p = w["old"] / "records.tsv"
    p.write_bytes(p.read_bytes().replace(b"\t3000\t", b"\t2999\t"))
    rep = _restage(m, w)
    assert rep["preset_rows_first_difference"]["row"] == 1
    rc = m.main(["restage-lr2", str(w["dl"]), str(w["meta"]), str(w["old"]), str(tmp_path / "cli"),
                 str(tmp_path / "cli" / "x.fna.gz"), "/tmp/if2453/s2/inputs/LR2", "--compressor", "gzip"])
    assert rc == 1


# --- card R-20: NX1 re-staged without the three 16S alignment files of the Lost City archive -------------------------
#
# The fixture stage has card C-66's layout (records.tsv and the other tables with CRLF, genomes.json compact with one
# newline at the end, contig-uid headers wrapped at 80). Truth by hand: which genome is an alignment file (its records hold gap
# bytes `-`), which contig a kept genome shares with one of them, and the bytes the kept records must keep.

pytestmark_r20 = pytest.mark.covers("R-20")

ALN = ("lost_city:Bipolaricaulota-16S-phylogeny", "lost_city:Methanosarcinales-16S-phylogeny",
       "lost_city:Thermodesulfovibrionales-16S-phylogeny")


def gapped(tag: str, n: int, gaps: int) -> bytes:
    """An aligned row: `gaps` gap bytes spread through n letters (as the Lost City alignments hold them)."""
    s = bytearray(dna(tag, n))
    for i in range(gaps):
        s[(i * 37 + 5) % n] = 0x2D
    return bytes(s)


NX_SEQ = {"a": dna("nx-a", 300), "b": dna("nx-b", 170), "c": dna("nx-c", 81), "x": dna("nx-x", 160),
          "e": dna("nx-e", 240), "f": dna("nx-f", 95),
          "al1": gapped("al1", 200, 7), "al2": gapped("al2", 200, 3), "al3": gapped("al3", 120, 11),
          "al4": gapped("al4", 90, 2)}
# genome, source, contigs (in C-66's order); contig x is in a MAG and in an alignment "genome" (it stays: a kept genome
# holds it); the three alignment genomes are interleaved with real MAGs, as in C-66's records.tsv
NX_GENOMES = [
    ("ehomd:GCA_000000011.1", "ehomd_v11_03_genomes", ["a", "b"]),
    ("lost_city:BinA", "lost_city_serpentinizing_mags", ["c", "x"]),
    (ALN[0], "lost_city_serpentinizing_mags", ["al1", "al2", "x"]),
    (ALN[1], "lost_city_serpentinizing_mags", ["al3"]),
    ("samail:GCF_000000033.2", "ena_c66_nx1_holobiont_subsurface", ["e"]),
    (ALN[2], "lost_city_serpentinizing_mags", ["al4"]),
    ("rhizo:MAG1", "rhizosphere_colombia_mags", ["f"]),
]
GC_HEAD = ["genome_id", "source_name", "genome_bp", "n_contigs", "assembly_quality", "genetic_code", "reason",
           "tested_alt", "cd_11", "mean_orf_11", "cd_25", "cd_4"]


def _crlf(rows) -> bytes:
    return b"".join(("\t".join(r) + "\r\n").encode() for r in rows)


def _fa80(tags) -> bytes:
    out = b""
    for t in tags:
        s = NX_SEQ[t]
        out += b">" + uid(s).encode() + b"\n" + b"".join(s[i:i + 80] + b"\n" for i in range(0, len(s), 80))
    return out


def _rec_rows(genomes):
    return [[g, src, str(i), uid(NX_SEQ[t]), str(len(NX_SEQ[t])), "0", "prokaryote", "mag_hq", "11", "tx"]
            for g, src, tags in genomes for i, t in enumerate(tags)]


def _gc_rows(genomes):
    return [[g, src, str(sum(len(NX_SEQ[t]) for t in tags)), str(len(tags)), "mag_hq", "11", "below_100kb", "False",
             "0.0", "0.0", "", ""] for g, src, tags in genomes]


def _gj(genomes) -> bytes:
    return json.dumps([{"contig_uids": [uid(NX_SEQ[t]) for t in tags], "genetic_code": 11, "kind": "prokaryote",
                        "circular": [], "assembly_quality": "mag_hq"} for _, _, tags in genomes],
                      separators=(",", ":")).encode() + b"\n"  # compact, one newline at the end, as C-66 wrote it


def _first_seen(genomes):
    seen = []
    for _, _, tags in genomes:
        for t in tags:
            if t not in seen:
                seen.append(t)
    return seen


def nx1_world(tmp: Path, *, genomes=NX_GENOMES, skani_rows=None) -> dict:
    old = tmp / "c66_NX1"
    old.mkdir(parents=True)
    (old / "records.tsv").write_bytes(_crlf([RECORDS_HEAD, *_rec_rows(genomes)]))
    (old / "genetic_codes.tsv").write_bytes(_crlf([GC_HEAD, *_gc_rows(genomes)]))
    (old / "genomes.json").write_bytes(_gj(genomes))
    (old / "shared_records.tsv").write_bytes(_crlf([["contig_uid", "owner_stage", "length"],
                                                    [uid(NX_SEQ["a"]), "PL1", "300"],
                                                    [uid(NX_SEQ["e"]), "MV1", "240"]]))
    (old / "trims.tsv").write_bytes(_crlf([["genome_id", "contig_name", "raw_len", "trimmed_k", "final_len"]]))
    (old / "repeated_records.tsv").write_bytes(_crlf([["genome_id", "source_name", "reason", "first_contig_uid"]]))
    skani = skani_rows if skani_rows is not None else [
        ["rhizo:MAG9", "rhizosphere_colombia_mags", "skani_within_stage_ge_96.0ani_rep:rhizo:MAG1"]]
    (old / "skani_within_dropped.tsv").write_bytes(_crlf([["genome_id", "source_name", "reason"], *skani]))
    obj = {"path": "/tmp/if2453/s2/inputs/NX1/nx1_compendia.fna.gz", "source_name": "c66_nx1_compendia",
           "source_release": "2026-11-24", "license": "CC BY 4.0 / other-open / INSDC", "redistributable": True,
           "retrieved_at": "2026-11-24T00:00:00Z"}
    (old / "objects.json").write_text(json.dumps([obj], indent=2) + "\n")
    (old / "objects_staging.json").write_text(json.dumps(
        [{**obj, "path": "/stor/x/NX1/nx1_compendia.fna.gz"}], indent=2) + "\n")
    fa = tmp / "nx1_compendia.fna.gz"
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(fa, "wb"), mtime=0) as g:
        g.write(_fa80(_first_seen(genomes)))
    return {"old": old, "fasta": fa, "tmp": tmp}


def _restage_nx1(m, w, name="out", **kw):
    return m.restage_nx1(w["old"], w["fasta"], w["tmp"] / name, w["tmp"] / name / "fa" / "nx1_compendia_r20.fna.gz",
                         "/tmp/if2453/s2/inputs/NX1", compressor="gzip", **kw)


KEPT = [g for g in NX_GENOMES if g[0] not in ALN]


@pytestmark_r20
def test_r20_restage_drops_the_alignment_files_and_keeps_every_other_byte(tmp_path: Path):
    m = _mod()
    w = nx1_world(tmp_path)
    rep = _restage_nx1(m, w)
    out = tmp_path / "out"
    new_fa = Path(rep["fasta_path"])
    # the kept records, byte for byte as the old file holds them, in its order; x stays (lost_city:BinA holds it)
    with gzip.open(new_fa, "rb") as f:
        assert f.read() == _fa80(["a", "b", "c", "x", "e", "f"])
    scan = _load("r19_fasta_scan", "r19_fasta_scan.py").scan_file(new_fa)
    assert scan["status"] == "clean" and scan["records"] == 6 and scan["non_iupac_bytes"] == 0
    # every table: C-66's rows less those of the three genomes, bytes and line ends kept
    assert (out / "records.tsv").read_bytes() == _crlf([RECORDS_HEAD, *_rec_rows(KEPT)])
    assert (out / "genetic_codes.tsv").read_bytes() == _crlf([GC_HEAD, *_gc_rows(KEPT)])
    assert (out / "genomes.json").read_bytes() == _gj(KEPT)
    for n in ("shared_records.tsv", "trims.tsv", "repeated_records.tsv", "skani_within_dropped.tsv"):
        assert (out / n).read_bytes() == (w["old"] / n).read_bytes(), n
    o = json.loads((out / "objects.json").read_text())
    assert o[0]["path"] == "/tmp/if2453/s2/inputs/NX1/nx1_compendia_r20.fna.gz"
    assert o[0]["source_release"] == "2026-11-27" and o[0]["source_name"] == "c66_nx1_compendia"
    assert json.loads((out / "objects_staging.json").read_text())[0]["path"] == str(new_fa)
    # what was removed, and why
    rm = (out / "removed_genomes.tsv").read_text().splitlines()
    assert rm[0].split("\t") == ["genome_id", "source_name", "records", "bases", "gap_bytes", "reason"]
    assert [r.split("\t")[:5] for r in rm[1:]] == [
        [ALN[0], "lost_city_serpentinizing_mags", "3", str(200 + 200 + 160), "10"],
        [ALN[1], "lost_city_serpentinizing_mags", "1", "120", "11"],
        [ALN[2], "lost_city_serpentinizing_mags", "1", "90", "2"]]
    kept_bp = sum(len(NX_SEQ[t]) for t in "abcxef")
    assert (rep["total_genomes"], rep["removed_genomes"], rep["removed_records_tsv_rows"]) == (4, 3, 5)
    assert (rep["distinct_contigs"], rep["total_bp"]) == (6, kept_bp)
    assert (rep["removed_contigs"], rep["removed_bases"], rep["removed_gap_bytes"]) == (4, 610, 23)
    assert rep["removed_contigs_kept_for_another_genome"] == 1
    assert (rep["old_distinct_contigs"], rep["old_total_bp"]) == (10, kept_bp + 610)
    assert rep["kept_records_byte_identical"] is True and rep["uid_mismatches"] == 0
    assert (rep["shared_contigs"], rep["shared_bp"]) == (2, 540)
    assert (rep["net_new_contigs"], rep["net_new_bp"]) == (4, kept_bp - 540)
    assert rep["genomes_by_source"] == {"ehomd_v11_03_genomes": 1, "ena_c66_nx1_holobiont_subsurface": 1,
                                        "lost_city_serpentinizing_mags": 1, "rhizosphere_colombia_mags": 1}
    assert rep["genetic_codes"] == {"11": 4}
    assert sha(w["fasta"]) == rep["old_fasta_sha256"]  # the staged file is only read


@pytestmark_r20
def test_r20_restage_is_byte_identical_when_run_again_and_leaves_the_old_stage_alone(tmp_path: Path):
    m = _mod()
    w = nx1_world(tmp_path)
    before = {p.name: sha(p) for p in w["old"].iterdir()} | {"fasta": sha(w["fasta"])}
    a, b = _restage_nx1(m, w, "a"), _restage_nx1(m, w, "b")
    assert before == {p.name: sha(p) for p in w["old"].iterdir()} | {"fasta": sha(w["fasta"])}
    assert a["fasta_sha256"] == b["fasta_sha256"]
    for p in sorted((tmp_path / "a").iterdir()):
        if p.is_file() and p.name not in ("objects_staging.json", "stage_report.json"):  # these name their own paths
            assert sha(p) == sha(tmp_path / "b" / p.name), p.name


@pytestmark_r20
def test_r20_restage_refusals(tmp_path: Path):
    m = _mod()
    w = nx1_world(tmp_path / "w")
    _restage_nx1(m, w)
    with pytest.raises(m.Refused, match="exists"):
        _restage_nx1(m, w)
    with pytest.raises(m.Refused, match="not in .*records.tsv"):
        _restage_nx1(m, w, "o2", drop=(*ALN, "lost_city:NoSuchBin"))
    # a kept record holding a gap byte is refused, not written
    bad = [*NX_GENOMES[:-1], ("rhizo:MAG1", "rhizosphere_colombia_mags", ["al4"])]
    w2 = nx1_world(tmp_path / "w2", genomes=[g for g in bad if g[0] != ALN[2]])
    with pytest.raises(m.Refused, match="IUPAC"):
        _restage_nx1(m, w2, drop=ALN[:2])
    assert not list((tmp_path / "w2" / "out" / "fa").glob("*"))  # nothing left behind, not even the .tmp
    # a genome C-66's skani step dropped for an alignment "representative" would be lost: refused
    w3 = nx1_world(tmp_path / "w3", skani_rows=[["lost_city:BinZ", "lost_city_serpentinizing_mags",
                                                 f"skani_within_stage_ge_96.0ani_rep:{ALN[1]}"]])
    with pytest.raises(m.Refused, match="skani"):
        _restage_nx1(m, w3)
    # genomes.json that does not follow records.tsv genome for genome
    w4 = nx1_world(tmp_path / "w4")
    (w4["old"] / "genomes.json").write_bytes(_gj(list(reversed(NX_GENOMES))))
    with pytest.raises(m.Refused, match="genomes.json"):
        _restage_nx1(m, w4)


@pytestmark_r20
def test_r20_the_command_line_restages_and_exits_zero(tmp_path: Path, capsys):
    m = _mod()
    w = nx1_world(tmp_path)
    rc = m.main(["restage-nx1", str(w["old"]), str(w["fasta"]), str(tmp_path / "cli"),
                 str(tmp_path / "cli" / "nx1_compendia_r20.fna.gz"), "/tmp/if2453/s2/inputs/NX1",
                 "--compressor", "gzip"])
    assert rc == 0
    rep = json.loads(capsys.readouterr().out)
    assert rep["removed_genomes"] == 3 and rep["total_genomes"] == 4 and rep["card"] == "R-20"
    assert rep["removed_genome_ids"] == sorted(ALN)  # the default list is the three Lost City alignment files


@pytestmark_r20
def test_r20_ingest_stores_a_gap_byte_as_sequence(tmp_path: Path):
    """Behaviour as it is (no change to src/): `biobox.fasta.iter_records` has no alphabet check and ingest hashes and
    stores the record's bytes, so a `-` is neither refused nor dropped. It is stored, counted in the contig's length,
    and part of its contig_uid; the reverse strand keeps it; the gene caller takes the sequence without an error."""
    import pyrodigal

    from biobox import fasta, ids
    from biobox.l0seq import SourceObject, store

    tail = dna("gv-tail", 400) + b"--"           # as GV1's `Hoko_1` ends
    inner = dna("gv-in1", 300) + b"-" + dna("gv-in2", 299)
    clean = dna("gv-clean", 500)
    src = tmp_path / "src"
    src.mkdir()
    p = src / "gv.fna"
    p.write_bytes(fasta_bytes([("Hoko_1 x", tail), ("inner", inner), ("clean", clean)]))
    assert [r.seq for r in fasta.iter_records(p)] == [tail, inner, clean]  # the reader keeps the bytes
    root = tmp_path / "root"
    store.ingest([SourceObject(p, "synth", "1", "public", True)], root, "S0")
    s = store.SeqStore(root, ["S0"])
    try:
        for seq in (tail, inner, clean):
            u = ids.contig_uid(seq)
            assert u == blake3(seq.upper()).digest(length=16)  # the uid of the bytes, gap bytes included
            assert s.get(u) == seq
            assert s.length_of(s.default_digest(u)) == len(seq)
        assert s.get(ids.contig_uid(tail), 398, 402) == tail[398:]
        with pytest.raises(KeyError):  # the same letters without the gaps are another contig, not stored
            s.get(ids.contig_uid(tail.replace(b"-", b"")))
        comp = bytes.maketrans(b"ACGT", b"TGCA")
        assert s.get(ids.contig_uid(inner), strand=-1) == inner.translate(comp)[::-1]  # '-' kept on the - strand
        assert s.get(ids.contig_uid(inner), strand=-1)[299:300] == b"-"
    finally:
        s.close()
    pyrodigal.GeneFinder(meta=True).find_genes(tail + inner)  # the gene caller does not refuse it either

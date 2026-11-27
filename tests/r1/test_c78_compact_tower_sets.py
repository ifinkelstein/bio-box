"""Card C-78 (DEC-79): stage the compact `tower` sets `LRP1`, `MIR2`, `SYM2` and `DPN1` (`MGE4` leaves `r26.11.5`).

Oracles are written here separately from the script: contig_uid by blake3 directly, terminal repeats planted by hand,
the archive members and their kinds written by the fixtures, skani's output tables written by hand (skani is not run:
GitHub Actions has no pod tool). Each test runs twice, with the normal PATH and with PATH reduced to `/usr/bin:/bin`
(as tests/r1/test_r232_index.py does).
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import json
import random
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-78")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-26_c78_compact_inputs.md"
INPUTS_HOME = Path("/tmp/if2453/s2/inputs")


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
    return _load("c78_compact_tower_sets", "c78_compact_tower_sets.py")


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def cuid(seq: str) -> str:
    return blake3(seq.upper().encode()).digest(length=16).hex()


def seq_of(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def fa_text(records: list[tuple[str, str]], width: int = 60) -> str:
    out = io.StringIO()
    for hdr, s in records:
        out.write(f">{hdr}\n")
        for i in range(0, len(s), width):
            out.write(s[i:i + width] + "\n")
    return out.getvalue()


def write_tgz(path: Path, members: dict[str, bytes]) -> None:
    """A deterministic tar.gz: members in the given order, mtime 0."""
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w", format=tarfile.PAX_FORMAT) as tf:
        dirs = sorted({str(Path(n).parent) for n in members if str(Path(n).parent) not in (".", "")})
        for d in dirs:
            ti = tarfile.TarInfo(d + "/")
            ti.type, ti.mtime, ti.mode = tarfile.DIRTYPE, 0, 0o755
            tf.addfile(ti)
        for name, data in members.items():
            ti = tarfile.TarInfo(name)
            ti.size, ti.mtime, ti.mode = len(data), 0, 0o644
            tf.addfile(ti, io.BytesIO(data))
    path.write_bytes(gzip.compress(raw.getvalue(), mtime=0))


def write_zip(path: Path, members: dict[str, bytes]) -> None:
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in members.items():
            zi = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
            z.writestr(zi, data)


def dtr(body: str, k: int) -> str:
    """A circular record as an assembler leaves it: the first k bases repeated at the end."""
    return body + body[:k]


def common(tmp: Path, tag: str, held: Path | None = None, later: Path | None = None) -> dict:
    return dict(root=None, root_stages=(), owner_inputs={"PL2": held} if held else {},
                later_inputs={"ECE2": later} if later else {}, inputs_home=INPUTS_HOME,
                retrieved_at="2026-11-25T22:00:00Z", work=tmp / f"work_{tag}", threads=1, memory_limit="1GB")


def tsv_rows(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check_inputs(stage: str, meta: Path) -> dict:
    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), stage, str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    doc = json.loads(chk.stdout)
    assert doc["problems"] == []
    return doc


# --- census --------------------------------------------------------------------------------------------------------

def test_census_names_every_member_and_tells_genome_dna_from_protein_and_alignments(tmp_path):
    mod = _mod()
    t = tmp_path / "a.tar.gz"
    write_tgz(t, {
        "x/genomes.fasta": fa_text([("g1", seq_of(1, 500)), ("g2", seq_of(2, 300))]).encode(),
        "x/prot.faa": fa_text([("p1", "MKVLAAGIRRW" * 5)]).encode(),
        "x/aln.fasta": fa_text([("a1", "ACGT--ACGT" * 5)]).encode(),
        "x/protein_named_fasta.fasta": fa_text([("p2", "MSTNPKPQRKTKRNTNRRPQDVKFPGG" * 3)]).encode(),
        "x/table.tsv": b"a\tb\n",
    })
    rep = mod.census([t])[str(t)]
    by = {m["name"]: m for m in rep["members"]}
    assert set(by) == {"x/", "x/genomes.fasta", "x/prot.faa", "x/aln.fasta", "x/protein_named_fasta.fasta",
                       "x/table.tsv"}
    assert by["x/genomes.fasta"]["alphabet"] == "dna"
    assert (by["x/genomes.fasta"]["records"], by["x/genomes.fasta"]["bases"]) == (2, 800)
    assert by["x/prot.faa"]["alphabet"] == "protein"
    assert by["x/aln.fasta"]["alphabet"] == "alignment"
    assert by["x/protein_named_fasta.fasta"]["alphabet"] == "protein"
    assert by["x/table.tsv"]["alphabet"] == "not_fasta" and by["x/"]["alphabet"] == "directory"


def test_census_reads_zips_and_plain_fasta(tmp_path):
    mod = _mod()
    z = tmp_path / "g.zip"
    write_zip(z, {"genomes_polinton_virophages.fasta": fa_text([("v1", seq_of(3, 400))]).encode()})
    f = tmp_path / "plain.fasta"
    f.write_text(fa_text([("c1", seq_of(4, 100)), ("c2", seq_of(5, 50))]))
    rep = mod.census([z, f])
    assert rep[str(z)]["members"][0]["alphabet"] == "dna" and rep[str(z)]["members"][0]["records"] == 1
    assert rep[str(f)]["members"] == [{"name": f.name, "bytes": f.stat().st_size, "alphabet": "dna",
                                       "records": 2, "bases": 150, "min_len": 50, "max_len": 100,
                                       "circular_in_header": 0, "first_headers": ["c1", "c2"]}]


def test_extract_copies_named_members_byte_for_byte_and_refuses_a_missing_one(tmp_path):
    mod = _mod()
    t = tmp_path / "s.tar.gz"
    write_tgz(t, {"source_data/data/metadata_mags.tab": b"bin\ttax\nA\tB\n", "source_data/x.ipynb": b"{}"})
    got = mod.extract(t, tmp_path / "out", ["source_data/data/metadata_mags.tab"])
    assert (tmp_path / "out" / "metadata_mags.tab").read_bytes() == b"bin\ttax\nA\tB\n"
    assert got["source_data/data/metadata_mags.tab"]["bytes"] == 12
    with pytest.raises(mod.Refused, match="no member"):
        mod.extract(t, tmp_path / "out2", ["source_data/nope.tab"])
    assert not (tmp_path / "out2").exists() or list((tmp_path / "out2").iterdir()) == []
    with pytest.raises(mod.Refused, match="exists"):
        mod.extract(t, tmp_path / "out", ["source_data/data/metadata_mags.tab"])


# --- DEC-65 (1) on a metagenome-assembled record -------------------------------------------------------------------

def test_dec65_trims_an_exact_terminal_repeat_once_and_keeps_others_whole():
    mod = _mod()
    body = seq_of(7, 3000).encode()
    rep = dtr(body.decode(), 55).encode()
    out, k, cls, circ = mod.dec65(rep, flagged=False)
    assert (out, k, cls, circ) == (body, 55, "exact", True)
    out, k, cls, circ = mod.dec65(body, flagged=False)
    assert (out, k, cls, circ) == (body, 0, "none", False)
    # a record its source flags circular with no repeat left in it is circular as it stands
    out, k, cls, circ = mod.dec65(body, flagged=True)
    assert (out, k, cls, circ) == (body, 0, "none", True)
    # a terminal repeat longer than 1,000 bp is not trimmed and the record stays linear
    long = dtr(seq_of(8, 5000), 1500).encode()
    out, k, cls, circ = mod.dec65(long, flagged=True)
    assert (k, cls, circ) == (0, "long", False) and out == long


# --- LRP1 ----------------------------------------------------------------------------------------------------------

LRP1_RECS = [
    ("S1_NODE_1_length_3055_cov_6.2_cutoff_0_type_circular", dtr(seq_of(11, 3000), 55)),   # trimmed, circular
    ("S1_NODE_2_length_3000_cov_2.8_cutoff_0_type_circular", seq_of(12, 3000)),             # circular as it stands
    ("S2_RNODE_7_length_4045_cov_47.1", seq_of(13, 4045)),                                   # linear
    ("S3_NODE_9_length_1030_cov_2.2", seq_of(14, 1030)),                                     # linear; held by PL2
    ("S3_NODE_10_length_1030_cov_2.2", seq_of(14, 1030)),                                    # exact repeat of NODE_9
    ("S4_NODE_11_length_2000_cov_1.0", seq_of(15, 2000)),                                    # held by ECE2 (later)
    ("S5_NODE_12_length_900_cov_1.0", "ACGTNNXX" * 100),                                     # not DNA: skipped
]


@pytest.fixture()
def lrp1_world(tmp_path: Path) -> dict:
    fa = tmp_path / "metagenome_plasmidome_plasmid2.fasta"
    fa.write_text(fa_text(LRP1_RECS))
    held = tmp_path / "chain" / "PL2"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(LRP1_RECS[3][1])]}]))
    later = tmp_path / "later" / "ECE2"
    later.mkdir(parents=True)
    (later / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(LRP1_RECS[5][1])]}]))
    return {"fa": fa, "held": held, "later": later, "tmp": tmp_path}


def run_lrp1(mod, w: dict, tag: str):
    t = w["tmp"]
    return mod.stage_lrp1(w["fa"], t / f"meta_{tag}", t / f"fasta_{tag}", **common(t, tag, w["held"], w["later"]))


def test_lrp1_stages_every_dna_record_as_a_plasmid_with_dec65_and_is_deterministic(lrp1_world):
    mod = _mod()
    rep = run_lrp1(mod, lrp1_world, "a")
    rep_b = run_lrp1(mod, lrp1_world, "b")
    t = lrp1_world["tmp"]
    meta, fa = t / "meta_a" / "LRP1", t / "fasta_a" / "LRP1" / "lrp1_urban_water_plasmidome.fna"
    genomes = json.loads((meta / "genomes.json").read_text())
    trimmed = seq_of(11, 3000)
    want = [cuid(trimmed), cuid(LRP1_RECS[1][1]), cuid(LRP1_RECS[2][1]), cuid(LRP1_RECS[3][1]),
            cuid(LRP1_RECS[5][1])]
    assert [g["contig_uids"] for g in genomes] == [[u] for u in want]
    assert {g["kind"] for g in genomes} == {"plasmid"} and {g["genetic_code"] for g in genomes} == {11}
    assert [g["circular"] for g in genomes] == [[want[0]], [want[1]], [], [], []]
    assert [g["assembly_quality"] for g in genomes] == ["complete", "complete", "contig", "contig", "contig"]

    text = fa.read_text()
    names = [ln[1:] for ln in text.splitlines() if ln.startswith(">")]
    assert names == [LRP1_RECS[i][0] for i in (0, 1, 2, 3, 5)]  # headers as deposited
    seqs = text.split(">")[1:]
    assert "".join(seqs[0].splitlines()[1:]) == trimmed  # the second copy of the 55 bp repeat is gone

    trims = tsv_rows(meta / "trims.tsv")
    assert [(r["genome"], r["repeat_length"]) for r in trims] == [(LRP1_RECS[0][0], "55")]
    skipped = tsv_rows(meta / "skipped_records.tsv")
    assert [(r["name"], r["reason"]) for r in skipped] == [
        (LRP1_RECS[4][0], f"exact_repeat_of {LRP1_RECS[3][0]}"), (LRP1_RECS[6][0], "not_dna")]
    shared = tsv_rows(meta / "shared_records.tsv")
    assert [(r["genome"], r["owner_stage"]) for r in shared] == [(LRP1_RECS[3][0], "PL2")]
    later = tsv_rows(meta / "later_records.tsv")
    assert [(r["genome"], r["holder"]) for r in later] == [(LRP1_RECS[5][0], "ECE2")]

    assert rep["records_read"] == 7 and rep["genomes"] == 5 and rep["trims"] == 1
    assert rep["owner_contigs"] == {"PL2": 1} and rep["later_contigs"] == {"ECE2": 1}
    assert rep["circular_flagged_in_header"] == 2

    objs = json.loads((meta / "objects.json").read_text())
    assert [o["path"] for o in objs] == [str(INPUTS_HOME / "LRP1" / "lrp1_urban_water_plasmidome.fna")]
    assert objs[0]["redistributable"] is False and "CC BY 4.0" in objs[0]["license"]
    assert json.loads((meta / "objects_staging.json").read_text())[0]["path"] == str(fa)

    skip = "objects_staging.json"
    assert {k: v for k, v in rep["sha256"].items() if k != skip} == {
        k: v for k, v in rep_b["sha256"].items() if k != skip}
    doc = check_inputs("LRP1", meta)
    assert doc["kinds"] == ["plasmid"] and doc["circular_rows"] == 2


def test_lrp1_refuses_an_existing_stage_and_leaves_nothing_when_it_refuses(lrp1_world, tmp_path):
    mod = _mod()
    run_lrp1(mod, lrp1_world, "a")
    with pytest.raises(mod.Refused, match="exists"):
        run_lrp1(mod, lrp1_world, "a")
    bad = tmp_path / "bad.fasta"
    bad.write_text("ACGT\n>x\nACGT\n")
    with pytest.raises(mod.Refused):
        mod.stage_lrp1(bad, tmp_path / "meta_x", tmp_path / "fasta_x", **common(tmp_path, "x"))
    assert not (tmp_path / "meta_x" / "LRP1").exists()
    d = tmp_path / "fasta_x" / "LRP1"
    assert not d.exists() or [p.name for p in d.iterdir()] == []


# --- MIR2 ----------------------------------------------------------------------------------------------------------

PLV1, PLV2, PLV3 = seq_of(21, 20000), seq_of(22, 15000), seq_of(23, 18000)
VIRO = dtr(seq_of(24, 19000), 40)


@pytest.fixture()
def mir2_world(tmp_path: Path) -> dict:
    z = tmp_path / "genomes_polinton_virophages.zip"
    write_zip(z, {"genomes_polinton_virophages.fasta":
                  fa_text([("N-HFG070301-C1323_2898-26282", PLV1), ("R-LAC-C9_1-19040", VIRO)]).encode()})
    t = tmp_path / "PLV_framework_supplementary_data.tar.gz"
    write_tgz(t, {
        "./fasta_unsorted_PLV/IMGVR_PLV_all.fasta": fa_text([("IMGVR_UViG_1_000001|1|2", PLV2),
                                                             ("IMGVR_UViG_1_000002|1|3", PLV3)]).encode(),
        "./fasta_unsorted_PLV/Austrian_lake_PLV_MCP_unsorted.faa": fa_text([("p", "MKLLIEPQ" * 20)]).encode(),
        # the annotated copy of a UViG the unsorted file holds: the same name and bases, listed once
        "./annotated_PLVs/Polintoviricetes_GKS1_20K.fasta": fa_text([("IMGVR_UViG_1_000001|1|2", PLV2)]).encode(),
        "./annotated_PLVs/Polintoviricetes_GKS1_20K.gff": b"##gff-version 3\n",
        "./Virus_cluster_MCP_profiles/VC_MCP_fasta_align/VC_1.fasta": fa_text([("a", "MKL--QP" * 10)]).encode(),
        "./Supplementary_Tables.xlsx": b"PK\x03\x04",
    })
    held = tmp_path / "chain" / "V2"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(PLV3)]}]))
    return {"zip": z, "tgz": t, "held": held, "tmp": tmp_path}


def test_mir2_stages_the_plv_and_virophage_genomes_of_both_archives_as_viruses(mir2_world):
    mod = _mod()
    w, t = mir2_world, mir2_world["tmp"]
    kw = common(t, "a")
    kw["owner_inputs"] = {"V2": w["held"]}
    rep = mod.stage_mir2(w["zip"], w["tgz"], t / "meta_a", t / "fasta_a", **kw)
    meta = t / "meta_a" / "MIR2"
    genomes = json.loads((meta / "genomes.json").read_text())
    assert [g["contig_uids"] for g in genomes] == [[cuid(PLV1)], [cuid(VIRO[:-40])], [cuid(PLV2)], [cuid(PLV3)]]
    assert {g["kind"] for g in genomes} == {"virus"}
    assert [g["circular"] for g in genomes] == [[], [cuid(VIRO[:-40])], [], []]
    fa_dir = t / "fasta_a" / "MIR2"
    assert sorted(p.name for p in fa_dir.iterdir()) == ["mir2_bellas_2026_plv_framework.fna",
                                                         "mir2_bulzu_2025_plv_virophage.fna"]
    members = {r["member"]: r for r in tsv_rows(meta / "skipped_members.tsv")}
    assert members["./fasta_unsorted_PLV/Austrian_lake_PLV_MCP_unsorted.faa"]["reason"] == "protein"
    assert members["./Virus_cluster_MCP_profiles/VC_MCP_fasta_align/VC_1.fasta"]["reason"] == "alignment"
    assert members["./annotated_PLVs/Polintoviricetes_GKS1_20K.gff"]["reason"] == "not_fasta"
    assert members["./Supplementary_Tables.xlsx"]["reason"] == "not_fasta"
    skipped = tsv_rows(meta / "skipped_records.tsv")
    assert [(r["member"], r["reason"]) for r in skipped] == [
        ("./annotated_PLVs/Polintoviricetes_GKS1_20K.fasta", "exact_repeat_of IMGVR_UViG_1_000001|1|2")]
    src = tsv_rows(meta / "record_sources.tsv")
    assert [(r["genome"], r["member"]) for r in src] == [
        ("N-HFG070301-C1323_2898-26282", "genomes_polinton_virophages.fasta"),
        ("R-LAC-C9_1-19040", "genomes_polinton_virophages.fasta"),
        ("IMGVR_UViG_1_000001|1|2", "./fasta_unsorted_PLV/IMGVR_PLV_all.fasta"),
        ("IMGVR_UViG_1_000002|1|3", "./fasta_unsorted_PLV/IMGVR_PLV_all.fasta")]
    assert rep["owner_contigs"] == {"V2": 1} and rep["by_category"]["bellas"]["genomes"] == 2
    doc = check_inputs("MIR2", meta)
    assert doc["kinds"] == ["virus"]


# --- SYM2 ----------------------------------------------------------------------------------------------------------

BIN_A = [seq_of(31, 30000), seq_of(32, 25000)]
BIN_B = [seq_of(33, 22000)]
VIR1 = dtr(seq_of(34, 5000), 33)       # a complete circular virus of the authors' final set
VIR2 = seq_of(35, 1200)                 # a fragment outside the final set
VIR_IN_BIN = BIN_B[0]                   # a viral contig the authors also binned into a MAG


@pytest.fixture()
def sym2_world(tmp_path: Path) -> dict:
    bins = tmp_path / "bins_fna.tar.gz"
    write_tgz(bins, {"bins_fna/LKH041_P01_Ci-25.fna": fa_text([("NODE_1", BIN_A[0]), ("NODE_2", BIN_A[1])]).encode(),
                     "bins_fna/LKH050_P02_Am-3.fna": fa_text([("NODE_9", BIN_B[0])]).encode()})
    vir = tmp_path / "viruses.fna.tar.gz"
    write_tgz(vir, {"viruses.fna": fa_text([("LKH041_P01_Ci_virus|NODE_57_length_5033_cov_23.7", VIR1),
                                            ("LKH041_P01_Ci_virus|NODE_88_length_1200_cov_2.0", VIR2),
                                            ("LKH050_P02_Am_virus|NODE_9_length_22000_cov_5.0", VIR_IN_BIN)]
                                           ).encode()})
    tables = tmp_path / "tables"
    tables.mkdir()
    (tables / "metadata_mags.tab").write_text(
        "MAGid\tLibrary\tHost\tMIMAG\tCompleteness_CheckM2\tContamination_CheckM2\tclassification\tsymbclass\n"
        "LKH041_P01_Ci-25\tLKH041\tGeleia sp.\tMQ\t60.1\t1.2\td__Bacteria;p__Pseudomonadota;c__Gamma\tsymbiont\n"
        "LKH050_P02_Am-3\tLKH050\tAmoeba sp.\tHQ\t95.0\t0.5\td__Bacteria;p__Bacillota;o__Mycoplasmatales\tsymbiont\n")
    (tables / "final_virus_df.tsv").write_text(
        "\tseq_name\tlength\ttopology\tcontig_id\tcheckv_quality\ttaxonomy\thost\n"
        "5\tNODE_57_length_5033_cov_23.7\t5033\tDTR\tLKH041_P01_Ci_virus|NODE_57_length_5033_cov_23.7\tComplete\t"
        "Viruses;Monodnaviria\tGeleia sp.\n")
    held = tmp_path / "chain" / "SYM1"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(VIR2)]}]))
    return {"bins": bins, "vir": vir, "tables": tables, "held": held, "tmp": tmp_path}


def fake_screen(seqs, organism_name=""):
    """The code screen's interface (card C-58's), without pyrodigal: code 4 for a Mycoplasma lineage."""
    return {"genetic_code": 4 if "Mycoplasma" in organism_name else 11, "reason": "fixture", "tested_alt": False,
            "cd_11": None, "mean_orf_11": None, "cd_25": None, "cd_4": None}


def test_sym2_stages_the_bins_as_prokaryotes_and_the_viral_contigs_as_viruses(sym2_world):
    mod = _mod()
    w, t = sym2_world, sym2_world["tmp"]
    kw = common(t, "a")
    kw["owner_inputs"] = {"SYM1": w["held"]}
    rep = mod.stage_sym2(w["bins"], w["vir"], w["tables"], t / "meta_a", t / "fasta_a", screen=fake_screen,
                         workers=0, **kw)
    meta = t / "meta_a" / "SYM2"
    genomes = json.loads((meta / "genomes.json").read_text())
    by_first = {g["contig_uids"][0]: g for g in genomes}
    a, b = by_first[cuid(BIN_A[0])], by_first[cuid(BIN_B[0])]
    assert a["contig_uids"] == [cuid(s) for s in BIN_A] and a["kind"] == "prokaryote" and a["genetic_code"] == 11
    assert a["assembly_quality"] == "mag_mq" and b["assembly_quality"] == "mag_hq" and b["genetic_code"] == 4
    v1 = by_first[cuid(VIR1[:-33])]
    assert v1["kind"] == "virus" and v1["circular"] == [cuid(VIR1[:-33])]
    assert by_first[cuid(VIR2)]["kind"] == "virus"
    assert len(genomes) == 4  # the binned viral contig is not listed a second time as a virus
    skipped = tsv_rows(meta / "skipped_records.tsv")
    assert [(r["name"], r["reason"]) for r in skipped] == [
        ("LKH050_P02_Am_virus|NODE_9_length_22000_cov_5.0", "contig_in_bin LKH050_P02_Am-3")]
    vt = {r["genome"]: r for r in tsv_rows(meta / "sym2_viruses.tsv")}
    assert vt["LKH041_P01_Ci_virus|NODE_57_length_5033_cov_23.7"]["in_final_set"] == "yes"
    assert vt["LKH041_P01_Ci_virus|NODE_88_length_1200_cov_2.0"]["in_final_set"] == "no"
    st = {r["genome"]: r["r232_status"] for r in tsv_rows(meta / "r232_status.tsv")}
    assert st == {"LKH041_P01_Ci-25": "not_visible", "LKH050_P02_Am-3": "not_visible"}
    codes = {r["genome"]: r["genetic_code"] for r in tsv_rows(meta / "genetic_codes.tsv")}
    assert codes == {"LKH041_P01_Ci-25": "11", "LKH050_P02_Am-3": "4"}
    assert rep["owner_contigs"] == {"SYM1": 1}
    doc = check_inputs("SYM2", meta)
    assert doc["kinds"] == ["prokaryote", "virus"] and doc["codes"] == [4, 11]


# --- the conf templates --------------------------------------------------------------------------------------------

# stage: (depth, library, AFTER, the staged FASTA names)
STAGED = {
    "LRP1": ("annotated", ("pfam", "38.2", "gathering"), ("SYM2",), {"lrp1_urban_water_plasmidome.fna"}),
}


@pytest.mark.parametrize("stage", sorted(STAGED))
def test_template_passes_stage_est_and_carries_the_cards_values(stage):
    est = _load("stage_est", "stage_est.py")
    depth, library, after, fastas = STAGED[stage]
    c = est.read_conf(TEMPLATES / f"{stage}.conf", stage, tree=REPO)
    assert (c.depth, c.crispr, c.library, c.after) == (depth, True, library, after)
    assert str(c.inputs) == f"/tmp/if2453/s2/inputs/{stage}" and c.min_contig_len is None
    assert NOTE in c.notes
    assert {n for _, n in c.sums} == {"genomes.json", "objects.json", *fastas}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    for h, _ in c.sums:
        assert h in note, f"SUM {h} not in the note"
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    allowed = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}
    allowed |= {s for s, _ in a26.R26_11_4_CHAIN}
    assert set(c.expect) <= allowed
    for owner, (n, _) in c.expect.items():
        assert n > 0


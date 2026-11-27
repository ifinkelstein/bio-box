"""Card C-81 (DEC-79): stage the Groundwater Virome Catalogue (`GWV1`) and the OMZ single cells (`SAG3`).

Oracles are written here separately from the script: contig_uid by blake3 directly, the DEC-65 trim by building each
fixture record from a known core and repeat, membership of GTDB R232 and GlobDB r232 by hand-written fixture tables,
the matching of archive members to Table S1 by hand, and Table S1 itself as an xlsx written here from raw XML. Each
test runs twice, with the normal PATH and with PATH reduced to `/usr/bin:/bin` (as tests/r1/test_r232_index.py does):
the stager needs only Python and its libraries, and GitHub Actions has none of the pod's tools; the one test that runs
skani skips cleanly without it.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import json
import random
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-81")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-26_gwv1_sag3_wet1_inputs.md"
ROOT13 = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}


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
    return _load("c81_gwv1_sag3_wet1", "c81_gwv1_sag3_wet1.py")


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


def _fa(records: list[tuple[str, str]], width: int = 60) -> str:
    out = io.StringIO()
    for hdr, s in records:
        out.write(f">{hdr}\n")
        for i in range(0, len(s), width):
            out.write(s[i:i + width] + "\n")
    return out.getvalue()


def _gz_lines(path: Path) -> list[str]:
    return gzip.decompress(path.read_bytes()).decode().splitlines()


def _rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ==================================================================================================================
# GWV1
# ==================================================================================================================

CORE_A = seq_of(101, 3000)
REP_A = CORE_A[:30]                       # an exact 30 bp direct terminal repeat: trimmed, circular
SEQ_A = CORE_A + REP_A
SEQ_B = seq_of(102, 2500)                 # flagged DTR, no repeat at the ends: circular, untrimmed (`none`)
CORE_C = seq_of(103, 2800)
SEQ_C = CORE_C + CORE_C[:30]              # the same kind of repeat, but CheckV says AAI-based: linear, untrimmed
SEQ_D = seq_of(104, 2600)                 # ITR: linear
SEQ_E = seq_of(105, 2400)                 # E and F carry the same bytes (DEC-27: listed once)
SEQ_G = seq_of(107, 2200)                 # held by V2 (an r26.11.4 stage)
SEQ_H = seq_of(108, 2100)                 # held by MV2 only (an r26.11.5 stage before GWV1)
GW = [("k141_1_length_3030_cov_9.0_B10", SEQ_A, "Complete", "DTR (high-confidence)"),
      ("k141_2_length_2500_cov_8.0_B10", SEQ_B, "High-quality", "DTR (medium-confidence)"),
      ("k141_3_length_2830_cov_7.0_B12", SEQ_C, "Medium-quality", "AAI-based (high-confidence)"),
      ("k141_4_length_2600_cov_6.0_B12", SEQ_D, "Complete", "ITR (high-confidence)"),
      ("k141_5_length_2400_cov_5.0_SD47", SEQ_E, "Low-quality", "HMM-based (lower-bound)"),
      ("k141_6_length_2400_cov_5.0_SD48", SEQ_E, "Low-quality", "HMM-based (lower-bound)"),
      ("k141_7_length_2200_cov_4.0_JX3", SEQ_G, "Not-determined", "NA"),
      ("k141_8_length_2100_cov_3.0_JX3", SEQ_H, "Low-quality", "AAI-based (medium-confidence)")]
INFO_HEAD = ["contig_id", "contig_length", "checkv_quality", "miuvig_quality", "completeness", "completeness_method",
             "gene_count", "viral_genes", "taxonomy", "lifestyle", "sample", "zone", "bioproject_accession",
             "biosample_accession", "SRA_accession"]


def write_gw(tmp: Path, recs=GW, *, info_recs=None, lengths=None) -> tuple[Path, Path]:
    fna = tmp / "GWVC.vOTUs.fna.gz"
    fna.write_bytes(gzip.compress(_fa([(n, s) for n, s, _, _ in recs], width=70).encode(), mtime=0))
    info = tmp / "GWVC.vOTUs.information.txt"
    rows = []
    for n, s, q, m in (recs if info_recs is None else info_recs):
        ln = (lengths or {}).get(n, len(s))
        rows.append([n, str(ln), q, "High-quality" if q == "Complete" else "Genome-fragment", "100", m, "9", "2",
                     "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes", "", n.rsplit("_", 1)[1],
                     "Zone_II", "PRJNA858913", "SAMN30550725", "SRR22388083"])
    info.write_text("\n".join("\t".join(r) for r in [INFO_HEAD, *rows]) + "\n")
    return fna, info


@pytest.fixture()
def gw(tmp_path: Path) -> dict:
    fna, info = write_gw(tmp_path)
    held = tmp_path / "chain" / "V2"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(SEQ_G)]}]))
    later = tmp_path / "later" / "MV2"
    later.mkdir(parents=True)
    (later / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(SEQ_H)]}]))
    return {"fna": fna, "info": info, "held": held, "later": later, "tmp": tmp_path}


def run_gwv1(mod, w: dict, tag: str, **kw):
    t = w["tmp"]
    return mod.stage_gwv1(
        w.get("fna_x", w["fna"]), w.get("info_x", w["info"]), t / f"meta_{tag}", t / f"fasta_{tag}", root=None,
        root_stages=(), owner_inputs={"V2": w["held"]}, later_inputs={"MV2": w["later"]},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=t / f"work_{tag}", **kw)


def test_gwv1_stage_applies_dec65_lists_repeats_once_and_is_deterministic(gw):
    mod = _mod()
    rep1 = run_gwv1(mod, gw, "a")
    rep2 = run_gwv1(mod, gw, "b")
    t = gw["tmp"]
    meta, fa = t / "meta_a" / "GWV1", t / "fasta_a" / "GWV1" / "gwv1_gwvc.fna.gz"

    # the FASTA: BGZF, every record in the file's order, A trimmed by its 30 bp repeat, C (not flagged DTR) untrimmed
    raw = fa.read_bytes()
    assert raw[:4] == bytes.fromhex("1f8b0804") and raw[12:14] == b"BC"
    lines = _gz_lines(fa)
    names = [ln[1:] for ln in lines if ln.startswith(">")]
    assert names == [n for n, *_ in GW]
    seqs, cur = [], []
    for ln in lines:
        if ln.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(ln)
    seqs.append("".join(cur))
    assert seqs[0] == CORE_A and seqs[2] == SEQ_C and seqs[1] == SEQ_B and seqs[3] == SEQ_D

    genomes = json.loads((meta / "genomes.json").read_text())
    want = [cuid(CORE_A), cuid(SEQ_B), cuid(SEQ_C), cuid(SEQ_D), cuid(SEQ_E), cuid(SEQ_G), cuid(SEQ_H)]
    assert [g["contig_uids"] for g in genomes] == [[u] for u in want]
    assert {g["kind"] for g in genomes} == {"virus"} and {g["assembly_quality"] for g in genomes} == {"vOTU"}
    assert [bool(g["circular"]) for g in genomes] == [True, True, False, False, False, False, False]
    assert genomes[0]["circular"] == [cuid(CORE_A)]

    trims = _rows(meta / "trims.tsv")
    assert [(r["genome"], r["source_contig_uid"], r["contig_uid"], r["bases_removed"]) for r in trims] == [
        (GW[0][0], cuid(SEQ_A), cuid(CORE_A), "30")]
    recs = {r["genome"]: r for r in _rows(meta / "records.tsv")}
    assert (recs[GW[0][0]]["dtr_class"], recs[GW[1][0]]["dtr_class"], recs[GW[2][0]]["dtr_class"]) == (
        "exact", "none", "")
    assert recs[GW[0][0]]["source_length"] == "3030" and recs[GW[0][0]]["length"] == "3000"
    rep = _rows(meta / "repeated_records.tsv")
    assert [(r["genome"], r["first_genome"]) for r in rep] == [(GW[5][0], GW[4][0])]
    assert [(r["genome"], r["owner_stage"]) for r in _rows(meta / "shared_records.tsv")] == [(GW[6][0], "V2")]
    assert [(r["genome"], r["holder"]) for r in _rows(meta / "later_records.tsv")] == [(GW[7][0], "MV2")]

    objs = json.loads((meta / "objects.json").read_text())
    assert objs[0]["path"] == "/tmp/if2453/s2/inputs/GWV1/gwv1_gwvc.fna.gz"
    assert objs[0]["redistributable"] is True and "cc-by-4.0" in objs[0]["license"]
    assert json.loads((meta / "objects_staging.json").read_text())[0]["path"] == str(fa)

    r = rep1
    assert r["records"] == 8 and r["genomes"] == 7 and r["genomes_repeated"] == 1
    assert r["dtr_records"] == 2 and r["dtr_classes"] == {"exact": 1, "none": 1} and r["trimmed_bases"] == 30
    assert r["source_bases"] == sum(len(s) for _, s, _, _ in GW) and r["staged_bases"] == r["source_bases"] - 30
    assert r["held_by"] == {"V2": 1} and r["later_held_by"] == {"MV2": 1}
    assert r["bases_new"] == sum(len(s) for s in (CORE_A, SEQ_B, SEQ_C, SEQ_D, SEQ_E, SEQ_H))
    assert r["bases_new_after_later"] == r["bases_new"] - len(SEQ_H)

    for f in ("genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
              "later_records.tsv", "trims.tsv"):
        assert (meta / f).read_bytes() == (t / "meta_b" / "GWV1" / f).read_bytes(), f
    assert raw == (t / "fasta_b" / "GWV1" / "gwv1_gwvc.fna.gz").read_bytes()
    assert {k: v for k, v in rep1["sha256"].items() if k != "objects_staging.json"} == {
        k: v for k, v in rep2["sha256"].items() if k != "objects_staging.json"}

    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), "GWV1", str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    doc = json.loads(chk.stdout)
    assert doc["problems"] == [] and doc["gene_call_mode"] == {"virus": "gv-meta"}


@pytest.mark.parametrize("case", ["not_in_table", "no_record", "length"])
def test_gwv1_refuses_a_table_that_does_not_describe_the_fasta_and_leaves_nothing(gw, case):
    mod = _mod()
    t = gw["tmp"] / case
    t.mkdir()
    if case == "not_in_table":
        fna, info = write_gw(t, info_recs=GW[:-1])
        match = "not in the information table"
    elif case == "no_record":
        fna, info = write_gw(t, GW[:-1], info_recs=GW)
        match = "no FASTA record"
    else:
        fna, info = write_gw(t, lengths={GW[3][0]: 2601})
        match = "2601 in the information table"
    with pytest.raises(mod.Refused, match=match):
        run_gwv1(mod, {**gw, "fna_x": fna, "info_x": info}, case)
    assert not (gw["tmp"] / f"meta_{case}" / "GWV1").exists()
    d = gw["tmp"] / f"fasta_{case}" / "GWV1"
    assert not d.exists() or list(d.iterdir()) == []


def test_gwv1_refuses_an_existing_output(gw):
    mod = _mod()
    run_gwv1(mod, gw, "z")
    with pytest.raises(mod.Refused, match="exists"):
        run_gwv1(mod, gw, "z")


def test_gwv1_refuses_a_later_stage_without_genomes_json(gw, tmp_path):
    mod = _mod()
    with pytest.raises(mod.Refused, match="no genomes.json"):
        mod.stage_gwv1(gw["fna"], gw["info"], tmp_path / "m", tmp_path / "f", root=None, root_stages=(),
                       owner_inputs={}, later_inputs={"MV2": tmp_path / "nowhere"}, inputs_home=Path("/x"))


def test_the_votu_threshold_rules_and_the_sample_hits_on_a_hand_written_skani_table(tmp_path):
    mod = _mod()
    assert mod.pair_meets(95.0, 30.0, 85.0) == {"shorter": True, "both": False, "near_identical": False}
    assert mod.pair_meets(99.5, 92.0, 95.0) == {"shorter": True, "both": True, "near_identical": True}
    assert mod.pair_meets(94.99, 100.0, 100.0)["shorter"] is False
    assert mod.pair_meets(97.0, 84.9, 84.0)["shorter"] is False
    tab = tmp_path / "a.tsv"
    head = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query", "Ref_name", "Query_name"]
    rows = [["s.fna", "mv1.fna", "96.0", "90.0", "20.0", "gw1", "mv1_a"],     # gw1 meets `shorter` only
            ["s.fna", "mv1.fna", "99.9", "99.0", "98.0", "gw2", "mv1_b"],     # gw2 meets all three
            ["s.fna", "mv1.fna", "90.0", "99.0", "99.0", "gw3", "mv1_c"],     # below 95
            ["s.fna", "s.fna", "99.9", "99.0", "99.0", "gw3", "gw4"]]         # two sample records: not counted
    tab.write_text("\n".join("\t".join(r) for r in [head, *rows]) + "\n")
    h = mod.sample_hits(tab, {"gw1", "gw2", "gw3", "gw4"})
    assert h == {"shorter": {"gw1", "gw2"}, "both": {"gw2"}, "near_identical": {"gw2"}}
    s = mod.ani_summary(h, {"gw1": 5000, "gw2": 12000, "gw3": 40000, "gw4": 6000},
                        {"gw1": "Low-quality", "gw2": "Complete", "gw3": "Low-quality", "gw4": "Complete"})
    assert s["queries"] == 4 and s["shorter"]["meet"] == 2
    assert s["shorter"]["by_checkv_quality"] == {"Complete": {"queries": 2, "meet": 1},
                                                 "Low-quality": {"queries": 2, "meet": 1}}


def test_gwv1_ani_with_skani_when_the_image_has_it(tmp_path):
    if shutil.which("skani") is None:
        pytest.skip("skani is not on PATH (GitHub Actions, or the tools-hidden run)")
    mod = _mod()
    big = [seq_of(400 + i, 20_000) for i in range(3)]
    recs = [(f"k141_{90 + i}_length_20000_cov_9.0_B10", s, "Low-quality", "AAI-based (high-confidence)")
            for i, s in enumerate(big)]
    fna, info = write_gw(tmp_path, recs)
    mod.stage_gwv1(fna, info, tmp_path / "meta_s", tmp_path / "fasta_s", root=None, root_stages=(), owner_inputs={},
                   later_inputs={}, inputs_home=Path("/x"), work=tmp_path / "work_s")
    meta, fa = tmp_path / "meta_s" / "GWV1", tmp_path / "fasta_s" / "GWV1" / "gwv1_gwvc.fna.gz"
    ref = tmp_path / "held.fna"
    ref.write_text(_fa([("held_copy_of_1", big[1]), ("held_other", seq_of(999, 20_000))]))
    # a chunk of 15 kb puts the two held records in two chunk files: two skani runs, the hit found in the first
    res = mod.gwv1_ani(meta, fa, tmp_path / "ani", {"V3": [ref]}, every=1, threads=1, chunk_bytes=15_000)
    assert res["samples"] == 3 and res["by_stage"]["V3"]["skani_runs"] == 2
    assert res["by_stage"]["V3"]["near_identical"]["meet"] == 1
    assert res["any_held_stage"]["shorter"]["meet"] == 1
    assert not list((tmp_path / "ani" / "work").glob("chunk_*"))  # every chunk file is removed


def test_fasta_chunks_cut_between_records_and_remove_each_file(tmp_path):
    mod = _mod()
    a, b = tmp_path / "a.fna", tmp_path / "b.fna"
    a.write_text(_fa([("r1", "A" * 10), ("r2", "C" * 10)]))
    b.write_text(_fa([("r3", "G" * 25)]))
    seen = []
    for p in mod.fasta_chunks([a, b], 15, tmp_path):
        seen.append([ln for ln in p.read_text().splitlines() if ln.startswith(">")])
    assert seen == [[">r1", ">r2"], [">r3"]]
    assert not list(tmp_path.glob("chunk_*"))


# ==================================================================================================================
# SAG3
# ==================================================================================================================

# Archive members (stem -> records) and what Table S1 says of them:
#   AAA001_A19_output_clean  -> AAA001_A19 (prefix); no BioSample, no WGS -> not_visible
#   AAA001_C10_contigs       -> AAA001_C10 (prefix); WGS AFHZ00000000, which GTDB R232 lists (GCA_000000101, a
#                               representative, also in GlobDB r232) -> in_gtdb_and_globdb_r232 via wgs
#   JGI_01_D11               -> JGI_01_D11 (exact); BioSample SAMN09212608 -> GCF_000000102 in GTDB only ->
#                               in_gtdb_r232 via biosample; its contig 2 is held by V2
#   JGI_01_D11_4000m         -> JGI_01_D11_4000m (exact, the longer name wins); BioSample not in GTDB ->
#                               not_in_gtdb_r232; one contig is held by SAG1 (later in the chain)
#   X_01_A1, X_01_A1_rep     -> both would take X_01_A1: tie, neither is matched -> not_visible
#   Unknown_9                -> no name -> not_visible
#   AG-311-D23               -> exact; GTDB-Tk calls it Mycoplasma -> code 4; WGS not in GTDB -> not_in_gtdb_r232
SAG = {
    "AAA001_A19_output_clean": [("AAA001_A19_contig00001_NODE_1", seq_of(201, 3000)),
                                ("AAA001_A19_contig00002_NODE_2", seq_of(202, 1500))],
    "AAA001_C10_contigs": [("AAA001_C10_NODE_1_length_4000", seq_of(203, 4000))],
    "JGI_01_D11": [("JGI01_D11_NODE_1", seq_of(204, 2200)), ("JGI01_D11_NODE_2", seq_of(205, 1300))],
    "JGI_01_D11_4000m": [("JGI_01_D11_4000m_NODE_1", seq_of(206, 2000))],
    "X_01_A1": [("X_01_A1_NODE_1", seq_of(207, 1800))],
    "X_01_A1_rep": [("X_01_A1_rep_NODE_1", seq_of(208, 1700))],
    "Unknown_9": [("Unknown_9_NODE_1", seq_of(209, 1600))],
    "AG-311-D23": [("2634166711_Ga0006891_101", seq_of(210, 2500))],
    # the real Table S1's GTDB-Tk sheet names these two cells `AB-754_C22` twice; its Sequencing sheet in full
    "AB-754_C22_AB-905_clean": [("AB-754_C22_AB-905_NODE_1", seq_of(211, 1900))],
    "AB-754_C22_AB-908_clean": [("AB-754_C22_AB-908_NODE_1", seq_of(212, 1950))],
}
EXTRA_SEQ = {"AB-754_C22_AB-905": ("SAMN32251362", ""), "AB-754_C22_AB-908": ("SAMN32251363", "")}
THIO = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__PS1;f__Thioglobaceae;g__Thioglobus;s__"
EXTRA_GT = [("AB-754_C22", THIO), ("AB-754_C22", THIO)]
TABLE = {  # name -> (GTDB-Tk classification, BioSample, WGS, quality)
    "AAA001_A19": ("d__Archaea;p__Thermoproteota;c__Nitrososphaeria;g__Nitrosopelagicus;s__", "", "", "Low"),
    "AAA001_C10": ("d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;g__Pelagibacter;s__", "", "AFHZ00000000",
                   "Medium"),
    "JGI_01_D11": ("d__Bacteria;p__Actinobacteriota;g__Actinomarina;s__", "SAMN09212608", "", "Medium"),
    "JGI_01_D11_4000m": ("Unclassified Bacteria", "SAMN13170097", "", "Low"),
    "X_01_A1": ("d__Bacteria;p__Marinisomatota;g__Marinisoma;s__", "SAMN00000001", "", "Low"),
    "AG-311-D23": ("d__Bacteria;p__Bacillota;c__Bacilli;o__Mycoplasmatales;g__Mycoplasma;s__", "SAMEA6076706",
                   "CABWYC000000000", "Low"),
    "NOT_IN_ARCHIVE": ("d__Bacteria;s__", "SAMN00000002", "", "Low"),
}
SEQ_HEAD = (["Sample Name", "Sequencing Location", "Sequencing Technology", "Assembly Tool", "Decontamination Tool",
             "Method", "Assay Name", "Raw Data File", "Data Repository", "Data Record (Bioproject)", "Data Repository",
             "Data Record Study ID", "Data Record Sample Accession", "Data Record Experiment ID",
             "Data Record Run Accession", "Data Repository", "Data Record Accession", "Method", "Derived Data File",
             "Data Repository", "Data Record Accession", "Method", "Derived Data File", "Data Repository",
             "Data Record Accession", "Method"])


def _colname(i: int) -> str:
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def _sheet_xml(rows: list[list[str]], shared: dict[str, int]) -> str:
    out = ['<?xml version="1.0" encoding="UTF-8"?>'
           '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>']
    for r, row in enumerate(rows, start=1):
        out.append(f'<row r="{r}">')
        for c, v in enumerate(row):
            if v == "":
                continue
            ref = f"{_colname(c)}{r}"
            v = v.replace("&", "&amp;").replace("<", "&lt;")
            if v in shared:
                out.append(f'<c r="{ref}" t="s"><v>{shared[v]}</v></c>')
            else:
                out.append(f'<c r="{ref}" t="inlineStr"><is><t>{v}</t></is></c>')
        out.append("</row>")
    out.append("</sheetData></worksheet>")
    return "".join(out)


def write_xlsx(path: Path, table=TABLE) -> None:
    """A minimal xlsx: the three sheets the stager reads, in a workbook order that differs from the files' numbering
    (the relationships map them), shared strings for two values and inline strings for the rest."""
    seq_rows = [SEQ_HEAD]
    for n, (bs, wgs) in [*((n, (bs, wgs)) for n, (_, bs, wgs, _) in table.items()), *EXTRA_SEQ.items()]:
        r = [""] * len(SEQ_HEAD)
        r[0], r[9], r[16], r[20], r[24] = n, "PRJNA1", bs, wgs, "2513237067"
        seq_rows.append(r)
    seq_rows.append([""] * 3)
    gt_rows = [["user_genome", "Method", "classification", "fastani_reference"]] + [
        [n, "gtdb-tk v2.1.0", c, "N/A"] for n, (c, _, _, _) in table.items()] + [
        [n, "gtdb-tk v2.1.0", c, "GCA_018644675.1"] for n, c in EXTRA_GT]
    q_rows = [["Sample ID", "Method", "Completeness", "Contamination", "Quality Classification"]] + [
        [n, "CheckM v1.2.1", "50", "1", q] for n, (_, _, _, q) in table.items()]
    shared_list = ["gtdb-tk v2.1.0", "CheckM v1.2.1"]
    shared = {s: i for i, s in enumerate(shared_list)}
    sst = ('<?xml version="1.0" encoding="UTF-8"?><sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/'
           'main">' + "".join(f"<si><t>{s}</t></si>" for s in shared_list) + "</sst>")
    wb = ('<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/'
          '2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
          '<sheet name="Inputs" sheetId="1" r:id="rId1"/>'
          '<sheet name="Outputs-Genome_Sequencing" sheetId="2" r:id="rId3"/>'
          '<sheet name="Outputs-Genome_Quality_Classifi" sheetId="3" r:id="rId4"/>'
          '<sheet name="Outputs_GTDB-TK" sheetId="4" r:id="rId2"/></sheets></workbook>')
    rels = ('<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/'
            '2006/relationships">'
            '<Relationship Id="rId1" Type="x" Target="worksheets/sheet1.xml"/>'
            '<Relationship Id="rId2" Type="x" Target="worksheets/sheet2.xml"/>'
            '<Relationship Id="rId3" Type="x" Target="worksheets/sheet3.xml"/>'
            '<Relationship Id="rId4" Type="x" Target="worksheets/sheet4.xml"/></Relationships>')
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("[Content_Types].xml", "<Types/>")
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", rels)
        z.writestr("xl/sharedStrings.xml", sst)
        z.writestr("xl/worksheets/sheet1.xml", _sheet_xml([["Source Name"], ["x"]], shared))
        z.writestr("xl/worksheets/sheet2.xml", _sheet_xml(gt_rows, shared))
        z.writestr("xl/worksheets/sheet3.xml", _sheet_xml(seq_rows, shared))
        z.writestr("xl/worksheets/sheet4.xml", _sheet_xml(q_rows, shared))


def write_tar(path: Path, members: dict[str, list[tuple[str, str]]], *, extra: dict[str, bytes] | None = None) -> None:
    with tarfile.open(path, "w:gz") as tf:
        d = tarfile.TarInfo("Final_Catalog")
        d.type = tarfile.DIRTYPE
        tf.addfile(d)
        # archive order is not name order: the FASTA follows the archive, the metadata the names
        for stem in sorted(members, reverse=True):
            data = _fa(members[stem]).encode()
            ti = tarfile.TarInfo(f"Final_Catalog/{stem}.fna")
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))
        for name, data in (extra or {}).items():
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))


def write_gtdb_meta(d: Path) -> None:
    d.mkdir()
    head = ["accession", "checkm2_completeness", "ncbi_biosample", "ncbi_wgs_master"]
    ar = [["GB_GCA_000000099.1", "90", "SAMN99999999", "AAAA00000000.1"]]
    bac = [["GB_GCA_000000101.1", "95", "SAMN00260143", "AFHZ00000000.1"],
           ["RS_GCF_000000102.1", "97", "SAMN09212608", "none"]]
    for name, rows in (("ar53_metadata_r232.tsv.gz", ar), ("bac120_metadata_r232.tsv.gz", bac)):
        (d / name).write_bytes(gzip.compress(("\n".join("\t".join(r) for r in [head, *rows]) + "\n").encode(),
                                             mtime=0))


def write_index(d: Path) -> None:
    d.mkdir()
    gt = [("GCA_000000099", "GCA_000000099", True, "GCA_000000099", "Archaeon one"),
          ("GCA_000000101", "GCA_000000101", True, "GCA_000000101", "Pelagibacter sp1"),
          ("GCF_000000102", "GCA_000000102", False, "GCA_000000777", "Actinomarina sp2")]
    pq.write_table(pa.table({
        "accession": [g[0] for g in gt], "genbank_accession": [g[1] for g in gt],
        "is_representative": [g[2] for g in gt], "representative": [g[3] for g in gt],
        "species": [g[4] for g in gt]}), d / "gtdb_r232.parquet")
    pq.write_table(pa.table({
        "globdb_id": ["GCA_000000101", "SPIRE_000001"], "insdc_accession": ["GCA_000000101", None],
        "original_id": [None, None], "cluster_id": [None, None]}), d / "globdb_r232.parquet")


@pytest.fixture()
def sg(tmp_path: Path) -> dict:
    tar = tmp_path / "OMZ_Compendium_SAGs_Full.tar.gz"
    write_tar(tar, SAG)
    xlsx = tmp_path / "Table_S1.xlsx"
    write_xlsx(xlsx)
    gtdb = tmp_path / "gtdb"
    write_gtdb_meta(gtdb)
    idx = tmp_path / "index"
    write_index(idx)
    held = tmp_path / "chain" / "V2"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(SAG["JGI_01_D11"][1][1])]}]))
    later = tmp_path / "later" / "SAG1"
    later.mkdir(parents=True)
    (later / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(SAG["JGI_01_D11_4000m"][0][1])]}]))
    return {"tar": tar, "xlsx": xlsx, "gtdb": gtdb, "index": idx, "held": held, "later": later, "tmp": tmp_path}


def run_sag3(mod, w: dict, tag: str, **kw):
    t = w["tmp"]
    return mod.stage_sag3(
        w["tar"], w["xlsx"], w["index"], w["gtdb"], t / f"meta_{tag}", t / f"fasta_{tag}", root=None, root_stages=(),
        owner_inputs={"V2": w["held"]}, later_inputs={"SAG1": w["later"]}, inputs_home=Path("/tmp/if2453/s2/inputs"),
        work=t / f"work_{tag}", **kw)


GENOME = {"AAA001_A19_output_clean": "AAA001_A19", "AAA001_C10_contigs": "AAA001_C10", "JGI_01_D11": "JGI_01_D11",
          "JGI_01_D11_4000m": "JGI_01_D11_4000m", "X_01_A1": "X_01_A1", "X_01_A1_rep": "X_01_A1_rep",
          "Unknown_9": "Unknown_9", "AG-311-D23": "AG-311-D23",
          "AB-754_C22_AB-905_clean": "AB-754_C22_AB-905", "AB-754_C22_AB-908_clean": "AB-754_C22_AB-908"}
NAMES = set(TABLE) | set(EXTRA_SEQ) | {"AB-754_C22"}


def test_table_s1_is_read_from_the_xlsx_by_sheet_name(sg):
    mod = _mod()
    s1 = mod.table_s1(sg["xlsx"])
    assert set(s1) == NAMES
    for n, (c, bs, wgs, q) in TABLE.items():
        assert (s1[n]["classification"], s1[n]["biosample"], s1[n]["wgs"], s1[n]["quality"]) == (c, bs, wgs, q)
    # a name the GTDB-Tk sheet gives twice has no classification; the Sequencing sheet's full names have none either
    assert (s1["AB-754_C22"]["classification"], s1["AB-754_C22"]["gtdbtk_rows"]) == ("", 2)
    assert (s1["AB-754_C22_AB-905"]["classification"], s1["AB-754_C22_AB-905"]["biosample"],
            s1["AB-754_C22_AB-905"]["gtdbtk_rows"]) == ("", "SAMN32251362", 0)
    assert mod.wgs_project("AFHZ00000000") == "AFHZ" and mod.wgs_project("CABWYC000000000.1") == "CABWYC"
    assert mod.wgs_project("AQSI01000000") == "AQSI" and mod.wgs_project("none") == ""


def test_a_sequencing_sheet_whose_columns_moved_is_refused(sg, tmp_path):
    mod = _mod()
    global SEQ_HEAD
    keep = SEQ_HEAD
    try:
        SEQ_HEAD = keep[:16] + ["Something else"] + keep[17:]
        x = tmp_path / "moved.xlsx"
        write_xlsx(x)
    finally:
        SEQ_HEAD = keep
    with pytest.raises(mod.Refused, match="column 16"):
        mod.table_s1(x)


def test_members_match_the_longest_table_name_and_a_shared_name_matches_neither():
    mod = _mod()
    m = mod.match_members(sorted(SAG), NAMES)
    assert m == {"AB-754_C22_AB-905_clean": ("AB-754_C22_AB-905", "prefix"),
                 "AB-754_C22_AB-908_clean": ("AB-754_C22_AB-908", "prefix"),
                 "AAA001_A19_output_clean": ("AAA001_A19", "prefix"),
                 "AAA001_C10_contigs": ("AAA001_C10", "prefix"),
                 "JGI_01_D11": ("JGI_01_D11", "exact"), "JGI_01_D11_4000m": ("JGI_01_D11_4000m", "exact"),
                 "X_01_A1": (None, "tie"), "X_01_A1_rep": (None, "tie"), "Unknown_9": (None, "none"),
                 "AG-311-D23": ("AG-311-D23", "exact")}


def test_census_names_every_entry_and_refuses_nothing_silently(sg, tmp_path):
    mod = _mod()
    t = tmp_path / "odd.tar.gz"
    write_tar(t, {"JGI_01_D11": SAG["JGI_01_D11"]}, extra={"README.txt": b"hi", "Final_Catalog/sub/x.fna": b">a\nA\n"})
    mem, census = mod.sag3_members(t)
    assert [s for _, s, _ in mem] == ["JGI_01_D11"]
    assert census["cells"] == 1 and census["directories"] == 1
    assert census["other"] == ["Final_Catalog/sub/x.fna", "README.txt"]


def test_sag3_stage_records_r232_status_screens_codes_and_is_deterministic(sg):
    mod = _mod()
    rep1 = run_sag3(mod, sg, "a")
    rep2 = run_sag3(mod, sg, "b")
    t = sg["tmp"]
    meta, fa = t / "meta_a" / "SAG3", t / "fasta_a" / "SAG3" / "sag3_omz.fna.gz"

    assert not (meta / "removed.tsv").exists()
    st = {r["genome"]: r for r in _rows(meta / "r232_status.tsv")}
    assert {g: r["r232_status"] for g, r in st.items()} == {
        "AAA001_A19": "not_visible", "AAA001_C10": "in_gtdb_and_globdb_r232", "AG-311-D23": "not_in_gtdb_r232",
        "JGI_01_D11": "in_gtdb_r232", "JGI_01_D11_4000m": "not_in_gtdb_r232", "Unknown_9": "not_visible",
        "X_01_A1": "not_visible", "X_01_A1_rep": "not_visible", "AB-754_C22_AB-905": "not_in_gtdb_r232",
        "AB-754_C22_AB-908": "not_in_gtdb_r232"}
    c10, d11 = st["AAA001_C10"], st["JGI_01_D11"]
    assert (c10["gtdb_r232_accession"], c10["via"], c10["is_gtdb_representative"], c10["globdb_r232_id"]) == (
        "GCA_000000101", "wgs", "True", "GCA_000000101")
    assert (d11["gtdb_r232_accession"], d11["via"], d11["in_globdb_r232"], d11["gtdb_species"]) == (
        "GCF_000000102", "biosample", "no", "Actinomarina sp2")

    kept = sorted(GENOME.values())
    genomes = json.loads((meta / "genomes.json").read_text())
    stem = {g: s for s, g in GENOME.items()}
    assert [g["contig_uids"] for g in genomes] == [[cuid(q) for _, q in SAG[stem[g]]] for g in kept]
    assert {g["kind"] for g in genomes} == {"prokaryote"}
    codes = dict(zip(kept, (g["genetic_code"] for g in genomes)))
    assert codes["AG-311-D23"] == 4 and {c for g, c in codes.items() if g != "AG-311-D23"} == {11}

    # the FASTA follows the archive (reverse name order here); a name not starting with the genome's gets its name
    names = [ln[1:] for ln in _gz_lines(fa) if ln.startswith(">")]
    want = []
    for s in sorted(SAG, reverse=True):
        g = GENOME[s]
        want += [h if h.startswith(g + "_") else f"{g}|{h}" for h, _ in SAG[s]]
    assert names == want
    assert "JGI_01_D11|JGI01_D11_NODE_1" in names and "AG-311-D23|2634166711_Ga0006891_101" in names

    mt = {r["genome"]: r for r in _rows(meta / "sag_match.tsv")}
    assert (mt["X_01_A1"]["basis"], mt["AAA001_A19"]["basis"], mt["Unknown_9"]["basis"]) == ("tie", "prefix", "none")
    assert [(r["genome"], r["owner_stage"]) for r in _rows(meta / "shared_records.tsv")] == [("JGI_01_D11", "V2")]
    assert [(r["genome"], r["holder"]) for r in _rows(meta / "later_records.tsv")] == [("JGI_01_D11_4000m", "SAG1")]
    objs = json.loads((meta / "objects.json").read_text())
    assert objs[0]["path"] == "/tmp/if2453/s2/inputs/SAG3/sag3_omz.fna.gz" and objs[0]["redistributable"] is True

    r = rep1
    # unmatched names: NOT_IN_ARCHIVE, X_01_A1 (tie), AB-754_C22 (the GTDB-Tk sheet's short name)
    assert r["cells"] == 10 and r["genomes"] == 10 and r["table_s1_unmatched"] == 3
    assert r["contigs_renamed"] == 3 and r["without_classification"] == 5  # two AB cells, two tied, Unknown_9
    assert r["match_basis"] == {"exact": 3, "none": 1, "prefix": 4, "tie": 2}
    assert r["via"] == {"biosample": 1, "wgs": 1}
    assert r["r232_status"] == {"in_gtdb_and_globdb_r232": 1, "in_gtdb_r232": 1, "not_in_gtdb_r232": 4,
                                "not_visible": 4}
    assert r["held_by"] == {"V2": 1} and r["later_held_by"] == {"SAG1": 1}
    assert r["bases"] == sum(len(q) for recs in SAG.values() for _, q in recs)

    for f in ("genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
              "later_records.tsv", "genetic_codes.tsv", "sag_match.tsv", "r232_status.tsv"):
        assert (meta / f).read_bytes() == (t / "meta_b" / "SAG3" / f).read_bytes(), f
    assert fa.read_bytes() == (t / "fasta_b" / "SAG3" / "sag3_omz.fna.gz").read_bytes()
    assert {k: v for k, v in rep1["sha256"].items() if k != "objects_staging.json"} == {
        k: v for k, v in rep2["sha256"].items() if k != "objects_staging.json"}

    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), "SAG3", str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


def test_sag3_lists_a_repeated_contig_set_once(sg, tmp_path):
    mod = _mod()
    same = [("JGI_01_D11_NODE_1", seq_of(301, 2000))]
    t = tmp_path / "dup.tar.gz"
    write_tar(t, {"JGI_01_D11": same, "Unknown_9": [("JGI_01_D11_NODE_1x", seq_of(301, 2000))]})
    rep = run_sag3(mod, {**sg, "tar": t}, "d")
    meta = sg["tmp"] / "meta_d" / "SAG3"
    assert len(json.loads((meta / "genomes.json").read_text())) == 1 and rep["genomes_repeated"] == 1
    assert [(r["genome"], r["first_genome"]) for r in _rows(meta / "repeated_records.tsv")] == [
        ("Unknown_9", "JGI_01_D11")]


def test_sag3_refuses_a_name_twice_an_empty_archive_and_a_bad_index_and_writes_nothing(sg, tmp_path):
    mod = _mod()
    t = tmp_path / "twice.tar.gz"
    write_tar(t, {"JGI_01_D11": [("JGI_01_D11_NODE_1", seq_of(1, 900)), ("JGI_01_D11_NODE_1", seq_of(2, 900))]})
    with pytest.raises(mod.Refused, match="occurs twice"):
        run_sag3(mod, {**sg, "tar": t}, "t")
    assert not (sg["tmp"] / "meta_t" / "SAG3").exists()
    d = sg["tmp"] / "fasta_t" / "SAG3"
    assert not d.exists() or list(d.iterdir()) == []
    e = tmp_path / "empty.tar.gz"
    write_tar(e, {})
    with pytest.raises(mod.Refused, match="no single-cell member"):
        run_sag3(mod, {**sg, "tar": e}, "e")
    with pytest.raises(mod.Refused, match="C-74's index"):
        run_sag3(mod, {**sg, "index": tmp_path / "nope"}, "i")
    run_sag3(mod, sg, "a")
    with pytest.raises(mod.Refused, match="exists"):
        run_sag3(mod, sg, "a")


def test_sag3_best_hits_keep_the_highest_ani_at_the_thresholds(tmp_path):
    mod = _mod()
    head = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query", "Ref_name", "Query_name"]
    rows = [["r/SAG1/a.fna", "q/g1.fna", "97.0", "60.0", "20.0", "x", "y"],
            ["r/SAG2/b.fna", "q/g1.fna", "98.5", "10.0", "55.0", "x", "y"],
            ["r/SAG2/c.fna", "q/g2.fna", "95.9", "99.0", "99.0", "x", "y"],
            ["r/SAG1/d.fna", "q/g3.fna", "99.0", "40.0", "49.9", "x", "y"]]
    tab = tmp_path / "t.tsv"
    tab.write_text("\n".join("\t".join(r) for r in [head, *rows]) + "\n")
    h = mod.best_hits(tab, {"q/g1.fna": "g1", "q/g2.fna": "g2", "q/g3.fna": "g3"},
                      {"r/SAG1/a.fna": ("SAG1", "a"), "r/SAG2/b.fna": ("SAG2", "b"), "r/SAG2/c.fna": ("SAG2", "c"),
                       "r/SAG1/d.fna": ("SAG1", "d")})
    assert h == {"g1": {"ani": 98.5, "af": 55.0, "stage": "SAG2", "genome": "b"}}


def test_archaea_refs_take_the_r226_archaea_of_the_asked_genera_only(tmp_path):
    mod = _mod()
    assert mod.genus_of("d__Archaea;p__Thermoproteota;g__Nitrosopelagicus;s__") == "Nitrosopelagicus"
    assert mod.genus_of("Unclassified Bacteria") == ""
    fna = tmp_path / "arc.fna.gz"
    fna.write_bytes(gzip.compress(_fa([("A1_c1", "ACGT" * 5), ("A1_c2", "TTTT" * 5), ("B2_c1", "GGGG" * 5),
                                       ("C3_c1", "CCCC" * 5)]).encode(), mtime=0))
    lst = tmp_path / "list.tsv"
    lst.write_text("genome\tgtdb\tfirst_record\trecords\tbp\theaders_start_with_genome\ttaxonomy\n"
                   "A1\ttrue\t0\t2\t40\ttrue\td__Archaea;p__Thermoproteota;g__Nitrosopelagicus;s__x\n"
                   "B2\tfalse\t2\t1\t20\ttrue\td__Archaea;p__Thermoproteota;g__Nitrosopumilus;s__y\n"
                   "C3\tfalse\t3\t1\t20\ttrue\td__Archaea;p__Halobacteriota;g__Nitrosopelagicus_A;s__z\n")
    files = mod.archaea_refs(fna, lst, {"Nitrosopelagicus"}, tmp_path / "out")
    assert sorted(files) == ["A1"]
    assert [ln for ln in files["A1"].read_text().splitlines() if ln.startswith(">")] == [">A1_c1", ">A1_c2"]


def test_split_genomes_writes_one_file_per_genome(sg, tmp_path):
    mod = _mod()
    run_sag3(mod, sg, "p")
    meta, fa = sg["tmp"] / "meta_p" / "SAG3", sg["tmp"] / "fasta_p" / "SAG3" / "sag3_omz.fna.gz"
    files = mod.split_genomes(meta, fa, tmp_path / "split")
    assert sorted(files) == sorted(GENOME.values())
    got = [ln for ln in files["JGI_01_D11"].read_text().splitlines() if not ln.startswith(">")]
    assert got == [q for _, q in SAG["JGI_01_D11"]]


# ==================================================================================================================
# WET1
# ==================================================================================================================

# The archive holds MAG00001..3 (MAG00004 is in the table, dereplicated away); every contig name is a sample's
# assembly name, which the stager prefixes with the MAG id. MAG00002 is in GTDB R232 under PRJNA1430164 by its bin name;
# MAG00003 is called Mycoplasma by GTDB-Tk (code 4); one MAG00001 contig is held by V2, one MAG00002 contig by SAG3.
WET = {
    "MAG00001": [("BYD1_k127_1", seq_of(501, 2200)), ("BYD1_k127_2", seq_of(502, 1300))],
    "MAG00002": [("BYD1_k127_3", seq_of(503, 2400)), ("BYD1_k127_4", seq_of(504, 900))],
    "MAG00003": [("ZJ2_k127_1", seq_of(505, 2500))],
}
WET_META = [
    ("MAG00001", "BYD1_bin.10", "SMP0001", "d__Bacteria;p__Acidobacteriota;g__QHXN01;s__", "95.97", "0.93", "High-quality"),
    ("MAG00002", "BYD1_bin.11", "SMP0001", "d__Archaea;p__Thermoproteota;g__Nitrosotalea;s__", "91.2", "8.08", "Medium-quality"),
    ("MAG00003", "ZJ2_bin.1", "SMP0002", "d__Bacteria;p__Bacillota;g__Mycoplasma;s__", "61.8", "1.59", "Medium-quality"),
    ("MAG00004", "ZJ2_bin.2", "SMP0002", "d__Bacteria;p__Bacillota;g__Bacillus;s__", "71.0", "1.0", "Medium-quality"),
]
WET_HEAD = ["ID", "Name", "Sample ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
            "GTDB_classification", "Completeness (%)", "Contamination (%)", "GC (%)", "N50 (bp)", " Genome size (bp)",
            "Quality score", "Num_5S_rRNA", "Num_16S_rRNA", "Num_23S_rRNA", "Num_tRNA", "Quality_tier (MIMAG)"]


def write_wet_xlsx(path: Path) -> None:
    rows = [["Table S2. The detailed information of all MAGs used in the study"], WET_HEAD]
    for i, n, s, c, comp, cont, q in WET_META:
        rows.append([i, n, s, "Bacteria", "", "", "", "", "", "", c, comp, cont, "50", "1000", "2000", "80", "0", "0",
                     "0", "10", q])
    wb = ('<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/'
          '2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
          '<sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>')
    rels = ('<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/'
            '2006/relationships"><Relationship Id="rId1" Type="x" Target="worksheets/sheet1.xml"/></Relationships>')
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", rels)
        z.writestr("xl/worksheets/sheet1.xml", _sheet_xml(rows, {}))


def write_wet_zip(path: Path, members=WET, *, extra: dict[str, bytes] | None = None) -> None:
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("dereplicated_genomes/", b"")
        for mag in sorted(members, reverse=True):  # archive order is not name order
            z.writestr(f"dereplicated_genomes/{mag}.fa", _fa(members[mag]))
        for n, d in (extra or {}).items():
            z.writestr(n, d)


def write_wet_gtdb(d: Path) -> None:
    d.mkdir()
    head = ["accession", "ncbi_bioproject", "ncbi_isolate", "ncbi_strain_identifiers", "ncbi_biosample",
            "ncbi_wgs_master"]
    ar = [["GB_GCA_000000201.1", "PRJNA1430164", "BYD1_bin.11", "none", "SAMN1", "JBAA00000000.1"]]
    bac = [["GB_GCA_000000202.1", "PRJNA1430164", "none", "ZJ9_bin.9", "SAMN2", "JBAB00000000.1"],
           ["GB_GCA_000000203.1", "PRJNA999", "BYD1_bin.10", "none", "SAMN3", "JBAC00000000.1"]]
    for name, rows in (("ar53_metadata_r232.tsv.gz", ar), ("bac120_metadata_r232.tsv.gz", bac)):
        (d / name).write_bytes(gzip.compress(("\n".join("\t".join(r) for r in [head, *rows]) + "\n").encode(),
                                             mtime=0))


def write_wet_index(d: Path) -> None:
    d.mkdir()
    gt = [("GCA_000000201", "GCA_000000201", True, "GCA_000000201", "Nitrosotalea sp1"),
          ("GCA_000000202", "GCA_000000202", False, "GCA_000000201", "Nitrosotalea sp1"),
          ("GCA_000000203", "GCA_000000203", True, "GCA_000000203", "Other sp")]
    pq.write_table(pa.table({
        "accession": [g[0] for g in gt], "genbank_accession": [g[1] for g in gt],
        "is_representative": [g[2] for g in gt], "representative": [g[3] for g in gt],
        "species": [g[4] for g in gt]}), d / "gtdb_r232.parquet")
    pq.write_table(pa.table({
        "globdb_id": ["GCA_000000201"], "insdc_accession": ["GCA_000000201"], "original_id": [None],
        "cluster_id": [None]}), d / "globdb_r232.parquet")


@pytest.fixture()
def wt(tmp_path: Path) -> dict:
    z = tmp_path / "Dereplicated_genomes.zip"
    write_wet_zip(z)
    x = tmp_path / "MAG_meta.xlsx"
    write_wet_xlsx(x)
    g = tmp_path / "gtdb"
    write_wet_gtdb(g)
    i = tmp_path / "index"
    write_wet_index(i)
    held = tmp_path / "chain" / "V2"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(WET["MAG00001"][1][1])]}]))
    later = tmp_path / "later" / "SAG3"
    later.mkdir(parents=True)
    (later / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(WET["MAG00002"][0][1])]}]))
    return {"zip": z, "xlsx": x, "gtdb": g, "index": i, "held": held, "later": later, "tmp": tmp_path}


def run_wet1(mod, w: dict, tag: str, **kw):
    t = w["tmp"]
    return mod.stage_wet1(
        w["zip"], w["xlsx"], w["index"], w["gtdb"], t / f"meta_{tag}", t / f"fasta_{tag}", root=None, root_stages=(),
        owner_inputs={"V2": w["held"]}, later_inputs={"SAG3": w["later"]}, inputs_home=Path("/tmp/if2453/s2/inputs"),
        work=t / f"work_{tag}", **kw)


def test_mag_meta_is_read_under_its_title_row_and_the_zip_census_names_every_entry(wt, tmp_path):
    mod = _mod()
    m = mod.mag_meta(wt["xlsx"])
    assert sorted(m) == ["MAG00001", "MAG00002", "MAG00003", "MAG00004"]
    assert (m["MAG00002"]["name"], m["MAG00002"]["sample"], m["MAG00002"]["classification"]) == (
        "BYD1_bin.11", "SMP0001", "d__Archaea;p__Thermoproteota;g__Nitrosotalea;s__")
    mem, census = mod.wet1_members(wt["zip"])
    assert [g for _, g in mem] == ["MAG00003", "MAG00002", "MAG00001"]  # archive order
    assert census == {"entries": 4, "directories": 1, "mags": 3, "other": []}
    odd = tmp_path / "odd.zip"
    write_wet_zip(odd, {"MAG00001": WET["MAG00001"]}, extra={"README.txt": b"x", "dereplicated_genomes/a/b.fa": b">a\nA\n"})
    assert mod.wet1_members(odd)[1]["other"] == ["README.txt", "dereplicated_genomes/a/b.fa"]


def test_project_crosswalk_matches_gtdb_rows_of_the_bioproject_by_bin_name(wt):
    mod = _mod()
    by_name, n_rows = mod.gtdb_project_crosswalk(wt["gtdb"], "PRJNA1430164", {"BYD1_bin.10", "BYD1_bin.11",
                                                                             "ZJ9_bin.9"})
    # BYD1_bin.10 is in R232 only under another BioProject, so it is not this MAG
    assert by_name == {"BYD1_bin.11": ["GCA_000000201"], "ZJ9_bin.9": ["GCA_000000202"]} and n_rows == 2


def test_wet1_stage_records_r232_status_prefixes_names_and_is_deterministic(wt):
    mod = _mod()
    rep1 = run_wet1(mod, wt, "a")
    rep2 = run_wet1(mod, wt, "b")
    t = wt["tmp"]
    meta, fa = t / "meta_a" / "WET1", t / "fasta_a" / "WET1" / "wet1_gwm.fna.gz"
    assert not (meta / "removed.tsv").exists()
    st = {r["genome"]: r for r in _rows(meta / "r232_status.tsv")}
    assert {g: r["r232_status"] for g, r in st.items()} == {
        "MAG00001": "not_in_gtdb_r232", "MAG00002": "in_gtdb_and_globdb_r232", "MAG00003": "not_in_gtdb_r232"}
    assert (st["MAG00002"]["gtdb_r232_accession"], st["MAG00002"]["name"]) == ("GCA_000000201", "BYD1_bin.11")
    genomes = json.loads((meta / "genomes.json").read_text())
    assert [g["contig_uids"] for g in genomes] == [[cuid(q) for _, q in WET[m]] for m in sorted(WET)]
    assert [g["genetic_code"] for g in genomes] == [11, 11, 4]
    assert {g["kind"] for g in genomes} == {"prokaryote"}
    assert [g["assembly_quality"] for g in genomes] == ["mag_hq", "mag_mq", "mag_mq"]  # MAG_meta's MIMAG tier
    names = [ln[1:] for ln in _gz_lines(fa) if ln.startswith(">")]
    assert names == [f"{m}|{h}" for m in sorted(WET, reverse=True) for h, _ in WET[m]]
    assert [(r["genome"], r["owner_stage"]) for r in _rows(meta / "shared_records.tsv")] == [("MAG00001", "V2")]
    assert [(r["genome"], r["holder"]) for r in _rows(meta / "later_records.tsv")] == [("MAG00002", "SAG3")]
    objs = json.loads((meta / "objects.json").read_text())
    assert objs[0]["path"] == "/tmp/if2453/s2/inputs/WET1/wet1_gwm.fna.gz" and objs[0]["redistributable"] is True
    r = rep1
    assert (r["mags"], r["genomes"], r["table_rows"], r["table_not_in_archive"]) == (3, 3, 4, 1)
    assert r["r232_status"] == {"in_gtdb_and_globdb_r232": 1, "not_in_gtdb_r232": 2}
    assert r["gtdb_rows_of_bioproject"] == 2 and r["held_by"] == {"V2": 1} and r["later_held_by"] == {"SAG3": 1}
    for f in ("genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
              "later_records.tsv", "genetic_codes.tsv", "r232_status.tsv"):
        assert (meta / f).read_bytes() == (t / "meta_b" / "WET1" / f).read_bytes(), f
    assert fa.read_bytes() == (t / "fasta_b" / "WET1" / "wet1_gwm.fna.gz").read_bytes()
    assert {k: v for k, v in rep1["sha256"].items() if k != "objects_staging.json"} == {
        k: v for k, v in rep2["sha256"].items() if k != "objects_staging.json"}
    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), "WET1", str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


def test_wet1_refuses_a_quality_tier_it_cannot_name():
    mod = _mod()
    assert (mod.mimag_quality("High-quality"), mod.mimag_quality("Medium-quality")) == ("mag_hq", "mag_mq")
    with pytest.raises(mod.Refused, match="Low-quality"):
        mod.mimag_quality("Low-quality")


def test_wet1_refuses_a_member_the_table_lacks_and_writes_nothing(wt, tmp_path):
    mod = _mod()
    z = tmp_path / "x.zip"
    write_wet_zip(z, {**WET, "MAG09999": [("Q_k127_1", seq_of(9, 900))]})
    with pytest.raises(mod.Refused, match="MAG09999"):
        run_wet1(mod, {**wt, "zip": z}, "x")
    assert not (wt["tmp"] / "meta_x" / "WET1").exists()
    d = wt["tmp"] / "fasta_x" / "WET1"
    assert not d.exists() or list(d.iterdir()) == []


def test_genus_genomes_and_split_genomes_keep_only_what_is_asked(wt, tmp_path):
    mod = _mod()
    tax = tmp_path / "non_gtdb_stages.tsv"
    tax.write_text("ID\tstage\ttaxonomy\nG_1\tGART1\td__Bacteria;p__Acidobacteriota;g__QHXN01;s__\n"
                   "G_2\tGART1\td__Bacteria;p__X;g__Other;s__\n")
    assert mod.genus_genomes(tax, {"QHXN01", "Nitrosotalea"}) == {"G_1"}
    run_wet1(mod, wt, "k")
    meta, fa = wt["tmp"] / "meta_k" / "WET1", wt["tmp"] / "fasta_k" / "WET1" / "wet1_gwm.fna.gz"
    files = mod.split_genomes(meta, [fa], tmp_path / "split", keep={"MAG00002"})
    assert sorted(files) == ["MAG00002"]


def test_split_genomes_reads_a_records_table_whose_column_is_genome_id(tmp_path):
    mod = _mod()
    meta = tmp_path / "NX1"
    meta.mkdir()
    a, b = seq_of(601, 900), seq_of(602, 800)
    (meta / "records.tsv").write_text(f"genome_id\tsource_name\trecord_ordinal\tcontig_uid\nX:1\ts\t0\t{cuid(a)}\n"
                                      f"X:2\ts\t0\t{cuid(b)}\n")
    fa = tmp_path / "nx1.fna"
    fa.write_text(_fa([("c1", a), ("c2", b)]))
    files = mod.split_genomes(meta, fa, tmp_path / "out")
    assert sorted(files) == ["X:1", "X:2"] and files["X:1"].name == "X_1.fna"


def test_split_genomes_writes_bytes_as_they_are_even_when_they_are_not_utf8(tmp_path):
    mod = _mod()
    meta = tmp_path / "S"
    meta.mkdir()
    a = b"ACGT\xe5ACGT"
    (meta / "records.tsv").write_text(
        f"genome\tcontig_uid\nG1\t{blake3(a.upper()).digest(length=16).hex()}\n")
    fa = tmp_path / "s.fna"
    fa.write_bytes(b">c1\n" + a + b"\n")
    files = mod.split_genomes(meta, fa, tmp_path / "out")
    assert files["G1"].read_bytes() == b">c1\n" + a + b"\n"


def test_split_genomes_refuses_an_unlisted_record_unless_asked_to_count_it(tmp_path):
    mod = _mod()
    meta = tmp_path / "S"
    meta.mkdir()
    a, b = seq_of(701, 900), seq_of(702, 800)
    (meta / "records.tsv").write_text(f"genome\tcontig_uid\nG1\t{cuid(a)}\n")
    fa = tmp_path / "s.fna"
    fa.write_text(_fa([("c1", a), ("c2", b)]))
    with pytest.raises(mod.Refused, match="in no genome"):
        mod.split_genomes(meta, fa, tmp_path / "o1")
    unlisted = []
    files = mod.split_genomes(meta, fa, tmp_path / "o2", unlisted=unlisted)
    assert sorted(files) == ["G1"] and unlisted == ["c2"]


def test_merge_best_keeps_each_genomes_best_hit_over_groups():
    mod = _mod()
    a = {"g1": {"ani": 97.0, "af": 60.0, "stage": "LR2", "genome": "x"},
         "g2": {"ani": 99.0, "af": 70.0, "stage": "LR2", "genome": "y"}}
    b = {"g1": {"ani": 98.0, "af": 55.0, "stage": "GART1", "genome": "z"}}
    assert mod.merge_best([a, b]) == {"g1": b["g1"], "g2": a["g2"]}


# ==================================================================================================================
# the conf templates
# ==================================================================================================================

@pytest.mark.parametrize("stage,fasta_name", [("GWV1", "gwv1_gwvc.fna.gz"), ("SAG3", "sag3_omz.fna.gz"),
                                              ("WET1", "wet1_gwm.fna.gz")])
def test_templates_pass_stage_est_and_carry_the_cards_values(stage, fasta_name):
    est = _load("stage_est", "stage_est.py")
    c = est.read_conf(TEMPLATES / f"{stage}.conf", stage, tree=REPO)
    assert (c.depth, c.crispr, c.library) == ("called", True, None)
    assert c.after == {"GWV1": ("INPH2604",), "SAG3": ("GWV1",), "WET1": ("SAG3",)}[stage]
    assert str(c.inputs) == f"/tmp/if2453/s2/inputs/{stage}" and c.min_contig_len is None
    assert c.genecall_mode == "default"  # kind="virus" is called by pyrodigal-gv in any case (genecall.mode_of)
    assert NOTE in c.notes
    assert {n for _, n in c.sums} == {"genomes.json", "objects.json", fasta_name}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    for h, _ in c.sums:
        assert h in note, f"SUM {h} not in the note"
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    allowed = ROOT13 | {s for s, _ in a26.R26_11_4_CHAIN} | {"MV2"}
    assert set(c.expect) <= allowed
    for owner, (n, _) in c.expect.items():
        assert n > 0

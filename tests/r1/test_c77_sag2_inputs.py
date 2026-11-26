"""Card C-77 (DEC-79): stage the GORG-Dark single cells (`SAG2`) from the figshare zip.

Oracles are written here separately from the script: contig_uid by blake3 directly, membership of GTDB R232 and
GlobDB r232 by hand-written fixture tables, the matching of single cells to the ENA table by hand. Each test runs
twice, with the normal PATH and with PATH reduced to `/usr/bin:/bin` (as tests/r1/test_r232_index.py does): the
stager needs only Python and its libraries, and GitHub Actions has none of the pod's tools.
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
import zipfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-77")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-26_sag2_inputs.md"


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
    return _load("c77_sag2_inputs", "c77_sag2_inputs.py")


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


def _fa_gz(records: list[tuple[str, str]], width: int = 60) -> bytes:
    out = io.StringIO()
    for hdr, s in records:
        out.write(f">{hdr}\n")
        for i in range(0, len(s), width):
            out.write(s[i:i + width] + "\n")
    return gzip.compress(out.getvalue().encode(), mtime=0)


# The fixture world. Seven single cells in the zip (plus the macOS resource forks and the directory entry a Mac zip
# carries); the ENA table lists six assemblies. The real zip's members carry no accession and no assembly name, so
# the only link to the ENA table is the total of bases (the table's base_count).
#   SCGC_AM-100-A01  4,500 bp, one ENA row of that size; GCA_900000001 is in GTDB R232 -> removed
#   SCGC_AM-100-A02  2,999 bp, one ENA row; GCA_900000002 is in GlobDB r232 only -> removed
#   SCGC_AM-100-A03  4,000 bp, one ENA row; in neither -> kept
#   SCGC_AM-100-A04  2,500 bp, one ENA row; in neither -> kept; the row's name is "Mycoplasma" -> code 4
#   SCGC_AM-100-A05  1,200 bp, as A06 and two ENA rows of 1,200 bp: no telling which is which -> kept, not seen
#   SCGC_AM-100-A06  1,200 bp
#   SCGC_AM-100-A07  2,501 bp, no ENA row -> kept, not seen; contig names that do not carry the cell's id; one
#                    contig is held by an earlier stage
SAGS: dict[str, list[tuple[str, str]]] = {
    "SCGC_AM-100-A01": [("AM-100-A01_NODE_1", seq_of(1, 3000)), ("AM-100-A01_NODE_2", seq_of(2, 1500))],
    "SCGC_AM-100-A02": [("AM-100-A02_NODE_1", seq_of(3, 2000)), ("AM-100-A02_NODE_2", seq_of(4, 999))],
    "SCGC_AM-100-A03": [("AM-100-A03_NODE_1", seq_of(5, 4000))],
    "SCGC_AM-100-A04": [("AM-100-A04_NODE_1", seq_of(6, 2500))],
    "SCGC_AM-100-A05": [("AM-100-A05_NODE_1", seq_of(7, 1200))],
    "SCGC_AM-100-A06": [("AM-100-A06_NODE_1", seq_of(8, 1200))],
    "SCGC_AM-100-A07": [("contig1", seq_of(9, 1800)), ("contig2", seq_of(10, 701))],
}
MEMBER = {s: f"GORG-Dark SAGs/{s}_contigs.fasta.gz" for s in SAGS}

ENA_HEAD = ["accession", "assembly_name", "base_count", "scientific_name", "tax_id", "wgs_set", "version"]
ENA = [
    ("GCA_900000001", "ASM90000001v1", "4500", "Alphaproteobacteria bacterium", "1", "JBAAAA01", "1"),
    ("GCA_900000002", "ASM90000002v1", "2999", "Gammaproteobacteria bacterium", "2", "JBAAAB01", "1"),
    ("GCA_900000003", "ASM90000003v1", "4000", "Candidatus Pelagibacter sp.", "3", "JBAAAC01", "1"),
    ("GCA_900000004", "ASM90000004v1", "2500", "Mycoplasma sp.", "4", "JBAAAD01", "1"),
    ("GCA_900000005", "ASM90000005v1", "1200", "Bacteria bacterium", "5", "JBAAAE01", "1"),
    ("GCA_900000006", "ASM90000006v1", "1200", "Bacteria bacterium", "6", "JBAAAF01", "1"),
]


def write_zip(path: Path, members: dict[str, bytes], *, junk: bool = True) -> None:
    with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as z:
        z.writestr("GORG-Dark SAGs/", b"")
        for name, data in members.items():
            z.writestr(name, data)
            if junk:
                z.writestr("__MACOSX/" + name.replace("GORG-Dark SAGs/", "GORG-Dark SAGs/._"), b"\x00" * 441)


def write_ena(path: Path, rows=ENA) -> None:
    path.write_text("\n".join("\t".join(r) for r in [ENA_HEAD, *rows]) + "\n")


def write_index(d: Path) -> None:
    """The columns of C-74's index that `lookup` reads, by hand: GTDB R232 holds A01's accession (and a second
    genome that nothing here has); GlobDB r232 holds A02's accession as a resource's INSDC accession."""
    d.mkdir()
    gt = [("GCA_900000001", "GCA_900000001", True, "GCA_900000001", "Species one"),
          ("GCA_900000099", "GCA_900000099", True, "GCA_900000099", "Species ninety nine")]
    pq.write_table(pa.table({
        "accession": [g[0] for g in gt], "genbank_accession": [g[1] for g in gt],
        "is_representative": [g[2] for g in gt], "representative": [g[3] for g in gt],
        "species": [g[4] for g in gt]}), d / "gtdb_r232.parquet")
    pq.write_table(pa.table({
        "globdb_id": ["MOTU40_000001", "GEM_000002"], "insdc_accession": [None, "GCA_900000002"],
        "original_id": [None, None], "cluster_id": [None, None]}), d / "globdb_r232.parquet")


@pytest.fixture()
def world(tmp_path: Path) -> dict:
    zp = tmp_path / "GORG-Dark_SAGs.zip"
    write_zip(zp, {MEMBER[s]: _fa_gz(recs) for s, recs in SAGS.items()})
    ena = tmp_path / "ena.tsv"
    write_ena(ena)
    idx = tmp_path / "index"
    write_index(idx)
    held = tmp_path / "chain" / "V2"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(SAGS["SCGC_AM-100-A07"][1][1])]}]))
    later = tmp_path / "later" / "SAG1"
    later.mkdir(parents=True)
    (later / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(SAGS["SCGC_AM-100-A03"][0][1])]}]))
    return {"zip": zp, "ena": ena, "index": idx, "held": held, "later": later, "tmp": tmp_path}


def run_stage(mod, w: dict, tag: str, **kw):
    t = w["tmp"]
    return mod.stage(
        w["zip"], w["ena"], w["index"], t / f"meta_{tag}", t / f"fasta_{tag}",
        root=None, root_stages=(), owner_inputs={"V2": w["held"]}, later_inputs={"SAG1": w["later"]},
        inputs_home=Path("/tmp/if2453/s2/inputs"), retrieved_at="2025-05-15T00:00:00Z", work=t / f"work_{tag}", **kw)


# --- the zip -------------------------------------------------------------------------------------------------------

def test_members_are_the_single_cells_and_not_the_macos_forks(world):
    mod = _mod()
    ms = mod.members(world["zip"])
    assert [m.sag for m in ms] == sorted(SAGS)
    assert {m.sag: m.name for m in ms} == MEMBER
    census = mod.census(world["zip"])
    assert census == {"entries": 1 + 2 * len(SAGS), "directories": 1, "macosx": len(SAGS), "sags": len(SAGS),
                      "other": []}


def test_a_zip_with_a_member_that_is_not_a_fasta_is_reported_not_hidden(world, tmp_path):
    mod = _mod()
    z = tmp_path / "odd.zip"
    write_zip(z, {MEMBER["SCGC_AM-100-A01"]: _fa_gz(SAGS["SCGC_AM-100-A01"]), "GORG-Dark SAGs/README.txt": b"hi"})
    assert mod.census(z)["other"] == ["GORG-Dark SAGs/README.txt"]
    assert [m.sag for m in mod.members(z)] == ["SCGC_AM-100-A01"]


# --- matching to the ENA table -------------------------------------------------------------------------------------

def test_match_by_total_bases_when_one_cell_and_one_row_have_them_and_never_by_a_tie(world):
    mod = _mod()
    scan = mod.scan(world["zip"])
    assert {s: (d["bp"], d["contigs"]) for s, d in scan.items()} == {
        s: (sum(len(q) for _, q in recs), len(recs)) for s, recs in SAGS.items()}
    ena = mod.read_ena(world["ena"])
    m = mod.match(scan, ena)
    want = {
        "SCGC_AM-100-A01": ("GCA_900000001", "base_count"),
        "SCGC_AM-100-A02": ("GCA_900000002", "base_count"),
        "SCGC_AM-100-A03": ("GCA_900000003", "base_count"),
        "SCGC_AM-100-A04": ("GCA_900000004", "base_count"),
        "SCGC_AM-100-A05": (None, "ambiguous_base_count"),
        "SCGC_AM-100-A06": (None, "ambiguous_base_count"),
        "SCGC_AM-100-A07": (None, "no_ena_row"),
    }
    assert {s: (r["accession"], r["basis"]) for s, r in m.items()} == want
    got = [r["accession"] for r in m.values() if r["accession"]]
    assert len(got) == len(set(got))  # two cells never take the same ENA row


def test_two_cells_of_the_same_size_and_one_row_match_nothing(world, tmp_path):
    mod = _mod()
    z = tmp_path / "tie.zip"
    write_zip(z, {"GORG-Dark SAGs/SCGC_AM-100-B01_contigs.fasta.gz": _fa_gz([("b1", seq_of(21, 4000))]),
                  "GORG-Dark SAGs/SCGC_AM-100-B02_contigs.fasta.gz": _fa_gz([("b2", seq_of(22, 4000))])})
    m = mod.match(mod.scan(z), mod.read_ena(world["ena"]))
    assert {s: (r["accession"], r["basis"]) for s, r in m.items()} == {
        "SCGC_AM-100-B01": (None, "ambiguous_base_count"), "SCGC_AM-100-B02": (None, "ambiguous_base_count")}


# --- staging, end to end -------------------------------------------------------------------------------------------

def test_stage_keeps_cells_r232_lists_and_records_their_status_deterministically(world):
    mod = _mod()
    rep1 = run_stage(mod, world, "a")
    rep2 = run_stage(mod, world, "b")
    t = world["tmp"]
    meta, fa = t / "meta_a" / "SAG2", t / "fasta_a" / "SAG2" / "sag2_gorg_dark.fna.gz"

    # nothing is removed because r232 lists it: the status is recorded for every genome (coordinator, 2026-11-26)
    assert not (meta / "removed.tsv").exists()
    st = {r["genome"]: r for r in csv.DictReader(open(meta / "r232_status.tsv"), delimiter="\t")}
    assert {g: r["r232_status"] for g, r in st.items()} == {
        "SCGC_AM-100-A01": "in_gtdb_r232", "SCGC_AM-100-A02": "in_globdb_r232", "SCGC_AM-100-A03": "not_in_r232",
        "SCGC_AM-100-A04": "not_in_r232", "SCGC_AM-100-A05": "not_visible", "SCGC_AM-100-A06": "not_visible",
        "SCGC_AM-100-A07": "not_visible"}
    a1 = st["SCGC_AM-100-A01"]
    assert (a1["in_gtdb_r232"], a1["is_gtdb_representative"], a1["gtdb_species"], a1["in_globdb_r232"]) == (
        "yes", "True", "Species one", "no")
    a2 = st["SCGC_AM-100-A02"]
    assert (a2["in_gtdb_r232"], a2["in_globdb_r232"], a2["globdb_r232_id"]) == ("no", "yes", "GEM_000002")

    kept = sorted(SAGS)
    genomes = json.loads((meta / "genomes.json").read_text())
    assert [g["contig_uids"] for g in genomes] == [[cuid(q) for _, q in SAGS[s]] for s in kept]
    assert {g["kind"] for g in genomes} == {"prokaryote"} and {g["assembly_quality"] for g in genomes} == {"contig"}
    assert all(g["circular"] == [] for g in genomes)
    codes = {s: g["genetic_code"] for s, g in zip(kept, genomes)}
    assert codes == {"SCGC_AM-100-A01": 11, "SCGC_AM-100-A02": 11, "SCGC_AM-100-A03": 11, "SCGC_AM-100-A04": 4,
                     "SCGC_AM-100-A05": 11, "SCGC_AM-100-A06": 11, "SCGC_AM-100-A07": 11}

    # the FASTA: BGZF (gzip members with the BC field, EOF block last), the kept cells only, in SAG order
    raw = fa.read_bytes()
    assert raw[:4] == bytes.fromhex("1f8b0804") and raw[12:14] == b"BC"
    assert raw.endswith(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))
    text = gzip.decompress(raw).decode()
    names = [ln[1:] for ln in text.splitlines() if ln.startswith(">")]
    assert len(names) == sum(len(SAGS[s]) for s in kept)
    # a contig name that carries the cell's id stays as it is; one that does not is given the member's stem
    assert names == ["AM-100-A01_NODE_1", "AM-100-A01_NODE_2", "AM-100-A02_NODE_1", "AM-100-A02_NODE_2",
                     "AM-100-A03_NODE_1", "AM-100-A04_NODE_1", "AM-100-A05_NODE_1", "AM-100-A06_NODE_1",
                     "SCGC_AM-100-A07|contig1", "SCGC_AM-100-A07|contig2"]
    assert len(set(n.split()[0] for n in names)) == len(names)  # one name per record in the file

    # match table, ownership, later holders, codes, objects
    mt = {r["sag"]: r for r in csv.DictReader(open(meta / "sag_match.tsv"), delimiter="\t")}
    assert mt["SCGC_AM-100-A05"]["seen_by_accession_test"] == "no"
    assert mt["SCGC_AM-100-A03"]["seen_by_accession_test"] == "yes" and mt["SCGC_AM-100-A03"]["status"] == "kept"
    shared = list(csv.DictReader(open(meta / "shared_records.tsv"), delimiter="\t"))
    assert [(r["genome"], r["owner_stage"]) for r in shared] == [("SCGC_AM-100-A07", "V2")]
    later = list(csv.DictReader(open(meta / "later_records.tsv"), delimiter="\t"))
    assert [(r["genome"], r["holder"]) for r in later] == [("SCGC_AM-100-A03", "SAG1")]
    gcs = list(csv.DictReader(open(meta / "genetic_codes.tsv"), delimiter="\t"))
    assert [r["genome"] for r in gcs] == kept
    objs = json.loads((meta / "objects.json").read_text())
    assert objs[0]["path"] == "/tmp/if2453/s2/inputs/SAG2/sag2_gorg_dark.fna.gz"
    assert objs[0]["redistributable"] is False and objs[0]["source_name"] == "gorg_dark_sags"
    assert "CC BY 4.0" in objs[0]["license"]
    assert json.loads((meta / "objects_staging.json").read_text())[0]["path"] == str(fa)

    r = rep1
    assert r["sags_in_zip"] == 7 and "removed" not in r
    assert r["r232_status"] == {"in_gtdb_r232": 1, "in_globdb_r232": 1, "not_in_r232": 2, "not_visible": 3}
    assert r["sags_staged"] == 7 and r["genomes"] == 7
    assert r["seen"] == {"matched_to_ena": 4, "not_matched": 3, "kept_after_test": 4}
    assert r["held_by"] == {"V2": 1} and r["later_held_by"] == {"SAG1": 1}

    for f in ("genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
              "later_records.tsv", "genetic_codes.tsv", "sag_match.tsv", "r232_status.tsv"):
        assert (meta / f).read_bytes() == (t / "meta_b" / "SAG2" / f).read_bytes(), f
    assert raw == (t / "fasta_b" / "SAG2" / "sag2_gorg_dark.fna.gz").read_bytes()
    skip = "objects_staging.json"  # it names the output directory, which differs between the two runs
    assert {k: v for k, v in rep1["sha256"].items() if k != skip} == {
        k: v for k, v in rep2["sha256"].items() if k != skip}

    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), "SAG2", str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


def test_stage_lists_a_repeated_contig_set_once(world, tmp_path):
    mod = _mod()
    z = tmp_path / "dup.zip"
    same = SAGS["SCGC_AM-100-A07"]
    write_zip(z, {"GORG-Dark SAGs/SCGC_AM-100-B01_contigs.fasta.gz": _fa_gz(same),
                  "GORG-Dark SAGs/SCGC_AM-100-B02_contigs.fasta.gz": _fa_gz(same)})
    w = {**world, "zip": z}
    rep = run_stage(mod, w, "d")
    meta = tmp_path / "meta_d" / "SAG2"
    assert len(json.loads((meta / "genomes.json").read_text())) == 1
    rep_rows = list(csv.DictReader(open(meta / "repeated_records.tsv"), delimiter="\t"))
    assert [(r["genome"], r["first_genome"]) for r in rep_rows] == [("SCGC_AM-100-B02", "SCGC_AM-100-B01")]
    assert rep["genomes_repeated"] == 1


def test_stage_refuses_an_existing_output_a_missing_ena_table_and_an_empty_zip_and_writes_nothing(world, tmp_path):
    mod = _mod()
    run_stage(mod, world, "a")
    with pytest.raises(mod.Refused, match="exists"):
        run_stage(mod, world, "a")
    w = {**world, "ena": tmp_path / "nope.tsv"}
    with pytest.raises(mod.Refused, match="ENA table"):
        run_stage(mod, w, "m")
    assert not (tmp_path / "meta_m" / "SAG2").exists()
    z = tmp_path / "empty.zip"
    write_zip(z, {})
    with pytest.raises(mod.Refused, match="no single cell"):
        run_stage(mod, {**world, "zip": z}, "e")
    assert not (tmp_path / "fasta_e" / "SAG2" / "sag2_gorg_dark.fna.gz").exists()


def test_stage_refuses_a_cell_whose_member_is_not_gzip_and_leaves_no_partial_file(world, tmp_path):
    mod = _mod()
    z = tmp_path / "bad.zip"
    write_zip(z, {MEMBER["SCGC_AM-100-A03"]: b"not gzip at all"})
    with pytest.raises(mod.Refused, match="SCGC_AM-100-A03"):
        run_stage(mod, {**world, "zip": z}, "x")
    d = tmp_path / "fasta_x" / "SAG2"
    assert not d.exists() or [p.name for p in d.iterdir()] == []
    assert not (tmp_path / "meta_x" / "SAG2").exists()


# --- the conf template ---------------------------------------------------------------------------------------------

def test_template_passes_stage_est_and_carries_the_cards_values():
    est = _load("stage_est", "stage_est.py")
    c = est.read_conf(TEMPLATES / "SAG2.conf", "SAG2", tree=REPO)
    assert (c.depth, c.crispr, c.library, c.after) == ("called", True, None, ("ASG1",))
    assert str(c.inputs) == "/tmp/if2453/s2/inputs/SAG2" and c.min_contig_len is None
    assert getattr(c, "genecall_mode", "default") == "default"
    assert NOTE in c.notes
    assert {n for _, n in c.sums} == {"genomes.json", "objects.json", "sag2_gorg_dark.fna.gz"}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    for h, _ in c.sums:
        assert h in note, f"SUM {h} not in the note"
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    allowed = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}
    allowed |= {s for s, _ in a26.R26_11_4_CHAIN}
    assert set(c.expect) <= allowed
    for owner, (n, _) in c.expect.items():
        assert n > 0

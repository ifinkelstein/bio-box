"""Card C-80 (DEC-79, DEC-80): Kavagutti's HiFi MAGs and NCLDVs (`KAVHIFI`) and INPHARED's 7Apr2026 delta (`INPH2604`).

Every expectation is written out by hand from the fixture (constructed truth): the parts taken from the archive by
member name, the INPHARED diff by accession and by `contig_uid`, the kinds and codes, the byte identity of database
records, the ownership lists, the accession test. contig_uid is computed here with blake3 directly. Each test runs
twice, with the normal PATH and with PATH reduced to `/usr/bin:/bin`: nothing here needs a tool of the pod's image
(GitHub Actions has none of them).
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import random
import subprocess
import sys
import tarfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-80")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-26_kavhifi_inph2604_inputs.md"


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


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
    return _load("c80_kavagutti_inphared", "c80_kavagutti_inphared.py")


def cuid(seq: str) -> str:
    return blake3(seq.upper().encode()).digest(length=16).hex()


def seq_of(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def fa(records: list[tuple[str, str]], width: int = 60) -> str:
    out = io.StringIO()
    for hdr, s in records:
        out.write(f">{hdr}\n")
        for i in range(0, len(s), width):
            out.write(s[i:i + width] + "\n")
    return out.getvalue()


def tar_gz(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as t:
        for name, data in members.items():
            if name.endswith("/"):
                info = tarfile.TarInfo(name.rstrip("/"))
                info.type = tarfile.DIRTYPE
                t.addfile(info)
                continue
            info = tarfile.TarInfo(name)
            info.size = len(data)
            t.addfile(info, io.BytesIO(data))
    return buf.getvalue()


# --- the Kavagutti archive fixture ---------------------------------------------------------------------------------

UNIT = seq_of(1, 4_960)
S = {
    "M1": seq_of(11, 3_000), "M1b": seq_of(12, 2_000),   # MAG M1: two contigs
    "M2": seq_of(13, 4_000),                              # MAG M2: the screen calls code 25
    "M3": seq_of(14, 3_500),                              # MAG M3: its name holds a GCA_ in GTDB R232
    "M4": seq_of(15, 2_500),                              # MAG M4: accession in neither release
    "N1": seq_of(21, 6_000),                              # NCLDV, also held by an earlier stage (V2)
    "N2": UNIT + UNIT[:40],                               # NCLDV with an exact 40 bp terminal repeat
    "N3": seq_of(22, 5_000),                              # NCLDV, held by a stage of r26.11.5 before this one
    "PHAGE": seq_of(23, 4_000), "EUK": seq_of(24, 3_000), "MITO": seq_of(25, 1_500), "TREE": seq_of(26, 100),
}
MAG_FILES = {
    "MAGs-fnas/HRS-A-bin-1.fasta": fa([("ctgA1 len=3000", S["M1"]), ("ctgA2", S["M1b"])], 70),
    "MAGs-fnas/HRS-A-bin-2.fasta": fa([("ctgB1", S["M2"])], 80),
    "MAGs-fnas/HRS-B-bin-3_GCA_900000001.1.fasta": fa([("ctgC1", S["M3"])]),
    "MAGs-fnas/HRS-C-bin-4.fasta": fa([("ctgD1 binned from GCF_900000009.1", S["M4"])]),
}
NCLDV_FILES = {
    "NCLDV/": b"", "NCLDV/NCLDV-tree/": b"", "NCLDV/NCLDV-fna/": b"",
    "NCLDV/NCLDV-tree/concat.aln.fna": fa([("aln", S["TREE"])]).encode(),
    "NCLDV/NCLDV-fna/HRS-EL6-C8.fna": fa([("HRS-EL6-C8", S["N1"])]).encode(),
    "NCLDV/NCLDV-fna/HRS-EL7-C102.fna": fa([("HRS-EL7-C102", S["N2"])]).encode(),
    "NCLDV/NCLDV-fna/HRS-EL7-C103.fna": fa([("HRS-EL7-C103", S["N3"])]).encode(),
}
EUK_TAR = tar_gz({"euk/E1.fna": fa([("e1", S["EUK"]), ("e2", S["EUK"][:1000])]).encode()})
MITO_TAR = tar_gz({"mito/": b"", "mito/HRS-EL3-C80_Synura_synuroidea.fna": fa([("mt", S["MITO"])]).encode(),
                   "mito/HRS-HL1-C288_Pleodorina_starrii.fna": fa([("mt2", S["MITO"][:900])]).encode()})
BAC_TAR = tar_gz({"MAGs-fnas/": b"", **{k: v.encode() for k, v in MAG_FILES.items()}})
NCLDV_TAR = tar_gz({k: v for k, v in NCLDV_FILES.items()})
PHAGE_TAR = tar_gz({"phages-fna/P1.fna": fa([("p1", S["PHAGE"])]).encode()})
SUPPORT = tar_gz({"figshare/": b"", "figshare/phages-fna.tar.gz": PHAGE_TAR, "figshare/euk-MAGs-fna.tar.gz": EUK_TAR,
                  "figshare/bac-MAGs-fnas.tar.gz": BAC_TAR, "figshare/mitochondrias-fna.tar.gz": MITO_TAR,
                  "figshare/NCLDV-data.tar.gz": NCLDV_TAR, "figshare/README": b"readme\n"})


def fake_screen(seqs):
    """The screen's result by the content: the MAG whose first contig is M2's is of code 25."""
    code = 25 if seqs[0].decode() == S["M2"] else 11
    return {"genetic_code": code, "reason": "stub", "cd_11": 0.5, "mean_orf_11": 300.0}


def fake_lookup(accs):
    held = {"GCA_900000001": (True, False)}
    return [{"accession": a, "in_gtdb_r232": held.get(a, (False, False))[0],
             "in_globdb_r232": held.get(a, (False, False))[1]} for a in accs]


def owner_dir(tmp_path: Path, name: str, seqs: list[str]) -> Path:
    d = tmp_path / "owners" / name
    d.mkdir(parents=True)
    (d / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s) for s in seqs]}]))
    return d


def run_kav(mod, tmp_path: Path, tag: str):
    support = tmp_path / "Support-data.gz"
    if not support.exists():
        support.write_bytes(SUPPORT)
        owner_dir(tmp_path, "V2", [S["N1"]])    # an r26.11.4 stage holds N1
        owner_dir(tmp_path, "GV1", [S["N3"]])   # a stage of r26.11.5 before KAVHIFI holds N3
    return mod.stage_kavagutti(
        support, tmp_path / f"meta_{tag}", tmp_path / f"fasta_{tag}", root=None, root_stages=(),
        owner_inputs={"V2": tmp_path / "owners" / "V2"}, later_inputs={"GV1": tmp_path / "owners" / "GV1"},
        inputs_home=Path("/tmp/if2453/s2/inputs"), lookup=fake_lookup, work=tmp_path / f"work_{tag}", workers=0,
        screen=fake_screen)


def test_kav_takes_the_mags_and_ncldv_by_member_name_and_sets_kinds_and_codes(tmp_path):
    mod = _mod()
    rep = run_kav(mod, tmp_path, "a")
    meta, fdir = tmp_path / "meta_a" / "KAVHIFI", tmp_path / "fasta_a" / "KAVHIFI"
    genomes = json.loads((meta / "genomes.json").read_text())
    by_first = {g["contig_uids"][0]: g for g in genomes}
    # M3 (accession in GTDB R232) is not staged; M1, M2, M4 are; the three NCLDV files are; nothing of the phages
    # (LRV1 took them), the trees, the eukaryotic MAGs or the mitochondria
    assert sorted(by_first) == sorted([cuid(S["M1"]), cuid(S["M2"]), cuid(S["M4"]), cuid(S["N1"]), cuid(S["N3"]),
                                       cuid(S["N2"][:-40])])
    assert by_first[cuid(S["M1"])]["contig_uids"] == [cuid(S["M1"]), cuid(S["M1b"])]
    assert by_first[cuid(S["M1"])]["kind"] == "prokaryote" and by_first[cuid(S["M1"])]["genetic_code"] == 11
    assert by_first[cuid(S["M2"])]["genetic_code"] == 25
    assert by_first[cuid(S["M4"])]["kind"] == "prokaryote"
    assert by_first[cuid(S["N1"])]["kind"] == "virus" and by_first[cuid(S["N1"])]["genetic_code"] == 11
    assert {g["assembly_quality"] for g in genomes} == {"contig"}
    # DEC-65 (1): N2's exact 40 bp terminal repeat is trimmed once and the record staged circular
    n2 = by_first[cuid(S["N2"][:-40])]
    assert n2["circular"] == [cuid(S["N2"][:-40])]
    assert all(g["circular"] == [] for g in genomes if g is not n2)
    assert rep["trims"] == 1 and rep["trimmed_bases"] == 40
    text = (fdir / "kavagutti_hifi_ncldv.fna").read_text()
    assert S["N2"][:-40] in text.replace("\n", "") and S["N2"] not in text.replace("\n", "")
    mags = (fdir / "kavagutti_hifi_bacterial_mags.fna").read_text()
    assert ">ctgA1 len=3000" in mags and S["M3"] not in mags.replace("\n", "")
    for other in ("PHAGE", "EUK", "MITO", "TREE"):
        assert S[other] not in text.replace("\n", "") + mags.replace("\n", "")
    skipped = (meta / "skipped_members.tsv").read_text().splitlines()
    assert skipped[1:] == ["HRS-B-bin-3_GCA_900000001.1\tbac_mag\taccession_in_gtdb_r232_or_globdb_r232\tGCA_900000001"]
    objs = json.loads((meta / "objects.json").read_text())
    assert [o["path"] for o in objs] == ["/tmp/if2453/s2/inputs/KAVHIFI/kavagutti_hifi_bacterial_mags.fna",
                                         "/tmp/if2453/s2/inputs/KAVHIFI/kavagutti_hifi_ncldv.fna"]
    assert all(o["redistributable"] is True and o["license"].startswith("CC BY 4.0")
               and o["source_name"] == "kavagutti_freshwater_hifi_phages" for o in objs)
    staging = json.loads((meta / "objects_staging.json").read_text())
    assert [o["path"] for o in staging] == [str(fdir / "kavagutti_hifi_bacterial_mags.fna"),
                                            str(fdir / "kavagutti_hifi_ncldv.fna")]


def test_kav_defers_the_eukaryotic_mags_and_the_mitochondria_without_staging_them(tmp_path):
    mod = _mod()
    rep = run_kav(mod, tmp_path, "a")
    meta, fdir = tmp_path / "meta_a" / "KAVHIFI", tmp_path / "fasta_a" / "KAVHIFI"
    assert (fdir / "deferred" / "euk-MAGs-fna.tar.gz").read_bytes() == EUK_TAR  # byte for byte
    assert (meta / "deferred" / "mitochondrias-fna.tar.gz").read_bytes() == MITO_TAR
    assert (meta / "euk_deferred.tsv").read_text().splitlines() == ["file\trecords\tbases", "euk/E1.fna\t2\t4000"]
    assert (meta / "mito_deferred.tsv").read_text().splitlines() == [
        "file\trecords\tbases", "mito/HRS-EL3-C80_Synura_synuroidea.fna\t1\t1500",
        "mito/HRS-HL1-C288_Pleodorina_starrii.fna\t1\t900"]
    assert rep["deferred"]["euk"]["files"] == 1 and rep["deferred"]["mito"]["bases"] == 2400
    assert rep["deferred"]["mito"]["path"] == str(meta / "deferred" / "mitochondrias-fna.tar.gz")


def test_kav_accession_test_counts_what_the_index_could_and_could_not_see(tmp_path):
    mod = _mod()
    rep = run_kav(mod, tmp_path, "a")
    t = rep["accession_test"]
    # four MAGs: M3 and M4 carry an INSDC accession (the index sees them; M3 is in GTDB R232), M1 and M2 carry none
    assert t["mags_staged_or_tested"] == 4
    assert t["mags_with_an_insdc_accession"] == 2
    assert t["mags_without_one_not_visible_to_the_index"] == 2
    assert t["found_in_gtdb_r232"] == 1 and t["found_in_globdb_r232"] == 0
    rows = (tmp_path / "meta_a" / "KAVHIFI" / "mag_accessions.tsv").read_text().splitlines()
    assert rows == ["genome\taccession\tin_gtdb_r232\tin_globdb_r232",
                    "HRS-A-bin-1\t\t\t", "HRS-A-bin-2\t\t\t",
                    "HRS-B-bin-3_GCA_900000001.1\tGCA_900000001\tTrue\tFalse",
                    "HRS-C-bin-4\tGCF_900000009\tFalse\tFalse"]


def test_the_accession_lookup_reads_card_c74s_index(tmp_path):
    mod = _mod()
    pq.write_table(pa.table({
        "accession": ["GCA_900000001", "GCF_900000002"], "genbank_accession": ["GCA_900000001", "GCA_900000002"],
        "is_representative": [True, False], "representative": ["GCA_900000001", "GCA_900000001"],
        "species": ["s__A a", "s__A a"]}), tmp_path / "gtdb_r232.parquet")
    pq.write_table(pa.table({"globdb_id": ["TPMC_1"], "insdc_accession": ["GCA_900000077"]}),
                   tmp_path / "globdb_r232.parquet")
    look = mod.accession_lookup(tmp_path)
    rows = look(["GCA_900000001", "GCA_900000002", "GCA_900000077", "GCF_900000009"])
    assert [(r["in_gtdb_r232"], r["in_globdb_r232"]) for r in rows] == [
        (True, False), (True, False), (False, True), (False, False)]
    assert mod.accessions_in("x_GCA_900000001.1_y", "GCF_900000009 and GCA_900000001") == [
        "GCA_900000001", "GCF_900000009"]


def test_kav_lists_ownership_by_exact_contig_uid_and_keeps_the_genomes(tmp_path):
    mod = _mod()
    rep = run_kav(mod, tmp_path, "a")
    meta = tmp_path / "meta_a" / "KAVHIFI"
    shared = [ln.split("\t") for ln in (meta / "shared_records.tsv").read_text().splitlines()[1:]]
    later = [ln.split("\t") for ln in (meta / "later_records.tsv").read_text().splitlines()[1:]]
    assert [(r[1], r[2], r[3]) for r in shared] == [("HRS-EL6-C8", cuid(S["N1"]), "V2")]
    assert [(r[1], r[2], r[3]) for r in later] == [("HRS-EL7-C103", cuid(S["N3"]), "GV1")]
    assert rep["owner_contigs"] == {"V2": 1} and rep["later_contigs"] == {"GV1": 1}
    genomes = json.loads((meta / "genomes.json").read_text())
    assert len(genomes) == 6  # the held contigs stay in their genomes (DEC-26)
    assert rep["new_contigs"] == rep["contigs"] - 1


def test_kav_is_byte_identical_between_runs_and_passes_the_stage_checker(tmp_path):
    mod = _mod()
    run_kav(mod, tmp_path, "a")
    run_kav(mod, tmp_path, "b")
    for rel in ("genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
                "later_records.tsv", "trims.tsv", "genetic_codes.tsv", "mag_accessions.tsv", "skipped_members.tsv",
                "euk_deferred.tsv", "mito_deferred.tsv", "deferred/mitochondrias-fna.tar.gz"):
        assert (tmp_path / "meta_a" / "KAVHIFI" / rel).read_bytes() == (tmp_path / "meta_b" / "KAVHIFI" / rel).read_bytes(), rel
    for name in ("kavagutti_hifi_bacterial_mags.fna", "kavagutti_hifi_ncldv.fna"):
        assert (tmp_path / "fasta_a" / "KAVHIFI" / name).read_bytes() == (tmp_path / "fasta_b" / "KAVHIFI" / name).read_bytes()
    meta = tmp_path / "meta_a" / "KAVHIFI"
    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), "KAVHIFI", str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    doc = json.loads(chk.stdout)
    assert doc["problems"] == [] and doc["kinds"] == ["prokaryote", "virus"] and doc["codes"] == [11, 25]
    assert doc["gene_call_mode"]["virus"] == "gv-meta"


def test_kav_refuses_an_archive_without_the_mag_member_and_writes_nothing(tmp_path):
    mod = _mod()
    bad = tmp_path / "bad.gz"
    bad.write_bytes(tar_gz({"figshare/NCLDV-data.tar.gz": NCLDV_TAR}))
    with pytest.raises(mod.Refused, match="bac-MAGs-fnas"):
        mod.stage_kavagutti(bad, tmp_path / "m", tmp_path / "f", root=None, root_stages=(), owner_inputs={},
                            later_inputs={}, inputs_home=Path("/x"), work=tmp_path / "w", workers=0,
                            screen=fake_screen)
    assert not (tmp_path / "m" / "KAVHIFI").exists()
    assert not (tmp_path / "f" / "KAVHIFI" / "kavagutti_hifi_ncldv.fna").exists()


def test_the_default_screen_calls_code_11_for_a_genome_under_100_kb():
    mod = _mod()
    r = mod.screen_codes([seq_of(5, 3_000).encode()])
    assert r["genetic_code"] == 11 and r["reason"] == "below_100kb"


# --- INPHARED -----------------------------------------------------------------------------------------------------

def _root_with(tmp_path: Path, labels: dict[str, list[tuple[str, str]]], index: dict[str, list[str]]) -> Path:
    """A root holding contig_label (uid, accession) of the stages in LABELS and contig_index of those in INDEX."""
    root = tmp_path / "root"
    for stage, rows in labels.items():
        d = root / "l0seq" / "contig_label" / f"stage={stage}"
        d.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.table({"contig_uid": pa.array([bytes.fromhex(cuid(s)) for _, s in rows], pa.binary()),
                                 "accession": [a for a, _ in rows]}), d / "part-0.parquet")
    for stage, seqs in index.items():
        d = root / "l0seq" / "contig_index" / f"stage={stage}"
        d.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.table({"contig_uid": pa.array([bytes.fromhex(cuid(s)) for s in seqs], pa.binary()),
                                 "is_default_variant": [True] * len(seqs)}), d / "part-0.parquet")
    return root


I = {k: seq_of(100 + i, 1_200 + 100 * i) for i, k in enumerate("ABCDEFGH")}
# the 7Apr2026 archive, record by record (accession, sequence):
#   A  in S1 under the same accession and sequence            -> not staged
#   B  in SAT1 under the same accession and sequence          -> not staged
#   C  new accession, new sequence                            -> staged
#   D  new accession, but S1 holds the same sequence          -> skipped (a renamed copy)
#   E  accession in S1, sequence changed (not held)           -> skipped (crossed)
#   F  new accession, sequence held by H1 only                -> staged, listed as owned by H1
#   G  new accession, sequence held by a later stage (GV1)    -> staged, listed as later
#   H  new, 200 kb long: counted as jumbo (header names Wolbachia)
ARCHIVE = [("AB000001", I["A"]), ("AB000002", I["B"]), ("AB000003", I["C"]), ("AB000004", I["D"]),
           ("AB000005", seq_of(900, 1_250)), ("AB000006", I["F"]), ("AB000007", I["G"]),
           ("AB000008", seq_of(901, 200_000))]
HEADERS = {a: f"{a} Test phage {a}, complete genome" if a != "AB000008" else f"{a} Wolbachia phage WO, complete genome"
           for a, _ in ARCHIVE}


def inph_fixture(tmp_path: Path):
    text = "".join(f">{HEADERS[a]}\n" + "\n".join(s[i:i + 70] for i in range(0, len(s), 70)) + "\n" for a, s in ARCHIVE)
    path = tmp_path / "7Apr2026_genomes.fa.gz"
    path.write_bytes(gzip.compress(text.encode(), mtime=0))
    root = _root_with(
        tmp_path,
        labels={"S1": [("AB000001.1", I["A"]), ("AB000005", seq_of(902, 1_300)), ("XX000001", seq_of(903, 1_000))]},
        index={"S1": [I["A"], I["D"], seq_of(902, 1_300), seq_of(903, 1_000)], "H1": [I["F"]]})
    sat1 = tmp_path / "owners" / "SAT1"
    sat1.mkdir(parents=True)
    (sat1 / "records.tsv").write_text(
        "ordinal\tsource\tcategory\taccession\tcontig_uid\tlength\tcircular\towner_stage\n"
        f"0\tx\ty\tAB000002.1\t{cuid(I['B'])}\t{len(I['B'])}\tfalse\t\n")
    (sat1 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(I["B"])]}]))
    gv1 = tmp_path / "owners" / "GV1"
    gv1.mkdir(parents=True)
    (gv1 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(I["G"])]}]))
    return path, root, sat1, gv1


def run_inph(mod, tmp_path: Path, tag: str):
    if not (tmp_path / "root").exists():
        inph_fixture(tmp_path)
    path, root = tmp_path / "7Apr2026_genomes.fa.gz", tmp_path / "root"
    sat1, gv1 = tmp_path / "owners" / "SAT1", tmp_path / "owners" / "GV1"
    return mod.stage_inphared(
        path, tmp_path / f"meta_{tag}", tmp_path / f"fasta_{tag}", root=root, root_stages=("S1", "H1"),
        owner_inputs={"SAT1": sat1}, later_inputs={"GV1": gv1}, known_root_stages=("S1",),
        known_inputs={"SAT1": sat1}, inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / f"work_{tag}")


def test_inphared_delta_is_diffed_by_accession_and_by_contig_uid_against_s1_and_sat1(tmp_path):
    mod = _mod()
    rep = run_inph(mod, tmp_path, "a")
    meta = tmp_path / "meta_a" / "INPH2604"
    genomes = json.loads((meta / "genomes.json").read_text())
    assert [g["contig_uids"] for g in genomes] == [[cuid(I["C"])], [cuid(I["F"])], [cuid(I["G"])],
                                                   [cuid(ARCHIVE[7][1])]]
    assert {g["kind"] for g in genomes} == {"virus"} and {g["genetic_code"] for g in genomes} == {11}
    assert {g["assembly_quality"] for g in genomes} == {"complete"} and all(g["circular"] == [] for g in genomes)
    assert rep["cells"] == {"accession_known+sequence_held": 2, "accession_known+sequence_new": 1,
                            "accession_new+sequence_held": 1, "accession_new+sequence_new": 4}
    skipped = [ln.split("\t") for ln in (meta / "skipped_records.tsv").read_text().splitlines()[1:]]
    assert [(r[0], r[3]) for r in skipped] == [
        ("AB000004", "accession_new_sequence_held_by_S1_or_SAT1"), ("AB000005", "accession_known_sequence_not_held")]
    assert rep["archive_records"] == 8 and rep["genomes"] == 4
    assert rep["jumbo_ge_200kb"] == 1 and rep["wolbachia"] == 1
    shared = [ln.split("\t") for ln in (meta / "shared_records.tsv").read_text().splitlines()[1:]]
    later = [ln.split("\t") for ln in (meta / "later_records.tsv").read_text().splitlines()[1:]]
    assert [(r[1], r[3]) for r in shared] == [("AB000006", "H1")]
    assert [(r[1], r[3]) for r in later] == [("AB000007", "GV1")]


def test_inphared_database_records_are_byte_identical_to_the_archive(tmp_path):
    mod = _mod()
    run_inph(mod, tmp_path, "a")
    want = "".join(f">{HEADERS[a]}\n" + "\n".join(s[i:i + 70] for i in range(0, len(s), 70)) + "\n"
                   for a, s in ARCHIVE if a in ("AB000003", "AB000006", "AB000007", "AB000008"))
    got = (tmp_path / "fasta_a" / "INPH2604" / "inphared_7apr2026_delta.fna").read_text()
    assert got == want
    objs = json.loads((tmp_path / "meta_a" / "INPH2604" / "objects.json").read_text())
    assert objs[0]["path"] == "/tmp/if2453/s2/inputs/INPH2604/inphared_7apr2026_delta.fna"
    assert objs[0]["source_name"] == "inphared_7apr2026" and objs[0]["redistributable"] is True


def test_inphared_is_deterministic_and_passes_the_stage_checker(tmp_path):
    mod = _mod()
    run_inph(mod, tmp_path, "a")
    run_inph(mod, tmp_path, "b")
    for rel in ("genomes.json", "objects.json", "records.tsv", "skipped_records.tsv", "shared_records.tsv",
                "later_records.tsv"):
        assert (tmp_path / "meta_a" / "INPH2604" / rel).read_bytes() == (tmp_path / "meta_b" / "INPH2604" / rel).read_bytes()
    assert (tmp_path / "fasta_a" / "INPH2604" / "inphared_7apr2026_delta.fna").read_bytes() == (
        tmp_path / "fasta_b" / "INPH2604" / "inphared_7apr2026_delta.fna").read_bytes()
    meta = tmp_path / "meta_a" / "INPH2604"
    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), "INPH2604", str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


def test_inphared_refuses_when_nothing_is_new_and_writes_nothing(tmp_path):
    mod = _mod()
    path, root, sat1, gv1 = inph_fixture(tmp_path)
    only_known = tmp_path / "known.fa.gz"
    only_known.write_bytes(gzip.compress(f">AB000001 x\n{I['A']}\n".encode(), mtime=0))
    with pytest.raises(mod.Refused, match="no INPHARED record is new"):
        mod.stage_inphared(only_known, tmp_path / "m", tmp_path / "f", root=root, root_stages=("S1",),
                           owner_inputs={}, later_inputs={}, known_root_stages=("S1",), known_inputs={"SAT1": sat1},
                           inputs_home=Path("/x"), work=tmp_path / "w")
    assert not (tmp_path / "m" / "INPH2604").exists()
    assert not (tmp_path / "f" / "INPH2604" / "inphared_7apr2026_delta.fna").exists()


# --- the conf templates -------------------------------------------------------------------------------------------

def test_templates_pass_stage_est_and_carry_the_cards_values():
    est = _load("stage_est", "stage_est.py")
    ck = est.read_conf(TEMPLATES / "KAVHIFI.conf", "KAVHIFI", tree=REPO)
    ci = est.read_conf(TEMPLATES / "INPH2604.conf", "INPH2604", tree=REPO)
    assert (ck.depth, ck.crispr, ck.library, ck.after) == ("called", True, None, ("MVEUKV",))
    assert (ci.depth, ci.crispr, ci.library, ci.after) == ("annotated", True, ("pfam", "38.2", "gathering"),
                                                           ("KAVHIFI",))
    assert str(ck.inputs) == "/tmp/if2453/s2/inputs/KAVHIFI" and str(ci.inputs) == "/tmp/if2453/s2/inputs/INPH2604"
    assert ck.min_contig_len is None and ci.min_contig_len is None
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    allowed = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}
    allowed |= {s for s, _ in a26.R26_11_4_CHAIN}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    want = {"KAVHIFI": {"genomes.json", "objects.json", "kavagutti_hifi_bacterial_mags.fna", "kavagutti_hifi_ncldv.fna"},
            "INPH2604": {"genomes.json", "objects.json", "inphared_7apr2026_delta.fna"}}
    for c in (ck, ci):
        assert NOTE in c.notes
        assert {n for _, n in c.sums} == want[c.stage]
        for h, _ in c.sums:
            assert h in note, f"{c.stage}: SUM {h} not in the note"
        assert set(c.expect) <= allowed
        for owner, (n, _) in c.expect.items():
            assert n > 0

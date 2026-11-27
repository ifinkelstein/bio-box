"""Card C-83 (DEC-85): stage the two coral-reef sets, Tara Pacific's Reef Microbiomics Database (`RMD1`) and the
Great Barrier Reef genomes (`GBR_MGD`, stage `GBRMGD`).

Oracles are written here separately from the script: contig_uid by blake3 directly; each fixture member's
compression chosen here and hidden behind names that say otherwise (a gzip member named `.fa`, a plain member named
`.fa.gz`), so that only the magic bytes can tell; record and base counts by hand from the fixture records; GTDB R232
membership and card C-74's index as hand-written fixture tables. Each test runs twice, with the normal PATH and with
PATH reduced to `/usr/bin:/bin` (as tests/r1/test_r232_index.py does): the stager needs only Python and its
libraries, and GitHub Actions has none of the pod's tools.
"""

from __future__ import annotations

import bz2
import csv
import gzip
import importlib.util
import io
import json
import lzma
import random
import subprocess
import sys
import tarfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-83")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-27_coral_inputs.md"
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
    return _load("c83_coral_inputs", "c83_coral_inputs.py")


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def cuid(seq: str) -> str:
    return blake3(seq.upper().encode()).digest(length=16).hex()


def seq_of(seed: int, n: int, alphabet: str = "ACGT") -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(alphabet) for _ in range(n))


def _fa(records: list[tuple[str, str]], width: int = 60) -> bytes:
    out = io.StringIO()
    for hdr, s in records:
        out.write(f">{hdr}\n")
        for i in range(0, len(s), width):
            out.write(s[i:i + width] + "\n")
    return out.getvalue().encode()


def _gz_lines(path: Path) -> list[str]:
    return gzip.decompress(path.read_bytes()).decode().splitlines()


def _rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tar(path: Path, members: list[tuple[str, bytes]], *, outer: str = "bz2", dirs: list[str] = ()) -> None:
    mode = {"bz2": "w:bz2", "gz": "w:gz", "none": "w"}[outer]
    with tarfile.open(path, mode) as tf:
        for d in dirs:
            ti = tarfile.TarInfo(d)
            ti.type = tarfile.DIRTYPE
            ti.mtime = 0
            tf.addfile(ti)
        for name, data in members:
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            ti.mtime = 0
            tf.addfile(ti, io.BytesIO(data))


# ==================================================================================================================
# the census: every member named, its compression by magic bytes, its content by its bytes
# ==================================================================================================================

A1, A2 = seq_of(1, 2500), seq_of(2, 1200)
B1 = seq_of(3, 3100)
C1 = seq_of(4, 2000)
D1 = seq_of(5, 1800)
PROT = seq_of(6, 300, "ACDEFGHIKLMNPQRSTVWY")
ALN = seq_of(7, 400)[:200] + "-" * 30 + seq_of(8, 170)
CENSUS_MEMBERS = [
    ("bins/A.fa", _fa([("A_1", A1), ("A_2", A2)])),                  # plain, named plain
    ("bins/B.fa", gzip.compress(_fa([("B_1", B1)]), mtime=0)),        # gzip, named plain
    ("bins/C.fna.bz2", bz2.compress(_fa([("C_1", C1)]))),             # bzip2, named so
    ("bins/D.fa.gz", _fa([("D_1", D1)])),                             # plain, named gzip
    ("bins/E.fa.xz", lzma.compress(_fa([("E_1", seq_of(9, 900))]))),  # xz, named so
    ("prot/A.faa", _fa([("p1", PROT)])),
    ("aln/x.fasta", _fa([("s1", ALN), ("s2", ALN)])),
    ("meta/table.tsv", b"genome\tcompleteness\nA\t95.1\n"),
    ("bins/empty.fa", b""),
]


def test_open_member_tells_compression_by_magic_bytes_not_by_name():
    mod = _mod()
    for raw, want in ((b">a\nACGT\n", "none"), (gzip.compress(b">a\nACGT\n"), "gzip"),
                      (bz2.compress(b">a\nACGT\n"), "bzip2"), (lzma.compress(b">a\nACGT\n"), "xz")):
        comp, fh = mod.open_member(io.BytesIO(raw))
        assert comp == want and fh.read() == b">a\nACGT\n"
    comp, fh = mod.open_member(io.BytesIO(b"\x28\xb5\x2f\xfd" + b"rest"))
    assert comp == "zstd"  # named, never read (no zstd reader in the standard library)


@pytest.mark.parametrize("outer", ["bz2", "gz"])
def test_census_names_every_member_with_its_compression_and_kind(tmp_path, outer):
    mod = _mod()
    tar = tmp_path / f"x.tar.{outer}"
    write_tar(tar, CENSUS_MEMBERS, outer=outer, dirs=["bins/"])
    out = tmp_path / "members.tsv"
    summary = mod.census(tar, out)
    rows = {r["member"]: r for r in _rows(out)}
    assert list(rows) == [n for n, _ in CENSUS_MEMBERS]  # archive order; the directory is counted, not listed
    comp = {n: r["compression"] for n, r in rows.items()}
    assert comp == {"bins/A.fa": "none", "bins/B.fa": "gzip", "bins/C.fna.bz2": "bzip2", "bins/D.fa.gz": "none",
                    "bins/E.fa.xz": "xz", "prot/A.faa": "none", "aln/x.fasta": "none", "meta/table.tsv": "none",
                    "bins/empty.fa": "none"}
    kind = {n: r["kind"] for n, r in rows.items()}
    assert kind == {"bins/A.fa": "dna", "bins/B.fa": "dna", "bins/C.fna.bz2": "dna", "bins/D.fa.gz": "dna",
                    "bins/E.fa.xz": "dna", "prot/A.faa": "protein", "aln/x.fasta": "dna_with_gaps",
                    "meta/table.tsv": "text", "bins/empty.fa": "empty"}
    assert (rows["bins/A.fa"]["records"], rows["bins/A.fa"]["bases"]) == ("2", str(len(A1) + len(A2)))
    assert (rows["bins/A.fa"]["min_record"], rows["bins/A.fa"]["max_record"]) == (str(len(A2)), str(len(A1)))
    assert rows["bins/B.fa"]["bases"] == str(len(B1)) and rows["aln/x.fasta"]["gap_bytes"] == "60"
    assert {n: r["stem"] for n, r in rows.items() if r["kind"] == "dna"} == {
        "bins/A.fa": "A", "bins/B.fa": "B", "bins/C.fna.bz2": "C", "bins/D.fa.gz": "D", "bins/E.fa.xz": "E"}
    assert summary["members"] == 9 and summary["directories"] == 1
    assert summary["by_kind"]["dna"] == {"members": 5, "records": 6,
                                         "bases": len(A1) + len(A2) + len(B1) + len(C1) + len(D1) + 900}
    assert summary["by_compression"] == {"bzip2": 1, "gzip": 1, "none": 6, "xz": 1}
    # deterministic
    out2 = tmp_path / "members2.tsv"
    mod.census(tar, out2)
    assert out.read_bytes() == out2.read_bytes()


@pytest.mark.parametrize("block", [5, 7, 64, 1 << 22])
def test_parse_fasta_in_blocks_agrees_with_biobox_fasta(tmp_path, monkeypatch, block):
    mod = _mod()
    from biobox import fasta

    recs = [("r1 first", seq_of(11, 130)), ("r2", seq_of(12, 1)), ("r3 x", seq_of(13, 77)), ("r4", seq_of(14, 200))]
    raw = _fa(recs, width=13).replace(b"\n", b"\r\n", 9)  # some CRLF line ends
    raw = raw.rstrip(b"\n")  # the last line without its newline
    p = tmp_path / "x.fa"
    p.write_bytes(raw)
    want = [(r.header.encode() if isinstance(r.header, str) else r.header, r.seq) for r in fasta.iter_records(p)]
    monkeypatch.setattr(mod, "BLOCK", block)
    got: list[list] = []
    mod.parse_fasta(io.BytesIO(raw), lambda h: got.append([h, b""]), lambda s: got[-1].__setitem__(1, got[-1][1] + s))
    assert [(h, s) for h, s in got] == want
    assert [h for h, _ in got] == [h.encode() for h, _ in recs]


def test_census_refuses_an_unreadable_archive_and_leaves_no_table(tmp_path):
    mod = _mod()
    bad = tmp_path / "bad.tar.bz2"
    bad.write_bytes(b"BZh9" + b"\x00" * 100)
    out = tmp_path / "members.tsv"
    with pytest.raises(mod.Refused, match="not a readable"):
        mod.census(bad, out)
    assert not out.exists()


def test_read_records_refuses_what_ingest_would_refuse():
    mod = _mod()
    assert mod.read_records(io.BytesIO(b">a x\nACGTN\nacg\n>b\nRYKM\n"), "m") == [("a x", b"ACGTNacg"), ("b", b"RYKM")]
    for raw, what in ((b">a\nAC-GT\n", "0x2d"), (b">a\nACGTX\n", "0x58"), (b">a\n>b\nACGT\n", "empty"),
                      (b"ACGT\n>a\nACGT\n", "before the first header"), (b">\xe5\nACGT\n", "not UTF-8")):
        with pytest.raises(mod.Refused, match=what):
            mod.read_records(io.BytesIO(raw), "m")


def test_quality_of_mag_uses_mimag_numbers():
    mod = _mod()
    assert [mod.quality_of_mag(c, k) for c, k in (("95.1", "1.2"), ("90.0", "1"), ("99", "5.0"), ("60", "3"),
                                                  ("NA", "1"))] == ["mag_hq", "mag_mq", "mag_mq", "mag_mq", "mag_mq"]


# ==================================================================================================================
# GBRMGD
# ==================================================================================================================

# Bins of two sites; contig names are the assembly's (`IMOS__P__<site>__contig_<n>_pilon:p [bin:...]`), which the
# stager prefixes with the bin's id. s2's second contig is held by V2 (an r26.11.4 stage), s3's by WET1 (an r26.11.5
# stage before GBRMGD); s4 repeats s1's contig set (DEC-27: listed once). The protein file, the table and a DNA file
# outside prokaryote_bins/ are left out.
G1A, G1B = seq_of(201, 2600), seq_of(202, 1400)
G2A, G2B = seq_of(203, 3000), seq_of(204, 1100)
G3A = seq_of(205, 2700)
GBR = {
    "hoskyn_2__rosella_refine_806_0": [("IMOS__P__hoskyn_2__contig_3485_pilon:p [group=P]", G1A),
                                       ("IMOS__P__hoskyn_2__contig_11_pilon:p [group=P]", G1B)],
    "erskine_4__metabat2_refined_691": [("IMOS__P__erskine_4__contig_16618_pilon:p", G2A),
                                        ("IMOS__P__erskine_4__contig_7_pilon:p", G2B)],
    "21580_2__metabat2_refined_129": [("IMOS__P__21580_2__contig_33182_pilon:p", G3A)],
    "hoskyn_2__semibin_9": [("IMOS__P__hoskyn_2__contig_3485_pilon:p again", G1A),
                            ("IMOS__P__hoskyn_2__contig_11_pilon:p again", G1B)],
}


def gbr_members() -> list[tuple[str, bytes]]:
    out = []
    for i, (g, recs) in enumerate(GBR.items()):
        raw = _fa(recs)
        out.append((f"prokaryote_bins/{g}.fna", gzip.compress(raw, mtime=0) if i == 1 else raw))
    out += [("prokaryote_bins/x.faa", _fa([("p", PROT)])), ("README.tsv", b"bin\tsite\n"),
            ("viruses/v1.fna", _fa([("v1", seq_of(206, 900))]))]
    return out


@pytest.fixture()
def gb(tmp_path: Path) -> dict:
    tar = tmp_path / "prokaryote_bins_5283.tar.bz2"
    write_tar(tar, gbr_members(), dirs=["prokaryote_bins/"])
    census = tmp_path / "members.tsv"
    _mod().census(tar, census)
    held = tmp_path / "chain" / "V2"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(G2B)]}]))
    later = tmp_path / "later" / "WET1"
    later.mkdir(parents=True)
    (later / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(G3A)]}]))
    return {"tar": tar, "census": census, "held": held, "later": later, "tmp": tmp_path}


def run_gbr(mod, w: dict, tag: str, **kw):
    t = w["tmp"]
    plan, left = mod.plan_gbrmgd(mod.read_census(w["census"]))
    return mod.stage_set(
        "GBRMGD", w["tar"], plan, t / f"meta_{tag}", t / f"fasta_{tag}", fasta_name="gbrmgd_bins.fna.gz",
        source=mod.GBRMGD_SOURCE, retrieved_at=mod.GBRMGD_RETRIEVED, status_head=mod.GBRMGD_STATUS_HEAD,
        record_extra="site", left_out=left, root=None, root_stages=(), owner_inputs={"V2": w["held"]},
        later_inputs={"WET1": w["later"]}, inputs_home=Path("/tmp/if2453/s2/inputs"), work=t / f"work_{tag}", **kw)


def test_gbrmgd_plan_takes_the_bins_and_names_what_it_leaves(gb):
    mod = _mod()
    plan, left = mod.plan_gbrmgd(mod.read_census(gb["census"]))
    assert sorted(v["genome"] for v in plan.values()) == sorted(GBR)
    assert {v["kind"] for v in plan.values()} == {"prokaryote"} and {v["quality"] for v in plan.values()} == {"mag_mq"}
    assert plan["prokaryote_bins/erskine_4__metabat2_refined_691.fna"]["status"] == (
        "erskine_4__metabat2_refined_691", "erskine_4", "metabat2_refined_691", "not_visible", "")
    assert {r[0]: r[-1] for r in left} == {
        "prokaryote_bins/x.faa": "protein FASTA (not genome sequence)", "README.tsv": "a table or text (not sequence)",
        "viruses/v1.fna": "DNA outside prokaryote_bins/ (not one of the set's prokaryote bins)"}


def test_gbrmgd_stage_prefixes_names_lists_repeats_once_and_is_deterministic(gb):
    mod = _mod()
    rep1 = run_gbr(mod, gb, "a")
    rep2 = run_gbr(mod, gb, "b")
    t = gb["tmp"]
    meta, fa = t / "meta_a" / "GBRMGD", t / "fasta_a" / "GBRMGD" / "gbrmgd_bins.fna.gz"
    genomes = json.loads((meta / "genomes.json").read_text())
    order = sorted(GBR)
    # one genome per bin in id order, less the bin that repeats an earlier one's contig set
    assert [g["contig_uids"] for g in genomes] == [[cuid(s) for _, s in GBR[g]] for g in order
                                                   if g != "hoskyn_2__semibin_9"]
    assert {(g["kind"], g["genetic_code"], g["assembly_quality"]) for g in genomes} == {("prokaryote", 11, "mag_mq")}
    assert [r["genome"] for r in _rows(meta / "repeated_records.tsv")] == ["hoskyn_2__semibin_9"]
    names = [ln[1:] for ln in _gz_lines(fa) if ln.startswith(">")]
    assert names == [f"{g}|{h}" for g in GBR for h, _ in GBR[g]]  # archive order, every name prefixed
    seqs = "".join(ln for ln in _gz_lines(fa) if not ln.startswith(">"))
    assert seqs == "".join(s for g in GBR for _, s in GBR[g])
    assert [(r["genome"], r["owner_stage"]) for r in _rows(meta / "shared_records.tsv")] == [
        ("erskine_4__metabat2_refined_691", "V2")]
    assert [(r["genome"], r["holder"]) for r in _rows(meta / "later_records.tsv")] == [
        ("21580_2__metabat2_refined_129", "WET1")]
    assert {r["genome"]: r["r232_status"] for r in _rows(meta / "r232_status.tsv")} == {g: "not_visible" for g in GBR}
    assert [r["member"] for r in _rows(meta / "left_out.tsv")] == ["prokaryote_bins/x.faa", "README.tsv",
                                                                    "viruses/v1.fna"]
    objs = json.loads((meta / "objects.json").read_text())
    assert objs[0]["path"] == "/tmp/if2453/s2/inputs/GBRMGD/gbrmgd_bins.fna.gz"
    assert objs[0]["redistributable"] is False and objs[0]["source_name"] == "tower_gbr_mgd"
    r = rep1
    assert (r["genomes_read"], r["genomes"], r["genomes_repeated"], r["contigs_distinct"]) == (4, 3, 1, 5)
    assert r["member_compression"] == {"gzip": 1, "none": 3}
    assert r["held_by"] == {"V2": 1} and r["later_held_by"] == {"WET1": 1}
    assert r["bases_new"] == len(G1A) + len(G1B) + len(G2A) + len(G3A)
    assert r["bases_new_after_later"] == len(G1A) + len(G1B) + len(G2A)
    for f in ("genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
              "later_records.tsv", "genetic_codes.tsv", "r232_status.tsv", "left_out.tsv"):
        assert (meta / f).read_bytes() == (t / "meta_b" / "GBRMGD" / f).read_bytes(), f
    assert fa.read_bytes() == (t / "fasta_b" / "GBRMGD" / "gbrmgd_bins.fna.gz").read_bytes()
    assert {k: v for k, v in rep1["sha256"].items() if k != "objects_staging.json"} == {
        k: v for k, v in rep2["sha256"].items() if k != "objects_staging.json"}
    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), "GBRMGD", str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


def test_stage_set_refuses_an_existing_output_and_writes_nothing_when_it_refuses(gb, tmp_path):
    mod = _mod()
    run_gbr(mod, gb, "e")
    with pytest.raises(mod.Refused, match="exists"):
        run_gbr(mod, gb, "e")
    bad = tmp_path / "bad.tar.bz2"
    write_tar(bad, [("prokaryote_bins/a__b.fna", _fa([("c1", seq_of(301, 900))]))])
    census = tmp_path / "bad.tsv"
    mod.census(bad, census)
    write_tar(bad, [("prokaryote_bins/a__b.fna", b">c1\nACGTXACGT\n")])  # the bytes changed after the census
    with pytest.raises(mod.Refused, match="0x58"):
        run_gbr(mod, {**gb, "tar": bad, "census": census}, "x")
    assert not (gb["tmp"] / "meta_x" / "GBRMGD").exists()
    d = gb["tmp"] / "fasta_x" / "GBRMGD"
    assert not d.exists() or list(d.iterdir()) == []


# ==================================================================================================================
# ANI helpers and the overlap with a stage staged later
# ==================================================================================================================

def test_genomes_where_picks_a_sets_mags_and_refuses_a_missing_column(tmp_path):
    mod = _mod()
    meta = tmp_path / "LR2"
    meta.mkdir()
    (meta / "records.tsv").write_text("genome_id\tcontig_uid\tassembly_quality\nm1\taa\tmag_hq\nm1\tbb\tmag_hq\n"
                                      "u1\tcc\tcontig\nm2\tdd\tmag_mq\n")
    assert mod.genomes_where(meta, "assembly_quality", {"mag_hq", "mag_mq"}) == {"m1", "m2"}
    with pytest.raises(mod.Refused, match="no column category"):
        mod.genomes_where(meta, "category", {"bin"})


def test_all_archaea_refs_writes_every_listed_genome_once(tmp_path):
    mod = _mod()
    a1, a2, b1 = seq_of(401, 900), seq_of(402, 700), seq_of(403, 800)
    fna = tmp_path / "arc.fna.gz"
    fna.write_bytes(gzip.compress(_fa([("GCA_1_c1", a1), ("GCA_1_c2", a2), ("GCA_2_c1", b1)]), mtime=0))
    lst = tmp_path / "list.tsv"
    lst.write_text("genome\tgtdb\tfirst_record\trecords\tbp\theaders_start_with_genome\ttaxonomy\n"
                   "GCA_1\ttrue\t0\t2\t1600\ttrue\td__Archaea;g__X\nGCA_2\ttrue\t2\t1\t800\ttrue\td__Archaea;g__Y\n")
    files = mod.all_archaea_refs(fna, lst, tmp_path / "out")
    assert sorted(files) == ["GCA_1", "GCA_2"]
    assert files["GCA_1"].read_bytes() == f">GCA_1_c1\n{a1}\n>GCA_1_c2\n{a2}\n".encode()


def test_overlap_counts_contigs_a_later_staged_set_also_lists(gb, tmp_path):
    mod = _mod()
    run_gbr(mod, gb, "o")
    meta = gb["tmp"] / "meta_o" / "GBRMGD"
    other = tmp_path / "RMD1"
    other.mkdir()
    (other / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(G1A), cuid(G3A), cuid(seq_of(9, 50))]}]))
    assert mod.overlap(meta, {"RMD1": other}, tmp_path / "w") == {"RMD1": {"contigs": 2,
                                                                          "bases": len(G1A) + len(G3A)}}


def test_ani_with_skani_when_the_image_has_it(gb, tmp_path):
    import shutil as _sh

    if _sh.which("skani") is None:
        pytest.skip("skani is not on PATH (GitHub Actions, or the tools-hidden PATH)")
    mod = _mod()
    big = {f"s{i}__b": [(f"c{i}", seq_of(500 + i, 60_000))] for i in range(3)}
    tar = tmp_path / "q.tar.bz2"
    write_tar(tar, [(f"prokaryote_bins/{g}.fna", _fa(r)) for g, r in big.items()])
    census = tmp_path / "q.tsv"
    mod.census(tar, census)
    w = {**gb, "tar": tar, "census": census}
    run_gbr(mod, w, "q")
    meta = gb["tmp"] / "meta_q" / "GBRMGD"
    fa = gb["tmp"] / "fasta_q" / "GBRMGD" / "gbrmgd_bins.fna.gz"
    ref = tmp_path / "REF"
    ref.mkdir()
    (ref / "records.tsv").write_text(f"genome\tcontig_uid\tcategory\nr0\t{cuid(big['s0__b'][0][1])}\tbin\n"
                                     f"r9\t{cuid(seq_of(999, 60_000))}\tvirus\n")
    rfa = tmp_path / "ref.fna"
    rfa.write_bytes(_fa([("r0c", big["s0__b"][0][1]), ("r9c", seq_of(999, 60_000))]))
    res = mod.ani(meta, [fa], tmp_path / "ani", refs={"REF": (ref, [rfa], ("category", ("bin",)))}, archaea=None,
                  threads=1)
    assert res["groups"]["REF"]["references"] == 1 and res["meet"] == 1
    assert [(r["genome"], r["stage"], r["ref_genome"]) for r in _rows(tmp_path / "ani" / "ani_hits.tsv")] == [
        ("s0__b", "REF", "r0")]


# ==================================================================================================================
# RMD1
# ==================================================================================================================

RM = {  # genome -> (data type, BioSample, GTDB-Tk taxonomy, completeness, contamination, records)
    "TARA_SAMEA6017985_MAG_00000015": ("TPAC-CORAL-MAGS", "SAMEA6017985", "d__Bacteria;p__Acidobacteriota;g__;s__",
                                       "98.59", "1.2", [("TARA_SAMEA6017985_MAG_00000015_000001", seq_of(601, 2600)),
                                                        ("TARA_SAMEA6017985_MAG_00000015_000002", seq_of(602, 1300))]),
    "CARD22-1_SAMN24292821_MAG_00000009": ("EXTERNAL-MAGS", "SAMN24292821", "d__Bacteria;p__Bacillota;g__Mycoplasma;s__",
                                           "61.0", "1.0", [("k141_77", seq_of(603, 2500))]),
    "ALME22-2_SAMEA112151714_REFG_10000001": ("EXTERNAL-ISOG", "SAMEA112151714",
                                              "d__Bacteria;p__Bacteroidota;g__Aquimarina;s__", "99.57", "0.95",
                                              [("NZ_X1.1", seq_of(604, 3000))]),
    "ASAH15-1_SAMN02947394_REFG_10000001": ("EXTERNAL-ISOG", "SAMN02947394",
                                            "d__Bacteria;p__Proteobacteria;g__Pseudoalteromonas;s__", "100", "0.51",
                                            [("NZ_Y1.1", seq_of(605, 2800))]),
}
RM_HEAD = ["data_type", "biosample", "genome", "GTDBTK_VERSION", "GTDBTK_DATABASE", "GTDBTK_TAXONOMY", "n_bgcs",
           "n_bgcs_complete", "primary_cluster", "secondary_cluster", "cluster_threshold",
           "secondary_cluster_representative", "QSCORE", "COMPLETENESS", "CONTAMINATION", "TOOL", "N50", "#SCAFFOLDS",
           "GC", "GENOME_SIZE"]


def write_rm_table(path: Path, extra: dict | None = None) -> None:
    lines = ["\t".join(RM_HEAD)]
    for g, (dt, bs, tax, comp, cont, recs) in {**RM, **(extra or {})}.items():
        vals = [dt, bs, g, "2.1.0", "R207", tax, "1", "0", "1", "1_0", "0.05", g, "50", comp, cont, "checkm", "1000",
                str(len(recs)), "0.5", str(sum(len(s) for _, s in recs))]
        lines.append("\t".join(f'"{v}"' if not v.replace(".", "").isdigit() else v for v in vals))
    path.write_text("\n".join(lines) + "\n")


def write_rm_gtdb(d: Path) -> None:
    d.mkdir()
    head = ["accession", "ncbi_biosample", "ncbi_wgs_master", "ncbi_bioproject", "ncbi_isolate",
            "ncbi_strain_identifiers"]
    bac = [["RS_GCF_000000301.1", "SAMEA112151714", "X", "P", "none", "none"],   # the Aquimarina isolate
           ["GB_GCA_000000302.1", "SAMN24292821", "Y", "P", "none", "none"],     # a MAG of the MAG's metagenome
           ["GB_GCA_000000303.1", "SAMN24292821", "Z", "P", "none", "none"]]
    for name, rows in (("ar53_metadata_r232.tsv.gz", []), ("bac120_metadata_r232.tsv.gz", bac)):
        (d / name).write_bytes(gzip.compress(("\n".join("\t".join(r) for r in [head, *rows]) + "\n").encode(),
                                             mtime=0))


def write_rm_index(d: Path) -> None:
    d.mkdir()
    pq.write_table(pa.table({
        "accession": ["GCF_000000301", "GCA_000000302", "GCA_000000303"],
        "genbank_accession": ["GCA_000000301", "GCA_000000302", "GCA_000000303"],
        "is_representative": [True, False, True], "representative": ["GCF_000000301", "GCA_000000303",
                                                                      "GCA_000000303"],
        "species": ["Aquimarina sp1", "Mycoplasma sp2", "Mycoplasma sp2"]}), d / "gtdb_r232.parquet")
    pq.write_table(pa.table({"globdb_id": ["GCF_000000301"], "insdc_accession": ["GCF_000000301"],
                             "original_id": [None], "cluster_id": [None]}), d / "globdb_r232.parquet")


@pytest.fixture()
def rm(tmp_path: Path) -> dict:
    members = [(f"RMDv1-genomes/{g}.fa.gz", gzip.compress(_fa(RM[g][5]), mtime=0)) for g in RM]
    members += [("RMDv1-genomes/NOT_IN_TABLE_MAG_1.fa", _fa([("x1", seq_of(606, 900))])),
                ("RMDv1-genomes/TARA_SAMEA6017985_MAG_00000015.faa", _fa([("p", PROT)]))]
    tar = tmp_path / "RMDv1-genomes.tar.gz"
    write_tar(tar, members, outer="gz", dirs=["RMDv1-genomes/"])
    census = tmp_path / "members.tsv"
    _mod().census(tar, census)
    table = tmp_path / "genomes-aggregated-summary.tsv"
    write_rm_table(table)
    g, i = tmp_path / "gtdb", tmp_path / "index"
    write_rm_gtdb(g)
    write_rm_index(i)
    held = tmp_path / "chain" / "V2"
    held.mkdir(parents=True)
    (held / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(RM["ALME22-2_SAMEA112151714_REFG_10000001"]
                                                                      [5][0][1])]}]))
    return {"tar": tar, "census": census, "table": table, "gtdb": g, "index": i, "held": held, "tmp": tmp_path}


def test_rmd1_table_is_read_quoted_and_refuses_an_unknown_data_type(rm, tmp_path):
    mod = _mod()
    t = mod.rmd1_table(rm["table"])
    assert sorted(t) == sorted(RM)
    assert t["CARD22-1_SAMN24292821_MAG_00000009"]["data_type"] == "EXTERNAL-MAGS"
    assert t["CARD22-1_SAMN24292821_MAG_00000009"]["GTDBTK_TAXONOMY"].endswith("g__Mycoplasma;s__")
    bad = tmp_path / "bad.tsv"
    write_rm_table(bad, {"X_EUK_1": ("EXTERNAL-EUK", "SAMN1", "", "0", "0", [("e", "ACGT")])})
    with pytest.raises(mod.Refused, match="EXTERNAL-EUK"):
        mod.rmd1_table(bad)


def test_rmd1_stage_records_r232_status_by_kind_and_is_deterministic(rm):
    mod = _mod()
    t = rm["tmp"]

    def run(tag):
        plan, left, info = mod.plan_rmd1(mod.read_census(rm["census"]), mod.rmd1_table(rm["table"]), rm["index"],
                                         rm["gtdb"])
        rep = mod.stage_set("RMD1", rm["tar"], plan, t / f"meta_{tag}", t / f"fasta_{tag}",
                            fasta_name="rmd1_reef_genomes.fna.gz", source=mod.RMD1_SOURCE,
                            retrieved_at=mod.RMD1_RETRIEVED, status_head=mod.RMD1_STATUS_HEAD,
                            record_extra="data_type", left_out=left, root=None, root_stages=(),
                            owner_inputs={"V2": rm["held"]}, later_inputs={},
                            inputs_home=Path("/tmp/if2453/s2/inputs"), work=t / f"work_{tag}")
        return rep, info

    (rep, info), (rep2, _) = run("a"), run("b")
    meta, fa = t / "meta_a" / "RMD1", t / "fasta_a" / "RMD1" / "rmd1_reef_genomes.fna.gz"
    st = {r["genome"]: r for r in _rows(meta / "r232_status.tsv")}
    assert {g: r["r232_status"] for g, r in st.items()} == {
        "TARA_SAMEA6017985_MAG_00000015": "not_visible", "CARD22-1_SAMN24292821_MAG_00000009": "not_visible",
        "ALME22-2_SAMEA112151714_REFG_10000001": "in_gtdb_and_globdb_r232",
        "ASAH15-1_SAMN02947394_REFG_10000001": "not_in_gtdb_r232"}
    assert st["ALME22-2_SAMEA112151714_REFG_10000001"]["gtdb_r232_accession"] == "GCF_000000301"
    assert st["CARD22-1_SAMN24292821_MAG_00000009"]["gtdb_r232_genomes_of_biosample"] == "2"
    assert {g: r["genome_kind"] for g, r in st.items()} == {g: ("isolate" if "REFG" in g else "mag") for g in RM}
    genomes = json.loads((meta / "genomes.json").read_text())
    order = sorted(RM)
    assert [g["contig_uids"] for g in genomes] == [[cuid(s) for _, s in RM[g][5]] for g in order]
    # in id order: two isolates (scaffold), the Mycoplasma MAG (61 % complete), the Tara MAG (98.6 % / 1.2 %)
    assert [g["assembly_quality"] for g in genomes] == ["scaffold", "scaffold", "mag_mq", "mag_hq"]
    assert [g["genetic_code"] for g in genomes] == [4 if "Mycoplasma" in RM[g][2] else 11 for g in order]
    names = [ln[1:] for ln in _gz_lines(fa) if ln.startswith(">")]
    assert names == [h if h.startswith(g + "_") else f"{g}|{h}" for g in RM for h, _ in RM[g][5]]
    assert {r["member"]: r["reason"] for r in _rows(meta / "left_out.tsv")} == {
        "RMDv1-genomes/NOT_IN_TABLE_MAG_1.fa": "DNA whose name is no genome of genomes-aggregated-summary.tsv",
        "RMDv1-genomes/TARA_SAMEA6017985_MAG_00000015.faa": "protein FASTA (not genome sequence)"}
    assert rep["member_compression"] == {"gzip": 4} and rep["held_by"] == {"V2": 1}
    assert info["isolates_found_by_biosample"] == 1 and info["mags_whose_biosample_gtdb_r232_holds"] == 1
    assert rep["r232_status"] == {"in_gtdb_and_globdb_r232": 1, "not_in_gtdb_r232": 1, "not_visible": 2}
    for f in ("genomes.json", "records.tsv", "r232_status.tsv", "genetic_codes.tsv", "left_out.tsv"):
        assert (meta / f).read_bytes() == (t / "meta_b" / "RMD1" / f).read_bytes(), f
    assert fa.read_bytes() == (t / "fasta_b" / "RMD1" / "rmd1_reef_genomes.fna.gz").read_bytes()
    chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), "RMD1", str(meta / "genomes.json"),
                          str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


# ==================================================================================================================
# the conf templates
# ==================================================================================================================

# RMD1's template is written when its staging ends; add ("RMD1", "rmd1_reef_genomes.fna.gz", "WET1") here with it.
@pytest.mark.parametrize("stage,fasta_name,after", [("GBRMGD", "gbrmgd_bins.fna.gz", "RMD1")])
def test_templates_pass_stage_est_and_carry_the_cards_values(stage, fasta_name, after):
    est = _load("stage_est", "stage_est.py")
    path = TEMPLATES / f"{stage}.conf"
    if not path.exists():
        pytest.fail(f"no template {path.name}")
    c = est.read_conf(path, stage, tree=REPO)
    assert (c.depth, c.crispr, c.library) == ("called", True, None)  # both stages bring well over 1.5 Gbp new
    assert c.after == (after,)  # DEC-85 (1): RMD1 after WET1, then GBRMGD
    assert str(c.inputs) == f"/tmp/if2453/s2/inputs/{stage}" and c.min_contig_len is None
    assert c.genecall_mode == "default"
    assert NOTE in c.notes
    assert {n for _, n in c.sums} == {"genomes.json", "objects.json", fasta_name}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    for h, _ in c.sums:
        assert h in note, f"SUM {h} not in the note"
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    allowed = ROOT13 | {s for s, _ in a26.R26_11_4_CHAIN}  # the template's scope; r26.11.5 holders go to A-29
    assert set(c.expect) <= allowed
    for owner, (n, _) in c.expect.items():
        assert n > 0


def test_stage_set_with_recycled_screen_workers_writes_the_same_bytes(gb, monkeypatch):
    mod = _mod()
    monkeypatch.setattr(mod, "TASKS_PER_CHILD", 1)  # every worker replaced after one genome
    run_gbr(mod, gb, "w1")
    run_gbr(mod, gb, "w2", workers=2)
    t = gb["tmp"]
    for f in ("genomes.json", "records.tsv", "genetic_codes.tsv", "r232_status.tsv"):
        assert (t / "meta_w1" / "GBRMGD" / f).read_bytes() == (t / "meta_w2" / "GBRMGD" / f).read_bytes(), f
    assert (t / "fasta_w1" / "GBRMGD" / "gbrmgd_bins.fna.gz").read_bytes() == (
        t / "fasta_w2" / "GBRMGD" / "gbrmgd_bins.fna.gz").read_bytes()

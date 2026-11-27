"""Card C-75 (DEC-79): GlobDB r226's rare phyla, unclassified domain and remaining archaea as GRARE1..3.

The truth of every fixture is written out by hand in this file: the selection rule by a brute-force filter over the
rows, the stage cut by sums worked out below, contig_uid by blake3 directly, the FASTA contents by the sequences the
test made. Input validation is `s2_check_inputs.py`. Each test runs twice, with the normal PATH and with PATH reduced
to `/usr/bin:/bin` (as tests/r1/test_r232_index.py does): nothing here needs a pod tool (bgzip, when absent, gives way
to Python's gzip, and the genetic-code screen is pyrodigal).
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
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-75")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-26_grare_inputs.md"


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
    return _load("c75_grare_inputs", "c75_grare_inputs.py")


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def cuid(seq: bytes) -> str:
    return blake3(seq.upper()).digest(length=16).hex()


# Sense codons without a stop codon (TAA, TAG, TGA)
_SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]


def synth(seed: int, *, n_genes: int = 115, codons: int = 300, uga_per_gene: int = 0) -> bytes:
    """A prokaryote-like contig of 80 + n_genes * (38 + 3 * codons) bp: RBS, ATG, sense codons, TAA. With
    uga_per_gene > 0 every ORF carries internal TGA codons (code 11 truncates them, codes 25 and 4 read through)."""
    rng = random.Random(seed)
    parts = ["A" * 40]
    for _ in range(n_genes):
        spacer = "".join(rng.choice("ACGT") for _ in range(7))
        cods = [rng.choice(_SENSE) for _ in range(codons - 2)]
        if uga_per_gene:
            step = max(1, len(cods) // (uga_per_gene + 1))
            for k in range(1, uga_per_gene + 1):
                cods[min(len(cods) - 1, k * step)] = "TGA"
        parts.append("".join(rng.choice("AT") for _ in range(25)) + "AGGAGG" + spacer + "ATG" + "".join(cods) + "TAA")
    parts.append("T" * 40)
    return "".join(parts).encode()


L = 80 + 115 * (38 + 3 * 300)  # the length synth() gives with its defaults


# ---- the tables ------------------------------------------------------------------------------------------------

NG_HEAD = ["ID", "stage", "position_in_stage", "resource", "genome_size", "contig_number", "checkm2_completeness",
           "checkm2_contamination", "taxonomy"]
ST_HEAD = ["ID", "domain", "phylum", "class", "order", "family", "genus", "species", "checkm2_completeness",
           "checkm2_contamination", "genome_size", "contig_number", "largest_contig", "contig_N50", "GC_content",
           "culture_collection"]


def tax(d: str, p: str, g: str = "g1") -> str:
    if d == "Unclassified":
        return ("d__Unclassified;p_Unclassified;c__Unclassified;o__Unclassified;f__Unclassified;g__Unclassified;"
                "s__Unclassified")
    return f"d__{d};p__{p};c__c1;o__o1;f__f1;g__{g};s__{g} sp"


def write_ng(path: Path, rows: list[tuple]) -> None:
    """rows: (ID, c41 stage, resource, genome_size, contig_number, completeness, contamination, taxonomy)"""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(NG_HEAD)
        for i, (gid, st, res, size, n, comp, cont, t) in enumerate(rows):
            w.writerow([gid, st, i, res, size, n, comp, cont, t])


def write_stats(path: Path, rows: list[tuple]) -> None:
    """rows: (ID, domain, phylum, genome_size, contig_number)"""
    buf = io.StringIO()
    buf.write("\t".join(ST_HEAD) + "\n")
    for gid, d, p, size, n in rows:
        buf.write("\t".join([gid, d, p, "c1", "o1", "f1", "g1", "s1", "95.0", "1.0", str(size), str(n), "1", "1",
                             "0.5", "0"]) + "\n")
    with gzip.open(path, "wt") as fh:
        fh.write(buf.getvalue())


def write_rekey(path: Path, rows: list[tuple]) -> None:
    path.write_text("r226_id\tstatus\tr232_id\treason\n" + "".join("\t".join(r) + "\n" for r in rows))


# ---- the selection rule ------------------------------------------------------------------------------------------

def selection_world(tmp: Path) -> dict:
    ng = [
        ("A_G1", "G1", "RBG", 100, 1, 90.0, 1.0, tax("Archaea", "Aenigmatarchaeota")),
        ("A_TP2", "G2", "mOTU", 100, 1, 95.0, 1.0, tax("Archaea", "Thermoproteota", "g2")),
        ("A_TP1", "G2", "RBG", 100, 1, 60.0, 2.0, tax("Archaea", "Thermoproteota", "g1")),
        ("U1", "G41", "SPIRE", 100, 1, 45.0, 1.0, tax("Unclassified", "")),
        ("ZB1", "G40", "SPIRE", 100, 1, 95.0, 1.0, tax("Bacteria", "ZeroBig")),
        ("ZB2", "G40", "SPIRE", 100, 1, 95.0, 1.0, tax("Bacteria", "ZeroBig")),
        ("ZB3", "G40", "SPIRE", 100, 1, 95.0, 12.0, tax("Bacteria", "ZeroBig")),
        ("ZP1", "G40", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "ZeroPhy")),
        ("T1", "G39", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Tinyota")),
        ("M2", "G38", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Minorota", "g1")),
        ("M1", "G38", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Minorota", "g1")),
        ("BP1", "G33", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Pseudomonadota")),
        ("BP2", "G33", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Pseudomonadota")),
        ("BP3", "G33", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Pseudomonadota")),
        ("PL1", "G29", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Planctomycetota")),
        ("PA1", "G27", "RBG", 100, 1, 95.0, 1.0, tax("Bacteria", "Patescibacteriota")),
    ]
    stats = [
        ("GCA_000000001", "Archaea", "Thermoproteota", 50, 1),
        ("GCF_000000002", "Bacteria", "Pseudomonadota", 60, 1),
        ("GCA_000000003", "Bacteria", "Minorota", 70, 1),
        ("GCA_000000004", "Bacteria", "Tinyota", 80, 1),
        ("A_G1", "Archaea", "Aenigmatarchaeota", 100, 1),
        ("A_TP1", "Archaea", "Thermoproteota", 100, 1),
        ("A_TP2", "Archaea", "Thermoproteota", 100, 1),
        ("U1", "Unclassified;p_Unclassified", "Unclassified", 100, 1),
    ]
    rekey = [
        ("A_TP1", "dropped", "", "dictionary_lists_it_but_not_in_r232"),
        ("A_TP2", "kept", "A_TP2", "same_original_id"),
        ("U1", "replaced", "GCA_000000099", "insdc_accession_in_gtdb_r232"),
        ("ZB1", "kept", "ZB1", "same_original_id"), ("ZB2", "kept", "ZB2", "same_original_id"),
        ("ZB3", "dropped", "", "resource_not_in_r232"), ("ZP1", "kept", "ZP1", "same_original_id"),
        ("T1", "kept", "T1", "same_id_no_dictionary"), ("M1", "kept", "M1", "same_original_id"),
        ("M2", "kept", "M2", "same_original_id"), ("BP1", "kept", "BP1", "same_original_id"),
        ("BP2", "kept", "BP2", "same_original_id"), ("BP3", "dropped", "", "not_in_r232"),
    ]
    p = {"ng": tmp / "non_gtdb_stages.tsv", "stats": tmp / "globdb_r226_tax_plus_stats.tsv.gz",
         "rekey": tmp / "rekey_unbuilt.tsv"}
    write_ng(p["ng"], ng)
    write_stats(p["stats"], stats)
    write_rekey(p["rekey"], rekey)
    return p


def _read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_the_selection_rule_and_its_parameters(tmp_path):
    m = _mod()
    w = selection_world(tmp_path)
    rep = m.select(w["ng"], w["stats"], w["rekey"], tmp_path / "sel", minor_max_genomes=3, n_stages=3,
                   max_stage_bp=600)
    rows = _read_tsv(tmp_path / "sel" / "selection.tsv")
    # Brute-force oracle: not G1, not a GART phylum; archaea, unclassified, a phylum with no GTDB genome (ZeroBig has 3
    # non-GTDB genomes, as many as minor_max, but no GTDB genome), or a bacterial phylum with fewer than 3 non-GTDB
    # genomes (Tinyota 1, Minorota 2; Pseudomonadota has 3 and is left out)
    assert {r["genome"] for r in rows} == {"A_TP1", "A_TP2", "U1", "ZB1", "ZB2", "ZB3", "ZP1", "T1", "M1", "M2"}
    by = {r["genome"]: r for r in rows}
    assert {g: by[g]["group"] for g in by} == {
        "A_TP1": "archaea", "A_TP2": "archaea", "U1": "unclassified_domain", "ZB1": "no_gtdb_phylum",
        "ZB2": "no_gtdb_phylum", "ZB3": "no_gtdb_phylum", "ZP1": "no_gtdb_phylum", "T1": "minor_phylum",
        "M1": "minor_phylum", "M2": "minor_phylum"}
    # order: archaea, unclassified, no-GTDB phyla by name, minor phyla by their non-GTDB count then name; inside a
    # phylum by (taxonomy, ID)
    assert [r["genome"] for r in rows] == ["A_TP1", "A_TP2", "U1", "ZB1", "ZB2", "ZB3", "ZP1", "T1", "M1", "M2"]
    # the cut: 1 + 3 * (bp before the phylum) // 1000: Thermoproteota 0, Unclassified 200, ZeroBig 300 -> GRARE1;
    # ZeroPhy 600 -> GRARE2; Tinyota 700, Minorota 800 -> GRARE3
    assert [r["stage"] for r in rows] == ["GRARE1"] * 6 + ["GRARE2"] + ["GRARE3"] * 3
    # the r232 key beside every genome
    assert (by["A_TP1"]["r232_status"], by["A_TP1"]["r232_id"], by["A_TP1"]["r232_reason"]) == (
        "dropped", "", "dictionary_lists_it_but_not_in_r232")
    assert (by["U1"]["r232_status"], by["U1"]["r232_id"]) == ("replaced_by GCA_000000099", "GCA_000000099")
    assert (by["M1"]["r232_status"], by["M1"]["r232_id"]) == ("kept", "M1")
    assert by["ZB3"]["assembly_quality"] == "contig" and by["ZB1"]["assembly_quality"] == "mag_hq"
    assert rep["parameters"]["minor_max_genomes"] == 3 and rep["parameters"]["n_stages"] == 3
    assert {s: (v["genomes"], v["bp"]) for s, v in rep["stages"].items()} == {
        "GRARE1": (6, 600), "GRARE2": (1, 100), "GRARE3": (3, 300)}
    assert rep["stages"]["GRARE1"]["r232_status"] == {"dropped": 2, "kept": 3, "replaced": 1}
    # the remainder of the unbuilt genomes, with its r232 key: the three Pseudomonadota genomes
    rem = _read_tsv(tmp_path / "sel" / "remainder.tsv")
    assert [r["genome"] for r in rem] == ["BP1", "BP2", "BP3"]
    assert rem[2]["r232_status"] == "dropped"
    assert rep["remainder"]["genomes"] == 3 and rep["remainder"]["bp"] == 300
    # the archaeal references: every archaeon of r226, GTDB and not, G1's too, from the stats table
    arch = _read_tsv(tmp_path / "sel" / "archaea.tsv")
    assert [(r["genome"], r["gtdb"], r["genome_size"]) for r in arch] == [
        ("A_G1", "false", "100"), ("A_TP1", "false", "100"), ("A_TP2", "false", "100"),
        ("GCA_000000001", "true", "50")]


def test_a_stage_over_the_cap_is_refused(tmp_path):
    m = _mod()
    w = selection_world(tmp_path)
    with pytest.raises(m.Refused, match="GRARE1"):
        m.select(w["ng"], w["stats"], w["rekey"], tmp_path / "sel", minor_max_genomes=3, n_stages=3,
                 max_stage_bp=500)
    assert not (tmp_path / "sel").exists()


def test_a_selected_genome_without_a_rekey_row_is_refused(tmp_path):
    m = _mod()
    w = selection_world(tmp_path)
    text = w["rekey"].read_text().replace("T1\tkept\tT1\tsame_id_no_dictionary\n", "")
    w["rekey"].write_text(text)
    with pytest.raises(m.Refused, match="T1"):
        m.select(w["ng"], w["stats"], w["rekey"], tmp_path / "sel", minor_max_genomes=3, n_stages=3,
                 max_stage_bp=600)


def test_the_command_line_selects_and_refuses(tmp_path):
    w = selection_world(tmp_path)
    script = BSP / "c75_grare_inputs.py"
    base = [sys.executable, str(script), "select", str(w["ng"]), str(w["stats"]), str(w["rekey"])]
    r = subprocess.run(base + [str(tmp_path / "sel"), "--minor-max-genomes", "3", "--max-stage-bp", "600"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert json.loads(r.stdout)["stages"]["GRARE2"]["genomes"] == 1
    r = subprocess.run(base + [str(tmp_path / "sel2"), "--minor-max-genomes", "3", "--max-stage-bp", "500"],
                       capture_output=True, text=True)
    assert r.returncode == 2 and "refused" in r.stderr


# ---- one stream, the code screen, and the stage files --------------------------------------------------------------

def _fa_gz(records: list[tuple[str, bytes]], width: int = 80) -> bytes:
    out = []
    for h, s in records:
        out.append(b">" + h.encode() + b"\n")
        out.extend(s[i:i + width] + b"\n" for i in range(0, len(s), width))
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0) as g:
        g.write(b"".join(out))
    return buf.getvalue()


def _tar(path: Path, members: list[tuple[str, list[tuple[str, bytes]]]]) -> None:
    with tarfile.open(path, "w:gz") as tf:
        info = tarfile.TarInfo("globdb_r226_genome_fasta/R226CHUNK000")
        info.type = tarfile.DIRTYPE
        tf.addfile(info)
        readme = b"not a genome\n"
        info = tarfile.TarInfo("globdb_r226_genome_fasta/README.txt")
        info.size = len(readme)
        tf.addfile(info, io.BytesIO(readme))
        for gid, recs in members:
            data = _fa_gz(recs)
            info = tarfile.TarInfo(f"globdb_r226_genome_fasta/R226CHUNK000/{gid}.fa.gz")
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))


def _hdr(gid: str, k: int) -> str:
    return f"{gid}_{k:012d}"


@pytest.fixture()
def world(tmp_path: Path) -> dict:
    s1 = synth(1)                      # A_TP1, also the second record of M1
    s2a, s2b = synth(2), synth(3)      # U1's two records; V2 (a chain stage) holds s2b
    s3 = synth(4, uga_per_gene=8)      # Z1, UGA read through
    s4 = synth(5)                      # M1's first record, held by the root's P2
    arch_g1 = synth(6, n_genes=3)      # a G1 archaeon: a reference only
    arch_gtdb = synth(7, n_genes=4)    # a GTDB archaeon: a reference only
    other = synth(8, n_genes=2)        # a bacterium nothing selects
    assert {len(s) for s in (s1, s2a, s2b, s3, s4)} == {L}
    x = 1000                           # A_TP0's table size: a selected archaeon the archive lacks
    ng = [
        ("A_G1", "G1", "RBG", len(arch_g1), 1, 90.0, 1.0, tax("Archaea", "Aenigmatarchaeota")),
        ("A_TP0", "G2", "RBG", x, 1, 95.0, 1.0, tax("Archaea", "Thermoproteota", "g0")),
        ("A_TP1", "G2", "RBG", L, 1, 95.0, 1.0, tax("Archaea", "Thermoproteota", "g1")),
        ("U1", "G41", "SPIRE", 2 * L, 2, 75.0, 2.0, tax("Unclassified", "")),
        ("Z1", "G40", "mOTU", L, 1, 92.0, 1.0, tax("Bacteria", "ZeroPhy")),
        ("M1", "G38", "mOTU", 2 * L, 2, 91.0, 0.5, tax("Bacteria", "Minorota")),
        ("BP1", "G33", "mOTU", len(other), 1, 95.0, 1.0, tax("Bacteria", "Pseudomonadota")),
        ("BP2", "G33", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Pseudomonadota")),
        ("BP3", "G33", "mOTU", 100, 1, 95.0, 1.0, tax("Bacteria", "Pseudomonadota")),
    ]
    stats = [
        ("GCA_000000001", "Archaea", "Thermoproteota", len(arch_gtdb), 1),
        ("GCF_000000002", "Bacteria", "Pseudomonadota", 60, 1),
        ("GCA_000000003", "Bacteria", "Minorota", 70, 1),
        ("A_G1", "Archaea", "Aenigmatarchaeota", len(arch_g1), 1),
        ("A_TP0", "Archaea", "Thermoproteota", x, 1),
        ("A_TP1", "Archaea", "Thermoproteota", L, 1),
    ]
    rekey = [("A_TP0", "kept", "A_TP0", "same_original_id"),
             ("A_TP1", "dropped", "", "dictionary_lists_it_but_not_in_r232"),
             ("U1", "replaced", "GCA_000000099", "insdc_accession_in_gtdb_r232"),
             ("Z1", "kept", "Z1", "same_original_id"), ("M1", "kept", "M1", "same_original_id"),
             ("BP1", "kept", "BP1", "same_original_id"), ("BP2", "kept", "BP2", "same_original_id"),
             ("BP3", "kept", "BP3", "same_original_id")]
    p = {"ng": tmp_path / "non_gtdb_stages.tsv", "stats": tmp_path / "stats.tsv.gz", "rekey": tmp_path / "rekey.tsv",
         "tar": tmp_path / "globdb_r226_genome_fasta.tar.gz"}
    write_ng(p["ng"], ng)
    write_stats(p["stats"], stats)
    write_rekey(p["rekey"], rekey)
    # archive order is not the selection order
    _tar(p["tar"], [
        ("M1", [(_hdr("M1", 1), s4), (_hdr("M1", 2), s1)]),
        ("BP1", [(_hdr("BP1", 1), other)]),
        ("GCA_000000001", [(_hdr("GCA_000000001", 1), arch_gtdb)]),
        ("Z1", [(_hdr("Z1", 1), s3)]),
        ("A_TP1", [(_hdr("A_TP1", 1), s1)]),
        ("U1", [(_hdr("U1", 1), s2a), (_hdr("U1", 2), s2b)]),
        ("A_G1", [(_hdr("A_G1", 1), arch_g1)]),
    ])
    root = tmp_path / "root"
    part = root / "l0seq" / "contig_index" / "stage=P2"
    part.mkdir(parents=True)
    pq.write_table(pa.table({"contig_uid": pa.array([bytes.fromhex(cuid(s4))], type=pa.binary(16)),
                             "is_default_variant": pa.array([True])}), part / "part-0.parquet")
    v2 = tmp_path / "chain" / "V2"
    v2.mkdir(parents=True)
    (v2 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s2b)], "genetic_code": 11, "kind": "virus",
                                                  "circular": [], "assembly_quality": "vOTU"}]))
    p.update(root=root, v2=v2, seqs={"s1": s1, "s2a": s2a, "s2b": s2b, "s3": s3, "s4": s4, "arch_g1": arch_g1,
                                     "arch_gtdb": arch_gtdb}, tmp=tmp_path)
    return p


def _run(m, w, name: str, *, workers: int = 1, cpu_guard_min: float = 280.0) -> dict:
    t = w["tmp"] / name
    sel = t / "sel"
    meta = t / "meta"
    fa = t / "fasta"
    m.select(w["ng"], w["stats"], w["rekey"], sel, minor_max_genomes=3, n_stages=3, max_stage_bp=4 * L)
    st = m.stream(sel, w["tar"], meta, fa, cpu_guard_min=cpu_guard_min)
    sc = m.screen(sel, meta, fa, workers=workers)
    fin = m.finish(sel, meta, fa, inputs_home=Path("/tmp/if2453/s2/inputs"), root=w["root"], root_stages=["P2"],
                   chain={"V2": w["v2"]}, work=t / "work")
    return {"sel": sel, "meta": meta, "fasta": fa, "stream": st, "screen": sc, "finish": fin}


def _fasta(path: Path) -> list[tuple[str, bytes]]:
    out, h, buf = [], None, []
    with gzip.open(path, "rb") as fh:
        for line in fh.read().splitlines():
            if line.startswith(b">"):
                if h is not None:
                    out.append((h, b"".join(buf)))
                h, buf = line[1:].decode(), []
            else:
                buf.append(line)
    if h is not None:
        out.append((h, b"".join(buf)))
    return out


def test_one_stream_writes_the_three_stages_and_the_archaeal_references(world):
    m = _mod()
    r = _run(m, world, "a")
    s = world["seqs"]
    fa = r["fasta"]
    # members taken by name, in archive order; the selected archaeon the archive lacks is reported, not invented
    assert _fasta(fa / "grare1.fna.gz") == [(_hdr("A_TP1", 1), s["s1"]), (_hdr("U1", 1), s["s2a"]),
                                            (_hdr("U1", 2), s["s2b"])]
    assert _fasta(fa / "grare2.fna.gz") == [(_hdr("Z1", 1), s["s3"])]
    assert _fasta(fa / "grare3.fna.gz") == [(_hdr("M1", 1), s["s4"]), (_hdr("M1", 2), s["s1"])]
    assert _fasta(fa / "globdb_r226_archaea.fna.gz") == [
        (_hdr("GCA_000000001", 1), s["arch_gtdb"]), (_hdr("A_TP1", 1), s["s1"]), (_hdr("A_G1", 1), s["arch_g1"])]
    st = r["stream"]
    assert st["streams"] == 1 and st["missing"] == {"GRARE1": ["A_TP0"], "archaea": ["A_TP0"]}
    assert st["members_seen"] == 7
    arch = _read_tsv(r["meta"] / "archaea" / "archaea_list.tsv")
    assert [(a["genome"], a["gtdb"], a["first_record"], a["records"], a["bp"]) for a in arch] == [
        ("GCA_000000001", "true", "0", "1", str(len(s["arch_gtdb"]))), ("A_TP1", "false", "1", "1", str(L)),
        ("A_G1", "false", "2", "1", str(len(s["arch_g1"])))]


def test_codes_records_ownership_and_s2_check(world):
    m = _mod()
    r = _run(m, world, "a")
    s = world["seqs"]
    meta = r["meta"]
    codes = {row["genome"]: row for st in ("GRARE1", "GRARE2", "GRARE3")
             for row in _read_tsv(meta / st / "genetic_codes.tsv")}
    assert set(codes) == {"A_TP1", "U1", "Z1", "M1"}
    assert codes["Z1"]["genetic_code"] in ("25", "4") and codes["Z1"]["tested_alt"] == "true"
    assert {codes[g]["genetic_code"] for g in ("A_TP1", "U1", "M1")} == {"11"}
    assert codes["U1"]["r232_status"] == "replaced_by GCA_000000099"
    rec1 = _read_tsv(meta / "GRARE1" / "records.tsv")
    assert [(x["genome"], x["record_ordinal"], x["contig_uid"], x["length"]) for x in rec1] == [
        ("A_TP1", "0", cuid(s["s1"]), str(L)), ("U1", "0", cuid(s["s2a"]), str(L)),
        ("U1", "1", cuid(s["s2b"]), str(L))]
    assert {x["r232_status"] for x in rec1 if x["genome"] == "A_TP1"} == {"dropped"}
    assert {x["group"] for x in rec1} == {"archaea", "unclassified_domain"}
    rep = r["finish"]
    assert rep["GRARE1"]["owner_contigs"] == {"V2": 1}
    assert rep["GRARE2"]["owner_contigs"] == {}
    # M1's first record is the root's (P2), its second is A_TP1's, which GRARE1 holds
    assert rep["GRARE3"]["owner_contigs"] == {"GRARE1": 1, "P2": 1}
    assert rep["GRARE3"]["new_contigs"] == 0 and rep["GRARE2"]["new_bases"] == L
    assert rep["GRARE1"]["missing_genomes"] == ["A_TP0"] and rep["GRARE1"]["genomes"] == 2
    assert rep["GRARE1"]["r232_status"] == {"dropped": 1, "replaced": 1}
    sh = _read_tsv(meta / "GRARE3" / "shared_records.tsv")
    assert sorted((x["genome"], x["contig_uid"], x["owner_stage"]) for x in sh) == sorted(
        [("M1", cuid(s["s4"]), "P2"), ("M1", cuid(s["s1"]), "GRARE1")])
    g3 = json.loads((meta / "GRARE3" / "genomes.json").read_text())
    assert g3 == [{"contig_uids": [cuid(s["s4"]), cuid(s["s1"])], "genetic_code": 11, "kind": "prokaryote",
                   "circular": [], "assembly_quality": "mag_hq"}]
    obj = json.loads((meta / "GRARE2" / "objects.json").read_text())
    assert obj[0]["path"] == "/tmp/if2453/s2/inputs/GRARE2/grare2.fna.gz"
    assert obj[0]["source_name"] == "globdb" and obj[0]["source_release"] == "r226"
    stg = json.loads((meta / "GRARE2" / "objects_staging.json").read_text())
    assert stg[0]["path"] == str(r["fasta"] / "grare2.fna.gz")
    for st in ("GRARE1", "GRARE2", "GRARE3"):
        assert rep[st]["s2_check"]["exit"] == 0 and rep[st]["s2_check"]["problems"] == [], rep[st]["s2_check"]


def test_every_file_is_byte_identical_when_made_again_and_with_more_workers(world):
    m = _mod()
    a = _run(m, world, "a")
    b = _run(m, world, "b", workers=2)
    names = ["genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
             "genetic_codes.tsv"]
    for st in ("GRARE1", "GRARE2", "GRARE3"):
        for n in names:
            assert (a["meta"] / st / n).read_bytes() == (b["meta"] / st / n).read_bytes(), (st, n)
    for f in ("grare1.fna.gz", "grare2.fna.gz", "grare3.fna.gz", "globdb_r226_archaea.fna.gz"):
        assert (a["fasta"] / f).read_bytes() == (b["fasta"] / f).read_bytes(), f
    for f in ("selection.tsv", "remainder.tsv", "archaea.tsv"):
        assert (a["sel"] / f).read_bytes() == (b["sel"] / f).read_bytes(), f
    assert (a["meta"] / "archaea" / "archaea_list.tsv").read_bytes() == (
        b["meta"] / "archaea" / "archaea_list.tsv").read_bytes()


def test_the_cpu_guard_stops_the_stream_and_leaves_no_file(world):
    m = _mod()
    t = world["tmp"] / "g"
    m.select(world["ng"], world["stats"], world["rekey"], t / "sel", minor_max_genomes=3, n_stages=3,
             max_stage_bp=4 * L)
    with pytest.raises(m.Refused, match="CPU"):
        m.stream(t / "sel", world["tar"], t / "meta", t / "fasta", cpu_guard_min=0.0)
    assert not list((t / "fasta").glob("*.fna.gz"))
    assert not (t / "meta" / "stream").exists()


def test_a_member_whose_records_differ_from_the_table_is_refused(world):
    m = _mod()
    t = world["tmp"] / "d"
    text = world["ng"].read_text().replace(f"Z1\tG40\t4\tmOTU\t{L}\t1", f"Z1\tG40\t4\tmOTU\t{L + 1}\t1")
    assert text != world["ng"].read_text()
    world["ng"].write_text(text)
    m.select(world["ng"], world["stats"], world["rekey"], t / "sel", minor_max_genomes=3, n_stages=3,
             max_stage_bp=4 * L)
    with pytest.raises(m.Refused, match="Z1"):
        m.stream(t / "sel", world["tar"], t / "meta", t / "fasta")


# ---- the conf templates -------------------------------------------------------------------------------------------

def test_templates_pass_stage_est_and_carry_the_cards_values():
    est = _load("stage_est", "stage_est.py")
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    allowed = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}
    allowed |= {s for s, _ in a26.R26_11_4_CHAIN}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    for n, after in ((1, "GART2"), (2, "GRARE1"), (3, "GRARE2")):
        st = f"GRARE{n}"
        c = est.read_conf(TEMPLATES / f"{st}.conf", st, tree=REPO)
        assert (c.depth, c.crispr, c.library, c.after, c.genecall_mode, c.min_contig_len) == (
            "called", True, None, (after,), "default", None)
        assert str(c.inputs) == f"/tmp/if2453/s2/inputs/{st}"
        assert NOTE in c.notes
        assert {name for _, name in c.sums} == {"genomes.json", "objects.json", f"grare{n}.fna.gz"}
        for h, _ in c.sums:
            assert h in note, f"{st}: SUM {h} not in the note"
        assert set(c.expect) <= allowed | {f"GRARE{k}" for k in range(1, n)}
        for _owner, (k, _) in c.expect.items():
            assert k > 0

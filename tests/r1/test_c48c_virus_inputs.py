"""Card C-48c (part 2): the selection rule and the staging of stage V2 (IMG/VR v4.1 high-confidence, high-quality vOTU
representatives, plus the Inovirus genomes), on small synthetic archives built in `tmp_path`.

Oracles are written here, separately from the script: contig uids are BLAKE3-128 of the upper-cased bases computed with
`blake3` directly, md5s with `hashlib`, and the representatives a brute-force sort of each vOTU's rows.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import io
import json
import random
import sys
import zipfile
from pathlib import Path

import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-48c")  # ring R1 by directory (tests/conftest.py)

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod"
HC = "Cus_VR/IMG_VR_2022-12-19_7.1/"
SEQ_TSV = HC + "IMGVR_all_Sequence_information-high_confidence.tsv"
FNA = HC + "IMGVR_all_nucleotides-high_confidence.fna.gz"
GB_ZIP = "Cusrus/data/Gb_files_inoviruses.zip"
HEAD = ["UVIG", "Taxon_oid", "Scaffold_oid", "Coordinates ('whole' if the UViG is the entire contig)",
        "Ecosystem classification", "vOTU", "Length", "Topology", "geNomad score", "Confidence",
        "Estimated completeness", "Estimated contamination", "MIUViG quality",
        "Gene content (total genes;cds;tRNA;geNomad marker)", "Taxonomic classification",
        "Taxonomic classification method", "Host taxonomy prediction", "Host prediction method",
        "Sequence origin (doi)"]


def _load():
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location("c48c_virus_inputs", SCRIPTS / "c48c_virus_inputs.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c48c_virus_inputs"] = mod
    spec.loader.exec_module(mod)
    return mod


vi = _load()


def cuid(seq: str) -> str:
    return blake3(seq.upper().encode()).digest(length=16).hex()


def row(uvig, votu, miuvig="High-quality", comp="100.00", cont="0.00", length=1000, topo="Linear",
        conf="High-confidence", eco="Environmental;Aquatic;Marine;Oceanic", host="d__Bacteria;p__X"):
    return {"UVIG": uvig, "vOTU": votu, "MIUViG quality": miuvig, "Estimated completeness": comp,
            "Estimated contamination": cont, "Length": str(length), "Topology": topo, "Confidence": conf,
            "Ecosystem classification": eco, "Host taxonomy prediction": host, "Host prediction method": "k-mer match",
            "Sequence origin (doi)": "geNomad v1.1.0 (10.5281/zenodo.7015982)"}


def brute(rows):
    """The rule written out by hand: per vOTU, its High-quality rows sorted by completeness capped at 100 (high
    first), contamination (low first), length (long first), UViG id; the first one."""
    out = {}
    for v in sorted({r["vOTU"] for r in rows}):
        hq = [r for r in rows if r["vOTU"] == v and r["MIUViG quality"] == "High-quality"]
        if hq:
            hq.sort(key=lambda r: (-min(float(r["Estimated completeness"]), 100.0),
                                   float(r["Estimated contamination"]), -int(r["Length"]), r["UVIG"]))
            out[v] = hq[0]["UVIG"]
    return out


RULE_ROWS = [
    # vOTU A: three HQ at 100% complete; the least contaminated wins over the longest
    row("U_A1", "vOTU_A", comp="100.00", cont="2.00", length=9000),
    row("U_A2", "vOTU_A", comp="100.00", cont="0.00", length=5000),
    row("U_A3", "vOTU_A", comp="95.00", cont="0.00", length=20000),
    row("U_A4", "vOTU_A", miuvig="Genome fragment", comp="50.00", length=30000),
    # vOTU B: fragments only, not taken
    row("U_B1", "vOTU_B", miuvig="Genome fragment", comp="80.00"),
    row("U_B2", "vOTU_B", miuvig="Unsure", comp="NA", cont="NA"),
    # vOTU C: a RefSeq genome only, not taken (MIUViG 'Reference' is not a quality tier)
    row("U_C1", "vOTU_C", miuvig="Reference", comp="100.00"),
    # vOTU D: 130% complete is no better than 100%: then contamination, then length, then the id
    row("U_D2", "vOTU_D", comp="130.00", cont="0.00", length=7000),
    row("U_D1", "vOTU_D", comp="100.00", cont="0.00", length=7000),
    row("U_D3", "vOTU_D", comp="99.00", cont="0.00", length=9000),
    # vOTU E: tie on everything but length
    row("U_E1", "vOTU_E", length=4000),
    row("U_E2", "vOTU_E", length=4100),
]


def test_imgvr_representative_rule_matches_brute_force_in_any_order():
    want = brute(RULE_ROWS)
    assert want == {"vOTU_A": "U_A2", "vOTU_D": "U_D1", "vOTU_E": "U_E2"}  # the rule, as the card states it
    rng = random.Random(4823)
    for _ in range(20):
        rows = RULE_ROWS[:]
        rng.shuffle(rows)
        reps, stats = vi.imgvr_representatives(rows)
        assert {v: r["UVIG"] for v, r in reps.items()} == want
        assert stats["rows"] == len(RULE_ROWS)
        assert stats["hq_rows"] == 8
        assert stats["hq_votus"] == 3


@pytest.mark.parametrize("bad", [row("U_X", "vOTU_X", conf="Putative"),
                                 row("U_Y", "vOTU_Y", comp="NA"),
                                 row("U_Z", "vOTU_Z", cont="")])
def test_imgvr_representative_rule_refuses_rows_it_cannot_judge(bad):
    with pytest.raises(vi.Refused):
        vi.imgvr_representatives(RULE_ROWS + [bad])


def test_circular_only_for_a_terminal_repeat_already_trimmed():
    unit = seq_of(11, 1200)  # random bases: the first 20 are not found again near the end
    with_repeat = unit + unit[:25]  # DEC-65: exact 20..1000 bp DTR is trimmed and staged circular
    with_long_repeat = unit + unit[:1001]  # DEC-65: > 1000 bp DTR stays linear and untrimmed
    assert vi.is_circular("Direct terminal repeat", unit.encode()) is True
    assert vi.is_circular("Direct terminal repeat", with_repeat.encode()) is True
    assert vi.is_circular("Direct terminal repeat", with_long_repeat.encode()) is False
    assert vi.is_circular("circular", unit.encode()) is True  # GenBank's word for the Inovirus records
    assert vi.is_circular("circular", with_repeat.encode()) is True  # DEC-65: Inovirus never trimmed, always circular
    for topo in ("Linear", "Provirus", "Inverted terminal repeat", "GVMAG", "linear", ""):
        assert vi.is_circular(topo, unit.encode()) is False


GB = """LOCUS       {name}            {n} bp    DNA     {topo}   PHG 01-JAN-2019
DEFINITION  Filamentous phage genome (fragment) identified in {name}
            second line of the definition.
ACCESSION   {name}
FEATURES             Location/Qualifiers
     source          1..{n}
ORIGIN
{origin}
//
"""


def genbank(name: str, seq: str, topo: str = "linear") -> str:
    lines = []
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60].lower()
        lines.append(f"{i + 1:>9} " + " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10)))
    return GB.format(name=name, n=len(seq), topo=topo, origin="\n".join(lines))


def test_genbank_records_parse_with_their_sequence_and_topology():
    s1 = "ACGT" * 40 + "A"
    s2 = "GGCCTTAA" * 9
    recs = list(vi.iter_genbank(io.BytesIO((genbank("Sp_one", s1) + genbank("Sp_two", s2, "circular")).encode())))
    assert [r["locus"] for r in recs] == ["Sp_one", "Sp_two"]
    assert [r["seq"] for r in recs] == [s1.lower().encode(), s2.lower().encode()]
    assert [r["topology"] for r in recs] == ["linear", "circular"]
    assert recs[0]["locus_length"] == len(s1)
    assert recs[0]["definition"].endswith("second line of the definition.")
    with pytest.raises(ValueError):
        list(vi.iter_genbank(io.BytesIO(genbank("Sp_one", s1).replace("//\n", "").encode())))


# --- the whole staging, on synthetic archives ------------------------------------------------------------------------

def seq_of(k: int, n: int) -> str:
    rng = random.Random(k)
    return "".join(rng.choice("ACGT") for _ in range(n))


def build_world(tmp: Path):
    """An HC zip (the table deflated, the FASTA stored gzip as JGI's), a metadata zip holding the Inovirus GenBank zip,
    and V1's records.tsv; returns the paths and what the oracle expects."""
    unit = seq_of(1, 3000)
    seqs = {
        "IMGVR_UViG_1_000001": unit,                    # vOTU_1 rep, DTR, repeat trimmed: circular
        "IMGVR_UViG_1_000002": seq_of(2, 2500),         # vOTU_1, shorter, not chosen
        "IMGVR_UViG_2_000001": seq_of(3, 4000) + seq_of(3, 4000)[:30],  # vOTU_2 rep, DTR still carrying its repeat
        "IMGVR_UViG_3_000001": seq_of(4, 3500),         # vOTU_3 fragment only: not taken
        "IMGVR_UViG_4_000001": seq_of(5, 5000),         # vOTU_4 rep, provirus; V1 holds the same bytes
        "IMGVR_UViG_GVMAG-1_000001": None,              # vOTU_5 rep, a bin of two records
    }
    gv_parts = [seq_of(6, 6000), seq_of(7, 4000)]
    rows = [
        row("IMGVR_UViG_1_000001", "vOTU_1", length=3000, topo="Direct terminal repeat"),
        row("IMGVR_UViG_1_000002", "vOTU_1", length=2500, cont="1.00"),
        row("IMGVR_UViG_2_000001", "vOTU_2", length=4030, topo="Direct terminal repeat"),
        row("IMGVR_UViG_3_000001", "vOTU_3", miuvig="Genome fragment", comp="40.00", length=3500),
        row("IMGVR_UViG_4_000001", "vOTU_4", length=5000, topo="Provirus"),
        row("IMGVR_UViG_GVMAG-1_000001", "vOTU_5", length=10000, topo="GVMAG"),
    ]
    tsv = "\t".join(HEAD) + "\n" + "".join("\t".join(r.get(h, "x") for h in HEAD) + "\n" for r in rows)
    fa = io.StringIO()
    for name, s in seqs.items():
        if s is None:
            for j, p in enumerate(gv_parts):
                fa.write(f">{name}|GVMAG-1|contig{j}\n" + "\n".join(p[i:i + 70] for i in range(0, len(p), 70)) + "\n")
            continue
        fa.write(f">{name}|1|2|1-{len(s)}\n" + "\n".join(s[i:i + 70] for i in range(0, len(s), 70)) + "\n")
    fna_gz = gzip.compress(fa.getvalue().encode(), mtime=0)
    hc = tmp / "hc.zip"
    with zipfile.ZipFile(hc, "w") as z:
        z.writestr(zipfile.ZipInfo(SEQ_TSV), tsv.encode(), compress_type=zipfile.ZIP_DEFLATED)
        z.writestr(zipfile.ZipInfo(FNA), fna_gz, compress_type=zipfile.ZIP_STORED)
    ino = {"Protoinoviridae/Protoinoviridae_Sf_1/Sp_b.gb": genbank("Sp_b", seq_of(8, 6000), "circular"),
           "Protoinoviridae/Protoinoviridae_Sf_1/Sp_a.gb": genbank("Sp_a", seq_of(9, 7000)),
           "Tandem/Tandem_Sf_2/Sp_c.gb": genbank("Sp_c", seq_of(5, 5000))}  # same bytes as IMGVR_UViG_4
    gb_buf = io.BytesIO()
    with zipfile.ZipFile(gb_buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in ino.items():
            z.writestr(zipfile.ZipInfo(name), text.encode())
    meta = tmp / "meta.zip"
    with zipfile.ZipFile(meta, "w") as z:
        z.writestr(zipfile.ZipInfo(GB_ZIP), gb_buf.getvalue(), compress_type=zipfile.ZIP_DEFLATED)
    v1 = tmp / "V1"
    v1.mkdir()
    (v1 / "records.tsv").write_text("ordinal\tgenome\tcontig_uid\tlength\n" + f"0\tuhgv-x\t{cuid(seq_of(5, 5000))}\t5000\n")
    md5 = {"seq_tsv": hashlib.md5(tsv.encode()).hexdigest(), "fna": hashlib.md5(fna_gz).hexdigest()}
    want = {
        "genomes": [[cuid(seqs["IMGVR_UViG_1_000001"])], [cuid(seq_of(3, 4000))],  # DEC-65: 30 bp exact DTR trimmed
                    [cuid(seqs["IMGVR_UViG_4_000001"])], [cuid(p) for p in gv_parts],
                    [cuid(seq_of(9, 7000))], [cuid(seq_of(8, 6000))]],
        "circular": [cuid(seqs["IMGVR_UViG_1_000001"]), cuid(seq_of(3, 4000)), cuid(seq_of(8, 6000))],
        "shared": [cuid(seqs["IMGVR_UViG_4_000001"])],
        "repeated": ["Sp_c"],
        "fasta_records": ["IMGVR_UViG_1_000001|1|2|1-3000", "IMGVR_UViG_2_000001|1|2|1-4030",
                          "IMGVR_UViG_4_000001|1|2|1-5000", "IMGVR_UViG_GVMAG-1_000001|GVMAG-1|contig0",
                          "IMGVR_UViG_GVMAG-1_000001|GVMAG-1|contig1"],
    }
    return hc, meta, v1, md5, want


def run_v2(hc, meta, v1, md5, out):
    return vi.stage_v2(hc, meta, out, root=None, root_stages=(), owner_inputs={"V1": v1},
                       expected_md5=md5, retrieved_imgvr="2026-11-18T23:15:37Z",
                       retrieved_inovirus="2026-11-18T23:16:08Z")


def test_stage_v2_inputs_match_the_oracle_and_repeat_byte_for_byte(tmp_path):
    hc, meta, v1, md5, want = build_world(tmp_path)
    out = tmp_path / "V2"
    rep = run_v2(hc, meta, v1, md5, out)
    genomes = json.loads((out / "genomes.json").read_text())
    assert [g["contig_uids"] for g in genomes] == want["genomes"]  # contigs in record order, genomes in FASTA order
    assert sorted(u for g in genomes for u in g["circular"]) == sorted(want["circular"])
    assert {g["kind"] for g in genomes} == {"virus"}
    assert {g["assembly_quality"] for g in genomes} == {"vOTU"}
    assert {g["genetic_code"] for g in genomes} == {11}
    shared = (out / "shared_records.tsv").read_text().splitlines()
    assert [ln.split("\t")[2] for ln in shared[1:]] == want["shared"] * 2  # the IMG/VR record and the Inovirus one
    assert {ln.split("\t")[3] for ln in shared[1:]} == {"V1"}
    repeated = (out / "repeated_records.tsv").read_text().splitlines()
    assert [ln.split("\t")[1] for ln in repeated[1:]] == want["repeated"]
    objects = json.loads((out / "objects.json").read_text())
    assert [o["source_name"] for o in objects] == ["jgi_img_vr", "jgi_inovirus"]
    assert all(o["redistributable"] is False for o in objects)
    assert all(Path(o["path"]).parent == out for o in objects)
    heads = [ln[1:] for ln in Path(objects[0]["path"]).read_text().splitlines() if ln.startswith(">")]
    assert heads == want["fasta_records"]
    ino_heads = [ln[1:].split()[0] for ln in Path(objects[1]["path"]).read_text().splitlines() if ln.startswith(">")]
    assert ino_heads == ["Sp_a", "Sp_b", "Sp_c"]  # file name order
    assert rep["hq_votus"] == 4 and rep["genomes"] == 6 and rep["imgvr_records"] == 5
    assert rep["md5"]["fna"]["equal"] and rep["md5"]["seq_tsv"]["equal"]
    assert not (tmp_path / "V2.tmp").exists()

    again = tmp_path / "again" / "V2"
    again.parent.mkdir()
    run_v2(hc, meta, v1, md5, again)
    for f in sorted(p.name for p in out.iterdir()):
        a, b = (out / f).read_bytes(), (again / f).read_bytes()
        if f == "objects.json":
            b = b.replace(str(again).encode(), str(out).encode())
        assert a == b, f


def test_stage_v2_refuses_a_member_whose_md5_differs_and_writes_nothing(tmp_path):
    hc, meta, v1, md5, _ = build_world(tmp_path)
    out = tmp_path / "V2"
    with pytest.raises(vi.Refused, match="md5"):
        run_v2(hc, meta, v1, {**md5, "fna": "0" * 32}, out)
    assert not out.exists()
    assert not (tmp_path / "V2.tmp").exists()


def test_stage_v2_refuses_an_existing_output(tmp_path):
    hc, meta, v1, md5, _ = build_world(tmp_path)
    out = tmp_path / "V2"
    out.mkdir()
    with pytest.raises(vi.Refused, match="exists"):
        run_v2(hc, meta, v1, md5, out)


# --- card C-48c part 3: stage GV1, every genome of GVMAGs V2's all_18kgenomes_fna.tar.gz ------------------------------

GV_MEMBER = "gvmags_v2/all_18kgenomes_fna.tar.gz"


def _tar_bytes(members: dict[str, bytes], mode: str) -> bytes:
    import tarfile

    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode=mode, format=tarfile.GNU_FORMAT) as tf:
        for name, data in members.items():
            ti = tarfile.TarInfo(name)
            ti.size, ti.mtime = len(data), 0
            tf.addfile(ti, io.BytesIO(data))
    return buf.getvalue()


def fna(genome: str, seqs: list[str]) -> bytes:
    return "".join(f">{genome}|c{j}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n"
                   for j, s in enumerate(seqs)).encode()


def build_gv(tmp: Path):
    """A GVMAGs-like tar: genome A of three contigs (one of them held by V2's inputs), B both of whose contigs V2 holds,
    C with A's bytes under another name, D sharing one contig with A, and a non-FASTA member."""
    a = [seq_of(21, 3000), seq_of(22, 2000), seq_of(23, 1500)]
    b = [seq_of(24, 2500), seq_of(25, 1800)]
    d = [seq_of(22, 2000), seq_of(26, 900)]  # A's second contig, and one of its own
    inner = _tar_bytes({"confirmed/B.fna": fna("B", b), "confirmed/A.fna": fna("A", a),
                        "confirmed/chunk_maker.sh~": b"#!/bin/bash\necho not data\n",
                        "confirmed/C.fna": fna("C", a), "confirmed/D.fna": fna("D", d)}, "w:gz")
    tar = tmp / "gvmags.tar"
    tar.write_bytes(_tar_bytes({GV_MEMBER: inner, "gvmags_v2/long_ortho_8000.tsv": b"ortho\tgenome\tgeneid\n"}, "w"))
    v2 = tmp / "V2"
    v2.mkdir()
    (v2 / "records.tsv").write_text("ordinal\tsource\tgenome\tcontig_uid\tlength\n" + "".join(
        f"{k}\timgvr\tU{k}\t{cuid(s)}\t{len(s)}\n" for k, s in enumerate([a[2], *b])))
    want = {"genomes": [[cuid(s) for s in b], [cuid(s) for s in a], [cuid(s) for s in d]],  # archive order, C repeats A
            "shared": sorted([cuid(a[2])] * 2 + [cuid(s) for s in b]), "repeated": ["C"],
            "several": [cuid(a[1])], "md5": hashlib.md5(inner).hexdigest(),
            "bases": sum(map(len, a + b + d)) - len(a[1])}
    return tar, v2, want


def run_gv1(tar, v2, out, md5):
    return vi.stage_gv1(tar, out, root=None, root_stages=(), owner_inputs={"V2": v2}, expected_md5=md5,
                        retrieved="2026-11-19T10:42:14Z")


def test_stage_gv1_lists_every_genome_and_repeats_byte_for_byte(tmp_path):
    tar, v2, want = build_gv(tmp_path)
    out = tmp_path / "GV1"
    rep = run_gv1(tar, v2, out, want["md5"])
    genomes = json.loads((out / "genomes.json").read_text())
    assert [g["contig_uids"] for g in genomes] == want["genomes"]
    assert {(g["kind"], g["genetic_code"], g["assembly_quality"]) for g in genomes} == {("virus", 11, "contig")}
    assert all(g["circular"] == [] for g in genomes)
    shared = [ln.split("\t") for ln in (out / "shared_records.tsv").read_text().splitlines()[1:]]
    assert sorted(r[3] for r in shared) == want["shared"]  # A's held contig in A and in C, and both of B's
    assert {r[4] for r in shared} == {"V2"}
    repeated = [ln.split("\t") for ln in (out / "repeated_records.tsv").read_text().splitlines()[1:]]
    assert [r[1] for r in repeated] == want["repeated"] and repeated[0][3] == "A"
    skipped = [ln.split("\t") for ln in (out / "skipped_members.tsv").read_text().splitlines()[1:]]
    assert [(r[0], r[2]) for r in skipped] == [("confirmed/chunk_maker.sh~", "not_fasta")]
    assert rep["genomes_all_contigs_held_earlier"] == 1  # B (C repeats A, and A has contigs of its own)
    assert rep["contigs_in_several_genomes"] == 1 and rep["fna_members"] == 4
    assert rep["fna_members"] == rep["genomes"] + rep["repeated_genomes"]
    assert rep["new_bases"] == want["bases"] - sum(len(seq_of(k, n)) for k, n in ((23, 1500), (24, 2500), (25, 1800)))
    objects = json.loads((out / "objects.json").read_text())
    assert [(o["source_name"], o["redistributable"]) for o in objects] == [("gvmags_v2", True)]
    heads = [ln[1:] for ln in Path(objects[0]["path"]).read_text().splitlines() if ln.startswith(">")]
    assert heads == ["B|c0", "B|c1", "A|c0", "A|c1", "A|c2", "C|c0", "C|c1", "C|c2", "D|c0", "D|c1"]
    again = tmp_path / "again" / "GV1"
    again.parent.mkdir()
    run_gv1(tar, v2, again, want["md5"])
    for f in sorted(p.name for p in out.iterdir()):
        a, b = (out / f).read_bytes(), (again / f).read_bytes()
        if f == "objects.json":
            b = b.replace(str(again).encode(), str(out).encode())
        assert a == b, f


def test_stage_gv1_refuses_a_member_whose_md5_differs(tmp_path):
    tar, v2, want = build_gv(tmp_path)
    out = tmp_path / "GV1"
    with pytest.raises(vi.Refused, match="md5"):
        run_gv1(tar, v2, out, "0" * 32)
    assert not out.exists() and not (tmp_path / "GV1.tmp").exists()

"""Card C-53 (DEC-65): terminal-repeat classification, trimming, and circularity for stage V2
(IMG/VR v4.1 HC vOTU representatives + Inovirus Roux 2019 GenBank records).

Oracles are written directly in this test (independent BLAKE3-128 digests and explicit constructed sequences for
every DEC-65 class).
"""

from __future__ import annotations

import csv
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

pytestmark = pytest.mark.covers("C-53", "DEC-65")  # ring R1 by directory (tests/conftest.py)

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod"
HC = "Cus_VR/IMG_VR_2022-12-19_7.1/"
SEQ_TSV = HC + "IMGVR_all_Sequence_information-high_confidence.tsv"
FNA = HC + "IMGVR_all_nucleotides-high_confidence.fna.gz"
GB_ZIP = "Cusrus/data/Gb_files_inoviruses.zip"
HEAD = [
    "UVIG", "Taxon_oid", "Scaffold_oid", "Coordinates ('whole' if the UViG is the entire contig)",
    "Ecosystem classification", "vOTU", "Length", "Topology", "geNomad score", "Confidence",
    "Estimated completeness", "Estimated contamination", "MIUViG quality",
    "Gene content (total genes;cds;tRNA;geNomad marker)", "Taxonomic classification",
    "Taxonomic classification method", "Host taxonomy prediction", "Host prediction method",
    "Sequence origin (doi)",
]


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


def seq_of(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s: str, pos: int) -> str:
    other = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[pos]]
    return s[:pos] + other + s[pos + 1:]


def row(uvig: str, votu: str, length: int, topo: str) -> dict[str, str]:
    return {
        "UVIG": uvig, "vOTU": votu, "MIUViG quality": "High-quality", "Estimated completeness": "100.00",
        "Estimated contamination": "0.00", "Length": str(length), "Topology": topo, "Confidence": "High-confidence",
        "Ecosystem classification": "Environmental;Aquatic;Marine", "Host taxonomy prediction": "d__Bacteria;p__X",
        "Host prediction method": "k-mer match", "Sequence origin (doi)": "geNomad v1.1.0 (10.5281/zenodo.7015982)",
    }


GB = """LOCUS       {name}            {n} bp    DNA     {topo}   PHG 01-JAN-2019
DEFINITION  Filamentous phage genome identified in {name}.
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


def test_ends_and_classify_cover_all_dec65_classes():
    # exact at lower bound (20), typical (75), upper bound (1000)
    for seed, k in ((101, 20), (102, 75), (103, 1000)):
        body = seq_of(seed, 3000)
        e = vi.ends((body + body[:k]).encode())
        assert e.exact == k
        assert vi.classify(e) == "exact"
        assert vi.trim((body + body[:k]).encode(), e.exact) == body.encode()

    # none: no exact repeat >= 20 bp and no near agreement
    body = seq_of(104, 3000)
    e = vi.ends(body.encode())
    assert e.exact == 0 and e.near == 0
    assert vi.classify(e) == "none"

    # partial: 60 bp near repeat with 1 substitution at index 35
    body = seq_of(105, 3000)
    rep = mutate(body[:60], 35)
    e = vi.ends((body + rep).encode())
    assert e.exact == 0 and e.near >= 1
    assert vi.classify(e) == "partial"

    # long: exact repeat of 1001 bp (> 1000)
    body = seq_of(106, 4000)
    e = vi.ends((body + body[:1001]).encode())
    assert e.exact == 1001
    assert vi.classify(e) == "long"

    # low_complexity: period <= 6, homopolymer/80% one base, or contains N
    mid = "GGTT" + seq_of(107, 2000) + "TTGG"
    for bad_rep in ("AC" * 20, "A" * 30, seq_of(108, 15) + "N" + seq_of(109, 15)):
        e = vi.ends((bad_rep + mid + bad_rep).encode())
        assert e.exact == len(bad_rep)
        assert vi.classify(e) == "low_complexity"

    # unresolved: > 1000 occurrences of the 20-bp head seed in the second half without an exact match
    head = seq_of(110, 20)
    filler = seq_of(111, 10)
    many = head + (filler + head) * 2100 + "TTTTTTTTTT"
    e = vi.ends(many.encode())
    assert e.capped is True
    assert vi.classify(e) == "unresolved"


def test_stage_v2_dec65_trims_only_exact_single_contig_imgvr_dtr_and_keeps_inovirus_byte_identical(tmp_path):
    # Constructed IMG/VR representatives covering every DEC-65 branch:
    # 1. DTR exact 20 bp -> trimmed, circular, in trims.tsv
    b_exact20 = seq_of(201, 2000)
    s_exact20 = b_exact20 + b_exact20[:20]
    # 2. DTR exact 1000 bp -> trimmed, circular, in trims.tsv
    b_exact1000 = seq_of(202, 3000)
    s_exact1000 = b_exact1000 + b_exact1000[:1000]
    # 3. DTR none -> untrimmed, circular, not in trims.tsv
    s_none = seq_of(203, 2500)
    # 4. DTR partial -> untrimmed, linear, not in trims.tsv
    b_partial = seq_of(204, 2500)
    s_partial = b_partial + mutate(b_partial[:60], 35)
    # 5. DTR long (1001 bp) -> untrimmed, linear, not in trims.tsv
    b_long = seq_of(205, 3000)
    s_long = b_long + b_long[:1001]
    # 6. DTR low_complexity -> untrimmed, linear, not in trims.tsv
    s_lowc = ("AC" * 15) + "GGTT" + seq_of(206, 2000) + "TTGG" + ("AC" * 15)
    # 7. DTR unresolved -> untrimmed, linear, not in trims.tsv
    h_unres = seq_of(207, 20)
    s_unres = h_unres + (seq_of(208, 10) + h_unres) * 2100 + "TTTTTTTTTT"
    # 8. Non-DTR Linear with exact 40 bp repeat -> NEVER trimmed, linear
    b_linrep = seq_of(209, 2000)
    s_linrep = b_linrep + b_linrep[:40]
    # 9. Multi-contig UViG flagged DTR where each contig has exact 40 bp repeat -> NEVER trimmed, linear
    b_mc0 = seq_of(210, 1500)
    b_mc1 = seq_of(211, 1600)
    s_mc = [b_mc0 + b_mc0[:40], b_mc1 + b_mc1[:40]]

    imgvr_specs = [
        ("IMGVR_UViG_01", "vOTU_01", "Direct terminal repeat", [s_exact20]),
        ("IMGVR_UViG_02", "vOTU_02", "Direct terminal repeat", [s_exact1000]),
        ("IMGVR_UViG_03", "vOTU_03", "Direct terminal repeat", [s_none]),
        ("IMGVR_UViG_04", "vOTU_04", "Direct terminal repeat", [s_partial]),
        ("IMGVR_UViG_05", "vOTU_05", "Direct terminal repeat", [s_long]),
        ("IMGVR_UViG_06", "vOTU_06", "Direct terminal repeat", [s_lowc]),
        ("IMGVR_UViG_07", "vOTU_07", "Direct terminal repeat", [s_unres]),
        ("IMGVR_UViG_08", "vOTU_08", "Linear", [s_linrep]),
        ("IMGVR_UViG_09", "vOTU_09", "Direct terminal repeat", s_mc),
    ]
    rows = [row(u, v, sum(len(x) for x in xs), topo) for u, v, topo, xs in imgvr_specs]
    tsv = "\t".join(HEAD) + "\n" + "".join("\t".join(r.get(h, "x") for h in HEAD) + "\n" for r in rows)
    fa = io.StringIO()
    for u, _v, _topo, xs in imgvr_specs:
        for j, s in enumerate(xs):
            fa.write(f">{u}|1|c{j}\n" + "\n".join(s[i:i + 70] for i in range(0, len(s), 70)) + "\n")
    fna_gz = gzip.compress(fa.getvalue().encode(), mtime=0)
    hc = tmp_path / "hc.zip"
    with zipfile.ZipFile(hc, "w") as z:
        z.writestr(zipfile.ZipInfo(SEQ_TSV), tsv.encode(), compress_type=zipfile.ZIP_DEFLATED)
        z.writestr(zipfile.ZipInfo(FNA), fna_gz, compress_type=zipfile.ZIP_STORED)

    # Constructed Inovirus GenBank records (public-database set: NEVER trimmed; circular iff LOCUS circular)
    ino_circ_rep_body = seq_of(301, 5000)
    ino_circ_rep = ino_circ_rep_body + ino_circ_rep_body[:45]  # exact 45 bp repeat, MUST NOT trim, MUST be circular!
    ino_circ_none = seq_of(302, 5200)
    ino_lin_rep_body = seq_of(303, 4800)
    ino_lin_rep = ino_lin_rep_body + ino_lin_rep_body[:45]  # linear fragment with repeat: untrimmed, linear
    ino_files = {
        "Fam/Sub/Sp_circ_rep.gb": genbank("Sp_circ_rep", ino_circ_rep, "circular"),
        "Fam/Sub/Sp_circ_none.gb": genbank("Sp_circ_none", ino_circ_none, "circular"),
        "Fam/Sub/Sp_lin_rep.gb": genbank("Sp_lin_rep", ino_lin_rep, "linear"),
    }
    gb_buf = io.BytesIO()
    with zipfile.ZipFile(gb_buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in ino_files.items():
            z.writestr(zipfile.ZipInfo(name), text.encode())
    meta = tmp_path / "meta.zip"
    with zipfile.ZipFile(meta, "w") as z:
        z.writestr(zipfile.ZipInfo(GB_ZIP), gb_buf.getvalue(), compress_type=zipfile.ZIP_DEFLATED)

    md5 = {
        "seq_tsv": hashlib.md5(tsv.encode()).hexdigest(),
        "fna": hashlib.md5(fna_gz).hexdigest(),
        "gb_zip": hashlib.md5(gb_buf.getvalue()).hexdigest(),
    }
    out = tmp_path / "V2"
    final_dir = tmp_path / "FINAL_V2"
    rep = vi.stage_v2(
        hc, meta, out, root=None, root_stages=(), owner_inputs={}, expected_md5=md5,
        retrieved_imgvr="2026-11-18T23:15:37Z", retrieved_inovirus="2026-11-18T23:16:08Z",
        final_dir=final_dir,
    )

    genomes = json.loads((out / "genomes.json").read_text())
    expected_genomes = [
        [cuid(b_exact20)],       # trimmed!
        [cuid(b_exact1000)],     # trimmed!
        [cuid(s_none)],          # untrimmed
        [cuid(s_partial)],       # untrimmed
        [cuid(s_long)],          # untrimmed
        [cuid(s_lowc)],          # untrimmed
        [cuid(s_unres)],         # untrimmed
        [cuid(s_linrep)],        # untrimmed
        [cuid(s_mc[0]), cuid(s_mc[1])],  # multi-contig untrimmed
        [cuid(ino_circ_none)],   # sorted filename order: Sp_circ_none, Sp_circ_rep, Sp_lin_rep
        [cuid(ino_circ_rep)],    # UNTRIMMED even with exact 45 bp repeat!
        [cuid(ino_lin_rep)],     # UNTRIMMED linear
    ]
    assert [g["contig_uids"] for g in genomes] == expected_genomes

    expected_circular = sorted([
        cuid(b_exact20),
        cuid(b_exact1000),
        cuid(s_none),
        cuid(ino_circ_none),
        cuid(ino_circ_rep),
    ])
    assert sorted(u for g in genomes for u in g["circular"]) == expected_circular

    # trims.tsv contains ONLY the two exact single-contig IMG/VR DTR trims, deterministically sorted
    with open(out / "trims.tsv", newline="") as fh:
        trims = list(csv.DictReader(fh, delimiter="\t"))
    assert len(trims) == 2
    assert [(t["genome"], int(t["repeat_length"]), int(t["bases_removed"]), int(t["length"])) for t in trims] == [
        ("IMGVR_UViG_01", 20, 20, 2000),
        ("IMGVR_UViG_02", 1000, 1000, 3000),
    ]
    assert trims[0]["source_contig_uid"] == cuid(s_exact20)
    assert trims[0]["contig_uid"] == cuid(b_exact20)
    assert trims[1]["source_contig_uid"] == cuid(s_exact1000)
    assert trims[1]["contig_uid"] == cuid(b_exact1000)

    # objects.json uses final_dir
    objects = json.loads((out / "objects.json").read_text())
    assert [Path(o["path"]).parent for o in objects] == [final_dir, final_dir]

    # report tallies
    assert rep["imgvr_dtr_classes"] == {
        "exact": 2, "none": 1, "partial": 1, "long": 1, "low_complexity": 1, "unresolved": 1,
    }
    assert rep["imgvr_trims"] == 2
    assert rep["imgvr_trimmed_bases"] == 1020
    assert rep["inovirus_circular"] == 2
    assert "trims.tsv" in rep["sha256"]

    # Determinism: a second run with the same final_dir is byte-identical across every file
    again = tmp_path / "V2_again"
    vi.stage_v2(
        hc, meta, again, root=None, root_stages=(), owner_inputs={}, expected_md5=md5,
        retrieved_imgvr="2026-11-18T23:15:37Z", retrieved_inovirus="2026-11-18T23:16:08Z",
        final_dir=final_dir,
    )
    for f in sorted(p.name for p in out.iterdir()):
        assert (out / f).read_bytes() == (again / f).read_bytes(), f

"""Card C-61 (DEC-71, DEC-65, DEC-27): staging CGVR viral genomes as stage CG1.

Oracles are independent BLAKE3-128 digests and explicit constructed sequences covering:
- header / bin_id parsing (single-contig CGVR headers and multi-contig `bin|contig` bins);
- DEC-65 terminal-repeat classification & trimming (exact 20..1,000 bp trimmed once & circular;
  none/partial/long/low_complexity/unresolved & multi-contig bins untrimmed & linear);
- DEC-27 repeated-record handling (kept in FASTA, genome listed once, logged in repeated_records.tsv);
- build-order contig_uid ownership across root stages and staged owner inputs (V2, PL2, PL3);
- deterministic byte-identical output across repeated runs and source sha256 verification.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import random
import sys
from pathlib import Path

import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-61", "DEC-65", "DEC-27")  # ring R1 by directory (tests/conftest.py)

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod"


def _load():
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location("c61_cgvr_inputs", SCRIPTS / "c61_cgvr_inputs.py")
    if spec is None or spec.loader is None:
        raise FileNotFoundError(SCRIPTS / "c61_cgvr_inputs.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c61_cgvr_inputs"] = mod
    spec.loader.exec_module(mod)
    return mod


def cuid(seq: str) -> str:
    return blake3(seq.upper().encode()).digest(length=16).hex()


def seq_of(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s: str, pos: int) -> str:
    other = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[pos]]
    return s[:pos] + other + s[pos + 1:]


def _write_gz_fasta(path: Path, records: list[tuple[str, str]]) -> str:
    body = "".join(f">{h}\n{s}\n" for h, s in records).encode()
    raw = gzip.compress(body, mtime=0)
    path.write_bytes(raw)
    return hashlib.sha256(raw).hexdigest()


def test_parse_bin_id_single_and_multi_contig():
    cg = _load()
    assert cg.parse_bin_id("AHFY1_vae_11") == "AHFY1_vae_11"
    assert cg.parse_bin_id("FAHBZ10_FDME220050692-1a_vae_233 extra info") == "FAHBZ10_FDME220050692-1a_vae_233"
    assert cg.parse_bin_id("BIN_01|contig_1 desc") == "BIN_01"
    assert cg.parse_bin_id("BIN_01|contig_2") == "BIN_01"


def test_stage_cg1_dec65_trimming_multi_contig_dec27_owners_and_determinism(tmp_path):
    cg = _load()

    # 1. Exact 141 bp MEGAHIT DTR -> trimmed, circular, in trims.tsv
    b_exact141 = seq_of(6101, 2500)
    s_exact141 = b_exact141 + b_exact141[:141]
    # 2. Exact 20 bp (lower bound) -> trimmed, circular, in trims.tsv
    b_exact20 = seq_of(6102, 2000)
    s_exact20 = b_exact20 + b_exact20[:20]
    # 3. Exact 1000 bp (upper bound) -> trimmed, circular, in trims.tsv
    b_exact1000 = seq_of(6103, 3000)
    s_exact1000 = b_exact1000 + b_exact1000[:1000]
    # 4. None (no repeat, no topology flag) -> untrimmed, linear
    s_none = seq_of(6104, 2400)
    # 5. Partial (near repeat with 1 substitution) -> untrimmed, linear
    b_partial = seq_of(6105, 2400)
    s_partial = b_partial + mutate(b_partial[:60], 30)
    # 6. Long (1001 bp exact repeat > 1000) -> untrimmed, linear
    b_long = seq_of(6106, 3000)
    s_long = b_long + b_long[:1001]
    # 7. Low-complexity exact repeat -> untrimmed, linear
    s_lowc = ("AC" * 15) + "GGTT" + seq_of(6107, 2000) + "TTGG" + ("AC" * 15)
    # 8. DEC-27 repeated record: identical after trimming to #1 (same body, 25 bp exact DTR)
    s_repeat_trimmed = b_exact141 + b_exact141[:25]
    # 9. Multi-contig bin (2 contigs sharing bin_id, even with 40 bp exact repeat) -> untrimmed, linear
    b_mc0 = seq_of(6108, 1500)
    b_mc1 = seq_of(6109, 1600)
    s_mc0 = b_mc0 + b_mc0[:40]
    s_mc1 = b_mc1 + b_mc1[:40]

    records = [
        ("S1_vae_1", s_exact141),
        ("S1_vae_2", s_exact20),
        ("S1_vae_3", s_exact1000),
        ("S1_vae_4", s_none),
        ("S1_vae_5", s_partial),
        ("S1_vae_6", s_long),
        ("S1_vae_7", s_lowc),
        ("S2_vae_1", s_repeat_trimmed),
        ("S3_vae_1|c1", s_mc0),
        ("S3_vae_1|c2", s_mc1),
    ]
    src_fa = tmp_path / "cgvr_src.fa.gz"
    src_sha = _write_gz_fasta(src_fa, records)

    # Staged owner inputs: V2 owns s_none, PL3 owns b_exact20 (trimmed!)
    v2_dir = tmp_path / "owners" / "V2"
    pl3_dir = tmp_path / "owners" / "PL3"
    v2_dir.mkdir(parents=True)
    pl3_dir.mkdir(parents=True)
    (v2_dir / "genomes.json").write_text(json.dumps([
        {"contig_uids": [cuid(s_none)], "genetic_code": 11, "kind": "virus", "circular": [], "assembly_quality": "vOTU"}
    ]))
    (pl3_dir / "genomes.json").write_text(json.dumps([
        {"contig_uids": [cuid(b_exact20)], "genetic_code": 11, "kind": "plasmid", "circular": [cuid(b_exact20)],
         "assembly_quality": "complete"}
    ]))

    out = tmp_path / "CG1"
    final_dir = tmp_path / "FINAL_CG1"
    rep = cg.stage_cg1(
        src_fa, out, root=None, root_stages=(), owner_inputs={"V2": v2_dir, "PL3": pl3_dir},
        expected_sha256=src_sha, retrieved="2026-10-23T14:49:00Z", final_dir=final_dir,
    )

    genomes = json.loads((out / "genomes.json").read_text())
    expected_genomes = [
        [cuid(b_exact141)],         # S1_vae_1 (trimmed 141 bp)
        [cuid(b_exact20)],          # S1_vae_2 (trimmed 20 bp)
        [cuid(b_exact1000)],        # S1_vae_3 (trimmed 1000 bp)
        [cuid(s_none)],             # S1_vae_4 (untrimmed linear)
        [cuid(s_partial)],          # S1_vae_5 (untrimmed linear)
        [cuid(s_long)],             # S1_vae_6 (untrimmed linear)
        [cuid(s_lowc)],             # S1_vae_7 (untrimmed linear)
        # S2_vae_1 repeats S1_vae_1 after trimming -> listed once (DEC-27)
        [cuid(s_mc0), cuid(s_mc1)],  # S3_vae_1 multi-contig untrimmed linear
    ]
    assert [g["contig_uids"] for g in genomes] == expected_genomes
    assert all(g["kind"] == "virus" and g["genetic_code"] == 11 and g["assembly_quality"] == "vOTU" for g in genomes)

    expected_circular = sorted([cuid(b_exact141), cuid(b_exact20), cuid(b_exact1000)])
    assert sorted(u for g in genomes for u in g["circular"]) == expected_circular

    # All 10 records are preserved in the staged FASTA (DEC-27: never drop)
    with gzip.open(out / cg.CGVR_FASTA, "rt") as fh:
        staged_headers = [ln[1:].strip() for ln in fh if ln.startswith(">")]
    assert len(staged_headers) == 10

    # trims.tsv has all 4 exact single-contig trims (including S2_vae_1 before DEC-27 dedup)
    with open(out / "trims.tsv", newline="") as fh:
        trims = list(csv.DictReader(fh, delimiter="\t"))
    assert [(t["genome"], int(t["repeat_length"]), int(t["bases_removed"]), int(t["length"])) for t in trims] == [
        ("S1_vae_1", 141, 141, 2500),
        ("S1_vae_2", 20, 20, 2000),
        ("S1_vae_3", 1000, 1000, 3000),
        ("S2_vae_1", 25, 25, 2500),
    ]

    # repeated_records.tsv logs S2_vae_1 -> S1_vae_1
    with open(out / "repeated_records.tsv", newline="") as fh:
        reps = list(csv.DictReader(fh, delimiter="\t"))
    assert len(reps) == 1 and reps[0]["genome"] == "S2_vae_1" and reps[0]["first_genome"] == "S1_vae_1"

    # objects.json points to final_dir and is redistributable
    objects = json.loads((out / "objects.json").read_text())
    assert len(objects) == 1
    assert Path(objects[0]["path"]) == final_dir / cg.CGVR_FASTA
    assert objects[0]["redistributable"] is True
    assert objects[0]["source_name"] == "cgvr"

    # ownership counts in report
    assert rep["owner_contigs"] == {"PL3": 1, "V2": 1}
    assert rep["owner_bases"] == {"PL3": 2000, "V2": 2400}
    assert rep["trims"] == 4
    assert rep["trimmed_bases"] == 141 + 20 + 1000 + 25
    assert rep["dtr_classes"] == {
        "exact": 4, "none": 1, "partial": 1, "long": 1, "low_complexity": 1, "multi_contig": 2,
    }

    # Determinism: second run is byte-identical across every file
    again = tmp_path / "CG1_again"
    cg.stage_cg1(
        src_fa, again, root=None, root_stages=(), owner_inputs={"V2": v2_dir, "PL3": pl3_dir},
        expected_sha256=src_sha, retrieved="2026-10-23T14:49:00Z", final_dir=final_dir,
    )
    for f in sorted(p.name for p in out.iterdir()):
        assert (out / f).read_bytes() == (again / f).read_bytes(), f


def test_stage_cg1_refuses_wrong_sha256(tmp_path):
    cg = _load()
    src_fa = tmp_path / "bad.fa.gz"
    _write_gz_fasta(src_fa, [("S1_vae_1", seq_of(6110, 500))])
    with pytest.raises(cg.Refused):
        cg.stage_cg1(
            src_fa, tmp_path / "out", root=None, root_stages=(), owner_inputs={},
            expected_sha256="0" * 64, retrieved="2026-10-23T14:49:00Z",
        )
    assert not (tmp_path / "out").exists()

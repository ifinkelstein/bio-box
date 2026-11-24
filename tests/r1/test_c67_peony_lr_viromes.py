"""Ring R1: Card C-67 (`LRV1` long-read aquatic viromes and `V4` PEONY 1.3 / PIGEON 3.0 staging)."""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import tarfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-67")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _load_c67():
    spec = importlib.util.spec_from_file_location("c67_peony_lr_viromes", BSP / "c67_peony_lr_viromes.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c67_peony_lr_viromes"] = mod
    spec.loader.exec_module(mod)
    return mod


def _make_seq(seed_int: int, length: int) -> bytes:
    bases = b"ACGT"
    x = seed_int
    out = bytearray(length)
    for i in range(length):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        out[i] = bases[(x >> 16) & 3]
    return bytes(out)


def test_dec65_dtr_trim_metagenomic_exact_vs_linear():
    c67 = _load_c67()
    body = _make_seq(42, 4000)
    dtr = _make_seq(99, 120)
    exact_seq = dtr + body + dtr
    trimmed, k, dtr_cls, is_circ = c67.apply_dec65_metagenomic(exact_seq)
    assert dtr_cls == "exact"
    assert k == 120
    assert is_circ is True
    assert trimmed == dtr + body

    # Low-complexity exact repeat (N-containing repeat of 30 bp) must NOT be trimmed and must be linear
    lc_rep = _make_seq(77, 15) + b"N" + _make_seq(78, 14)
    lc_seq = lc_rep + body + lc_rep
    trimmed_lc, k_lc, cls_lc, circ_lc = c67.apply_dec65_metagenomic(lc_seq)
    assert cls_lc == "low_complexity"
    assert k_lc == 0
    assert circ_lc is False
    assert trimmed_lc == lc_seq

    # Plain linear sequence without DTR
    lin_seq = _make_seq(7, 3500)
    trimmed_lin, k_lin, cls_lin, circ_lin = c67.apply_dec65_metagenomic(lin_seq)
    assert cls_lin == "none"
    assert k_lin == 0
    assert circ_lin is False
    assert trimmed_lin == lin_seq

    # Exact repeat > 1,000 bp (long) must NOT be trimmed and must be linear
    long_rep = _make_seq(123, 1200)
    long_seq = long_rep + body + long_rep
    trimmed_long, k_long, cls_long, circ_long = c67.apply_dec65_metagenomic(long_seq)
    assert cls_long == "long"
    assert k_long == 0
    assert circ_long is False
    assert trimmed_long == long_seq


def test_stage_lrv1_trims_dtrs_dedups_owners_and_writes_metadata(tmp_path):
    from biobox import ids

    c67 = _load_c67()
    src_dir = tmp_path / "src"
    src_dir.mkdir()

    # 1. Beaulaurier AFVG with exact 85 bp DTR
    dtr1 = _make_seq(11, 85)
    body1 = _make_seq(12, 3000)
    seq_beau_dtr = dtr1 + body1 + dtr1
    seq_beau_lin = _make_seq(13, 2800)
    f_beau = src_dir / "beaulaurier_afvgs.fna"
    f_beau.write_bytes(
        b">MN694825.1 AFPP_250M12\n" + seq_beau_dtr + b"\n"
        b">MN694826.1 AFVG_1\n" + seq_beau_lin + b"\n"
    )

    # 2. Kavagutti phage (one unique linear, one identical to trimmed Beaulaurier -> within-stage repeat!)
    seq_kava_new = _make_seq(21, 5000)
    f_kava = src_dir / "kavagutti_phages.fna"
    f_kava.write_bytes(
        b">kava_phage_1 complete\n" + seq_kava_new + b"\n"
        b">kava_dup_of_beau\n" + (dtr1 + body1) + b"\n"
    )

    # 3. VirION2 contig (already owned by earlier stage V2)
    seq_owned = _make_seq(31, 4200)
    f_vir2 = src_dir / "virion2_wec.fna"
    f_vir2.write_bytes(b">wec_virion2_1\n" + seq_owned + b"\n")

    # Fake earlier stage V2 owner_inputs directory
    v2_dir = tmp_path / "owner_V2"
    v2_dir.mkdir()
    owned_uid = ids.contig_uid(seq_owned).hex()
    (v2_dir / "genomes.json").write_text(json.dumps([{"contig_uids": [owned_uid]}]))

    sources_spec = [
        {
            "name": "beaulaurier_aloha_afvgs",
            "category": "ont_marine_afvg",
            "fasta": str(f_beau),
            "out_name": "beaulaurier_aloha_afvgs.fna",
            "source_release": "Beaulaurier 2020 PRJNA529454",
            "license": "INSDC (public domain)",
            "redistributable": True,
        },
        {
            "name": "kavagutti_freshwater_hifi_phages",
            "category": "hifi_freshwater_phage",
            "fasta": str(f_kava),
            "out_name": "kavagutti_freshwater_hifi_phages.fna",
            "source_release": "Kavagutti 2023 Figshare 20431332",
            "license": "CC BY 4.0",
            "redistributable": True,
        },
        {
            "name": "virion2_wec_marine_virome",
            "category": "ont_marine_virion2",
            "fasta": str(f_vir2),
            "out_name": "virion2_wec_marine_virome.fna",
            "source_release": "Zablocki 2021 PeerJ e11088",
            "license": "CC BY 4.0",
            "redistributable": True,
        },
    ]

    out = tmp_path / "LRV1"
    rep = c67.stage_lrv1(
        sources_spec,
        out,
        root=None,
        root_stages=[],
        owner_inputs={"V2": v2_dir},
        retrieved="2026-11-24T00:00:00Z",
    )

    assert rep["total_records"] == 5
    assert rep["trims"] == 1
    assert rep["trimmed_bases"] == 85
    assert rep["genomes"] == 4  # 1 within-stage repeat (kava_dup_of_beau matches trimmed MN694825.1)
    assert rep["repeated_genomes"] == 1
    assert rep["circular_genomes"] == 1
    assert rep["held_contigs"] == 1
    assert rep["owner_contigs"] == {"V2": 1}
    assert rep["new_contigs"] == 3

    genomes = json.loads((out / "genomes.json").read_text())
    assert len(genomes) == 4
    assert all(g["kind"] == "virus" and g["genetic_code"] == 11 for g in genomes)
    objects = json.loads((out / "objects.json").read_text())
    assert len(objects) == 3
    assert all(o["redistributable"] is True for o in objects)
    assert (out / "trims.tsv").is_file()
    assert (out / "lrv1_metadata.tsv").is_file()
    assert (out / "records.tsv").is_file()
    assert (out / "shared_records.tsv").is_file()
    assert (out / "repeated_records.tsv").is_file()


def test_stage_lrv1_refuses_non_redistributable(tmp_path):
    c67 = _load_c67()
    fa = tmp_path / "x.fna"
    fa.write_bytes(b">a\n" + _make_seq(1, 2000) + b"\n")
    with pytest.raises(c67.Refused):
        c67.stage_lrv1(
            [{
                "name": "bad",
                "category": "test",
                "fasta": str(fa),
                "out_name": "bad.fna",
                "source_release": "test",
                "license": "restricted",
                "redistributable": False,
            }],
            tmp_path / "LRV1_bad",
            root=None,
            root_stages=[],
            owner_inputs={},
            retrieved="2026-11-24T00:00:00Z",
        )


def test_extract_kavagutti_phages_from_nested_support_tar(tmp_path):
    c67 = _load_c67()
    inner_tar = tmp_path / "phages-fna.tar.gz"
    with tarfile.open(inner_tar, "w:gz") as tf:
        for name, seq in (("phage_b.fna", _make_seq(2, 1500)), ("phage_a.fna", _make_seq(1, 1500))):
            data = f">{name[:-4]}\n".encode() + seq + b"\n"
            info = tarfile.TarInfo(name=f"phages/{name}")
            info.size = len(data)
            tf.addfile(info, fileobj=c67.io.BytesIO(data))

    outer_tar = tmp_path / "Support-data.gz"
    with tarfile.open(outer_tar, "w:gz") as tf:
        tf.add(inner_tar, arcname="Support-data/phages-fna.tar.gz")

    out_fa = tmp_path / "kavagutti_extracted.fna"
    res = c67.extract_kavagutti_phages(outer_tar, out_fa)
    assert res["records"] == 2
    assert res["bases"] == 3000
    # Deterministic sorted order by member name: phage_a before phage_b
    lines = out_fa.read_text().splitlines()
    assert lines[0] == ">phage_a"
    assert lines[2] == ">phage_b"


def test_stage_v4_streams_peony_and_pigeon_dedups_and_splits_layout(tmp_path):
    from biobox import ids

    c67 = _load_c67()
    dl = tmp_path / "dl"
    dl.mkdir()

    dtr = _make_seq(55, 60)
    body = _make_seq(56, 12000)
    seq_shared_peony_pigeon = dtr + body + dtr
    seq_peony_only = _make_seq(57, 11000)
    seq_pigeon_only = _make_seq(58, 10500)
    seq_held_v3 = _make_seq(59, 10200)

    peony_gz = dl / "PEONY1.3_Nov26.fa.gz"
    with gzip.open(peony_gz, "wb") as fh:
        fh.write(b">peony_votu_1\n" + seq_shared_peony_pigeon + b"\n")
        fh.write(b">peony_votu_2\n" + seq_peony_only + b"\n")

    pigeon_gz = dl / "PIGEON3.0.fa.gz"
    with gzip.open(pigeon_gz, "wb") as fh:
        fh.write(b">pigeon_dup_1\n" + seq_shared_peony_pigeon + b"\n")
        fh.write(b">pigeon_votu_3\n" + seq_pigeon_only + b"\n")
        fh.write(b">pigeon_held_v3\n" + seq_held_v3 + b"\n")

    v3_dir = tmp_path / "owner_V3"
    v3_dir.mkdir()
    (v3_dir / "genomes.json").write_text(json.dumps([{"contig_uids": [ids.contig_uid(seq_held_v3).hex()]}]))

    fasta_out = tmp_path / "stor_V4" / "v4_peony_pigeon3.fna"
    meta_out = tmp_path / "tmp_V4"
    final_dir = tmp_path / "final_V4"

    rep = c67.stage_v4(
        [("peony1.3", peony_gz), ("pigeon3.0", pigeon_gz)],
        fasta_out=fasta_out,
        meta_out=meta_out,
        final_dir=final_dir,
        root=None,
        root_stages=[],
        owner_inputs={"V3": v3_dir},
        retrieved="2026-11-24T00:00:00Z",
    )

    assert rep["total_records"] == 5
    assert rep["genomes"] == 4  # pigeon_dup_1 collapsed into peony_votu_1
    assert rep["repeated_genomes"] == 1
    assert rep["trims"] == 2  # both input records had exact DTR before dedup (or 1 if deduped on write)
    assert rep["circular_genomes"] == 1
    assert rep["held_contigs"] == 1
    assert rep["owner_contigs"] == {"V3": 1}
    assert rep["new_contigs"] == 3

    assert fasta_out.is_file()
    assert (meta_out / "genomes.json").is_file()
    assert (meta_out / "objects.json").is_file()
    assert (meta_out / "objects_staging.json").is_file()
    assert (meta_out / "V4_left_for_later.tsv.gz").is_file()

    objs_final = json.loads((meta_out / "objects.json").read_text())
    objs_stage = json.loads((meta_out / "objects_staging.json").read_text())
    assert objs_final[0]["path"] == str(final_dir / fasta_out.name)
    assert objs_stage[0]["path"] == str(fasta_out)
    assert objs_final[0]["redistributable"] is True

    # SHA-256 verification on cached file succeeds on match and refuses on mismatch
    n, md5_h, sha_h = c67._download_file(
        "https://example.invalid/unused", peony_gz, state={"last": 0.0, "refusals": 0},
        expected_sha256=hashlib.sha256(peony_gz.read_bytes()).hexdigest(),
    )
    assert n == peony_gz.stat().st_size
    assert sha_h == hashlib.sha256(peony_gz.read_bytes()).hexdigest()
    with pytest.raises(c67.Refused):
        c67._download_file(
            "https://example.invalid/unused", peony_gz, state={"last": 0.0, "refusals": 0},
            expected_sha256="0" * 64,
        )

"""Ring R1: Card C-69 — staging Egovirales gut giant-virus genomes (EGV1).

Tests:
- Exact read-only ITR detection (including 55 bp ITRs) without ever trimming a single base (DEC-65).
- Byte-identical staging of linear untrimmed `kind="virus"` (`genetic_code=11`) genomes with no 5 kb floor.
- Byte-for-byte determinism across archive member ordering and zip timestamps (AGENTS.md §4.3).
- Within-stage `genome_uid` deduplication (DEC-27) and cross-stage prefix ownership (DEC-26).
- Refusal on non-redistributable sources (DEC-51), MD5 mismatch, and existing output directory.
- `EGV1.conf` validation via `stage_est.read_conf`.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import zipfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-69")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
STAGES_TEMPLATES = BSP / "stages" / "templates"


def _load(name: str):
    path = BSP / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"missing {path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _revcomp(seq: bytes) -> bytes:
    return seq.translate(bytes.maketrans(b"ACGTacgt", b"TGCAtgca"))[::-1]


def _make_zip(path: Path, members: list[tuple[str, bytes]], *, date_time=(2026, 8, 25, 8, 20, 26)) -> Path:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, data in members:
            zi = zipfile.ZipInfo(name, date_time=date_time)
            zf.writestr(zi, data)
    return path


def test_detect_itr_finds_exact_55bp_itr_without_modifying_sequence():
    m = _load("c69_egovirales_inputs")
    itr55 = b"ACGTTGCAAACCGGTTACGTTGCAAACCGGTTACGTTGCAAACCGGTTACGTTGC"
    assert len(itr55) == 55
    body = b"ATGCGATCGATCGATC" * 40
    seq = itr55 + body + _revcomp(itr55)
    orig = bytes(seq)

    assert m.detect_itr(seq) == 55
    assert seq == orig  # read-only, never modified

    # Mixed-case sequence still detects exact 55 bp ITR
    mixed = itr55.lower() + body + _revcomp(itr55).upper()
    assert m.detect_itr(mixed) == 55

    # Non-ITR sequence returns 0
    assert m.detect_itr(b"AAAAAAAAAA" + body + b"CCCCCCCCCC") == 0

    # Short (< min_len) inverted repeat returns 0
    short = b"ACGTACGTACGTACG"  # 15 bp
    assert m.detect_itr(short + body + _revcomp(short), min_len=20) == 0


def test_stage_egv1_preserves_55bp_itrs_untrimmed_no_5kb_floor_and_linear_topology(tmp_path):
    from biobox import ids

    m = _load("c69_egovirales_inputs")
    itr55 = b"ACGTTGCAAACCGGTTACGTTGCAAACCGGTTACGTTGCAAACCGGTTACGTTGC"
    complete_seq = itr55 + (b"GCTAGCTAGGATCCAT" * 400) + _revcomp(itr55)  # 6,510 bp with 55 bp ITR
    mag_c1 = b"TTGGCCAATTGGCCAA" * 350  # 5,600 bp
    mag_c2_small = b"ACGTACGTGGGCCCAA" * 75  # 1,200 bp (< 5 kb; must be kept)

    zpath = _make_zip(
        tmp_path / "ITEM_08_Egovirales_Genomes.zip",
        [
            ("ITEM_08_Egovirales_Genomes/README.txt", b"Egovirales 224 genomes\n"),
            (
                "ITEM_08_Egovirales_Genomes/EGO_complete_01.fa",
                b">EGO_complete_01_c1 complete linear genome\n" + complete_seq + b"\n",
            ),
            (
                "ITEM_08_Egovirales_Genomes/EGO_mag_02.fna",
                b">EGO_mag_02_c1\n" + mag_c1 + b"\n>EGO_mag_02_c2\n" + mag_c2_small + b"\n",
            ),
        ],
    )
    z_md5 = hashlib.md5(zpath.read_bytes()).hexdigest()

    out = tmp_path / "EGV1"
    rep = m.stage_egv1(
        zpath,
        out,
        root=None,
        root_stages=(),
        owner_inputs={},
        retrieved="2026-11-24T00:00:00Z",
        expected_md5=z_md5,
    )

    assert rep["stage"] == "EGV1"
    assert rep["redistributable"] is True
    assert rep["trims"] == 0
    assert rep["circular_genomes"] == 0
    assert rep["linear_genomes"] == 2
    assert rep["genomes"] == 2
    assert rep["complete_itr_genomes"] == 1
    assert rep["contigs"] == 3
    assert rep["lt_5kb_contigs"] == 1
    assert rep["skipped_members"] == 1
    assert rep["bases"] == len(complete_seq) + len(mag_c1) + len(mag_c2_small)

    # Exact byte-identical contig_uids (ITR preserved untrimmed, < 5 kb contig kept)
    genomes = json.loads((out / "genomes.json").read_text())
    assert len(genomes) == 2
    by_first = {g["contig_uids"][0]: g for g in genomes}
    cu_complete = ids.contig_uid(complete_seq).hex()
    cu_mag1 = ids.contig_uid(mag_c1).hex()
    cu_mag2 = ids.contig_uid(mag_c2_small).hex()
    assert cu_complete in by_first
    assert by_first[cu_complete]["kind"] == "virus"
    assert by_first[cu_complete]["genetic_code"] == 11
    assert by_first[cu_complete]["circular"] == []
    assert by_first[cu_complete]["assembly_quality"] == "complete"
    assert by_first[cu_mag1]["contig_uids"] == [cu_mag1, cu_mag2]
    assert by_first[cu_mag1]["circular"] == []
    assert by_first[cu_mag1]["assembly_quality"] == "contig"

    # objects.json carries CC BY 4.0 and redistributable: True
    objs = json.loads((out / "objects.json").read_text())
    assert len(objs) == 1
    assert objs[0]["source_name"] == "egovirales_gaia_2026"
    assert objs[0]["license"] == "CC BY 4.0"
    assert objs[0]["redistributable"] is True
    assert objs[0]["path"] == str(out / "egovirales_gaia_2026.fna")

    # trims.tsv always has 0 data rows
    with open(out / "trims.tsv", newline="") as fh:
        assert list(csv.DictReader(fh, delimiter="\t")) == []

    # genome_input.tsv and egv1_metadata.tsv record ITR=55 and linear topology
    with open(out / "genome_input.tsv", newline="") as fh:
        gi_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(gi_rows) == 3
    by_cu = {r["contig_uid"]: r for r in gi_rows}
    assert by_cu[cu_complete]["itr_length"] == "55"
    assert by_cu[cu_complete]["topology"] == "linear"
    assert by_cu[cu_complete]["circular"] == "false"
    assert by_cu[cu_mag2]["length"] == str(len(mag_c2_small))

    with open(out / "egv1_metadata.tsv", newline="") as fh:
        meta_rows = list(csv.DictReader(fh, delimiter="\t"))
    by_name = {r["genome_name"]: r for r in meta_rows}
    assert by_name["EGO_complete_01"]["max_itr_length"] == "55"
    assert by_name["EGO_complete_01"]["assembly_quality"] == "complete"
    assert by_name["EGO_mag_02"]["contig_count"] == "2"

    # inputs_SUMS covers core files
    sums_text = (out / "inputs_SUMS").read_text()
    for fname in ("genomes.json", "objects.json", "egovirales_gaia_2026.fna", "genome_input.tsv", "egv1_metadata.tsv"):
        assert fname in sums_text


def test_stage_egv1_deterministic_byte_identical_across_zip_member_orders(tmp_path):
    m = _load("c69_egovirales_inputs")
    itr = b"CGTACGTTGCAAACCGGTTACGTACGTTGCAAACCGGTTACGTACGTTGCAAACC"
    fa_a = b">ego_a_1\n" + itr + (b"ACGT" * 200) + _revcomp(itr) + b"\n"
    fa_b = b">ego_b_1\n" + (b"GCTA" * 180) + b"\n>ego_b_2\n" + (b"TTAA" * 120) + b"\n"

    z1 = _make_zip(
        tmp_path / "z1.zip",
        [("genomes/ego_a.fa", fa_a), ("genomes/ego_b.fa", fa_b)],
        date_time=(2026, 1, 1, 0, 0, 0),
    )
    z2 = _make_zip(
        tmp_path / "z2.zip",
        [("genomes/ego_b.fa", fa_b), ("genomes/ego_a.fa", fa_a)],
        date_time=(2026, 11, 24, 12, 30, 0),
    )

    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"
    final = Path("/tmp/if2453/s2/inputs/EGV1")
    m.stage_egv1(z1, out1, root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z", final_dir=final)
    m.stage_egv1(z2, out2, root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z", final_dir=final)

    for fname in (
        "egovirales_gaia_2026.fna",
        "genomes.json",
        "objects.json",
        "genome_input.tsv",
        "records.tsv",
        "egv1_metadata.tsv",
        "trims.tsv",
        "inputs_SUMS",
    ):
        assert (out1 / fname).read_bytes() == (out2 / fname).read_bytes(), f"{fname} differed across member orders"


def test_stage_egv1_dec27_dedup_and_dec26_prefix_overlap(tmp_path):
    from biobox import ids

    m = _load("c69_egovirales_inputs")
    seq_shared_gv1 = b"AACCGGTTAACCGGTT" * 100
    seq_unique = b"CGTACGTAGCTAGCTA" * 100

    # Three genomes: g1 (seq_shared_gv1), g2 (exact duplicate of g1 -> DEC-27 repeated),
    # g3 (multi-contig sharing seq_shared_gv1 + seq_unique -> distinct genome_uid, kept under DEC-27)
    zpath = _make_zip(
        tmp_path / "egv_dedup.zip",
        [
            ("g1.fna", b">g1_c1\n" + seq_shared_gv1 + b"\n"),
            ("g2.fna", b">g2_c1\n" + seq_shared_gv1 + b"\n"),
            ("g3.fna", b">g3_c1\n" + seq_shared_gv1 + b"\n>g3_c2\n" + seq_unique + b"\n"),
        ],
    )

    # Simulate earlier GV1 stage owning seq_shared_gv1 via records.tsv (DEC-26)
    gv1_dir = tmp_path / "GV1"
    gv1_dir.mkdir()
    cu_gv1 = ids.contig_uid(seq_shared_gv1).hex()
    (gv1_dir / "records.tsv").write_text(
        f"ordinal\tgenome\tcontig_uid\tlength\n0\tgv_ref\t{cu_gv1}\t{len(seq_shared_gv1)}\n"
    )

    out = tmp_path / "EGV1"
    rep = m.stage_egv1(
        zpath,
        out,
        root=None,
        root_stages=(),
        owner_inputs={"GV1": gv1_dir},
        retrieved="2026-11-24T00:00:00Z",
    )

    assert rep["genomes"] == 2
    assert rep["repeated_genomes"] == 1
    assert rep["contigs"] == 2
    assert rep["held_contigs"] == 1
    assert rep["new_contigs"] == 1
    assert rep["owned_by_stage"] == {"GV1": 1}
    assert rep["genomes_all_contigs_held_earlier"] == 1

    with open(out / "repeated_records.tsv", newline="") as fh:
        rep_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rep_rows) == 1
    assert rep_rows[0]["genome"] == "g2"
    assert rep_rows[0]["first_genome"] == "g1"

    with open(out / "shared_records.tsv", newline="") as fh:
        sh_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(sh_rows) == 3  # g1_c1, g2_c1, g3_c1 all carry seq_shared_gv1
    assert all(r["owner_stage"] == "GV1" for r in sh_rows)


def test_stage_egv1_refuses_non_redistributable_md5_mismatch_and_existing_out(tmp_path):
    m = _load("c69_egovirales_inputs")
    zpath = _make_zip(tmp_path / "egv.zip", [("g1.fa", b">g1\n" + b"ACGT" * 50 + b"\n")])

    with pytest.raises(m.Refused, match="non-redistributable"):
        m.stage_egv1(
            zpath,
            tmp_path / "out_nr",
            root=None,
            root_stages=(),
            owner_inputs={},
            retrieved="2026-11-24T00:00:00Z",
            redistributable=False,
        )

    with pytest.raises(m.Refused, match="md5 mismatch"):
        m.stage_egv1(
            zpath,
            tmp_path / "out_md5",
            root=None,
            root_stages=(),
            owner_inputs={},
            retrieved="2026-11-24T00:00:00Z",
            expected_md5="0" * 32,
        )

    existing = tmp_path / "existing"
    existing.mkdir()
    with pytest.raises(m.Refused, match="exists"):
        m.stage_egv1(
            zpath,
            existing,
            root=None,
            root_stages=(),
            owner_inputs={},
            retrieved="2026-11-24T00:00:00Z",
        )


def test_egv1_conf_parses_with_stage_est_read_conf():
    est = _load("stage_est")
    c = est.read_conf(STAGES_TEMPLATES / "EGV1.conf", "EGV1", tree=REPO)
    assert c.stage == "EGV1"
    assert c.inputs == Path("/tmp/if2453/s2/inputs/EGV1")
    assert c.depth == "annotated"
    assert c.genecall_mode == "gv-meta"
    assert c.library == ("pfam", "38.2", "gathering")
    assert c.crispr is True
    assert c.min_contig_len is None
    assert c.after == ("GV1",)
    assert "planning/measurements/2026-11-24_egv1_prep.md" in c.notes


# --- additions of the staging session (2026-11-26): the zip's real layout, the r26.11.4 chain, the record -----------

INPUTS_NOTE = REPO / "planning" / "measurements" / "2026-11-26_egv1_inputs.md"
ZIP_ROOT = "ITEM_08_Egovirales_Genomes/"
FOLDER_01 = ZIP_ROOT + "01_Egovirales_224_MAGs/"


def test_stage_egv1_members_prefix_takes_the_224_folder_and_skips_its_byte_identical_copies(tmp_path):
    m = _load("c69_egovirales_inputs")
    fa_a = b">a_scaffold_1\n" + b"GATTACAGATTACAGG" * 60 + b"\n"
    fa_b = b">b_scaffold_1\n" + b"TTGACCAATGGCCAAT" * 70 + b"\n"
    members = [
        (ZIP_ROOT + "READ_ME.txt", b"This item holds folders\n"),
        (ZIP_ROOT + ".DS_Store", b"\x00\x00"),
        ("__MACOSX/" + ZIP_ROOT + "01_Egovirales_224_MAGs/._a-contigs.fa", b"\x00junk"),
        (FOLDER_01 + "a-contigs.fa", fa_a),
        (FOLDER_01 + "b-contigs.fa", fa_b),
        (ZIP_ROOT + "02_Egovirales_116_nrMAgs-indidual_genomes/a-contigs.fa", fa_a),  # same name, same bytes
        (ZIP_ROOT + "03_Egovirales_6_nrMAGs-human-indidual_genomes/all_six.fa", fa_a + fa_b),
        (ZIP_ROOT + "04_ZEEV/DNApol_B", b">p\nMKV\n"),
    ]
    zpath = _make_zip(tmp_path / "item08.zip", members)
    kw = dict(root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-26T00:00:00Z",
              final_dir=Path("/tmp/if2453/s2/inputs/EGV1"))

    with pytest.raises(m.Refused, match="duplicate genome name"):
        m.stage_egv1(zpath, tmp_path / "all", **kw)

    rep = m.stage_egv1(zpath, tmp_path / "only01", members_prefix=FOLDER_01, **kw)
    assert rep["genomes"] == 2
    assert rep["fasta_members"] == 2
    assert rep["contigs"] == 2
    assert rep["bases"] == 16 * 60 + 16 * 70
    with open(tmp_path / "only01" / "skipped_members.tsv", newline="") as fh:
        reasons = {r["member"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")}
    assert reasons[ZIP_ROOT + "02_Egovirales_116_nrMAgs-indidual_genomes/a-contigs.fa"] == "outside_members_prefix"
    assert reasons[ZIP_ROOT + "03_Egovirales_6_nrMAGs-human-indidual_genomes/all_six.fa"] == "outside_members_prefix"
    assert all(r in ("outside_members_prefix", "ignored_system_entry", "not_fasta") for r in reasons.values())

    # the folder's own order and the zip's member order do not matter
    zrev = _make_zip(tmp_path / "item08_rev.zip", list(reversed(members)), date_time=(2026, 11, 1, 1, 1, 2))
    m.stage_egv1(zrev, tmp_path / "rev", members_prefix=FOLDER_01, **kw)
    for fname in ("egovirales_gaia_2026.fna", "genomes.json", "records.tsv", "skipped_members.tsv", "inputs_SUMS"):
        assert (tmp_path / "only01" / fname).read_bytes() == (tmp_path / "rev" / fname).read_bytes(), fname


def test_chain_owner_inputs_follow_the_r26_11_4_order_and_end_with_gv1(tmp_path):
    m = _load("c69_egovirales_inputs")
    a26 = _load("a26_r26_11_4_chain")
    owners = m.chain_owner_inputs(tmp_path)
    assert list(owners) == [s for s, _ in a26.R26_11_4_CHAIN] + ["GV1"]
    assert len(owners) == 20
    assert owners["V2"] == tmp_path / "V2"
    assert owners["GV1"] == tmp_path / "GV1"


def test_stage_egv1_owner_inputs_in_chain_order_first_holder_owns(tmp_path):
    from biobox import ids

    m = _load("c69_egovirales_inputs")
    seq_both = b"AACCGGTTAACCGGTT" * 100
    seq_new = b"CGTACGTAGCTAGCTA" * 100
    zpath = _make_zip(tmp_path / "z.zip", [("g.fna", b">c1\n" + seq_both + b"\n>c2\n" + seq_new + b"\n")])
    for name in ("MV1", "GV1"):  # both hold seq_both; MV1 is earlier in the chain, so it owns it
        d = tmp_path / "held" / name
        d.mkdir(parents=True)
        (d / "records.tsv").write_text(f"ordinal\tgenome\tcontig_uid\tlength\n0\tx\t{ids.contig_uid(seq_both).hex()}\t1600\n")
    rep = m.stage_egv1(zpath, tmp_path / "o", root=None, root_stages=(),
                       owner_inputs={"MV1": tmp_path / "held" / "MV1", "GV1": tmp_path / "held" / "GV1"},
                       retrieved="2026-11-26T00:00:00Z")
    assert rep["owned_by_stage"] == {"MV1": 1}
    assert rep["new_contigs"] == 1


def _conf_lines(name: str) -> list[str]:
    return (STAGES_TEMPLATES / name).read_text().splitlines()


def test_egv1_template_has_real_sums_found_in_the_inputs_note_and_the_chain_expectations():
    est = _load("stage_est")
    c = est.read_conf(STAGES_TEMPLATES / "EGV1.conf", "EGV1", tree=REPO)
    assert "planning/measurements/2026-11-26_egv1_inputs.md" in c.notes
    note = INPUTS_NOTE.read_text()
    sums = [ln.split()[0] for ln in _conf_lines("EGV1.conf") if ln.startswith("SUM=")]
    assert len(sums) == 3
    for s in (x.removeprefix("SUM=") for x in sums):
        assert len(s) == 64 and set(s) != {"0"}, "placeholder checksum left in the template"
        assert s in note, f"{s} is not recorded in the note"
    assert c.after == ("GV1",)
    assert c.depth == "annotated" and c.crispr is True and c.min_contig_len is None


def test_stage_egv1_writes_objects_staging_json_naming_the_file_where_it_is_outside_the_sums(tmp_path):
    m = _load("c69_egovirales_inputs")
    zpath = _make_zip(tmp_path / "z.zip", [("g.fna", b">c1\n" + b"ACGGTCAA" * 100 + b"\n")])
    final = Path("/tmp/if2453/s2/inputs/EGV1")
    out = tmp_path / "EGV1"
    m.stage_egv1(zpath, out, root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-26T00:00:00Z",
                 final_dir=final)
    staged = json.loads((out / "objects_staging.json").read_text())
    final_objs = json.loads((out / "objects.json").read_text())
    assert staged[0]["path"] == str(out / "egovirales_gaia_2026.fna")
    assert final_objs[0]["path"] == str(final / "egovirales_gaia_2026.fna")
    assert {k: v for k, v in staged[0].items() if k != "path"} == {k: v for k, v in final_objs[0].items() if k != "path"}
    assert "objects_staging.json" not in (out / "inputs_SUMS").read_text()
    # staged in place (no separate final directory): no second objects file
    out2 = tmp_path / "EGV1b"
    m.stage_egv1(zpath, out2, root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-26T00:00:00Z")
    assert not (out2 / "objects_staging.json").exists()

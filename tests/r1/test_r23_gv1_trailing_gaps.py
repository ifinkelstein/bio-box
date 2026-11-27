"""Card R-23: `GV1` re-staged with the trailing gap bytes of three contigs removed.

Card R-19's scan found the last two bytes of `GV1`'s contigs `Hoko_1`, `Hoko_3` and `Hoko_4` to be `--`, and card R-20
read the ingest path: a `-` is stored as a base, counted in the length and hashed into the `contig_uid`. The re-stage
reads card C-48c's staged `GV1` (its FASTA, `records.tsv`, `genomes.json` and the other tables) and writes it again under
a new FASTA name with the trailing `-` of those records removed and their uids and lengths recomputed; every other byte
of the FASTA and every other table row is copied as it is, and the staged files are only read.

The fixture is a small stage in C-48c's layout (headers and line widths as in the archive, `records.tsv` one row per
FASTA record in file order, `genomes.json` in `json.dumps` form without a final newline, a genome that repeats another's
contig set, a contig two genomes share). Every expected byte is written out by hand in this file, and the uids are
computed here with `blake3` directly (the definition of `contig_uid`), not with the function the script uses. Each test
runs twice, with the normal PATH and with PATH reduced to `/usr/bin:/bin`: the re-stage needs no pod tool, and GitHub
Actions has none of them.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("R-23")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
POD = REPO / "scripts" / "stress" / "build_stage_pod"
SCRIPT = POD / "r23_gv1_trailing_gaps.py"
TEMPLATE = POD / "stages" / "templates" / "GV1.conf"
LIVE_CONF = POD / "stages" / "GV1.conf"
NOTE = "planning/measurements/2026-11-27_non_base_bytes.md"
NEW_FASTA = "gvmags_v2_all_18kgenomes_r23.fna"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _mod():
    return _load("r23_gv1_trailing_gaps", SCRIPT)


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def uid(seq: bytes) -> str:
    return blake3(seq.upper()).digest(length=16).hex()


def guid(uids: list[str]) -> str:
    return blake3(b"".join(sorted(bytes.fromhex(u) for u in uids))).digest(length=16).hex()


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# The fixture's records, in archive order: (genome, member, accession, line layout as written). `Hoko_1` ends in `--`
# inside its last line, `Hoko_3`'s last line is `--` alone, `Hoko_4`'s gaps are split over two lines (`-` ends a full
# line, `-` alone after it). `Rep` repeats `Alpha`'s contig set (DEC-27) and `Beta` shares `a2` with `Alpha`.
SEQ = {
    "a1": b"ACGTACGTAC" b"GTTGCA",
    "a2": b"TTGACCAGTA" b"acgtn",
    "h1": b"GATTACAGAT" b"TACA--",
    "h2": b"CCCCGGGGAA" b"TT",
    "h3": b"AAACCCGGGT" b"--",
    "h4": b"GGGTTTAAA-" b"-",
    "b1": b"TTTTTTTTTT" b"GA",
}
OLD_FASTA = (
    b">a1 alpha contig one\n"
    b"ACGTACGTAC\n"
    b"GTTGCA\n"
    b">a2\n"
    b"TTGACCAGTA\n"
    b"acgtn\n"
    b">Hoko_1\n"
    b"GATTACAGAT\n"
    b"TACA--\n"
    b">Hoko_2\n"
    b"CCCCGGGGAA\n"
    b"TT\n"
    b">Hoko_3\n"
    b"AAACCCGGGT\n"
    b"--\n"
    b">Hoko_4\n"
    b"GGGTTTAAA-\n"
    b"-\n"
    b">a1 alpha contig one\n"
    b"ACGTACGTAC\n"
    b"GTTGCA\n"
    b">a2\n"
    b"TTGACCAGTA\n"
    b"acgtn\n"
    b">b1\n"
    b"TTTTTTTTTT\n"
    b"GA\n"
    b">a2\n"
    b"TTGACCAGTA\n"
    b"acgtn\n"
)
# the same file with the six gap bytes gone, written out by hand: `Hoko_3`'s line `--` and `Hoko_4`'s line `-` go with
# their line ends, and `Hoko_4`'s first line loses its last byte
NEW_FASTA_BYTES = (
    b">a1 alpha contig one\n"
    b"ACGTACGTAC\n"
    b"GTTGCA\n"
    b">a2\n"
    b"TTGACCAGTA\n"
    b"acgtn\n"
    b">Hoko_1\n"
    b"GATTACAGAT\n"
    b"TACA\n"
    b">Hoko_2\n"
    b"CCCCGGGGAA\n"
    b"TT\n"
    b">Hoko_3\n"
    b"AAACCCGGGT\n"
    b">Hoko_4\n"
    b"GGGTTTAAA\n"
    b">a1 alpha contig one\n"
    b"ACGTACGTAC\n"
    b"GTTGCA\n"
    b">a2\n"
    b"TTGACCAGTA\n"
    b"acgtn\n"
    b">b1\n"
    b"TTTTTTTTTT\n"
    b"GA\n"
    b">a2\n"
    b"TTGACCAGTA\n"
    b"acgtn\n"
)
TRIMMED = {"h1": b"GATTACAGATTACA", "h3": b"AAACCCGGGT", "h4": b"GGGTTTAAA"}
ROWS = [  # (genome, record, key, member)
    ("Alpha", 0, "a1", "confirmed/Alpha.fna"),
    ("Alpha", 1, "a2", "confirmed/Alpha.fna"),
    ("Hoko", 0, "h1", "confirmed/Hoko.fna"),
    ("Hoko", 1, "h2", "confirmed/Hoko.fna"),
    ("Hoko", 2, "h3", "confirmed/Hoko.fna"),
    ("Hoko", 3, "h4", "confirmed/Hoko.fna"),
    ("Rep", 0, "a1", "confirmed/Rep.fna"),
    ("Rep", 1, "a2", "confirmed/Rep.fna"),
    ("Beta", 0, "b1", "uncertain/Beta.fna"),
    ("Beta", 1, "a2", "uncertain/Beta.fna"),
]
ACC = {"a1": "a1", "a2": "a2", "h1": "Hoko_1", "h2": "Hoko_2", "h3": "Hoko_3", "h4": "Hoko_4", "b1": "b1"}
GENOMES = [("Alpha", ["a1", "a2"]), ("Hoko", ["h1", "h2", "h3", "h4"]), ("Beta", ["b1", "a2"])]


def _records_tsv(seqs: dict[str, bytes]) -> str:
    lines = ["ordinal\tgenome\trecord\taccession\tcontig_uid\tlength\tmember"]
    for k, (g, r, key, mem) in enumerate(ROWS):
        lines.append(f"{k}\t{g}\t{r}\t{ACC[key]}\t{uid(seqs[key])}\t{len(seqs[key])}\t{mem}")
    return "".join(ln + "\n" for ln in lines)


def _genomes_json(seqs: dict[str, bytes]) -> str:
    return json.dumps([{"contig_uids": [uid(seqs[k]) for k in ks], "genetic_code": 11, "kind": "virus",
                        "circular": [], "assembly_quality": "contig"} for _, ks in GENOMES])


def _write_stage(d: Path, fasta: bytes = OLD_FASTA, seqs: dict[str, bytes] | None = None,
                 repeated_extra: str = "", shared_extra: str = "") -> Path:
    seqs = dict(SEQ if seqs is None else seqs)
    d.mkdir(parents=True)
    (d / "gvmags_v2_all_18kgenomes.fna").write_bytes(fasta)
    (d / "records.tsv").write_text(_records_tsv(seqs))
    (d / "genomes.json").write_text(_genomes_json(seqs))
    (d / "objects.json").write_text(json.dumps([{
        "path": str(d / "gvmags_v2_all_18kgenomes.fna"), "source_name": "gvmags_v2",
        "source_release": "zenodo:10.5281/zenodo.20126600 GVMAGs V2, all_18kgenomes_fna.tar.gz (every genome)",
        "license": "CC BY-SA 4.0", "redistributable": True, "retrieved_at": "2026-11-19T10:42:14Z"}], indent=1))
    rep_guid = guid([uid(seqs["a1"]), uid(seqs["a2"])])
    (d / "repeated_records.tsv").write_text(
        "ordinal\tgenome\tgenome_uid\tfirst_genome\n" f"6\tRep\t{rep_guid}\tAlpha\n" + repeated_extra)
    (d / "shared_records.tsv").write_text(
        "ordinal\tgenome\trecord\tcontig_uid\towner_stage\n" f"8\tBeta\t0\t{uid(seqs['b1'])}\tV2\n" + shared_extra)
    (d / "skipped_members.tsv").write_text("member\tbytes\treason\nconfirmed/chunk_maker.sh~\t1001\tnot_fasta\n")
    return d


@pytest.fixture
def stage(tmp_path):
    old = _write_stage(tmp_path / "inputs" / "GV1")
    return {"old": old, "meta": tmp_path / "r23" / "GV1", "fasta": tmp_path / "stor" / "GV1" / NEW_FASTA,
            "final": Path("/tmp/if2453/s2/inputs/GV1"), "tmp": tmp_path}


def _run(m, st, **kw):
    return m.restage(st["old"], st["meta"], st["fasta"], st["final"], **kw)


# --- the re-stage ---------------------------------------------------------------------------------------------------

def test_the_trailing_gaps_go_and_every_other_byte_stays(stage):
    m = _mod()
    rep = _run(m, stage)
    assert stage["fasta"].read_bytes() == NEW_FASTA_BYTES
    new = dict(SEQ, **TRIMMED)
    assert (stage["meta"] / "records.tsv").read_text() == _records_tsv(new)
    assert (stage["meta"] / "genomes.json").read_text() == _genomes_json(new)  # same form, no final newline
    for n in ("repeated_records.tsv", "shared_records.tsv", "skipped_members.tsv"):
        assert (stage["meta"] / n).read_bytes() == (stage["old"] / n).read_bytes(), n
    objs = json.loads((stage["meta"] / "objects.json").read_text())
    stg = json.loads((stage["meta"] / "objects_staging.json").read_text())
    old_obj = json.loads((stage["old"] / "objects.json").read_text())[0]
    assert objs == [{**old_obj, "path": f"/tmp/if2453/s2/inputs/GV1/{NEW_FASTA}"}]
    assert stg == [{**old_obj, "path": str(stage["fasta"])}]
    assert (stage["meta"] / "objects.json").read_text() == json.dumps(objs, indent=1)  # C-48c's form
    gaps = (stage["meta"] / "gap_trims.tsv").read_text().splitlines()
    assert gaps == [
        "ordinal\tgenome\taccession\told_contig_uid\told_length\tremoved\tnew_contig_uid\tnew_length",
        f"2\tHoko\tHoko_1\t{uid(SEQ['h1'])}\t16\t--\t{uid(TRIMMED['h1'])}\t14",
        f"4\tHoko\tHoko_3\t{uid(SEQ['h3'])}\t12\t--\t{uid(TRIMMED['h3'])}\t10",
        f"5\tHoko\tHoko_4\t{uid(SEQ['h4'])}\t11\t--\t{uid(TRIMMED['h4'])}\t9",
    ]
    # counts, by hand: 10 records, 7 distinct contigs, 3 genomes listed (Rep is a repeat of Alpha)
    assert rep["old"]["records"] == rep["new"]["records"] == 10
    old_bases = sum(len(SEQ[k]) for _, _, k, _ in ROWS)
    assert rep["old"]["record_bases"] == old_bases and rep["new"]["record_bases"] == old_bases - 6
    assert rep["old"]["contigs"] == rep["new"]["contigs"] == 7
    distinct_old = sum(len(s) for s in SEQ.values())
    assert rep["old"]["contig_bases"] == distinct_old and rep["new"]["contig_bases"] == distinct_old - 6
    assert rep["old"]["genomes"] == rep["new"]["genomes"] == 3
    assert rep["old"]["memberships"] == rep["new"]["memberships"] == 8
    assert rep["trimmed_records"] == 3 and rep["removed_gap_bytes"] == 6
    assert rep["fasta_bytes_removed"] == len(OLD_FASTA) - len(NEW_FASTA_BYTES) == 8  # 2 + 3 + 1 + 2
    assert rep["other_records_byte_identical"] is True and rep["uid_mismatches"] == 0
    assert rep["new_uids_already_in_stage"] == []
    assert rep["old_fasta_sha256"] == hashlib.sha256(OLD_FASTA).hexdigest()
    assert rep["fasta_sha256"] == hashlib.sha256(NEW_FASTA_BYTES).hexdigest()
    # the removed bytes' places in the old file: the two bytes before `Hoko_1`'s last line end; `Hoko_3`'s line `--\n`;
    # `Hoko_4`'s last byte of its first line, then its line `-\n`
    o_h1 = OLD_FASTA.index(b"TACA--\n") + 4
    o_h3 = OLD_FASTA.index(b"AAACCCGGGT\n--\n") + 11
    o_h4 = OLD_FASTA.index(b"GGGTTTAAA-\n") + 9
    assert rep["removed_ranges"] == [[o_h1, 2], [o_h3, 3], [o_h4, 1], [o_h4 + 2, 2]]
    keep, at = bytearray(), 0
    for off, n in rep["removed_ranges"]:
        keep += OLD_FASTA[at:off]
        at = off + n
    keep += OLD_FASTA[at:]
    assert bytes(keep) == NEW_FASTA_BYTES


def test_the_new_file_scans_clean_and_the_old_one_does_not(stage):
    scan = _load("r19_fasta_scan", POD / "r19_fasta_scan.py").scan_file
    _run(_mod(), stage)
    old = scan(stage["old"] / "gvmags_v2_all_18kgenomes.fna")
    assert old["status"] == "damaged" and old["records_with_non_iupac"] == 3 and old["non_iupac_bytes"] == 6
    new = scan(stage["fasta"])
    assert new["status"] == "clean" and new["records"] == 10 and new["non_iupac_bytes"] == 0
    assert new["bases"] == sum(len(SEQ[k]) for _, _, k, _ in ROWS) - 6


def test_the_staged_files_are_only_read_and_a_second_run_is_byte_identical(stage):
    m = _mod()
    before = {p.name: sha(p) for p in stage["old"].iterdir()}
    _run(m, stage)
    second = {"old": stage["old"], "meta": stage["tmp"] / "again" / "GV1",
              "fasta": stage["tmp"] / "again" / NEW_FASTA, "final": stage["final"]}
    _run(m, second)
    assert {p.name: sha(p) for p in stage["old"].iterdir()} == before
    assert stage["fasta"].read_bytes() == second["fasta"].read_bytes()
    for p in stage["meta"].iterdir():
        if p.name in ("objects_staging.json", "stage_report.json"):
            continue  # they name the run's own paths
        assert p.read_bytes() == (second["meta"] / p.name).read_bytes(), p.name
    assert not list(stage["meta"].glob(".*")) and not list(stage["fasta"].parent.glob("*.tmp"))


def test_a_new_uid_an_earlier_stage_holds_is_reported(stage):
    m = _mod()
    chain = stage["tmp"] / "chain"
    (chain / "V2").mkdir(parents=True)
    (chain / "V2" / "genomes.json").write_text(json.dumps([{"contig_uids": [uid(TRIMMED["h3"]), uid(SEQ["b1"])]}]))
    rep = _run(m, stage, chain_home=chain, work=stage["tmp"] / "work")
    assert rep["ownership_old"]["held_contigs"] == {"V2": 1}  # b1
    assert rep["ownership_new"]["held_contigs"] == {"V2": 2}  # b1 and the trimmed Hoko_3 (DEC-26)
    assert rep["trimmed_uids_held"] == {uid(TRIMMED["h3"]): "V2"}


def test_a_trimmed_contig_equal_to_another_contig_of_its_genome_is_refused(stage, tmp_path):
    m = _mod()
    seqs = dict(SEQ, h3=b"CCCCGGGGAA" b"TT--")  # Hoko_3 without its gaps is Hoko_2
    fasta = OLD_FASTA.replace(b">Hoko_3\nAAACCCGGGT\n--\n", b">Hoko_3\nCCCCGGGGAA\nTT--\n")
    old = _write_stage(tmp_path / "x" / "GV1", fasta=fasta, seqs=seqs)
    with pytest.raises(m.Refused, match="Hoko_3.*twice in genome Hoko"):
        m.restage(old, tmp_path / "x" / "meta", tmp_path / "x" / "new.fna", stage["final"])


def test_a_trimmed_contig_equal_to_a_contig_of_another_genome_is_reported(stage, tmp_path):
    m = _mod()
    seqs = dict(SEQ, h1=b"TTTTTTTTTT" b"GA--")  # Hoko_1 without its gaps is Beta's b1 (DEC-27: listed in both)
    fasta = OLD_FASTA.replace(b">Hoko_1\nGATTACAGAT\nTACA--\n", b">Hoko_1\nTTTTTTTTTT\nGA--\n")
    old = _write_stage(tmp_path / "y" / "GV1", fasta=fasta, seqs=seqs)
    rep = m.restage(old, tmp_path / "y" / "meta", tmp_path / "y" / "new.fna", stage["final"])
    assert rep["new_uids_already_in_stage"] == [uid(SEQ["b1"])]
    assert rep["old"]["contigs"] == 7 and rep["new"]["contigs"] == 6
    assert rep["new"]["memberships"] == 8 and rep["new"]["genomes"] == 3


def test_a_trimmed_genome_equal_to_an_earlier_genome_becomes_its_repeat(tmp_path):
    """GV1's case: `Hoko` without its gaps is the GenBank assembly staged earlier in the archive. By C-48c's DEC-27
    rule the later genome is the earlier one: out of genomes.json, into repeated_records.tsv; its records stay."""
    m = _mod()
    gb1, gb2, hk1 = b"ACGTTGCAAC" b"GT", b"TTTTGGGGCC" b"A", b"ACGTTGCAAC" b"GT--"
    fasta = (b">Gb|K1__1\nACGTTGCAAC\nGT\n>Gb|K2__1\nTTTTGGGGCC\nA\n>x1\nCCCCCCCCCC\n"
             b">Hk_1\nACGTTGCAAC\nGT--\n>Hk_2\nTTTTGGGGCC\nA\n")
    rows = [("Gb", 0, "Gb|K1__1", gb1), ("Gb", 1, "Gb|K2__1", gb2), ("X", 0, "x1", b"CCCCCCCCCC"),
            ("Hk", 0, "Hk_1", hk1), ("Hk", 1, "Hk_2", gb2)]
    d = tmp_path / "in" / "GV1"
    d.mkdir(parents=True)
    (d / "g.fna").write_bytes(fasta)
    (d / "records.tsv").write_text("ordinal\tgenome\trecord\taccession\tcontig_uid\tlength\tmember\n" + "".join(
        f"{k}\t{g}\t{r}\t{a}\t{uid(s)}\t{len(s)}\tconfirmed/{g}.fna\n" for k, (g, r, a, s) in enumerate(rows)))
    ent = lambda ss: {"contig_uids": [uid(s) for s in ss], "genetic_code": 11, "kind": "virus",  # noqa: E731
                      "circular": [], "assembly_quality": "contig"}
    (d / "genomes.json").write_text(json.dumps([ent([gb1, gb2]), ent([b"CCCCCCCCCC"]), ent([hk1, gb2])]))
    (d / "objects.json").write_text(json.dumps([{"path": str(d / "g.fna"), "source_name": "gvmags_v2"}], indent=1))
    rep_old = "ordinal\tgenome\tgenome_uid\tfirst_genome\n" "9\tLater\tab\tX\n"
    (d / "repeated_records.tsv").write_text(rep_old)
    (d / "shared_records.tsv").write_text("ordinal\tgenome\trecord\tcontig_uid\towner_stage\n")
    (d / "skipped_members.tsv").write_text("member\tbytes\treason\n")
    meta, fa = tmp_path / "meta", tmp_path / "new.fna"
    rep = m.restage(d, meta, fa, Path("/tmp/if2453/s2/inputs/GV1"), expect=("Hk_1",))
    assert fa.read_bytes() == fasta.replace(b"GT--\n", b"GT\n")
    assert json.loads((meta / "genomes.json").read_text()) == [ent([gb1, gb2]), ent([b"CCCCCCCCCC"])]
    assert (meta / "repeated_records.tsv").read_text() == (
        "ordinal\tgenome\tgenome_uid\tfirst_genome\n" f"3\tHk\t{guid([uid(gb1), uid(gb2)])}\tGb\n" "9\tLater\tab\tX\n")
    assert rep["genomes_now_repeats"] == [{"ordinal": 3, "genome": "Hk", "genome_uid": guid([uid(gb1), uid(gb2)]),
                                           "first_genome": "Gb"}]
    assert rep["new_uids_already_in_stage"] == [uid(gb1)]
    assert (rep["old"]["genomes"], rep["new"]["genomes"]) == (3, 2)
    assert (rep["old"]["contigs"], rep["new"]["contigs"]) == (4, 3)
    assert rep["new"]["records"] == 5 and rep["new"]["memberships"] == 3
    assert rep["other_records_byte_identical"] is True and rep["uid_mismatches"] == 0


# --- refusals ------------------------------------------------------------------------------------------------------

def _refused(m, tmp_path, match, **kw):
    old = _write_stage(tmp_path / "bad" / "GV1", **kw)
    meta, fa = tmp_path / "bad" / "meta", tmp_path / "bad" / "new.fna"
    with pytest.raises(m.Refused, match=match):
        m.restage(old, meta, fa, Path("/tmp/if2453/s2/inputs/GV1"))
    assert not fa.exists() and not fa.with_name(fa.name + ".tmp").exists()


def test_an_inner_gap_is_refused(tmp_path):
    seqs = dict(SEQ, a1=b"ACGTACG-AC" b"GTTGCA")
    _refused(_mod(), tmp_path, r"record 1 \(a1\).*not only trailing", seqs=seqs,
             fasta=OLD_FASTA.replace(b"ACGTACGTAC\nGTTGCA", b"ACGTACG-AC\nGTTGCA"))


def test_another_byte_that_is_not_a_base_is_refused(tmp_path):
    seqs = dict(SEQ, b1=b"TTTTTTTTTT" b"G*")
    _refused(_mod(), tmp_path, r"record 9 \(b1\)", seqs=seqs, fasta=OLD_FASTA.replace(b"TTTTTTTTTT\nGA", b"TTTTTTTTTT\nG*"))


def test_a_records_table_out_of_step_with_the_fasta_is_refused(tmp_path):
    seqs = dict(SEQ, h2=b"CCCCGGGGAA" b"TA")  # records.tsv says TA, the FASTA holds TT
    _refused(_mod(), tmp_path, r"record 4 \(Hoko_2\).*records.tsv", seqs=seqs)


def test_a_trimmed_genome_in_the_repeats_table_is_refused(tmp_path):
    _refused(_mod(), tmp_path, "Hoko.*repeated_records.tsv", repeated_extra="99\tHokoCopy\tff\tHoko\n")


def test_a_trimmed_contig_in_the_shared_table_is_refused(tmp_path):
    _refused(_mod(), tmp_path, "shared_records.tsv", shared_extra=f"2\tHoko\t0\t{uid(SEQ['h1'])}\tV2\n")


def test_unexpected_trimmed_records_are_refused(stage):
    m = _mod()
    with pytest.raises(m.Refused, match="expected.*Hoko_9"):
        _run(m, stage, expect=("Hoko_1", "Hoko_3", "Hoko_9"))
    assert not stage["fasta"].exists()


def test_existing_targets_are_refused(stage):
    m = _mod()
    _run(m, stage)
    with pytest.raises(m.Refused, match="exists"):
        _run(m, stage)
    with pytest.raises(m.Refused, match="not empty"):
        m.restage(stage["old"], stage["meta"], stage["tmp"] / "other.fna", stage["final"])


def test_the_command_line(stage):
    p = subprocess.run([sys.executable, str(SCRIPT), "restage", str(stage["old"]), str(stage["meta"]),
                        str(stage["fasta"]), str(stage["final"])], capture_output=True, text=True, check=False)
    assert p.returncode == 0, p.stderr
    rep = json.loads(p.stdout)
    assert rep["card"] == "R-23" and rep["trimmed_records"] == 3
    p2 = subprocess.run([sys.executable, str(SCRIPT), "restage", str(stage["old"]), str(stage["meta"]),
                         str(stage["fasta"]), str(stage["final"])], capture_output=True, text=True, check=False)
    assert p2.returncode == 2 and "refused" in p2.stderr


# --- the template ---------------------------------------------------------------------------------------------------

def test_the_template_names_the_new_files_and_keeps_the_stages_values():
    est = _load("stage_est", POD / "stage_est.py")
    c = est.read_conf(TEMPLATE, "GV1", tree=REPO)
    live = est.read_conf(LIVE_CONF, "GV1", tree=REPO)
    assert str(c.inputs) == "/tmp/if2453/s2/inputs/GV1"
    assert (c.depth, c.library, c.crispr, c.min_contig_len, c.genecall_mode) == \
        (live.depth, live.library, live.crispr, live.min_contig_len, live.genecall_mode)
    assert c.after == live.after == ("V2", "PL3", "GART2")
    assert {k: v[0] for k, v in c.judge.items()} == {k: v[0] for k, v in live.judge.items()}
    assert {k: v[0] for k, v in c.expect.items()} == {k: v[0] for k, v in live.expect.items()}
    assert NOTE in c.notes
    assert {n for _, n in c.sums} == {"genomes.json", "objects.json", NEW_FASTA}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    for h, n in c.sums:
        assert h in note, f"SUM of {n} is not in the note"
    assert not set(c.sums) & set(live.sums), "the template's files are new; the staged ones stay as they are"
    for st, (n, src) in c.expect.items():
        assert src == NOTE and (f"{n:,}" in note or f" {n} " in note), st

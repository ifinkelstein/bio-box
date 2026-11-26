"""Card C-79 (DEC-79): MetaVR's `>= 30 kb` spillover (`MV2`) and its eukaryotic-virus slice (`MVEUKV`).

Oracles are written here separately from the script: the selection by a brute-force loop over shuffled rows,
contig_uid by blake3 directly. The selection logic needs no native tool; only the skani test skips without skani.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-79")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-26_mv2_mveukv_inputs.md"


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
    return _load("c79_metavr_mv2_eukv", "c79_metavr_mv2_eukv.py")


def cuid(seq: str) -> str:
    return blake3(seq.upper().encode()).digest(length=16).hex()


def seq_of(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


PRE = "r__Varidnaviria;k__Bamfordvirae;"
NC = PRE + "p__Nucleocytoviricota;c__Megaviricetes;o__Imitervirales;f__Mimiviridae;g__;s__"
MAV = PRE + "p__Preplasmiviricota;c__Maveriviricetes;o__Priklausovirales;f__;g__;s__"
POL = PRE + "p__Preplasmiviricota;c__Polintoviricetes;o__Orthopolintovirales;f__Adintoviridae;g__;s__"
TEC = PRE + "p__Preplasmiviricota;c__Tectiliviricetes;o__Kalamavirales;f__Tectiviridae;g__;s__"
PHAGE = "r__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;c__Caudoviricetes;o__;f__;g__;s__"
DECOY = "r__;k__;p__;c__;o__;f__;g__;s__Nucleocytoviricota-like virus"
QRANK = {"Complete": 0, "High-quality": 1, "Medium-quality": 2, "Low-quality": 3, "Not-determined": 4}
TAXON_OF = (("Nucleocytoviricota", "p__Nucleocytoviricota"), ("Maveriviricetes", "c__Maveriviricetes"),
            ("Polintoviricetes", "c__Polintoviricetes"))


def row(uvig, votu, tax, *, length, quality="Low-quality", conf="High", comp=20.0, contam=0.0,
        topology="Linear", concat=False):
    return {"uvig": uvig, "taxon_oid": 1, "votu": votu, "quality": quality, "completeness": comp,
            "checkv_contamination": contam, "length": length, "uvig_topology": topology, "is_concatemer": concat,
            "viral_confidence": conf, "genetic_code": "11", "genome_type": "dsDNA", "source": "Metagenome",
            "Ecosystem": "Environmental", "host_taxonomy": "", "ictv_taxonomy": tax}


SCHEMA = pa.schema([
    ("uvig", pa.string()), ("taxon_oid", pa.int64()), ("votu", pa.string()), ("quality", pa.string()),
    ("completeness", pa.float64()), ("checkv_contamination", pa.float64()), ("length", pa.int64()),
    ("uvig_topology", pa.string()), ("is_concatemer", pa.bool_()), ("viral_confidence", pa.string()),
    ("genetic_code", pa.string()), ("genome_type", pa.string()), ("source", pa.string()), ("Ecosystem", pa.string()),
    ("host_taxonomy", pa.string()), ("ictv_taxonomy", pa.string()),
])
LEFT_HEAD = ["reason", "mv2_rank", "uvig", "votu", "quality", "viral_confidence", "length", "completeness",
             "checkv_contamination", "uvig_topology", "is_concatemer", "genetic_code"]


def write_parquet(path: Path, rows: list[dict]) -> None:
    pq.write_table(pa.Table.from_pylist(rows, schema=SCHEMA), path)


def write_left(path: Path, left: list[tuple[str, int, dict]]) -> None:
    """C-62's MV1_left_for_later.tsv.gz form: (reason, mv2_rank, row)."""
    buf = io.StringIO()
    buf.write("\t".join(LEFT_HEAD) + "\n")
    for reason, rank, r in left:
        buf.write("\t".join(str(x) for x in (
            reason, rank, r["uvig"], r["votu"], r["quality"], r["viral_confidence"], r["length"], r["completeness"],
            r["checkv_contamination"], r["uvig_topology"], str(r["is_concatemer"]).lower(), r["genetic_code"])) + "\n")
    path.write_bytes(gzip.compress(buf.getvalue().encode(), mtime=0))


# The rule fixture. MV1 took M_JUMBO (absent from the left list); MV2 takes S1..S3 by rank.
R = {r["uvig"]: r for r in [
    row("M_JUMBO", "vA", NC, length=250_000, quality="High-quality"),
    row("A2", "vA", NC, length=8_000),  # vOTU touched by MV1
    row("S1", "vS1", PHAGE, length=40_000),
    row("S2", "vS2", NC, length=35_000),
    row("S3", "vS3", MAV, length=31_000),
    row("S2b", "vS2", NC, length=9_000),  # vOTU touched by MV2
    row("C1", "vC", NC, length=12_000),
    row("D1", "vD", NC, length=7_000, comp=30.0),
    row("D2", "vD", NC, length=9_000, comp=20.0),
    row("E1", "vE", NC, length=20_000, contam=15.0),  # high contamination: out; E2 represents vE
    row("E2", "vE", NC, length=6_000, contam=2.0),
    row("F1", "vF", NC, length=8_000, conf="Low"),
    row("G1", "vG", NC, length=3_000),
    row("H1", "vH", NC, length=9_000, concat=True),
    row("I1", "vI", MAV, length=6_000),
    row("J1", "vJ", POL, length=5_000),  # the floor is inclusive
    row("K1", "vK", POL, length=4_999),
    row("L1", "vL", TEC, length=10_000),  # Tectiliviricetes: phages and adenoviruses, not taken
    row("N1", "vN", PHAGE, length=10_000),
    row("O1", "vO", DECOY, length=10_000),
    row("P1", "vP", None, length=10_000),
    row("Q1", "vQ", NC, length=6_000, comp=10.0),  # qualifies; Q2 is the better row but of low confidence
    row("Q2", "vQ", NC, length=20_000, comp=50.0, conf="Low"),
    row("T1", "vT", NC, length=6_500, quality="Not-determined"),
]}
MV2_RANK = {"S2": 1, "S1": 2, "S3": 3}
LEFT = [("mv2_spillover" if u in MV2_RANK else "lq_nd_lt30k", MV2_RANK.get(u, 0), r)
        for u, r in R.items() if u != "M_JUMBO"]


def taxon_of(tax):
    toks = (tax or "").split(";")
    for name, tok in TAXON_OF:
        if tok in toks:
            return name
    return None


def brute(rows: list[dict], left: list[tuple[str, int, dict]]):
    """Independent oracle: MV2 in rank order; MVEUKV one best row per untouched vOTU among the qualifying rows."""
    in_left = {r["uvig"] for _, _, r in left}
    mv2 = sorted(((rank, r["uvig"]) for reason, rank, r in left if reason == "mv2_spillover"))
    mv2_ids = [u for _, u in mv2]
    touched = {r["votu"] for r in rows if r["uvig"] not in in_left or r["uvig"] in mv2_ids}

    def key(r):
        return (QRANK.get(r["quality"], 5), -min(r["completeness"] if r["completeness"] is not None else -1.0, 100.0),
                r["checkv_contamination"] if r["checkv_contamination"] is not None else 999.0, -r["length"], r["uvig"])

    groups: dict[str, list[dict]] = {}
    for r in rows:
        t = taxon_of(r["ictv_taxonomy"])
        if t is None or r["votu"] in touched:
            continue
        groups.setdefault(r["votu"], [])
        if r["is_concatemer"] or (r["checkv_contamination"] or 0.0) > 10.0:
            continue
        groups[r["votu"]].append(r)
    sel, left_out = [], {}
    for v, grp in groups.items():
        ok = sorted((r for r in grp if r["viral_confidence"] == "High" and r["length"] >= 5000), key=key)
        if ok:
            sel.append(ok[0])
        elif not grp:
            left_out[v] = "no_eligible_member"
        else:
            best = sorted(grp, key=key)[0]
            left_out[v] = "lt5kb" if best["length"] < 5000 else "low_confidence_ge5kb"
    order = [n for n, _ in TAXON_OF]
    sel.sort(key=lambda r: (order.index(taxon_of(r["ictv_taxonomy"])),) + key(r))
    return mv2_ids, [r["uvig"] for r in sel], left_out


def test_selection_mv2_in_rank_order_and_mveukv_rule_independent_of_row_order(tmp_path):
    mv = _mod()
    rows = list(R.values())
    want_mv2, want_eukv, want_left = brute(rows, LEFT)
    assert want_mv2 == ["S2", "S1", "S3"]
    assert sorted(want_eukv) == ["C1", "D1", "E2", "I1", "J1", "Q1", "T1"]
    assert want_left == {"vF": "low_confidence_ge5kb", "vG": "lt5kb", "vH": "no_eligible_member", "vK": "lt5kb"}
    left_gz = tmp_path / "left.tsv.gz"
    write_left(left_gz, LEFT)
    rng = random.Random(7901)
    for i in range(6):
        shuffled = rows[:]
        rng.shuffle(shuffled)
        p = tmp_path / f"t{i}.parquet"
        write_parquet(p, shuffled)
        mv2, eukv, left, stats = mv.select_lists(p, left_gz, work=tmp_path / f"w{i}")
        assert [r["uvig"] for r in mv2] == want_mv2
        assert [r["mv2_rank"] for r in mv2] == [1, 2, 3]
        assert [r["uvig"] for r in eukv] == want_eukv
        assert [r["rank"] for r in eukv] == list(range(1, len(want_eukv) + 1))
        assert {r["votu"]: r["reason"] for r in left} == want_left
        assert {r["uvig"]: r["taxon"] for r in eukv}["I1"] == "Maveriviricetes"
        assert stats["mv2_uvigs"] == 3 and stats["mv2_bases"] == 106_000
        assert stats["mveukv_uvigs"] == 7
        assert stats["mveukv_bases"] == sum(R[u]["length"] for u in want_eukv)


def test_selection_refuses_a_left_list_that_does_not_match_the_table(tmp_path):
    mv = _mod()
    p = tmp_path / "t.parquet"
    write_parquet(p, list(R.values()))
    stray = tmp_path / "stray.tsv.gz"
    write_left(stray, LEFT + [("lq_nd_lt30k", 0, row("ZZ", "vZ", NC, length=6_000))])
    with pytest.raises(mv.Refused, match="not in the table"):
        mv.select_lists(p, stray, work=tmp_path / "w1")
    gap = tmp_path / "gap.tsv.gz"
    write_left(gap, [(reason, 4 if r["uvig"] == "S3" else rank, r) for reason, rank, r in LEFT])
    with pytest.raises(mv.Refused, match="mv2_rank"):
        mv.select_lists(p, gap, work=tmp_path / "w2")


# --- staging end to end -------------------------------------------------------------------------------------------

def _fa(records: list[tuple[str, str]]) -> bytes:
    out = io.StringIO()
    for hdr, s in records:
        out.write(f">{hdr}\n")
        for i in range(0, len(s), 70):
            out.write(s[i:i + 70] + "\n")
    return gzip.compress(out.getvalue().encode(), mtime=0)


def _stage_fixture(tmp_path: Path):
    unit = seq_of(11, 34_960)
    seqs = {
        "S2": unit + unit[:40],  # MV2, exact 40 bp DTR
        "S1": seq_of(12, 40_000),
        "S3": seq_of(13, 31_000),
        "C1a": seq_of(14, 7_000), "C1b": seq_of(15, 5_000),  # MVEUKV Virus MAG of two records
        "D1": seq_of(16, 7_000),
        "I1": seq_of(17, 6_000),
        "J1": seq_of(18, 5_000),
        "X": seq_of(19, 3_000),  # not selected
    }
    eukv_gv1 = seq_of(20, 8_000)  # also a GV1 contig
    eukv_v2 = seq_of(21, 9_000)  # also a contig of an r26.11.4 stage (V2)
    rows = [
        row("S2", "vS2", NC, length=len(seqs["S2"]), topology="Direct terminal repeat"),
        row("S1", "vS1", PHAGE, length=40_000),
        row("S3", "vS3", MAV, length=31_000),
        row("C1", "vC", NC, length=12_000, topology="Virus MAG"),
        row("D1", "vD", NC, length=7_000),
        row("I1", "vI", MAV, length=6_000),
        row("J1", "vJ", POL, length=5_000),
        row("GV", "vGV", NC, length=8_000),
        row("V2C", "vV2", NC, length=9_000),
        row("DUP", "vDUP", NC, length=6_000),  # its sequence is I1's: one genome listed, the other repeated (DEC-27)
        row("X1", "vX", PHAGE, length=3_000),
        row("MJ", "vMJ", NC, length=3_000),  # MV1 took it
    ]
    seqs["DUP"] = seqs["I1"]
    left = []
    for r in rows:
        if r["uvig"] == "MJ":
            continue
        rank = {"S2": 1, "S1": 2, "S3": 3}.get(r["uvig"], 0)
        left.append(("mv2_spillover" if rank else "lq_nd_lt30k", rank, r))
    pq_path = tmp_path / "METAVR_main_table.parquet"
    write_parquet(pq_path, rows)
    left_gz = tmp_path / "left.tsv.gz"
    write_left(left_gz, left)
    fna = tmp_path / "METAVR.fna.bgz"
    fna.write_bytes(_fa([
        ("MJ|1|s0", seq_of(22, 3_000)),
        ("S1|1|s1", seqs["S1"]), ("C1|1|c0", seqs["C1a"]), ("C1|1|c1", seqs["C1b"]), ("S2|1|s2", seqs["S2"]),
        ("D1|1|s3", seqs["D1"]), ("GV|1|s4", eukv_gv1), ("V2C|1|s5", eukv_v2), ("I1|1|s6", seqs["I1"]),
        ("DUP|1|s7", seqs["DUP"]), ("J1|1|s8", seqs["J1"]), ("S3|1|s9", seqs["S3"]), ("X1|1|s10", seqs["X"]),
    ]))
    v2 = tmp_path / "V2"
    v2.mkdir()
    (v2 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(eukv_v2)]}]))
    gv1 = tmp_path / "GV1"
    gv1.mkdir()
    (gv1 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(eukv_gv1), cuid(seqs["S1"])]}]))
    return pq_path, left_gz, fna, v2, gv1, seqs, unit, eukv_gv1, eukv_v2


def _run_stage(mv, tmp_path, tag, pq_path, left_gz, fna, v2, gv1):
    return mv.stage(
        pq_path, left_gz, fna, tmp_path / f"meta_{tag}", tmp_path / f"fasta_{tag}",
        root=None, root_stages=(), owner_inputs={"V2": v2}, later_inputs={"GV1": gv1},
        inputs_home=Path("/tmp/if2453/s2/inputs"), retrieved_at="2026-10-23T15:22:00Z", work=tmp_path / f"work_{tag}",
    )


def test_stage_writes_both_stages_from_one_stream_deterministically(tmp_path):
    mv = _mod()
    pq_path, left_gz, fna, v2, gv1, seqs, unit, eukv_gv1, eukv_v2 = _stage_fixture(tmp_path)
    rep1 = _run_stage(mv, tmp_path, "a", pq_path, left_gz, fna, v2, gv1)
    rep2 = _run_stage(mv, tmp_path, "b", pq_path, left_gz, fna, v2, gv1)
    m2, me = tmp_path / "meta_a" / "MV2", tmp_path / "meta_a" / "MVEUKV"
    f2, fe = tmp_path / "fasta_a" / "MV2" / "mv2.fna", tmp_path / "fasta_a" / "MVEUKV" / "mveukv.fna"
    assert rep1["streams"] == 1

    g2 = json.loads((m2 / "genomes.json").read_text())
    assert [g["contig_uids"] for g in g2] == [[cuid(unit)], [cuid(seqs["S1"])], [cuid(seqs["S3"])]]  # mv2_rank order
    assert g2[0]["circular"] == [cuid(unit)]
    assert {g["kind"] for g in g2} == {"virus"} and {g["genetic_code"] for g in g2} == {11}
    names2 = [ln[1:].split("|")[0] for ln in f2.read_text().splitlines() if ln.startswith(">")]
    assert sorted(names2) == ["S1", "S2", "S3"]

    ge = json.loads((me / "genomes.json").read_text())
    # rank order: Nucleocytoviricota (C1 12 kb, V2C 9 kb, GV 8 kb, D1 7 kb, DUP 6 kb by length), then
    # Maveriviricetes (I1, which repeats DUP's contig and is not listed again), then Polintoviricetes (J1)
    assert [g["contig_uids"] for g in ge] == [
        [cuid(seqs["C1a"]), cuid(seqs["C1b"])], [cuid(eukv_v2)], [cuid(eukv_gv1)], [cuid(seqs["D1"])],
        [cuid(seqs["I1"])], [cuid(seqs["J1"])],
    ]
    rep_e = rep1["stages"]["MVEUKV"]
    assert rep_e["repeated_genomes"] == 1
    assert rep_e["owner_contigs"] == {"V2": 1}
    assert rep_e["later_contigs"] == {"GV1": 1}
    rep_2 = rep1["stages"]["MV2"]
    assert rep_2["later_contigs"] == {"GV1": 1} and rep_2["trims"] == 1 and rep_2["trimmed_bases"] == 40
    shared = list(csv.DictReader(open(me / "shared_records.tsv"), delimiter="\t"))
    assert [(r["genome"], r["owner_stage"]) for r in shared] == [("V2C", "V2")]
    later = list(csv.DictReader(open(me / "later_records.tsv"), delimiter="\t"))
    assert [(r["genome"], r["holder"]) for r in later] == [("GV", "GV1")]

    objs = json.loads((me / "objects.json").read_text())
    assert objs[0]["path"] == "/tmp/if2453/s2/inputs/MVEUKV/mveukv.fna"
    assert objs[0]["redistributable"] is False and objs[0]["source_name"] == "metavr"
    assert json.loads((me / "objects_staging.json").read_text())[0]["path"] == str(fe)
    left = gzip.decompress((me / "MVEUKV_left_for_later.tsv.gz").read_bytes()).decode().splitlines()
    assert len(left) == 1  # header only: every eukaryotic-virus vOTU of the fixture is taken

    for stage, files in (("MV2", ("genomes.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
                                  "later_records.tsv", "trims.tsv", "objects.json")),
                         ("MVEUKV", ("genomes.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv",
                                     "later_records.tsv", "trims.tsv", "objects.json",
                                     "MVEUKV_left_for_later.tsv.gz"))):
        for f in files:
            a = (tmp_path / "meta_a" / stage / f).read_bytes()
            b = (tmp_path / "meta_b" / stage / f).read_bytes()
            assert a == b, f"{stage}/{f} differs between runs"
    assert f2.read_bytes() == (tmp_path / "fasta_b" / "MV2" / "mv2.fna").read_bytes()
    assert fe.read_bytes() == (tmp_path / "fasta_b" / "MVEUKV" / "mveukv.fna").read_bytes()
    assert rep1["stages"]["MVEUKV"]["sha256"]["genomes.json"] == rep2["stages"]["MVEUKV"]["sha256"]["genomes.json"]

    for stage, meta in (("MV2", m2), ("MVEUKV", me)):
        chk = subprocess.run([sys.executable, str(BSP / "s2_check_inputs.py"), stage, str(meta / "genomes.json"),
                              str(meta / "objects_staging.json")], capture_output=True, text=True, check=False)
        assert chk.returncode == 0, chk.stdout + chk.stderr
        assert json.loads(chk.stdout)["problems"] == []


def test_stage_refuses_a_selected_uvig_missing_from_the_fasta_and_writes_nothing(tmp_path):
    mv = _mod()
    pq_path, left_gz, fna, v2, gv1, *_ = _stage_fixture(tmp_path)
    lines = gzip.decompress(fna.read_bytes()).decode().split(">")
    fna.write_bytes(gzip.compress(">".join(x for x in lines if not x.startswith("J1|")).encode(), mtime=0))
    with pytest.raises(mv.Refused, match="no FASTA record"):
        _run_stage(mv, tmp_path, "a", pq_path, left_gz, fna, v2, gv1)
    assert not (tmp_path / "meta_a" / "MVEUKV").exists()
    assert not (tmp_path / "fasta_a" / "MVEUKV" / "mveukv.fna").exists()


# --- the ANI sample -----------------------------------------------------------------------------------------------

def test_ani_summary_counts_the_vOTU_threshold_by_taxon(tmp_path):
    mv = _mod()
    table = tmp_path / "a.tsv"
    table.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "r\tq\t99.0\t10.0\t95.0\tgv_big\tq1\n"   # q1 (8 kb) inside a 300 kb GV1 contig: meets `shorter`
        "r\tq\t94.0\t90.0\t90.0\tgv_a\tq2\n"     # ANI under 95: no
        "r\tq\t97.0\t40.0\t50.0\tgv_a\tq3\n"     # AF of the shorter under 85: no
    )
    sample = {"q1": 8_000, "q2": 9_000, "q3": 20_000, "q4": 6_000}
    taxon = {"q1": "Nucleocytoviricota", "q2": "Nucleocytoviricota", "q3": "Maveriviricetes",
             "q4": "Polintoviricetes"}
    refs = {"gv_big": 300_000, "gv_a": 30_000}
    s = mv.ani_summary(table, sample, taxon, refs)
    assert s["queries"] == 4
    assert s["shorter"]["meet"] == 1
    assert s["shorter"]["by_taxon"]["Nucleocytoviricota"] == {"queries": 2, "meet": 1}
    assert s["both"]["meet"] == 0


@pytest.mark.skipif(shutil.which("skani") is None, reason="skani is not on PATH (the pod's image has it)")
def test_ani_sample_finds_a_copy_of_a_reference_with_skani(tmp_path):
    mv = _mod()
    ref = seq_of(31, 30_000)
    q_same, q_other = ref[2_000:22_000], seq_of(32, 20_000)
    meta = tmp_path / "MVEUKV"
    meta.mkdir()
    fa = tmp_path / "mveukv.fna"
    fa.write_text(f">QS|1|a\n{q_same}\n>QO|1|b\n{q_other}\n")
    (meta / "records.tsv").write_text("ordinal\tgenome\ttaxon\n0\tQS\tNucleocytoviricota\n1\tQO\tMaveriviricetes\n")
    (meta / "repeated_records.tsv").write_text("ordinal\tgenome\tgenome_uid\tfirst_genome\n")
    gv = tmp_path / "gv1.fna"
    gv.write_text(f">gvref\n{ref}\n")
    s = mv.ani_sample(meta, fa, gv, tmp_path / "ani", threads=1, every=1)
    assert s["samples"] == 2
    assert s["vs_gv1"]["shorter"]["meet"] == 1


# --- the conf templates -------------------------------------------------------------------------------------------

def test_templates_pass_stage_est_and_carry_the_cards_values():
    est = _load("stage_est", "stage_est.py")
    c2 = est.read_conf(TEMPLATES / "MV2.conf", "MV2", tree=REPO)
    ce = est.read_conf(TEMPLATES / "MVEUKV.conf", "MVEUKV", tree=REPO)
    assert (c2.depth, c2.crispr, c2.library, c2.after) == ("called", True, None, ("MGE4",))
    assert (ce.depth, ce.crispr, ce.library, ce.after) == ("annotated", True, ("pfam", "38.2", "gathering"), ("MV2",))
    assert str(c2.inputs) == "/tmp/if2453/s2/inputs/MV2" and str(ce.inputs) == "/tmp/if2453/s2/inputs/MVEUKV"
    assert c2.min_contig_len is None and ce.min_contig_len is None
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    allowed = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}
    allowed |= {s for s, _ in a26.R26_11_4_CHAIN}
    note = (REPO / NOTE).read_text(encoding="utf-8")
    for c, fasta_name in ((c2, "mv2.fna"), (ce, "mveukv.fna")):
        assert NOTE in c.notes
        assert {n for _, n in c.sums} == {"genomes.json", "objects.json", fasta_name}
        for h, _ in c.sums:
            assert h in note, f"{c.stage}: SUM {h} not in the note"
        assert set(c.expect) <= allowed
        for owner, (n, _) in c.expect.items():
            assert n > 0

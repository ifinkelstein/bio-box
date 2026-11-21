"""Ring R1: card C-47 part (2), the breakdown of the NR90 evaluation (`c47_nr90_pfam_eval.py breakdown`).

On C5 the evaluation tripped DEC-35 (6) on missed pairs (2026-11-21). The breakdown says where the missed and wrong
pairs come from and measures alternative rules for the PI. Every expected number below is worked out by hand from the
fixture (AGENTS §5.2): six proteins in two clusters, the in-memory test's scenario with p2 shortened to 92 residues
(8% shorter than its representative) and given identity 0.96, so that the length, identity and f sweeps separate.

    p1 rep (100, copied)  PF00001 [10,50) 40/38, PF00002 [60,90) 30/29
    p2 member of p1 (92, identity 0.96, f 1.0): own PF00001 41/39       -> inferred: PF00002 wrong
    p3 member of p1 (70, 30% shorter)                                    -> scanned for length
    p4 rep (200)          PF00010 [50,150) 60/59
    p5 member of p4 (195, identity 0.98, f 0.5)                          -> scanned (unmapped)
    p6 member of p4 (200, identity 0.995, f 0.9): own PF00010, PF00012 21.5/19.0 -> inferred: PF00012 missed
    GA (seq/dom): PF00001 25/25, PF00002 25/24, PF00003 20/20, PF00010 30/30, PF00011 25/25, PF00012 20/18
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-47")  # ring R1 by directory (tests/conftest.py)

EVAL_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "c47_nr90_pfam_eval.py"


def _load():
    spec = importlib.util.spec_from_file_location("c47_eval_bd", EVAL_SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c47_eval_bd"] = mod
    spec.loader.exec_module(mod)
    return mod


HMM = """HMMER3/f [3.4 | Aug 2023]
NAME  A
ACC   PF00001.26
GA    25.00 25.00;
//
HMMER3/f [3.4 | Aug 2023]
NAME  B
ACC   PF00002.31
GA    25.00 24.00;
//
NAME  C
ACC   PF00003.9
GA    20.00 20.00;
//
NAME  D
ACC   PF00010.12
GA    30.00 30.00;
//
NAME  E
ACC   PF00011.4
GA    25.00 25.00;
//
NAME  F
ACC   PF00012.21
GA    20.00 18.00;
//
"""


def test_parse_ga_keys_models_by_accession_without_version():
    m = _load()
    got = m.parse_ga(HMM.splitlines())
    assert got["PF00002"] == (25.0, 24.0)
    assert got["PF00012"] == (20.0, 18.0)
    assert len(got) == 6


def test_margin_bins():
    m = _load()
    assert m.margin_bin(None) == "no_score"
    assert m.margin_bin(-0.1) == "<0"
    assert m.margin_bin(0.0) == "[0,1)"
    assert m.margin_bin(1.5) == "[1,2)"
    assert m.margin_bin(5.0) == "[5,10)"
    assert m.margin_bin(30.0) == "[25,inf)"


def test_match_instance_takes_the_member_row_that_overlaps_most():
    m = _load()
    assert m.match_instance((100, 150), [(10, 60), (95, 150)]) == (95, 150)
    assert m.match_instance((100, 150), [(10, 60)]) is None
    assert m.match_instance((100, 150), []) is None


def _uid(name: str) -> bytes:
    return name.encode().ljust(16, b"\0")


def _fixture(tmp_path: Path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from biobox.l1.genecall import GENE_CALL_VERSION as GCV

    root = tmp_path / "root"
    U = pa.binary(16)

    def put(rel, table):
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, p)
        return rel

    def record(stage, step, rels):
        p = root / "runs" / f"stage={stage}" / f"{step}.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"outputs": [{"path": r} for r in rels]}))

    lens = {"p1": 100, "p2": 92, "p3": 70, "p4": 200, "p5": 195, "p6": 200}
    names = sorted(lens)
    ps = put(f"l0seq/protein_seq/gene_call_version={GCV}/stage=C5/protein_prefix=00/part-0.parquet", pa.table({
        "protein_uid": pa.array([_uid(n) for n in names], U), "length": pa.array([lens[n] for n in names], pa.uint32())}))
    contig = {"p1": "cA", "p2": "cA", "p3": "cA", "p4": "cB", "p5": "cB", "p6": "cB"}
    fp = put(f"l1/feature_protein/gene_call_version={GCV}/stage=C5/part-0.parquet", pa.table({
        "contig_uid": pa.array([_uid(contig[n]) for n in names], U), "protein_uid": pa.array([_uid(n) for n in names], U)}))
    record("C5", "genecall", [ps, fp])
    sha = {"GCA_A": b"\x01" * 32, "GCA_B": b"\x02" * 32}
    man = put("l0/l0_manifest/stage=C5/part-0.parquet", pa.table({
        "object_path": ["GCA_A.fna.gz", "GCA_B.fna.gz"], "object_sha256": pa.array(list(sha.values()), pa.binary(32))}))
    occ = put("l0seq/contig_occurrence/stage=C5/source=ncbi_genomes/part-0.parquet", pa.table({
        "l0_object_sha256": pa.array(list(sha.values()), pa.binary(32)), "contig_uid": pa.array([_uid("cA"), _uid("cB")], U)}))
    record("C5", "ingest", [man, occ])
    doms = [("p1", "PF00001.26", 10, 50, 40.0, 38.0), ("p1", "PF00002.30", 60, 90, 30.0, 29.0),
            ("p2", "PF00001.26", 10, 50, 41.0, 39.0), ("p3", "PF00001.26", 10, 50, 35.0, 34.0),
            ("p3", "PF00003.9", 55, 68, 22.0, 21.0), ("p4", "PF00010.12", 50, 150, 60.0, 59.0),
            ("p5", "PF00010.12", 50, 150, 61.0, 60.0), ("p5", "PF00011.4", 155, 190, 33.0, 30.0),
            ("p6", "PF00010.12", 50, 150, 62.0, 61.0), ("p6", "PF00012.21", 155, 195, 21.5, 19.0)]
    dom = put("l1/domains/library=pfam-38.2/stage=C5/part-0.parquet", pa.table({
        "protein_uid": pa.array([_uid(d[0]) for d in doms], U), "model_acc": [d[1] for d in doms],
        "env_start": pa.array([d[2] for d in doms], pa.uint32()), "env_end": pa.array([d[3] for d in doms], pa.uint32()),
        "seq_bitscore": pa.array([d[4] for d in doms], pa.float32()),
        "dom_bitscore": pa.array([d[5] for d in doms], pa.float32())}))
    record("C5", "domains-pfam-38.2", [dom])
    src = put(f"l0seq/protein_seq/gene_call_version={GCV}/stage=P1/protein_prefix=00/part-0.parquet",
              pa.table({"protein_uid": pa.array([_uid("p1")], U)}))
    plan = root / "l1/domain_scan/library=pfam-38.2/stage=C5/plan.json"
    plan.parent.mkdir(parents=True, exist_ok=True)
    plan.write_text(json.dumps({"copy": {"sources": [{"stage": "P1", "protein_files": [{"path": src}]}]}}))
    rep = {"p1": "p1", "p2": "p1", "p3": "p1", "p4": "p4", "p5": "p4", "p6": "p4"}
    cm = tmp_path / "cm.parquet"
    pq.write_table(pa.table({"member_uid": pa.array([_uid(n) for n in names], U), "level": ["nr90"] * 6,
                             "rep_uid": pa.array([_uid(rep[n]) for n in names], U)}), cm)
    # q* member, t* representative; the CIGAR's query is the representative (I: representative only, D: member only)
    rows = [("p1", 0, 100, 0, 100, b"100M", 1.0), ("p2", 0, 92, 0, 100, b"50M8I42M", 0.96),
            ("p3", 0, 70, 0, 70, b"70M", 0.97), ("p4", 0, 200, 0, 200, b"200M", 1.0),
            ("p5", 0, 100, 0, 100, b"100M", 0.98), ("p6", 60, 200, 60, 200, b"140M", 0.995)]
    aln = tmp_path / "aln.parquet"
    pq.write_table(pa.table({
        "member_uid": pa.array([_uid(r[0]) for r in rows], U), "level": ["nr90"] * 6,
        "rep_uid": pa.array([_uid(rep[r[0]]) for r in rows], U),
        **{k: pa.array([r[i] for r in rows], pa.uint32()) for i, k in enumerate(("qstart", "qend", "tstart", "tend"), 1)},
        "pident": pa.array([r[6] for r in rows], pa.float32()), "cigar": pa.array([r[5] for r in rows], pa.binary())}), aln)
    sel = tmp_path / "selection.tsv"
    sel.write_text("accession\tphylum\nGCA_A\tp__A\nGCA_B\tp__B\n")
    objs = tmp_path / "objects.json"
    objs.write_text("[]")
    hmm = tmp_path / "Pfam-A.hmm.gz"
    with gzip.open(hmm, "wt") as fh:
        fh.write(HMM)
    return root, cm, aln, sel, objs, hmm


def test_breakdown_on_a_constructed_fixture(tmp_path):
    m = _load()
    root, cm, aln, sel, objs, hmm = _fixture(tmp_path)
    out = tmp_path / "breakdown.json"
    assert m.main(["breakdown", "--root", str(root), "--stage", "C5", "--stage-cm", str(cm), "--stage-aln", str(aln),
                   "--selection", str(sel), "--objects", str(objs), "--hmm", str(hmm), "--out", str(out),
                   "--tmp", str(tmp_path / "t")]) == 0
    d = json.loads(out.read_text())

    base = d["base"]
    assert (base["true_pairs"], base["pairs_right"], base["pairs_false"], base["pairs_missed"]) == (10, 9, 1, 1)
    assert (base["scanned_residues"], base["total_residues"], base["scanned_residues_after_b47b"]) == (565, 857, 465)
    assert d["consistent_with_eval"] is True

    assert d["missed_by_cause"] == {"rep_has_no_row_for_the_model": 1, "rep_row_not_mapped": 0,
                                    "of_which_rep_has_no_domain_at_all": 0}
    assert d["wrong_by_cause"] == {"member_scan_lacks_the_model": 1}

    mg = d["margins"]
    assert mg["missed"]["seq"] == {"[1,2)": 1} and mg["missed"]["dom"] == {"[1,2)": 1}  # p6 PF00012: 21.5-20, 19-18
    assert mg["wrong"]["seq"] == {"[5,10)": 1} and mg["wrong"]["dom"] == {"[5,10)": 1}  # p1 PF00002: 30-25, 29-24
    assert mg["right_inferred"]["seq"] == {"[10,25)": 1, "[25,inf)": 1}  # p2 PF00001 41-25; p6 PF00010 62-30

    idb = d["by_identity"]
    assert idb["[0.95,0.97)"] == {"members": 1, "true_pairs": 1, "reported_pairs": 2, "pairs_false": 1, "pairs_missed": 0}
    assert idb["[0.99,1)"] == {"members": 1, "true_pairs": 2, "reported_pairs": 1, "pairs_false": 0, "pairs_missed": 1}
    lb = d["by_length"]
    assert lb["(0.05,0.10]"]["pairs_false"] == 1 and lb["0"]["pairs_missed"] == 1
    assert d["by_model"]["missed"] == [{"model": "PF00012", "pairs": 1, "share": 1.0}]
    assert d["by_model"]["wrong"] == [{"model": "PF00002", "pairs": 1, "share": 1.0}]

    sweep = {(r["length_cut"], r["identity_min"], r["f_cut"]): r for r in d["sweep"]}
    assert len(sweep) == 4 * 4 * 3

    def rates(key):
        r = sweep[key]
        return round(r["wrong_rate"], 4), round(r["missed_rate"], 4), r["scanned_residues"], r["scanned_residues_after_b47b"]

    assert rates((0.2, None, 0.8)) == (0.1, 0.1, 565, 465)
    assert rates((0.1, None, 0.8)) == (0.1, 0.1, 565, 465)  # p2 is 8% shorter
    assert rates((0.05, None, 0.8)) == (0.0, 0.1, 657, 557)  # p2 scanned
    assert rates((0.2, 0.95, 0.8)) == (0.1, 0.1, 565, 465)  # p2's identity 0.96 is not under 0.95
    assert rates((0.2, 0.97, 0.8)) == (0.0, 0.1, 657, 557)
    assert rates((0.2, None, 0.9)) == (0.1, 0.1, 565, 465)  # p6's f 0.9 is not under 0.9
    assert rates((0.2, None, 1.0)) == (round(1 / 11, 4), 0.0, 765, 665)  # p6 scanned: 11 reported, 1 wrong
    for level in ("0.01", "0.005"):
        c = d["cheapest"][level]
        assert (c["length_cut"], c["identity_min"], c["f_cut"]) == (0.2, 0.97, 1.0)
        assert (c["wrong_rate"], c["missed_rate"], c["scanned_residues"]) == (0.0, 0.0, 857)

    co = d["coordinates"]  # p2's two rep rows and p6's one; p2 has no PF00002 of its own
    assert (co["rep_rows"], co["matched"], co["no_member_instance"]) == (3, 2, 1)
    assert co["both_within_5"] == 0.5 and co["both_within_10"] == 1.0

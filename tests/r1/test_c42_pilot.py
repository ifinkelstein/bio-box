"""Card C-42: unit tests for the SPIRE ~1,000-sample unbinned-contig pilot scripts.

Covers:
1. `c42_select.py`: microntology habitat assignment, square-root stratified allocation (Hamilton
   largest-remainder with floor >= 1 and cap <= N_h), deterministic seeded sample selection
   independent of input row order, and per-habitat assembly base summaries.
2. `c42_contigs.py`: streaming FASTA / `.fa.gz` contig length shares (>= 1.5, 5, 10 kb), binned
   base share, and exact `contig_uid` (strand-canonical BLAKE3-128) repeat counting across samples,
   checked against an independent brute-force oracle.
3. `c42_inputs.py`: direct generation of `genomes.json` and `objects.json` (one prokaryote
   `contig`-quality `GenomeInput` per distinct `contig_uid`, optional `--min-contig-len`), byte
   determinism across input sample orders, refusal on bad inputs, and round-trip compatibility with
   `biobox.l0seq.inputs.write_genome_input` and `biobox.l0seq.ingest`.
4. `c42_novelty.py`: novel NR50 family counting against a reference release's protein UIDs and
   exact sample-by-sample accumulation curves per habitat, checked against a brute-force set oracle.
5. `c42_tniq.py`: strand-oriented +-10-gene neighborhood window extraction from `contig_tapes`,
   contig-edge truncation flags, nearby `l1_crispr` array matching, and CAST locus classification.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest

from biobox import ids
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import read_genome_input, write_genome_input

pytestmark = pytest.mark.covers("C-42", "MR-01")

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod"


def _load(name: str):
    path = SCRIPTS / f"{name}.py"
    if not path.is_file():
        raise FileNotFoundError(f"missing script {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# --- 1. c42_select.py -------------------------------------------------------------------------


def test_c42_allocate_sqrt_preserves_rare_strata_and_hits_exact_target():
    m = _load("c42_select")
    counts = {
        "host-associated:animal host": 62392,
        "aquatic:marine": 5340,
        "terrestrial:soil": 1874,
        "anthropogenic:wastewater": 591,
        "air": 17,
        "terrestrial:tundra": 1,
    }
    alloc = m.allocate_sqrt(counts, target=100)
    assert sum(alloc.values()) == 100
    assert set(alloc) == set(counts)
    for k, n in counts.items():
        assert 1 <= alloc[k] <= n
    # Rare strata are preserved: tundra (N=1) gets 1, air (N=17) gets >= 1, and animal host
    # gets far less than proportional-to-N (62392/70215 = 88.9% -> sqrt gives ~58).
    assert alloc["terrestrial:tundra"] == 1
    assert alloc["air"] >= 1
    assert alloc["host-associated:animal host"] < 70
    assert alloc["host-associated:animal host"] > alloc["aquatic:marine"] > alloc["terrestrial:soil"]

    # Order independence: reversed dict insertion order produces identical allocation.
    rev = dict(reversed(list(counts.items())))
    assert m.allocate_sqrt(rev, target=100) == alloc

    # Target exceeding total available samples caps every stratum at N_h.
    capped = m.allocate_sqrt({"a": 3, "b": 1}, target=10)
    assert capped == {"a": 3, "b": 1}


def test_c42_select_pilot_is_deterministic_and_summarizes_habitats():
    m = _load("c42_select")
    study_rows = [
        {"study": "S_gut", "sample_id": f"SAM_G{i:02d}", "microntology": "['host-associated', 'host-associated:animal host']"}
        for i in range(16)
    ] + [
        {"study": "S_mar", "sample_id": f"SAM_M{i:02d}", "microntology": "['aquatic', 'aquatic:marine']"}
        for i in range(4)
    ] + [
        {"study": "S_air", "sample_id": "SAM_A00", "microntology": "['air']"}
    ]
    api_rows = [
        {"sample_name": r["sample_id"], "study": r["study"], "assembly_bp": 10_000_000 * (idx + 1),
         "assembly_n50": 2500, "binned_fraction": 0.25}
        for idx, r in enumerate(study_rows)
    ]
    mag_bp = {"SAM_G00": (2, 3_000_000), "SAM_M00": (1, 4_000_000)}

    merged_a = m.merge_sample_tables(study_rows, api_rows, mag_bp)
    merged_b = m.merge_sample_tables(list(reversed(study_rows)), list(reversed(api_rows)), mag_bp)
    assert merged_a == merged_b

    hab_summary = {r["habitat"]: r for r in m.summarize_by_habitat(merged_a, key="habitat")}
    assert hab_summary["host-associated"]["samples"] == 16
    assert hab_summary["aquatic"]["samples"] == 4
    assert hab_summary["air"]["samples"] == 1
    assert hab_summary["host-associated"]["assembly_bp"] == sum(10_000_000 * (i + 1) for i in range(16))
    assert hab_summary["host-associated"]["mag_bp"] == 3_000_000

    pick_1 = m.select_pilot(merged_a, target=7, seed=20261118, stratum_key="habitat2")
    pick_2 = m.select_pilot(list(reversed(merged_a)), target=7, seed=20261118, stratum_key="habitat2")
    assert pick_1 == pick_2
    assert len(pick_1) == 7
    assert any(r["sample_id"] == "SAM_A00" for r in pick_1)


# --- 2. c42_contigs.py ------------------------------------------------------------------------


def _revcomp(seq: bytes) -> bytes:
    tbl = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
    return seq.translate(tbl)[::-1]


def test_c42_contigs_length_shares_bins_and_exact_repeats(tmp_path: Path):
    m = _load("c42_contigs")
    c_short = b"AAACGT" * 133 + b"AA"  # 800 bp (< 1.5 kb, non-palindromic)
    c_mid = b"TTTGCA" * 333 + b"TT"    # 2,000 bp (>= 1.5 kb, < 5 kb, non-palindromic)
    c_5k = b"AAACCG" * 1000            # 6,000 bp (>= 5 kb, < 10 kb)
    c_10k = b"ATGCCG" * 2000           # 12,000 bp (>= 10 kb)

    # Sample 1 (aquatic): c_short, c_mid, c_10k (total = 14,800 bp)
    s1 = tmp_path / "S1.fa.gz"
    with gzip.open(s1, "wb") as fh:
        fh.write(b">k1\n" + c_short + b"\n>k2\n" + c_mid + b"\n>k3\n" + c_10k + b"\n")

    # Sample 2 (terrestrial): lowercase c_short (exact contig_uid match!), reverse-complement of c_mid
    # (strand-canonical match!), plus c_5k (total = 8,800 bp)
    s2 = tmp_path / "S2.fa.gz"
    with gzip.open(s2, "wb") as fh:
        fh.write(b">x0\n" + c_short.lower() + b"\n>x1\n" + _revcomp(c_mid) + b"\n>x2\n" + c_5k + b"\n")

    meta = [
        {"sample_id": "S1", "habitat": "aquatic", "habitat2": "aquatic:marine", "mag_bp": 6000, "binned_fraction": 0.5},
        {"sample_id": "S2", "habitat": "terrestrial", "habitat2": "terrestrial:soil", "mag_bp": 2000, "binned_fraction": 0.25},
    ]
    res = m.analyze_pilot_contigs(meta, {"S1": s1, "S2": s2}, tmp_path / "work", keep_chunks=True)
    res_reuse = m.analyze_pilot_contigs(meta, {}, tmp_path / "work", reuse_chunks=True)
    assert res == res_reuse

    # Brute-force oracle:
    # S1: total 14,800; >=1.5k: 14,000; >=5k: 12,000; >=10k: 12,000; mag_bp = 6,000
    # S2: total 8,800; >=1.5k: 8,000; >=5k: 6,000; >=10k: 0; mag_bp = 2,000
    by_hab = {r["habitat"]: r for r in res["by_habitat"]}
    assert by_hab["aquatic"]["total_bp"] == 14_800
    assert by_hab["aquatic"]["bp_ge_1500"] == 14_000
    assert by_hab["aquatic"]["bp_ge_5000"] == 12_000
    assert by_hab["aquatic"]["bp_ge_10000"] == 12_000
    assert by_hab["aquatic"]["mag_bp"] == 6_000
    assert abs(by_hab["aquatic"]["share_ge_1500"] - round(14_000 / 14_800, 6)) < 1e-9

    assert by_hab["terrestrial"]["total_bp"] == 8_800
    assert by_hab["terrestrial"]["bp_ge_1500"] == 8_000
    assert by_hab["terrestrial"]["bp_ge_5000"] == 6_000
    assert by_hab["terrestrial"]["bp_ge_10000"] == 0

    # Exact contig_uid repeats: 5 distinct contig_uids (c_short shared across S1 & S2, revcomp(c_mid) distinct),
    # while strand-canonical repeats has 4 distinct canon_uids and 2 shared across S1 & S2 (c_short + c_mid)!
    rep = res["repeats"]
    assert rep["total_records"] == 6
    assert rep["distinct_contig_uids"] == 5
    assert rep["contigs_in_multiple_samples"] == 1
    assert rep["contigs_across_habitats"] == 1
    assert rep["repeated_bp_between_samples"] == 800
    assert rep["distinct_strand_canonical_uids"] == 4
    assert rep["strand_canonical_contigs_in_multiple_samples"] == 2
    assert rep["strand_canonical_repeated_bp_between_samples"] == 2_800


# --- 3. c42_inputs.py -------------------------------------------------------------------------


def test_c42_inputs_deterministic_and_roundtrips_with_biobox_ingest(tmp_path: Path):
    m = _load("c42_inputs")
    seq_a = b"ATGAAATTTCCCGGG" * 120   # 1,800 bp
    seq_b = b"ATGCCCTTTAAAGGG" * 40    # 600 bp
    seq_c = b"ATGGGGCCCTTTAAA" * 150   # 2,250 bp

    f1 = tmp_path / "raw" / "SAM1.fa.gz"
    f2 = tmp_path / "raw" / "SAM2.fa.gz"
    f1.parent.mkdir(parents=True)
    with gzip.open(f1, "wb") as fh:
        fh.write(b">c1\n" + seq_a + b"\n>c2\n" + seq_b + b"\n")
    with gzip.open(f2, "wb") as fh:
        # seq_a repeated in SAM2 (lowercase -> identical contig_uid), plus seq_c
        fh.write(b">d1\n" + seq_a.lower() + b"\n>d2\n" + seq_c + b"\n")

    samples = [
        {"sample_id": "SAM1", "retrieved_at": "2026-11-18T20:00:00Z"},
        {"sample_id": "SAM2", "retrieved_at": "2026-11-18T20:01:00Z"},
    ]
    fastas = {"SAM1": f1, "SAM2": f2}

    out1 = tmp_path / "in_order1"
    out2 = tmp_path / "in_order2"
    rep1 = m.build_stage_inputs(samples, fastas, out1, min_contig_len=0)
    rep2 = m.build_stage_inputs(list(reversed(samples)), fastas, out2, min_contig_len=0)

    assert (out1 / "genomes.json").read_bytes() == (out2 / "genomes.json").read_bytes()
    assert (out1 / "objects.json").read_bytes() == (out2 / "objects.json").read_bytes()
    assert rep1["distinct_contigs"] == 3
    assert rep1["total_records"] == 4

    # Round-trip through biobox write_genome_input and ingest!
    root = tmp_path / "root"
    genomes_raw = json.loads((out1 / "genomes.json").read_text())
    from biobox.l0seq.inputs import GenomeInput

    gis = [
        GenomeInput(
            contig_uids=tuple(bytes.fromhex(u) for u in g["contig_uids"]),
            genetic_code=g["genetic_code"],
            kind=g["kind"],
            circular=tuple(bytes.fromhex(u) for u in g["circular"]),
            assembly_quality=g["assembly_quality"],
        )
        for g in genomes_raw
    ]
    gi_res = write_genome_input(root, "SP1", gis)
    assert gi_res == {"genomes": 3, "contigs": 3}
    assert len(read_genome_input(root, "SP1")) == 3

    objs_raw = json.loads((out1 / "objects.json").read_text())
    from datetime import datetime

    objs = [
        SourceObject(
            path=Path(o["path"]),
            source_name=o["source_name"],
            source_release=o["source_release"],
            license=o["license"],
            redistributable=o["redistributable"],
            retrieved_at=datetime.fromisoformat(o["retrieved_at"]),
        )
        for o in objs_raw
    ]
    ing_res = ingest(objs, root, "SP1")
    assert ing_res["objects"] == 2
    assert ing_res["records"] == 4
    import pyarrow.parquet as pq
    ci_tbl = pq.read_table(root / "l0seq" / "contig_index" / "stage=SP1")
    assert sum(1 for v in ci_tbl["is_default_variant"].to_pylist() if v) == 3

    # Test --min-contig-len 1500 drops seq_b (600 bp) consistently from both FASTAs and genomes.json
    out_filt = tmp_path / "in_filt"
    rep_f = m.build_stage_inputs(samples, fastas, out_filt, min_contig_len=1500)
    assert rep_f["distinct_contigs"] == 2
    assert rep_f["dropped_records"] == 1
    assert rep_f["dropped_bp"] == 600
    assert {ids.contig_uid(seq_a).hex(), ids.contig_uid(seq_c).hex()} == {
        g["contig_uids"][0] for g in json.loads((out_filt / "genomes.json").read_text())
    }


# --- 4. c42_novelty.py ------------------------------------------------------------------------


def test_c42_novelty_matches_bruteforce_set_oracle(tmp_path: Path):
    m = _load("c42_novelty")
    p_ref1 = ids.protein_uid("MKKKLLLAAA")
    p_ref2 = ids.protein_uid("MGGGSSSTTT")
    p_sh = ids.protein_uid("MKKKLLLAAB")    # clusters with p_ref1 -> NOT novel
    p_n1 = ids.protein_uid("MWWWYYYFFF")    # novel cluster C1 (singleton in S1, seen again via p_n1b in S2)
    p_n1b = ids.protein_uid("MWWWYYYFFE")   # same novel cluster C1 (>=2 members)
    p_n2 = ids.protein_uid("MCCCQQQNNN")    # novel singleton cluster C2 in S1
    p_n3 = ids.protein_uid("MHHHIIIVVV")    # novel singleton cluster C3 in S3

    clu_tsv = tmp_path / "clu50.tsv"
    rows = [
        (p_ref1.hex(), p_ref1.hex()),
        (p_ref1.hex(), p_sh.hex()),
        (p_ref2.hex(), p_ref2.hex()),
        (p_n1.hex(), p_n1.hex()),
        (p_n1.hex(), p_n1b.hex()),
        (p_n2.hex(), p_n2.hex()),
        (p_n3.hex(), p_n3.hex()),
    ]
    clu_tsv.write_text("".join(f"{r}\t{mem}\n" for r, mem in rows))

    samples = [
        {"sample_id": "S1", "habitat": "aquatic", "bp": 5_000_000, "proteins": [p_sh, p_n1, p_n2]},
        {"sample_id": "S2", "habitat": "aquatic", "bp": 5_000_000, "proteins": [p_n1b]},
        {"sample_id": "S3", "habitat": "terrestrial", "bp": 2_000_000, "proteins": [p_n1, p_n3]},
    ]
    out = m.compute_novelty(clu_tsv, {p_ref1, p_ref2}, samples, tmp_path / "nov_work")

    # Brute-force truth:
    # Total clusters = 5: {p_ref1, p_sh} (ref+pilot), {p_ref2} (ref-only), {p_n1, p_n1b} (novel, size 2),
    # {p_n2} (novel, size 1), {p_n3} (novel, size 1).
    assert out["summary"]["total_clusters"] == 5
    assert out["summary"]["ref_clusters"] == 2
    assert out["summary"]["shared_clusters"] == 1
    assert out["summary"]["novel_clusters"] == 3
    assert out["summary"]["novel_nonsingleton_clusters"] == 1

    by_hab = {r["habitat"]: r for r in out["by_habitat"]}
    # aquatic sees {p_n1, p_n1b} and {p_n2} -> 2 novel clusters over 0.01 Gbp = 200 / Gbp
    assert by_hab["aquatic"]["novel_clusters_in_habitat"] == 2
    assert abs(by_hab["aquatic"]["novel_per_gbp"] - 200.0) < 1e-6
    # terrestrial sees {p_n1, p_n1b} and {p_n3} -> 2 in habitat (1 marginal after aquatic) over 0.002 Gbp = 1000 / Gbp
    assert by_hab["terrestrial"]["novel_clusters_in_habitat"] == 2
    assert abs(by_hab["terrestrial"]["novel_per_gbp"] - 1000.0) < 1e-6

    # Global curve across S1 -> S2 -> S3: +2 at S1, +0 at S2 (p_n1b joins C1 already found in S1!), +1 at S3 (p_n3)
    curve = out["global_curve"]
    assert [(r["sample_id"], r["new_clusters"], r["cum_new_clusters"]) for r in curve] == [
        ("S1", 2, 2),
        ("S2", 0, 2),
        ("S3", 1, 3),
    ]


def test_c42_novelty_linclust_split_align_and_direct_createdb(tmp_path: Path):
    import shutil
    import pyarrow as pa
    import pyarrow.parquet as pq
    from biobox.l1.genecall import GENE_CALL_VERSION

    if shutil.which("mmseqs") is None:
        pytest.skip("mmseqs not on PATH")

    m = _load("c42_novelty")

    def write_stage_proteins(root: Path, stage: str, seqs: list[str]) -> None:
        by_px: dict[str, list[tuple[bytes, str]]] = {}
        for s in seqs:
            u = ids.protein_uid(s)
            by_px.setdefault(u[:1].hex(), []).append((u, s))
        for px, items in sorted(by_px.items()):
            items.sort(key=lambda x: x[0])
            uniq = []
            for u, s in items:
                if not uniq or uniq[-1][0] != u:
                    uniq.append((u, s))
            out_dir = root / "l0seq" / "protein_seq" / f"gene_call_version={GENE_CALL_VERSION}" / f"stage={stage}" / f"protein_prefix={px}"
            out_dir.mkdir(parents=True, exist_ok=True)
            pq.write_table(
                pa.table({
                    "protein_uid": pa.array([u for u, _ in uniq], type=pa.binary(16)),
                    "aa_len": pa.array([len(s) for _, s in uniq], type=pa.int32()),
                    "seq": pa.array([s for _, s in uniq], type=pa.string()),
                }),
                out_dir / "part-00000.parquet",
            )

    base_a = "MKKLLVVAAIIGGSSTTDDNNQQEERRHHYYWWFFCCMM" * 3
    base_b = "MPPQQNNSSTTVVLLAAIIGGFFEEDDKKRRHHYYWWCC" * 3
    mut_a1 = base_a[:80] + "A" * 10 + base_a[90:]
    mut_a2 = base_a[:60] + "V" * 12 + base_a[72:]
    novel_c = "MWWYYFFHHCCMMNNQQEERRKKDDSSGGTTIIVVLLAA" * 3
    novel_c2 = novel_c[:70] + "L" * 10 + novel_c[80:]
    novel_d = "MCCAALLVVIIGGSSTTDDNNQQEERRKKHHYYWWFFMM" * 3
    novel_e = "MGGSSTTAALLVVIIPPQQNNDDEERRKKHHYYWWFFCC" * 3

    ref_root = tmp_path / "ref_root"
    pilot_root = tmp_path / "pilot_root"
    write_stage_proteins(ref_root, "S1", [base_a, base_b])
    write_stage_proteins(pilot_root, "SP1", [mut_a1, mut_a2, novel_c, novel_c2, novel_d, novel_e])

    # Oracle: single-part linclust-50 via combined.faa
    w1 = tmp_path / "w_one"
    assert m.main([
        "export-faa", "--ref-root", str(ref_root), "--ref-stages", "S1",
        "--pilot-root", str(pilot_root), "--pilot-stages", "SP1", "--work-dir", str(w1),
    ]) == 0
    assert m.main([
        "linclust-50", "--work-dir", str(w1), "--threads", "1", "--align-parts", "1", "--split-mem", "1G",
    ]) == 0

    # Under test: direct stdin createdb (no combined.faa written) + 3-part split align + cleanup-db
    w3 = tmp_path / "w_three"
    assert m.main([
        "export-faa", "--ref-root", str(ref_root), "--ref-stages", "S1",
        "--pilot-root", str(pilot_root), "--pilot-stages", "SP1", "--work-dir", str(w3),
        "--direct-createdb",
    ]) == 0
    assert not (w3 / "combined.faa").exists()
    assert (w3 / "db.dbtype").is_file()
    assert m.main([
        "linclust-50", "--work-dir", str(w3), "--threads", "1", "--align-parts", "3", "--split-mem", "1G",
        "--cleanup-db",
    ]) == 0
    assert not (w3 / "db.dbtype").exists()
    assert not (w3 / "lin50").exists()

    def parse_clusters(tsv: Path) -> set[frozenset[str]]:
        by_rep: dict[str, set[str]] = {}
        for ln in tsv.read_text().splitlines():
            if ln.strip():
                r, mem = ln.split("\t")
                by_rep.setdefault(r, set()).add(mem)
        return {frozenset(v) for v in by_rep.values()}

    assert parse_clusters(w3 / "clu50.tsv") == parse_clusters(w1 / "clu50.tsv")


def test_c42_per_stage_roots_exclude_earlier_owned_contigs_and_summarize(tmp_path: Path):
    m = _load("c42_inputs")
    seq_shared = b"ATGAAATTTCCCGGG" * 120   # 1,800 bp (shared across stage SPA1 and SPQ1)
    seq_only_q = b"ATGGGGCCCTTTAAA" * 150   # 2,250 bp (unique to SPQ1)
    cu_shared = ids.contig_uid(seq_shared)
    cu_only_q = ids.contig_uid(seq_only_q)

    fa1 = tmp_path / "raw" / "S_A1.fa.gz"
    fa2 = tmp_path / "raw" / "S_Q1.fa.gz"
    fa1.parent.mkdir(parents=True)
    with gzip.open(fa1, "wb") as fh:
        fh.write(b">a1\n" + seq_shared + b"\n")
    with gzip.open(fa2, "wb") as fh:
        fh.write(b">q_dup\n" + seq_shared + b"\n>q_new\n" + seq_only_q + b"\n")

    in_a = tmp_path / "inputs" / "SPA1"
    rep_a = m.build_stage_inputs(
        [{"sample_id": "S_A1", "retrieved_at": "2026-11-18T20:00:00Z"}],
        {"S_A1": fa1},
        in_a,
        min_contig_len=1500,
    )
    assert rep_a["distinct_contigs"] == 1

    # Ingest SPA1 into its own isolated per-stage root `roots/SPA1`
    roots = tmp_path / "roots"
    from datetime import datetime
    objs_a = [
        SourceObject(
            path=Path(o["path"]),
            source_name=o["source_name"],
            source_release=o["source_release"],
            license=o["license"],
            redistributable=o["redistributable"],
            retrieved_at=datetime.fromisoformat(o["retrieved_at"]),
        )
        for o in json.loads((in_a / "objects.json").read_text())
    ]
    ingest(objs_a, roots / "SPA1", "SPA1")

    # Export owned UIDs from SPA1 and build SPQ1 inputs excluding them
    owned_bin = tmp_path / "owned.bin"
    assert m.write_owned_uids_bin(roots, ["SPA1"], owned_bin) == 1

    in_q = tmp_path / "inputs" / "SPQ1"
    rep_q = m.build_stage_inputs(
        [{"sample_id": "S_Q1", "retrieved_at": "2026-11-18T20:01:00Z"}],
        {"S_Q1": fa2},
        in_q,
        min_contig_len=1500,
        exclude_uids_bin=owned_bin,
    )
    assert rep_q["distinct_contigs"] == 1
    assert rep_q["held_by_earlier_stage_records"] == 1
    assert rep_q["held_by_earlier_stage_bp"] == 1800
    assert [g["contig_uids"][0] for g in json.loads((in_q / "genomes.json").read_text())] == [cu_only_q.hex()]
    cross_lines = (in_q / "cross_stage_occurrences.tsv").read_text().splitlines()
    assert len(cross_lines) == 2
    assert f"S_Q1\t{cu_shared.hex()}\tq_dup\t1800" == cross_lines[1]
    assert m.stage_root(roots, "SPA1") == roots / "SPA1"




# --- 5. c42_tniq.py ---------------------------------------------------------------------------


def test_c42_tniq_strand_oriented_windows_and_cast_classification():
    m = _load("c42_tniq")
    cu = ids.contig_uid(b"ACGT" * 5000)
    # 5 genes on minus strand: in genomic order (0..4): cas6f, cas7f, cas5f, tniQ (anchor at idx=3, strand=-1), tnsC
    # Along anchor's strand (-1), +1 downstream of tniQ is idx=2 (cas5f), +2 is idx=1 (cas7f), +3 is idx=0 (cas6f),
    # and -1 upstream of tniQ is idx=4 (tnsC)!
    p_cas6f = ids.protein_uid("MCAS6F")
    p_cas7f = ids.protein_uid("MCAS7F")
    p_cas5f = ids.protein_uid("MCAS5F")
    p_tniq = ids.protein_uid("MTNIQ")
    p_tnsc = ids.protein_uid("MTNSC")
    p_other = ids.protein_uid("MOTHER")
    p_uids = [p_cas6f, p_cas7f, p_cas5f, p_tniq, p_tnsc]
    starts = [500, 1500, 2600, 4800, 6200]
    ends = [1200, 2500, 4700, 6000, 7500]
    strands = [-1, -1, -1, -1, -1]
    f_uids = [ids.feature_uid(cu, s, e, st, ids.KIND["cds"]) for s, e, st in zip(starts, ends, strands)]

    tapes_rows = [{
        "contig_uid": cu,
        "contig_len": 8000,
        "feature_uids": f_uids,
        "protein_uids": p_uids,
        "starts": starts,
        "ends": ends,
        "strands": strands,
        "kinds": [ids.KIND["cds"]] * 5,
    }]
    crispr_rows = [{"contig_uid": cu, "start": 100, "end": 350, "n_repeats": 4, "spacer_count": 3}]

    windows, req_uids = m.extract_hit_windows(tapes_rows, {p_tniq}, crispr_rows=crispr_rows, window=10)
    assert len(windows) == 1
    assert req_uids == sorted(set(p_uids))
    assert p_other not in req_uids

    w = windows[0]
    assert w["truncated_5p"] is True   # only 1 upstream gene (tnsC) < 10
    assert w["truncated_3p"] is True   # only 3 downstream genes < 10
    assert len(w["crispr_arrays"]) == 1
    by_off = {s["offset_rel"]: s["protein_uid"] for s in w["slots"]}
    assert by_off == {-1: p_tnsc, 0: p_tniq, 1: p_cas5f, 2: p_cas7f, 3: p_cas6f}

    pfam = {
        p_tniq: [("PF06527.17", "TniQ", 65.0)],
        p_cas5f: [("PF09614.15", "Cas_Csy2", 80.0)],
        p_cas7f: [("PF09615.16", "Cas_Csy3", 120.0)],
        p_cas6f: [("PF09618.15", "Cas_Csy4", 55.0)],
        p_tnsc: [("PF05621.16", "TniB", 70.0)],
    }
    cls = m.classify_tniq_locus(w["slots"], pfam, w["crispr_arrays"])
    assert cls["locus_class"] == "I-F3_CAST"
    assert cls["has_cas1_cas2_cas3"] is False
    assert cls["has_crispr_array"] is True
    assert cls["cas5f_at_plus1"] is True

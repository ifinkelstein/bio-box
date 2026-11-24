"""Card C-62 (DEC-71): deterministic MetaVR vOTU representative + >= 30 kb jumbo/large UViG selection,
DEC-65 DTR trimming, <= 3.0M contig / <= 40 Gbp stage-cap MV1/MV2 splitting, and stage MV1 inputs.

Oracles are written here separately from the script: contig_uid via blake3 directly, representative
selection via a brute-force Python sort over shuffled rows.
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import random
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-62")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts" / "stress" / "build_stage_pod"


def _load():
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location("c62_metavr_inputs", SCRIPTS / "c62_metavr_inputs.py")
    if spec is None or spec.loader is None:
        raise FileNotFoundError(SCRIPTS / "c62_metavr_inputs.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c62_metavr_inputs"] = mod
    spec.loader.exec_module(mod)
    return mod


def cuid(seq: str) -> str:
    return blake3(seq.upper().encode()).digest(length=16).hex()


def seq_of(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_row(
    uvig: str,
    votu: str,
    *,
    quality: str = "High-quality",
    completeness: float | None = 100.0,
    contamination: float | None = 0.0,
    length: int = 40_000,
    topology: str = "Linear",
    is_concatemer: bool = False,
    viral_confidence: str = "High",
    genetic_code: str = "11",
    genome_type: str = "dsDNA",
    source: str = "Metagenome",
    taxon_oid: int = 3300000001,
    ecosystem: str = "Environmental",
    host_taxonomy: str = "d__Bacteria;p__Pseudomonadota",
) -> dict:
    return {
        "uvig": uvig,
        "taxon_oid": taxon_oid,
        "votu": votu,
        "quality": quality,
        "completeness": completeness,
        "checkv_contamination": contamination,
        "length": length,
        "uvig_topology": topology,
        "is_concatemer": is_concatemer,
        "viral_confidence": viral_confidence,
        "genetic_code": genetic_code,
        "genome_type": genome_type,
        "source": source,
        "Ecosystem": ecosystem,
        "host_taxonomy": host_taxonomy,
    }


QUALITY_RANK = {
    "Complete": 0,
    "High-quality": 1,
    "Medium-quality": 2,
    "Low-quality": 3,
    "Not-determined": 4,
}


def brute_select(rows: list[dict], *, max_contigs: int = 3_000_000, max_bases: int = 40_000_000_000):
    """Independent brute-force oracle for C-62 selection and MV1/MV2 splitting."""
    non_concat = [r for r in rows if not r["is_concatemer"]]
    by_votu: dict[str, list[dict]] = {}
    for r in non_concat:
        by_votu.setdefault(r["votu"], []).append(r)

    reps: dict[str, dict] = {}
    for v, grp in by_votu.items():
        grp.sort(
            key=lambda r: (
                QUALITY_RANK.get(r["quality"], 5),
                -min(r["completeness"] if r["completeness"] is not None else -1.0, 100.0),
                r["checkv_contamination"] if r["checkv_contamination"] is not None else 999.0,
                -r["length"],
                r["uvig"],
            )
        )
        reps[v] = grp[0]
    rep_uvigs = {r["uvig"] for r in reps.values()}

    candidates = []
    deferred = []
    for r in rows:
        if r["is_concatemer"]:
            deferred.append((r["uvig"], "concatemer"))
            continue
        is_rep = r["uvig"] in rep_uvigs
        contam = r["checkv_contamination"] if r["checkv_contamination"] is not None else 0.0
        comp = min(r["completeness"] if r["completeness"] is not None else -1.0, 100.0)
        qrank = QUALITY_RANK.get(r["quality"], 5)
        if r["length"] >= 200_000:
            key = (1, 0 if is_rep else 1, qrank, -r["length"], r["uvig"])
            candidates.append((key, "jumbo_ge200k", r))
        elif not is_rep:
            deferred.append((r["uvig"], "not_representative"))
        elif contam > 10.0:
            deferred.append((r["uvig"], "high_contamination"))
        elif r["quality"] in ("Complete", "High-quality"):
            key = (2, qrank, -comp, contam, -r["length"], r["uvig"])
            candidates.append((key, "rep_complete_hq", r))
        elif r["quality"] == "Medium-quality" and r["length"] >= 30_000:
            key = (3, -r["length"], -comp, contam, r["uvig"])
            candidates.append((key, "rep_mq_ge30k", r))
        elif r["quality"] == "Medium-quality":
            key = (4, -comp, -r["length"], contam, r["uvig"])
            candidates.append((key, "rep_mq_lt30k", r))
        elif r["length"] >= 30_000:
            conf_rank = 0 if r["viral_confidence"] == "High" else 1
            key = (5, conf_rank, -r["length"], contam, r["uvig"])
            candidates.append((key, "rep_large_ge30k_lq_nd", r))
        else:
            deferred.append((r["uvig"], "lq_nd_lt30k"))

    candidates.sort(key=lambda x: x[0])
    selected = []
    tot_c = tot_b = 0
    capped = False
    for _, tier, r in candidates:
        if not capped and (tot_c + 1 <= max_contigs) and (tot_b + r["length"] <= max_bases):
            selected.append((r["uvig"], tier))
            tot_c += 1
            tot_b += r["length"]
        else:
            capped = True
            deferred.append((r["uvig"], "mv2_spillover"))
    return selected, deferred


def _write_parquet(path: Path, rows: list[dict]) -> None:
    schema = pa.schema([
        ("uvig", pa.string()),
        ("taxon_oid", pa.int64()),
        ("votu", pa.string()),
        ("quality", pa.string()),
        ("completeness", pa.float64()),
        ("checkv_contamination", pa.float64()),
        ("length", pa.int64()),
        ("uvig_topology", pa.string()),
        ("is_concatemer", pa.bool_()),
        ("viral_confidence", pa.string()),
        ("genetic_code", pa.string()),
        ("genome_type", pa.string()),
        ("source", pa.string()),
        ("Ecosystem", pa.string()),
        ("host_taxonomy", pa.string()),
    ])
    table = pa.Table.from_pylist(rows, schema=schema)
    pq.write_table(table, path)


RULE_ROWS = [
    # Concatemer jumbo: MUST be excluded despite Complete + >=200kb
    make_row("U_CAT", "vOTU_1", quality="Complete", length=250_000, is_concatemer=True),
    # vOTU_1: Complete beats HQ; 130% completeness capped at 100% so lower contam wins
    make_row("U_1A", "vOTU_1", quality="Complete", completeness=130.0, contamination=2.0, length=60_000),
    make_row("U_1B", "vOTU_1", quality="Complete", completeness=100.0, contamination=0.0, length=45_000),
    make_row("U_1C", "vOTU_1", quality="High-quality", completeness=100.0, contamination=0.0, length=90_000),
    # vOTU_2: Jumbo rep + jumbo non-rep (BOTH selected in Priority 1) + non-rep <200kb (deferred)
    make_row("U_2A", "vOTU_2", quality="High-quality", length=240_000),
    make_row("U_2B", "vOTU_2", quality="Medium-quality", length=210_000),
    make_row("U_2C", "vOTU_2", quality="High-quality", length=80_000),
    # vOTU_3: MQ >=30kb low contam (Priority 3)
    make_row("U_3A", "vOTU_3", quality="Medium-quality", completeness=70.0, contamination=3.0, length=35_000),
    # vOTU_4: MQ <30kb low contam (Priority 4)
    make_row("U_4A", "vOTU_4", quality="Medium-quality", completeness=65.0, contamination=0.0, length=18_000),
    # vOTU_5: MQ high contam >10% (deferred: high_contamination)
    make_row("U_5A", "vOTU_5", quality="Medium-quality", completeness=80.0, contamination=15.0, length=40_000),
    # vOTU_6: LQ >=30kb High confidence (Priority 5)
    make_row("U_6A", "vOTU_6", quality="Low-quality", completeness=40.0, contamination=0.0, length=55_000, viral_confidence="High"),
    # vOTU_7: LQ >=30kb Low confidence (Priority 5, after High confidence)
    make_row("U_7A", "vOTU_7", quality="Low-quality", completeness=40.0, contamination=0.0, length=75_000, viral_confidence="Low"),
    # vOTU_8: LQ <30kb (deferred: lq_nd_lt30k)
    make_row("U_8A", "vOTU_8", quality="Low-quality", completeness=25.0, contamination=0.0, length=12_000),
]


def test_select_mv1_deterministic_and_excludes_concatemers(tmp_path):
    mv = _load()
    want_sel, want_def = brute_select(RULE_ROWS)
    assert [u for u, _ in want_sel] == ["U_2A", "U_2B", "U_1B", "U_3A", "U_4A", "U_6A", "U_7A"]
    assert dict(want_def)["U_CAT"] == "concatemer"
    assert dict(want_def)["U_5A"] == "high_contamination"
    assert dict(want_def)["U_8A"] == "lq_nd_lt30k"
    assert dict(want_def)["U_2C"] == "not_representative"

    rng = random.Random(6201)
    for i in range(10):
        shuffled = RULE_ROWS[:]
        rng.shuffle(shuffled)
        pq_path = tmp_path / f"shuf_{i}.parquet"
        _write_parquet(pq_path, shuffled)
        sel, deferred, stats = mv.select_mv1(pq_path, work=tmp_path / f"w_{i}")
        assert [(r["uvig"], r["tier"]) for r in sel] == want_sel
        assert sorted((r["uvig"], r["reason"]) for r in deferred) == sorted(want_def)
        assert len(sel) + len(deferred) == len(RULE_ROWS)
        assert stats["total_rows"] == len(RULE_ROWS)
        assert stats["selected_uvigs"] == len(want_sel)


def test_select_mv1_enforces_stage_caps_and_orders_mv2_remainder(tmp_path):
    mv = _load()
    pq_path = tmp_path / "cap.parquet"
    _write_parquet(pq_path, RULE_ROWS)

    # Cap at 550,000 bp: fits U_2A (240k) + U_2B (210k) + U_1B (45k) + U_3A (35k) + U_4A (18k) = 548,000 bp;
    # U_6A (55k) and U_7A (75k) spill over in exact priority order into MV2!
    sel, deferred, _ = mv.select_mv1(pq_path, max_bases=550_000, max_contigs=100, work=tmp_path / "w_bp")
    want_sel, want_def = brute_select(RULE_ROWS, max_bases=550_000, max_contigs=100)
    assert [(r["uvig"], r["tier"]) for r in sel] == want_sel
    assert [r["uvig"] for r in sel] == ["U_2A", "U_2B", "U_1B", "U_3A", "U_4A"]
    spill = [r for r in deferred if r["reason"] == "mv2_spillover"]
    assert [r["uvig"] for r in spill] == ["U_6A", "U_7A"]

    # Cap at 3 contigs: fits U_2A, U_2B, U_1B; rest spill over in exact priority order
    sel_c, def_c, _ = mv.select_mv1(pq_path, max_bases=40_000_000_000, max_contigs=3, work=tmp_path / "w_c")
    assert [r["uvig"] for r in sel_c] == ["U_2A", "U_2B", "U_1B"]
    spill_c = [r for r in def_c if r["reason"] == "mv2_spillover"]
    assert [r["uvig"] for r in spill_c] == ["U_3A", "U_4A", "U_6A", "U_7A"]


def test_dec65_dtr_trimming_and_circularity():
    mv = _load()
    unit = seq_of(6211, 2000).encode()
    exact_dtr = unit + unit[:35]
    long_dtr = unit + unit[:1005]
    assert mv.apply_dec65("Direct terminal repeat", exact_dtr) == (unit, True, "exact", 35)
    assert mv.apply_dec65("Direct terminal repeat", unit) == (unit, True, "none", 0)
    assert mv.apply_dec65("Direct terminal repeat", long_dtr) == (long_dtr, False, "long", 0)
    for topo in ("Linear", "Provirus", "Inverted terminal repeat", "Virus MAG", ""):
        assert mv.apply_dec65(topo, exact_dtr) == (exact_dtr, False, "", 0)


def test_stage_mv1_end_to_end_and_determinism(tmp_path):
    mv = _load()
    s_dtr_unit = seq_of(101, 4000)
    s_dtr_exact = s_dtr_unit + s_dtr_unit[:30]
    s_dtr_none = seq_of(102, 4500)
    s_shared = seq_of(103, 5000)
    s_vmag_0 = seq_of(104, 3000)
    s_vmag_1 = seq_of(105, 2500)
    s_jumbo_rep = seq_of(106, 205_000)
    s_jumbo_nonrep = seq_of(107, 201_000)
    s_concat = seq_of(108, 210_000)
    s_lq_short = seq_of(109, 2000)

    rows = [
        make_row("MV_DTR_EXACT", "v1", quality="Complete", length=len(s_dtr_exact), topology="Direct terminal repeat"),
        make_row("MV_DTR_NONE", "v2", quality="Complete", length=len(s_dtr_none), topology="Direct terminal repeat"),
        make_row("MV_SHARED", "v3", quality="High-quality", length=len(s_shared), topology="Provirus"),
        make_row("MV_DUP", "v4", quality="High-quality", length=len(s_shared), topology="Linear"),  # repeats MV_SHARED
        make_row("MV_VMAG", "v5", quality="High-quality", length=len(s_vmag_0) + len(s_vmag_1), topology="Virus MAG"),
        make_row("MV_JUMBO_REP", "v6", quality="High-quality", length=len(s_jumbo_rep), topology="Linear", genetic_code="15"),
        make_row("MV_JUMBO_NONREP", "v6", quality="Medium-quality", length=len(s_jumbo_nonrep), topology="Linear", genetic_code="4"),
        make_row("MV_CONCAT", "v7", quality="Complete", length=len(s_concat), is_concatemer=True),
        make_row("MV_LQ_SHORT", "v8", quality="Low-quality", length=len(s_lq_short)),
    ]
    pq_path = tmp_path / "METAVR_main_table.parquet"
    _write_parquet(pq_path, rows)

    fa = io.StringIO()
    def _w(hdr, seq):
        fa.write(f">{hdr}\n")
        for i in range(0, len(seq), 80):
            fa.write(seq[i:i + 80] + "\n")

    _w("MV_DTR_EXACT|1|s1", s_dtr_exact)
    _w("MV_DTR_NONE|1|s2", s_dtr_none)
    _w("MV_SHARED|1|s3", s_shared)
    _w("MV_DUP|1|s4", s_shared)
    _w("MV_VMAG|1|c0", s_vmag_0)
    _w("MV_VMAG|1|c1", s_vmag_1)
    _w("MV_JUMBO_REP|1|s5", s_jumbo_rep)
    _w("MV_JUMBO_NONREP|1|s6", s_jumbo_nonrep)
    _w("MV_CONCAT|1|s7", s_concat)
    _w("MV_LQ_SHORT|1|s8", s_lq_short)

    fna_bgz = tmp_path / "METAVR.fna.bgz"
    fna_bgz.write_bytes(gzip.compress(fa.getvalue().encode(), mtime=0))

    # Fake V2 owner holding MV_SHARED
    v2_dir = tmp_path / "V2"
    v2_dir.mkdir()
    (v2_dir / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s_shared)]}]))

    out1 = tmp_path / "MV1_run1"
    out2 = tmp_path / "MV1_run2"
    rep1 = mv.stage_mv1(
        pq_path, fna_bgz, out1,
        root=None, root_stages=(), owner_inputs={"V2": v2_dir},
        retrieved_at="2026-10-23T15:22:00Z",
        work=tmp_path / "work1",
    )
    rep2 = mv.stage_mv1(
        pq_path, fna_bgz, out2,
        root=None, root_stages=(), owner_inputs={"V2": v2_dir},
        retrieved_at="2026-10-23T15:22:00Z",
        work=tmp_path / "work2",
    )

    genomes = json.loads((out1 / "genomes.json").read_text())
    assert {g["kind"] for g in genomes} == {"virus"}
    assert {g["genetic_code"] for g in genomes} == {11}
    assert {g["assembly_quality"] for g in genomes} <= {"vOTU", "contig"}
    # MV_DUP repeats MV_SHARED, so 6 unique genomes listed out of 7 selected UViGs
    assert len(genomes) == 6
    assert [g["contig_uids"] for g in genomes] == [
        [cuid(s_dtr_unit)],
        [cuid(s_dtr_none)],
        [cuid(s_shared)],
        [cuid(s_vmag_0), cuid(s_vmag_1)],
        [cuid(s_jumbo_rep)],
        [cuid(s_jumbo_nonrep)],
    ]
    assert sorted(u for g in genomes for u in g["circular"]) == sorted([cuid(s_dtr_unit), cuid(s_dtr_none)])

    objs = json.loads((out1 / "objects.json").read_text())
    assert len(objs) == 1
    assert objs[0]["redistributable"] is False
    assert objs[0]["source_name"] == "metavr"

    trims = (out1 / "trims.tsv").read_text().splitlines()
    assert len(trims) == 2
    assert "MV_DTR_EXACT" in trims[1]
    assert cuid(s_dtr_unit) in trims[1]

    shared = (out1 / "shared_records.tsv").read_text().splitlines()
    assert len(shared) == 3  # header + MV_SHARED + MV_DUP
    assert all(ln.endswith("\tV2") for ln in shared[1:])

    repeated = (out1 / "repeated_records.tsv").read_text().splitlines()
    assert len(repeated) == 2
    assert "MV_DUP" in repeated[1]

    left_lines = gzip.decompress((out1 / "MV1_left_for_later.tsv.gz").read_bytes()).decode().splitlines()
    assert len(left_lines) == 3  # header + MV_CONCAT + MV_LQ_SHORT
    assert rep1["trims"] == 1
    assert rep1["trimmed_bases"] == 30
    assert rep1["length_mismatches_count"] == 0

    # Byte-for-byte determinism across runs (except objects.json path)
    for fname in ("genomes.json", "metavr_mv1.fna", "records.tsv", "trims.tsv",
                  "repeated_records.tsv", "shared_records.tsv", "MV1_left_for_later.tsv.gz"):
        assert (out1 / fname).read_bytes() == (out2 / fname).read_bytes(), f"{fname} differs across runs"
    assert rep1["sha256"]["genomes.json"] == rep2["sha256"]["genomes.json"]

    # s2_check_inputs.py passes on staged MV1
    chk = subprocess.run(
        [sys.executable, str(SCRIPTS / "s2_check_inputs.py"), "MV1",
         str(out1 / "genomes.json"), str(out1 / "objects.json")],
        capture_output=True, text=True, check=False,
    )
    assert chk.returncode == 0, chk.stdout + chk.stderr

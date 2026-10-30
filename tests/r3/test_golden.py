"""Ring R3: the golden set through O-INGEST and O-CALL.

Runs only when BIOBOX_GOLDEN_DIR points at a fetched golden set
(`python -m biobox.fetch golden.yaml DIR --lock golden.lock.yaml`); skipped in CI.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import shutil
import time
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, contracts, fasta, ids, ledger, provenance
from biobox.l0seq import SeqStore, SourceObject, ingest
from biobox.l1 import igr
from biobox.l3 import pairs, tapes, tokens
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, call_genome, run

pytestmark = pytest.mark.r3

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = yaml.safe_load((ROOT / "golden.yaml").read_text())
LOCK = yaml.safe_load((ROOT / "golden.lock.yaml").read_text()) if (ROOT / "golden.lock.yaml").exists() else None
DIR = os.environ.get("BIOBOX_GOLDEN_DIR")
if not DIR or LOCK is None:
    pytest.skip("golden set not fetched (set BIOBOX_GOLDEN_DIR)", allow_module_level=True)
DIR = Path(DIR)

_BASES = "TCAG"
_STD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_AA = {
    11: _STD,
    4: _STD[:14] + "W" + _STD[15:],
    15: _STD[:11] + "Q" + _STD[12:],
    25: _STD[:14] + "G" + _STD[15:],
}


def _translate(nt: bytes, code: int) -> str:
    aa = _AA[code]
    out = []
    for i in range(0, len(nt) - 2, 3):
        c = nt[i : i + 3].decode().upper()
        out.append(
            "X" if any(b not in _BASES for b in c) else aa[16 * _BASES.index(c[0]) + 4 * _BASES.index(c[1]) + _BASES.index(c[2])]
        )
    return "".join(out)


def _records(gid):
    return fasta.read(DIR / LOCK["files"][gid]["path"])


def _input(entry, code=None, kind=None):
    recs = _records(entry["id"])
    uids = tuple(ids.contig_uid(r.seq) for r in recs)
    circ = tuple(u for u, spec in zip(uids, entry["records"]) if spec["topology"] == "circular")
    return GenomeInput(uids, code or entry["genetic_code"], kind=kind or entry["kind"], circular=circ)


@pytest.fixture(scope="module")
def store_root(tmp_path_factory):
    root = tmp_path_factory.mktemp("golden")
    objs = [SourceObject(DIR / LOCK["files"][e["id"]]["path"], "golden", str(GOLDEN["version"])) for e in GOLDEN["genomes"]]
    t0 = time.perf_counter()
    ingest(objs, root, "S0")
    t1 = time.perf_counter()
    inputs = [_input(e) for e in GOLDEN["genomes"]]
    report = run(root, "S0", inputs)
    t2 = time.perf_counter()
    igrs = igr.run(root, "S0", {u for gi in inputs for u in gi.circular})
    tp = tapes.run(root, "S0", {u for gi in inputs for u in gi.circular})
    print(f"\ngolden: ingest {t1 - t0:.1f}s, gene calling {t2 - t1:.1f}s, {report}, igr {igrs}, tapes {tp}")
    if shutil.which("mmseqs"):
        from biobox.l2 import cluster
        cluster.run(root, "rS0.golden", GENE_CALL_VERSION, threads=4, tmp=root.parent)
        tokens.run(root, "rS0.golden")
        print("pairs", pairs.run(root, "rS0.golden"))
    return root


def test_files_match_lock():
    for gid, f in LOCK["files"].items():
        with gzip.open(DIR / f["path"], "rb") as fh:
            assert hashlib.sha256(fh.read()).hexdigest() == f["sha256_fna"], gid


def test_every_record_has_genes(store_root):
    occ = pq.read_table(
        store_root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet"
    ).to_pylist()
    per_contig = Counter(r["contig_uid"] for r in occ)
    for e in GOLDEN["genomes"]:
        for r in _records(e["id"]):
            assert per_contig[ids.contig_uid(r.seq)] >= 1, (e["id"], r.header[:40])


def test_prokaryote_gene_density_is_plausible(store_root):
    store = SeqStore(store_root)
    for e in GOLDEN["genomes"]:
        if e["kind"] != "prokaryote":
            continue
        calls = call_genome(store, _input(e))
        mb = sum(len(r.seq) for r in _records(e["id"])) / 1e6
        assert 700 <= len(calls) / mb <= 1250, (e["id"], len(calls) / mb)


@pytest.mark.parametrize("gid,expected", [("lak_a1", 15), ("lak_b8", 15), ("lambda", 11), ("t4", 11)])
def test_gv_picks_the_phage_genetic_code(store_root, gid, expected):
    # lak_b8: NCBI says code 11, the sequence says 15 (golden.yaml note); lambda and t4 are code-11 controls.
    e = next(e for e in GOLDEN["genomes"] if e["id"] == gid)
    calls = call_genome(SeqStore(store_root), _input(e))
    codes = Counter(c["genetic_code"] for c in calls)
    assert codes.most_common(1)[0][0] == expected, codes
    assert codes[expected] / sum(codes.values()) >= 0.95, codes


@pytest.mark.parametrize("gid", ["lak_a1", "lak_b8", "lambda", "t4"])
def test_phage_coding_density_under_called_code(store_root, gid):
    e = next(e for e in GOLDEN["genomes"] if e["id"] == gid)
    calls = call_genome(SeqStore(store_root), _input(e))
    n = sum(len(r.seq) for r in _records(gid))
    covered = sum(c["end"] - c["start"] for c in calls)
    assert covered / n >= 0.85, (gid, covered / n)


def test_sr1_needs_code_25(store_root):
    e = next(e for e in GOLDEN["genomes"] if e["id"] == "sr1_raac1")
    store = SeqStore(store_root)
    c25 = call_genome(store, _input(e, code=25))
    c11 = call_genome(store, _input(e, code=11))
    mean = lambda cs: sum(len(c["aa"]) for c in cs) / len(cs)  # noqa: E731
    assert mean(c25) > 1.3 * mean(c11), (mean(c25), mean(c11))


def test_obelisk_oblin1_is_called(store_root):
    e = next(e for e in GOLDEN["genomes"] if e["id"] == "obelisk_ss")
    calls = call_genome(SeqStore(store_root), _input(e))
    recs = _records("obelisk_ss")
    for truth in e["truth"]:
        rec = next(r for r in recs if truth["record"] in r.header)
        cu = ids.contig_uid(rec.seq)
        lo, hi = truth["coords_1based_inclusive"][0] - 1, truth["coords_1based_inclusive"][1]
        best = max(
            (min(hi, c["end"]) - max(lo, c["start"]) for c in calls if c["contig_uid"] == cu and c["strand"] == 1),
            default=0,
        )
        assert best >= 0.9 * (hi - lo), (truth, best)


def test_r04_independent_translation_on_golden(store_root):
    store = SeqStore(store_root)
    checked = 0
    for e in GOLDEN["genomes"]:
        gi = _input(e)
        for c in call_genome(store, gi):
            if c["dist_to_contig_end_bp"] == 0:
                continue
            nt = store.get(c["contig_uid"], c["start"], c["end"], c["strand"], circular=c["contig_uid"] in gi.circular)
            aa = _translate(nt, c["genetic_code"])
            assert aa.endswith("*"), (e["id"], c["start"])
            assert c["aa"] == "M" + aa[1:-1], (e["id"], c["start"])
            checked += 1
    assert checked > 30_000


@pytest.fixture(scope="module")
def golden_catalog(store_root, tmp_path_factory):
    man = catalog.build_manifest(store_root, "rS0.golden")
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path_factory.mktemp("cat") / "c.duckdb"))
    yield man, con
    con.close()


@pytest.mark.covers("LY-02", "LY-03")
def test_golden_layout_sort_and_statistics(golden_catalog):
    man, _ = golden_catalog
    assert contracts.check_sort_order(man) == []
    assert contracts.check_statistics(man) == []


@pytest.mark.covers("C-03", "C-05", "R-06")
def test_golden_contracts(golden_catalog):
    _, con = golden_catalog
    assert contracts.check_foreign_keys(con) == []
    assert contracts.check_coordinates(con) == []
    assert contracts.check_lineage(con) == []
    wraps = con.execute('SELECT count(*) FROM protein_occurrence o JOIN contig_index ci USING (contig_uid) '
                        'WHERE o."end" > ci.contig_len').fetchall()[0][0]
    assert wraps >= 1  # origin-wrapping genes exist, so C-05's wrap rule is exercised


@pytest.mark.covers("LY-04")
def test_golden_ly04(golden_catalog):
    man, _ = golden_catalog
    problems = contracts.check_file_sizes(man)
    print("\nLY-04:", problems)
    assert problems == []  # one shard and one protein_prefix=all file for the golden stage


@pytest.mark.covers("C-08")
def test_golden_provenance(golden_catalog):
    man, _ = golden_catalog
    assert provenance.check_provenance(man) == []


@pytest.mark.covers("C-03", "LY-02")
def test_golden_shards(golden_catalog):
    man, _ = golden_catalog
    assert contracts.check_shards(man) == []
    shards = [a for a in man["artifacts"] if a["path"].endswith(".fa.gz")]
    assert len(shards) == 1


@pytest.mark.covers("L-01", "L-02", "L-03")
def test_golden_ledger_reconciles(store_root):
    assert ledger.check_balance(store_root, "S0", ["ingest", "genecall"]) == []
    expected = {ids.sha256_file(DIR / LOCK["files"][e["id"]]["path"]) for e in GOLDEN["genomes"]}
    assert ledger.check_inputs(store_root, "S0", expected) == []
    assert ledger.check_chain(store_root, "S0") == []
    gc = ledger.read(store_root, "S0", "genecall")
    assert sum(r["n_in"] for r in gc) == 33_716 and sum(r["n_kept"] for r in gc) == 33_442


@pytest.mark.covers("C-04", "C-05", "C-07", "Q2")
def test_golden_tapes(golden_catalog, store_root):
    man, con = golden_catalog
    assert contracts.check_list_lengths(con) == []
    assert contracts.check_saturation(con) == []
    assert contracts.check_coordinates(con) == []
    n_slots, wraps = con.execute("SELECT sum(n_slots), count(*) FILTER (WHERE list_max(end_bps) > contig_len) "
                                 "FROM contig_tapes").fetchall()[0]
    assert n_slots == 33_716 and wraps >= 1
    occ = {}
    for f in store_root.glob(f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet"):
        for r in pq.read_table(f).to_pylist():
            occ.setdefault(r["contig_uid"], []).append(r)
    cuids = sorted(occ)[::3]
    for g in con.sql("SELECT contig_uid, start_bps, protein_uids FROM tapes_of(?)", params=[cuids]).to_arrow_table().to_pylist():
        slots = sorted(occ[g["contig_uid"]], key=lambda r: r["slot_ordinal"])
        assert g["start_bps"] == [r["start"] for r in slots]
        assert g["protein_uids"] == [r["protein_uid"] for r in slots]


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="MMseqs2 not on PATH")
@pytest.mark.covers("Q1", "C-03")
def test_golden_pairs(golden_catalog):
    _, con = golden_catalog
    assert contracts.check_foreign_keys(con) == []
    a, b, n = con.execute("SELECT anchor_token, partner_token, count(DISTINCT (contig_uid, anchor_ordinal)) "
                          "FROM slot_pairs WHERE abs(orf_offset) <= 3 GROUP BY 1, 2 ORDER BY 3 DESC, 1, 2 LIMIT 1").fetchall()[0]
    assert con.execute("SELECT q1_count(?, ?, 3, NULL)", [a, b]).fetchall()[0][0] == n

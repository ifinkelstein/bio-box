"""Card C-48e part 1: the first SPIRE tranche (selection, fetch, stage inputs).

Covers:
1. `c48e_select.py`: a tranche chosen as a pure function of the tables: habitats filtered, level-2 strata
   interleaved in proportion to sqrt(N_h) by a divisor (Webster) sequence checked against an independent sort of
   all quotients, studies taken in turn within a stratum, the two caps respected by the longest prefix of the
   order, a later tranche excluding the samples an earlier one took and continuing the same order, and the
   nearest-N50 estimate of contigs and bases at 5 kb checked against a hand computation.
2. `c48e_fetch.py` against a local HTTP server: a log line per sample with sizes and hashes; a plain FASTA stored
   gzipped; multi-member gzip read whole; HTTP 404 recorded and not retried; a transient error retried and a
   sample that fails three tries listed, not replaced; HTTP 429 stops the whole fetch; a rerun skips what is
   complete and verified and fetches again a file whose bytes no longer match the log; the disk guard stops
   before any request.
3. `c48e_inputs.py`: the stage's inputs as the pilot's depth-1 stages had them (one prokaryote `contig` genome per
   distinct contig, here only those of at least `--min-contig-len`, the assemblies staged whole by hard link) checked
   against sets computed in the test from the sequences; byte-identical files whatever the selection's order; a
   file whose sha256 is not the log's refused; and a round trip through `write_genome_input` and `ingest` with
   `min_contig_len`, whose ledger counts the dropped records and bases the test computes.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import math
import random
import sys
import threading
from collections import Counter
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-48e")

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


# --- 1. c48e_select.py -----------------------------------------------------------------------------


def _row(sid: str, study: str, hab2: str, bp: int | None, n50: int | None = 1000) -> dict:
    return {"sample_id": sid, "study": study, "habitat": hab2.split(":", 1)[0], "habitat2": hab2,
            "assembly_bp": bp, "assembly_n50": n50 if bp is not None else None}


def _table() -> list[dict]:
    """10 marine (1 without a published size), 5 benthic, 17 soil, 1 level-1-only terrestrial, 10 animal host."""
    rows = []
    marine_studies = ["Sm1"] * 5 + ["Sm2"] * 3 + ["Sm3"] * 2
    for i, st in enumerate(marine_studies):
        rows.append(_row(f"MAR{i:02d}", st, "aquatic:marine", None if i == 9 else 100_000_000 + i * 7_000_000, 800 + 90 * i))
    for i in range(5):
        rows.append(_row(f"BEN{i:02d}", "Sb1" if i < 3 else "Sb2", "aquatic:benthic", 50_000_000 + i * 3_000_000, 1500 + 10 * i))
    for i in range(17):
        rows.append(_row(f"SOI{i:02d}", f"Ss{i % 4}", "terrestrial:soil", 200_000_000 + i * 11_000_000, 600 + 40 * i))
    rows.append(_row("TER00", "St1", "terrestrial", 30_000_000, 2000))
    for i in range(10):
        rows.append(_row(f"GUT{i:02d}", "Sg1", "host-associated:animal host", 90_000_000, 3000))
    return rows


def _webster_oracle(avail: dict[str, int], weights: dict[str, float], k: int) -> Counter:
    """Independent check: the k largest quotients w_h / (2j + 1), j < avail_h, by one global sort."""
    q = sorted(((weights[h] / (2 * j + 1), h) for h in avail for j in range(avail[h])), reverse=True)
    return Counter(h for _, h in q[:k])


def test_c48e_order_interleaves_strata_by_sqrt_webster_and_studies_in_turn():
    m = _load("c48e_select")
    rows = _table()
    order, unsized = m.tranche_order(rows, habitats=("aquatic", "terrestrial"), seed=7)

    assert [r["sample_id"] for r in unsized] == ["MAR09"]
    assert all(r["habitat"] in ("aquatic", "terrestrial") for r in order)
    assert len(order) == 9 + 5 + 17 + 1
    assert [r["order_rank"] for r in order] == list(range(1, len(order) + 1))

    # N_h counts every sample of the stratum (as the pilot did), the sized ones are what can be taken.
    n_h = {"aquatic:marine": 10, "aquatic:benthic": 5, "terrestrial:soil": 17, "terrestrial": 1}
    avail = {"aquatic:marine": 9, "aquatic:benthic": 5, "terrestrial:soil": 17, "terrestrial": 1}
    weights = {h: math.sqrt(n) for h, n in n_h.items()}
    for r in order:
        assert r["stratum_n"] == n_h[r["habitat2"]]
    for k in range(1, len(order) + 1):
        got = Counter(r["habitat2"] for r in order[:k])
        assert got == _webster_oracle(avail, weights, k), k

    # Within a stratum, studies are taken in turn: the first three marine picks are three different studies,
    # and the first four soil picks are four different studies.
    marine = [r for r in order if r["habitat2"] == "aquatic:marine"]
    soil = [r for r in order if r["habitat2"] == "terrestrial:soil"]
    assert len({r["study"] for r in marine[:3]}) == 3
    assert len({r["study"] for r in soil[:4]}) == 4
    assert [r["stratum_rank"] for r in marine] == list(range(1, 10))


def test_c48e_order_is_a_pure_function_of_the_table_and_seed():
    m = _load("c48e_select")
    rows = _table()
    a, ua = m.tranche_order(rows, habitats=("aquatic", "terrestrial"), seed=7)
    shuffled = list(rows)
    random.Random(3).shuffle(shuffled)
    b, ub = m.tranche_order(shuffled, habitats=("terrestrial", "aquatic"), seed=7)
    assert a == b and ua == ub
    c, _ = m.tranche_order(rows, habitats=("aquatic", "terrestrial"), seed=8)
    assert [r["sample_id"] for r in a] != [r["sample_id"] for r in c]


def test_c48e_nearest_n50_estimate_matches_hand_computation_and_keeps_measured_samples():
    m = _load("c48e_select")
    train = [
        {"sample_id": "T1", "assembly_n50": 500, "total_bp": 100, "bp_ge_5000": 10, "n_ge_5000": 1},
        {"sample_id": "T2", "assembly_n50": 1000, "total_bp": 300, "bp_ge_5000": 60, "n_ge_5000": 4},
        {"sample_id": "T3", "assembly_n50": 4000, "total_bp": 200, "bp_ge_5000": 120, "n_ge_5000": 9},
        {"sample_id": "T4", "assembly_n50": 900, "total_bp": 50, "bp_ge_5000": 5, "n_ge_5000": 1},
    ]
    rows = [
        {"sample_id": "X", "assembly_bp": 1_000_000, "assembly_n50": 950},
        {"sample_id": "T3", "assembly_bp": 200, "assembly_n50": 4000},
    ]
    out = {r["sample_id"]: r for r in m.estimate_5k(rows, train, k=2)}
    # X: |ln 950 - ln n50| is least for T4 (900) then T2 (1000): ratios over their summed bases (350).
    assert out["X"]["est_source"] == "nearest_n50"
    assert out["X"]["est_bp_ge_5000"] == round(1_000_000 * (60 + 5) / 350)
    assert out["X"]["est_n_ge_5000"] == round(1_000_000 * (4 + 1) / 350)
    # A sample the pilot measured keeps its measured values.
    assert out["T3"]["est_source"] == "measured_pilot"
    assert (out["T3"]["est_bp_ge_5000"], out["T3"]["est_n_ge_5000"]) == (120, 9)


def _with_estimates(m, order):
    return [{**r, "est_n_ge_5000": r["assembly_bp"] // 10_000_000, "est_bp_ge_5000": r["assembly_bp"] // 8}
            for r in order]


def test_c48e_take_is_the_longest_prefix_within_both_caps():
    m = _load("c48e_select")
    order, _ = m.tranche_order(_table(), habitats=("aquatic", "terrestrial"), seed=7)
    est = _with_estimates(m, order)
    for max_n, max_bp in ((100, 10**12), (10**9, 400_000_000), (60, 300_000_000), (0, 0), (10**9, 10**12)):
        taken, rest = m.take(est, max_contigs=max_n, max_bp=max_bp)
        assert [r["sample_id"] for r in taken + rest] == [r["sample_id"] for r in est]
        n = sum(r["est_n_ge_5000"] for r in taken)
        bp = sum(r["est_bp_ge_5000"] for r in taken)
        assert n <= max_n and bp <= max_bp
        if rest:
            nxt = rest[0]
            assert n + nxt["est_n_ge_5000"] > max_n or bp + nxt["est_bp_ge_5000"] > max_bp
        if taken:
            assert (taken[-1]["cum_est_n_ge_5000"], taken[-1]["cum_est_bp_ge_5000"]) == (n, bp)


def test_c48e_a_later_tranche_excludes_the_earlier_one_and_continues_its_order():
    m = _load("c48e_select")
    order, _ = m.tranche_order(_table(), habitats=("aquatic", "terrestrial"), seed=7)
    est = _with_estimates(m, order)
    t1, rest1 = m.take(est, max_contigs=80, max_bp=500_000_000)
    assert t1 and rest1
    t2, rest2 = m.take(est, max_contigs=80, max_bp=500_000_000, exclude={r["sample_id"] for r in t1})
    ids1 = {r["sample_id"] for r in t1}
    assert t2 and not ids1 & {r["sample_id"] for r in t2 + rest2}
    assert [r["sample_id"] for r in t2] == [r["sample_id"] for r in rest1[: len(t2)]]
    assert [r["sample_id"] for r in t2 + rest2] == [r["sample_id"] for r in rest1]


# --- 2. c48e_fetch.py ------------------------------------------------------------------------------

_FA1 = b">k1 flag=1\n" + b"ACGT" * 3000 + b"\n"
_FA2 = b">k2\n" + b"TTGCA" * 2000 + b"\n"
_PLAIN = b">p1\n" + b"GATTACA" * 700 + b"\n"


class _Server:
    """A local stand-in for the provider: per-path behaviours, every request counted."""

    def __init__(self):
        self.hits: Counter = Counter()
        self.flaky_left = 2
        srv = self

        class H(BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def _send(self, code, body=b"", headers=()):
                self.send_response(code)
                for k, v in headers:
                    self.send_header(k, v)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def do_GET(self):
                p = self.path
                srv.hits[p] += 1
                assert self.headers.get("User-Agent") == "curl"
                if p == "/download_assembly/OK1":
                    return self._send(302, headers=[("Location", "/files/OK1.fa.gz")])
                if p == "/files/OK1.fa.gz":  # two gzip members, as some servers concatenate
                    return self._send(200, gzip.compress(_FA1, mtime=0) + gzip.compress(_FA2, mtime=0))
                if p == "/download_assembly/PLAIN":
                    return self._send(200, _PLAIN)
                if p == "/download_assembly/FLAKY":
                    if srv.flaky_left > 0:
                        srv.flaky_left -= 1
                        return self._send(503)
                    return self._send(200, gzip.compress(_FA2, mtime=0))
                if p == "/download_assembly/BAD":
                    return self._send(200, gzip.compress(_FA1, mtime=0)[:-20])
                if p == "/download_assembly/RATE":
                    return self._send(429)
                return self._send(404)

        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), H)
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}/download_assembly/{{sample_id}}"
        self.t = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self.t.start()

    def close(self):
        self.httpd.shutdown()
        self.httpd.server_close()


@pytest.fixture
def server():
    s = _Server()
    yield s
    s.close()


def _fetch(m, server, tmp_path, sids, **kw):
    args = dict(url_template=server.url, interval=0.5, tries=3, retry_wait=0.0, max_dl_bytes=10**9,
                min_free_bytes=0, stop_files=[tmp_path / "STOP"])
    args.update(kw)
    return m.fetch_all(sids, tmp_path / "dl", **args)


def _log(tmp_path):
    return [json.loads(x) for x in (tmp_path / "dl" / "fetch_log.jsonl").read_text().splitlines()]


def test_c48e_fetch_logs_hashes_retries_and_lists_failures(tmp_path, server):
    m = _load("c48e_fetch")
    res = _fetch(m, server, tmp_path, ["OK1", "PLAIN", "MISS", "FLAKY", "BAD"])
    assert res["ok"] == ["FLAKY", "OK1", "PLAIN"]
    assert res["not_found"] == ["MISS"] and res["failed"] == ["BAD"]
    final = {r["sample_id"]: r for r in _log(tmp_path) if r["status"] != "retry"}
    for sid in ("OK1", "PLAIN", "FLAKY"):
        f = tmp_path / "dl" / f"{sid}.fa.gz"
        assert final[sid]["disk_sha256"] == hashlib.sha256(f.read_bytes()).hexdigest()
        assert final[sid]["disk_bytes"] == f.stat().st_size
    with gzip.open(tmp_path / "dl" / "OK1.fa.gz") as fh:  # both members kept
        assert fh.read() == _FA1 + _FA2
    assert final["OK1"]["uncompressed_bytes"] == len(_FA1) + len(_FA2)
    assert final["OK1"]["final_url"].endswith("/files/OK1.fa.gz")
    assert gzip.decompress((tmp_path / "dl" / "PLAIN.fa.gz").read_bytes()) == _PLAIN
    assert final["PLAIN"]["wire_gzipped"] is False
    assert final["FLAKY"]["attempt"] == 3 and server.hits["/download_assembly/FLAKY"] == 3
    assert final["BAD"]["status"] == "failed" and server.hits["/download_assembly/BAD"] == 3
    assert server.hits["/download_assembly/MISS"] == 1
    assert not (tmp_path / "dl" / "BAD.fa.gz").exists()
    assert not list((tmp_path / "dl").glob(".*part"))


def test_c48e_fetch_rerun_skips_verified_and_refetches_a_changed_file(tmp_path, server):
    m = _load("c48e_fetch")
    _fetch(m, server, tmp_path, ["OK1", "MISS", "FLAKY"])
    before = Counter(server.hits)
    res = _fetch(m, server, tmp_path, ["OK1", "MISS", "FLAKY"])
    assert server.hits == before
    assert res["skipped_verified"] == ["FLAKY", "OK1"] and res["not_found"] == ["MISS"]
    f = tmp_path / "dl" / "OK1.fa.gz"
    f.write_bytes(gzip.compress(b">x\nAC\n", mtime=0))
    res = _fetch(m, server, tmp_path, ["OK1", "MISS", "FLAKY"])
    assert res["ok"] == ["OK1"] and server.hits["/files/OK1.fa.gz"] == before["/files/OK1.fa.gz"] + 1
    with gzip.open(f) as fh:
        assert fh.read() == _FA1 + _FA2


def test_c48e_fetch_stops_on_rate_limit_and_on_the_disk_guard(tmp_path, server):
    m = _load("c48e_fetch")
    with pytest.raises(m.Stop) as e:
        _fetch(m, server, tmp_path, ["RATE", "OK1"])
    assert e.value.code == 4
    assert server.hits["/download_assembly/OK1"] == 0
    assert [r["status"] for r in _log(tmp_path)] == ["refused"]
    with pytest.raises(m.Stop) as e:
        _fetch(m, server, tmp_path / "g", ["OK1"], min_free_bytes=10**18)
    assert e.value.code == 3 and server.hits["/download_assembly/OK1"] == 0
    (tmp_path / "STOP").write_text("")
    with pytest.raises(m.Stop) as e:
        _fetch(m, server, tmp_path, ["OK1"])
    assert e.value.code == 75


# --- 3. c48e_inputs.py -----------------------------------------------------------------------------


def _seq(n: int, seed: int) -> bytes:
    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(n))


def _inputs_fixture(tmp_path: Path):
    from biobox import ids

    x, y = _seq(6000, 1), _seq(5000, 2)  # at or over the length
    short = [_seq(1200, 3), _seq(4999, 4), _seq(800, 5)]
    files = {
        "SAM_A": b">a1 len=6000\n" + x + b"\n>a2\n" + short[0] + b"\n>a3\n" + y[:2500] + b"\n" + y[2500:] + b"\n",
        "SAM_B": b">b1\n" + x.lower() + b"\n>b2\n" + short[1] + b"\n",
        "SAM_C": b">c1\n" + short[2] + b"\n",
    }
    dl = tmp_path / "dl"
    dl.mkdir()
    log = []
    for sid, data in files.items():
        f = dl / f"{sid}.fa.gz"
        f.write_bytes(gzip.compress(data, mtime=0))
        log.append({"sample_id": sid, "status": "ok", "disk_bytes": f.stat().st_size,
                    "disk_sha256": hashlib.sha256(f.read_bytes()).hexdigest(), "fetched_utc": "2026-11-20T17:00:00Z"})
    log.append({"sample_id": "SAM_D", "status": "not_found", "fetched_utc": "2026-11-20T17:00:05Z"})
    (dl / "fetch_log.jsonl").write_text("".join(json.dumps(r) + "\n" for r in log))
    sel = [{"order_rank": i + 1, "sample_id": s, "study": f"St{i}", "habitat": "aquatic", "habitat2": "aquatic:marine"}
           for i, s in enumerate(["SAM_C", "SAM_A", "SAM_D", "SAM_B"])]
    truth = {
        "genomes": sorted({ids.contig_uid(x).hex(), ids.contig_uid(y).hex()}),
        "records": 6,
        "bp": 6000 + 1200 + 5000 + 6000 + 4999 + 800,
        "below_records": 3,
        "below_bp": 1200 + 4999 + 800,
    }
    return dl, sel, truth


def _write_sel(path: Path, rows: list[dict]) -> Path:
    cols = list(rows[0].keys())
    path.write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))
    return path


def test_c48e_inputs_one_genome_per_long_contig_whole_assemblies_and_byte_identical(tmp_path):
    m = _load("c48e_inputs")
    dl, sel, truth = _inputs_fixture(tmp_path)
    final = tmp_path / "final" / "SPT"
    outs = []
    for i, rows in enumerate((sel, list(reversed(sel)))):
        s = _write_sel(tmp_path / f"sel{i}.tsv", rows)
        out = tmp_path / f"run{i}" / "SPT"
        rep = m.build(s, dl / "fetch_log.jsonl", dl, "SPT", final, out, min_contig_len=5000,
                      uids_out=tmp_path / f"run{i}" / "uids.parquet", tmp=tmp_path / "tmp")
        outs.append(out)
    a, b = outs
    for name in ("genomes.json", "objects.json", "samples.tsv", "input_report.json", "assemblies.sha256"):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name
    assert (tmp_path / "run0" / "uids.parquet").read_bytes() == (tmp_path / "run1" / "uids.parquet").read_bytes()

    genomes = json.loads((a / "genomes.json").read_text())
    assert [g["contig_uids"][0] for g in genomes] == truth["genomes"]
    assert all(g == {"contig_uids": g["contig_uids"], "genetic_code": 11, "kind": "prokaryote", "circular": [],
                     "assembly_quality": "contig"} for g in genomes)
    objs = json.loads((a / "objects.json").read_text())
    assert [o["path"] for o in objs] == [str(final / "assemblies" / f"{s}.fa.gz") for s in ("SAM_A", "SAM_B", "SAM_C")]
    assert {o["source_name"] for o in objs} == {"spire"} and all(o["redistributable"] for o in objs)
    for s in ("SAM_A", "SAM_B", "SAM_C"):
        assert (a / "assemblies" / f"{s}.fa.gz").stat().st_ino == (dl / f"{s}.fa.gz").stat().st_ino
    assert rep["records"] == truth["records"] and rep["bp"] == truth["bp"]
    assert (rep["records_below"], rep["bp_below"]) == (truth["below_records"], truth["below_bp"])
    assert rep["records_at_least"] == 3 and rep["distinct_contigs"] == 2 and rep["contigs_in_2plus_samples"] == 1
    assert rep["samples_staged"] == 3 and rep["not_staged"] == [{"sample_id": "SAM_D", "status": "not_found"}]
    assert not (tmp_path / "run0" / ".SPT.tmp").exists()
    with pytest.raises(m.Refused):
        m.build(_write_sel(tmp_path / "sel0.tsv", sel), dl / "fetch_log.jsonl", dl, "SPT", final, a,
                min_contig_len=5000, uids_out=tmp_path / "u.parquet", tmp=tmp_path / "tmp")


def test_c48e_inputs_refuse_a_file_that_differs_from_the_log(tmp_path):
    m = _load("c48e_inputs")
    dl, sel, _ = _inputs_fixture(tmp_path)
    (dl / "SAM_B.fa.gz").write_bytes(gzip.compress(b">b1\nACGT\n", mtime=0))
    with pytest.raises(m.Refused, match="SAM_B"):
        m.build(_write_sel(tmp_path / "sel.tsv", sel), dl / "fetch_log.jsonl", dl, "SPT", tmp_path / "f" / "SPT",
                tmp_path / "o" / "SPT", min_contig_len=5000, uids_out=tmp_path / "u.parquet", tmp=tmp_path / "tmp")
    assert not (tmp_path / "o" / "SPT").exists()


def test_c48e_inputs_roundtrip_through_write_genome_input_and_ingest_with_min_contig_len(tmp_path):
    from datetime import datetime

    import pyarrow.parquet as pq

    from biobox.l0seq import SourceObject, ingest
    from biobox.l0seq.inputs import GenomeInput, write_genome_input

    m = _load("c48e_inputs")
    dl, sel, truth = _inputs_fixture(tmp_path)
    final = tmp_path / "in" / "SPT"
    m.build(_write_sel(tmp_path / "sel.tsv", sel), dl / "fetch_log.jsonl", dl, "SPT", final, final,
            min_contig_len=5000, uids_out=tmp_path / "u.parquet", tmp=tmp_path / "tmp")
    root = tmp_path / "root"
    gis = [GenomeInput(contig_uids=tuple(bytes.fromhex(u) for u in g["contig_uids"]), genetic_code=g["genetic_code"],
                       kind=g["kind"], circular=(), assembly_quality=g["assembly_quality"])
           for g in json.loads((final / "genomes.json").read_text())]
    assert write_genome_input(root, "SPT", gis) == {"genomes": 2, "contigs": 2}
    objs = [SourceObject(path=Path(o["path"]), source_name=o["source_name"], source_release=o["source_release"],
                         license=o["license"], redistributable=o["redistributable"],
                         retrieved_at=datetime.fromisoformat(o["retrieved_at"]))
            for o in json.loads((final / "objects.json").read_text())]
    res = ingest(objs, root, "SPT", min_contig_len=5000)
    assert res["records"] == truth["records"]
    assert (res["rejected_records"], res["rejected_bases"]) == (truth["below_records"], truth["below_bp"])
    ci = pq.read_table(root / "l0seq" / "contig_index" / "stage=SPT")
    kept = sorted(bytes(u).hex() for u, d in zip(ci["contig_uid"].to_pylist(), ci["is_default_variant"].to_pylist()) if d)
    assert kept == truth["genomes"]


# --- 4. c48e_est.py --------------------------------------------------------------------------------------------------


def test_c48e_est_drops_samples_from_the_end_of_the_order_until_both_caps_hold():
    m = _load("c48e_est")
    sel = [{"order_rank": str(i), "sample_id": f"S{i}", "habitat": "aquatic", "habitat2": "aquatic:marine",
            "est_source": "nearest_n50", "assembly_bp": "1000", "est_n_ge_5000": "10", "est_bp_ge_5000": "100"}
           for i in (3, 1, 2, 4)]
    meas = {"S1": (12, 90), "S2": (9, 140), "S3": (11, 100)}  # S4 not staged
    samples = [{"sample_id": s, "bp": "1000", "records_at_least": str(n), "bp_at_least": str(b)}
               for s, (n, b) in meas.items()]
    got = m.measured_vs_estimate(sel, samples, max_n=25, max_bp=1000)
    assert got["drop_from_end"] == ["S3"]  # 32 contigs > 25: drop the last staged in order; 21 <= 25
    assert (got["kept_sum_n_ge_5000"], got["kept_sum_bp_ge_5000"]) == (21, 230)
    got = m.measured_vs_estimate(sel, samples, max_n=100, max_bp=200)
    assert got["drop_from_end"] == ["S3", "S2"]  # bases: 330 > 200, 230 > 200, 90 <= 200
    allrow = [g for g in got["groups"] if g["level"] == "all"][0]
    assert (allrow["samples"], allrow["est_n_ge_5000"], allrow["meas_n_ge_5000"]) == (3, 30, 32)
    assert allrow["ratio_n_meas_over_est"] == round(32 / 30, 4)

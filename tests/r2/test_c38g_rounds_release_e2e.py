"""Ring R2: card C-38g — §3.8 tests (5) & (6), `nbr summarize` round source, micro-S0 hash, and 100x memory.

Covers:
  - §3.8 test (5): `catalog.profile_query` on a release with 10 rounds merged by `rounds.merge_release`
    reads O(1) = 1 row group of the merged `release=` files for both literal and joined keys (for both
    `domain_state(puid, lib)` and `domains_of(puid)`), reading 0 row groups from the 10 per-round files.
  - §3.8 test (6): family tokens (`token_dict`, `tape_tokens`) and `token_postings` of a mixed release
    are byte-identical with and without rounds; `nbr summarize` shows round rows with `source='round'`;
    `hits_tsv` and `locus_tables` reflect round `measured` and round `none`; `provenance.check_provenance`
    and `catalog.verify_manifest` pass cleanly; and a release without rounds keeps its exact unchanged
    `manifest_hash` (`5fa9e49761393e3468ebc3c1d7d7224a`).
  - `AGENTS.md` §4.2: 100x synthetic store (`scripts/stress/synth_scale.py`) streaming k-way merge across
    5 overlapping rounds stays under 512 MiB peak RSS and <= 64 MiB growth from 10x to 100x.
"""

from __future__ import annotations

import contextlib
import hashlib
import io as stdio
import json
import runpy
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, dim, ids, provenance, runstep
from biobox import io as bio
from biobox.l0seq import inputs as l0inputs
from biobox.l0seq import store
from biobox.l1 import domains, genecall, igr
from biobox.l1 import rounds as l1rounds
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l2 import cluster as l2cluster
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import TABLES, table
from biobox.search import export, nbr

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402
import test_c38a_mixed_release as c38a  # noqa: E402

direct = c08.direct

pytestmark = [pytest.mark.r2, pytest.mark.covers("C-38g", "LY-06", "LY-10", "MR-08", "W-05", "W-06")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
LIB = "pfam-38.2"
SMALL, BIG = 10, 100
GROWTH_BOUND = 64 * 2**20
PEAK_BOUND = 512 * 2**20

CHILD_MERGE = r"""
import json, sys
from pathlib import Path
from biobox.l1 import rounds
root, rel, out, tmp = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
res = rounds.merge_release(root, rel, out=out, tmp=tmp)
print(json.dumps(res, sort_keys=True))
"""


def _puid(i: int) -> bytes:
    return ids.protein_uid(f"M{'Q' * 12}{i:07d}")


def _write_minimal_stage(root: Path, stage: str = "S0") -> None:
    spec = TABLES["contig_index"]
    cu = ids.contig_uid(f"GATTACA{stage}".encode() * 30)
    rd = ids.raw_digest(f"GATTACA{stage}".encode() * 30)
    bio.write_table(
        pa.table({
            "contig_uid": pa.array([cu], pa.binary(16)),
            "raw_digest": pa.array([rd], pa.binary(32)),
            "contig_len": pa.array([len(f"GATTACA{stage}") * 30], pa.int64()),
            "is_default_variant": pa.array([True], pa.bool_()),
        }, schema=spec.schema),
        spec,
        root / spec.location / f"stage={stage}" / "part-0.parquet",
    )


def _dom_table(uids: list[bytes]) -> pa.Table:
    spec = TABLES["round_domains"]
    n = len(uids)
    return pa.table({
        "protein_uid": pa.array(uids, pa.binary(16)),
        "model_acc": pa.array([f"PF{i % 50:05d}.1" for i in range(n)], pa.string()),
        "model_name": pa.array([f"fam_{i % 50}" for i in range(n)], pa.string()),
        "model_len": pa.array([120] * n, pa.uint32()),
        "seq_bitscore": pa.array([60.0] * n, pa.float32()),
        "seq_evalue_zref": pa.array([1e-12] * n, pa.float64()),
        "dom_bitscore": pa.array([55.0] * n, pa.float32()),
        "dom_ievalue_zref": pa.array([1e-11] * n, pa.float64()),
        "env_start": pa.array([2] * n, pa.uint32()),
        "env_end": pa.array([15] * n, pa.uint32()),
        "ali_start": pa.array([2] * n, pa.uint32()),
        "ali_end": pa.array([15] * n, pa.uint32()),
        "hmm_from": pa.array([0] * n, pa.uint32()),
        "hmm_to": pa.array([13] * n, pa.uint32()),
        "cutoff": pa.array(["gathering"] * n, pa.string()),
        "z": pa.array([1e8] * n, pa.float64()),
        "dom_z": pa.array([1e8] * n, pa.float64()),
    }, schema=spec.schema)


def _cov_table(uids: list[bytes], n_rows: list[int]) -> pa.Table:
    return pa.table({
        "protein_uid": pa.array(uids, pa.binary(16)),
        "n_rows": pa.array(n_rows, pa.uint32()),
    }, schema=TABLES["round_coverage"].schema)


def test_c38g_ten_rounds_profile_query_reads_one_row_group(tmp_path):
    """§3.8 test (5): 10 rounds merged into release= files -> O(1) = 1 row group for literal and joined keys."""
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")

    n_rounds, per_round, rg_size = 10, 200, 100
    all_uids: list[bytes] = []
    hit_uids: set[bytes] = set()
    dspec, cspec = TABLES["round_domains"], TABLES["round_coverage"]
    for r in range(n_rounds):
        rid = f"r-{r:016x}"
        uids = sorted(_puid(r * per_round + j) for j in range(per_round))
        all_uids.extend(uids)
        # Even-index proteins get a hit; odd-index get n_rows=0 (round none)
        r_hits = [u for idx, u in enumerate(uids) if idx % 2 == 0]
        hit_uids.update(r_hits)
        bio.write_table(_dom_table(r_hits), dspec,
                        root / dspec.location / f"round={rid}" / f"library={LIB}" / "part-0.parquet")
        bio.write_table(_cov_table(uids, [1 if idx % 2 == 0 else 0 for idx in range(len(uids))]), cspec,
                        root / cspec.location / f"round={rid}" / f"library={LIB}" / "part-0.parquet")

    catalog.define_release(root, "r10", ["S0"])
    res = l1rounds.merge_release(root, "r10", row_group_rows=rg_size)
    assert res["rounds"] == n_rounds
    assert res["libraries"][LIB]["proteins"] == n_rounds * per_round
    assert res["libraries"][LIB]["domains"] == len(hit_uids)

    merged_cov = bio.part_files(root / f"l1/round_coverage/release=r10/library={LIB}")[0]
    merged_dom = bio.part_files(root / f"l1/round_domains/release=r10/library={LIB}")[0]
    assert pq.read_metadata(merged_cov).num_row_groups == 20
    assert pq.read_metadata(merged_dom).num_row_groups == 10

    man = catalog.build_manifest(root, "r10", ["S0"])
    assert catalog.verify_manifest(man) == []
    cat_path = catalog.generate_catalog(man, root / "catalog/release=r10.duckdb")

    sorted_uids = sorted(all_uids)
    probe_hit = next(u for u in sorted_uids[900:1100] if u in hit_uids)
    probe_none = next(u for u in sorted_uids[900:1100] if u not in hit_uids)

    con = catalog.open_catalog(cat_path)
    try:
        # 1. Literal key via domain_state UDF -> 1 row group of merged round_coverage
        q1 = catalog.profile_query(con, "SELECT domain_state(?, ?)", [probe_hit, LIB])
        assert q1["result"] == [("measured",)]
        ds1 = [s for s in q1["scans"] if s.get("lookup") == "bb_domain_state"]
        assert ds1 and sum(s["row_groups_scanned"] for s in ds1) == 1, q1["scans"]
        assert all("release=r10" in p for s in ds1 for p in s["files"]), ds1

        q2 = catalog.profile_query(con, "SELECT domain_state(?, ?)", [probe_none, LIB])
        assert q2["result"] == [("none",)]
        ds2 = [s for s in q2["scans"] if s.get("lookup") == "bb_domain_state"]
        assert ds2 and sum(s["row_groups_scanned"] for s in ds2) == 1, q2["scans"]
        assert all("release=r10" in p for s in ds2 for p in s["files"]), ds2

        # 2. Joined key via domain_state UDF -> 1 row group of merged round_coverage
        con.execute("CREATE TEMP TABLE probes (protein_uid BLOB)")
        con.execute("INSERT INTO probes VALUES (?)", [probe_hit])
        q3 = catalog.profile_query(con, "SELECT protein_uid, domain_state(protein_uid, ?) FROM probes", [LIB])
        assert q3["result"] == [(probe_hit, "measured")]
        ds3 = [s for s in q3["scans"] if s.get("lookup") == "bb_domain_state"]
        assert ds3 and sum(s["row_groups_scanned"] for s in ds3) == 1, q3["scans"]
        assert all("release=r10" in p for s in ds3 for p in s["files"]), ds3

        # 3. domains_of(puid) reads 1 file (the merged release= file) and 1 of its 10 row groups
        q4 = catalog.profile_query(con, "SELECT model_acc, source FROM domains_of(?)", [probe_hit])
        assert len(q4["result"]) == 1 and q4["result"][0][1] == "round"
        rd_scans = [s for s in q4["scans"] if any("round_domains" in p for p in s.get("files", ()))]
        assert rd_scans, q4["scans"]
        assert sum(s["files_read"] for s in rd_scans) == 1, rd_scans
        assert sum(s["row_groups_scanned"] for s in rd_scans) == 1, rd_scans
        assert sum(s["row_groups_total"] for s in rd_scans) == 10, rd_scans
    finally:
        con.close()


def _build_two_stage_root(base: Path, direct) -> tuple[Path, Path]:
    """Build S0 (annotated for test-1) and S1 (called, not annotated) in `base/root`."""
    from biobox import cli

    root = base / "root"
    work = direct["work"]
    # S0: micro-S0 annotated
    assert cli.main(["stage", "inputs", "--root", str(root), "--stage", "S0",
                     "--genomes", str(work.genomes("S0"))]) == 0
    assert cli.main(["stage", "ingest", "--root", str(root), "--stage", "S0",
                     "--objects", str(work.objects("S0"))]) == 0
    assert cli.main(["stage", "genecall", "--root", str(root), "--stage", "S0", "--workers", "1"]) == 0
    assert cli.main(["stage", "igr", "--root", str(root), "--stage", "S0"]) == 0
    assert cli.main(["stage", "tapes", "--root", str(root), "--stage", "S0"]) == 0
    domains.run(root, "S0", V, "test", "1", work.lib, "gathering", cpus=1)

    # S1: a second stage with mutated contigs at depth called (no domains)
    import gzip

    fa = base / "s1.fna"
    s0_fa = Path(json.loads(work.objects("S0").read_text())[0]["path"])
    raw_text = gzip.open(s0_fa, "rt").read() if str(s0_fa).endswith(".gz") else s0_fa.read_text()
    lines = raw_text.splitlines()
    s1_lines = []
    for ln in lines:
        if ln.startswith(">"):
            s1_lines.append(ln + "_s1")
        else:
            chars = list(ln)
            for idx in range(0, len(chars), 30):
                if chars[idx] == "A":
                    chars[idx] = "T"
                elif chars[idx] == "T":
                    chars[idx] = "A"
            s1_lines.append("".join(chars))
    fa.write_text("\n".join(s1_lines) + "\n")
    # Parse mutated sequences to compute contig_uids for GenomeInput
    s1_seqs: list[str] = []
    cur_seq: list[str] = []
    for ln in s1_lines:
        if ln.startswith(">"):
            if cur_seq:
                s1_seqs.append("".join(cur_seq))
                cur_seq = []
        else:
            cur_seq.append(ln.strip())
    if cur_seq:
        s1_seqs.append("".join(cur_seq))
    s1_cuids = [ids.contig_uid(seq.encode()).hex() for seq in s1_seqs]
    obj = base / "s1_objects.json"
    obj.write_text(json.dumps([{"path": str(fa), "source_name": "s1", "source_release": "1"}]))
    gfile = base / "s1_genomes.json"
    gfile.write_text(json.dumps([{"contig_uids": s1_cuids, "kind": "prokaryote"}]))
    assert cli.main(["stage", "inputs", "--root", str(root), "--stage", "S1", "--genomes", str(gfile)]) == 0
    assert cli.main(["stage", "ingest", "--root", str(root), "--stage", "S1", "--objects", str(obj)]) == 0
    assert cli.main(["stage", "genecall", "--root", str(root), "--stage", "S1", "--workers", "1"]) == 0
    assert cli.main(["stage", "igr", "--root", str(root), "--stage", "S1"]) == 0
    assert cli.main(["stage", "tapes", "--root", str(root), "--stage", "S1"]) == 0
    return root, work.lib


def test_c38g_tokens_postings_byte_identical_and_nbr_round_source(direct, tmp_path):
    """§3.8 test (6): tokens & postings byte-identical with/without rounds; nbr summarize shows source='round'."""
    root_plain, _ = _build_two_stage_root(tmp_path / "plain", direct)
    root_round = tmp_path / "with_rounds" / "root"
    shutil.copytree(root_plain, root_round)

    s1_prots = sorted({
        u for f in sorted((root_round / f"l0seq/protein_seq/gene_call_version={V}/stage=S1").rglob("*.parquet"))
        for u in pq.read_table(f, columns=["protein_uid"]).column(0).to_pylist()
    })
    assert len(s1_prots) >= 5
    hit_u = s1_prots[:2]
    none_u = s1_prots[2:4]
    unscan_u = s1_prots[4]
    covered = hit_u + none_u
    dspec, cspec = TABLES["round_domains"], TABLES["round_coverage"]
    rid = "r-00000000000000c1"
    dom_rows = [
        {
            "protein_uid": u, "model_acc": "PF99901.1", "model_name": "RoundPlanted", "model_len": 100,
            "seq_bitscore": 60.0, "seq_evalue_zref": 1e-12, "dom_bitscore": 55.0, "dom_ievalue_zref": 1e-11,
            "env_start": 5, "env_end": 55, "ali_start": 5, "ali_end": 55, "hmm_from": 0, "hmm_to": 50,
            "cutoff": "gathering", "z": 1e8, "dom_z": 1e8,
        }
        for u in sorted(hit_u)
    ]
    bio.write_table(pa.Table.from_pylist(dom_rows, schema=dspec.schema), dspec,
                    root_round / dspec.location / f"round={rid}" / "library=test-1" / "part-0.parquet")
    bio.write_table(
        pa.Table.from_pylist([{"protein_uid": u, "n_rows": 1 if u in hit_u else 0} for u in sorted(covered)],
                             schema=cspec.schema),
        cspec,
        root_round / cspec.location / f"round={rid}" / "library=test-1" / "part-0.parquet",
    )

    rel = "rC38g.1"
    for r in (root_plain, root_round):
        catalog.define_release(r, rel, ["S0", "S1"])
        l1rounds.merge_release(r, rel)
        tokens.run(r, rel, ["S0", "S1"], without_clusters=True)
        postings.run(r, rel, ["S0", "S1"])
        dim.run(r, rel, ["S0", "S1"])

    # §3.8 test (6): token_dict, tape_tokens, token_postings are BYTE-IDENTICAL with and without rounds!
    for sub in (f"l3/token_dict/release={rel}", f"l3/tape_tokens", f"l3/token_postings/release={rel}"):
        f_plain = sorted((root_plain / sub).rglob("*.parquet"))
        f_round = sorted((root_round / sub).rglob("*.parquet"))
        assert [p.relative_to(root_plain) for p in f_plain] == [p.relative_to(root_round) for p in f_round]
        for a, b in zip(f_plain, f_round):
            assert a.read_bytes() == b.read_bytes(), (a, b)

    # Manifest, verify_manifest, provenance all pass on root_round
    man = catalog.build_manifest(root_round, rel, ["S0", "S1"])
    assert catalog.verify_manifest(man) == []
    assert provenance.check_provenance(man) == []

    # Check hits_tsv reflects round measured ('RoundPlanted'/'measured'), round none ('-'/'none'), unscanned ('?'/'not scanned')
    tv = tapes.tape_version([V])
    fake_hits = [
        {"target_kind": "protein", "target_uid": u, "seq_bitscore": 50.0, "seq_evalue_zref": 1e-10, "bitscore": 50.0}
        for u in (hit_u[0], none_u[0], unscan_u)
    ]
    tsv_lines = list(export.hits_tsv(root_round, fake_hits, release=rel, tape_version=tv, gene_call_version=V,
                                     manifest=man, library="test-1"))
    hdr = tsv_lines[0].split("\t")
    pu_i, arch_i, ds_i = hdr.index("protein_uid"), hdr.index("pfam_architecture"), hdr.index("domain_state")
    got_hits = {ln.split("\t")[pu_i]: (ln.split("\t")[arch_i], ln.split("\t")[ds_i]) for ln in tsv_lines[1:]}
    assert got_hits[hit_u[0].hex()] == ("RoundPlanted", "measured")
    assert got_hits[none_u[0].hex()] == ("-", "none")
    assert got_hits[unscan_u.hex()] == ("?", "not scanned")

    # Check nbr.summarize on an anchor in S1 near hit_u[0] -> sees RoundPlanted with source='round'
    occ_files = bio.part_files(root_round / f"l1/protein_occurrence/gene_call_version={V}/stage=S1")
    occ_rows = [r for f in occ_files for r in pq.read_table(f).to_pylist()]
    target_contig = next(r["contig_uid"] for r in occ_rows if r["protein_uid"] == hit_u[0])
    anchor_cand = next((r["protein_uid"] for r in occ_rows
                        if r["contig_uid"] == target_contig and r["protein_uid"] != hit_u[0]), None)
    assert anchor_cand is not None
    summ = nbr.summarize(root_round, {anchor_cand}, nbr.Policy.parse("20g"), by=("pfam",),
                         release=rel, manifest=man)
    round_rows = [r for r in summ["rows"] if r.get("model_acc") == "PF99901.1"]
    assert round_rows and all(r.get("source") == "round" for r in round_rows), summ["rows"]


def test_c38g_micro_s0_without_rounds_keeps_exact_manifest_hash(tmp_path):
    """Accept: fixture manifest hash of a release without rounds is unchanged (5fa9e49761393e3468ebc3c1d7d7224a)."""
    if shutil.which("mmseqs") is None:
        pytest.skip("build_fixture needs mmseqs")
    out = tmp_path / "fixture"
    p = subprocess.run([sys.executable, str(REPO / "scripts/stress/build_fixture.py"),
                        str(REPO), str(out), "rX.fixture"],
                       capture_output=True, text=True, check=True)
    report = json.loads(p.stdout)
    assert report["manifest_hash"] == c38a.FIXTURE_HASH, report
    man = catalog.build_manifest(out, "rX.fixture", ["S0"])
    assert catalog.manifest_hash(man) == c38a.FIXTURE_HASH


def test_c38g_build_stage_runs_rounds_step_between_cluster_and_tokens(direct, tmp_path):
    """stage.build_stage runs `rounds` after `cluster` and before `tokens` when `release_round` is non-empty."""
    if shutil.which("mmseqs") is None:
        pytest.skip("build-stage clusters with MMseqs2")
    from biobox import stage as bs

    prof, w = direct["prof"], direct["work"]
    root = tmp_path / "root"
    # Pre-seed one complete round on disk before build_stage cuts the release
    u0 = ids.protein_uid("M" + "A" * 40)
    rid = "r-00000000000000aa"
    bio.write_table(_dom_table([u0]), TABLES["round_domains"],
                    root / TABLES["round_domains"].location / f"round={rid}" / "library=test-1" / "part-0.parquet")
    bio.write_table(_cov_table([u0], [1]), TABLES["round_coverage"],
                    root / TABLES["round_coverage"].location / f"round={rid}" / "library=test-1" / "part-0.parquet")

    seen: list[str] = []
    b = bs.Build(
        root=root, stage="S0", release="rBS.1", stages=("S0",),
        genomes=w.genomes("S0"), objects=w.objects("S0"), depth="called", threads=1,
    )
    rc = bs.build_stage(prof, b, run=c38a.runner, _hook=seen.append)
    assert rc == 0
    assert seen.index("after:cluster") < seen.index("after:rounds") < seen.index("after:tokens"), seen
    assert (root / "runs/stage=rBS.1/rounds.json").is_file()
    assert bio.part_files(root / "l1/round_domains/release=rBS.1/library=test-1")
    assert bio.part_files(root / "l1/round_coverage/release=rBS.1/library=test-1")
    man = catalog.read_manifest(root / "manifests/release=rBS.1.json")
    assert catalog.verify_manifest(man) == []
    assert provenance.check_provenance(man) == []


def _synth_store(out: Path, copies: int) -> Path:
    argv = sys.argv
    sys.argv = [str(SCRIPT), str(REPO), str(out), str(copies), "--shard-bases", str(2**26), "--seed", "1"]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()):
            g = runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
    g["con"].close()
    return Path(g["root"])


def _make_five_rounds(root: Path) -> int:
    """Split S0's synthetic proteins and domains into 5 overlapping rounds."""
    fp = bio.part_files(root / f"l1/feature_protein/gene_call_version={V}/stage=S0")
    uids = sorted({u for f in fp for u in pq.read_table(f, columns=["protein_uid"], partitioning=None)
                   .column(0).to_pylist()})
    dom_files = bio.part_files(root / "l1/domains/library=synfam-1/stage=S0")
    dspec, cspec = TABLES["round_domains"], TABLES["round_coverage"]
    all_dom = pa.concat_tables([pq.read_table(f, columns=dspec.schema.names, partitioning=None) for f in dom_files])
    n = len(uids)
    chunk = max(1, n // 5)
    for r in range(5):
        lo = max(0, r * chunk - chunk // 4)
        hi = n if r == 4 else min(n, (r + 1) * chunk + chunk // 4)
        sub_uids = uids[lo:hi]
        sub_set = set(sub_uids)
        rid = f"r-{r:016x}"
        mask = pa.array([u in sub_set for u in all_dom.column("protein_uid").to_pylist()])
        sub_dom = all_dom.filter(mask)
        bio.write_table(sub_dom, dspec,
                        root / dspec.location / f"round={rid}" / "library=synfam-1" / "part-0.parquet")
        hit_counts: dict[bytes, int] = {}
        for u in sub_dom.column("protein_uid").to_pylist():
            hit_counts[u] = hit_counts.get(u, 0) + 1
        bio.write_table(
            pa.table({"protein_uid": pa.array(sub_uids, pa.binary(16)),
                      "n_rows": pa.array([hit_counts.get(u, 0) for u in sub_uids], pa.uint32())},
                     schema=cspec.schema),
            cspec,
            root / cspec.location / f"round={rid}" / "library=synfam-1" / "part-0.parquet",
        )
    catalog.define_release(root, "rMem.1", ["S0"])
    return n


def test_c38g_merge_release_100x_memory_bounded(tmp_path, record_property):
    """AGENTS §4.2: streaming k-way merge of 5 overlapping rounds stays bounded on 10x -> 100x synthetic store."""
    got = {}
    for n in (SMALL, BIG):
        base = tmp_path / f"x{n}"
        root = _synth_store(base, n)
        n_prots = _make_five_rounds(root)
        peak, out, scratch = tmp_path / f"peak-{n}.json", tmp_path / f"out-{n}", tmp_path / f"scratch-{n}"
        scratch.mkdir()
        cmd = runstep.measured([sys.executable, "-c", CHILD_MERGE, str(root), "rMem.1", str(out), str(scratch)],
                               peak, False)
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        assert p.returncode == 0, p.stderr[-4000:]
        rss, source = runstep.step_peak(peak)
        assert rss is not None, source
        res = json.loads(p.stdout.strip().splitlines()[-1])
        assert res["rounds"] == 5 and res["libraries"]["synfam-1"]["proteins"] == n_prots
        got[n] = (rss, source, res)
    record_property("c38g_peak_bytes", {n: g[0] for n, g in got.items()})
    growth = got[BIG][0] - got[SMALL][0]
    print(f"C-38g merge_release peaks: 10x={got[SMALL][0] / 2**20:.1f} MiB, "
          f"100x={got[BIG][0] / 2**20:.1f} MiB, growth={growth / 2**20:.1f} MiB")
    assert growth <= GROWTH_BOUND, f"grew {growth / 2**20:.1f} MiB (> {GROWTH_BOUND / 2**20:.0f} MiB)"
    assert max(g[0] for g in got.values()) <= PEAK_BOUND, got

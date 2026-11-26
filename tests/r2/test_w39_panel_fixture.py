"""Ring R2: end-to-end fixture, external hmmsearch oracle, --carry-from, sql --panel, and 100x RSS tests (card W-39)."""

from __future__ import annotations

import hashlib
import io
import json
import os
import random
import resource
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, cli, dim, fasta, ids
from biobox import io as bio
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import write_genome_input
from biobox.l1 import domains, genecall, igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import GenomeInput
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import table
from biobox.search import engine, pack, panel, results

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-39", "S-03", "S-04", "MR-12", "W-28", "W-29")]

# Deterministic synthetic protein templates (distinct, non-homologous amino acid sequences)
def _rand_aa(seed: int, length: int) -> str:
    rng = random.Random(seed)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    return "M" + "".join(rng.choice(alphabet) for _ in range(length - 1))


# Codon table (standard code 11, no internal stops) to back-translate proteins into valid bacterial CDSs
_CODON = {
    "A": "GCT", "C": "TGC", "D": "GAT", "E": "GAA", "F": "TTC", "G": "GGT", "H": "CAT", "I": "ATC",
    "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAC", "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT",
    "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAC",
}


def _to_dna(aa: str) -> str:
    return "".join(_CODON[c] for c in aa) + "TAA"


def _hmm_bytes(name: str, seq: str, *, acc: str | None = None, ga: tuple[float, float] | None = None) -> bytes:
    hmm = engine.build_from_sequence(name, seq)
    if acc is not None:
        hmm.accession = acc.encode()
    if ga is not None:
        hmm.cutoffs.gathering = ga
    buf = io.BytesIO()
    hmm.write(buf)
    return buf.getvalue()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _build_stage(root: Path, stage: str, contigs: list[tuple[str, str, bool, str]], tmp: Path) -> None:
    """Ingest contigs (acc, dna, is_circular, kind), call genes, build IGRs and tapes."""
    fa = tmp / f"{stage}.fna"
    lines = []
    for acc, dna, _, _ in contigs:
        lines += [f">{acc}", dna]
    fa.write_text("\n".join(lines) + "\n", encoding="utf-8")
    ingest([SourceObject(fa, "w39syn", "1")], root, stage)
    genomes = []
    for acc, dna, circ, kind in contigs:
        cu = ids.contig_uid(dna.encode("ascii"))
        genomes.append(GenomeInput((cu,), 11, kind=kind, circular=(cu,) if circ else ()))
    write_genome_input(root, stage, genomes)
    genecall.run(root, stage)
    igr.run(root, stage)
    tapes.run(root, stage)


def _finalize_release(root: Path, release: str, stages: list[str], packs: Path) -> None:
    catalog.define_release(root, release, stages)
    dim.run(root, release, stages)
    tokens.run(root, release, stages, V, without_clusters=True)
    postings.run(root, release, stages, V)
    man = catalog.build_manifest(root, release, stages)
    catalog.write_manifest(man, catalog.manifest_path(root, release))
    catalog.generate_catalog(man, catalog.catalog_path(root, release))
    pack.build(root, release, packs=packs)


def _spacer_dna(seed: int, bp: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(bp))


def _make_contig_dna(proteins: list[str], *, gap_bp: int = 120, pad_bp: int = 300, seed: int = 1) -> str:
    """Assemble a synthetic DNA contig carrying CDSs with RBS and spacers so pyrodigal calls each cleanly."""
    parts = [_spacer_dna(seed, pad_bp)]
    for i, aa in enumerate(proteins):
        # Shine-Dalgarno RBS + 7 bp spacer before ATG
        parts.append("AGGAGGTAAAAAAA" + _to_dna(aa))
        parts.append(_spacer_dna(seed + 100 + i, gap_bp))
    parts.append(_spacer_dna(seed + 999, pad_bp))
    return "".join(parts)


@pytest.fixture(scope="module")
def w39_world(tmp_path_factory):
    """Two-stage synthetic store (S0, S1) with planted multi-subunit, fusion, decoy-vetoed, and 0-hit proteins."""
    base = tmp_path_factory.mktemp("w39_world")
    root = base / "root"
    packs = base / "packs"
    runs = base / "runs"
    work = base / "work"
    hmms_dir = base / "hmms"
    hmms_dir.mkdir()

    # Core protein motifs (150 aa each)
    seq_m1 = _rand_aa(101, 150)
    seq_m2 = _rand_aa(202, 150)
    seq_m3 = _rand_aa(303, 150)
    seq_decoy = _rand_aa(404, 150)
    seq_bg1 = _rand_aa(505, 160)
    seq_bg2 = _rand_aa(606, 160)
    seq_bg3 = _rand_aa(707, 160)
    seq_neg = _rand_aa(808, 160)

    # Fusion protein M1 + linker + M2
    seq_fusion = seq_m1 + _rand_aa(901, 20) + seq_m2[1:]
    # Contig 1 (S0, linear): M1 + M2 co-occurring operon + background
    c1_dna = _make_contig_dna([seq_bg1, seq_m1, seq_m2, seq_bg2], seed=11)
    # Contig 2 (S0, circular): fusion M1+M2 + negative control + decoy protein
    c2_dna = _make_contig_dna([seq_fusion, seq_neg, seq_decoy, seq_bg3], seed=22)
    # Contig 3 (S1, linear): M1 + M3 operon + 0-hit background
    seq_bg4 = _rand_aa(909, 170)
    c3_dna = _make_contig_dna([seq_m1, seq_m3, seq_bg4], seed=33)

    _build_stage(root, "S0", [("C1_OPERON", c1_dna, False, "prokaryote"), ("C2_CIRC", c2_dna, True, "plasmid")], base)
    _finalize_release(root, "rW39.1", ["S0"], packs)

    _build_stage(root, "S1", [("C3_S1", c3_dna, False, "virus")], base)
    _finalize_release(root, "rW39.2", ["S0", "S1"], packs)

    # Write HMM files
    hmm_m12 = hmms_dir / "panel_m12.hmm"
    hmm_m12.write_bytes(
        _hmm_bytes("Model1", seq_m1, acc="M1", ga=(30.0, 25.0))
        + _hmm_bytes("Model2", seq_m2, acc="M2", ga=(30.0, 25.0))
    )
    hmm_m3 = hmms_dir / "panel_m3.hmm"
    hmm_m3.write_bytes(_hmm_bytes("Model3", seq_m3, acc="M3", ga=(30.0, 25.0)))
    # Decoy HMM matching seq_decoy strongly, plus a weak cross-hit target HMM matching seq_decoy weakly
    hmm_decoy = hmms_dir / "decoy.hmm"
    hmm_decoy.write_bytes(_hmm_bytes("DecoyRdRP", seq_decoy, acc="DECOY1", ga=(30.0, 25.0)))

    # Resolve protein_uid for M1-only positive control and seq_neg negative control from S0 protein_seq
    prots_s0 = domains.load_proteins(root, V, "S0")
    pos_m1_uid = next(u for u, s in prots_s0.items() if seq_m1[5:-5] in s and seq_m2[5:-5] not in s)
    neg_uid = next(u for u, s in prots_s0.items() if seq_neg[5:-5] in s)
    occ_s0 = pq.read_table(next((root / f"l1/protein_occurrence/gene_call_version={V}/stage=S0").glob("*.parquet"))).to_pylist()
    pos_m1_occ = next(r for r in occ_s0 if r["protein_uid"] == pos_m1_uid)
    pos_m1_span = f"C1_OPERON:{pos_m1_occ['start']}-{pos_m1_occ['end']}"

    return {
        "base": base, "root": root, "packs": packs, "runs": runs, "work": work,
        "hmm_m12": hmm_m12, "hmm_m3": hmm_m3, "hmm_decoy": hmm_decoy,
        "pos_m1_uid": pos_m1_uid.hex(), "pos_m1_span": pos_m1_span, "neg_uid": neg_uid.hex(),
    }


def _external_hmmsearch_oracle(hmm_files: list[Path], faa: Path, out_domtbl: Path) -> set[tuple[str, str, int, int]]:
    """Run external HMMER `hmmsearch` binary at Z=1e8, domZ=1e8 with --cut_ga."""
    hmmsearch = shutil.which("hmmsearch")
    if not hmmsearch:
        pytest.skip("hmmsearch binary not on PATH")
    combined = out_domtbl.parent / "oracle_combined.hmm"
    combined.write_bytes(b"".join(p.read_bytes() for p in hmm_files))
    subprocess.run(
        [hmmsearch, "--cpu", "1", "-Z", "100000000", "--domZ", "100000000", "--cut_ga",
         "--domtblout", str(out_domtbl), "-o", os.devnull, str(combined), str(faa)],
        check=True,
    )
    hits = set()
    for line in out_domtbl.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split()
        # domtblout columns: 0=target_name, 3=query_name, 4=query_acc, 19=env_from(1-based), 20=env_to
        tname, qname, qacc = f[0], f[3], f[4]
        qkey = qacc if qacc != "-" else qname
        hits.add((tname, qkey, int(f[19]) - 1, int(f[20])))
    return hits


def test_w39_end_to_end_oracle_hmmsearch_and_results_roundtrip(w39_world, tmp_path):
    """End-to-end panel run matches external hmmsearch 1:1, produces all tables/diagnostics, and archives cleanly."""
    out = tmp_path / "panel_out"
    rc = cli.main([
        "search", "panel",
        "--root", str(w39_world["root"]),
        "--release", "rW39.1",
        "--pack-dir", str(w39_world["packs"]),
        "--runs", str(w39_world["runs"]),
        "--work", str(w39_world["work"]),
        "--hmm", str(w39_world["hmm_m12"]),
        "--decoy-hmm", f"decoy:{w39_world['hmm_decoy']}",
        "--cutoff", "gathering",
        "--positive", f"M1={w39_world['pos_m1_span']}",
        "--negative-control", w39_world["neg_uid"],
        "--negative", "shuffled_seed",
        "--min-models", "2",
        "--max-gap-genes", "5",
        "--max-gap-bp", "10000",
        "--neighbors", "5",
        "--cpus", "1",
        "--jobs", "1",
        "--out", str(out),
    ])
    assert rc == 0

    for fname in ("models.parquet", "models.tsv", "protein_hits.parquet", "orf_hits.parquet", "orf_hits.tsv",
                  "loci_summary.parquet", "loci_summary.tsv", "controls.tsv", "coverage.json",
                  "unannotated_proteins.txt", "RUN"):
        assert (out / fname).is_file(), f"missing {fname}"
    assert (out / "loci").is_dir() and len(list((out / "loci").glob("locus_*.tsv"))) >= 1

    # Verify external hmmsearch oracle 1:1 concordance on target models
    faa = tmp_path / "all_s0.faa"
    prots = domains.load_proteins(w39_world["root"], V, "S0")
    faa.write_text("".join(f">{u.hex()}\n{s}\n" for u, s in sorted(prots.items())), encoding="utf-8")
    oracle_hits = _external_hmmsearch_oracle([w39_world["hmm_m12"]], faa, tmp_path / "oracle.domtbl")

    ph = pq.read_table(out / "protein_hits.parquet").to_pylist()
    target_hits = {(r["protein_uid"].hex(), r["model_key"], r["env_start"], r["env_end"])
                   for r in ph if not r["is_decoy"]}
    assert target_hits == oracle_hits

    # Verify loci_summary recovered both C1_OPERON (2-ORF M1+M2) and C2_CIRC (1-ORF fusion M1+M2)
    loci = pq.read_table(out / "loci_summary.parquet").to_pylist()
    assert len(loci) == 2
    assert {r["contig_accession"] for r in loci} == {"C1_OPERON", "C2_CIRC"}

    # Verify results archive + flush + get round-trip (W-28 / W-29)
    run_id, rel = results.read_run_file(out)
    assert rel == "rW39.1"
    rdir = w39_world["runs"] / f"run={run_id}"
    assert (rdir / "run.parquet").is_file()
    assert (rdir / "panel_hits.parquet").is_file()
    run_row = pq.read_table(rdir / "run.parquet").to_pylist()[0]
    assert run_row["status"] == "ok"
    assert "panel.v1" in run_row["recipe_yaml"]
    assert "command" in json.loads(run_row["limits_json"])

    bulk = tmp_path / "bulk"
    bulk.mkdir()
    where = results.Places(spool=tmp_path / "spool", area=bulk / "results", bulk=bulk, floor=64 * 1024 * 1024)
    line, state = results.archive(rdir, [out], where)
    assert state == "spooled"
    restored = tmp_path / "restored"
    _, n_files = results.get(where, run_id, restored)
    assert n_files > 0
    assert (restored / "outputs" / out.name / "loci_summary.parquet").is_file()


def test_w39_positive_control_miss_fails_run(w39_world, tmp_path):
    """A missed positive control exits non-zero and records status='failed'."""
    out = tmp_path / "fail_out"
    rc = cli.main([
        "search", "panel",
        "--root", str(w39_world["root"]),
        "--release", "rW39.1",
        "--pack-dir", str(w39_world["packs"]),
        "--runs", str(w39_world["runs"]),
        "--work", str(w39_world["work"]),
        "--hmm", str(w39_world["hmm_m12"]),
        "--cutoff", "gathering",
        "--positive", f"M2={w39_world['pos_m1_uid']}",  # M2 will NOT hit M1's protein
        "--negative", "none",
        "--out", str(out),
    ])
    assert rc == 1
    run_id, _ = results.read_run_file(out)
    run_row = pq.read_table(w39_world["runs"] / f"run={run_id}" / "run.parquet").to_pylist()[0]
    assert run_row["status"] == "failed"


def test_w39_byte_identity_jobs_and_carry_from(w39_world, tmp_path):
    """--jobs 1 vs --jobs 2 and --carry-from (new stage S0->S0,S1 and new HMM M1,M2->M1,M2,M3) are byte-identical."""
    out_j1 = tmp_path / "j1"
    out_j2 = tmp_path / "j2"
    common = [
        "search", "panel",
        "--root", str(w39_world["root"]),
        "--pack-dir", str(w39_world["packs"]),
        "--runs", str(w39_world["runs"]),
        "--work", str(w39_world["work"]),
        "--cutoff", "gathering",
        "--positive", f"M1={w39_world['pos_m1_uid']}",
        "--negative", "none",
        "--min-models", "1",
    ]
    assert cli.main([*common, "--release", "rW39.1", "--hmm", str(w39_world["hmm_m12"]),
                     "--jobs", "1", "--cpus", "1", "--out", str(out_j1)]) == 0
    assert cli.main([*common, "--release", "rW39.1", "--hmm", str(w39_world["hmm_m12"]),
                     "--jobs", "2", "--cpus", "1", "--out", str(out_j2)]) == 0

    for fname in ("models.parquet", "protein_hits.parquet", "orf_hits.parquet", "loci_summary.parquet"):
        assert _sha256(out_j1 / fname) == _sha256(out_j2 / fname), f"jobs 1 vs 2 mismatch on {fname}"

    # 1. Fresh full scan on rW39.2 (S0, S1) with M1, M2, M3
    out_full = tmp_path / "full_r2_m123"
    assert cli.main([*common, "--release", "rW39.2",
                     "--hmm", str(w39_world["hmm_m12"]), "--hmm", str(w39_world["hmm_m3"]),
                     "--jobs", "1", "--cpus", "1", "--out", str(out_full)]) == 0

    # 2. Carry-forward from out_j1 (rW39.1, S0, M1+M2) -> rW39.2 (S0+S1, M1+M2+M3): BOTH new stage AND new HMM!
    out_carry = tmp_path / "carry_r2_m123"
    assert cli.main([*common, "--release", "rW39.2",
                     "--hmm", str(w39_world["hmm_m12"]), "--hmm", str(w39_world["hmm_m3"]),
                     "--carry-from", str(out_j1),
                     "--jobs", "1", "--cpus", "1", "--out", str(out_carry)]) == 0

    for fname in ("models.parquet", "protein_hits.parquet", "orf_hits.parquet", "loci_summary.parquet"):
        assert _sha256(out_full / fname) == _sha256(out_carry / fname), f"--carry-from mismatch on {fname}"


def test_w39_sql_panel_mounting_and_row_group_pruning(w39_world, tmp_path, capsys):
    """`bio-box sql --panel rt=DIR` mounts views, warns on release mismatch, refuses collisions, and prunes row groups."""
    out = tmp_path / "panel_sql"
    assert cli.main([
        "search", "panel",
        "--root", str(w39_world["root"]),
        "--release", "rW39.1",
        "--pack-dir", str(w39_world["packs"]),
        "--runs", str(w39_world["runs"]),
        "--work", str(w39_world["work"]),
        "--hmm", str(w39_world["hmm_m12"]),
        "--cutoff", "gathering",
        "--positive", f"M1={w39_world['pos_m1_uid']}",
        "--negative", "none",
        "--min-models", "1",
        "--out", str(out),
    ]) == 0

    # Query via bio-box sql --panel rt=DIR
    q_out = tmp_path / "q.tsv"
    rc = cli.main([
        "sql", "--root", str(w39_world["root"]), "--release", "rW39.1",
        "--panel", f"rt={out}",
        "-e", "SELECT m.model_key, count(*) AS n FROM rt_models m JOIN rt_orf_hits h USING (model_key) GROUP BY m.model_key ORDER BY m.model_key",
        "--out", str(q_out),
    ])
    assert rc == 0
    lines = q_out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "model_key\tn"
    assert len(lines) == 3

    # Release mismatch prints 1-line stderr note
    capsys.readouterr()
    rc2 = cli.main([
        "sql", "--root", str(w39_world["root"]), "--release", "rW39.2",
        "--panel", f"rt={out}",
        "-e", "SELECT count(*) FROM rt_loci",
    ])
    assert rc2 == 0
    err = capsys.readouterr().err
    assert "rW39.1" in err and "rW39.2" in err

    # Colliding panel name refused with exit 2
    rc3 = cli.main([
        "sql", "--root", str(w39_world["root"]), "--release", "rW39.1",
        "--panel", f"l1_domains={out}",
        "-e", "SELECT 1",
    ])
    assert rc3 == 2

    # Multi-row-group pruning check via catalog.profile_query
    rg_dir = tmp_path / "panel_rg"
    shutil.copytree(out, rg_dir)
    # Write a 5-row-group sorted protein_hits.parquet with distinct protein_uids
    ph_schema = pq.read_schema(rg_dir / "protein_hits.parquet")
    base_row = pq.read_table(rg_dir / "protein_hits.parquet").to_pylist()[0]
    with pq.ParquetWriter(rg_dir / "protein_hits.parquet", ph_schema, compression="zstd") as w:
        for g in range(5):
            rows = [{**base_row, "protein_uid": bytes([g * 10 + i]) + b"\x00" * 15} for i in range(4)]
            w.write_table(pa.Table.from_pylist(rows, schema=ph_schema))
    con, scratch = cli.sql_connection(w39_world["root"], "rW39.1", panels=[("rt", rg_dir)])
    try:
        target_hex = (bytes([20]) + b"\x00" * 15).hex()
        prof = catalog.profile_query(
            con, f"SELECT count(*) FROM rt_protein_hits WHERE protein_uid = unhex('{target_hex}')",
            profile_dir=scratch / "duckdb" / "profile",
        )
        assert prof["result"][0][0] == 1
        ph_scans = [s for s in prof["scans"] if any("protein_hits.parquet" in f for f in s["files"])]
        assert ph_scans and ph_scans[0]["row_groups_scanned"] == 1
        assert ph_scans[0]["row_groups_total"] == 5
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)


def test_w39_100x_synthetic_memory_bound(w39_world, tmp_path):
    """Streaming panel search over a 100x scaled protein stream stays well under 2 GB peak RSS."""
    # Replicate S0 proteins 100x into a synthetic fastz shard pack and stream-search it in a child process
    from biobox.search import fastz

    prots = domains.load_proteins(w39_world["root"], V, "S0")
    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()
    records = []
    for copy_idx in range(100):
        for u, s in sorted(prots.items()):
            uid = hashlib.blake2b(u + copy_idx.to_bytes(4, "big"), digest_size=16).digest()
            records.append((uid.hex(), s))
    records.sort(key=lambda x: x[0])
    written = fastz.write_records(shard_dir, ".faa.zst", records, 2**30)
    shard_path = written[0]["path"]

    # The child's own peak is VmHWM of /proc/self/status, the high-water mark of the memory it has had since its exec.
    # Not `ru_maxrss`: on Linux a child's `ru_maxrss` carries the peak of the parent it was started from (exec keeps
    # the larger of the old and new images' peaks), so under one pytest process over the whole suite (`uv run pytest`,
    # GitHub Actions) it reports the suite's peak, 3,982,860,288 B on run 45711591139. Card R-01 (measured, wilkcomp02,
    # 2026-11-25, Python 3.14.10 in the test image): a parent wrote 3 GB and freed it, then `subprocess.run` started a
    # child that wrote 64 MiB; the child read ru_maxrss 3,011,248,128 B and VmHWM 78,209,024 B. This scan's VmHWM was
    # 346,992,640 B under a parent that had peaked at 3 GB (its ru_maxrss 3,036,643,328 B). Card B-24's test measures
    # its child the same way (tests/r2/test_micro_s0_cluster.py). `ru_maxrss` is used only where there is no /proc
    # (macOS), and is printed beside VmHWM.
    child = (
        "import json, os, resource, sys\n"
        "from pathlib import Path\n"
        "from biobox.search import panel\n"
        "hmms = panel.load_panel_specs([sys.argv[1]], [], [], [], cutoff='gathering')\n"
        "hits = panel.scan_fastz_shards([Path(sys.argv[2])], hmms, cpus=1, block_residues=65536)\n"
        "rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss\n"
        "maxrss = rss if sys.platform == 'darwin' else rss * 1024\n"
        "st = '/proc/self/status'\n"
        "hwm = (next(int(x.split()[1]) * 1024 for x in open(st) if x.startswith('VmHWM:'))\n"
        "       if os.path.exists(st) else None)\n"
        "print(json.dumps({'n_hits': len(hits), 'hwm_bytes': hwm, 'maxrss_bytes': maxrss}))\n"
    )
    p = subprocess.run(
        [sys.executable, "-c", child, str(w39_world["hmm_m12"]), str(shard_path)],
        capture_output=True, text=True, check=True, timeout=180,
    )
    info = json.loads(p.stdout.strip().splitlines()[-1])
    peak = info["hwm_bytes"] if info["hwm_bytes"] is not None else info["maxrss_bytes"]
    print(f"\nw39 100x scan in a child: {info['n_hits']} hits, VmHWM {info['hwm_bytes']} B, "
          f"ru_maxrss {info['maxrss_bytes']} B (with its parent's peak)")
    assert info["n_hits"] >= 100
    assert peak < 2 * 1024**3, (f"peak RSS {peak} >= 2 GiB (VmHWM {info['hwm_bytes']}, "
                                f"ru_maxrss {info['maxrss_bytes']})")

"""Ring R2: card C-38b — Domain state: "none" is not "not scanned" (`planning/design/C-38_import_depths.md` §3.4).

Covers:
  (1) A mixed fixture release over `S0` (`annotated` for `test-1`, carrying `measured`, `inferred`, mixed
      `measured`+`inferred`, and `none` proteins) and `S1` (`called`, carrying T4 `not scanned` proteins plus a
      mutated lambda contig sharing exact `protein_uid`s with `S0`'s `measured`, `inferred`, and `none` proteins
      alongside an `S1`-only mutated protein `Q'`): every protein's state in the catalog view `domain_state` and
      scalar macro `domain_state(protein_uid, library)` (under both `open_catalog` and plain `duckdb.connect`,
      and on a pre-C-38b catalog) equals an independent brute-force reader over the Parquet files written
      separately in this test (AGENTS §5.2).
  (2) `catalog.profile_query` proves O(1) row-group reads for single-protein lookups across multi-row-group
      `protein_seq` and `l1_domains` files, both for literal keys on `domain_state` view / `domain_state(p, l)`
      and for joined-table keys on `domain_state(q.protein_uid, 'test-1')` (AGENTS §4.4).
  (3) Exports on the mixed release (`hits_tsv`, `locus_tables` including a non-CDS slot, `nbr.summarize`,
      `flanks.run`, and `search hmm` `README.txt`) print `-` for `none` and `?` for `not scanned` in
      `pfam_architecture` and populate `domain_state` against the brute-force Parquet oracle.
  (4) Exports on an all-annotated release (`S0` and `S1` both annotated) have zero `?` or `not scanned` cells,
      `hits.faa`, `neighbors.tsv` and `nbr flanks` are byte-identical to pre-C-38b, and `hits.tsv` and `loci/*.tsv`
      are byte-identical once the single `domain_state` column is removed.
"""

from __future__ import annotations

import dataclasses
import json
import re
import shutil
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from biobox import catalog, cli, dim, fasta, ids
from biobox.l0seq import SourceObject, ingest, revcomp
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import domains, igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.l2 import cluster as l2_cluster
from biobox.l3 import postings, tapes, tokens
from biobox.registry import schemas
from biobox.search import export, nbr, pack

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_search as ms  # noqa: E402

world = ms.world
FIX = ms.REPO / "tests" / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
ALA = {b"GCT", b"GCC", b"GCA", b"GCG"}

pytestmark = [
    pytest.mark.r2,
    pytest.mark.covers("C-38b", "LY-06", "LY-09", "LY-10", "W-05", "W-06", "W-07"),
]


def _mutate_middle_codon(seq: bytes, start: int, end: int, strand: int) -> bytes:
    """Change one interior codon so the contig gets a new contig_uid and this one gene gets a new protein_uid,
    leaving every other gene on the contig with its exact original protein_uid."""
    k = (end - start) // 3 // 2
    pos = start + 3 * k if strand == 1 else end - 3 * (k + 1)
    codon = seq[pos:pos + 3] if strand == 1 else revcomp(seq[pos:pos + 3])
    new = b"TGG" if codon.upper() in ALA else b"GCT"
    return seq[:pos] + (new if strand == 1 else revcomp(new)) + seq[pos + 3:]


def _standin_clusters(root: Path, rel: str, stages: list[str]) -> None:
    """Deterministic stand-in NR90/NR50 clusters and `member_cluster` lookup over `stages`."""
    uids = sorted({
        u
        for s in stages
        for f in sorted((root / f"l0seq/protein_seq/gene_call_version={V}/stage={s}").rglob("*.parquet"))
        for u in pq.read_table(f, columns=["protein_uid"], partitioning=None).column(0).to_pylist()
    })
    rows = [
        {"member_uid": p, "level": lvl, "cluster_uid": bytes([p[0] % g, len(lvl)]) * 8, "rep_uid": p,
         "cluster_params_hash": "stand-in"}
        for p in uids for lvl, g in (("nr90", 97), ("nr50", 23))
    ]
    rows.sort(key=lambda r: (r["level"], r["cluster_uid"], r["member_uid"]))
    d = root / f"l2/cluster_member/release={rel}"
    d.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=schemas.table("cluster_member").schema),
                   d / "part-0.parquet", compression="zstd", write_statistics=True)
    l2_cluster.write_member_cluster([d / "part-0.parquet"], root / f"l2/member_cluster/release={rel}")


def _rewrite_small_row_groups(path: Path, schema: pa.Schema, row_group_size: int = 16) -> None:
    t = pq.read_table(path, partitioning=None)
    if t.schema != schema:
        t = t.cast(schema)
    pq.write_table(t, path, row_group_size=row_group_size, compression="zstd", write_statistics=True)


@pytest.fixture()
def l1_domains_with_inferred(monkeypatch):
    """Enable the nullable `inferred_from` column on `TABLES['l1_domains']` (card B-47c / LY-09 compatibility)."""
    old = schemas.TABLES["l1_domains"]
    if "inferred_from" not in old.schema.names:
        ext = dataclasses.replace(old, schema=old.schema.append(pa.field("inferred_from", pa.binary(16), nullable=True)))
        monkeypatch.setitem(schemas.TABLES, "l1_domains", ext)
    return schemas.TABLES["l1_domains"]


def _build_mixed_fixture(base: Path, dom_spec) -> dict:
    """Build a two-stage store under `base / 'root'`:
      - `S0` (`annotated test-1`): phage lambda (`NC_001416.1`) + *M. genitalium* (`NC_000908.2`), with
        `measured`, `inferred`, mixed `measured`+`inferred`, and `none` proteins in `l1_domains`.
      - `S1` (`called`): phage T4 (`NC_000866.4`) + `lambda_q_mut` (a mutated copy of lambda whose gene `Q`
        became an `S1`-only protein `Q'` while 60+ other lambda genes share exact `protein_uid`s with `S0`'s
        `measured`, `inferred`, and `none` proteins).
    Also rewrites `protein_seq` and `l1_domains` into 16-row row groups so `profile_query` tests real pruning.
    """
    root = base / "root"
    rel = "r26.11.2"
    objs = CONTRACT["objects"]
    s0_objs = [o for o in objs if o["accession"] != "NC_000866.4"]
    t4_obj = next(o for o in objs if o["accession"] == "NC_000866.4")

    ingest([SourceObject(FIX / o["file"], "micro", "1", "PD", True) for o in s0_objs], root, "S0", shard_bases=200_000)
    s0_genomes, s0_circ = [], set()
    lam_seq = b""
    lam_cu = b""
    for o in s0_objs:
        seq = fasta.read(FIX / o["file"])[0].seq
        cu = ids.contig_uid(seq)
        if o["accession"] == "NC_001416.1":
            lam_seq, lam_cu = seq, cu
        circ = (cu,) if o["topology"] == "circular" else ()
        s0_circ |= set(circ)
        s0_genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=circ))
    write_genome_input(root, "S0", s0_genomes)
    call_genes(root, "S0")
    igr.run(root, "S0", s0_circ)
    tapes.run(root, "S0", s0_circ)

    # Pick lambda genes in S0 so that measured, mixed measured+inferred, inferred, and none proteins all sit
    # on lambda (and therefore will also appear on S1's `lambda_q_mut`!), and mutate one interior gene Q -> Q'.
    occ_s0 = [
        r for f in sorted((root / f"l1/protein_occurrence/gene_call_version={V}/stage=S0").glob("*.parquet"))
        for r in pq.read_table(f, partitioning=None).to_pylist()
    ]
    prots_s0 = domains.load_proteins(root, V, "S0")
    lam_genes = sorted([r for r in occ_s0 if r["contig_uid"] == lam_cu], key=lambda r: r["start"])
    assert len(lam_genes) >= 20
    # Mutate gene at index 10 on lambda_q_mut so genes 8, 9, 11, 12 share protein_uids with S0 and flank Q' in S1!
    q_gene = lam_genes[10]
    lam_mut_seq = _mutate_middle_codon(lam_seq, q_gene["start"], q_gene["end"], q_gene["strand"])
    lam_mut_cu = ids.contig_uid(lam_mut_seq)
    assert lam_mut_cu != lam_cu
    mut_fna = base / "lambda_q_mut.fna"
    mut_fna.write_bytes(b"\n".join(fasta.format_record("LAM_MUT.1", lam_mut_seq, 70)) + b"\n")

    # Ingest and call S1 (`called` depth: NO domain scan!)
    ingest([
        SourceObject(FIX / t4_obj["file"], "micro", "1", "PD", True),
        SourceObject(mut_fna, "micro", "1", "PD", True),
    ], root, "S1", shard_bases=200_000)
    t4_cu = ids.contig_uid(fasta.read(FIX / t4_obj["file"])[0].seq)
    s1_genomes = [
        GenomeInput((t4_cu,), t4_obj["genetic_code"], kind=t4_obj["kind"], circular=()),
        GenomeInput((lam_mut_cu,), 11, kind="virus", circular=()),
    ]
    write_genome_input(root, "S1", s1_genomes)
    call_genes(root, "S1")
    igr.run(root, "S1", set())
    tapes.run(root, "S1", set())

    occ_s1 = [
        r for f in sorted((root / f"l1/protein_occurrence/gene_call_version={V}/stage=S1").glob("*.parquet"))
        for r in pq.read_table(f, partitioning=None).to_pylist()
    ]
    s0_puids = {r["protein_uid"] for r in occ_s0}
    s1_puids = {r["protein_uid"] for r in occ_s1}
    shared_lam_genes = [g for g in lam_genes if g["protein_uid"] in s1_puids]
    assert len(shared_lam_genes) >= 15
    s1_mut_only = [r["protein_uid"] for r in occ_s1 if r["contig_uid"] == lam_mut_cu and r["protein_uid"] not in s0_puids]
    assert len(s1_mut_only) == 1
    q_prime_uid = s1_mut_only[0]

    # Build HMM library: 3 lambda shared proteins + 1 T4 protein (so we can also search T4 hits!)
    prots_s1 = domains.load_proteins(root, V, "S1")
    t4_seed_uid = max(
        (r["protein_uid"] for r in occ_s1 if r["contig_uid"] == t4_cu and r["protein_uid"] not in s0_puids),
        key=lambda u: (len(prots_s1[u]), u),
    )
    # Shared lambda targets around index 10 (indices 8, 9, 11 are right next to Q' on lambda_q_mut!)
    p_meas_shared = lam_genes[9]["protein_uid"]
    p_both_shared = lam_genes[11]["protein_uid"]
    p_inf_shared = lam_genes[8]["protein_uid"]
    p_none_shared = lam_genes[12]["protein_uid"]
    assert {p_meas_shared, p_both_shared, p_inf_shared, p_none_shared} <= s1_puids

    # Also pick S0-only (*M. genitalium*) proteins for measured, inferred, and none
    mg_puids = sorted(
        {r["protein_uid"] for r in occ_s0 if r["contig_uid"] != lam_cu and r["protein_uid"] not in s1_puids},
        key=lambda u: (-len(prots_s0[u]), u),
    )
    p_meas_s0, p_inf_s0, p_none_s0 = mg_puids[0], mg_puids[1], mg_puids[2]

    hmm_path = base / "c38b.hmm"
    lam_query_hmm = base / "lam_query.hmm"
    t4_query_hmm = base / "t4_query.hmm"
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder, bg = pyhmmer.plan7.Builder(alphabet, seed=42), pyhmmer.plan7.Background(alphabet)

    def build_hmm(uid: bytes, seq_map: dict[bytes, str], acc: str, name: str, cutoff: float):
        seq = pyhmmer.easel.TextSequence(name=b"s", sequence=seq_map[uid].replace("*", "")).digitize(alphabet)
        hmm, _, _ = builder.build(seq, bg)
        hmm.name, hmm.accession = name.encode(), acc.encode()
        hmm.cutoffs.gathering = (cutoff, cutoff)
        return hmm

    h_lam0 = build_hmm(p_meas_shared, prots_s0, "TEST0001", "LamMeas", 120.0)
    h_lam1 = build_hmm(p_both_shared, prots_s0, "TEST0002", "LamBoth", 120.0)
    h_mg0 = build_hmm(p_meas_s0, prots_s0, "TEST0003", "MgMeas", 120.0)
    h_t4 = build_hmm(t4_seed_uid, prots_s1, "T4M0001", "T4Head", 200.0)
    with open(hmm_path, "wb") as fh:
        for h in (h_lam0, h_lam1, h_mg0, h_t4):
            h.write(fh)
    with open(lam_query_hmm, "wb") as fh:
        h_lam0.write(fh)
    with open(t4_query_hmm, "wb") as fh:
        h_t4.write(fh)

    # Scan ONLY S0 with `test-1`
    domains.run(root, "S0", V, "test", "1", hmm_path, "gathering")

    # Augment S0's l1_domains with `inferred_from` rows for `p_both_shared`, `p_inf_shared`, `p_inf_s0` and pad
    # additional measured rows on M. genitalium proteins so l1_domains has >= 6 row groups at row_group_size=16,
    # keeping `p_none_shared` and `p_none_s0` free of any domain row.
    dom_file = root / "l1/domains/library=test-1/stage=S0/part-0.parquet"
    existing = pq.read_table(dom_file, partitioning=None).to_pylist()
    existing = [
        {**r, "inferred_from": None}
        for r in existing
        if r["protein_uid"] not in {p_inf_shared, p_inf_s0, p_none_shared, p_none_s0}
    ]
    template = dict(existing[0])

    def synth_dom(uid: bytes, acc: str, name: str, env_s: int, inferred_from: bytes | None) -> dict:
        return {**template, "protein_uid": uid, "model_acc": acc, "model_name": name,
                "env_start": env_s, "env_end": env_s + 50, "ali_start": env_s, "ali_end": env_s + 50,
                "inferred_from": inferred_from}

    existing.append(synth_dom(p_both_shared, "TEST0099", "InfExtra", 10, p_meas_shared))
    existing.append(synth_dom(p_inf_shared, "TEST0098", "LamInf", 0, p_meas_shared))
    existing.append(synth_dom(p_inf_s0, "TEST0097", "MgInf", 0, p_meas_s0))
    for extra_uid in mg_puids[10:95]:
        if extra_uid not in {p_inf_shared, p_inf_s0, p_none_shared, p_none_s0}:
            existing.append(synth_dom(extra_uid, "TEST0050", "MgPad", 0, None))
    existing.sort(key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
    pq.write_table(pa.Table.from_pylist(existing, schema=dom_spec.schema), dom_file,
                   row_group_size=16, compression="zstd", write_statistics=True)

    for s in ("S0", "S1"):
        for pf in sorted((root / f"l0seq/protein_seq/gene_call_version={V}/stage={s}").rglob("*.parquet")):
            _rewrite_small_row_groups(pf, schemas.table("protein_seq").schema, row_group_size=16)

    catalog.define_release(root, rel, ["S0", "S1"])
    _standin_clusters(root, rel, ["S0", "S1"])
    dim.run(root, rel, ["S0", "S1"])
    tokens.run(root, rel, ["S0", "S1"])
    postings.run(root, rel, ["S0", "S1"])
    man = catalog.build_manifest(root, rel, ["S0", "S1"])
    catalog.write_manifest(man, catalog.manifest_path(root, rel))
    cat = catalog.generate_catalog(man, root / f"catalog/release={rel}.duckdb")
    packs = base / "packs"
    pack.build(root, rel, models={"test-1": hmm_path}, packs=packs)

    return {
        "root": root,
        "rel": rel,
        "man": man,
        "cat": cat,
        "packs": packs,
        "hmm": hmm_path,
        "lam_query_hmm": lam_query_hmm,
        "t4_query_hmm": t4_query_hmm,
        "p_meas_shared": p_meas_shared,
        "p_both_shared": p_both_shared,
        "p_inf_shared": p_inf_shared,
        "p_none_shared": p_none_shared,
        "p_meas_s0": p_meas_s0,
        "p_inf_s0": p_inf_s0,
        "p_none_s0": p_none_s0,
        "q_prime_uid": q_prime_uid,
        "t4_seed_uid": t4_seed_uid,
        "lam_cu": lam_cu,
        "lam_mut_cu": lam_mut_cu,
        "t4_cu": t4_cu,
    }


def _brute_force_domain_oracle(root: Path, manifest: dict, library: str) -> tuple[dict[bytes, str], dict[bytes, str]]:
    """Independent brute-force reader over the manifest's Parquet files (AGENTS §5.2).
    Returns ({protein_uid: domain_state}, {protein_uid: pfam_architecture})."""
    annotated_stages = {
        f["partitions"]["stage"]
        for f in manifest["files"]
        if f["table"] == "l1_domains" and f["partitions"].get("library") == library
    }
    all_proteins: set[bytes] = set()
    annotated_proteins: set[bytes] = set()
    for f in manifest["files"]:
        if f["table"] == "protein_seq":
            st = f["partitions"]["stage"]
            uids = set(pq.read_table(root / f["path"], columns=["protein_uid"], partitioning=None).column(0).to_pylist())
            all_proteins |= uids
            if st in annotated_stages:
                annotated_proteins |= uids

    dom_by_prot: dict[bytes, list[dict]] = {}
    for f in manifest["files"]:
        if f["table"] == "l1_domains" and f["partitions"].get("library") == library:
            t = pq.read_table(root / f["path"], partitioning=None)
            has_inf = "inferred_from" in t.schema.names
            for r in t.to_pylist():
                dom_by_prot.setdefault(r["protein_uid"], []).append({
                    "env_start": r["env_start"],
                    "model_acc": r["model_acc"],
                    "model_name": r["model_name"],
                    "inferred_from": r["inferred_from"] if has_inf else None,
                })

    states: dict[bytes, str] = {}
    archs: dict[bytes, str] = {}
    for pu in all_proteins | set(dom_by_prot):
        rows = dom_by_prot.get(pu, [])
        if rows:
            if any(r["inferred_from"] is None for r in rows):
                states[pu] = "measured"
            else:
                states[pu] = "inferred"
            ordered = sorted({(r["env_start"], r["model_acc"], r["model_name"]) for r in rows})
            names: list[str] = []
            for _, _, nm in ordered:
                if not names or names[-1] != nm:
                    names.append(nm)
            archs[pu] = "+".join(names)
        elif pu in annotated_proteins:
            states[pu] = "none"
            archs[pu] = "-"
        else:
            states[pu] = "not scanned"
            archs[pu] = "?"
    return states, archs


@pytest.mark.covers("C-38b", "LY-06", "LY-09", "LY-10", "W-05", "W-06", "W-07")
def test_c38b_mixed_release_catalog_pruning_and_exports_match_brute_force_oracle(
    tmp_path, capsys, l1_domains_with_inferred
):
    fx = _build_mixed_fixture(tmp_path / "fx", l1_domains_with_inferred)
    root, rel, man, cat = fx["root"], fx["rel"], fx["man"], fx["cat"]

    want_state, want_arch = _brute_force_domain_oracle(root, man, "test-1")
    # Every one of the four states (plus shared annotated+called proteins!) is present in constructed truth
    assert set(want_state.values()) == {"measured", "inferred", "none", "not scanned"}
    assert want_state[fx["p_meas_shared"]] == "measured"
    assert want_state[fx["p_both_shared"]] == "measured"
    assert want_state[fx["p_inf_shared"]] == "inferred"
    assert want_state[fx["p_none_shared"]] == "none"  # in BOTH S0 (annotated) and S1 (called) -> 'none', NOT 'not scanned'!
    assert want_state[fx["p_meas_s0"]] == "measured"
    assert want_state[fx["p_inf_s0"]] == "inferred"
    assert want_state[fx["p_none_s0"]] == "none"
    assert want_state[fx["q_prime_uid"]] == "not scanned"
    assert want_state[fx["t4_seed_uid"]] == "not scanned"

    # --- (1) Catalog view `domain_state` and scalar macro `domain_state(protein_uid, library)` ---
    plain = duckdb.connect(str(cat), read_only=True)
    try:
        view_rows = plain.execute("SELECT protein_uid, library, domain_state FROM domain_state").fetchall()
        assert {(u, lib): st for u, lib, st in view_rows} == {(u, "test-1"): st for u, st in want_state.items()}
        for probe in (fx["p_meas_shared"], fx["p_both_shared"], fx["p_inf_shared"], fx["p_none_shared"],
                      fx["q_prime_uid"], fx["t4_seed_uid"]):
            assert plain.execute("SELECT domain_state(?, 'test-1')", [probe]).fetchone()[0] == want_state[probe]
        assert plain.execute("SELECT domain_state(?, 'other-lib')", [fx["p_meas_shared"]]).fetchone()[0] == "not scanned"
        assert plain.execute("SELECT domain_state(NULL, 'test-1')").fetchone()[0] is None
    finally:
        plain.close()

    con = catalog.open_catalog(cat, threads=1)
    try:
        got_view = con.execute("SELECT protein_uid, library, domain_state FROM domain_state").fetchall()
        assert {(u, lib): st for u, lib, st in got_view} == {(u, "test-1"): st for u, st in want_state.items()}
        con.register("_all_p", pa.table({"protein_uid": pa.array(sorted(want_state), pa.binary(16))}))
        got_macro = dict(con.execute("SELECT protein_uid, domain_state(protein_uid, 'test-1') FROM _all_p").fetchall())
        assert got_macro == want_state
        assert con.execute("SELECT domain_state(?, 'other-lib')", [fx["p_meas_shared"]]).fetchone()[0] == "not scanned"
        assert con.execute("SELECT domain_state(NULL, 'test-1')").fetchone()[0] is None

        # --- (2) O(1) row-group pruning via `catalog.profile_query` (AGENTS §4.4) ---
        p_rgs = sum(f["row_groups"] for f in man["files"] if f["table"] == "protein_seq")
        d_rgs = sum(f["row_groups"] for f in man["files"] if f["table"] == "l1_domains")
        assert p_rgs >= 40 and d_rgs >= 6, (p_rgs, d_rgs)

        for probe in (fx["p_meas_shared"], fx["p_inf_shared"], fx["p_none_shared"], fx["t4_seed_uid"]):
            lo = 1 if want_state[probe] != "not scanned" else 0
            # Literal lookup through `domain_state(?, 'test-1')`
            prof_lit = catalog.profile_query(con, "SELECT domain_state(?, 'test-1')", [probe])
            assert prof_lit["result"] == [(want_state[probe],)]
            scans_lit = [s for s in prof_lit["scans"] if s.get("lookup") == "bb_domain_state"]
            rg_lit = sum(s["row_groups_scanned"] for s in scans_lit)
            assert lo <= rg_lit <= 3, prof_lit["scans"]

            # Joined-table lookup through `domain_state(q.protein_uid, 'test-1')`
            con.register("_one_q", pa.table({"protein_uid": pa.array([probe], pa.binary(16))}))
            prof_join = catalog.profile_query(con, "SELECT domain_state(q.protein_uid, 'test-1') FROM _one_q q")
            assert prof_join["result"] == [(want_state[probe],)]
            scans_join = [s for s in prof_join["scans"] if s.get("lookup") == "bb_domain_state"]
            rg_join = sum(s["row_groups_scanned"] for s in scans_join)
            assert lo <= rg_join <= 3, prof_join["scans"]

            # Literal lookup on the `domain_state` SQL view
            prof_view = catalog.profile_query(
                con, "SELECT domain_state FROM domain_state WHERE protein_uid = ? AND library = 'test-1'", [probe]
            )
            assert prof_view["result"] == [(want_state[probe],)]
            assert prof_view["scans"]
            for sc in prof_view["scans"]:
                assert sc["row_groups_scanned"] <= 4, prof_view["scans"]
    finally:
        con.close()

    # Pre-C-38b catalog compatibility: drop `domain_state` view from a copy of the catalog file and verify `open_catalog`
    old_cat = tmp_path / "pre_c38b.duckdb"
    shutil.copyfile(cat, old_cat)
    rw = duckdb.connect(str(old_cat))
    try:
        rw.execute("DROP VIEW domain_state")
        rw.execute("DROP MACRO domain_state")
    finally:
        rw.close()
    con_old = catalog.open_catalog(old_cat, threads=1)
    try:
        assert dict(con_old.execute("SELECT protein_uid, domain_state FROM domain_state WHERE library = 'test-1'").fetchall()) == want_state
        assert con_old.execute("SELECT domain_state(?, 'test-1')", [fx["t4_seed_uid"]]).fetchone()[0] == "not scanned"
    finally:
        con_old.close()

    # --- (3) Exports on the mixed release (`hits_tsv`, `locus_tables`, `nbr.summarize`, `nbr flanks`, `search hmm`) ---
    tv = tapes.tape_version([V])
    probe_uids = [
        fx["p_meas_shared"], fx["p_both_shared"], fx["p_inf_shared"], fx["p_none_shared"],
        fx["p_meas_s0"], fx["p_inf_s0"], fx["p_none_s0"], fx["q_prime_uid"], fx["t4_seed_uid"],
    ]
    made_hits = [
        {"target_kind": "protein", "target_uid": u, "seq_bitscore": 100.0 - k, "seq_evalue_zref": 1e-20}
        for k, u in enumerate(probe_uids)
    ]
    hits_lines = [ln.split("\t") for ln in export.hits_tsv(
        root, made_hits, release=rel, tape_version=tv, gene_call_version=V
    )]
    h_cols = hits_lines[0]
    assert h_cols[h_cols.index("pfam_architecture") + 1] == "domain_state"
    seen_stages_for_none_shared = set()
    for row in hits_lines[1:]:
        rec = dict(zip(h_cols, row))
        pu = bytes.fromhex(rec["protein_uid"])
        assert rec["domain_state"] == want_state[pu], (rec["stage"], pu.hex(), rec)
        assert rec["pfam_architecture"] == want_arch[pu], (rec["stage"], pu.hex(), rec)
        if pu == fx["p_none_shared"]:
            seen_stages_for_none_shared.add(rec["stage"])
    # `p_none_shared` sits on BOTH S0 and S1 contigs, and in BOTH rows its state is 'none' (`-`), never '?'!
    assert seen_stages_for_none_shared == {"S0", "S1"}

    # Locus tables: check all 4 states + a non-CDS slot (`protein_uid=None` -> empty `pfam_architecture` and `domain_state`)
    nb_rows, _ = nbr.fetch(root, {fx["p_meas_shared"], fx["t4_seed_uid"]}, nbr.Policy.parse("3g"),
                           files=catalog.ReleaseFiles.of(root, rel))
    # Inject one synthetic non-CDS slot row in one locus to verify non-CDS stays empty in both columns
    non_cds = {**nb_rows[0], "offset": -99, "offset_rel": -99, "protein_uid": None,
               "feature_uid": b"\xff" * 16, "delta_bp": 50}
    loci_out = list(export.locus_tables(root, nb_rows + [non_cds], release=rel, gene_call_version=V))
    assert loci_out
    saw_question = saw_dash = saw_inferred = saw_measured = saw_non_cds = False
    for _, lines in loci_out:
        l_cols = lines[2].split("\t")
        assert l_cols == export.LOCUS_TSV_COLUMNS
        for ln in lines[3:]:
            rec = dict(zip(l_cols, ln.split("\t")))
            if not rec["protein_uid"]:
                assert rec["pfam_architecture"] == "" and rec["domain_state"] == ""
                saw_non_cds = True
                continue
            pu = bytes.fromhex(rec["protein_uid"])
            assert rec["domain_state"] == want_state[pu]
            assert rec["pfam_architecture"] == want_arch[pu]
            saw_question |= rec["pfam_architecture"] == "?" and rec["domain_state"] == "not scanned"
            saw_dash |= rec["pfam_architecture"] == "-" and rec["domain_state"] == "none"
            saw_inferred |= rec["domain_state"] == "inferred"
            saw_measured |= rec["domain_state"] == "measured"
    assert saw_question and saw_dash and saw_inferred and saw_measured and saw_non_cds

    # `nbr.summarize` on the mixed release separates `none` ("") and `not scanned` ("?")
    summary = nbr.summarize(root, {fx["p_meas_shared"], fx["t4_seed_uid"]}, nbr.Policy.parse("3g"),
                            by=("pfam", "nr50"), release=rel)
    pfam_unann = [r for r in summary["rows"] if r["level"] == "pfam" and r["unannotated"]]
    assert [(r["key"], r["model_name"]) for r in pfam_unann] == [("", None), ("?", "not scanned")]
    # Brute-force check of `not scanned` neighbour count in `summary`: on `lambda_q_mut` Q' is 1 unscanned neighbour,
    # and on T4 all 6 flanking genes of `t4_seed_uid` (3g each side, or truncated) are unscanned!
    nb_only = [r for r in nb_rows if r["offset"] != 0]
    want_ns_nb = sum(1 for r in nb_only if r["protein_uid"] is not None and want_state[r["protein_uid"]] == "not scanned")
    want_none_nb = sum(1 for r in nb_only if r["protein_uid"] is None or want_state[r["protein_uid"]] == "none")
    assert pfam_unann[0]["n_neighbors"] == want_none_nb
    assert pfam_unann[1]["n_neighbors"] == want_ns_nb and want_ns_nb >= 2

    # `nbr flanks` on the mixed release: unscanned and none genes both keep unique `~<feature_uid>` tokens
    flanks_dir = tmp_path / "flanks_mixed"
    assert cli.main([
        "nbr", "flanks", f"{fx['p_meas_shared'].hex()},{fx['t4_seed_uid'].hex()}",
        "--release", rel, "--root", str(root), "--window", "3g", "--level", "pfam", "--out", str(flanks_dir),
    ]) == 0
    capsys.readouterr()
    leaves_txts = [p.read_text() for p in sorted(flanks_dir.glob("groups/*/leaves.tsv"))]
    assert leaves_txts and all("?" not in t.split("\n", 1)[0] for t in leaves_txts)
    # Neither `-` nor `?` is ever treated as a shared family token in `leaves.tsv`!
    for txt in leaves_txts:
        for ln in txt.splitlines()[1:]:
            toks = ln.split("\t")[9:11]
            assert not any(tok.startswith("-:") or tok.startswith("?:") for side in toks for tok in side.split(",") if tok)

    # `search hmm` end-to-end on the mixed release (using T4 query HMM -> hits T4 in called stage S1!)
    pos_spec = ms._w07_spec(root, fx["t4_seed_uid"])
    hmm_out = tmp_path / "hmm_mixed"
    runs_dir = tmp_path / "runs"
    code = cli.main([
        "search", "hmm", "--root", str(root), "--release", rel, "--hmm", str(fx["t4_query_hmm"]),
        "--bits", "200", "--neighbors", "3", "--positive", pos_spec, "--out", str(hmm_out), "--cpus", "1",
        "--pack-dir", str(fx["packs"]), "--runs", str(runs_dir), "--work", str(tmp_path / "work"),
    ])
    capsys.readouterr()
    assert code == 0
    readme = (hmm_out / "README.txt").read_text()
    assert "domain_state" in readme and "not scanned" in readme
    htsv = [ln.split("\t") for ln in (hmm_out / "hits.tsv").read_text().splitlines()]
    assert {r[htsv[0].index("pfam_architecture")] for r in htsv[1:]} == {"?"}
    assert {r[htsv[0].index("domain_state")] for r in htsv[1:]} == {"not scanned"}
    ntsv = (hmm_out / "neighbors.tsv").read_text()
    assert "\npfam\t?\ttrue\t\t\tnot scanned\t" in "\n" + ntsv


def _strip_column(tsv_text: str, col_name: str, header_idx: int = 0) -> str:
    lines = tsv_text.splitlines()
    cols = lines[header_idx].split("\t")
    idx = cols.index(col_name)
    out = list(lines[:header_idx])
    for ln in lines[header_idx:]:
        parts = ln.split("\t")
        out.append("\t".join(parts[:idx] + parts[idx + 1:]))
    return "\n".join(out) + ("\n" if tsv_text.endswith("\n") else "")


@pytest.mark.covers("C-38b", "W-05", "W-06", "W-07")
def test_c38b_all_annotated_release_exports_are_byte_identical_except_for_domain_state_column(
    world, tmp_path, capsys
):
    """On a release whose stages are all annotated (`ms.world`), zero rows are `?` / `not scanned`, `hits.faa`,
    `neighbors.tsv` and `nbr flanks` have zero `?` or `not scanned` rows, and `hits.tsv` and `loci/*.tsv` stripped of
    `domain_state` equal the pre-C-38b TSV columns and values byte-for-byte."""
    root = ms._copy(world, tmp_path)
    man = catalog.build_manifest(root, ms.REL, ["S0"])
    catalog.write_manifest(man, catalog.manifest_path(root, ms.REL))
    packs = ms._packs_of(root)

    hmm, pos = ms._w07_hmm(world, tmp_path)
    spec = ms._w07_spec(root, pos)
    out = tmp_path / "out"
    runs = root.with_name(root.name + ".runs")
    assert cli.main([
        "search", "hmm", "--root", str(root), "--release", ms.REL, "--hmm", str(hmm), "--bits", "25",
        "--neighbors", "5", "--positive", spec, "--out", str(out), "--cpus", "1",
        "--pack-dir", str(packs), "--runs", str(runs), "--work", str(tmp_path / "work"),
    ]) == 0
    rid = capsys.readouterr().out.splitlines()[0].split("\t")[0]

    # 1. `neighbors.tsv` has NO `?` / `not scanned` row on an all-annotated release
    neighbors_txt = (out / "neighbors.tsv").read_text()
    assert "\t?\t" not in neighbors_txt and "not scanned" not in neighbors_txt

    # 2. `hits.tsv` has only `measured` and `none` (`-`), never `?` or `not scanned`, and stripping `domain_state`
    # matches the pre-C-38b column layout byte-for-byte against `_w05_oracle_genes` without `domain_state`
    hits_txt = (out / "hits.tsv").read_text()
    h_rows = [ln.split("\t") for ln in hits_txt.splitlines()]
    ds_idx = h_rows[0].index("domain_state")
    pa_idx = h_rows[0].index("pfam_architecture")
    assert {r[ds_idx] for r in h_rows[1:]} <= {"measured", "none"}
    assert "?" not in {r[pa_idx] for r in h_rows[1:]}

    old_hit_cols = [c for c in export.HITS_TSV_COLUMNS if c not in ("tier", "domain_state")]
    final_hits = ms._w05_final_hits(root, rid)
    assert final_hits
    oracle_genes = ms._w05_oracle_genes(root, final_hits, {})
    want_old_hits = "\n".join(ms._w05_tsv(old_hit_cols, oracle_genes)) + "\n"
    assert _strip_column(hits_txt, "domain_state", header_idx=0) == want_old_hits

    # 3. `loci/*.tsv` have no `?` or `not scanned`, and stripping `domain_state` matches the pre-C-38b oracle
    want_loci = dict(ms._w05_oracle_loci(root, final_hits, "5g"))
    for p in sorted((out / "loci").glob("locus_*.tsv")):
        txt = p.read_text()
        l_rows = [ln.split("\t") for ln in txt.splitlines()[2:]]
        lds = l_rows[0].index("domain_state")
        assert {r[lds] for r in l_rows[1:] if r[l_rows[0].index("protein_uid")]} <= {"measured", "none"}
        want_old_locus_txt = "\n".join(want_loci[p.name]) + "\n"
        if "domain_state" in want_loci[p.name][2].split("\t"):
            want_old_locus_txt = _strip_column(want_old_locus_txt, "domain_state", header_idx=2)
        assert _strip_column(txt, "domain_state", header_idx=2) == want_old_locus_txt

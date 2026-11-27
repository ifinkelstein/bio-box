"""Card A-33: release-level protein-family and syntenic-pair saturation tracker (`DEC-86 (1)`).

Tests cover:
1. Manifest-only reading (unmanifested stray Parquet files ignored; missing manifested files or archived releases refused).
2. Linear-contig edge-truncated CDS (`dist_to_contig_end_bp == 0`) excluded from `complete_orf` while circular-contig CDS
   (`dist_to_contig_end_bp == contig_len > 0`) are included.
3. Multi-copy identical or homologous ORFs within one `contig_dim.genome_uid` counted as 1 genome incidence (`g_c = 1`).
4. Good-Turing sample coverage `c_hat = 1 - f1 / n_incidences` (`1.0` at `0`) and bias-corrected Chao1
   `s_obs + f1 * (f1 - 1) / (2 * (f2 + 1))` (safe at `f2 == 0`) against an independent brute-force Python oracle, plus
   thread/bucket byte-identity.
5. Stage-prefix collector curve and marginal `k >= 3` / `k >= 5` decomposition into `brand_new_ge_k` + `promoted_to_ge_k`
   (both per-stage and `--earlier` multi-stage batch).
6. Strand-invariant adjacent CDS NR50 pairs from `contig_tapes` + `tape_tokens` + `contig_dim` (reverse-complement
   invariance, circular `(n_slots - 1, 0)` wrap iff `is_circular` and `n_slots >= 2`, skipping null non-CDS tokens while
   preserving full 64-bit `uint64` token precision).
7. `train_check.py` announcement formatting with and without `--saturation-summary`.
8. `history` number/label/formula validator, `render` deterministic `.png`/`.svg` output, and subprocess execution with
   native tools hidden from `PATH` (`PATH=/usr/bin:/bin`).
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import random
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import ids
from biobox.l3.tapes import build_tape, tape_version
from biobox.l3.tokens import token_id
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("A-33")

REPO = Path(__file__).resolve().parents[2]
POD_DIR = REPO / "scripts" / "stress" / "build_stage_pod"
SAT_SCRIPT = POD_DIR / "release_saturation.py"
CHK_SCRIPT = POD_DIR / "train_check.py"
HISTORY_JSON = REPO / "planning" / "measurements" / "release_saturation_history.json"

GCV = "pyrodigal-3.7.1-gv-0.3.2-p2"
TV = tape_version([GCV])


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _uid(tag: str) -> bytes:
    """Deterministic 16-byte digest for test fixtures."""
    return hashlib.sha256(tag.encode("utf-8")).digest()[:16]


def _uid_shared_prefix(low: int) -> bytes:
    """Adversarial 16-byte UID with upper 64 bits all zero (tests full 128-bit comparison)."""
    return b"\x00" * 8 + int(low).to_bytes(8, "big")


def _write_pq(
    root: Path,
    rel_path: str,
    spec_name: str,
    rows: list[dict],
    sort_keys: tuple[str, ...],
    *,
    row_group_size: int | None = None,
) -> dict:
    spec = table(spec_name)
    rows_sorted = sorted(rows, key=lambda r: tuple(r[k] for k in sort_keys))
    tbl = pa.Table.from_pylist(rows_sorted, schema=spec.schema)
    out_path = root / rel_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(tbl, out_path, compression="zstd", write_statistics=True, row_group_size=row_group_size)
    raw = out_path.read_bytes()
    pf = pq.ParquetFile(out_path)
    return {
        "path": rel_path,
        "table": spec_name,
        "bytes": len(raw),
        "rows": tbl.num_rows,
        "row_groups": pf.metadata.num_row_groups,
        "sha256": hashlib.sha256(raw).hexdigest(),
    }


def _build_synthetic_release(
    root: Path,
    release: str,
    stage_defs: list[dict],
    clusters_nr50: dict[bytes, bytes],
    *,
    token_map: dict[bytes, int] | None = None,
    row_group_size: int | None = None,
) -> dict:
    """Write a minimal manifest-complete synthetic release under `root` and return its manifest dict.

    Each entry in `stage_defs`:
      {
        "stage": "S0",
        "contigs": [
          {
            "contig_uid": bytes(16),
            "genome_uid": bytes(16),
            "contig_len": int,
            "is_circular": bool,
            "is_default_variant": bool (default True),
            "slots": [
              # CDS slot: {"kind": "cds", "start": int, "end": int, "strand": +1|-1,
              #            "protein_uid": bytes(16), "dist_to_contig_end_bp": int}
              # non-CDS slot: {"kind": "ncrna", "start": int, "end": int, "strand": +1|-1}
            ],
          },
          ...
        ],
      }
    """
    if token_map is None:
        distinct_clusters = sorted({c for c in clusters_nr50.values()})
        token_map = {c: token_id("nr50", 0, idx) for idx, c in enumerate(distinct_clusters)}

    files: list[dict] = []
    stages = [sd["stage"] for sd in stage_defs]
    cdim_rows: list[dict] = []
    tt_rows: list[dict] = []

    for sd in stage_defs:
        st = sd["stage"]
        ci_rows: list[dict] = []
        po_rows: list[dict] = []
        ct_rows: list[dict] = []

        for c in sd["contigs"]:
            cu = c["contig_uid"]
            gu = c["genome_uid"]
            clen = int(c["contig_len"])
            circ = bool(c["is_circular"])
            is_def = bool(c.get("is_default_variant", True))
            rd = hashlib.sha256(cu + clen.to_bytes(8, "big")).digest()
            ci_rows.append({
                "contig_uid": cu,
                "raw_digest": rd,
                "contig_len": clen,
                "is_default_variant": is_def,
            })
            if not is_def:
                continue

            cdim_rows.append({
                "contig_uid": cu,
                "genome_uid": gu,
                "contig_len": clen,
                "is_circular": circ,
                "tier": None,
                "assembly_quality": "contig",
                "tax_order": None,
                "is_species_rep": True,
                "n_genomes": 1,
            })

            slots = sorted(
                c.get("slots", []),
                key=lambda s: (int(s["start"]), int(s["end"]), int(s["strand"]), ids.KIND[s.get("kind", "cds")]),
            )
            coords = [(int(s["start"]), int(s["end"]), int(s["strand"])) for s in slots]
            kinds = [ids.KIND[s.get("kind", "cds")] for s in slots]
            tape_cols = build_tape(coords, clen, circ, kinds)
            pus: list[bytes | None] = []
            fus: list[bytes] = []
            nr50_toks: list[int | None] = []
            nr90_toks: list[int | None] = []
            fam_toks: list[int | None] = []
            dom_toks: list[list[int]] = []

            for s in slots:
                kname = s.get("kind", "cds")
                fu = ids.feature_uid(cu, int(s["start"]), int(s["end"]), int(s["strand"]), ids.KIND[kname])
                fus.append(fu)
                if kname == "cds":
                    pu = s["protein_uid"]
                    pus.append(pu)
                    cu_nr50 = clusters_nr50[pu]
                    tok50 = token_map[cu_nr50]
                    nr50_toks.append(tok50)
                    nr90_toks.append(token_id("nr90", 0, 0))
                    fam_toks.append(None)
                    dom_toks.append([])
                    po_rows.append({
                        "protein_uid": pu,
                        "contig_uid": cu,
                        "start": int(s["start"]),
                        "end": int(s["end"]),
                        "strand": int(s["strand"]),
                        "dist_to_contig_end_bp": int(s["dist_to_contig_end_bp"]),
                        "feature_uid": fu,
                    })
                else:
                    pus.append(None)
                    nr50_toks.append(None)
                    nr90_toks.append(None)
                    fam_toks.append(None)
                    dom_toks.append([])

            ct_rows.append({
                "contig_uid": cu,
                "contig_len": clen,
                "is_circular": circ,
                "n_slots": len(slots),
                **tape_cols,
                "protein_uids": pus,
                "feature_uids": fus,
            })
            tt_rows.append({
                "contig_uid": cu,
                "family_tokens": fam_toks,
                "nr50_tokens": nr50_toks,
                "nr90_tokens": nr90_toks,
                "domain_family_tokens": dom_toks,
            })

        e_ci = _write_pq(root, f"l0seq/contig_index/stage={st}/part-0.parquet", "contig_index", ci_rows,
                         ("contig_uid", "raw_digest"))
        e_ci["partitions"] = {"stage": st}
        files.append(e_ci)

        e_po = _write_pq(
            root,
            f"l1/protein_occurrence/gene_call_version={GCV}/stage={st}/part-0.parquet",
            "protein_occurrence",
            po_rows,
            ("protein_uid", "contig_uid", "start", "end", "strand"),
            row_group_size=row_group_size,
        )
        e_po["partitions"] = {"gene_call_version": GCV, "stage": st}
        files.append(e_po)

        e_ct = _write_pq(
            root,
            f"l3/tapes/tape_version={TV}/stage={st}/part-0.parquet",
            "contig_tapes",
            ct_rows,
            ("contig_uid",),
            row_group_size=row_group_size,
        )
        e_ct["partitions"] = {"tape_version": TV, "stage": st}
        files.append(e_ct)

    # member_cluster holds both nr50 and nr90 sorted by (level, member_uid)
    mc_rows: list[dict] = []
    for pu, cu_id in clusters_nr50.items():
        mc_rows.append({"member_uid": pu, "level": "nr50", "cluster_uid": cu_id, "rep_uid": cu_id})
        mc_rows.append({"member_uid": pu, "level": "nr90", "cluster_uid": pu, "rep_uid": pu})
    e_mc = _write_pq(
        root,
        f"l2/member_cluster/release={release}/part-0.parquet",
        "member_cluster",
        mc_rows,
        ("level", "member_uid"),
        row_group_size=row_group_size,
    )
    e_mc["partitions"] = {"release": release}
    files.append(e_mc)

    e_cd = _write_pq(
        root,
        f"dim/contig_dim/release={release}/part-0.parquet",
        "contig_dim",
        cdim_rows,
        ("contig_uid",),
        row_group_size=row_group_size,
    )
    e_cd["partitions"] = {"release": release}
    files.append(e_cd)

    e_tt = _write_pq(
        root,
        f"l3/tape_tokens/tape_version={TV}/release={release}/part-0.parquet",
        "tape_tokens",
        tt_rows,
        ("contig_uid",),
        row_group_size=row_group_size,
    )
    e_tt["partitions"] = {"tape_version": TV, "release": release}
    files.append(e_tt)

    manifest = {
        "manifest_version": 3,
        "release_id": release,
        "stages": stages,
        "files": files,
        "artifacts": [],
    }
    mpath = root / "manifests" / f"release={release}.json"
    mpath.parent.mkdir(parents=True, exist_ok=True)
    mpath.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


# --- Independent brute-force Python oracle (AGENTS.md §5.2) --------------------------------------------------


def _oracle_stats_from_counts(counts: dict[object, int]) -> dict:
    pos = [g for g in counts.values() if g >= 1]
    s_obs = len(pos)
    f1 = sum(1 for g in pos if g == 1)
    f2 = sum(1 for g in pos if g == 2)
    ge3 = sum(1 for g in pos if g >= 3)
    ge5 = sum(1 for g in pos if g >= 5)
    n_inc = sum(pos)
    c_hat = 1.0 if n_inc == 0 else 1.0 - (f1 / n_inc)
    chao1 = float(s_obs) + (f1 * (f1 - 1)) / (2.0 * (f2 + 1))
    return {
        "s_obs": s_obs,
        "f1": f1,
        "f2": f2,
        "ge3": ge3,
        "ge5": ge5,
        "n_incidences": n_inc,
        "c_hat": c_hat,
        "chao1": chao1,
    }


def _oracle_canonical_pair(ta: int, sa: int, tb: int, sb: int) -> tuple[tuple[int, int], tuple[int, int]]:
    fwd = ((int(ta), int(sa)), (int(tb), int(sb)))
    rev = ((int(tb), -int(sb)), (int(ta), -int(sa)))
    return fwd if fwd <= rev else rev


def _oracle_saturation(
    stage_defs: list[dict],
    clusters_nr50: dict[bytes, bytes],
    token_map: dict[bytes, int],
    *,
    earlier_stages: list[str] | None = None,
) -> dict:
    """Brute-force pure-Python oracle computing exact prefix and marginal saturation metrics."""
    cum_all_genomes: dict[bytes, set[bytes]] = {}
    cum_comp_genomes: dict[bytes, set[bytes]] = {}
    cum_pair_genomes: dict[tuple, set[bytes]] = {}

    prev_comp_counts: dict[bytes, int] = {}
    prev_pair_counts: dict[tuple, int] = {}
    prev_all_s_obs = 0
    prev_comp_s_obs = 0
    prev_pair_s_obs = 0

    cum_bases = 0
    cum_contigs = 0
    stage_curve = []

    snap_at_earlier_comp: dict[bytes, int] | None = None
    snap_at_earlier_pair: dict[tuple, int] | None = None
    snap_at_earlier_bases = 0
    if earlier_stages is not None and len(earlier_stages) == 0:
        snap_at_earlier_comp = {}
        snap_at_earlier_pair = {}

    for idx, sd in enumerate(stage_defs):
        st = sd["stage"]
        st_bases = 0
        st_contigs = 0
        st_occ_all = 0
        st_occ_comp = 0

        for c in sd["contigs"]:
            if not c.get("is_default_variant", True):
                continue
            st_contigs += 1
            st_bases += int(c["contig_len"])
            gu = c["genome_uid"]
            circ = bool(c["is_circular"])
            slots = sorted(
                c.get("slots", []),
                key=lambda s: (int(s["start"]), int(s["end"]), int(s["strand"]), ids.KIND[s.get("kind", "cds")]),
            )
            toks: list[int | None] = []
            strands: list[int] = []
            for s in slots:
                strands.append(int(s["strand"]))
                if s.get("kind", "cds") == "cds":
                    pu = s["protein_uid"]
                    cu_id = clusters_nr50[pu]
                    toks.append(token_map[cu_id])
                    st_occ_all += 1
                    cum_all_genomes.setdefault(cu_id, set()).add(gu)
                    if int(s["dist_to_contig_end_bp"]) > 0:
                        st_occ_comp += 1
                        cum_comp_genomes.setdefault(cu_id, set()).add(gu)
                else:
                    toks.append(None)

            n = len(slots)
            if n >= 2:
                adj_indices = [(i, i + 1) for i in range(n - 1)]
                if circ:
                    adj_indices.append((n - 1, 0))
                for ia, ib in adj_indices:
                    if toks[ia] is None or toks[ib] is None:
                        continue
                    cp = _oracle_canonical_pair(toks[ia], strands[ia], toks[ib], strands[ib])
                    cum_pair_genomes.setdefault(cp, set()).add(gu)

        cum_bases += st_bases
        cum_contigs += st_contigs

        cur_all_counts = {k: len(v) for k, v in cum_all_genomes.items()}
        cur_comp_counts = {k: len(v) for k, v in cum_comp_genomes.items()}
        cur_pair_counts = {k: len(v) for k, v in cum_pair_genomes.items()}

        s_all = _oracle_stats_from_counts(cur_all_counts)
        s_comp = _oracle_stats_from_counts(cur_comp_counts)
        s_pair_raw = _oracle_stats_from_counts(cur_pair_counts)

        def _decomp(prev_d: dict, cur_d: dict, k: int) -> tuple[int, int, int]:
            brand = sum(1 for item, cnt in cur_d.items() if cnt >= k and prev_d.get(item, 0) == 0)
            prom = sum(1 for item, cnt in cur_d.items() if cnt >= k and 1 <= prev_d.get(item, 0) < k)
            return brand + prom, brand, prom

        net3, b3, p3 = _decomp(prev_comp_counts, cur_comp_counts, 3)
        net5, b5, p5 = _decomp(prev_comp_counts, cur_comp_counts, 5)
        pnet3, pb3, pp3 = _decomp(prev_pair_counts, cur_pair_counts, 3)
        pnet5, pb5, pp5 = _decomp(prev_pair_counts, cur_pair_counts, 5)

        gbp = st_bases / 1e9
        stage_curve.append({
            "stage": st,
            "stage_index": idx,
            "stored_contigs": st_contigs,
            "stored_bases": st_bases,
            "cum_stored_contigs": cum_contigs,
            "cum_stored_bases": cum_bases,
            "protein_occurrences_all": st_occ_all,
            "protein_occurrences_complete": st_occ_comp,
            "protein_occurrences_edge_truncated": st_occ_all - st_occ_comp,
            "all_orf": s_all,
            "complete_orf": s_comp,
            "syntenic_pairs": {
                "pairs_ge1": s_pair_raw["s_obs"],
                "pairs_f1": s_pair_raw["f1"],
                "pairs_f2": s_pair_raw["f2"],
                "pairs_ge3": s_pair_raw["ge3"],
                "pairs_ge5": s_pair_raw["ge5"],
                "n_incidences": s_pair_raw["n_incidences"],
                "c_hat_pairs": s_pair_raw["c_hat"],
                "chao1_pairs": s_pair_raw["chao1"],
            },
            "marginal": {
                "new_all_orf_ge1": s_all["s_obs"] - prev_all_s_obs,
                "new_complete_orf_ge1": s_comp["s_obs"] - prev_comp_s_obs,
                "net_new_ge3": net3,
                "brand_new_ge3": b3,
                "promoted_to_ge3": p3,
                "net_new_ge3_per_gbp": None if st_bases == 0 else net3 / gbp,
                "net_new_ge5": net5,
                "brand_new_ge5": b5,
                "promoted_to_ge5": p5,
                "net_new_ge5_per_gbp": None if st_bases == 0 else net5 / gbp,
                "new_pairs_ge1": s_pair_raw["s_obs"] - prev_pair_s_obs,
                "net_new_pairs_ge3": pnet3,
                "brand_new_pairs_ge3": pb3,
                "promoted_to_pairs_ge3": pp3,
                "net_new_pairs_ge3_per_gbp": None if st_bases == 0 else pnet3 / gbp,
                "net_new_pairs_ge5": pnet5,
                "brand_new_pairs_ge5": pb5,
                "promoted_to_pairs_ge5": pp5,
            },
        })

        prev_comp_counts = dict(cur_comp_counts)
        prev_pair_counts = dict(cur_pair_counts)
        prev_all_s_obs = s_all["s_obs"]
        prev_comp_s_obs = s_comp["s_obs"]
        prev_pair_s_obs = s_pair_raw["s_obs"]

        if earlier_stages is not None and idx + 1 == len(earlier_stages):
            snap_at_earlier_comp = dict(cur_comp_counts)
            snap_at_earlier_pair = dict(cur_pair_counts)
            snap_at_earlier_bases = cum_bases

    since_earlier = None
    if snap_at_earlier_comp is not None and snap_at_earlier_pair is not None:
        final_comp = {k: len(v) for k, v in cum_comp_genomes.items()}
        final_pair = {k: len(v) for k, v in cum_pair_genomes.items()}
        new_b = cum_bases - snap_at_earlier_bases
        e_net3 = sum(1 for item, cnt in final_comp.items() if cnt >= 3 and snap_at_earlier_comp.get(item, 0) < 3)
        e_b3 = sum(1 for item, cnt in final_comp.items() if cnt >= 3 and snap_at_earlier_comp.get(item, 0) == 0)
        e_p3 = sum(1 for item, cnt in final_comp.items() if cnt >= 3 and 1 <= snap_at_earlier_comp.get(item, 0) < 3)
        e_net5 = sum(1 for item, cnt in final_comp.items() if cnt >= 5 and snap_at_earlier_comp.get(item, 0) < 5)
        e_b5 = sum(1 for item, cnt in final_comp.items() if cnt >= 5 and snap_at_earlier_comp.get(item, 0) == 0)
        e_p5 = sum(1 for item, cnt in final_comp.items() if cnt >= 5 and 1 <= snap_at_earlier_comp.get(item, 0) < 5)
        ep_net3 = sum(1 for item, cnt in final_pair.items() if cnt >= 3 and snap_at_earlier_pair.get(item, 0) < 3)
        ep_b3 = sum(1 for item, cnt in final_pair.items() if cnt >= 3 and snap_at_earlier_pair.get(item, 0) == 0)
        ep_p3 = sum(1 for item, cnt in final_pair.items() if cnt >= 3 and 1 <= snap_at_earlier_pair.get(item, 0) < 3)
        since_earlier = {
            "new_stored_bases": new_b,
            "net_new_ge3": e_net3,
            "brand_new_ge3": e_b3,
            "promoted_to_ge3": e_p3,
            "net_new_ge5": e_net5,
            "brand_new_ge5": e_b5,
            "promoted_to_ge5": e_p5,
            "net_new_pairs_ge3": ep_net3,
            "brand_new_pairs_ge3": ep_b3,
            "promoted_to_pairs_ge3": ep_p3,
        }

    return {"stage_curve": stage_curve, "since_earlier": since_earlier}


def _assert_matches_oracle(got: dict, oracle: dict) -> None:
    assert len(got["stage_curve"]) == len(oracle["stage_curve"])
    for gs, os_ in zip(got["stage_curve"], oracle["stage_curve"]):
        assert gs["stage"] == os_["stage"]
        assert gs["stored_contigs"] == os_["stored_contigs"]
        assert gs["stored_bases"] == os_["stored_bases"]
        assert gs["cum_stored_contigs"] == os_["cum_stored_contigs"]
        assert gs["cum_stored_bases"] == os_["cum_stored_bases"]
        assert gs["protein_occurrences_all"] == os_["protein_occurrences_all"]
        assert gs["protein_occurrences_complete"] == os_["protein_occurrences_complete"]
        assert gs["protein_occurrences_edge_truncated"] == os_["protein_occurrences_edge_truncated"]

        for section in ("all_orf", "complete_orf"):
            for k in ("s_obs", "f1", "f2", "ge3", "ge5", "n_incidences"):
                assert gs[section][k] == os_[section][k], (gs["stage"], section, k)
            assert gs[section]["c_hat"] == pytest.approx(os_[section]["c_hat"], abs=1e-12)
            assert gs[section]["chao1"] == pytest.approx(os_[section]["chao1"], abs=1e-9)

        for k in ("pairs_ge1", "pairs_f1", "pairs_f2", "pairs_ge3", "pairs_ge5", "n_incidences"):
            assert gs["syntenic_pairs"][k] == os_["syntenic_pairs"][k], (gs["stage"], k)
        assert gs["syntenic_pairs"]["c_hat_pairs"] == pytest.approx(os_["syntenic_pairs"]["c_hat_pairs"], abs=1e-12)
        assert gs["syntenic_pairs"]["chao1_pairs"] == pytest.approx(os_["syntenic_pairs"]["chao1_pairs"], abs=1e-9)

        gm, om = gs["marginal"], os_["marginal"]
        for k in (
            "new_all_orf_ge1", "new_complete_orf_ge1",
            "net_new_ge3", "brand_new_ge3", "promoted_to_ge3",
            "net_new_ge5", "brand_new_ge5", "promoted_to_ge5",
            "new_pairs_ge1", "net_new_pairs_ge3", "brand_new_pairs_ge3", "promoted_to_pairs_ge3",
            "net_new_pairs_ge5", "brand_new_pairs_ge5", "promoted_to_pairs_ge5",
        ):
            assert gm[k] == om[k], (gs["stage"], "marginal", k)
        for k in ("net_new_ge3_per_gbp", "net_new_ge5_per_gbp", "net_new_pairs_ge3_per_gbp"):
            if om[k] is None:
                assert gm[k] is None
            else:
                assert gm[k] == pytest.approx(om[k], rel=1e-9, abs=1e-9)

    if oracle["since_earlier"] is not None:
        se = got["since_earlier"]
        assert se is not None
        for k, v in oracle["since_earlier"].items():
            assert se[k] == v, ("since_earlier", k)


# --- Tests ---------------------------------------------------------------------------------------------------


def test_a33_manifest_only_reading_ignores_stray_parquet_and_refuses_missing(tmp_path):
    rs = _load_module("release_saturation", SAT_SCRIPT)
    root = tmp_path / "root"
    p1, p2, p_stray = _uid("p1"), _uid("p2"), _uid("p_stray")
    c1, c_stray = _uid("c1"), _uid("c_stray")
    clusters = {p1: c1, p2: c1}
    tok_map = {c1: token_id("nr50", 0, 0)}

    stages = [{
        "stage": "S0",
        "contigs": [{
            "contig_uid": _uid("cu0"),
            "genome_uid": _uid("g0"),
            "contig_len": 5000,
            "is_circular": False,
            "slots": [
                {"start": 100, "end": 400, "strand": 1, "protein_uid": p1, "dist_to_contig_end_bp": 100},
                {"start": 500, "end": 900, "strand": -1, "protein_uid": p2, "dist_to_contig_end_bp": 500},
            ],
        }],
    }]
    _build_synthetic_release(root, "r26.11.1", stages, clusters, token_map=tok_map)

    # Plant unmanifested stray part-99.parquet files in every partition directory (AGENTS.md §3.9)
    cu_stray = _uid("cu_stray")
    gu_stray = _uid("g_stray")
    _write_pq(root, "l0seq/contig_index/stage=S0/part-99.parquet", "contig_index", [{
        "contig_uid": cu_stray, "raw_digest": b"\xaa" * 32, "contig_len": 999999, "is_default_variant": True,
    }], ("contig_uid", "raw_digest"))
    _write_pq(root, f"l1/protein_occurrence/gene_call_version={GCV}/stage=S0/part-99.parquet", "protein_occurrence", [{
        "protein_uid": p_stray, "contig_uid": cu_stray, "start": 10, "end": 100, "strand": 1,
        "dist_to_contig_end_bp": 10, "feature_uid": _uid("fu_stray"),
    }], ("protein_uid",))
    _write_pq(root, "l2/member_cluster/release=r26.11.1/part-99.parquet", "member_cluster", [{
        "member_uid": p_stray, "level": "nr50", "cluster_uid": c_stray, "rep_uid": c_stray,
    }], ("level", "member_uid"))
    _write_pq(root, "dim/contig_dim/release=r26.11.1/part-99.parquet", "contig_dim", [{
        "contig_uid": cu_stray, "genome_uid": gu_stray, "contig_len": 999999, "is_circular": False,
        "tier": None, "assembly_quality": "contig", "tax_order": None, "is_species_rep": True, "n_genomes": 1,
    }], ("contig_uid",))

    out_dir = tmp_path / "out"
    tmp_dir = tmp_path / "tmp"
    rc = rs.main(["measure", "--root", str(root), "--release", "r26.11.1", "--out", str(out_dir),
                  "--tmp", str(tmp_dir), "--threads", "2"])
    assert rc == 0
    got = json.loads((out_dir / "saturation.json").read_text(encoding="utf-8"))
    oracle = _oracle_saturation(stages, clusters, tok_map)
    _assert_matches_oracle(got, oracle)
    assert got["totals"]["stored_bases"] == 5000  # stray 999,999 bp ignored

    # Refuses --tmp inside --root
    with pytest.raises((ValueError, SystemExit)):
        rs.measure_release(root=root, release="r26.11.1", out=tmp_path / "out2", tmp=root / "bad_tmp")

    # Refuses retired release
    arch = root / "archived" / "release=r26.11.1.json"
    arch.parent.mkdir(parents=True, exist_ok=True)
    arch.write_text('{"release": "r26.11.1"}\n', encoding="utf-8")
    with pytest.raises(FileNotFoundError):
        rs.measure_release(root=root, release="r26.11.1", out=tmp_path / "out3", tmp=tmp_path / "tmp3")
    arch.unlink()

    # Refuses missing manifested file
    (root / "dim/contig_dim/release=r26.11.1/part-0.parquet").unlink()
    with pytest.raises(FileNotFoundError):
        rs.measure_release(root=root, release="r26.11.1", out=tmp_path / "out4", tmp=tmp_path / "tmp4")


def test_a33_linear_edge_truncated_vs_circular_complete_orfs(tmp_path):
    rs = _load_module("release_saturation", SAT_SCRIPT)
    root = tmp_path / "root"

    p_edge_5p = _uid("p_edge_5p")
    p_edge_3p = _uid("p_edge_3p")
    p_lin_ok = _uid("p_lin_ok")
    p_circ_wrap = _uid("p_circ_wrap")
    clusters = {
        p_edge_5p: _uid("c_edge_5p"),
        p_edge_3p: _uid("c_edge_3p"),
        p_lin_ok: _uid("c_lin_ok"),
        p_circ_wrap: _uid("c_circ_wrap"),
    }
    tok_map = {cid: token_id("nr50", 0, i) for i, cid in enumerate(sorted(clusters.values()))}

    stages = [{
        "stage": "S0",
        "contigs": [
            {
                "contig_uid": _uid("lin_contig"),
                "genome_uid": _uid("g_lin"),
                "contig_len": 3000,
                "is_circular": False,
                "slots": [
                    {"start": 0, "end": 450, "strand": 1, "protein_uid": p_edge_5p, "dist_to_contig_end_bp": 0},
                    {"start": 600, "end": 1200, "strand": -1, "protein_uid": p_lin_ok, "dist_to_contig_end_bp": 600},
                    {"start": 2500, "end": 3000, "strand": 1, "protein_uid": p_edge_3p, "dist_to_contig_end_bp": 0},
                ],
            },
            {
                "contig_uid": _uid("circ_contig"),
                "genome_uid": _uid("g_circ"),
                "contig_len": 4000,
                "is_circular": True,
                "slots": [
                    # Origin-crossing circular CDS: end > contig_len, dist_to_contig_end_bp == contig_len > 0
                    {"start": 3800, "end": 4300, "strand": 1, "protein_uid": p_circ_wrap, "dist_to_contig_end_bp": 4000},
                ],
            },
        ],
    }]
    _build_synthetic_release(root, "r26.11.1", stages, clusters, token_map=tok_map)
    got = rs.measure_release(root=root, release="r26.11.1", out=tmp_path / "out", tmp=tmp_path / "tmp")
    oracle = _oracle_saturation(stages, clusters, tok_map)
    _assert_matches_oracle(got, oracle)

    assert got["totals"]["all_orf"]["s_obs"] == 4
    assert got["totals"]["complete_orf"]["s_obs"] == 2  # only p_lin_ok + p_circ_wrap
    assert got["totals"]["protein_occurrences_edge_truncated"] == 2
    assert got["totals"]["protein_occurrences_complete"] == 2


def test_a33_multicopy_identical_orfs_in_same_genome_count_as_one_incidence(tmp_path):
    rs = _load_module("release_saturation", SAT_SCRIPT)
    root = tmp_path / "root"

    # Two homologs p_a, p_b in the SAME NR50 cluster c_dup, repeated 5 times across 2 contigs and 2 stages
    # of the SAME genome g_same -> must still have g_c = 1 (singleton f1 = 1, f2 = 0, ge3 = 0)!
    p_a, p_b = _uid("p_a"), _uid("p_b")
    c_dup = _uid("c_dup")
    g_same = _uid("g_same")
    clusters = {p_a: c_dup, p_b: c_dup}
    tok_map = {c_dup: token_id("nr50", 0, 1)}

    stages = [
        {
            "stage": "S0",
            "contigs": [
                {
                    "contig_uid": _uid("cu_1"),
                    "genome_uid": g_same,
                    "contig_len": 5000,
                    "is_circular": False,
                    "slots": [
                        {"start": 100, "end": 400, "strand": 1, "protein_uid": p_a, "dist_to_contig_end_bp": 100},
                        {"start": 500, "end": 800, "strand": 1, "protein_uid": p_a, "dist_to_contig_end_bp": 500},
                        {"start": 900, "end": 1200, "strand": -1, "protein_uid": p_b, "dist_to_contig_end_bp": 900},
                    ],
                },
                {
                    "contig_uid": _uid("cu_2"),
                    "genome_uid": g_same,
                    "contig_len": 3000,
                    "is_circular": True,
                    "slots": [
                        {"start": 200, "end": 600, "strand": 1, "protein_uid": p_a, "dist_to_contig_end_bp": 3000},
                    ],
                },
            ],
        },
        {
            "stage": "S1",
            "contigs": [
                {
                    "contig_uid": _uid("cu_3"),
                    "genome_uid": g_same,
                    "contig_len": 4000,
                    "is_circular": False,
                    "slots": [
                        {"start": 300, "end": 700, "strand": 1, "protein_uid": p_b, "dist_to_contig_end_bp": 300},
                    ],
                },
            ],
        },
    ]
    _build_synthetic_release(root, "r26.11.2", stages, clusters, token_map=tok_map)
    got = rs.measure_release(root=root, release="r26.11.2", out=tmp_path / "out", tmp=tmp_path / "tmp")
    oracle = _oracle_saturation(stages, clusters, tok_map)
    _assert_matches_oracle(got, oracle)

    comp = got["totals"]["complete_orf"]
    assert comp["s_obs"] == 1
    assert comp["f1"] == 1
    assert comp["f2"] == 0
    assert comp["ge3"] == 0
    assert comp["n_incidences"] == 1
    assert comp["c_hat"] == pytest.approx(0.0)


def test_a33_good_turing_coverage_and_chao1_edge_cases_and_seeded_oracle(tmp_path):
    rs = _load_module("release_saturation", SAT_SCRIPT)

    # Edge case 1: empty incidences -> c_hat == 1.0, chao1 == 0.0
    root_empty = tmp_path / "empty"
    p_only_edge = _uid("p_only_edge")
    c_only = _uid("c_only")
    stages_empty = [{
        "stage": "S0",
        "contigs": [{
            "contig_uid": _uid("cu_e"),
            "genome_uid": _uid("g_e"),
            "contig_len": 1000,
            "is_circular": False,
            "slots": [
                {"start": 0, "end": 300, "strand": 1, "protein_uid": p_only_edge, "dist_to_contig_end_bp": 0},
            ],
        }],
    }]
    _build_synthetic_release(root_empty, "r26.11.1", stages_empty, {p_only_edge: c_only})
    got_empty = rs.measure_release(root=root_empty, release="r26.11.1", out=tmp_path / "out_e", tmp=tmp_path / "tmp_e")
    assert got_empty["totals"]["complete_orf"]["n_incidences"] == 0
    assert got_empty["totals"]["complete_orf"]["c_hat"] == 1.0
    assert got_empty["totals"]["complete_orf"]["chao1"] == 0.0
    assert got_empty["totals"]["syntenic_pairs"]["n_incidences"] == 0
    assert got_empty["totals"]["syntenic_pairs"]["c_hat_pairs"] == 1.0

    # Seeded multi-stage world (seed = 20261127) including adversarial shared-64bit-prefix UIDs
    rng = random.Random(20261127)
    n_clusters = 45
    cluster_uids = [_uid_shared_prefix(i + 1) if i < 15 else _uid(f"cl_{i}") for i in range(n_clusters)]
    proteins: list[bytes] = []
    clusters_map: dict[bytes, bytes] = {}
    for i in range(120):
        pu = _uid_shared_prefix(10_000 + i) if i < 30 else _uid(f"pr_{i}")
        proteins.append(pu)
        clusters_map[pu] = cluster_uids[rng.randrange(n_clusters)]
    tok_map = {cid: token_id("nr50", 0, idx) for idx, cid in enumerate(sorted(set(clusters_map.values())))}

    genome_uids = [_uid_shared_prefix(50_000 + g) if g < 8 else _uid(f"gn_{g}") for g in range(24)]
    stages_seeded = []
    for s_idx in range(4):
        contigs = []
        for c_idx in range(12):
            cu = _uid_shared_prefix(90_000 + s_idx * 100 + c_idx) if c_idx < 4 else _uid(f"cu_{s_idx}_{c_idx}")
            gu = genome_uids[rng.randrange(len(genome_uids))]
            circ = rng.random() < 0.3
            clen = 20_000
            n_sl = rng.randint(0, 7)
            slots = []
            for sl in range(n_sl):
                start = sl * 2000 + (0 if (not circ and sl == 0 and rng.random() < 0.35) else 100)
                end = start + 900
                strand = 1 if rng.random() < 0.5 else -1
                if rng.random() < 0.15:
                    slots.append({"kind": "ncrna", "start": start, "end": end, "strand": strand})
                else:
                    dist = clen if circ else (0 if start == 0 else min(start, clen - end))
                    slots.append({
                        "kind": "cds",
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "protein_uid": proteins[rng.randrange(len(proteins))],
                        "dist_to_contig_end_bp": dist,
                    })
            contigs.append({"contig_uid": cu, "genome_uid": gu, "contig_len": clen, "is_circular": circ, "slots": slots})
        stages_seeded.append({"stage": f"S{s_idx}", "contigs": contigs})

    root_s = tmp_path / "seeded"
    # row_group_size=12 forces multiple Parquet row groups on member_cluster (240 rows -> 20 RGs),
    # tape_tokens (48 rows -> 4 RGs), protein_occurrence, and contig_tapes so threads=4 exercises
    # multi-slice ProcessPoolExecutor row-group partitioning and min/max pruning.
    _build_synthetic_release(root_s, "r26.11.2", stages_seeded, clusters_map, token_map=tok_map, row_group_size=12)
    # Also write an earlier manifest over S0,S1
    (root_s / "manifests" / "release=r26.11.1.json").write_text(
        json.dumps({"manifest_version": 3, "release_id": "r26.11.1", "stages": ["S0", "S1"], "files": [], "artifacts": []}),
        encoding="utf-8",
    )

    out1, out4 = tmp_path / "out_t1", tmp_path / "out_t4"
    got1 = rs.measure_release(
        root=root_s, release="r26.11.2", earlier="r26.11.1", out=out1, tmp=tmp_path / "t1", threads=1, buckets=1
    )
    got4 = rs.measure_release(
        root=root_s, release="r26.11.2", earlier="r26.11.1", out=out4, tmp=tmp_path / "t4", threads=4, buckets=8
    )
    oracle = _oracle_saturation(stages_seeded, clusters_map, tok_map, earlier_stages=["S0", "S1"])
    _assert_matches_oracle(got1, oracle)
    _assert_matches_oracle(got4, oracle)

    # Byte-identical outputs across thread and bucket counts (AGENTS.md §4.3)
    assert (out1 / "saturation.json").read_bytes() == (out4 / "saturation.json").read_bytes()
    assert (out1 / "saturation.txt").read_bytes() == (out4 / "saturation.txt").read_bytes()


def test_a33_marginal_family_yield_decomposition_and_earlier_prefix(tmp_path):
    rs = _load_module("release_saturation", SAT_SCRIPT)
    root = tmp_path / "root"

    # Planted clusters A..E across S0, S1, S2:
    # A: S0=0, S1=+3 (brand_new_ge3 at S1), S2=+2 (promoted_to_ge5 at S2)
    # B: S0=1, S1=+2 (promoted_to_ge3 at S1), S2=0
    # C: S0=2, S1=+3 (promoted_to_ge3 and promoted_to_ge5 at S1), S2=0
    # D: S0=0, S1=0, S2=+5 (brand_new_ge3 and brand_new_ge5 at S2)
    # E: S0=0, S1=+1, S2=+2 (promoted_to_ge3 at S2 step-by-step; brand_new_ge3 vs earlier S0!)
    names = ("A", "B", "C", "D", "E")
    clusters = {_uid(f"p_{nm}"): _uid(f"c_{nm}") for nm in names}
    tok_map = {cid: token_id("nr50", 0, i) for i, cid in enumerate(sorted(clusters.values()))}
    by_stage = {
        "S0": {"A": 0, "B": 1, "C": 2, "D": 0, "E": 0},
        "S1": {"A": 3, "B": 2, "C": 3, "D": 0, "E": 1},
        "S2": {"A": 2, "B": 0, "C": 0, "D": 5, "E": 2},
    }

    g_seq = 0
    stage_defs = []
    for st in ("S0", "S1", "S2"):
        contigs = []
        for nm, n_gen in by_stage[st].items():
            pu = _uid(f"p_{nm}")
            for _ in range(n_gen):
                g_seq += 1
                contigs.append({
                    "contig_uid": _uid(f"cu_{st}_{g_seq}"),
                    "genome_uid": _uid(f"g_{g_seq}"),
                    "contig_len": 10_000,
                    "is_circular": False,
                    "slots": [{"start": 200, "end": 800, "strand": 1, "protein_uid": pu, "dist_to_contig_end_bp": 200}],
                })
        stage_defs.append({"stage": st, "contigs": contigs})

    _build_synthetic_release(root, "r26.11.2", stage_defs, clusters, token_map=tok_map)
    (root / "manifests" / "release=r26.11.1.json").write_text(
        json.dumps({"manifest_version": 3, "release_id": "r26.11.1", "stages": ["S0"], "files": [], "artifacts": []}),
        encoding="utf-8",
    )

    got = rs.measure_release(
        root=root, release="r26.11.2", earlier="r26.11.1", out=tmp_path / "out", tmp=tmp_path / "tmp"
    )
    oracle = _oracle_saturation(stage_defs, clusters, tok_map, earlier_stages=["S0"])
    _assert_matches_oracle(got, oracle)

    m_s1 = got["stage_curve"][1]["marginal"]
    assert (m_s1["brand_new_ge3"], m_s1["promoted_to_ge3"], m_s1["net_new_ge3"]) == (1, 2, 3)  # A brand_new; B,C promoted
    assert (m_s1["brand_new_ge5"], m_s1["promoted_to_ge5"], m_s1["net_new_ge5"]) == (0, 1, 1)  # C promoted

    m_s2 = got["stage_curve"][2]["marginal"]
    assert (m_s2["brand_new_ge3"], m_s2["promoted_to_ge3"], m_s2["net_new_ge3"]) == (1, 1, 2)  # D brand_new; E promoted
    assert (m_s2["brand_new_ge5"], m_s2["promoted_to_ge5"], m_s2["net_new_ge5"]) == (1, 1, 2)  # D brand_new; A promoted

    # Relative to r26.11.1 (S0 only): A, D, and E all had 0 genomes in S0 -> brand_new_ge3 == 3, promoted_to_ge3 == 2 (B, C)
    se = got["since_earlier"]
    assert (se["brand_new_ge3"], se["promoted_to_ge3"], se["net_new_ge3"]) == (3, 2, 5)
    assert (se["brand_new_ge5"], se["promoted_to_ge5"], se["net_new_ge5"]) == (2, 1, 3)  # A and D brand_new vs S0; C promoted


def test_a33_syntenic_pairs_strand_invariance_circular_wrap_and_null_skip(tmp_path):
    rs = _load_module("release_saturation", SAT_SCRIPT)

    # Use large uint64 local ids (> 2^53) so any float64 conversion on nullable arrays would corrupt tokens!
    p1, p2, p3, p4 = _uid("p1"), _uid("p2"), _uid("p3"), _uid("p4")
    c1, c2, c3, c4 = _uid("c1"), _uid("c2"), _uid("c3"), _uid("c4")
    clusters = {p1: c1, p2: c2, p3: c3, p4: c4}
    tok_map = {
        c1: token_id("nr50", 0, (1 << 54) + 1),
        c2: token_id("nr50", 0, (1 << 54) + 2),
        c3: token_id("nr50", 0, (1 << 54) + 3),
        c4: token_id("nr50", 0, (1 << 54) + 4),
    }

    def _rc_contig(c: dict) -> dict:
        L = c["contig_len"]
        rc_slots = []
        for s in c["slots"]:
            ns = dict(s)
            ns["start"] = L - int(s["end"])
            ns["end"] = L - int(s["start"])
            ns["strand"] = -int(s["strand"])
            rc_slots.append(ns)
        out = dict(c)
        out["slots"] = rc_slots
        return out

    fwd_contigs = [
        # Linear contig with interspersed non-CDS slot: (p1, p2) adjacent, then ncrna blocks (p2, p3), then (p3, p4) adjacent
        {
            "contig_uid": _uid("lin_1"),
            "genome_uid": _uid("g1"),
            "contig_len": 10_000,
            "is_circular": False,
            "slots": [
                {"kind": "cds", "start": 200, "end": 700, "strand": 1, "protein_uid": p1, "dist_to_contig_end_bp": 200},
                {"kind": "cds", "start": 800, "end": 1300, "strand": -1, "protein_uid": p2, "dist_to_contig_end_bp": 800},
                {"kind": "ncrna", "start": 1500, "end": 1600, "strand": 1},
                {"kind": "cds", "start": 2000, "end": 2500, "strand": 1, "protein_uid": p3, "dist_to_contig_end_bp": 2000},
                {"kind": "cds", "start": 2700, "end": 3200, "strand": 1, "protein_uid": p4, "dist_to_contig_end_bp": 2700},
            ],
        },
        # Circular 2-slot contig: wraps (1, 0) as well as (0, 1) -> 2 distinct junctions when same strand & distinct tokens
        {
            "contig_uid": _uid("circ_2"),
            "genome_uid": _uid("g2"),
            "contig_len": 6_000,
            "is_circular": True,
            "slots": [
                {"kind": "cds", "start": 300, "end": 900, "strand": 1, "protein_uid": p1, "dist_to_contig_end_bp": 6000},
                {"kind": "cds", "start": 2000, "end": 2600, "strand": 1, "protein_uid": p3, "dist_to_contig_end_bp": 6000},
            ],
        },
        # Circular 1-slot contig: n_slots < 2 -> must NOT self-wrap!
        {
            "contig_uid": _uid("circ_1"),
            "genome_uid": _uid("g3"),
            "contig_len": 3_000,
            "is_circular": True,
            "slots": [
                {"kind": "cds", "start": 300, "end": 900, "strand": 1, "protein_uid": p4, "dist_to_contig_end_bp": 3000},
            ],
        },
    ]
    rc_contigs = [_rc_contig(c) for c in fwd_contigs]

    stages_fwd = [{"stage": "S0", "contigs": fwd_contigs}]
    stages_rc = [{"stage": "S0", "contigs": rc_contigs}]

    _build_synthetic_release(tmp_path / "fwd", "r26.11.1", stages_fwd, clusters, token_map=tok_map)
    _build_synthetic_release(tmp_path / "rc", "r26.11.1", stages_rc, clusters, token_map=tok_map)

    got_fwd = rs.measure_release(root=tmp_path / "fwd", release="r26.11.1", out=tmp_path / "of", tmp=tmp_path / "tf")
    got_rc = rs.measure_release(root=tmp_path / "rc", release="r26.11.1", out=tmp_path / "or", tmp=tmp_path / "tr")

    oracle_fwd = _oracle_saturation(stages_fwd, clusters, tok_map)
    _assert_matches_oracle(got_fwd, oracle_fwd)
    # Forward vs reverse-complemented contigs yield identical saturation summaries!
    assert got_fwd == got_rc
    # 2 pairs from lin_1 ((p1+, p2-) and (p3+, p4+)) + 2 pairs from circ_2 ((p1+, p3+) and wrap (p3+, p1+)) + 0 from circ_1 = 4
    assert got_fwd["totals"]["syntenic_pairs"]["pairs_ge1"] == 4


def test_a33_train_check_with_and_without_saturation_summary(tmp_path):
    tc = _load_module("train_check", CHK_SCRIPT)
    root = tmp_path / "root"
    (root / "manifests").mkdir(parents=True)
    m1 = '{"release":"r26.11.1","stages":["S0"],"files":[],"artifacts":[]}'
    m2 = '{"release":"r26.11.2","stages":["S0","S1"],"files":[],"artifacts":[]}'
    (root / "manifests" / "release=r26.11.1.json").write_text(m1, encoding="utf-8")
    (root / "manifests" / "release=r26.11.2.json").write_text(m2, encoding="utf-8")

    before = tmp_path / "before.tsv"
    after = tmp_path / "after.tsv"
    before.write_text(f"manifests/release=r26.11.1.json\t{len(m1)}\t1700000000.0\n", encoding="utf-8")
    after.write_text(
        f"manifests/release=r26.11.1.json\t{len(m1)}\t1700000000.0\n"
        f"manifests/release=r26.11.2.json\t{len(m2)}\t1700000100.0\n",
        encoding="utf-8",
    )

    out_plain = tmp_path / "out_plain"
    assert tc.main([
        "--root", str(root), "--before", str(before), "--after", str(after),
        "--release", "r26.11.2", "--stages", "S0,S1", "--earlier", "r26.11.1", "--out", str(out_plain),
    ]) == 0
    ann_plain = (out_plain / "announcement.md").read_text(encoding="utf-8")
    assert "Saturation" not in ann_plain

    sat_path = tmp_path / "saturation.json"
    sat_doc = {
        "release": "r26.11.2",
        "earlier": "r26.11.1",
        "stages": ["S0", "S1"],
        "totals": {
            "stored_bases": 15_000_000_000,
            "stored_gbp": 15.0,
            "complete_orf": {
                "s_obs": 1_200_000, "f1": 600_000, "f2": 200_000, "ge3": 400_000, "ge5": 180_000,
                "n_incidences": 3_000_000, "c_hat": 0.80, "chao1": 2_099_995.5,
            },
            "syntenic_pairs": {
                "pairs_ge1": 2_500_000, "pairs_f1": 1_500_000, "pairs_f2": 400_000,
                "pairs_ge3": 600_000, "pairs_ge5": 250_000, "n_incidences": 5_000_000,
                "c_hat_pairs": 0.70, "chao1_pairs": 5_312_491.0,
            },
        },
        "since_earlier": {
            "earlier": "r26.11.1",
            "new_stages": ["S1"],
            "new_stored_bases": 5_000_000_000,
            "new_stored_gbp": 5.0,
            "net_new_ge3": 100_000,
            "brand_new_ge3": 65_000,
            "promoted_to_ge3": 35_000,
            "net_new_ge3_per_gbp": 20_000.0,
            "net_new_pairs_ge3": 150_000,
            "brand_new_pairs_ge3": 90_000,
            "promoted_to_pairs_ge3": 60_000,
            "net_new_pairs_ge3_per_gbp": 30_000.0,
        },
    }
    sat_path.write_text(json.dumps(sat_doc, indent=2), encoding="utf-8")

    out_sat = tmp_path / "out_sat"
    assert tc.main([
        "--root", str(root), "--before", str(before), "--after", str(after),
        "--release", "r26.11.2", "--stages", "S0,S1", "--earlier", "r26.11.1",
        "--saturation-summary", str(sat_path), "--out", str(out_sat),
    ]) == 0
    ann_sat = (out_sat / "announcement.md").read_text(encoding="utf-8")
    assert "- **Saturation (`release_saturation.py`):**" in ann_sat
    assert "400,000" in ann_sat
    assert "0.8000" in ann_sat or "0.8" in ann_sat
    assert "65,000" in ann_sat and "35,000" in ann_sat


def test_a33_history_render_and_minimal_path_subprocess(tmp_path):
    rs = _load_module("release_saturation", SAT_SCRIPT)
    root = tmp_path / "root"
    p1, p2, p3 = _uid("hp1"), _uid("hp2"), _uid("hp3")
    c1, c2 = _uid("hc1"), _uid("hc2")
    clusters = {p1: c1, p2: c1, p3: c2}
    tok_map = {c1: token_id("nr50", 0, 0), c2: token_id("nr50", 0, 1)}
    stages = [
        {
            "stage": "S1",
            "contigs": [
                {"contig_uid": _uid("hcu1"), "genome_uid": _uid("hg1"), "contig_len": 10_000, "is_circular": False,
                 "slots": [{"start": 100, "end": 400, "strand": 1, "protein_uid": p1, "dist_to_contig_end_bp": 100},
                           {"start": 500, "end": 800, "strand": 1, "protein_uid": p3, "dist_to_contig_end_bp": 500}]},
                {"contig_uid": _uid("hcu2"), "genome_uid": _uid("hg2"), "contig_len": 10_000, "is_circular": False,
                 "slots": [{"start": 100, "end": 400, "strand": 1, "protein_uid": p1, "dist_to_contig_end_bp": 100},
                           {"start": 500, "end": 800, "strand": 1, "protein_uid": p3, "dist_to_contig_end_bp": 500}]},
            ],
        },
        {
            "stage": "G1",
            "contigs": [
                {"contig_uid": _uid("hcu3"), "genome_uid": _uid("hg3"), "contig_len": 15_000, "is_circular": True,
                 "slots": [{"start": 100, "end": 400, "strand": 1, "protein_uid": p2, "dist_to_contig_end_bp": 15_000},
                           {"start": 500, "end": 800, "strand": 1, "protein_uid": p3, "dist_to_contig_end_bp": 15_000}]},
            ],
        },
    ]
    _build_synthetic_release(root, "r26.11.3", stages, clusters, token_map=tok_map)

    # Run `measure` via subprocess with minimal PATH=/usr/bin:/bin (hiding pod native tools)
    out_dir = tmp_path / "sat_out"
    stats_file = tmp_path / "stats.json"
    env_min = {
        "PATH": "/usr/bin:/bin",
        "PYTHONPATH": str(REPO / "src"),
        "HOME": str(tmp_path),
    }
    proc = subprocess.run(
        [sys.executable, str(SAT_SCRIPT), "measure", "--root", str(root), "--release", "r26.11.3",
         "--out", str(out_dir), "--tmp", str(tmp_path / "sat_tmp"), "--threads", "2", "--stats-out", str(stats_file)],
        env=env_min, capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert stats_file.is_file()
    stats = json.loads(stats_file.read_text(encoding="utf-8"))
    assert stats["wall_s"] >= 0.0 and stats["vmhwm_bytes"] >= 0

    # Build and check history JSON
    hist_path = tmp_path / "history.json"
    assert rs.main([
        "history", "--data", str(hist_path),
        "--add", str(out_dir / "saturation.json"),
        "--host", "wilkcomp02", "--date", "2026-11-27",
        "--command", "release_saturation.py measure --release r26.11.3",
        "--set-collector",
    ]) == 0
    assert rs.main(["history", "--data", str(hist_path), "--check"]) == 0

    # Tampering with a derived value or stripping a label must fail --check
    bad = json.loads(hist_path.read_text(encoding="utf-8"))
    bad["releases"][0]["totals"]["complete_orf"]["c_hat"]["value"] = 0.123456
    bad_path = tmp_path / "bad_history.json"
    bad_path.write_text(json.dumps(bad), encoding="utf-8")
    assert rs.main(["history", "--data", str(bad_path), "--check"]) == 1

    # Render deterministic PNG and SVG (no <dc:date>, no prose em/en-dashes)
    pytest.importorskip("matplotlib")
    fig_prefix = tmp_path / "fig" / "sat"
    assert rs.main(["render", "--data", str(hist_path), "--out", str(fig_prefix)]) == 0
    png, svg = Path(f"{fig_prefix}.png"), Path(f"{fig_prefix}.svg")
    assert png.is_file() and png.stat().st_size > 20_000
    svg_text = svg.read_text(encoding="utf-8")
    assert "<dc:date>" not in svg_text
    assert "—" not in svg_text and "–" not in svg_text
    fig_prefix2 = tmp_path / "fig2" / "sat"
    assert rs.main(["render", "--data", str(hist_path), "--out", str(fig_prefix2)]) == 0
    assert Path(f"{fig_prefix2}.svg").read_text(encoding="utf-8") == svg_text

    # If the repository's release_saturation_history.json is present, verify it passes --check too
    if HISTORY_JSON.is_file():
        assert rs.main(["history", "--data", str(HISTORY_JSON), "--check"]) == 0


def test_a33_release_train_invokes_release_saturation_and_passes_summary_to_train_check(tmp_path):
    """Verify release_train.sh step 8 runs release_saturation.py (clamped <= 4 threads), copies saturation.json,
    cleans $WORK/tmp/saturation, and passes --saturation-summary to train_check.py."""
    sys.path.insert(0, str(REPO / "tests" / "r1"))
    import test_c46_release_train as tc46

    drv, capture, work, env = tc46._a25_make_train_world(tmp_path, mem_scope="750G", max_day=40, max_night=80)
    env["BIOBOX_TRAIN_NOW"] = "2026-11-24T22:00:00"  # night -> eff_ex_threads=24, sat_threads must clamp to 4!

    fake = tmp_path / "fakebin"
    # Extend singularity stub to capture release_saturation.py invocation and write a synthetic saturation.json
    sing = fake / "singularity"
    sing_text = sing.read_text(encoding="utf-8")
    inject = (
        '  *release_saturation.py*measure*)\n'
        '    echo "release_saturation args=$*" >> "$CAPTURE_LOG"\n'
        '    out=""; prev=""\n'
        '    for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done\n'
        '    mkdir -p "$out"\n'
        '    printf \'{"release":"r26.11.2","totals":{"complete_orf":{"s_obs":10,"ge3":3,"ge5":1,"c_hat":0.75,"chao1":12.0},"syntenic_pairs":{"pairs_ge3":2,"c_hat_pairs":0.5}}}\' > "$out/saturation.json"\n'
        '    exit 0 ;;\n'
    )
    sing.write_text(sing_text.replace('case "$*" in\n', 'case "$*" in\n' + inject, 1), encoding="utf-8")

    # Extend python3 stub to log train_check.py args
    py3 = fake / "python3"
    py3_text = py3.read_text(encoding="utf-8")
    py3.write_text(
        py3_text.replace(
            "*train_check.py*)\n",
            '*train_check.py*)\n    echo "train_check args=$*" >> "$CAPTURE_LOG"\n',
            1,
        ),
        encoding="utf-8",
    )

    got = subprocess.run(["bash", str(drv), "run", "--no-live-cpu"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    lines = capture.read_text(encoding="utf-8").splitlines()
    sat_lines = [ln for ln in lines if ln.startswith("release_saturation ")]
    assert len(sat_lines) == 1, lines
    assert "--release r26.11.2" in sat_lines[0]
    assert "--earlier r26.11.1" in sat_lines[0]
    assert "--threads 4" in sat_lines[0]  # clamped from 24 to 4!
    assert (work / "saturation.json").is_file()
    assert not (work / "tmp" / "saturation").exists()  # cleaned!

    chk_lines = [ln for ln in lines if ln.startswith("train_check ")]
    assert len(chk_lines) == 1
    assert f"--saturation-summary {work}/saturation/saturation.json" in chk_lines[0]


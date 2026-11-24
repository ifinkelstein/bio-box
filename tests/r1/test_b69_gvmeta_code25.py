"""Ring R1: Rescue `code 25` (`UGA -> Gly`) CPR contigs in `--genecall-mode gv-meta-c25` (card B-69, DEC-71).

Tests:
1. `GENECALL_MODES == ("default", "gv-meta", "gv-meta-c25")`, `mode_of`, and pricing:
   `gv-meta-c25` uses `gv-meta` rates without altering `CALL_S_PER_MBP` or `GENES_PER_MBP`,
   keeping `"default"` and `"gv-meta"` plans and outputs 100% byte-identical.
2. Constructed `>= 20 kb` `UGA -> Gly` (`code 25`) CPR contig carrying an RT ORF with conserved
   palm/fingers glycines encoded by `TGA`:
   - `gv-meta` calls it under `genetic_code == 4`, translating internal `TGA` as `W`;
   - `gv-meta-c25` rescues it to `genetic_code == 25`, translating internal `TGA` as `G`
     (recovering intact RT palm glycine motifs and matching an independent `translation_table=25`
     re-translation oracle 1:1).
3. True `code 4` (*Mycoplasma genitalium* `TGA -> Trp`), `code 15` (`TAG -> Gln` jumbophage),
   and `code 11` standard prokaryotic contigs are untouched by `gv-meta-c25` (identical
   coordinates, genetic codes, and protein sequences to `gv-meta`).
4. Determinism & equivalence: `gv-meta-c25` 1-process (`workers=1` and `workers=2`) and batched
   (`plan` -> `publish_plan` -> `run_batch` -> `merge`) produce byte-identical Parquet outputs
   and identical result dicts.
5. `stage.py` / `cli.py` `--genecall-mode gv-meta-c25` and `stage_est.py` / `stage_driver.sh`
   `GENECALL_MODE=gv-meta-c25` parsing, estimation, and rate guard.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import random
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pyrodigal
import pyrodigal_gv
import pytest

from biobox import cli, fasta, ids, profile
from biobox import stage as bs
from biobox.l0seq import inputs as l0inputs
from biobox.l0seq import store as l0store
from biobox.l0seq.inputs import GenomeInput
from biobox.l0seq.store import SourceObject
from biobox.l1 import genecall

pytestmark = pytest.mark.covers("B-69")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"

# Canonical RT palm/fingers motif tokens with conserved glycines (G) that will be encoded by TGA
_RT_MOTIF_AA = (
    "MKTQILVDEK"
    "GLP"        # fingers glycine
    "VKRIASTNDLFEYH"
    "QGF"        # palm motif B glycine
    "PLSPILANVFTKDR"
    "YADD"       # palm catalytic motif C
    "LVIAKSEWTP" # real TGG Trp alongside TGA Gly
    "GFLG"       # conserved double-glycine motif
    "HTRNQVYICW"
)


def _load_stage_est():
    spec = importlib.util.spec_from_file_location("stage_est_b69", BSP / "stage_est.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage_est_b69"] = mod
    spec.loader.exec_module(mod)
    return mod


def _encode_aa_code25(aa: str, seed: int = 69) -> str:
    """Encode an amino acid sequence as DNA using code 25 (`G` -> `TGA`, `W` -> `TGG`)."""
    table = {
        "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
        "G": "TGA", "H": "CAT", "I": "ATT", "K": "AAA", "L": "TTA",
        "M": "ATG", "N": "AAT", "P": "CCA", "Q": "CAA", "R": "AGA",
        "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
    }
    rng = random.Random(seed)
    codons = [table[ch] for ch in aa]
    # Pad to ~400 aa so Prodigal calls it cleanly as a full-length RT ORF
    filler = [
        "GCT", "GAT", "GAA", "TTT", "CAT", "ATT", "AAA", "TTA", "AAT",
        "CCA", "CAA", "AGA", "TCT", "ACT", "GTT", "TAT", "TGA", "GGT",
    ]
    for _ in range(320):
        codons.append(rng.choice(filler))
    return "".join(codons) + "TAA"


def _make_code25_cpr_contig() -> tuple[bytes, str]:
    """Construct a >= 25 kb `UGA -> Gly` (`code 25`) CPR contig from a 30 kb slice of T4
    (AT-rich ~35% GC, similar to Absconditabacteria/Gracilibacteria) by recoding ~65% of
    internal in-frame `Gly` codons (`GGN` on +1, `NCC` on -1) to `TGA` (`TCA` on -1) while
    preserving all `TGG` (`Trp`) codons, plus an embedded RT ORF with conserved palm glycines
    encoded by `TGA` (`GLP`, `QGF`, `YADD`, `GFLG`)."""
    t4_path = REPO / "tests" / "fixtures" / "micro_s0" / "NC_000866.4.fna.gz"
    raw = fasta.read(t4_path)[0].seq.upper().decode("ascii")[:30_000]
    buf = list(raw)
    finder = pyrodigal.GeneFinder(meta=True)
    last_end = 0
    gly_idx = 0
    for g in finder.find_genes(raw):
        start, end = g.begin - 1, g.end
        if start < last_end or (end - start) < 240:
            continue
        last_end = end
        if g.strand == 1:
            for pos in range(start + 3, end - 3, 3):
                if raw[pos : pos + 2] == "GG":
                    gly_idx += 1
                    if gly_idx % 3 != 0:  # recode 2/3 of Gly codons to TGA (~4-6% of codons)
                        buf[pos : pos + 3] = list("TGA")
        else:
            for pos in range(end - 6, start, -3):
                if raw[pos + 1 : pos + 3] == "CC":
                    gly_idx += 1
                    if gly_idx % 3 != 0:
                        buf[pos : pos + 3] = list("TCA")
    rt_orf = "AAAGGAGGTTATAA" + _encode_aa_code25(_RT_MOTIF_AA) + "TTTAAATATAAACC"
    seq = ("".join(buf) + "TTATTTAA" + rt_orf + "TTATTTAA").encode("ascii")
    assert len(seq) >= 25_000
    return seq, _RT_MOTIF_AA


def _make_code4_mycoplasma_contig() -> bytes:
    """Real >= 25 kb `code 4` (`TGA -> Trp`) contig from Mycoplasma genitalium (micro_s0)."""
    mg_path = REPO / "tests" / "fixtures" / "micro_s0" / "NC_000908.2.fna.gz"
    seq = fasta.read(mg_path)[0].seq.upper()[:30_000]
    assert len(seq) == 30_000
    return seq


def _make_code15_jumbophage_contig() -> bytes:
    """Construct a >= 35 kb TAG-recoded (`code 15`) jumbophage contig from T4 (micro_s0)."""
    t4_path = REPO / "tests" / "fixtures" / "micro_s0" / "NC_000866.4.fna.gz"
    raw = fasta.read(t4_path)[0].seq.upper().decode("ascii")[40_000:80_000]
    buf = list(raw)
    finder = pyrodigal.GeneFinder(meta=True)
    last_end = 0
    for g in finder.find_genes(raw):
        start, end = g.begin - 1, g.end
        if start < last_end or (end - start) < 300:
            continue
        last_end = end
        if g.strand == 1:
            for pos in range(start + 3, end - 3, 3):
                if raw[pos : pos + 3] in ("CAG", "CAA"):
                    buf[pos : pos + 3] = list("TAG")
        else:
            for pos in range(end - 6, start, -3):
                if raw[pos : pos + 3] in ("CTG", "TTG"):
                    buf[pos : pos + 3] = list("CTA")
    return "".join(buf).encode("ascii")


def _make_code11_contig() -> bytes:
    """Real standard `code 11` contig from phage lambda (micro_s0)."""
    lam_path = REPO / "tests" / "fixtures" / "micro_s0" / "NC_001416.1.fna.gz"
    return fasta.read(lam_path)[0].seq.upper()[:15_000]


def _stage_fixture(tmp_path: Path, stage: str = "C25"):
    root = tmp_path / "root"
    inp = tmp_path / "inputs"
    inp.mkdir(parents=True)

    cpr_seq, rt_motif = _make_code25_cpr_contig()
    mg_seq = _make_code4_mycoplasma_contig()
    c15_seq = _make_code15_jumbophage_contig()
    c11_seq = _make_code11_contig()

    fna = inp / "stage.fna"
    with open(fna, "wb") as fh:
        for name, seq in (
            ("cpr_code25", cpr_seq),
            ("mycoplasma_code4", mg_seq),
            ("jumbo_code15", c15_seq),
            ("standard_code11", c11_seq),
        ):
            fh.write(b"\n".join(fasta.format_record(name, seq, 70)) + b"\n")

    cu_cpr = ids.contig_uid(cpr_seq)
    cu_mg = ids.contig_uid(mg_seq)
    cu_c15 = ids.contig_uid(c15_seq)
    cu_c11 = ids.contig_uid(c11_seq)

    # All staged as unbinned metagenome contigs (`kind="prokaryote"`, `genetic_code=11`)
    genomes = [
        GenomeInput((cu_cpr,), 11, kind="prokaryote"),
        GenomeInput((cu_mg,), 11, kind="prokaryote"),
        GenomeInput((cu_c15,), 11, kind="prokaryote"),
        GenomeInput((cu_c11,), 11, kind="prokaryote"),
    ]
    objs = [SourceObject(fna, "test", "1", "CC0", True)]
    l0inputs.write_genome_input(root, stage, genomes)
    l0store.ingest(objs, root, stage, tmp=tmp_path)

    (inp / "genomes.json").write_text(json.dumps([
        {"contig_uids": [c.hex() for c in g.contig_uids], "genetic_code": g.genetic_code,
         "kind": g.kind, "circular": [], "assembly_quality": "contig"}
        for g in genomes
    ], indent=1))
    (inp / "objects.json").write_text(json.dumps([
        {"path": str(fna), "source_name": "test", "source_release": "1",
         "license": "CC0", "redistributable": True}
    ], indent=1))
    return {
        "root": root, "inp": inp, "stage": stage,
        "cpr_seq": cpr_seq, "rt_motif": rt_motif,
        "mg_seq": mg_seq, "c15_seq": c15_seq, "c11_seq": c11_seq,
        "cu_cpr": cu_cpr, "cu_mg": cu_mg, "cu_c15": cu_c15, "cu_c11": cu_c11,
        "genomes": genomes,
    }


def _tree_hashes(base: Path) -> dict[str, str]:
    return {
        str(p.relative_to(base)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(base.rglob("*"))
        if p.is_file()
    }


def _read_calls_by_contig(out: Path, stage: str) -> dict[bytes, list[tuple]]:
    fp = pq.read_table(
        out / f"l1/feature_protein/gene_call_version={genecall.GENE_CALL_VERSION}/stage={stage}"
    ).to_pylist()
    pseq = {
        r["protein_uid"]: r["seq"]
        for r in pq.read_table(
            out / f"l0seq/protein_seq/gene_call_version={genecall.GENE_CALL_VERSION}/stage={stage}"
        ).to_pylist()
    }
    by_cu: dict[bytes, list[tuple]] = {}
    for r in fp:
        by_cu.setdefault(r["contig_uid"], []).append(
            (r["start"], r["end"], r["strand"], r["genetic_code"], pseq[r["protein_uid"]])
        )
    for v in by_cu.values():
        v.sort()
    return by_cu


def test_b69_modes_and_pricing():
    """`GENECALL_MODES` includes `gv-meta-c25`, prices with `gv-meta` rates, and keeps
    `CALL_S_PER_MBP` / `GENES_PER_MBP` untouched."""
    assert genecall.GENECALL_MODES == ("default", "gv-meta", "gv-meta-c25")
    for kind in ("prokaryote", "plasmid", "virus"):
        for bases in (50_000, 120_000):
            assert genecall.mode_of(kind, bases, genecall_mode="gv-meta-c25") == "gv-meta"
            for circ in (False, True):
                assert genecall.genome_cost_us(kind, [(bases, circ)], genecall_mode="gv-meta-c25") == (
                    genecall.genome_cost_us(kind, [(bases, circ)], genecall_mode="gv-meta")
                )
            assert genecall.genome_genes(kind, bases, genecall_mode="gv-meta-c25") == (
                genecall.genome_genes(kind, bases, genecall_mode="gv-meta")
            )
    assert set(genecall.CALL_S_PER_MBP) == {
        ("single", False), ("single", True),
        ("gv-meta", False), ("gv-meta", True),
        ("meta", False), ("meta", True),
    }


def test_b69_code25_llr_discriminator_calibration():
    """Trp-only binomial LLR separates empirical CPR code-25 profiles from code-4 profiles."""
    # Empirical CPR code-25 profiles (20 kb, ~6,000 codons; planning/measurements/2026-11-23_globdb_art_phyla.md §2)
    assert genecall._code25_llr(6000, n_tga=123, n_tgg=60) >= genecall.CODE25_MIN_LLR   # Absconditabacterales (6.83 TGA/kb)
    assert genecall._code25_llr(6000, n_tga=356, n_tgg=60) >= genecall.CODE25_MIN_LLR   # Gracilibacteria (19.79 TGA/kb)
    assert genecall._code25_llr(6000, n_tga=40, n_tgg=63) >= genecall.CODE25_MIN_LLR    # Minisyncoccia (2.2 TGA/kb)

    # True code-4 profiles (Mycoplasma genitalium + high-Trp code-4 phage)
    assert genecall._code25_llr(6000, n_tga=38, n_tgg=20) < 0.0   # M. genitalium (~0.97% total Trp)
    assert genecall._code25_llr(6000, n_tga=60, n_tgg=24) < 0.0   # High-Trp code-4 phage (~1.40% total Trp)


def test_b69_rescues_code25_cpr_and_leaves_code4_code15_code11_unchanged(tmp_path: Path):
    """On a mixed stage with code-25 CPR (with RT palm glycines at TGA), true code-4 Mycoplasma,
    code-15 jumbophage, and code-11 contigs:
    - `gv-meta` calls CPR as code 4 (`TGA -> W`, corrupting RT palm glycines);
    - `gv-meta-c25` calls CPR as code 25 (`TGA -> G`, recovering intact RT palm glycines and
      matching an independent `translation_table=25` oracle 1:1);
    - code 4, code 15, and code 11 contigs are 100% identical between `gv-meta` and `gv-meta-c25`;
    - `workers=1`, `workers=2`, and batched (`plan`/`run_batch`/`merge`) are byte-identical."""
    fx = _stage_fixture(tmp_path)
    root, st = fx["root"], fx["stage"]

    # Independent oracle for the CPR contig: ViralGeneFinder(meta=True) with table-25 translation
    oracle_vf = pyrodigal_gv.ViralGeneFinder(meta=True)
    oracle_cpr = sorted(
        (g.begin - 1, g.end, g.strand,
         25 if g.translation_table == 4 else g.translation_table,
         g.translate(translation_table=25 if g.translation_table == 4 else None, include_stop=False))
        for g in oracle_vf.find_genes(fx["cpr_seq"].decode("ascii"))
    )
    assert any(g[3] == 25 for g in oracle_cpr)

    # Run gv-meta (must stay untouched: CPR called as code 4 with W at TGA)
    out_gv = tmp_path / "out_gv"
    res_gv = genecall.run(root, st, out=out_gv, tmp=tmp_path, workers=1, genecall_mode="gv-meta")
    assert res_gv["modes"] == ["gv-meta"]
    gv_calls = _read_calls_by_contig(out_gv, st)

    cpr_gv = gv_calls[fx["cu_cpr"]]
    assert any(g[3] == 4 for g in cpr_gv), "gv-meta must call code-25 CPR contig under code 4"
    assert not any(g[3] == 25 for g in cpr_gv)
    corrupted_motif = fx["rt_motif"].replace("G", "W")
    assert any(corrupted_motif in g[4] for g in cpr_gv), "gv-meta must translate TGA palm glycines as W"
    assert not any(fx["rt_motif"] in g[4] for g in cpr_gv)

    # Run gv-meta-c25 (workers=1 and workers=2)
    out_c25_1 = tmp_path / "out_c25_1"
    out_c25_2 = tmp_path / "out_c25_2"
    res_c25_1 = genecall.run(root, st, out=out_c25_1, tmp=tmp_path, workers=1, genecall_mode="gv-meta-c25")
    res_c25_2 = genecall.run(root, st, out=out_c25_2, tmp=tmp_path, workers=2, genecall_mode="gv-meta-c25")
    assert res_c25_1["modes"] == ["gv-meta-c25"]
    assert res_c25_1 == res_c25_2
    assert _tree_hashes(out_c25_1) == _tree_hashes(out_c25_2)

    c25_calls = _read_calls_by_contig(out_c25_1, st)

    # 1. CPR contig rescued to code 25 with intact Gly RT palm motifs, matching oracle 1:1
    cpr_c25 = c25_calls[fx["cu_cpr"]]
    assert any(g[3] == 25 for g in cpr_c25)
    assert not any(g[3] == 4 for g in cpr_c25)
    assert any(fx["rt_motif"] in g[4] for g in cpr_c25), "gv-meta-c25 must recover intact Gly RT palm motifs"
    assert cpr_c25 == oracle_cpr

    # 2. True code-4 Mycoplasma contig stays code 4 and 100% identical to gv-meta
    assert any(g[3] == 4 for g in c25_calls[fx["cu_mg"]])
    assert not any(g[3] == 25 for g in c25_calls[fx["cu_mg"]])
    assert c25_calls[fx["cu_mg"]] == gv_calls[fx["cu_mg"]]

    # 3. Code-15 jumbophage contig stays code 15 and 100% identical to gv-meta
    assert any(g[3] == 15 for g in c25_calls[fx["cu_c15"]])
    assert c25_calls[fx["cu_c15"]] == gv_calls[fx["cu_c15"]]

    # 4. Code-11 standard contig stays code 11 and 100% identical to gv-meta
    assert all(g[3] == 11 for g in c25_calls[fx["cu_c11"]])
    assert c25_calls[fx["cu_c11"]] == gv_calls[fx["cu_c11"]]

    # 5. Batched plan -> publish_plan -> run_batch -> merge is byte-identical to 1-process gv-meta-c25
    p_c25 = genecall.plan(root, st, batch_cpu_min=240.0, batches=2, tmp=tmp_path, genecall_mode="gv-meta-c25")
    assert p_c25["genecall_mode"] == "gv-meta-c25"
    work = genecall.work_dir(tmp_path / "runs_tmp", root, st, p_c25)
    plan_path = genecall.write_plan(p_c25, work)
    genecall.publish_plan(root, st, plan_path, root, tmp=tmp_path)
    for j in range(len(p_c25["batches"])):
        rec = genecall.run_batch(root, st, plan_path, j, workers=2)
        assert rec["modes"] == ["gv-meta-c25"]
        genecall.publish_batch(root, st, plan_path, rec, root)
    out_batched = tmp_path / "out_batched"
    res_batched = genecall.merge(root, st, plan_path, out=out_batched, tmp=tmp_path)
    assert res_batched == res_c25_1
    assert _tree_hashes(out_batched) == _tree_hashes(out_c25_1)


def test_b69_cli_stage_and_stage_est_gvmeta_c25(tmp_path: Path):
    """`stage.Build`, `cli.main`, and `stage_est.py` accept `gv-meta-c25` and price it with
    `gv-meta` rates."""
    est = _load_stage_est()
    assert est.GENECALL_MODES == ("default", "gv-meta", "gv-meta-c25")

    prof = profile.Profile(
        max_threads_day=2, max_threads_night=2, mem_scope="4G", tmp=str(tmp_path / "ptmp"),
        durable=str(tmp_path / "dur"), bulk=str(tmp_path / "bulk"), nice_monolithic=0, nice_fanout=0,
        day_start_hour=8, day_end_hour=20, cpu_min_per_process=300, small_file_floor_bytes=1000,
    )
    b = bs.Build(root=tmp_path / "r", stage="S0", genomes=tmp_path / "g.json", objects=tmp_path / "o.json",
                 genecall_mode="gv-meta-c25")
    gc = next(s for s in bs.steps(prof, b) if s.name == "genecall")
    idx = gc.argv.index("--genecall-mode")
    assert gc.argv[idx + 1] == "gv-meta-c25"

    fx = _stage_fixture(tmp_path / "cli_fx")
    out_cli = tmp_path / "cli_out"
    cli_tmp = tmp_path / "cli_tmp"
    cli_tmp.mkdir()
    rc = cli.main(["stage", "genecall", "--root", str(fx["root"]), "--stage", fx["stage"],
                   "--genecall-mode", "gv-meta-c25", "--out", str(out_cli), "--tmp", str(cli_tmp)])
    assert rc == 0
    cli_calls = _read_calls_by_contig(out_cli, fx["stage"])
    assert any(g[3] == 25 for g in cli_calls[fx["cu_cpr"]])
    assert any(g[3] == 4 for g in cli_calls[fx["cu_mg"]])

    # stage_est.py conf + before + rates_doc
    inp = fx["inp"]
    sums = {
        n: hashlib.sha256((inp / n).read_bytes()).hexdigest()
        for n in ("genomes.json", "objects.json", "stage.fna")
    }
    conf_p = tmp_path / f"{fx['stage']}.conf"
    conf_p.write_text("\n".join([
        f"STAGE={fx['stage']}",
        "NOTE=planning/measurements/2026-11-20_stage_V1.md",
        f"INPUTS={inp}",
        *(f"SUM={h} {n}" for n, h in sums.items()),
        "DEPTH=called",
        "CRISPR=no",
        "GENECALL_MODE=gv-meta-c25",
    ]) + "\n")
    parsed = est.read_conf(conf_p, fx["stage"])
    assert parsed.genecall_mode == "gv-meta-c25"

    prof_path = tmp_path / "profile.yaml"
    prof_path.write_text(
        f"max_threads_day: 2\nmax_threads_night: 2\nmem_scope: 4G\ntmp: {tmp_path / 'pt'}\n"
        f"durable: {tmp_path / 'pd'}\nbulk: {tmp_path / 'pb'}\nnice_monolithic: 0\nnice_fanout: 0\n"
        "day_start_hour: 8\nday_end_hour: 20\ncpu_min_per_process: 300\nsmall_file_floor_bytes: 1000\n"
    )
    plan_p = tmp_path / "plan_c25.json"
    p = subprocess.run(
        [sys.executable, str(BSP / "stage_est.py"), "before", "--conf", str(conf_p),
         "--stage", fx["stage"], "--root", str(fx["root"]), "--profile", str(prof_path),
         "--inputs", str(inp), "--tmp", str(tmp_path / "bt_c25"), "--out-plan", str(plan_p)],
        capture_output=True, text=True, check=True,
    )
    doc = json.loads(p.stdout)
    assert doc["genecall_mode"] == "gv-meta-c25"
    assert json.loads(plan_p.read_text())["genecall_mode"] == "gv-meta-c25"
    want = genecall.plan(fx["root"], fx["stage"], batch_cpu_min=240.0, tmp=tmp_path, genecall_mode="gv-meta-c25")
    assert doc["genecall"]["estimate_cpu_min"] == want["estimate_cpu_min"]

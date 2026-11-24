"""Ring R1: `--genecall-mode gv-meta` for unbinned metagenome stages (card B-67, DEC-71).

Tests:
1. `genecall_mode="default"` (and omitted) leaves `mode_of`, `genome_cost_us`, `genome_genes`,
   `plan` (omitting `"genecall_mode"`), `run`, `run_batch` + `merge`, `stage.steps`, and
   `stage_est.py` 100% byte-identical.
2. On a constructed `>= 120 kb` TAG-recoded (genetic code 15) jumbophage contig staged as
   `GenomeInput(kind="prokaryote")` alongside a `< 10 kb` prokaryote contig:
   - `default` trains code-11 single mode on the `>= 120 kb` contig and fragments the TAG-recoded ORFs;
   - `gv-meta` calls intact full-length `genetic_code == 15` ORFs (`mode == "gv-meta"`) matching
     a direct `pyrodigal_gv.ViralGeneFinder(meta=True)` oracle 1:1 across 1-process (`workers=1`
     and `workers=2`) and batched (`plan` -> `publish_plan` -> `run_batch` -> `merge`) runs.
3. Pricing (`CALL_S_PER_MBP[("gv-meta", is_circular)]` without `META_LINEAR_LONG_FACTOR`,
   `GENES_PER_MBP["gv-meta"]`), `stage.py` / `cli.py` `--genecall-mode`, and `stage_est.py` /
   `stage_driver.sh` `GENECALL_MODE=gv-meta` parsing, validation, estimation, and rate guard.
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

pytestmark = pytest.mark.covers("B-67")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _load_stage_est():
    spec = importlib.util.spec_from_file_location("stage_est_b67", BSP / "stage_est.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage_est_b67"] = mod
    spec.loader.exec_module(mod)
    return mod


def _make_code15_jumbophage_contig() -> tuple[bytes, list[tuple[int, int]]]:
    """Construct a >= 120 kb TAG-recoded (code 15) jumbophage contig from the micro_s0 T4 genome
    (168,903 bp) by recoding internal in-frame Gln codons (`CAG`/`CAA` on +1, `CTG`/`TTG` on -1)
    to `TAG` (`CTA` on -1) inside T4's real CDSs.

    Authentic phage RBSs, start/stop codons, and dicodon usage are preserved: under code-11 single
    training (`default`), every internal `TAG` prematurely terminates the ORF; under `gv-meta`
    (`pyrodigal_gv.ViralGeneFinder(meta=True)`), the viral code-15 model recovers the intact
    full-length `genetic_code == 15` ORFs.
    """
    t4_path = REPO / "tests" / "fixtures" / "micro_s0" / "NC_000866.4.fna.gz"
    raw = fasta.read(t4_path)[0].seq.upper().decode("ascii")
    buf = list(raw)
    finder = pyrodigal.GeneFinder(meta=True)
    recoded_long: list[tuple[int, int]] = []
    last_end = 0
    for g in finder.find_genes(raw):
        start, end = g.begin - 1, g.end
        if start < last_end or (end - start) < 300:
            continue
        last_end = end
        n_mut = 0
        if g.strand == 1:
            for pos in range(start + 3, end - 3, 3):
                if raw[pos : pos + 3] in ("CAG", "CAA"):
                    buf[pos : pos + 3] = list("TAG")
                    n_mut += 1
        else:
            for pos in range(end - 6, start, -3):
                if raw[pos : pos + 3] in ("CTG", "TTG"):
                    buf[pos : pos + 3] = list("CTA")
                    n_mut += 1
        if n_mut >= 5 and (end - start) >= 1200:
            recoded_long.append((start, end))
    seq = "".join(buf).encode("ascii")
    assert len(seq) >= 120_000
    assert len(recoded_long) >= 20
    return seq, recoded_long


def _make_standard_contig(n_genes: int = 6, seed: int = 167) -> bytes:
    """Deterministic < 10 kb standard code-11 prokaryotic contig."""
    rng = random.Random(seed)
    sense = [
        "GCT", "GCC", "CGT", "AAT", "GAT", "CAA", "GAA", "GGT", "CAT", "ATT",
        "CTG", "AAA", "TTT", "CCG", "TCT", "ACT", "TGG", "TAT", "GTT", "GTG",
    ]
    parts = ["T" * 60]
    for _ in range(n_genes):
        parts.append("AAAGGAGGTTATAA")
        parts.append("ATG" + "".join(rng.choice(sense) for _ in range(300)) + "TAA")
        parts.append("TTTAAATATAAACC")
    seq = "".join(parts).encode("ascii")
    assert len(seq) < 10_000
    return seq


def _stage_fixture(tmp_path: Path, stage: str = "M1"):
    """Create a staged root with a >=120 kb code-15 contig (prokaryote), a <10 kb prokaryote
    contig, a circular plasmid, and a virus contig."""
    root = tmp_path / "root"
    inp = tmp_path / "inputs"
    inp.mkdir(parents=True)
    jumbo_seq, spans = _make_code15_jumbophage_contig()
    short_seq = _make_standard_contig(6, seed=167)
    plasmid_seq = _make_standard_contig(8, seed=267)
    virus_seq = _make_standard_contig(7, seed=367)

    fna = inp / "stage.fna"
    with open(fna, "wb") as fh:
        for name, seq in (
            ("jumbo_code15", jumbo_seq),
            ("short_prok", short_seq),
            ("plasmid_circ", plasmid_seq),
            ("virus_lin", virus_seq),
        ):
            fh.write(b"\n".join(fasta.format_record(name, seq, 70)) + b"\n")

    cu_jumbo = ids.contig_uid(jumbo_seq)
    cu_short = ids.contig_uid(short_seq)
    cu_plasmid = ids.contig_uid(plasmid_seq)
    cu_virus = ids.contig_uid(virus_seq)

    genomes = [
        GenomeInput((cu_jumbo,), 11, kind="prokaryote"),
        GenomeInput((cu_short,), 11, kind="prokaryote"),
        GenomeInput((cu_plasmid,), 11, kind="plasmid", circular=frozenset({cu_plasmid})),
        GenomeInput((cu_virus,), 11, kind="virus"),
    ]
    objs = [SourceObject(fna, "test", "1", "CC0", True)]
    l0inputs.write_genome_input(root, stage, genomes)
    l0store.ingest(objs, root, stage, tmp=tmp_path)

    (inp / "genomes.json").write_text(json.dumps([
        {"contig_uids": [c.hex() for c in g.contig_uids], "genetic_code": g.genetic_code,
         "kind": g.kind, "circular": [c.hex() for c in sorted(g.circular)], "assembly_quality": "complete"}
        for g in genomes
    ], indent=1))
    (inp / "objects.json").write_text(json.dumps([
        {"path": str(fna), "source_name": "test", "source_release": "1",
         "license": "CC0", "redistributable": True}
    ], indent=1))
    return {
        "root": root,
        "inp": inp,
        "fna": fna,
        "stage": stage,
        "jumbo_seq": jumbo_seq,
        "spans": spans,
        "cu_jumbo": cu_jumbo,
        "cu_short": cu_short,
        "cu_plasmid": cu_plasmid,
        "cu_virus": cu_virus,
        "genomes": genomes,
    }


def _tree_hashes(base: Path) -> dict[str, str]:
    return {
        str(p.relative_to(base)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(base.rglob("*"))
        if p.is_file()
    }


def test_b67_mode_of_and_pricing_default_vs_gvmeta():
    """`default` preserves policy p2; `gv-meta` forces `gv-meta` and prices with `gv-meta` rates
    (no `META_LINEAR_LONG_FACTOR`)."""
    assert genecall.GENECALL_MODES == ("default", "gv-meta", "gv-meta-c25")
    for kind, bases, want in (
        ("prokaryote", 50_000, "meta"),
        ("prokaryote", 120_000, "single"),
        ("plasmid", 120_000, "meta"),
        ("virus", 50_000, "gv-meta"),
        ("virus", 120_000, "gv-meta"),
    ):
        assert genecall.mode_of(kind, bases) == want
        assert genecall.mode_of(kind, bases, genecall_mode="default") == want
        assert genecall.mode_of(kind, bases, genecall_mode="gv-meta") == "gv-meta"

    with pytest.raises(ValueError, match="genecall_mode"):
        genecall.mode_of("prokaryote", 50_000, genecall_mode="bogus")

    # Pricing: gv-meta uses CALL_S_PER_MBP[("gv-meta", circular)] and GENES_PER_MBP["gv-meta"],
    # with NO META_LINEAR_LONG_FACTOR on >= 10 kb linear contigs.
    for circ in (False, True):
        n = 120_000
        want_us = int(-(-genecall.STEP_OVERHEAD * genecall.CALL_S_PER_MBP[("gv-meta", circ)]
                        * genecall.called_bases(n, circ) // 1))
        assert genecall.genome_cost_us("prokaryote", [(n, circ)], genecall_mode="gv-meta") == want_us
        assert genecall.genome_cost_us("plasmid", [(n, circ)], genecall_mode="gv-meta") == want_us
    assert genecall.genome_genes("prokaryote", 120_000, genecall_mode="gv-meta") == pytest.approx(
        genecall.GENES_PER_MBP["gv-meta"] * 120_000 / 1e6
    )
    # Default pricing unchanged
    assert genecall.genome_cost_us("prokaryote", [(120_000, False)]) == genecall.genome_cost_us(
        "prokaryote", [(120_000, False)], genecall_mode="default"
    )
    assert genecall.genome_genes("prokaryote", 120_000) == genecall.genome_genes(
        "prokaryote", 120_000, genecall_mode="default"
    )


def test_b67_default_mode_byte_identical_plan_and_outputs(tmp_path: Path):
    """Omitting `genecall_mode` and passing `genecall_mode='default'` produce byte-identical plans
    (without a `genecall_mode` key) and Parquet outputs."""
    fx = _stage_fixture(tmp_path)
    root, st = fx["root"], fx["stage"]

    p_omitted = genecall.plan(root, st, batch_cpu_min=240.0, tmp=tmp_path)
    p_default = genecall.plan(root, st, batch_cpu_min=240.0, tmp=tmp_path, genecall_mode="default")
    assert "genecall_mode" not in p_omitted
    assert "genecall_mode" not in p_default
    assert json.dumps(p_omitted, sort_keys=True) == json.dumps(p_default, sort_keys=True)

    out1, out2 = tmp_path / "out_omitted", tmp_path / "out_default"
    r1 = genecall.run(root, st, out=out1, tmp=tmp_path)
    r2 = genecall.run(root, st, out=out2, tmp=tmp_path, genecall_mode="default")
    assert r1 == r2
    assert _tree_hashes(out1) == _tree_hashes(out2)


def test_b67_gvmeta_recovers_code15_jumbophage_matching_oracle_one_proc_and_batched(tmp_path: Path):
    """On a >= 120 kb TAG-recoded (code 15) jumbophage staged as `kind='prokaryote'`,
    `default` fragments ORFs under code-11 single training, while `gv-meta` recovers the full-length
    code-15 ORFs matching `ViralGeneFinder(meta=True)` 1:1 in 1-process (`workers=1`, `workers=2`)
    and batched (`plan` -> `publish_plan` -> `run_batch` -> `merge`) runs."""
    fx = _stage_fixture(tmp_path)
    root, st = fx["root"], fx["stage"]

    # Independent oracle: direct pyrodigal_gv.ViralGeneFinder(meta=True)
    oracle_finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    oracle_genes = sorted(
        (g.begin - 1, g.end, g.strand, g.translation_table, g.translate(include_stop=False))
        for g in oracle_finder.find_genes(fx["jumbo_seq"].decode("ascii"))
    )
    oracle_code15 = [g for g in oracle_genes if g[3] == 15]
    assert len(oracle_code15) >= 80, f"expected >=80 code-15 ORFs from oracle, got {len(oracle_code15)}"

    # Default mode fragments the recoded ORFs under code-11 single training
    out_def = tmp_path / "out_def"
    res_def = genecall.run(root, st, out=out_def, tmp=tmp_path)
    assert "single" in res_def["modes"]
    fp_def = pq.read_table(
        out_def / f"l1/feature_protein/gene_call_version={genecall.GENE_CALL_VERSION}/stage={st}"
    ).to_pylist()
    jumbo_def = [r for r in fp_def if r["contig_uid"] == fx["cu_jumbo"]]
    assert all(r["genetic_code"] == 11 for r in jumbo_def)
    def_spans = {(r["start"], r["end"]) for r in jumbo_def}
    assert not any(sp in def_spans for sp in fx["spans"]), "code-11 single mode must fragment recoded TAG ORFs"

    # gv-meta 1-process (workers=1 and workers=2: exercises both worker threads and subprocesses)
    out_gv1 = tmp_path / "out_gv1"
    out_gv2 = tmp_path / "out_gv2"
    res_gv1 = genecall.run(root, st, out=out_gv1, tmp=tmp_path, workers=1, genecall_mode="gv-meta")
    res_gv2 = genecall.run(root, st, out=out_gv2, tmp=tmp_path, workers=2, genecall_mode="gv-meta")
    assert res_gv1["modes"] == ["gv-meta"]
    assert res_gv1 == res_gv2
    assert _tree_hashes(out_gv1) == _tree_hashes(out_gv2)

    # Compare gv-meta calls on the jumbophage contig 1:1 against the ViralGeneFinder(meta=True) oracle
    fp_gv = pq.read_table(
        out_gv1 / f"l1/feature_protein/gene_call_version={genecall.GENE_CALL_VERSION}/stage={st}"
    ).to_pylist()
    pseq_gv = {
        r["protein_uid"]: r["seq"]
        for r in pq.read_table(
            out_gv1 / f"l0seq/protein_seq/gene_call_version={genecall.GENE_CALL_VERSION}/stage={st}"
        ).to_pylist()
    }
    got_jumbo = sorted(
        (r["start"], r["end"], r["strand"], r["genetic_code"], pseq_gv[r["protein_uid"]])
        for r in fp_gv
        if r["contig_uid"] == fx["cu_jumbo"]
    )
    assert got_jumbo == oracle_genes

    # Batched: plan -> publish_plan -> run_batch -> merge must be byte-identical to 1-process gv-meta
    p_gv = genecall.plan(root, st, batch_cpu_min=240.0, batches=2, tmp=tmp_path, genecall_mode="gv-meta")
    assert p_gv["genecall_mode"] == "gv-meta"
    work = genecall.work_dir(tmp_path / "runs_tmp", root, st, p_gv)
    plan_path = genecall.write_plan(p_gv, work)
    genecall.publish_plan(root, st, plan_path, root, tmp=tmp_path)
    for j in range(len(p_gv["batches"])):
        rec = genecall.run_batch(root, st, plan_path, j, workers=2)
        assert rec["modes"] == ["gv-meta"]
        genecall.publish_batch(root, st, plan_path, rec, root)
    out_batched = tmp_path / "out_batched"
    res_batched = genecall.merge(root, st, plan_path, out=out_batched, tmp=tmp_path)
    assert res_batched == res_gv1
    assert _tree_hashes(out_batched) == _tree_hashes(out_gv1)


def test_b67_stage_and_cli_options(tmp_path: Path):
    """`stage.Build` / `stage.steps` and `cli.main` accept `--genecall-mode default|gv-meta`,
    omit the flag on `default`, and refuse invalid combinations."""
    prof = profile.Profile(
        max_threads_day=2, max_threads_night=2, mem_scope="4G", tmp=str(tmp_path / "ptmp"),
        durable=str(tmp_path / "dur"), bulk=str(tmp_path / "bulk"), nice_monolithic=0, nice_fanout=0,
        day_start_hour=8, day_end_hour=20, cpu_min_per_process=300, small_file_floor_bytes=1000,
    )
    b_def = bs.Build(root=tmp_path / "r", stage="S0", genomes=tmp_path / "g.json", objects=tmp_path / "o.json")
    gc_def = next(s for s in bs.steps(prof, b_def) if s.name == "genecall")
    assert "--genecall-mode" not in gc_def.argv

    b_gv = bs.Build(root=tmp_path / "r", stage="S0", genomes=tmp_path / "g.json", objects=tmp_path / "o.json",
                    genecall_mode="gv-meta")
    gc_gv = next(s for s in bs.steps(prof, b_gv) if s.name == "genecall")
    idx = gc_gv.argv.index("--genecall-mode")
    assert gc_gv.argv[idx + 1] == "gv-meta"

    with pytest.raises(ValueError, match="genecall_mode"):
        bs.steps(prof, bs.Build(root=tmp_path / "r", stage="S0", genomes=tmp_path / "g.json",
                                objects=tmp_path / "o.json", genecall_mode="bogus"))
    with pytest.raises(ValueError, match="stored"):
        bs.steps(prof, bs.Build(root=tmp_path / "r", stage="S0", genomes=tmp_path / "g.json",
                                objects=tmp_path / "o.json", depth="stored", genecall_mode="gv-meta"))

    # CLI `stage genecall --genecall-mode gv-meta` runs end-to-end
    fx = _stage_fixture(tmp_path / "cli_fx")
    out_cli = tmp_path / "cli_out"
    cli_tmp = tmp_path / "cli_tmp"
    cli_tmp.mkdir()
    rc = cli.main(["stage", "genecall", "--root", str(fx["root"]), "--stage", fx["stage"],
                   "--genecall-mode", "gv-meta", "--out", str(out_cli), "--tmp", str(cli_tmp)])
    assert rc == 0
    fp = pq.read_table(
        out_cli / f"l1/feature_protein/gene_call_version={genecall.GENE_CALL_VERSION}/stage={fx['stage']}"
    ).to_pylist()
    assert any(r["contig_uid"] == fx["cu_jumbo"] and r["genetic_code"] == 15 for r in fp)


def test_b67_stage_est_conf_before_and_rates(tmp_path: Path):
    """`stage_est.py` parses optional `GENECALL_MODE=gv-meta`, refuses unknown values, and prices
    `before` and `rates` with `gv-meta` rates while keeping default byte-identical."""
    est = _load_stage_est()
    fx = _stage_fixture(tmp_path / "est_fx")
    inp = fx["inp"]
    sums = {
        n: hashlib.sha256((inp / n).read_bytes()).hexdigest()
        for n in ("genomes.json", "objects.json", "stage.fna")
    }
    base_lines = [
        f"STAGE={fx['stage']}",
        "NOTE=planning/measurements/2026-11-20_stage_V1.md",
        f"INPUTS={inp}",
        *(f"SUM={h} {n}" for n, h in sums.items()),
        "DEPTH=called",
        "CRISPR=no",
    ]
    c_def = tmp_path / f"{fx['stage']}.conf"
    c_def.write_text("\n".join(base_lines) + "\n")
    parsed_def = est.read_conf(c_def, fx["stage"])
    assert parsed_def.genecall_mode == "default"

    c_gv_dir = tmp_path / "gv_conf"
    c_gv_dir.mkdir()
    c_gv = c_gv_dir / f"{fx['stage']}.conf"
    c_gv.write_text("\n".join([*base_lines, "GENECALL_MODE=gv-meta"]) + "\n")
    parsed_gv = est.read_conf(c_gv, fx["stage"])
    assert parsed_gv.genecall_mode == "gv-meta"

    for bad_line in ("GENECALL_MODE=bogus", "GENECALL_MODE=", "GENECALL_MODE=gv-meta\nGENECALL_MODE=gv-meta"):
        c_bad_dir = tmp_path / f"bad_{abs(hash(bad_line))}"
        c_bad_dir.mkdir()
        c_bad = c_bad_dir / f"{fx['stage']}.conf"
        c_bad.write_text("\n".join([*base_lines, bad_line]) + "\n")
        with pytest.raises(est.ConfError):
            est.read_conf(c_bad, fx["stage"])

    # `before` prices gv-meta higher than default on >=120 kb prokaryote and matches genecall.plan
    prof_path = tmp_path / "profile.yaml"
    prof_path.write_text(
        f"max_threads_day: 2\nmax_threads_night: 2\nmem_scope: 4G\ntmp: {tmp_path / 'pt'}\n"
        f"durable: {tmp_path / 'pd'}\nbulk: {tmp_path / 'pb'}\nnice_monolithic: 0\nnice_fanout: 0\n"
        "day_start_hour: 8\nday_end_hour: 20\ncpu_min_per_process: 300\nsmall_file_floor_bytes: 1000\n"
    )
    def _run_before(conf_path: Path, plan_out: Path) -> dict:
        p = subprocess.run(
            [sys.executable, str(BSP / "stage_est.py"), "before", "--conf", str(conf_path),
             "--stage", fx["stage"], "--root", str(fx["root"]), "--profile", str(prof_path),
             "--inputs", str(inp), "--tmp", str(tmp_path / f"bt_{plan_out.stem}"),
             "--out-plan", str(plan_out)],
            capture_output=True, text=True, check=True,
        )
        return json.loads(p.stdout)

    plan_def_p, plan_gv_p = tmp_path / "plan_def.json", tmp_path / "plan_gv.json"
    doc_def = _run_before(c_def, plan_def_p)
    doc_gv = _run_before(c_gv, plan_gv_p)
    assert "genecall_mode" not in doc_def
    assert "genecall_mode" not in json.loads(plan_def_p.read_text())
    assert doc_gv["genecall_mode"] == "gv-meta"
    assert json.loads(plan_gv_p.read_text())["genecall_mode"] == "gv-meta"

    want_gv = genecall.plan(fx["root"], fx["stage"], batch_cpu_min=240.0, tmp=tmp_path, genecall_mode="gv-meta")
    assert doc_gv["genecall"]["estimate_cpu_min"] == want_gv["estimate_cpu_min"]
    assert doc_gv["genecall"]["genes_expected"] == want_gv["genes_expected"]
    assert doc_gv["genecall"]["estimate_cpu_min"] > doc_def["genecall"]["estimate_cpu_min"]

    # `rates_doc` for a 1-process genecall run uses gv-meta's estimate when conf sets GENECALL_MODE=gv-meta
    runs_dir = fx["root"] / "runs" / f"stage={fx['stage']}"
    runs_dir.mkdir(parents=True, exist_ok=True)
    (runs_dir / "genecall.json").write_text(json.dumps({"cpu_min": 0.5, "exit_code": 0}))
    prof = profile.load(prof_path)
    rdoc_gv = est.rates_doc(parsed_gv, fx["root"], prof, inp / "objects.json", None, tmp_path / "rt_gv")
    gc_row = next(r for r in rdoc_gv["rows"] if r["step"] == "genecall")
    assert gc_row["estimate"] == round(want_gv["estimate_cpu_min"], 3)


def test_b67_stage_driver_passes_genecall_mode_flag(tmp_path: Path):
    """`stage_driver.sh` parses `GENECALL_MODE` from `conf.env`, appends `--genecall-mode gv-meta`
    only when set, and forwards it in `rehearse`."""
    from tests.r1.test_stage_driver import World

    w_def = World(tmp_path / "w_def", stage="BM1", depth="called", library="", crispr="no")
    for m in ("estimate", "rehearse", "before", "run"):
        assert w_def.drive(m).returncode == 0
    pod_def = [ln for ln in w_def.calls().splitlines() if ln.startswith("pod_run bm1-1 ")]
    assert len(pod_def) == 1
    assert "--genecall-mode" not in pod_def[0]

    w_gv = World(tmp_path / "w_gv", stage="BM2", depth="called", library="", crispr="no",
                 extra="GENECALL_MODE=gv-meta\n")
    for m in ("estimate", "rehearse", "before", "run"):
        assert w_gv.drive(m).returncode == 0
    pod_gv = [ln for ln in w_gv.calls().splitlines() if ln.startswith("pod_run bm2-1 ")]
    assert len(pod_gv) == 1
    assert "--genecall-mode gv-meta" in pod_gv[0]
    rh_conf = (w_gv.c / "rh" / "BM2R.conf").read_text()
    assert "GENECALL_MODE=gv-meta" in rh_conf
    pod_again = [ln for ln in w_gv.calls().splitlines() if "bm2-rh-again" in ln]
    assert any("--genecall-mode gv-meta" in ln for ln in pod_again)

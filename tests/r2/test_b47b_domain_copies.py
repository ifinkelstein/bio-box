"""Ring R2, card B-47b (DEC-35 (2)): a stage's Pfam-style scan copies the rows of proteins an earlier stage measured.

Fixture: the micro-S0 proteins (lambda, T4, M. genitalium, gene-called by the repo) split into four stages by their
index in protein_uid order (`MEMBERS`): B shares proteins with A and with C (some with both), D shares none. The HMM
library is the 15 lambda models of `test_micro_s0_domains` (gathering 25 bits). Every domain step runs as build-stage
runs a scan in shards, each part a real `run-step` (plan, shards, merge), so the sources have their run records.

- Either order: B scanned after A (copying from A) and B scanned alone give byte-identical domain files and ledger;
  so do A after B and A alone (MR-12, MR-08: accretion leaves the hits of existing targets as they were).
- The copied rows equal the `hmmsearch` binary's scan of those proteins (the oracle of `test_micro_s0_domains`).
- What is copied is decided from the records: every protein of B that A measured, with or without a hit, is left out
  of B's shards (sets and lengths computed here from the sequences, not from the plan), and the plan's counts are
  those sets'.
- Several sources: the first in the order given is the source, and B's files are the same in either order.
- A named stage sharing no protein, and a plan naming no other stage, give today's files and today's plan.
- A source scanned with another HMM file, another cutoff, another pyhmmer (its ledger rewritten), another Z (its plan
  rewritten), in one process (no plan) or not finished (no merge record) is passed over, and the plan says why.
- The merge refuses a source file whose sha256 is not the plan's.
- build-stage: a second stage built after a first whose scan has a plan copies from it through the plan step (given
  `--seq-stages`), and its files equal a direct full scan of its proteins.
"""

from __future__ import annotations

import dataclasses
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import fasta, ids, ledger, runstep
from biobox import io as bio
from biobox import stage as bs
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.registry.schemas import table

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402
import test_micro_s0_domains as m0  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("MR-12", "MR-08", "D-01")]

LIB, VER = "test", "1"
SCAN = f"domains-{LIB}-{VER}"
MEMBERS = {"A": (0, 3), "B": (0, 1, 2), "C": (0, 1, 4), "D": (5,)}  # index of a protein in uid order, mod 6
SPEC = table("protein_seq")


def _protein_rows(root: Path, stage: str) -> pa.Table:
    t = pa.concat_tables(pq.read_table(f, partitioning=None) for f in domains._protein_files(root, V, stage))
    return t.sort_by("protein_uid")


def _write_stage(root: Path, stage: str, rows: pa.Table) -> None:
    bio.write_table(rows, SPEC, root / f"l0seq/protein_seq/gene_call_version={V}/stage={stage}/protein_prefix=all"
                                       "/part-0.parquet")


def _run(prof, root: Path, stage: str, step: str, argv: list[str]) -> int:
    cmd = [sys.executable, "-m", "biobox.cli", *argv]
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0)


def _scan(w: dict, root: Path, stage: str, *, seq: tuple[str, ...] = (), hmm: Path | None = None,
          cutoff: str = "gathering", upto: str = "merge") -> None:
    """The stage's scan as build-stage runs it in shards, each part a run-step; `seq` is the plan's --seq-stages."""
    hmm = w["lib"] if hmm is None else hmm
    com = ["--root", str(root), "--stage", stage, "--library", LIB, "--library-version", VER]
    extra = ["--seq-stages", ",".join(seq)] if seq else []
    assert _run(w["prof"], root, stage, f"{SCAN}.plan", ["stage", "domains", "plan", *com, "--hmm", str(hmm),
                                                         "--cutoff", cutoff, "--shard-residues", str(w["size"]),
                                                         *extra]) == 0
    if upto == "plan":
        return
    for k in range(len(domains.read_plan(root, stage, LIB, VER)["shards"])):
        assert _run(w["prof"], root, stage, f"{SCAN}.shard-{k:05d}",
                    ["stage", "domains", "shard", *com, "--hmm", str(hmm), "--shard", str(k), "--cpus", "1",
                     "--tmp", str(w["scratch"])]) == 0
    if upto == "shards":
        return
    assert _run(w["prof"], root, stage, SCAN, ["stage", "domains", "merge", *com, "--tmp", str(w["scratch"])]) == 0


def _outputs(root: Path, stage: str) -> dict[str, bytes]:
    dirs = (f"l1/domains/library={LIB}-{VER}/stage={stage}", f"ledger/stage={stage}/step={SCAN}")
    return {p.relative_to(root).as_posix(): p.read_bytes() for d in dirs for p in sorted((root / d).rglob("*"))
            if p.is_file()}


def _copy(src: Path, dst: Path) -> Path:
    shutil.copytree(src, dst, symlinks=True)
    return dst


def _domain_rows(root: Path, stage: str) -> list[dict]:
    files = bio.part_files(root / f"l1/domains/library={LIB}-{VER}/stage={stage}")
    return [r for f in files for r in pq.read_table(f, partitioning=None).to_pylist()]


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    base = tmp_path_factory.mktemp("b47b")
    src = base / "src"
    m0._build_store(src)
    lib = base / "lib.hmm"
    m0._make_library(m0._lambda_proteins(src), lib)
    lib2 = base / "lib2.hmm"
    lib2.write_bytes(lib.read_bytes().split(b"//\n", 1)[1])  # every model but the first: another HMM file
    allp = _protein_rows(src, "S0")
    r0 = base / "r0"
    sets, lengths = {}, dict(zip(allp.column("protein_uid").to_pylist(), allp.column("length").to_pylist()))
    for s, keep in MEMBERS.items():
        mask = pa.array([i % 6 in keep for i in range(allp.num_rows)])
        rows = allp.filter(mask)
        _write_stage(r0, s, rows)
        sets[s] = set(rows.column("protein_uid").to_pylist())
    size = sum(lengths[u] for u in sets["B"]) // 2 + max(lengths.values())
    w = {"base": base, "src": src, "lib": lib, "lib2": lib2, "r0": r0, "sets": sets, "lengths": lengths,
         "size": size, "scratch": base / "scratch", "prof": c08.profile_at(base / "prof")}
    rA = _copy(r0, base / "rA")
    _scan(w, rA, "A")
    rAC = _copy(rA, base / "rAC")
    _scan(w, rAC, "C")
    ab = _copy(rA, base / "ab")
    _scan(w, ab, "B", seq=("B", "A"))
    ba = _copy(r0, base / "ba")
    _scan(w, ba, "B")
    _scan(w, ba, "A", seq=("A", "B"))
    acb = _copy(rAC, base / "acb")
    _scan(w, acb, "B", seq=("B", "A", "C"))
    cab = _copy(rAC, base / "cab")
    _scan(w, cab, "B", seq=("B", "C", "A"))
    rD = _copy(r0, base / "rD")
    _scan(w, rD, "D")
    _scan(w, rD, "B", seq=("B", "D"))
    w.update(rA=rA, rAC=rAC, ab=ab, ba=ba, acb=acb, cab=cab, rD=rD)
    return w


def _shard_proteins(root: Path, stage: str) -> int:
    recs = sorted((root / domains.scan_dir(stage, LIB, VER)).glob("shard-*.json"))
    return sum(json.loads(p.read_text())["proteins"] for p in recs)


def test_b47b_two_stages_built_in_either_order_give_the_full_scans_files(world):
    w = world
    assert _outputs(w["ab"], "B") and _outputs(w["ab"], "B") == _outputs(w["ba"], "B")  # B after A == B alone
    assert _outputs(w["ba"], "A") == _outputs(w["rA"], "A")  # A after B == A alone
    # the plan of B after A: what A measured is left out of B's shards, with or without a hit (sets from here)
    shared = w["sets"]["B"] & w["sets"]["A"]
    a_hit = {r["protein_uid"] for r in _domain_rows(w["rA"], "A")}
    assert shared & a_hit and shared - a_hit  # copied with rows, and copied as no rows
    plan = domains.read_plan(w["ab"], "B", LIB, VER)
    cp = plan["copy"]
    assert [s["stage"] for s in cp["sources"]] == ["A"] and cp["passed_over"] == []
    assert cp["proteins"] == len(shared) and cp["residues"] == sum(w["lengths"][u] for u in shared)
    assert cp["rows"] == sum(1 for r in _domain_rows(w["rA"], "A") if r["protein_uid"] in shared)
    assert cp["proteins_with_hit"] == len(shared & a_hit)
    scanned = w["sets"]["B"] - shared
    assert sum(s["proteins"] for s in plan["shards"]) == len(scanned) == _shard_proteins(w["ab"], "B")
    assert sum(s["residues"] for s in plan["shards"]) == sum(w["lengths"][u] for u in scanned)
    assert plan["proteins"] == len(w["sets"]["B"])  # the ledger's n_in stays the stage's
    # the plan names the files copied from, with the sha256 of the bytes on disk
    for f in cp["sources"][0]["files"]:
        assert f["sha256"] == ids.sha256_file(w["ab"] / f["path"]).hex() and "stage=A/" in f["path"]


@pytest.mark.skipif(shutil.which("hmmsearch") is None, reason="hmmsearch binary not on PATH")
def test_b47b_copied_rows_equal_the_hmmsearch_binary(world, tmp_path):
    w = world
    shared = w["sets"]["B"] & w["sets"]["A"]
    prots = domains.load_proteins(w["ab"], V, "B")
    faa = tmp_path / "p.faa"
    faa.write_text("".join(f">{u.hex()}\n{s.replace('*', '')}\n" for u, s in sorted(prots.items()) if u in shared))
    dom = tmp_path / "dom.tsv"
    subprocess.run(["hmmsearch", "--cut_ga", "-Z", "1e8", "--domZ", "1e8", "--cpu", "1", "--noali",
                    "--domtblout", str(dom), "-o", str(tmp_path / "out.txt"), str(w["lib"]), str(faa)], check=True)
    ref = {}
    for line in dom.read_text().splitlines():
        if not line.startswith("#"):
            f = line.split()
            ref[(f[0], f[4], int(f[19]) - 1, int(f[20]))] = (float(f[7]), float(f[13]), float(f[12]))
    ours = {(r["protein_uid"].hex(), r["model_acc"], r["env_start"], r["env_end"]):
            (r["seq_bitscore"], r["dom_bitscore"], r["dom_ievalue_zref"])
            for r in _domain_rows(w["ab"], "B") if r["protein_uid"] in shared}
    assert ours and set(ours) == set(ref)
    for k, (s, d, e) in ours.items():
        rs, rd, re_ = ref[k]
        assert s == pytest.approx(rs, abs=0.05) and d == pytest.approx(rd, abs=0.05)
        assert e == pytest.approx(re_, rel=0.02)


def test_b47b_the_first_named_source_is_taken_and_either_order_gives_the_same_files(world):
    w = world
    sb, sa, sc = w["sets"]["B"], w["sets"]["A"], w["sets"]["C"]
    assert sb & sa & sc and (sb & sc) - sa  # some proteins in both sources, some in C only
    for root, order in ((w["acb"], ("A", "C")), (w["cab"], ("C", "A"))):
        cp = domains.read_plan(root, "B", LIB, VER)["copy"]
        assert [s["stage"] for s in cp["sources"]] == list(order)
        first, second = (sa, sc) if order == ("A", "C") else (sc, sa)
        assert [s["proteins"] for s in cp["sources"]] == [len(sb & first), len(sb & second - first)]
        assert _outputs(root, "B") == _outputs(w["ba"], "B")


def test_b47b_a_stage_sharing_nothing_and_a_plan_naming_no_other_stage_give_todays_files(world, tmp_path):
    w = world
    assert not w["sets"]["B"] & w["sets"]["D"]
    cp = domains.read_plan(w["rD"], "B", LIB, VER)["copy"]
    assert [s["stage"] for s in cp["sources"]] == ["D"] and cp["proteins"] == 0 and cp["rows"] == 0
    assert _outputs(w["rD"], "B") == _outputs(w["ba"], "B")
    # naming only the stage itself is naming no other stage: the plan is today's, byte for byte
    alone, itself = _copy(w["r0"], tmp_path / "alone"), _copy(w["r0"], tmp_path / "itself")
    _scan(w, alone, "B", upto="plan")
    _scan(w, itself, "B", seq=("B",), upto="plan")
    rel = domains.scan_dir("B", LIB, VER) + "/plan.json"
    assert (itself / rel).read_bytes() == (alone / rel).read_bytes()
    assert "copy" not in json.loads((alone / rel).read_text())


def _rewrite_ledger_pyhmmer(root: Path) -> None:
    """A's ledger as if another pyhmmer had written it, its run record updated to the new bytes."""
    rec_path = root / f"runs/stage=A/{SCAN}.json"
    rec = json.loads(rec_path.read_text())
    rows = ledger.read(root, "A", SCAN)
    for r in rows:
        t = json.loads(r["tool_versions_json"])
        t["pyhmmer"] = "0.0.0"
        r["tool_versions_json"] = json.dumps(t, sort_keys=True)
    d = root / f"ledger/stage=A/step={SCAN}"
    shutil.rmtree(d)
    ledger.write(root, "A", SCAN, rows)
    for o in rec["outputs"]:
        if o["path"].startswith(f"ledger/stage=A/step={SCAN}/"):
            o["bytes"] = (root / o["path"]).stat().st_size
            o["sha256"] = ids.sha256_file(root / o["path"]).hex()
    rec_path.write_text(json.dumps(rec, sort_keys=True, indent=1) + "\n")


@pytest.mark.parametrize("kind", ["hmm", "cutoff", "pyhmmer", "z", "one-process", "incomplete"])
def test_b47b_a_source_scanned_otherwise_is_passed_over(world, tmp_path, kind):
    w = world
    root = tmp_path / "root"
    if kind in ("pyhmmer", "z"):
        _copy(w["rA"], root)
        if kind == "pyhmmer":
            _rewrite_ledger_pyhmmer(root)
        else:
            p = root / domains.scan_dir("A", LIB, VER) / "plan.json"
            doc = json.loads(p.read_text())
            doc["z"] = 1e9
            p.write_text(json.dumps(doc, sort_keys=True, indent=1) + "\n")
    else:
        _copy(w["r0"], root)
        if kind == "hmm":
            _scan(w, root, "A", hmm=w["lib2"])
        elif kind == "cutoff":
            _scan(w, root, "A", cutoff="none")
        elif kind == "incomplete":
            _scan(w, root, "A", upto="shards")
        else:  # in one process: no plan
            assert _run(w["prof"], root, "A", SCAN, ["stage", "domains", "--root", str(root), "--stage", "A",
                                                      "--library", LIB, "--library-version", VER, "--hmm",
                                                      str(w["lib"]), "--cutoff", "gathering"]) == 0
    _scan(w, root, "B", seq=("B", "A"), upto="plan")
    plan = domains.read_plan(root, "B", LIB, VER)
    cp = plan["copy"]
    assert cp["sources"] == [] and cp["proteins"] == 0
    assert [p["stage"] for p in cp["passed_over"]] == ["A"]
    why = cp["passed_over"][0]["why"]
    want = {"hmm": "hmm_sha256", "cutoff": "cutoff", "pyhmmer": "pyhmmer", "z": "z", "one-process": "plan",
            "incomplete": "run record"}[kind]
    assert want in why, why
    assert sum(s["proteins"] for s in plan["shards"]) == len(w["sets"]["B"])  # every protein is scanned


def test_b47b_the_merge_refuses_a_source_file_whose_sha256_is_not_the_plans(world, tmp_path):
    w = world
    root = _copy(w["rA"], tmp_path / "root")
    _scan(w, root, "B", seq=("B", "A"), upto="shards")
    f = bio.part_files(root / f"l1/domains/library={LIB}-{VER}/stage=A")[0]
    data = bytearray(f.read_bytes())
    data[len(data) // 2] ^= 0xFF  # the same size, other bytes
    f.unlink()
    f.write_bytes(bytes(data))
    with pytest.raises(ValueError, match="sha256"):
        domains.merge(root, "B", LIB, VER, tmp=w["scratch"])
    assert _outputs(root, "B") == {}


# --- build-stage: the plan step is given the root's other stages (card B-47b) ----------------------------------------


class Killed(Exception):
    pass


def _lambda_variant(work: Path) -> tuple[Path, Path]:
    """Genomes and objects files of a stage S2: lambda with three substitutions, a contig of its own whose proteins
    are mostly lambda's."""
    rec = fasta.read(m0.FIX / "NC_001416.1.fna.gz")[0]
    seq = bytearray(rec.seq.upper())
    for p in (1_000, 20_000, 40_000):
        seq[p] = ord("A") if seq[p] != ord("A") else ord("C")
    fa = work / "lambda_variant.fna.gz"
    fa.write_bytes(gzip.compress(b">lambda_variant\n" + bytes(seq) + b"\n", mtime=0))
    genomes, objects = work / "genomes-S2.json", work / "objects-S2.json"
    genomes.write_text(json.dumps([{"contig_uids": [ids.contig_uid(bytes(seq)).hex()], "genetic_code": 11,
                                    "kind": "virus", "circular": []}]))
    objects.write_text(json.dumps([{"path": str(fa), "source_name": "fixture", "source_release": "1",
                                    "license": "INSDC (public)", "redistributable": True}]))
    return genomes, objects


def test_b47b_build_stage_copies_from_an_earlier_stage_through_the_plan_step(world, tmp_path, capsys):
    w = c08.Work(tmp_path / "work")
    w.prepare()
    genomes2, objects2 = _lambda_variant(w.base)
    prof = dataclasses.replace(c08.profile_at(tmp_path / "prof"), cpu_min_per_process=1)  # S0's scan in shards
    root = tmp_path / "root"
    lib = (bs.Library(LIB, VER, world["lib"], "gathering"),)

    def runner(prof_, root_, stage, step, cmd, **kw):
        return runstep.run_step(prof_, root_, stage, step, cmd, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0, **kw)

    def kill(point: str) -> None:
        if point == f"after:{SCAN}":
            raise Killed(point)

    b0 = bs.Build(root=root, stage="S0", release="rB.1", stages=("S0",), genomes=w.genomes("S0"),
                  objects=w.objects("S0"), libraries=lib, threads=1)
    with pytest.raises(Killed):
        bs.build_stage(prof, b0, run=runner, _hook=kill)
    assert (root / domains.scan_dir("S0", LIB, VER) / "plan.json").exists()
    capsys.readouterr()
    b2 = dataclasses.replace(b0, stage="S2", release="rB.2", stages=("S0", "S2"), genomes=genomes2, objects=objects2)
    with pytest.raises(Killed):
        bs.build_stage(prof, b2, run=runner, _hook=kill)
    err = capsys.readouterr().err
    plan = domains.read_plan(root, "S2", LIB, VER)
    rec = json.loads((root / f"runs/stage=S2/{SCAN}.plan.json").read_text())
    assert "--seq-stages" in rec["argv"] and rec["argv"][rec["argv"].index("--seq-stages") + 1] == "S2,S0"
    s0 = set(_protein_rows(root, "S0").column("protein_uid").to_pylist())
    s2 = set(_protein_rows(root, "S2").column("protein_uid").to_pylist())
    assert s0 & s2 and s2 - s0
    assert [s["stage"] for s in plan["copy"]["sources"]] == ["S0"] and plan["copy"]["proteins"] == len(s0 & s2)
    assert "S0" in err and "copies" in err
    # the oracle: a full scan of S2's proteins in one process, in a root that holds nothing else
    fresh = tmp_path / "fresh"
    shutil.copytree(root / f"l0seq/protein_seq/gene_call_version={V}/stage=S2",
                    fresh / f"l0seq/protein_seq/gene_call_version={V}/stage=S2")
    domains.run(fresh, "S2", V, LIB, VER, world["lib"], "gathering")
    assert _outputs(root, "S2") == _outputs(fresh, "S2") != {}

"""Ring R2: a contig shared by genomes in two stages belongs to the stage that stored it first (card C-20, DEC-26).

The store is content-addressed: ingest stores a sequence once, in the first stage that sees it. When a later
stage's genome holds the same bytes, that stage reads the earlier stages' sequences (so single-mode gene calling
still trains on its whole genome), writes no feature, IGR or tape rows for the shared contig, and records the skip
in its result and ledger (reason `owned_by:<stage>`). `contig_dim` keeps the owning stage's genome for the contig
and lists the later row as skipped. A release that lists the later stage without the owning one is refused with a
message naming the owning stage.

The fixture (built from the micro-S0 genomes, written as new FASTA objects):

    S0: G0 = (M. genitalium M, plasmid P)        prokaryote, code 4, single mode
        GL = (phage lambda L)                    virus, gv-meta
    S1: G1 = (reverse complement of M, P)        prokaryote, code 4, single mode: P in the same mode and code
        G2 = (M rotated by 100 kb, L)            prokaryote, code 4, single mode: the S0 phage inside a prokaryote

P is the first 8,000 bases of phage T4. Both stages are built with `build-stage` (every step a real run-step
subprocess), each with a release: in order A, S0 then S1; in order B, S1 then S0. The result depends on that order
for the shared contigs only, and the tests say so: their calls, IGRs, tapes and `contig_dim` genome come from the
stage built first; every other row of a release over both stages is the same in both orders.
"""

from __future__ import annotations

import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, dim, fasta, ids, io, ledger, runstep
from biobox import stage as bs
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l0seq.store import revcomp
from biobox.l1 import genecall
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l3 import tapes
from biobox.registry.schemas import table

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

pytestmark = [
    pytest.mark.r2,
    pytest.mark.covers("D-01", "R-09"),
    pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)"),
]

FIX = c08.FIX
TAPE_V = tapes.tape_version([V])
# protein_seq is written per stage with the stage's distinct proteins, so a protein made by genomes of two stages
# has a row in each (stage.py: 10,634,028 rows for 10,620,684 proteins in rN.1, measured). That is the table's
# design, not a shared contig: the reverse-complemented chromosome of S1 makes every S0 protein again.
PER_STAGE_BY_DESIGN = {"protein_seq"}


def _seq(acc: str) -> bytes:
    return fasta.read(FIX / f"{acc}.fna.gz")[0].seq


M, L = _seq("NC_000908.2"), _seq("NC_001416.1")
P = _seq("NC_000866.4")[:8000]
M_RC = revcomp(M)
M_ROT = M[100_000:] + M[:100_000]
CU = {name: ids.contig_uid(s) for name, s in
      {"M": M, "L": L, "P": P, "M_RC": M_RC, "M_ROT": M_ROT}.items()}
SHARED = {CU["P"], CU["L"]}

# per stage: objects (file name, [(header, sequence)]) and genomes (contig names, code, kind, circular names)
STAGES = {
    "S0": {"objects": [("s0_host.fna", [("hostA chromosome", M), ("hostA plasmid", P)]),
                       ("s0_lambda.fna", [("phage lambda", L)])],
           "genomes": [(("M", "P"), 4, "prokaryote", ("M",)), (("L",), 11, "virus", ())]},
    "S1": {"objects": [("s1_hostB.fna", [("hostB chromosome", M_RC), ("hostB plasmid", P)]),
                       ("s1_hostC.fna", [("hostC chromosome", M_ROT), ("hostC prophage", L)])],
           "genomes": [(("M_RC", "P"), 4, "prokaryote", ("M_RC",)), (("M_ROT", "L"), 4, "prokaryote", ("M_ROT",))]},
}


def genome_inputs(stage: str) -> list[GenomeInput]:
    return [GenomeInput(tuple(CU[n] for n in names), code, kind=kind, circular=tuple(CU[n] for n in circ))
            for names, code, kind, circ in STAGES[stage]["genomes"]]


def genome_uid(names: tuple[str, ...]) -> bytes:
    return ids.genome_uid([CU[n] for n in names])


def write_inputs(base: Path) -> dict[str, tuple[Path, Path]]:
    """The FASTA objects and the CLI's genomes and objects JSON files of each stage."""
    base.mkdir(parents=True, exist_ok=True)
    out = {}
    for stage, spec in STAGES.items():
        objs = []
        for name, recs in spec["objects"]:
            lines = [ln for h, s in recs for ln in fasta.format_record(h, s, 80)]
            (base / name).write_bytes(b"\n".join(lines) + b"\n")
            objs.append({"path": str(base / name), "source_name": "c20_fixture", "source_release": "1",
                         "license": "INSDC (public)", "redistributable": True})
        g = base / f"genomes-{stage}.json"
        g.write_text(json.dumps([{"contig_uids": [u.hex() for u in x.contig_uids], "genetic_code": x.genetic_code,
                                  "kind": x.kind, "circular": [u.hex() for u in x.circular],
                                  "assembly_quality": x.assembly_quality} for x in genome_inputs(stage)]))
        o = base / f"objects-{stage}.json"
        o.write_text(json.dumps(objs))
        out[stage] = (g, o)
    return out


def source_objects(base: Path, stage: str) -> list[SourceObject]:
    return [SourceObject(base / name, "c20_fixture", "1", "INSDC (public)", True) for name, _ in STAGES[stage]["objects"]]


def runner(prof, root, stage, step, cmd, **kw) -> int:
    """run-step at a fixed hour and load, as `test_build_stage.runner`."""
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0, **kw)


def build(prof, root: Path, inputs, stage: str, release: str, stages: tuple[str, ...]) -> int:
    g, o = inputs[stage]
    b = bs.Build(root=root, stage=stage, release=release, stages=stages, genomes=g, objects=o, threads=1)
    return bs.build_stage(prof, b, run=runner)


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("c20")
    prof = c08.profile_at(d / "prof")
    inputs = write_inputs(d / "inputs")
    a, b = d / "a" / "root", d / "b" / "root"
    codes = {}
    codes["A1"] = build(prof, a, inputs, "S0", "rC20.a1", ("S0",))
    s0_only = d / "a-after-s0"
    if codes["A1"] == 0:
        shutil.copytree(a, s0_only, symlinks=True)
    codes["A2"] = build(prof, a, inputs, "S1", "rC20.a2", ("S0", "S1"))
    codes["B1"] = build(prof, b, inputs, "S1", "rC20.b1", ("S1",))
    codes["B2"] = build(prof, b, inputs, "S0", "rC20.b2", ("S1", "S0"))
    return {"d": d, "prof": prof, "inputs": inputs, "a": a, "b": b, "s0_only": s0_only, "codes": codes}


def _rows(root: Path, rel: str, columns=None) -> list[dict]:
    files = io.part_files(root / rel)
    if not files:
        return []
    return pa.concat_tables(pq.read_table(f, columns=columns, partitioning=None) for f in files).to_pylist()


def _release_rows(root: Path, release: str, name: str) -> list[dict]:
    """Every row of table `name` in the files of the release's manifest, with the file's stage."""
    man = catalog.read_manifest(catalog.manifest_path(root, release))
    out = []
    for f in man["files"]:
        if f["table"] == name:
            for r in pq.read_table(root / f["path"], partitioning=None).to_pylist():
                out.append({**r, "_stage": f["partitions"].get("stage")})
    return out


def test_both_orders_build_through_build_stage(world):
    assert world["codes"] == {"A1": 0, "A2": 0, "B1": 0, "B2": 0}


@pytest.mark.covers("L-01")
@pytest.mark.parametrize("order, later, owner", [("a", "S1", "S0"), ("b", "S0", "S1")])
def test_the_later_stage_writes_no_gene_igr_or_tape_row_for_a_shared_contig_and_ledgers_the_skip(
        world, order, later, owner):
    root = world[order]
    for rel, col in ((f"l1/feature/producer={V}/stage={later}", "contig_uid"),
                     (f"l1/feature_protein/gene_call_version={V}/stage={later}", "contig_uid"),
                     (f"l1/protein_occurrence/gene_call_version={V}/stage={later}", "contig_uid"),
                     (f"l0seq/igr/igr_version={TAPE_V}-igr1/stage={later}", "contig_uid"),
                     (f"l3/tapes/tape_version={TAPE_V}/stage={later}", "contig_uid")):
        got = {r[col] for r in _rows(root, rel, [col])}
        assert got, rel
        assert not got & SHARED, rel
    # the owning stage holds the shared contigs' genes, IGRs and tapes
    for rel in (f"l1/feature/producer={V}/stage={owner}", f"l0seq/igr/igr_version={TAPE_V}-igr1/stage={owner}",
                f"l3/tapes/tape_version={TAPE_V}/stage={owner}"):
        assert SHARED <= {r["contig_uid"] for r in _rows(root, rel, ["contig_uid"])}, rel
    rows = ledger.read(root, later, "genecall")
    skipped = {r["object_key"]: r for r in rows if r["reason_code"] != "ok"}
    assert set(skipped) == SHARED
    for r in skipped.values():
        assert r["reason_code"] == f"owned_by:{owner}"
        assert (r["object_kind"], r["n_in"], r["n_kept"], r["n_merged"], r["n_rejected"]) == ("contig", 0, 0, 0, 0)
    assert ledger.check_balance(root, later, ["genecall"]) == []
    # the stage's own contigs are all in its ledger, each once with reason ok
    own = {r["object_key"] for r in rows if r["reason_code"] == "ok"}
    assert own == {cu for g in genome_inputs(later) for cu in g.contig_uids} - SHARED


@pytest.mark.parametrize("order", ["a", "b"])
def test_every_contig_of_every_genome_has_its_tape_and_igrs_in_a_release_over_both_stages(world, order):
    root = world[order]
    rel = {"a": "rC20.a2", "b": "rC20.b2"}[order]
    tape = {r["contig_uid"] for r in _release_rows(root, rel, "contig_tapes")}
    igr = {r["contig_uid"] for r in _release_rows(root, rel, "igr")}
    want = {cu for s in STAGES for g in genome_inputs(s) for cu in g.contig_uids}
    assert tape == want and igr == want
    # one tape per contig, and each tape's slots are the features of its owning stage
    tapes_ = _release_rows(root, rel, "contig_tapes")
    assert len(tapes_) == len(want)
    feats = defaultdict(set)
    for r in _release_rows(root, rel, "feature"):
        feats[r["contig_uid"]].add(r["feature_uid"])
    for t in tapes_:
        assert set(t["feature_uids"]) == feats[t["contig_uid"]]


@pytest.mark.parametrize("order", ["a", "b"])
def test_no_key_is_duplicated_in_any_table_a_release_over_both_stages_reads(world, order):
    root = world[order]
    rel = {"a": "rC20.a2", "b": "rC20.b2"}[order]
    man = catalog.read_manifest(catalog.manifest_path(root, rel))
    by = defaultdict(list)
    for f in man["files"]:
        by[f["table"]].append(f)
    assert {"feature", "feature_protein", "protein_occurrence", "igr", "contig_tapes", "contig_dim", "contig_index",
            "seq_locator", "tape_tokens", "token_postings"} <= set(by)
    con = duckdb.connect()
    dup = {}
    for name, files in sorted(by.items()):
        if name in PER_STAGE_BY_DESIGN:
            continue
        spec = table(name)
        parts = {k for f in files for k in f["partitions"]}
        pk = [c for c in spec.primary_key if c in spec.schema.names or c in parts]

        def lit(v) -> str:
            return "NULL::VARCHAR" if v is None else "'" + str(v).replace("'", "''") + "'"

        # a key column is the file's column, or its partition value from the manifest (NULL where it has none)
        sel = " UNION ALL ".join(
            "SELECT " + ", ".join(f'"{c}"' if c in spec.schema.names else f'{lit(f["partitions"].get(c))} AS "{c}"'
                                  for c in pk)
            + f" FROM read_parquet({lit(root / f['path'])}, hive_partitioning = false)" for f in files)
        cols = ", ".join(f'"{c}"' for c in pk)
        (n,), = con.execute(f"SELECT count(*) FROM (SELECT {cols} FROM ({sel}) GROUP BY {cols} "
                            "HAVING count(*) > 1)").fetchall()
        if n:
            dup[name] = n
    assert dup == {}


def test_contig_dim_keeps_the_owning_stages_genome_and_the_order_decides_which_that_is(world):
    """DEC-26's accepted consequence, stated: the genome of a shared contig in `contig_dim` is the genome of the
    stage built first. Every other contig has the same row in both orders."""
    a = {r["contig_uid"]: r for r in _release_rows(world["a"], "rC20.a2", "contig_dim")}
    b = {r["contig_uid"]: r for r in _release_rows(world["b"], "rC20.b2", "contig_dim")}
    want = {cu for s in STAGES for g in genome_inputs(s) for cu in g.contig_uids}
    assert set(a) == set(b) == want
    for cu in want - SHARED:
        assert {k: v for k, v in a[cu].items() if k != "_stage"} == {k: v for k, v in b[cu].items() if k != "_stage"}
    assert a[CU["P"]]["genome_uid"] == genome_uid(("M", "P")) and b[CU["P"]]["genome_uid"] == genome_uid(("M_RC", "P"))
    assert a[CU["L"]]["genome_uid"] == genome_uid(("L",)) and b[CU["L"]]["genome_uid"] == genome_uid(("M_ROT", "L"))


def test_the_genes_of_every_unshared_contig_are_the_same_in_both_orders(world):
    """Gene calling trains on the whole genome, the shared contig included, whichever stage owns it, so the calls of
    the other contigs do not depend on the order. The shared contigs' calls come from the owning stage."""
    def by_contig(root, rel):
        out = defaultdict(set)
        for r in _release_rows(root, rel, "feature_protein"):
            out[r["contig_uid"]].add((r["feature_uid"], r["start"], r["end"], r["strand"], r["protein_uid"],
                                      r["genetic_code"], r["_stage"]))
        return out

    a, b = by_contig(world["a"], "rC20.a2"), by_contig(world["b"], "rC20.b2")
    for cu in {CU["M"], CU["M_RC"], CU["M_ROT"]}:
        assert a[cu] and a[cu] == b[cu]
    for cu in SHARED:
        assert {x[-1] for x in a[cu]} == {"S0"} and {x[-1] for x in b[cu]} == {"S1"}


def test_a_release_that_lists_the_later_stage_without_the_owning_stage_is_refused_naming_it(world, tmp_path, capfd):
    root = tmp_path / "root"
    shutil.copytree(world["a"], root, symlinks=True)
    assert build(world["prof"], root, world["inputs"], "S1", "rC20.late", ("S1",)) == runstep.EXIT_REFUSED
    err = capfd.readouterr().err
    assert "stage S0" in err and "DEC-26" in err
    assert not (root / "dim/contig_dim/release=rC20.late").exists()
    with pytest.raises(ValueError, match=r"stage S0 .*DEC-26"):
        dim.run(root, "rC20.late2", ["S1"])


@pytest.fixture(scope="module")
def direct_s1(world):
    """S1's inputs and ingest called directly on a copy of order A's root as it stood after S0."""
    if world["codes"]["A1"] != 0:
        pytest.fail("order A's first build failed")
    root = world["d"] / "direct-s1"
    shutil.copytree(world["s0_only"], root, symlinks=True)
    write_genome_input(root, "S1", genome_inputs("S1"))
    ingest(source_objects(world["d"] / "inputs", "S1"), root, "S1")
    return root


def test_gene_calling_names_the_owning_stage_when_it_is_not_read_and_reports_its_skips(direct_s1, tmp_path):
    with pytest.raises(KeyError, match=r"stage S0 .*DEC-26"):
        genecall.run(direct_s1, "S1", out=tmp_path / "o1", tmp=tmp_path)
    res = genecall.run(direct_s1, "S1", seq_stages=["S1", "S0"], out=tmp_path / "o2", tmp=tmp_path)
    assert res["skipped_contigs"] == [{"contig_uid": cu.hex(), "owner_stage": "S0"} for cu in sorted(SHARED)]
    assert res["modes"] == ["single"]


def test_the_ncrna_producer_skips_a_shared_contig_and_ledgers_it(direct_s1, tmp_path):
    if shutil.which("cmsearch") is None:
        pytest.skip("cmsearch (Infernal) is in the image")
    from biobox.l1 import ncrna

    lib = FIX / "rfam15.1_micro_s0.cm.gz"
    out = tmp_path / "o"
    res = ncrna.run(direct_s1, "S1", lib, 2, out=out, tmp=tmp_path, seq_stages=["S1", "S0"])
    assert res["skipped_contigs"] == [{"contig_uid": cu.hex(), "owner_stage": "S0"} for cu in sorted(SHARED)]
    prod = ncrna.producer(lib)
    feats = {r["contig_uid"] for r in _rows(out, f"l1/feature/producer={prod}/stage=S1", ["contig_uid"])}
    assert feats and not feats & SHARED
    rows = ledger.read(out, "S1", f"ncrna-{prod}")
    assert {r["object_key"] for r in rows if r["reason_code"] == "owned_by:S0"} == SHARED
    assert ledger.check_balance(out, "S1", [f"ncrna-{prod}"]) == []

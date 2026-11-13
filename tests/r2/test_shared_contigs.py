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
for the shared contigs only, and the tests say so: their calls, IGRs and tapes come from the stage built first;
every other row of a release over both stages is the same in both orders.

Card C-21a (DEC-27) adds a stage whose own genomes share a contig (section near the end): `genome_input` lists it in
each, the genome with the least `genome_uid` writes its calls, and the stage does not depend on the genomes' order.

Card C-21b (DEC-27, last section) reads a release's membership through the catalog view `genome_member` and writes
`genome_dim` and `contig_dim.n_genomes`; `contig_dim`'s genome of a shared contig is the least of its genomes, so
membership and the dimension tables do not depend on the build order or on the order `release_def` lists the stages.
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
from biobox.registry.schemas import conformance_problems, table

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
    assert duplicated_keys(root, rel) == {}


def duplicated_keys(root: Path, rel: str) -> dict[str, int]:
    """Per table of the release's manifest (but `PER_STAGE_BY_DESIGN`), the number of primary-key values held by more
    than one row; a table without any is not listed."""
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
    return dup


def test_contig_dim_keeps_the_least_genome_of_a_shared_contig_whatever_the_build_order(world):
    """Card C-21b (DEC-27) migrates C-20's `test_contig_dim_keeps_the_owning_stages_genome_and_the_order_decides_which_
    that_is`: DEC-26 gave a shared contig the genome of the stage built first, so its `contig_dim` row depended on the
    build order. Now its genome is the least `genome_uid` of all its genomes, a function of the release's set of
    stages, and every row is the same in both orders, the shared contigs' too."""
    a = {r["contig_uid"]: r for r in _release_rows(world["a"], "rC20.a2", "contig_dim")}
    b = {r["contig_uid"]: r for r in _release_rows(world["b"], "rC20.b2", "contig_dim")}
    want = {cu for s in STAGES for g in genome_inputs(s) for cu in g.contig_uids}
    assert set(a) == set(b) == want
    for cu in want:
        assert {k: v for k, v in a[cu].items() if k != "_stage"} == {k: v for k, v in b[cu].items() if k != "_stage"}
    assert a[CU["P"]]["genome_uid"] == min(genome_uid(("M", "P")), genome_uid(("M_RC", "P")))
    assert a[CU["L"]]["genome_uid"] == min(genome_uid(("L",)), genome_uid(("M_ROT", "L")))


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
    # card E-07: the fixture library's models on whole contigs (S1 is not gene-called here, so no intergenic models)
    model_set = FIX / "ncrna_set_whole.tsv"
    out = tmp_path / "o"
    res = ncrna.run(direct_s1, "S1", lib, 2, out=out, tmp=tmp_path, seq_stages=["S1", "S0"], model_set=model_set)
    assert res["skipped_contigs"] == [{"contig_uid": cu.hex(), "owner_stage": "S0"} for cu in sorted(SHARED)]
    prod = ncrna.producer(lib, model_set)
    feats = {r["contig_uid"] for r in _rows(out, f"l1/feature/producer={prod}/stage=S1", ["contig_uid"])}
    assert feats and not feats & SHARED
    rows = ledger.read(out, "S1", f"ncrna-{prod}")
    assert {r["object_key"] for r in rows if r["reason_code"] == "owned_by:S0"} == SHARED
    assert ledger.check_balance(out, "S1", [f"ncrna-{prod}"]) == []


# --- card C-21a (DEC-27): a contig in several genomes of one stage --------------------------------------------------
#
# Stage S2 holds three genomes that share the plasmid P:
#
#     GA = (M, P)        prokaryote, code 4, M circular: single mode, trained on M and P
#     GB = (M_ROT, P)    prokaryote, code 4, M_ROT circular: single mode, P in the same mode and code as GA
#     GC = (L, P)        prokaryote, code 4, 56.5 kb: under 100 kb, so meta mode: P in another mode
#
# `genome_input` lists P in all three. P's calls are written once, by the genome with the least genome_uid (DEC-27),
# whichever that is; the others' calls of P are dropped and the skip is recorded. A genome whose code differs from
# another's for P would be refused at `stage inputs` (the design note's rule, test_micro_s0_inputs.py). The stage is
# built twice with build-stage, in two fresh roots, with the genomes and the contigs of each listed in opposite orders.

W_STAGE, W_REL = "S2", "rC21.w"
W_OBJECTS = [("w_hostA.fna", [("hostA chromosome", M), ("hostA plasmid", P)]),
             ("w_hostB.fna", [("hostB chromosome", M_ROT), ("hostB plasmid", P)]),
             ("w_small.fna", [("small host phage", L), ("small host plasmid", P)])]
W_GENOMES = [(("M", "P"), 4, "prokaryote", ("M",)), (("M_ROT", "P"), 4, "prokaryote", ("M_ROT",)),
             (("L", "P"), 4, "prokaryote", ())]
W_GUID = {names: genome_uid(names) for names, *_ in W_GENOMES}
W_OWNER = min(W_GUID, key=W_GUID.get)  # the contig names of the genome that writes P's calls
SEQ = {"M": M, "L": L, "P": P, "M_RC": M_RC, "M_ROT": M_ROT}


def w_genomes(order: int = 1) -> list[GenomeInput]:
    """The stage's genomes; order -1 lists the genomes, and the contigs of each, the other way round."""
    return [GenomeInput(tuple(CU[n] for n in names[::order]), code, kind=kind, circular=tuple(CU[n] for n in circ))
            for names, code, kind, circ in W_GENOMES[::order]]


def _genomes_json(genomes: list[GenomeInput]) -> str:
    return json.dumps([{"contig_uids": [u.hex() for u in x.contig_uids], "genetic_code": x.genetic_code,
                        "kind": x.kind, "circular": [u.hex() for u in x.circular],
                        "assembly_quality": x.assembly_quality} for x in genomes])


def w_inputs(base: Path) -> Path:
    """The FASTA objects of stage S2 and its objects JSON file."""
    base.mkdir(parents=True, exist_ok=True)
    objs = []
    for name, recs in W_OBJECTS:
        lines = [ln for h, s in recs for ln in fasta.format_record(h, s, 80)]
        (base / name).write_bytes(b"\n".join(lines) + b"\n")
        objs.append({"path": str(base / name), "source_name": "c21a_fixture", "source_release": "1",
                     "license": "INSDC (public)", "redistributable": True})
    (base / "objects.json").write_text(json.dumps(objs))
    return base / "objects.json"


@pytest.fixture(scope="module")
def within(tmp_path_factory):
    d = tmp_path_factory.mktemp("c21a")
    prof = c08.profile_at(d / "prof")
    objects = w_inputs(d / "inputs")
    codes, roots = {}, {}
    for order in (1, -1):
        g = d / "inputs" / f"genomes{order}.json"
        g.write_text(_genomes_json(w_genomes(order)))
        roots[order] = d / f"order{order}" / "root"
        b = bs.Build(root=roots[order], stage=W_STAGE, release=W_REL, stages=(W_STAGE,), genomes=g, objects=objects,
                     threads=1)
        codes[order] = bs.build_stage(prof, b, run=runner)
    return {"d": d, "codes": codes, "roots": roots, "root": roots[1]}


def _oracle_calls(names: tuple[str, ...], code: int) -> set[tuple[int, int, int]]:
    """P's genes as pyrodigal itself gives them for the genome `names` (independent of `genecall`): meta mode below
    100 kb of genome, else single mode trained on the genome's contigs in contig_uid order; 0-based half-open."""
    import pyrodigal

    seqs = {CU[n]: SEQ[n].upper().decode() for n in names}
    if sum(len(s) for s in seqs.values()) < 100_000:
        finder = pyrodigal.GeneFinder(meta=True)
    else:
        finder = pyrodigal.GeneFinder(meta=False)
        finder.train(*(seqs[cu] for cu in sorted(seqs)), translation_table=code)
    return {(g.begin - 1, g.end, g.strand) for g in finder.find_genes(P.upper().decode())}


def test_c21a_a_stage_whose_genomes_share_a_contig_builds_through_build_stage(within):
    assert within["codes"] == {1: 0, -1: 0}


@pytest.mark.covers("C-01")
def test_c21a_genome_input_lists_the_shared_contig_in_every_genome(within):
    rows = _rows(within["root"], f"l0seq/genome_input/stage={W_STAGE}")
    got = [(r["contig_uid"], r["genome_uid"]) for r in rows]
    want = {(CU[n], W_GUID[names]) for names, *_ in W_GENOMES for n in names}
    assert sorted(got) == got and set(got) == want and len(got) == len(want)
    assert {gu for cu, gu in got if cu == CU["P"]} == set(W_GUID.values())


@pytest.mark.covers("D-01")
def test_c21a_the_shared_contig_has_one_call_set_the_least_genomes(within):
    """P's genes are those of the genome with the least genome_uid, by pyrodigal run directly as the oracle; the three
    genomes call P differently, so the rule decides something."""
    root = within["root"]
    by_genome = {names: _oracle_calls(names, code) for names, code, *_ in W_GENOMES}
    assert len({frozenset(v) for v in by_genome.values()}) > 1
    feats = [r for r in _rows(root, f"l1/feature/producer={V}/stage={W_STAGE}") if r["contig_uid"] == CU["P"]]
    assert {(r["start"], r["end"], r["strand"]) for r in feats} == by_genome[W_OWNER]
    assert len(feats) == len(by_genome[W_OWNER])
    for rel in (f"l1/feature_protein/gene_call_version={V}/stage={W_STAGE}",
                f"l1/protein_occurrence/gene_call_version={V}/stage={W_STAGE}"):
        rows = [r for r in _rows(root, rel) if r["contig_uid"] == CU["P"]]
        assert sorted(r["feature_uid"] for r in rows) == sorted(r["feature_uid"] for r in feats), rel


@pytest.mark.covers("L-01")
def test_c21a_the_ledger_records_the_skip_naming_the_owning_genome(within):
    root = within["root"]
    rows = [r for r in ledger.read(root, W_STAGE, "genecall") if r["object_key"] == CU["P"]]
    assert sorted(r["reason_code"] for r in rows) == ["ok", f"owned_by:{W_STAGE}/{W_GUID[W_OWNER].hex()}"]
    assert ledger.check_balance(root, W_STAGE, ["genecall"]) == []


def test_c21a_contig_dim_holds_one_row_per_contig_with_the_owning_genome(within):
    rows = _release_rows(within["root"], W_REL, "contig_dim")
    by = {r["contig_uid"]: r["genome_uid"] for r in rows}
    assert len(rows) == len(by) == 4
    assert by[CU["P"]] == W_GUID[W_OWNER]
    for names, *_ in W_GENOMES:
        assert by[CU[names[0]]] == W_GUID[names]


def test_c21a_no_key_is_duplicated_and_the_release_verifies(within):
    for root in within["roots"].values():
        assert duplicated_keys(root, W_REL) == {}
        man = catalog.read_manifest(catalog.manifest_path(root, W_REL))
        assert catalog.verify_manifest(man, root) == []


def test_c21a_every_tape_holds_each_feature_of_its_contig_once(within):
    """IGRs and tapes are per contig: one tape for P, whose slots are its features, each once and in order."""
    root = within["root"]
    feats = defaultdict(list)
    rows = _release_rows(root, W_REL, "feature")
    for r in sorted(rows, key=lambda r: (r["start"], r["end"], r["strand"], r["kind"])):
        feats[r["contig_uid"]].append(r["feature_uid"])
    tapes_ = _release_rows(root, W_REL, "contig_tapes")
    assert sorted(t["contig_uid"] for t in tapes_) == sorted(feats)
    for t in tapes_:
        assert t["feature_uids"] == feats[t["contig_uid"]], t["contig_uid"].hex()
        assert t["n_slots"] == len(feats[t["contig_uid"]])
    igr = defaultdict(int)
    for r in _release_rows(root, W_REL, "igr"):
        igr[r["contig_uid"]] += 1
    assert set(igr) == set(feats)


@pytest.mark.covers("D-01")
def test_c21a_the_stage_does_not_depend_on_the_order_of_its_genomes(within):
    """AGENTS §4.3: every file of the two roots (but run records, manifests and the catalog, which hold the root's
    path) is the same, byte for byte."""
    a, b = (c08.data(c08.tree(within["roots"][k])) for k in (1, -1))
    assert a and a == b


@pytest.fixture(scope="module")
def direct_w(within):
    """Stage S2's inputs and ingest called directly in a fresh root (no gene calls yet)."""
    root = within["d"] / "direct-w"
    ingest([SourceObject(within["d"] / "inputs" / name, "c21a_fixture", "1", "INSDC (public)", True)
            for name, _ in W_OBJECTS], root, W_STAGE)
    write_genome_input(root, W_STAGE, w_genomes())
    return root


def test_c21a_gene_calling_reports_the_skip_and_does_not_depend_on_genome_order(direct_w, tmp_path):
    want = [{"contig_uid": CU["P"].hex(), "owner_stage": W_STAGE, "owner_genome": W_GUID[W_OWNER].hex()}]
    res = genecall.run(direct_w, W_STAGE, out=tmp_path / "o1", tmp=tmp_path)
    assert res["skipped_contigs"] == want
    assert res["genes"] == sum(1 for _ in _rows(tmp_path / "o1", f"l1/feature/producer={V}/stage={W_STAGE}"))
    # the genomes given explicitly, in either order, and called two at a time: the same files. Each genome's contigs
    # are given in contig_uid order, as `genome_input` gives them: single mode trains on them in the order given, and
    # that order is not this card's question
    for k, (order, workers) in enumerate(((1, 2), (-1, 1), (-1, 2))):
        o = tmp_path / f"e{k}"
        gs = [GenomeInput(tuple(sorted(g.contig_uids)), g.genetic_code, kind=g.kind, circular=g.circular)
              for g in w_genomes(order)]
        r = genecall.run(direct_w, W_STAGE, genomes=gs, out=o, tmp=tmp_path, workers=workers)
        assert r == res
        assert c08.tree(o) == c08.tree(tmp_path / "o1"), (order, workers)


@pytest.mark.covers("L-01", "D-01")
def test_c21a_the_ncrna_producer_searches_a_shared_contig_once_and_ledgers_the_skip(direct_w, tmp_path):
    """Under card E-07's scan (a model set of whole-genome and intergenic models, on whole contigs and on the regions
    between the stage's gene calls), P is searched once and ledgered for the genome that owns it; the files are the
    same with the genomes given in the other order, and scanned in shards."""
    if shutil.which("cmsearch") is None or shutil.which("cmfetch") is None:
        pytest.skip("cmsearch and cmfetch (Infernal) are in the image")
    from biobox.l1 import ncrna

    import test_ncrna_regions as rg

    root, scratch = tmp_path / "root", tmp_path / "scratch"
    scratch.mkdir()
    shutil.copytree(direct_w, root, symlinks=True)
    genecall.run(root, W_STAGE, tmp=scratch)  # the intergenic models search between these genes
    ms = rg.write_set(tmp_path / "mixed.tsv", rg.MIXED)
    lib = rg.LIB
    res = ncrna.run(root, W_STAGE, lib, 2, model_set=ms, out=tmp_path / "o", tmp=scratch)
    owner = f"owned_by:{W_STAGE}/{W_GUID[W_OWNER].hex()}"
    assert res["skipped_contigs"] == [{"contig_uid": CU["P"].hex(), "owner_stage": W_STAGE,
                                       "owner_genome": W_GUID[W_OWNER].hex()}]
    assert res["contigs"] == 4 and res["bases"] == sum(len(SEQ[n]) for n in ("M", "M_ROT", "L", "P"))
    assert res["regions"] > 0
    prod = ncrna.producer(lib, ms)
    rows = ledger.read(tmp_path / "o", W_STAGE, f"ncrna-{prod}")
    assert sorted(r["reason_code"] for r in rows if r["object_key"] == CU["P"]) == ["ok", owner]
    assert len({r["object_key"] for r in rows}) == 4
    assert ledger.check_balance(tmp_path / "o", W_STAGE, [f"ncrna-{prod}"]) == []
    # the genomes given explicitly, the other way round: the same files
    res2 = ncrna.run(root, W_STAGE, lib, 2, genomes=w_genomes(-1), model_set=ms, out=tmp_path / "o2", tmp=scratch)
    assert res2 == res
    assert c08.tree(tmp_path / "o2") == c08.tree(tmp_path / "o")
    # in shards (each big contig a shard of its own): the files `run` writes
    o3 = tmp_path / "o3"
    p = ncrna.plan(root, W_STAGE, lib, len(M), model_set=ms, tmp=scratch, out=o3)
    assert p["shards"] >= 2 and p["contigs"] == 4
    for k in range(p["shards"]):
        ncrna.scan_shard(root, W_STAGE, lib, k, 1, model_set=ms, tmp=scratch, out=o3)
    ncrna.merge(root, W_STAGE, lib, model_set=ms, tmp=scratch, out=o3)
    produced = {k: v for k, v in c08.tree(o3).items() if not k.startswith(("l1/ncrna_scan/", "l1/domain_scan/"))}
    assert produced == c08.tree(tmp_path / "o")


# --- card C-21b (DEC-27): genome membership, `genome_dim` and `contig_dim.n_genomes` ------------------------------
#
# A release's membership is the union of its stages' `genome_input` rows, read through the catalog view `genome_member`
# (contig, genome, stage) and the macro `genomes_of`. `contig_dim` keeps one row per contig, its genome the least of
# the contig's genomes and `n_genomes` their number; `genome_dim` holds every genome whole. The oracle is the fixture's
# own genome lists (`W_GENOMES`, `STAGES`) and sequence lengths, not the code under test.

LEN = {CU[n]: len(s) for n, s in SEQ.items()}


def _con(root: Path, rel: str):
    return catalog.open_catalog(catalog.catalog_path(root, rel))


def _member(root: Path, rel: str) -> list[tuple]:
    con = _con(root, rel)
    try:
        return con.execute("SELECT contig_uid, genome_uid, stage FROM genome_member ORDER BY ALL").fetchall()
    finally:
        con.close()


def _strip(rows: list[dict]) -> list[dict]:
    return [{k: v for k, v in r.items() if k != "_stage"} for r in rows]


@pytest.mark.covers("C-01")
def test_c21b_genome_member_lists_every_genome_of_a_contig_shared_within_a_stage(within):
    want = sorted((CU[n], W_GUID[names], W_STAGE) for names, *_ in W_GENOMES for n in names)
    assert _member(within["root"], W_REL) == want
    con = _con(within["root"], W_REL)
    try:
        least_first = [(g,) for g in sorted(W_GUID.values())]
        assert con.execute("SELECT genome_uid FROM genomes_of(?)", [[CU["P"]]]).fetchall() == least_first
        lit = f"[from_hex('{CU['P'].hex()}')]"
        assert con.execute(f"SELECT genome_uid FROM genomes_of({lit})").fetchall() == least_first
        two = con.execute("SELECT contig_uid, genome_uid FROM genomes_of(?)", [[CU["M"], CU["L"]]]).fetchall()
        assert two == sorted([(CU["M"], W_GUID[("M", "P")]), (CU["L"], W_GUID[("L", "P")])])
        assert con.execute("SELECT count(*) FROM genomes_of(?)", [[CU["M_RC"]]]).fetchone()[0] == 0  # not in it
    finally:
        con.close()


@pytest.mark.covers("C-01")
def test_c21b_contig_dim_counts_a_shared_contigs_genomes_and_keeps_the_least(within):
    rows = {r["contig_uid"]: r for r in _release_rows(within["root"], W_REL, "contig_dim")}
    assert len(rows) == 4
    assert (rows[CU["P"]]["genome_uid"], rows[CU["P"]]["n_genomes"]) == (min(W_GUID.values()), 3)
    for names, *_ in W_GENOMES:
        assert (rows[CU[names[0]]]["genome_uid"], rows[CU[names[0]]]["n_genomes"]) == (W_GUID[names], 1)
    assert {r["assembly_quality"] for r in rows.values()} == {"complete"}  # the genomes agree


@pytest.mark.covers("C-01")
def test_c21b_genome_dim_holds_every_genome_whole(within):
    root = within["root"]
    rows = _strip(_release_rows(root, W_REL, "genome_dim"))
    want = [{"genome_uid": W_GUID[names], "n_contigs": len(names), "genome_len": sum(LEN[CU[n]] for n in names),
             "kind": kind, "genetic_code": code, "assembly_quality": "complete", "tier": None, "tax_order": None,
             "is_species_rep": False} for names, code, kind, _ in W_GENOMES]
    assert rows == sorted(want, key=lambda r: r["genome_uid"])
    man = catalog.read_manifest(catalog.manifest_path(root, W_REL))
    files = [root / f["path"] for f in man["files"] if f["table"] == "genome_dim"]
    assert files and all(f"release={W_REL}/" in p.as_posix() for p in files)
    for p in files:
        assert conformance_problems(pq.read_table(p, partitioning=None), table("genome_dim")) == []
    con = _con(root, W_REL)
    try:
        cols = ", ".join(table("genome_dim").schema.names)
        got = con.execute(f"SELECT {cols} FROM genome_dim_now() ORDER BY genome_uid").fetchall()
        assert got == [tuple(r.values()) for r in rows]
    finally:
        con.close()


@pytest.mark.covers("D-01")
def test_c21b_membership_and_dims_across_stages_do_not_depend_on_the_build_order(world):
    """Orders A and B build S0 and S1 the other way round, so each owns the shared contigs' calls in one of them
    (DEC-26). Membership, `contig_dim` and `genome_dim` are the same in both, byte for byte, and S1's genome that
    shares P with S0 is whole in `genome_dim` in order A, where S0 owns P."""
    a_rel, b_rel = "rC20.a2", "rC20.b2"
    want = sorted((cu, ids.genome_uid(list(g.contig_uids)), s) for s in STAGES for g in genome_inputs(s)
                  for cu in g.contig_uids)
    assert _member(world["a"], a_rel) == _member(world["b"], b_rel) == want
    for name in ("contig_dim", "genome_dim"):
        got = []
        for root, rel in ((world["a"], a_rel), (world["b"], b_rel)):
            man = catalog.read_manifest(catalog.manifest_path(root, rel))
            got.append([(root / f["path"]).read_bytes() for f in man["files"] if f["table"] == name])
        assert got[0] and got[0] == got[1], name
    dims = {r["contig_uid"]: r for r in _release_rows(world["a"], a_rel, "contig_dim")}
    assert {cu for cu, r in dims.items() if r["n_genomes"] == 2} == SHARED
    assert {r["n_genomes"] for cu, r in dims.items() if cu not in SHARED} == {1}
    genomes = {r["genome_uid"]: r for r in _release_rows(world["a"], a_rel, "genome_dim")}
    assert len(genomes) == sum(len(STAGES[s]["genomes"]) for s in STAGES)
    for s in STAGES:
        for g in genome_inputs(s):
            r = genomes[ids.genome_uid(list(g.contig_uids))]
            assert (r["n_contigs"], r["genome_len"], r["kind"], r["genetic_code"]) == (
                len(g.contig_uids), sum(LEN[cu] for cu in g.contig_uids), g.kind, g.genetic_code)


@pytest.mark.covers("D-01")
def test_c21b_dims_and_membership_do_not_depend_on_the_order_release_def_lists_the_stages(world, tmp_path):
    got = []
    for k, stages in enumerate((["S0", "S1"], ["S1", "S0"])):
        root, rel = tmp_path / f"r{k}", f"rC21b.o{k}"
        shutil.copytree(world["a"], root, symlinks=True)
        dim.run(root, rel, stages, tmp=tmp_path)
        man = catalog.build_manifest(root, rel, stages)
        con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / f"c{k}.duckdb"))
        try:
            mem = con.execute("SELECT contig_uid, genome_uid, stage FROM genome_member ORDER BY ALL").fetchall()
        finally:
            con.close()
        files = {name: [(root / f["path"]).read_bytes() for f in man["files"] if f["table"] == name]
                 for name in ("contig_dim", "genome_dim")}
        assert mem and files["contig_dim"] and files["genome_dim"]
        got.append((mem, files))
    assert got[0] == got[1]

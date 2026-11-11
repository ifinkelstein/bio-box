"""Ring R2: clustering in parts (card B-41, design `planning/design/B-37_split_clustering.md`).

The oracle is `cluster.run`, the builder's two `cluster` workflow calls in one process, kept as it is for this.
The split path replays the module calls those workflow calls make, one process each, with every per-query
module (`prefilter`, `rescorediagonal`, `align`, and the builder's member `align`) run as N processes over
consecutive key ranges of its input and joined in key order; prefilter parts ask for the self hit
(`--add-self-matches 1`). On two inputs, the micro-S0 fixture corpus (plus a 2% mutated lambda, as
`test_micro_s0_cluster`) and a constructed set with dense families and many length ties (the B-37 design's
`synth`, smaller), the split build at 1, 2 and 4 parts must give `cluster.run`'s files byte for byte, and every
intermediate database must be equal entry for entry at any part count.

Also: the hard-coded module list equals what the image's MMseqs2 logs for the two workflow calls at `-v 3`
(an upgrade that changes the workflow fails here); without `--add-self-matches 1` the prefilter parts give other
intermediate databases (the reason the flag is there); a deleted scratch database is found, rebuilt to its
recorded digest, and the build ends byte-identical; a step that cannot reproduce its record is refused, naming it.

Needs `mmseqs` on PATH (the bio-box image); skipped otherwise.
"""

from __future__ import annotations

import hashlib
import json
import random
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pytest

from biobox import fasta, ids, io
from biobox.l0seq import SourceObject, ingest
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import GenomeInput
from biobox.l1.genecall import run as call_genes
from biobox.l2 import cluster
from biobox.registry.schemas import table

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_cluster as c24  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("MR-16", "C-11"),
              pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")]

REL = "rB41.1"
THREADS = 2
PLAN_RECORDS = "l2/cluster_plan/"  # the split build's own records: never in a manifest


# --- inputs ---------------------------------------------------------------------------------------------


def fixture_root(root: Path) -> Path:
    """The micro-S0 corpus and a 2% mutated lambda, ingested and gene-called as `test_micro_s0_cluster` does."""
    lam = fasta.read(c24.FIX / "NC_001416.1.fna.gz")[0].seq
    src = root.parent / "src"
    src.mkdir(parents=True, exist_ok=True)
    mut = src / "lambda_mut.fna"
    mut.write_bytes(b"\n".join(fasta.format_record("lambda_mut", c24._mutate(lam, 0.02, 11), 70)) + b"\n")
    objs = [SourceObject(c24.FIX / o["file"], "fixture", "1") for o in c24.CONTRACT["objects"]] + [
        SourceObject(mut, "mut", "1")]
    ingest(objs, root, "S0")
    genomes = []
    for o in c24.CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(c24.FIX / o["file"])[0].seq)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                                   circular=(cu,) if o["topology"] == "circular" else ()))
    genomes.append(GenomeInput((ids.contig_uid(fasta.read(mut)[0].seq),), 11, kind="virus"))
    call_genes(root, "S0", genomes)
    return root


def synth_root(root: Path, bases: int = 300, seed: int = 37) -> Path:
    """The B-37 design's constructed set at a fifth of its size: `bases` random proteins of 80 to 299 residues,
    each with 5 variants at 25% substitution (about 75% identity: one NR50 cluster, several NR90 ones) and each
    variant with 4 copies at 4% (one NR90 cluster). Substitutions only, so lengths tie within a family, where a
    greedy clustering by length could depend on order. Written as stage S0's `protein_seq`."""
    rng = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def mutate(s: str, rate: float) -> str:
        return "".join(rng.choice(aa.replace(c, "")) if rng.random() < rate else c for c in s)

    seqs = set()
    for _ in range(bases):
        base = "M" + "".join(rng.choice(aa) for _ in range(rng.randint(79, 298)))
        for _ in range(5):
            var = mutate(base, 0.25)
            for _ in range(4):
                seqs.add(mutate(var, 0.04))
    rows = sorted((ids.protein_uid(s), s) for s in seqs)
    t = pa.table({"protein_uid": [u for u, _ in rows],
                  "protein_digest": [hashlib.sha256(s.encode()).digest() for _, s in rows],
                  "length": [len(s) for _, s in rows], "seq": [s for _, s in rows]})
    io.write_table(t, table("protein_seq"), root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0/"
                                                  "protein_prefix=all/part-0.parquet")
    return root


def tree(root: Path, skip: str | None = PLAN_RECORDS) -> dict[str, str]:
    root = Path(root)
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file() and not (skip and p.relative_to(root).as_posix().startswith(skip))}


def records(out: Path) -> dict[str, dict]:
    base = Path(out) / PLAN_RECORDS / f"release={REL}"
    return {p.name.removesuffix(".json"): json.loads(p.read_bytes()) for p in sorted(base.glob("*.json"))}


def digests(out: Path) -> dict[str, dict]:
    """Every database each step wrote, by step: its outputs' digests. Part steps (whose outputs depend on the cut)
    are left out; a joined call's record holds the joined database."""
    return {n: r["outputs"] for n, r in records(out).items() if ".part-" not in n and not n.endswith(".cut")
            and n != cluster.PLAN_STEP}


INPUTS = {"fixture": fixture_root, "synth": synth_root}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """For each input: its root, `cluster.run`'s output root, and the split builds' output roots at 1, 2 and 4
    parts (and at 4 without the self hit on the fixture)."""
    base = tmp_path_factory.mktemp("b41")
    out = {}
    for name, make in INPUTS.items():
        root = make(base / name / "root")
        oracle = base / name / "oracle"
        (base / name / "tmp-oracle").mkdir()
        cluster.run(root, REL, V, ["S0"], threads=THREADS, tmp=base / name / "tmp-oracle", out=oracle)
        split = {}
        for n in (1, 2, 4):
            split[n] = base / name / f"split-{n}"
            cluster.run_split(root, REL, V, ["S0"], parts=n, threads=THREADS, tmp=base / name / f"tmp-{n}",
                              out=split[n])
        out[name] = {"root": root, "oracle": oracle, "split": split, "base": base / name}
    fx = out["fixture"]
    fx["naive"] = base / "fixture" / "naive-4"
    cluster.run_split(fx["root"], REL, V, ["S0"], parts=4, threads=THREADS, tmp=base / "fixture" / "tmp-naive",
                      out=fx["naive"], _self_matches=False)
    return out


# --- the module list ----------------------------------------------------------------------------------------

_LINCLUST = re.compile(r"\btmp(90|50)/\d+/linclust/\d+")
_CASCADE = re.compile(r"\btmp(90|50)/\d+")


def logged_calls(log: str) -> dict[str, list[list[str]]]:
    """The module calls each `cluster` workflow call made, from its `-v 3` log: every line that starts with a
    module's name between that call's header and the next one (`kmermatcher` prints its own command line twice:
    a line equal to the one before is one call), without `--threads` and `-v`, with the workflow's temporary
    directories renamed as the split path names them."""
    modules = {"kmermatcher", "rescorediagonal", "clust", "createsubdb", "filterdb", "prefilter", "align",
               "mergeclusters"}
    out: dict[str, list[list[str]]] = {}
    cur, last = None, None
    for line in log.splitlines():
        m = re.match(r"^### (nr90|nr50)$", line)
        if m:
            cur, last = m.group(1), None
            out[cur] = []
            continue
        word = line.split(" ", 1)[0]
        if cur is None or word not in modules:
            continue
        if line != last:
            line_ = _LINCLUST.sub(lambda x: f"tmp{x.group(1)}/w/linclust/l", line)
            line_ = _CASCADE.sub(lambda x: f"tmp{x.group(1)}/w", line_)
            argv, keep, i = shlex.split(line_), [], 0
            while i < len(argv):
                if argv[i] in ("--threads", "-v"):
                    i += 2
                    continue
                keep.append(argv[i])
                i += 1
            out[cur].append(keep)
        last = line
    return out


def test_b41_the_hard_coded_module_list_is_the_one_the_images_mmseqs2_logs(built, tmp_path):
    """The two `cluster` calls of `cluster.run`, with its arguments, at `-v 3` on the fixture: the module calls they
    log are `cluster.workflow(level)`, call for call and argument for argument (15 at NR90, 23 at NR50)."""
    w = tmp_path
    prot = built["fixture"]["root"]
    files = cluster._protein_files(prot, V, ["S0"])
    assert cluster._write_fasta(files, w / "in.faa", w / "prot.parquet") > 800
    log = w / "v3.log"

    def mm(*args: str, head: str | None = None) -> None:
        with open(log, "a") as fh:
            if head:
                fh.write(f"### {head}\n")
                fh.flush()
            subprocess.run(["mmseqs", *args, "-v", "3", "--threads", "1"], check=True, cwd=w, stdout=fh,
                           stderr=subprocess.STDOUT)

    subprocess.run(["mmseqs", "createdb", "in.faa", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "1",
                    "-v", "1"], check=True, cwd=w)
    mm("cluster", "db", "clu90", "tmp90", "--min-seq-id", "0.9", "-c", "0.8", *cluster.COMMON, head="nr90")
    subprocess.run(["mmseqs", "createsubdb", "clu90", "db", "rep90", "-v", "1"], check=True, cwd=w)
    mm("cluster", "rep90", "clu50r", "tmp50", "--min-seq-id", "0.5", "-c", "0.8", *cluster.COMMON, head="nr50")
    got = logged_calls(log.read_text())
    for level, n in (("nr90", 15), ("nr50", 23)):
        want = [list(c) for c in cluster.workflow(level)]
        assert len(got[level]) == n == len(want), level
        for i, (g, x) in enumerate(zip(got[level], want)):
            assert g == x, f"{level} call {i}: the image logs {shlex.join(g)}\nthe split path runs {shlex.join(x)}"


# --- equality with cluster.run --------------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(INPUTS))
@pytest.mark.parametrize("n", [1, 2, 4])
def test_b41_the_split_build_equals_cluster_run_file_for_file(built, name, n):
    b = built[name]
    want = tree(b["oracle"])
    assert {k.split("/")[1] for k in want if k.startswith("l2/")} == {"cluster_member", "cluster_member_aln",
                                                                     "member_cluster"}
    assert tree(b["split"][n]) == want
    recs = records(b["split"][n])
    per_query = [c for c in cluster.calls() if c.split]
    assert len(per_query) == 16  # 4 prefilter, 6 align and 4 rescorediagonal calls, and the two member aligns
    for c in per_query:
        parts = recs[f"{c.name}.cut"]["parts"]
        if n == 1 or name == "synth":  # every per-query input of the constructed set has 4 keys or more
            assert len(parts) == n, c.name
        assert [f"{c.name}.part-{k:05d}" in recs for k in range(len(parts))] == [len(parts) > 1] * len(parts)
    # nothing of MMseqs2's scratch is in the output root
    assert not [k for k in tree(b["split"][n], skip=None) if not k.startswith(("l2/", "ledger/", "release_def/"))]


@pytest.mark.parametrize("name", sorted(INPUTS))
def test_b41_every_intermediate_database_is_the_same_at_any_part_count(built, name):
    b = built[name]
    one = digests(b["split"][1])
    assert len(one) == len(cluster.calls()) and all(one.values())
    for n in (2, 4):
        assert digests(b["split"][n]) == one, n


def test_b41_without_the_self_hit_the_prefilter_parts_give_other_intermediates(built):
    """The B-37 design measured 7 databases that differ on the fixture at 4 parts (every prefilter step and 3
    align steps) when the parts do not ask for the self hit; the final clusters happened to be equal there."""
    fx = built["fixture"]
    one, naive = digests(fx["split"][1]), digests(fx["naive"])
    differ = sorted(n for n in one if one[n] != naive[n])
    prefilters = [c.name for c in cluster.calls() if c.module == "prefilter"]
    assert set(prefilters) <= set(differ), differ
    assert all(digests(fx["split"][4])[n] == one[n] for n in differ)  # with the flag, the same parts are equal


# --- lost scratch -----------------------------------------------------------------------------------------------


def test_b41_a_deleted_scratch_database_is_rebuilt_to_its_digest_and_the_build_ends_byte_identical(built):
    fx = built["fixture"]
    base = fx["base"] / "lost"
    tmp, out = base / "tmp", base / "out"
    cluster.run_split(fx["root"], REL, V, ["S0"], parts=2, threads=THREADS, tmp=tmp, out=out,
                      until=cluster.FINAL_STEP)
    w = cluster.work_dir(tmp, fx["root"], REL, out=out)
    assert cluster.lost(fx["root"], REL, tmp, out=out) == []
    # the joined NR90 cascade prefilter (a per-query step in parts) and NR50's linclust redundancy (one process)
    pref = next(c for c in cluster.calls() if c.name.startswith("cluster.nr90.") and c.module == "prefilter")
    clust = next(c for c in cluster.calls() if c.name == "cluster.nr50.02.clust")
    for c in (pref, clust):
        for p in (w / c.output).parent.glob(Path(c.output).name + "*"):
            p.unlink()
    lost = cluster.lost(fx["root"], REL, tmp, out=out)
    parts = [k for k in records(out) if k.startswith(f"{pref.name}.part-")]
    assert len(parts) == 2 and lost == [*parts, pref.name, clust.name]
    for name in lost:
        cluster.redo(fx["root"], REL, name, tmp, threads=1, out=out)
    assert cluster.lost(fx["root"], REL, tmp, out=out) == []
    cluster.run_split(fx["root"], REL, V, ["S0"], parts=2, threads=THREADS, tmp=tmp, out=out)
    assert tree(out) == tree(fx["oracle"])


def test_b41_a_step_that_cannot_reproduce_its_record_is_refused_naming_it(built, tmp_path):
    fx = built["fixture"]
    tmp, out = tmp_path / "tmp", tmp_path / "out"
    cluster.run_split(fx["root"], REL, V, ["S0"], parts=2, threads=THREADS, tmp=tmp, out=out,
                      until="cluster.nr90.01.rescorediagonal.cut")
    name = "cluster.nr90.00.kmermatcher"
    rec = out / PLAN_RECORDS / f"release={REL}" / f"{name}.json"
    doc = json.loads(rec.read_bytes())
    path = next(iter(doc["outputs"]))
    rec.write_text(json.dumps(doc | {"outputs": {path: "0" * 64}}, sort_keys=True, indent=1) + "\n")
    assert cluster.lost(fx["root"], REL, tmp, out=out) == [name]
    with pytest.raises(ValueError, match=re.escape(name)):
        cluster.redo(fx["root"], REL, name, tmp, threads=1, out=out)

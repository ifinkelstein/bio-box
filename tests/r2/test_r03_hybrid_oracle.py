"""Ring R2: an independent oracle for the hybrid clustering modes (card R-03; AGENTS §5.2).

For `hybrid-a`, `hybrid-b1` and `linclust`, `cluster.run` replays `cluster.workflow(level, mode)`, the call list
`run_split` replays too, so the byte identity of tests/r1/test_cluster_hybrid.py and tests/r1/test_cluster_pruning.py
compares two drivers of one decomposition. Here each level is compared with MMseqs2 run as a user runs it, never
through `workflow()`, as sets of member sets (representatives and row order are not compared). There are three
inputs: the micro-S0 fixture's proteins (as tests/r2/test_cluster_split.py builds them, 924 proteins), a planted set
(`planted_root`, seed 20261125) and a set near the NR50 threshold where the cascade's second step changes the
clusters (`remote_root`, seed 2). Every MMseqs2 call uses 2 threads; the database is made as the builder makes its
own (FASTA sorted by uid, `createdb --shuffle 0`), from the proteins read here.

Open question (for the PI, card R-03): plain `mmseqs linclust` with the level's parameters differs from the build in
exactly the three flags `DECOMPOSITION_FLAGS["nr90"]`: the build has 857 NR90 clusters against plain linclust's 859
on the fixture and 240 against 239 on the planted set (measured on wilkcomp02, 2026-11-25, MMseqs2 18.8cc5c, through
`scripts/pod_pytest.sh`). The tests of that difference are in
tests/r2/test_r03_nr90_plain_linclust.py on branch task/R-03, pending the PI's decision. Card C-38c's measurement,
on which DEC-77 rests, ran `l2/cluster.workflow`'s own ten calls (`scripts/stress/c38c_linclust_vs_cascade.py`,
`linclust_redundancy_calls`), not plain `mmseqs linclust`.

What each check here is:

- NR90, all three modes (one workflow): independent oracle, one `mmseqs linclust` call with the level's parameters
  and the three flags the decomposition uses, stated explicitly (`--spaced-kmer-mode 1 --alignment-mode 3
  --comp-bias-corr 0`, the values `mmseqs cluster` gives its own linclust step). Equal on all three inputs: the 10
  calls are linclust's call structure with those flags. On the near-threshold set, where no two proteins are within
  90%, plain `mmseqs linclust` agrees as well.
- NR50, `hybrid-b1`: a single MMseqs2 command reproduces the decomposition, `mmseqs cluster rep90 clu50r tmp
  --min-seq-id 0.5 -c 0.8 --cov-mode 1 --cluster-mode 2 -s 2 --cluster-steps 2`: its `-v 3` log is the 19 calls of
  `workflow("nr50", "hybrid-b1")` argument for argument, and its clusters, on the build's own rep90, are the build's.
- NR50, `hybrid-a`: the default `mmseqs cluster` with the level's parameters (MMseqs2 picks `-s 3` in 3 steps at
  50% identity), the same two checks, all 23 calls. On the fixture and the planted set 1, 2 and 3 cascade steps
  give the same NR50 clusters; on the near-threshold set one step gives other clusters than two, so a dropped or
  added step shows there; two against three steps is told apart only by the call checks.
- NR50, `linclust` mode: one `mmseqs linclust` with the level's parameters and `--spaced-kmer-mode 1
  --alignment-mode 3` (the decomposition's flags at NR50).
- Planted truth (constructed, not an MMseqs2 output): NR90 groups whole and never two in one cluster, NR50 never
  joins unrelated bases, and the cascaded NR50 of `hybrid-a` and `hybrid-b1` keeps each planted family whole.

Not checked: that every member meets the level's identity and coverage against its final representative. MMseqs2's
own cascade does not promise it (`mmseqs cluster -h`: "Cascaded clustering can cluster sequence that do not fulfill
the clustering criteria", unless `--cluster-reassign`), and linclust's `mergeclusters` of its pre-clustering and its
clustering is transitive too; since a single command reproduces each NR50 decomposition, equality with that command
is the stronger check.

Needs `mmseqs` on PATH (the bio-box image); skipped otherwise (GitHub Actions).
"""

from __future__ import annotations

import hashlib
import itertools
import random
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import ids, io
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l2 import cluster
from biobox.registry.schemas import table

sys.path.insert(0, str(Path(__file__).parent))
import test_cluster_split as b41  # noqa: E402  (the micro-S0 fixture's proteins, as card B-41 builds them)

pytestmark = [pytest.mark.r2, pytest.mark.covers("R-03", "B-71", "DEC-77"),
              pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")]

REL = "rR03"
THREADS = 2
MODES = ("hybrid-b1", "hybrid-a", "linclust")
INPUTS = ("fixture", "planted", "remote")
# a level's parameters as `cluster.LEVELS` and `cluster.COMMON` name them, and nothing else
LEVEL_ARGS = {lvl: ["--min-seq-id", str(p["min_seq_id"]), "-c", str(p["coverage"]), *cluster.COMMON]
              for lvl, p in cluster.LEVELS.items()}
# the one MMseqs2 command each mode's NR50 step is checked against, on the same rep90 database
NR50_ORACLE = {
    "hybrid-b1": ("cluster", "-s", "2", "--cluster-steps", "2"),
    "hybrid-a": ("cluster",),
}
# The flags whose values `workflow()` takes from the linclust step inside `mmseqs cluster` (its `-v 3` log) where
# plain `mmseqs linclust` with the same level parameters runs `--spaced-kmer-mode 0 --alignment-mode 2
# --comp-bias-corr 1` (NR90) and `--spaced-kmer-mode 0 --alignment-mode 2` (NR50).
DECOMPOSITION_FLAGS = {"nr90": ["--spaced-kmer-mode", "1", "--alignment-mode", "3", "--comp-bias-corr", "0"],
                       "nr50": ["--spaced-kmer-mode", "1", "--alignment-mode", "3"]}
AA = "ACDEFGHIKLMNPQRSTVWY"


# --- inputs ---------------------------------------------------------------------------------------------


def planted_root(root: Path, bases: int = 60, seed: int = 20261125) -> tuple[Path, dict[str, tuple[int, int]]]:
    """A planted set, written as stage S0's `protein_seq`: `bases` unrelated random proteins (hubs) of 120 to 300
    residues; for each, 3 variants at 20% substitution (80% identity to the hub, about 64% to each other) with their
    last 3 residues cut (so the hub is the longest of its family); the hub and each variant are a group, with 3 copies
    each at 2% substitution (about 96 to 98% identity inside a group). Returns the root and the truth,
    protein uid hex -> (base, group)."""
    rng = random.Random(seed)

    def mutate(s: str, rate: float) -> str:
        return "".join(rng.choice(AA.replace(c, "")) if rng.random() < rate else c for c in s)

    truth: dict[str, tuple[int, int]] = {}
    seqs: dict[str, str] = {}
    for b in range(bases):
        hub = "M" + "".join(rng.choice(AA) for _ in range(rng.randint(119, 299)))
        heads = [hub] + [mutate(hub, 0.20)[:-3] for _ in range(3)]
        for g, head in enumerate(heads):
            for s in [head] + [mutate(head, 0.02) for _ in range(3)]:
                h = ids.protein_uid(s).hex()
                truth.setdefault(h, (b, g))
                seqs[h] = s
    return write_proteins(root, seqs.values()), truth


def remote_root(root: Path, bases: int = 100, seed: int = 2) -> Path:
    """A set near the NR50 threshold, written as stage S0's `protein_seq`: `bases` random hubs of 150 to 300
    residues, each with 6 variants at 40 to 50% substitution (50 to 60% identity to the hub) cut by 3 residues. Here
    the cascade's second step changes NR50 (measured on wilkcomp02, 2026-11-25, MMseqs2 18.8cc5c, through
    `scripts/pod_pytest.sh`: 700 proteins; `mmseqs cluster` at the level's parameters gives 166 clusters with `-s 1
    --cluster-steps 1`, 163 with `-s 2 --cluster-steps 2` and 163 at its default), so the NR50 checks see a step
    dropped or added."""
    rng = random.Random(seed)

    def mutate(s: str, rate: float) -> str:
        return "".join(rng.choice(AA.replace(c, "")) if rng.random() < rate else c for c in s)

    seqs = []
    for _ in range(bases):
        hub = "M" + "".join(rng.choice(AA) for _ in range(rng.randint(149, 299)))
        seqs += [hub] + [mutate(hub, r)[:-3] for r in (0.40, 0.42, 0.44, 0.46, 0.48, 0.50)]
    return write_proteins(root, seqs)


def write_proteins(root: Path, seqs) -> Path:
    """`seqs` as stage S0's `protein_seq` under `root`."""
    rows = sorted({(ids.protein_uid(s), s) for s in seqs})
    t = pa.table({"protein_uid": [u for u, _ in rows],
                  "protein_digest": [hashlib.sha256(s.encode()).digest() for _, s in rows],
                  "length": [len(s) for _, s in rows], "seq": [s for _, s in rows]})
    io.write_table(t, table("protein_seq"), root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0/"
                                                  "protein_prefix=all/part-0.parquet")
    return root


def proteins(root: Path) -> dict[str, str]:
    """The distinct proteins of stage S0, uid hex -> sequence, read here from the files (not by `cluster`)."""
    out: dict[str, str] = {}
    for f in sorted(root.glob(f"l0seq/protein_seq/gene_call_version={V}/stage=S0/protein_prefix=*/part-*.parquet")):
        t = pq.read_table(f, columns=["protein_uid", "seq"])
        for u, s in zip(t.column(0).to_pylist(), t.column(1).to_pylist()):
            out[u.hex()] = s
    return out


# --- MMseqs2 as a user runs it --------------------------------------------------------------------------


def mm(w: Path, *args: str, threads: bool = True) -> None:
    cmd = ["mmseqs", *args, "-v", "1", *(["--threads", str(THREADS)] if threads else [])]
    r = subprocess.run(cmd, cwd=w, capture_output=True, text=True)
    assert r.returncode == 0, f"{shlex.join(cmd)}\n{r.stdout}\n{r.stderr}"


def mm_log(w: Path, *args: str) -> str:
    """One command at `-v 3` on one thread: its log, in which MMseqs2 prints every module call it makes."""
    cmd = ["mmseqs", *args, "-v", "3", "--threads", "1"]
    r = subprocess.run(cmd, cwd=w, capture_output=True, text=True)
    assert r.returncode == 0, f"{shlex.join(cmd)}\n{r.stdout}\n{r.stderr}"
    return r.stdout + r.stderr


def oracle_db(w: Path, prot: dict[str, str]) -> Path:
    """The proteins as a FASTA sorted by uid, `>{uid hex}`, and its database made as the builder makes its own (the
    input, not the clustering: the same database is what makes the two clusterings comparable)."""
    w.mkdir(parents=True, exist_ok=True)
    with open(w / "in.faa", "w") as fh:
        for h in sorted(prot):
            fh.write(f">{h}\n{prot[h]}\n")
    mm(w, "createdb", "in.faa", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "1", threads=False)
    return w


def tsv_sets(w: Path, db: str, clu: str) -> dict[str, frozenset[str]]:
    """A clustering database as {representative: member set}, from `createtsv`."""
    mm(w, "createtsv", db, db, clu, f"{clu}.tsv")
    out: dict[str, set[str]] = {}
    for line in (w / f"{clu}.tsv").read_text().splitlines():
        if line:
            rep, mem = line.split("\t")
            out.setdefault(rep, set()).add(mem)
    return {r: frozenset(m) for r, m in out.items()}


def built_sets(out: Path) -> dict[str, dict[str, frozenset[str]]]:
    """The build's `cluster_member` as {level: {representative: member set}}."""
    by: dict[str, dict[str, set[str]]] = {}
    for f in sorted((out / f"l2/cluster_member/release={REL}").glob("part-*.parquet")):
        t = pq.read_table(f, columns=["member_uid", "level", "rep_uid"])
        for m, lv, r in zip(*(t.column(i).to_pylist() for i in range(3))):
            by.setdefault(lv, {}).setdefault(r.hex(), set()).add(m.hex())
    return {lv: {r: frozenset(m) for r, m in d.items()} for lv, d in by.items()}


def same_clustering(got: dict[str, frozenset[str]], want: dict[str, frozenset[str]], what: str) -> None:
    """Fail unless the two clusterings have the same member sets (representatives and order are not compared)."""
    g, x = set(got.values()), set(want.values())
    if g == x:
        return
    only_g, only_x = g - x, x - g
    moved = set().union(*only_g) if only_g else set()
    where = {m: c for c in x for m in c}
    lines = [f"{what}: the build has {len(g)} clusters, the oracle {len(x)}; {len(only_g)} of the build's member "
             f"sets are not the oracle's and {len(only_x)} of the oracle's are not the build's, covering "
             f"{len(moved)} proteins"]
    for c in sorted(only_g, key=lambda c: (-len(c), min(c)))[:6]:
        parts = sorted({where[m] for m in c}, key=lambda s: (-len(s), min(s)))
        lines.append(f"  build cluster of {len(c)} (e.g. {min(c)}) -> oracle clusters "
                     f"{[(len(p & c), len(p)) for p in parts]} (shared, size)")
    pytest.fail("\n".join(lines))


def linclust_nr90(prot: dict[str, str], w: Path, *flags: str) -> dict[str, frozenset[str]]:
    """One `mmseqs linclust` call at NR90 with the level's parameters and `flags`, on the proteins' database."""
    w = oracle_db(w, prot)
    mm(w, "linclust", "db", "clu90", "tmp90", *LEVEL_ARGS["nr90"], *flags)
    return tsv_sets(w, "db", "clu90")


# --- the builds and their oracles -------------------------------------------------------------------------


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """For each input: its proteins and the build of each mode (`cluster.run`, one process, `THREADS`) as member
    sets per level."""
    base = tmp_path_factory.mktemp("r03")
    roots = {"fixture": b41.fixture_root(base / "fixture" / "root")}
    roots["planted"], truth = planted_root(base / "planted" / "root")
    roots["remote"] = remote_root(base / "remote" / "root")
    out = {}
    for name, root in roots.items():
        d = {"prot": proteins(root), "runs": {}, "base": base / name}
        for mode in MODES:
            (base / name / f"tmp-{mode}").mkdir()
            cluster.run(root, REL, V, ["S0"], threads=THREADS, tmp=base / name / f"tmp-{mode}",
                        out=base / name / f"out-{mode}", mode=mode)
            d["runs"][mode] = built_sets(base / name / f"out-{mode}")
        out[name] = d
    out["planted"]["truth"] = truth
    return out


# --- NR90: one `mmseqs linclust` call with the decomposition's flags ------------------------------------------


@pytest.mark.parametrize("name", INPUTS)
def test_r03_nr90_decomposition_equals_mmseqs_linclust_with_its_three_flags(built, name, tmp_path):
    """Oracle (independent implementation): `mmseqs linclust db clu90 tmp --min-seq-id 0.9 -c 0.8 --cov-mode 1
    --cluster-mode 2 --spaced-kmer-mode 1 --alignment-mode 3 --comp-bias-corr 0`, the level's parameters with the
    three flags the decomposition uses stated explicitly, on the same database. Equal in all three modes (one
    workflow): the 10 calls are linclust's call structure with those flags."""
    b = built[name]
    want = linclust_nr90(b["prot"], tmp_path, *DECOMPOSITION_FLAGS["nr90"])
    for mode in MODES:
        assert sum(map(len, b["runs"][mode]["nr90"].values())) == len(b["prot"])
        same_clustering(b["runs"][mode]["nr90"], want,
                        f"{name} {mode} NR90 against mmseqs linclust {shlex.join(DECOMPOSITION_FLAGS['nr90'])}")


@pytest.mark.parametrize("mode", MODES)
def test_r03_nr90_on_the_near_threshold_set_equals_plain_mmseqs_linclust(built, mode, tmp_path):
    """Oracle (independent implementation): plain `mmseqs linclust` with the level's parameters only, on the
    near-threshold set, where no two proteins are within 90% identity: there the three flags change nothing."""
    b = built["remote"]
    same_clustering(b["runs"][mode]["nr90"], linclust_nr90(b["prot"], tmp_path),
                    f"remote {mode} NR90 against plain mmseqs linclust")


# --- NR50: one MMseqs2 command on the same rep90 ----------------------------------------------------------


@pytest.mark.parametrize("mode", sorted(NR50_ORACLE))
@pytest.mark.parametrize("name", INPUTS)
def test_r03_nr50_equals_one_mmseqs_command_on_the_builds_rep90(built, name, mode, tmp_path):
    """Oracle (independent implementation): `NR50_ORACLE[mode]` with `LEVELS["nr50"]` and `COMMON`, run on the
    build's own NR90 representatives (a `createsubdb` of the same database by their keys), so NR50 is checked apart
    from NR90. Compared twice: on the representatives (the NR50 step's own output), and on every protein, each NR50
    cluster of the oracle taken with the build's NR90 clusters of its members (the build's merge of the levels)."""
    nr50_check(built[name], built[name]["runs"][mode], NR50_ORACLE[mode], tmp_path, f"{name} {mode}")


@pytest.mark.parametrize("name", INPUTS)
def test_r03_nr50_of_linclust_mode_equals_mmseqs_linclust_with_its_two_flags(built, name, tmp_path):
    """Oracle (independent implementation): the `linclust` mode's NR50 against one `mmseqs linclust` call with the
    level's parameters and `--spaced-kmer-mode 1 --alignment-mode 3`, the decomposition's flags at NR50."""
    nr50_check(built[name], built[name]["runs"]["linclust"], ("linclust", *DECOMPOSITION_FLAGS["nr50"]), tmp_path,
               f"{name} linclust")


def nr50_check(b: dict, run: dict, command: tuple[str, ...], w: Path, label: str) -> None:
    w = oracle_db(w, b["prot"])
    key = {}
    for line in (w / "db.lookup").read_text().splitlines():
        k, h, _ = line.split("\t")
        key[h] = int(k)
    reps90 = set(run["nr90"])
    (w / "reps.txt").write_text("".join(f"{k}\n" for k in sorted(key[h] for h in reps90)))
    mm(w, "createsubdb", "reps.txt", "db", "rep90", threads=False)
    cmd, *extra = command
    mm(w, cmd, "rep90", "clu50r", "tmp50", *LEVEL_ARGS["nr50"], *extra)
    want = tsv_sets(w, "rep90", "clu50r")
    got_reps = {r: c & reps90 for r, c in run["nr50"].items()}
    what = f"{label} NR50 against mmseqs {shlex.join(command)}"
    same_clustering(got_reps, want, what + " (rep90)")
    full = {r: frozenset().union(*(run["nr90"][m] for m in c)) for r, c in want.items()}
    same_clustering(run["nr50"], full, what + " (every protein)")


# --- the module calls MMseqs2 itself makes -----------------------------------------------------------------

MODULES = {"kmermatcher", "rescorediagonal", "clust", "createsubdb", "filterdb", "prefilter", "align",
           "mergeclusters"}


def logged_calls(log: str, tmp: str, linclust_only: bool) -> list[list[str]]:
    """The module calls in a `-v 3` log (a line equal to the one before is one call: `kmermatcher` prints its own
    twice), without `--threads` and `-v`, with MMseqs2's temporary directories renamed as `workflow()` names them:
    `<tmp>/<n>/linclust/<m>` (inside `cluster`) and `<tmp>/<n>` (`linclust` alone) as `<tmp>/w/linclust/l`, and
    `<tmp>/<n>` inside `cluster` as `<tmp>/w`. Only paths are renamed; every parameter is compared as logged."""
    out, last = [], None
    for line in log.splitlines():
        if line.split(" ", 1)[0] not in MODULES:
            continue
        if line != last:
            s = re.sub(rf"\b{tmp}/\d+/linclust/\d+", f"{tmp}/w/linclust/l", line)
            s = re.sub(rf"\b{tmp}/\d+", f"{tmp}/w/linclust/l" if linclust_only else f"{tmp}/w", s)
            argv, keep, i = shlex.split(s), [], 0
            while i < len(argv):
                if argv[i] in ("--threads", "-v"):
                    i += 2
                    continue
                keep.append(argv[i])
                i += 1
            out.append(keep)
        last = line
    return out


def split_argv(argv: list[str]) -> tuple[list[str], dict[str, str]]:
    """(the module and its positional arguments, {flag: value}): every flag of these calls takes one value, and a
    value may itself start with '-' (`--realign-score-bias -0.2`)."""
    pos = [argv[0], *itertools.takewhile(lambda a: not a.startswith("-"), argv[1:])]
    rest = argv[len(pos):]
    return pos, dict(zip(rest[::2], rest[1::2]))


def call_check(prot: dict[str, str], w: Path, level: str, mode: str, command: tuple[str, ...]) -> None:
    """`workflow(level, mode)` against the module calls `command` logs at `-v 3` on `prot`, call for call and
    argument for argument; a difference names the flags."""
    w = oracle_db(w, prot)
    mm(w, "linclust", "db", "clu90", "tmp90", *LEVEL_ARGS["nr90"])
    mm(w, "createsubdb", "clu90", "db", "rep90", threads=False)
    cmd, *extra = command
    db, out, tmp = ("db", "clu90", "tmp90") if level == "nr90" else ("rep90", "clu50r", "tmp50")
    log = mm_log(w, cmd, db, f"{out}-log", tmp, *LEVEL_ARGS[level], *extra)
    got = [[a if a != f"{out}-log" else out for a in c] for c in logged_calls(log, tmp, cmd == "linclust")]
    want = [list(c) for c in cluster.workflow(level, mode)]
    diffs = []
    for i, (g, x) in enumerate(zip(got, want)):
        if g != x:
            (pg, fg), (px, fx) = split_argv(g), split_argv(x)
            flags = [f"{f} {fg.get(f, 'absent')} (workflow: {fx.get(f, 'absent')})"
                     for f in sorted(set(fg) | set(fx)) if fg.get(f) != fx.get(f)]
            diffs.append(f"call {i} ({g[0]}): {'positional arguments differ; ' if pg != px else ''}"
                         f"flags that differ: {', '.join(flags)}\n  mmseqs {cmd} logs  {shlex.join(g)}\n"
                         f"  workflow() runs {shlex.join(x)}")
    assert len(got) == len(want), f"mmseqs {cmd} logs {len(got)} module calls, workflow({level!r}, {mode!r}) {len(want)}"
    assert not diffs, f"{len(diffs)} of {len(want)} calls differ\n" + "\n".join(diffs)


@pytest.mark.parametrize("mode", sorted(NR50_ORACLE))
def test_r03_the_nr50_decomposition_is_the_module_calls_of_one_mmseqs_cluster_call(built, tmp_path, mode):
    """`workflow("nr50", mode)` against the module calls `NR50_ORACLE[mode]` logs at `-v 3` on the fixture, call for
    call and argument for argument (19 calls for `hybrid-b1`, 23 for `hybrid-a`)."""
    call_check(built["fixture"]["prot"], tmp_path, "nr50", mode, NR50_ORACLE[mode])


# --- planted truth (constructed, not an MMseqs2 output) ------------------------------------------------------


@pytest.mark.parametrize("mode", MODES)
def test_r03_planted_families_cluster_as_constructed(built, mode):
    """Constructed truth on the planted set: at NR90 each group (96 to 98% identity) lies in one cluster and no
    cluster holds two groups (at most 80% apart, under the 90% threshold, with nothing between them); at NR50 no
    cluster holds two bases (unrelated sequences); under the cascaded NR50 of `hybrid-a` and `hybrid-b1` every base's
    family (each variant 80% identical to the hub, the longest) lies in one cluster."""
    b = built["planted"]
    truth, run = b["truth"], b["runs"][mode]
    assert set(truth) == set(b["prot"])
    groups: dict[tuple[int, int], set[str]] = {}
    for h, t in truth.items():
        groups.setdefault(t, set()).add(h)
    where90 = {m: r for r, c in run["nr90"].items() for m in c}
    split90 = [t for t, g in groups.items() if len({where90[m] for m in g}) > 1]
    mixed90 = [r for r, c in run["nr90"].items() if len({truth[m] for m in c}) > 1]
    mixed50 = [r for r, c in run["nr50"].items() if len({truth[m][0] for m in c}) > 1]
    assert not split90, f"planted NR90 groups split: {split90[:10]}"
    assert not mixed90, f"NR90 clusters holding two planted groups: {len(mixed90)}"
    assert not mixed50, f"NR50 clusters holding two unrelated bases: {len(mixed50)}"
    if mode != "linclust":
        where50 = {m: r for r, c in run["nr50"].items() for m in c}
        bases = {t[0] for t in truth.values()}
        split50 = [x for x in bases if len({where50[m] for m, t in truth.items() if t[0] == x}) > 1]
        assert not split50, f"planted NR50 families split: {split50[:10]}"

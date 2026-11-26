"""Ring R2: an independent oracle for the hybrid clustering modes (card R-03; AGENTS §5.2).

For `hybrid-a`, `hybrid-b1` and `linclust`, `cluster.run` replays `cluster.workflow(level, mode)`, the call list
`run_split` replays too, so the byte identity of tests/r1/test_cluster_hybrid.py and tests/r1/test_cluster_pruning.py
compares two drivers of one decomposition. Here each level is compared with MMseqs2 run as a user runs it, never
through `workflow()`, as sets of member sets (representatives and row order are not compared). There are three
inputs: the micro-S0 fixture's proteins (as tests/r2/test_cluster_split.py builds them, 924 proteins), a planted set
(`planted_root`, seed 20261125) and a set near the NR50 threshold where the cascade's second step changes the
clusters (`remote_root`, seed 2). Every MMseqs2 call uses 2 threads; the database is made as the builder makes its
own (FASTA sorted by uid, `createdb --shuffle 0`), from the proteins read here.

What each check is:

- NR90, all three modes (one workflow): independent oracle, one `mmseqs linclust` call with the parameters
  `LEVELS["nr90"]` and `COMMON` name. **It differs** (measured on wilkcomp02, 2026-11-25, this file through
  `scripts/pod_pytest.sh`): on the fixture the build has 857 clusters and linclust 859 (6 proteins placed otherwise),
  on the planted set 240 and 239 (8 proteins); on the near-threshold set, where no two proteins are within 90%,
  they agree. The cause is three flags: `workflow()` copied the linclust step that
  `mmseqs cluster` runs inside itself, where MMseqs2 sets `kmermatcher --spaced-kmer-mode 1` and `align
  --alignment-mode 3 --comp-bias-corr 0`, while `mmseqs linclust` with the same level parameters runs
  `--spaced-kmer-mode 0`, `--alignment-mode 2` and `--comp-bias-corr 1` (the module-call test names them). The
  smallest inputs that show it are two proteins each (`SMALLEST`): on the fixture's pair `--spaced-kmer-mode` alone
  decides (the build joins a 77-residue fragment that linclust leaves apart); on the planted pair `--alignment-mode`
  alone decides (linclust, judging identity from the score, joins two proteins 85.6% identical at NR90; the build
  keeps them apart). These tests assert equality and are left red: a finding for the PI, not a test to bend.
- NR90 diagnostic: `mmseqs linclust` with those three flags at `workflow()`'s values (`PRESTEP`) gives the build's
  clusters on both inputs: the 10 calls are linclust's own call structure; only the flags differ.
- NR50, `hybrid-b1`: a single MMseqs2 command reproduces the decomposition, `mmseqs cluster rep90 clu50r tmp
  --min-seq-id 0.5 -c 0.8 --cov-mode 1 --cluster-mode 2 -s 2 --cluster-steps 2`: its `-v 3` log is the 19 calls of
  `workflow("nr50", "hybrid-b1")` argument for argument, and its clusters, on the build's own rep90, are the build's.
- NR50, `hybrid-a`: the default `mmseqs cluster` with the level's parameters (MMseqs2 picks `-s 3` in 3 steps at
  50% identity), the same two checks, all 23 calls.
- NR50, `linclust` mode: one `mmseqs linclust` at the level's parameters; it differs by `--spaced-kmer-mode` and
  `--alignment-mode` as NR90 does, and equals the build with those two flags set (diagnostic).
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
    "linclust": ("linclust",),
}
# The flags whose values `workflow()` takes from the linclust step inside `mmseqs cluster` (its `-v 3` log) and
# `mmseqs linclust` sets otherwise (`test_r03_the_decomposition_is_the_module_calls_mmseqs2_makes`). Used only by the
# diagnostic oracles that say which flag a difference comes from.
PRESTEP = {"nr90": ["--spaced-kmer-mode", "1", "--alignment-mode", "3", "--comp-bias-corr", "0"],
           "nr50": ["--spaced-kmer-mode", "1", "--alignment-mode", "3"]}
PRESTEP_VALUE = dict(zip(PRESTEP["nr90"][::2], PRESTEP["nr90"][1::2]))
# The smallest inputs found that show the NR90 difference, and the flag that alone decides each (see the test).
SMALLEST = {
    "fragment": ((
        "MKTPTIPTLLGPDGMTSLREYAGYHGGGSGFGGQLRSWNPPSESVDAALLPNFTRGNARADDLVRNNGYAANAIQLHQDHIVGSFFRLSHRPSWRYLGIGEEEARA"
        "FSREVEAAWKEFAEDDCCCIDVERKRTFTMMIREGVAMHAFNGELFVQATWDTSSSRLFRTQFRMVSPKRISNPNNTGDSRNCRAGVQINDSGAALGYYVSEDGYPG"
        "WMPQKWTWIPRELPGGRASFIHVFEPVEDGQTRGANVFYSVMEQMKMLDTLQNTQLQSAIVKAMYAATIESELDTQSAMDFILGANSQEQRERLTGWIGEIAAYYAAA"
        "PVRLGGAKVPHLMPGDSLNLQTAQDTDNGYSVFEQSLLRYIAAGLGVSYEQLSRNYAQMSYSTARASANESWAYFMGRRKFVASRQASQMFLCWLEEAIVRRVVTLP"
        "SKARFSFQEARSAWGNCDWIGSGRMAIDGLKEVQEAVMLIEAGLSTYEKECAKRGDDYQEIFAQQVRETMERRAAGLKPPAWAAAAFESGLRQSTEEEKSDSRAA",
        "MFLCWQEGAIVRRLVALPSKARFSLQEARSAWGNCDWIGSGRMGIDGLKEVQEAVMLIEAGLSTYEKECAKRGDDYQ",
    ), "--spaced-kmer-mode"),
    "estimated-identity": ((
        "MLFGDNFLMPHRWEAGKGNVSMALGPQFNFRYHKKNPNDFNYVIAYNFMVSQIVTCRQGNLINCWMDCIASKWRHPMVLKCGHGGNSCTFRNNQYARGTMWPWEII"
        "IGMKIDSYSMMINCVRENPHLYDPWWSWESIMIHIIMGFR",
        "MLFGGNFLMPHNWEAGKGNMGMALGPNFNFRTHKKYPNDFNYVIAYNFMVMRIVTCRQCNTINCSHDCIASPWPHPMVLKCGHGGHSCVFRNNQYARGTMWPWEII"
        "IGGKIESYSMMIMCVRENPHLYDPWWSWESIMIHIIMGWRPHG",
    ), "--alignment-mode"),
}
SMALLEST_UIDS = {"fragment": ["265373fae12d768f7aabeccf9b45d634", "a024462b0e6be4a34fe85622fcbf57a8"],
                 "estimated-identity": ["7f01129eb212025368fa9abccbc3a635", "d5b91b9961a08b00513b3e7714bb2ef3"]}
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
    the cascade's second step changes NR50 (measured on wilkcomp02, 2026-11-25, through `scripts/pod_pytest.sh`: 700
    proteins; `mmseqs cluster` at the level's parameters gives 166 clusters with `-s 1 --cluster-steps 1`, 163 with
    `-s 2 --cluster-steps 2` and 163 at its default), so the NR50 checks see a step dropped or added."""
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


# --- the builds and their oracles -------------------------------------------------------------------------


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """For each input: its proteins, the build of each mode (`cluster.run`, one process, `THREADS`) as member sets
    per level, and the oracle's NR90: one `mmseqs linclust` call with the level's parameters."""
    base = tmp_path_factory.mktemp("r03")
    roots = {"fixture": b41.fixture_root(base / "fixture" / "root")}
    roots["planted"], truth = planted_root(base / "planted" / "root")
    roots["remote"] = remote_root(base / "remote" / "root")
    out = {}
    for name, root in roots.items():
        d = {"prot": proteins(root), "runs": {}}
        for mode in MODES:
            (base / name / f"tmp-{mode}").mkdir()
            cluster.run(root, REL, V, ["S0"], threads=THREADS, tmp=base / name / f"tmp-{mode}",
                        out=base / name / f"out-{mode}", mode=mode)
            d["runs"][mode] = built_sets(base / name / f"out-{mode}")
        w = oracle_db(base / name / "oracle", d["prot"])
        mm(w, "linclust", "db", "clu90", "tmp90", *LEVEL_ARGS["nr90"])
        d["nr90"] = tsv_sets(w, "db", "clu90")
        mm(w, "linclust", "db", "clu90p", "tmp90p", *LEVEL_ARGS["nr90"], *PRESTEP["nr90"])
        d["nr90_prestep"] = tsv_sets(w, "db", "clu90p")
        out[name] = d
    out["planted"]["truth"] = truth
    return out


# --- NR90: one `mmseqs linclust` call ---------------------------------------------------------------------


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("name", INPUTS)
def test_r03_nr90_equals_one_mmseqs_linclust_call_with_the_levels_parameters(built, name, mode):
    """Oracle (independent implementation): `mmseqs linclust db clu90 tmp --min-seq-id 0.9 -c 0.8 --cov-mode 1
    --cluster-mode 2`, the parameters `LEVELS["nr90"]` and `COMMON` name, on the same database."""
    b = built[name]
    assert sum(map(len, b["runs"][mode]["nr90"].values())) == len(b["prot"])
    same_clustering(b["runs"][mode]["nr90"], b["nr90"], f"{name} {mode} NR90 against mmseqs linclust")


@pytest.mark.parametrize("name", INPUTS)
def test_r03_nr90_equals_mmseqs_linclust_with_the_flags_mmseqs_cluster_gives_its_prestep(built, name):
    """Diagnostic oracle (one parameter set changed on purpose): `mmseqs linclust` with the level's parameters and
    `PRESTEP["nr90"]`, the three flags whose values `mmseqs cluster` gives its own linclust step. Equal: the 10 calls
    are linclust's call structure, and the difference of the test above is those three flags and nothing else."""
    b = built[name]
    for mode in MODES:
        same_clustering(b["runs"][mode]["nr90"], b["nr90_prestep"],
                        f"{name} {mode} NR90 against mmseqs linclust {shlex.join(PRESTEP['nr90'])}")


@pytest.mark.parametrize("pair", sorted(SMALLEST))
def test_r03_smallest_input_that_shows_the_nr90_difference(pair, tmp_path):
    """The smallest inputs found (greedy deletion from the differing clusters of each set, on the pod, 2026-11-25):
    two proteins each. `fragment` (the fixture): the build joins a 77-residue fragment to a 533-residue protein,
    `mmseqs linclust` leaves them apart; `--spaced-kmer-mode` alone decides it. `estimated-identity` (the planted
    set): a hub copy and a variant of one planted family, substitutions only, 146 aligned positions of which 85.6%
    are identical (computed here); `mmseqs linclust` puts them in one NR90 cluster, its `--alignment-mode 2` judging
    identity from the score, and the build, at `--alignment-mode 3`, keeps them apart; `--alignment-mode` alone
    decides it. Oracle: one `mmseqs linclust` call with the level's parameters."""
    seqs, flag = SMALLEST[pair]
    assert [ids.protein_uid(s).hex() for s in seqs] == SMALLEST_UIDS[pair]  # as found in the sets above
    if pair == "estimated-identity":
        a, c = seqs
        assert sum(x == y for x, y in zip(a, c)) / min(len(a), len(c)) < cluster.LEVELS["nr90"]["min_seq_id"]
    root = write_proteins(tmp_path / "root", seqs)
    (tmp_path / "tmp").mkdir()
    cluster.run(root, REL, V, ["S0"], threads=THREADS, tmp=tmp_path / "tmp", out=tmp_path / "out", mode="hybrid-b1")
    got = built_sets(tmp_path / "out")["nr90"]
    w = oracle_db(tmp_path / "oracle", proteins(root))
    mm(w, "linclust", "db", "clu90", "tmp90", *LEVEL_ARGS["nr90"])
    mm(w, "linclust", "db", "clu90f", "tmp90f", *LEVEL_ARGS["nr90"], flag, PRESTEP_VALUE[flag])
    # the one flag at the value `workflow()` gives it makes linclust agree: it is the flag that decides
    same_clustering(got, tsv_sets(w, "db", "clu90f"), f"{pair}: NR90 against mmseqs linclust {flag} "
                                                      f"{PRESTEP_VALUE[flag]}")
    same_clustering(got, tsv_sets(w, "db", "clu90"), f"{pair}: NR90 against mmseqs linclust")


# --- NR50: one MMseqs2 command on the same rep90 ----------------------------------------------------------


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("name", INPUTS)
def test_r03_nr50_equals_one_mmseqs_command_on_the_builds_rep90(built, name, mode, tmp_path):
    """Oracle (independent implementation): `NR50_ORACLE[mode]` with `LEVELS["nr50"]` and `COMMON`, run on the
    build's own NR90 representatives (a `createsubdb` of the same database by their keys), so NR50 is checked apart
    from NR90. Compared twice: on the representatives (the NR50 step's own output), and on every protein, each NR50
    cluster of the oracle taken with the build's NR90 clusters of its members (the build's merge of the levels)."""
    nr50_check(built[name], built[name]["runs"][mode], NR50_ORACLE[mode], tmp_path, f"{name} {mode}")


@pytest.mark.parametrize("name", INPUTS)
def test_r03_nr50_of_linclust_mode_equals_mmseqs_linclust_with_the_prestep_flags(built, name, tmp_path):
    """Diagnostic oracle (one parameter set changed on purpose): the `linclust` mode's NR50 against `mmseqs linclust`
    with `PRESTEP["nr50"]` added. Equal: the 10 calls are linclust's own call structure, and the difference the test
    above reports is those flags alone."""
    nr50_check(built[name], built[name]["runs"]["linclust"], ("linclust", *PRESTEP["nr50"]), tmp_path,
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


CHECKS = {  # (level, mode): the command whose module calls workflow(level, mode) replays
    ("nr90", "hybrid-b1"): ("linclust",),
    ("nr50", "linclust"): NR50_ORACLE["linclust"],
    ("nr50", "hybrid-b1"): NR50_ORACLE["hybrid-b1"],
    ("nr50", "hybrid-a"): NR50_ORACLE["hybrid-a"],
}


@pytest.mark.parametrize("level,mode", sorted(CHECKS))
def test_r03_the_decomposition_is_the_module_calls_mmseqs2_makes(built, tmp_path, level, mode):
    """`workflow(level, mode)` against the module calls the one command (`CHECKS`) logs at `-v 3` on the fixture,
    call for call and argument for argument. NR90 is the same 10 calls in all three modes, so `hybrid-b1` stands for
    them. This names the parameter where a clustering differs; the clusterings above are the evidence."""
    w = oracle_db(tmp_path, built["fixture"]["prot"])
    mm(w, "linclust", "db", "clu90", "tmp90", *LEVEL_ARGS["nr90"])
    mm(w, "createsubdb", "clu90", "db", "rep90", threads=False)
    cmd, *extra = CHECKS[(level, mode)]
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


def split_argv(argv: list[str]) -> tuple[list[str], dict[str, str]]:
    """(the module and its positional arguments, {flag: value}): every flag of these calls takes one value, and a
    value may itself start with '-' (`--realign-score-bias -0.2`)."""
    pos = [argv[0], *itertools.takewhile(lambda a: not a.startswith("-"), argv[1:])]
    rest = argv[len(pos):]
    return pos, dict(zip(rest[::2], rest[1::2]))


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

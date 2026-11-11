"""Card B-37 experiment (design only): MMseqs2 clustering as the release builder runs it, and the same
clustering decomposed into its modules with the per-query modules split into independent processes.

    python b37_split_cluster.py export ROOT STAGES OUT.faa [--limit N]
    python b37_split_cluster.py synth OUT.faa [--seed 37]
    python b37_split_cluster.py workflow IN.faa WORK --threads T
    python b37_split_cluster.py split IN.faa WORK --splits N --threads T --plan-from WORKFLOW_WORK [--self-matches]
    python b37_split_cluster.py compare WORK_A WORK_B

`export` writes the distinct proteins of the stages' `protein_seq` files as `biobox.l2.cluster` writes its
FASTA (`>{uid hex}`, sorted by protein_uid), the first N only with `--limit`. `synth` writes a constructed set
with many clusters and many length ties (below).

`workflow` runs exactly the MMseqs2 calls of `biobox.l2.cluster.run` (createdb, `cluster` at NR90,
createsubdb, `cluster` at NR50 over the NR90 representatives, mergeclusters, then per level createtsv, align,
convertalis) at `-v 3`, each top-level command its own process timed with `os.wait4`.

`split` replays the module calls that the two `cluster` workflow calls made, read from that run's own `-v 3`
log (`--plan-from`), with the same arguments and paths, one process per module call, plus the one step the
workflow's `linclust.sh` does in shell (`order_redundancy`, the first column of `pre_clust.index`). Every
per-query module (`prefilter`, `rescorediagonal`, `align`, and the builder's own member `align`) runs as N
processes: `prefilter` over N consecutive key ranges of its query database (target: the whole database),
`rescorediagonal` and `align` over N consecutive key ranges of their input result database (query and
target: the whole databases). Ranges are a pure function of the sorted key list and N. The parts are joined
in key order with `concatdbs --preserve-keys 1`. `kmermatcher`, `clust`, `filterdb`, `createsubdb` and
`mergeclusters` run once, as in the workflow.

`compare` checks every MMseqs2 database the two runs share (same relative path) entry by entry (key and
bytes; the data files' layout may differ) and the final TSVs as sorted rows and as bytes.

Every number this script prints is a measurement of the host it runs on. It writes only under WORK. It is
not a test and is not run by pytest.

`synth` (seed fixed): 1,500 random base proteins of 80 to 299 residues; each gives 5 variants at 25%
substitution (about 75% identity: one NR50 cluster, several NR90 clusters) and each variant 4 copies at 4%
substitution (one NR90 cluster). 30,000 proteins; substitutions only, so lengths tie within a family, which is
where a greedy clustering by length could depend on order. Records are named by the hex of the BLAKE3-128 of
the sequence and sorted by it, as the builder names proteins.
"""

from __future__ import annotations

import argparse
import heapq
import json
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

MMSEQS = "mmseqs"
COMMON = ["--cov-mode", "1", "--cluster-mode", "2"]
ALN_FORMAT = "query,target,qstart,qend,tstart,tend,fident,alnlen,cigar"
# the modules the `cluster` workflow ran (its -v 3 log on the fixture, MMseqs2 18.8cc5c)
MODULES = {"kmermatcher", "rescorediagonal", "clust", "createsubdb", "filterdb", "prefilter", "align", "mergeclusters"}
# per-query modules: which positional argument is split (the query db of prefilter, the result db of the others)
# and which is the output
SPLIT = {"prefilter": (0, 2), "rescorediagonal": (2, 3), "align": (2, 3)}
SELF_MATCHES = False  # `split --self-matches`: a prefilter part runs with --add-self-matches 1


class Runner:
    """Runs one MMseqs2 command per process in `work`, logging its command line and output to
    `mmseqs.log`, and its wall and CPU seconds (`os.wait4` of the child, which includes the children
    it reaped) to `times.jsonl`."""

    def __init__(self, work: Path, verbosity: int = 3):
        self.work, self.v = work, verbosity
        self.log = open(work / "mmseqs.log", "a")
        self.times = open(work / "times.jsonl", "a")

    def __call__(self, label: str, *args: str, threads: int | None = None) -> dict:
        cmd = [MMSEQS, *args, "-v", str(self.v)]
        if threads is not None:
            cmd += ["--threads", str(threads)]
        self.log.write("### " + label + ": " + shlex.join(cmd) + "\n")
        self.log.flush()
        t = time.perf_counter()
        p = subprocess.Popen(cmd, stdout=self.log, stderr=subprocess.STDOUT, cwd=self.work)
        _, status, ru = os.wait4(p.pid, 0)
        wall = time.perf_counter() - t
        code = os.waitstatus_to_exitcode(status)
        rec = {"label": label, "module": args[0], "wall_s": round(wall, 3),
               "cpu_s": round(ru.ru_utime + ru.ru_stime, 3), "maxrss_kb": ru.ru_maxrss, "exit": code}
        self.times.write(json.dumps(rec) + "\n")
        self.times.flush()
        if code != 0:
            raise RuntimeError(f"{label}: {shlex.join(cmd)} exited {code}")
        return rec


def _fresh(work: Path) -> Path:
    if work.exists():
        sys.exit(f"{work}: exists; each run writes into a fresh directory")
    work.mkdir(parents=True)
    return work


# --- inputs ---------------------------------------------------------------------------------------------


def export(root: Path, stages: list[str], out: Path, limit: int | None) -> int:
    from biobox.l1.genecall import GENE_CALL_VERSION
    from biobox.l2 import cluster

    files = cluster._protein_files(root, GENE_CALL_VERSION, stages)
    n, last = 0, None
    with open(out, "w") as fh:
        for u, s in heapq.merge(*(cluster._stage_proteins(f) for f in files), key=lambda r: r[0]):
            if u == last:
                continue
            last = u
            fh.write(f">{u.hex()}\n{s}\n")
            n += 1
            if limit is not None and n >= limit:
                break
    return n


def synth(out: Path, seed: int) -> int:
    from blake3 import blake3

    rng = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def mutate(s: str, rate: float) -> str:
        return "".join(rng.choice(aa.replace(c, "")) if rng.random() < rate else c for c in s)

    recs = {}
    for _ in range(1500):
        base = "M" + "".join(rng.choice(aa) for _ in range(rng.randint(79, 298)))
        for _ in range(5):
            var = mutate(base, 0.25)
            for _ in range(4):
                s = mutate(var, 0.04)
                recs[blake3(s.encode()).digest(length=16).hex()] = s
    with open(out, "w") as fh:
        for h in sorted(recs):
            fh.write(f">{h}\n{recs[h]}\n")
    return len(recs)


# --- the builder's calls, as biobox.l2.cluster.run makes them --------------------------------------------


def workflow(faa: Path, work: Path, threads: int) -> None:
    w = _fresh(work)
    shutil.copy(faa, w / "in.faa")
    mm = Runner(w)
    mm("createdb", "createdb", "in.faa", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "1")
    mm("cluster90", "cluster", "db", "clu90", "tmp90", "--min-seq-id", "0.9", "-c", "0.8", *COMMON, threads=threads)
    mm("createsubdb90", "createsubdb", "clu90", "db", "rep90")
    mm("cluster50", "cluster", "rep90", "clu50r", "tmp50", "--min-seq-id", "0.5", "-c", "0.8", *COMMON,
       threads=threads)
    mm("mergeclusters", "mergeclusters", "db", "clu50", "clu90", "clu50r", threads=threads)
    _finish(mm, threads, 1)


def _member_align(level: str, clu: str) -> list[str]:
    return ["align", "db", "db", clu, f"aln_{level}", "-a", "-e", "inf", "--min-seq-id", "0", "-c", "0",
            "--alignment-mode", "3"]


def _finish(mm: Runner, threads: int, splits: int) -> None:
    """Per level: the members' TSV, the members aligned to their representatives (the builder's `align` over
    the clustering; split by representative range when `splits` > 1), and that alignment's TSV."""
    for level, clu in (("nr90", "clu90"), ("nr50", "clu50")):
        mm(f"createtsv_{level}", "createtsv", "db", "db", clu, f"{clu}.tsv", threads=threads)
        _module(mm, f"align_{level}", _member_align(level, clu), splits, threads)
        mm(f"convertalis_{level}", "convertalis", "db", "db", f"aln_{level}", f"aln_{level}.tsv",
           "--format-output", ALN_FORMAT, threads=threads)


# --- the same clustering, module by module, with splits ---------------------------------------------------


def _index_keys(index: Path) -> list[int]:
    return sorted(int(line.split("\t", 1)[0]) for line in index.read_text().splitlines() if line)


def ranges(keys: list[int], splits: int) -> list[list[int]]:
    """`splits` consecutive runs of the sorted keys, as equal in count as integer division makes them: a pure
    function of the key list and `splits`."""
    n = len(keys)
    return [keys[k * n // splits:(k + 1) * n // splits] for k in range(splits)]


def _module(mm: Runner, label: str, argv: list[str], splits: int, threads: int) -> None:
    """One module call; a per-query module with `splits` > 1 as that many processes over key ranges, joined
    in key order."""
    mod = argv[0]
    if mod not in SPLIT or splits == 1:
        mm(label, *argv, threads=threads)
        return
    pos, out_pos = SPLIT[mod]
    args = argv[1:]
    w = mm.work
    src, out = args[pos], args[out_pos]
    parts = []
    for k, sub in enumerate(ranges(_index_keys(w / f"{src}.index"), splits)):
        ids = Path(f"{out}.b37ids{k}")
        (w / ids).write_text("".join(f"{x}\n" for x in sub))
        part_in = f"{src}.b37in{k}_{Path(out).name}"
        # a sequence db (prefilter's query) is linked, as the workflow links its subsets; a result db is copied
        mm(f"{label}.in{k}", "createsubdb", str(ids), src, part_in, "--subdb-mode", "1" if mod == "prefilter" else "0")
        part_out = f"{out}.b37part{k}"
        a = list(args)
        a[pos], a[out_pos] = part_in, part_out
        if mod == "prefilter" and SELF_MATCHES:
            # the part's query db is not the target db by path, so prefilter does not treat it as a search of
            # a db against itself (measured on the fixture: the self hit's score and, for short proteins, its
            # presence change); ask for the self hit explicitly instead
            i = a.index("--add-self-matches")
            a[i + 1] = "1"
        mm(f"{label}.part{k}", mod, *a, threads=threads)
        parts.append(part_out)
    cur = parts[0]
    for k, p in enumerate(parts[1:], 1):
        nxt = f"{out}.b37cat{k}"
        mm(f"{label}.cat{k}", "concatdbs", cur, p, nxt, "--preserve-keys", "1", threads=1)
        cur = nxt
    mm(f"{label}.mv", "mvdb", cur, out)


_HEAD = re.compile(r"^### (cluster90|cluster50): ")


def plan_from(log: Path) -> dict[str, list[list[str]]]:
    """The module calls each `cluster` workflow call made, in order, from its run's -v 3 log: every line that
    starts with a module name between that call's header and the next header (a module prints its own command
    line; `kmermatcher` prints it twice, so a line equal to the one before is one call)."""
    plans: dict[str, list[list[str]]] = {}
    cur, last = None, None
    for line in log.read_text().splitlines():
        m = _HEAD.match(line)
        if m:
            cur, last = m.group(1), None
            plans[cur] = []
            continue
        if line.startswith("### "):
            cur = None
            continue
        if cur is None:
            continue
        word = line.split(" ", 1)[0]
        if word in MODULES and line != last:
            argv = shlex.split(line)
            # threads and verbosity are set by the replay
            out, i = [], 0
            while i < len(argv):
                if argv[i] in ("--threads", "-v"):
                    i += 2
                    continue
                out.append(argv[i])
                i += 1
            plans[cur].append(out)
        if word in MODULES:
            last = line
    return plans


def _dirs(w: Path, argv: list[str]) -> None:
    for a in argv[1:]:
        if not a.startswith("-") and "/" in a:
            (w / a).parent.mkdir(parents=True, exist_ok=True)


def _replay(mm: Runner, name: str, calls: list[list[str]], splits: int, threads: int) -> None:
    for i, argv in enumerate(calls):
        _dirs(mm.work, argv)
        _module(mm, f"{name}.{i:02d}.{argv[0]}", argv, splits, threads if argv[0] != "createsubdb" else None)
        if argv[0] == "clust" and argv[3].endswith("/pre_clust"):
            # linclust.sh: awk '{ print $1 }' pre_clust.index > order_redundancy
            idx = mm.work / f"{argv[3]}.index"
            (idx.parent / "order_redundancy").write_text(
                "".join(line.split("\t", 1)[0] + "\n" for line in idx.read_text().splitlines()))


def split(faa: Path, work: Path, splits: int, threads: int, plan_dir: Path) -> None:
    plans = plan_from(plan_dir / "mmseqs.log")
    w = _fresh(work)
    (w / "plan.json").write_text(json.dumps(plans, indent=1))
    shutil.copy(faa, w / "in.faa")
    mm = Runner(w)
    mm("createdb", "createdb", "in.faa", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "1")
    _replay(mm, "cluster90", plans["cluster90"], splits, threads)
    mm("createsubdb90", "createsubdb", "clu90", "db", "rep90")
    _replay(mm, "cluster50", plans["cluster50"], splits, threads)
    mm("mergeclusters", "mergeclusters", "db", "clu50", "clu90", "clu50r", threads=threads)
    _finish(mm, threads, splits)


# --- comparison --------------------------------------------------------------------------------------------


def read_db(base: Path) -> dict[int, bytes]:
    """An MMseqs2 database as {key: entry bytes}: the index's (key, offset, length) over the data file, or over
    its numbered parts `.0`, `.1`, ... taken as one run in that order. Each entry must end in a NUL."""
    if base.exists() and base.is_file():
        data = base.read_bytes()
    else:
        parts, k = [], 0
        while (base.parent / f"{base.name}.{k}").exists():
            parts.append((base.parent / f"{base.name}.{k}").read_bytes())
            k += 1
        data = b"".join(parts)
    out = {}
    for line in (base.parent / f"{base.name}.index").read_text().splitlines():
        key, off, ln = (int(x) for x in line.split("\t"))
        e = data[off:off + ln]
        if len(e) != ln or not e.endswith(b"\0"):
            raise RuntimeError(f"{base}: entry {key} is not a NUL-terminated entry at {off}+{ln}")
        out[key] = e
    return out


def _pairs(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text().splitlines():
        rep, mem = line.split("\t")
        if mem in out:
            raise RuntimeError(f"{path}: member {mem} listed twice")
        out[mem] = rep
    return out


def compare(a: Path, b: Path) -> dict:
    dbs = {}
    for t in sorted(a.rglob("*.dbtype")):
        rel = t.relative_to(a).as_posix()[:-len(".dbtype")]
        if ".b37" in rel or "/latest/" in rel or rel.endswith("_h"):
            continue
        if not (b / f"{rel}.dbtype").exists():
            dbs[rel] = "missing in B"
            continue
        da, db_ = read_db(a / rel), read_db(b / rel)
        if da == db_:
            dbs[rel] = f"equal ({len(da)} entries)"
        else:
            diff = sorted(k for k in set(da) | set(db_) if da.get(k) != db_.get(k))
            dbs[rel] = f"DIFFERENT: {len(diff)} of {len(da)} entries, first keys {diff[:5]}"
    res: dict = {"databases": dbs, "all_databases_equal": all(v.startswith("equal") for v in dbs.values())}
    for level, clu in (("nr90", "clu90"), ("nr50", "clu50")):
        pa_, pb = _pairs(a / f"{clu}.tsv"), _pairs(b / f"{clu}.tsv")
        moved = sorted(m for m in pa_ if pb.get(m) != pa_[m])
        sizes: dict[str, int] = {}
        for r in pa_.values():
            sizes[r] = sizes.get(r, 0) + 1
        aln_a = sorted((a / f"aln_{level}.tsv").read_text().splitlines())
        aln_b = sorted((b / f"aln_{level}.tsv").read_text().splitlines())
        res[level] = {"members": len(pa_), "clusters": len(sizes), "clusters_of_2_or_more": sum(v > 1 for v in sizes.values()),
                      "largest": max(sizes.values()), "same_members": set(pa_) == set(pb),
                      "member_rep_pairs_equal": not moved and set(pa_) == set(pb), "moved": len(moved),
                      "moved_examples": [(m, pa_[m], pb.get(m)) for m in moved[:10]],
                      "cluster_tsv_bytes_equal": (a / f"{clu}.tsv").read_bytes() == (b / f"{clu}.tsv").read_bytes(),
                      "aln_rows": len(aln_a), "aln_rows_equal_sorted": aln_a == aln_b,
                      "aln_tsv_bytes_equal": (a / f"aln_{level}.tsv").read_bytes() == (b / f"aln_{level}.tsv").read_bytes()}
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("export")
    e.add_argument("root", type=Path)
    e.add_argument("stages")
    e.add_argument("out", type=Path)
    e.add_argument("--limit", type=int)
    s = sub.add_parser("synth")
    s.add_argument("out", type=Path)
    s.add_argument("--seed", type=int, default=37)
    wf = sub.add_parser("workflow")
    wf.add_argument("faa", type=Path)
    wf.add_argument("work", type=Path)
    wf.add_argument("--threads", type=int, required=True)
    sp = sub.add_parser("split")
    sp.add_argument("faa", type=Path)
    sp.add_argument("work", type=Path)
    sp.add_argument("--splits", type=int, required=True)
    sp.add_argument("--threads", type=int, required=True)
    sp.add_argument("--plan-from", type=Path, required=True)
    sp.add_argument("--self-matches", action="store_true")
    c = sub.add_parser("compare")
    c.add_argument("a", type=Path)
    c.add_argument("b", type=Path)
    a = ap.parse_args()
    if a.cmd == "export":
        print(json.dumps({"proteins": export(a.root, a.stages.split(","), a.out, a.limit)}))
    elif a.cmd == "synth":
        print(json.dumps({"proteins": synth(a.out, a.seed)}))
    elif a.cmd == "workflow":
        workflow(a.faa, a.work, a.threads)
    elif a.cmd == "split":
        global SELF_MATCHES
        SELF_MATCHES = a.self_matches
        split(a.faa, a.work, a.splits, a.threads, a.plan_from)
    else:
        print(json.dumps(compare(a.a, a.b), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

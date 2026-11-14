#!/usr/bin/env python3
"""Card E-13: CRISPR arrays at a minimum of two repeats against three (MinCED's `-minNR`), and how much is false.

Usage:
    python scripts/stress/e13_crispr_min_repeats.py view SRC DST STAGES
    python scripts/stress/e13_crispr_min_repeats.py produce --root VIEW --stage S --cpus N --min-repeats K --tmp DIR
        [--seq-stages S1,P1,P2]
    python scripts/stress/e13_crispr_min_repeats.py null --root ROOT --stages S1,P1,P2 --cpus N --tmp DIR --out DIR
        [--window 1000] [--seed 13] [--batch-bases 25000000]
    python scripts/stress/e13_crispr_min_repeats.py analyse --root ROOT --release rN.3 --view VIEW --null DIR
        --out DIR --tmp DIR [--threads 4] [--library LIB]

Every subcommand runs on the pod, in the runtime image, with the tree's `src` on PYTHONPATH (AGENTS §5.5, §7); none
writes under ROOT, which is bound read-only.

- `view`: a scratch root DST holding one symlink per stage directory of the sequence tables the producer reads
  (`contig_index`, `seq_locator`, `seq_shards`, `genome_input`) and `store_params.json`, as card E-12 made it (the
  pattern of card B-42). The producer's files, ledger and run records land in DST's own directories.
- `produce`: the producer (`biobox.l1.crispr.run`) with a chosen minimum of repeats, as a `bio-box run-step` command:
  it writes into `$BIOBOX_STAGING`, which run-step publishes into the view. A minimum other than 3 is its own producer
  (`minced-<v>-w8-nr<k>-r23x47-s26x50`), so nothing of the default producer's changes.
- `null`: the first null. Every contig of the stages, its sequence shuffled within non-overlapping windows of
  `--window` bases (the last one shorter) by the Altschul-Erickson doublet shuffle (`dishuffle`): within each window
  the multiset of dinucleotides and the first and last base are kept exactly, so the window's base composition and,
  since the bases at every window join are kept too, the contig's dinucleotide counts are kept exactly; repeats
  longer than chance are destroyed. Checked for every window of every contig as it is made (`shuffle_windows`).
  Seeded per contig from `--seed` and the contig_uid (BLAKE2b), so the result depends neither on the order of the
  contigs nor on the processes. MinCED is then run on the shuffled records at a minimum of 2 and of 3, with the
  producer's JVM settings and arguments (`crispr._argv`), each record linear (a circular contig's flanks are not
  added: a chance array across an origin is not counted). Writes `null_summary.json` (calls per Mbp) and
  `null_arrays.tsv`; the coordinates of a shuffled contig are the real contig's, so the analysis can place them
  beside the real genes.
- `analyse`: the arrays found at 2 and not at 3 (`extra`: an array of the nr2 producer that overlaps no array of the
  nr3 producer on its contig), on the release's own tables of `rN.3` (gene calls from `protein_occurrence`, Pfam
  domains from `l1_domains`, contigs from `contig_dim` and `genome_input`; the catalog opened read-only from a copy in
  `--tmp`): within 10 genes and within 10 kb of a gene with a Cas domain (the Cas families of the CRISPR example,
  `CAS_FAMILIES`, taken from `examples/crispr_cas_neighbourhoods/sql/q01_loci_by_domains.sql`), wholly inside a
  called gene, in a genome with no Cas domain at all; their repeat length, spacer length and repeat identity against
  the arrays found at 3; the example's Tn7-like loci (a TniQ gene with a Cas7f gene within 10 genes) that gain an
  array within 10 genes; and the second null: the share of extra calls near a Cas gene expected if each were placed
  at a uniformly random position of its own genome (`placement_share`, formula there). Writes `summary.json` and
  TSV tables into `--out`.

Definitions (both the observed arrays and the nulls use these, `Windows`):
- "within K genes" of a gene C: at most K - 1 called genes lie wholly between the array and C (the gene next to an
  array is 1 gene away; a gene that overlaps the array is 0 away). With the genes of the contig, b_K is the K-th
  largest start among the genes that end at or before C's start and e_K the K-th smallest end among the genes that
  start at or after C's end (none: the contig's end, unbounded); an array [a, b) is within K genes of C exactly when
  a < e_K and b > b_K.
- "within 10 kb": the gap between the array and the gene, max(0, C.start - b, a - C.end), is at most 10,000 bases.
- The tighter windows `1_gene` (K = 1: the gene next to the array, or one it overlaps) and `1_kb` (a gap of at
  most 1,000 bases) are measured too, because the minimal arrays of the example lie about 130 bases from Cas6f.
- On a circular contig every gene is also taken shifted by -n and +n, so windows run across the origin.
- Repeat identity of an array: of its first two repeats, matching bases over the length when they have one length,
  else 1 - (Levenshtein distance / longer length). For arrays of 3 or more repeats also the mean identity of every
  repeat to the array's column-majority consensus.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

STAGES = ("S1", "P1", "P2")
PROK = ("P1", "P2")
K_GENES = 10
BP_WINDOW = 10_000
NEG, POS = -(10 ** 15), 10 ** 15
REPEAT_MATCH = 0.8  # identity at which a call's repeat counts as one of its genome's arrays of 3 (`mark`)
_RC = str.maketrans("ACGTN", "TGCAN")
# The Cas families of the CRISPR worked example: every role of `examples/crispr_cas_neighbourhoods/sql/
# q01_loci_by_domains.sql` (and `checks_bp.sql`) whose name starts with `cas` (31 Pfam families); TniQ is the anchor of
# its Tn7-like loci. Copied from those files on 2026-11-14.
CAS_FAMILIES = {
    "PF09611": "cas8f", "PF09614": "cas5f", "PF09615": "cas7f", "PF09618": "cas6f",
    "PF24121": "cas8b", "PF24122": "cas8b", "PF05107": "cas7", "PF09704": "cas5", "PF09709": "cas8c",
    "PF09559": "cas6", "PF01881": "cas6", "PF10040": "cas6",
    "PF01867": "cas1", "PF09827": "cas2", "PF09707": "cas2", "PF01930": "cas4",
    "PF21384": "cas3", "PF18019": "cas3", "PF21802": "cas3", "PF18395": "cas3", "PF22590": "cas3",
    "PF30455": "cas3", "PF09481": "cas8e", "PF09485": "cas11e",
    "PF16593": "cas9", "PF16592": "cas9", "PF16595": "cas9", "PF22702": "cas9",
    "PF22335": "cas10", "PF12469": "cas10", "PF21918": "cas12a",
}
TNIQ = "PF06527"


def _cpu_s() -> float:
    """CPU-seconds of this process (DuckDB's threads included) and its waited-for children."""
    s, c = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return s.ru_utime + s.ru_stime + c.ru_utime + c.ru_stime


# --- view --------------------------------------------------------------------------------------------------------


def cmd_view(a) -> int:
    src, dst, stages = Path(a.src), Path(a.dst), a.stages.split(",")
    dst.mkdir(parents=True)
    n = 0
    for t in ("contig_index", "seq_locator", "seq_shards", "genome_input"):
        (dst / "l0seq" / t).mkdir(parents=True)
        for s in stages:
            p = src / "l0seq" / t / f"stage={s}"
            if not p.is_dir():
                sys.exit(f"{p}: no such directory")
            os.symlink(p, dst / "l0seq" / t / f"stage={s}")
            n += 1
    os.symlink(src / "l0seq" / "store_params.json", dst / "l0seq" / "store_params.json")
    print(f"{n} stage directories and store_params.json linked into {dst} (stages {stages})")
    return 0


# --- produce -------------------------------------------------------------------------------------------------------


def cmd_produce(a) -> int:
    from biobox.l1 import crispr

    staging = os.environ.get("BIOBOX_STAGING")
    if not staging:
        sys.exit("produce runs under bio-box run-step, which sets BIOBOX_STAGING")
    seq = a.seq_stages.split(",") if a.seq_stages else None
    res = crispr.run(Path(a.root), a.stage, a.cpus, tmp=Path(a.tmp), out=Path(staging), seq_stages=seq,
                     min_repeats=a.min_repeats)
    print(json.dumps(res, sort_keys=True, default=str))
    return 0


# --- the shuffle ----------------------------------------------------------------------------------------------------


def _reaches(v: int, exit_to: dict[int, int], last: int) -> bool:
    seen = set()
    while v != last:
        if v in seen or v not in exit_to:
            return False
        seen.add(v)
        v = exit_to[v]
    return True


def dishuffle(seq: bytes, rng: random.Random) -> bytes:
    """The Altschul-Erickson doublet shuffle (Altschul and Erickson 1985; Kandel et al. 1996): a sequence with exactly
    the dinucleotides of `seq` (as a multiset) and its first and last letter, drawn uniformly among all such sequences.
    The dinucleotides are the edges of a multigraph on the letters; each letter but the last picks one of its out-edges
    uniformly as its last exit, kept only when the last exits lead every letter to the last letter (a tree; rejection
    sampling gives a tree with probability proportional to the product of its edges' multiplicities, which is what
    uniformity over sequences needs); the other out-edges of each letter are shuffled, its last exit appended, and the
    walk from the first letter, taking each letter's edges in order, is the result (an Eulerian trail). Any byte is a
    letter (N and IUPAC codes too)."""
    n = len(seq)
    if n < 3:
        return seq
    a = np.frombuffer(seq, dtype=np.uint8)
    src, dst = a[:-1], a[1:]
    last = int(a[-1])
    edges = {int(v): dst[src == v].tolist() for v in np.unique(src)}
    others = sorted(v for v in edges if v != last)
    while True:
        pick = {v: rng.randrange(len(edges[v])) for v in others}
        exit_to = {v: edges[v][pick[v]] for v in others}
        if all(_reaches(v, exit_to, last) for v in others):
            break
    walk: list = [None] * 256
    for v in sorted(edges):
        e = edges[v]
        if v != last:
            k = pick[v]
            tail = e[k]
            e[k] = e[-1]
            e.pop()
            rng.shuffle(e)
            e.append(tail)
        else:
            rng.shuffle(e)
        walk[v] = e
    pos = [0] * 256
    v = int(a[0])
    out = bytearray(n)
    out[0] = v
    for i in range(1, n):
        e = walk[v]
        j = pos[v]
        pos[v] = j + 1
        v = e[j]
        out[i] = v
    return bytes(out)


def _dinucleotides(a: np.ndarray) -> np.ndarray:
    return np.sort(a[:-1].astype(np.int32) * 256 + a[1:])


def shuffle_windows(seq: bytes, rng: random.Random, window: int) -> bytes:
    """`seq` shuffled by `dishuffle` in non-overlapping windows of `window` bases (the last one shorter). Each window
    is checked as it is made: same length, same first and last base, same multiset of dinucleotides
    (AssertionError otherwise)."""
    out = []
    for i in range(0, len(seq), window):
        w = seq[i:i + window]
        s = dishuffle(w, rng)
        a, b = np.frombuffer(w, np.uint8), np.frombuffer(s, np.uint8)
        if len(s) != len(w) or s[:1] != w[:1] or s[-1:] != w[-1:] or not np.array_equal(_dinucleotides(a),
                                                                                          _dinucleotides(b)):
            raise AssertionError(f"the shuffle of window {i // window} does not keep its dinucleotides")
        out.append(s)
    return b"".join(out)


def contig_rng(seed: int, contig_uid: bytes) -> random.Random:
    return random.Random(int.from_bytes(hashlib.blake2b(seed.to_bytes(8, "big") + contig_uid, digest_size=8).digest(),
                                        "big"))


# --- identity -----------------------------------------------------------------------------------------------------


def _levenshtein(x: str, y: str) -> int:
    prev = list(range(len(y) + 1))
    for i, cx in enumerate(x, 1):
        cur = [i]
        for j, cy in enumerate(y, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (cx != cy)))
        prev = cur
    return prev[-1]


def pair_identity(x: str, y: str) -> float:
    if len(x) == len(y):
        return sum(p == q for p, q in zip(x, y)) / len(x)
    return 1 - _levenshtein(x, y) / max(len(x), len(y))


def consensus_identity(reps: list[str]) -> float | None:
    """Mean identity of the repeats to their column-majority consensus (repeats of one length only; None otherwise
    or for fewer than 3)."""
    if len(reps) < 3 or len({len(r) for r in reps}) != 1:
        return None
    cons = [Counter(col).most_common(1)[0][0] for col in zip(*reps)]
    return sum(sum(p == q for p, q in zip(r, cons)) / len(r) for r in reps) / len(reps)


# --- null 1: shuffled contigs ---------------------------------------------------------------------------------------

_W: dict = {}


def _worker_init(root: str, stages: list[str], scratch: str, window: int, seed: int) -> None:
    from biobox.l0seq.store import SeqStore
    from biobox.l1 import crispr

    _W.update(store=SeqStore(Path(root), stages), tool=crispr._tool(), scratch=Path(scratch), window=window,
              seed=seed, crispr=crispr)


def _null_task(task: tuple[int, str, list[bytes]]) -> dict:
    """Shuffle one batch of contigs into a FASTA file, run MinCED at minimum 2 and 3 on it, and return its arrays."""
    k, stage, contigs = task
    cr, (java, jar), scratch = _W["crispr"], _W["tool"], _W["scratch"]
    d = Path(tempfile.mkdtemp(prefix=f"null-{k:05d}-", dir=scratch))
    t0 = time.process_time()
    try:
        fna = d / "batch.fna"
        where: dict[str, tuple[int, int]] = {}
        bases = 0
        with open(fna, "wb") as fh:
            at = 0
            for cu in contigs:
                seq = _W["store"].get(cu).upper()
                shuf = shuffle_windows(seq, contig_rng(_W["seed"], cu), _W["window"])
                head = b">" + cu.hex().encode() + b"\n"
                fh.write(head + shuf + b"\n")
                where[cu.hex()] = (at + len(head), len(shuf))
                at += len(head) + len(shuf) + 1
                bases += len(shuf)
        shuffle_cpu = time.process_time() - t0
        rows = []
        with open(fna, "rb") as fh:
            for mr in (2, 3):
                gff, txt = d / f"nr{mr}.gff", d / f"nr{mr}.txt"
                done = subprocess.run(cr._argv(java, jar, fna, txt, gff, d, mr), stdin=subprocess.DEVNULL,
                                      stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, cwd=d)
                if done.returncode != 0:
                    raise RuntimeError(f"minced exited {done.returncode}: {done.stderr.decode(errors='replace')[-800:]}")
                with open(gff) as g:
                    for arr in cr.parse_gff(g):
                        off, _ = where[arr.record]
                        reps = []
                        for s, e in arr.units:
                            fh.seek(off + s)
                            reps.append(fh.read(e - s).decode("ascii"))
                        rows.append({"stage": stage, "min_repeats": mr, "contig_uid": arr.record, "start": arr.start,
                                     "end": arr.end, "n_repeats": len(arr.units),
                                     "repeat_len": arr.units[0][1] - arr.units[0][0],
                                     "spacer_len": arr.units[1][0] - arr.units[0][1],
                                     "identity": round(pair_identity(reps[0], reps[1]), 4),
                                     "repeat": reps[0]})
        return {"stage": stage, "contigs": len(contigs), "bases": bases, "rows": rows, "shuffle_cpu_s": shuffle_cpu}
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _stage_contigs(root: Path, stage: str) -> list[tuple[bytes, int]]:
    """(contig_uid, length) of the stage's contigs, once each, sorted by contig_uid: genome_input's contigs with
    contig_index's lengths."""
    import pyarrow.parquet as pq

    from biobox import io
    from biobox.registry.schemas import table

    def rows(t: str, cols: list[str]):
        for f in io.part_files(root / table(t).location / f"stage={stage}"):
            tb = pq.read_table(f, columns=cols, partitioning=None)
            yield from zip(*(tb.column(c).to_pylist() for c in cols))

    want = {cu for (cu,) in rows("genome_input", ["contig_uid"])}
    lens = {cu: n for cu, n in rows("contig_index", ["contig_uid", "contig_len"]) if cu in want}
    missing = want - set(lens)
    if missing:
        raise RuntimeError(f"stage {stage}: {len(missing)} contigs of genome_input have no contig_index row")
    return sorted(lens.items())


def cmd_null(a) -> int:
    root, out = Path(a.root), Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    Path(a.tmp).mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix="e13-null-", dir=a.tmp))
    stages = a.stages.split(",")
    t0 = time.time()
    tasks = []
    for s in stages:
        cur, size = [], 0
        for cu, n in _stage_contigs(root, s):
            cur.append(cu)
            size += n
            if size >= a.batch_bases:
                tasks.append((len(tasks), s, cur))
                cur, size = [], 0
        if cur:
            tasks.append((len(tasks), s, cur))
    per = {s: {"contigs": 0, "bases": 0, "arrays": {2: 0, 3: 0}} for s in stages}
    rows, shuffle_cpu = [], 0.0
    try:
        with ProcessPoolExecutor(a.cpus, initializer=_worker_init,
                                 initargs=(str(root), stages, str(scratch), a.window, a.seed)) as ex:
            for r in ex.map(_null_task, tasks):
                p = per[r["stage"]]
                p["contigs"] += r["contigs"]
                p["bases"] += r["bases"]
                shuffle_cpu += r["shuffle_cpu_s"]
                for row in r["rows"]:
                    p["arrays"][row["min_repeats"]] += 1
                rows += r["rows"]
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    rows.sort(key=lambda r: (r["stage"], r["min_repeats"], r["contig_uid"], r["start"]))
    cols = ["stage", "min_repeats", "contig_uid", "start", "end", "n_repeats", "repeat_len", "spacer_len", "identity",
            "repeat"]
    with open(out / "null_arrays.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    # the pool's workers may be children of a fork server, not of this process, so their CPU is not in this
    # process's rusage: the run's CPU is the unit's (systemd's "CPU time consumed"); the shuffle's is summed here
    summary = {"window": a.window, "seed": a.seed, "stages": {}, "shuffle_cpu_s": round(shuffle_cpu, 1),
               "wall_s": round(time.time() - t0, 1), "tasks": len(tasks), "processes": a.cpus}
    for s, p in per.items():
        mbp = p["bases"] / 1e6
        summary["stages"][s] = {"contigs": p["contigs"], "bases": p["bases"],
                                "arrays_nr2": p["arrays"][2], "arrays_nr3": p["arrays"][3],
                                "per_mbp_nr2": p["arrays"][2] / mbp if mbp else None,
                                "per_mbp_nr3": p["arrays"][3] / mbp if mbp else None}
    (out / "null_summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))
    return 0


# --- windows and the second null ---------------------------------------------------------------------------------


class Windows:
    """The windows of a contig's anchor genes (Cas genes, or one TniQ): per anchor copy, (b_K, e_K) for "within K
    genes" and (start, end) for "within B bases" (module docstring)."""

    def __init__(self, starts: np.ndarray, ends: np.ndarray, n: int, circular: bool, anchors: np.ndarray,
                 k: int = K_GENES):
        shifts = (-n, 0, n) if circular else (0,)
        s_ext = np.concatenate([starts + sh for sh in shifts]) if len(starts) else starts
        e_ext = np.concatenate([ends + sh for sh in shifts]) if len(ends) else ends
        b, e, cs, ce = [], [], [], []
        for i in anchors:
            for sh in shifts:
                c0, c1 = int(starts[i]) + sh, int(ends[i]) + sh
                left = s_ext[e_ext <= c0]
                right = e_ext[s_ext >= c1]
                b.append(int(np.partition(left, -k)[-k]) if len(left) >= k else NEG)
                e.append(int(np.partition(right, k - 1)[k - 1]) if len(right) >= k else POS)
                cs.append(c0)
                ce.append(c1)
        self.b, self.e = np.array(b, dtype=np.int64), np.array(e, dtype=np.int64)
        self.cs, self.ce = np.array(cs, dtype=np.int64), np.array(ce, dtype=np.int64)

    def near_genes(self, a: int, b: int) -> bool:
        return bool(len(self.b)) and bool(np.any((a < self.e) & (b > self.b)))

    def near_bp(self, a: int, b: int, bp: int = BP_WINDOW) -> bool:
        return bool(len(self.cs)) and bool(np.any((self.cs - b <= bp) & (a - self.ce <= bp)))

    def gene_intervals(self, length: int, n: int) -> list:
        """The array starts s in [0, n - length] whose array [s, s + length) is within K genes of an anchor, as
        closed integer intervals (b_K - length + 1, e_K - 1), clipped."""
        hi = n - length
        g = [(max(0, int(x) - length + 1), min(hi, int(y) - 1)) for x, y in zip(self.b, self.e)]
        return [iv for iv in g if iv[0] <= iv[1]]

    def bp_intervals(self, length: int, n: int, bp: int) -> list:
        """The same for "within `bp` bases": (C.start - bp - length, C.end + bp), clipped."""
        hi = n - length
        p = [(max(0, int(x) - bp - length), min(hi, int(y) + bp)) for x, y in zip(self.cs, self.ce)]
        return [iv for iv in p if iv[0] <= iv[1]]


def _union_len(ivs: list[tuple[int, int]]) -> int:
    tot, cur = 0, None
    for lo, hi in sorted(ivs):
        if cur is None or lo > cur[1] + 1:
            if cur is not None:
                tot += cur[1] - cur[0] + 1
            cur = [lo, hi]
        else:
            cur[1] = max(cur[1], hi)
    return tot + (cur[1] - cur[0] + 1 if cur else 0)


NEAR = ("10_genes", "10_kb", "1_gene", "1_kb")  # the windows beside a Cas gene that are measured


def placement_share(length: int, contigs: list[tuple[int, Windows | None, Windows | None]]) -> dict[str, float]:
    """The second null: the probability that an array of `length` bases, its start drawn uniformly from every start
    at which it fits on a contig of its genome, lies in each window of `NEAR` beside a Cas gene:
    f = sum_c |N_c(length)| / sum_c (n_c - length + 1), over the genome's contigs c with n_c >= length, where N_c is
    the union of the start intervals of the window (`Windows.gene_intervals` of the contig's 10-gene and 1-gene
    windows, `bp_intervals` at 10,000 and 1,000 bases; a contig without a Cas gene: empty). The expected number of
    calls in a window, were they placed at random in their own genomes, is the sum of f over them."""
    num = dict.fromkeys(NEAR, 0)
    den = 0
    for n, w10, w1 in contigs:
        if n < length:
            continue
        den += n - length + 1
        if w10 is not None:
            num["10_genes"] += _union_len(w10.gene_intervals(length, n))
            num["1_gene"] += _union_len(w1.gene_intervals(length, n))
            num["10_kb"] += _union_len(w10.bp_intervals(length, n, BP_WINDOW))
            num["1_kb"] += _union_len(w10.bp_intervals(length, n, 1000))
    return {k: (v / den if den else 0.0) for k, v in num.items()}


# --- analyse ---------------------------------------------------------------------------------------------------------


def _q(con, sql: str, params=None):
    cur = con.execute(sql, params or [])
    return cur.to_arrow_table() if hasattr(cur, "to_arrow_table") else cur.fetch_arrow_table()  # DuckDB 2 and 1


def _files_sql(files: list[Path]) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def _arrays(con, view: Path, prod: str, stage: str) -> list[dict]:
    """The producer's arrays of the stage in the view: coordinates, repeats (in order), first spacer length."""
    from biobox import io

    ff = io.part_files(view / f"l1/feature/producer={prod}/stage={stage}")
    ef = io.part_files(view / f"l1/crispr/producer={prod}/stage={stage}")
    if not ff:
        return []
    t = _q(con, f"""
        SELECT f.feature_uid, f.contig_uid, f.start, f."end", e.contig_len, e.n_rep, e.sp_len, e.reps
        FROM read_parquet({_files_sql(ff)}) f JOIN (
          SELECT feature_uid, any_value(contig_len) AS contig_len,
                 count(*) FILTER (WHERE element = 'repeat') AS n_rep,
                 min(CASE WHEN ordinal = 1 THEN "end" - start END) AS sp_len,
                 list(seq ORDER BY ordinal) FILTER (WHERE element = 'repeat') AS reps
          FROM read_parquet({_files_sql(ef)}) GROUP BY feature_uid) e USING (feature_uid)
        ORDER BY f.contig_uid, f.start, f."end" """)
    out = []
    for r in t.to_pylist():
        reps = r["reps"]
        out.append({"stage": stage, "feature_uid": r["feature_uid"], "contig_uid": r["contig_uid"], "start": r["start"],
                    "end": r["end"], "contig_len": r["contig_len"], "n_repeats": r["n_rep"],
                    "repeat_len": len(reps[0]), "spacer_len": r["sp_len"],
                    "identity": round(pair_identity(reps[0], reps[1]), 4),
                    "consensus_identity": consensus_identity(reps), "repeat": reps[0]})
    return out


def _quantiles(xs: list[float]) -> dict:
    if not xs:
        return {"n": 0}
    q = np.percentile(np.array(xs, dtype=float), [0, 10, 25, 50, 75, 90, 100], method="lower")
    return {"n": len(xs), **{k: round(float(v), 4) for k, v in zip(("min", "p10", "p25", "median", "p75", "p90",
                                                                     "max"), q)}}


def _identity_bins(xs: list[float]) -> dict:
    return {"1.0": sum(x == 1.0 for x in xs), "[0.9,1)": sum(0.9 <= x < 1.0 for x in xs),
            "[0.8,0.9)": sum(0.8 <= x < 0.9 for x in xs), "<0.8": sum(x < 0.8 for x in xs)}


def _overlaps(arr: dict, others: list[dict]) -> bool:
    return any(o["start"] < arr["end"] and arr["start"] < o["end"] for o in others)


def cmd_analyse(a) -> int:
    from biobox import catalog
    from biobox.l1 import crispr

    root, view, out = Path(a.root), Path(a.view), Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    Path(a.tmp).mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix="e13-analyse-", dir=a.tmp))
    t0, c0 = time.time(), _cpu_s()
    try:
        (scratch / "duckdb").mkdir()
        con = catalog.open_catalog(catalog.catalog_path(root, a.release), local_dir=scratch / "catalog",
                                   threads=a.threads, memory_limit="8GB", temp_directory=scratch / "duckdb")
        res = _analyse(con, view, Path(a.null), a, crispr, out)
        con.close()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    res["cost"] = {"cpu_s": round(_cpu_s() - c0, 1), "wall_s": round(time.time() - t0, 1)}
    (out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True, default=str) + "\n")
    print(json.dumps(res["cost"]))
    return 0


def _write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    def cell(v):
        if isinstance(v, (bytes, bytearray)):
            return bytes(v).hex()
        if isinstance(v, bool):
            return "1" if v else "0"
        return "" if v is None else str(v)

    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(cell(r.get(c)) for c in cols) + "\n")


def _analyse(con, view: Path, null_dir: Path, a, crispr, out: Path) -> dict:
    import pyarrow as pa

    v = crispr.minced_version()
    p2, p3 = f"minced-{v}-{crispr.param_tag(2)}", f"minced-{v}-{crispr.param_tag(3)}"
    res: dict = {"producers": {"nr2": p2, "nr3": p3}, "release": a.release, "cas_families": CAS_FAMILIES,
                 "k_genes": K_GENES, "bp_window": BP_WINDOW}
    libs = [r["library"] for r in _q(con, "SELECT DISTINCT library FROM l1_domains ORDER BY 1").to_pylist()]
    lib = a.library or (libs[0] if len(libs) == 1 else None)
    if lib is None:
        raise SystemExit(f"l1_domains holds libraries {libs}; name one with --library")
    res["library"] = lib

    # contigs of the stages: genome, length, topology, stage
    ct = _q(con, f"""
        SELECT d.contig_uid, d.genome_uid, d.contig_len, d.is_circular, min(gi.stage) AS stage,
               count(DISTINCT gi.stage) AS n_stages
        FROM contig_dim d JOIN genome_input gi USING (contig_uid)
        WHERE gi.stage IN {tuple(STAGES)} GROUP BY ALL ORDER BY d.contig_uid""").to_pylist()
    res["contigs_in_two_stages"] = sum(r["n_stages"] > 1 for r in ct)
    cidx = {r["contig_uid"]: i for i, r in enumerate(ct)}
    bases = Counter()
    for r in ct:
        bases[r["stage"]] += r["contig_len"]
    res["bases"] = dict(bases)
    gbp_prok = sum(bases[s] for s in PROK) / 1e9
    res["gbp_p1_p2"] = gbp_prok

    # genes with their Cas / TniQ roles, in contig order (genes sorted by start, end: the tape's order)
    con.register("ctab", pa.table({"contig_uid": pa.array([r["contig_uid"] for r in ct], pa.binary(16)),
                                   "ci": pa.array(range(len(ct)), pa.int64())}))
    roles = dict(CAS_FAMILIES, **{TNIQ: "tniQ"})
    values = ", ".join(f"('{pf}', '{rl}')" for pf, rl in sorted(roles.items()))
    g = _q(con, f"""
        WITH roles(pf, role) AS (VALUES {values}),
        dom AS (SELECT DISTINCT d.protein_uid, r.role FROM l1_domains d
                JOIN roles r ON r.pf = split_part(d.model_acc, '.', 1) WHERE d.library = ?),
        pr AS (SELECT protein_uid, bool_or(role <> 'tniQ') AS is_cas, bool_or(role = 'tniQ') AS is_tniq,
                      bool_or(role = 'cas7f') AS is_cas7f, bool_or(role = 'cas6f') AS is_cas6f
               FROM dom GROUP BY protein_uid),
        g AS (SELECT DISTINCT stage, contig_uid, start, "end", protein_uid FROM protein_occurrence
              WHERE stage IN {tuple(STAGES)})
        SELECT c.ci, g.start, g."end", coalesce(pr.is_cas, false) AS is_cas, coalesce(pr.is_tniq, false) AS is_tniq,
               coalesce(pr.is_cas7f, false) AS is_cas7f, coalesce(pr.is_cas6f, false) AS is_cas6f
        FROM g JOIN ctab c USING (contig_uid) LEFT JOIN pr USING (protein_uid)
        ORDER BY c.ci, g.start, g."end" """, [lib])
    ci = g.column("ci").to_numpy()
    gs, ge = g.column("start").to_numpy().astype(np.int64), g.column("end").to_numpy().astype(np.int64)
    flags = {k: np.asarray(g.column(k).to_numpy(), dtype=bool) for k in ("is_cas", "is_tniq", "is_cas7f", "is_cas6f")}
    bounds = np.searchsorted(ci, np.arange(len(ct) + 1))
    res["genes"] = int(len(ci))
    res["genes_with_cas_domain"] = int(flags["is_cas"].sum())

    def genes_of(i: int):
        lo, hi = bounds[i], bounds[i + 1]
        return gs[lo:hi], ge[lo:hi], {k: np.flatnonzero(f[lo:hi]) for k, f in flags.items()}

    cas_win: dict[int, Windows] = {}
    cas_win1: dict[int, Windows] = {}  # the 1-gene windows (an adjacent gene, or one the array overlaps)
    genome_cas: Counter = Counter()
    for i in np.unique(ci[flags["is_cas"]]):
        s, e, fl = genes_of(int(i))
        r = ct[int(i)]
        cas_win[int(i)] = Windows(s, e, r["contig_len"], r["is_circular"], fl["is_cas"])
        cas_win1[int(i)] = Windows(s, e, r["contig_len"], r["is_circular"], fl["is_cas"], k=1)
        genome_cas[r["genome_uid"]] += len(fl["is_cas"])
    genome_contigs: dict[bytes, list[int]] = {}
    for i, r in enumerate(ct):
        genome_contigs.setdefault(r["genome_uid"], []).append(i)
    res["genomes"] = {s: len({r["genome_uid"] for r in ct if r["stage"] == s}) for s in STAGES}
    res["genomes_with_cas_domain"] = {s: len({r["genome_uid"] for r in ct if r["stage"] == s
                                               and genome_cas[r["genome_uid"]]}) for s in STAGES}

    def annotate(arr: dict) -> dict:
        i = cidx[arr["contig_uid"]]
        r = ct[i]
        w = cas_win.get(i)
        s, e, _ = genes_of(i)
        n = r["contig_len"]
        a_, b_ = arr["start"], arr["end"]
        if r["is_circular"]:
            s3, e3 = np.concatenate([s - n, s, s + n]), np.concatenate([e - n, e, e + n])
        else:
            s3, e3 = s, e
        f = placement_share(b_ - a_, [(ct[j]["contig_len"], cas_win.get(j), cas_win1.get(j))
                                      for j in genome_contigs[r["genome_uid"]]])
        near = {"near_10_genes": w is not None and w.near_genes(a_, b_),
                "near_10_kb": w is not None and w.near_bp(a_, b_, BP_WINDOW),
                "near_1_gene": w is not None and cas_win1[i].near_genes(a_, b_),
                "near_1_kb": w is not None and w.near_bp(a_, b_, 1000)}
        return {**arr, "genome_uid": r["genome_uid"], **near,
                "inside_gene": bool(np.any((s3 <= a_) & (e3 >= b_))),
                "overlaps_gene": bool(np.any((s3 < b_) & (e3 > a_))),
                "genome_has_cas": genome_cas[r["genome_uid"]] > 0, **{f"f_{k}": round(v, 6) for k, v in f.items()}}

    # the producers' arrays
    arrays = {(k, s): _arrays(con, view, p, s) for k, p in (("nr2", p2), ("nr3", p3)) for s in STAGES}
    per_stage, extras, nr3_all = {}, [], []
    for s in STAGES:
        a2, a3 = arrays[("nr2", s)], arrays[("nr3", s)]
        by3: dict[bytes, list[dict]] = {}
        for x in a3:
            by3.setdefault(x["contig_uid"], []).append(x)
        fu2 = {x["feature_uid"]: x for x in a2}
        fu3 = {x["feature_uid"] for x in a3}
        same = sum(1 for x in a3 if x["feature_uid"] in fu2 and fu2[x["feature_uid"]]["repeat"] == x["repeat"]
                   and fu2[x["feature_uid"]]["n_repeats"] == x["n_repeats"])
        ex = [x for x in a2 if not _overlaps(x, by3.get(x["contig_uid"], []))]
        by2: dict[bytes, list[dict]] = {}
        for x in a2:
            by2.setdefault(x["contig_uid"], []).append(x)
        per_stage[s] = {"nr3": len(a3), "nr2": len(a2), "nr3_found_identically_at_nr2": same,
                        "nr3_missing_at_nr2": sum(1 for x in a3 if x["feature_uid"] not in fu2),
                        "nr3_overlapped_by_no_nr2": sum(1 for x in a3 if not _overlaps(x, by2.get(x["contig_uid"], []))),
                        "nr2_of_2_repeats": sum(x["n_repeats"] == 2 for x in a2),
                        "nr2_overlapping_nr3_not_equal": sum(1 for x in a2 if x["feature_uid"] not in fu3
                                                             and _overlaps(x, by3.get(x["contig_uid"], []))),
                        "extra": len(ex), "extra_repeats_hist": dict(sorted(Counter(x["n_repeats"] for x in ex).items())),
                        "genomes_with_extra": len({ct[cidx[x["contig_uid"]]]["genome_uid"] for x in ex}),
                        "bases": bases[s]}
        if s in PROK:
            extras += [annotate(x) for x in ex]
            nr3_all += [annotate(x) for x in a3]
    res["per_stage"] = per_stage

    # an independent marker of a real CRISPR repeat: the call's repeat matches (identity >= REPEAT_MATCH, either
    # strand) the first repeat of an array found at 3 in the same genome. It misses true arrays whose repeat family
    # has no array of 3 in the genome (the minimal I-F3 arrays of the example differ from their host's I-F repeat)
    known: dict[bytes, set[str]] = {}
    for x in nr3_all:
        known.setdefault(x["genome_uid"], set()).add(x["repeat"])

    def mark(x: dict) -> dict:
        reps = known.get(x["genome_uid"], ())
        rc = x["repeat"].translate(_RC)[::-1]
        best = max((max(pair_identity(x["repeat"], r), pair_identity(rc, r)) for r in reps), default=0.0)
        x["best_nr3_repeat_identity"] = round(best, 4)
        x["matches_nr3_repeat"] = best >= REPEAT_MATCH
        return x

    extras = [mark(x) for x in extras]

    def tally(xs: list[dict]) -> dict:
        n = len(xs)
        d = {"n": n, "of_2_repeats": sum(x["n_repeats"] == 2 for x in xs)}
        for k in ("inside_gene", "overlaps_gene"):
            d[k] = sum(x[k] for x in xs)
        for k in NEAR:
            d[f"near_{k}"] = sum(x[f"near_{k}"] for x in xs)
            d[f"near_{k}_not_inside_gene"] = sum(x[f"near_{k}"] and not x["inside_gene"] for x in xs)
            d[f"expected_near_{k}_random_placement"] = round(sum(x[f"f_{k}"] for x in xs), 2)
            d[f"expected_near_{k}_not_inside_gene_random_placement"] = round(
                sum(x[f"f_{k}"] for x in xs if not x["inside_gene"]), 2)
        d["in_genome_without_cas"] = sum(not x["genome_has_cas"] for x in xs)
        if xs and "matches_nr3_repeat" in xs[0]:
            d["matches_nr3_repeat"] = sum(x["matches_nr3_repeat"] for x in xs)
        d["near_10_genes_inside_gene"] = sum(x["near_10_genes"] and x["inside_gene"] for x in xs)
        d["far_not_inside_gene_genome_with_cas"] = sum(not x["near_10_kb"] and not x["near_10_genes"]
                                                        and not x["inside_gene"] and x["genome_has_cas"] for x in xs)
        return d

    res["extra_p1_p2"] = tally(extras)
    res["extra_by_stage"] = {s: tally([x for x in extras if x["stage"] == s]) for s in PROK}
    res["extra_of_2_repeats_p1_p2"] = tally([x for x in extras if x["n_repeats"] == 2])
    res["extra_of_3_or_more_p1_p2"] = tally([x for x in extras if x["n_repeats"] > 2])
    res["nr3_p1_p2"] = tally(nr3_all)
    res["extra_cross"] = dict(sorted(Counter(
        f"near10g={int(x['near_10_genes'])} inside_gene={int(x['inside_gene'])} genome_cas={int(x['genome_has_cas'])}"
        for x in extras).items()))
    res["distributions"] = {
        name: {"repeat_len": _quantiles([x["repeat_len"] for x in xs]),
               "spacer_len": _quantiles([x["spacer_len"] for x in xs]),
               "identity_first_two": _quantiles([x["identity"] for x in xs]),
               "identity_bins": _identity_bins([x["identity"] for x in xs]),
               "consensus_identity": _quantiles([x["consensus_identity"] for x in xs
                                                 if x["consensus_identity"] is not None])}
        for name, xs in (("extra_p1_p2", extras), ("extra_near_10_genes", [x for x in extras if x["near_10_genes"]]),
                         ("extra_far", [x for x in extras if not x["near_10_genes"] and not x["near_10_kb"]]),
                         ("nr3_p1_p2", nr3_all))}

    # the shuffled null, placed beside the real genes (a shuffled contig keeps its coordinates)
    ns = json.loads((null_dir / "null_summary.json").read_text())
    nrows = []
    with open(null_dir / "null_arrays.tsv") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            for k in ("min_repeats", "start", "end", "n_repeats", "repeat_len", "spacer_len"):
                r[k] = int(r[k])
            r["identity"] = float(r["identity"])
            r["contig_uid"] = bytes.fromhex(r["contig_uid"])
            r["consensus_identity"] = None
            nrows.append(r)
    null_extra = []
    for s in PROK:
        n3: dict[bytes, list[dict]] = {}
        for r in nrows:
            if r["stage"] == s and r["min_repeats"] == 3:
                n3.setdefault(r["contig_uid"], []).append(r)
        null_extra += [annotate(r) for r in nrows if r["stage"] == s and r["min_repeats"] == 2
                       and not _overlaps(r, n3.get(r["contig_uid"], []))]
    null_extra = [mark(x) for x in null_extra]
    null_bases_prok = sum(ns["stages"][s]["bases"] for s in PROK if s in ns["stages"])
    if not null_bases_prok:
        raise SystemExit(f"{null_dir}: the null holds neither P1 nor P2")
    res["null_shuffle"] = {"summary": ns, "extra_p1_p2": tally(null_extra),
                           "null_bases_p1_p2": null_bases_prok,
                           "identity_bins_extra": _identity_bins([x["identity"] for x in null_extra]),
                           "distributions_extra": {"repeat_len": _quantiles([x["repeat_len"] for x in null_extra]),
                                                   "spacer_len": _quantiles([x["spacer_len"] for x in null_extra]),
                                                   "identity_first_two": _quantiles([x["identity"]
                                                                                     for x in null_extra])}}
    null_gbp = null_bases_prok / 1e9

    # rules: what each adds and what it wrongly adds, per Gbp of P1 + P2. A rule is a filter `pre` on the extra
    # calls and, for the A rules, a window beside a Cas gene. For a windowed rule, the false calls it keeps are at
    # most the second null's sum of f over every call that passes `pre` (were every such call false and placed at
    # random in its genome, that many would land in the window by chance). For a rule without a window the estimate
    # of what it wrongly adds is the kept calls far from any Cas gene (neither within 10 genes nor 10 kb), which
    # assumes true two-repeat arrays sit beside Cas genes. Both carry the shuffled null's chance calls the rule keeps.
    def rule(name: str, desc: str, pre, window: str | None) -> dict:
        passed = [x for x in extras if pre(x)]
        keep = (lambda x: x[f"near_{window}"]) if window else (lambda x: True)
        kept = [x for x in passed if keep(x)]
        kn = [x for x in null_extra if pre(x) and keep(x)]
        far = sum(not x["near_10_genes"] and not x["near_10_kb"] for x in kept)
        wrong = round(sum(x[f"f_{window}"] for x in passed), 1) if window else far
        return {"rule": name, "desc": desc, "adds": len(kept), "adds_per_gbp": round(len(kept) / gbp_prok, 2),
                "wrong_estimate": wrong, "wrong_per_gbp": round(wrong / gbp_prok, 2),
                "wrong_basis": "random placement, upper bound" if window else "far from any Cas gene",
                "chance_shuffle": len(kn), "chance_shuffle_per_gbp": round(len(kn) / null_gbp, 3),
                "inside_gene": sum(x["inside_gene"] for x in kept), "far_from_cas": far,
                "matches_nr3_repeat": sum(x["matches_nr3_repeat"] for x in kept),
                "null_shuffle_matches_nr3_repeat": sum(x["matches_nr3_repeat"] for x in kn),
                "in_genome_without_cas": sum(not x["genome_has_cas"] for x in kept)}

    def everything(x):
        return True

    def not_inside(x):
        return not x["inside_gene"]

    def identical(x):
        return x["identity"] == 1.0

    rules = [
        rule("A", "nr3 everywhere, plus nr2 calls within 10 genes of a Cas gene", everything, "10_genes"),
        rule("A1", "as A, not wholly inside a called gene", not_inside, "10_genes"),
        rule("A2", "nr3 everywhere, plus nr2 calls within 10 kb of a Cas gene", everything, "10_kb"),
        rule("A3", "as A, the two repeats identical", identical, "10_genes"),
        rule("A4", "nr3 everywhere, plus nr2 calls next to a Cas gene (1 gene), not wholly inside a called gene",
             not_inside, "1_gene"),
        rule("A5", "nr3 everywhere, plus nr2 calls within 1 kb of a Cas gene, not wholly inside a called gene",
             not_inside, "1_kb"),
        rule("B", "nr2 everywhere", everything, None),
        rule("B1", "nr2 everywhere, the two repeats identical", identical, None),
        rule("B2", "nr2 everywhere, not wholly inside a called gene", not_inside, None),
    ]
    res["rules"] = rules
    _write_tsv(out / "rules.tsv", rules, list(rules[0]))

    # the worked example's Tn7-like loci: a TniQ gene with a Cas7f gene within 10 genes (by rank, on its contig)
    by_contig: dict[tuple[str, bytes], list[dict]] = {}
    for (k, s), xs in arrays.items():
        for x in xs:
            by_contig.setdefault((k, x["contig_uid"]), []).append(x)
    loci = []
    for i in np.unique(ci[flags["is_tniq"]]):
        i = int(i)
        s, e, fl = genes_of(i)
        r = ct[i]
        m = len(s)
        for t in fl["is_tniq"]:
            d = np.abs(fl["is_cas7f"] - t)
            if r["is_circular"]:
                d = np.minimum(d, m - d)
            if not len(d) or d.min() > K_GENES:
                continue
            w = Windows(s, e, r["contig_len"], r["is_circular"], np.array([t]))
            near = {k: [x for x in by_contig.get((k, r["contig_uid"]), []) if w.near_genes(x["start"], x["end"])]
                    for k in ("nr2", "nr3")}
            c6 = fl["is_cas6f"]
            gained = [x for x in near["nr2"] if not _overlaps(x, near["nr3"])]

            def to_cas6f(x):
                if not len(c6):
                    return None
                return int(min(max(0, int(s[j]) - x["end"], x["start"] - int(e[j])) for j in c6))

            loci.append({"stage": r["stage"], "contig_uid": r["contig_uid"], "tniq_start": int(s[t]),
                         "tniq_end": int(e[t]), "cas7f_rank_distance": int(d.min()), "arrays_nr3": len(near["nr3"]),
                         "arrays_nr2": len(near["nr2"]), "gains": bool(gained and not near["nr3"]),
                         "gained_arrays": ";".join(f"{x['start']}-{x['end']}:{x['n_repeats']}r:{x['repeat']}"
                                                   f":to_cas6f={to_cas6f(x)}" for x in gained)})
    if loci:
        con.register("lc", pa.table({"contig_uid": pa.array(sorted({x["contig_uid"] for x in loci}), pa.binary(16))}))
        acc = {r["contig_uid"]: r["accession"] for r in _q(con, """
            SELECT contig_uid, first(accession ORDER BY l0_object_sha256, record_ordinal) AS accession
            FROM contig_occurrence WHERE contig_uid IN (SELECT contig_uid FROM lc) GROUP BY contig_uid""").to_pylist()}
        for x in loci:
            x["accession"] = acc.get(x["contig_uid"])
    loci.sort(key=lambda x: (x["stage"], x["accession"] or "", x["tniq_start"]))
    res["tn7_loci"] = {"loci": len(loci), "with_array_nr3": sum(x["arrays_nr3"] > 0 for x in loci),
                       "with_array_nr2": sum(x["arrays_nr2"] > 0 for x in loci),
                       "gain_at_nr2": sum(x["gains"] for x in loci)}
    _write_tsv(out / "tn7_loci.tsv", loci, ["stage", "accession", "contig_uid", "tniq_start", "tniq_end",
                                            "cas7f_rank_distance", "arrays_nr3", "arrays_nr2", "gains",
                                            "gained_arrays"])
    cols = ["stage", "contig_uid", "genome_uid", "start", "end", "n_repeats", "repeat_len", "spacer_len", "identity",
            "near_10_genes", "near_10_kb", "near_1_gene", "near_1_kb", "inside_gene", "overlaps_gene", "genome_has_cas",
            "f_10_genes", "f_10_kb", "f_1_gene", "f_1_kb", "best_nr3_repeat_identity", "matches_nr3_repeat", "repeat"]
    _write_tsv(out / "extra_arrays.tsv", extras, cols)
    _write_tsv(out / "null_extra_arrays.tsv", null_extra, cols)
    return res


# --- command line --------------------------------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("view", help="the scratch view of the stages' sequence store")
    v.add_argument("src")
    v.add_argument("dst")
    v.add_argument("stages")
    p = sub.add_parser("produce", help="the producer at a chosen minimum, as a run-step command")
    p.add_argument("--root", required=True)
    p.add_argument("--stage", required=True)
    p.add_argument("--cpus", type=int, required=True)
    p.add_argument("--min-repeats", type=int, required=True)
    p.add_argument("--tmp", required=True)
    p.add_argument("--seq-stages")
    n = sub.add_parser("null", help="MinCED at 2 and 3 on shuffled contigs")
    n.add_argument("--root", required=True)
    n.add_argument("--stages", required=True)
    n.add_argument("--cpus", type=int, required=True)
    n.add_argument("--tmp", required=True)
    n.add_argument("--out", required=True)
    n.add_argument("--window", type=int, default=1000)
    n.add_argument("--seed", type=int, default=13)
    n.add_argument("--batch-bases", type=int, default=25_000_000)
    s = sub.add_parser("analyse", help="the extra calls on the release's genes and domains, and the nulls")
    s.add_argument("--root", required=True)
    s.add_argument("--release", required=True)
    s.add_argument("--view", required=True)
    s.add_argument("--null", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--tmp", required=True)
    s.add_argument("--threads", type=int, default=4)
    s.add_argument("--library")
    a = ap.parse_args(argv)
    return {"view": cmd_view, "produce": cmd_produce, "null": cmd_null, "analyse": cmd_analyse}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())

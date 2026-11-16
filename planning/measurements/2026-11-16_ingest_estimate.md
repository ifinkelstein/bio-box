# The ingest estimate counts what a stage stores (card B-52, DEC-47 (2))

Date: 2026-11-16, wilkcomp02 (the pod's clock, CST). Every number carries a label: *measured* (host, date, command or
record), *derived* (formula and inputs) or *estimate*. Read-only on the data root `/tmp/if2453/r10/root` (bound `:ro`
into every command that read it); nothing written under `/stor`; scratch `/tmp/if2453/b52/`. The build of `PL1`
(unit `c26-pl1-2`) was running throughout; only the published `ingest` record of `PL1` was read, no file of a step
still being written.

## Result

- **What the cost follows**: what a stage stores, more than what it reads. After hashing, 0.84 CPU-minutes per Gbp
  read plus 3.41 per Gbp stored hold all nine records within 6.2% (*derived*, below); but which contigs a stage
  stores is known only after hashing.
- **What is known before launch**: the stage's bytes of FASTA (as before) and the bases of the contigs its
  `genome_input` lists that another stage's `contig_index` already holds (the inputs step publishes `genome_input`
  before ingest; ingest reads those contigs and stores nothing of them, DEC-26). For `PL1` that is 1,177,838,866
  bases (7,046 contigs), for H1 141,683,767 (514), for P3 to P6 7 to 186 kbp (one contig each), for S1, P1, P2 none.
- **The model** (`stage.ingest_cpu_min(n_bytes, held_bases)`): 0.07 CPU-minutes + 4.11 per Gbyte of FASTA not held
  (`INGEST_CPU_MIN = 4.11` per `INGEST_BP = 10**9`; the 0.07 of card C-13 unchanged). **Every built stage's record is
  within 20% of its estimate, eight of nine within 10%**: from P1 7.2% over to S1 14.2% under; `PL1` 3.4% over
  (26.612 against 25.749), H1 1.2% over (*derived*, table). Before the card: `PL1` 20.0% under, S1 21.1% under, H1
  13.5% under.
- `build-stage` computes the estimate when it reaches the ingest step (`held_bases`, one DuckDB query on one thread
  over the stage's `genome_input` and the other stages' `contig_index` files, those ingest's own anti-join reads) and
  says it in one line, within the margin or over it. The code's `held_bases` on the real root gave exactly the
  table's held bases for all nine stages (*measured*, unit `b52-held`, below).
- The batch cut is unchanged in form: a pass-1 shard and a batch are estimated with every byte counted as stored
  (a pass-1 shard stores nothing; a batch stores at most what it reads); `pass1_bytes` grows by 8.8% with the lower
  rate (*derived*: 4.4705 / 4.11). Shown on the fixture (test below).

## 1. Setup

- Branch `task/B-52` from `origin/plan/task-cards` at `d033220`.
- Reads in the runtime image `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, each in a transient unit
  (`systemd-run --user -p MemoryMax=8G -p CPUQuota=100% --nice=19`, one thread), the data root bound read-only:
  - unit `b52-table`, 08:35:58 to 08:36:08 CST (file times of its output): `b52_ingest_table.py` (Appendix A),
    `/tmp/if2453/b52/table.json` (sha256 `46911d9c87571b20d6379fe46719d5ae3fc3a4792371ae14dc7e1e11716b2461`).
  - unit `b52-held`, ended 08:48:03 CST (file time; a first launch at 08:47 failed on `fetchone`, fixed in the code
    before this run): `b52_held_check.py` (Appendix B) with this card's `src`, the code's
    `held_bases` and `ingest_cpu_min` on the real root, `/tmp/if2453/b52/held.jsonl`.
  - The fit is arithmetic on `table.json` (`b52_fit.py`, Appendix C).
- At 08:34 the account was at its weekday cap (the `PL1` build's 9 scan shards of 4 and another agent's 4); these
  reads used one CPU each, one at a time, as the coordinator asked.

## 2. What each built stage read and stored (*measured*)

Per stage, in ingest order: its ingest run record (`runs/stage=<S>/ingest.json`: `cpu_min`, every one a single
process, none in batches: no `l0seq/ingest_plan` in the root); the bytes of FASTA of the record's `--objects` (every
object a plain file, so bytes on disk); the records of its `contig_occurrence` and the bases they hold (repeats
counted); the bases its `genome_input` lists that a stage ingested before it held as a default variant ("held");
the default variants and bases of its own `contig_index` ("stored").

| Stage | Run record (start, UTC) | FASTA bytes read | Records | Bases read | Held before | Stored contigs | Stored bases | CPU-min |
|---|---|---|---|---|---|---|---|---|
| S1 | 2026-11-10 16:57 | 2,121,124,779 | 34,076 | 2,084,350,440 | 0 | 27,258 | 1,664,109,938 | 7.540 |
| P1 | 2026-11-10 18:13 | 1,029,384,352 | 539 | 1,016,625,803 | 0 | 539 | 1,016,625,803 | 4.611 |
| P2 | 2026-11-10 18:43 | 9,716,599,273 | 5,018 | 9,596,180,328 | 0 | 5,018 | 9,596,180,328 | 38.335 |
| P3 | 2026-11-11 19:22 | 10,074,413,081 | 9,108 | 9,949,084,485 | 6,963 | 9,107 | 9,949,077,522 | 43.889 |
| P4 | 2026-11-12 02:26 | 10,101,497,437 | 4,006 | 9,976,464,003 | 7,053 | 4,005 | 9,976,456,950 | 42.979 |
| P5 | 2026-11-12 03:27 | 10,092,827,773 | 4,112 | 9,967,897,315 | 185,683 | 4,111 | 9,967,711,632 | 43.855 |
| P6 | 2026-11-12 04:29 | 10,087,447,605 | 4,228 | 9,962,553,301 | 17,198 | 4,227 | 9,962,536,103 | 41.486 |
| H1 | 2026-11-15 13:44 | 1,962,116,403 | 7,230 | 1,961,987,211 | 141,683,767 | 6,714 | 1,819,831,338 | 7.645 |
| PL1 | 2026-11-16 03:06 | 7,425,668,755 | 72,556 | 7,297,669,596 | 1,177,838,866 | 61,647 | 5,924,552,875 | 26.612 |

- Consistency (*derived*): for every stage, listed contigs = held + stored (H1 7,228 = 514 + 6,714; `PL1` 68,693 =
  7,046 + 61,647), and listed bases = held + stored bases (`PL1` 7,102,391,741 = 1,177,838,866 + 5,924,552,875, the
  card's figures). `PL1`'s held set is the inputs note's 7,046 contigs, 1,177,838,866 bases.
- S1's FASTA holds 34,076 records of 27,258 distinct contigs: 6,818 records, 420,240,502 bases, repeat a sequence
  already read (*derived*: 34,076 − 27,258; 2,084,350,440 − 1,664,109,938). `PL1` repeats 3,863 records (72,556 −
  68,693; the inputs note's repeats), H1 2.
- The run records were written by two images (`c30be585…` for S1 to P6, `d67300aa…` for H1 and `PL1`) and at 16 to
  36 `--threads`; ingest is one process at one thread (wall ≈ CPU in every record: P3 2,730 s wall, 2,633 s CPU).

## 3. The residual of the present estimate, and the fits tried (*derived*, `b52_fit.py`)

The measure is the guard's: (record − estimate) / estimate. Each fit below minimises the squared relative error
over the nine records, with the fixed 0.07 CPU-minutes (card C-13) kept.

| Model | Coefficients (CPU-min per 10⁹) | Worst | Within 10% | Within 20% | Residuals S1, P1, P2, P3, P4, P5, P6, H1, PL1 |
|---|---|---|---|---|---|
| Before B-52: per byte of FASTA read | 4.4705 (42.9 / 9.596) | −21.1% | 5 | 7 | −21.1, −1.3, −11.9, −2.7, −5.0, −3.0, −8.1, −13.5, −20.0 |
| One rate per byte read, refitted | 3.987 | −11.6% | 6 | 9 | −11.6, +10.5, −1.2, +9.1, +6.5, +8.8, +3.0, −3.1, −10.3 |
| **One rate per byte not held (chosen)** | **4.11** (fit 4.1100) | −14.2% | **8** | **9** | −14.2, +7.2, −4.2, +5.8, +3.3, +5.5, −0.1, +1.2, +3.4 |
| Per byte read + per byte not held | 0.881 + 3.205 | −13.7% | 8 | 9 | −13.7, +7.9, −3.6, +6.5, +4.0, +6.2, +0.5, +0.2, −0.1 |
| After hashing: per base stored | 4.302 | −7.3% | 9 | 9 | +4.3, +3.8, −7.3, +2.4, −0.0, +2.1, −3.4, −3.2, +4.1 |
| After hashing: per base read + per base stored | 0.838 + 3.415 | −6.2% | 9 | 9 | +0.6, +5.0, −6.2, +3.6, +1.1, +3.3, −2.2, −3.6, +0.7 |

- **Rejected, the rate per byte read alone** (refitted): within 20% everywhere, but S1 and `PL1` 10 to 12% under and
  P1, P3, P5 9 to 11% over: it cannot tell a stage that stores what it reads from one that does not.
- **Rejected, a rate per byte read plus a rate per byte not held**: no record moves into 10% and the worst improves
  by half a point (−13.7 against −14.2); the read coefficient is not identified by the known-before-launch figures
  (H1 and `PL1` are the only stages with bases held, and the two-term fit gives them +0.2% and −0.1% by spending a
  second constant on them). Not smaller, no better.
- **Not usable, the two models after hashing**: they say what the cost follows (stored bases carry about 80% of it:
  3.415 of 4.253 per Gbp, *derived*), but the bases a stage stores, and its repeated records, are known only after
  ingest's pass over the sequences. S1 is the stage this costs: its 6,818 repeated records are counted as stored
  before launch, so its estimate is 14.2% high.
- **Chosen**: one rate per byte of FASTA not held, the least-squares fit 4.1100 (CPU-minutes − 0.07 per Gbyte not
  held: S1 3.522, P1 4.412, P2 3.938, P3 4.350, P4 4.248, P5 4.338, P6 4.106, H1 4.161, `PL1` 4.248). The interval of
  rates that keeps all nine within 20% is 3.665 to 4.411 (*derived*: P1's (4.611 / 1.2 − 0.07) / 1.0294 and S1's
  (7.540 / 0.8 − 0.07) / 2.1211). The rates that put eight within 10% are 4.004 to 4.376 (*derived*: P1's lower,
  P2's upper bound); S1 can be within 10% only below 3.917, where P1 is not, so eight is the most any one rate gives.

### The estimate before and after, per stage (`ingest_cpu_min`'s ratio, record / estimate)

| Stage | Read bytes | Held bases | Stored bases | CPU-min (record) | Old estimate | Old ratio | New estimate | **New ratio** |
|---|---|---|---|---|---|---|---|---|
| S1 | 2,121,124,779 | 0 | 1,664,109,938 | 7.540 | 9.553 | 0.789 | 8.788 | **0.858** |
| P1 | 1,029,384,352 | 0 | 1,016,625,803 | 4.611 | 4.672 | 0.987 | 4.301 | **1.072** |
| P2 | 9,716,599,273 | 0 | 9,596,180,328 | 38.335 | 43.508 | 0.881 | 40.005 | **0.958** |
| P3 | 10,074,413,081 | 6,963 | 9,949,077,522 | 43.889 | 45.108 | 0.973 | 41.476 | **1.058** |
| P4 | 10,101,497,437 | 7,053 | 9,976,456,950 | 42.979 | 45.229 | 0.950 | 41.587 | **1.033** |
| P5 | 10,092,827,773 | 185,683 | 9,967,711,632 | 43.855 | 45.190 | 0.970 | 41.551 | **1.055** |
| P6 | 10,087,447,605 | 17,198 | 9,962,536,103 | 41.486 | 45.166 | 0.919 | 41.529 | **0.999** |
| H1 | 1,962,116,403 | 141,683,767 | 1,819,831,338 | 7.645 | 8.842 | 0.865 | 7.552 | **1.012** |
| PL1 | 7,425,668,755 | 1,177,838,866 | 5,924,552,875 | 26.612 | 33.267 | 0.800 | 25.749 | **1.034** |

The new estimates are those the code printed on the real root (unit `b52-held`, *measured*: `held_bases` equal to the
table's for all nine). Records above their estimate: at most 7.2% (P1), so a process estimated at the 240 CPU-minute
margin would use at most 257 of 300 (*derived*: 240 × 1.072).

One more record, not a built stage: card C-36's rehearsal ingest of 300 phase-2 objects (`.gz`, 1,463,528,965 bytes
of FASTA, 3 contigs of 968 held, in a scratch root) took 6.48 CPU-minutes (*measured*, wilkcomp02, 2026-11-15,
`planning/measurements/2026-11-15_s2c_phase2_inputs.md`); the new estimate is at most 6.085 (*derived*: 0.07 + 4.11
× 1.4635, less the 3 contigs' bases), the record at least 6.5% over it. It was 6.61 before (2.0% under).

## 4. Before launch, in the code

- `held_bases(root, stage, tmp)`: one DuckDB query (one thread, the gene-call plan's 2 GB memory limit, spilling
  under the scratch directory) over the stage's `genome_input` files and every other stage's `contig_index` files on
  disk (`store.index_files`, the files ingest's anti-join reads): the sum of `contig_len` of the distinct listed
  contigs held as a default variant elsewhere. 0 before the stage has a `genome_input`.
- `_batched` (the step's decision, after the inputs step): `ingest_cpu_min(FASTA bytes, held_bases)`, said in one
  line either way, e.g. for `PL1` "7,425,668,755 bytes of FASTA, 1,177,838,866 bases of them in contigs another stage
  holds, are estimated at 25.7 CPU-minutes in one process (...)". `PL1` would still run in one process, as it did.
- `_refusal` (before any step, before `genome_input` exists): one object alone, every byte counted as stored, as
  before; the message names the new rate.
- What is not known before launch, and so not in the estimate: repeated records within the stage's own FASTA
  (S1: 20% of its bases read), and contigs repeated between objects of the same stage. Both lower the cost, so the
  estimate errs high for them; S1 is the worst case seen (14.2% high).

## 5. The batch cut

No built stage was ingested in batches, so there is no batch record to hold. The form is unchanged: `pass1_bytes`
cuts pass-1 shards to fit `SHARD_MARGIN` of the budget at the rate, every byte counted as stored; a batch is warned
when its bases read, every base counted as stored, pass the budget. A pass-1 shard reads and hashes and stores
nothing, so its estimate errs high (after hashing a base read costs about 0.84 of the 4.25 CPU-minutes per Gbp of a
base stored, *derived* above); a batch stores at most what it reads, so its estimate is within the fit's error.
`tests/r2/test_ingest_batches.py::test_b52_every_pass1_shard_and_batch_of_an_ingest_in_batches_fits_the_margin` sets
the rate so that S0 is over the margin in one process and each object under it, runs build-stage to the end of
ingest, and checks every pass-1 shard and both batches within the margin.

## 6. Tests

All on the pod, test image `bio-box-test_faec1f4.sif` (`scripts/pod_pytest.sh`), after the `PL1` build's unit had
ended (10:08 CST).

- Test first: the new tests on the commit that adds them alone (`8232fee`, the present `stage.py`), 10:10 CST: `13
  failed, 1 passed, 15 deselected`; every estimate test failed with `TypeError: ingest_cpu_min() got an unexpected
  keyword argument 'held_bases'` and the two build-stage tests with no `held_bases`; the batch-cut test passed (its
  form is unchanged). On the present function `PL1`'s ratio is 0.800 and S1's 0.789 (section 3).
- `tests/r1/test_b52_ingest_estimate.py`: each of the nine records within 20% of its estimate; exactly S1 outside 10%;
  the form without held bases. `tests/r2/test_ingest_batches.py`: `held_bases` on the fixture (lambda, which S0 holds,
  counted once though two genomes of S2 list it: 48,502 bases, the fixture's table); build-stage runs S2's ingest in
  one process in a root that holds S0 and in batches in one that does not, at a rate between the two estimates; every
  pass-1 shard and both batches of S0's ingest in batches within the margin.
- With the code: those files, `tests/r1/test_stage_estimates.py` (card C-13's three ingest points still within 20%)
  and `tests/r0`: `256 passed in 219.39s`. The whole suite in 2 shards with `BIOBOX_PEER_MANIFEST`: `1057 passed, 12
  skipped` and `1347 passed, 11 skipped`; the cross-host test passed against the peer's manifest hash
  `28ffa9808ecd4829e52f6892f94f7b85`. No existing test was changed.

## Appendix A: `b52_ingest_table.py` (run as `python -I b52_ingest_table.py /tmp/if2453/r10/root S1 P1 P2 P3 P4 P5 P6 H1 PL1`)

```python
"""Card B-52: what each built stage's ingest read and stored, beside its run record (read-only on the root).

python b52_ingest_table.py ROOT STAGE... > table.json

Per stage, in the order given (the order the stages were ingested):
  cpu_min, wall_s, started_at       the run record runs/stage=S/ingest.json
  objects, disk_bytes, fasta_bytes  the --objects of that record: files, bytes on disk, bytes of FASTA (decompressed)
  records, bases_read               contig_occurrence rows of S (the files the record lists) and the sum of their contig
                                    lengths (each record counted, repeats too)
  contigs_listed, bases_listed      distinct contig_uids of S's genome_input (the files inputs.json lists) and their
                                    lengths
  shared_contigs, shared_bases      of those, the ones a stage given earlier holds as a default variant in its
                                    contig_index (known before S's ingest starts: the inputs step and the root)
  stored_contigs, stored_bases      the default variants of S's own contig_index (the files the record lists)
One thread, no writes.
"""
from __future__ import annotations

import json
import sys
import zlib
from pathlib import Path

import pyarrow.parquet as pq


def fasta_bytes(p: Path) -> int:
    if not str(p).endswith(".gz"):
        return p.stat().st_size
    n = 0
    with open(p, "rb") as fh:
        d = zlib.decompressobj(31)
        while True:
            chunk = fh.read(1 << 22)
            if not chunk:
                n += len(d.flush())
                break
            while chunk:
                n += len(d.decompress(chunk))
                if d.eof:  # a multi-member gzip (bgzip): start the next member
                    chunk = d.unused_data
                    d = zlib.decompressobj(31)
                else:
                    chunk = b""
    return n


def listed(root: Path, rec: dict, prefix: str) -> list[Path]:
    return [root / o["path"] for o in rec["outputs"] if o["path"].startswith(prefix)]


def main() -> None:
    root = Path(sys.argv[1])
    lens: dict[bytes, int] = {}  # contig_uid -> length, default variants of the stages read so far
    out = []
    for s in sys.argv[2:]:
        rec = json.loads((root / "runs" / f"stage={s}" / "ingest.json").read_text())
        inp = json.loads((root / "runs" / f"stage={s}" / "inputs.json").read_text())
        objs = Path(rec["argv"][rec["argv"].index("--objects") + 1])
        disk = fa = n_obj = 0
        missing = 0
        for o in json.loads(objs.read_text()):
            p = Path(o["path"])
            p = p if p.is_absolute() else objs.parent / p
            n_obj += 1
            if not p.exists():
                missing += 1
                continue
            disk += p.stat().st_size
            fa += fasta_bytes(p)
        own: dict[bytes, int] = {}
        for f in listed(root, rec, f"l0seq/contig_index/stage={s}/"):
            t = pq.read_table(f, columns=["contig_uid", "contig_len", "is_default_variant"], use_threads=False)
            for u, n, d in zip(t["contig_uid"].to_pylist(), t["contig_len"].to_pylist(),
                               t["is_default_variant"].to_pylist()):
                if d:
                    own[u] = n
        listed_uids: set[bytes] = set()
        for f in listed(root, inp, f"l0seq/genome_input/stage={s}/"):
            listed_uids.update(pq.read_table(f, columns=["contig_uid"], use_threads=False)["contig_uid"].to_pylist())
        shared = [u for u in listed_uids if u in lens]
        allx = {**lens, **own}
        records = bases_read = 0
        unknown = 0
        for f in listed(root, rec, f"l0seq/contig_occurrence/stage={s}/"):
            for u in pq.read_table(f, columns=["contig_uid"], use_threads=False)["contig_uid"].to_pylist():
                records += 1
                if u in allx:
                    bases_read += allx[u]
                else:
                    unknown += 1
        out.append({
            "stage": s, "cpu_min": rec["cpu_min"], "wall_s": rec["wall_s"], "started_at": rec["started_at"],
            "threads": rec.get("threads"), "image_digest": rec.get("image_digest"),
            "objects": n_obj, "objects_missing": missing, "disk_bytes": disk, "fasta_bytes": fa,
            "records": records, "records_unknown_len": unknown, "bases_read": bases_read,
            "contigs_listed": len(listed_uids), "bases_listed": sum(allx.get(u, 0) for u in listed_uids),
            "shared_contigs": len(shared), "shared_bases": sum(lens[u] for u in shared),
            "stored_contigs": len(own), "stored_bases": sum(own.values()),
        })
        print(json.dumps(out[-1]), file=sys.stderr, flush=True)
        for u, n in own.items():
            lens.setdefault(u, n)
    json.dump(out, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
```

## Appendix B: `b52_held_check.py` (run as `python b52_held_check.py /tmp/if2453/r10/root table.json /tmp/if2453/b52/tmp`, `PYTHONPATH` the card's `src`)

```python
"""Card B-52: `stage.held_bases` of the code on the real root, read-only, one thread (an oracle check against the
table of b52_ingest_table.py, which counts only stages ingested earlier; PL1, now on disk, holds no contig of the
others' inputs that they do not hold first)."""
import json
import sys
from pathlib import Path

from biobox import stage as bs

root = Path(sys.argv[1])
table = {r["stage"]: r for r in json.load(open(sys.argv[2]))}
for s, r in table.items():
    held = bs.held_bases(root, s, Path(sys.argv[3]))
    est = bs.ingest_cpu_min(r["fasta_bytes"], held)
    print(json.dumps({"stage": s, "held_bases": held, "table_shared_bases": r["shared_bases"],
                      "equal": held == r["shared_bases"], "estimate": round(est, 3), "cpu_min": r["cpu_min"],
                      "ratio": round(r["cpu_min"] / est, 3)}), flush=True)
```

## Appendix C: `b52_fit.py` (run as `python3 -I b52_fit.py table.json`; arithmetic only)

```python
"""Card B-52: the fits of the ingest estimate over the nine records (arithmetic on the table read by b52_ingest_table.py).

python b52_fit.py table.json
"""
import json
import sys

F0 = 0.07  # INGEST_FIXED_CPU_MIN
T = {r["stage"]: r for r in json.load(open(sys.argv[1]))}


def fit_rel(xs):
    """Least squares of the relative error: min sum(((C - F0) - sum_i r_i x_i) / C)^2; x in Gbytes or Gbp."""
    rows = [[x / 1e9 / T[s]["cpu_min"] for x in xs[s]] for s in T]
    ys = [(T[s]["cpu_min"] - F0) / T[s]["cpu_min"] for s in T]
    k = len(rows[0])
    m = [[sum(r[i] * r[j] for r in rows) for j in range(k)] + [sum(r[i] * y for r, y in zip(rows, ys))]
         for i in range(k)]
    for i in range(k):
        p = m[i][i]
        m[i] = [a / p for a in m[i]]
        for j in range(k):
            if j != i:
                q = m[j][i]
                m[j] = [a - q * c for a, c in zip(m[j], m[i])]
    return [m[i][k] for i in range(k)]


def show(name, xs, coef):
    est = {s: F0 + sum(c * x / 1e9 for c, x in zip(coef, xs[s])) for s in T}
    off = {s: (T[s]["cpu_min"] - est[s]) / est[s] for s in T}
    print(f"{name}: coef {[round(c, 4) for c in coef]}; max |off| {max(map(abs, off.values())):.3f}; "
          f"within 10%: {sum(abs(o) <= 0.10 for o in off.values())}; within 20%: "
          f"{sum(abs(o) <= 0.20 for o in off.values())}")
    print("   " + " ".join(f"{s} {o:+.3f}" for s, o in off.items()))


F = {s: [r["fasta_bytes"]] for s, r in T.items()}
show("before B-52: 42.9 per 9,596,180,328 bytes read", F, [42.9 / 9.596180328])
show("one rate per byte read, refitted", F, fit_rel(F))
K = {s: [r["fasta_bytes"] - r["shared_bases"]] for s, r in T.items()}
show("one rate per byte not held (chosen form)", K, fit_rel(K))
show("one rate per byte not held, 4.11 (the code)", K, [4.11])
FK = {s: [r["fasta_bytes"], r["fasta_bytes"] - r["shared_bases"]] for s, r in T.items()}
show("per byte read plus per byte not held", FK, fit_rel(FK))
N = {s: [r["stored_bases"]] for s, r in T.items()}
show("after hashing: one rate per base stored", N, fit_rel(N))
RN = {s: [r["bases_read"], r["stored_bases"]] for s, r in T.items()}
show("after hashing: per base read plus per base stored", RN, fit_rel(RN))
print()
print("stage | read bytes | records | bases read | held before | stored contigs | stored bases | cpu_min | "
      "old est | old ratio | new est | new ratio | CPU-min per Gbyte not held")
for s, r in T.items():
    c = r["cpu_min"]
    old = F0 + r["fasta_bytes"] * 42.9 / 9596180328
    k = r["fasta_bytes"] - r["shared_bases"]
    new = F0 + k * 4.11 / 1e9
    print(f"{s} | {r['fasta_bytes']:,} | {r['records']:,} | {r['bases_read']:,} | {r['shared_bases']:,} | "
          f"{r['stored_contigs']:,} | {r['stored_bases']:,} | {c:.3f} | {old:.3f} | {c / old:.3f} | {new:.3f} | "
          f"{c / new:.3f} | {(c - F0) / (k / 1e9):.3f}")
```

# The ncRNA producer of card E-07 on the 8 genomes of 0-09c, 2026-11-11

Card E-07's pod check (Pod: yes). Every number is *measured* on wilkcomp02 on 2026-11-11 by the commands in
"Commands", unless labelled *derived* (formula given). Image `bio-box-test_74b3e05.sif` (sha256
`c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898`, Infernal 1.1.5); code: branch `task/E-07` at
`5c4da05` (its `src/` copied to `/tmp/if2453/e07/tree`); library `/tmp/if2453/golden/rfam/Rfam.cm.gz` (sha256
`e2e636e9…8168`, refdb.yaml's Rfam 15.1); model set `prok-1` (588 models, 23 whole-genome).

## Result

| quantity | value | label |
|---|---|---|
| producer | `infernal-1.1.5-rfam-15.1-prok1-igr200nt` | measured (the plan's and merge's output) |
| genomes, contigs, genome bases | 8, 8, 18,106,231 | measured (plan) |
| shards (plan at `--scan-bases 9500000`) | 2: 9,200,471 bases in 5 contigs, 8,905,760 in 3 | measured |
| intergenic records (regions) | 5,707 + 5,875 = 11,582 | measured (shard results); 0-09c's margin-200 count was 11,582 |
| hits (`l1_ncrna` rows) | 1,027 in 254 families | measured |
| truth: 0-09c margin-200 `--notrunc` hits of the set's intergenic models, plus its whole-genome producer-settings hits of the 23 whole-genome models | 1,027 (1,004 + 23) in 254 families | measured (`/tmp/if2453/m09c/out/merged/*.hits.jsonl`, read only) |
| truth hits found at the same contig, start, end, strand and model | 1,027 of 1,027 (recall 1.000); 0 producer hits outside the truth | measured (`e07_compare.py`) |
| CPU of the two shard steps (run records, user + system) | 11.41 + 16.25 = 27.65 CPU-minutes | measured |
| CPU of the plan and merge steps | 0.10 + 0.04 CPU-minutes | measured |
| CPU of the two shard units (systemd `CPUUsageNSec`, the whole unit with singularity and run-step) | 688.5 + 979.2 = 1,667.7 CPU-s | measured |
| **CPU-hours per 10 Gbp, shard steps** | **254.5** | derived: 27.65 × 60 / 18,106,231 × 1e10 / 3,600 |
| the same, every step (plan, shards, merge) | 255.8 | derived: 27.79 × 60 / 18,106,231 × 1e10 / 3,600 |
| against the design note's 247.4 | +2.9% (shard steps), +3.4% (every step) | derived; within the 20% rule (AGENTS §6) |
| per shard | 206.6 (shard 0), 304.0 (shard 1) CPU-hours per 10 Gbp | derived: each shard's CPU-s / its bases × 1e10 / 3,600 |
| wall time of the shards (two at once, `--cpu 2`) | 371.1 s and 526.1 s | measured (run records) |
| peak RSS of a shard step | 1.96 GB and 2.30 GB (run records); unit `MemoryPeak` 2.52 GB and 2.86 GB | measured; below the 6.4 GiB per two threads that the shard pool counts |

Hit recall against the whole-genome scan of the set's models (1,033 hits in 254 families, *measured* count from
0-09c's `wg_default.hits.jsonl`) is the design note's 1,027 of 1,033 (99.4%, *derived* there); this check does not
re-derive it, it shows that the producer reproduces the configuration's hits exactly.

## Reading

The producer, run as plan, two shards and merge under `run-step`, gives exactly the hits of 0-09c's margin-200
`--notrunc` scan restricted to the set plus the whole-genome hits of the in-CDS class, at 254.5 CPU-hours per
10 Gbp, 2.9% above the 247.4 the note derived from per-model CPU lines (the note's derivation summed per-model lines,
which ran 0.4% below process CPU in 0-09c). The two shards differ by a factor of 1.5 per base (206.6 and 304.0):
cost follows the number and length of intergenic records and the genomes' content (shard 1 holds the *E. coli*
and *B. subtilis* contigs and one more, by the plan's contig_uid ranges), so a shard cut by genome bases at 247.4
can run about 1.23 times its estimate (derived: 304.0 / 247.4); with the 80% margin (240 of 300 CPU-minutes) a shard at shard 1's rate would use about 295 CPU-minutes
(derived: 240 × 304.0 / 247.4), inside the budget but with little room. Memory stays far below the 6.4 GiB per two
threads the pool counts.

## Commands

On the laptop: the worktree's tree (without `.git` and `planning/`) was copied to `/tmp/if2453/e07/tree` with
`rsync -a --delete`; `/tmp/if2453/m09c/root` was copied to `/tmp/if2453/e07/root` (`cp -a`) so the steps could write.
The profile `/tmp/if2453/e07/e07-profile.yaml` is `profiles/wilkcomp02.yaml` with `tmp` under `/tmp/if2453/e07`,
`mem_scope: 16G` and both thread caps 4. Before each launch the account's `ps -u if2453 -o pcpu=` sum was 205% to
209%. Each step ran in its own unit, `MemoryMax=16G`:

```bash
P=infernal-1.1.5-rfam-15.1-prok1-igr200nt
RUN="bash /tmp/if2453/e07/tree/scripts/pod_run.sh --profile /tmp/if2453/e07/e07-profile.yaml --image /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif --bind /tmp/if2453 --env PYTHONPATH=/tmp/if2453/e07/tree/src:/tmp/if2453/pilot_s1/pydeps --mem 16G"
STEP="run-step --profile /tmp/if2453/e07/e07-profile.yaml --root /tmp/if2453/e07/root --stage M9C"
NC="python -m biobox.cli stage ncrna --root /tmp/if2453/e07/root --stage M9C --cm /tmp/if2453/golden/rfam/Rfam.cm.gz"
$RUN --threads 2 e07-plan2 $STEP --step ncrna-$P.plan --threads 1 -- $NC plan --scan-bases 9500000 --tmp /tmp/if2453/e07/tmp/plan
$RUN --threads 2 e07-shard-0 $STEP --step ncrna-$P.shard-00000 --threads 2 -- $NC shard --shard 0 --cpus 2 --tmp /tmp/if2453/e07/tmp/scan --seq-stages M9C
$RUN --threads 2 e07-shard-1 $STEP --step ncrna-$P.shard-00001 --threads 2 -- $NC shard --shard 1 --cpus 2 --tmp /tmp/if2453/e07/tmp/scan --seq-stages M9C
$RUN --threads 1 e07-merge $STEP --step ncrna-$P --threads 1 -- $NC merge --tmp /tmp/if2453/e07/tmp/scan
```

(`$NC plan` is `stage ncrna plan ...`: the action follows `stage ncrna`; the variables are shorthand for this note.)
The two shard units ran at the same time: two `cmsearch` at most, two threads each. A first plan unit, `e07-plan`,
was refused because its `--tmp` directory did not exist; the producer now creates it (`ncrna._mkscratch`, in the
code above). The comparison, in the image in a `MemoryMax=4G`, `CPUQuota=100%` scope at nice 19:
`python -I e07_compare.py /tmp/if2453/e07/root M9C $P prok-1.tsv /tmp/if2453/m09c/out /tmp/if2453/e07/compare.json`
(the script is in `/tmp/if2453/e07/`). Outputs stay on the pod under `/tmp/if2453/e07`; nothing was written under
`/stor` or the m09c, r10 and r50 directories.

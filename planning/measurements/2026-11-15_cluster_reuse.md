# A release takes an earlier release's clustering (card B-49, DEC-45 (2))

Pod wilkcomp02, 2026-11-15 (times CST, from the pod's clock), runtime image `bio-box_faec1f4.sif`
(`sha256:d67300aa…`), the tree of branch `task/B-49` (synced to `/tmp/if2453/b49/tree`, its `src` first on the path).
The data root `/tmp/if2453/r10/root` was only read: it was bound read-only into every command, and the scratch view
`/tmp/if2453/b49/view` held one symlink per file of the root's stages S1, P1 to P6 and of releases `rM.1`, `rM.2`
(`scripts/stress/build_stage_pod/b49_view.py`; 3,199 links). Driver: `scripts/stress/build_stage_pod/b49_driver.sh`.

## The rule

What is compared, what is written, what names the source: see the module docstring of `src/biobox/l2/cluster.py`
("Reuse"). In short: R1's plan record must equal the plan record R2's clustering would write in every field but
`release` (stages, `protein_seq` files by path and size, the distinct protein uids by the sha256 of `prot.parquet`,
MMseqs2, cut, format, counts); R1's clustering must be finished and its ledger must name today's MMseqs2 and pyarrow.
R2 gets copies of R1's tables (sha256 checked against R1's run record), R1's step records with `release` and
`plan_sha256` made R2's, a ledger with R1's rows plus one `reused:<R1>:<inputs_sha256>` row per level, and the record
`l2/cluster_plan/release=<R2>/cluster.reuse.json`; the run record of R2's step `cluster` names R1 in its argv.

## Plans of `rM.1` and `rM.2` (all 443 records)

*Measured* (`b49_compare_plans.py /tmp/if2453/r10/root rM.1 rM.2`, in the image, read-only): both have 443 records
with the same names. Beyond `release` and `plan_sha256`, 17 records differ: the plan's `cut` (`rM.1`: `{"parts": 24}`;
`rM.2`: `{"one_process_cpu_min_max": 20000.0, "parts": 24}`) and the 16 cut records' `cpu_min` estimates (the same
part boundaries, keys and units; `rM.1`'s were made by the cost model before card C-22's refit). Both clustering ledgers
name `mmseqs2 18.8cc5c` and `pyarrow 26.0.0`, as the runtime image of `faec1f4` does. So the earlier statement that the
plans differ "only in the release's name and the plan hash" holds for 426 of 443 records; under the rule `rM.2` is a
source for a release planned today and `rM.1` is not (its cut differs).

## The scratch release `rB49.1` over the seven stages of `rM.2`

*Measured* (unit `b49-build`, `b49_driver.sh build`, 14:37:47 to 15:40:40; `--threads 4`, unit quota 6 CPUs, 32G):
`build-stage --stage S1 --release rB49.1 --stages S1,P1,P2,P3,P4,P5,P6 --cluster-parts 24` with the wilkcomp02 profile
(so the plan's cut is `rM.2`'s, `{"one_process_cpu_min_max": 20000.0, "parts": 24}`). The six stage steps were
passed over (their run records are linked). Before the cluster step build-stage compared the plans (14:37:57 to
14:40:13, with the DEC-26 owner check and the protein uids of 45,815,843 proteins): `rM.1` refused (its cut), `rM.2`
taken. Then:

| Step | CPU-minutes | wall s | peak RSS | note |
|---|---|---|---|---|
| cluster (`release cluster reuse --from rM.2`) | 3.29 | 219.1 | 0.70 GB | 11 table files, 10,587,564,216 bytes copied; 443 records; published at 14:47:38 (rsync to the view) |
| tokens | 22.01 | 1,395.2 | 2.48 GB | `nr50` 8,192,065 and `nr90` 31,183,876 tokens, 19,036 family |
| postings | 21.65 | 1,381.5 | 2.63 GB | |
| dim | 0.06 | 3.1 | 0.34 GB | |
| manifest | 4.21 | 286.7 | 0.87 GB | |
| catalog | 0.09 | 15.1 | 0.34 GB | |

(run records `runs/stage=rB49.1/*.json`, `b49_check.py`). The whole unit: 3,317.1 CPU-seconds, 55.3 CPU-minutes
(*measured*, `CPUUsageNSec`), of which the release steps' records hold 51.3 (*derived*: the sum above) and build-stage's
own comparison and checks the rest, about 4.0 (*derived*). The unit's `MemoryPeak` read 32 GiB, its cap: it counts page
cache (card C-15's caveat), and no step's own peak passed 2.63 GB. Against `rM.2`'s clustering, 14,086.8 CPU-minutes over
444 run records (*measured*, the sum of `runs/stage=rM.2/cluster*.json`), the cluster step cost 3.29 CPU-minutes (*derived*:
0.023%) and 219 s plus its publish, about 7 minutes from launch to published (14:40:24 to 14:47:38).

Checks (*measured*, `b49_check.py`, `release verify`, in the image, the root read-only):
- every one of the 11 copied files has the sha256 that `rM.2`'s run record lists, and is a real file of the view with one
  link (not a symlink, not a hard link);
- tokens (5 files), postings (16) and dim (2) of `rB49.1` have the same sha256 as `rM.2`'s files of the same names
  (their run records): the release steps on top of the copied clusters give `rM.2`'s bytes (no code of `l3/`, `dim.py`
  or `io.py` changed since `rM.2`'s tree `5ba68a6`);
- the ledger holds `rM.2`'s two rows and two `reused:rM.2:1d0d49fc…` rows; the reuse record names `rM.2`, its plan
  `27cfd07b…` and `rB49.1`'s `18a71fe3…`, both with inputs hash `1d0d49fc…`;
- `release verify`: "release rB49.1: verified (stages S1,P1,P2,P3,P4,P5,P6, 206 files, 154 artifacts, manifest hash
  b727a588617df111cf8e7cbcc8a91816)".

Disk: free on `/tmp` 635,499,229,184 bytes before (*measured*, `df`); the view's real files reached 18,672,170,438 bytes
(*measured*); its peak was not sampled: at most about 21.2 GB while the cluster step's 10.6 GB were in staging and in the
view at once, and at most 26.8 GB during the postings publish (*derived*: 18.7 GB at the end plus at most the 8.1 GB the
later steps wrote, held twice). The view and scratch were removed by 16:17; `/tmp/if2453/b49` keeps 34 MB (logs, tree,
scripts).

## Partial reuse (a release over a superset of the stages)

Not done: such a release is clustered. To take R1's clusters for its proteins and cluster only the new ones would need
(1) MMseqs2's `clusterupdate`-like assignment, which DEC-24 rejected because the result depends on the order stages
arrived; or (2) evidence that a from-scratch clustering of the superset keeps R1's clusters of R1's proteins: a
measurement of how many members move (card C-35 asks it for the nine stages against `rM.2`), and a rule for which
cluster uids may then be copied (a cluster with no new member and no member moved), checked against a full clustering.

## The tokens refusal (card item 5): not done, stopped

Seen building silently on `c9e336a` (*measured*, pod, `release define` then `release tokens` on a release never
clustered: exit 0, no message, tokens without cluster levels). A refusal in `tokens.run` or `release tokens` would break
existing tests on the GitHub `ci` job, whose runner has no MMseqs2: `tests/r2/test_builders_out.py` (`RELEASE_STEPS`
without `cluster` when MMseqs2 is absent, and its `tokens` step) and `tests/r2/test_micro_s0_release_def.py`
(`build_release`) build tokens with no clustering by design there. The card names no migration, so the refusal waits.

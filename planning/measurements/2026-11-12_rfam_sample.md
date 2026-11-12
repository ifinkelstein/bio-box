# Which long-tail Rfam models hit on a sample of the corpus, 2026-11-12 (card 0-08a)

Card 0-08a, M1 of `planning/design/0-08_ncrna_models.md` (DEC-28's follow-up): full Rfam 15.1 on whole contigs of
a 100-genome prokaryote sample of release `rN.3` (stages P1 and P2, 2,497 genomes, stratified by GTDB R226
phylum) plus 300 phages of stage S1, with the producer's `cmsearch` settings, to decide which long-tail models join
the frozen set. Every number is *measured* on wilkcomp02 in the image `bio-box_80d4c55.sif` (Infernal 1.1.5,
Rfam 15.1 `Rfam.cm.gz`, sha256 `e2e636e9…8168`, equal to `refdb.yaml`; checked on the pod, `lib/Rfam.cm.gz.sha256`)
by the commands below, unless labelled *derived* (formula given) or *estimate*. The scan started 2026-11-11
21:58 CST (2026-11-12 UTC, the run records' dates).

**Status when this file was first committed: the scan is running** (an unattended driver on the pod; 4 of 27
processes done). The results section is written by `rfam_sample_report.py` from the scan's outputs; see "How to
finish" at the end.

**Measured against the estimate (AGENTS §6).** The note's M1 cost is 81 CPU-hours (*derived* there from 1,903
CPU-hours per 10 Gbp × an *estimated* 425 Mbp). The sample is 412,032,974 bases (*measured*), and the first 4
processes ran at 1,702 to 2,108 bases per CPU-second (*measured*, `runs/s00.json` to `s03.json`), against 1,454 in
0-09c: projected over the whole sample about 60 CPU-hours, about 25% under the estimate (*derived*: 9.21 CPU-hours
/ 62,675,593 bases × 412,032,974). The deviation is in the cheap direction and per-genome rates differ by 24% among
these four processes; the final figure is in the summary table.

<!-- results:begin -->

(Written by `rfam_sample_report.py` when the scan is done.)

<!-- results:end -->

## Method

**The sample** (`rfam_sample.py sample`; reproducible from the seed 20261112 and the sorted accession lists):

- Prokaryotes: the 2,497 genomes of stages P1 (250) and P2 (2,247) of `/tmp/if2453/r10/root`, read through the
  current code (`genome_input`, `contig_index`). A genome's accession is the file name of the source object at
  the same index of its stage's `objects.json` as the genome's entry in `genomes.json` (the two lists were made
  together by `make_inputs.py` and have equal lengths); for every sampled genome this pairing was checked by
  hashing its FASTA records to the genome's contig uids. Each genome's GTDB R226 phylum is read from
  `bac120_metadata_r226.tsv.gz` and `ar53_metadata_r226.tsv.gz` (`gtdb_taxonomy`), matched by accession (all
  2,497 matched exactly). Allocation of the 100 over the 55 phyla present (8 archaeal): one each, the remaining
  45 in proportion to the phylum's genome count (largest remainder, ties by phylum name, never more than a phylum
  has); within a phylum `random.Random("0-08a:20261112:<domain>:<phylum>").sample(<sorted accessions>, k)`.
- Phages: the 27,258 genomes of S1 sorted by `genome_uid` (hex), every 90th from the first (`27,258 // 300`).
- The sample's 676 contigs (none circular) were exported to `/tmp/if2453/m08a/fa/s<k>.fna` (`rfam_sample.py
  fasta`, through `SeqStore`; every record's BLAKE3 contig uid rechecked against its name: 676 records,
  412,032,974 bases, 0 mismatches).
- Gene calls on those contigs (`l1/feature`, producer `pyrodigal-3.7.1-gv-0.3.2-p2`, stages S1, P1, P2) and the
  Pfam 38.2 domains of their proteins (`l1/feature_protein`, `l1/domains`) were exported to
  `sample/genes.jsonl` (`rfam_sample.py genes`): 379,579 genes on all 676 contigs, 289,697 with a Pfam domain
  (*measured*).

**The plan.** Genomes in the order prokaryotes by accession, then phages by genome_uid, packed greedily into
processes of at most 18,000,000 bases (a genome never split): 27 processes of 8,981,933 to 17,990,498 bases
(`sample/plan.json`). At 1,450 bases per CPU-second (0-09c, *measured*) the largest is 207 CPU-minutes, and 227
at E-01's slower 1,319 (*derived*), under the 300 CPU-minute limit (AGENTS §7.2).

**The scan** is the producer's whole-contig argv (`biobox.l1.ncrna._cmsearch_argv`: `cmsearch --cpu 2 -Z 100000
--cut_ga --noali --tformat fasta --tblout <tbl>`; not `--notrunc`, not intergenic) with `-o <file>` instead of
`/dev/null`, so that each model's own `# CPU time` line is kept, over all 4,227 models of `Rfam.cm`, one process
per planned part. Hits are mapped to contigs by `infernal_igr.py collect` (the producer's rules) and per-model CPU
read by `infernal_igr.py permodel` (both from card 0-09c). CPU is user + system of the `cmsearch` process
(`timed.py`).

**Pod limits kept.** Each process is a unit `m08a-s<k>` with `MemoryMax=16G`, `MemorySwapMax=0`,
`CPUQuota=200%`, `nice 19`, in the image; at most 2 at a time (4 threads), beside the 50 Gbp build. Before every
launch the driver sums `ps -u if2453 -o pcpu=` (the ps line left out) and launches only when that plus 400% is at
most 4000% (the account's day cap, and its cap whenever the load average is above 8); otherwise it polls every 4
minutes. Nothing is written under `/tmp/if2453/r10/root` or `/stor`; the driver only reads for a `STOP` file there.

## Commands

On the laptop: `git -C <worktree> checkout -b measure/0-08a origin/plan/task-cards` (at `4648286`); the tree was
mirrored to `/tmp/if2453/ci/m08a` by `sh scripts/pod_pytest.sh m08a tests/r0/test_runbook.py` (34 passed), and on
the pod copied to the frozen `/tmp/if2453/m08a/tree` (`rsync -a --delete`), from which every step runs. `IMG` is
`/tmp/if2453/images/bio-box_80d4c55.sif`, `BIN` is `/tmp/if2453/m08a/tree/scripts/stress/annotation_cost`, and
each python step runs as
`systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% nice -n 19 singularity exec
--cleanenv --bind /tmp/if2453 --env PYTHONPATH=/tmp/if2453/m08a/tree/src $IMG python $BIN/timed.py <record> -- python ...`:

```bash
sha256sum /tmp/if2453/golden/rfam/Rfam.cm.gz > /tmp/if2453/m08a/lib/Rfam.cm.gz.sha256
gunzip -c /tmp/if2453/golden/rfam/Rfam.cm.gz > /tmp/if2453/m08a/lib/Rfam.cm
python $BIN/rfam_sample.py sample /tmp/if2453/r10/root /tmp/if2453/r10/inputs /tmp/if2453/b18/gtdb /tmp/if2453/m08a   # 0.66 CPU-min
python $BIN/rfam_sample.py fasta  /tmp/if2453/r10/root /tmp/if2453/m08a /tmp/if2453/m08a/lib/Rfam.cm                   # 0.16 CPU-min
python $BIN/rfam_sample.py genes  /tmp/if2453/r10/root /tmp/if2453/m08a                                                 # 0.23 CPU-min
systemd-run --user --unit=m08a-driver --collect --quiet -p MemoryMax=4G -p MemorySwapMax=0 -p CPUQuota=100% \
  -p StandardOutput=append:/tmp/if2453/m08a/driver.log -p StandardError=append:/tmp/if2453/m08a/driver.log \
  nice -n 19 bash $BIN/rfam_sample_driver.sh
```

The driver (`rfam_sample_driver.sh`) launches each process as `rfam_sample_shard.sh <k>` in unit `m08a-s<k>`:
`timed.py runs/s<k>.json -- cmsearch --cpu 2 -Z 100000 --cut_ga --noali --tformat fasta --tblout out/s<k>.tbl -o
out/s<k>.out lib/Rfam.cm fa/s<k>.fna`, then `infernal_igr.py collect out/s<k>.tbl fa/s<k>.map.json
out/s<k>.hits.jsonl` and `infernal_igr.py permodel out/s<k>.out out/s<k>.permodel.json`, then `out/s<k>.ok`. When
all 27 are done it runs `rfam_sample_report.py /tmp/if2453/m08a /tmp/if2453/m08a/tree /tmp/if2453/m08a/report`
and writes `DONE`.

## How to finish (if the scan was still running when this was committed)

1. Is it done? On the pod: `ls /tmp/if2453/m08a/DONE /tmp/if2453/m08a/FAILED`; progress is
   `tail /tmp/if2453/m08a/driver.log` and `ls /tmp/if2453/m08a/out/*.ok | wc -l` (27 when complete);
   `systemctl --user list-units 'm08a-*'` shows the driver and the running processes.
2. `FAILED` names the failed processes (their logs: `logs/s<k>.log`; `out/s<k>.fail` holds the exit code). Run the
   driver again (the `systemd-run` line above): done processes are skipped, failed ones retried. To pause: `touch
   /tmp/if2453/m08a/HALT` (running processes finish; remove the file and relaunch the driver to resume). To use 4
   processes once the build has finished and the pod is otherwise idle: `echo 4 > /tmp/if2453/m08a/MAXP`.
3. When `DONE` exists, `/tmp/if2453/m08a/report/` holds `report.json`, `results.md`, `0-08a_prok-2_proposed.tsv`
   and this file with its results filled in (`2026-11-12_rfam_sample.md`). Copy the last two into the repository
   (`planning/measurements/` and `planning/design/`), write the reading paragraph below, and commit. To rebuild
   the report by hand (also mid-scan with `--partial`), run the scope line above with `python
   $BIN/rfam_sample_report.py /tmp/if2453/m08a /tmp/if2453/m08a/tree /tmp/if2453/m08a/report [--partial]`.

## Reading

(To be written from the final results, and marked as the writer's interpretation.)

## What was not measured, and why

- Intergenic cost of the added models: the scan is whole genome (the truth scan); the cost of prok-2 is *derived*
  from the 8-genome per-model lines of 0-09c, as prok-1's was.
- Whether each hit of an excluded or long-tail model is a real RNA: the report gives coordinates, score against
  GA, overlapping gene calls with their Pfam domains and overlapping prok-1 hits; the judgement is the PI's.
- Clan competition is not applied: every hit at or above GA is counted (card E-08 applies it). The report says
  which candidate hits overlap a prok-1 hit of the same clan.
- Genomes outside P1 and P2, metagenome contigs, and stage S1's phages beyond the 300 sampled.

## Files

Scripts in `scripts/stress/annotation_cost/`: `rfam_sample.py` (sample, export, gene calls),
`rfam_sample_driver.sh` and `rfam_sample_shard.sh` (the scan), `rfam_sample_report.py` (every table), reusing
`infernal_igr.py` and `timed.py`. Outputs stay on the pod under `/tmp/if2453/m08a` (`sample/`, `fa/` the exported
FASTA, `runs/`, `out/` tables, gzipped `-o` outputs, hits and per-model CPU, `logs/`, `report/`).

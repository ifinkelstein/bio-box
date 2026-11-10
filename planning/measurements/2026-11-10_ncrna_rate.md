# ncRNA scan rate with full Rfam 15.1, 2026-11-10 (card E-01)

Every number here is *measured* on wilkcomp02 (112 cores; load average 13 to 28 from other users' jobs
during the runs) on 2026-11-10, in the pinned image `bio-box-test_74b3e05.sif` (Infernal 1.1.5, as
`container/environment.yml` pins it), unless it is labelled *derived* or *estimate*. Each run is its own
`systemd-run --user --unit=e01-<name> --collect -p MemoryMax=32G -p MemorySwapMax=0` unit at `nice -n 19`, with
`singularity exec --cleanenv --env PYTHONPATH=<tree>/src:/tmp/if2453/pilot_s1/pydeps --bind /tmp/if2453`, timed
by the host's `/usr/bin/time -v` around `singularity exec` (the image has no `time`; user and system time
include every waited-for child, so cmsearch's CPU is counted; peak RSS is the largest single process).
The scripts are in `/tmp/if2453/e01/` on the pod and are reproduced in the appendix; they are not in the
repository (the card lists no script file). Nothing was written under `/stor`.

## Library

Rfam 15.1 `Rfam.cm.gz` from the data root's reference area
(`/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/rfam/15.1/`, read only), copied to `/tmp/if2453/e01/rfam/`:
45,753,387 bytes, sha256 `e2e636e9ce138dae506769aa74d067ce35466e94002818cfe125ea2a3cab8168`, equal to
`refdb.yaml`. 4,227 models, every one with a GA line; the largest window `W` is 4,139 and the largest
filter-HMM `MAXL` 5,746 (both RF02543, LSU_rRNA_eukarya), so the builder's circular flank is 5,746 bases.
So the producer is `infernal-1.1.5-rfam-15.1`.

## Probe: the micro-S0 fixture (one cmsearch call)

`cmsearch --cpu 8 -Z 100000 --cut_ga --noali --tblout ... Rfam.cm all.fna` over lambda, T4 and
*M. genitalium* (797,481 bp, 3 contigs), the binary alone (no builder):

| bases | user | system | CPU | wall | CPU/wall | peak RSS |
|---|---|---|---|---|---|---|
| 797,481 | 1,010.80 s | 19.81 s | 1,030.61 s | 446.91 s | 2.31 | 2,144,164 kB |

Rate 773.8 bases per CPU-second (*derived*: 797,481 / 1,030.61). 58 hits at or above GA in 16 families;
the fixture library of the tests (`tests/fixtures/micro_s0/rfam15.1_micro_s0.cm.gz`, nine families) was
chosen from them. Overlaps seen: the four SSU and three LSU rRNA models of the rRNA clans all hit the one
rRNA operon of *M. genitalium* (each above its GA), and tRNA-Sec (RF01852) hits a tRNA locus that tRNA
(RF00005) also hits (seq 257,558-257,647 and 257,559-257,647, 1-based); the builder keeps every one (no clan
rule in the spec or the decision rows).

## The golden set through the builder

The 15 golden genomes (`/tmp/if2453/golden`, `golden.lock.yaml`; 36,475,789 bp in 21 contigs) were ingested
into one root as six stages, so that no single process passes 300 CPU-minutes (AGENTS §7.2; the probe's rate
put the whole set near 785 CPU-minutes, *derived*). Then each stage was scanned by one builder process,
`python -m biobox.cli stage ncrna --root R --stage G<k> --cm Rfam.cm.gz --cpus 8 --tmp <scratch>` (tree
`task/E-01` at 5ed1781), in one unit `e01-rate`, one stage after the other, 09:47 to 11:15 CST. Each stage is
one cmsearch call (`batches` 1). CPU is user + system of the process and its children.

| stage | genomes | contigs | bases | hits | user s | system s | CPU s | wall | CPU/wall | peak RSS kB |
|---|---|---|---|---|---|---|---|---|---|---|
| G1 | btheta_vpi5482, lambda, t4, obelisk_ss | 6 | 6,513,108 | 176 | 4,513.78 | 87.38 | 4,601.16 | 14:47.86 | 5.18 | 9,766,768 |
| G2 | macetivorans_c2a, mgenitalium, lak_a1, lak_b8 | 4 | 7,424,839 | 293 | 6,240.43 | 112.25 | 6,352.68 | 18:44.90 | 5.65 | 9,194,132 |
| G3 | ecoli_mg1655, sr1_raac1, spyogenes_sf370 | 3 | 7,671,845 | 663 | 5,975.92 | 131.37 | 6,107.29 | 18:42.53 | 5.44 | 12,903,916 |
| G4 | bsub_168, sacidocaldarius_dsm639 | 2 | 6,441,565 | 368 | 5,107.69 | 116.04 | 5,223.73 | 16:48.96 | 5.18 | 11,808,300 |
| G5 | mtb_h37rv | 1 | 4,411,532 | 92 | 2,858.98 | 53.09 | 2,912.07 | 10:06.67 | 4.80 | 1,892,712 |
| G6 | hvolcanii_ds2 | 5 | 4,012,900 | 147 | 2,410.46 | 40.99 | 2,451.45 | 8:25.69 | 4.85 | 2,793,484 |
| all | 15 | 21 | 36,475,789 | 1,739 | 27,107.26 | 541.12 | 27,648.38 | 1:27:36.61 | 5.26 | 12,903,916 (largest) |

- **Rate: 1,319 bases per CPU-second** over the golden set (*derived*: 36,475,789 / 27,648.38 = 1,319.3;
  per stage 1,169 (G2) to 1,637 (G6), *derived*). Equivalently 758 CPU-seconds per Mbp.
- **Peak RSS: 12.3 GiB** (12,903,916 kB, G3), the largest single process, which is cmsearch: the builder's
  Python holds one contig (at most 6.26 Mbp here) and DuckDB is capped at 2 GB. It does not follow contig
  length: the single 4.41 Mbp contig of G5 peaked at 1.9 GB, while the stages with the most rRNA operons
  (G3 with *E. coli*, G4 with *B. subtilis*) peaked highest. That is an *estimate* of the cause (the
  alignment of long rRNA hits by up to 8 worker threads at once), not measured.
- **Threads used**: 5.26 CPU-seconds per wall second at `--cpu 8` (*derived*); on the three-contig probe,
  2.31.
- **The largest process** used 105.9 CPU-minutes (G2, *derived*: 6,352.68 / 60), under the 300-minute limit.
- The probe's lower rate (774 against 1,319) is consistent with a fixed cost per cmsearch call of reading
  and configuring all 4,227 models: fitting `CPU = a + b x Mbp` to the probe and G5 gives a = 615 CPU-s and
  b = 521 CPU-s per Mbp (*derived*, two points of different content, the probe being the bare binary; the
  linear model is an *estimate*). The builder therefore searches
  large batches (`BATCH_BASES` = 50 Mbp, chosen, not tuned).

## Projection (*derived* from the golden rate of 1,319.3 bases per CPU-second)

| bases | CPU-hours | wall at 40 threads, as used here (5.26 of 8) | processes under 300 CPU-min |
|---|---|---|---|
| 10 Gbp | 2,105 (1e10 / 1,319.3 / 3,600) | 80 h (2,105 / (40 x 5.26 / 8)) | at least 422 (23.7 Mbp each) |
| 50 Gbp | 10,527 (5e10 / 1,319.3 / 3,600) | 400 h | at least 2,106 |

The 300 CPU-minute limit allows about 23.7 Mbp per process (*derived*: 18,000 x 1,319.3). The builder scans
a stage in one process, so a stage above about 24 Mbp needs the scan split into shard processes as card B-24
did for the domain scan; that is not in card E-01 (see its report). Golden genomes are not a random sample:
the rate on metagenome or phage contigs (more, shorter contigs; fewer rRNA operons) is not measured.

## Appendix: scripts (on the pod, `/tmp/if2453/e01/`)

- `probe.sh`: unzips the three fixture FASTA files into `fx/all.fna` and runs the probe's cmsearch under
  `/usr/bin/time -v`.
- `cmstats.py`: reads every CM header of `Rfam.cm` (models, largest `W` and `MAXL`, models without GA).
- `mklib.sh`: `cmfetch -o lib/rfam15.1_micro_s0.cm -f rfam/Rfam.cm lib/keys.txt` for the nine accessions,
  then `gzip -n -9`: the test fixture.
- `measure/prepare_golden.py`: ingests the golden genomes as stages G1 to G6 (`ingest`, `write_genome_input`,
  genetic codes, kinds and topologies from `golden.yaml`), printing `measure/prepared.json`.
- `measure/rate.sh`: for each stage, `/usr/bin/time -v -o logs/G<k>.time singularity exec ... python -m
  biobox.cli stage ncrna --root measure/root --stage G<k> --cm rfam/Rfam.cm.gz --cpus 8 --tmp measure/scratch`,
  skipping a stage whose `logs/G<k>.json` exists and stopping at a `measure/STOP` file.
- `measure/collect.py`: the table above from `logs/`.

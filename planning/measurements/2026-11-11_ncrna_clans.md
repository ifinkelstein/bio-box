# Clan competition for ncRNA hits (card E-08) on the 8 genomes of 0-09c, 2026-11-11

Card E-08's pod check (Pod: yes) and what the Infernal binaries say about clans. Every number is *measured* on
wilkcomp02 on 2026-11-11 by the commands in "Commands", unless labelled *derived* (formula given). Image
`bio-box-test_74b3e05.sif` (sha256 `c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898`, Infernal
1.1.5); code: branch `task/E-08` (its tree copied to `/tmp/if2453/e08/tree`); library
`/tmp/if2453/golden/rfam/Rfam.cm.gz` (refdb.yaml's Rfam 15.1); model set `prok-1`; clan file
`src/biobox/l1/rfam_sets/Rfam-15.1.clanin` (refdb.yaml `rfam_clanin`, 151 clans).

## Result

| quantity | value | label |
|---|---|---|
| producer | `infernal-1.1.5-rfam-15.1-prok1-igr200nt-clan1` | measured (plan and merge output) |
| shards (plan at `--scan-bases 9500000`) | 2: 9,200,471 bases in 5 contigs, 8,905,760 in 3 (as card E-07's) | measured |
| hits (`l1_ncrna` rows) | 1,027 | measured |
| truth hits (0-09c, as card E-07's check) found at the same contig, start, end, strand and model | 1,027 of 1,027; no other row | measured (`e08_compare.py`) |
| hits that lose clan competition | **52**: 25 in CL00111 (SSU), 25 in CL00112 (LSU), 1 in CL00003 (SRP), 1 in CL00011 (Glm) | measured |
| hits that win, or whose model is in no clan | 975 (465 of them in no clan) | measured |
| features (`feature` rows) | 975, the distinct loci of the winners; card E-07's producer wrote 1,027 on the same stage | measured |
| clan flags against the rule's loop written in the compare script (base sets, rank order) | 0 differences | measured |
| competing groups (one clan, strand and connected overlap) with more than two hits | 0: no chain | measured |
| SSU loci (connected overlapping CL00111 hits on one strand) | 25, each 2 hits (SSU_rRNA_bacteria, SSU_rRNA_archaea), each with exactly 1 winner | measured |
| LSU loci (CL00112) | 25, each 2 hits (LSU_rRNA_bacteria, LSU_rRNA_archaea), each with exactly 1 winner | measured |
| rRNA operons by genome (SSU loci = LSU loci) | B. subtilis 10, E. coli 7, S. pyogenes 6, M. tuberculosis 1, S. acidocaldarius 1 | measured |
| losers outside rRNA | B. subtilis Bacteria_small_SRP 26533-26633 (+, 76.7 bits) under Bacteria_large_SRP (250.4); E. coli GlmY_tke1 3986433-3986582 (+, 66.8) under GlmZ_SraJ (218.4) | measured |
| ledger, summed over the 8 contigs | n_in 1,027, n_kept 975, n_merged 52 (the losers); `gene_call_version` `pyrodigal-3.7.1-gv-0.3.2-p2`, `rfam_clanin` 15.1 among the tool versions | measured |
| CPU of the two shard steps (run records) | 11.53 + 16.83 = 28.36 CPU-minutes | measured |
| CPU-hours per 10 Gbp, shard steps | 261.0 (shard 0: 208.9, shard 1: 314.9) | derived: CPU-minutes × 60 / bases × 1e10 / 3,600 |
| the same in card E-07's run of the same shards | 254.5 (206.6, 304.0) | derived there |
| CPU of the plan and merge steps | 0.08 + 0.05 CPU-minutes | measured |
| peak RSS of a shard step | 2.01 GB and 2.29 GB (run records); unit `MemoryPeak` 2.57 GB and 2.86 GB | measured |

The merge, which now decides clan competition, costs what card E-07's merge did (0.05 against 0.04 CPU-minutes). The
shard steps ran 2.6% above card E-07's (28.36 against 27.65 CPU-minutes, *derived*); the scan is the same code, so
this is the host's variation between two runs. Shard 1 ran at 314.9 CPU-hours per 10 Gbp, 3.6% above the 304.0 this
card cuts shards at (*derived*: 314.9 / 304.0); a shard cut at 80% of 300 CPU-minutes at 304.0 would have used 248.6
CPU-minutes at that rate (*derived*: 240 × 314.9 / 304.0), inside the budget.

## Reading

**One kept SSU and one kept LSU feature per operon: yes**, on all 25 operons. With the set `prok-1` each rRNA gene is
hit by two models (bacterial and archaeal; the eukaryotic models are out of the set, DEC-28 (1)), not four and three
as with the full library (design note §5); the better-scoring model wins every time (the bacterial one on bacteria,
the archaeal one on S. acidocaldarius). 5S rRNA is in its own clan (CL00113), so it never competes with LSU.

**Every other truth hit kept: all but two, by the rule's design.** Rfam puts Bacteria_large_SRP and
Bacteria_small_SRP in one clan (CL00003) and GlmZ_SraJ and GlmY_tke1 in another (CL00011); on these genomes each pair
hits one locus on one strand, so the lower-scoring model loses, as clan competition intends (the design note §5
expected them kept only "unless Rfam puts them in one clan with their neighbour"; it does). The four tRNA plus RtT
loci of E. coli are kept (RtT is in no clan), and so are the 12 antisense overlaps (opposite strands never compete),
including SdsR_RyeB across SraC_RyeA, which share clan CL00105. The tRNA and tRNA-Sec models never overlapped on
these genomes.

## What the binaries offer for clans (Infernal 1.1.5, help text in the image)

- `cmsearch -h`: no clan option.
- `cmscan -h`: `--clanin <f>  : read clan information from file <f>`; `--oclan : w/'--fmt 2' and '--tblout', only
  mark overlaps within clans`; `--oskip : w/'--fmt 2' and '--tblout', do not output lower scoring overlaps`; `--fmt
  <n> : set hit table format to <n> (1<=n<=3)`.
- `cmscan --clanin` refuses a clan file that names a model its database lacks (`Error: unable to find model
  alpha_tmRNA listed in ... in SSI index file`), so a test that gives it a subset of Rfam cuts the clan file to the
  subset's models.
- `cmpress` and `cmscan -T` accept uncalibrated models, so single-sequence models built with `cmbuild --enone` from
  segments of a random sequence give hits exactly where the segments were cut: constructed cases in seconds
  (`cmcalibrate` of 11 such models ran for over 50 minutes and was stopped; it is not needed).
- In `--fmt 2` output, column `olp` (with `--oclan`, only hits of the same clan count) is one of four marks, read
  from the binary's output on constructed cases (`olp_probe.py`, `olp_probe3.py` below):
  - `*`: overlaps no hit;
  - `^`: overlaps hits, none better;
  - `=`: overlaps a better hit marked `^` (`anyidx` names it); `--oskip` drops exactly these;
  - `$`: overlaps better hits, none marked `^` (each of them overlaps a better one itself); `--oskip` keeps it.

What the probes show (hits as cmscan printed them, 1-based inclusive):

| case | cmscan | meaning |
|---|---|---|
| D 721-770 (40.7 bits), E 770-799 (21.0), one clan | D `^`, E `=` (afrct1 0.033: 1 of 30 bases) | one shared base is an overlap |
| K and L, the same segment, both 30.9 bits | K `^`, L `=` | a tie goes to the first in the database, which is also the smaller name here; cmscan's tie order is not otherwise established |
| I 1451-1500 (+), J 1461-1490 (-), one clan | both `*` | opposite strands never compete |
| F 1101-1150 (clan), G 1111-1140 (no clan) | both `*` | a model in no clan never competes |
| chain of three: A 301-360 (51.1), B 341-380 (31.3), C 371-400 (20.5); A and C apart | A `^`, B `=`, C `$`; `--oskip` keeps A, C | a hit overlapped only by a loser stays |
| chain of four: P1 301-370 (61.3), P2 361-420 (50.6), P3 411-460 (40.1), P4 451-490 (31.3), each overlapping the next by 10 | P1 `^`, P2 `=`, P3 `$`, P4 `$`; `--oskip` keeps P1, P3, P4 | P3 and P4 overlap and both stay |
| chain of five: R1 > R2 > ... > R5 likewise | R1 `^`, R2 `=`, R3, R4, R5 `$`; `--oskip` keeps R1, R3, R4, R5 | the same, one level deeper |

So Rfam's own rule (`cmscan --oclan --oskip`) is one level deep: a hit is dropped when it overlaps a better hit of its
clan that overlaps no better hit; past the first loser in a chain, overlapping hits all stay. The producer's rule
(`ncrna.compete`) is greedy: in rank order (bit score, then accession, then coordinates), a hit is kept when it
shares no base with a hit of its clan and strand kept before it. The overlap is the binary's (one shared base, same
strand, same clan). Every hit cmscan drops the producer drops (a `^` hit is first in its neighbourhood, so it is kept
by both), every `*` and `^` hit it keeps the producer keeps, and the chain of three gives the same answer. They
differ only from a chain of four on: cmscan keeps P4 beside the kept P3, the producer drops it, so no two kept hits of
a clan ever overlap on a strand (the design note's §5 test assertion). The design note's original wording ("among
overlapping hits of models in one clan, keep the best") recalled `cmscan --oclan` as the reference; the binary shows
that neither reading of it (any better overlap loses; greedy) is what it does, and greedy is the one of the two that
agrees with it on every case above but the chains of four and five.

On the G-04b panel (`tests/r3/test_g04b_operon.py`: E. coli rrnB and B. subtilis rrnO operons, an M. genitalium tRNA
cluster with a tRNA-Sec hit on a tRNA, GlmZ/GlmY, the two SRPs, SdsR/SraC, tRNA/HPnc0260, tRNA/RtT, cut from the
golden genomes into one 17,050-base contig), `cmscan` and `cmsearch` report the same 33 hits, cmscan marks no hit `$`
(no chain), and it marks `=` exactly the 7 hits the producer flags as losers (both SSU_rRNA_archaea, both
LSU_rRNA_archaea, tRNA-Sec, Bacteria_small_SRP, GlmY_tke1). `tests/r2/test_ncrna_clans.py` runs the constructed cases
above through `cmscan` in every CI run and checks the producer against its marks, the chain of four included. The
8 genomes have no chain either: every competing group there is a pair (the compare script's groups), so the two
rules agree on them, and a merge rerun with the final (greedy) code into a copy of the root gave files byte-identical
to the pod run's (sha256 of every feature, `l1_ncrna` and ledger file).

## Commands

On the laptop: the worktree (without `.git` and `planning/`) was copied to `/tmp/if2453/e08/tree` with `rsync -a
--delete`; `/tmp/if2453/m09c/root` was copied to `/tmp/if2453/e08/root` (`cp -a`) so the steps could write. The
profile `/tmp/if2453/e08/e08-profile.yaml` is card E-07's with `tmp` under `/tmp/if2453/e08`. Before each launch the
account's `ps -u if2453 -o pcpu=` sum was 940% to 987% (night: the cap is 5600%). Each step ran in its own unit,
`MemoryMax=16G`, the two shards at once with two `cmsearch` threads each:

```bash
P=infernal-1.1.5-rfam-15.1-prok1-igr200nt-clan1
RUN="bash /tmp/if2453/e08/tree/scripts/pod_run.sh --profile /tmp/if2453/e08/e08-profile.yaml --image /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif --bind /tmp/if2453 --env PYTHONPATH=/tmp/if2453/e08/tree/src:/tmp/if2453/pilot_s1/pydeps --mem 16G"
STEP="run-step --profile /tmp/if2453/e08/e08-profile.yaml --root /tmp/if2453/e08/root --stage M9C"
NC="python -m biobox.cli stage ncrna --root /tmp/if2453/e08/root --stage M9C --cm /tmp/if2453/golden/rfam/Rfam.cm.gz"
$RUN --threads 2 e08-plan $STEP --step ncrna-$P.plan --threads 1 -- $NC plan --scan-bases 9500000 --tmp /tmp/if2453/e08/tmp/plan
$RUN --threads 2 e08-shard-0 $STEP --step ncrna-$P.shard-00000 --threads 2 -- $NC shard --shard 0 --cpus 2 --tmp /tmp/if2453/e08/tmp/scan --seq-stages M9C
$RUN --threads 2 e08-shard-1 $STEP --step ncrna-$P.shard-00001 --threads 2 -- $NC shard --shard 1 --cpus 2 --tmp /tmp/if2453/e08/tmp/scan --seq-stages M9C
$RUN --threads 1 e08-merge $STEP --step ncrna-$P --threads 1 -- $NC merge --tmp /tmp/if2453/e08/tmp/scan
```

The comparison and the derived numbers, in the image in a `MemoryMax=4G`, `CPUQuota=100%` scope at nice 19:
`python -I e08_compare.py /tmp/if2453/e08/root M9C $P prok-1.tsv Rfam-15.1.clanin /tmp/if2453/m09c/out
/tmp/if2453/m09c/fa/wg_all.map.json /tmp/if2453/e08/compare.json` and `python -I e08_numbers.py /tmp/if2453/e07/root
/tmp/if2453/e08/root` (both scripts in `/tmp/if2453/e08/`). The help texts: `cmscan -h` and `cmsearch -h` in the image.
The probes, in the image in `MemoryMax=2G` scopes at nice 19: `olp_probe.py` built the models of the first table's
first five rows (`cmbuild -F --enone -n <name>`, a Stockholm file of one sequence with an all-`.` `SS_cons`) and
`cmpress`ed them in `/tmp/if2453/e08/probe2` (its `cmcalibrate` step was stopped), then `cmscan --cpu 2 -T 5 --incT 5
--notrunc --fmt 2 --clanin clanin --oclan [--oskip] --noali --tblout ... db.cm contig.fa`; `python -I olp_probe3.py
/tmp/if2453/e08/probe3` did the same for the chains of four and five. The pod run's steps ran with the branch's first
rule (a hit loses to any better overlapping hit); after the switch to the greedy rule, `ncrna.merge` was run again on a
copy of the root (`/tmp/if2453/e08/root2`, the producer's outputs removed, the shards' scratch runs reused) and its
files compared by sha256 with the published ones: identical. Outputs stay on the pod under `/tmp/if2453/e08`; nothing
was written under `/stor` or the m09c, e07, r10 and r50 directories.

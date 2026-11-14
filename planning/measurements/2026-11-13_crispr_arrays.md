# CRISPR arrays by MinCED on stages P1, P2 and S1 (card E-12)

Date: 2026-11-13, wilkcomp02. Every number carries a label: *measured* (host, date, command or record), *derived*
(formula and inputs) or *estimate*.

## Setup

- Code: branch `task/E-12`, synced to `/tmp/if2453/ci/E-12` by `scripts/pod_pytest.sh E-12 ...`.
- Image: the runtime image built from `main` at `faec1f4`, `/tmp/if2453/images/bio-box_faec1f4.sif`,
  sha256 `d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a` (`sha256sum` on the pod, the `--digest`
  of every launch and the `image_digest` of every run record). MinCED 0.4.2, OpenJDK 25.0.2-internal (`java -version`
  in the image).
- Data: read-only. The data root `/tmp/if2453/r10/root` was bound read-only into every unit
  (`--bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro`). The scratch root `/tmp/if2453/e12/view` was made by
  `/tmp/if2453/e12/bin/mkview_e12.py /tmp/if2453/r10/root /tmp/if2453/e12/view S1,P1,P2` (the pattern of card B-42):
  real directories holding one symlink per stage directory of `contig_index`, `seq_locator`, `seq_shards` and
  `genome_input` of S1, P1 and P2, and `store_params.json`; no run record, ledger, `l1` file or `STOP` file. The
  producer's files, its ledger and run-step's record went into the view's own directories. Nothing was written under
  the real root.
- Profile: `profiles/wilkcomp02.yaml` with `tmp: /tmp/if2453/e12/tmp` (`/tmp/if2453/e12/profile.yaml`).
- Launch, per stage S (units `e12-p1`, `e12-p2`, `e12-s1`; the other two stages after `--seq-stages`):

```text
bash /tmp/if2453/ci/E-12/scripts/pod_run.sh --profile /tmp/if2453/e12/profile.yaml \
  --image /tmp/if2453/images/bio-box_faec1f4.sif --digest sha256:d67300aa...f098a \
  --bind /tmp/if2453 --bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro --env PYTHONPATH=/tmp/if2453/ci/E-12/src \
  --mem 16G --threads 4 e12-p1 run-step --profile /tmp/if2453/e12/profile.yaml --root /tmp/if2453/e12/view \
  --stage P1 --step crispr-minced-0.4.2-w8-nr3-r23x47-s26x50 --threads 4 -- \
  python -m biobox.cli stage crispr --root /tmp/if2453/e12/view --stage P1 --cpus 4 \
  --tmp /tmp/if2453/e12/scratch-p1 --seq-stages P1,S1,P2
```

  4 MinCED processes at once, each a JVM with a fixed 512 MiB heap, nice 19, `MemoryMax=16G`. The account ran a
  40-thread Pfam scan beside them (`top`, about 41 to 46 CPUs before each launch, under the night cap of 56).

## Results

Run records `runs/stage=<S>/crispr-minced-0.4.2-w8-nr3-r23x47-s26x50.json` of the view and the step's JSON result
(*measured*, wilkcomp02, 2026-11-13). Peak RSS is the record's `peak_rss_bytes` (`child_maxrss`: the largest process
of the step, the producer or one JVM); the unit's peak is `pod_run`'s `MemoryPeak` (all four JVMs and the producer).

| | P1 | P2 | S1 |
|---|---|---|---|
| genomes | 250 | 2,247 | 27,258 |
| contigs searched | 539 | 5,018 | 27,258 |
| genome bases | 1,016,625,803 | 9,596,180,328 | 1,664,109,938 |
| arrays | 286 | 2,291 | 213 |
| repeats | 9,287 | 63,317 | 1,365 |
| spacers | 9,001 | 61,026 | 1,152 |
| distinct spacer sequences | 8,638 | 59,424 | 339 |
| genomes with an array | 106 | 915 | 184 |
| arrays per genome (*derived*: arrays / genomes) | 1.14 | 1.02 | 0.0078 |
| repeats per array, median / max | 19 / 357 | 16 / 400 | 7 / 36 |
| arrays of 3 or 4 repeats | 26 | 214 | 43 |
| median repeat / spacer length (bases) | 32 / 35 | 30 / 35 | 30 / 36 |
| longest array (bases) | 26,487 | 28,431 | 2,297 |
| arrays across a circular origin | 0 | 0 | 0 |
| arrays rejected (longer than the circle) | 0 | 0 | 0 |
| skipped contigs (another stage's) | 0 | 0 | 0 |
| CPU-minutes (record `cpu_min`) | 3.083 | 26.623 | 5.096 |
| wall seconds (record `wall_s`) | 50.6 | 405.3 | 87.5 |
| peak RSS of one process (MB) | 337.6 | 435.5 | 531.3 |
| unit MemoryPeak (GB) | 1.77 | 1.85 | 1.62 |
| CPU-seconds per Gbp (*derived*: cpu_min × 60 / Gbp) | 182.0 | 166.5 | 183.7 |
| bases per CPU-second (*derived*: bases / (cpu_min × 60)) | 5.50 M | 6.01 M | 5.44 M |

The table's distributions are from `/tmp/if2453/e12/bin/stats.py` over the view's files (*measured*, same date).

The CPU budget (*derived*). The step's planning rate is P1's, 182.0 CPU-seconds per Gbp (`CRISPR_CPU_S_PER_GBP` in
`src/biobox/stage.py`); S1 ran 0.9% slower per genome base (its circular phage records carry two flanks of their own
length), P2 8.5% faster. At 182.0 a stage reaches 80% of the 300 CPU-minute budget at 79.1 Gbp (0.8 × 300 × 60 /
182.0). The largest stage built so far, P2, used 26.6 CPU-minutes, 8.9% of the budget, so the step has no shards; above
79.1 Gbp `build-stage` refuses it (exit 2) and the stage must be split (or the step sharded, a later card). The
earlier probe of one 8-genome batch (11.8 CPU-s for about 0.1 Gbp, *measured*, same host and date, `/usr/bin/time`)
gave about 118 CPU-s per Gbp without the producer's reading of the sequence store; the producer's own process reads
and decompresses every contig, which is the difference (*inferred*, not profiled).

Memory (*measured*, as above). One MinCED JVM on the largest P2 genome (GCA_037157765.1, 13.4 Mbp in one record)
peaked at 147 MB with a 256 MiB heap and 351 MB with 1 GiB, with byte-identical GFF (`/usr/bin/time -f %M`). The
test of the 10x and 100x synthetic stores (`tests/r2/test_crispr_memory.py`) measured 339.6 MB and 339.5 MB (bound
768 MiB, growth 32 MiB).

## The worked example's control: *Photobacterium damselae* 04Ya311 (GCF_026001825.1, stage P2)

`/tmp/if2453/e12/bin/control.py` (in the image, read-only): the producer's P2 arrays on the genome's three records,
and the stage's gene calls (`pyrodigal-3.7.1-gv-0.3.2-p2`) with their Pfam 38.2 domains from the real root. Offsets
count genes from the array: -1 is the gene just before it, +1 the gene just after (*measured*, wilkcomp02,
2026-11-13).

- **Chromosome 2 (`NZ_AP026781.1`), the Tn7-like element of the control** (TniQ 426340-427525, Cas7f 423377-424433,
  Cas5f-8f 424438-426340, Cas6f 422775-423375, all on the minus strand; TnsA-TnsB at 440353-442837): **no array within
  10 genes** at MinCED's default `-minNR 3`. The nearest arrays are those of the complete system below, 209 genes
  away. MinCED with `-minNR 2` on the same genome (a probe, not the producer; `minced -minNR 2 -gffFull`) reports a
  two-repeat array at 422558-422647 (1-based; repeat `TTATTATCAGCTACCTATTCGGCAGTTCAC`, one spacer), 128 bases from
  Cas6f's 3' end: the paper's "short CRISPR array" of a minimal I-F3 system is here a single spacer, which the
  producer's three-repeat minimum does not call.
- **Chromosome 2, the complete type I-F system** (Cas1 681243-682230, Cas2-Cas3 682226-685571, Cas8f 685586-686888,
  Cas5f 686887-687844, Cas7f 687857-688841, Cas6f 688842-689397, all on the plus strand): two arrays, one on each
  side. [679543, 681134): 27 repeats, 26 spacers, ending 109 bases before Cas1 (Cas1 at offset +2: one called ORF lies
  inside the array). [689534, 690762): 21 repeats, 20 spacers, starting 137 bases after Cas6f (offset -1). Repeat
  consensus of both `CTTCACTGCCGAGTAGGCAGCTTAGAAA` (MinCED's `rpt_unit_seq`).
- **Chromosome 1 (`NZ_AP026780.1`), the second Tn7-like I-F3 element**: one array, [2224313, 2224582), 5 repeats and
  4 spacers, 128 bases from Cas6f (offset +1), then Cas7f (+2), Cas5f-8f (+3, 2,088 bases), TniQ (+4); TnsB and TnsA
  are 13 and 14 genes away. Its repeat (consensus `ATTTTTAGCTACCTATGCGGCAGGACAC`) differs from the complete system's
  and resembles the chromosome 2 element's two-repeat array (`...AGCTACCTAT.CGGCAG...`).
- The plasmid (`NZ_AP026782.1`): no array.

No array check was added to the example (`examples/**` is not this card's), and no release holds arrays.

## Later cards (what this does not do)

- Orientation (CRISPRDirection, `PLAN.md` step 4): every array is stored on strand +1, the strand its sequences read
  on. An orientation would be a new nullable column of `l1_crispr` (LY-09), so `feature_uid` does not change; it needs
  CRISPRDirection (or an equivalent rule from the repeat's conserved 3' end and AT-rich leader) pinned in the image,
  and a test against arrays of known orientation (the golden genomes' Cas operons).
- The spacer corpus of the search pack (spec invariant 16, `planning/data_layer_and_tests.md` §2B.2): spacers keyed by
  `region_uid` (stored in `l1_crispr`), with their array (`feature_uid`), repeat and contig; the pack builder
  (`biobox.search.pack`) would write `spacers/` from `l1_crispr` rows of element `spacer` (both strands searchable),
  remove it from `missing_targets`, and C-12 would count `spacers` = spacer calls. Distinct spacers are many fewer than
  spacer calls in phage collections (S1: 339 of 1,152).
- Tape membership: arrays are not tape slots yet. `stage tapes --producers <gene caller>,<this producer>` builds a
  tape with them (tested), but the C-08 check (`biobox.provenance._producer_problems`) knows only gene-call and ncRNA
  producer names and would reject the tape version; `STEP_TOOLS` has no `crispr-*` entry. An array slot also shifts
  `nbr` gene windows (`--window Ng` would count it as a slot).
- Two-repeat arrays (`-minNR 2`): the control's minimal I-F3 array has one spacer. Lowering the minimum is a new
  producer and needs its false-positive rate measured (for example against CRISPRCasFinder's evidence levels).

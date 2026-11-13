# T7-like RNA polymerases in phages and prokaryotes: the guide's worked example as a release check

The first worked example of `GUIDE.md` (sections 3 to 7 teach its commands one at a time), in the form of the
CRISPR-Cas example (`examples/crispr_cas_neighbourhoods/`): one script that runs it on the pod, the facts its answer
must show on any release (`invariant`), and its counts per release (`recorded/`). Card E-11; rerun on every milestone
release (AGENTS §5.7, `examples/README.md`). State of 2026-11-13, release `rN.3`. Numbers carry a label: *measured*
(host, date, command), *derived* (formula and inputs) or *estimate*.

## 1. The plan (written before the run)

**The question.** Where are the single-subunit, T7-like RNA polymerases in a release, how many are in phages and how
many in prokaryotes, and which genes sit next to them? It is the guide's example because its answer is known twice
over: T7's genome map puts gene 0.7 (the protein kinase) just upstream of gene 1 (the polymerase) and gene 1.3 (the DNA
ligase) three genes downstream, and two trials answered it on `rN.1` by hand and with `search hmm`
(`planning/measurements/2026-11-10_trial_search_rN1.md`, `2026-11-11_trial_search_rN1_second.md`).

**The model and cutoff.** Pfam `RNA_pol` (PF00940.25, 449 match states), fetched by name from the pod's Pfam 38.2
file, the model of both trials (sha256 `32e6f0f0...0e87e`, *measured*: `sha256sum`, wilkcomp02, 2026-11-13). Searched
at 30 bits on the sequence score, the guide's cutoff, above the model's gathering cutoff (GA 22.00, the model's `GA`
line), so the release's own Pfam scan (at GA) holds proteins the search does not, and the comparison with it is made
at 30 bits.

**The control.** T7 gene 1, GenBank `NC_001604` 3171..5822 (`NC_001604:3170-5822` in bio-box's 0-based form), feature
`39f9ef263e234b012d2671634468546e`, protein `2eef240afefddcd2c8ac7660eb50e384` on contig
`719747a8fef97b2380955d92845861db` (*measured*, `bio-box find`, second trial file, section 1). Its window by
coordinates on `rN.3`, from `protein_occurrence` and the release's Pfam scan (*measured*, wilkcomp02, 2026-11-13,
`bio-box sql` before the run, the query of `sql/checks_coord.sql` widened to 8 genes): `ocr` (gene 0.3) at -5,
`Phage_T7_GP04`, `Gp0_5`, `Phage_T7_0_6`, `Phage_T7_kinase` (gene 0.7) at -1, gene 1 (`RPOL_N+RNA_pol`), then
`Phage_T7_1_1`, `T7-like_gp12`, `DNA_ligase_A_M+DNA_ligase_C` (gene 1.3) at +3, `DUF5464`, `BC10` at +5.

**What would count as agreement.** In kind, on any release that holds T7: the control found at or above 30 bits with
that architecture, its kinase at -1 and its ligase at +3, the run's window equal to the contig's own gene order, the
shuffled-model control empty, hits in phages and in prokaryotes, and `hmmsearch` run by hand agreeing with the search.
In number, on `rN.3`: the trials' values on `rN.1` (the same genomes): 1,957 proteins, 1,931 of them in phages and 26
in prokaryotes, on 2,581 genes (second trial, section 3); 2,581 loci on 2,526 contigs, 463 cut by a contig end, 24,473
neighbours (same); 432 flank groups, 1,558 leaves, 598 context types (`planning/measurements/2026-11-12_flanks_rN3.md`);
2,569 hit genes with a Pfam architecture and `Phage_T7_1_1` next to 1,032 loci (`2026-11-11_build_r50.md` §A3).

**What bio-box cannot do for this question today.** No taxonomy (`contig_dim.tax_order` is null), so "prokaryote" is
the genome kind of the stage's inputs, not a lineage; no enrichment statistic (counts are counts; cards F-01 to F-03);
no ncRNA in `rN.3`; neighbours are genes in a window, not operons.

## 2. Running it

On the pod host, in `bash`, as `GUIDE.md` section 2 sets it up; the laptop runs nothing. The environment is every
example's (`examples/README.md`):

```bash
export TREE=/tmp/if2453/e11/tree
export PROFILE="$TREE/profiles/wilkcomp02.yaml"
export IMAGE=/tmp/if2453/images/bio-box_64f1bf0.sif
export DIGEST=sha256:cb7353f85316981be5268c2212100322b5fbf1086c32b07e2f9b47802b6afba3
export ROOT=/tmp/if2453/r10/root
export RELEASE=rN.3
export OUT=/tmp/if2453/<you>/rnap-rN.3
export EX="$TREE/examples/rnap"
bash "$EX/run.sh"
tail -n 1 "$OUT/check.txt"
```

The image is the runtime image built from `main` (RUNBOOK section 0), which holds matplotlib 3.11.2 for `--figure`
(*measured*: `import matplotlib` in it, wilkcomp02, 2026-11-13). The steps, as `run.sh` runs them (`bb` launches one
`bio-box` command through `scripts/pod_run.sh`, `tool` one HMMER or Python command in a capped scope; both wait for
room under the account's thread cap: `examples/_shared/launch.sh`; `$S`, `$F`, `$Q` are the shared flags written out
in `run.sh`):

1. The model: `hmmfetch -o "$OUT/hmm/RNA_pol.hmm" "$PFAM" RNA_pol`, from the Pfam file on the pod (read-only).
2. The control must name one gene: `bio-box find --root "$ROOT" --release "$RELEASE" --accession NC_001604:3170-5822`.
3. The search, as `GUIDE.md` section 4: `bio-box search hmm --root "$ROOT" --release "$RELEASE" --profile "$PROFILE"
   --runs "$OUT/runs" --work "$OUT/work" --neighbors 5 --cpus 2 --mem-gb 14 --threads 2 --hmm "$OUT/hmm/RNA_pol.hmm"
   --bits 30 --positive NC_001604:3170-5822 --out "$OUT/rnap"`; the run id is read from its log.
4. Flanks, as `GUIDE.md` section 6, with figures: `bio-box nbr flanks "$RUN" --root "$ROOT" --profile "$PROFILE" --runs
   "$OUT/runs" --threads 2 --window 5g --level nr50 --figure --out "$OUT/flanks"`.
5. The SQL checks with the run attached, `bio-box sql --root "$ROOT" --release "$RELEASE" --runs "$OUT/runs" --threads
   2 "$EX/sql/checks.sql" --out "$OUT/observed_sql.tsv"`, and without it, by coordinates, `bio-box sql --root "$ROOT"
   --release "$RELEASE" --threads 2 "$EX/sql/checks_coord.sql" --out "$OUT/observed_coord.tsv"`.
6. The oracle: `sql/oracle_uids.sql` (every hit and every protein the Pfam scan gave `RNA_pol`), `bio-box seq
   --proteins-file`, `hmmsearch --cpu 2 -Z 1e8 --domZ 1e8 -T 0 --noali --tblout "$OUT/oracle.tbl" ...`, then
   `python -I "$TREE/examples/_shared/oracle.py" --hmm ... --tblout ... --prefix rnap --model RNA_pol --bits 30
   "$OUT/rnap/hits.tsv"`.
7. The values read from the run's files: `python -I "$EX/observe.py"` on `hits.tsv`, `neighbors.tsv` and the flanks.
8. `observed.tsv` (the four files' rows), then `python -I "$TREE/examples/_shared/check.py" --expected
   "$EX/expected.tsv" --recorded "$EX/recorded/$RELEASE.tsv" --write-recorded "$OUT/recorded_new.tsv"
   "$OUT/observed.tsv"`.

## 3. What it found on `rN.3`

All *measured* on wilkcomp02, 2026-11-13, by `run.sh` (first run `OUT=/tmp/if2453/e11/rnap1`, 13:47 to 14:05 CST,
run `45dd1fdcbf56eedaf627c0dd37d2a7a9`, pack `a334c89ccf048904f297de3b2597af22`; the code of commit `8ea227f` before
its `recorded/` and `results/` files existed, runtime image `bio-box_64f1bf0.sif`), unless another source is named. `recorded/rN.3.tsv` is that run's `recorded_new.tsv`;
`results/` holds the release check's run (`scripts/release_examples.sh rN.3`, `OUT=/tmp/if2453/e11/acc/rnap`),
whose values equal it (section 5).

| | `rN.3` (this example) | Measurement notes | Same? |
|---|---|---|---|
| proteins at 30 bits | 1,957 | 1,957 (second trial on `rN.1`; §A3 on `rN.3`) | yes |
| in phages / in prokaryotes (proteins) | 1,931 / 26 | 1,931 / 26 (second trial, section 3) | yes |
| hit genes (loci) | 2,581: 2,555 on phage contigs, 26 on prokaryote or plasmid contigs | 2,581; 2,555 and 26 (same) | yes |
| contigs, genomes | 2,526, 2,526 | 2,526 contigs and genomes (same) | yes |
| T7 gene 1 | 367.7 bits, `RPOL_N+RNA_pol`, recovered | 367.7 bits (second trial; §A3) | yes |
| negative control (shuffled model) | 0 hits | 0 hits (second trial; §A3) | yes |
| `hmmsearch` by hand | 0 hits only in the search, 0 only in `hmmsearch`, largest score difference 0.00 bit, over the hits and the Pfam-annotated proteins | 0 and 0, 0.000 bit, over every protein of `rN.1` (second trial, section 3) | yes; this oracle reads fewer proteins (below) |
| hit genes with a Pfam architecture | 2,569 | 2,569 (§A3) | yes |
| Pfam scan: proteins with `RNA_pol` (GA 22) | 1,957; 11 of them below 30 bits; 11 hits have no `RNA_pol` row | not in the notes | new |
| `neighbors.tsv` header | 2,581 loci, 2,526 contigs, 463 censored, 24,473 neighbours | the same four (second trial; §A3) | yes |
| top Pfam neighbour | `Phage_T7_1_1`, 1,032 loci | `Phage_T7_1_1` (PF08200.18), 1,032 loci (§A3) | yes |
| `nbr flanks` (nr50, 5g, anchor NR90 groups, cut 0.5) | 432 groups, 1,558 leaves, 598 context types, 0 types with phage and non-phage loci | 432, 1,558, 598, 0 (flanks note) | yes |
| T7 gene 1's flank group | `c1bf639f...`, 162 loci, 64 leaves, 10 context types; its type holds phage loci only | the note's "T7-like group" is another NR90 cluster, `4a0062a1...` (306 loci, 202 leaves, 11 types); the 162-locus, 10-type group is in its list of largest groups | yes (different groups, both present: 306 and 11 also here, `groups.tsv`) |

So every number of the notes is reproduced on `rN.3` with the current code, none differs. The search hits and the
Pfam scan's `RNA_pol` set are the same size and not the same set: 11 proteins score 30 bits or more with no domain row
at GA (the sequence passes, no single domain does), and 11 others have a domain row at GA with a sequence score under
30. Every protein the scan gives `RNA_pol` with 30 bits or more is a hit, and the scores of the 1,946 common proteins
are equal (largest difference 0.0 bit; *derived*: 1,957 - 11).

**The control's window**, in the run (`neighborhoods`, 5 genes on each side) and by coordinates alone: the same 10
genes at the same offsets, gene 0.7's kinase at -1 and gene 1.3's ligase at +3, as in section 1.

**Its flank group.** In T7 gene 1's group (`c1bf639f92687b483b659b4dadf7c52c`, the anchor's NR90 cluster) the first
context type holds 132 of the 162 loci. Its conserved neighbours (`conservation.tsv`, NR50 clusters along the anchor's
strand): `916b7b66...` at +1 in 129 of 131 loci that have it, `eb4f50cb...` at +2 in 128 of 130, `fa4cb451...` (the
DNA ligase's cluster, `DNA_ligase_A_M` by its genes' Pfam names) at +3 in 126 of 131, and `5fe1410c...` (the
kinase's, `Phage_T7_kinase`) at -1 in 113 of 113. The figure `figures/t7_group_c1bf639f.svg` (and `.pdf`) draws it:
64 leaves in tree order beside the UPGMA tree, genes coloured by NR50 family, grey below 10% of the group's loci.

**The oracle reads fewer proteins than the trials.** `hmmsearch` here scores the hits and the proteins the release's
Pfam scan gives `RNA_pol` (1,968 proteins, *derived*: 1,957 + 11), not all 10,620,684 proteins of the release as the
trials did (88.2 s and 3.7 GB of export, second trial). A protein that scores 30 bits but is neither a hit nor in the
Pfam scan would escape it; the trials found none on `rN.1`, and the Pfam scan holds every protein with a domain at GA
22, eight bits under the cutoff.

## 4. Cost on `rN.3`

*Measured*, wilkcomp02, 2026-11-13, the first run (`OUT=/tmp/if2453/e11/rnap1`), 2 threads, the account at 38 to 41 CPUs
of its day cap of 40 for most of it (the `waiting for room` lines):

| Step | Wall | CPU | Memory |
|---|---|---|---|
| `hmmfetch` of `RNA_pol` from `Pfam-A.hmm.gz` | 42.9 s | 67.3 CPU-s | 45 MB |
| `find` | | 5.8 CPU-s | 0.57 GB (unit peak) |
| `search hmm` | 453 s for the unit (13:48:55 to 13:56:28); 394.1 s the run (`wall_s`) | 675.2 CPU-s for the unit (11.25 CPU-min); worker 9.33 CPU-min | worker 2.62 GB (`peak_rss_bytes`); unit 8.89 GB with page cache |
| `nbr flanks --figure` (432 figures) | 268 s | 279.4 CPU-s (4.66 CPU-min) | 4.03 GB (unit peak) |
| three `sql`, one `seq` | | 5.8 to 16.0 CPU-s each, 37.0 in all | 0.44 to 1.16 GB (unit peak) |
| `hmmsearch` oracle on 1,968 proteins | 15.1 s | 28.6 CPU-s | 46 MB |
| `oracle.py`, `observe.py`, `check.py` | 1.5 s | 1.5 CPU-s | 47 MB |

In all 1,094.8 CPU-seconds, 18.2 CPU-minutes (*derived*: the units' `CPUUsageNSec` from the `pod_run: ... usage`
lines, 997.4 s, plus `/usr/bin/time -v` user and system time of the scopes, 97.4 s); 18 min 10 s from start to end, of
which at least 4 minutes waiting for room (*derived*: 4 `waiting for room` lines of 60 s or more). Against the notes:
the search's unit CPU is 11.25 CPU-minutes here and 10.78 in §A3 (+4%), the worker's 9.33 against 9.18 (+2%) and its
peak memory 2.62 GB against 2.69 GB (-2%); the run's wall time is 394.1 s against 280.5 s (+40%) because this run gave
the search 2 CPUs and §A3 gave it 4 (the same CPU work on half the CPUs; not a like-for-like number). `nbr flanks
--figure` took 4.66 CPU-minutes against 4.57 (W-27, flanks note; +2%).

## 5. Running this as a release check

`scripts/release_examples.sh RELEASE` runs this example and every other one (`examples/README.md`); alone, it is
section 2's command with `RELEASE` and `OUT` changed. `run.sh` exits with `check.py`'s code: 0 when every invariant
holds, 1 when one fails or is missing, 2 on a malformed input or a step that failed. The release needs a search pack in
the profile's pack directory (`bio-box releases` names it), and the control genome, T7 (`NC_001604`, stage S1), which
every `rN.*` and `rM.*` release holds.

The release check of `rN.3` (*measured*, wilkcomp02, 2026-11-13, `scripts/release_examples.sh rN.3` with
`OUT=/tmp/if2453/e11/acc`, 14:37 to 15:03 for this example, 18.4 CPU-minutes): every invariant of this example held
and every recorded value equalled `recorded/rN.3.tsv` (`results/check.txt`: 42 ok); its outputs were byte-identical
to the first run's but for the run id (`planning/measurements/2026-11-13_release_examples_rN3.md`).

**What a failed invariant would mean.**

| Check | If it fails |
|---|---|
| `t7_control_bits`, `t7_positive_recovered` | The model no longer scores T7 gene 1 at 30 bits: the model file (`rnap_hmm_sha256` differs: another Pfam version), the pack or the scoring changed. (A lost control already stops the search, `positive control not recovered`.) |
| `t7_control_architecture` | The release's Pfam scan of the control changed: another library, or the scan's cutoffs |
| `t7_control_order` | Gene calls, the window's orientation (`offset_rel`) or the Pfam names of T7's genes changed. If `t7_window_is_the_contig_order` holds, it is the annotation, not the window |
| `t7_window_is_the_contig_order` | The run's window of the control is not the contig's gene order by coordinates: a defect of `nbr` or of the tapes. Treat every window-based check as unproven until it is explained |
| `coord_t7_kinase_offset`, `coord_t7_ligase_offset` | By coordinates and the Pfam scan alone, without the run: T7's gene calls or their domains changed. If these hold and `t7_control_order` fails, suspect the run's windows |
| `rnap_negative_control_hits` | The shuffled model finds proteins: a scoring or composition-bias change (the run itself fails above its allowance) |
| `rnap_loci_virus`, `rnap_loci_not_virus` | No hit in phages, or none in prokaryotes: the genome kinds of the stage inputs (`genome_input.kind`) were lost, or a release without one of the two |
| `rnap_loci_missing_windows` | A gene of a hit protein has no window in the run: `nbr`'s anchors and `protein_occurrence` disagree |
| `rnap_hits_below_cutoff`, `rnap_pfam_proteins_at_cutoff_not_hit`, `max_bits_diff_search_vs_pfam` | The search and the release's own Pfam scan disagree: a pack that misses proteins, or a scoring difference |
| `oracle_rnap_*` | The search and `hmmsearch` run by hand disagree on the hits or their scores |
| `rnap_flank_loci_missing` | `nbr flanks` dropped a hit locus from `context_types.tsv` |

## 6. Files here

| File | What it is |
|---|---|
| `README.md` | this example |
| `run.sh` | the whole example as one command (section 2); it sources `examples/_shared/launch.sh` |
| `expected.tsv` | the checks: id, assertion, kind, expectation |
| `recorded/rN.3.tsv` | the recorded values of `rN.3` (the first run's `recorded_new.tsv`) |
| `observe.py` | the values read from the run's files (`neighbors.tsv`, `nbr flanks`, the model's sha256); standard library only; tested by `tests/r1/test_examples_check.py` |
| `sql/checks.sql`, `sql/checks_coord.sql`, `sql/oracle_uids.sql` | what `run.sh` runs: with the run, by coordinates without it, and the oracle's proteins |
| `results/` | the release check's run: `check.txt`, `observed.tsv`, `SHA256SUMS` of what stays on the pod |
| `figures/t7_group_c1bf639f.svg`, `.pdf` | the `nbr flanks` figure of T7 gene 1's group |

The comparison (`check.py`) and the oracle (`oracle.py`) are shared by every example: `examples/_shared/`. Kept on the
pod, not here (local disk, wilkcomp02): `/tmp/if2453/e11/rnap1` and `/tmp/if2453/e11/acc/rnap`, each with the run
(`runs/`), `rnap/` (`hits.tsv`, `hits.faa`, `neighbors.tsv`, 2,581 locus tables) and `flanks/` (432 groups with their
figures).

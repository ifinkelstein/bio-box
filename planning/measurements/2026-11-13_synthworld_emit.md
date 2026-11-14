# synthworld S11 (emit) on `micro`, and the emitted world through the normal builders

Card F-02, session 8 (branch `task/F-02-8`). Code: `src/biobox/synthworld/emit.py` (its docstring states every
rule), preset section `emit` of `presets/micro.yaml`, driver `scripts/stress/f028_synthworld_builders.py`. Tests:
`tests/r1/test_synthworld_emit_units.py` (the contaminant window and the lift-and-place chain against a
provenance oracle), `tests/r2/test_synthworld_emit.py` (one fixture builds S1-S11). Small tables beside this note in
`2026-11-13_synthworld_emit/`.

Host and commands for every *measured* number: wilkcomp02, 2026-11-13, one process at nice 19 in a systemd scope,
the account at 33 to 37 threads of its 40-thread day cap (a 32-thread Pfam scan of another card beside it).
**Run 1** and **run 2**: `sh scripts/pod_pytest.sh F-02-8 -s tests/r2/test_synthworld_emit.py` (16:28 and 16:35;
default test image `bio-box-test_64f1bf0.sif`). **World**: unit `f028-world`, `python
scripts/stress/f028_synthworld_builders.py world /tmp/if2453/f028/world` in the runtime image `bio-box_64f1bf0.sif`
(16:33). `micro` default seed 20261112.

## Part A: S11 emit

### Time and memory

| quantity | run 1 | run 2 | world | label |
|---|---|---|---|---|
| S1-S2 backbone | 37.03 s | 36.38 s | 36.25 s | measured |
| S3-S4 sequences | 82.58 s | 79.51 s | 81.35 s | measured |
| S5 replicons | 43.41 s | 43.25 s | 43.98 s | measured |
| S6 mobile | 11.36 s | 11.13 s | 10.75 s | measured |
| S7 fusions | 5.01 s | 5.73 s | 5.52 s | measured |
| S8-S9 assembly | 18.61 s | 18.02 s | 17.83 s | measured |
| S10 sources | 6.63 s | 7.27 s | 6.94 s | measured |
| **S11 emit** | 27.04 s | 25.69 s | 26.78 s | measured |
| S11 in a child (`PYTHONHASHSEED=919`, one thread) | 25.55 s | 28.59 s | - | measured |
| S11 peak RSS in that child (VmHWM) | 167 MB | 165 MB | - | measured |
| S1-S11 in one process | 231.67 s | 226.98 s | 229.39 s | run 1 and 2 derived (sum of the stage times, one process); world measured (the driver's sum, same) |
| S1-S11 peak RSS of the whole process (`ru_maxrss`) | - | - | 415 MB | measured |

S11 is about 27 s, against the session-7 note's *estimate* of about 5 s ("one write of every source at S7's cost
per rewrite"): more than 20% off (AGENTS §6), reported here. The estimate priced the FASTA only; most of S11 is the
truth: every gene piece, junction, element, motif, spacer, fused ORF and call error lifted through three frames
and placed on every occurrence, and the hashing of every contig and IGR (*estimate*, not profiled). S1-S11 is
227 to 232 s, so the 3-minute target of spec §2.6 is missed by 47 to 52 s (*derived*: 226.98 - 180 and
231.67 - 180); the world was not thinned.

### Realized values (measured, runs 1 and 2 identical)

| quantity | value |
|---|---|
| digests | emit `a2f0e448484dab63d41281a5ef6a00c0`, world `f0262027e4ccdcbc72cd7d90f23fbb5c`, release_tag `sw-51319743b457c684`; identical in run 1, run 2, the child process and the world build in the runtime image |
| earlier stage digests after S11 | unchanged: backbone `025fdd0d27fdc65e678782216ecd21fb`, sequences `34ebce1324413a261acd874d7ee4fe3e`, replicons `55bda97155d75148e5cdca7b90a1d7be`, mobile `0bd7f0db5691de58ea5e1a0df1e6e570`, fusions `db09abf3fe79ae34a73e3ab0664a9b13`, assembly `2f425158dc38f2302c1d02e9ebed30da`, sources `770007bde8615f41d9dd5f2eb912e753`; every file of every earlier stage byte-identical before and after |
| sources | 62: 46 original, 2 exact and 2 near duplicates, 12 phages; 2,327 contig occurrences, 2,217 distinct contigs (exact duplicates share their uid; 2 near-duplicate contigs with no site to change keep the original's uid, as session 7 recorded) |
| build-stage genomes | 60 (the 2 exact duplicates are the same genome by `genome_uid` and are listed once); 2 contigs in two genomes (those 2 unchanged near-duplicate contigs, DEC-27); accepted by the builders' `--genomes`/`--objects` readers and `genome_input`'s refusals |
| labels | 1 of 50 host sources mislabelled (rate 0.05); 7 classes, one biome each |
| truth_gene rows | 86,194 whole and 3,812 cut on the assembly contigs (every emitted gene piece of `gene_contigs.tsv`), 2,876 whole and 155 cut on near-duplicate copies, 98 whole and 18 cut on contaminants, 730 phage genes; every whole CDS translates to its protein_uid under its genome's code (frameshifted genes excepted, whose uid is the intended protein), every row's bases are its S8 gene's (near copies: synonymous changes only) |
| truth_igr | 86,350 whole and 309 cut on assembly contigs, 2,731 and 13 on near copies, 106 and 2 on contaminants; not placed: 2,126 dropped, 52 touched by an S7 event, 401 by an S8 edit (none by S6, whose insertions sit at gene ends) |
| truth_rna | 1,058 whole on assembly contigs, 30 on near copies, 2 on contaminants, 20 dropped |
| truth_spacer_target | 3,711 rows with array and target whole, 20 decoy spacers (no target), 59 whose prophage target is dropped, 95 whose array is dropped (their phage target whole), 2 both dropped |
| truth_fusion | 58 whole and 4 cut on assembly contigs, 1 on a near copy, 3 dropped, 3 touched by an S8 edit (recoding in the code-4 clade) |
| truth_callerr | 6,080 whole on assembly contigs, 1,758 on near copies (the code-4 near duplicate's alt_code genes), 268 partial, 403 dropped |
| contigs | 2,327 emitted rows (2,217 edge rows, one per distinct uid), 321 dropped rows, 38 circular, 30 prophage-bearing |
| families, modules, traps | 2,810 families (host, NOGD, Cas, phage); 175 module rows; traps: 186 promiscuous-domain rows, 16 cargo, 5 hubs, 5 decoys, 1 clonal burst |
| sizes | `emit/` 169 MB on disk, of which the truth TSVs 56.3 MB (truth_gene 26.7 MB, truth_igr 24.6 MB) |

The lift-and-place chain was checked against an independent oracle (R1): a record carried through random S6
insertions, S7 and S8 edits and S9 cuts as an array of each base's S5 position, 2,333 random features: the Placer
agrees on every one (687 whole, 133 cut, 190 dropped, 1,323 touched by a later edit).

One finding about the world, not the emit: the code-4 near duplicate (g031's src_c) carries 1,758 alt_code
call-error rows, so `truth_callerr` counts a near copy's errors again; a scorer that counts errors per genome
should key on `true_gene_id` and `occurrence`.

## Part B: the emitted world through the normal builders

Commands (wilkcomp02, 2026-11-13, runtime image `/tmp/if2453/images/bio-box_64f1bf0.sif`, sha256
`cb7353f8…afba3`, tree of this branch at `/tmp/if2453/f028/tree`, scratch profile `/tmp/if2453/f028/profile.yaml`:
RUNBOOK Appendix B's, with the real 40/56 thread caps so that build-stage's allowance subtracts the account's
other work, `mem_scope` 16G, root `/tmp/if2453/f028/root`):

1. `f028_synthworld_builders.py inputs /tmp/if2453/f028/world /tmp/if2453/f028/inputs` (capped scope): stage
   **F0** = the original sources (46 host genomes with their 10 contaminant contigs, and the 12 phages: 58
   genomes, 2,132 contigs); stage **F1** = the duplicate sources (2 exact, 2 near: 4 genomes, 195 contigs).
2. `bash scripts/pod_run.sh … --threads 4 f028-F0 build-stage --stage F0 --release rF8.1 --stages F0 … --threads 4`
   (16:37-16:42, all 12 steps, unit totals 9.0 CPU-minutes, MemoryPeak 1.75 GB), then `f028-F1 build-stage
   --stage F1 --release rF8.2 --stages F0,F1` (16:43-16:47, 12 steps, 7.2 CPU-minutes, 1.69 GB); per-step wall,
   CPU and peak RSS in `steps.tsv` (clustering 106-107 s and 6.3 CPU-minutes, every other step under 33 s).
3. `f028_synthworld_builders.py compare /tmp/if2453/f028/world /tmp/if2453/f028/root rF8.2 /tmp/if2453/f028/out`
   (capped scope): reads the run records' file lists and `ledger.read`, never a glob.

No domain step: the world's panel is the `synthetic` one (Dirichlet rows fitted to Pfam statistics, spec §2.3 as
built), which exists only as `sequences/panel.tsv`, not as an HMMER file, so the normal domain step has nothing to
scan; a real Pfam scan would test nothing the world's truth knows (the families come from the synthetic panel).
No ncRNA step (card). Gene-call version `pyrodigal-3.7.1-gv-0.3.2-p2`.

### Contigs, duplicates, shared contigs (measured)

| check | result |
|---|---|
| every emitted contig stored once | 2,217 distinct uids, 2,217 default variants (F0 2,132, F1 85), none missing or extra |
| lengths and digests | all 2,327 occurrences: `contig_len` and `raw_digest` (sha256 of the FASTA bases) equal the emitted contig's; all 2,327 `contig_occurrence` rows map their accession to the right uid |
| exact duplicates (MR-06) | 108 member contigs, all with the original's uid, stored once by F0; their genome is the same `genome_uid`, so `contig_dim.n_genomes` is 1 for all 108 and `genome_dim` has 60 rows, equal to the truth's 60 distinct genomes, with kind, genetic code and assembly quality as the truth gives them for every genome |
| near duplicates | 87 member contigs: 85 new uids, both copies stored; 2 unchanged contigs (no site, above) share the original's uid and are the release's only contigs in two genomes (`n_genomes` 2; dim `shared_contigs` 2, DEC-27) |
| DEC-26 | F1's gene calling skipped exactly the 110 contigs F0 had stored (108 exact, 2 unchanged near), each with ledger reason `owned_by:F0`; their calls are F0's |

### Gene calls against `truth_gene` (measured; `recall_by_class.tsv`, `precision_by_origin.tsv`)

Stop-anchored match (spec §3.1): a truth gene is found when a call has its `contig_uid`, strand and `stop_bp`; a
piece whose 3' end a contig edge cut is found by a call in its frame (anchored at its 5' end) covering at least
half of it (pieces cut at both ends are never found by this rule). Start agreement: the call equals the whole
gene's start and end.

| class | n | found | recall | exact start of found |
|---|---|---|---|---|
| untouched (host CDS, code 11, whole) | 77,938 | 77,456 | 0.994 | 0.825 |
| call error: frameshift insertion | 29 | 23 | 0.793 | 0 (the call is the 3' piece) |
| call error: frameshift deletion | 28 | 23 | 0.821 | 0 |
| call error: internal stop | 34 | 27 | 0.794 | 0 |
| call error: upstream alternative start | 29 | 29 | 1.000 | 0.414 (the rest start at the planted ATG or elsewhere) |
| call error: downstream alternative start | 30 | 28 | 0.933 | 0.107 (most start at the competing ATG) |
| call error: smORF (10-50 aa) | 39 | 16 | 0.410 | 0.938 |
| code-4 genomes, gene with internal TGA | 4,921 | 4,904 | 0.997 | 0.824 |
| code-4 genomes, gene without | 684 | 673 | 0.984 | 0.834 |
| code-25 genome, gene with internal TGA | 970 | 962 | 0.992 | 0.835 |
| code-25 genome, gene without | 285 | 281 | 0.986 | 0.840 |
| contig edge, 3' end kept (stop match) | 1,906 | 1,749 | 0.918 | - |
| contig edge, 3' end cut (frame match) | 1,906 | 1,803 | 0.946 | - |
| across a circular origin | 36 | 36 | 1.000 | 0.861 |
| prophage genes (code 11, whole) | 1,027 | 1,005 | 0.979 | 0.823 |
| Cas genes | 96 | 96 | 1.000 | 0.885 |
| fused, split, NOGD genes | 44, 28, 12 | 44, 28, 12 | 1.000 | 0.841, 0.893, 0.917 |
| phage genomes (meta mode, pyrodigal-gv) | 730 | 704 | 0.964 | 0.824 |
| near-duplicate copies | 3,031 | 3,008 | 0.992 | 0.784 |
| contaminant contigs | 116 | 112 | 0.966 | 0.616 |
| **all truth genes** | 93,883 | 92,983 | 0.990 | 0.791 (partial pieces counted as not exact) |

Precision (calls matching a truth gene by the same rules, over 95,492 calls): assembly contigs 0.974 (2,410
unmatched, median 61 codons), near copies 0.973, contaminants 0.974, phages 0.982.

The code-11 pieces that a code-11 reading of a broken gene gives (`truth_callerr.expected_calls`) were found by
stop: frameshift deletions 52/57, insertions 49/57, internal stops 59/70, alternative starts 30/30 and 29/31,
smORFs 16/39. For `alt_code` rows 4,557/9,963 are found: the builders call code-4 and code-25 genomes with their
own code (`genome_input.genetic_code`), so a code-4 gene is called whole and only its 3' code-11 piece (which ends
at the true stop) matches; that is the intended behaviour, not a miss.

What this says (not tuned; the world and the builders are as they were): with the right genetic code in
`genome_input` the caller recovers recoded genes as well as untouched ones; a frameshift or internal stop is never
called whole (start agreement 0) and its 3' piece is called four times in five; planted smORFs are mostly missed,
as Prodigal's minimum gene length predicts for 10-30 aa (*estimate*); edge-cut genes are found about 92-95% of
the time. No builder defect was found: every check of the contigs, the duplicates and the shared contigs came out
as DEC-26 and DEC-27 say.

## Proposal: truth tables as registered Parquet (for the PI)

This session leaves the truth as TSV, because the registry (`src/biobox/registry/schemas.py`) and its locks are
outside the card's files. Registering them would need:

1. **Where they live.** Not in a release's data root. Truth describes a test world, not the store, and a store
   reader must never see it (a join from a release table to truth would make a test pass by construction). Proposed
   location: the world directory, `worlds/<world_id>/truth/<table>/part-0.parquet`, keyed by `world_digest`
   (spec §6.2, §8.2 risk 6: worlds are cached by digest and regenerated, not archived). The test data roots of R2
   and R5 would read them by path from the world's `manifest.json`, never by glob (AGENTS §3.9).
2. **Which registry.** A second registry module beside the store's, for example `biobox.synthworld.truth_schemas`
   with `TableSpec`s of layer `truth`, and its own lock file under `tests/registry/` (`truth_schema_lock.json`)
   updated by `scripts/update_schema_lock.py`, so a truth column cannot change silently and the store's locks are
   not touched. Until the first signed release the same LY-09 rule applies (append nullable columns only after).
3. **Tables and keys.** The 14 tables of `emit.TABLES`, with `world_id string` and `release_tag string` first.
   Coordinate tables (`truth_gene`, `truth_igr`, `truth_rna`, `truth_spacer_target`, `truth_fusion`,
   `truth_callerr`) sorted by (`contig_uid`, `start`, `end`, `strand`) with min/max statistics on those columns;
   `truth_contig` and `truth_edges` by `contig_uid`; `truth_genome` by (`genome_id`, `source_label`); `truth_dup`
   by `dup_group_id`; `truth_contam` by `contig_uid`; `truth_family` by `family_id`; `truth_module` by
   (`module_id`, `order_index`); `truth_traps` by `trap_id`. Row groups of 122,880 rows (AGENTS §4.6). Primary keys:
   `feature_uid` (genes, RNA, call errors, fused ORFs), `region_uid` (IGRs), (`array_feature_uid`, `spacer_id`,
   `target_contig_uid`, `proto_start`) for spacer targets, `contig_uid` for contigs (per occurrence:
   `accession`), never an ordinal.
4. **Types.** `contig_uid`, `feature_uid`, `region_uid`, `igr_uid`, `protein_uid`, `genome_uid`,
   `array_feature_uid`, `target_contig_uid`, `rep_contig_uid`, `member_contig_uid`: `FIXED_LEN_BYTE_ARRAY(16)`
   (hex in the TSV). Coordinates `uint32` (start, end, stop_bp, error_pos, alt_start, spacer and protospacer
   bounds, replicon_start/end, rotation_offset, length), nullable where the TSV cell is empty (an unplaced row);
   strand `int8` (+1/-1); `overlap_bp`, `genetic_code`, `n_mismatches`, `n_partial_genes` `uint8`;
   `split_domain_index`, `order_index`, `spacer_index` `uint16`; counts `uint32`; identities, GC, completeness,
   rates `float32`; flags `bool`; the rest `string`. List columns become nested types:
   `domain_arch list<struct<hmm_acc string, aa_start uint16, aa_end uint16>>`,
   `elements list<struct<kind string, start uint32, end uint32, strand int8>>`,
   `expected_calls list<struct<start uint32, end uint32, strand int8>>`, `gene_ids list<string>`,
   `partner_family_ids list<string>`, `domain_hmms list<string>`.
5. **The ingest manifest** (`ingest.tsv`) and the build-stage inputs stay outside the registry: they are inputs,
   which the builders read through `genome_input` and the ingest's source objects.
6. **Cost.** At `micro` the 14 TSVs total 56.3 MB (measured, part A); Parquet with dictionary encoding would be
   smaller. Nothing here needs a new dependency (pyarrow is locked).

## What remains of the card

- §3: the truth as registered Parquet (the proposal above, a PI decision); columns left empty (`ref_bitscore`,
  `ref_bitscore_median`, `gene_tree_extant`, `origin_node`, `strength_class`, which needs the binning of the
  strength grid, `realized_independent_transfers`, `realized_vehicle_independent_loci`, `powered_expected`);
  `truth_igr` rows for the junctions S6-S8 made (around prophages, Cas loci, smORFs and fused genes); domain
  architectures of Cas genes and smORFs; `tandem_paralog` traps; the `truth_module_pairs` view and the replicon-to-
  contig helper view of §3.3.
- §4: null worlds (N0-full, N0-strict, N1, N2) and the difficulty audit of §4.3.
- §6.2-§6.4: the macOS-vs-Linux digest; a `synthworld build` command line; `synthworld_version` and the git
  commit in `manifest.json` (left out: the tree in the image is not a git checkout, and a commit would make the
  digest depend on it); the `test_run` columns.
- §7: realism against the golden genomes (§7.1), tool parity beyond pyrodigal (§7.2: hmmsearch with a Pfam panel,
  MinCED, cmsearch), the rest of §7.3's internal checks (module contiguity, `ani_realized` within 0.02%, already in
  the S10 test).
- Left over from earlier sessions and still open: an outgroup world and vector-like contaminants; smORFs of
  uniform residues; call errors per genome, not inherited; the 3-minute target (S1-S11 227 to 232 s, above).

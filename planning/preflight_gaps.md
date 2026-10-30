# Pre-mortem and preflight gaps (before S0)

Status: plan only. Date: 2026-10-28. Search gaps #23–#28 and §1.1 added 2026-10-29 from `planning/search/` (`SEARCH_PLAN.md` calls them G-23..G-28). Test gap #29, K7 and the §7 test-harness group added 2026-10-29 from `planning/tests/TEST_ARCHITECTURE.md`. PI decisions DEC-01..DEC-12 recorded 2026-10-29 in `planning/DECISIONS.md`.

Premise: it is April 2027 and S0–S3 failed or stalled. This file lists the most likely causes that the current plan documents do not cover, and the preflight item that closes each one. Inputs read: `STAGES.md`, `00-spec.md` (v3.0), `PLAN.md` (v3.0), `DEPLOYMENT-AND-MALLEABILITY.md` (v3.0), `PRIOR-ART.md`, `planning/sources_and_sizing.md`, `planning/data_layer_and_tests.md`, `audits/SYNTHESIS_2.md`, and the pod notes (`WILKECOMP.org`, `fragment-defense/old-plan.2026-10-18.bak` §9).

Facts marked *(verify)* were not checked against a primary source for this file.

---

## 1. Ranked gap table

Severity: **B** means it blocks S0 from starting; **H** means S0–S3 likely fail without it; **M** means it causes rework or delay.

| # | Gap | Why it bites | Preflight action | Owner role | Blocking for |
|---|---|---|---|---|---|
| 1 | **No agreed first users or questions.** No document says which queries the lab needs first. | Six months of plumbing that answers nobody's question. The S0–S3 science checks are generic, not the lab's. | Write and sign off the "first 10 questions" (§3). Map each to the earliest stage that can answer it, and make two of them part of the S1 and S3 gates. | PI | S0 |
| 2 | **No code repository, CI or release process.** The environment notes say this is not a git repo; there is no test runner and no image build pipeline. | Nothing in `data_layer_and_tests.md` can run. Determinism and accretion tests need committed code and pinned images. | Create the repo (code, `golden.yaml`, `sources.yaml`, ADRs, tests). Set up CI: unit and property tests plus a *micro-S0* fixture (3 genomes, under 10 min) on the R2 CI host (GitHub Actions or another, chosen at S−1; TEST_ARCHITECTURE §11 item 2); the full micro-corpus (~50 synthetic + 3 real genomes, under 15 min) follows at S0; image build to the registry; the `.sif` pulled on the pod by digest. Choose a code license (see #18). | Data engineer | S0 |
| 3 | **Cross-document contradictions** (§2). Stage labels in `PLAN.md` still follow the pre-STAGES order. ID hash width and algorithm differ across three files. The sequence-store format is stated three ways. | The implementer builds the wrong schema or gates the wrong stage. IDs are the one decision that can't be changed cheaply. | Resolve every row of §2 in one ADR pass. Add a CI doc-lint that greps stage labels, hash algorithm and ID width against `STAGES.md` and the ID ADR. | Data engineer + PI | S0 |
| 4 | **FCS-GX memory vs the 200 GB cgroup.** `00-spec.md:47` and `PLAN.md:114` put FCS-GX screening in Stage 0b. The FCS-GX database is about 470 GB and is meant to be held in RAM/shared memory *(verify)*. | Stage 0b can't run on wilkcomp02 under the house rules. It is noted as "check in S0" with no fallback. | **Decided 2026-10-29: deferred, accreted later.** FCS-GX runs when a host with ≥512 GB RAM is available (rented box or later hardware) and lands as a per-genome QC sidecar (flagged contigs and regions with taxon and reason code). The next release recomputes statistics with those contigs masked; sequences, tapes and IDs are unchanged. Until then: CheckM2 + GUNC for quality and chimerism, plus a low-RAM human/vertebrate host screen from S0 (minimap2 or sketch containment against GRCh38/T2T and common hosts), so gap #16 holds. Every genome carries `fcs_gx_status` (`ncbi_screened`, `screened_here`, `not_screened`), and results from releases with unscreened genomes carry the flag. One-off check: FCS-GX on a few hundred older NCBI assemblies on a rented box, to measure what deferral misses. Test MR-19 (TEST_ARCHITECTURE §4). | PI + pod owner | S0 (S3 for non-NCBI sources) |
| 5 | **Pod approvals not requested.** The plan assumes multi-week runs, TB-scale scratch, a possible memory exception (S5), and that wilkcomp02 is free. Under the house rules, anything beyond the agreed footprint is "a request to Claus, not a notice". | A stop request mid-S5, or watchdog mail to rctf-support every hour, which damages the relationship. | Send one written footprint request covering S0–S5: threads, memory, scratch growth to about 2 TB, run lengths, a STOP-file pause guarantee, and a README per large directory. Also ask about the wilkcomp03 GPU for Stage 2c (#9) and a second account for a postdoc. File the reply in the repo. | PI | S0 (memory exception: S5) |
| 6 | **L0 snapshots live only on unbacked scratch.** `DEPLOYMENT-AND-MALLEABILITY.md:30` treats source snapshots and the sequence store as "regenerable" on `/stor/scratch`. But upstreams change: INPHARED is rolling, IMG/VR needs a login, and GTDB retires releases. | A scratch loss or upstream change makes a stage unreproducible, and the accretion-equivalence test can no longer be run. | Classify L0 as **durable**; only L1+ is "regenerable". **Decided 2026-10-29 (DEC-11):** manifests, checksums and `contig_occurrence` on `/stor/work`; the PI backs L0 snapshots up to the lab NAS (date to set), sources that version away first (INPHARED is rolling, IMG/VR needs a login, GTDB retires releases). The earlier proposal (second copy on `/stor/work` while under ~50 GB, object storage from S5) remains available if the NAS copy is late. | Data engineer | S0 |
| 7 | **Reference databases not inventoried or pinned** (§4). Sizes, licenses, versions and RAM needs for Pfam, NCBIfam, PHROGs, Rfam, eggNOG, the CheckM2 DB, the geNomad DB, iPHoP, Foldseek AFDB/ESMAtlas30 and FCS-GX are not in any file. | Surprise downloads in the hundreds of GB. A database version bump mid-stage changes `annot_version` and breaks accretion equivalence. | Fill `refdb.yaml` (§4) with version, URL, sha256, size, RAM, license and redistribution flag. Download and checksum into `l0/refdb/` during preflight. Freeze for S0–S3. | Bioinformatician | S0 |
| 8 | **No owners, calendar or time box.** STAGES gives core-hours but no weeks, people or decision rights. | Stages drift; nobody can say "stop". | Use the role table and calendar in §5. Each stage gets a time box; overrunning it by 50% triggers the kill review (#12). | PI | S0 |
| 9 | **The GPU on the pod is ignored.** `00-spec.md:31` says "No GPU: ProstT5 and gLM stages are out of scope here", but wilkcomp03 has 4× RTX 6000 Ada (WILKECOMP §2). | Renting a GEX130 for Stage 2c is unnecessary cost, and dark-fold assignment is delayed to S6+. | Add wilkcomp03 to the target table as "GPU, by request". Ask in the #5 request. Probe its disk and load before relying on it. | PI | S3 (first dark-segment set) |
| 10 | **Hetzner rehearsal (P1) has no budget or billing model.** STAGES §5 says "rent for a few days". Hetzner dedicated servers bill monthly and some carry setup fees *(verify)*. Object storage, egress and the AWS Logan budget have no numbers. | P1 slips, or a surprise invoice arrives. | Get quotes for: an AX162-class box for one month, or Hetzner Cloud CCX (hourly) for the rehearsal; object storage for 1–2 TB; R2/S3 egress. Set budget caps and alerts. Record the figures in `planning/costs.md`. AWS Logan costing is deferred to an S7 gate. | PI | P1 (after S2) |
| 11 | **Tool and database upgrade policy is undefined.** Pfam, Rfam, geNomad and MMseqs2 release on their own schedules. The pod still has MMseqs2 15, and the plan needs ≥18. | A mid-stage upgrade silently changes outputs; accretion equivalence fails for reasons nobody can explain. | Policy: versions are frozen per release train. Bumps happen only at stage boundaries, through a "re-annotation release": the L1 append under a new `annot_version`, plus golden regression, a drift report (G-05) and a changelog entry. Pod-installed tools are never used; only the `.sif` is. | Bioinformatician | S0 |
| 12 | **No exit or kill criteria.** The plan defines pass gates but not when to stop the project or a stage. | Sunk-cost continuation. | Use §6. | PI | S0 |
| 13 | **Golden set not fully pinned.** SR1, Lak and Obelisk accessions are marked "pin". STAGES says ~14 genomes, while `data_layer_and_tests.md:152,303` says ~16 including CAMI and spike-ins. Truth tables (RegulonDB, SubtiWiki, M-CSA extracts) are not yet downloaded or versioned. | B-G can't be scored, and S0 can't exit. | Write `golden.yaml` with accessions, truth sources and versions, and a per-genome expected-result table. Each truth extract is an L0 source with its own contract. | Bioinformatician | S0 |
| 14 | **Volatile `/tmp` and long runs.** Whether `/tmp` is cleaned on reboot is unknown (WILKECOMP "Could not verify"). S5 runs for weeks. | A reboot wipes MMseqs2 working state mid-pass and the whole pass reruns. | Every monolithic step checkpoints per shard to NFS through the one-rsync rule. Ask rctf-support about the tmpfiles policy (#5). D-02 kill-resume tests must include "`/tmp` wiped". | Data engineer | S1 |
| 15 | **No postdoc onboarding.** The account `if2453` must not be shared. There is no runbook or data dictionary. | Only one person can run anything, and the plan has a bus factor of 1. | Request a second pod account (#5). Write `RUNBOOK.md` (launch, pause, resume, STOP, cleanup, who to email), a generated schema/data dictionary, a query cookbook notebook answering the first-10 questions, and a glossary. | Data engineer | S1 |
| 16 | **Human-derived sequence exposure** is only addressed for Logan and GSA-Human. Host contamination in MAGs and SAGs is not, and screening rejects are retained in "quarantine". | Retaining and serving human-read-derived contigs creates a privacy and policy problem, even from public SRA. | Rule: contigs classified as human or vertebrate host are rejected and **not** retained in quarantine. Only their digest and reason code are kept. No SRA-derived inputs before S7; any later use passes STAT preselection plus FCS-GX. Write this into `sources.yaml`. | PI | S0 (rule); S7 (enforcement) |
| 17 | **Build vs reuse is not decided per question.** PRIOR-ART now names tools that already answer slices: webFlaGs, GCsnap2, Gaia, IMG/VR's own neighborhoods. | S0–S3 rebuild what a lab member could query today, so nobody uses it. | For each of the first-10 questions, record the best existing tool and its answer as the external baseline. A stage counts as "useful" only when bio-box answers at least one question that existing tools can't, or answers it materially faster. | PI | S1 |
| 18 | **Code license and the name.** `bio-box` collides with the existing "Bioboxes" container standard (bioboxes.org, Belmann et al. 2015) *(verify)* and possibly with registry or PyPI names. INPHARED scripts are AGPL; vendoring them would put AGPL terms on bio-box. | Rename churn after release, or license contamination. | Check the registry, PyPI and GitHub names. Pick a code license (for example MIT or Apache-2.0) and call INPHARED as an external tool instead of vendoring it. | PI | S0 (license); before any public release (name) |
| 19 | **Reuse of the existing pod harness is not stated.** `homology_pod.sh` and the fragment-defense pod scripts already implement nice, cgroup, STOP and rsync. The fragment-defense project already holds DefenseFinder and PADLOC runs on phages. | The wrapper and its pod bugs get reinvented; free validation data goes unused. | Port `homology_pod.sh` behaviors into `bio-box run-step` and list its known pitfalls (§9.3a/b). Use fragment-defense's phage DefenseFinder/PADLOC outputs as an S1 cross-check sidecar. | Data engineer | S0 |
| 20 | **Oracle and test compute on the pod.** Brute-force oracles at S3+ and the 1–16-client DuckDB tests consume pod CPU and NFS bandwidth. | Tests trip etiquette or the watchdog, or get skipped. | Run fan-out test shards at nice 15 and monolithic R3 steps at nice 19, each under 300 CPU-min (TEST_ARCHITECTURE §8). Run multi-client tests at night only. Budget test core-hours per stage in `sizing.yaml`. | Data engineer | S1 |
| 21 | **No data management plan.** If grant-funded, the NIH DMS policy may apply. Citations and attribution obligations (CC BY-SA for GTDB, CC-BY for CAMI) accumulate. | Packs ship without required attribution. | Generate `CITATION.cff` and an attribution file per release from `sources.yaml`. Record the funding source and DMS obligations. | PI | S2 (first pack) |
| 22 | **Documentation drift is likely.** There are about 3,400 lines of plan across 14 files. | Readers act on stale text, as happened with the round-1 PRIOR-ART. | One authoritative file per topic: STAGES for stages, the ID ADR, the `sizing.yaml` ledger. Other documents link rather than restate. The doc-lint from #3 enforces this. | Data engineer | S0 |
| 23 | **No search command or runner.** The CLI and both query-class lists take tokens, UIDs or SQL, never a sequence, HMM, MSA, CM or structure (`planning/search/06_infra_inventory.md` §3 item 1; cited as `06` below). | 348/365 synteny hypotheses need a neighborhood step and 89 need custom profiles; minute-to-hour runs cannot live under `serve` timeouts. | Resolved in plan by `planning/search/SEARCH_PLAN.md` SD-01 (with SD-08, SD-10, SD-15), pending implementation. | Data engineer | S0 (ST0 recipe runner) |
| 24 | **No search pack.** FASTA exists only on demand; MMseqs2 DBs are transient; no NR90/NR50, IGR or HH targets ship (`06` §3 item 2). | Extracting all S5 proteins at the R-08 rate takes ~16 h (*derived*) before any search starts. | Resolved in plan by `planning/search/SEARCH_PLAN.md` SD-06, pending implementation. | Data engineer | S1 (ST1 pack) |
| 25 | **Member alignments not stored.** PLAN cites a `protein_member` alignment table; the L2 layout has none, and 16 B/protein sizing is too small for it (`06` §3 items 10–14). | Hits on representatives cannot be projected to members, so `scan-new-models` and NR90 final rounds have no exact path. | Resolved in plan by `planning/search/SEARCH_PLAN.md` SD-03 and SD-05 (`cluster_member_aln`, ~19 GB at S5, *estimate*), pending implementation. | Data engineer | S1 (ST1 projection) |
| 26 | **Z undefined.** "Fixed per release" with no value or definition; a per-stage Z changes stored E-values and breaks A-02 (`06` §3 items 16–17). | E-values are not comparable across stages or between library and user searches. | Resolved in plan by `planning/search/SEARCH_PLAN.md` SD-04, pending implementation. Constants open in §1.1. | Bioinformatician | S0 |
| 27 | **Custom hits cannot reach statistics.** Tokens come from frozen dictionaries; sidecars join after the filter; the sidecar coverage rule could refuse sparse hit lists (`06` §3 items 4–5). | A user family never reaches `slot_pairs`, `pair_stats`, planes or `rescore`. | Resolved in plan by `planning/search/SEARCH_PLAN.md` SD-09 (coverage rule clarified in its §3.1), pending implementation. | Data engineer | S2C (ST2) |
| 28 | **Nucleotide search outside gene calls absent.** No six-frame or unannotated-ORF search and no persistent nucleotide DB (`06` §3 item 6). | ~85% of RODEO lasso precursors were not GenBank CDS; absence claims without dark-ORF checks are unsafe. | Resolved in plan by `planning/search/SEARCH_PLAN.md` SD-07 (IGR, dark-ORF and spacer corpora; `bio-box region translate`), pending implementation. | Bioinformatician | S0 (ST0 corpora) |
| 29 | **Tests are organized by layer and stage, not by operation.** Families C, L, R, Q, A, D, G, K and S exist, but nothing makes a new workflow inherit coverage; there is no planted truth or null world for the statistics, no metamorphic relations, no mutation testing, and no results table. | New workflows ship untested; S1-severity failures (strand errors, ID leakage, miscalibrated nulls) pass golden snapshots. | Resolved in plan by `planning/tests/TEST_ARCHITECTURE.md` (operation registry O-*, F-01..F-15, MR-01..MR-19, workflow manifests and coverage matrix, micro-corpus, synthetic and null worlds, LY-01..LY-10, rings R0–R5, `test_run`, mutation testing, fault injection for every S1 failure), pending implementation. Spec items in §7; open items in TEST_ARCHITECTURE §11. | Data engineer + bioinformatician | S−1 (spec), S0 (micro-corpus and R0–R3) |

### 1.1 Open items from the search round (2026-10-29)

| Item | What is open | Decides | Needed by |
|---|---|---|---|
| S2C approval | **Approved 2026-10-29:** stage S2C (GTDB R226 `Complete Genome` assemblies, ≤20 per species by CheckM2 quality then ANI diversity; *estimate* 15–30k assemblies, 60–120 Gbp) between S2 and S3 (`SEARCH_PLAN.md` §6.1). Optional follow-on S2C-b (the remaining complete genomes, uncapped) is added after the S2C gate only if the signed-off first-10 questions include strain-level work (defense hotspots, pangenome turnover). The exact count is still an S−1 task. | PI | S−1 (count) |
| Z constants | **Decided 2026-10-29 (DEC-04):** `Z_prot = domZ = 1e8` sequences and `Z_nt = 1e5` Mb, one constant across the release line (SD-04) | PI + bioinformatician | done |
| ST4 search host | Size and timing of a Hetzner search host for the S6 pack: ≥512 GB RAM, ≥4 TB NVMe (*estimate*). Quote it with P1 (#10). | PI | S4 / S6 (ST4) |
| Synteny corpus repair | About 46 `lit/` files hold a different paper than their name (9 in 2022–2026, ~25 in 2017–2021, 12 in 1993–2016), and the synteny audit counts are zero for all 365 rows. Replays use supplements and are unaffected; hypotheses drawn from mislabeled texts need re-extraction before they drive experiments. | Bioinformatician | before any hypothesis drives experiments |
| Replay supplement fetch list | Doron 2018 S1–S3, Gao 2020 S1–S3, Altae-Tran Data S1–S4, Millman Table S1, Brewer 2021 files 1–2, CRISPR-Cas 2025 Zenodo PSSMs, Domainator Zenodo; each pinned as an L0 source | Bioinformatician | S−1 (each tier's gate uses them) |

---

## 2. Cross-document contradictions found

| # | Topic | Statement A | Statement B | Resolve to |
|---|---|---|---|---|
| C1 | Stage labels | `STAGES.md:91-94`: S4 = IMG/VR, S5 = GTDB species reps, S6 = all GTDB + RefSeq flag | `PLAN.md:156-157`: GTDB reps (S4), RefSeq + multi-strain (S5), IMG/VR (S6); `PLAN.md:48`: same old order; `00-spec.md:41`: same old order | STAGES labels; fix PLAN §0 and §2 tables |
| C2 | Benchmark stages | `STAGES.md:245` B1 at S5; `STAGES.md:276` B3 at S6 | `PLAN.md:513` B1 at S4; `PLAN.md:515` B3 at S5; `PLAN.md:516,518,519` B4/B6/B-EXT at S4 | Shift PLAN §5 by STAGES numbering (B1, B4, B6, B-EXT → S5; B3 → S6) |
| C3 | Content-UID hash | `data_layer_and_tests.md:54-59`: blake3 truncated to 64 bits (`[:8]`, uint64); `raw_digest` = sha256 | `DEPLOYMENT-AND-MALLEABILITY.md:127`: BLAKE3 truncated to **128 bits**; `sources_and_sizing.md:48`: "xxh128 or SHA-256"; `critique_7_database.md:154`: xxh128 | One ADR. At 64 bits the birthday bound is ~1 collision per 6×10⁹ items (derived: n²/2⁶⁵ ≈ 1 at n ≈ 6×10⁹). That is reachable at S7, so use 128-bit BLAKE3 for UIDs, uint64 dense IDs per release, and keep sha256 only for file/object integrity |
| C4 | Sequence-store primary format | `sources_and_sizing.md:47`: bgzip+faidx primary, zstd-seekable optional | `data_layer_and_tests.md:80`: Parquet `seq_chunks` primary, zstd-seekable **rejected**; `DEPLOYMENT-AND-MALLEABILITY.md:193,207`: Parquet primary plus optional bgzip or zstd-seekable mirror; `PLAN.md:351`: three options, decided by "Q7" | Decided by K-13 at S0 (`STAGES.md:122`). PLAN's "Q7" should read K-13. All docs say "decided at S0 by K-13" until then |
| C5 | IGR storage | `00-spec.md` invariant 2, `sources_and_sizing.md:46,58`: IGRs stored (0.032 GB/Gbp budgeted) | `data_layer_and_tests.md:84`: IGRs are coordinates only, sliced on demand | Coordinates only (no second copy of bases); sizing drops the 0.032 term; spec wording "retrievable" not "stored" |
| C6 | `lite` edition size | `DEPLOYMENT-AND-MALLEABILITY.md:250`: `lite` ~85 GB with download times | `DEPLOYMENT-AND-MALLEABILITY.md:305` and `STAGES.md:331`: `lite` = S0–S5, whose sequence store alone is ~207 GB (`STAGES.md:93`) | Replace the size row with "estimate after S5"; drop the download-time row or derive it from `sizing.yaml` |
| C7 | Release naming | `STAGES.md:62,142`: `rS<k>.<n>` (e.g. `rS0.1`) | `DEPLOYMENT-AND-MALLEABILITY.md:53`: `--release r0.1` | `rS<k>.<n>` |
| C8 | Golden set size | `STAGES.md:88,113`: ~14 genomes + CAMI subset | `data_layer_and_tests.md:152,303`: ~16 incl. CAMI + spike-in panel | `golden.yaml` is authoritative; STAGES lists 14 genomes + CAMI subset + spike-in panel |
| C9 | GPU on the pod | `00-spec.md:31`: "No GPU" | WILKECOMP §2: wilkcomp03 has 4× RTX 6000 Ada | Add wilkcomp03 as a by-request GPU host (gap #9) |
| C10 | Invariant count | `00-spec.md:39`: fourteen | `PLAN.md:44`: "The spec lists twelve" (historical v2.0 correction, reads as current) | Reword PLAN:44 to "v2.0 said ten but listed twelve; v3.0 has fourteen" |

None of these require new design; they are edit-level fixes. C1–C4 must be fixed before S0, because they determine code.

**Status 2026-10-28:** C1–C10 were resolved in the documents. They now use the STAGES labels and benchmark stages, 128-bit BLAKE3 UIDs (sha256 for exact bytes and objects), Parquet `seq_chunks` as primary storage pending K-13, IGRs as coordinates, `lite` resized, `rS<k>.<n>` release names, 16 golden entries, wilkcomp03 as a GPU option, and the invariant count reworded. The remaining preflight item is the CI doc-lint that keeps them resolved.

---

## 3. First 10 questions (draft, for PI sign-off)

These are written to exercise the plan's distinct capabilities and to match the lab's current work (defense and anti-defense systems, retrons, DNA repair). The PI should replace or rank them. **Earliest stage** is the first stage whose data can answer the question at all; **useful at** is when the answer is likely better than existing tools.

| # | Question | Capability exercised | Earliest stage | Useful at | Existing-tool baseline |
|---|---|---|---|---|---|
| Q1 | For a bait family (e.g. a retron RT), return every neighborhood ±10 genes with raw contig FASTA and coordinates | Neighborhood fetch + sequence store | S1 | S4 | webFlaGs, GCsnap2, IMG/VR |
| Q2 | Which phage-encoded genes sit within ±5 ORFs of known anti-defense genes (Acr, Ocr-like, anti-RecBCD) more often than passenger expectation? | Score test + Tier 2 geometry test | S1 | S4 | Manual DefenseFinder/PADLOC runs (fragment-defense) |
| Q3 | Which uncharacterized families are conserved neighbors of RecBCD/AddAB across orders? | Anchor scan + dispersion | S3 | S5 | STRING neighborhood channel |
| Q4 | On plasmids, which pairs are coupled modules (TA, partition) and which are AMR–transposon passengers? | Vehicle dedup + passenger traps | S2 | S2 | Coinfinder (pangenome-level) |
| Q5 | Export all members of a subfamily with neighborhoods, sequences and taxonomy, for SSN/GNN or structure prediction | Subfamily tokens + export | S1 | S5 | EFI-GNT |
| Q6 | Which defense-system cassettes vary between strains of one species, and with what flanking vehicles? | Pangenome planes + vehicles | S6 | S6 | PADLOC/DefenseFinder per genome |
| Q7 | In reduced genomes, which pathways are lost entirely vs holed with an unexplained replacement? | Missing-step search with masks | S3 | S6 | GapMind |
| Q8 | Which riboswitch subtypes sit upstream of a given transporter family? | ncRNA per-member subtyping + tape | S3 | S5 | Rfam + manual |
| Q9 | Which intron/HEG insertions occur in phage core genes, and in which lineages? | Split-gene reconstruction | S1 | S4 | Manual |
| Q10 | Taxonomic distribution (tree-ordered) of a family or module, with closed-genome support | Presence bitmaps + quality mask | S3 | S5 | AnnoTree |

Gate additions: Q1 and Q4 become part of the S1 and S2 gates (answerable, with sequences, by `bio-box` CLI and a cookbook notebook). Q3 and Q10 become part of the S3 gate. If no question beats its baseline by S4, trigger kill review K3 (§6).

---

## 4. Reference databases to inventory (`refdb.yaml`)

All sizes, RAM needs and licenses are *(verify)*; the purpose of this table is the checklist, not the values.

| Database | Used for | Version to pin | Approx. size | RAM concern | License / redistribution |
|---|---|---|---|---|---|
| Pfam-A HMMs | domains, Codetta | current release at freeze | ~1–2 GB | no | CC0 |
| NCBIfam | domains | current | ~1–3 GB | no | public domain (NCBI) |
| PHROGs HMMs | phage domains | v4+ | ~1 GB | no | check |
| arCOG | archaeal domains | 2018/2020 release | small | no | check |
| eggNOG (for COG) | COG assignment | 5.x/6 | ~50 GB | moderate | check |
| Rfam CMs | ncRNA | 15.1 (4,227 families) | <1 GB | no | CC0 |
| Orphan CMs (Weinberg/Breaker etc.) | ncRNA | per paper | small | no | per paper; cite |
| MacSyFinder models (DefenseFinder, CONJScan, TXSScan) | system sidecars | pinned model versions | small | no | check per model package |
| PADLOC DB | system sidecars | pinned | small | no | check |
| geNomad DB | regime, viral taxonomy | pinned | ~few GB | low | check |
| CheckV DB | viral completeness (S1, S4 calibration) | pinned | ~few GB | low | check |
| CheckM2 DB | MAG/SAG QC (S7; skip if GTDB metadata suffices) | pinned | ~3 GB | low | check |
| GUNC DB | chimera screening | pinned | ~10–15 GB | moderate | check |
| FCS-GX DB | contamination | pinned | ~470 GB | **high (gap #4)** | public (NCBI) |
| GTDB R226 metadata, trees, protein reps (for taxonomy) | taxonomy, tree order | R226 | tens of GB | MMseqs2 taxonomy DB build | CC BY-SA (verify) |
| GTDB-Tk R226 package | only for genomes without GTDB taxonomy (S7) | R226 | ~100+ GB | moderate | GPL tool; data per GTDB |
| iPHoP DB | virus hosts (S4) | pinned | ~100+ GB | moderate | check |
| Foldseek AFDB/ESMAtlas30 3Di DBs | fold transfer (D-10) | pinned | hundreds of GB | MMseqs2 search memory | CC-BY 4.0 |
| Truth sets (EcoCyc/RegulonDB, SubtiWiki, M-CSA, Rfam seeds) | golden regression | pinned extracts | small | no | **EcoCyc/BioCyc licensing (verify)**; RegulonDB academic terms |

Preflight: every row has version, URL, sha256, measured size, and license filled; rows needed for S0–S3 are downloaded and frozen.

---

## 5. Roles, decision rights and calendar

| Role | Holder (to assign) | Decides |
|---|---|---|
| PI / product owner | Ilya | Scope, first-10 questions, kill reviews, budget, licenses, pod requests |
| Data engineer | to assign (postdoc or staff) | Repo, CI, wrapper, data layer, tests, releases |
| Bioinformatician | to assign (may be the same person) | Tool versions, refdb, golden set, annotation rules |
| Independent curator | someone not tuning scores | Trap panels, B5/B6 sets (pre-registered) |
| Test harness owners (gap #29) | data engineer (schema and operation registries, manifests, rings R0–R4, `test_run`, CI host, mutation testing) + bioinformatician (`synthworld` spec and the ZOMBI2 fit spike, truth tables, oracles, null-world calibration of O-STATS) | Harness spec sign-off at S−1; S0 harness gate (R0–R3, LY-06 layout calls with the data engineer) |
| Flexibility-check lab member | to assign; not a harness author (TEST_ARCHITECTURE §11 item 4) | None; adds one workflow from a manifest and recipe only |
| Pod owner | Claus Wilke | Footprint exceptions, wilkcomp03 access, second account |
| Facility | rctf-support | Watchdog, `/tmp` policy, quotas |

Calendar (time boxes, assuming one person at ~60% time; revise after S0):

| Phase | Time box | Ends with |
|---|---|---|
| Stage −1 preflight | 3–4 weeks (*estimate*; +1 week for the harness spec) | checklist below green, test-harness group included |
| S0 | 6–8 weeks (*estimate*: 4 for the S0 build plus 2–4 for the TEST_ARCHITECTURE §10 harness build) | `rS0.1`, `sizing.yaml` refit, ADRs frozen, harness gate (R0–R3, LY-06, mutation score, null-world calibration) as a `test_run` query |
| S1 | 2 weeks | `rS1.1`, Q1 answered |
| S2 | 2–3 weeks | `rS2.1`, Q4 answered, DuckLake/LIST decisions |
| S2C | 2–5 weeks (*estimate*: 3,000–22,000 core-h at ~800 core-h/day, plus gates) | `rS2C.1`, ST2 search replays, first stage-scale null-world calibration (K7) |
| P1 | 1 week | cross-host hash equality |
| S3 | 3 weeks | `rS3.1`, Q3/Q10 answered, S5 go/no-go from scaling fit |
| S4 | 3 weeks | `rS4.1` |
| S5 | 6–14 weeks (CPU-bound) | `rS5.1`, B1 with CI, `lite` edition |
| S6 | not boxed yet (split pod/Hetzner; set at the S5 gate from the scaling fit) | `rS6.1`, `core` edition |

Overrunning any time box by 50% triggers a kill review.

---

## 6. Exit and kill criteria

| ID | Trigger | Action |
|---|---|---|
| K1 | S0 measured rates predict S5 > 150 days on the pod **and** no Hetzner budget approved | Stop at S4; ship `lite` without S5; re-plan |
| K2 | B-G gene-calling or ncRNA targets fail after two fix iterations | Stop; adopt an existing annotation (e.g. PGAP/Bakta outputs) for Tier 1 and re-plan the annotation layer |
| K3 | By the end of S4, no first-10 question is answered better than its baseline | Stop building; keep the sequence store and per-genome annotations as a lab resource |
| K4 | S2 trap FPR 95% upper bound > 10% after the D-13 statistics are implemented | Pause statistics work; the engine ships as a retrieval tool only |
| K5 | Accretion equivalence (A-tests) cannot be made to pass by S2 | Drop incremental L2; rebuild L2 per release and re-budget |
| K6 | Pod owner withdraws or restricts the footprint | Move to P1 host early; costs come from `planning/costs.md` |
| K7 | (Evaluated after K4: if K4 has paused statistics at S2, K7 is evaluated only once the K4 fix lands.) O-STATS still fails null-world calibration (MR-18: p-values not uniform by KS, or FDR above nominal at the binomial 95% bound) at S2C scale after two fix iterations | Redesign the statistics layer (nulls, dependence model) before S5; until then statistics ship flagged uncalibrated and retrieval and search continue |

---

## 7. Stage −1 preflight checklist (all green before S0 starts)

Governance
- [ ] First-10 questions signed off by PI, each mapped to a stage and an existing-tool baseline (§3). **Deferred** (DEC-05): ask again before the S1 gate; the draft Q1/Q4 stand as provisional S1/S2 gate items until then.
- [x] Calendar and time boxes accepted (DEC-08). Role holders still to name in §5.
- [x] Kill criteria K1–K7 accepted (DEC-09).
- [ ] Footprint request approved (DEC-10); draft in `planning/footprint_request.md`. Still to do: send it to Claus Wilke (threads, memory, scratch to ~2 TB, run lengths, STOP guarantee, wilkcomp03 GPU, second account) and reply filed; `/tmp` cleanup policy asked of rctf-support.

Decisions (ADRs committed)
- [ ] Contradictions C1–C10 resolved; doc-lint added to CI.
- [ ] ID ADR: 128-bit BLAKE3 content UIDs, uint64 per-release dense IDs, sha256 for object integrity, `id_history`.
- [x] FCS-GX decision: deferred, accreted later as a QC sidecar; human/vertebrate host screen from S0 (gap #4, 2026-10-29).
- [ ] Human-data rule written into `sources.yaml` (no retention of human-classified contigs).
- [ ] Upgrade policy (frozen per release train; bumps only at stage boundaries).
- [x] L0 durability rule (DEC-11): manifests and checksums on `/stor/work`; snapshots to the lab NAS by the PI (date to set), sources that version away first.
- [x] Code license: none for now, private repo (DEC-07).
- [ ] Name collision checked (`bio-box` vs Bioboxes, PyPI, registries); needed before any public release.
- [x] S2C approved, ≤20 per species (2026-10-29).
- [ ] Complete-genome count from GTDB R226 metadata recorded (§1.1).
- [x] Z constants fixed for the release line: `Z_prot = domZ = 1e8`, `Z_nt = 1e5` Mb (DEC-04).
- [x] Reference databases frozen for S0–S3, upgrades only at release boundaries (DEC-12); `refdb.yaml` still to fill.

Infrastructure
- [ ] Git repo created with `golden.yaml`, `sources.yaml`, `refdb.yaml`, `params.yaml` stubs, ADR folder, `RUNBOOK.md` skeleton.
- [ ] CI: unit and property tests plus micro-S0 fixture green; image built and pushed by digest.
- [ ] `.sif` pulled on wilkcomp02 by digest; `bio-box doctor` output recorded (CPU flags, cgroup delegation, `ulimit`, filesystems, singularity/Apptainer version, DuckDB version).
- [ ] `bio-box run-step` wrapper ported from `homology_pod.sh` (nice 19/15, systemd-run MemoryMax, thread env vars, STOP file, `/tmp` staging, one rsync, small-file guard, ledger row), with a kill-resume test including a wiped `/tmp`.
- [ ] Directories created with setgid and README (what, who, until when, safe to delete): `/stor/work/Wilke/if2453/all-the-data`, `/stor/scratch/Wilke/if2453/all-the-data`.

Test harness, on paper (gap #29; `STAGES.md` S−1 item 15; § numbers in this group are `planning/tests/TEST_ARCHITECTURE.md`'s)
- [ ] Schema registry specified: pyarrow schema plus semantic constraints for every L0–L4 and search table, 16-byte UIDs as FIXED_LEN_BYTE_ARRAY (R0, LY-01). **v0 done** (`src/biobox/registry/schemas.py`, 17 tables); the rest are listed in `tests/registry/pending_tables.yaml`.
- [x] Operation registry (`tests/registry/operations.yaml`): every O-* row with inputs and outputs, oracle, MRs and existing tests (§2); R0 lint green.
- [ ] Synthetic-world generator (`synthworld`) spec: planted features and truth tables (§6.2), null worlds with phylogeny and HGT. **Spec written** (`planning/tests/synthworld_spec.md`; ZOMBI2 backbone plus custom code, Zombi 1 and ALF not used); S−1 fit spike done 2026-10-29: GO with conditions (`planning/tests/synthworld_spike.md`, DEC-15).
- [x] CI host for R2 chosen: GitHub Actions on private `ifinkelstein/bio-box` (DEC-06).
- [x] `test_run` schema registered (fields in §8).
- [ ] S0 gate written as a SQL query over `test_run`.
- [ ] S-12 manifests on paper: the three Domainator case studies and five synteny hypotheses as recipes, each with a workflow manifest (§5.1).

Data
- [x] `golden.yaml` pinned (2026-10-29): 15 genomes and phiX verified against NCBI; SR1, Lak (A1 plus a B8 contrast) and Obelisk chosen; checksums in `golden.lock.yaml`.
- [x] CAMI II marine sample chosen: short-read sample 0 gold-standard contigs (`golden.yaml`, fetched on the pod at S0).
- [ ] Truth extracts (EcoCyc, SubtiWiki, RegulonDB, M-CSA) downloaded as L0 sources with contracts.
- [ ] `refdb.yaml` filled for every S0–S3 database; downloads checksummed and frozen (§4).
- [ ] `sources.yaml` contracts drafted for S0–S3 sources, with license and redistribution flags (GTDB taxonomy, PLSDB terms checked).
- [ ] Replay supplements on the fetch list (§1.1) downloaded and pinned as L0 sources.

Budget
- [ ] `planning/costs.md` with quotes for P1 (dedicated monthly vs cloud hourly), object storage, egress; caps and alerts defined.

# B-47 design: Pfam by representatives, with inferred rows (DEC-32)

Card B-47 (design half, tier H). DEC-32 decides *that* S2C phase 2 is annotated by scanning NR90 representatives
(plus members whose length differs by more than 20%) and inferring the rest. This note decides nothing; it lays out
*how*, so that the PI can choose in one reading. It reads the code on `plan/task-cards` at `466d609`
(`schemas.py`, `l1/domains.py`, `l2/cluster.py`, `stage.py`, `l3/tokens.py`, `l3/postings.py`, `catalog.py`) and
the records of measurements 0-09b and C-25 step 1, and adds one read-only measurement (§7). Labels: *measured*
(host, date, command), *derived* (formula and inputs), *estimate* (with its basis).

## 0. The problem in ten lines

1. Domains are a property of a **stage**: `l1_domains` is partitioned by (`library`, `stage`), written once by the
   stage's `domains-<lib>` step, before any release exists (`stage.py` step list; `domains.py`).
2. Clusters are a property of a **release**: `cluster_member` and `cluster_member_aln` are partitioned by `release`
   and recomputed from scratch for each one (DEC-24).
3. DEC-32's rule needs clusters to pick what to scan, and alignments to place inferred envelopes.
4. If those are the release's, a stage's annotation depends on a release: the scan cannot run until a release
   including the stage is clustered, and inferred rows change whenever clusters move.
5. Inferred rows that depend on a release cannot live in a stage partition without one stage file serving two
   releases with different contents (AGENTS §3.10; card B-33's scope rule).
6. Representatives move: from `rN.1b` to `rN.3`, 107,870 of 2,053,075 proteins stopped being NR90 representatives
   and 35,639 members changed representative (*measured*, §7).
7. A protein that becomes a representative and was never scanned must be scanned then, and that late scan has no
   stage partition left to go into.
8. The alternative is to make the clustering that drives the scan a property of the stage too: cluster the
   stage's own proteins, scan its representatives, infer within it.
9. That keeps every rule of AGENTS §3 as it is, but gives up the redundancy that lies only across stages.
10. The choice is between those two, and C-25 step (3) measures the size of what the second gives up.

## 1. Where inferred rows live

Three options from the card, each paired with the scan scope it can serve (§2 defines 2a and 2c).

**(a) In `l1_domains`, marked.** Append one nullable column, `inferred_from` (`UID`): the representative's
`protein_uid`, null for a measured row. Inferred rows sit beside measured ones in
`l1/domains/library=<L>/stage=<S>/`, under the same key (`protein_uid`, `model_acc`, `env_start`, `library`,
`stage`) and sort key. Under LY-09 this is an allowed append; files written before it read as null (measured)
through `union_by_name` (`catalog._view_sql`).
- *Written when:* by the stage's domain step, once.
- *Strains:* nothing, **if** the clustering that decides inference is the stage's own (2c). With the release's
  clustering (2a) it breaks §3.10: release R2's different clusters would need different inferred rows in a file R1
  already lists. Writing a second part there is what `verify_manifest` reports as "not in manifest" for R1.
- *Tokens and postings:* read as today (`tokens._domain_files`, `postings._family_run`), so they include inferred
  rows with no code change. One change is needed: one source per protein (§1.1).
- *"Measured only":* `WHERE inferred_from IS NULL` on the same files. It prunes nothing more but reads nothing
  more either (the files are sorted by `protein_uid`).
- *Storage:* once per stage. At 10 Gbp, 3,600,029 inferred rows (*measured*, `rN.3`, §7) beside 16,749,026 measured
  ones, +21%. At about 70 bytes per row, about 250 MB (*estimate*: 51.9 B per measured row, *measured* on the P2
  file, 757,921,910 B / 14,608,733 rows, plus 16 B of uid that compresses poorly).
- *R2 when clusters move:* nothing changes. Release clusters do not drive domains under 2c.

**(b) A per-release table, `l1_domains_inferred`, partitioned (`library`, `release`).** It has the columns of
`l1_domains` plus `inferred_from`. `l1_domains` keeps measured rows only.
- *Partitioning:* it must be partitioned by `release`, not by stage *and* release. `catalog._scope` reads `stage`
  first, so a (`library`, `stage`, `release`) path would enter every release that lists the stage.
- *Written when:* by a new release step, `infer`, after `cluster` and before `tokens`.
- *Strains:* none of §3's rules for the inferred rows themselves (B-33 rule (b): a `release=` partition is that
  release's only). The strain is on the measured side (§2, 2a): late scans have no home.
- *Tokens and postings:* both must read a second source (the release's inferred files) beside the stages' domain
  files. The two paths would be duplicated in `tokens.run` and `postings._family_run`.
- *"Measured only":* free, since the inferred table is not read.
- *Storage:* one copy per release, about 250 MB per 10 Gbp (*estimate*, as above). At 95 M proteins and a release a
  month, it grows with every release.
- *R2 when clusters move:* R2 writes its own inferred rows from its own alignments. That costs no scan, only a
  CIGAR walk over about 2.2 M members per 10 Gbp (*measured* count, 2,188,469 in `rN.3`). It also needs the late
  scans of §4, which (b) does not house.

**(c) No stored rows: a view.** The catalog joins the representative's measured rows through `member_cluster` and
`cluster_member_aln` at query time.
- *Lookup:* a point lookup prunes. `member_cluster` is sorted by (`level`, `member_uid`), giving the
  representative; then `cluster_member_aln`, sorted by (`level`, `rep_uid`, `member_uid`), gives the alignment;
  then `l1_domains`, sorted by `protein_uid`. That is three joined keys, each pruned only with the join filter
  `open_catalog` raises (AGENTS §4.4, to be checked with `profile_query`).
- *Mapping:* the envelope mapping is a CIGAR walk, which is a Python UDF in the catalog.
- *Bulk readers:* `tokens`, `postings`, the search pack and the exports need every inferred row. Each would run the
  same join and walk at build time, so the logic exists in four places, or in one place that is (b) without the
  file.
- *"Measured only":* free.
- *Storage:* none.
- *R2:* recomputed. The late scans (§4) are still needed and still have no home.

### 1.1 One source per protein (needed by every option)

A protein can be listed in two stages: 13,302 of `rN.3`'s 10,620,684 proteins are (*measured*, §7). Today both
stages scan it and the rows agree (0-09b), so `DISTINCT` in the readers hides the duplicate. With inference, one
stage may infer a protein that another stage measured, or two stages may infer it from different
representatives, with envelopes a residue apart. `postings._family_run` would then write two family postings for
one domain.

The rule I propose for every reader of a release:

- a protein's rows are its measured rows if any stage of the release holds measured rows for it;
- otherwise they are the inferred rows of the stage that comes first in the release's `release_def` order.

The second clause depends on declared order, not on build order. The fact table and the catalog view stay
unchanged; the rule is applied where domains are read: `tokens`, `postings`, and a `domains_now` view used by
`domains_of` and the exports.

## 2. What gets scanned and when

**2a. Release first.** The release is clustered before its new stage is scanned. The scan becomes a stage step that
takes a release: `domains-<L> --release R`, where R lists the stage. Its plan is the stage's proteins that are R's
NR90 representatives or length-mismatched members, minus every protein already measured in a stage of R.
- *Step order:* `inputs, ingest, genecall, features, igr, tapes, cluster (R), domains-<L> (S, from R), infer-<L> (R),
  tokens, postings, dim, manifest, catalog`.
- *Shards:* cut as today (`domains._cut`, by residues, `protein_uid` order) from that list. The plan names R's
  `member_cluster` files and their sha256.
- *A stage added later:* it waits for a release that contains it. Its members of old, already measured
  representatives are inferred, which is the saving 2a has over 2c. Old proteins whose role changed must be scanned
  (§4). Those proteins belong to other stages, so their measured rows would be written into the new stage's
  partition: a stage file holding rows for proteins not in its `protein_seq`.
- *A release over existing stages only* (as `rN.3` was over `rN.1`'s stages): its late scans have no new stage to
  go into. It needs a new partition kind (a scan "round" owned by the release that made it, and readable by later
  releases, which is a new scope rule beside B-33's), or it scans them again in every release.
- *Order:* the measured set then depends on which releases were built before. Two roots with the same stages but
  different release histories hold different measured and inferred rows for the same release, so a release stops
  being a function of its stages (AGENTS §4.3, card B-25).
- *Precondition:* under 2a a stage cannot be scanned before the whole release is clustered. At S2C scale that is
  about 36 h per release (DEC-24, *derived* there) on the critical path of every new stage.

**2b. Copy, never rescan (no inference).** A stage's scan plan leaves out every protein already measured, with the
same HMM file sha256, cutoff and tool versions, in an earlier stage of the root. The merge copies those rows into
the stage's file.
- *Content:* at fixed Z a protein's hits do not depend on the other proteins (MR-12; `domains.py` docstring), so a
  copied row is byte-identical to the row a scan would write. The stage's `l1_domains` file is then the same
  whichever stages were built first.
- *Records:* only the split between scanned and copied depends on build order, and it goes in the scan plan and
  shard records (`domain_scan`, never in a manifest), not in the ledger.
- *Error:* no error at all.
- *Saving:* the fraction of a stage's distinct proteins identical to one already scanned. For one-genome-per-species
  stages it is tiny: 13,302 proteins in two stages at 10 Gbp, 0.13% (*measured*). For phase 2 against phase 1 (the
  same species' best genome is in P1 to P6) it is unmeasured and probably large (§6).
- *Step order:* unchanged. 2b is a change inside `domains.plan` and `domains.merge`.

**2c. Stage-local representatives (with 2b).** A new stage step, `represent`, clusters the stage's own distinct
proteins at NR90 with the release clustering's MMseqs2 parameters and module calls (`l2/cluster.py` `LEVELS`,
`COMMON`, the NR90 calls of `calls()`, in parts under B-41's rules). It then aligns each member to its
representative as `cluster.nr90.members` does.
- *Output:* a new table, `stage_rep`, partitioned (`rep_version`, `stage`). `rep_version` is `nr90-<params_hash>`,
  so C-19's owner rule spares releases built before it. It has one row per protein: `member_uid`, `rep_uid`,
  alignment columns as `cluster_member_aln`, and `role` (`rep`, `length`, `inferred`).
- *Scan set:* the domain step's plan is the `rep` and `length` proteins minus those copied (2b), cut into shards as
  today.
- *Merge:* writes scanned, copied and inferred rows into the stage's one file, sorted as today. The inferred rows
  are each `inferred` member's representative rows mapped through its `stage_rep` alignment (§3).
- *Step order:* `inputs, ingest, genecall, features, igr, tapes, represent, domains-<L>, cluster, tokens, postings,
  dim, manifest, catalog`. The release steps are unchanged.
- *A stage added later:* it is complete when built. Joining old clusters at the next release changes its NR90 and
  NR50 tokens, never its domains.
- *Order:* the stage's file is a function of the stage's proteins and the library alone. Clustering is
  deterministic (`--shuffle 0`, explicit threads, B-41's equality), and copies equal scans.
- *Cost:* one NR90 clustering per stage, about 990 CPU-minutes at 10.6 M proteins (*derived*: the sum of the
  `cluster.plan`, `cluster.nr90.*` and `cluster.nr90.members` run records of `rN.3`, 992.2 CPU-minutes, *measured*
  in `runs/stage=rN.3/`), against about 22,100 for that stage's full Pfam scan (*derived*: 2,776.8 M residues of
  P2 × 7.96).
- *Loss:* redundancy that lies only across stages. A phase-2 protein whose only NR90 relative is in another stage
  (typically the phase-1 genome of its species) is scanned. §7 measures the case at hand: S1 clustered alone scans
  45.12% of its residues, against 45.08% when clustered with P1 and 44.83% with P1 and P2, a loss of 0.3 points.
  For P1 (one genome per species) the loss is large: 94.3% alone with S1 against 83.1% once P2 joins, because
  131,072 of P1's 149,381 inferred members have their representative in P2 (*measured*). Phase 2 is between these.

**2d (an extension of 2c, not proposed now).** Before scanning, search the stage's representatives against
earlier stages' measured proteins (`mmseqs search`, at 90% identity and 80% member coverage). A representative
that hits one is inferred from it. This recovers 2a's cross-stage saving. The cost is that the stage's content then
depends on which stages existed when it was built, as DEC-26 accepted for contigs. Worth it only if C-25 step (3)
finds 2c's cross-stage loss large (Q7).

## 3. Mapping an envelope through the stored alignment

The alignment is `cluster_member_aln` under 2a and `stage_rep` under 2c, with the same columns. q* are the member's
coordinates and t* the representative's, 0-based half-open. The CIGAR has query = representative: `M` consumes
both, `I` the representative only, `D` the member only. 0-09b checked this convention on 5,000 of 5,000 rows per
level, and every non-representative CIGAR spans exactly its stored ends. Strands do not apply.

1. **Walk** the CIGAR from (`tstart`, `qstart`) and build the map from representative residue to member residue
   for `M` columns only.
2. **Start:** `env_start` maps to the member position of the first `M`-column representative residue at or after
   `env_start`.
3. **End:** `env_end` maps to one past the member position of the last `M`-column residue before `env_end`.
4. **`ali_start`, `ali_end`:** the same way. An envelope end inside an `I` run (representative-only residues) moves
   inward to the next aligned residue. A `D` run (member-only residues) inside the envelope is included.
5. **Model coordinates:** `hmm_from` and `hmm_to` are the representative's. They are model coordinates and are not
   mapped.
6. **Domains partly outside the aligned span.** Let f be the fraction of the representative's envelope that lies
   in [`tstart`, `tend`). On `rN.3` (*measured*, §7, 3,600,029 inferred rows):

   | f | rows | member's own scan has the model |
   |---|---|---|
   | 1 | 3,525,802 (97.93%) | 98.46% |
   | [0.8, 1) | 65,488 (1.82%) | 96.53% |
   | [0.5, 0.8) | 5,164 | 71.86% |
   | (0, 0.5) | 1,807 | 25.12% |
   | 0 | 1,768 | 15.55% |

   Proposed rule: a row with f ≥ 0.8 is inferred, mapped with its clipped end. If any of a member's rows has
   f < 0.8, the member is **scanned** (role `length` → `unmapped`), never written with a guessed envelope. At most
   8,739 members per 10 Gbp (*measured*: rows below 0.8, so members ≤ rows), 0.40% of inferred members. That is at
   most 2.8 M more residues, 0.1% of the release's scan (*derived*: 8,739 × 323.6 residues, P2's mean length). The
   alternative (Q4) is to drop such a domain; at these rates, dropping removes about 4,400 true pairs and about
   4,300 false ones (*derived* from the table).
7. **Length rule** (DEC-32): a member with |L_m − L_r| > 0.2 L_r is scanned, decided when the plan is made.
   Lengths are `protein_seq.length`, which equals the stripped sequence length on every row (0-09b).
8. **Scores:** an inferred row carries the representative's `seq_bitscore`, `seq_evalue_zref`, `dom_bitscore`,
   `dom_ievalue_zref`, `cutoff`, `z` and `dom_z` unchanged, marked by `inferred_from`. They are the
   representative's search results, never re-scored or estimated for the member, and the table doc must say so.
9. **No domains:** a member whose representative has no domain gets no rows. "Inferred to have none" and
   "scanned, none" are told apart by `stage_rep.role` and by the ledger, never by the absence of rows.

## 4. Stability between releases

- **Measured rows are kept forever and never rescanned.** Under 2c this follows from §3.10 (written once) and from
  2b (copied to a later stage, never rescanned). Under 2a it needs the late-scan home of §2.
- **Under 2c nothing moves when release clusters move.** The inference source is the stage's clustering, which no
  release changes. A protein's role is fixed when its stage is built.
- **Under 2a, from `rN.1b` to `rN.3`** (adding P2 to S1 and P1; *measured*, §7), over the 2,053,075 old proteins:

  | Transition at NR90 | proteins |
  |---|---|
  | representative → representative | 1,287,842 |
  | representative → member | 107,870 |
  | member → representative | 1,417 |
  | member → member, same representative | 620,307 |
  | member → member, other representative | 35,639 |

  - Scanned in `rN.1b` but inferable in `rN.3`: 103,000 proteins. They keep their measured rows.
  - Inferred in `rN.1b` but to be scanned in `rN.3`: **2,423 proteins, 663,720 residues**. These are the late scans,
    0.03% of the 2,269.5 M residues `rN.3` scans when `rN.1b`'s are kept.
  - Inferred rows to rewrite: every member whose representative changed, plus the late-scanned ones, about 38,000.
  - From `rN.1a` to `rN.1b` (adding P1 to S1): 140 late scans.

  So churn is cheap in CPU but structurally awkward: the late scans are few, and they still need a home (§2a).
- **New stage under 2a:** of P2's 8,567,609 new proteins, 131,812 were inferable from a protein `rN.1b` had
  measured, 1.5% (*measured*). For one-genome-per-species stages, measured representatives from old stages save
  little. Phase 2 against phase 1 is the case where they would save more (§6).

## 5. Verification and accounting

- **Ledger** of the domain step (`domains._write`), one row per reason, for object `library`:
  - `n_in` = the stage's proteins;
  - `no_hit` = measured with no domain;
  - `inferred` = members whose rows are inferred;
  - `inferred_none` = inferred, with a representative with no domain;
  - `unmapped_scanned` = members scanned under §3.6;
  - `length_scanned` = members scanned under the length rule.

  These counts are a function of the stage alone. The scanned-or-copied split goes in the plan and shard records,
  which are not in a manifest, because it depends on build order. The tool versions include MMseqs2's for
  `represent`.
- **`release verify`:**
  - (a) with 2c: no change. `l1_domains` keeps its paths. `stage_rep` is a stage table with a second partition key,
    so C-19's `_owner_key` spares releases built before it.
  - (b): `l1_domains_inferred` under `release=<R>` is R's only (B-33). The 2a late-scan rounds need a new scope
    rule.
  - (c): nothing stored.
- **Provenance of an inferred row:**
  - which clustering: `stage_rep` of its stage, at `rep_version`, whose params hash names the MMseqs2 parameters;
  - which representative: `inferred_from`;
  - which alignment: the `stage_rep` row of the member;
  - which scan: the representative's measured row in the same file, from a shard whose record names the plan, and
    the plan names the HMM sha256.

  Under (b) it is the same with "release R's clustering".
- **Tests** (AGENTS §5; each fails first):
  1. *Oracle on a fixture:* a constructed stage with planted cases: point substitutions of a representative, a
     fragment shorter by more than 20%, a fusion (representative with an extra domain), a representative with no
     domain, a domain at the representative's end outside the alignment, and a protein shared with an earlier
     stage. The oracle is a full `hmmsearch` scan (the binary in the image) of every protein, compared pair by pair.
     The envelope mapping is tested against a brute-force walk written separately: residue by residue over the
     two sequences, not over the CIGAR.
  2. *Determinism:* the same files at 1 and 4 threads, with the stage's input order shuffled, and with the two
     stages built in either order (2b's copies).
  3. *Memory:* the 100× synthetic store's `represent` and merge under the card's bound.
  4. *Real stage:* one phase-2 stage scanned in full and by representatives. Accept if the wrong-pair rate
     (false / reported, 0-09b's definition) is ≤ 1% (DEC-32: stop above 1%). Report missed / true, and the
     95th-percentile envelope difference (≤ 5 residues at NR90 proposed; 2 *measured* on `rN.1b`).

  The full scan of tests 1 and 4 is not wasted: those measured rows are the stage's result.
- **What C-25 step (3) must measure** on one fetched phase-2 stage, for 2c to be chosen:
  1. distinct proteins over genes, within the stage;
  2. the share of its distinct proteins identical to a protein of P1 to P6 (2b's saving);
  3. the hybrid's scanned residue fraction with the stage clustered alone (2c), and clustered with P1 to P6 (the 2a
     counterfactual: one clustering of those stages, run read-only into scratch);
  4. the wrong-pair and missed rates of stage-local inference against the full scan;
  5. §3.6's f distribution;
  6. the CPU-minutes of the stage's NR90 clustering.

  2c is chosen if the gap in (3) is under 10 points of the stage's residues (Q7) and (4) is ≤ 1%. A warning from
  §7: S1, the phage set with many isolates per species, already sits at **0.98% to 1.00% wrong** (*measured*:
  0.983% clustered alone, 1.004% in `rN.3`). A strain stage may land close to the threshold.

## 6. Cost for S2C phase 2

Inputs:
- 51,597,124 GTDB proteins (gene calls) in the 13,485 remainder genomes (*measured*, metadata,
  `2026-11-12_s2c_selection.md`).
- 323.6 residues per protein (*measured* on P2: 2,776,784,358 / 8,580,713).
- 7.96 CPU-minutes per M residues (*measured*, C-10).

**Today's rule** (every stage's distinct proteins scanned) costs at most **2,215 CPU-hours** (*derived*:
51.6 M × 323.6 × 7.96 / 60 / 1e6, assuming no protein repeats within a stage; DEC-32's 2,500 is the per-Gbp
version of the same estimate). Call this T.

Unknowns, each an *estimate* until C-25 step (3) replaces it:
- **x**, the share of a phase-2 stage's distinct proteins identical to an already-scanned protein. Assumed 0.3,
  since strains of one species share many identical core proteins with the phase-1 genome. Replaced by (3)(2).
- **h**, the hybrid's scanned residue fraction for a strain stage clustered alone. Assumed 0.45: S1's *measured*
  45.1%, the one redundant set we have. Replaced by (3)(3).
- **g**, the extra points 2a gains from cross-stage relatives. Assumed 0.05: S1 *measured* 0.003, P1 0.11; a
  phase-2 stage cut by species is closer to S1. Replaced by (3)(3).

Copies and inference are taken as independent (*estimate*).

| Option | Scanned share of T (*derived*) | CPU-hours scanned | Saved | Extra |
|---|---|---|---|---|
| today | 1 | 2,215 | 0 | none |
| 2b alone | 1 − x = 0.70 | 1,550 | 665 | none; no error |
| **2c + 2b** | h (1 − x) = 0.315 | 700 | 1,515 | NR90 clustering of each phase-2 stage, 50 to 80 CPU-h (*derived*: 992 CPU-minutes per 10.6 M proteins, scaled by the B-36 exponent 1.244 to 5 or 6 stages of 51.6 M × 0.7 to 1.0 proteins) |
| 2a | (h − g)(1 − x) = 0.28 | 620 | 1,595 | late scans about 0.03% (*measured* churn ratio); the release clustering before every new stage's scan |
| 2a without copies | h − g = 0.40 | 885 | 1,330 | as 2a |

The difference between 2a and 2c is about 80 CPU-hours of scan plus 50 to 80 of clustering, roughly 6% to 7% of
T. That is the price of keeping domains a stage property, at these assumptions. It scales with g; if C-25 (3)
finds g = 0.15, the scan difference becomes about 230 CPU-hours (*derived*: 0.15 × 0.7 × 2,215).

## 7. Measurement made for this note

*Measured* on wilkcomp02, 2026-11-12 23:32 to 23:36 CST. Command:

```
systemd-run --user --unit b47-churn --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
  -p StandardOutput=append:/tmp/if2453/b47/logs/b47-churn.log -p StandardError=append:/tmp/if2453/b47/logs/b47-churn.log \
  nice -n 19 singularity exec --cleanenv --bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro --bind /tmp/if2453/b47 \
  /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif python /tmp/if2453/b47/bin/timed.py \
  /tmp/if2453/b47/runs/b47-churn.json -- python /tmp/if2453/b47/bin/b47_churn.py /tmp/if2453/r10/root \
  /tmp/if2453/b47/work/churn.json --memory 5GB --tmp /tmp/if2453/b47/work
```

- **Script:** `scripts/stress/annotation_cost/b47_churn.py`. It runs in DuckDB, one thread.
- **Inputs:** only files the manifests of `rN.1a`, `rN.1b` and `rN.3` list, each checked against its sha256. These
  are `protein_seq`, `cluster_member` and (for `rN.3`) `cluster_member_aln`, at level `nr90`, and the three
  `pfam-38.2` `l1_domains` files `rN.3` lists (S1 605,564, P1 1,552,653, P2 14,608,733 rows).
- **Cost:** 3.75 CPU-minutes, 240 s, peak RSS 4,681 MiB. The account was at about 45 threads before launch,
  under the night cap of 56. Nothing was written outside `/tmp/if2453/b47`.
- **Check:** `rN.1b`'s decision row reproduces 0-09b's exactly (right 2,024,291, false 7,429, missed 6,853,
  400,152,460 residues scanned).
- **Note on `rN.3`:** its clusters equal `rN.1`'s byte for byte (equal sha256 of `cluster_member` and
  `cluster_member_aln` in both manifests).

Hybrid rule (DEC-32) at NR90; wrong = false / reported, missed = missed / true:

| Release (stages) | Scope | Residues scanned | Wrong | Missed | Inferred members | Inferred rows | Inferred members with representative in another stage |
|---|---|---|---|---|---|---|---|
| `rN.1a` (S1) | S1, stage-local | 45.12% | 0.983% | 0.907% | 545,002 | 352,798 | 0 |
| `rN.1b` (S1, P1) | S1 | 45.08% | 0.984% | 0.909% | 545,317 | 353,261 | 444 |
| `rN.1b` | P1 | 94.27% | 0.116% | 0.107% | 51,339 | 101,719 | 558 |
| `rN.1b` | all | 70.66% | 0.366% | 0.337% | 596,590 | 454,925 | |
| `rN.3` (S1, P1, P2; 10 Gbp) | S1 | 44.83% | 1.004% | 0.917% | 547,871 | 356,871 | 4,389 |
| `rN.3` | P1 | 83.11% | 0.340% | 0.346% | 149,381 | 294,730 | 131,072 |
| `rN.3` | P2 | 81.74% | 0.363% | 0.363% | 1,499,477 | 2,960,886 | 128,908 |
| `rN.3` | all | **78.89%** | **0.385%** | **0.382%** | 2,188,469 | 3,600,029 | |

New against 0-09b:
- the error at 10 Gbp, which 0-09b could not measure because P2's scan was not yet published;
- the stage-local case (`rN.1a`);
- the churn of §4;
- the envelope coverage of §3.6.

## 8. Recommendation

**2c with 2b, rows in `l1_domains` (option a); release steps unchanged.**

Reasons:
- Domains stay a property of the stage. Every rule of AGENTS §3 holds as written: §3.7 (one `TableSpec` each, one
  new), §3.8 (an appended nullable column), §3.9 and §3.10 (written once, in its stage), and B-33/C-19's scope.
- DEC-24 needs no amendment. A release stays a function of its stages, independent of the history of earlier
  releases, which 2a cannot claim (§2a, "Order").
- Measured rows are never rescanned, by 2b, and nothing ever needs a late scan.
- `tokens` and `postings` change only by the one-source rule (§1.1).
- The scan of a new stage does not wait for a 36-hour release clustering.
- The price is the cross-stage redundancy, 6% to 7% of phase 2's scan at the assumptions of §6, which C-25 step (3)
  measures. 2d recovers it later if it proves large, without changing the storage.

Risks:
1. *g large.* If phase 2's strains are spread across stages, or most of their relatives are the phase-1 genome,
   2c's loss is large. Mitigation: cut phase-2 stages by species (Q3); 2d.
2. *Wrong pairs near 1%.* The redundant set measured, S1, is at 0.98% to 1.00%. A strain stage could exceed DEC-32's
   threshold, under any option. Then DEC-32 says stop and return to the PI.
3. *Extra clustering per stage* (50 to 80 CPU-hours for phase 2). Card C-22 found the per-call estimates 25% to 36%
   low at 10.6 M proteins, so `represent` must use B-41's parts and C-22's refit.
4. *The one-source rule* is a reader change in three places. A reader missed would double family postings on
   proteins in two stages, 0.13% of proteins at 10 Gbp (*measured*). The test is a protein measured in one stage
   and inferred in another.
5. *A later full scan* of a stage scanned by representatives cannot rewrite its file. It is possible only as a new
   library version, or as a later card adding a "completion" partition and its scope rule. That is not proposed
   now.

## 9. Questions only the PI can answer

1. Does a stage's own NR90 clustering count as DEC-32's "NR90 cluster" (option 2c), rather than the release's?
   *yes / no*
2. Copy measured rows of a protein already scanned in an earlier stage, instead of scanning it again, in every new
   stage (2b)? *yes / no*
3. Cut S2C phase-2 stages so that all of a species' phase-2 genomes fall in one stage? *yes / no*
4. A domain whose envelope lies less than 80% inside the alignment: scan the member, or drop the domain? *scan /
   drop*
5. Do family tokens and postings include inferred rows by default? *yes / no* (DEC-32 implies yes)
6. On the strain stage, stop above 1% wrong pairs only, or above 1% wrong or 1% missed? *wrong / either*
7. Revisit (2d or 2a) if C-25 step (3) finds 2c scanning more than 10 points more residues than 2a? *yes / other
   threshold*
8. A protein inferred in two stages takes the rows of the stage first in `release_def` order? *yes / no*

## 10. Cards that follow from the recommendation (text for `TASKS.md`)

```
### B-47a Stage-local NR90 representatives: the `represent` step
**M** · Depends: B-47 (PI answers Q1, Q3), C-25 step (3), B-41, C-22 · Design: `planning/design/B-47_inferred_domains.md` §2c, §5
Files: `src/biobox/l1/represent.py` (new), `src/biobox/l2/cluster.py` (expose the NR90 call list and member alignment for reuse; no change to release clustering), `src/biobox/registry/schemas.py` (`stage_rep`, partitioned (`rep_version`, `stage`)), the locks by `scripts/update_schema_lock.py`, `src/biobox/stage.py` (step `represent` before `domains-<lib>` when `--library ... --by-representatives`), `src/biobox/cli.py` (`stage represent`), `RUNBOOK.md`, tests
Do, test first: cluster the stage's distinct proteins at NR90 with the release clustering's parameters and module calls, in parts under B-41's rules and C-22's constants; align each member to its representative; write one `stage_rep` row per protein with its role (`rep`; `length` when |L_m − L_r| > 0.2 L_r; `inferred`). Byte-identical at 1 and 4 threads and with input order shuffled; on the fixture, the representatives equal `cluster.run`'s NR90 representatives over the same proteins (the release clustering as an independent path). A release built before the step still verifies (C-19). Memory test on the 100× store.
Accept: the tests; on the pod, `represent` on the C-25 (3) stage with its CPU-minutes against §6's derivation.

### B-47b Copy measured rows instead of scanning again
**M** · Depends: B-47 (Q2)
Files: `src/biobox/l1/domains.py` (plan leaves out proteins measured in an earlier stage with the same HMM sha256, cutoff and tool versions; merge copies their rows; the plan names the files copied from, with sha256), `src/biobox/stage.py` (passes the root's other stages, as `--seq-stages`), tests
Do, test first: two stages sharing proteins, built in either order, give byte-identical `l1_domains` files; a copied row equals the row of a full scan of that protein (`hmmsearch` binary as oracle); a source scanned with another HMM file or pyhmmer version is not copied. The scanned and copied counts go in the plan and shard records only.
Accept: the tests.

### B-47c Inferred rows: envelope mapping and the merge
**M** · Depends: B-47a, B-47b, B-47 (Q4)
Files: `src/biobox/l1/infer.py` (new: the CIGAR walk of §3), `src/biobox/l1/domains.py` (merge writes scanned, copied and inferred rows in one sorted stream, bounded memory), `src/biobox/registry/schemas.py` (`l1_domains.inferred_from`, appended, nullable; table doc says inferred rows carry the representative's scores), the locks by `scripts/update_schema_lock.py`, `src/biobox/ledger.py` reasons if needed, tests
Do, test first: §3's mapping against a brute-force residue walk written separately; §3.6's rule (f < 0.8 → the member is scanned, ledger `unmapped_scanned`); ledger rows of §5; on the planted fixture every inferred pair checked against a full `hmmsearch` scan; files written before the column read as measured.
Accept: the tests; the fixture's wrong-pair rate reported.

### B-47d Readers: one source per protein, and "measured only"
**M** · Depends: B-47c, B-47 (Q5, Q8)
Files: `src/biobox/l3/tokens.py`, `src/biobox/l3/postings.py`, `src/biobox/catalog.py` (`domains_now` view and `domains_of` with `inferred_from`; a `measured_domains` view), `src/biobox/search/export.py` (an `inferred_from` column, README line), the W-07 byte-identity test (migrated for the appended column only), tests
Do, test first: a protein measured in one stage and inferred in another gets the measured rows only, one family posting per domain; inferred in two stages, the first stage in `release_def` order; tokens and postings of a release without inferred rows are byte-identical to today's.
Accept: the tests.

### B-47e The first phase-2 stage by representatives, against its full scan
**M** · Pod: yes · Depends: B-47a to B-47d, C-25 steps (2) and (3)
Do: build the C-25 (3) stage with `--by-representatives` into a scratch root that reads the full scan's stage files without changing them; compare pair by pair with the full scan; report wrong / reported, missed / true, envelope differences (median, p95), the f distribution, residues scanned, CPU-minutes of `represent` and of the scan, against §6. Stop and report to the PI if wrong > 1% (DEC-32).
Accept: the measurement file.
```

## 11. Not determined here, and what settles each

1. x, h, g for a strain stage. Settled by C-25 step (3), measurements (1) to (3) of §5.
2. Wrong-pair rate on a strain stage. Settled by C-25 step (3)(4); S1's 0.98% to 1.00% is the nearest evidence.
3. The cost of `represent` at a phase-2 stage's size. Settled by B-47a's pod run.
4. Whether linclust alone (the NR90 calls `00` to `09`, 73.0 CPU-minutes of the 992, *derived* from `rN.3`'s run
   records) would do for stage-local representatives, at a tenth of the cost. Its clusters differ from the
   cascade's, and that difference is not measured. A follow-up for B-47e, not a reason to wait.

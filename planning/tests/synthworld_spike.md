# `synthworld` S−1 fit spike: ZOMBI2 0.54.0

Date: 2026-10-29. Answers the open questions in `synthworld_spec.md` §1.4 and §8.1 by running ZOMBI2. All runs, scripts and analysis JSON are in the session scratchpad (`…/scratchpad/zombi2_spike/{env,runs,scripts,analysis}`). They are not copied into the repo.

**Verdict: GO, with conditions.** ZOMBI2 0.54.0 works as the backbone. It is deterministic, fast (13 s at 4,000 genes per genome) and its copy-to-leaf mapping is exact. Group integrity holds without exception, but only with `fates="group"`. Seven spec statements need correcting (§6).

---

## 1. Setup

| item | value |
|---|---|
| host | Apple M1 Pro, 10 cores, 32 GB, macOS (Darwin 24.6.0) |
| Python | CPython 3.12.12, venv made by uv 0.10.7 |
| packages | `zombi2==0.54.0`, `numpy==2.5.4`, `tqdm==4.70.1` (the full resolved set) |
| isolation | Separate venv in the scratchpad. Our own scripts sit in `scripts/` and each run writes to a new empty folder in `runs/`. Every Python call used `-I`. Nothing went into the repo `.venv`. |

```
uv venv <scratch>/zombi2_spike/env --python 3.12
uv pip install --python <scratch>/zombi2_spike/env/bin/python zombi2==0.54.0
/usr/bin/time -l env/bin/python -I scripts/run_ordered.py runs/<new> --genes N --groups --group-reps 3 \
    --group-loss 0.03 --group-transfer 0.15 --chromosomes 2 --topology circular,linear --translocation 0.01
env/bin/python -I scripts/analyze.py runs/<run> analysis/<run>.json      # Q1, Q2, schema
env/bin/python -I scripts/attribute.py runs/<run>                         # cause of each broken unit
env/bin/python -I scripts/digest.py runs/<a> runs/<b> [--exclude NAME]    # Q3
env/bin/python -I env/bin/zombi2 species|genomes runs/<new> …             # CLI runs (Q3)
```

**Species tree.** All runs use the same tree: `species.simulate_species_tree(birth=1.0, death=0.3, n_extant=46, seed=20261029)`. It has 46 extant and 9 extinct tips, 109 nodes, height 5.02 and total branch length 49.1.

**Genome rates.** These are spike rates, because §5.1 does not give DTL rates. All rates are per copy:
- duplication 0.02, transfer 0.05, loss 0.12;
- inversion 0.05, transposition 0.02, translocation 0.01;
- origination 0.05 × G per lineage, where G is `initial_families`.

Extents are geometric means: duplication 2, loss 2, transfer 3, inversion 5, transposition 3, translocation 3. Other settings: `inversion_probability` 0.5, `transfer_to="uniform"`, `replacement=0`, two replicons (one circular, one linear).

**Groups.** Each group run declares 30 groups of 2–6 families. They are 3 replicates of the full dial grid:
- `together` ∈ {1, 0.7, 0.3, 0} × `fates` ∈ {family, group}, with `adjacent=True`;
- plus two `adjacent=False` groups (`together` 1 and 0).

Member families have their own loss rate (0.03) and transfer rate (0.15), so groups survive and travel. A second Q1 run, `dupext1`, uses extent 1 for duplication, loss and transfer.

**Analysis rules.**
- A *unit* is one row set of `groups.tsv` sharing a `unit` id: one copy of a group in one genome.
- A unit counts as **intact** when all four hold: its members present in an extant genome are on one chromosome, form a contiguous run (ring wrap allowed), share one strand, and appear in declared order (reversed on −).
- Units with fewer than two surviving members are not scored.

---

## 2. Q1: do group members stay co-oriented and adjacent?

**Answer: yes for `fates="group"`, without exception. No for `fates="family"`: about 10–14% of multi-member units break.** The break comes from tandem duplication, not from rearrangements.

### 2.1 Integrity

| run | fates | units with ≥2 members | intact | co-oriented | on − strand (flipped since root) |
|---|---|---|---|---|---|
| `q1_groups_2000_a` (5,193 inversions, 2,039 transpositions, 1,006 translocations) | group | 1,905 | **1,905 (100%)** | 1,905 (100%) | 953 (50%) |
| | family | 1,185 | 1,068 (90.1%) | 1,132 (95.5%) | 482 |
| `q1_groups_2000_dupext1` (5,405 / 2,148 / 1,062) | group | 995 | **995 (100%)** | 995 (100%) | 502 (50%) |
| | family | 1,404 | 1,205 (85.8%) | 1,334 (95.0%) | 645 |
| `q4_groups_4000` (9,175 / 3,659 / 1,835) | group | 1,625 | **1,625 (100%)** | 1,625 (100%) | 909 (56%) |
| | family | 1,398 | 1,199 (85.8%) | 1,322 (94.6%) | 597 |

Every gene in `initial_genome.tsv` is on + (2,120 of 2,120 in run a). Half the intact `fates="group"` units sit on −, so inversions and inverted moves flipped them as whole blocks thousands of times. They were never cut.

### 2.2 Why `fates="family"` units break

`scripts/attribute.py` walks each broken unit's ancestry up the tree. For each one it finds the first branch on which the block was broken and lists the events on that branch that touched the members' ancestors.

| run | broken units | a member was duplicated on that branch | transfer only | no member event (root stem) |
|---|---|---|---|---|
| a | 117 | 106 | 7 | 4 |
| dupext1 | 199 | 182 | 10 | 7 |

Two traced examples:
- **Run dupext1, lineage `n103`, chromosome 206.** Positions 126–128 hold `g251590` (family 2000, grp0 unit 31), then `g251591` (family 2000, *not* a member), then `g251592` (family 2001, unit 31). The non-member tandem copy of the first member sits *inside* the unit.
  - The `group()` docstring says "a tandem copy of a member is placed beside the group, not inside it". This contradicts it, even at duplication extent 1.
- **Run a, lineage `n101`, unit 109 (grp12).** The four members read strands −1, +1, +1, +1. The block first broke on branch `n21`, where duplication event 9739 re-minted the unit's copy of family 2048.

Once a unit is no longer contiguous, the engine stops protecting it. Later translocations then carry members onto the other replicon (31 such cases in run a).

### 2.3 Dials

Each transfer event is read through the `event` id in `genome_events.tsv`, and classified by how many group members it moved.

**`together`.** Share of member transfers that moved the whole group, for `fates="group"`:

| `together` | extent 1 | extent mean 3 (run a) | extent mean 3 (4,000-gene run) |
|---|---|---|---|
| 1.0 | 87/87 (100%) | 224/224 (100%) | 279/279 (100%) |
| 0.7 | 160/217 (74%) | 277/468 (59%) | 312/519 (60%) |
| 0.3 | 22/73 (30%) | 36/114 (32%) | 14/67 (21%) |
| 0.0 | 0/81 (0%) | 8/303 (3%) | 4/117 (3%) |

The dial behaves as documented at extent 1. With longer extents the realized share drifts from the dial in both directions, so truth must count events, not trust the dial.

**`adjacent=False`.** These groups write **no units** to `groups.tsv`. Their co-transfer is visible only through `event` ids. With `together=1`:
- extent 1: 348/348 transfers moved the whole group;
- extent mean 3: only 148/410 did (36%; 44% at 4,000 genes). A run that starts on a neighbour and covers one member carries that member alone.

**`fates`.**
- `fates="group"` with `together=1`: all 88 member losses took the whole unit (run a), and duplications copy the whole unit.
- `fates="family"`: losses are mostly of single members (65 single-member vs 18 whole-group in the same cell).

**Refusals** (`scripts/check_refusals.py`):
- `replacement` is refused with groups at *any* value above 0, not only `True`: 1.0 and 0.5 both raise `ValueError: replacement above 0 is not built for a run with groups`. Without groups, 0.5 is accepted.
- An `adjacent=True` group cannot have members with `origin=`: `ValueError … an adjacent group starts in the initial genome, whole, so its members take no origin` (`scripts/check_origin.py`). `adjacent=False` accepts it.

**Root layout.** Declared (named) families are laid out as **one consecutive block at the end of each initial chromosome**, in declaration order, dealt round-robin across replicons. In run a they occupy positions 1000–1059 on both replicons. This holds for `adjacent=False` groups too.
- So every module starts next to other modules, and the members of "independent" `adjacent=False` groups start next to each other.
- Extant genomes still hold 431 neighbouring pairs inside `adjacent=False` groups, and 304 distinct cross-group neighbour pairs. Examples: grp18–grp28 sit next to each other 75 times, grp18–grp19 66 times.

---

## 3. Q2: do gene copies map one-to-one to gene-tree leaves?

**Yes, exactly, against the *extant* gene trees.** The join key is `f"{lineage}_{copy}"`, because `gene_order.tsv` writes `copy` as `g<id>`.

| run | extant `gene_order` rows | extant-tree leaves | rows without a leaf | leaves without a row | duplicate keys | family mismatch |
|---|---|---|---|---|---|---|
| probe (200 genes) | 8,654 | 8,654 | 0 | 0 | 0 | 0 |
| a (2,000) | 102,951 | 102,951 | 0 | 0 | 0 | 0 |
| dupext1 (2,000) | 105,656 | 105,656 | 0 | 0 | 0 | 0 |
| 4,000 | 183,292 | 183,292 | 0 | 0 | 0 | 0 |

**Trap: the leaf prefix does not mean extant.** Complete gene trees contain leaves named `n<extant tip>_g<id>` that are **lost copies**: 8,801 in run a and 17,104 at 4,000 genes.
- Example: `n100_g448668` is born at a speciation and then lost on branch `n100` (`genome_events.tsv`, lines 265131 and 268855).
- So the adapter must take extant copies from `gene_order.tsv` or the `extant` column, never from the `n` prefix in `complete`.
- Every extant row does appear in its complete tree (0 missing).

Other facts the adapter needs:
- **Copy ids are re-minted at every speciation.** For example, `n33_g7562` becomes `n83_g18224` and `n84_g18340`. A copy id names one node's gene, not a lineage.
- **Internal node labels are `<event>_n<species>`**, for example `speciation_n69` and `transfer_n36`, as documented.
- 7,092 families were born in run a; 5,696 have an extant tree. The `extant` cell is empty for the rest.

---

## 4. Q3: determinism; serial vs parallel

**The serial engine is byte-identical across repeated runs. The ordered resolution has no parallel engine at all.** Digests below are sha256: run digests are given in full, file digests as the first 16 hex characters.

| comparison | result |
|---|---|
| ordered, Python API with groups: run a vs run b (same args) | All 18 ZOMBI2 files identical. Run digest (all files except our `spike_info.json`) `b3ac7652509f9d83052ae8aa07281636afcbefd11fefa251554dac7dfa9d130c` for both. Files: `gene_order` `106e4c06929aa2f2`, `gene_trees` `e1a5f390a3b53fb3`, `genome_events` `afdcc12633d8cda4`, `groups` `06d4dc750a9ef6b2`, `rearrangement_events` `fcad8d15ec799c71`, `species_complete.nwk` `b477c09e5194d94d`. |
| ordered, CLI without groups: `q3_cli_a` vs `q3_cli_b` | All data files identical (`gene_order` `73193b0769fbf642`). **`run.zombi2`, `genomes/genomes.log` and `species/species.log` differ** between the two runs. They embed wall-clock timestamps (`seed 20261029 · 2026-10-29 18:13:31`) and the absolute run path. |
| ordered, CLI vs Python API (same params) | `genomes/` byte-identical, apart from `species_complete.nwk`, which the API's `write()` adds there. Digest of the CLI `genomes/` without its log: `bcb5c3e4cf28823b6a5689ca246f9268feac400ffd02513e939b4236296c90fe`. |
| ordered, parallel | Python: `TypeError: simulate_genomes_ordered() got an unexpected keyword argument 'parallel'`. CLI: `--parallel applies to --resolution family only … the ordered resolution couples families by position … so it has no per-family engine`. |
| family resolution: serial ×2 vs `parallel=2` vs `parallel=4` | Serial a = serial b: `1512d6e88c5562d8e602c6808e7d0c14cef564e0a5d56264b9fe115d25dd1e3d`. `parallel=2` = `parallel=4`: `9592424ec618159d82b72b4fd89b05018ab68d9229dbbdcb61f00ff18c244418`. Serial and parallel differ in 6 files (129,249 vs 128,407 events). The species tree is identical. So the parallel engine is worker-count-invariant but is a different draw from the serial one, as documented. |

macOS vs Linux was **not tested**; only one host was available. That check stays a CI item.

---

## 5. Q4: runtime and memory (46 tips, ordered, single process)

Configuration as in §1, with 30 groups and two replicons. "Wall" is the whole process from `/usr/bin/time -l`: interpreter start, imports, species tree, simulation and writing every output with `compact=True`.

| initial families | gene events | rearrangements | extant genes/genome, mean (min–max) | simulate s | write s | wall s | max RSS MiB | peak footprint MiB |
|---|---|---|---|---|---|---|---|---|
| 1,000 | 103,559 | 4,984 | 1,454 (1,357–1,581) | 0.95 | 2.65 | 4.16 | 290 | 211 |
| 2,000 (run a) | 166,086 | 8,238 | 2,238 (2,119–2,398) | 1.50 | 3.98 | 6.56 | 427 | 322 |
| 2,000 (run b, repeat) | same | same | same | 1.52 | 4.00 | 6.58 | 432 | 321 |
| 2,000 (extents 1) | 148,909 | 8,615 | 2,297 (2,206–2,366) | 1.55 | 3.84 | 6.18 | 399 | 306 |
| 4,000 | 300,207 | 14,669 | 3,985 (3,756–4,228) | 2.80 | 8.22 | 12.94 | 681 | 511 |
| 2,000, CLI, no groups, one replicon | — | — | 5,520 families | — | — | 3.51 | 276 | 195 |

§1.4 item 1 passes: 6.6 s against the 60 s limit. Time and memory grow roughly linearly with genome size. Writing the TSVs takes about 70% of the time.

At 4,000 genes the run directory is 59 MB: `genome_events` 22 MB, `gene_trees.tsv` 20 MB, `gene_order` 14 MB.

The 5,000-tip × 3,000-gene stage world (§5.2) was not run. Extrapolating linearly in branch length × genes gives roughly 100× the 4,000-gene run in events, and memory in tens of GB unless `stream_to=` is used. Treat that as an estimate to measure on the stage host.

---

## 6. Output files the adapter reads (§2.2)

All five exist. Headers as written by 0.54.0:

| file | columns | difference from the docs or spec |
|---|---|---|
| `species/species_complete.nwk` (also `species_extant.nwk`, `species_fates.tsv`) | Newick, nodes `n<id>`; extinct tips `e<id>` | none |
| `genomes/gene_order.tsv` | `lineage chromosome topology position strand family copy` | `copy` is a string `g<id>`; `family` is an integer id |
| `genomes/gene_trees.tsv` (from `write(dir, compact=True)`) | `family complete extant` | In Python, `compact=` on `simulate_genomes_ordered` needs `stream_to=`; pass it to `write()` instead. |
| `genomes/genome_events.tsv` | `time kind family parents children chromosome start length dest_chromosome dest_position event` | **11 columns**. The docs header omits `event`. Transfer `kind` is `transfer_additive` or `transfer_replacing`, never `transfer`. |
| `genomes/groups.tsv` | `group unit family gene` | `gene` is a bare integer (no `g`); `family` is an id. Names map through `genome_summary.json` → `families.names`. No rows for `adjacent=False` groups. |
| also present | `rearrangement_events.tsv`: `time kind lineage chromosome start length dest_chromosome dest_position flipped cuts`; `chromosome_events.tsv`; `initial_genome.tsv`; `profiles.tsv`; `genome_summary.json`; `links.tsv`; `family_multipliers.tsv`; `lineage_multipliers.tsv` | `rearrangement_events.tsv` has a `cuts` column the docs do not list |

---

## 7. Spec corrections needed (`synthworld_spec.md`)

1. **§1.3 reason 2, §2.1 S2, §3.3 `truth_module`, §7.3: group adjacency is guaranteed only for `fates="group"`.** With `fates="family"`, 10–14% of multi-member units are broken in extant genomes (§2.2), so the §7.3 invariant "members with `adjacent=True` are contiguous in every extant carrier" fails. Choose one:
   - plant all adjacency-scored modules with `fates="group"`, and reinterpret the §5.3 "per-family loss" strength axis;
   - or keep `fates="family"`, compute `realized_loci` from `gene_order.tsv` rather than assuming contiguity, and file an upstream issue: tandem copies land inside the unit, against the docstring.
2. **§5.3 "replacement probability" knob and §4.1 N0: `replacement` must be 0 in any run with groups.** Vehicles and Cas operons are groups, so only N0-strict could use replacing HGT. Either drop the knob or note that N0-strict then differs from N0-full in more than the planted groups.
3. **§5.3 "origin depth ∈ {young, old}": not available for `adjacent=True` groups.** They always start whole in the root genome. Young modules need another mechanism, such as clade-restricted presence through per-family loss or transfer rules, or an upstream feature.
4. **§2.2 and §6.2 (`run.zombi2` and digests):**
   - Groups can be declared only from the Python API (no CLI flag or TOML key). The API does not write `run.zombi2`, which is built by the CLI from per-level `.log` records. `synthworld` must write its own ZOMBI2 provenance: version, numpy version, the full kwargs and the seed.
   - Where a `run.zombi2` or `.log` exists, it holds timestamps and absolute paths and is not byte-stable. Exclude it from `world_digest`, or normalize it first.
5. **§2.2 schema pins:**
   - `genome_events.tsv` has 11 columns including `event`;
   - transfer kinds are `transfer_additive` and `transfer_replacing`;
   - `rearrangement_events.tsv` has a `cuts` column;
   - `copy` is `g<id>` but `groups.tsv.gene` is a bare integer;
   - the leaf label is `<lineage>_<copy>`.
6. **§2.3 and §3.2 (`copy_label`):**
   - Complete gene trees label lost copies with the extant tip's `n` prefix. Extancy comes from `extant` trees or `gene_order.tsv` only.
   - Copy ids change at every speciation.
7. **§5.1 genome size "range 1,000–4,000 (tip-level origination/loss balance)": one root genome with uniform rates gives only about ±9% spread across tips.** Measured ranges: 2,119–2,398 and 3,756–4,228.
   - A 4× spread needs lineage-varying rates (`varying_among('lineages', Drift(...))`, which the ordered engine accepts; not tested here) or a different sizing approach.
   - The mean also drifts from target with the rate mix: −19% at 2,000 genes without groups, +32% at 1,000 genes with groups. The §1.4 item 5 ±30% band therefore needs fitted rates. §5.1 also needs actual DTL and rearrangement rates.
8. **§4.1, §4.2 and §2.4 confound: all declared families start in one consecutive block at the end of the root chromosomes.** Modules therefore start adjacent to other modules, and `adjacent=False` members start adjacent to each other. This plants inherited adjacency between planted families that the null does not have.
   - Mitigation to test at S0: a high transposition rate on the root stem (`changing_at` skyline), or upstream random placement.
   - Related: every root gene is on +, so strand switches come only from inversions. This matters for the §7.1 directon and strand-switch statistics.
9. **§1.4 item 4 and §6.1: resolved.** `simulate_genomes_ordered` has no `parallel` argument, so ordered runs are always serial. Parallel vs serial differences apply only at the family resolution. macOS vs Linux digest equality is still untested.
10. **§2.1 S2: "`together`" is a per-event probability for transfers that start on a member.** Realized co-transfer share differs at extent > 1 (§2.3: 59–60% at `together=0.7`; 36–44% for `adjacent=False`, `together=1`). `truth_module.realized_independent_transfers` must count `event` ids, never derive them from the dial.

---

## 8. Go / no-go

**GO for ZOMBI2 0.54.0 as the S2 backbone, under three conditions:**
- (a) planted adjacency modules use `fates="group"`, or truth is computed from realized order;
- (b) `replacement=0` in every world with groups;
- (c) root-layout scrambling is added and checked at S0.

| criterion | result |
|---|---|
| copy ↔ leaf mapping (Q2) | exact |
| serial determinism (Q3) | byte-identical data files |
| speed (Q4) | 13 s and 0.7 GB at 4,000 genes |
| group integrity under inversions, transpositions, translocations and transfers (Q1) | perfect with `fates="group"` |

The failures found are semantics to design around, plus one docstring-contradicting bug to report upstream (tandem copies inside `fates="family"` units). None calls for the §1.1 fallback.

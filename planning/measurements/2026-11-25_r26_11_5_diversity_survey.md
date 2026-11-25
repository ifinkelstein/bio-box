# `r26.11.5+` Diversity Expansion Survey: Cleanest × Highest-Diversity Sources, OpenAlex Citation Chains & `OpenGenome2` Audit

- **Date**: 2026-11-25
- **Method**: Six parallel read-only reconnaissance audits (`r5-inventory`, `r5-phages`, `r5-longread`, `r5-extremophiles`, `r5-uncultured`, `r5-euk-rna`) covering: (1) full audit of built (`r26.11.3`), queued (`r26.11.4`), staged (`/stor/scratch`, `/stor/work`, `tower`), and deferred (`CATALOG_GOAL.md`) sources; (2) anonymous OpenAlex citation-chain traversal across 2022–2026 viral, long-read, extremophile, uncultured, organellar, RNA-virus, and microbial-eukaryote literature; (3) subset-by-subset audit of `arcinstitute/opengenome2`; and (4) ranking strictly by **Cleanliness × Expected Diversity Yield** (complete/circular/HQ or curated-GFF3 assemblies with clean metadata, low contamination, and highest net-new NR50/Gbp or new biological-realm breadth first).
- **Number labels (`AGENTS.md` §6)**: Every number carries `*measured*`, `*derived*`, `*published*`, or `*estimate*`.

---

## 1. Where We Stand After `r26.11.4` & Our Biggest Remaining Gaps

From [`release_diversity.json`](release_diversity.json) and [`2026-11-24_a26_r26_11_4_confs.md`](2026-11-24_a26_11_4_confs.md):

| Source Class | `r26.11.3` (*measured*) | `r26.11.4` queued (*measured* staged) | `r26.11.4` Share (*derived*) | Biggest Post-`r26.11.4` Blind Spots |
|---|---:|---:|---:|---|
| Viruses & phages | `8.662 Gbp` | `83.200 Gbp` | `29.0%` | Eukaryotic giant/divergent dsDNA viruses (`GV1`, `EGV1`, `MIR1`), non-gut extreme/soil/rumen/cold-seep/vent viromes |
| Uncultured archaea & MAGs | `10.043 Gbp` | `60.891 Gbp` | `21.2%` | `193` unbuilt `GlobDB` bacterial phyla (`349 Gbp` on `/stor`), deep-ocean SAGs (`GORG-Dark`), post-`GlobDB` extremophiles (hot springs, vents, cold seeps, caves, closed ONT DPANN, protist endosymbionts), `OMD v2` marine |
| Complete prokaryote isolates | `59.578 Gbp` | `59.578 Gbp` | `20.8%` | Intentionally held (`C1..C4, C6` deferred under `DEC-70`; `10.5k` NR50/Gbp) |
| Unbinned & long-read metagenomes | `35.689 Gbp` | `56.295 Gbp` | `19.6%` | Dereplicated long-read cMAG catalogs (`MFD-LR`, `MiDAS Global`, `ONT Compost`, `LR2`), `SPX2` (`plant host` & `anthropogenic`), `Sandpiper × Logan` dark matter |
| Plasmids & MGEs | `5.925 Gbp` | `24.417 Gbp` | `8.5%` | Cryptic environmental circular plasmids (`Urban Water` 10.9k), `Inocles`/`Tycheposons` (`ECE2`), `PIMEs`, `proMGE2` |
| Organelles | `0.000 Gbp` | `2.653 Gbp` | `0.9%` | Complete fungal/giant-plant mitogenomes (`O3`), cryptophyte nucleomorphs, nitroplasts, apicoplasts, dinoflagellate/kinetoplastid organelles |
| **RNA viruses, Obelisks, viroids** | `0.000 Gbp` | **`0.000 Gbp`** | **`0.0%`** | Entire realm absent (blocked on 4 small, zero-schema-break `C-50` edits) |
| **Nuclear microbial eukaryotes** | `0.000 Gbp` | **`0.000 Gbp`** | **`0.0%`** | Entire domain absent (`ET1` spliced cDNAs ready now; `C-60` Wave 1 `Y1000+` `CC0` yeasts + `VEuPathDB` intron-poor protists) |
| **Total** | **`119.897 Gbp`** (`22.16M` NR50) | **`287.034 Gbp`** (`~64M` NR50 *est.*) | **`100.0%`** ($D = 4.9 / 6$) | — |

---

## 2. `arcinstitute/opengenome2` Audit: Raw Nucleotides Are Available, But Skip >85%

- **Availability (*published*, Hugging Face `arcinstitute/opengenome2`, Brixi et al. *Nature* 2026):** **Yes** — raw compressed DNA FASTAs (`fasta/*.fasta.gz`), pre-tokenized Evo 2 `jsonl/` shards, and `species_metadata.csv` (~15.3 MB) are publicly downloadable under **Apache 2.0** (~8.84 Tbp across ~5.52 TB).
- **Buy-vs-Skip Verdict per Subset:**
  - **SKIP `fasta/prokaryotes`** (~351 B bp): Older `GTDB v214.1/v220.0` + `IMG/PR` (we hold `GTDB R226`, `GlobDB r226`, and `PL3` with `DEC-65` DTR trims & circularity).
  - **SKIP `fasta/viral`** (~34 B+ bp): Older `IMG/VR` with **eukaryotic viruses stripped for biosecurity** and CheckV/geNomad/DTR metadata lost (we hold `V1`, `V2`, `V3`, `MV1`, `CG1`, `H1`, `SAT1`).
  - **SKIP `fasta/organelle`** (~2.8 B bp): Superseded by our fresher RefSeq r238 `O1`+`O2`+`LPL1` (`2.92 Gbp`), which preserve per-record `transl_table` (`1,3,4,11,16,22`) and circularity flags.
  - **SKIP `fasta/eukaryotes`** (~7.0 Tbp): ~95% macroscopic `Animalia`/`Plantae` + synthetic mRNA/window/promoter views; its `Protista`+`Fungi` (~300 Gbp) strip the curated GFF3 exon/intron models we get natively from `Y1000+`, `VEuPathDB`, and NCBI RefSeq under `C-60`.
  - **SKIP `fasta/rna`**: Short isolated ncRNAs without genomic neighborhood context.
  - **AUDIT-ONLY (Low Priority): `fasta/metagenomics` (`MGD DB`, ~854 Gbp):** Only potentially net-new prokaryotic slice (unauthenticated bundle of `JGI IMG/M`, `MG-RAST`, `MGnify`, `Tara`, `NCBI`), **but** Evo 2 preprocessing **splits contigs at `>= 5` consecutive `N`s** (breaking `±10`-gene synteny), strips per-contig SRA/sample/habitat/quality metadata (breaking `C-43` run deduplication), and carries JGI upstream attribution caveats (`DEC-51`). Do **not** bulk-ingest ahead of clean upstream catalogs (`MFD-LR`, `MiDAS Global`, `OMD v2`, `SPX2`, `Sandpiper × Logan`).

---

## 3. Ranked `r26.11.5+` Recommendations (Strictly Ordered by Cleanliness × Diversity Yield)

### Tier 1 — Immediate Zero-Fetch & Compact Jewel-Box Stages (`~38–42 Gbp` total, highest cleanliness × NR50/Gbp)

All datasets below are either **already verified on `/stor/scratch`**, **already code-ready & actively staging to `tower` (`Batch 1`)**, or **sub-Gbp closed/circular/curated open-access archives** that fill major taxonomic/functional blind spots with minimal noise:

| Rank | Stage | Dataset & Primary Source / DOI | Size & Counts | Cleanliness & Diversity Payoff | License / `redist.` | Status |
|---:|---|---|---|---|---|---|
| **1a** | **`LR2` + `SAG1` + `NX1` + `GV1`** | Verified zero-fetch quartet on `/stor/scratch` & `/tmp` (`2026-11-24_lr2_nx1_inputs.md`, `c58`, `a26`) | **`31.512 Gbp`** net new (`74,137` genomes, `1.32M` contigs, *measured*); `~6.9–8.8M` NR50 (*est.*) | **0 network fetch, 0 blockers**: doubles long-read sequence (`LR2` `21.55 Gbp`, `1,863` circular), opens marine SAGs (`SAG1` `4.92 Gbp`), adds extremophile/holobiont MAGs (`NX1` `3.40 Gbp`), brings in giant viruses (`GV1` `1.64 Gbp`) | `CC BY 4.0` / `INSDC` (`SAG1` `CC BY-NC` `false`) | **Ready on pod `/stor/scratch`** |
| **1b** | **`EGV1` + `ECE2` + `MIR1` + `ASG1`** | Pre-build ready `C-69`/`C-70` compacts (*Egovirales* `figshare.25398895`, `Inocles` `zenodo.16656823` + `Tycheposons`, *Mirusviricota* `figshare.28955240` + `PLVs`, Appler 2026 *Asgardarchaeota* `figshare.26426707`) | **`< 2.1 Gbp`** combined (*est.*) | Ultra-clean, `DEPTH=annotated`; unlocks *Egovirales*, giant oral `Inocles`, marine `Tycheposons`, *Mirusviricota* viral phylum, and `404` post-`GlobDB` Asgard MAGs | `CC BY 4.0` / `MIT` (`true`) | **Staging to `tower` (`Batch 1`)** |
| **1c** | **`DPN1`** | ***Nanobdellota* Complete Closed ONT Genomes** (Nielsen & Lui May 2026, `10.64898/2026.05.13.725050`, Zenodo `10.5281/zenodo.20174424`) + **Wu 2024 `DPANN-MGE`** (Zenodo `10.5281/zenodo.11215303`) | **`256` DPANN (`208` complete closed circular!)** (`~0.3 Gbp`) + `3,662` DPANN/MGEs (`~1–2 Gbp` post-dedup *est.*) | Closed circular ultra-small DPANN chromosomes rotated to `ORC1/Cdc6` + DPANN MGEs systematically missed by CheckM | `CC BY 4.0` (`true`; verify `20174424`) | **Top `Batch 2` `tower` pull** |
| **1d** | **`SYM2`** | **Protist Single-Cell Endosymbionts, Host-CPR & Giant Viruses** (Schulz/Woyke *Nat. Commun.* Nov 2025, `10.1038/s41467-025-65263-4`, `https://portal.nersc.gov/cfs/nelli/ms_protistsymb/`) | `117` ciliate/amoeba endosymbionts + `258` host-CPR + `>80` giant viruses (`~0.5–0.8 Gbp` *est.*) | CheckM-excluded intracellular clades (*Holosporales*/*Chlamydiae*/*Babelota*/code-4 & code-25 CPR); `DEPTH=annotated` | `CC BY 4.0` (`true`) | **Top `Batch 2` `tower` pull** |
| **1e** | **`LRP1` (`PL4`)** | **European Urban Water Cryptic Circular Plasmidome** (He et al. *Environ. Int.* 2024 / *Microbiome* 2026, Zenodo `10.5281/zenodo.13143815` / `PRJEB85938`) | **`10,942` true circular plasmids** (`9,538` = `87%` novel, **`< 0.3 Gbp`** *est.*) | Purified circular environmental plasmidome; ~0% overlap with `PL1..PL3`; defense-island & TA-system hotspot (`DEPTH=annotated`) | `CC BY 4.0` (`true`) | **Top `Batch 2` `tower` pull** |
| **1f** | **`MIR2` + `MGE4`** | **Divergent Viral/MGE Micro-Pack**: **Bulzu 2025 Long-Read PLVs** (`zenodo.14540677`, `>1,000` complete PLVs) + **Bellas 2026 PLV Framework** (`figshare.30112102`) + **`PIMEs` 2026** (`zenodo.22794128`, immune-rich IMEs) + **Capsid Satellites** (`zenodo.24320115`, `916` seqs) + **`AVGD` 2026 Archaeal Viruses** (*Microbiome* `10.1186/s40168-026-02445-2`, `3,708` genomes) | **`< 0.5 Gbp`** combined (*est.*) | Extreme NR50/Gbp and anti-phage/defense density; full-length long-read PLVs & archaeal viruses (`DEPTH=annotated`) | `CC BY 4.0` (`true`; verify per deposit) | **Top `Batch 2` `tower` pull** |
| **1g** | **`O3`** | **Organelle Expansion (>2.5× breadth, zero intron caller)**: **Sarhan 2025/2026** (`2,695` complete circular fungal mtDNAs, `figshare.28750034`, `41.3 MB`, `CC BY-NC-ND` -> split `O3b`), **`PMGD`/GenBank giant plant mtDNAs** (`1,829` complete), **cryptophyte nucleomorphs** (`0–17` introns), **nitroplasts** (`UCYN-A`), **apicoplasts**, **dinoflagellate plastid minicircles**, **kinetoplastid maxicircles**, & on-pod **Kavagutti HiFi mtDNAs** | **`~1.0–1.5 Gbp`** (*est.*) | Triples fungal mtDNAs, doubles plant mtDNAs, adds vestigial eukaryotic nucleomorphs & nitroplasts using existing `O2` code routing (`DEPTH=annotated`) | `INSDC`/`CC BY 4.0` (`O3a`) + `CC BY-NC-ND` (`O3b`) | **Top `Batch 2` `tower` pull** |

---

### Tier 2 — Clean, Dereplicated Flagship Compendia (`~110–150 Gbp`, massive NR50 + biome leap)

Every dataset in Tier 2 is **already species/vOTU-dereplicated or long-read cMAG/MQ-HQ filtered** (clean contiguity, minimal strain bloat, verified open-access bulk download):

| Rank | Stage | Dataset & Primary Source / DOI | Size & Counts | Cleanliness & Diversity Payoff | License / `redist.` |
|---:|---|---|---|---|---|
| **2a** | **`SAG2`** | **`GORG-Dark` Aphotic Ocean SAGs** (Chang/Pachiadaki *Cell* 2026, `PRJNA1259050` + `figshare.29087090` + `zenodo.22086263`) | `9,698` SAGs (`8,088` ENA assemblies, **`6.466 Gbp`** *measured*, **`100.0%` novel** vs GTDB/GlobDB/IUIG) | Deep-ocean single-cell genomes excluded by MAG pipelines; unlike `SAG1`, `PRJNA1259050` is public INSDC (`redistributable: true`) | `INSDC` / `CC BY 4.0` (`true`) |
| **2b** | **`LR-COMP`** | **`ONT Mature Compost` (`compost_mgx_2026.04`)** (`s3://ont-open-data/compost_mgx_2026.04`, Apr 2026) | **`5,598` MQ/HQ MIMAGs incl. `1,353` single-contig circular cMAGs** (`~18–22 Gbp` *est.*) | Public AWS S3 `--no-sign-request`; `1,353` closed circular R10.4.1 (`metaMDBG`) chromosomes; `0%` overlap with `r26.11.4` | ONT Open Data (verify bucket `LICENSE`) |
| **2c** | **`MFD-LR`** | **Microflora Danica Long-Read MAGs** (Sereika et al. *Nat. Microbiol.* 2025, `10.1038/s41564-025-02062-z`, Zenodo `10.5281/zenodo.14537680` `MFD-LR-drep.tar.gz`, `20.8 GB`) | **`15,314` novel terrestrial soil/sediment ONT species MAGs** (`+8%` prokaryotic tree, `~55–65 Gbp` *est.*) | **Replaces pulling raw `254.7 Gbp` `PRJEB58634` from `LR1_left_for_later`!** Distinct from Illumina `MFD.ilm_mag.tar.gz` | `CC BY 4.0` (`true`) |
| **2d** | **`MIDAS-LR`** | **MiDAS Global Genome Catalog** (Liu/Albertsen/Dueholm *bioRxiv* Jul 2026, `10.64898/2026.07.10.737647`, Zenodo `10.5281/zenodo.19066579` / `PRJEB83983`) | **`53,501` ONT R10.4.1 MAGs -> `12,047` global WWTP species** (`~42–55 Gbp` species-drep *est.*) | `14×` global expansion of `Singleton 2021` (`LR2`); **replaces pulling raw `2.55 Tbp` `PRJEB122849`** first; rich in *Myxococcota*/CPR | `CC BY 4.0` (`true`) |
| **2e** | **`GHM1` + `VHM1`** | **Global Hot Spring Microbiomes & Virome** (Chen/Tiedje Nov 2026, `figshare.35073963`, **`12.67 GB`** archive *published*) | **`7,885` species MAGs (`>70%` novel)** + **`36,936` viral spp** (`~21–27 Gbp` *est.*) | Single clean `CC BY 4.0` deposit covering global terrestrial geothermal archaea, bacteria, and thermophilic viruses | `CC BY 4.0` (`true`) |
| **2f** | **`NX2`** | **Global Deep-Sea Extremophiles**: **`DSV70`** Hydrothermal Vents (St. John/Reysenbach *Sci. Data* Jan 2026, BioStudies `S-BSST222717` / `PRJNA1244896`, `7,422` MQ/HQ MAGs, `18.3%` Archaea, `~14–18 Gbp`) + **`CSMC`** Cold Seeps (`figshare.22568107`, `3,164` spp `~94%` novel, `~6–8 Gbp`) + **`Hatiba Mons`** Red Sea vents (`zenodo.19235833`, `1,217` MAGs `~70%` novel) + **Global Caves** (*Sci. Data* Oct 2025, `figshare.29955365`, `1,979` cave spp, `~5–7 Gbp`) | **`~13,782` extreme-biosphere MAGs** (`~28–36 Gbp` pre-dedup *est.*) | All `CC BY 4.0` / `INSDC`, post-`GlobDB r226` chemosynthetic & subterranean dark matter (`skani` 96% dedup vs `G1`/`GlobDB`) | `CC BY 4.0` / `INSDC` (`true`) |
| **2g** | **`OMD2`** | **Ocean Microbiomics Database v2 (`OMDB`)** (Sunagawa Lab 2025–2026, `https://omdb.microbiomics.io/`) | `274,282` marine genomes -> **`31,058` species reps** (`3.7×` OMD v1, `2×` GOMC) -> **`~8k–15k` net-new species (`~20–35 Gbp` *est.*)** after `skani` vs `GTDB+GlobDB` | Cleanest, best-annotated global marine prokaryotic expansion post-`GlobDB r226` | `CC BY 4.0` (`true`) |
| **2h** | **`DSV1` + `TSV1` + `RV1` + `WWV1`** | **Clean Non-Human-Gut Virome Wave** (all Zenodo/Figshare `CC BY 4.0` / `CC0`):<br>• **`DSV1`** (`~5.5–7.5 Gbp`): `GCSV` Cold Seep (`193,465` vOTUs, Asgard/DPANN viruses) + `Langwig 2024` Vents (`49,962` genomes, `370 MB`, `figshare.25968037`) + `DSVv1` Deep Subsurface to 3,136 m (`10,422` vOTUs `>86%` novel, `figshare.34745214`) + `GCVD` Coral (`>20k`, `figshare.24968805.v1`)<br>• **`TSV1`** (`~3.5–5.5 Gbp`): `GSV Atlas` (`916 MiB` already in `sources.yaml`!) + `TPLV` Tibetan Lakes (`>71k`, `zenodo.21773290`) + `UK Peatland` (`figshare.28143446`) + `GSAPD` Active Soil Prophages (`21,397`, `zenodo.19510112`) + `GtRVC`/`GFHPD` glaciers<br>• **`RV1`** (`~6.9 Gbp`): `RVD` Rumen Virome (`397,180` vOTUs, `91%` unclassified, `zenodo.7412085`)<br>• **`WWV1`** (`~2.0–3.5 Gbp`): `ASVC` Sludge (`85,114` vOTUs, `zenodo.7847962`) + `WVDB` Ultra-Deep Wastewater (`>21k` near-complete, `zenodo.20276352`) | **`~18–23 Gbp`** combined (*est.*) | Massive non-human-gut viral diversity (`gv-meta` + `DEC-65` DTR trim) filling soil, peatland, glacier, alpine lake, cold seep, hydrothermal vent, subsurface, coral, rumen, and wastewater virome gaps | `CC BY 4.0` / `CC0` (`true`) |
| **2i** | **`LR-EXT` + `LRV2`** | **Compact Weird-Habitat HiFi & Long-Read Virome Bundle**: **Aquatic Symbiosis Genomics (`PRJEB43743`)** marine sponge/coral/clam HiFi symbiont MAGs + **Arctic Permafrost HiFi** (`figshare.32054802`, `~500` MAGs) + **Katase Hot Spring HiFi** (`14` cMAGs) + **Tibetan Saline Lake HiFi** (`426` MAGs + `29` circular viruses) + **Termite HiFi cMAGs** (`PRJNA1126444`, `14` cMAGs) + **PacBio Revio Sludge** (`zenodo.23755631`, `447` MAGs) + **Andaman Sea ONT Virome** (`1,604` HQ) + **Cook 2024 (`zenodo.10650983`) / Wirbel 2026 (`figshare.30900317`)** | **`< 6 Gbp`** combined (*est.*) | High-accuracy HiFi/ONT closed chromosomes & viral genomes across holobionts, permafrost, hot springs, hypersaline lakes, and deep ocean (`DEPTH=annotated`) | `CC BY 4.0` / `INSDC` (`true`) |

---

### Tier 3 — Unlocking the Two Missing Realms (`0.0%` -> Live in `r26.11.5–r26.11.6`, `~22–28 Gbp`, near-100% net-new NR50)

1. **Immediate Bridge (No `C-60` needed): `ET1` Spliced Protist & Algal Transcriptomes (`C-68`, `~8–15 Gbp`, `MMETSP` `zenodo.3247846` + `1KP` `GigaDB 100910`)** — code, tests, and template already ready; unlocks nuclear eukaryotic proteins (`POLRMT`/`RpoT`, cGAS/viperin/Argonaute) across `~1,821+` protists & algae. *(Staging in `Batch 1` now.)*
2. **Fast-Track `C-50` — Global RNA Virosphere, Obelisks & Viroids (`~0.35 Gbp` Phase 1 -> `~2.5–3.5 Gbp` Phase 2, `>95%` net-new NR50):**
   - Requires only **4 small, zero-schema-break edits** (`U -> T` at staging, `"rna_virus"`/`"cccrna"` enums, `rna-orf6`/`cccrna-orf6` on existing `six_frame_orfs`, `ribozyme-1` ~22 CMs from pinned `Rfam 15.1`).
   - Phase 1 (`~0.35 Gbp`, `<2 h` on pod): RefSeq `Riboviria` + `VIRE Riboviria` (`8,056`) + `Obelisks` (`29,959` + `Logan` `26,263` 3.7× expansion) + `ViroidDB` + `Ambiviricota` + `Lee/Neri 2023` (`11,378` cccRNAs).
   - Phase 2 (`~2.5 Gbp`): `LucaProt` (`161,979` spp, `figshare.26298802`) + **`UPVAtlas` 2026** (`54,945` urban RvOTUs, `77%` novel, `10.1038/s41467-026-73605-z`) + **Global WWTP RNA Virome** (`11,414` vOTUs, prokaryotic RNA phages) + `RVMT v3` (`zenodo.7368133`) + `Serratus`.
3. **Cleanest `C-60` Wave 1 — Intron-Poor Curated-GFF3 Microbial Eukaryotes (`~18–20 Gbp`, `<4%` clustering load):**
   - **`Y1000+` Budding Yeasts** (Opulente et al. 2024 *Science*, Figshare+ `10.25452/figshare.plus.22802147`, **`CC0` Public Domain**, **`1,154` genomes / `~1,090` spp**, **~95% intronless genes** with curated GFF3s, `~14 Gbp`) — cleanest possible first nuclear-eukaryote wave!
   - **`VEuPathDB` Intron-Poor Protists** (`TriTrypDB` `80+` kinetoplastids, `MicrosporidiaDB` `50+` compact parasites, `GiardiaDB`, `CC BY 4.0`, `~3.5 Gbp`, virtually 0 cis-introns + curated GFF3s) + **Gaïa *Parabasalia* hosts** (`~1 Gbp`) + **`starbase` 2025/2026** (`zenodo.17533381`, `5,493` giant fungal `Starship` transposons, `~0.6 Gbp`).
   - Followed in Wave 2/3 by on-pod `IMG-UMAG` `E1`/`E2` (`14.44 Gbp`), `Tagirdzhanova` lichen eMAGs (`10.16 Gbp`), `OMD` eMAGs, `Tara SMAGs`, and `TOPAZ`.

---

### Tier 4 — High-Yield Unbinned & On-Disk Scale Expansion (`r26.11.6+`)

1. **Zero-Fetch On-Disk Scale (`wilkcomp02` `/stor`):** **`MV2`** (`3.91 Gbp`, `98,920` `>= 30 kb` `MetaVR` reps) + **targeted `GlobDB r226` rare-phylum waves** (`G2..G41` holds `130,960` remaining non-GTDB species / `349 Gbp`, `~200k–266k` NR50/Gbp — already absorbs `SMAG` Ma 2023!) + **`J1` `IUIG`** (`~9.65 Gbp`, `~4,235` novel reps).
2. **`SPX2` Expanded-Habitat SPIRE Tranche 2 (`35–40 Gbp >= 5 kb`, `gv-meta`):** Pull **untouched high-yield L2 strata** from `2026-11-18_spire_pilot.md` (`host-associated:plant host` `568k` NR50/Gbp @ `98.1%` retention; `anthropogenic:{waste,agriculture,wastewater}` `492k–650k`/Gbp; `aquatic:briny water` `696k`/Gbp) while skipping `air` (`91k`), `artificial` (`66k`), and `animal host` (`230k`).
3. **`Sandpiper × Logan` Dark-Matter Tranche (`LGDM1`, `35–40 Gbp >= 5 kb`) + Full Multi-Panel `logan50` (`LGRT2` + `LGDEF1`):** Anti-join `Sandpiper` (`zenodo.10547493`, `248k+` SRA SingleM profiles) against `C-43` (`182k` runs) to stream `s3://logan-pub/c/` contigs from un-assembled SRA runs with lowest `known_species_fraction`; finish the remaining `99.03%` of `logan50` (`~2.6 h` at 40 overnight threads) in one pass with rare-RT + CRISPR/defense HMMs.
4. **`LR3` Non-Aalborg Open Unbinned Sweep (`LR1_left_for_later.tsv.gz`, `~45–55 Gbp`):** All `17` freshwater (`11.2 Gbp`) + all `28` sludge (`9.1 Gbp`) + top `40` submitter `marine_sediment` (`~25–30 Gbp`) + `15` `LR1` retries (skipping raw `PRJEB58634`/`PRJEB122849` until after `MFD-LR`/`MiDAS Global`).
5. **Next PI JGI Portal Batch:** `532` unrestricted JGI long-read assemblies (`LR_jgi_batch.tsv.gz`, `257.48 Gbp`) ahead of bulk short-read `IMG/M`.

---

## 4. Explicit Skips, Embargoes & Traps Identified

- **Skip `OpenGenome2`** (`>85%` superseded/stripped; `MGD DB` splits scaffolds at `>=5 N`s & drops metadata).
- **Skip standalone `SMAG` (Ma 2023):** Already absorbed into our on-disk `GlobDB r226` (`6,216` reps, `18.70 Gbp`, `*measured*`).
- **Skip `AllTheBacteria` (ATB):** `2.44M` assemblies collapse to only `11,273` GTDB species (`216×` clinical strain redundancy, `*published*`). Also skip `BV-BRC`, `BiG-FAM`/`MIBiG` (use as oracles, not sequence stages), bulk `MetaSUB`/`ISS`/`Plastisphere`, and raw-read-only `NEON`/`AusMicrobiome`/`EMP500`/`TREC` (harvest assembled outputs via `NMDC`/`MGnify`/`Logan` instead).
- **Do NOT bulk-pull raw `PRJEB58634` (`254.7 Gbp`) or `PRJEB122849` (`2.55 Tbp`)** before `MFD-LR` (`20.8 GB`) & `MiDAS Global` (`12,047` species).
- **Embargoed / portal-gated (do not schedule yet):** Chinese Saline Lakes (`zenodo.22827828`, `48,830` MAGs — **embargoed until 2027-12-31**; replaced by open **`SODA1` `PRJNA857294`** below); `MEER` Hadal (`OEP004067` / `CNP0004890` — JS-portal-only; replaced by open **`ATAC1` Atacama Trench `PRJEB57745`** below); `V4` `PEONY`/`PIGEON 3.0` (`dryad.prr4xgz35` — needs PI Dryad ORCID token); `Chlorarachniophyte` nucleomorphs (`>850` pygmy introns -> defer to `C-60`).

---

## 5. Addendum — Aquatic & Unknown-Unknown Blind-Spot Scouts (`r5-aquatic-scout` & `r5-blindspots`, 2026-11-25)

Follow-up scouts expanded **clean aquatic compendia** (16 net-new sets) and **untapped ecological niches & public databases** (15 net-new sets) ranked strictly by **cleanliness × expected diversity yield**:

### 5.1 Top Clean Aquatic Picks (All Open `CC BY 4.0` / `INSDC`, Zero Overlap with `r26.11.4` or §3)

| Stage | Dataset & Primary Source / DOI | Counts & Size (*published* / *estimate*) | Why High Cleanliness × Diversity |
|---|---|---|---|
| **`FUTM1` + `GXM1`** | **Futian Mangrove 7-Yr (`figshare.30883646`, *Sci. Data* 2026)** + **Guangxi Mangrove (`figshare.29320385` / `PRJEB96880`, 2025/2026)** | `3,404` spp reps (`336` Arch, **82 phyla**) + `2,965` MAGs (**582 Arch = 19.6%**, **78 phyla**); `~15–19 Gbp` pre-joint-dedup | Intertidal anoxic mangrove sediments; `<10%` overlap with `GlobDB`/`OMD v2`; dedup jointly at 96% ANI |
| **`FLHC1`** | **Freshwater Lake metaHi-C Genomes & Virome (`zenodo.15856862`, 2025/2026)** | **`>900` Hi-C-resolved MAGs + `>33,000` viral genomes** with physical virus-host links (`~2.5–4.0 Gbp`) | Physical Hi-C linkage resolves freshwater phage-host pairs (*Synechococcus*, *Limnohabitans*, *Planktophila*) |
| **`GFS1`** | **Vanishing Glaciers — Glacier-Fed Streams (`PRJNA781406`, *Nat. Microbiol.* Jan 2025)** | **`2,855` bacterial MAGs** (+ archaea/fungi/algae/viruses) across 9 mountain ranges (`~7–9 Gbp`) | High-altitude glacial-stream cryosphere dark matter |
| **`MOS1`** | **MOSAiC Central Arctic Ocean & Sea-Ice (`figshare.25746222` / `PRJNA1106878`, *Sci. Data* Feb 2025)** | **`2,407` prokaryotic + `56` eukaryotic MAGs** (`~5–7 Gbp`) | Multi-year Arctic sea-ice, melt-pond, & ice-ridge microbiome |
| **`SODA1`** | **Global Soda Lake Biogeography (`PRJNA857294` + `figshare.28152113`, *Front. Microbiol.* Jul 2025)** | `2,227` MAGs -> **`1,330` species** across 4 continents (`~3.5–5.5 Gbp`) | **Open non-embargoed replacement** for embargoed Chinese Saline Lakes (`zenodo.22827828`) |
| **`MACA1`** | **Macroalgal Epiphytic Bacterial Atlas (`zenodo.21563086`, *Phycology* Sep 2026)** | **`1,915` QC MAGs** from paired long/short-read kelp/macroalgae (`~5–7 Gbp`) | Under-sampled marine algal epibiont niche (`<10%` overlap) |
| **`NGOM1` + `ATAC1`** | **Northern Gulf of Mexico (`PRJNA1171802`, `1,313` strict-QC MAGs)** + **Atacama Trench Hadal (`PRJEB57745`, `1,357` MAGs to 8,085 m)** | `2,670` MAGs (`~5.5–7.5 Gbp`) | **`ATAC1` is open ENA replacement for JS-gated `MEER` Hadal**; `NGOM1` is `>=75%` comp / `<=5%` cont incl. *Asgardarchaeota* |
| **Wave B Aquatic** | **`GROWdb`** (`2,093` unique river MAGs), **`DWGC`** (`1,141` drinking water spp), **`LakePulse`** (`1,008` Canadian lake spp), **`GPDC`/`PRJNA858913`** aquifer defensome (`141` phyla), **`HYP1`** Great Salt Lake + Cryopeg brines | `~15–20 Gbp` combined | Clean freshwater/river/aquifer/brine depth for `r26.11.6` |

### 5.2 Top Unknown-Unknown Blind Spots (Untapped Biomes & Open Portals)

| Stage / Portal | Dataset & Primary Source / DOI | Counts & Size (*published* / *estimate*) | Why High Cleanliness × Diversity |
|---|---|---|---|
| **`CANDI1`** | **`CAnDi` — Catalogue of Anaerobic Digestion (`zenodo.23322224`)** | **`36,480` MAGs -> `11,149` prok spp + `~1.9M` viral seqs (`~536k` vOTUs)** across `1,256` global AD/biogas reactors | Massive syntrophic methanogen/anaerobic-bacteria + virome reservoir; direct open Zenodo |
| **`SEGC1`** | **`SEGC` — Standardized Extremophilic Genomic Catalogue (*bioRxiv* Apr 2026, `10.64898/2026.04.17.719132`)** | **`54,661` MQ/HQ MAGs -> `21,805` spp (`66.1%` novel!) + `162,855` BGCs** across 7 extreme niches incl. Acid Mine Drainage (AMD) | Fills AMD, cryogenic, and saline-alkaline extremophile gaps in one standardized catalog (verify bulk FASTA mirror first) |
| **`OIL1`** | **Chinese Oilfield / Deep Petroleum Reservoirs (Zhou *AEM* 2026, `figshare.31337758`)** | **`3,057` MQ/HQ MAGs (`73.77%` novel species!)** across 13 deep oilfields (`~6–8 Gbp`) | Open `CC BY 4.0` Figshare; deep hydrocarbon-biosphere dark matter |
| **`DTOL1`** | **Sanger Tree of Life / DToL / EBP HiFi Cobionts (`MarkerScan`/`ASCC`, ENA REST API)** | Closed circular PacBio HiFi intracellular endosymbionts (`110` complete *Wolbachia* from 368 UK insects alone + ongoing `PRJEB*` cobionts) | Systematically discarded by CheckM/GTDB/GlobDB; unauthenticated ENA REST API |
| **`CFMD1`** | **`cFMD` — Curated Food Metagenomic Data (Carlino *Cell* 2024, `zenodo.10891046` / `13285428`)** | `10,899` MAGs -> **`1,144` SGBs (`1,036` prok + `108` euk, `320` novel)** (`~3 Gbp`) | Clean open Zenodo fermented/food microbiome catalog (rich in lactic acid bacteria, phages & yeasts) |
| **`PLNT1`** | **Non-Soil Plant Micro-Pack**: Forest Phyllosphere (`PRJEB93989`, `423` MAGs `>400` novel spp) + Boreal Moss (`PRJNA1189183`) + *Sphagnum* (`PRJNA1048504`) + Deadwood (`figshare.14821752`) | `~1.5–2.5 Gbp` combined | Targets `SPIRE`’s highest-retention host niche (`plant host` `568k` NR50/Gbp @ `98.1%` retention) |
| **Open Portals** | **EMBL-EBI `MGnify Genomes`** (`ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/` — `non-model-fish-gut`, `honeybee-gut`, `maize-rhizosphere`, `HoloFood`), **`DOE NMDC` open REST API** (`api.microbiomedata.org`), **`CNCB-NGDC GWH`** (`ftp://download.big.ac.cn/gwh/`, unauthenticated FTP unlike JS-gated `NODE`/`OEP`) | Modular per-biome pulls | Clean unauthenticated programmatic access for `r26.11.6+` |


# `r26.11.5+` Deep Non-Duplicate Collections Addendum (Giant Viruses, Phages, Novel Prokaryotes & Login-Gated Portals)

- **Date:** 2026-11-25
- **Companion to:** [`2026-11-25_r26_11_5_diversity_survey.md`](2026-11-25_r26_11_5_diversity_survey.md)
- **Scope:** Synthesis of 4 dedicated deep reconnaissance scouts (`r5-giant-viruses`, `r5-phages-deep`, `r5-prok-portals`, `r5-prok-novel`) hunting clean, high-quality, non-duplicate collections across **both open-access and login-gated/national portals** (JGI Data Portal / MycoCosm / PhycoCosm, CNCB / Biosino NODE / eLMSG / CNGBdb, Dryad ORCID, CyVerse / iVirus, EMBL-EBI MGnify, DOE NMDC, NERSC, Zenodo, Figshare, INSDC).
- **Prioritization Rubric:** **Cleanliness** (closed/HiFi/HQ/species-dereplicated, low host contamination) × **Expected Non-Duplicate Diversity Yield** beyond `GTDB R226`, `GlobDB r226`, `S1..SPX1`, `r26.11.4` (`V2..GART2`), `GV1`, and `tower` Batches 1–4.
- **Number Labels (`AGENTS.md` §6):** Every count and size is labeled `*measured*`, `*published*`, `*derived*`, or `*estimate*`.

---

## 0. Track 0 — Zero-Fetch On-Pod Goldmines (0 Downloads, 0 Logins)

Three high-yield collections are **already sitting on disk** on `wilkcomp02` / `/stor` and were left behind by earlier stage slices:

| Priority | Stage Idea | On-Disk Location | What Was Left Behind | Net-New Scale (*published* / *measured* / *estimate*) |
|---:|---|---|---|---|
| **0.1** | **`MV-EUKV`** | `/stor/scratch/Wilke/if2453/fragdef/db/metavr/METAVR.fna.bgz` | `MV1` capped at top `1M` reps and `MV2 >= 30 kb` misses `15–30 kb` virophages/PLVs. Leaves `>15,000` *Nucleocytoviricota* (incl. novel order *Skiavirales*), *Mirusviricota*, and `15–30 kb` *Maveriviricetes* (virophage) + *Preplasmiviricota* (PLV) vOTUs un-ingested. | `~8,900` NCLDV vOTUs + virophage/PLV slice (`~0.8–1.8 Gbp` net new `*estimate*`; `redistributable: false`) |
| **0.2** | **`GV-HIFI`** | `/tmp/if2453/c67/dl/kavagutti/Support-data.gz` (`NCLDV-data.tar.gz`) | `LRV1` (`C-67`) extracted only `phages-fna.tar.gz`, leaving Kavagutti et al. 2023 PacBio HiFi freshwater spring-bloom NCLDV contigs untouched on disk. | **`3,309` PacBio HiFi NCLDV contigs** (`*published*`, `~0.15–0.35 Gbp` `*estimate*`; `CC BY 4.0`) |
| **0.3** | **`V2_EXT`** | `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/jgi_img_vr_v4{.1_hc,_full}_20261118.zip` | `IMG/VR v5` has **not** released as of Nov 2026 (`v4.1` remains latest). `V2` only took `224,520` HQ HC reps, leaving **`2,693,001` HC vOTUs** and **`10,101,426` putative rows** (`*measured*`) on disk. Taking **HC MQ (`>=50%`) + v4.0 Putative Complete/DTR/HQ** rescues divergent archaeal viruses & phages lacking canonical geNomad hallmarks. | `~400k–600k` MQ + Putative DTR/HQ vOTUs (`~10–15 Gbp` `*estimate*`; `0%` overlap vs `V2`; `redistributable: false`) |

---

## 1. Track A — Open-Access Batch 5 (Auto-Staging Detached on `tower`)

All collections below are **open-access** (no login wall) and can be pulled directly onto `tower:/mnt/user/research-bulk/bio-box/sources_incoming/` in detached Batch 5 (chained behind Batches 3–4, zero laptop streaming).

### 1.1 Giant Viruses, Virophages & PLVs (`~1.5–2.5 Gbp` total, `DEPTH=annotated`)

| Rank | Stage | Collection & Exact Source / DOI | Counts & Size (*published* / *estimate*) | Why High Cleanliness × Non-Duplicate Yield |
|---:|---|---|---|---|
| **G1** | **`GV-POLAR`** | Polar Cryosphere Quartet (4 open Figshare records, 2024–2026):<br>• Pitot 2024 Arctic/Antarctic lakes (`10.6084/m9.figshare.25000571`)<br>• Polar Ecosystems 2026 (`10.6084/m9.figshare.32316834`)<br>• Low-Temp Marine (`10.6084/m9.figshare.26314543`)<br>• Polar Microalgal GEVEs (`10.6084/m9.figshare.29546846`) | `3,304` + `1,573` + marine/algal = **`~4,900+` GVMAGs** (`*published*`); `~0.6–1.1 Gbp` net new (`*estimate*`) | Largest biome-specific giant virus expansion beyond `GV1` (`GVMAGs V2`); cryosphere lakes & sea-ice brines are heavily under-sampled (`>60–75%` net new `*estimate*`). |
| **G2** | **`GV-FUNG`** | ***Mycodnaviridae*** (*PLOS Biol.* Aug 2026, UMich Deep Blue Data `10.7302/p5ay-ng61`, `PRJNA1257431`) + **`GKS`** Alpine Lake Single-Cell GV MAGs (Sep 2026, `zenodo.22834836`) + **Danish WWTP** NCLDV/Virophage/PLV (*ISME J* 2026, `zenodo.17735769`) + **Siberian Permafrost** (`10.1038/s41467-022-33633-x`) + **Yellowstone Geothermal** (`10.1038/s42003-024-05931-1`) | Fungal giant viruses (up to `350 kb`) + `GKS` single-cell + `61` WWTP NCLDVs/`15` virophages/`14` PLVs + permafrost/geothermal (`~0.15–0.35 Gbp` `*estimate*`) | Extreme host-kingdom novelty: unlocks **Kingdom Fungi giant viruses** (`0%` in `GV1`), single-cell host-paired alpine infections, Pleistocene permafrost *Pithoviridae*, and geothermal algal-mat NCLDVs. |
| **G3** | **`GV-ISO`** | Complete Cultured Giant Virus Isolates in NCBI GenBank absent from RefSeq `V1` (NCBI E-utilities/Datasets query on `"Nucleocytoviricota"[Organism]` & `"Maveriviricetes"[Organism]`) | `~250–450` closed isolate genomes (`~0.15–0.30 Gbp` `*estimate*`) | Zero binning chimerism; captures proposed ***Manesviridae*** (*Ushikuvirus*, *Usurpativirus*, *Furtivovirus*, *Clandestinovirus*), *Medusavirus euryale*, newer *Pandoravirus*/*Pithovirus*/*Cedratvirus*/*Yaravirus*/*Orpheovirus*. |
| **G4** | **`GV-GEVE`** | Curated Eukaryotic GEVEs: Moniruzzaman 2020 chlorophytes + *Euglena* *Asfuvirales* (`zenodo.15238723`) + Sarre 2024 *Amoebidium* holozoan GEVEs (`PRJEB68378`) | `~100–250` integrated GEVE loci (`~0.08–0.20 Gbp` `*estimate*`) | Invisible to viral-fraction metagenomes (`GV1`/`V2`/`MV1`). *Note:* route intron-rich loci via `C-60`. |
| **G5** | **`VIR-PLV`** | *Preplasmiviricota* / *Bamfordvirae* backbone beyond `MIR1`/`MIR2`: Roux 2023 ICTV *Maveriviricetes* + Hackl 2021 EMALEs (`33` provirophages) + Starrett/Buck 2021 Adintoviruses (`figshare.14161250`) + *Preplasmiviricota* ref (`zenodo.6787989`) + *Oxyrrhis* fused-PolB PLVs | `~500–1,500` curated virophages, EMALEs, adintoviruses, dinoflagellate PLVs (`<0.05 Gbp` `*estimate*`) | Ultra-compact eukaryotic satellite/virophage/polinton reference backbone. |

### 1.2 Open Phages & Archaeal Viruses (`~10–15 Gbp` total)

| Rank | Stage | Collection & Exact Source / DOI | Counts & Size (*published* / *estimate*) | Why High Cleanliness × Non-Duplicate Yield |
|---:|---|---|---|---|
| **P1** | **`GSVC`** | Global Swine Gut Virome Catalog (Quan et al., late 2026; `zenodo.22159188`; supersedes `PVD`) | **`224,572` MQ/HQ vOTUs** from `5,064` swine metagenomes (`*published*`; `~6–9 Gbp` `*estimate*`) | 95% ANI CheckV-curated non-ruminant mammalian gut phage atlas; `<5–10%` overlap (`*estimate*`) with `V1`/`RV1`. |
| **P2** | **`GHV`** | Global Halophilic Virome (*ISME Commun.* 2026, `zenodo.21450953`, `PRJEB45291`) | **`10,705` halophilic vOTUs** across 20 hypersaline sites on 4 continents (`*published*`; `~0.3–0.5 Gbp` `*estimate*`) | Top open haloarchaeal & halobacterial virus expansion beyond `AVGD`/`SAT1` (`DEPTH=annotated`). |
| **P3** | **`CGVD`** | Chicken Gastrointestinal Virus Genome Catalog (`zenodo.14684892`) | **`39,380` MQ/HQ species-level vOTUs** (`*published*`; `~1.2–1.8 Gbp` `*estimate*`) | **`84.9%` novel**, `<3.1%` overlap (`*published*`) with human, ruminant, or swine viromes. |
| **P4** | **`TYMEFLIES-V`** | Lake Mendota 20-Yr Freshwater vMAGs (`figshare.24915750`) + `GtRVC` Glacier-to-River (`36,358` vOTUs) | `~1.5–2.7 Gbp` combined (`*estimate*`) | Eutrophic freshwater cyanophages/heterotroph phages + proglacial cold-adapted viral dark matter. |
| **P5** | **`BPHAGE`** | Honeybee *Apis mellifera* Gut Phageome (*Nat. Commun.* 2026, `zenodo.23266178`) + Non-gut Human (`HOVD` `zenodo.17143670` + `SVD` skin `2,873` vOTUs, `90.85%` novel) | `~1.1–2.8 Gbp` combined (`*estimate*`) | Fills social-insect pollinator phageome (`<2%` overlap) and oral/skin/airway human phages left out of gut-only `V1`. |
| **P6** | **`ARC-V2`** | Archaeal & Cultured Jewel-Box Pack: Methanogen Virus Compendium (`10.17632/wdyjgcdgbm.1`, `282` HQ + `205` VLP across `87` families) + Meng Li Jul 2026 Asgard Viruses (`10.1186/s40168-026-02481-y`, `23` viruses + first Asgard satellite) + `INPHARED2`/`PhagesDB` post-`9Aug2025` isolate delta (`~2.4k–5k` complete isolates) + `Grandevirales` Lak megaphages (`540–735 kb`, alt-code) | `~0.3–0.5 Gbp` combined (`*estimate*`; `DEPTH=annotated`) | Highest evolutionary novelty per Mbp in archaeal virology + 100% complete isolated phage delta. |

### 1.3 Open Novel Prokaryotes & Open National/Consortium Portals (`~85–135 Gbp` net new after `skani` 96% dedup)

| Rank | Stage | Collection & Exact Source / DOI | Counts & Net-New Size (*published* / *estimate*) | Cleanliness, Genetic Code & Non-Duplication Notes |
|---:|---|---|---|---|
| **K1** | **`SYM3`** | **`SymGs` / `NeLLi2023`** (Villada/Schulz/Woyke *Nat. Biotechnol.* 2026; `portal.nersc.gov/cfs/nelli/symgs/` + `zenodo.20399071`) | **`14,070` bacterial/archaeal symbionts** (`*published*`; `~6–10 Gbp` net new `*estimate*`) | CheckM-bypassed symbionts across half of all prok phyla; distinct from `SYM1`/`SYM2`. Screen codes `4/11/25`. |
| **K2** | **`MGNIFY2`** | **EMBL-EBI `MGnify Genomes` 2024–2026 Non-Human Biomes** (`ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/`): `BlueRemediomics` `marine-v2.0` (`13,223`), `soil-v1.0` (`19,472`), rumen/pig/chicken/`HoloFood`/rhizosphere/fish/honeybee | **`41,041` spp reps** (`*published*`) -> `~12k–18k` net new (`~30–48 Gbp` `*estimate*`) | Uniform `Comp - 5*Cont >= 50`, `Cont <= 5%`, 95% ANI reps. `GlobDB r226` only holds 2023 `MGnify`. |
| **K3** | **`EHI1`** | **Earth Hologenome Initiative DR1** (*GigaScience* 2025; `zenodo.16757580` + metadata `zenodo.16784684`) | **`17,666` MQ/HQ MAGs** (`*published*`) -> `~12–18 Gbp` species-drep net new (`*estimate*`) | Global wild vertebrates (amphibians, reptiles, birds, fish, wild mammals) absent from human/livestock catalogs. |
| **K4** | **`EMERGE`** | **Stordalen Mire Permafrost Thaw v2.0.1** (`zenodo.17468704`, Oct 2025) | `36,419` strict-QC MAGs (`>=70%` comp) -> **`2,048` 95% ANI spp reps** (`*published*`; `~5–7 Gbp` `*estimate*`) | Cleanest Arctic permafrost thaw-gradient catalog; species-rep subset avoids strain bloat. |
| **K5** | **`QXLSG1`** | **Qinghai-Xizang Lake Sediment Catalog** (*Cell Reports* 2025; Figshare + `PRJCA051988`) | `5,866` MAGs -> **`2,742` species** (`*published*`; `~7–9 Gbp` `*estimate*`) | Open complement to embargoed saline lakes (`zenodo.22827828`); `>4,000 m` Tibetan salinity gradient. |
| **K6** | **`ISYM1`** | **`iSymBase`** (*ISME Comm.* 2026; `http://symbiont.insect-genome.com/`) | **`3,648` HQ insect symbionts** incl. `1,570` new MAGs (`*published*`; `~3.5–5.0 Gbp` net new `*estimate*`) | CheckM-bypassed & code-4 obligate insect endosymbionts (`skani` vs `SYM1`/`SYM3`/`DTOL1`). |
| **K7** | **`ASG2`** | **Deep-Origin *Asgardarchaeota* + Asgard MGEs** (Zhang *Nature* 2025 `figshare.34064073` + Valentin-Alvarado 2024 `zenodo.12617226`) | **`223` HQ Asgard MAGs** + closed Asgard genomes/MGEs (`*published*`; `~0.8–1.1 Gbp` `*estimate*`; `DEPTH=annotated`) | Non-Heimdallarchaeia clades independent of `ASG1` (Appler 2026). |
| **K8** | **`MYXO1`** | **Expanded *Myxococcota*** (Agrawal/Müller Feb 2026 + Khanal Pokharel *AEM* 2026) | **`154` giant (`9–16 Mbp`) HQ/complete genomes** (`*published*`; `~1.5–1.9 Gbp` `*estimate*`; `DEPTH=annotated`) | Expands phylum `11 -> 28` families and `32 -> 90` genera; extraordinary BGC & predatory machinery density. |
| **K9** | **`DESERT1`** | **Global Arid + Qaidam Mars-Analog + Atacama + McMurdo Cryptoendoliths** (`10.26180/24302833` + `zenodo.15743210` `1,773` mMAGs `>94%` novel + `zenodo.14537000` + `zenodo.7313591`) | **`4,000+` arid/desert/cryptoendolith MAGs** (`*published*`; `~9–12 Gbp` net new `*estimate*`) | Closes hyperarid desert, Mars-analog hypersaline basin, and Antarctic rock-endolith blind spots. |
| **K10** | **`NMDC-NEON`** | **DOE `NMDC` NEON Continental Soil (`nmdc:sty-11-34xj1150`, open REST API `api.microbiomedata.org`) + `ESS-DIVE` East River Watershed SFA** | `6,000+` assemblies -> `~3k–8k` MQ/HQ MAGs (`~8–20 Gbp` `*estimate*`) | Clean programmatic harvest of continental US `NEON` soils post-`IMG-UMAG` (`40–65%` novel `*estimate*`). |
| **K11** | **`NICHE-PK2`** | **Compact Open High-Novelty Pack:** `MBGG1` wild migratory birds (`5,823` MAGs, `figshare.com/s/8699370b85d1875209b2` + `zenodo.13643296`) + `GBR-SED` coral-reef sediment (`2,698` MAGs -> `1,201` spp, `zenodo.24187765`) + `CODE4-SYM` (`GLMC` leafhoppers `zenodo.20042172` + coral *Mollicutes* + fungal `MRE1`) + `ACIDO-SPG` (`figshare.33269484` + `zenodo.13266492` + `zenodo.21135542`) + `BATHY1` (first circular Bathy-6 cMAG + `119+` MAGs) + `SUBSURF1` (`zenodo.22053915`) + `PARAMO-CPR` (`zenodo.21344952`) + `4GDB` Global Glacier (`4,501` genomes, `nmdc.cn/4gdb/`) + `gcType` non-GenBank type strains (`gctype.wdcm.org`) + Ruminant Archaeome (`PRJEB81441`, `998` archaea) | `~18–30 Gbp` net new combined (`*estimate*`) | Covers wild avian flyways, carbonate reef sediments, code-4 fungal/insect/coral symbionts, sponges, circular *Bathyarchaeia*, deep crystalline bedrock, Andean páramo, Three Poles glaciers, and CNGB/NBRC-only type strains. |

---

## 2. Track B — Login-Gated PI Pull Playbook (Grouped by Portal)

Four logins unlock the largest, cleanest, highest-novelty collections unavailable to anonymous downloaders:

### 2.1 JGI Data Portal / SSO (PI Already Has Account; `redistributable: false`)
1. **`LR_jgi_batch.tsv.gz` (`532` unrestricted JGI long-read metagenomes, `257.48 Gbp` `*measured*`)** — already queued; piggyback `ViralRecall`/`GVClass` giant-virus extraction on the same assemblies!
2. **`IMG/M` Post-Cutoff Deltas (`img.jgi.doe.gov/m` -> `data.jgi.doe.gov` / Globus):**
   - **Post-`2024-02-06` Unrestricted MQ/HQ MAGs** (delta beyond our on-disk `IMG-UMAG` `193,015`-bin snapshot): `~40k–70k` bins -> **`~10k–20k` net-new species (`~25–50 Gbp` `*estimate*`)**.
   - **Post-`2023-11-10` Unrestricted SAGs** (delta beyond our on-disk `IUIG` snapshot): **`~1.5k–3k` SAGs (`~2–4 Gbp` `*estimate*`)** — captures CPR/DPANN/ultra-small lineages missed by binning.
3. **`JGI-GEVE` (`MycoCosm` zoosporic fungi + `PhycoCosm` unmasked algae):** `~150–300` genomes (`~0.1–0.4 Gbp` extracted viral loci `*estimate*`) — primary source of *Mycodnaviridae* and green/stramenopile algal GEVEs.

### 2.2 CNCB / Biosino NODE / CNGBdb (One Free Biosino Account Unlocks 4 Flagship Catalogs!)
Register one free account at `https://www.biosino.org/node` (pull via `lftp -e "set ftp:ssl-force true; set ssl:verify-certificate no; open ftps://fms.biosino.org:2122; user <email> <password>"`) + optional free `CNGBdb` account (`https://db.cngb.org/`):
1. **`MEER` Mariana Trench Hadal (`OEP004067` / `eLMSG LMSG_G000037053.1–44616.1` / `CNP0004890`, *Cell* Mar 2025):** **`7,564` SRGs** (`~18–24 Gbp` `*estimate*`, **`89.4%` novel** `*published*`) — single cleanest hadal prokaryote atlas.
2. **`CPR1` (`PBGD` + `PPGD`, `OEZ00021651`, *Nat. Microbiol.* Nov 2026):** **`24,523` *Patescibacteria* (CPR) genomes (`13,242` HQ)** + **`18,876` patesciphages (`12,867` vOTUs)** (`*published*`; `~6–10 Gbp` net-new CPR + `~0.8 Gbp` phages `*estimate*`; screen code `25`/`11`).
3. **`GDSV 2.0` Global Deep-Sea Sediment Virome (`OEP003844`):** **`554,664` benthic deep-sea DNA vOTUs** from 159 cores (`*published*`; `~3–6 Gbp` MQ/HQ `*estimate*`).
4. **`GCSV` Global Cold Seep Virome (`PRJCA048260`):** **`193,465` species-level vOTUs** (`*published*`, **`>99%` novel**, `~4–6 Gbp` `*estimate*`) — unmatched reservoir of cold-seep archaeal viruses (*Asgard*, *DPANN*, *Halobacteriota*).
- *Same-login bonuses:* `CnGVR` (`OMIX007646`), `GPB` (`CNGBdb GDS0000055`, `104` isolated obligate-anaerobe gut phages), `DeepseaDB` SCS/Cold Seep (`CNGBdb`, `2,655+` MAGs), and `GCM-AMD` (`gcMeta` hyperacidic AMD after checking `SEGC1`).

### 2.3 Dryad (ORCID Login)
1. **`PEONY 1.3` / `PIGEON 3.0`** (`10.5061/dryad.prr4xgz35`, `PEONY1.3.fa.gz` `8.55 GB`, `PIGEON3.0.fa.gz` `9.13 GB` `*published*`): **`>1.35M` encapsidated VLP vOTUs** (`>=10 kb`, 95% ANI; `~28–35 Gbp` `*estimate*`) — **#1 overall portal-gated phage target** (VLP-purified free virions eliminate cellular/prophage contamination).
2. *Optional same-token niche MAGs:* Boreal river Hg-methylators (`dryad.9zw3r22mz`, `850` MAGs), Kelp collapse (`dryad.79cnp5j78`).

### 2.4 CyVerse / iVirus (Free Account or Public `/iplant/home/shared/iVirus/`)
1. **Stordalen Mire Permafrost Virome + `Guliya` 41 kyr Ice-Core Virome** (`10.25739/*`): `~15k–35k` permafrost/polar vOTUs + **`1,705` ancient 41,000-yr glacier vOTUs** (`*published*`, `~0.8–1.8 Gbp` `*estimate*`).

---

## 3. Verified Duplicate Traps & Explicit Skips (Do Not Fetch)

| Source | Reason to Skip |
|---|---|
| **`Jug` Phages** (*Sci. Adv.* Feb 2026, `10.1126/sciadv.aeb6265`) | **Already ingested in `H1`** (`HPGC.7230.genomes.final.fasta.gz`, `figshare.30506864.v1`). |
| **JGI `GVMAGs v1`** (Schulz 2020, tape restore) & bulk **`GVDB`** / *Tara* NCLDVs | Superseded by **`GV1`** (`GVMAGs V2`, `18,499` genomes) and **`MIR1`** (`GOEV`). Do not file a JGI tape-restore ticket. |
| **`IMG/VR v5`** | Not released as of Nov 2026 (`v4.1` remains latest; both `v4.1 HC` and `v4.0 full` are already on `/stor`). |
| **Bulk `PhageScope 2.0`, `GPD`, `Prophage-DB`** | Heavy overlap with `V1` (`UHGV`), `H1` (`HPGC`), `MV1`, `V3`, `P1..P6`/`G1`/`SPX1`. Only surgically slice `PhageScope 2.0` (`VMGC`/`ELGV`/`BAPS`) later if needed. |
| **Old `GEBA` waves (`Pilot`, `KMG`, `RNB`, `MDM`)** | Already in `GTDB R226` and on-disk `IUIG` (`2023-11-10`). |
| **`TPMC`, `SMAG`, `GOMC`, `SHGO`, `RGM`, `GFS`** (Chinese portals) | Already inside **`GlobDB r226`** on `/stor`. |
| **`Bioplatforms Australia`, `EMBL TREC`** | Raw FASTQs/ASVs only (no centralized MAG catalog yet; monitor MGnify API for future `TREC` MAGs). |
| **`DDBJ/JAMSTEC`, `DOE KBase`, `DSMZ/BacDive`** | Fragmented per-study only (`DDBJ`), mirror of JGI/ESS-DIVE (`KBase`), or INSDC metadata hub (`BacDive`). |

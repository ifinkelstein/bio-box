# PI Login-Gated & Manual-Action Data Pull Checklist (`r26.11.5+`)

*Created: 2026-11-25 | Scope: Every high-value `bio-box` source that requires a portal login, API token, ORCID OAuth, or author access-request click. Everything open-access (including **CNGBdb** `ftp.cngb.org/pub/CNSA/` and **CyVerse `iVirus`** `data.cyverse.org/dav-anon/`, both confirmed anonymous HTTPS!) is already auto-chained detached onto `tower` across Batches 1–7.*

---

## 1. JGI Data Portal / IMG/M / MycoCosm / PhycoCosm *(PI already has JGI SSO login)*

- **Start here:** `https://data.jgi.doe.gov/` (Profile icon → **Copy My Session Token**, or Globus)
- **License note:** JGI Data Utilization Policy (`redistributable: false`, `DEC-51`).
- **Targets (in priority order):**
  1. **`LR_jgi_batch.tsv.gz`** — `532` unrestricted JGI long-read metagenome assemblies (`257.48 Gbp` `*measured*`; piggyback `ViralRecall`/`GVClass` giant-virus extraction on the same assemblies).
  2. **`IMG/M` Post-`2024-02-06` Unrestricted MQ/HQ MAG Delta** (`img.jgi.doe.gov/m` → `Find Genomes` → `Metagenome Bins` → `Public=Yes`, `Unrestricted`, `HQ/MQ`, `Add Date > 2024-02-06`): `~10k–20k` net-new species (`~25–50 Gbp` `*estimate*`).
  3. **`IMG/M` Post-`2023-11-10` Unrestricted SAG Delta** (`Genome Search` → `Single Cell (SAG)`, `Unrestricted`, `Add Date > 2023-11-10`): `~1.5k–3k` SAGs (`~2–4 Gbp` `*estimate*`).
  4. **`MycoCosm` (`https://mycocosm.jgi.doe.gov/mycocosm/home`)** — early-diverging non-Dikarya fungi (`1KFG`/`ZyGoLife`/*Neocallimastigomycota* anaerobic fungi + unmined *Mycodnaviridae* GEVEs; GFF3 + FASTA for `C-60`).
  5. **`PhycoCosm` (`https://phycocosm.jgi.doe.gov/phycocosm/home`)** — `~294–300` microalgal genomes (incl. *100 Diatom Genomes* + algal GEVEs; GFF3 + FASTA for `C-60`).
  6. **JGI Proposal `505025`** — `815` chemosynthetic marine nematode symbiont MAGs.
- *(Do NOT file a JGI tape-restore for `GVMAGs v1` — superseded by `GV1` `GVMAGs V2`.)*

---

## 2. Biosino NODE / `eLMSG` *(1 free account unlocks 3–4 flagship catalogs via `lftp` FTPS)*

- **Register / login:** `https://www.biosino.org/node`
- **Bulk FTPS pull directly on `tower` once registered:**
  ```bash
  lftp -e "set ftp:ssl-force true; set ssl:verify-certificate no; open ftps://fms.biosino.org:2122; user <email> <password>; ls /Public/"
  ```
- **Targets:**
  1. **`CPR1` (`PBGD` + `PPGD`, *Nat. Microbiol.* Nov 2026):** `https://www.biosino.org/node/project/detail/OEZ00021651` — **`24,523` *Patescibacteria* (CPR) genomes (`13,242` HQ) + `18,876` patesciphages (`12,867` vOTUs)** (`*published*`).
  2. **`GDSV 2.0` (Global Deep-Sea Sediment Virome):** `https://www.biosino.org/node/project/detail/OEP003844` — **`554,664` deep-sea sediment DNA vOTUs** from 159 cores (`*published*`).
  3. **`DSR-RNA` (Deep-Sea Sediment RNA Virome):** `https://www.biosino.org/node/project/detail/OEP002537` — **`85,059` VLP-purified benthic RNA vOTUs (`98.28%` novel) + `1,463` complete benthic RNA genomes** (`*published*`, for `C-50`).
  4. **`MEER` Mariana Trench Hadal (`OEP004067` / `eLMSG LMSG_G000037053.1–44616.1`):** `https://www.biosino.org/node/project/detail/OEP004067` — **`7,564` hadal SRGs (`89.4%` novel)** *(only needed if CNGBdb mirror `CNP0004890` in Batch 7 lacks SRG FASTAs)*.
  5. *Bonus on same `eLMSG` login:* **`MBAMD-ISO`** (`652` acid-mine-drainage isolates).

---

## 3. Dryad *(ORCID login → API token unlocks all Dryad v2 datasets directly from `tower`)*

- **Login / token:** `https://datadryad.org/` (Sign in with ORCID → Account → API Application / Token — OR click Download in browser and paste presigned S3 link).
- **Targets:**
  1. **`PEONY 1.3` + `PIGEON 3.0`** (`https://datadryad.org/dataset/doi:10.5061/dryad.prr4xgz35`): `PEONY1.3_Nov26.fa.gz` (`8.05 GB`, SHA-256 `0107a228...`) + `PIGEON3.0.fa.gz` (`9.13 GB`, SHA-256 `f12f09a5...`) + metadata CSV — **`>1.35M` encapsidated VLP vOTUs** (`#1` phage portal target, `CC0-1.0`).
  2. **Shipworm Symbiont MAGs:** `10.5061/dryad.ksn02v7jd` (`90` MAGs).
  3. **Targeted Niche MAGs:** Boreal river Hg-methylators (`10.5061/dryad.9zw3r22mz`, `850` MAGs) & Kelp collapse (`10.5061/dryad.79cnp5j78`).

---

## 4. CNCB / NGDC Portal-Only Items *(only if Batch 7 open FTP probe misses them)*

- **Portal:** `https://ngdc.cncb.ac.cn/`
- **Targets:**
  1. **`GCSV` (Global Cold Seep Virome):** `https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA048260` — **`193,465` cold-seep vOTUs (`>99%` novel, huge archaeal virus reservoir)**.
  2. **`CnGVR`:** `OMIX007646` (Chinese gut viral ref).
  3. **`FUTM1`:** Check `Dataset 2.xlsx` from Batch 3 `figshare:30883646` for CNCB/NODE accession.

---

## 5. One-Click Author Access Request & WAF-Blocked Browser Grabs

- **Zenodo `20336193` (`EUK-DSDNA-EVE`):** `https://zenodo.org/records/20336193` — `781,111` eukaryotic Maverick/Polinton/EVE loci (`CC BY 4.0`, `access_right: restricted`). Log into Zenodo, click **Request access**, and forward the approval `?token=...` URL when emailed.
- **Browser-only WAF sites (small files, zero login, just AWS/Cloudflare WAF blocking `curl`):**
  - `GV-MYCODNA` (*Mycodnaviridae*, first fungal giant virus family): `https://doi.org/10.7302/p5ay-ng61` (UMich Deep Blue Data, Cloudflare WAF)
  - `MBGG1` (Migratory bird `5,823` MAGs): `https://figshare.com/s/8699370b85d1875209b2` (AWS WAF)
  - `ET1_1KP` (1,000 Plant/Algal Transcriptomes): CyVerse WebDAV (Cloudflare Turnstile)
  - `GCHVD1` (`36,860` coral holobiont vOTUs): Mendeley Data `rkchknnsw5` (wait for public publish / check `PRJNA1310261` & `PRJNA1313875`).

---

## 6. Embargo Calendar Watchlist (Auto-Open Later — No Login Needed Once Unembargoed)

- **2027-02-01:** `GBR-SED` (`zenodo:24187765`, `2,698` Great Barrier Reef coral-reef sediment MAGs)
- **2027-09-03:** PFAS bioremediation MAGs (`zenodo:23030372`)
- **2027-12-04:** Karst cave MAGs (`figshare:32209623`)
- **Monitor `TEG` (`51,921` Tibetan MAGs):** `zenodo:17089467` / `24492917` (waiting for authors to upload the MAG FASTA tarball alongside existing tables/PDBs).

# W-40: Unified Single-Subunit & Single-Polypeptide Viral RNAP Super-Panel (`ssrnap`)

Date: 2026-11-24  
Host: `wilkcomp02.ccbb.utexas.edu` (read-only against `/tmp/if2453/r10/root`, release `r26.11.2`, `--stages S1,H1,PL1,P1,C5,V1`)  
Task card: `W-40` (`DEC-04`, `DEC-29`, `DEC-50`, `DEC-51`, `DEC-69`, `DEC-73`, `DEC-74`)

## 1. Two-Namespace Panel Composition (`refdb.yaml` & `examples/rnap/`)

The `ssrnap` super-panel keeps the two evolutionarily distinct catalytic folds strictly separated into two target namespaces (`30` target profiles = `25` `palm:` + `5` `dpbb_single:`) guarded by `16` outgroup decoy profiles (`46` profiles total), with per-model thresholds in [`examples/rnap/ssrnap_cutoffs.tsv`](../../examples/rnap/ssrnap_cutoffs.tsv) and a runnable companion check in [`examples/rnap/ssrnap_panel.sh`](../../examples/rnap/ssrnap_panel.sh) (leaving the locked `PF00940` release check `examples/rnap/run.sh`, `expected.tsv`, and `recorded/*.tsv` untouched).

| Namespace / Role | Source Pack (pinned in `refdb.yaml`) | Models | Cutoff | SHA256 / Provenance |
|---|---|---:|---|---|
| `palm:` (Fold A — right-handed T7 / Pol I RRM) | `pfam` `38.2` (`CC0`, redistributable) | `7` | Pfam `GA` | `PF00940.25` (`RNA_pol`, GA `22.0`), `PF14700.12` (`RPOL_N`, `30.1`), `PF29426.1` (`Phage_RNApol`, `27.0`), N4 `vRNAP` quartet `PF21769.5` (`29.1`), `PF21867.4` (`28.5`), `PF21894.4` (`28.1`), `PF21624.3` (`29.4`) |
| `palm:` (Fold A) | `vogdb_rnap` (`vog227`, `CC BY 4.0`, redistributable) | `7` | `30.0` bits, dom `20.0`, `iE <= 1e-5` | `vogdb_rnap_palm_clean.hmm` (`5,907,002 B`, sha256 `69dacd8555b6d2bcb63f68bd692cce420a01d3bc85b60ddd5c41ef99c5b61c25`; `VOG10353`, `VOG27844`, `VOG44798`, `VOG57737`, `VOG67840`, `VOG68249`, `VOG72184`) |
| `palm:` (Fold A) | `phrogs_rnap` (`v4`, Zenodo `10.5281/zenodo.23342382`, `CC BY 4.0`, redistributable) | `6` | `30.0` bits, dom `20.0`, `iE <= 1e-5` | `phrogs_rnap_palm_clean.hmm` (`2,608,689 B`, sha256 `5a25ed171396f4e6a05fbe4698b154837010ce7c83aedc92f8a0cc022076a40f`; `phrog_414`, `2051`, `2166`, `2206`, `20599`, `36028`, built via `hmmbuild --seed 42` from `MSA_phrogs.tar.gz` sha256 `1ca216aa...`) |
| `palm:` (Fold A) | `ssrnap_custom_seeds` (`CC0`, redistributable) | `5` | `30.0` bits, dom `20.0`, `iE <= 1e-5` | `seeds_palm.hmm` (`1,242,938 B`, sha256 `81ade5eaf5599ee11fd8d550d717e34e7f9e87feb70beaee82221af555ff0055`; `hmmbuild --seed 42`: `N4_RNAP2_gp15` `60` seqs, `N4_RNAP2_gp16` `60` seqs, `mito_invertron_rpo` `35` seqs, `psychro_autographi_rnap` `40` seqs, `basal_protist_rpot` `15` seqs from `Pfam-A.seed.gz`) |
| `dpbb_single:` (Fold B — single-polypeptide / split DPBB) | `pfam` `38.2` (`CC0`, redistributable) | `1` | Pfam `GA` | `PF05183.18` (`RdRP` / `QDE-1`-like, GA `23.8`) |
| `dpbb_single:` (Fold B) | `ssrnap_custom_seeds` (`CC0`, redistributable) | `4` | `30.0` bits, dom `20.0`, `iE <= 1e-5` | `seeds_dpbb.hmm` (`1,911,040 B`, sha256 `4ff13b69ca22fe3c33b6a364fbbabca58bbfadd932a50992d2f751889530e8eb`; `hmmbuild --seed 42`: `crass_vrnap_dpbb` `60` seqs, `yono_dpbb` `9` seqs, `pgkl2_orf6_dpbb` `17` seqs, `pgkl2_orf7_dpbb` `12` seqs) |
| `decoy_pol1:` (Pol I outgroup) | `pfam` `38.2` | `1` | Pfam `GA` | `PF00476.27` (`DNA_pol_A`, GA `34.3`) |
| `decoy_cell_rnap:` (cellular DPBB outgroup) | `pfam` `38.2` | `15` | Pfam `GA` | `RNA_pol_Rpb1_1..7` (`PF04997`, `PF00623`, `PF04983`, `PF05000`, `PF04998`, `PF04992`, `PF04990`) + `RNA_pol_Rpb2_1..7,45` (`PF04563`, `PF04561`, `PF04565`, `PF04566`, `PF10385`, `PF04567`, `PF00562`, `PF04560`) |

Note on `NCBIfam 20.0`: scanning all `18,950` `hmm_PGAP.HMM.tgz` profiles confirmed `0` phage/organellar single-subunit `ssRNAP` models (all `RNA polymerase` models in `NCBIfam 20.0` are cellular `rpoA/B/C/D/E/K/L/N/P/omega` subunits or sigma factors), so `0` NCBIfam models were added. Note on `PHROGs v4`: `HMM_phrog.tar.gz` contains HHsearch `.hhm` files rather than HMMER3 profiles; HMMER3 profiles were therefore built deterministically (`hmmbuild --seed 42`) from the sibling `MSA_phrogs.tar.gz` (`md5: ba9f677b7d56c53b49e297405d34a921`, `sha256: 1ca216aa830ebb94c925d2733c14d5870d3922bbe7a6eddffc66122cfc945923`) in Zenodo record `10.5281/zenodo.23342382`.

## 2. Benchmark on `r26.11.2` (`--stages S1,H1,PL1,P1,C5,V1`)

All scans ran read-only (`--bind /tmp/if2453/r10:/tmp/if2453/r10:ro`) on `wilkcomp02` (2026-11-24) inside `bio-box_faec1f4.sif` capped to `<= 4` CPUs (`--cpus 4 --threads 4`), `nice 19`, `MemoryMax=32G`, `MemorySwapMax=0`, across `14` deterministic sub-step shards (`full_task_0000`..`full_task_0013`, each `< 200` CPU-min; unit `w40-ssrnap-1` shards `0..7` stopped cleanly at `261.5` CPU-min / peak `8.5 GB`; unit `w40-ssrnap-2` per-task subprocesses shards `8..13` `131.8` CPU-min total / peak `4.2 GB`; merge + loci + oracle + `LPL1`/`O2` verification `435.7 s` wall; total `/tmp/if2453/w40` disk footprint `44 MB` `< 2 GB`).

| Metric | Value | Label |
|---|---:|---|
| Panel models (`targets` / `decoys`) | `46` (`25` `palm:` + `5` `dpbb_single:` + `1` `decoy_pol1:` + `15` `decoy_cell_rnap:`) | *measured* |
| Scan sub-steps | `14` (`full_task_0000`..`full_task_0013`) | *measured* |
| Total `protein_hits.parquet` rows / distinct proteins | `77,291` / `29,321` | *measured* |
| `vetoed_by_decoy = True` rows (`decoy_cell_rnap` vetoing cellular `rpoB`/`rpoC` cross-hits) | `1,525` (`1,465` distinct proteins on `dpbb_single:pgkl2_orf6_dpbb`) | *measured* |
| Clean target proteins (`NOT is_decoy AND best_on_overlap AND NOT vetoed_by_decoy`) | `4,922` `palm:` + `4,168` `dpbb_single:` (`9,090` total) | *measured* |
| Total `orf_hits.parquet` rows / distinct ORFs | `57,548` / `14,095` | *measured* |
| Total `loci_summary.parquet` loci (`--max-gap-genes 2`) | `10,279` (`V1: 5,086; S1: 3,555; H1: 1,561; C5: 38; PL1: 32; P1: 7`) | *measured* |
| `split_model_pair` loci | `363` (`V1: 271; S1: 67; H1: 24; PL1: 1`) | *measured* |
| N4 `vRNAP` (*Schitoviridae*) recovery (`S1` / all 6 stages) | `S1: 400` proteins / `496` contigs (`V1: 126/123`, `C5: 28/30`, `PL1: 3/3`, `H1: 3/1`) — exceeds `253+` phages missed in `2026-11-09_ssrnap_workflow` | *measured* |
| N4 heterodimeric `RNAP2` (`gp15`+`gp16`) paired loci (`--max-gap-genes 2`) | `267` loci (`S1: 211`, `V1: 55`, `PL1: 1`) | *measured* |
| `DNA_pol_A` (`PF00476.27`) false positives (`16,183` decoy hits) | `0` | *measured* |
| Cellular `rpoB`/`rpoC` false positives in `P1`/`C5` and `S1` | `0` (`1,465` cellular proteins vetoed by `decoy_cell_rnap`; remaining `H1`/`V1`/`PL1` multi-domain overlaps are viral jumbo-phage RNAPs) | *measured* |
| Twilight-zone viral/plasmid dark ORFs shortlisted (`twilight_dark_orfs.tsv`, no Pfam architecture) | `1,872` (`1,506` in `700–1,300` aa; `366` in `2,000–4,000` aa) | *measured* |

### Stage × Namespace Breakdown (Clean Non-Decoy `best_on_overlap AND NOT vetoed_by_decoy` ORFs / Contigs)

| Stage | `palm:` ORFs / Contigs | `dpbb_single:` ORFs / Contigs |
|---|---:|---:|
| `S1` (phages) | `3,746` / `2,844` | `155` / `135` |
| `H1` (huge phages) | `1,301` / `1,277` | `285` / `272` |
| `V1` (viromes) | `1,185` / `1,011` | `4,294` / `3,979` |
| `PL1` (plasmids) | `9` / `7` | `25` / `24` |
| `C5` (complete prokaryotes) | `42` / `38` | `0` / `0` |
| `P1` (prokaryote species reps) | `8` / `6` | `0` / `0` |

## 3. Canonical Controls Recovery (`100%`)

### 3.1 `r26.11.2` `S1` Controls (all `vetoed_by_decoy = False`, *measured*)

| Control | Contig | `protein_uid` | Namespace | Top Model & Sequence Bitscore |
|---|---|---|---|---|
| T7 `gp1` | `NC_001604` | `2eef240afefddcd2c8ac7660eb50e384` | `palm` | `palm:phrog_414` (`1,138.9`), `palm:VOG10353` (`1,065.4`), `palm:PF00940.25` (`367.7`) |
| `phiKMV` | `NC_005045` | `eb8377fda717b96aac266fda63a2b096` | `palm` | `palm:phrog_414` (`884.8`), `palm:VOG10353` (`747.0`), `palm:PF00940.25` (`340.5`) |
| `VSW-3` | `KX066068` | `42c5546422cdf8424b5a0db0d0bdc6f8` | `palm` | `palm:psychro_autographi_rnap` (`1,333.5`), `palm:phrog_414` (`905.4`), `palm:PF00940.25` (`304.2`) |
| `Syn5` | `NC_009531` | `bf8803f4abe8a136c0d132ad2f31882b` | `palm` | `palm:VOG10353` (`799.3`), `palm:phrog_414` (`778.0`), `palm:PF00940.25` (`252.5`) |
| `KP34` | `NC_013649` | `46ecf9a9d959fa66eca4cff738575523` | `palm` | `palm:phrog_414` (`909.0`), `palm:VOG10353` (`795.0`), `palm:PF00940.25` (`295.3`) |
| N4 `vRNAP` (`gp50`) | `NC_008720` | `1e396b2d4e3d7a5280a6835bcf18e10d` | `palm` | `palm:VOG27844` (`4,076.9`), `palm:phrog_2206` (`3,964.3`), `palm:PF21894.4` (`215.9`), `palm:PF21769.5` (`174.7`) |
| N4 `RNAP2` (`gp15`) | `NC_008720` | `a17b8fbf735e80abf6546c9a641ee2fb` | `palm` | `palm:N4_RNAP2_gp15` (`509.3`), `palm:phrog_2166` (`397.9`) |
| N4 `RNAP2` (`gp16`) | `NC_008720` | `eda1349c1f642a9294c36fda6fd0793a` | `palm` | `palm:N4_RNAP2_gp16` (`791.2`), `palm:phrog_2051` (`573.4`) |
| `phi14:2` `vRNAP` (`6VR4`) | `NC_021806` | `835ad638bcf19577465c5c2792dcb759` | `dpbb_single` | `dpbb_single:crass_vrnap_dpbb` (`3,231.2`) |
| `SPβ` `YonO` | `NC_001884` | `82bbc9c36fa5629026f5ea0a668ba1be` | `dpbb_single` | `dpbb_single:yono_dpbb` (`1,581.5`) |

On `NC_008720`, `bio-box search panel` (`--max-gap-genes 2`) automatically paired N4 `gp15`+`gp16` into a single `2,072` bp heterodimeric locus (`palm:N4_RNAP2_gp15(+)>palm:N4_RNAP2_gp16(+)`) alongside the separate `10,503` bp `palm:VOG27844(-)` `vRNAP` locus.

### 3.2 Staged `LPL1` & `O2` FASTA Verification (*measured*, `hmmsearch` `-Z 1e8 --domZ 1e8 -T 30`)

Translating `/tmp/if2453/s2/inputs/LPL1/ncbi_eutils_c59a_lpl1.fna.gz` (`2,709` ORFs >= 80 aa -> `112` hit contigs) and `/tmp/if2453/s2/inputs/O2/ncbi_refseq_plastid_r238_o2.fna.gz` (`134` ORFs -> `7` hits) confirmed 100% recovery of every mitochondrial invertron and cytoplasmic killer-plasmid control:

| Control | Accession | Namespace | Top Model & Sequence Bitscore |
|---|---|---|---|
| *Neurospora* `kalilo` | `X52106.1` | `palm` | `mito_invertron_rpo` (`1,058.0`; Pfam `RNA_pol` `220.1`) |
| *Podospora* `pAL2-1` | `X60707.1` | `palm` | `mito_invertron_rpo` (`637.8`; Pfam `RNA_pol` `204.4`) |
| *Neurospora* `maranhar` | `X55361.1` | `palm` | `mito_invertron_rpo` (`620.5`; Pfam `RNA_pol` `262.7`) |
| *Claviceps* `pClK1` | `X15648.1` | `palm` | `mito_invertron_rpo` (`225.8`; Pfam `RNA_pol` `72.6`) |
| Maize `S-2` | `J01426.1` | `palm` | `mito_invertron_rpo` (`113.2`; Pfam `RNA_pol` `55.5`) |
| *K. lactis* `pGKL2` `ORF6` + `ORF7` | `X07776.1` | `dpbb_single` | `pgkl2_orf6_dpbb` (`1,714.7`) + `pgkl2_orf7_dpbb` (`195.9`) |
| `pPE1B` / `pPac1-1` / *K. phaffii* | `AJ278986.2`, `AM180622.1`, `MG491503.1` | `dpbb_single` | `pgkl2_orf6_dpbb` (`1,516.1`–`1,624.9`) + `pgkl2_orf7_dpbb` (`169.5`–`188.6`) |

## 4. 1:1 `hmmsearch` Binary Oracle Concordance (*measured*)

Standalone `hmmsearch` (`-Z 1e8 --domZ 1e8`) was run on `50` sampled `S1`/`H1`/`V1` proteins plus all `10` canonical controls across both custom seed packs and the Pfam GA subset (`PF00940.25`, `PF21769.5`, `PF05183.18`):

| Oracle Subset | Standalone `hmmsearch` Domains | Panel Domains | Exact 1:1 Match (`protein_uid`, `model`, `env_start`, `env_end`) |
|---|---:|---:|---|
| `seeds_palm.hmm` (`5` custom Fold A profiles) | `87` | `87` | `True` |
| `seeds_dpbb.hmm` (`4` custom Fold B profiles) | `29` | `29` | `True` |
| Pfam GA subset (`PF00940`, `PF21769`, `PF05183`) | `39` | `39` | `True` |
| **Total** | **`155`** | **`155`** | **`True` (`155 / 155`)** |

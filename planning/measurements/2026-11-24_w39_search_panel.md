# W-39: One-Pass Multi-HMM Panel Search (`bio-box search panel`) & `prok-rt` Super-Panel Benchmark

Date: 2026-11-24  
Host: `wilkcomp02.ccbb.utexas.edu` (read-only against `/tmp/if2453/r10/root`, release `r26.11.2`)  
Task card: `W-39` (`DEC-04`, `DEC-29`, `DEC-51`, `DEC-71`, `DEC-73`, `DEC-74`)

## 1. `anchors.hmm` (`TniQ + Cas_Csy3`) 1:1 Parity Check on `r26.11.2`

Command (*measured*, `wilkcomp02`, 2026-11-24, unit `w39-anchors-7`, 4 CPUs, `nice 19`, `MemoryMax=32G`):
```bash
python -m biobox.cli search panel --root /tmp/if2453/r10/root --release r26.11.2 \
  --profile profiles/wilkcomp02.yaml --hmm /tmp/if2453/w39/anchors/anchors.hmm \
  --cutoff gathering --positive PF06527.17=NZ_AP026781.1:426340-427525 \
  --positive PF09615.16=NZ_AP026781.1:423377-424433 --negative none \
  --min-models 2 --max-gap-genes 9 --max-gap-bp 100000 --neighbors 10 \
  --cpus 4 --threads 4 --out /tmp/if2453/w39/anchors/out
```

| Metric | Value | Label |
|---|---:|---|
| Run ID | `ac558dba3280ba8e578a1cf812638017` | *measured* |
| Wall time (full `r26.11.2` pack, 2 HMMs, 4 CPUs) | `651.4 s` (`10.9 min`) | *measured* |
| `TniQ` (`PF06527.17`) hit proteins / ORFs | `3,368` / `4,980` | *measured* |
| `Cas_Csy3` (`PF09615.16`) hit proteins / ORFs | `828` / `1,110` | *measured* |
| Called co-occurring `TniQ + Cas_Csy3` loci (`loci_summary.parquet` & `loci/*.tsv`) | `152` loci across `149` contigs | *measured* |
| 1:1 SQL parity against `r26.11.2` `l1_domains` `TniQ` + `Cas_Csy3` (`<= 10` slots) | `149 / 149` contigs (`0` only-SQL, `0` only-panel) | *measured* |
| Stage breakdown of `152` loci | `C5: 87, P2: 17, P4: 13, P5: 13, P6: 9, PL1: 7, P3: 5, P1: 1` | *measured* |
| `results archive` + `get` round-trip (`W-28`/`W-29`) | `spooled`, `5,027,840 B`, `167 / 167` files restored | *measured* |

## 2. Unified 9-Pack `prok-rt` Super-Panel Assembly (`refdb.yaml`)

All 9 prokaryotic/viral/ART-RT HMM packs were pinned in `refdb.yaml` with exact version/DOI/URL/local path, license, `redistributable` flag, and `sha256` checksums (`388` total profiles = `375` targets + `13` outgroup/composition decoys):

| Pack | Models | Namespace(s) | License / `redistributable` | Source / Notes |
|---|---:|---|---|---|
| 1. `lab_art_rt` | `15` (`4` targets + `8` decoys + `3` components) | `art`, `decoy_art`, `art_comp` | unpublished / `false` | `both.hmm` (`min_dom_frac=0.45` on targets; `*_REV`/`*_PERM` routed to `--decoy-hmm`) + `our_components.hmm` |
| 2. `spire_retrons_2026` | `30` (`1` general + `29` types) | `spire` | MIT / `true` | Zenodo `10.5281/zenodo.20409741` (`toro_general.hmm`, `toro_types.hmm`, `32.0`-bit cutoff) |
| 3. `toro_ug_abi_retrons` | `0` standalone `.hmm` (`DRT2` zip pinned) | — | CC BY 4.0 / `true` | Zenodo `10.5281/zenodo.13257156` (`DRT2-analysis-v1.0.0.zip`, `sha256 c98d18d8...`; `UG1..28` & `DRT_1..9` covered by `padloc`/`df`/`myrt`) |
| 4. `padloc_rt` | `172` | `padloc` | MIT / `true` | `padloc-db v2.0.0` (`sha256 29c3ada5...`) RT/retron/DRT/Abi HMMs + `hmm_meta.txt` thresholds |
| 5. `defensefinder_rt` | `66` | `df` | GPL-3.0 / `false` | `defense-finder-models v3.1.0` `Retron_*`, `DRT_*`, `All_UG`, `AbiA/K/P2` |
| 6. `myrt` | `45` | `myrt` | unlicensed / `false` | `RVT-All.hmm` (`sha256 74556bd2...`, `30.0` bits, `i_evalue <= 1e-5`) |
| 7. `roux_dgr` | `14` | `dgr` | BSD-3-Clause-LBNL / `true` | `dgr_scripts v1.0` (`sha256 634bbd80...`, `Custom_DGR_RT.hmm`) |
| 8. `vogdb_rt` | `10` | `vog` | CC BY 4.0 / `true` | `VOGdb vog227` (`sha256 d3f9d239...`, 10 RT-annotated VOGs) |
| 9. `pfam 38.2` + `ncbifam 20.0` | `36` (`15` Pfam + `5` `decoy_rdrp` + `16` NCBIfam) | `pfam`, `decoy_rdrp`, `ncbifam` | CC0 / PD / `true` | Pfam `GA` (`PF00078`, ART partners, `RdRP_1..4`/`Flavi_NS5` decoys) + NCBIfam `TC` (`TIGR04416`, `NF038233`, RT models) |

## 3. `prok-rt` Super-Panel Benchmark on `r26.11.2` (`--stages S1,H1,PL1,P1,C5`)

Run on `wilkcomp02` (2026-11-24, unit `w39-prokrt-1`, 4 CPUs, `nice 19`, `MemoryMax=32G`, `39` sub-step shards each `< 200` CPU-min, scan peak RSS `2.89 GB` *measured*, merge+verify `263.9 s` *measured*):

| Metric | Value | Label |
|---|---:|---|
| Scan sub-steps | `39` (`full_task_0000`..`full_task_0038`), each `< 200` CPU-min | *measured* |
| Scan peak RSS (`w39-prokrt-1` across all 39 scan sub-steps) | `2.89 GB` (`2,885,419,008 B`, `< 4 GB` target) | *measured* |
| Total `protein_hits.parquet` rows / distinct proteins | `418,307` / `50,478` | *measured* |
| `vetoed_by_decoy = True` hits | `2` | *measured* |
| Total `orf_hits.parquet` rows / distinct ORFs | `1,735,589` / `111,183` | *measured* |
| Total `loci_summary.parquet` clusters | `98,744` (`C5: 45,225; PL1: 33,002; S1: 9,135; H1: 6,663; P1: 4,719`) | *measured* |
| Global diagnostics (`cas_less_crispr` / `split_model_pair`) | `288` / `616` | *measured* |
| `H1` huge-phage diagnostics (`loci` / `cas_less` / `igr>=500` / `bridged>=500` / `split_pair`) | `6,663` / `5` / `928` / `1,638` / `7` | *measured* |
| `PL1` plasmid diagnostics (`loci` / `cas_less` / `igr>=500` / `bridged>=500` / `split_pair`) | `33,002` / `234` / `17,805` / `22,112` / `396` | *measured* |
| `S1` phage diagnostics (`loci` / `cas_less` / `igr>=500` / `bridged>=500` / `split_pair`) | `9,135` / `1` / `825` / `2,022` / `19` | *measured* |
| Decoy-filtered novel/unclassified non-`PF00078` RT proteins (`best_on_overlap AND NOT vetoed_by_decoy`) | `45,489` | *measured* |
| Standalone `hmmsearch` binary oracle concordance (25 sampled `S1` `PF00078.33` proteins) | `25 / 25` (`exact_match: True`) | *measured* |

### Recovery of Known `S1` ART-RT Operons

`bio-box search panel` recovered intact 6-ORF huge-phage ART-RT operons in `S1` (`huge_phage_500`, `huge_phage_666`, `huge_phage_375`, `huge_phage_694`, `huge_phage_697`, `huge_phage_831`, `huge_phage_305`, `huge_phage_371`, `huge_phage_630`, `huge_phage_920`) with exact syntenic order `art_comp:DUF8567(+)>art_comp:PG_binding_1(+)>art:DISTORF(+)>art:NTERM+spire:retron_reference+art:CTERM(+)>art:PROXORF(+)>art_comp:DUF8549(+)` and `max_adjacent_igr_bp = 1,096–1,285 bp`, plus *Staphylococcus* phage retron loci (`OR845005`, `OR844994`, `OM439675`, `OM439669`, `MK075001`). Note: `LPJP1` (`PX530545`, post-14Apr2025 INPHARED delta) is not yet in `r26.11.2` `S1` (staged for `V3`/`SAT1` in `r26.11.4` under `DEC-71`).

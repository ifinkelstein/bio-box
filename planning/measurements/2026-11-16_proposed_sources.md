# Proposed source entries for cards C-39, C-40 and C-43

Date: 2026-11-16. Step (1) of the three cards, done by the orchestrator at the PI's word ("do your suggestion re
c-39, c-40, c-43"): the entries proposed for `sources.yaml`. **Nothing here is approved and nothing was fetched**:
the pages below were read through web searches and page reads on 2026-11-16, and every size, count, checksum and
licence is the provider's own statement (*published*), to be checked against the files when they arrive. An entry
goes into `sources.yaml` only with the PI's approval (AGENTS §8).

## A. Ready to approve: exact files, open licences

```yaml
  - name: uhgv
    what: Unified Human Gut Virome v1.0, the vOTU representatives and their tables; card C-39
    url: https://zenodo.org/records/17402089/files/{file}?download=1   # Zenodo's file form; confirm at the first request
    doi: 10.5281/zenodo.17402089   # published 2025-10-24
    files:
      - {file: votus_full.fna.zst, size: 964.7 MB, md5: 9ad996fea08c6087a452194f0c2f9d33}   # 168,536 vOTU representatives
      - {file: votus_metadata_extended.tsv.zst, size: 14.0 MB, md5: 6702355349e933fb605e62a8bbaa50a9}
      - {file: source_biosample_metadata.tsv.zst, size: 152.1 kB, md5: 9934204a79cdac60ce9681a6233282b6}
      - {file: host_range_breadth.tsv.zst, size: 490.6 kB, md5: 25960f854fa411956ef1353213133bfc}
    not_taken: uhgv_full.fna.zst (1.7 GB, all 873,995 genomes) and the protein files (we call our own genes)
    license: CC BY 4.0 (the Zenodo record)
  - name: gsv_atlas
    what: Global Soil Virus Atlas (Graham et al. 2024), one archive of six files; card C-39
    url: https://g-83fdd0.1beed.03c0.data.globus.org/cartuids/35846_b8b9149f2f59ff950549b5c89b8b35e2.tar
    readme: https://datahub.pnnl.gov/sites/default/files/2023-12/README_Global%20Soil%20Virus%20%28GSV%29%20Atlas.txt
    doi: 10.25584/2229733   # PNNL DataHub
    contents: 616,935 uncultivated viral genomes; the 49,649 that passed quality control (38,508 vOTUs); proteins; gene records; the table of 2,953 soil samples
    size: not stated on the page; read from the response header before the download
    license: CC0 (the DataHub page)
  - name: spire_tables
    what: SPIRE v1 tables only (no sequence); card C-43, and the sample choice of card C-42
    url: https://black.embl.de/~fullam/spire/metadata/{file}
    files:
      - {file: README, size: 5.7 kB, md5: 51152f40061eeab04fabbe18d37c22ac}
      - {file: spire_v1_genome_metadata.tsv.gz, size: 92 MB, md5: a4ba393df92e9d1a35b8b996a69b39fd}
      - {file: spire_v1_cluster_metadata.tsv.gz, size: 2.1 MB, md5: 3f572296ed0d5918b18916b6a9099618}
    also: https://black.embl.de/~fullam/spire/microntology/ena_drop.v_2_0.microntology.v_1_0.annotations.v_0_2.tsv.gz   # habitat labels, 1.7 MB, md5 62056caa199fa3aca9b019f400149527
    license: none stated on the downloads page; the paper is open access. Tables only until the terms are found
  - name: ncbi_metagenome_assembly_summary
    what: NCBI's table of GenBank metagenome assemblies (no sequence); card C-43
    url: https://ftp.ncbi.nlm.nih.gov/genomes/genbank/metagenomes/assembly_summary.txt   # the path NCBI's FTP guide describes; confirm at the first request
    readme: https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
    license: INSDC (public)
```

## B. Not ready: what is missing for each

| Source | Card | What is missing |
|---|---|---|
| IMG/VR v4.1, high-confidence genomes | C-39 | JGI's portal (`genome.jgi.doe.gov/portal/IMG_VR`, the folder whose name ends `_7.1`) answered 503 to an anonymous read today and a download needs a JGI login and acceptance of JGI's data usage policy. This is the PI's account: either the PI downloads the nucleotide FASTA and the table to the pod, or says how an agent may. Local only, `redistributable = false` |
| IMG/PR | C-40 | the same, at `genome.jgi.doe.gov/portal/IMG_PR` |
| Gut plasmids of Yu, Fogarty and Eren 2024 | C-40 | the record found (`doi:10.6082/cx5ay-1ev28`, CC BY 4.0) holds the article and a 35.8 MB supplement; the plasmid sequences are in Zenodo records named in the article's data availability statement, which were not read today. One page read of the open-access article gives them |
| PlasmidScope | C-40 | its site (`plasmid.deepomics.org`) was not read: the download page, the file per source repository and the licence |
| VIRE | C-39 | only its paper was found; its tables' location was not looked for (overlap question only) |
| MGnify assemblies table | C-43 | the listing route (MGnify's API or ENA's portal API) was not read today |
| IMG/M and GOLD sample tables | C-43 | JGI again: unreachable anonymously today |
| SPIRE per-sample assemblies | C-42 | the downloads page offers assemblies by study; the URL form is in its README, which is in block A |

## C. What the agents of the three cards do once block A is approved

- C-39: fetch `uhgv` and `gsv_atlas`, check the checksums, make the inputs of the virus stage and the overlap
  tables; IMG/VR waits for the PI.
- C-43: fetch the two tables, read SPIRE's README, and come back with the entries for MGnify's and JGI's tables.
- C-40: nothing in block A is a plasmid set, so it starts with the two page reads of block B (the article's data
  availability; PlasmidScope's download page) and comes back with exact entries; IMG/PR waits for the PI.

# CRISPR-Cas systems carried by Tn7-like transposons: a worked example and a release check

The second worked example of `GUIDE.md` (the first is the T7-like RNA polymerase, `GUIDE.md` sections 3 to 7). It
runs the guide's commands (`search hmm`, `nbr summarize`, `nbr flanks --figure`, `find`, `seq`, `sql`) on a system
whose answer is published, says where bio-box reproduces it and where it cannot, and is kept as a check that an
agent reruns on every milestone release (section 6). State of 2026-11-13, release `rN.3`. Numbers carry a label:
*measured* (host, date, command), *derived* (formula and inputs) or *estimate*; a number quoted from the paper says
so and names its section, and is never one of ours.

## 1. The plan (written before the run)

**The paper.** Peters JE, Makarova KS, Shmakov S, Koonin EV. Recruitment of CRISPR-Cas systems by Tn7-like
transposons. *PNAS* 2017;114(35):E7358-E7366, doi:10.1073/pnas.1709035114, open access as PMC5584455 (read there,
2026-11-13). Chosen by the PI. It suits a test because its anchors are Pfam families, its main claims are about
which genes sit beside the anchor and which are missing, and its scale (thousands of loci) is within what one
complete-genome release holds.

**The claim reproduced.** A minimal type I-F CRISPR-Cas system travels with a distinct family of Tn7-like
transposons. Quoted from the paper:

- the minimal systems "consist of fused *cas8f* and *cas5f*, *cas7f*, and *cas6f* genes and a short CRISPR array"
  (Abstract), and lack the adaptation module and *cas3*: both the I-F and the I-B variants have "lost the
  adaptation module (*cas1* and *cas2*) and the *cas3* gene" (section "Evolution of the Association Between
  CRISPR-Cas Variants and Tn7-Like Elements");
- "In the vast majority of the loci corresponding to the cas1-less branch, a tnsD(tniQ) gene is located next to
  the cas genes" (section "A Variant of the Type I-F CRISPR-Cas System Is Specifically Associated with a Distinct
  Family of Tn7-Like Elements"; "the I-F section" below);
- the hosts are Gammaproteobacteria: "the branch encompassed 423 sequences from 19 genera of Gammaproteobacteria",
  out of "2,905 Cas7f protein sequences" (the I-F section); organisms named in the Results include *Vibrio
  parahaemolyticus*, *V. fluvialis*, *Photobacterium kishitanii*, *P. leiognathi*, *Pseudoalteromonas rubra*,
  *Shewanella baltica*;
- two further, independent captures of type I-B systems, in cyanobacteria, which "lack adaptation genes and Cas3"
  (section "The cas1-Less Type I-F CRISPR-Cas System Is Mobilized Together with Conserved Transposition Genes");
- "14 complete bacterial genomes containing Tn7-like elements with the I-F CRISPR-Cas" and four complete genomes
  with the I-B elements (section "Analysis of CRISPR Arrays Associated with the cas1-Less I-F Systems");
- three new attachment sites for the I-F elements: downstream of an inosine-5'-monophosphate dehydrogenase (IMPDH)
  gene, upstream of *yciA*, and inside the gene of the signal-recognition-particle RNA; with Tn7's *glmS* and
  Tn6022's *yifB* as the known ones (section "Chromosomal Insertions of the I-F CRISPR-Cas-Associated Elements Show
  Three Recognizable Attachment Sites ...").

How the paper defined its loci (Methods, "Prokaryotic Genome Database and ORF Annotation"): PSI-BLAST from the CDD
profiles of Cas7f (cd09737, pfam09615), TnsA (pfam08722, pfam08721) and TnsD/TniQ (pfam06527) at E 1e-4 against
4,961 complete and 43,599 draft genomes of March 2016; then every ORF "within 10-kb regions up- and downstream of
cas7f genes ... and 20-kb regions up- and downstream of tnsD(tniQ) and tnsA" annotated with 30,953 CDD profiles and
217 Cas profiles. The association itself was judged on protein trees (FastTree, bootstrap): the cas1-less I-F
variants form "a single, monophyletic, strongly supported branch". 10,349 loci held tniQ/tnsD or tnsA (the I-F
section).

**The anchors.** Both of the paper's probes that Pfam 38.2 holds as one model each, so the test runs in both
directions, as the paper's did:

| Anchor | Pfam 38.2 | Gathering cutoff (GA, bits) | Role |
|---|---|---|---|
| TniQ (TnsD) | `TniQ` PF06527.17 | 27.0 | the transposon's target-selection protein, next to the cas genes in the I-F3 elements |
| Cas7f | `Cas_Csy3` PF09615.16 | 25.0 | the I-F Cascade backbone; the paper's I-F probe ("the most conserved component in all systems of this subtype including the minimal variant lacking cas1, cas2, and cas3 genes") |

Searched with `--bits` equal to the model's GA, so a hit is what the release's own Pfam scan (gathering cutoffs)
would call that family; the scan is then an independent check of the search (section 6). TnsA was not used as an
anchor: `Tn7_TnsA-like_N` (PF08722) hits 2,017 phage genes in `rN.3` (*measured*, `sql/q00_models.sql`, below),
a family of phage endonucleases far wider than Tn7.

The other Pfam families read in the windows ("roles" in the tables): Cas8f `Cas_Csy1` PF09611, Cas5f `Cas_Csy2`
PF09614, Cas6f `Cas_Csy4` PF09618; types I-B and I-C: `Cas8b_N`/`Cas8b_C` PF24121/PF24122, `Cas_Csd1` (Cas8c) PF09709, `Cas_Cas7` PF05107, `Cas_Cas5d`
PF09704; Cas6 PF09559, PF01881, PF10040; Cas1 `Cas_Cas1` PF01867; Cas2 PF09827, PF09707; Cas3 PF18019, PF18395,
PF21384 (the I-F Cas2-Cas3 fusion), PF21802, PF22590, PF30455; Cas4 PF01930; Cse1, Cse2, Cas9, Cas10, Cas12a
domains; the transposition proteins TnsA (PF08722, PF08721), TnsB (PF27165, PF29687, PF29688, PF27135), TnsC
(`Tn7_TnsC_Int` PF11426, and `TniB` PF05621, the AAA+ family of TnsC in these elements), TnsD (PF15978), TnsE
(PF18623); and for the attachment sites IMPDH PF00478, 4HBT PF03061 (the thioesterase family of YciA), SIS PF01380
and GATase_6 PF13522 (GlmS), Mg_chelatase PF01078 (YifB). The list is in every SQL file (`roles`).

**The release.** `rN.3`: 27,258 phage genomes (stage S1) and 2,497 complete prokaryote genomes (P1, P2), one per
GTDB species, 10,620,684 distinct proteins, Pfam 38.2 at gathering cutoffs (*measured*, `bio-box releases`,
2026-11-13). The paper's database held about 20 times as many genomes (*derived*: (4,961 + 43,599) / 2,497 = 19.4).

**Windows and levels.** `nbr summarize` at the paper's windows, 20 kb on each side of TniQ and 10 kb of Cas7f, at
the `pfam` and `nr50` levels. The checks and `nbr flanks` use 10 genes on each side (`--neighbors 10`, `--window
10g`): a minimal I-F3 element is about seven genes (*tnsA*, *tnsB*, *tnsC*, *tniQ*, *cas8f-cas5f*, *cas7f*,
*cas6f*), and `nbr flanks` compares positions, so it takes gene windows only (`planning/design/W-25_flanks.md`).
`nbr flanks` at `pfam` (these are prokaryote loci, where most genes carry a domain; on phages `nr50` is the level,
`planning/measurements/2026-11-12_flanks_rN3.md`), loci grouped by the anchor's NR50 cluster. (Changed after the first run, which made one context type per locus: the final runs are in section 2, step 5, and the reason in section 3.)

**What would count as agreement.** In kind: (1) TniQ loci that hold Cas7f hold no Cas1, Cas2 or Cas3 within the
window, and Cas7f loci that hold TniQ likewise; (2) in most such loci the gene next to *tniQ* is the Cas8f-Cas5f
fusion, followed by *cas7f* and *cas6f*, with *tnsA*, *tnsB*, *tnsC* on the other side; (3) the hosts are
Gammaproteobacteria (*Vibrio*, *Photobacterium*, *Pseudoalteromonas*, *Shewanella* and kin); (4) a cyanobacterial
I-B association may be present, without Cas1, Cas2, Cas3; (5) some elements lie beside IMPDH, YciA-like, GlmS or
YifB genes. In number: only as a share. The paper's counts (423 sequences, 19 genera, 14 complete genomes) are of a
database 20 times larger with many genomes per species; ours can only be smaller, and is not expected to match.

**What bio-box cannot do for this analysis today** (stated before the run):

- No CRISPR arrays: no stage of any release has repeat or array annotations (`feature` kind `crispr_array` exists
  in `ids.KIND` but no producer writes it; `PLAN.md` step 4 names MinCED and CRISPRDirection, and no card builds
  it). So no array downstream of *cas6f*, no spacers, no protospacers, no self-targeting; the paper's "more than
  800 spacers" cannot be approached.
- No non-coding RNA in `rN.3` (`releases` prints `-` for ncrna): the SRP-RNA (*ffs*) attachment site cannot be
  seen. Card E-07 (ncRNA scan in shards) and E-09 (the `prok-2` model set) build it.
- No transposon ends: no search for the 5-bp target-site duplication, the 8-bp terminal sequence or the 22-bp
  TnsB-binding sites, so element boundaries are not found; "beside an IMPDH gene" means within the gene window.
- No taxonomy: `contig_dim.tax_order` is null. Hosts are named from the source records' FASTA headers, which is a
  name, not a GTDB call.
- No phylogenies: the paper's evidence is a tree on which the cas1-less branch is monophyletic; bio-box has protein
  clusters (NR90, NR50) and flank trees (UPGMA on gene-order distances), not protein trees with support.
- No enrichment statistic: the counts below are counts (cards F-01 to F-03). Genomes are one per species, so
  redundancy is low, but nothing is down-weighted.
- One complete genome per species, no drafts, no metagenomes: the association may be found in a handful of genomes
  or none.

## 2. Running it

Everything runs on the pod host, in `bash`, as `GUIDE.md` section 2 sets it up; the laptop runs nothing. One script
runs the whole example (section 6 explains it as a check); the steps it runs are listed after it, each one command,
with the placeholders of the runbook: `$TREE` a checkout of current code (here `/tmp/if2453/ci/ex2`, mirrored from the
laptop as `scripts/pod_pytest.sh` does), `$PROFILE` the host profile, `$IMAGE` the test image with matplotlib
(`/tmp/if2453/images/bio-box-test_64f1bf0.sif`), `$DIGEST` its sha256, `$ROOT` the data root (bound read-only),
`$OUT` a new scratch directory, and `$EX` this directory in `$TREE`.

```bash
export TREE=/tmp/if2453/ci/ex2
export PROFILE="$TREE/profiles/wilkcomp02.yaml"
export IMAGE=/tmp/if2453/images/bio-box-test_64f1bf0.sif
export DIGEST=sha256:3d77df6b7f0b894b5371e42eaf077dfda09ebc97b1f78f2ed1d2cfc5fa9fbda4
export ROOT=/tmp/if2453/r10/root
export RELEASE=rN.3
export OUT=/tmp/if2453/ex2/rN.3-b
export EX="$TREE/examples/crispr_cas_neighbourhoods"
bash "$EX/run.sh"
```

So that it outlives the shell, the run of section 3 was started detached with the same values, as a small user unit
of its own (the script itself sleeps and launches; each step is its own capped unit):

```bash
systemd-run --user --unit ex2-driver-b --collect --quiet -p MemoryMax=2G -p MemorySwapMax=0 -p StandardOutput=append:/tmp/if2453/ex2/driver-b.log -p StandardError=append:/tmp/if2453/ex2/driver-b.log --setenv=TREE="$TREE" --setenv=PROFILE="$PROFILE" --setenv=IMAGE="$IMAGE" --setenv=DIGEST="$DIGEST" --setenv=ROOT="$ROOT" --setenv=RELEASE="$RELEASE" --setenv=OUT="$OUT" --setenv=UNIT=ex2b --setenv=THREADS=2 --setenv=USER="$USER" /usr/bin/nice -n 19 bash "$EX/run.sh"
```

`run.sh` waits before each launch until the account's CPU (the sum of `ps -u $USER -o pcpu=`) plus its own threads
(`THREADS`, default 2) fits under the hour's cap (40 by day, 56 at night; AGENTS 7.2), and stops if a `STOP` file
appears in the root. The steps, as it runs them (`$S`, `$N`, `$F`, `$Q` are the shared flags written out in
`run.sh`; `bb` launches one `bio-box` command through `scripts/pod_run.sh`, `tool` one HMMER or Python command in a
capped scope):

1. The anchors, by name, from the Pfam file on the pod (read-only; nothing is downloaded): `hmmfetch -o
   "$OUT/hmm/anchors.hmm" -f "$PFAM" "$OUT/hmm/anchors.keys"` with the keys `TniQ` and `Cas_Csy3`, then one file
   per model. Their GA lines give `--bits`.
2. The positive controls, the TniQ and the Cas7f of one transposon-associated element of *Photobacterium damselae*
   subsp. *damselae* 04Ya311, chromosome 2 (a genus the paper names; section 3): `bio-box find --root "$ROOT"
   --release "$RELEASE" --accession NZ_AP026781.1:426340-427525` (TniQ) and `... NZ_AP026781.1:423377-424433`
   (Cas7f).
3. The two searches: `bio-box search hmm --root "$ROOT" --release "$RELEASE" --profile "$PROFILE" --runs
   "$OUT/runs" --work "$OUT/work" --neighbors 10 --cpus 2 --mem-gb 14 --threads 2 --hmm "$OUT/hmm/TniQ.hmm" --bits
   27.00 --positive NZ_AP026781.1:426340-427525 --out "$OUT/tniq"`, and the same with `Cas_Csy3.hmm`, `--bits
   25.00`, `--positive NZ_AP026781.1:423377-424433`, `--out "$OUT/cas7f"`.
4. Neighbours at the paper's windows: `bio-box nbr summarize "$RUN_TNIQ" --root "$ROOT" --profile "$PROFILE" --runs
   "$OUT/runs" --by pfam,nr50 --top 0 --format tsv --window 20kb --out "$OUT/summary/tniq_20kb.tsv"`, and `--window
   10kb` for the Cas7f run.
5. Loci grouped by their flanks: (a) every TniQ locus as the guide does it, `bio-box nbr flanks "$RUN_TNIQ" --root
   "$ROOT" --profile "$PROFILE" --runs "$OUT/runs" --threads 2 --window 10g --level nr50 --out
   "$OUT/flanks/tniq_nr50"` (one group per anchor NR90 cluster, no figures); (b) every Cas7f locus in one group,
   `nbr flanks "$RUN_CAS7F" ... --window 6g --level nr50 --group-by none --figure --out "$OUT/flanks/cas7f_all_nr50"`;
   (c) the TniQ loci with a Cas domain only: `bio-box sql --root "$ROOT" --release "$RELEASE" --runs "$OUT/runs"
   --threads 2 "$EX/sql/uids_tniq_cas.sql" --out "$OUT/uids_tniq_cas.txt"`, then `bio-box nbr flanks --anchors-file
   "$OUT/uids_tniq_cas.txt" --release "$RELEASE" --tape-version "$TV" --root "$ROOT" --profile "$PROFILE" --threads 2
   --window 10g --level pfam --group-by none --figure --out "$OUT/flanks/tniq_cas_pfam"` (`$TV` is the run's tape
   version, from the `# tape_version` line of `$OUT/tniq/neighbors.tsv`); and the locus table, `bio-box sql ...
   "$EX/sql/q03_tniq_cas_loci.sql" --out "$OUT/tniq_cas_loci.tsv"`.
6. The checks: `bio-box sql --root "$ROOT" --release "$RELEASE" --runs "$OUT/runs" --threads 2 "$EX/sql/checks.sql"
   --out "$OUT/observed_sql.tsv"`, and by base-pair windows without the runs, `bio-box sql --root "$ROOT" --release
   "$RELEASE" --threads 2 "$EX/sql/checks_bp.sql" --out "$OUT/observed_bp.tsv"`.
7. The oracle: `bio-box sql ... "$EX/sql/oracle_uids.sql" --out "$OUT/oracle_uids.txt"`, `bio-box seq --root
   "$ROOT" --release "$RELEASE" --proteins-file "$OUT/oracle_uids.txt" --out "$OUT/oracle.faa"`, `hmmsearch --cpu 2
   -Z 1e8 --domZ 1e8 -T 0 --noali --tblout "$OUT/oracle.tbl" -o "$OUT/logs/hmmsearch.out" "$OUT/hmm/anchors.hmm"
   "$OUT/oracle.faa"`, `python -I "$EX/oracle.py" "$OUT/hmm/anchors.hmm" "$OUT/oracle.tbl" tniq="$OUT/tniq/hits.tsv"
   cas7f="$OUT/cas7f/hits.tsv" > "$OUT/observed_oracle.tsv"`.
8. `observed.tsv` (the three files' rows), then `python -I "$EX/check.py" --expected "$EX/expected.tsv" --recorded
   "$EX/recorded/$RELEASE.tsv" --write-recorded "$OUT/recorded_new.tsv" "$OUT/observed.tsv"`.

The run ids (`$RUN_TNIQ`, `$RUN_CAS7F`) are read from the searches' logs (`search run: <id>`). Three tables of
section 3 and one figure are not written by the script; they came from these commands, each run once in a capped
scope on the pod (`sql/q01_loci_by_domains.sql`, every TniQ and Cas7f gene with the roles in its base-pair window,
was used while writing the checks and is kept for reading):

```bash
bio-box sql --root "$ROOT" --release "$RELEASE" --threads 1 "$EX/sql/q00_models.sql" --out "$OUT/extra/q00_models.tsv"
bio-box sql --root "$ROOT" --release "$RELEASE" --threads 1 "$EX/sql/q02_paper_organisms.sql" --out "$OUT/extra/q02_paper_organisms.tsv"
bio-box sql --root "$ROOT" --release "$RELEASE" --threads 1 "$EX/sql/q04_cas5f_lengths.sql" --out "$OUT/extra/q04_cas5f_lengths.tsv"
bio-box nbr flanks "$RUN_CAS7F" --root "$ROOT" --profile "$PROFILE" --runs "$OUT/runs" --threads 2 --window 6g --level pfam --group-by none --figure --out "$OUT/flanks/cas7f_all_pfam"
```

## 3. What it found on `rN.3`

All *measured* on wilkcomp02, 2026-11-13, by `run.sh` with `OUT=/tmp/if2453/ex2/rN.3-b` (tree `$TREE` mirrored from
the commit that adds this directory, image `bio-box-test_64f1bf0.sif`, sha256 `3d77df6b...fbda4`), unless another
command is named. Run ids: TniQ `ef8d86b8e0e4b756a8b1d63d9eb92aab`, Cas7f `47db59881f7779cb86d7e35984b06705`, pack
`a334c89ccf048904f297de3b2597af22`. `results/check.txt` is the check's own output and `results/observed.tsv` every
value it read. An earlier run (`/tmp/if2453/ex2/rN.3`, other positive controls, an earlier `run.sh`) wrote
byte-identical `hits.tsv`, `hits.faa`, locus tables and 10 kb and 20 kb summaries (`sha256sum`, `diff -r -q`).

**The anchors' families in the release** (`sql/q00_models.sql`, no run; `results/q00_models.tsv`): the release's own
Pfam scan gives a TniQ domain to proteins on 622 genes, all in prokaryote stages (51 proteins in P1, 524 in P2; 574
distinct, *derived*: the 587 TniQ hits below less the 13 without a Pfam row), none in phages;
Cas_Csy3 (Cas7f) 114 proteins on 157 genes, 46 of them in phages (3 distinct proteins); Cas_Csy1 (Cas8f) 78
proteins and Cas_Cas1 1,200 (both summed over the three stages). `Tn7_TnsA-like_N` (TnsA) 509 proteins on 2,017
phage genes and 275 on 284 prokaryote genes: the reason TnsA is not an anchor.

**The searches.** TniQ at 27.0 bits: 587 proteins on 635 genes, all prokaryote (57 genes in P1, 578 in P2). Cas7f at
25.0 bits: 114 proteins on 157 genes (6 in P1, 105 in P2, 46 phage genes in S1). Every protein the Pfam scan gave
either family is a hit; 13 TniQ hits have no TniQ row in the scan (their sequence score passes GA but no single
domain passes the domain GA); the scores of the shared proteins equal the scan's (largest difference 0.0 bit), and
`hmmsearch` run by hand on the hits and the Pfam-annotated proteins finds the same proteins at GA with the same
scores (0 only in either, largest difference 0.00 bit; `oracle.py`). The controls score 61.3 bits (TniQ) and 138.5
bits (Cas7f).

**The association.** 27 TniQ loci have a Cas domain within 10 genes (`results/tniq_cas_loci.tsv`, from
`sql/q03_tniq_cas_loci.sql`, one row per locus with the roles of its 21 genes along the TniQ's strand):

| Class | TniQ loci | What they hold | Hosts (source record names) |
|---|---|---|---|
| I-F3: Cas7f, no Cas1, Cas2 or Cas3 | 20, on 18 contigs of 16 genomes | TniQ, then a Cas5f-domain protein, Cas7f, Cas6f; TnsA, TnsB and TnsC upstream | *Vibrio crassostreae*, *V. pomeroyi*, *V. toranzoniae*, *V. apostichopi*, *V. hyugaensis*, *Vibrio* sp. ED004 (both chromosomes), *Vibrio* sp. 16; *Photobacterium damselae* subsp. *damselae* (both chromosomes); *Pseudoalteromonas translucida*, *Pseudoalteromonas* sp. JSTW; *Ferrimonas lipolytica*; *Planctobacterium marinum*; *Halomonas* sp. 'Soap Lake #7'; *Lelliottia steviae*; *Endozoicomonas euniceicola* (2 loci) and *Endozoicomonas* sp. GU-1 (2 loci): 8 genera, all Gammaproteobacteria |
| I-C-like: Cas8c (Csd1), Cas5, Cas7, no Cas7f, no Cas1, Cas2 or Cas3 | 3, two of them one element (a TniQ on either side) | a Cas5d-Cas8c-Cas7 module next to TniQ (one TniQ also has a TnsD domain), TnsB in all three, TnsA in two | *Endozoicomonas euniceicola* EF212 |
| I-B (Cas8b), no Cas1, Cas2, Cas3 | 0 | | |
| other | 4 | Cas1, Cas3 or Cas10: other systems near a transposon | *Allocoleopsis franciscana* PCC 7113 plasmid (2 loci, a cyanobacterium: one with Cas10, one with Cas1, Cas3, Cas4, Cas6), *Ideonella dechloratans* and *Paracoccus qomolangmaensis* (type I-E: Cse2, Cas3, Cas5; Cse1 in *Paracoccus*) |

The two Endozoicomonas pairs are one element each, with a TniQ on either side of the Cascade genes, so the 20 I-F3
TniQ loci are 18 elements, one per contig. By the paper's base-pair windows, built in SQL from coordinates and the
Pfam scan alone (`sql/checks_bp.sql`), 21 TniQ genes have a Cas7f within 20 kb and 19 Cas7f genes a TniQ within
10 kb (the 21st TniQ is 19 kb from its Cas7f, beyond 10 genes).

- **No Cas1, Cas2 or Cas3** in any of the 20 I-F3 windows, nor in the 19 Cas7f windows that hold a TniQ, by gene
  windows (`nbr`) and by base-pair windows (SQL): 0 of each, as the paper says.
- **The order.** In 16 of the 20 TniQ loci (0.8) the next gene along TniQ's strand carries the Cas5f domain. From
  the Cas7f side (`results/cas7f_10kb.tsv`, `nbr summarize` at 10 kb): Cas6f (`Cas_Csy4`) at +1 in 153 of 156 loci
  that have it, the Cas5f-domain gene at -1 in 151 of 154, and TniQ at -2 in 17 of its 19 loci. So: *tniQ*,
  *cas8f-cas5f*, *cas7f*, *cas6f*, all on one strand (same-strand share 1.0), the paper's order. TnsA, TnsB and TnsC
  lie upstream of *tniQ*, often with cargo between: all three within 10 genes in only 6 of the 20 loci; within 20
  kb, TnsB in 18 of the 21 TniQ genes with a Cas7f and all three in 12 (`sql/checks_bp.sql`).
- **The fusion.** No I-F3 window holds a Cas8f (`Cas_Csy1`) domain. The gene next to TniQ carries only Cas5f's
  domain but is about twice as long: 17 Cas5f-domain genes within 3 kb of a TniQ gene have proteins of 390 to 718
  amino acids (median 695), the 134 others 248 to 399 (median 300); none of the 151 carries a Cas8f domain as well
  (`sql/q04_cas5f_lengths.sql`, `results/q04_cas5f_lengths.tsv`). That is the paper's Cas8f-Cas5f fusion, whose
  Cas8f part Pfam's `Cas_Csy1` model does not see.
- **The other direction.** Of the 111 prokaryote Cas7f loci, 22 have neither Cas1 nor Cas3 within 10 genes, and 19
  of those 22 have a TniQ. The other 89 hold Cas1 or Cas3, with a separate Cas8f (`Cas_Csy1` at -2 in 77 of 79):
  ordinary type I-F systems. All 46 phage Cas7f loci (*Vibrio* phage ICP1 isolates and the JSF phages) hold Cas1 and
  Cas3 and no TniQ. *Photobacterium damselae* 04Ya311 carries both kinds: a transposon-borne I-F3 element on each
  chromosome and a complete I-F system on chromosome 2.
- **Attachment sites.** One I-F3 element has an IMPDH gene (PF00478) as the gene next to its TnsA-TnsB end (*P.
  damselae*, chromosome 2, the control): the paper's IMPDH site. None of the 20 has a GlmS (SIS), YciA-family (4HBT)
  or YifB (Mg_chelatase) domain within 10 genes; across all 635 TniQ loci, SIS is within 20 kb of 85 and
  Mg_chelatase of 34 (`nbr summarize` at 20 kb): these are the Tn7 and Tn6022-like elements of the release. The
  SRP-RNA site cannot be seen (no ncRNA).
- **Cargo.** Between the transposition genes and *tniQ* several *Vibrio* elements carry restriction-modification
  genes (`Methylase_S` in 29% of the 28 loci of the figure below, `HsdM_N+N6_Mtase` in 18%), as in figure 1.

**The paper's organisms in `rN.3`** (`sql/q02_paper_organisms.sql`, header names; `results/q02_paper_organisms.tsv`).
None of the species the Results name has a prokaryote genome in the release: *V. parahaemolyticus*, *V. fluvialis*,
*Photobacterium kishitanii*, *P. leiognathi*, *Pseudoalteromonas rubra*, *Shewanella baltica*, *S. putrefaciens*,
*Thiomicrospira* (*Hydrogenovibrio*) *crunogena*, *Cyanothece*, nor *V. cholerae* (header matches for some of these
are phages). The genera are: *Vibrio* 41 genomes (6 in P1, 35 in P2), *Pseudoalteromonas* 10, *Shewanella* 11,
*Aliivibrio* 4, *Photobacterium* 1. So the I-F3 element is in 7 of 41 *Vibrio* genomes, the one *Photobacterium*, 2
of 10 *Pseudoalteromonas*, and no *Shewanella* or *Aliivibrio* (*derived* from the two tables; genus by header).

**Figures** (`nbr flanks --figure`, rows in UPGMA order of the flank distance; colours are Pfam architectures,
grey a family in under 10% of the group's loci, white a gene with no domain):

- `figures/tniq_cas_pfam.svg` (and `.pdf`): the 27 TniQ loci with a Cas domain plus one locus brought in by a shared
  anchor protein (28 loci; `nbr flanks --anchors-file results/uids_tniq_cas.txt --level pfam --window 10g
  --group-by none`). The first 16 rows of the tree are the I-F3 elements of *Vibrio*, *Photobacterium*,
  *Pseudoalteromonas*, *Ferrimonas*, *Planctobacterium*, *Halomonas* and *Lelliottia*, each with the Cas5f, Cas7f,
  Cas6f run (green, orange, light blue) right of TniQ and, in most, TnsA-TnsB-TniB on the left; then the four
  *Endozoicomonas* I-F3 loci (a TniQ on either side of the Cascade genes), the *Endozoicomonas* I-C-like loci, and
  the complete systems near transposons.
- `figures/cas7f_all_pfam.svg` (and `.pdf`): every Cas7f locus (157) in one group at the Pfam level, 6 genes a side
  (the extra command above): the complete I-F systems with Cas1 and the Cas2-Cas3 fusion upstream and Cas8f
  (`Cas_Csy1`) at -2, the I-F3 elements with TniQ at -2 instead, and the ICP1 phages. 60 of 157 leaves are drawn
  (the largest leaf of each context type, then the largest leaves; 152 types). Read top to bottom: complete I-F
  systems (Cas1 in green, the Cas2-Cas3 fusion in pink, olive or teal by its domain architecture, Cas8f in yellow,
  the Cas5f-domain gene in light blue, Cas7f black, Cas6f orange); the ICP1 phages (Cas1, a Cas3 with only the
  `Cas3-like_C_2` domain, no Cas8f domain); and at the bottom the I-F3 elements, TniQ (brown) at -2, a longer
  Cas5f-domain gene at -1, no Cas1, Cas3 or Cas8f. Took 24.5 s and 33.4 CPU-seconds, at most 1.39 GB
  (*measured*, `/usr/bin/time -v`).
- `figures/cas7f_all_nr50.svg` (and `.pdf`): the same loci at the NR50 level, as `run.sh` draws them. It shows a
  limit rather than the biology: across 157 loci of many genera each NR50 cluster of a Cas protein is under 10% of
  the loci, so every prokaryote gene is grey, and the only coloured families are those of the 42 near-identical
  ICP1 phage loci that form one leaf (section 5).

**Context types.** `nbr flanks` cuts each tree at a distance of 0.5. At the Pfam level with 10 genes a side, each of
the 28 loci of figure 1 is its own context type (28 types), though the tree puts the I-F3 loci together. At the
NR50 level the 635 TniQ loci fall into 543 groups (one per anchor NR90 cluster), and the 157 Cas7f loci in one group
into 116 leaves and 108 types; in the first run, at the Pfam level grouped by the anchor's NR50 cluster, 442 groups
and one type per locus (a 16-locus group, 16 types). Variable cargo and genes without an id at the level (each a
token of its own) push nearly every pair above the cut. The tree is the useful output here, not the types (section
5).

## 4. Against the paper

The paper's numbers are quoted with their section; ours are *measured* on `rN.3` (section 3) unless labelled.

| The paper (section) | `rN.3` | Reading |
|---|---|---|
| Minimal I-F systems lack *cas1*, *cas2*, *cas3* ("Evolution of the Association ...") | 0 of 20 TniQ-Cas7f loci and 0 of 19 Cas7f-TniQ loci have a Cas1, Cas2 or Cas3 domain, by gene and by base-pair windows | Agrees |
| "In the vast majority of the loci corresponding to the cas1-less branch, a tnsD(tniQ) gene is located next to the cas genes" (I-F section) | 19 of the 22 prokaryote Cas7f loci without Cas1 and Cas3 have TniQ within 10 genes; TniQ at -2 of Cas7f in 17 of 19 | Agrees. We define the cas1-less set by the window, not by a tree branch |
| Systems "consist of fused *cas8f* and *cas5f*, *cas7f*, and *cas6f* genes" (Abstract) | tniQ, cas5f-domain gene, cas7f, cas6f on one strand; the cas5f-domain gene next to TniQ has a median 695 aa against 300 aa elsewhere, and no Cas8f domain | Agrees in order and in kind; Pfam names the fusion by its Cas5f half only |
| Hosts: "423 sequences from 19 genera of Gammaproteobacteria", of 2,905 Cas7f (I-F section) | 20 TniQ loci (19 Cas7f loci) in 16 genomes of 8 genera, all Gammaproteobacteria, of 157 Cas7f genes | Agrees in kind. A database of 4,961 complete and 43,599 draft genomes against 2,497 complete ones, one per species: the counts are not comparable |
| "14 complete bacterial genomes containing Tn7-like elements with the I-F CRISPR-Cas" ("Analysis of CRISPR Arrays ...") | 16 complete genomes | Similar by chance more than by design: 14 of 4,961 complete genomes is 0.28%, 16 of 2,497 is 0.64% (*derived*); we keep one genome per species, they kept every genome |
| Named hosts: *V. parahaemolyticus*, *V. fluvialis*, *Photobacterium kishitanii*, *P. leiognathi*, *Pseudoalteromonas rubra*, *Shewanella baltica* (Results) | None of these species is in the release; the element is in other species of *Vibrio*, *Photobacterium*, *Pseudoalteromonas* | Cannot be tested on `rN.3`; see section 6 for the 50 Gbp release |
| Two independent I-B captures, in cyanobacteria, without adaptation genes and Cas3 ("The cas1-Less Type I-F ...") | No TniQ locus with Cas8b; the two cyanobacterial TniQ-Cas loci (*Allocoleopsis franciscana* PCC 7113, one plasmid) have Cas10, or Cas1 and Cas3, not a minimal system. `rN.3` has 6 proteins with a Cas8b_N and 7 with a Cas8b_C domain | Not reproduced; absence is expected with so few cyanobacteria and Cas8b genes, so it is not evidence against the paper |
| (not in the paper) | 3 TniQ loci in *Endozoicomonas euniceicola* with a Cas5d-Cas8c-Cas7 (type I-C-like) module and no Cas1, Cas2, Cas3 | An observation beyond the paper, from Pfam domains only; it needs a tree and the arrays before anyone calls it a CAST |
| Attachment sites: downstream of an IMPDH gene, upstream of *yciA*, in the SRP-RNA gene; Tn7's *glmS* ("Chromosomal Insertions ...") | IMPDH beside one element's TnsA-TnsB end (*P. damselae* chromosome 2); no 4HBT (YciA family), SIS (GlmS) or YifB within 10 genes of an I-F3 TniQ; SIS within 20 kb of 85 of all 635 TniQ loci | Partly: the IMPDH site is seen once. *ffs* cannot be seen (no ncRNA); without element ends "beside" is a gene window, not an insertion site |
| Short CRISPR arrays downstream of *cas6f*; "more than 800 spacers" (Analysis of CRISPR Arrays) | not available | No CRISPR arrays in any release |
| 10,349 loci with tniQ/tnsD or tnsA (I-F section) | 635 TniQ loci | Not comparable (database size, and TnsA not used) |

Where we differ and why: (1) one complete genome per GTDB species, no drafts, so most of the paper's 423 sequences
come from genomes we do not have, and many named species are absent; (2) Pfam gathering cutoffs instead of
PSI-BLAST at 1e-4 with CDD and 217 custom Cas profiles: the Cas8f half of the fusion is invisible to Pfam, and a
weak TniQ (13 hits) passes the sequence cutoff with no domain row; (3) windows are counted in genes for the checks
and the figures, in bases for the summaries; (4) no CRISPR arrays, ncRNA, element ends or trees, so the paper's
evidence that the I-F3 systems are one monophyletic branch, and everything about spacers, cannot be checked here.

## 5. What the example shows bio-box needs

Each is a candidate card; the card that covers it is named, or "none".

1. **CRISPR arrays.** No producer writes `feature` rows of kind `crispr_array` (`ids.KIND` has the kind; `PLAN.md`
   step 4 names MinCED and CRISPRDirection with spacers into a spacer table). Without it the example cannot show an
   array downstream of *cas6f*, count spacers, or look for self-targeting spacers. Card: none.
2. **ncRNA in a release.** The SRP-RNA (*ffs*) attachment site, and any tRNA site, need ncRNA features in the
   release the search runs on. Cards E-07 (ncRNA scan in shards) and E-09 (`prok-2`); `rN.3` has none.
3. **Taxonomy.** Hosts are named here from FASTA headers; `contig_dim.tax_order` is null. A genus or GTDB species
   per genome would make "which genera" and "is the paper's species present" a join instead of a text match.
   Card: none (`STAGES.md` lists "Taxonomy assignment wired into `tax_order_id`" for a later stage; card B-12
   built `contig_dim` with the taxonomy columns null "until a taxonomy step exists").
4. **Neighbourhood statistics.** Every number here is a count; "the I-F3 elements lack Cas1" is a structural
   statement that needs none, but "TniQ is enriched next to Cas7f" would. Cards F-01 to F-03.
5. **Context types, and figure colours, across genera.** `nbr flanks` made one context type per locus for these
   elements (28 of 28 at Pfam; 108 types for 116 leaves at NR50), because cargo varies and each gene without an id
   is a token of its own. The tree is right; the cut at 0.5 is not useful for mobile elements with cargo. At the
   NR50 level the figure greys out every prokaryote Cas gene, since each species' Cas proteins form their own NR50
   cluster, each under the 10% colour threshold, while 42 near-identical phage loci (one leaf) take all the colours.
   The question whether the Pfam level should compare annotated positions only is already open for the PI
   (`planning/measurements/2026-11-12_flanks_rN3.md`, card W-25); this example adds two cases: count leaves rather
   than loci for the colour threshold, and colour NR50 clusters by their Pfam architecture. Card: a W-25/W-27
   follow-up (none yet).
6. **`bio-box sql` counts one row too few when the first column's name starts with `#`.** The guide's idiom for a
   uid file (`AGENT_GUIDE.md` section 5, `AS "#protein_uid"`) makes the column line look like a comment to the
   counter, which then takes the first row for the column line: `sql: statement 1: 26 rows` and `wrote
   .../uids.txt: 25 rows` for the same file of 26 uids (`_Tally` in `src/biobox/cli.py`). The file is right; the
   count is wrong. Card: none (a W-23 follow-up).
7. **`bio-box sql` names multi-file tables `part-0.parquet` on stderr.** For `l1_domains`, `protein_occurrence` and
   `contig_occurrence` the row-group report reads `part-0.parquet 137 of 137 row groups` instead of the table's
   name (single-file tables such as `contig_dim` are named): `bio-box sql --root "$ROOT" --release rN.3 -e "SELECT
   count(*) FROM l1_domains"` prints `read part-0.parquet 137 of 137 row groups`. Card: none (a W-23 follow-up).
8. **No contig-name table in `rN.3`.** Every `find --accession` and `--positive` reads all contig names (the
   `note: ... slow path` line; 5.8 CPU-seconds and 0.54 GB for each `find` here). Known: `GUIDE.md` section 8, card
   B-27; the 50 Gbp release has the table.
9. **A launcher that waits for room.** `run.sh` has to sum `ps` itself and wait before every launch, because the
   account was at its day cap for most of the run (section 6). `scripts/pod_run.sh` could do that once for every
   caller. Card: C-28 (the night allowance) is related, not the same; none for the day.
10. **Protein trees.** The paper's evidence for one captured I-F variant is a tree with support; bio-box has NR90
    and NR50 clusters and flank trees, not protein phylogenies. Card: none.

## 6. Running this as a release check

Before a milestone release is signed, an agent reruns the example on it and reads the check's last line. The
command is section 2's with `RELEASE` and `OUT` changed; the release needs a search pack in the profile's pack
directory (`bio-box releases` names it).

```bash
export RELEASE=rN.3
export OUT=/tmp/if2453/<you>/crispr-$RELEASE
bash "$EX/run.sh"
tail -n 1 "$OUT/check.txt"
```

`run.sh` exits with `check.py`'s code: 0 when every invariant holds, 1 when one fails or is missing, 2 on a
malformed input. `$OUT/check.txt` has one line per check: `ok`, `FAIL` (an invariant), `differs` (a recorded value
of this release that changed), `new` (no recorded value for this release yet), `missing`, `unexpected`.

**Two kinds of check** (`expected.tsv`: check id, what it asserts, kind, expectation):

- `invariant`: must hold on every release that holds the control genome (*Photobacterium damselae* subsp.
  *damselae* 04Ya311, stage P2) and the phage ICP1_2011_A (S1). Written as `== V`, `>= N`, `<= N` or `in A,B`.
- `recorded`: a count of one release, kept in `recorded/<release>.tsv`. A new release adds its own file: copy
  `$OUT/recorded_new.tsv` (written by every run) to `recorded/<release>.tsv` in the commit that records the run;
  never edit another release's file. A `differs` line on a rerun of the same release means the release or the code
  changed: find out which before recording anything.

**What a failed invariant would mean.**

| Check | If it fails |
|---|---|
| `tniq_control_bits`, `cas7f_control_bits` | The model no longer scores the control above GA: the Pfam version, the pack or the scoring changed. (A lost control already stops the search, `positive control not recovered`.) |
| `tniq_control_order` | The control element's gene order, strands or annotations differ: gene calls, the window's orientation (`offset_rel`), the Pfam scan of those proteins, or the role list in the SQL. Compare the locus table in `$OUT/tniq/loci/` with the string in `expected.tsv` |
| `pdam_complete_if_window_has_cas1_and_cas3`, `icp1_cas7f_window_has_cas1_and_cas3` | The windows no longer see Cas1 or Cas3 where they are. Then the two "no Cas1, Cas2 or Cas3" checks pass for the wrong reason: treat them as failed too |
| `tniq_cas7f_loci_with_cas1_cas2_cas3`, `cas7f_tniq_loci_with_cas1_cas2_cas3` | A transposon-associated I-F locus now has adaptation or interference genes within 10 genes. A window that crosses into the next system, a false domain call, or biology the paper did not see in its genomes (new genomes). Read the locus before calling it a defect; if the `bp_` checks below agree, it is the data, not `nbr` |
| `bp_tniq_cas7f_loci_with_cas1_cas2_cas3`, `bp_cas7f_tniq_loci_with_cas1_cas2_cas3` | The same by base-pair windows from coordinates, without the runs or `nbr`. If only the gene-window check fails, suspect `nbr`; if only this one, suspect coordinates or `protein_occurrence` |
| `tniq_cas7f_share_cas5f_at_plus1` | Under half of the I-F3 TniQ genes have the Cas8f-Cas5f gene next to them: an ordering or strand error, or a release whose I-F3 elements are mostly of another arrangement |
| `tniq_pfam_proteins_not_hit`, `cas7f_pfam_proteins_not_hit`, `max_bits_diff_search_vs_pfam` | The search and the release's own Pfam scan disagree: a pack that misses proteins, or a scoring difference |
| `oracle_*` | The search and `hmmsearch` run by hand disagree on the hits or their scores |

**Cost on `rN.3`** (*measured*, wilkcomp02, 2026-11-13, the run of section 3, 2 threads, the account at its day
cap for most of it):

| Step | Wall | CPU | Memory |
|---|---|---|---|
| `hmmfetch` of the two anchors from `Pfam-A.hmm.gz` | 52.9 s | 82.4 CPU-s | 46 MB |
| two `find` | | 6.2 and 6.1 CPU-s | 0.55 GB (unit peak) |
| `search hmm` TniQ | 275.0 s (run record) | 469.5 CPU-s for the unit (7.8 CPU-min); worker 5.8 CPU-min | worker 3.01 GB; unit 8.13 GB with page cache |
| `search hmm` Cas7f | 297.0 s | 533.6 CPU-s (8.9 CPU-min); worker 7.1 CPU-min | worker 2.62 GB; unit 7.73 GB |
| two `nbr summarize` | | 43.5 and 34.7 CPU-s | 3.8 and 3.5 GB (unit peak) |
| three `nbr flanks` (two with a figure) | | 50.1, 46.6 and 29.9 CPU-s | 2.6 to 4.3 GB (unit peak) |
| five `sql`, one `seq` | | 6.8 to 15.6 CPU-s each, 60.2 in all | 0.44 to 1.11 GB (unit peak) |
| `hmmsearch` oracle, `oracle.py`, `check.py` | 2.9 s | 2.9 CPU-s | 46 MB |

In all 1,366 CPU-seconds, 22.8 CPU-minutes (*derived*: the sum of the units' `CPUUsageNSec` from the `pod_run:
... usage` lines and of `/usr/bin/time -v` user plus system for the rest); 29 min 46 s from start to end, of which
about 13 minutes waiting for room (*derived*: 13 `waiting for room` lines, one 60 s wait each). With room it is
about a quarter of an hour (*estimate*). The first example's search took 340.7 s and 10.9 CPU-minutes at 4 threads
(`planning/measurements/2026-11-12_flanks_rN3.md`). The largest unit, 8.13 GB, is mostly the pack's page cache;
`run.sh` gives each unit `MemoryMax=16G`.

**On the 50 Gbp release.** `rM.2` (11,966 genomes, being built; `planning/tasks/STATUS.md`) is the better test: more
genomes per genus, and perhaps the paper's own species. When it has a pack, run the same command with
`RELEASE=rM.2`, then see which of the paper's species it holds:

```bash
bio-box sql --root "$ROOT" --release rM.2 "$EX/sql/q02_paper_organisms.sql" --out "$OUT/q02_paper_organisms.tsv"
```

The controls are in stage P2, which `rM.*` includes. The paper's accessions for its 12 representative elements are
in its Table S1 and Dataset S1, which were not read for this example (the PNAS supplement was not fetched); the
strains named in the Results are *Vibrio parahaemolyticus* RIMD2210633, *V. fluvialis* ATCC 33809, *Pseudoalteromonas
rubra* SCSIO6842, *Shewanella baltica* OS678 and *Thiomicrospira crunogena* XCL_2. For each one `rM.2` holds, the
next agent should find its element with `sql/q03_tniq_cas_loci.sql` and add its TniQ as a second control (a new row
in `expected.tsv` with the gene's `feature_uid`, written into `sql/checks.sql` as the first control is).

## 7. Files here

| File | What it is |
|---|---|
| `README.md` | this example |
| `run.sh` | the whole example as one command (section 6) |
| `expected.tsv` | the checks: id, assertion, kind, expectation |
| `recorded/rN.3.tsv` | the recorded values of `rN.3` |
| `check.py`, `oracle.py` | the comparison, and the `hmmsearch` oracle (standard library only); tested by `tests/r1/test_example_crispr_check.py` |
| `sql/checks.sql`, `sql/checks_bp.sql`, `sql/oracle_uids.sql`, `sql/uids_tniq_cas.sql` | what `run.sh` runs |
| `sql/q00_models.sql` to `sql/q04_cas5f_lengths.sql` | the tables of section 3 that the check does not need (`q01` is for reading: every TniQ and Cas7f gene with the roles in its base-pair window) |
| `results/` | this run's small tables: `check.txt`, `observed.tsv`, `tniq_cas_loci.tsv`, `cas7f_10kb.tsv` (`nbr summarize`), `q00_models.tsv`, `q02_paper_organisms.tsv`, `q04_cas5f_lengths.tsv`, `uids_tniq_cas.txt`, `SHA256SUMS` of what stays on the pod |
| `figures/` | the three figures of section 3, SVG and PDF |

Kept on the pod, not here (`$OUT` = `/tmp/if2453/ex2/rN.3-b` on wilkcomp02, local disk; paths and sha256 in `results/SHA256SUMS`): the runs (`$OUT/runs`), the hit tables and
locus tables (`$OUT/tniq`, `$OUT/cas7f`), the 20 kb summary of the TniQ run (2,379,860 bytes), and the `nbr flanks`
outputs.

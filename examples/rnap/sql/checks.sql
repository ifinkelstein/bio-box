-- checks.sql: the RNA polymerase example's observed values from its run, one row per check (check_id, value), for
-- the shared check.py. Run with the example's run attached (`bio-box sql --runs OUT/runs`): Pfam RNA_pol (PF00940)
-- searched at 30 bits with --neighbors 5. Windows are the run's own neighbourhoods (`neighborhoods`, 5 genes on
-- each side of every hit gene); names are the release's Pfam domains (l1_domains). Keys are feature_uid and
-- protein_uid, never an ordinal (AGENTS 3.3).
-- The control: T7 gene 1, the RNA polymerase (GenBank NC_001604 3171..5822; 0-based half-open 3170-5822, +1),
--   feature 39f9ef263e234b012d2671634468546e, protein 2eef240afefddcd2c8ac7660eb50e384, contig
--   719747a8fef97b2380955d92845861db.
WITH runs AS (  -- the example's run: the one whose seed is the PF00940 model
  SELECT run, controls_json FROM search_run WHERE seed_json LIKE '%PF00940%'
),
hits AS (
  SELECT h.target_uid AS protein_uid, max(h.seq_bitscore) AS bits
  FROM search_hit h JOIN runs r ON r.run = h.run WHERE h.target_kind = 'protein' GROUP BY ALL
),
nb AS (  -- every window slot of the run; the anchor itself is offset 0
  SELECT n.* FROM neighborhoods n JOIN runs r ON r.run = n.run
),
pfam AS (  -- the release's own Pfam rows of RNA_pol (gathering cutoff, 22 bits in Pfam 38.2)
  SELECT protein_uid, max(seq_bitscore) AS bits FROM l1_domains WHERE split_part(model_acc, '.', 1) = 'PF00940'
  GROUP BY ALL
),
dom AS (  -- the Pfam architecture of every protein in the windows: names in order along it, each name once
  SELECT protein_uid, string_agg(model_name, '+' ORDER BY first_start, model_name) AS arch
  FROM (SELECT protein_uid, model_name, min(env_start) AS first_start FROM l1_domains
        WHERE protein_uid IN (SELECT protein_uid FROM nb) GROUP BY ALL)
  GROUP BY protein_uid
),
anchors AS (  -- one row per locus (a hit gene)
  SELECT DISTINCT anchor_feature_uid, anchor_protein_uid, contig_uid FROM nb
),
loc AS (  -- with the genome kind (virus, prokaryote, plasmid) of its contig
  SELECT a.*, coalesce(g.kind, '?') AS kind
  FROM anchors a
  LEFT JOIN (SELECT contig_uid, CASE WHEN count(DISTINCT kind) = 1 THEN any_value(kind) ELSE 'mixed' END AS kind
             FROM genome_input WHERE contig_uid IN (SELECT contig_uid FROM anchors) GROUP BY contig_uid) g
    ON g.contig_uid = a.contig_uid
),
genes AS (  -- every gene of a hit protein, from protein_occurrence (not from the run's windows)
  SELECT feature_uid FROM protein_occurrence WHERE protein_uid IN (SELECT protein_uid FROM hits)
),
ctl AS (  -- the control's window, slot by slot along its strand, with each slot's architecture ('.' for none)
  SELECT nb.offset_rel, nb.feature_uid, coalesce(dom.arch, '.') AS arch
  FROM nb LEFT JOIN dom ON dom.protein_uid = nb.protein_uid
  WHERE nb.anchor_feature_uid = unhex('39f9ef263e234b012d2671634468546e')
),
coord AS (  -- the control's contig in coordinate order, from protein_occurrence alone
  SELECT feature_uid, row_number() OVER (ORDER BY start, "end", feature_uid) AS k
  FROM protein_occurrence WHERE contig_uid = unhex('719747a8fef97b2380955d92845861db')
),
coord_win AS (  -- the 5 genes on each side of the control by coordinates (it is on +1: offset = offset_rel)
  SELECT c.feature_uid, c.k - a.k AS off
  FROM coord c, (SELECT k FROM coord WHERE feature_uid = unhex('39f9ef263e234b012d2671634468546e')) a
  WHERE c.k - a.k BETWEEN -5 AND 5
)
-- invariants
SELECT 't7_control_bits' AS check_id, CAST(round(max(bits), 1) AS VARCHAR) AS value FROM hits
  WHERE protein_uid = unhex('2eef240afefddcd2c8ac7660eb50e384')
UNION ALL SELECT 't7_control_architecture', max(arch) FROM ctl WHERE offset_rel = 0
UNION ALL SELECT 't7_control_order', string_agg(CASE WHEN offset_rel = 0 THEN '^' ELSE arch END, ' ' ORDER BY offset_rel)
  FROM ctl
UNION ALL SELECT 't7_window_is_the_contig_order',
  CAST((SELECT count(*) FROM ctl) = 11 AND (SELECT count(*) FROM coord_win) = 11
       AND NOT EXISTS (SELECT offset_rel, feature_uid FROM ctl EXCEPT SELECT off, feature_uid FROM coord_win) AS VARCHAR)
UNION ALL SELECT 't7_positive_recovered', max(regexp_extract(controls_json, '"recovered": (true|false)', 1)) FROM runs
UNION ALL SELECT 'rnap_negative_control_hits', max(regexp_extract(controls_json, '"n_hits": ([0-9]+)', 1)) FROM runs
UNION ALL SELECT 'rnap_loci_virus', CAST(count(*) AS VARCHAR) FROM loc WHERE kind = 'virus'
UNION ALL SELECT 'rnap_loci_not_virus', CAST(count(*) AS VARCHAR) FROM loc WHERE kind <> 'virus'
UNION ALL SELECT 'rnap_loci_missing_windows', CAST(count(*) AS VARCHAR) FROM genes
  WHERE feature_uid NOT IN (SELECT anchor_feature_uid FROM anchors)
UNION ALL SELECT 'rnap_hits_below_cutoff', CAST(count(*) AS VARCHAR) FROM hits WHERE bits < 30
UNION ALL SELECT 'rnap_pfam_proteins_at_cutoff_not_hit', CAST(count(*) AS VARCHAR) FROM pfam
  WHERE bits >= 30 AND protein_uid NOT IN (SELECT protein_uid FROM hits)
UNION ALL SELECT 'max_bits_diff_search_vs_pfam', CAST(round(max(abs(h.bits - p.bits)), 2) AS VARCHAR)
  FROM hits h JOIN pfam p ON p.protein_uid = h.protein_uid
-- recorded (release-specific counts)
UNION ALL SELECT 'rnap_hit_proteins', CAST(count(*) AS VARCHAR) FROM hits
UNION ALL SELECT 'rnap_hit_genes', CAST(count(*) AS VARCHAR) FROM loc
UNION ALL SELECT 'rnap_hit_contigs', CAST(count(DISTINCT contig_uid) AS VARCHAR) FROM loc
UNION ALL SELECT 'rnap_hit_genomes', CAST(count(DISTINCT d.genome_uid) AS VARCHAR)
  FROM loc JOIN contig_dim_now() d ON d.contig_uid = loc.contig_uid
UNION ALL SELECT 'rnap_proteins_virus', CAST(count(DISTINCT anchor_protein_uid) AS VARCHAR) FROM loc WHERE kind = 'virus'
UNION ALL SELECT 'rnap_proteins_not_virus', CAST(count(DISTINCT anchor_protein_uid) AS VARCHAR) FROM loc
  WHERE kind <> 'virus'
UNION ALL SELECT 'rnap_genes_with_architecture', CAST(count(*) AS VARCHAR) FROM loc
  WHERE anchor_protein_uid IN (SELECT protein_uid FROM dom)
UNION ALL SELECT 'rnap_hits_not_in_pfam', CAST(count(*) AS VARCHAR) FROM hits
  WHERE protein_uid NOT IN (SELECT protein_uid FROM pfam)
UNION ALL SELECT 'rnap_pfam_proteins', CAST(count(*) AS VARCHAR) FROM pfam
UNION ALL SELECT 'rnap_pfam_proteins_below_cutoff', CAST(count(*) AS VARCHAR) FROM pfam WHERE bits < 30
ORDER BY check_id;

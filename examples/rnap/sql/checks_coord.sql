-- checks_coord.sql: the control's known neighbours by coordinates, independent of the search run and of `nbr` (no
-- run attached): the genes of T7's contig (NC_001604, contig 719747a8fef97b2380955d92845861db) in coordinate order
-- from protein_occurrence, and the release's Pfam scan of their proteins. T7 gene 0.7, the protein kinase
-- (Phage_T7_kinase, PF27024), is the gene before gene 1 (feature 39f9ef263e234b012d2671634468546e, on +1); gene 1.3,
-- the DNA ligase (DNA_ligase_A_M, PF01068), is the third after it. Rows (check_id, value) for check.py; the offset
-- is the least offset within 5 genes that carries the family, empty when none does.
WITH coord AS (
  SELECT feature_uid, protein_uid, row_number() OVER (ORDER BY start, "end", feature_uid) AS k
  FROM protein_occurrence WHERE contig_uid = unhex('719747a8fef97b2380955d92845861db')
),
win AS (
  SELECT c.feature_uid, c.protein_uid, c.k - a.k AS off
  FROM coord c, (SELECT k FROM coord WHERE feature_uid = unhex('39f9ef263e234b012d2671634468546e')) a
  WHERE c.k - a.k BETWEEN -5 AND 5
),
fam AS (
  SELECT DISTINCT w.off, split_part(d.model_acc, '.', 1) AS pf
  FROM win w JOIN l1_domains d ON d.protein_uid = w.protein_uid
)
SELECT 'coord_t7_kinase_offset' AS check_id, CAST(min(off) AS VARCHAR) AS value FROM fam WHERE pf = 'PF27024'
UNION ALL SELECT 'coord_t7_ligase_offset', CAST(min(off) AS VARCHAR) FROM fam WHERE pf = 'PF01068'
ORDER BY check_id;

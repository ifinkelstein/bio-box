-- checks_bp: the structural checks again, by base-pair windows built here from gene coordinates, independent of
-- the search runs and of `nbr` (no run attached): the loci are the genes with a Cas7f (Cas_Csy3) or TniQ domain in
-- the release's own Pfam scan, the windows the paper's (Peters et al. 2017, Methods: 10 kb up- and downstream of
-- cas7f, 20 kb of tniQ; a gene counts when any of its bases lies within the window; no wrap round the origin of a
-- circular contig). Rows (check_id, value) for check.py. Full scans of l1_domains and protein_occurrence.
WITH roles(pf, role) AS (VALUES
  ('PF06527', 'tniQ'), ('PF15978', 'tnsD'), ('PF08722', 'tnsA'), ('PF08721', 'tnsA'),
  ('PF27165', 'tnsB'), ('PF29687', 'tnsB'), ('PF29688', 'tnsB'), ('PF27135', 'tnsB'),
  ('PF11426', 'tnsC'), ('PF05621', 'tniB'), ('PF18623', 'tnsE'),
  ('PF09611', 'cas8f'), ('PF09614', 'cas5f'), ('PF09615', 'cas7f'), ('PF09618', 'cas6f'),
  ('PF24121', 'cas8b'), ('PF24122', 'cas8b'), ('PF05107', 'cas7'), ('PF09704', 'cas5'), ('PF09709', 'cas8c'),
  ('PF09559', 'cas6'), ('PF01881', 'cas6'), ('PF10040', 'cas6'),
  ('PF01867', 'cas1'), ('PF09827', 'cas2'), ('PF09707', 'cas2'), ('PF01930', 'cas4'),
  ('PF21384', 'cas3'), ('PF18019', 'cas3'), ('PF21802', 'cas3'), ('PF18395', 'cas3'), ('PF22590', 'cas3'),
  ('PF30455', 'cas3'), ('PF09481', 'cas8e'), ('PF09485', 'cas11e'),
  ('PF16593', 'cas9'), ('PF16592', 'cas9'), ('PF16595', 'cas9'), ('PF22702', 'cas9'),
  ('PF22335', 'cas10'), ('PF12469', 'cas10'), ('PF21918', 'cas12a')
),
dom AS (
  SELECT DISTINCT d.protein_uid, r.role
  FROM l1_domains d JOIN roles r ON r.pf = split_part(d.model_acc, '.', 1)
),
anchor AS (
  SELECT DISTINCT o.feature_uid, o.contig_uid, o.start, o."end", o.strand, o.stage, dom.role AS anchor_role,
         CASE dom.role WHEN 'cas7f' THEN 10000 ELSE 20000 END AS win
  FROM protein_occurrence o JOIN dom ON dom.protein_uid = o.protein_uid
  WHERE dom.role IN ('cas7f', 'tniQ')
),
gene AS (
  SELECT o.feature_uid, o.contig_uid, o.start, o."end", o.protein_uid
  FROM protein_occurrence o WHERE o.contig_uid IN (SELECT contig_uid FROM anchor)
),
near AS (
  SELECT a.feature_uid AS anchor_uid, a.anchor_role, dom.role
  FROM anchor a
  JOIN gene g ON g.contig_uid = a.contig_uid AND g.feature_uid <> a.feature_uid
             AND g."end" > a.start - a.win AND g.start < a."end" + a.win
  JOIN dom ON dom.protein_uid = g.protein_uid
),
loc AS (
  SELECT a.anchor_role, a.feature_uid, a.contig_uid, list_distinct(list(n.role)) AS roles
  FROM anchor a LEFT JOIN near n ON n.anchor_uid = a.feature_uid AND n.anchor_role = a.anchor_role
  GROUP BY ALL
),
f AS (
  SELECT *, list_contains(roles, 'cas7f') AS has_cas7f, list_contains(roles, 'tniQ') AS has_tniq,
         list_has_any(roles, ['cas1', 'cas2', 'cas3']) AS has_adapt_or_cas3
  FROM loc
)
SELECT 'bp_tniq_cas7f_loci_with_cas1_cas2_cas3' AS check_id,
       CAST(count(*) FILTER (WHERE has_adapt_or_cas3) AS VARCHAR) AS value
  FROM f WHERE anchor_role = 'tniQ' AND has_cas7f
UNION ALL SELECT 'bp_cas7f_tniq_loci_with_cas1_cas2_cas3', CAST(count(*) FILTER (WHERE has_adapt_or_cas3) AS VARCHAR)
  FROM f WHERE anchor_role = 'cas7f' AND has_tniq
UNION ALL SELECT 'bp_tniq_cas7f_loci', CAST(count(*) AS VARCHAR) FROM f WHERE anchor_role = 'tniQ' AND has_cas7f
UNION ALL SELECT 'bp_cas7f_tniq_loci', CAST(count(*) AS VARCHAR) FROM f WHERE anchor_role = 'cas7f' AND has_tniq
UNION ALL SELECT 'bp_tniq_cas7f_loci_with_tnsA_tnsB_tnsC',
  CAST(count(*) FILTER (WHERE list_contains(roles, 'tnsA') AND list_contains(roles, 'tnsB')
                          AND list_has_any(roles, ['tniB', 'tnsC'])) AS VARCHAR)
  FROM f WHERE anchor_role = 'tniQ' AND has_cas7f
UNION ALL SELECT 'bp_tniq_cas7f_loci_with_tnsB', CAST(count(*) FILTER (WHERE list_contains(roles, 'tnsB')) AS VARCHAR)
  FROM f WHERE anchor_role = 'tniQ' AND has_cas7f
UNION ALL SELECT 'bp_tniq_loci',CAST(count(*) AS VARCHAR) FROM f WHERE anchor_role = 'tniQ'
UNION ALL SELECT 'bp_cas7f_loci', CAST(count(*) AS VARCHAR) FROM f WHERE anchor_role = 'cas7f'
ORDER BY check_id;

-- q01: the example's loci from the release's own Pfam scan alone (no search run): one row per gene that carries
-- a Cas7f (Cas_Csy3) or TniQ domain, with the roles found on other genes within the paper's windows on the same
-- contig (Peters et al. 2017, Methods: 10 kb up- and downstream of cas7f, 20 kb of tniQ; here a gene counts when
-- any of its bases lies within the window). Windows do not wrap round the origin of a circular contig.
-- Full scans of l1_domains, protein_occurrence and contig_occurrence.
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
lab AS (  -- the contig's first source record, as hits.tsv names it (rN.3 has no contig_label table)
  SELECT contig_uid, first(accession ORDER BY l0_object_sha256, record_ordinal) AS accession,
         first(original_header ORDER BY l0_object_sha256, record_ordinal) AS original_header
  FROM contig_occurrence WHERE contig_uid IN (SELECT contig_uid FROM anchor) GROUP BY contig_uid
),
near AS (
  SELECT a.feature_uid AS anchor_uid, a.anchor_role, dom.role
  FROM anchor a
  JOIN gene g ON g.contig_uid = a.contig_uid AND g.feature_uid <> a.feature_uid
             AND g."end" > a.start - a.win AND g.start < a."end" + a.win
  JOIN dom ON dom.protein_uid = g.protein_uid
)
SELECT a.anchor_role, a.stage, lower(hex(a.feature_uid)) AS feature_uid, l.accession, a.start, a."end", a.strand,
       string_agg(DISTINCT n.role, ',' ORDER BY n.role) AS roles_in_window,
       l.original_header AS header
FROM anchor a
LEFT JOIN near n ON n.anchor_uid = a.feature_uid AND n.anchor_role = a.anchor_role
LEFT JOIN lab l ON l.contig_uid = a.contig_uid
GROUP BY ALL
ORDER BY a.anchor_role, a.stage, l.accession, a.start;

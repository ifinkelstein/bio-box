-- q03: every TniQ hit locus with a Cas domain within 10 genes (the TniQ run attached: `--runs OUT/runs`), one row
-- per locus: the source record, the TniQ gene, a class, and the roles of the genes from 10 upstream to 10
-- downstream along the TniQ's strand ('^' the TniQ, '.' a gene with none of the example's domains).
-- Classes: I-F3 (Cas7f, no Cas1/Cas2/Cas3); I-C-like (Cas8c/Csd1, no Cas7f, no Cas1/Cas2/Cas3); I-B-like (Cas8b, the same);
-- Cas7-only (Cas_Cas7 but neither Cas8b nor Cas8c, the same); other.
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
  ('PF22335', 'cas10'), ('PF12469', 'cas10'), ('PF21918', 'cas12a'),
  ('PF00478', 'IMPDH'), ('PF03061', '4HBT'), ('PF01380', 'SIS'), ('PF13522', 'GATase_6'), ('PF01078', 'Mg_chelatase')
),
nb AS (
  SELECT n.* FROM neighborhoods n JOIN search_run r ON r.run = n.run WHERE r.seed_json LIKE '%PF06527%'
),
dom AS (
  SELECT DISTINCT d.protein_uid, r.role
  FROM l1_domains d JOIN roles r ON r.pf = split_part(d.model_acc, '.', 1)
  WHERE d.protein_uid IN (SELECT protein_uid FROM nb)
),
slot AS (
  SELECT nb.anchor_feature_uid, nb.offset_rel,
         CASE WHEN nb.offset_rel = 0 THEN '^' ELSE coalesce(string_agg(DISTINCT dom.role, '+' ORDER BY dom.role), '.') END AS r,
         list(dom.role) FILTER (WHERE dom.role IS NOT NULL) AS rl
  FROM nb LEFT JOIN dom ON dom.protein_uid = nb.protein_uid
  GROUP BY nb.anchor_feature_uid, nb.offset_rel
),
locus AS (
  SELECT anchor_feature_uid, string_agg(r, ' ' ORDER BY offset_rel) AS window_roles,
         flatten(list(rl) FILTER (WHERE offset_rel <> 0)) AS flat
  FROM slot GROUP BY anchor_feature_uid
),
cls AS (
  SELECT *, list_has_any(flat, ['cas1', 'cas2', 'cas3']) AS adapt3 FROM locus
),
lab AS (
  SELECT contig_uid, first(accession ORDER BY l0_object_sha256, record_ordinal) AS accession,
         first(original_header ORDER BY l0_object_sha256, record_ordinal) AS header
  FROM contig_occurrence
  WHERE contig_uid IN (SELECT contig_uid FROM nb) GROUP BY contig_uid
)
SELECT CASE WHEN list_contains(c.flat, 'cas7f') AND NOT c.adapt3 THEN 'I-F3'
            WHEN c.adapt3 OR list_contains(c.flat, 'cas7f') THEN 'other'
            WHEN list_contains(c.flat, 'cas8c') THEN 'I-C-like'
            WHEN list_contains(c.flat, 'cas8b') THEN 'I-B-like'
            WHEN list_contains(c.flat, 'cas7') THEN 'Cas7-only'
            ELSE 'other' END AS class,
       l.accession, o.start, o."end", o.strand, lower(hex(c.anchor_feature_uid)) AS tniq_feature_uid,
       lower(hex(o.protein_uid)) AS tniq_protein_uid, c.window_roles, l.header
FROM cls c
JOIN (SELECT DISTINCT feature_uid, contig_uid, start, "end", strand, protein_uid FROM protein_occurrence
      WHERE feature_uid IN (SELECT anchor_feature_uid FROM nb)) o ON o.feature_uid = c.anchor_feature_uid
JOIN lab l ON l.contig_uid = o.contig_uid
WHERE list_has_any(c.flat, ['cas1', 'cas2', 'cas3', 'cas5', 'cas5f', 'cas6', 'cas6f', 'cas7', 'cas7f', 'cas8b', 'cas8c',
                            'cas8e', 'cas8f', 'cas9', 'cas10', 'cas12a'])
ORDER BY class, l.header, o.start;

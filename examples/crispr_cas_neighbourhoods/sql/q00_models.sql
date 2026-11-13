-- q00: every Pfam model of the example, with the proteins and genes the release's own Pfam scan gave it
-- (l1_domains, gathering cutoffs), split by stage. A full scan of l1_domains and protein_occurrence.
WITH roles(pf, role) AS (VALUES
  ('PF06527', 'tniQ'), ('PF15978', 'tnsD'), ('PF08722', 'tnsA'), ('PF08721', 'tnsA'),
  ('PF27165', 'tnsB'), ('PF29687', 'tnsB'), ('PF29688', 'tnsB'), ('PF27135', 'tnsB'),
  ('PF11426', 'tnsC'), ('PF05621', 'tniB'), ('PF18623', 'tnsE'),
  ('PF09611', 'cas8f'), ('PF09614', 'cas5f'), ('PF09615', 'cas7f'), ('PF09618', 'cas6f'),
  ('PF24121', 'cas8b'), ('PF24122', 'cas8b'), ('PF05107', 'cas7'), ('PF09704', 'cas5'), ('PF09709', 'cas8c'),
  ('PF09559', 'cas6'), ('PF01881', 'cas6'), ('PF10040', 'cas6'),
  ('PF01867', 'cas1'), ('PF09827', 'cas2'), ('PF09707', 'cas2'), ('PF21384', 'cas3'), ('PF01930', 'cas4'),
  ('PF18019', 'cas3'), ('PF21802', 'cas3'), ('PF18395', 'cas3'), ('PF22590', 'cas3'), ('PF30455', 'cas3'),
  ('PF09481', 'cas8e'), ('PF09485', 'cas11e'),
  ('PF16593', 'cas9'), ('PF16592', 'cas9'), ('PF16595', 'cas9'), ('PF22702', 'cas9'),
  ('PF22335', 'cas10'), ('PF12469', 'cas10'), ('PF21918', 'cas12a')
),
dom AS (
  SELECT DISTINCT d.protein_uid, r.role, d.model_name, split_part(d.model_acc, '.', 1) AS pf
  FROM l1_domains d JOIN roles r ON r.pf = split_part(d.model_acc, '.', 1)
)
SELECT dom.role, dom.pf, dom.model_name, o.stage,
       count(DISTINCT dom.protein_uid) AS proteins, count(*) AS genes, count(DISTINCT o.contig_uid) AS contigs
FROM dom JOIN protein_occurrence o ON o.protein_uid = dom.protein_uid
GROUP BY ALL ORDER BY dom.role, dom.pf, o.stage;

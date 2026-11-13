-- checks.sql: the example's observed values, one row per check (check_id, value), for check.py.
-- Run with both search runs attached (`bio-box sql --runs OUT/runs`): the TniQ run (seed PF06527) and the Cas7f
-- run (seed PF09615), each searched with --neighbors 10. Windows here are the runs' own neighbourhoods
-- (`neighborhoods`, 10 genes on each side of every hit gene); roles are the release's Pfam domains of the window's
-- proteins (l1_domains). Keys are feature_uid and protein_uid, never an ordinal (AGENTS 3.3).
-- Control genes (Photobacterium damselae subsp. damselae 04Ya311, chromosome 2, NZ_AP026781.1; 0-based half-open):
--   TniQ  NZ_AP026781.1:426340-427525 (-1)  feature f76022e1d2b79430f7a0ff9febf01a47  protein 28e32e6533f7f55e5bca64f44c89b426
--   Cas7f NZ_AP026781.1:423377-424433 (-1)  feature b5ab7b2fd228711801fbdf56af1f5fc4  protein cf65ec9336feb980ebf09b593078bf7d
--   and, in the same chromosome, a complete type I-F system: Cas7f NZ_AP026781.1:687857-688841 (+1), feature
--   099ffda89602165c6368d6a601339856.
-- Phage control (Vibrio phage ICP1_2011_A, MH310933): Cas7f MH310933:48150-49071 (-1), feature
--   23032ea6b24e8cc12738c0ef18f3a114, a phage-borne type I-F system with its adaptation genes.
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
dom AS (
  SELECT DISTINCT d.protein_uid, r.role
  FROM l1_domains d JOIN roles r ON r.pf = split_part(d.model_acc, '.', 1)
  WHERE d.protein_uid IN (SELECT protein_uid FROM neighborhoods)
),
runs AS (
  SELECT run, CASE WHEN seed_json LIKE '%PF06527%' THEN 'tniq' WHEN seed_json LIKE '%PF09615%' THEN 'cas7f' END AS anchor
  FROM search_run
),
nb AS (  -- every window slot, with its run's anchor model; the anchor itself is offset 0
  SELECT r.anchor, n.anchor_feature_uid, n.anchor_protein_uid, n.contig_uid, n.feature_uid, n.protein_uid,
         n.offset_rel, n.flags
  FROM neighborhoods n JOIN runs r ON r.run = n.run
),
win AS (  -- one row per locus: the roles of its window (the anchor gene excluded)
  SELECT nb.anchor, nb.anchor_feature_uid, any_value(nb.contig_uid) AS contig_uid,
         list_distinct(list(dom.role) FILTER (WHERE nb.offset_rel <> 0 AND dom.role IS NOT NULL)) AS roles
  FROM nb LEFT JOIN dom ON dom.protein_uid = nb.protein_uid
  GROUP BY nb.anchor, nb.anchor_feature_uid
),
loc AS (
  SELECT w.*, coalesce(g.kind, '?') AS kind,  -- the genome kind (virus, prokaryote, plasmid) of the contig
         list_contains(roles, 'cas7f') AS has_cas7f, list_contains(roles, 'tniQ') AS has_tniq,
         list_contains(roles, 'cas1') OR list_contains(roles, 'cas2') OR list_contains(roles, 'cas3') AS has_adapt_or_cas3,
         list_contains(roles, 'cas7') OR list_contains(roles, 'cas8b') OR list_contains(roles, 'cas8c') AS has_cas7_ib,
         list_contains(roles, 'cas8b') AS has_cas8b, list_contains(roles, 'cas8c') AS has_cas8c,
         list_contains(roles, 'cas8f') AS has_cas8f,
         list_contains(roles, 'tnsA') AS has_tnsa, list_contains(roles, 'tnsB') AS has_tnsb,
         list_contains(roles, 'tniB') OR list_contains(roles, 'tnsC') AS has_tnsc
  FROM win w
  LEFT JOIN (SELECT contig_uid, CASE WHEN count(DISTINCT kind) = 1 THEN any_value(kind) ELSE 'mixed' END AS kind
             FROM genome_input WHERE contig_uid IN (SELECT contig_uid FROM win) GROUP BY contig_uid) g
    ON g.contig_uid = w.contig_uid
),
cast_if AS (  -- TniQ loci with a Cas7f domain within 10 genes: the candidate transposon-associated type I-F loci
  SELECT * FROM loc WHERE anchor = 'tniq' AND has_cas7f
),
plus1 AS (  -- the gene just downstream of each such TniQ, along its strand, and its roles
  SELECT c.anchor_feature_uid, list_distinct(list(dom.role) FILTER (WHERE dom.role IS NOT NULL)) AS roles
  FROM cast_if c JOIN nb ON nb.anchor = 'tniq' AND nb.anchor_feature_uid = c.anchor_feature_uid AND nb.offset_rel = 1
  LEFT JOIN dom ON dom.protein_uid = nb.protein_uid
  GROUP BY c.anchor_feature_uid
),
ctl AS (  -- the control locus, slot by slot from -10 to +10 along the TniQ's strand
  SELECT nb.offset_rel, coalesce(string_agg(DISTINCT dom.role, '+' ORDER BY dom.role), '.') AS r
  FROM nb LEFT JOIN dom ON dom.protein_uid = nb.protein_uid
  WHERE nb.anchor = 'tniq' AND nb.anchor_feature_uid = unhex('f76022e1d2b79430f7a0ff9febf01a47')
    AND nb.offset_rel BETWEEN -10 AND 10
  GROUP BY nb.offset_rel
),
hits AS (
  SELECT r.anchor, h.target_uid AS protein_uid, max(h.seq_bitscore) AS bits
  FROM search_hit h JOIN runs r ON r.run = h.run WHERE h.target_kind = 'protein' GROUP BY ALL
),
pfam AS (  -- the release's own Pfam rows of the two anchor models (gathering cutoffs)
  SELECT CASE split_part(model_acc, '.', 1) WHEN 'PF06527' THEN 'tniq' ELSE 'cas7f' END AS anchor,
         protein_uid, max(seq_bitscore) AS bits
  FROM l1_domains WHERE split_part(model_acc, '.', 1) IN ('PF06527', 'PF09615') GROUP BY ALL
)
-- invariants
SELECT 'tniq_control_bits' AS check_id, CAST(round(max(bits), 1) AS VARCHAR) AS value FROM hits
  WHERE anchor = 'tniq' AND protein_uid = unhex('28e32e6533f7f55e5bca64f44c89b426')
UNION ALL SELECT 'cas7f_control_bits', CAST(round(max(bits), 1) AS VARCHAR) FROM hits
  WHERE anchor = 'cas7f' AND protein_uid = unhex('cf65ec9336feb980ebf09b593078bf7d')
UNION ALL SELECT 'tniq_control_order', string_agg(CASE WHEN offset_rel = 0 THEN '^' ELSE r END, ' ' ORDER BY offset_rel)
  FROM ctl
UNION ALL SELECT 'pdam_complete_if_window_has_cas1_and_cas3',
  CAST(coalesce(bool_and(list_contains(roles, 'cas1') AND list_contains(roles, 'cas3')), false) AS VARCHAR)
  FROM loc WHERE anchor = 'cas7f' AND anchor_feature_uid = unhex('099ffda89602165c6368d6a601339856')
UNION ALL SELECT 'icp1_cas7f_window_has_cas1_and_cas3',
  CAST(coalesce(bool_and(list_contains(roles, 'cas1') AND list_contains(roles, 'cas3')), false) AS VARCHAR)
  FROM loc WHERE anchor = 'cas7f' AND anchor_feature_uid = unhex('23032ea6b24e8cc12738c0ef18f3a114')
UNION ALL SELECT 'tniq_cas7f_loci_with_cas1_cas2_cas3', CAST(count(*) FILTER (WHERE has_adapt_or_cas3) AS VARCHAR)
  FROM cast_if
UNION ALL SELECT 'cas7f_tniq_loci_with_cas1_cas2_cas3', CAST(count(*) FILTER (WHERE has_adapt_or_cas3) AS VARCHAR)
  FROM loc WHERE anchor = 'cas7f' AND has_tniq
UNION ALL SELECT 'tniq_cas7f_share_cas5f_at_plus1',
  CAST(round(avg(CASE WHEN list_contains(roles, 'cas5f') THEN 1.0 ELSE 0.0 END), 3) AS VARCHAR) FROM plus1
UNION ALL SELECT 'tniq_pfam_proteins_not_hit', CAST(count(*) AS VARCHAR) FROM pfam p
  WHERE p.anchor = 'tniq' AND p.protein_uid NOT IN (SELECT protein_uid FROM hits WHERE anchor = 'tniq')
UNION ALL SELECT 'cas7f_pfam_proteins_not_hit', CAST(count(*) AS VARCHAR) FROM pfam p
  WHERE p.anchor = 'cas7f' AND p.protein_uid NOT IN (SELECT protein_uid FROM hits WHERE anchor = 'cas7f')
UNION ALL SELECT 'max_bits_diff_search_vs_pfam', CAST(round(max(abs(h.bits - p.bits)), 2) AS VARCHAR)
  FROM hits h JOIN pfam p ON p.anchor = h.anchor AND p.protein_uid = h.protein_uid
-- recorded (release-specific counts)
UNION ALL SELECT 'tniq_hit_proteins', CAST(count(*) AS VARCHAR) FROM hits WHERE anchor = 'tniq'
UNION ALL SELECT 'cas7f_hit_proteins', CAST(count(*) AS VARCHAR) FROM hits WHERE anchor = 'cas7f'
UNION ALL SELECT 'tniq_hits_not_in_pfam', CAST(count(*) AS VARCHAR) FROM hits h
  WHERE h.anchor = 'tniq' AND h.protein_uid NOT IN (SELECT protein_uid FROM pfam WHERE anchor = 'tniq')
UNION ALL SELECT 'cas7f_hits_not_in_pfam', CAST(count(*) AS VARCHAR) FROM hits h
  WHERE h.anchor = 'cas7f' AND h.protein_uid NOT IN (SELECT protein_uid FROM pfam WHERE anchor = 'cas7f')
UNION ALL SELECT 'tniq_loci_virus', CAST(count(*) AS VARCHAR) FROM loc WHERE anchor = 'tniq' AND kind = 'virus'
UNION ALL SELECT 'tniq_loci_not_virus', CAST(count(*) AS VARCHAR) FROM loc WHERE anchor = 'tniq' AND kind <> 'virus'
UNION ALL SELECT 'cas7f_loci_virus', CAST(count(*) AS VARCHAR) FROM loc WHERE anchor = 'cas7f' AND kind = 'virus'
UNION ALL SELECT 'cas7f_loci_not_virus', CAST(count(*) AS VARCHAR) FROM loc WHERE anchor = 'cas7f' AND kind <> 'virus'
UNION ALL SELECT 'tniq_cas7f_loci_virus', CAST(count(*) AS VARCHAR) FROM cast_if WHERE kind = 'virus'
UNION ALL SELECT 'tniq_cas7f_loci', CAST(count(*) AS VARCHAR) FROM cast_if
UNION ALL SELECT 'tniq_cas7f_contigs', CAST(count(DISTINCT contig_uid) AS VARCHAR) FROM cast_if
UNION ALL SELECT 'tniq_cas7f_genomes', CAST(count(DISTINCT d.genome_uid) AS VARCHAR)
  FROM cast_if c JOIN contig_dim_now() d ON d.contig_uid = c.contig_uid
UNION ALL SELECT 'tniq_cas7f_loci_with_tnsA_tnsB_tnsC', CAST(count(*) FILTER (WHERE has_tnsa AND has_tnsb AND has_tnsc) AS VARCHAR)
  FROM cast_if
UNION ALL SELECT 'tniq_cas7f_loci_with_cas8f_csy1', CAST(count(*) FILTER (WHERE has_cas8f) AS VARCHAR) FROM cast_if
UNION ALL SELECT 'tniq_minimal_non_if_loci', CAST(count(*) AS VARCHAR) FROM loc
  WHERE anchor = 'tniq' AND has_cas7_ib AND NOT has_cas7f AND NOT has_adapt_or_cas3
UNION ALL SELECT 'tniq_minimal_ib_loci', CAST(count(*) AS VARCHAR) FROM loc
  WHERE anchor = 'tniq' AND has_cas8b AND NOT has_cas7f AND NOT has_adapt_or_cas3
UNION ALL SELECT 'tniq_minimal_ic_loci', CAST(count(*) AS VARCHAR) FROM loc
  WHERE anchor = 'tniq' AND has_cas8c AND NOT has_cas7f AND NOT has_adapt_or_cas3
UNION ALL SELECT 'tniq_other_cas_loci', CAST(count(*) AS VARCHAR) FROM loc
  WHERE anchor = 'tniq' AND NOT has_cas7f AND NOT has_cas7_ib
    AND list_has_any(roles, ['cas1', 'cas2', 'cas3', 'cas5', 'cas6', 'cas8e', 'cas10', 'cas9', 'cas12a'])
UNION ALL SELECT 'cas7f_not_virus_loci_without_cas1_cas3', CAST(count(*) AS VARCHAR) FROM loc
  WHERE anchor = 'cas7f' AND kind <> 'virus' AND NOT (list_contains(roles, 'cas1') OR list_contains(roles, 'cas3'))
UNION ALL SELECT 'cas7f_not_virus_loci_without_cas1_cas3_with_tniq', CAST(count(*) AS VARCHAR) FROM loc
  WHERE anchor = 'cas7f' AND kind <> 'virus' AND NOT (list_contains(roles, 'cas1') OR list_contains(roles, 'cas3'))
    AND has_tniq
UNION ALL SELECT 'cas7f_virus_loci_with_cas1_and_cas3', CAST(count(*) AS VARCHAR) FROM loc
  WHERE anchor = 'cas7f' AND kind = 'virus' AND list_contains(roles, 'cas1') AND list_contains(roles, 'cas3')
UNION ALL SELECT 'tniq_cas7f_loci_near_' || r.role, CAST(count(DISTINCT nb.anchor_feature_uid) AS VARCHAR)
  FROM (VALUES ('IMPDH'), ('4HBT'), ('SIS'), ('Mg_chelatase')) r(role)
  LEFT JOIN dom ON dom.role = r.role
  LEFT JOIN nb ON nb.anchor = 'tniq' AND nb.protein_uid = dom.protein_uid
    AND nb.anchor_feature_uid IN (SELECT anchor_feature_uid FROM cast_if)
  GROUP BY r.role
ORDER BY check_id;

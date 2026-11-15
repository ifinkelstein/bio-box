-- checks_arrays.sql: the example's CRISPR array checks (card E-14), one row per check (check_id, value), for check.py.
-- Run as checks.sql is, with both search runs attached (`bio-box sql --runs OUT/runs`; each searched with
-- --neighbors 10). Every value is 'not applicable' when the release has no l1_crispr file (no CRISPR array step:
-- every release up to rN.3); their rows in expected.tsv say `[requires l1_crispr]`, so check.py reports them n/a.
--
-- The window in bases of a locus, as `bio-box nbr` defines it (src/biobox/search/arrays.py), built here again in SQL
-- from the runs' own windows without `nbr`: from the least start to the greatest end of the locus's rows, a row
-- reached across a circular origin (flags bit 1) moved by one contig length to the anchor's side. An array (a
-- `feature` row of kind 4, crispr_array, of a `minced-` producer) is within the window when it shares at least one
-- base with that span, on a circular contig also when moved by one contig length either way. Arrays are counted
-- once per (feature_uid, producer), of every `minced-` producer the release holds: the default one (at least 3
-- repeats) and, when the release has it, the filtered one of card E-15 (2 repeats, beside a Cas gene, never
-- overlapping a default array; its name carries `-cas1kb-`). Keys are feature_uid, never an ordinal (AGENTS 3.3).
-- Control genes (as in checks.sql): TniQ feature f76022e1d2b79430f7a0ff9febf01a47 and Cas7f feature
-- b5ab7b2fd228711801fbdf56af1f5fc4 of the Tn7-like element of Photobacterium damselae subsp. damselae 04Ya311
-- chromosome 2 (NZ_AP026781.1); the complete type I-F system's Cas7f 099ffda89602165c6368d6a601339856 on the same
-- chromosome; the Cas7f of Vibrio phage ICP1_2011_A 23032ea6b24e8cc12738c0ef18f3a114. The element's Cas6f is
-- NZ_AP026781.1:422775-423375 (README section 3).
WITH has AS (SELECT count(*) > 0 AS yes FROM catalog_files WHERE "table" = 'l1_crispr'),
roles(pf, role) AS (VALUES
  ('PF06527', 'tniQ'), ('PF09615', 'cas7f'), ('PF01867', 'cas1'), ('PF09827', 'cas2'), ('PF09707', 'cas2'),
  ('PF21384', 'cas3'), ('PF18019', 'cas3'), ('PF21802', 'cas3'), ('PF18395', 'cas3'), ('PF22590', 'cas3'),
  ('PF30455', 'cas3')
),
runs AS (
  SELECT run, CASE WHEN seed_json LIKE '%PF06527%' THEN 'tniq' WHEN seed_json LIKE '%PF09615%' THEN 'cas7f' END AS anchor
  FROM search_run
),
nb AS (
  SELECT r.anchor, n.anchor_feature_uid, n.contig_uid, n.protein_uid, n."offset", n.start, n."end", n.flags
  FROM neighborhoods n JOIN runs r ON r.run = n.run
),
dims AS (
  SELECT contig_uid, max(contig_len) AS n, bool_or(is_circular) AS circ FROM contig_dim_now()
  WHERE contig_uid IN (SELECT contig_uid FROM nb) GROUP BY contig_uid
),
dom AS (
  SELECT DISTINCT d.protein_uid, r.role
  FROM l1_domains d JOIN roles r ON r.pf = split_part(d.model_acc, '.', 1)
  WHERE d.protein_uid IN (SELECT protein_uid FROM nb)
),
win AS (  -- one row per locus: its span in the anchor's frame, its contig's length and topology, its roles
  SELECT nb.anchor, nb.anchor_feature_uid, any_value(nb.contig_uid) AS contig_uid,
         min(nb.start + CASE WHEN (nb.flags::INTEGER & 2) = 2 AND nb."offset" <> 0
                             THEN sign(nb."offset") * d.n ELSE 0 END) AS lo,
         max(nb."end" + CASE WHEN (nb.flags::INTEGER & 2) = 2 AND nb."offset" <> 0
                             THEN sign(nb."offset") * d.n ELSE 0 END) AS hi,
         any_value(d.n) AS n, coalesce(any_value(d.circ), false) AS circ,
         list_distinct(list(dom.role) FILTER (WHERE nb."offset" <> 0 AND dom.role IS NOT NULL)) AS roles
  FROM nb LEFT JOIN dims d ON d.contig_uid = nb.contig_uid LEFT JOIN dom ON dom.protein_uid = nb.protein_uid
  GROUP BY nb.anchor, nb.anchor_feature_uid
),
loc AS (
  SELECT w.*, coalesce(g.kind, '?') AS kind,
         list_contains(roles, 'cas7f') AS has_cas7f, list_contains(roles, 'tniQ') AS has_tniq,
         list_contains(roles, 'cas1') OR list_contains(roles, 'cas2') OR list_contains(roles, 'cas3') AS has_adapt_or_cas3,
         list_contains(roles, 'cas1') AND list_contains(roles, 'cas3') AS has_cas1_and_cas3
  FROM win w
  LEFT JOIN (SELECT contig_uid, CASE WHEN count(DISTINCT kind) = 1 THEN any_value(kind) ELSE 'mixed' END AS kind
             FROM genome_input WHERE contig_uid IN (SELECT contig_uid FROM win) GROUP BY contig_uid) g
    ON g.contig_uid = w.contig_uid
),
arr AS (
  SELECT contig_uid, feature_uid, producer, start, "end" FROM feature
  WHERE kind = 4 AND producer LIKE 'minced-%' AND contig_uid IN (SELECT contig_uid FROM win)
),
inwin AS (  -- (locus, array) pairs: the array shares a base with the window, moved by -n, 0 or +n on a circle
  SELECT DISTINCT l.anchor, l.anchor_feature_uid, a.feature_uid, a.producer
  FROM loc l JOIN arr a ON a.contig_uid = l.contig_uid
  JOIN (VALUES (-1), (0), (1)) k(s) ON k.s = 0 OR l.circ
  WHERE a.start + k.s * l.n < l.hi AND a."end" + k.s * l.n > l.lo
),
per AS (  -- arrays within each locus's window; n_near_cas: those of the filtered producer of card E-15 (DEC-40),
          -- the two-repeat arrays beside Cas genes, whose name carries the filter (`-cas1kb-`)
  SELECT l.*, (SELECT count(*) FROM inwin i WHERE i.anchor = l.anchor AND i.anchor_feature_uid = l.anchor_feature_uid)
         AS n_arrays,
         (SELECT count(*) FROM inwin i WHERE i.anchor = l.anchor AND i.anchor_feature_uid = l.anchor_feature_uid
            AND regexp_matches(i.producer, '-cas[0-9]+kb-')) AS n_near_cas
  FROM loc l
),
ctl AS (  -- the control chromosome, from the control TniQ's locus in the run (no scan of `feature` by feature_uid)
  SELECT any_value(contig_uid) AS contig_uid FROM nb
  WHERE anchor = 'tniq' AND anchor_feature_uid = unhex('f76022e1d2b79430f7a0ff9febf01a47')
),
na AS (SELECT CASE WHEN (SELECT yes FROM has) THEN NULL ELSE 'not applicable' END AS v)
-- invariant: the complete type I-F system of the control chromosome has an array within 10 genes of its Cas7f
SELECT 'arrays_pdam_complete_if_window_has_array' AS check_id,
  coalesce((SELECT v FROM na), CAST(coalesce((SELECT bool_or(n_arrays > 0) FROM per WHERE anchor = 'cas7f'
    AND anchor_feature_uid = unhex('099ffda89602165c6368d6a601339856')), false) AS VARCHAR)) AS value
-- recorded (release- and setting-specific)
UNION ALL SELECT 'arrays_producers', coalesce((SELECT v FROM na),  -- the settings in force (from the file list)
  (SELECT string_agg(DISTINCT regexp_extract(path, 'producer=([^/]+)', 1), ',' ORDER BY regexp_extract(path,
   'producer=([^/]+)', 1)) FROM catalog_files WHERE "table" = 'l1_crispr'))
UNION ALL SELECT 'arrays_tniq_control_window', coalesce((SELECT v FROM na), CAST((SELECT coalesce(max(n_arrays), 0)
  FROM per WHERE anchor = 'tniq' AND anchor_feature_uid = unhex('f76022e1d2b79430f7a0ff9febf01a47')) AS VARCHAR))
UNION ALL SELECT 'arrays_cas7f_control_window', coalesce((SELECT v FROM na), CAST((SELECT coalesce(max(n_arrays), 0)
  FROM per WHERE anchor = 'cas7f' AND anchor_feature_uid = unhex('b5ab7b2fd228711801fbdf56af1f5fc4')) AS VARCHAR))
UNION ALL SELECT 'arrays_control_cas6f_nearest_bp', coalesce((SELECT v FROM na), CAST((
  SELECT min(greatest(a.start, 422775) - least(a."end", 423375)) FROM feature a
  WHERE a.kind = 4 AND a.producer LIKE 'minced-%' AND a.contig_uid = (SELECT contig_uid FROM ctl)) AS VARCHAR))
UNION ALL SELECT 'arrays_pdam_complete_if_window', coalesce((SELECT v FROM na), CAST((SELECT coalesce(max(n_arrays), 0)
  FROM per WHERE anchor = 'cas7f' AND anchor_feature_uid = unhex('099ffda89602165c6368d6a601339856')) AS VARCHAR))
UNION ALL SELECT 'arrays_icp1_cas7f_window', coalesce((SELECT v FROM na), CAST((SELECT coalesce(max(n_arrays), 0)
  FROM per WHERE anchor = 'cas7f' AND anchor_feature_uid = unhex('23032ea6b24e8cc12738c0ef18f3a114')) AS VARCHAR))
UNION ALL SELECT 'arrays_tniq_cas7f_loci_with_array', coalesce((SELECT v FROM na), CAST((SELECT count(*) FROM per
  WHERE anchor = 'tniq' AND has_cas7f AND n_arrays > 0) AS VARCHAR))
UNION ALL SELECT 'arrays_tniq_cas7f_loci_with_near_cas_array', coalesce((SELECT v FROM na), CAST((SELECT count(*)
  FROM per WHERE anchor = 'tniq' AND has_cas7f AND n_near_cas > 0) AS VARCHAR))
UNION ALL SELECT 'arrays_cas7f_tniq_loci_with_array', coalesce((SELECT v FROM na), CAST((SELECT count(*) FROM per
  WHERE anchor = 'cas7f' AND has_tniq AND NOT has_adapt_or_cas3 AND n_arrays > 0) AS VARCHAR))
UNION ALL SELECT 'arrays_cas7f_complete_not_virus_loci_with_array', coalesce((SELECT v FROM na), CAST((SELECT count(*)
  FROM per WHERE anchor = 'cas7f' AND kind <> 'virus' AND has_cas1_and_cas3 AND n_arrays > 0) AS VARCHAR))
UNION ALL SELECT 'arrays_cas7f_complete_not_virus_loci', coalesce((SELECT v FROM na), CAST((SELECT count(*)
  FROM per WHERE anchor = 'cas7f' AND kind <> 'virus' AND has_cas1_and_cas3) AS VARCHAR))
UNION ALL SELECT 'arrays_cas7f_virus_loci_with_array', coalesce((SELECT v FROM na), CAST((SELECT count(*) FROM per
  WHERE anchor = 'cas7f' AND kind = 'virus' AND n_arrays > 0) AS VARCHAR))
ORDER BY check_id;

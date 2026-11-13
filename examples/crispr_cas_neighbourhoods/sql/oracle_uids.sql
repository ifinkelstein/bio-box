-- oracle_uids.sql: the proteins the hmmsearch oracle scores: every hit of the two attached runs and every protein
-- the release's Pfam scan gave a TniQ (PF06527) or Cas_Csy3 (PF09615) domain. One uid per line; the column line
-- starts with '#', so `bio-box seq --proteins-file` skips it.
SELECT DISTINCT lower(hex(protein_uid)) AS "#protein_uid" FROM (
  SELECT target_uid AS protein_uid FROM search_hit WHERE target_kind = 'protein'
  UNION
  SELECT protein_uid FROM l1_domains WHERE split_part(model_acc, '.', 1) IN ('PF06527', 'PF09615')
) ORDER BY 1;

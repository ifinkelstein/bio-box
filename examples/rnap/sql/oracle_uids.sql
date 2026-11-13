-- oracle_uids.sql: the proteins the hmmsearch oracle scores: every hit of the attached run and every protein the
-- release's Pfam scan gave an RNA_pol (PF00940) domain at its gathering cutoff (22 bits, below the search's 30). One
-- uid per line; the column line starts with '#', so `bio-box seq --proteins-file` skips it.
SELECT DISTINCT lower(hex(protein_uid)) AS "#protein_uid" FROM (
  SELECT target_uid AS protein_uid FROM search_hit WHERE target_kind = 'protein'
  UNION
  SELECT protein_uid FROM l1_domains WHERE split_part(model_acc, '.', 1) = 'PF00940'
) ORDER BY 1;

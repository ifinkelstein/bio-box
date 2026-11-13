-- q04: the Cas8f-Cas5f fusion. Every gene with a Cas_Csy2 (Cas5f, PF09614) domain, split by whether a TniQ gene
-- (PF06527) starts within 3 kb of it on the same contig and whether the same protein also carries a Cas_Csy1
-- (Cas8f, PF09611) domain; with the protein lengths (amino acids) of each class. No run needed. Full scans of
-- l1_domains and protein_occurrence; protein_seq read by protein uid.
WITH dom AS (
  SELECT protein_uid, bool_or(split_part(model_acc, '.', 1) = 'PF09614') AS csy2,
         bool_or(split_part(model_acc, '.', 1) = 'PF09611') AS csy1,
         bool_or(split_part(model_acc, '.', 1) = 'PF06527') AS tniq
  FROM l1_domains WHERE split_part(model_acc, '.', 1) IN ('PF09614', 'PF09611', 'PF06527') GROUP BY protein_uid
),
g AS (
  SELECT o.feature_uid, o.contig_uid, o.start, o."end", o.protein_uid, d.csy2, d.csy1, d.tniq
  FROM protein_occurrence o JOIN dom d ON d.protein_uid = o.protein_uid
),
c AS (
  SELECT a.feature_uid, a.protein_uid, a.csy1,
         EXISTS (SELECT 1 FROM g t WHERE t.tniq AND t.contig_uid = a.contig_uid
                 AND t.start < a."end" + 3000 AND t."end" > a.start - 3000) AS tniq_near
  FROM g a WHERE a.csy2
)
SELECT tniq_near, csy1 AS also_csy1, count(*) AS genes, count(DISTINCT c.protein_uid) AS proteins,
       min(length(protein_of(c.protein_uid))) AS min_aa,
       median(length(protein_of(c.protein_uid))) AS median_aa,
       max(length(protein_of(c.protein_uid))) AS max_aa
FROM c GROUP BY ALL ORDER BY ALL;

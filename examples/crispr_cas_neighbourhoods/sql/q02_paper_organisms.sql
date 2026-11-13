-- q02: are the organisms Peters et al. 2017 name (Results) in the release? The release has no taxonomy table
-- (contig_dim.tax_order is null), so this matches the source records' headers (contig_occurrence.original_header,
-- a full scan): a header match is a name in a FASTA header, not a GTDB species call.
WITH names(name) AS (VALUES
  ('Vibrio parahaemolyticus'), ('Vibrio fluvialis'), ('Pseudoalteromonas rubra'), ('Photobacterium kishitanii'),
  ('Photobacterium leiognathi'), ('Shewanella baltica'), ('Thiomicrospira crunogena'), ('Hydrogenovibrio crunogenus'),
  ('Cyanothece'), ('Shewanella putrefaciens'), ('Vibrio cholerae'), ('Aliivibrio'), ('Photobacterium'),
  ('Pseudoalteromonas'), ('Shewanella'), ('Vibrio'), ('Scytonema'), ('Anabaena'), ('Nostoc')
),
c AS (
  SELECT DISTINCT o.contig_uid, o.stage, o.original_header FROM contig_occurrence o
)
SELECT n.name, c.stage, count(DISTINCT c.contig_uid) AS contigs, count(DISTINCT d.genome_uid) AS genomes,
       min(c.original_header) AS example_header
FROM names n LEFT JOIN c ON c.original_header LIKE '%' || n.name || '%'
LEFT JOIN contig_dim_now() d ON d.contig_uid = c.contig_uid
GROUP BY ALL ORDER BY n.name, c.stage;

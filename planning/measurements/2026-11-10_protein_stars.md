# Stop codons in `protein_seq` (card W-11)

Question: does `protein_seq` hold any `*`, so that FASTA for a person would need to drop or flag it?

Command (*measured*, host wilkcomp02, 2026-11-10, image `bio-box-test_74b3e05.sif`, one thread, read-only bind of the store,
`systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=0 nice -n 19`):

```
python scripts/stress/protein_stars.py /data rL.2
```

with `/data` bound read-only to `/stor/scratch/Wilke/if2453/all-the-data/provisional/rL`. The file list is the
`protein_seq` entries of `manifests/release=rL.2.json`.

Result (*measured*, same run):

| quantity | count |
|---|---|
| `protein_seq` files | 2 |
| sequences | 2,053,315 |
| ending in `*` | 0 |
| with an internal `*` | 0 |
| equal to `*` | 0 |

Consequence: the FASTA writers are unchanged; `tests/r2/test_micro_s0_search.py::test_w11_no_stored_protein_contains_a_star`
pins the fixture to the same property.

# Search pack size, and why the pack is kept off shared storage (card W-03, 2026-11-10)

## What was measured

*Measured* on wilkcomp02, 2026-11-10, in the pinned test image, by a probe test run once with
`sh scripts/pod_pytest.sh W-03 tests/r2/test_w03_probe.py -x` (the probe was not committed). It
copied the C-08 fixture root (`tests/r2/test_builders_out.py`, `direct`: stage S0 = lambda and
*M. genitalium*, 628,578 bp; release `rS0.c08`; the six-model test library `test-1`) twice, and built
the release's pack in each copy with the code at `dd6103a`. Then it applied `contracts.check_file_sizes`
to the pack's files, registered for the probe as one artifact partitioned by `release`, and compared the
two packs byte for byte.

| file | bytes |
|---|---:|
| `pack_manifest.json` | 3,111 |
| `unique/part-0.faa.zst` (589 proteins, 197,346 aa) | 125,394 |
| `igr/part-0.fna.zst` (189 IGRs of at least 30 nt, 36,004 nt) | 15,785 |
| `igr/igr_meta.parquet` | 17,470 |
| `dark_orf/part-0.faa.zst` (611 dark ORFs, 30,891 aa) | 33,893 |
| `dark_orf/dark_orf_meta.parquet` | 22,608 |
| `models/test-1.hmm` (6 models) | 1,924,240 |
| total, 7 files | 2,142,501 |

- The LY-04 rule at the test profiles' 1,000-byte floor: no problem. At the pod's 64 MiB floor
  (`small_file_floor_bytes: 67108864`): `search_pack {'release': 'rS0.c08'}: 7 of 7 files under 64 MB
  (2.0 MB in all); split only above the band`.
- Determinism: the two builds are byte-identical in every file, with the same `pack_id`
  (`885750263f48d6c68567b4796504c2dc`). The `DATE` lines of `models/test-1.hmm` are copied from the
  input library (`DATE  Tue Nov 10 12:24:56 2026` in both). pyhmmer does not stamp a build time when it
  re-writes a model it has read. The library's sha256 is among the pack manifest's `sources`, so the
  `pack_id` changes only when the library file does.

## What that implies at 10 Gbp

*Derived*: the fixture's rates × 10^10 / 628,578 = 15,909. FASTA text is the residues plus 35 bytes
per record (`>`, 32 hex digits, two newlines). The zstd ratio is the fixture's: 0.575 bytes per text
byte for `unique`, 0.370 for `igr`, 0.648 for `dark_orf`. Shards are cut at 2 GiB of text
(`pack.SHARD_TEXT_BYTES`).

| target | text | files | zstd bytes |
|---|---:|---|---:|
| `unique` | 3.47 GB | 2 shards | about 1.24 GB and 0.76 GB (an upper bound: distinct proteins grow sublinearly) |
| `igr` | 0.68 GB | 1 shard | about 0.25 GB |
| `dark_orf` | 0.83 GB | 1 shard | about 0.54 GB |
| `igr_meta.parquet` | | 1 | about 0.28 GB |
| `dark_orf_meta.parquet` | | 1 | about 0.36 GB |
| `pack_manifest.json` | | 1 | a few KB at any scale |
| `models/<library>.hmm` | | 1 per library | not measured for Pfam |

That is 8 files with one library, about 3.4 GB plus the models. The fixture's rate is about 3.4 MB of
pack per Mbp (*derived*, 2,142,501 / 0.6286). A release under about 20 Mbp therefore has a whole
pack under 64 MiB.

## Why the pack lives outside the data root

Every release's pack has files under 64 MiB, so the size rule refuses it on `/stor` at any scale,
however its files are grouped. Its own `pack_manifest.json` is always a few KB. A small library's
model file is small. The last shard of a target can be any size: at 6.5 Gbp the last `unique`
shard falls to about 61 MB (*derived*: 2.254 GB of text, less 2 GiB, × 0.575). In a small release every file is small. Exempting the pack
from the rule was not wanted.

The pack is a derived, rebuildable cache of a release: the same inputs give the same bytes and the
same `pack_id`. So since card W-03 it is built into a pack directory outside the data root
(`bio-box pack build --pack-dir DIR`; by default `<profile tmp>/packs`). It is not a release file and
not a `build-stage` step, and a search whose pack is missing says how to build it.

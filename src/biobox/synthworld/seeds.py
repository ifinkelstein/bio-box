"""Seeds for synthworld (planning/tests/synthworld_spec.md §6.1).

One master seed (u64) per world. Every stage and entity draws from its own stream,
``SeedSequence(seed, spawn_key=blake3("<stage>/<entity>")[:16] as four big-endian uint32)``, so a
draw never depends on worker count, iteration order, the clock or any global random state.
External tools (ZOMBI2 ``seed=``, ``cmemit``/``hmmemit --seed``) get a derived non-zero
31-bit seed: 0 means "arbitrary" in HMMER and Infernal and is never passed.

The component paths in use (part of the output contract; renaming one changes every world):

- ``S1/species_tree``: ZOMBI2 species seed;
- ``S1/clonal_burst``: the tip a burst is grafted onto, and the burst's split times;
- ``S2/genomes``: ZOMBI2 ordered-genomes seed;
- ``S3/panel``: the synthetic HMM panel; ``S3/panel/<NAME>``: a Pfam HMM's selection, dark and promiscuous
  draws;
- ``S3/family/<family_id>``: a family's architecture, target identity, root protein and rate categories;
- ``S4/lineages``: per species branch, its GC3 step and relaxed-clock multiplier;
- ``S4/family/<family_id>``: a family's evolution down its gene tree (pairs for the rate solve, root codons,
  substitutions, synonymous events, indels);
- ``S5/family/<family_id>``: a family's cassette layout, root cassettes and their evolution down its gene
  tree (with its RNA motif, when linked);
- ``S5/motif/panel``: the synthetic RNA motif models; ``S5/motif/links``: which families carry which model;
- ``S5/emit/<link>``: a link's clade emissions (synthetic draws, or the ``cmemit --seed`` tool seed);
- ``S5/origin/<genome_id>/<replicon>``: where a circular replicon is cut;
- ``S5/tight/<up>/<down>``: not a stream but a BLAKE3 hash of ``<seed>/S5/tight/<up>/<down>`` (families in
  transcription order), the background tight-junction draw of §2.4;
- ``S6/phages/world``: the phage set's master seed (its first two uint32 words); the phage set's own streams
  are these paths under that seed;
- ``S6/placement``: which phages are integrated, and per prophage and array its ranked branches, insertion
  points along them and strand; ``S6/anchor/<item>``: the item's anchor family;
- ``S6/prophage/<prophage_id>``: a prophage's evolution down its clade;
- ``S6/crispr/<array_id>``: an array's repeat, leader, spacers and their acquisition and deletion; its Cas
  families use ``S3/family/<id>``, ``S4/family/<id>`` and ``S5/family/<id>`` with family numbers 900000 + 10 a
  + k (`mobile`).
"""

from __future__ import annotations

import numpy as np
from blake3 import blake3

_U64 = 2**64
_TOOL_MAX = 2**31 - 1


def check_seed(seed) -> int:
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or not 0 <= int(seed) < _U64:
        raise ValueError(f"a world seed is an unsigned 64-bit integer, got {seed!r}")
    return int(seed)


def spawn_key(path: str) -> tuple[int, int, int, int]:
    """The spawn key of one component: blake3(path) truncated to 16 bytes, as four uint32 words."""
    d = blake3(path.encode("utf-8")).digest(length=16)
    return tuple(int.from_bytes(d[i:i + 4], "big") for i in (0, 4, 8, 12))  # type: ignore[return-value]


def seed_sequence(seed: int, path: str) -> np.random.SeedSequence:
    return np.random.SeedSequence(check_seed(seed), spawn_key=spawn_key(path))


def generator(seed: int, path: str) -> np.random.Generator:
    """The component's own generator (PCG64)."""
    return np.random.Generator(np.random.PCG64(seed_sequence(seed, path)))


def tool_seed(seed: int, path: str) -> int:
    """A derived seed for an external tool, in [1, 2**31 - 1] (never 0)."""
    word = int(seed_sequence(seed, path).generate_state(1, dtype=np.uint32)[0])
    return word % _TOOL_MAX + 1

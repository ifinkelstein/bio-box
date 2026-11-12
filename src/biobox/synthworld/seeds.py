"""Seeds for synthworld (planning/tests/synthworld_spec.md §6.1).

One master seed (u64) per world. Every stage and entity draws from its own stream,
``SeedSequence(seed, spawn_key=blake3("<stage>/<entity>")[:16] as four big-endian uint32)``, so a
draw never depends on worker count, iteration order, the clock or any global random state.
External tools (ZOMBI2 ``seed=``, later ``cmemit``/``hmmemit --seed``) get a derived non-zero
31-bit seed: 0 means "arbitrary" in HMMER and Infernal and is never passed.

The component paths in use (part of the output contract; renaming one changes every world):

- ``S1/species_tree``: ZOMBI2 species seed;
- ``S1/clonal_burst``: the tip a burst is grafted onto, and the burst's split times;
- ``S2/genomes``: ZOMBI2 ordered-genomes seed.
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

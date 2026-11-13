"""Genetic code and codon usage for synthworld's coding sequences (planning/tests/synthworld_spec.md §2.3
"Codons"; card F-02 session 3).

Codons are indexed 0..63 as ``16 i + 4 j + k`` over the bases ``TCAG`` (``CODONS``). Amino acids are
indexed 0..19 in HMMER's order ``ACDEFGHIKLMNPQRSTVWY`` (``AA``), the order of a Pfam HMM's emission rows;
two pseudo-residues complete a coding sequence: ``START`` (20, the initiation codon, translated as M) and
``STOP`` (21).

Codon usage model (`usage_table`), one table per lineage: an amino acid's synonymous codons under genetic
code 11 are weighted ``gc3`` if the codon ends in G or C and ``1 - gc3`` otherwise, then normalized (so a
lineage's GC3 sets its usage; Met and Trp have one codon each). Initiation and stop codons have fixed
weights from the preset (``codons.start_codons``, ``codons.stop_codons``). Genetic codes 4 and 25 (S8) are
not implemented here.
"""

from __future__ import annotations

import numpy as np

BASES = "TCAG"
CODONS: list[str] = [a + b + c for a in BASES for b in BASES for c in BASES]
#: NCBI transl_table=11 amino acids in CODONS order (the standard code's assignments)
_CODE11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE11: dict[str, str] = dict(zip(CODONS, _CODE11))
#: NCBI transl_table=11 initiation codons
STARTS11 = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
AA = "ACDEFGHIKLMNPQRSTVWY"
START, STOP = 20, 21
#: the protein letter of each residue code (START reads as M; STOP has none)
LETTERS = np.frombuffer((AA + "M*").encode("ascii"), dtype=np.uint8)
#: codon index -> its three bases, for joining a coding sequence with one fancy index
CODON_BYTES = np.frombuffer("".join(CODONS).encode("ascii"), dtype=np.uint8).reshape(64, 3)
_GC3 = np.array([c[2] in "GC" for c in CODONS])


def synonyms(code: dict[str, str] = CODE11) -> list[np.ndarray]:
    """For each amino acid of ``AA``, the indices of its codons under ``code``."""
    return [np.array([i for i, c in enumerate(CODONS) if code[c] == aa]) for aa in AA]


_SYN = synonyms()


def usage_table(gc3: float, start_weights: dict[str, float], stop_weights: dict[str, float]) -> np.ndarray:
    """(22, 64) codon probabilities given a residue code: rows 0..19 the amino acids, ``START`` and
    ``STOP``. Refuses an initiation codon outside code 11's or a stop codon that is not one."""
    g = float(gc3)
    if not 0.0 < g < 1.0:
        raise ValueError(f"gc3 must lie in (0, 1), got {gc3}")
    tab = np.zeros((22, 64))
    for i, idx in enumerate(_SYN):
        w = np.where(_GC3[idx], g, 1.0 - g)
        tab[i, idx] = w / w.sum()
    for row, weights, ok in ((START, start_weights, set(STARTS11)),
                             (STOP, stop_weights, {c for c, a in CODE11.items() if a == "*"})):
        bad = set(weights) - ok
        if bad or not weights:
            raise ValueError(f"codon weights {weights} include {sorted(bad)}, not allowed in that row")
        items = sorted((str(c), float(w)) for c, w in weights.items())   # key order never matters
        total = sum(w for _, w in items)
        for c, w in items:
            tab[row, CODONS.index(c)] = w / total
    return tab


def cumulative(tab: np.ndarray) -> np.ndarray:
    """Row-wise cumulative probabilities with the last column exactly 1, for inverse-CDF draws."""
    cum = np.cumsum(tab, axis=-1)
    cum[..., -1] = 1.0
    return cum


def draw(cum: np.ndarray, rows: np.ndarray, u: np.ndarray) -> np.ndarray:
    """One draw per element: the first column whose cumulative probability in row ``rows[i]`` reaches
    ``u[i]`` (``u`` uniform in [0, 1))."""
    idx = np.count_nonzero(cum[rows] <= u[:, None], axis=1)
    return np.minimum(idx, cum.shape[-1] - 1)


#: the residue code of each codon under code 11 (stops STOP), for checking a CDS against its residues
CODE11_RESIDUE = np.array([AA.index(a) if a != "*" else STOP for a in _CODE11], dtype=np.int8)


def translate(cds: str) -> str:
    """Translate a coding sequence under code 11, its first codon as the initiator M when it is one of
    code 11's initiation codons; stops read ``*``."""
    out = [CODE11[cds[i:i + 3]] for i in range(0, len(cds) - len(cds) % 3, 3)]
    if out and cds[:3] in STARTS11:
        out[0] = "M"
    return "".join(out)

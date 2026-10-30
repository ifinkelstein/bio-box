"""Stop-to-stop ORFs in six frames, and the dark-ORF rule for the search pack (SEARCH_PLAN SD-07).

An ORF runs from the codon after one stop to its next stop, inclusive (the span includes the
closing stop codon, as gene calls do; the amino-acid sequence does not). Coordinates are
0-based half-open on the + strand. On a linear contig the runs before the first stop and
after the last one are kept and marked partial (open at the contig end). On a circular
contig ORFs are read around the origin, so one may have `end > contig_len`, and every ORF is
bounded by two stops.

A dark ORF is a stop-to-stop ORF of at least `MIN_AA` residues that

    1. is not a called gene's own ORF (same strand and same stop codon), and
    2. has less than `MAX_GENE_OVERLAP` of its span inside called genes (either strand),
       unless it touches the `END_BP` at either end of a linear contig, where all frames are kept.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass

import pyhmmer

from biobox.l0seq.store import revcomp

MIN_AA = 30
END_BP = 3_000
MAX_GENE_OVERLAP = 0.5

_DNA = pyhmmer.easel.Alphabet.dna()
_CODES: dict[int, pyhmmer.easel.GeneticCode] = {}


@dataclass(frozen=True)
class Orf:
    start: int
    end: int
    strand: int
    aa: str
    partial: bool


def _translate(nt: str, code: int) -> str:
    if code not in _CODES:
        _CODES[code] = pyhmmer.easel.GeneticCode(code)
    n = len(nt) - len(nt) % 3
    if n == 0:
        return ""
    seq = pyhmmer.easel.TextSequence(name=b"f", sequence=nt[:n]).digitize(_DNA)
    return seq.translate(_CODES[code]).textize().sequence


def _runs(aa: str) -> list[tuple[int, int, bool, bool]]:
    """(first codon, stop codon index or len(aa), open left, open right) of every run between stops."""
    stops = [i for i, c in enumerate(aa) if c == "*"]
    out, prev = [], -1
    for s in stops:
        out.append((prev + 1, s, prev == -1, False))
        prev = s
    out.append((prev + 1, len(aa), prev == -1, True))
    return out


def six_frame_orfs(seq: bytes, circular: bool, code: int = 11, min_aa: int = MIN_AA) -> list[Orf]:
    s = seq.upper().decode("ascii")
    L = len(s)
    if L < 3:
        return []
    out: list[Orf] = []
    for strand in (1, -1):
        base = s if strand == 1 else revcomp(seq.upper()).decode("ascii")
        text = base * 3 if circular else base  # room for an ORF of up to L nt starting anywhere in [3, L + 3)
        for f in range(3):
            aa = _translate(text[f:], code)
            for a, b, open_l, open_r in _runs(aa):
                if b - a < min_aa:
                    continue
                lo, hi = f + 3 * a, f + 3 * b + (0 if open_r else 3)
                if circular:
                    # one representative per circular ORF: a start in [3, L + 3), both stops seen
                    if open_l or open_r or not 3 <= lo < L + 3 or hi - lo > L:
                        continue
                    lo, hi = (lo - L, hi - L) if lo >= L else (lo, hi)
                if strand == 1:
                    start, end = lo, hi
                elif circular:
                    start, end = L - hi, L - lo
                    if start < 0:
                        start, end = start + L, end + L
                else:
                    start, end = L - hi, L - lo
                out.append(Orf(start, end, strand, aa[a:b], not circular and (open_l or open_r)))
    out.sort(key=lambda o: (o.start, o.end, o.strand))
    return out


def _merged(genes: list[tuple[int, int, int]], length: int, circular: bool) -> list[tuple[int, int]]:
    spans = []
    for s, e, _ in genes:
        spans.append((s, e))
        if circular:
            spans += [(s - length, e - length), (s + length, e + length)]
    spans.sort()
    merged: list[list[int]] = []
    for s, e in spans:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def _covered(merged: list[tuple[int, int]], starts: list[int], lo: int, hi: int) -> int:
    i = max(bisect.bisect_right(starts, lo) - 1, 0)
    total = 0
    while i < len(merged) and merged[i][0] < hi:
        s, e = merged[i]
        total += max(0, min(e, hi) - max(s, lo))
        i += 1
    return total


def dark_orfs(seq: bytes, circular: bool, genes: list[tuple[int, int, int]], code: int = 11,
              min_aa: int = MIN_AA, end_bp: int = END_BP,
              max_gene_overlap: float = MAX_GENE_OVERLAP) -> list[tuple[Orf, float, bool]]:
    """(orf, fraction of its span in called genes, touches a linear contig end) for each dark ORF."""
    L = len(seq)
    own = set()
    for s, e, strand in genes:
        own.add((strand, e % L if strand == 1 else s % L))
    merged = _merged(genes, L, circular)
    starts = [s for s, _ in merged]
    out = []
    for o in six_frame_orfs(seq, circular, code, min_aa):
        if (o.strand, (o.end if o.strand == 1 else o.start) % L) in own:
            continue
        frac = _covered(merged, starts, o.start, o.end) / (o.end - o.start)
        at_end = not circular and (o.start < end_bp or o.end > L - end_bp)
        if frac >= max_gene_overlap and not at_end:
            continue
        out.append((o, frac, at_end))
    return out

"""RNA motif models for synthworld's IGRs (planning/tests/synthworld_spec.md §2.1 S5, §2.5; card F-02 session 4).

A motif *emission* is one DNA sequence with its secondary structure: base pairs (i, j), i < j, 0-based on the
emitted sequence. Two model sources (preset ``rna_motifs.source``):

- ``rfam`` (§2.5): the preset's ``rfam_models`` accessions are copied verbatim out of an Rfam CM library
  (``Rfam.cm`` or ``Rfam.cm.gz``, read line by line, `extract_rfam`) into one CM file each, and a clade's
  emissions come from Infernal's ``cmemit -N <clades> -a --dna --seed <s> [--exp x] [-l]`` (`cmemit`): the
  Stockholm alignment it writes gives each sequence (its letters, inserts included) and the consensus
  structure ``SS_cons`` (WUSS), from which a pair is kept where both of its columns hold a residue of that
  sequence (`emissions_from_stockholm`). One call per link with ``-N`` equal to its number of clades draws
  the clades' emissions independently from the model, the same law as one ``-N 1`` call per clade (§2.5
  step 2), with one process per link. ``--seed`` is a derived non-zero seed (§6.1).
- ``synthetic`` (the ``micro`` default, so a world builds where neither Infernal nor Rfam is present, as in
  GitHub Actions): stream ``S5/motif/panel`` draws ``n_models`` models. A model is a chain of hairpins
  (``synthetic.hairpins`` of them, each a 5' linker, a stem, a loop and the closing stem, then a 3' linker;
  lengths uniform on the preset ranges); an unpaired column's base distribution is Dirichlet(``unpaired_
  concentration`` x 1/4 each); a stem pair's is Dirichlet(``pair_concentration`` x the pair weights) over the
  six canonical pairs. A clade's emission (stream ``S5/emit/<link>``) draws each unpaired column and each
  pair from its distribution, so every synthetic pair is canonical.

Evolution within a clade (§2.5 step 3, `evolve_motif`): on a branch, an unpaired position has a substitution
event with probability 1 - exp(-rate t) and takes a base from the lineage's intergenic composition; a pair
has an event with the same probability and is replaced by a pair drawn from the canonical pair weights
(``intergenic.cassettes.pair_weights``): paired positions change only by swaps among the Watson-Crick and GU
pairs, so a pair that was canonical stays canonical (a non-canonical pair cmemit emitted stays until its
first event).
"""

from __future__ import annotations

import gzip
import math
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from biobox.synthworld import seeds

BASES = "ACGT"
#: the six canonical pairs (5' base, 3' base) in DNA letters: Watson-Crick and GU (GT)
CANONICAL = ("AT", "TA", "GC", "CG", "GT", "TG")
_OPEN, _CLOSE = "<([{", ">)]}"
_CODE = {b: i for i, b in enumerate(BASES)}


@dataclass(frozen=True)
class Emission:
    seq: str                                  # uppercase DNA
    pairs: tuple[tuple[int, int], ...]        # (i, j), i < j, 0-based in seq, sorted by i


@dataclass
class MotifModel:
    index: int
    accession: str
    name: str
    source: str                               # synthetic | rfam
    length: int                               # consensus length (Rfam CLEN) or the synthetic model's length
    ga: float                                 # Rfam gathering cutoff (bits); nan for a synthetic model
    cm_file: Path | None = None               # rfam: the model's own CM file in scratch
    structure: str = ""                       # synthetic: dot-bracket consensus
    unpaired: np.ndarray | None = field(default=None, repr=False)     # (L, 4) per column (rows of paired ones unused)
    pair_cum: np.ndarray | None = field(default=None, repr=False)     # (n_pairs, 6) cumulative pair distributions
    pairs: tuple[tuple[int, int], ...] = ()


# --- structures ------------------------------------------------------------------------------------------

def parse_wuss(ss: str) -> list[tuple[int, int]]:
    """Base pairs of a WUSS (or dot-bracket) string: ``<>``, ``()``, ``[]`` and ``{}`` nest; every other
    character (``.,:_-~`` and the pseudoknot letters) is unpaired. Sorted by the opening position."""
    stacks: dict[str, list[int]] = {c: [] for c in _OPEN}
    out = []
    for i, c in enumerate(ss):
        if c in _OPEN:
            stacks[c].append(i)
        elif c in _CLOSE:
            o = _OPEN[_CLOSE.index(c)]
            if not stacks[o]:
                raise ValueError(f"unbalanced {c!r} at {i} in {ss!r}")
            out.append((stacks[o].pop(), i))
    if any(stacks.values()):
        raise ValueError(f"unbalanced structure {ss!r}")
    return sorted(out)


def dot_bracket(n: int, pairs) -> str:
    s = ["."] * n
    for i, j in pairs:
        s[i], s[j] = "(", ")"
    return "".join(s)


def parse_stockholm(text: str) -> tuple[list[tuple[str, str]], str]:
    """(name, aligned sequence) in file order, and ``#=GC SS_cons``, of one Stockholm alignment; blocks of an
    interleaved file are joined per name."""
    seqs: dict[str, list[str]] = {}
    ss: list[str] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("# STOCKHOLM") or line.startswith("//"):
            continue
        if line.startswith("#=GC"):
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "SS_cons":
                ss.append(parts[2])
            continue
        if line.startswith("#"):
            continue
        name, aln = line.split()
        seqs.setdefault(name, []).append(aln)
    return [(k, "".join(v)) for k, v in seqs.items()], "".join(ss)


def emissions_from_stockholm(text: str) -> list[Emission]:
    """Each aligned sequence's letters (gaps ``-`` and ``.`` removed, uppercased, U read as T) and the
    consensus pairs whose two columns both hold one of its residues."""
    rows, ss = parse_stockholm(text)
    cons = parse_wuss(ss)
    out = []
    for _name, aln in rows:
        if len(aln) != len(ss):
            raise ValueError(f"aligned length {len(aln)} differs from SS_cons {len(ss)}")
        col2pos, seq = {}, []
        for c, ch in enumerate(aln):
            if ch in "-.~":
                continue
            col2pos[c] = len(seq)
            seq.append(ch.upper().replace("U", "T"))
        pairs = tuple(sorted((col2pos[a], col2pos[b]) for a, b in cons if a in col2pos and b in col2pos))
        s = "".join(seq)
        if set(s) - set(BASES):
            raise ValueError(f"an emission with letters outside ACGT: {sorted(set(s) - set(BASES))}")
        out.append(Emission(s, pairs))
    return out


# --- synthetic models ------------------------------------------------------------------------------------

def _pair_cum(weights: dict[str, float]) -> np.ndarray:
    w = np.array([float(weights[p]) for p in CANONICAL])
    c = np.cumsum(w / w.sum())
    c[-1] = 1.0
    return c


def synthetic_models(cfg: dict, pair_weights: dict[str, float], seed: int) -> list[MotifModel]:
    """``rna_motifs.n_models`` synthetic models from stream ``S5/motif/panel`` (module docstring)."""
    syn = cfg["synthetic"]
    rng = seeds.generator(seed, "S5/motif/panel")
    base_w = np.array([float(pair_weights[p]) for p in CANONICAL])
    base_w = base_w / base_w.sum()
    out = []
    for k in range(int(cfg["n_models"])):
        def u(lo_hi):
            lo, hi = (int(x) for x in lo_hi)
            return int(rng.integers(lo, hi + 1))
        n_h = u(syn["hairpins"])
        ss: list[str] = []
        for _ in range(n_h):
            ss.append("." * u(syn["linker"]))
            stem, loop = u(syn["stem"]), u(syn["loop"])
            ss.append("(" * stem + "." * loop + ")" * stem)
        ss.append("." * u(syn["linker"]))
        structure = "".join(ss)
        pairs = tuple(parse_wuss(structure))
        L = len(structure)
        unp = rng.dirichlet(np.full(4, float(syn["unpaired_concentration"]) / 4.0), size=L)
        pc = np.cumsum(rng.dirichlet(float(syn["pair_concentration"]) * base_w * 6.0, size=len(pairs)), axis=1)
        pc[:, -1] = 1.0
        out.append(MotifModel(k, f"SYNRF{k + 1:04d}", f"synRNA{k + 1:02d}", "synthetic", L, float("nan"),
                              None, structure, unp, pc, pairs))
    return out


def synthetic_emit(model: MotifModel, n: int, rng: np.random.Generator) -> list[Emission]:
    out = []
    paired = {i for p in model.pairs for i in p}
    for _ in range(int(n)):
        seq = [""] * model.length
        u = rng.random(model.length)
        for i in range(model.length):
            if i not in paired:
                c = np.cumsum(model.unpaired[i])
                c[-1] = 1.0
                seq[i] = BASES[int(np.searchsorted(c, u[i], side="right"))]
        up = rng.random(len(model.pairs))
        for k, (i, j) in enumerate(model.pairs):
            p = CANONICAL[min(int(np.searchsorted(model.pair_cum[k], up[k], side="right")), 5)]
            seq[i], seq[j] = p[0], p[1]
        out.append(Emission("".join(seq), model.pairs))
    return out


# --- Rfam models and cmemit -----------------------------------------------------------------------------

def _open_text(path: Path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="ascii") if path.name.endswith(".gz") else open(path, encoding="ascii")


def extract_rfam(library: Path, accessions, out_dir: Path) -> dict[str, dict]:
    """Copy each model of ``accessions`` verbatim (its CM block and its filter HMM block, from the
    ``INFERNAL1/`` line through the ``//`` that closes the HMM) out of the library into ``out_dir/<acc>.cm``,
    reading the library once, line by line. Returns per accession its NAME, CLEN and GA. Refuses a model
    the library lacks, holds twice, or that has no GA."""
    want = set(accessions)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    found: dict[str, dict] = {}
    buf: list[str] = []
    head: dict[str, str] = {}
    closes = 0
    with _open_text(library) as fh:
        for line in fh:
            if line.startswith("INFERNAL1/"):
                buf, head, closes = [line], {}, 0
                continue
            if not buf:
                continue
            buf.append(line)
            if closes == 0:
                parts = line.split(None, 1)
                if parts and parts[0] in ("NAME", "ACC", "CLEN", "GA") and parts[0] not in head:
                    head[parts[0]] = parts[1].split()[0]
            if line.startswith("//"):
                closes += 1
                if closes == 2:
                    acc = head.get("ACC", head.get("NAME", ""))
                    if acc in want:
                        if acc in found:
                            raise ValueError(f"{library}: model {acc} appears twice")
                        if "GA" not in head:
                            raise ValueError(f"{library}: model {acc} has no GA cutoff")
                        (out_dir / f"{acc}.cm").write_text("".join(buf), encoding="ascii")
                        found[acc] = {"name": head["NAME"], "clen": int(head["CLEN"]), "ga": float(head["GA"])}
                    buf, head, closes = [], {}, 0
    missing = sorted(want - set(found))
    if missing:
        raise ValueError(f"{library}: {len(missing)} model(s) not in the library, first {missing[0]}")
    return found


def infernal_version() -> str:
    exe = shutil.which("cmemit")
    if exe is None:
        raise FileNotFoundError("cmemit (Infernal) is not on PATH; rna_motifs.source 'rfam' needs it")
    done = subprocess.run([exe, "-h"], capture_output=True, text=True, check=True)
    m = re.search(r"INFERNAL (\S+)", done.stdout)
    return m.group(1) if m else "unknown"


def cmemit(cm_file: Path, n: int, seed: int, *, exp: float = 1.0, local: bool = False) -> list[Emission]:
    """``n`` emissions of the one model in ``cm_file`` by ``cmemit -N n -a --dna --seed seed``."""
    if not 1 <= int(seed) <= 2**31 - 1:
        raise ValueError(f"cmemit seed must be in [1, 2^31 - 1] (0 means arbitrary), got {seed}")
    exe = shutil.which("cmemit")
    if exe is None:
        raise FileNotFoundError("cmemit (Infernal) is not on PATH")
    argv = [exe, "-N", str(int(n)), "-a", "--dna", "--seed", str(int(seed))]
    if float(exp) != 1.0:
        argv += ["--exp", f"{float(exp):g}"]
    if local:
        argv.append("-l")
    done = subprocess.run(argv + [str(cm_file)], capture_output=True, text=True, stdin=subprocess.DEVNULL)
    if done.returncode != 0:
        raise RuntimeError(f"cmemit exited with code {done.returncode}: {done.stderr[-2000:]}")
    em = emissions_from_stockholm(done.stdout)
    if len(em) != int(n):
        raise RuntimeError(f"cmemit gave {len(em)} sequences, asked for {n}")
    return em


# --- evolution (§2.5 step 3) ----------------------------------------------------------------------------

_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(BASES):
    _LUT[ord(_b)] = _i
_PAIR_A = np.array([_CODE[p[0]] for p in CANONICAL], dtype=np.uint8)
_PAIR_B = np.array([_CODE[p[1]] for p in CANONICAL], dtype=np.uint8)


def to_codes(seq: str) -> np.ndarray:
    """ACGT -> 0..3."""
    c = _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    if (c == 255).any():
        raise ValueError("a sequence with letters outside ACGT")
    return c


def to_str(c: np.ndarray) -> str:
    return np.frombuffer(b"ACGT", dtype=np.uint8)[c].tobytes().decode("ascii")


def evolve_motif(motif: np.ndarray, pairs: np.ndarray, unpaired: np.ndarray, x: float, comp_cum: np.ndarray,
                 pair_cum: np.ndarray, rng: np.random.Generator) -> None:
    """One branch of expected events ``x`` per position (unpaired) or per pair, in place on the base codes
    ``motif``: an unpaired position's new base from ``comp_cum`` (cumulative over ACGT), a pair's new pair
    from ``pair_cum`` (cumulative over `CANONICAL`). Draw order: unpaired hits, their bases, pair hits, their
    pairs."""
    p = -math.expm1(-float(x))
    hit = unpaired[rng.random(len(unpaired)) < p]
    if hit.size:
        motif[hit] = np.minimum(np.count_nonzero(comp_cum[None, :] <= rng.random(hit.size)[:, None], axis=1), 3)
    if len(pairs):
        ph = np.nonzero(rng.random(len(pairs)) < p)[0]
        if ph.size:
            k = np.minimum(np.count_nonzero(pair_cum[None, :] <= rng.random(ph.size)[:, None], axis=1), 5)
            motif[pairs[ph, 0]] = _PAIR_A[k]
            motif[pairs[ph, 1]] = _PAIR_B[k]


def is_canonical(a: str, b: str) -> bool:
    return a + b in CANONICAL


__all__ = ["BASES", "CANONICAL", "Emission", "MotifModel", "parse_wuss", "dot_bracket", "parse_stockholm",
           "emissions_from_stockholm", "synthetic_models", "synthetic_emit", "extract_rfam", "infernal_version",
           "cmemit", "evolve_motif", "to_codes", "to_str", "is_canonical"]

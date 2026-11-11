"""Profile search rounds and the iteration controller (operations O-SEARCH, O-ITER; SEARCH_PLAN §4.1, SD-10).

Every round is an `hmmsearch` of one profile against the whole target at the fixed search
space (Z_prot = domZ = 1e8). Round 1 searches with the seed: an HMM as given (from a pack
library or a user's file), or a profile built from a sequence (what phmmer and the first round
of jackhmmer do), from an MSA, or from unaligned seed sequences (`build_from_sequences`).

The target is never held whole (card W-19): a search reads it from the pack's shards one block of
`BLOCK_RESIDUES` residues at a time (`stream_target`), each block searched across the run's threads, and
keeps only the hits (and, for a round that builds a next profile, its hits' sequences). At the fixed Z a
target's scores and E-values do not depend on what it is searched with, so this is one search of the target.

Per round k:

  1. Keep targets whose sequence bitscore is >= `inclusion_bits`; for each, every domain
     with bitscore >= `domain_bits`, and always its best domain.
  2. Curation, as labels with reasons (never numbers): `query_coverage` when the best
     domain covers less than `min_query_coverage` of the profile; `partial` under the
     `exclude` policy. Dropped targets are logged and never enter the next profile.
  3. From round 2, checks on the targets new in this round:
       drift_alarm   a library family not on the seed takes more than `alarm_fraction` of
                     at least `alarm_min` new targets (stored library hits, l1_domains)
       context_lost  fewer than `context_retention` of the new targets have a neighbor
                     family seen around the round-1 targets
     Either reverts: round k is logged with label `drop`, and round k-1 is final.
  4. `no_new_reps` stops when a round adds no target; `iterations.max` bounds the rounds.
  5. The next profile is built from the kept domain envelopes only (aligned to this round's
     profile with hmmalign, trimmed to match states), plus the seed sequence, so a fused
     partner domain never enters the alignment.

At S0 there are no NR50 representatives, so "new representatives" are new targets.
"""

from __future__ import annotations

import io
import math
import random
import re
from dataclasses import dataclass, field
from typing import Callable, Iterator

import pyhmmer
from pyhmmer.plan7 import Background, Builder, TraceAligner

from biobox.l1.domains import DOM_Z, Z_PROT, iter_hmms

AMINO = pyhmmer.easel.Alphabet.amino()
# residues per block of a target read from the pack (card W-19): about 100 to 150 MB of digitized sequences with
# their Easel records, for proteins of about 300 aa (estimate); a search holds one block at a time
BLOCK_RESIDUES = 2**25
BUILDER_SEED = 42
FALLBACK_FAMILY_BITS = 25.0  # family calls on the seed for models without GA cutoffs


@dataclass
class Target:
    """What a search reads: `size` sequences, given as blocks of digitized sequences, read afresh by each call of
    `blocks()` (one pass over the target); `seqs` holds the text of the sequences a search kept
    (`search(keep_sequences=True)`: the hits of a round that builds the next profile), or all of them for a target
    loaded whole (`load_target`). `blocks_read` counts the blocks of the last pass."""
    kind: str  # protein | dark_orf
    size: int
    _blocks: Callable[[], Iterator[pyhmmer.easel.DigitalSequenceBlock]]
    seqs: dict[bytes, str] = field(default_factory=dict)
    blocks_read: int = 0

    def blocks(self) -> Iterator[pyhmmer.easel.DigitalSequenceBlock]:
        self.blocks_read = 0
        for block in self._blocks():
            self.blocks_read += 1
            yield block


def load_target(kind: str, records: Iterator[tuple[str, str]]) -> Target:
    """A small target held whole, in one block (tests and small inputs; a search reads a pack's shards with
    `stream_target`)."""
    seqs = {bytes.fromhex(n): s for n, s in records}
    uids = sorted(seqs)
    block = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=u.hex().encode(), sequence=seqs[u]) for u in uids).digitize(AMINO)
    return Target(kind, len(uids), lambda: iter([block]), seqs)


def stream_target(kind: str, paths: list, size: int, block_residues: int | None = None) -> Target:
    """A target read from its pack shards (`paths`, in pack order) a block of about `block_residues` residues
    (default `BLOCK_RESIDUES`) at a time, never held whole (card W-19); `size` is its number of sequences (the pack
    manifest's `n_seq`: a searchable target holds each name once). A block is searched by pyhmmer across the
    search's threads, and at the fixed Z (`Z_PROT`, `DOM_Z`) a sequence's scores and E-values do not depend on the
    other sequences searched with it, so the hits are those of one search of the whole target."""
    from biobox.search import fastz

    n = int(block_residues or BLOCK_RESIDUES)

    def blocks():
        for p in paths:
            yield from fastz.blocks(p, AMINO, n)

    return Target(kind, size, blocks)


def text(x) -> str:
    """pyhmmer names and accessions as str (0.12 returns str, earlier releases bytes)."""
    return x.decode() if isinstance(x, bytes) else (x or "")


def _digital(name: str, seq: str) -> pyhmmer.easel.DigitalSequence:
    return pyhmmer.easel.TextSequence(name=name.encode(), sequence=seq).digitize(AMINO)


def _hmm_text(hmm) -> bytes:
    buf = io.BytesIO()
    hmm.write(buf)
    # HMMER writes the build time (DATE); drop it so a profile's bytes depend only on its content
    return b"".join(l for l in buf.getvalue().splitlines(keepends=True) if not l.startswith(b"DATE"))


def _msa_text(msa) -> bytes:
    buf = io.BytesIO()
    (msa.textize() if hasattr(msa, "textize") else msa).write(buf, "stockholm")
    return buf.getvalue()


def build_from_sequence(name: str, seq: str):
    hmm, _, _ = Builder(AMINO, seed=BUILDER_SEED).build(_digital(name, seq), Background(AMINO))
    hmm.name = name.encode()
    return hmm


def build_from_msa(name: str, msa):
    msa.name = name.encode()
    hmm, _, _ = Builder(AMINO, seed=BUILDER_SEED).build_msa(msa, Background(AMINO))
    hmm.name = name.encode()
    return hmm


def build_from_sequences(name: str, records: list[tuple[str, str]]):
    """A profile from unaligned seeds, as hmmalign + hmmbuild would make it with no multiple aligner:
    a single-sequence profile of the first seed, every seed aligned to it (trimmed to the model's
    match states) and the profile built from that alignment. Returns (profile, alignment)."""
    first = build_from_sequence(f"{name}_first", records[0][1])
    block = pyhmmer.easel.DigitalSequenceBlock(AMINO, [_digital(n, s) for n, s in records])
    msa = pyhmmer.hmmer.hmmalign(first, block, trim=True, digitize=True, cpus=1)
    return build_from_msa(name, msa), msa


def read_msa(path) -> pyhmmer.easel.DigitalMSA:
    with pyhmmer.easel.MSAFile(str(path), digital=True, alphabet=AMINO) as fh:
        return fh.read()


def _model_keys(hmm) -> set[str]:
    acc = text(hmm.accession)
    return {k for k in (acc, acc.split(".")[0], text(hmm.name)) if k}


def read_hmm_file(data: bytes, model: str | None = None):
    """The one HMMER3 profile in a user's file, or the one `model` names (name or accession, with or
    without version). ValueError says why when the file holds none, several, or a non-protein one."""
    try:
        with pyhmmer.plan7.HMMFile(io.BytesIO(data)) as fh:
            hmms = list(fh)
    except (ValueError, EOFError, OSError, RuntimeError) as e:
        raise ValueError(f"not an HMMER3 profile file ({e})") from None
    if not hmms:
        raise ValueError("not an HMMER3 profile file (no model in it)")
    names = [text(h.accession) or text(h.name) for h in hmms]
    if model is None:
        if len(hmms) > 1:
            raise ValueError(f"{len(hmms)} models in the file ({', '.join(names[:5])}"
                             f"{', ...' if len(names) > 5 else ''}); name one with `model:` (name or accession)")
        hmm = hmms[0]
    else:
        match = [h for h in hmms if model in _model_keys(h)]
        if not match:
            raise ValueError(f"no model {model!r} in the file (it holds {', '.join(names[:5])}"
                             f"{', ...' if len(names) > 5 else ''})")
        if len(match) > 1:
            raise ValueError(f"model {model!r} names {len(match)} models in the file; give a unique accession")
        hmm = match[0]
    if not hmm.alphabet.is_amino():
        raise ValueError(f"model {text(hmm.name)} is not an amino-acid profile; protein targets need an "
                         "amino-acid profile")
    return hmm


_NUCLEOTIDE = set("ACGTUN")


def read_seed_fasta(data: bytes) -> list[tuple[str, str]]:
    """(name, sequence) of a protein FASTA of seeds, in file order; ValueError says why it is unusable:
    fewer than 2 sequences, a duplicate name, an empty sequence, a character outside A-Z, or only
    nucleotide letters."""
    try:
        with pyhmmer.easel.SequenceFile(io.BytesIO(data), format="fasta") as fh:
            recs = [(text(s.name), s.sequence) for s in fh]
    except (ValueError, EOFError, OSError, RuntimeError) as e:
        raise ValueError(f"not a FASTA file ({e})") from None
    if len(recs) < 2:
        raise ValueError(f"{len(recs)} sequence(s); a profile from sequences needs at least 2 "
                         "(give one sequence as a `sequence` seed)")
    out, seen = [], set()
    for name, seq in recs:
        if name in seen:
            raise ValueError(f"duplicate sequence name {name!r}")
        seen.add(name)
        seq = re.sub(r"\s+", "", seq).upper().rstrip("*")
        if not seq:
            raise ValueError(f"sequence {name!r} is empty")
        if not re.fullmatch(r"[A-Z]+", seq):
            raise ValueError(f"sequence {name!r} has characters outside A-Z")
        out.append((name, seq))
    if set("".join(s for _, s in out)) <= _NUCLEOTIDE:
        raise ValueError("every sequence holds only nucleotide letters (ACGTUN); give protein sequences")
    return out


def find_model(paths: list, key: str):
    """A model by accession (with or without version) or name, from the pack's libraries."""
    for p in paths:
        # read each library to the end: leaving iter_hmms open while holding a model crashes pyhmmer
        for hmm in list(iter_hmms(p)):
            acc = text(hmm.accession)
            if key in (acc, acc.split(".")[0], text(hmm.name)):
                return hmm
    return None


def families_of(seq: str, model_paths: list, threads: int) -> set[str]:
    """Library families on one sequence, at each model's GA cutoff (25 bits without one)."""
    out = set()
    target = pyhmmer.easel.DigitalSequenceBlock(AMINO, [_digital("seed", seq)])
    for p in model_paths:
        for hits in pyhmmer.hmmsearch(list(iter_hmms(p)), target, cpus=threads, Z=Z_PROT, domZ=DOM_Z):
            q = hits.query
            ga = q.cutoffs.gathering if q.cutoffs.gathering_available() else None
            cut = ga[0] if ga else FALLBACK_FAMILY_BITS
            if any(h.score >= cut for h in hits):
                out.add(text(q.accession) or text(q.name))
    return out


def search(hmm, target: Target, inclusion_bits: float, domain_bits: float, threads: int, *,
           keep_sequences: bool = False) -> dict[bytes, dict]:
    """Included targets -> {seq_bitscore, seq_evalue_zref, domains: [...]} for one profile, one block of the
    target at a time (`Target.blocks`). With `keep_sequences`, the text of each included target's sequence is kept
    in `target.seqs` (for the next round's alignment), read from the block it was found in."""
    out = {}
    for block in target.blocks():
        found = []
        for hits in pyhmmer.hmmsearch([hmm], block, cpus=threads, Z=Z_PROT, domZ=DOM_Z,
                                      T=inclusion_bits, incT=inclusion_bits):
            for hit in hits:
                if hit.score < inclusion_bits:
                    continue
                doms = sorted(hit.domains, key=lambda d: -d.score)
                keep = [d for i, d in enumerate(doms) if i == 0 or d.score >= domain_bits]
                uid = bytes.fromhex(text(hit.name))
                found.append(uid)
                out[uid] = {
                    "seq_bitscore": hit.score, "seq_evalue_zref": hit.evalue,
                    "domains": sorted(({"bitscore": d.score, "evalue_zref": d.i_evalue, "env_start": d.env_from - 1,
                                        "env_end": d.env_to, "ali_start": d.alignment.target_from - 1,
                                        "ali_end": d.alignment.target_to, "hmm_from": d.alignment.hmm_from - 1,
                                        "hmm_to": d.alignment.hmm_to} for d in keep), key=lambda d: d["env_start"]),
                }
        if keep_sequences:
            _keep(block, {u.hex() for u in found if u not in target.seqs}, target.seqs)
    return out


def _keep(block, names: set[str], seqs: dict[bytes, str]) -> None:
    """The text of the sequences of `block` named in `names`, into `seqs`."""
    if not names:
        return
    for s in block:
        n = text(s.name)
        if n in names:
            seqs[bytes.fromhex(n)] = s.textize().sequence


def tail_expectation(hmm, bits: float, size: int) -> float:
    """Expected number of chance targets scoring >= bits: size * P(S >= bits) from the profile's
    forward-score exponential tail (HMMER's tau, lambda)."""
    ep = hmm.evalue_parameters
    return size * min(1.0, math.exp(-ep.f_lambda * (bits - ep.f_tau)))


@dataclass
class Settings:
    inclusion_bits: float
    domain_bits: float
    max_rounds: int
    stop: tuple[str, ...]
    min_query_coverage: float = 0.0
    alarm_fraction: float = 0.5
    alarm_min: int = 3
    context_retention: float = 0.5
    partial_policy: str = "flag"
    cpus: int = 1


@dataclass
class Round:
    k: int
    hmm_text: bytes
    msa_text: bytes | None
    n_seq: int
    n_eff: float | None
    hits: dict[bytes, dict]
    labels: dict[bytes, tuple[str, str | None]]
    new: list[bytes]
    checks: dict = field(default_factory=dict)


def _label(uid: bytes, hit: dict, model_len: int, s: Settings, partial: set[bytes]) -> tuple[str, str | None]:
    if s.partial_policy == "exclude" and uid in partial:
        return "drop", "partial"
    best = max(hit["domains"], key=lambda d: d["bitscore"])
    cov = (best["hmm_to"] - best["hmm_from"]) / model_len if model_len else 0
    if cov < s.min_query_coverage:
        return "drop", f"query_coverage {cov:.2f} < {s.min_query_coverage:.2f}"
    return "keep", None


def iterate(seed_hmm, seed_msa_text: bytes | None, seed_seq: str | None, target: Target, s: Settings,
            seed_families: set[str], families: Callable[[list[bytes]], dict[bytes, set[str]]] | None,
            context: Callable[[list[bytes]], dict[bytes, set[str]]] | None,
            partial: set[bytes], name: str) -> Iterator[Round]:
    """Yield each round; the caller reads `.checks` on the last one for the stop reason."""
    hmm, msa_text = seed_hmm, seed_msa_text
    seen: set[bytes] = set()
    round1_context: set[str] | None = None
    for k in range(1, s.max_rounds + 1):
        # a round that may build a next profile keeps its hits' sequences: the alignment reads them
        hits = search(hmm, target, s.inclusion_bits, s.domain_bits, threads=s.cpus, keep_sequences=k < s.max_rounds)
        labels = {u: _label(u, h, hmm.M, s, partial) for u, h in hits.items()}
        kept = sorted(u for u, (lab, _) in labels.items() if lab == "keep")
        new = [u for u in kept if u not in seen]
        rnd = Round(k, _hmm_text(hmm), msa_text, int(hmm.nseq or 1), float(hmm.nseq_effective or 1.0),
                    hits, labels, new)
        checks = rnd.checks
        if k == 1 and context is not None and "context_lost" in s.stop:
            round1_context = set().union(*context(kept).values()) if kept else set()
        if k > 1:
            if "drift_alarm" in s.stop and families is not None and len(new) >= s.alarm_min:
                fam = families(new)
                counts: dict[str, int] = {}
                for u in new:
                    for f in fam.get(u, ()):
                        counts[f] = counts.get(f, 0) + 1
                foreign = {f: c / len(new) for f, c in counts.items() if f not in seed_families}
                worst = max(foreign.items(), key=lambda x: (x[1], x[0]), default=None)
                checks["foreign_share"] = worst
                if worst and worst[1] > s.alarm_fraction:
                    checks["stop"] = f"drift_alarm: {worst[0]} in {worst[1]:.2f} of {len(new)} new targets"
            if "stop" not in checks and round1_context is not None and new:
                ctx = context(new)
                retained = sum(1 for u in new if ctx.get(u, set()) & round1_context) / len(new)
                checks["context_retention"] = retained
                if retained < s.context_retention:
                    checks["stop"] = f"context_lost: {retained:.2f} of {len(new)} new targets keep round-1 context"
            if "stop" in checks:
                checks["revert"] = True
                rnd.labels = {u: ("drop", checks["stop"].split(":")[0]) for u in hits}
                yield rnd
                return
            if not new and "no_new_reps" in s.stop:
                checks["stop"] = "no_new_reps"
                yield rnd
                return
        seen |= set(kept)
        if k == s.max_rounds:
            checks["stop"] = "max_rounds"
            yield rnd
            return
        if not kept:
            checks["stop"] = "no_hits"
            yield rnd
            return
        yield rnd
        # next profile: kept envelopes only, aligned to this round's profile
        seqs = []
        for u in kept:
            for d in hits[u]["domains"]:
                seqs.append(_digital(f"{u.hex()}/{d['env_start']}-{d['env_end']}",
                                     target.seqs[u][d["env_start"]:d["env_end"]]))
        if seed_seq is not None:
            seqs.append(_digital("seed", seed_seq))
        block = pyhmmer.easel.DigitalSequenceBlock(AMINO, seqs)
        aligner = TraceAligner()
        msa = aligner.align_traces(hmm, block, aligner.compute_traces(hmm, block), trim=True, digitize=True)
        msa_text = _msa_text(msa)
        hmm = build_from_msa(f"{name}_r{k + 1}", msa)


def shuffled(seq: str, seed: int) -> str:
    aa = list(seq)
    random.Random(seed).shuffle(aa)
    return "".join(aa)


def negative_control(seed_seq: str, target: Target, s: Settings, seed: int) -> dict:
    """MR-09: a composition-preserving shuffle of the seed, searched as round 1."""
    hmm = build_from_sequence("shuffled_seed", shuffled(seed_seq, seed))
    n = len(search(hmm, target, s.inclusion_bits, s.domain_bits, threads=s.cpus))
    exp = tail_expectation(hmm, s.inclusion_bits, target.size)
    return {"kind": "shuffled_seed", "n_hits": n, "expectation": exp, "allowed": math.ceil(exp),
            "pass": n <= math.ceil(exp)}

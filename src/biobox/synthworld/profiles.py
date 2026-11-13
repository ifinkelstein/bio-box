"""The world's HMM panel and the profile root sampler (planning/tests/synthworld_spec.md §2.1 S3, §2.3 "Sites",
"Root"; card F-02 session 3).

A profile is a Plan7 core model: per match state ``k = 1..M`` an emission row over ``AA`` (HMMER order), and
per node ``k = 0..M`` the transition probabilities M->I, I->I, M->D, D->D (node 0 is the begin state; node M
has none). Two panel sources (preset ``sequences.profiles.source``):

- ``pfam``: HMMs read from a Pfam HMM file (``Pfam-A.hmm`` or ``.gz``) with pyhmmer, in file order. An HMM
  is selected when its match length lies in ``pfam_match_range``, it has gathering cutoffs, and the first
  uniform of its stream ``S3/panel/<NAME>`` is below ``pfam_select_p``; reading stops when the panel holds
  ``n_scan + n_dark`` HMMs (so only the head of the file is read). Of those, the ``n_dark`` with the
  smallest second uniform are held out of the scan panel (dark), and of the rest the ``n_promiscuous``
  with the smallest third uniform are the promiscuous-domain traps.
- ``synthetic`` (the micro default, so the world builds where no Pfam file is present): stream
  ``S3/panel`` draws ``n_scan + n_dark`` profiles. Match length: log-normal (``synthetic_match``: median,
  sigma, clipped to [min, max]); each site's row ~ Dirichlet(A_s x background), A_s log-normal
  (``synthetic_concentration``); transitions constant (``synthetic_transitions``), none out of node 0
  into an insert and none out of node M. Profiles ``SYN0001``.. ; the last ``n_dark`` are dark and the first
  ``n_promiscuous`` the promiscuous traps.

`sample_domain` is the root sampler of §2.3: a walk of the core model from the begin state (inserts and
deletes by the transitions, as ``hmmemit`` samples the core model), match states emitting from their row
and insert states from the world background (§2.3: insert residues use the background composition).
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from biobox.synthworld import seeds

#: HMMER 3.4's amino-acid null model in ``AA`` order (pyhmmer 0.12.3
#: ``Background(Alphabet.amino()).residue_frequencies``, read on wilkcomp02 2026-11-12), normalized
_BG = np.array([0.0787945, 0.01516, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092,
                0.0594422, 0.0963728, 0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364,
                0.0540687, 0.0673417, 0.0114135, 0.0304133])
BACKGROUND = _BG / _BG.sum()


@dataclass(frozen=True, eq=False)
class Profile:
    index: int
    name: str
    accession: str
    source: str              # pfam | synthetic
    match: np.ndarray        # (M, 20) emission rows of match states 1..M
    t_mi: np.ndarray         # (M + 1,) per node 0..M
    t_ii: np.ndarray
    t_md: np.ndarray
    t_dd: np.ndarray
    ga: tuple[float, float] | None
    in_scan_panel: bool
    promiscuous: bool

    @property
    def M(self) -> int:
        return int(self.match.shape[0])

    @property
    def saturation(self) -> float:
        """Mean over match states of the sum of squared emission probabilities: the identity of two
        independent draws at match sites."""
        return float((self.match ** 2).sum(axis=1).mean())


def sample_domain(profile: Profile, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """One draw of the core model: residue codes (0..19) and, per residue, its match state (1..M) or 0 for an
    insert-state residue. Inserts emit from `BACKGROUND`."""
    M = profile.M
    cum_match = np.cumsum(profile.match, axis=1)
    cum_match[:, -1] = 1.0
    cum_bg = np.cumsum(BACKGROUND)
    cum_bg[-1] = 1.0
    nodes: list[int] = []
    k, state = 0, "M"                      # node 0, the begin state, behaves as a match state
    while True:
        if state == "M":
            u = rng.random()
            if u < profile.t_mi[k]:
                state = "I"
                nodes.append(0)
                continue
            k += 1
            state = "D" if u < profile.t_mi[k - 1] + profile.t_md[k - 1] else "M"
        elif state == "I":
            if rng.random() < profile.t_ii[k]:
                nodes.append(0)
                continue
            k += 1
            state = "M"
        else:  # D
            k += 1
            state = "D" if rng.random() < profile.t_dd[k - 1] else "M"
        if k > M:
            break
        if state == "M":
            nodes.append(k)
        if k == M and state == "D":
            break
    node = np.array(nodes, dtype=np.int32)
    u = rng.random(len(node))
    aa = np.empty(len(node), dtype=np.int8)
    is_match = node > 0
    aa[is_match] = (cum_match[node[is_match] - 1] <= u[is_match, None]).sum(axis=1).clip(max=19)
    aa[~is_match] = (cum_bg[None, :] <= u[~is_match, None]).sum(axis=1).clip(max=19)
    return aa, node


# --- panels --------------------------------------------------------------------------------------------

def _transitions(M: int, t: dict) -> tuple[np.ndarray, ...]:
    mi = np.full(M + 1, float(t["mi"]))
    ii = np.full(M + 1, float(t["ii"]))
    md = np.full(M + 1, float(t["md"]))
    dd = np.full(M + 1, float(t["dd"]))
    mi[0] = 0.0                            # no insert before the first match state: tails cover it
    mi[M] = md[M] = 0.0                    # node M ends the model
    dd[M] = 0.0
    return mi, ii, md, dd


def synthetic_panel(seed: int, cfg: dict) -> list[Profile]:
    """``n_scan + n_dark`` synthetic profiles from stream ``S3/panel`` (module docstring)."""
    rng = seeds.generator(seed, "S3/panel")
    n_scan, n_dark, n_prom = int(cfg["n_scan"]), int(cfg["n_dark"]), int(cfg["n_promiscuous"])
    sm, sc = cfg["synthetic_match"], cfg["synthetic_concentration"]
    out = []
    for i in range(n_scan + n_dark):
        M = int(np.clip(round(float(sm["median"]) * np.exp(float(sm["sigma"]) * rng.standard_normal())),
                        int(sm["min"]), int(sm["max"])))
        conc = float(sc["median"]) * np.exp(float(sc["sigma"]) * rng.standard_normal(M))
        g = rng.gamma(conc[:, None] * BACKGROUND[None, :]) + 1e-12
        match = g / g.sum(axis=1, keepdims=True)
        mi, ii, md, dd = _transitions(M, cfg["synthetic_transitions"])
        out.append(Profile(i, f"SYN{i + 1:04d}", f"SYN{i + 1:04d}.1", "synthetic", match, mi, ii, md, dd,
                           None, i < n_scan, i < n_prom))
    return out


def _text(v) -> str:
    return v.decode() if isinstance(v, bytes) else str(v)


def pfam_panel(path, seed: int, cfg: dict) -> list[Profile]:
    """The panel read from a Pfam HMM file (module docstring)."""
    from pyhmmer.plan7 import HMMFile

    n_scan, n_dark, n_prom = int(cfg["n_scan"]), int(cfg["n_dark"]), int(cfg["n_promiscuous"])
    lo, hi = (int(x) for x in cfg["pfam_match_range"])
    p_sel = float(cfg["pfam_select_p"])
    path = Path(path)
    picked = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as fh, HMMFile(fh) as hmms:
        for hmm in hmms:
            if not (lo <= hmm.M <= hi) or not hmm.cutoffs.gathering_available():
                continue
            name = _text(hmm.name)
            u = seeds.generator(seed, f"S3/panel/{name}").random(3)
            if u[0] >= p_sel:
                continue
            tp = np.asarray(hmm.transition_probabilities, dtype=np.float64)
            match = np.asarray(hmm.match_emissions, dtype=np.float64)[1:, :20].copy()
            match /= match.sum(axis=1, keepdims=True)
            picked.append((name, _text(hmm.accession), match, tp, tuple(float(x) for x in hmm.cutoffs.gathering),
                           float(u[1]), float(u[2])))
            if len(picked) == n_scan + n_dark:
                break
    if len(picked) < n_scan + n_dark:
        raise ValueError(f"{path.name}: only {len(picked)} HMMs pass the panel rule, need {n_scan + n_dark}")
    dark = {i for i, _ in sorted(enumerate(picked), key=lambda x: (x[1][5], x[1][0]))[:n_dark]}
    scan = [i for i in range(len(picked)) if i not in dark]
    prom = {i for i in sorted(scan, key=lambda i: (picked[i][6], picked[i][0]))[:n_prom]}
    out = []
    for i, (name, acc, match, tp, ga, _, _) in enumerate(picked):
        # pyhmmer columns: MM MI MD IM II DM DD, rows nodes 0..M
        out.append(Profile(i, name, acc, "pfam", match, tp[:, 1].copy(), tp[:, 4].copy(), tp[:, 2].copy(),
                           tp[:, 6].copy(), ga, i not in dark, i in prom))
    return out


def panel(cfg: dict, seed: int, pfam_hmm=None) -> list[Profile]:
    """The world's panel from the preset's ``sequences.profiles`` section."""
    src = cfg["source"]
    if src == "synthetic":
        return synthetic_panel(seed, cfg)
    if src == "pfam":
        if pfam_hmm is None:
            raise ValueError("sequences.profiles.source is 'pfam': pass the Pfam HMM file (pfam_hmm=...)")
        return pfam_panel(pfam_hmm, seed, cfg)
    raise ValueError(f"sequences.profiles.source is 'synthetic' or 'pfam', got {src!r}")

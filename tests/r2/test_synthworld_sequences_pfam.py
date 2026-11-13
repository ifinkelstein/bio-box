"""Ring R2: tool parity of synthworld's Pfam-seeded families against the HMMER binaries (synthworld_spec §2.3
cross-check, §7.2 hmmsearch; card F-02 session 3).

Needs the HMMER binaries (`hmmfetch`, `hmmemit`, `hmmsearch`, in the bio-box image), ZOMBI2, and Pfam's HMM
file named by ``BIOBOX_PFAM_HMM`` (on the pod, read-only:
``/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz``); skips otherwise, so it skips
in GitHub Actions. The models are fetched from the Pfam file by ``hmmfetch``, not by the generator's reader.

What the profile-F81 model of §2.3 predicts, and what these tests hold it to:
- at every site the equilibrium is the source HMM's match emission row, and the root is drawn from it, so
  each copy is, on its own, a draw of the same kind as the root (§2.3: "at T -> infinity the members become
  independent hmmemit-like draws"). A copy's score against the source HMM therefore does not drift with its
  divergence from the root; divergence lowers identity between copies, not detectability. The detection
  test checks both ends: copies close to the root are found above the gathering cutoff, and the score
  does not fall systematically with divergence.
- §2.3 cross-check: for 20 HMMs, hmmsearch scores of the numpy root sampler's draws and of `hmmemit`'s draws
  agree (two-sample KS p > 0.01 at fixed seeds).
"""

from __future__ import annotations

import collections
import csv
import importlib.util
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from biobox.synthworld import backbone, codons, profiles, seeds, sequences

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]

PFAM = os.environ.get("BIOBOX_PFAM_HMM", "")
_missing = [b for b in ("hmmfetch", "hmmemit", "hmmsearch") if shutil.which(b) is None]
needs_pfam = pytest.mark.skipif(
    bool(_missing) or not (PFAM and Path(PFAM).is_file()) or importlib.util.find_spec("zombi2") is None,
    reason=f"needs the HMMER binaries (missing: {_missing or 'none'}), Pfam's HMM file in BIOBOX_PFAM_HMM "
           f"({PFAM or 'unset'}) and zombi2; runs on the pod in the bio-box image")


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fasta(path: Path) -> dict[str, str]:
    out, name, chunks = {}, None, []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line)
    if name is not None:
        out[name] = "".join(chunks)
    return out


def _tblout(path: Path) -> dict[tuple[str, str], float]:
    """(target, query name) -> full-sequence bit score, from hmmsearch --tblout."""
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        out[(f[0], f[2])] = float(f[5])
    return out


def _ks_p(a, b) -> float:
    """Two-sample Kolmogorov-Smirnov p-value, asymptotic (Numerical Recipes kstwo/probks), written here."""
    a, b = np.sort(a), np.sort(b)
    allv = np.concatenate([a, b])
    d = float(np.max(np.abs(np.searchsorted(a, allv, "right") / len(a) - np.searchsorted(b, allv, "right") / len(b))))
    ne = len(a) * len(b) / (len(a) + len(b))
    lam = (math.sqrt(ne) + 0.12 + 0.11 / math.sqrt(ne)) * d
    s, sign = 0.0, 1.0
    for k in range(1, 101):
        term = sign * 2.0 * math.exp(-2.0 * k * k * lam * lam)
        s += term
        if abs(term) < 1e-10:
            break
        sign = -sign
    return min(1.0, max(0.0, s))


def _write_fasta(path: Path, recs) -> None:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in recs), encoding="ascii")


@pytest.fixture(scope="module")
def pfam_world(tmp_path_factory):
    out = tmp_path_factory.mktemp("micro_pfam") / "world"
    params = backbone.load_preset("micro", overrides={"sequences": {"profiles": {"source": "pfam"}}})
    backbone.build(params, params["default_seed"], out)
    sw = sequences.build(out, pfam_hmm=PFAM)
    return sw, out


@pytest.fixture(scope="module")
def fetched(pfam_world, tmp_path_factory):
    """Every panel HMM the tests use, fetched from the Pfam file by hmmfetch (one pass over the file)."""
    _, out = pfam_world
    panel = _tsv(out / "sequences/panel.tsv")
    d = tmp_path_factory.mktemp("hmms")
    names = d / "names.txt"
    names.write_text("".join(p["name"] + "\n" for p in panel if p["in_scan_panel"] == "true"))
    hmm = d / "panel.hmm"
    with open(hmm, "wb") as fh:
        subprocess.run(["hmmfetch", "-f", PFAM, str(names)], check=True, stdout=fh)
    return hmm, d


def _one(hmm_all: Path, name: str, dst: Path) -> Path:
    path = dst / f"{name}.hmm"
    with open(path, "wb") as fh:
        subprocess.run(["hmmfetch", str(hmm_all), name], check=True, stdout=fh)
    return path


@needs_pfam
def test_pfam_panel_matches_the_file(pfam_world, fetched):
    """The panel's accessions, lengths and gathering cutoffs are those of the HMMs hmmfetch returns."""
    _, out = pfam_world
    hmm, _ = fetched
    text = hmm.read_text()
    got = {}
    for block in text.split("//\n"):
        if "NAME" not in block:
            continue
        kv = {}
        for line in block.splitlines():
            parts = line.split(None, 1)
            if parts and parts[0] in ("NAME", "ACC", "LENG", "GA"):
                kv[parts[0]] = parts[1].strip()
        got[kv["NAME"]] = kv
    panel = [p for p in _tsv(out / "sequences/panel.tsv") if p["in_scan_panel"] == "true"]
    assert len(panel) == len(got) == int(backbone.load_preset("micro")["sequences"]["profiles"]["n_scan"])
    for p in panel:
        kv = got[p["name"]]
        assert kv["ACC"] == p["accession"] and int(kv["LENG"]) == int(p["match_length"])
        ga = [float(x) for x in kv["GA"].rstrip(";").split()]
        assert (float(p["ga_seq"]), float(p["ga_dom"])) == pytest.approx((ga[0], ga[1]))


@needs_pfam
@pytest.mark.covers("F-02", "B-SIM")
def test_hmmsearch_finds_copies_near_the_root_and_scores_do_not_fall_with_divergence(pfam_world, fetched,
                                                                                       tmp_path):
    """Families with a scan-panel domain and at least 6 copies (up to 40 families). Truth: each copy's
    expected substitutions per site from the family root (`genes.tsv` subst_per_site) and the root itself
    (`roots.faa`). Assertions:
    (1) a copy's result tracks its family root's at low divergence and loses it as divergence grows:
    over all copies sorted by divergence, the nearest third agrees with its root on detection above the
    gathering cutoff (`hmmsearch --cut_ga`) more often than the farthest third, and differs from the
    root's bit score by less (median absolute difference). The root is one profile draw and, like
    `hmmemit` output, need not pass GA, so "every near copy is found" is not what the model predicts
    (first run: 21 of 26 copies within 0.1 substitutions per site found; second: 377 of 441 copies within
    0.25 agreed with their root, both on wilkcomp02 2026-11-12);
    (2) the profile-F81 model predicts no decline of a copy's score with divergence: splitting each
    family's copies at their median divergence, the far half scores lower on average in at most
    0.5 + 3 sd (binomial) of families."""
    _, out = pfam_world
    hmm_all, _ = fetched
    prot = _fasta(out / "sequences/proteins.faa")
    genes = _tsv(out / "sequences/genes.tsv")
    panel = {p["accession"]: p for p in _tsv(out / "sequences/panel.tsv")}
    by_fam = collections.defaultdict(list)
    for g in genes:
        by_fam[int(g["family"])].append(g)
    chosen = []
    for f in _tsv(out / "sequences/families.tsv"):
        if f["dark"] == "true" or len(by_fam[int(f["family"])]) < 6:
            continue
        accs = [a for a in f["domain_hmms"].split(";") if panel[a]["promiscuous"] == "false"]
        if accs:
            chosen.append((int(f["family"]), panel[accs[0]]["name"]))
        if len(chosen) == 40:
            break
    assert len(chosen) >= 20, len(chosen)
    q = tmp_path / "queries.hmm"
    names = sorted({c[1] for c in chosen})
    (tmp_path / "names.txt").write_text("".join(n + "\n" for n in names))
    with open(q, "wb") as fh:
        subprocess.run(["hmmfetch", "-f", str(hmm_all), str(tmp_path / "names.txt")], check=True, stdout=fh)
    roots = _fasta(out / "sequences/roots.faa")
    recs = [(g["copy_label"], prot[g["copy_label"]]) for fam, _ in chosen for g in by_fam[fam]]
    recs += [(f"F{fam:06d}", roots[f"F{fam:06d}"]) for fam, _ in chosen]
    fa = tmp_path / "copies.faa"
    _write_fasta(fa, recs)
    ga_tbl, all_tbl = tmp_path / "ga.tbl", tmp_path / "all.tbl"
    subprocess.run(["hmmsearch", "--cut_ga", "--cpu", "1", "--noali", "--tblout", str(ga_tbl), str(q), str(fa)],
                   check=True, capture_output=True)
    subprocess.run(["hmmsearch", "--max", "-T", "-1000", "--domT", "-1000", "--cpu", "1", "--noali",
                    "--tblout", str(all_tbl), str(q), str(fa)], check=True, capture_output=True)
    hits, scores = _tblout(ga_tbl), _tblout(all_tbl)
    near, far_lower, n_fam = [], 0, 0
    bins = collections.defaultdict(lambda: [0, 0])
    roots_found = 0
    per_copy = []                          # (divergence, agrees with root on GA, |score - root score|)
    for fam, name in chosen:
        root_found = (f"F{fam:06d}", name) in hits
        root_score = scores.get((f"F{fam:06d}", name), -1000.0)
        roots_found += root_found
        cs = sorted(by_fam[fam], key=lambda g: float(g["subst_per_site"]))
        for g in cs:
            d = float(g["subst_per_site"])
            found = (g["copy_label"], name) in hits
            per_copy.append((d, found == root_found,
                             abs(scores.get((g["copy_label"], name), -1000.0) - root_score)))
            if d <= 0.25:
                near.append(found == root_found)
            b = min(int(d / 0.5), 6)
            bins[b][0] += found
            bins[b][1] += 1
        s = [scores.get((g["copy_label"], name), -1000.0) for g in cs]
        half = len(s) // 2
        if float(cs[half - 1]["subst_per_site"]) < float(cs[-half]["subst_per_site"]):
            n_fam += 1
            far_lower += np.mean(s[-half:]) < np.mean(s[:half])
    curve = {f"{b * 0.5:.1f}-{(b + 1) * 0.5:.1f}" if b < 6 else ">=3.0": f"{k}/{n}" for b, (k, n) in sorted(bins.items())}
    per_copy.sort()
    k = len(per_copy) // 3
    first, last = per_copy[:k], per_copy[-k:]
    agree = (np.mean([a for _, a, _ in first]), np.mean([a for _, a, _ in last]))
    delta = (np.median([s for _, _, s in first]), np.median([s for _, _, s in last]))
    print(f"\nhmmsearch --cut_ga detection by substitutions per site from the root: {curve}; roots found "
          f"{roots_found}/{len(chosen)}; copies within 0.25 agreeing with their root {sum(near)}/{len(near)}; "
          f"agreement with the root nearest/farthest third {agree[0]:.3f}/{agree[1]:.3f}, median |score - "
          f"root score| {delta[0]:.1f}/{delta[1]:.1f} bits; families whose far half scores lower: "
          f"{far_lower}/{n_fam}")
    assert agree[0] > agree[1]
    assert delta[0] < delta[1]
    assert far_lower <= 0.5 * n_fam + 3 * math.sqrt(0.25 * n_fam)


@needs_pfam
@pytest.mark.covers("F-02")
def test_root_sampler_scores_like_hmmemit(pfam_world, fetched, tmp_path):
    """§2.3 cross-check (R1 in the spec; here because it needs the binaries and the Pfam file): for 20
    scan-panel HMMs, 500 draws of the numpy root sampler (`profiles.sample_domain`, the t -> infinity
    distribution of every site) and `hmmemit -N 500 --seed <derived, non-zero>` give hmmsearch full-sequence
    scores whose two-sample KS p-value exceeds 0.01 for every HMM."""
    sw, out = pfam_world
    hmm_all, d = fetched
    model = sequences.model(out, pfam_hmm=PFAM)
    scan = [p for p in model.panel if p.in_scan_panel and not p.promiscuous][:20]
    assert len(scan) == 20
    report = []
    for prof in scan:
        hmm = _one(hmm_all, prof.name, tmp_path)
        rng = seeds.generator(sw.seed, f"test/crosscheck/{prof.name}")
        ours = []
        for i in range(500):
            aa, _ = profiles.sample_domain(prof, rng)
            ours.append((f"ours{i}", "".join(codons.AA[x] for x in aa)))
        theirs_fa = tmp_path / f"{prof.name}.emit.faa"
        seed = seeds.tool_seed(sw.seed, f"test/crosscheck/hmmemit/{prof.name}")
        with open(theirs_fa, "wb") as fh:
            subprocess.run(["hmmemit", "-N", "500", "--seed", str(seed), str(hmm)], check=True, stdout=fh)
        theirs = [(f"emit{i}", s) for i, s in enumerate(_fasta(theirs_fa).values())]
        both = tmp_path / f"{prof.name}.both.faa"
        _write_fasta(both, ours + theirs)
        tbl = tmp_path / f"{prof.name}.tbl"
        subprocess.run(["hmmsearch", "--max", "-T", "-1000", "--domT", "-1000", "--cpu", "1", "--noali",
                        "--tblout", str(tbl), str(hmm), str(both)], check=True, capture_output=True)
        sc = _tblout(tbl)
        a = [sc.get((k, prof.name), -1000.0) for k, _ in ours]
        b = [sc.get((k, prof.name), -1000.0) for k, _ in theirs]
        p = _ks_p(a, b)
        report.append((prof.name, round(float(np.median(a)), 1), round(float(np.median(b)), 1), round(p, 4)))
    print(f"\nroot sampler vs hmmemit (name, median ours, median hmmemit, KS p): {report}")
    assert all(p > 0.01 for *_, p in report), report

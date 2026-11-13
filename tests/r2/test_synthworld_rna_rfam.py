"""Ring R2: tool parity of synthworld's Rfam-emitted RNA motifs against the Infernal binaries (synthworld_spec
§2.5, §7.2 cmsearch; card F-02 session 4).

Needs Infernal (`cmfetch`, `cmemit`, `cmsearch`, in the bio-box image), ZOMBI2, and Rfam 15.1's
`Rfam.cm.gz` from ``$BIOBOX_RFAM_CM``, else ``$BIOBOX_GOLDEN_DIR/rfam/Rfam.cm.gz`` (on the pod
``scripts/pod_pytest.sh --golden`` sets the latter, read-only under /tmp/if2453); skips otherwise, so it skips in
GitHub Actions. The micro world is built with ``rna_motifs.source: rfam``. The models the oracle uses are
fetched from the library by `cmfetch -f` (Infernal's own reader), not by the generator's `rna.extract_rfam`.

What is checked: the world's clade emissions are exactly what `cmemit -N <clades> -a --dna --seed <the
derived seed>` gives on the fetched model (the seed plumbing of §6.1), and `cmsearch --cut_ga` with the
linked model finds each planted instance at its truth locus (same strand, at least half of it covered), for
nearly every instance at low divergence from its emission; recall by divergence bin is printed (§7.2
"recall of truth_rna by difficulty bin").
"""

from __future__ import annotations

import collections
import csv
import gzip
import importlib.util
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from biobox.synthworld import backbone, intergenic, rna, seeds, sequences

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]


def _rfam() -> str:
    for p in (os.environ.get("BIOBOX_RFAM_CM"),
              os.path.join(os.environ.get("BIOBOX_GOLDEN_DIR", ""), "rfam", "Rfam.cm.gz")):
        if p and Path(p).is_file():
            return p
    return ""


RFAM = _rfam()
_missing = [b for b in ("cmfetch", "cmemit", "cmsearch") if shutil.which(b) is None]
needs_rfam = pytest.mark.skipif(
    bool(_missing) or not RFAM or importlib.util.find_spec("zombi2") is None,
    reason=f"needs Infernal (missing: {_missing or 'none'}), Rfam.cm.gz in BIOBOX_RFAM_CM or "
           f"BIOBOX_GOLDEN_DIR/rfam ({RFAM or 'not found'}) and zombi2; runs on the pod in the bio-box image")

_RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s: str) -> str:
    return "".join(_RC[c] for c in reversed(s))


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


def _tblout(path: Path) -> list[dict]:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        out.append({"target": f[0], "query": f[2], "from": int(f[7]), "to": int(f[8]), "strand": f[9],
                    "score": float(f[14])})
    return out


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if _missing or not RFAM or importlib.util.find_spec("zombi2") is None:
        pytest.skip("needs Infernal, Rfam.cm.gz and zombi2")
    base = tmp_path_factory.mktemp("micro_rfam")
    out = base / "world"
    params = backbone.load_preset("micro", {"rna_motifs": {"source": "rfam"}})
    backbone.build(params, params["default_seed"], out)
    sequences.build(out)
    rw = intergenic.build(out, rfam_cm=RFAM)
    # the oracle's models: Infernal's cmfetch on the decompressed library, one CM file per model
    plain = base / "Rfam.cm"
    with gzip.open(RFAM, "rb") as src, open(plain, "wb") as dst:
        shutil.copyfileobj(src, dst, 1 << 20)
    cms = {}
    for m in rw.models:
        keys = base / f"{m.accession}.key"
        keys.write_text(m.accession + "\n")
        cm = base / f"{m.accession}.cm"
        subprocess.run(["cmfetch", "-f", "-o", str(cm), str(plain), str(keys)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        cms[m.accession] = cm
    plain.unlink()
    return rw, out, cms, base


@needs_rfam
def test_rfam_models_are_the_presets_and_carry_their_ga(world):
    rw, out, cms, _ = world
    want = backbone.load_preset("micro")["rna_motifs"]["rfam_models"]
    assert [m.accession for m in rw.models] == want
    assert all(m.source == "rfam" and m.ga > 0 for m in rw.models)
    assert rw.meta["motif_source"]["release"] == "15.1" and rw.meta["versions"]["infernal"].startswith("1.1")
    # the generator's verbatim copy and cmfetch's agree on each model's header
    for m in rw.models:
        ga = [line.split()[1] for line in cms[m.accession].read_text().splitlines() if line.startswith("GA ")]
        assert float(ga[0]) == m.ga


@needs_rfam
def test_clade_emissions_are_cmemit_at_the_derived_seed(world):
    rw, out, cms, _ = world
    em = collections.defaultdict(list)
    for r in _tsv(out / "replicons/emissions.tsv"):
        em[r["link"]].append(r)
    p = backbone.load_preset("micro")
    for lk in rw.links:
        rows = sorted(em[lk.link], key=lambda r: int(r["emission_id"].split(".")[1]))
        assert len(rows) == lk.n_clades > 0
        argv = ["cmemit", "-N", str(lk.n_clades), "-a", "--dna", "--seed",
                str(seeds.tool_seed(p["default_seed"], f"S5/emit/{lk.link}"))]
        if p["rna_motifs"]["exp"] != 1.0:
            argv += ["--exp", str(p["rna_motifs"]["exp"])]
        done = subprocess.run(argv + [str(cms[lk.model_acc])], check=True, capture_output=True, text=True)
        per: dict[str, str] = {}        # Stockholm blocks of a long alignment, joined per name
        for line in done.stdout.splitlines():
            if line.strip() and not line.startswith("#") and not line.startswith("//"):
                name, aln = line.split()
                per[name] = per.get(name, "") + aln
        seqs = [s.replace("-", "").replace(".", "").upper() for s in per.values()]
        assert seqs == [r["sequence"] for r in rows], lk.link


@needs_rfam
@pytest.mark.covers("F-02", "B-SIM")
def test_cmsearch_finds_planted_motifs_at_their_truth_loci(world, record_property):
    rw, out, cms, base = world
    seqs = {}
    for f in sorted((out / "replicons/fasta").glob("*.fna")):
        seqs.update(_fasta(f))
    rec = {(r.genome_id, r.replicon): r.record for r in rw.replicons}
    rnas = _tsv(out / "replicons/rna.tsv")
    ga = {m.accession: m.ga for m in rw.models}
    flank = 100
    by_model = collections.defaultdict(list)
    for r in rnas:
        by_model[r["model_acc"]].append(r)
    found, score = {}, {}
    for acc, rows in sorted(by_model.items()):
        fa = base / f"{acc}.targets.fna"
        with open(fa, "w") as fh:
            for r in rows:
                seq = seqs[rec[(r["genome_id"], int(r["replicon"]))]]
                N = len(seq)
                s = int(r["start"]) - flank
                ring = seq * 3
                cut = ring[N + s:N + int(r["end"]) + flank]
                fh.write(f">{r['rna_id']}\n{cut}\n")
        tbl = base / f"{acc}.tbl"
        subprocess.run(["cmsearch", "--cpu", "1", "--cut_ga", "--notrunc", "--noali", "--tblout", str(tbl),
                        "-o", os.devnull, str(cms[acc]), str(fa)], check=True, stdout=subprocess.DEVNULL)
        hits = collections.defaultdict(list)
        for h in _tblout(tbl):
            hits[h["target"]].append(h)
        for r in rows:
            L = int(r["length"])
            truth = (flank, flank + L)           # the instance's span in its target, on +
            want = "+" if int(r["strand"]) == 1 else "-"
            ok, best = False, float("nan")
            for h in hits[r["rna_id"]]:
                a, b = min(h["from"], h["to"]) - 1, max(h["from"], h["to"])
                cover = max(0, min(b, truth[1]) - max(a, truth[0]))
                if h["strand"] == want and cover >= 0.5 * L and h["score"] >= ga[acc]:
                    ok = True
                    best = h["score"] if np.isnan(best) else max(best, h["score"])
            found[r["rna_id"]], score[r["rna_id"]] = ok, best
    bins = [(0.0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 1.01)]
    table = []
    for lo, hi in bins:
        ids_ = [r["rna_id"] for r in rnas
                if (lo < int(r["n_subst_from_emission"]) / int(r["length"]) <= hi) or (lo == 0.0 and r["n_subst_from_emission"] == "0")]
        n = len(ids_)
        table.append((f"{lo}-{hi}", n, sum(found[i] for i in ids_)))
    per_model = {acc: (len(rows), sum(found[r["rna_id"]] for r in rows)) for acc, rows in sorted(by_model.items())}
    print(f"\ncmsearch --cut_ga recall by divergence from the emission (bin, instances, found): {table}; "
          f"by model (instances, found): {per_model}")
    record_property("rfam_recall_bins", str(table))
    low = [r["rna_id"] for r in rnas if int(r["n_subst_from_emission"]) / int(r["length"]) <= 0.05]
    assert len(low) >= 50
    assert sum(found[i] for i in low) >= 0.9 * len(low)

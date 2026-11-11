"""Card 0-09c, item 4: the E-02 spike-in panel under the scan variants (run ON THE POD, in the image).

    python infernal_spikein.py emit  <repo> <plain Rfam.cm> <work dir>
    python infernal_spikein.py score <repo> <work dir> <hits.jsonl>

`emit` builds the panel's planted contig with the test's own helpers, imported from
`tests/r3/test_g04_spikeins.py` (`PANEL`, `_emit`, `_plant`; same seeds, so the same contig and truth as the
test), fetches the panel's 20 models with `cmfetch` as the test does, and writes, in the work dir:
`panel.cm`, `truth.json`, and the contig as one intergenic-style region (`region.fna`, `region.map.json`:
the whole linear contig is the region, so a scan of it isolates the cmsearch settings). lane.sh then runs
cmsearch on it with the producer's argv (default filters) and with `--rfam`, and `infernal_igr.py collect`
maps the hits back.

`score` applies the test's rule (a hit of the planted family, on its strand, overlapping at least
`MIN_OVERLAP` of the planted length) and prints the recovery with the test's exact interval.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from infernal_igr import g04_module, read_hits  # noqa: E402


def emit(repo: Path, full_cm: Path, work: Path) -> None:
    g = g04_module(repo)
    work.mkdir(parents=True, exist_ok=True)
    keys = work / "panel.keys"
    keys.write_text("".join(a + "\n" for a in g.PANEL))
    panel_cm = work / "panel.cm"
    subprocess.run([g.CMFETCH, "-f", "-o", str(panel_cm), str(full_cm), str(keys)], check=True,
                   stdin=subprocess.DEVNULL, capture_output=True)
    got = sorted({ln.split()[1] for ln in panel_cm.read_text().splitlines() if ln.startswith("ACC ")})
    assert got == sorted(g.PANEL), "cmfetch did not return exactly the panel"
    draws_dir = work / "draws"
    draws_dir.mkdir(exist_ok=True)
    emitted, draws = g._emit(panel_cm, draws_dir)
    contig, truth = g._plant(emitted)
    (work / "truth.json").write_text(json.dumps({"truth": truth, "draws": draws, "contig_len": len(contig)},
                                                sort_keys=True) + "\n")
    name = "g04_0"
    (work / "region.fna").write_bytes(b">" + name.encode() + b"\n" + contig + b"\n")
    doc = {"mode": "igr", "batch": "g04", "margin": 0, "flank": 0, "records": 1, "bases_searched": len(contig),
           "genome_bases": len(contig), "genomes": ["g04_synthetic"],
           "map": {name: {"contig": "g04", "genome": "g04_synthetic", "len": len(contig), "circular": False,
                          "start": 0, "end": len(contig)}}}
    (work / "region.map.json").write_text(json.dumps(doc, sort_keys=True) + "\n")
    print(json.dumps({"planted": len(truth), "contig_len": len(contig), "draws": draws}, sort_keys=True))


def score(repo: Path, work: Path, hits_path: Path) -> None:
    g = g04_module(repo)
    t = json.loads((work / "truth.json").read_text())
    truth, hits = t["truth"], read_hits(hits_path)
    per, missed = {a: [] for a in g.PANEL}, []
    for p in truth:
        need = g.MIN_OVERLAP * (p["end"] - p["start"])
        ok = any(h["acc"] == p["family"] and h["strand"] == p["strand"]
                 and g._overlap(h["start"], h["end"], p["start"], p["end"]) >= need for h in hits)
        per[p["family"]].append(ok)
        if not ok:
            missed.append({k: p[k] for k in ("family", "start", "end", "strand", "alone")})
    k, n = sum(map(sum, per.values())), len(truth)
    lo, hi = g.clopper_pearson(k, n)
    extra = [h for h in hits if not any(g._overlap(h["start"], h["end"], p["start"], p["end"]) for p in truth)]
    print(json.dumps({"recovered": k, "planted": n, "rate": round(k / n, 4), "ci95": [round(lo, 3), round(hi, 3)],
                      "hits": len(hits), "on_background": len(extra), "missed": missed,
                      "per_family": {a: f"{sum(v)}/{len(v)}" for a, v in per.items()}}, sort_keys=True))


if __name__ == "__main__":
    cmd, a = sys.argv[1], sys.argv[2:]
    if cmd == "emit":
        emit(Path(a[0]), Path(a[1]), Path(a[2]))
    elif cmd == "score":
        score(Path(a[0]), Path(a[1]), Path(a[2]))
    else:
        print(__doc__, file=sys.stderr)
        sys.exit(2)

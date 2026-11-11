"""Ring R3, test G-04 (card E-02): Rfam spike-ins recovered by the ncRNA producer (card E-01).

The truth is constructed (AGENTS §5.6, §5.2): 20 Rfam 15.1 families, fixed below by accession; for each,
`cmemit -N 5 --seed <SEED + index>` draws 5 sequences from the family's own model (Infernal's emitter, not
the producer); the 100 sequences are planted, in a seeded random order and on a seeded random strand, into
one linear synthetic contig of seeded uniform-random background (A, C, G, T equally likely), at recorded
coordinates (0-based half-open, strand +1 or -1; a minus-strand sequence is planted as its reverse
complement). The producer (`biobox.l1.ncrna.run`) then searches that contig.

Panel rule (written before any result was seen): the families are the accessions of `PANEL`, chosen by
category to cover tRNA (RF00005, RF01852), rRNA fragments that are short (RF00001, RF00002), seven
riboswitches (RF00162, RF00167, RF00059, RF00050, RF00174, RF00504, RF00168), bacterial RNase P (RF00010,
RF00011), tmRNA (RF00023), 6S (RF00013), group I and group II introns (RF00028, RF00029), bacterial SRP
(RF00169) and two small structured RNAs (RF00230 T-box, RF00080 yybP-ykoY). Large rRNA models (SSU, LSU) are
left out for cost (E-01 measured 1,319 bases per CPU-second for the full library). The list is not
changed after results.

A planted sequence is recovered when the producer reports a hit of its own family (`model_acc`) on the same
strand whose overlap with the planted interval is at least half of the planted length. Anything else
(another family's hit on it, a hit on the wrong strand, a shorter overlap) is a miss; hits that overlap no
planted sequence are counted as extra and reported, not scored.

Search space: the producer runs `cmsearch -Z 1e5 --cut_ga` (DEC-04, `ncrna.Z_NT_MB`) whatever the library;
so searching the 20 chosen models alone (fetched verbatim from Rfam 15.1 with `cmfetch`) gives each of their
hits the E-value and the GA decision the full library gives it (the GA cutoff is a bit score, and Z is
fixed, not the library size or the contig length).

Where the models are: the full `Rfam.cm.gz` of Rfam 15.1, whose sha256 must be refdb.yaml's. The file is
read from `$BIOBOX_RFAM_CM`, else `$BIOBOX_GOLDEN_DIR/rfam/Rfam.cm.gz`; the test skips when neither exists,
when Infernal (`cmemit`, `cmfetch`, `cmsearch`) is missing or when BIOBOX_GOLDEN_DIR is unset. Nothing is
downloaded.

Recovery rate is reported overall and per family with the exact (Clopper-Pearson) two-sided 95% binomial
interval, computed here by bisection on the binomial tail (no new dependency; scipy is not one), and the
interval code is itself tested against two published values.
"""

from __future__ import annotations

import gzip
import hashlib
import math
import os
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from biobox import fasta, ids, io, provenance
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input

pytestmark = [pytest.mark.r3, pytest.mark.covers("G-04")]

PANEL = {  # accession: name (Rfam 15.1)
    "RF00005": "tRNA", "RF01852": "tRNA-Sec", "RF00001": "5S_rRNA", "RF00002": "5_8S_rRNA",
    "RF00162": "SAM", "RF00167": "Purine", "RF00059": "TPP", "RF00050": "FMN", "RF00174": "Cobalamin",
    "RF00504": "Glycine", "RF00168": "Lysine", "RF00010": "RNaseP_bact_a", "RF00011": "RNaseP_bact_b",
    "RF00023": "tmRNA", "RF00013": "6S", "RF00028": "Intron_gpI", "RF00029": "Intron_gpII",
    "RF00169": "Bacteria_small_SRP", "RF00230": "T-box", "RF00080": "yybP-ykoY",
}
PER_FAMILY = 5
SEED = 20261110  # background and layout; cmemit uses SEED + the family's index in PANEL
MIN_OVERLAP = 0.5  # of the planted length
SPACER = (200, 800)  # background bases between planted sequences, inclusive range
EDGE = (200, 800)  # background before the first and after the last
# Lower bound on overall recovery, set from the first measured run (pod wilkcomp02, 2026-11-10, Infernal 1.1.5,
# these seeds and this panel): 67 of 100 recovered (exact 95% interval 0.569 to 0.761), 67 hits and none on
# background, so every miss is a planted sequence that scored below its model's GA cutoff (cmemit draws from
# the model, so some draws are weak). Per family that run gave 5 of 5 for 8 families and 0 of 5 for tRNA-Sec,
# FMN and T-box; the panel was not changed after seeing this. The bound is 0.55, below the interval's lower
# end, so a change in Infernal or in the producer that costs 13 or more recoveries fails the test.
MEASURED_OVERALL = 0.67  # measured
MIN_OVERALL = 0.55

CMSEARCH, CMEMIT, CMFETCH = (shutil.which(t) for t in ("cmsearch", "cmemit", "cmfetch"))
RFAM_SHA256 = next((db.get("sha256") for db in provenance.refdb_libraries() if db["name"] == "rfam"), None)


# --- the exact binomial interval (no dependency) ---------------------------------------------------------


def _cdf(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p)."""
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k + 1))


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact two-sided (1 - alpha) interval for a binomial proportion, by bisection on the binomial tail:
    lower p with P(X >= k | p) = alpha/2, upper p with P(X <= k | p) = alpha/2."""
    if not 0 <= k <= n or n < 1:
        raise ValueError(f"k={k}, n={n}")

    def solve(f, target):  # f is decreasing in p; find p with f(p) = target
        lo, hi = 0.0, 1.0
        for _ in range(100):
            mid = (lo + hi) / 2
            if f(mid) > target:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    # P(X >= k | p) increases with p: solve on q = 1 - p instead
    lower = 0.0 if k == 0 else 1 - solve(lambda q: 1 - _cdf(k - 1, n, 1 - q), alpha / 2)
    upper = 1.0 if k == n else solve(lambda p: _cdf(k, n, p), alpha / 2)
    return lower, upper


def test_clopper_pearson_matches_published_values():
    """Published exact 95% limits: 5 of 10 gives (0.1871, 0.8129) and 0 of 10 gives (0, 0.3085) (Clopper and
    Pearson 1934 as tabulated in Brown, Cai and DasGupta 2001 and in R's `binom.test`)."""
    lo, hi = clopper_pearson(5, 10)
    assert (round(lo, 4), round(hi, 4)) == (0.1871, 0.8129)
    lo, hi = clopper_pearson(0, 10)
    assert lo == 0.0 and round(hi, 4) == 0.3085
    lo, hi = clopper_pearson(10, 10)
    assert hi == 1.0 and round(lo, 4) == 0.6915


# --- the panel -------------------------------------------------------------------------------------------


def _rfam_path() -> Path | None:
    for p in (os.environ.get("BIOBOX_RFAM_CM"),
              os.path.join(os.environ["BIOBOX_GOLDEN_DIR"], "rfam", "Rfam.cm.gz")
              if os.environ.get("BIOBOX_GOLDEN_DIR") else None):
        if p and Path(p).is_file():
            return Path(p)
    return None


def _run(argv: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(argv, check=True, stdin=subprocess.DEVNULL, capture_output=True, **kw)


def _revcomp(s: bytes) -> bytes:
    return s.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]


def _emit(plain_panel: Path, work: Path) -> dict[str, list[bytes]]:
    """5 sequences per family from `cmemit` (RNA letters turned to DNA), a family's model fetched on its own."""
    out = {}
    for i, acc in enumerate(PANEL):
        keys, one, fa = work / f"{acc}.keys", work / f"{acc}.cm", work / f"{acc}.fa"
        keys.write_text(acc + "\n")
        _run([CMFETCH, "-f", "-o", str(one), str(plain_panel), str(keys)])
        _run([CMEMIT, "-N", str(PER_FAMILY), "--seed", str(SEED + i), "-o", str(fa), str(one)])
        seqs = [r.seq.upper().replace(b"U", b"T") for r in fasta.read(fa)]
        assert len(seqs) == PER_FAMILY and all(set(s) <= set(b"ACGT") for s in seqs), acc
        out[acc] = seqs
    return out


def _plant(emitted: dict[str, list[bytes]]) -> tuple[bytes, list[dict]]:
    """The contig and the truth: one dict per planted sequence (family, start, end, strand), 0-based half-open."""
    rng = random.Random(SEED)

    def background(n: int) -> bytes:
        return bytes(rng.choices(b"ACGT", k=n))

    items = [(acc, s) for acc, seqs in emitted.items() for s in seqs]
    rng.shuffle(items)
    parts, truth, pos = [background(rng.randint(*EDGE))], [], 0
    pos = len(parts[0])
    for acc, s in items:
        strand = rng.choice((1, -1))
        parts.append(s if strand == 1 else _revcomp(s))
        truth.append({"family": acc, "start": pos, "end": pos + len(s), "strand": strand})
        pos += len(s)
        parts.append(background(rng.randint(*SPACER)))
        pos += len(parts[-1])
    parts.append(background(rng.randint(*EDGE)))
    return b"".join(parts), truth


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


@pytest.fixture(scope="module")
def panel_run(tmp_path_factory):
    if CMSEARCH is None or CMEMIT is None or CMFETCH is None:
        pytest.skip("Infernal (cmsearch, cmemit, cmfetch) not on PATH; it is in the test image")
    rfam = _rfam_path()
    if rfam is None:
        pytest.skip("Rfam 15.1 not found (set BIOBOX_RFAM_CM, or put Rfam.cm.gz in $BIOBOX_GOLDEN_DIR/rfam/)")
    digest = hashlib.sha256(rfam.read_bytes()).hexdigest()
    assert digest == RFAM_SHA256, f"{rfam} is not the Rfam 15.1 file refdb.yaml pins"
    from biobox.l1 import ncrna

    base = tmp_path_factory.mktemp("g04")
    work, scratch, root = base / "work", base / "scratch", base / "root"
    work.mkdir()
    scratch.mkdir()
    full = work / "Rfam.cm"
    with gzip.open(rfam, "rb") as src, open(full, "wb") as dst:
        shutil.copyfileobj(src, dst, 1 << 20)
    keys = work / "panel.keys"
    keys.write_text("".join(a + "\n" for a in PANEL))
    panel_cm = work / "panel.cm"
    _run([CMFETCH, "-f", "-o", str(panel_cm), str(full), str(keys)])
    full.unlink()
    got = [ln.split()[1] for ln in panel_cm.read_text().splitlines() if ln.startswith("ACC ")]
    assert sorted(set(got)) == sorted(PANEL), "cmfetch did not return exactly the panel"
    # a panel library of the plain file for the producer, gzipped like Rfam.cm.gz
    lib = base / "panel.cm.gz"
    lib.write_bytes(gzip.compress(panel_cm.read_bytes(), 9, mtime=0))

    emitted = _emit(panel_cm, work)
    contig, truth = _plant(emitted)
    fna = base / "synthetic.fna"
    fna.write_bytes(b">g04_synthetic\n" + contig + b"\n")
    ingest([SourceObject(fna, "g04", "1", "CC0-1.0", True)], root, "S0")
    cu = ids.contig_uid(contig)
    write_genome_input(root, "S0", [GenomeInput((cu,), 11, kind="bacteria", circular=())])
    report = ncrna.run(root, "S0", lib, cpus=2, tmp=scratch)
    import pyarrow.parquet as pq

    rel = f"l1/ncrna/producer={report['producer']}/stage=S0"
    hits = [r for f in io.part_files(root / rel) for r in pq.read_table(f, partitioning=None).to_pylist()]
    return {"contig_len": len(contig), "truth": truth, "hits": hits, "report": report}


def test_g04_rfam_spike_ins_are_recovered(panel_run):
    truth, hits = panel_run["truth"], panel_run["hits"]
    assert len(truth) == PER_FAMILY * len(PANEL)
    assert all(0 <= t["start"] < t["end"] <= panel_run["contig_len"] for t in truth)
    per: dict[str, list[bool]] = {a: [] for a in PANEL}
    for t in truth:
        need = MIN_OVERLAP * (t["end"] - t["start"])
        ok = any(h["model_acc"] == t["family"] and h["strand"] == t["strand"]
                 and _overlap(h["start"], h["end"], t["start"], t["end"]) >= need for h in hits)
        per[t["family"]].append(ok)
    extra = [h for h in hits if not any(_overlap(h["start"], h["end"], t["start"], t["end"]) for t in truth)]
    k, n = sum(map(sum, per.values())), len(truth)
    lo, hi = clopper_pearson(k, n)
    lines = [f"G-04 overall recovery {k}/{n} = {k / n:.3f}, exact 95% interval ({lo:.3f}, {hi:.3f}); "
             f"{len(hits)} hits, {len(extra)} on no planted sequence; contig {panel_run['contig_len']} bases"]
    for acc, flags in per.items():
        a, b = clopper_pearson(sum(flags), len(flags))
        lines.append(f"  {acc:8} {PANEL[acc]:20} {sum(flags)}/{len(flags)}  ({a:.3f}, {b:.3f})")
    print("\n" + "\n".join(lines))
    assert k / n >= MIN_OVERALL, "\n".join(lines)

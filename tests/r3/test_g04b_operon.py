"""Ring R3, test G-04b (card E-08, DEC-28): clan competition on planted loci where models overlap.

The E-02 panel (`test_g04_spikeins`) plants `cmemit` draws, which one model each finds; here the overlaps are the
point, so the planted loci are real ones. Windows of three golden genomes (sha256 checked against golden.lock.yaml)
are cut at the coordinates below and joined, in this order and on the + strand, by seeded uniform-random spacers
into one linear synthetic contig. The window coordinates (0-based half-open on the genome) were chosen from the hits
of measurement 0-09c's margin-200 `--notrunc` scan on the 8 golden genomes, which card E-07's producer reproduces
exactly (planning/measurements/2026-11-11_ncrna_producer.md), widened by a few hundred bases; the truth is not taken
from there, it is the binaries' own output on the planted contig, below.

    operon_ecoli    E. coli rrnB: 16S, the spacer, 23S, 5S. Two SSU models (bacteria, archaea; clan CL00111) hit
                    the 16S gene and two LSU models (CL00112) the 23S gene.
    operon_bsub     B. subtilis rrnO: the same, and Bacillaceae-1 (no clan) antisense to the 5' end of the 16S gene.
    trna_sec        M. genitalium tRNA cluster: tRNA-Sec (CL00001) on a tRNA of the same clan, and four more tRNAs.
    glm             E. coli GlmZ_SraJ and GlmY_tke1 (both CL00011) on one locus and strand.
    srp             B. subtilis Bacteria_large_SRP and Bacteria_small_SRP (both CL00003) on one locus and strand.
    sdsr_srac       E. coli SraC_RyeA and SdsR_RyeB (both CL00105) on opposite strands of one locus.
    trna_hpnc       E. coli tRNA (CL00001) and HPnc0260 (no clan) on opposite strands of one locus.
    trna_rtt        E. coli tRNA (CL00001) and RtT (no clan) on one locus and strand.

The model set is these 16 models, every one intergenic; the stage's gene calls are none, so the contig is one region
searched with `--notrunc` (as `test_g04_spikeins` does). Oracles (AGENTS §5.2), both Infernal 1.1.5 binaries:

- `cmsearch --cut_ga --notrunc -Z 1e5` on the contig: every hit it reports is one `l1_ncrna` row (all hits stay);
- `cmscan --cut_ga --notrunc -Z 1e5 --fmt 2 --clanin <Rfam.clanin cut to the panel's models> --oclan` on the contig
  with the same models (cmscan refuses a clan file that names a model its database lacks): its `olp` column marks a
  hit `=` when it overlaps a better-scoring hit of its clan (`--oclan`: only within clans). A producer hit is a clan
  loser exactly when cmscan marks it `=`.

Expectations, written before any result was seen: per operon exactly one kept SSU hit and one kept LSU hit on the
rRNA genes, the others of those clans present and flagged; 5S and the tRNAs kept; at the tRNA cluster one of tRNA and
tRNA-Sec kept on the shared locus; one of GlmZ/GlmY and one of the two SRPs kept; both antisense partners kept (one
clan on opposite strands, and a clan against no clan); tRNA and RtT both kept; no two kept hits of one clan overlap on
a strand; the features are the distinct loci of the kept hits.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import random
import shutil
import sys
from pathlib import Path

import pytest
import yaml

from biobox import fasta, ids, io
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parents[1] / "r2"))
import test_g04_spikeins as g04  # noqa: E402
import test_ncrna_clans as cl  # noqa: E402

pytestmark = [pytest.mark.r3, pytest.mark.covers("G-04", "E-08")]

REPO = Path(__file__).resolve().parents[2]
CMSCAN, CMPRESS = shutil.which("cmscan"), shutil.which("cmpress")
PANEL = {  # accession: name (Rfam 15.1); all in the set prok-1
    "RF00177": "SSU_rRNA_bacteria", "RF01959": "SSU_rRNA_archaea", "RF02541": "LSU_rRNA_bacteria",
    "RF02540": "LSU_rRNA_archaea", "RF00001": "5S_rRNA", "RF00005": "tRNA", "RF01852": "tRNA-Sec",
    "RF00391": "RtT", "RF00083": "GlmZ_SraJ", "RF00128": "GlmY_tke1", "RF00101": "SraC_RyeA",
    "RF00111": "SdsR_RyeB", "RF02194": "HPnc0260", "RF01690": "Bacillaceae-1", "RF01854": "Bacteria_large_SRP",
    "RF00169": "Bacteria_small_SRP",
}
WINDOWS = [  # name, golden genome, start, end (0-based half-open on the genome, + strand)
    ("operon_ecoli", "ecoli_mg1655", 4166300, 4172100),
    ("operon_bsub", "bsub_168", 9500, 15200),
    ("trna_sec", "mgenitalium", 257250, 257900),
    ("glm", "ecoli_mg1655", 3986300, 3986750),
    ("srp", "bsub_168", 26300, 26800),
    ("sdsr_srac", "ecoli_mg1655", 1922950, 1923400),
    ("trna_hpnc", "ecoli_mg1655", 1746400, 1746850),
    ("trna_rtt", "ecoli_mg1655", 780450, 780800),
]
SEED = 20261111
SPACER = 300  # random bases between windows and at both ends
Z_NT_MB = g04.Z_NT_MB


def _golden(name: str) -> bytes:
    """The golden genome's one sequence, its gzip file checked against golden.lock.yaml."""
    entry = yaml.safe_load((REPO / "golden.lock.yaml").read_text())["files"][name]
    path = Path(os.environ["BIOBOX_GOLDEN_DIR"]) / entry["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256_gz"], f"{path} is not golden.lock's"
    recs = fasta.read(path)
    assert len(recs) == 1
    return recs[0].seq.upper()


def _parse_fmt2(path: Path) -> list[dict]:
    """`cmscan --fmt 2 --tblout` rows, parsed here: idx, model name and accession, query, clan, mdl, model and
    sequence coordinates, strand, trunc, pass, gc, bias, score, E-value, inc, olp, ..."""
    out = []
    for ln in path.read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        lo, hi = sorted((int(f[9]), int(f[10])))
        out.append({"name": f[1], "acc": f[2], "query": f[3], "clan": None if f[5] == "-" else f[5],
                    "start": lo - 1, "end": hi, "strand": 1 if f[11] == "+" else -1, "score": float(f[16]),
                    "olp": f[19]})
    return out


def _overlap(a: dict, b: dict) -> int:
    return g04._overlap(a["start"], a["end"], b["start"], b["end"])


@pytest.fixture(scope="module")
def operon_run(tmp_path_factory):
    if g04.CMSEARCH is None or g04.CMFETCH is None or CMSCAN is None or CMPRESS is None:
        pytest.skip("Infernal (cmsearch, cmfetch, cmscan, cmpress) not on PATH; it is in the test image")
    rfam = g04._rfam_path()
    if rfam is None or not os.environ.get("BIOBOX_GOLDEN_DIR"):
        pytest.skip("Rfam 15.1 and the golden genomes not found (set BIOBOX_GOLDEN_DIR)")
    assert hashlib.sha256(rfam.read_bytes()).hexdigest() == g04.RFAM_SHA256
    from biobox.l1 import ncrna
    from biobox.l1.genecall import GENE_CALL_VERSION
    from biobox.registry.schemas import table

    base = tmp_path_factory.mktemp("g04b")
    work, scratch, root = base / "work", base / "scratch", base / "root"
    work.mkdir()
    scratch.mkdir()
    full = work / "Rfam.cm"
    with gzip.open(rfam, "rb") as src, open(full, "wb") as dst:
        shutil.copyfileobj(src, dst, 1 << 20)
    keys = work / "panel.keys"
    keys.write_text("".join(a + "\n" for a in PANEL))
    panel_cm = work / "panel.cm"
    g04._run([g04.CMFETCH, "-f", "-o", str(panel_cm), str(full), str(keys)])
    full.unlink()
    lib = base / "panel.cm.gz"
    lib.write_bytes(gzip.compress(panel_cm.read_bytes(), 9, mtime=0))

    rng = random.Random(SEED)
    genomes = {g: _golden(g) for g in sorted({w[1] for w in WINDOWS})}
    parts, where, pos = [], {}, 0
    for name, g, s, e in WINDOWS:
        parts.append(bytes(rng.choices(b"ACGT", k=SPACER)))
        pos += SPACER
        parts.append(genomes[g][s:e])
        where[name] = (pos, pos + e - s)
        pos += e - s
    parts.append(bytes(rng.choices(b"ACGT", k=SPACER)))
    contig = b"".join(parts)
    fna = base / "operons.fna"
    fna.write_bytes(b">g04b_operons\n" + contig + b"\n")
    ingest([SourceObject(fna, "g04b", "1", "INSDC (public)", True)], root, "S0")
    cu = ids.contig_uid(contig)
    write_genome_input(root, "S0", [GenomeInput((cu,), 11, kind="bacteria", circular=())])
    panel_set = base / "panel.tsv"
    panel_set.write_text("acc\tname\twhole_genome\n" + "".join(f"{a}\t{n}\t0\n" for a, n in PANEL.items()))
    io.write_parts([], table("feature"), root / f"l1/feature/producer={GENE_CALL_VERSION}/stage=S0", 24)
    report = ncrna.run(root, "S0", lib, cpus=2, tmp=scratch, model_set=panel_set)
    assert report["regions"] == 1 and report["region_bases"] == len(contig)
    prod = report["producer"]
    hits = cl_rows(root, f"l1/ncrna/producer={prod}/stage=S0")
    feats = cl_rows(root, f"l1/feature/producer={prod}/stage=S0")
    # the oracles: cmsearch for the hits, cmscan with the clan file for the overlap marks
    searched = g04._tblout(["--cut_ga", "--notrunc"], panel_cm, fna, work / "cmsearch.tbl")
    g04._run([CMPRESS, str(panel_cm)])
    # cmscan refuses a clan file naming a model its database lacks: the clan file's lines cut to the panel's models
    clanin = work / "panel.clanin"
    lines = [ln.split("\t") for ln in cl.CLANIN.read_text().splitlines()]
    clanin.write_text("".join("\t".join([f[0], *mine]) + "\n" for f in lines
                              if (mine := [m for m in f[1:] if m in PANEL.values()])))
    tbl = work / "cmscan.tbl"
    g04._run([CMSCAN, "--cpu", "2", "-Z", f"{Z_NT_MB:g}", "--cut_ga", "--notrunc", "--noali", "--fmt", "2",
              "--clanin", str(clanin), "--oclan", "--tblout", str(tbl), "-o", str(work / "cmscan.out"),
              str(panel_cm), str(fna)])
    scanned = _parse_fmt2(tbl)
    return {"hits": hits, "features": feats, "where": where, "searched": searched, "scanned": scanned,
            "report": report, "contig_len": len(contig)}


def cl_rows(root: Path, rel: str) -> list[dict]:
    import pyarrow.parquet as pq

    return [r for f in io.part_files(root / rel) for r in pq.read_table(f, partitioning=None).to_pylist()]


def _in(name: str, run: dict, rows: list[dict]) -> list[dict]:
    s, e = run["where"][name]
    return [r for r in rows if s <= r["start"] and r["end"] <= e]


def test_every_hit_of_the_binary_is_a_row_and_cmscan_reports_the_same_hits(operon_run):
    hits = operon_run["hits"]
    got = {(r["start"], r["end"], r["strand"], r["model_acc"], round(r["bitscore"], 1)) for r in hits}
    want = {(h["start"], h["end"], h["strand"], h["acc"], h["score"]) for h in operon_run["searched"]}
    assert want and got == want
    scanned = {(h["start"], h["end"], h["strand"], h["acc"], h["score"]) for h in operon_run["scanned"]}
    assert scanned == want
    # the planted loci are what the windows were chosen for: several SSU and LSU models on each operon, etc.
    for op in ("operon_ecoli", "operon_bsub"):
        accs = {r["model_acc"] for r in _in(op, operon_run, hits)}
        assert {"RF00177", "RF01959", "RF02541", "RF02540", "RF00001"} <= accs, (op, accs)
    assert {"RF01690"} <= {r["model_acc"] for r in _in("operon_bsub", operon_run, hits)}
    for name, accs in (("trna_sec", {"RF00005", "RF01852"}), ("glm", {"RF00083", "RF00128"}),
                       ("srp", {"RF01854", "RF00169"}), ("sdsr_srac", {"RF00101", "RF00111"}),
                       ("trna_hpnc", {"RF00005", "RF02194"}), ("trna_rtt", {"RF00005", "RF00391"})):
        assert accs <= {r["model_acc"] for r in _in(name, operon_run, hits)}, name


def test_a_hit_loses_exactly_when_cmscan_marks_it_below_a_better_hit_of_its_clan(operon_run):
    olp = {(h["start"], h["end"], h["strand"], h["acc"]): h for h in operon_run["scanned"]}
    clans = cl.clan_file()
    lost = 0
    for r in operon_run["hits"]:
        h = olp[(r["start"], r["end"], r["strand"], r["model_acc"])]
        assert r["clan"] == h["clan"] == clans.get(r["model_name"])
        assert r["clan_winner"] == (h["olp"] != "="), (r, h)
        lost += not r["clan_winner"]
    assert lost >= 5  # 2 operons x (SSU, LSU), tRNA-Sec, GlmY, the small SRP: at least one loser each


def test_one_kept_ssu_and_lsu_per_operon_and_every_other_planted_gene_kept(operon_run):
    hits, feats = operon_run["hits"], operon_run["features"]
    for op in ("operon_ecoli", "operon_bsub"):
        here = _in(op, operon_run, hits)
        for clan in ("CL00111", "CL00112"):
            mine = [r for r in here if r["clan"] == clan]
            kept = [r for r in mine if r["clan_winner"]]
            assert len(mine) >= 2 and len(kept) == 1, (op, clan, mine)
            assert kept[0]["bitscore"] == max(r["bitscore"] for r in mine)
            assert all(_overlap(kept[0], r) for r in mine)  # the losers lie on the kept hit's gene
        assert all(r["clan_winner"] for r in here if r["model_acc"] in ("RF00001", "RF00005", "RF01690"))
    tsec = _in("trna_sec", operon_run, hits)
    sec = [r for r in tsec if r["model_acc"] == "RF01852"]
    for s in sec:
        rivals = [t for t in tsec if t["model_acc"] == "RF00005" and t["strand"] == s["strand"] and _overlap(s, t)]
        assert rivals and sum(r["clan_winner"] for r in [s, *rivals]) == 1
    assert sum(r["model_acc"] == "RF00005" and r["clan_winner"] for r in tsec) >= 4
    for name in ("glm", "srp"):
        here = _in(name, operon_run, hits)
        assert len(here) == 2 and sum(r["clan_winner"] for r in here) == 1, here
        assert max(here, key=lambda r: r["bitscore"])["clan_winner"]
    for name in ("sdsr_srac", "trna_hpnc", "trna_rtt"):
        here = _in(name, operon_run, hits)
        assert here and all(r["clan_winner"] for r in here), here
        if name != "trna_rtt":  # the antisense partners: a hit on each strand, overlapping
            assert any(a["strand"] == 1 and b["strand"] == -1 and _overlap(a, b) for a in here for b in here), here
        else:  # tRNA and RtT nested on one strand
            assert any(a["model_acc"] == "RF00005" and b["model_acc"] == "RF00391" and a["strand"] == b["strand"]
                       and _overlap(a, b) for a in here for b in here), here
    # no two kept hits of one clan overlap on a strand
    kept = [r for r in hits if r["clan_winner"] and r["clan"] is not None]
    for a in kept:
        for b in kept:
            assert a is b or not (a["clan"] == b["clan"] and a["strand"] == b["strand"] and _overlap(a, b))
    # the features are the distinct loci of the kept hits; every losing hit is in l1_ncrna with its flag
    want = sorted({(r["feature_uid"], r["start"], r["end"], r["strand"]) for r in hits if r["clan_winner"]},
                  key=lambda t: t[1:])
    assert [(f["feature_uid"], f["start"], f["end"], f["strand"]) for f in feats] == want
    assert operon_run["report"]["clan_losers"] == sum(not r["clan_winner"] for r in hits)
    print(f"\nG-04b: {len(hits)} hits, {len(feats)} features, {operon_run['report']['clan_losers']} clan losers: "
          + ", ".join(sorted(f"{r['model_name']}@{r['start']}" for r in hits if not r["clan_winner"])))

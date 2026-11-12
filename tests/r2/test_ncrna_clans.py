"""Ring R2, card E-08 (DEC-28): clan competition among the ncRNA producer's hits, and the clan file.

The rule under test (`biobox.l1.ncrna.compete`): two hits compete when they lie on the same contig and strand, share
at least one base (on a circular contig, around the origin too: a hit across it has end > contig_len), and their
models are in one clan of `Rfam.clanin`. Rank is the bit score, higher first, then the model accession, then (start,
end), so no two hits tie. Taken in rank order, a hit is kept when it shares no base with a competitor already kept,
else it loses. A hit of a model in no clan never competes.

Oracles (AGENTS §5.2), written here, separately from the producer:

- a brute-force loop: overlap by intersecting the sets of bases each hit covers (positions modulo the contig length),
  the hits taken in rank order by a sort of this file's own key, each compared with every hit kept before it. It is
  compared with the producer's function on seeded random hit sets (ties, nested hits, opposite strands, hits across
  the origin, models in no clan) and on hand-made edge cases;
- the `cmscan` binary (Infernal 1.1.5, `--clanin --oclan --fmt 2`) on constructed cases: single-sequence models built
  with `cmbuild` from segments of a seeded random sequence, so their hits fall where the segments are (chains of
  three and four hits of one clan, a one-base overlap, a tie, opposite strands, a model in no clan nested in one of a
  clan). The producer's function, given cmscan's own hits, must drop every hit cmscan marks `=` and keep every hit it
  marks `*` or `^`; cmscan keeps every `$` hit (one whose better competitors all overlap better ones), the producer
  keeps it unless a kept better hit overlaps it, which happens in the chain of four only.

The clan mapping of the loop's cases is written by hand (the function takes a mapping); the repository's clan file is
read by `ncrna.read_clans`, checked here against its sha256 in refdb.yaml and a parse of this file.
"""

from __future__ import annotations

import hashlib
import random
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-08")]

REPO = Path(__file__).resolve().parents[2]
CLANIN = REPO / "src/biobox/l1/rfam_sets/Rfam-15.1.clanin"
# a hand-written clan mapping: model name -> clan; the names not listed are in no clan
CLANS = {"tRNA": "CL00001", "tRNA-Sec": "CL00001", "tmRNA": "CL00001",
         "SSU_rRNA_bacteria": "CL00111", "SSU_rRNA_archaea": "CL00111",
         "LSU_rRNA_bacteria": "CL00112", "LSU_rRNA_archaea": "CL00112"}
MODELS = {"tRNA": "RF00005", "tRNA-Sec": "RF01852", "tmRNA": "RF00023", "SSU_rRNA_bacteria": "RF00177",
          "SSU_rRNA_archaea": "RF01959", "LSU_rRNA_bacteria": "RF02541", "LSU_rRNA_archaea": "RF02540",
          "RtT": "RF00391", "Bacillaceae-1": "RF01690"}


def _ncrna():
    from biobox.l1 import ncrna  # imported here so that a missing function fails the tests, not the collection

    return ncrna


# --- the oracle --------------------------------------------------------------------------------------------------


def bases(h: dict) -> set[int]:
    n = h["contig_len"]
    return {x % n for x in range(h["start"], h["end"])}


def clan_file() -> dict[str, str]:
    """The repository's clan file, parsed here: model name -> clan accession (each line a clan, then its members)."""
    out = {}
    for ln in CLANIN.read_text().splitlines():
        acc, *names = ln.split("\t")
        for m in names:
            out[m] = acc
    return out


def better(a: dict, b: dict) -> bool:
    """a ranks above b: a higher bit score, then the smaller accession, then the smaller (start, end)."""
    if a["bitscore"] != b["bitscore"]:
        return a["bitscore"] > b["bitscore"]
    if a["model_acc"] != b["model_acc"]:
        return a["model_acc"] < b["model_acc"]
    return (a["start"], a["end"]) < (b["start"], b["end"])


def oracle(hits: list[dict], clans: dict[str, str]) -> list[tuple[str | None, bool]]:
    """Each hit's clan and whether it is kept: in rank order (a selection sort by `better`), kept when no hit kept
    before it has its clan and strand and shares a base with it."""
    order, left = [], list(range(len(hits)))
    while left:
        best = left[0]
        for j in left[1:]:
            if better(hits[j], hits[best]):
                best = j
        order.append(best)
        left.remove(best)
    kept: list[int] = []
    won = [True] * len(hits)
    for i in order:
        h, c = hits[i], clans.get(hits[i]["model_name"])
        if c is not None and any(clans.get(hits[j]["model_name"]) == c and hits[j]["strand"] == h["strand"]
                                 and bases(hits[j]) & bases(h) for j in kept):
            won[i] = False
        else:
            kept.append(i)
    return [(clans.get(h["model_name"]), w) for h, w in zip(hits, won)]


def hit(name: str, start: int, end: int, strand: int, score: float, n: int) -> dict:
    return {"model_name": name, "model_acc": MODELS[name], "start": start, "end": end, "strand": strand,
            "bitscore": score, "contig_len": n}


def random_hits(rng: random.Random, n: int, circular: bool) -> list[dict]:
    out, seen = [], set()
    for _ in range(rng.randint(0, 25)):
        name = rng.choice(sorted(MODELS))
        length = rng.randint(1, min(n - 1, rng.choice([30, 120, 400])))
        start = rng.randrange(n) if circular else rng.randrange(n - length + 1)
        strand = rng.choice((1, -1))
        score = rng.choice([25.0, 40.5, 40.5, 61.25, 90.0, round(rng.uniform(20, 200), 1)])  # ties on purpose
        key = (name, start, start + length, strand)
        if key in seen:  # one row per (locus, model): the table's key
            continue
        seen.add(key)
        out.append(hit(name, start, start + length, strand, score, n))
    return out


# --- the competition -----------------------------------------------------------------------------------------------


EDGES = [
    # ties: the same locus and score, two models of one clan; the smaller accession wins
    [hit("tRNA", 100, 172, 1, 50.0, 1000), hit("tRNA-Sec", 100, 172, 1, 50.0, 1000)],
    # nested, different clans (and one in no clan): all kept
    [hit("SSU_rRNA_bacteria", 100, 900, 1, 1500.0, 2000), hit("tRNA", 300, 372, 1, 60.0, 2000),
     hit("RtT", 310, 480, 1, 55.0, 2000)],
    # nested, one clan: the inner one loses to the better outer one
    [hit("SSU_rRNA_bacteria", 100, 900, 1, 1500.0, 2000), hit("SSU_rRNA_archaea", 105, 895, 1, 1000.0, 2000)],
    # opposite strands, one clan: both kept
    [hit("tRNA", 100, 172, 1, 50.0, 1000), hit("tRNA-Sec", 100, 172, -1, 70.0, 1000)],
    # one shared base is an overlap; touching is not
    [hit("tRNA", 100, 172, 1, 50.0, 1000), hit("tmRNA", 171, 500, 1, 90.0, 1000),
     hit("tRNA-Sec", 500, 590, 1, 30.0, 1000)],
    # across the origin of a circular contig (end > n) against a hit at its start, and one just short of it
    [hit("tRNA", 990, 1030, 1, 50.0, 1000), hit("tRNA-Sec", 29, 60, 1, 40.0, 1000),
     hit("tmRNA", 30, 90, 1, 45.0, 1000)],
    # a chain: A beats B, B overlaps C but A does not; C ranks below B, which loses, so C is kept
    [hit("tmRNA", 100, 400, 1, 300.0, 2000), hit("tRNA", 350, 450, 1, 200.0, 2000),
     hit("tRNA-Sec", 440, 600, 1, 100.0, 2000)],
    # a chain of four, each overlapping the next only: the first and third are kept
    [hit("SSU_rRNA_bacteria", 0, 70, -1, 70.0, 2000), hit("SSU_rRNA_archaea", 60, 120, -1, 60.0, 2000),
     hit("SSU_rRNA_bacteria", 110, 160, -1, 50.0, 2000), hit("SSU_rRNA_archaea", 150, 190, -1, 40.0, 2000)],
    # models in no clan never compete, even with the same locus
    [hit("RtT", 10, 80, 1, 50.0, 1000), hit("Bacillaceae-1", 10, 80, 1, 70.0, 1000)],
    [],
]


@pytest.mark.parametrize("case", range(len(EDGES)))
def test_competition_equals_the_brute_force_loop_on_edge_cases(case):
    hits = EDGES[case]
    assert _ncrna().compete(hits, CLANS) == oracle(hits, CLANS)


def test_competition_edge_cases_say_what_the_rule_is():
    """The expected winners of the edge cases, written out (the oracle's reading of the rule, made explicit)."""
    nc = _ncrna()
    won = [[w for _, w in nc.compete(h, CLANS)] for h in EDGES]
    assert won[0] == [True, False]  # RF00005 < RF01852
    assert won[1] == [True, True, True]
    assert won[2] == [True, False]
    assert won[3] == [True, True]
    assert won[4] == [False, True, True]  # tRNA shares base 171 with tmRNA; tRNA-Sec starts where tmRNA ends
    assert won[5] == [True, False, True]  # base 29 is in both; tmRNA starts at 30, after tRNA's end (30), beats tRNA-Sec
    assert won[6] == [True, False, True]  # C is overlapped only by B, which A beat
    assert won[7] == [True, False, True, False]
    assert won[8] == [True, True]
    assert [c for c, _ in nc.compete(EDGES[8], CLANS)] == [None, None]


@pytest.mark.covers("MR-04", "C-05")
def test_competition_equals_the_brute_force_loop_on_random_hit_sets():
    nc = _ncrna()
    rng = random.Random(20261111)
    checked = lost = cross = 0
    for _ in range(400):
        circular = rng.random() < 0.5
        n = rng.choice([150, 600, 2500])
        hits = random_hits(rng, n, circular)
        want = oracle(hits, CLANS)
        assert nc.compete(hits, CLANS) == want, hits
        # the answer does not depend on the input order
        order = list(range(len(hits)))
        rng.shuffle(order)
        got = nc.compete([hits[i] for i in order], CLANS)
        assert got == [want[i] for i in order]
        checked += len(hits)
        lost += sum(not w for _, w in want)
        cross += sum(h["end"] > n for h in hits)
    assert checked > 3000 and lost > 300 and cross > 100  # the sets exercise losing and the origin


# --- against cmscan on constructed cases ---------------------------------------------------------------------------


CMBUILD, CMPRESS, CMSCAN = (shutil.which(t) for t in ("cmbuild", "cmpress", "cmscan"))


def parse_fmt2(path: Path) -> list[dict]:
    """`cmscan --fmt 2 --tblout` rows: model name, clan, sequence coordinates (to 0-based half-open), strand, score,
    olp. A model without an accession is named by its name, as `ncrna.read_models` does."""
    out = []
    for ln in path.read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        lo, hi = sorted((int(f[9]), int(f[10])))
        out.append({"model_name": f[1], "model_acc": f[1] if f[2] == "-" else f[2],
                    "clan": None if f[5] == "-" else f[5], "start": lo - 1, "end": hi,
                    "strand": 1 if f[11] == "+" else -1, "bitscore": float(f[16]), "olp": f[19]})
    return out


@pytest.mark.skipif(None in (CMBUILD, CMPRESS, CMSCAN), reason="Infernal (cmbuild, cmpress, cmscan) not on PATH; it "
                                                               "is in the test image")
def test_competition_against_cmscan_on_constructed_cases(tmp_path):
    rng = random.Random(20261111)

    def rs(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    def rc(s: str) -> str:
        return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    x1, x3, y, z, w, t, s = rs(60), rs(60), rs(80), rs(50), rs(50), rs(40), rs(220)
    segments = {
        "cA": x1, "cB": x1[40:] + x3[:20], "cC": x3[10:40],  # chain of three: A > B > C, A and C apart
        "dD": y[:50], "dE": y[49:79],  # one shared base
        "nF": z, "nG": z[10:40],  # G, in no clan, nested in F
        "sI": w, "sJ": rc(w[10:40]),  # one clan, opposite strands
        "tK": t, "tL": t,  # a tie: the same segment
        "p1": s[0:70], "p2": s[60:120], "p3": s[110:160], "p4": s[150:190],  # chain of four
    }
    clans = {"cA": "CLX01", "cB": "CLX01", "cC": "CLX01", "dD": "CLX02", "dE": "CLX02", "nF": "CLX03",
             "sI": "CLX03", "sJ": "CLX03", "tK": "CLX04", "tL": "CLX04", "p1": "CLX05", "p2": "CLX05",
             "p3": "CLX05", "p4": "CLX05"}
    contig = "".join([rs(300), x1 + x3, rs(300), y, rs(300), z, rs(300), w, rs(300), t, rs(300), s, rs(300)])
    fa = tmp_path / "contig.fa"
    fa.write_text(">c\n" + contig + "\n")
    db = tmp_path / "db.cm"
    with open(db, "w") as out:
        for name, seq in segments.items():
            sto, cm = tmp_path / f"{name}.sto", tmp_path / f"{name}.cm"
            sto.write_text(f"# STOCKHOLM 1.0\n\ns1 {seq}\n#=GC SS_cons {'.' * len(seq)}\n//\n")
            subprocess.run([CMBUILD, "-F", "--enone", "-n", name, str(cm), str(sto)], check=True, capture_output=True)
            out.write(cm.read_text())
    subprocess.run([CMPRESS, "-F", str(db)], check=True, capture_output=True)
    clanin = tmp_path / "clanin"
    by: dict[str, list[str]] = {}
    for m, c in clans.items():
        by.setdefault(c, []).append(m)
    clanin.write_text("".join(c + "\t" + "\t".join(ms) + "\n" for c, ms in sorted(by.items())))
    tbl = tmp_path / "hits.tbl"
    subprocess.run([CMSCAN, "--cpu", "1", "-T", "10", "--incT", "10", "--notrunc", "--noali", "--fmt", "2", "--clanin",
                    str(clanin), "--oclan", "--tblout", str(tbl), "-o", str(tmp_path / "hits.out"), str(db), str(fa)],
                   check=True, capture_output=True)
    scanned = parse_fmt2(tbl)
    assert sorted(h["model_name"] for h in scanned) == sorted(segments)  # each model hits once, where it was cut
    hits = [{**h, "contig_len": len(contig)} for h in scanned]
    got = _ncrna().compete(hits, clans)
    assert got == oracle(hits, clans)
    olp = {h["model_name"]: h["olp"] for h in scanned}
    won = {h["model_name"]: w for h, (_, w) in zip(hits, got)}
    for name in segments:
        if olp[name] == "=":
            assert not won[name], name
        elif olp[name] in "*^":
            assert won[name], name
    # what the binary shows: one shared base is an overlap; the chain of three keeps its third; strands and a model
    # in no clan never compete; in the chain of four cmscan keeps p4 (`$`, overlapped only by p3, itself `$`) where
    # the producer drops it under the kept p3
    assert (olp["dE"], olp["cB"], olp["cC"], olp["tL"]) == ("=", "=", "$", "=")
    assert {olp[m] for m in ("nF", "nG", "sI", "sJ")} == {"*"}
    assert [olp[m] for m in ("p1", "p2", "p3", "p4")] == ["^", "=", "$", "$"]
    assert [won[m] for m in ("p1", "p2", "p3", "p4")] == [True, False, True, False]
    assert won["cC"] and won["tK"] and not won["tL"]


# --- the clan file -------------------------------------------------------------------------------------------------


def test_the_clan_file_is_refdbs_and_parses_to_its_lines():
    nc = _ncrna()
    entry = next(d for d in yaml.safe_load((REPO / "refdb.yaml").read_text())["databases"] if d["name"] == "rfam_clanin")
    assert hashlib.sha256(CLANIN.read_bytes()).hexdigest() == entry["sha256"]
    want = clan_file()
    assert len({a for a in want.values()}) == 151
    got = nc.read_clans()
    assert got == want
    assert got["tRNA"] == got["tRNA-Sec"] == "CL00001" and "RtT" not in got
    assert nc.clan_version() == "15.1"


def test_a_clan_file_that_differs_is_refused(tmp_path):
    nc = _ncrna()
    bad = tmp_path / "Rfam.clanin"
    shutil.copy(CLANIN, bad)
    bad.write_text(bad.read_text().replace("\ttRNA-Sec", ""))
    with pytest.raises(ValueError, match="sha256"):
        nc.read_clans(bad)


@pytest.mark.covers("C-08")
def test_provenance_accepts_the_producer_with_set_region_and_clan_tails():
    """E-07 report (a): `NCRNA_RE` names Infernal's and Rfam's versions with an optional `-<set>-<region>` tail
    (card E-07) and an optional `-<clan rule>` after it (this card)."""
    from biobox import provenance

    names = {"infernal-1.1.5-rfam-15.1": True, "infernal-1.1.5-rfam-15.1-prok1-igr200nt": True,
             "infernal-1.1.5-rfam-15.1-prok1-igr200nt-clan1": True,
             "infernal-1.1.5-rfam-15.1-sha0123456789ab-igr200nt-clan1": True,
             "infernal-1.1.5-rfam-15.1-prok1": False, "infernal-1.1.5-rfam-15.1-prok1-igr200nt-x1": False,
             "infernal-1.1.5-rfam-15.1-clan1": False}
    assert {n: provenance.NCRNA_RE.match(n) is not None for n in names} == names
    lock = provenance.locked_versions()
    for p in ("infernal-1.1.5-rfam-15.1-prok1-igr200nt", "infernal-1.1.5-rfam-15.1-prok1-igr200nt-clan1"):
        assert provenance._producer_problems("x", "p", p.replace("1.1.5", lock["infernal"]), lock) == []
    assert provenance._producer_problems("x", "p", "infernal-0.0-rfam-15.1-prok1-igr200nt-clan1", lock)
    assert provenance._producer_problems("x", "p", f"infernal-{lock['infernal']}-rfam-15.0-prok1-igr200nt-clan1", lock)


@pytest.mark.covers("C-08")
def test_provenance_checks_the_gene_call_version_an_ncrna_ledger_names(tmp_path):
    """E-07 report (c): the ncRNA ledger names the gene-call version whose intergenic regions were scanned
    (`gene_call_version` among its tool versions); C-08 checks it as it checks a gene-call partition (its pyrodigal
    versions equal the lock), not as a tool of the lock."""
    from biobox import ledger, provenance
    from biobox.l1.genecall import GENE_CALL_VERSION as V

    lock = provenance.locked_versions()
    tools = {"infernal": lock["infernal"], "pyarrow": lock["pyarrow"], "rfam": "15.1", "rfam_clanin": "15.1"}

    def manifest(step: str, gcv: str) -> dict:
        ledger.write(tmp_path, "S", step, [ledger.row(b"\x01" * 16, "contig", 1, 1, tools={**tools,
                                                                                         "gene_call_version": gcv})])
        rel = ledger.path(tmp_path, "S", step).relative_to(tmp_path).as_posix()
        return {"root": str(tmp_path), "files": [{"table": "ledger", "path": rel,
                                                  "partitions": {"stage": "S", "step": step}}]}

    assert provenance.check_provenance(manifest("ncrna-a", V), lock) == []
    probs = provenance.check_provenance(manifest("ncrna-b", "pyrodigal-0.0.1-gv-0.0.1-p1"), lock)
    assert len(probs) == 1 and "gene_call_version" in probs[0]
    assert provenance.check_provenance(manifest("ncrna-c", "not-a-version"), lock)


def test_the_producer_name_carries_the_clan_rule():
    nc = _ncrna()
    assert nc.CLAN_RULE == "clan1"
    lib = REPO / "tests/fixtures/micro_s0/rfam15.1_micro_s0.cm.gz"
    if shutil.which("cmsearch") is None:
        pytest.skip("cmsearch (Infernal) not on PATH; it is in the test image")
    assert nc.producer(lib).endswith("-prok1-igr200nt-clan1")

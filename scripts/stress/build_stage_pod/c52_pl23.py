"""Card C-52 (DEC-63): stages PL2 and PL3 from card C-40's staging of ten plasmid sets, with the sources' circular flags.

    python c52_pl23.py survey   --inputs IN --out DIR [--long-bound N] [--pl1 PL1_IN]
    python c52_pl23.py deferred --gut-table TSV --skipped TSV --out FILE
    python c52_pl23.py stage    --inputs IN --survey DIR --root ROOT --out-pl2 DIR --out-pl3 DIR --work DIR
                                [--final-pl2 DIR] [--final-pl3 DIR] [--long-bound N] [--fragments linear|circular]
    python c52_pl23.py overlap  --staging PL2=DIR --staging PL3=DIR --root ROOT [--other STAGE=DIR ...] --work DIR
                                --out FILE

IN is what `make_inputs.py plasmid-sources` wrote for card C-40 (`objects.json`, the FASTA objects, `records.tsv`, whose
rows are the objects' records in order). Nothing in IN is written.

The source's flag (DEC-63 (3): "use circular if that is the source flag"), by set: `gut_plasmids_yu2024`, its table's
column `circular` = True; PlasmidScope's `IMG-PR`, its column `Topology` = `direct terminal repeat` (IMG-PR never says
`circular`); every other PlasmidScope set (`mMGE`, `RefSeq`, `GenBank`, `ENA`, `DDBJ`, `COMPASS`, `TPA`, `Kraken2`),
`Topology` = `circular`. Any other value is not a circular flag; a source or value the rule does not know is refused.

The two ends of a flagged record (`ends`, on the sequence in upper case): `exact` is the longest k, from `SEED` (20,
DEC-63's "a repeat of 20 bp or more") up to half the record, with the record's first k bases equal to its last k; the
record then holds the terminal repeat twice, and `trim` removes the second copy (the last k bases), so that the circle
carries it once and the record's start is unchanged. A further agreement of the ends is any place within `WINDOW`
bases of the end where the record's first `SEED` bases occur again, or within `WINDOW` bases of the start where its
last `SEED` bases do, whose span is longer than the exact repeat and which is not itself an exact repeat of the ends.
It is `near` when it looks like a second copy of a repeat (identity at least `NEAR_IDENTITY` over its span: a repeat
with substitutions; or an unbroken run from its seed over half its span: a repeat broken by an indel), else `far` (a
shared stretch, then the ends differ: an internal repeat, counted and not acted on). `classify`: `exact` (one exact
repeat and no near agreement: the trim is mechanical), `none` (no repeat and no near agreement: the record is the
circle as it is), and the classes a trim cannot settle: `partial` (a near agreement: the ends agree in part, with or
without a shorter exact repeat), `low_complexity` (the exact repeat holds an N, or has a period of 6 or less, or one
base is 80% of it), `long` (the exact repeat is longer than `--long-bound`), `unresolved` (more than `MAX_HITS` places
to look at: a low-complexity start).

survey: every record its source flags circular, with its class, into DIR/survey.tsv, and DIR/summary.json (per source:
records, flagged, flagged and incomplete, classes, the lengths of the exact repeats, the first five accessions of each
class, the classes among the records that are a genome the stage would store). With --pl1, the same over the circular
genomes of stage PL1's inputs (card C-26, PLSDB), for comparison. Streams; refuses (exit 2, no summary) when
`records.tsv` does not line up with the objects' records.

deferred (DEC-59): the gut table's contigs that are circular and not among the authors' representatives
(`in_dereplicated_plasmids` False, `circular` True), sorted by contig, with the table's cluster, cluster type and the
cluster's representative, the contig's length and its `contig_uid` from the FASTA record (the staging's
`skipped_records.tsv` lists every gut record not taken); a contig with no record is refused (exit 2, nothing written).

stage (DEC-63): two inputs directories from IN's records, none left out and none added. `PL2` takes the objects whose
`redistributable` is true, `PL3` the others, each in IN's order. A record is staged circular when its source flags it
(and, for a record its source also calls incomplete, only under `--fragments circular`; under `--fragments linear` it
is linear, card C-40's rule; with neither, such a record is refused). A circular record of class `exact` is trimmed
(`trim`) and its `contig_uid` is taken from the trimmed bytes; one of class `none` is staged as it is; a record of any
other class is refused unless `--ambiguous` says how to stage it, untrimmed in every case (`linear`: card C-39's rule
for a terminal repeat left in the record; `circular`: the source's flag as it stands; `asis`: the topology card C-40's
staging gave the record, which is linear for IMG-PR and circular for the other sets). The
survey's classes are checked before anything is written and again record by record. An object
with no trim is a hard link to IN's file; one with a trim is written again, record by record, with every other record
byte for byte. Per stage: the FASTA objects; `objects.json` (IN's entries, each path in the stage's final directory,
`--final-*`, default the output directory: the orchestrator renames a staging into place); `genomes.json` (one genome
per distinct sequence over the stage's sources, the first record taken, as card C-40's staging: kind plasmid, code 11,
circular from that record, assembly_quality as C-40 mapped it); `records.tsv` (C-40's columns, the first record and the
owning stage within this stage and the root, and the record's source `contig_uid`, length and trimmed bases);
`trims.tsv` (every trim: record, repeat length, bases removed); `skipped_records.tsv` (IN's, for the stage's sources);
`stage_report.json`. Each stage is written into a sibling scratch directory and renamed into place only when both are
complete. Bounded memory: records stream, and the per-record rows go through Parquet and DuckDB.

overlap: what each new staging's contigs (as staged) share with the root's stages, with stage PL1's genome_input, with
staged inputs of stages not yet built (`--other`), and with each other; beside the same counts for the records' source
bytes, before any trim.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from biobox import fasta, ids

csv.field_size_limit(1 << 30)

SEED = 20
WINDOW = 10_000
MAX_HITS = 1_000
NEAR_IDENTITY = 0.9
LONG_BOUND = 1_000
LOW_PERIOD = 6
LOW_TOP = 0.8

GUT = "gut_plasmids_yu2024"
SCOPE = "plasmidscope_"
METAGENOMIC_SETS = frozenset({GUT, SCOPE + "IMG-PR", SCOPE + "mMGE"})
PUBLIC_SETS = frozenset({SCOPE + s for s in ("RefSeq", "GenBank", "ENA", "DDBJ", "COMPASS", "TPA", "Kraken2")})
SCOPE_VALUES = ("circular", "linear", "direct terminal repeat", "inverted terminal repeat", "concatemer")
CLASSES = ("exact", "none", "partial", "low_complexity", "long", "unresolved")
REPEAT_BINS = ((20, 100), (100, 1000), (1000, 5000), (5000, None))


class Refused(Exception):
    """Inputs the script will not act on (exit 2)."""


def source_flag(source: str, value: str) -> bool:
    """Whether the source marks the record circular (module docstring)."""
    if source == GUT:
        if value not in ("True", "False"):
            raise Refused(f"{source}: circular {value!r} is not True or False")
        return value == "True"
    if source.startswith(SCOPE):
        if value not in SCOPE_VALUES:
            raise Refused(f"{source}: Topology {value!r} is not one of {SCOPE_VALUES}")
        if source == SCOPE + "IMG-PR":
            return value == "direct terminal repeat"
        return value == "circular"
    raise Refused(f"source {source!r}: no rule for its circular flag")


def fragment(source: str, completeness: str) -> bool:
    """Whether the source calls the record incomplete (card C-40's `fragment`: PlasmidScope's `incomplete`)."""
    return source.startswith(SCOPE) and completeness == "incomplete"


@dataclass(frozen=True)
class Ends:
    n: int
    exact: int = 0  # the exact terminal repeat's length (0: none)
    near: int = 0  # other agreements of the ends, longer than the exact repeat, that look like a repeat
    far: int = 0  # other agreements: a shared stretch, then the ends differ (an internal repeat; not a terminal one)
    near_len: int = 0  # the longest near agreement: its span, mismatches and unbroken run from its seed
    near_mm: int = 0
    near_run: int = 0
    period: int = 0  # of the exact repeat: least period up to LOW_PERIOD (0: none), largest base share, N count
    top: float = 0.0
    n_bases: int = 0
    capped: bool = False


def _same(s: bytes, a: int, b: int, length: int) -> bool:
    """s[a:a+length] == s[b:b+length], in growing chunks so that a mismatch near the start costs little."""
    i, step = 0, 64
    while i < length:
        j = min(length, i + step)
        if s[a + i:a + j] != s[b + i:b + j]:
            return False
        i, step = j, step * 8
    return True


def _mismatches(x: bytes, y: bytes) -> int:
    return sum(1 for p, q in zip(x, y) if p != q)


def _run(x: bytes, y: bytes) -> int:
    """The length of the common prefix of x and y."""
    k = 0
    for p, q in zip(x, y):
        if p != q:
            break
        k += 1
    return k


def _is_near(span: int, mm: int, run: int) -> bool:
    """An agreement that looks like a second copy of a repeat: identity at least NEAR_IDENTITY over the span (a repeat
    with substitutions), or an unbroken run from the seed over half the span or more (a repeat broken by an indel)."""
    return 1 - mm / span >= NEAR_IDENTITY or 2 * run >= span


def _period(r: bytes) -> int:
    for p in range(1, LOW_PERIOD + 1):
        if len(r) > p and r[p:] == r[:-p]:
            return p
    return 0


def ends(seq: bytes) -> Ends:
    """The two ends of a record (module docstring)."""
    s = seq.upper()
    n = len(s)
    if n < 2 * SEED:
        return Ends(n)
    head, tail = s[:SEED], s[n - SEED:]
    lo = n - n // 2  # a start p >= lo is a repeat of k = n - p <= n // 2
    exact, tries, capped = 0, 0, False
    p = s.find(head, lo)
    while p != -1:
        tries += 1
        if tries > MAX_HITS:
            capped = True
            break
        if _same(s, p, 0, n - p):
            exact = n - p
            break
        p = s.find(head, p + 1)
    # other agreements of the ends: (span, mismatches, unbroken run from the seed). One inside the exact repeat's own
    # second copy (span shorter than it) is the repeat's inner structure and is not counted; nor is a shorter exact
    # repeat of the ends (it is part of the longest).
    found = []
    p = s.find(head, max(lo, n - WINDOW))
    k = 0
    while p != -1 and not capped:
        k += 1
        if k > MAX_HITS:
            capped = True
            break
        span = n - p
        if span > exact and not _same(s, p, 0, span):
            found.append((span, _mismatches(s[p:], s[:span]), _run(s[p:], s[:span])))
        p = s.find(head, p + 1)
    top_q = min(WINDOW, n // 2) - SEED  # an end-seed at q covers s[:q + SEED]
    q = s.find(tail, 0, top_q + SEED) if top_q >= 0 else -1
    k = 0
    while q != -1 and not capped:
        k += 1
        if k > MAX_HITS:
            capped = True
            break
        span = q + SEED
        if span > exact and not _same(s, 0, n - span, span):
            found.append((span, _mismatches(s[:span], s[n - span:]), _run(s[:span][::-1], s[n - span:][::-1])))
        q = s.find(tail, q + 1, top_q + SEED)
    near = [f for f in found if _is_near(*f)]
    near_len, near_mm, near_run = max(near, key=lambda t: (t[0], -t[1])) if near else (0, 0, 0)
    period, top, n_bases = 0, 0.0, 0
    if exact:
        r = s[:exact]
        period = _period(r)
        top = max(r.count(b) for b in (b"A", b"C", b"G", b"T")) / exact
        n_bases = r.count(b"N")
    return Ends(n, exact, len(near), len(found) - len(near), near_len, near_mm, near_run, period, round(top, 4),
                n_bases, capped)


def low_complexity(e: Ends) -> bool:
    return e.exact > 0 and (e.n_bases > 0 or e.period > 0 or e.top >= LOW_TOP)


def classify(e: Ends, long_bound: int = LONG_BOUND) -> str:
    if e.capped:
        return "unresolved"
    if e.near:
        return "partial"
    if e.exact:
        if low_complexity(e):
            return "low_complexity"
        if e.exact > long_bound:
            return "long"
        return "exact"
    return "none"


def trim(seq: bytes, k: int) -> bytes:
    """The record without the second copy of its terminal repeat of k bases (its last k)."""
    if not 0 < k <= len(seq) // 2:
        raise Refused(f"a trim of {k} bases from a record of {len(seq)}")
    return seq[:len(seq) - k]


def repeat_bin(k: int) -> str:
    for lo, hi in REPEAT_BINS:
        if k >= lo and (hi is None or k < hi):
            return f"{lo}-{hi - 1}" if hi else f">={lo}"
    return "<20"


# --------------------------------------------------------------------------------------------------------------
# survey


def _tsv_rows(path: Path):
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)


def records_with_rows(inputs: Path):
    """(object, record, its records.tsv row) for every record of the staging's objects, in order; Refused when the two
    do not line up (another accession or length, a row too many or too few)."""
    objs = json.loads((inputs / "objects.json").read_text())
    rows = _tsv_rows(inputs / "records.tsv")
    for o in objs:
        stem = Path(o["path"]).stem
        for r in fasta.iter_records(Path(o["path"])):
            row = next(rows, None)
            if row is None or row["source"] != stem or row["accession"] != r.accession or int(row["length"]) != len(r.seq):
                raise Refused(f"records.tsv does not line up with {o['path']} at {r.accession} "
                              f"(row {None if row is None else (row['source'], row['accession'], row['length'])})")
            yield o, r, row
    extra = next(rows, None)
    if extra is not None:
        raise Refused(f"records.tsv does not line up: row {extra['ordinal']} has no record in the objects")


SURVEY_HEAD = ["ordinal", "source", "accession", "contig_uid", "length", "topology_value", "completeness_value",
               "fragment", "genome_stored", "exact", "near", "far", "near_len", "near_mm", "near_run", "period", "top",
               "n_bases", "class"]


def _new_tally() -> dict:
    return {"records": 0, "flagged": 0, "flagged_and_fragment": 0, "classes": Counter(), "stored_classes": Counter(),
            "fragment_classes": Counter(), "repeat_bins": Counter(), "repeat_lengths": Counter(),
            "repeat_share_of_record": Counter(), "examples": {}}


def _share_bin(k: int, n: int) -> str:
    f = k / n
    for edge in (0.01, 0.05, 0.1, 0.25):
        if f < edge:
            return f"<{edge}"
    return ">=0.25"


def _tally(t: dict, acc: str, e: Ends, cls: str, frag: bool, stored: bool) -> None:
    t["flagged"] += 1
    t["flagged_and_fragment"] += frag
    t["classes"][cls] += 1
    if stored:
        t["stored_classes"][cls] += 1
    if frag:
        t["fragment_classes"][cls] += 1
    if e.exact:
        t["repeat_bins"][repeat_bin(e.exact)] += 1
        t["repeat_lengths"][e.exact] += 1
        t["repeat_share_of_record"][_share_bin(e.exact, e.n)] += 1
    ex = t["examples"].setdefault(cls, [])
    if len(ex) < 5:
        ex.append(acc)


def _done(t: dict) -> dict:
    out = dict(t)
    for k in ("classes", "stored_classes", "fragment_classes", "repeat_bins", "repeat_share_of_record"):
        out[k] = dict(sorted(t[k].items()))
    out["repeat_lengths_top"] = [[k, v] for k, v in sorted(t["repeat_lengths"].items(), key=lambda kv: (-kv[1], kv[0]))[:25]]
    del out["repeat_lengths"]
    return out


def survey(inputs: Path, out: Path, long_bound: int, pl1: Path | None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    tallies: dict[str, dict] = {}
    tmp = out / ".survey.tsv.tmp"
    with open(tmp, "w") as w:
        w.write("\t".join(SURVEY_HEAD) + "\n")
        for _o, r, row in records_with_rows(inputs):
            src = row["source"]
            t = tallies.setdefault(src, _new_tally())
            t["records"] += 1
            if not source_flag(src, row["topology_value"]):
                continue
            frag = fragment(src, row["completeness_value"])
            stored = not row["first_accession"] and not row["owner_stage"]
            e = ends(r.seq)
            cls = classify(e, long_bound)
            _tally(t, r.accession, e, cls, frag, stored)
            w.write("\t".join(map(str, [
                row["ordinal"], src, r.accession, row["contig_uid"], len(r.seq), row["topology_value"],
                row["completeness_value"], str(frag).lower(), str(stored).lower(), e.exact, e.near, e.far,
                e.near_len, e.near_mm, e.near_run, e.period, e.top, e.n_bases, cls])) + "\n")
    doc = {"long_bound": long_bound, "seed": SEED, "window": WINDOW, "near_identity": NEAR_IDENTITY,
           "sources": {s: _done(t) for s, t in tallies.items()}}
    if pl1 is not None:
        doc["pl1"] = survey_pl1(pl1, long_bound)
    tmp.replace(out / "survey.tsv")
    (out / "summary.json").write_text(json.dumps(doc, indent=1, sort_keys=True))
    return doc


def survey_pl1(pl1: Path, long_bound: int) -> dict:
    """The circular genomes of stage PL1's inputs (one genome per distinct sequence, `circular` lists its contig):
    their ends, as PL1 staged them (circular, untrimmed)."""
    circ = set()
    for g in json.loads((pl1 / "genomes.json").read_text()):
        circ.update(g.get("circular", ()))
    t = _new_tally()
    seen = set()
    for o in json.loads((pl1 / "objects.json").read_text()):
        for r in fasta.iter_records(Path(o["path"])):
            t["records"] += 1
            cu = ids.contig_uid(r.seq).hex()
            if cu not in circ or cu in seen:
                continue
            seen.add(cu)
            e = ends(r.seq)
            _tally(t, r.accession, e, classify(e, long_bound), False, True)
    return _done(t)


# --------------------------------------------------------------------------------------------------------------
# deferred

DEFERRED_HEAD = ["contig", "cluster_ID", "cluster_type", "cluster_dereplication_representative", "contig_length",
                 "contig_uid"]


def deferred(gut_table: Path, skipped: Path, out: Path) -> int:
    want: dict[str, list[str]] = {}
    with open(gut_table, newline="") as fh:
        lines = (ln for ln in fh if not ln.startswith("#"))
        for r in csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE):
            if r["in_dereplicated_plasmids"] == "False" and r["circular"] == "True":
                if r["contig"] in want:
                    raise Refused(f"{gut_table}: contig {r['contig']} has two rows")
                want[r["contig"]] = [r["cluster_ID"], r["cluster_type"], r["cluster_dereplication_representative"],
                                     r["contig_length"], ""]
    for r in _tsv_rows(skipped):
        if r["source"] == GUT and r["accession"] in want:
            v = want[r["accession"]]
            if v[3] != r["length"]:
                raise Refused(f"{r['accession']}: the table says {v[3]} bases, the FASTA record {r['length']}")
            v[4] = r["contig_uid"]
    missing = sorted(c for c, v in want.items() if not v[4])
    if missing:
        raise Refused(f"{len(missing)} deferred contigs have no FASTA record in {skipped}: {', '.join(missing[:5])}")
    tmp = out.with_name(f".{out.name}.tmp")
    with open(tmp, "w") as w:
        w.write("\t".join(DEFERRED_HEAD) + "\n")
        for c in sorted(want):
            w.write("\t".join([c, *want[c]]) + "\n")
    tmp.replace(out)
    return len(want)


# --------------------------------------------------------------------------------------------------------------
# stage

STAGES = ("PL2", "PL3")  # PL2: the objects card C-40 marked redistributable; PL3: the others (DEC-63 (1), DEC-51)
RECORD_HEAD = ["ordinal", "source", "accession", "contig_uid", "length", "topology_value", "completeness_value",
               "circular", "assembly_quality", "fragment", "first_source", "first_accession", "owner_stage",
               "source_contig_uid", "source_length", "trimmed_bases"]
TRIM_HEAD = ["source", "accession", "source_contig_uid", "contig_uid", "source_length", "repeat_length",
             "bases_removed", "length"]
FLUSH = 100_000


def _bad_line(by: dict, key: str, n: int, ex: list[str]) -> str:
    return f"{key}: {n} (e.g. {', '.join(ex)})"


def plan_from_survey(survey_dir: Path, long_bound: int, fragments: str | None, ambiguous: str | None,
                     rule: str | None = None) -> Counter:
    """Per source, the trims the survey found among the records that will be staged circular. Refused before anything
    is written when a class a trim cannot settle is among those records and no rule for them (`ambiguous`) was given,
    or when a record its source flags circular is also incomplete and no rule for it (`fragments`) was given.
    Under `rule='dec65'` (DEC-65): trim only exact 20..long_bound bp repeats in METAGENOMIC_SETS; PUBLIC_SETS are
    never trimmed and stay byte-identical as deposited; incomplete fragments and non-exact metagenomic repeats are
    staged linear and untrimmed."""
    doc = json.loads((survey_dir / "summary.json").read_text())
    if doc.get("long_bound") != long_bound:
        raise Refused(f"the survey was made with --long-bound {doc.get('long_bound')}, not {long_bound}")
    trims: Counter = Counter()
    bad: dict[str, list] = {}
    frag: dict[str, list] = {}
    for r in _tsv_rows(survey_dir / "survey.tsv"):
        src, cls = r["source"], r["class"]
        if rule == "dec65":
            if src not in METAGENOMIC_SETS and src not in PUBLIC_SETS:
                raise Refused(f"source {src!r}: not in DEC-65 metagenomic or public sets")
            if src in PUBLIC_SETS or r["fragment"] == "true":
                continue
            if cls == "exact":
                trims[src] += 1
            continue
        if r["fragment"] == "true":
            v = frag.setdefault(src, [0, []])
            v[0] += 1
            if len(v[1]) < 5:
                v[1].append(r["accession"])
            if fragments != "circular":
                continue
        if cls == "exact":
            trims[src] += 1
        elif cls != "none":
            v = bad.setdefault(f"{src} {cls}", [0, []])
            v[0] += 1
            if len(v[1]) < 5:
                v[1].append(r["accession"])
    problems = []
    if frag and fragments is None:
        problems.append("records flagged circular by their source and also incomplete, with no rule given "
                        "(--fragments linear, card C-40's, or circular): "
                        + "; ".join(_bad_line(frag, k, *v) for k, v in sorted(frag.items())))
    if bad and ambiguous is None:
        problems.append("records to be staged circular whose terminal repeat a trim cannot settle, with no rule given "
                        "(--ambiguous linear, circular or asis, untrimmed in every case): "
                        + "; ".join(_bad_line(bad, k, *v) for k, v in sorted(bad.items())))
    if problems:
        raise Refused(" | ".join(problems))
    return trims


def _owners_sql(root: Path, root_stages: list[str] | None = None) -> str:
    """Each contig a stage of the root stores as its default variant, with that stage (DEC-26), as stage_est.py's
    `before` reads them (the stage's contig_index). When `root_stages` is given, read only those explicit stages."""
    from biobox import io

    base = root / "l0seq" / "contig_index"
    parts = []
    if root_stages is not None:
        for st_name in root_stages:
            d = base / f"stage={st_name}"
            files = io.part_files(d) if d.is_dir() else []
            if not files:
                raise Refused(f"root stage {st_name!r} has no contig_index part files under {root}")
            st = st_name.replace("'", "''")
            parts.append(f"SELECT contig_uid, '{st}' AS st FROM read_parquet({_lit(files)}, hive_partitioning = false) "
                         "WHERE is_default_variant")
    else:
        for d in sorted(base.glob("stage=*")) if base.is_dir() else []:
            files = io.part_files(d)
            if files:
                st = d.name.removeprefix("stage=").replace("'", "''")
                parts.append(f"SELECT contig_uid, '{st}' AS st FROM read_parquet({_lit(files)}, hive_partitioning = false) "
                             "WHERE is_default_variant")
    return " UNION ALL ".join(parts) if parts else "SELECT NULL::BLOB AS contig_uid, NULL::VARCHAR AS st LIMIT 0"


def _lit(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def _duck(work: Path):
    import duckdb

    (work / "duck").mkdir(parents=True, exist_ok=True)
    return duckdb.connect(config={"threads": 2, "memory_limit": "3GB", "temp_directory": str(work / "duck")})


def stage(inputs: Path, survey_dir: Path, root: Path, outs: dict[str, Path], finals: dict[str, Path], work: Path,
          long_bound: int, fragments: str | None, ambiguous: str | None = None,
          rule: str | None = None, root_stages: list[str] | None = None) -> dict:
    """The stage mode (module docstring): checks everything it can before writing, writes each stage into a sibling
    scratch directory, and renames both into place only when nothing was refused."""
    import shutil

    objs = json.loads((inputs / "objects.json").read_text())
    which = {Path(o["path"]).stem: ("PL2" if o["redistributable"] else "PL3") for o in objs}
    for s in STAGES:
        if s not in which.values():
            raise Refused(f"no object of {inputs} goes to {s}")
        if outs[s].exists():
            raise Refused(f"{outs[s]} exists")
    trims = plan_from_survey(survey_dir, long_bound, fragments, ambiguous, rule=rule)
    tmps = {s: outs[s].with_name(f".{outs[s].name}.tmp") for s in STAGES}
    for t in tmps.values():
        shutil.rmtree(t, ignore_errors=True)
        t.mkdir(parents=True)
    work.mkdir(parents=True, exist_ok=True)
    try:
        rep = _stage(inputs, survey_dir, root, objs, which, trims, tmps, finals, work, long_bound, fragments,
                     ambiguous, rule, root_stages)
    except BaseException:
        for t in tmps.values():
            shutil.rmtree(t, ignore_errors=True)
        raise
    for s in STAGES:
        tmps[s].rename(outs[s])
    return rep


def _stage(inputs, survey_dir, root, objs, which, trims, tmps, finals, work, long_bound, fragments, ambiguous,
           rule=None, root_stages=None) -> dict:
    import os

    import pyarrow as pa
    import pyarrow.parquet as pq

    schema = pa.schema([("ord", pa.int64()), ("source", pa.string()), ("accession", pa.string()),
                        ("uid", pa.binary(16)), ("length", pa.int64()), ("tv", pa.string()), ("cv", pa.string()),
                        ("circular", pa.bool_()), ("quality", pa.string()), ("fragment", pa.bool_()),
                        ("src_uid", pa.binary(16)), ("src_length", pa.int64()), ("trimmed", pa.int64())])
    writers = {s: pq.ParquetWriter(work / f"rows-{s}.parquet", schema) for s in STAGES}
    buf: dict[str, list] = {s: [] for s in STAGES}
    ords = {s: 0 for s in STAGES}

    def flush(s):
        if buf[s]:
            writers[s].write_table(pa.table({n: list(c) for n, c in zip(schema.names, zip(*buf[s]))}, schema=schema))
            buf[s] = []

    trim_rows: dict[str, list[tuple]] = {s: [] for s in STAGES}
    srows = _tsv_rows(survey_dir / "survey.tsv")
    pending = next(srows, None)
    per: dict[str, dict] = {}
    cur, fa_w = None, None
    try:
        for o, r, row in records_with_rows(inputs):
            stem = Path(o["path"]).stem
            s = which[stem]
            if stem != cur:
                if fa_w is not None:
                    fa_w.close()
                cur = stem
                dst = tmps[s] / f"{stem}.fasta"
                if trims[stem]:
                    fa_w = open(dst, "wb")
                else:
                    os.link(o["path"], dst)  # staged unchanged: the same bytes, no copy
                    fa_w = None
                per[stem] = {"stage": s, "records": 0, "source_bases": 0, "bases": 0, "circular_records": 0,
                             "trims": 0, "trimmed_bases": 0, "ambiguous": Counter(),
                             "object": "rewritten" if fa_w else "hard link"}
            flagged = source_flag(stem, row["topology_value"])
            frag = fragment(stem, row["completeness_value"])
            sv = None
            if flagged:
                sv = pending
                if sv is None or sv["ordinal"] != row["ordinal"] or sv["accession"] != r.accession:
                    raise Refused(f"the survey does not line up with the records at {stem} {r.accession}")
                pending = next(srows, None)
            seq, k = r.seq, 0
            if rule == "dec65":
                if stem not in METAGENOMIC_SETS and stem not in PUBLIC_SETS:
                    raise Refused(f"source {stem!r}: not in DEC-65 metagenomic or public sets")
                circ = flagged and not frag
                if circ and stem in METAGENOMIC_SETS:
                    e = ends(seq)
                    cls = classify(e, long_bound)
                    if cls != sv["class"]:
                        raise Refused(f"{stem} {r.accession}: the survey says {sv['class']}, the record is {cls}")
                    if cls == "exact":
                        k = e.exact
                        seq = trim(seq, k)
                    elif cls != "none":
                        per[stem]["ambiguous"][cls] += 1
                        circ = False
            else:
                circ = flagged and (not frag or fragments == "circular")
                if circ:
                    e = ends(seq)
                    cls = classify(e, long_bound)
                    if cls != sv["class"]:
                        raise Refused(f"{stem} {r.accession}: the survey says {sv['class']}, the record is {cls}")
                    if cls == "exact":
                        k = e.exact
                        seq = trim(seq, k)
                    elif cls != "none":
                        if ambiguous is None:
                            raise Refused(f"{stem} {r.accession}: class {cls} and no --ambiguous rule")
                        per[stem]["ambiguous"][cls] += 1
                        circ = ambiguous == "circular" or (ambiguous == "asis" and row["circular"] == "true")
            if fa_w is not None:
                fa_w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
            elif k:
                raise Refused(f"{stem} {r.accession}: a trim in an object the survey said had none")
            src_uid = bytes.fromhex(row["contig_uid"])
            uid = ids.contig_uid(seq) if k else src_uid
            if k:
                trim_rows[s].append((stem, r.accession, src_uid.hex(), uid.hex(), len(r.seq), k, k, len(seq)))
            buf[s].append((ords[s], stem, r.accession, uid, len(seq), row["topology_value"],
                           row["completeness_value"], circ, row["assembly_quality"], frag, src_uid, len(r.seq), k))
            ords[s] += 1
            if len(buf[s]) >= FLUSH:
                flush(s)
            t = per[stem]
            t["records"] += 1
            t["source_bases"] += len(r.seq)
            t["bases"] += len(seq)
            t["circular_records"] += circ
            t["trims"] += k > 0
            t["trimmed_bases"] += k
        if pending is not None:
            raise Refused(f"the survey has a row with no record: {pending['source']} {pending['accession']}")
    finally:
        if fa_w is not None:
            fa_w.close()
        for s in STAGES:
            flush(s)
            writers[s].close()
    for s in STAGES:
        trim_rows[s].sort(key=lambda t: (t[0], t[1]))
        with open(tmps[s] / "trims.tsv", "w") as w:
            w.write("\t".join(TRIM_HEAD) + "\n")
            for tr in trim_rows[s]:
                w.write("\t".join(map(str, tr)) + "\n")
    for stem, n in trims.items():
        if per.get(stem, {}).get("trims") != n:
            raise Refused(f"{stem}: the survey has {n} trims, the staging made {per.get(stem, {}).get('trims')}")
    for t in per.values():
        t["ambiguous"] = dict(sorted(t["ambiguous"].items()))
    rep = {"sources": per, "rule": rule, "fragments": fragments, "ambiguous": ambiguous, "long_bound": long_bound,
           "stages": {}}
    con = _duck(work)
    con.execute(f"CREATE TABLE own AS SELECT contig_uid AS uid, min(st) AS st, count(DISTINCT st) AS m "
                f"FROM ({_owners_sql(root, root_stages)}) GROUP BY contig_uid")
    for s in STAGES:
        rep["stages"][s] = _write_stage(con, s, work / f"rows-{s}.parquet", tmps[s], finals[s], objs, which)
    con.close()
    for s in STAGES:
        with open(inputs / "skipped_records.tsv") as fh, open(tmps[s] / "skipped_records.tsv", "w") as w:
            for i, ln in enumerate(fh):
                if i == 0 or which.get(ln.split("\t", 1)[0]) == s:
                    w.write(ln)
        (tmps[s] / "stage_report.json").write_text(json.dumps(
            {"stage": s, "rule": rule, "fragments": fragments, "ambiguous": ambiguous, "long_bound": long_bound,
             "sources": {k: v for k, v in per.items() if v["stage"] == s}, **rep["stages"][s]}, indent=1,
            sort_keys=True))
    return rep


def _write_stage(con, s: str, rows: Path, tmp: Path, final: Path, objs: list[dict], which: dict) -> dict:
    con.execute(f"CREATE OR REPLACE TABLE r AS SELECT * FROM read_parquet('{rows}')")
    con.execute("CREATE OR REPLACE TABLE f AS SELECT uid, min(ord) AS ford FROM r GROUP BY uid")
    cur = con.execute(
        "SELECT r.ord, r.source, r.accession, lower(hex(r.uid)), r.length, r.tv, r.cv, r.circular, r.quality, "
        "r.fragment, CASE WHEN r.ord = f.ford THEN '' ELSE r0.source END, "
        "CASE WHEN r.ord = f.ford THEN '' ELSE r0.accession END, coalesce(o.st, ''), lower(hex(r.src_uid)), "
        "r.src_length, r.trimmed FROM r JOIN f USING (uid) JOIN r AS r0 ON r0.ord = f.ford "
        "LEFT JOIN own AS o USING (uid) ORDER BY r.ord")
    with open(tmp / "records.tsv", "w") as w:
        w.write("\t".join(RECORD_HEAD) + "\n")
        while rows_ := cur.fetchmany(FLUSH):
            for t in rows_:
                w.write("\t".join(str(v).lower() if isinstance(v, bool) else str(v) for v in t) + "\n")
    cur = con.execute("SELECT lower(hex(r.uid)), r.circular, r.quality FROM r JOIN f ON r.ord = f.ford ORDER BY r.ord")
    n = 0
    with open(tmp / "genomes.json", "w") as w:
        w.write("[")
        while rows_ := cur.fetchmany(FLUSH):
            for uid, circ, q in rows_:
                g = {"contig_uids": [uid], "genetic_code": 11, "kind": "plasmid", "circular": [uid] if circ else [],
                     "assembly_quality": q}
                w.write((", " if n else "") + json.dumps(g))
                n += 1
        w.write("]")
    objects = [dict(o, path=str(final / f"{Path(o['path']).stem}.fasta")) for o in objs
               if which[Path(o["path"]).stem] == s]
    (tmp / "objects.json").write_text(json.dumps(objects, indent=1))
    by_source = {}
    for src, recs, bases, gen, gbases, gcirc, gcomp, new, nbases, held in con.execute(
            "SELECT r.source, count(*), sum(r.length), count(*) FILTER (WHERE r.ord = f.ford), "
            "coalesce(sum(r.length) FILTER (WHERE r.ord = f.ford), 0), "
            "count(*) FILTER (WHERE r.ord = f.ford AND r.circular), "
            "count(*) FILTER (WHERE r.ord = f.ford AND r.quality = 'complete'), "
            "count(*) FILTER (WHERE r.ord = f.ford AND o.uid IS NULL), "
            "coalesce(sum(r.length) FILTER (WHERE r.ord = f.ford AND o.uid IS NULL), 0), "
            "count(*) FILTER (WHERE o.uid IS NOT NULL) "
            "FROM r JOIN f USING (uid) LEFT JOIN own AS o USING (uid) GROUP BY r.source ORDER BY min(r.ord)").fetchall():
        by_source[src] = {"records": recs, "bases": int(bases), "genomes": gen, "genome_bases": int(gbases),
                          "genomes_circular": gcirc, "genomes_complete": gcomp, "new_contigs": new,
                          "new_bases": int(nbases), "records_owned_by_earlier_stage": held}
    owned = dict(con.execute("SELECT o.st, count(*) FROM f JOIN own AS o USING (uid) GROUP BY o.st ORDER BY o.st"
                             ).fetchall())
    twice, = con.execute("SELECT count(*) FROM f JOIN own AS o USING (uid) WHERE o.m > 1").fetchone()
    return {"genomes": n, "contigs": n, "by_source": by_source, "owned_by_stage": owned,
            "contigs_held_by_two_stages": twice,
            "new_contigs": sum(v["new_contigs"] for v in by_source.values()),
            "new_bases": sum(v["new_bases"] for v in by_source.values())}


# --------------------------------------------------------------------------------------------------------------
# overlap


def _genome_uids_sql(genomes_json: Path) -> str:
    """The distinct contig uids a genomes.json lists (what build-stage and stage_est.py's `before` read)."""
    p = str(genomes_json).replace("'", "''")
    return (f"SELECT DISTINCT from_hex(u) AS uid FROM (SELECT unnest(contig_uids) AS u FROM read_json('{p}', "
            "format = 'array', columns = {'contig_uids': 'VARCHAR[]'}))")


def overlap(stagings: dict[str, Path], root: Path, others: dict[str, Path], work: Path,
            root_stages: list[str] | None = None) -> dict:
    """What each new staging's contigs (its genomes.json: the sequences as staged) share with the stages of the root
    (default variants, as `before` counts them: one stage per contig, the least name), with the staged inputs of
    stages not built yet (their genomes.json), and with each other; beside the same counts for the records' source
    sequences (records.tsv `source_contig_uid`: the bytes before any trim)."""
    con = _duck(work)
    con.execute(f"CREATE TABLE own AS SELECT contig_uid AS uid, min(st) AS st FROM ({_owners_sql(root, root_stages)}) "
                "GROUP BY contig_uid")
    gi = sorted((root / "l0seq" / "genome_input" / "stage=PL1").glob("part-*.parquet"))
    con.execute("CREATE TABLE pl1 AS " + (f"SELECT DISTINCT contig_uid AS uid FROM read_parquet({_lit(gi)}, "
                                           "hive_partitioning = false)" if gi else "SELECT NULL::BLOB AS uid LIMIT 0"))
    for name, d in others.items():
        con.execute(f"CREATE TABLE \"o_{name}\" AS {_genome_uids_sql(d / 'genomes.json')}")
    doc: dict = {"root_stages": sorted(r[0] for r in con.execute("SELECT DISTINCT st FROM own").fetchall()),
                 "others": {k: str(v) for k, v in others.items()}, "stages": {}}
    for s, d in stagings.items():
        p = str(d / "records.tsv").replace("'", "''")
        con.execute(f"CREATE TABLE \"n_{s}\" AS {_genome_uids_sql(d / 'genomes.json')}")
        con.execute(f"CREATE TABLE \"len_{s}\" AS SELECT from_hex(contig_uid) AS uid, min(length::BIGINT) AS n "
                    f"FROM read_csv('{p}', delim = '\t', header = true, quote = '', all_varchar = true) GROUP BY uid")
        con.execute(f"CREATE TABLE \"s_{s}\" AS SELECT DISTINCT from_hex(source_contig_uid) AS uid FROM read_csv('{p}', "
                    "delim = '\t', header = true, quote = '', all_varchar = true)")
        out = {}
        for kind, t in (("as_staged", f'"n_{s}"'), ("source_bytes", f'"s_{s}"')):
            k = {"contigs": con.execute(f"SELECT count(*) FROM {t}").fetchone()[0],
                 "held_by": dict(con.execute(f"SELECT st, count(*) FROM {t} JOIN own USING (uid) GROUP BY st "
                                             "ORDER BY st").fetchall()),
                 "in_pl1_genome_input": con.execute(f"SELECT count(*) FROM {t} SEMI JOIN pl1 USING (uid)").fetchone()[0],
                 "in_staged": {name: con.execute(f'SELECT count(*) FROM {t} SEMI JOIN "o_{name}" USING (uid)').fetchone()[0]
                               for name in others}}
            out[kind] = k
        doc["stages"][s] = out
    names = list(stagings)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            doc[f"shared_{a}_{b}"] = con.execute(f'SELECT count(*) FROM "n_{a}" SEMI JOIN "n_{b}" USING (uid)').fetchone()[0]
    # who would hold each contig if the stages not built yet were built in the order given (the others, then the
    # stagings): a root stage first, else the first of them that lists it (DEC-26); what `estimate` would then count
    doc["build_order"] = [*others, *names]
    for i, s in enumerate(names):
        whens = [f"WHEN uid IN (SELECT uid FROM \"o_{x}\") THEN '{x}'" for x in others]
        whens += [f"WHEN uid IN (SELECT uid FROM \"n_{t}\") THEN '{t}'" for t in names[:i]]
        case = "CASE WHEN own.st IS NOT NULL THEN own.st " + " ".join(whens) + " END"
        q = (f'SELECT {case} AS h, l.n FROM "n_{s}" LEFT JOIN own USING (uid) '
             f'LEFT JOIN "len_{s}" l USING (uid)')
        doc["stages"][s]["held_in_build_order"] = dict(con.execute(
            f"SELECT h, count(*) FROM ({q}) WHERE h IS NOT NULL GROUP BY h ORDER BY h").fetchall())
        doc["stages"][s]["held_in_build_order_bases"] = dict(con.execute(
            f"SELECT h, coalesce(sum(n), 0)::BIGINT FROM ({q}) WHERE h IS NOT NULL GROUP BY h ORDER BY h").fetchall())
        tot_c, tot_b, new_c, new_b = con.execute(
            f"SELECT count(*), coalesce(sum(n), 0)::BIGINT, count(*) FILTER (WHERE h IS NULL), "
            f"coalesce(sum(n) FILTER (WHERE h IS NULL), 0)::BIGINT FROM ({q})").fetchone()
        doc["stages"][s]["total_contigs"] = tot_c
        doc["stages"][s]["total_bases"] = int(tot_b)
        doc["stages"][s]["new_contigs"] = new_c
        doc["stages"][s]["new_bases"] = int(new_b)
    con.close()
    return doc


# --------------------------------------------------------------------------------------------------------------


def _parse_stages(val: str | None) -> list[str] | None:
    if val is None:
        return None
    out = [s.strip() for s in val.split(",") if s.strip()]
    if not out:
        raise Refused("--root-stages must name at least one stage")
    return out


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(prog="c52_pl23.py")
    sub = ap.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("survey")
    p.add_argument("--inputs", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--long-bound", type=int, default=LONG_BOUND)
    p.add_argument("--pl1", type=Path, default=None)
    p = sub.add_parser("deferred")
    p.add_argument("--gut-table", required=True, type=Path)
    p.add_argument("--skipped", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p = sub.add_parser("stage")
    p.add_argument("--inputs", required=True, type=Path)
    p.add_argument("--survey", required=True, type=Path)
    p.add_argument("--root", required=True, type=Path)
    for s in STAGES:
        p.add_argument(f"--out-{s.lower()}", required=True, type=Path)
        p.add_argument(f"--final-{s.lower()}", type=Path, default=None)
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--long-bound", type=int, default=LONG_BOUND)
    p.add_argument("--fragments", choices=("linear", "circular"), default=None)
    p.add_argument("--ambiguous", choices=("linear", "circular", "asis"), default=None)
    p.add_argument("--rule", choices=("dec65",), default=None)
    p.add_argument("--root-stages", default=None)
    p = sub.add_parser("overlap")
    p.add_argument("--staging", action="append", required=True, help="STAGE=DIR")
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--root-stages", default=None)
    p.add_argument("--other", action="append", default=[], help="STAGE=DIR, a staged inputs directory not built yet")
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    try:
        if a.mode == "survey":
            doc = survey(a.inputs, a.out, a.long_bound, a.pl1)
            print(json.dumps({s: {"flagged": v["flagged"], "classes": v["classes"]} for s, v in doc["sources"].items()}))
        elif a.mode == "deferred":
            print(json.dumps({"deferred": deferred(a.gut_table, a.skipped, a.out), "out": str(a.out)}))
        elif a.mode == "stage":
            outs = {s: getattr(a, f"out_{s.lower()}") for s in STAGES}
            finals = {s: getattr(a, f"final_{s.lower()}") or outs[s] for s in STAGES}
            print(json.dumps(stage(a.inputs, a.survey, a.root, outs, finals, a.work, a.long_bound, a.fragments,
                                   a.ambiguous, rule=a.rule, root_stages=_parse_stages(a.root_stages)), sort_keys=True))
        elif a.mode == "overlap":
            pairs = lambda xs: {k: Path(v) for k, _, v in (x.partition("=") for x in xs)}  # noqa: E731
            doc = overlap(pairs(a.staging), a.root, pairs(a.other), a.work,
                          root_stages=_parse_stages(a.root_stages))
            a.out.write_text(json.dumps(doc, indent=1, sort_keys=True))
            print(json.dumps(doc, sort_keys=True))
    except Refused as e:
        print(f"c52_pl23.py {a.mode}: refused: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

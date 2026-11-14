"""CRISPR arrays beside a locus (card E-14; arrays are card E-12's `feature` rows of kind `crispr_array` and its
`l1_crispr` rows, one per repeat and per spacer).

Arrays are not tape slots: a gene window counts genes, and an array on the tape would shift every window (card E-12).
So a locus's arrays are found by coordinates, after `nbr.fetch` has made the window, and shown as rows of their own
beside the genes, never as slots and never with an ordinal.

The window in bases. Each locus (one anchor gene, `(contig_uid, anchor_feature_uid)`) gets a base span `[lo, hi)`
in the anchor's frame:

- a gene window (`--window Ng`): from the least start to the greatest end of the locus's rows, the anchor and the
  genes the window shows (after `--strand` and `--operon-gap`, and cut where a linear contig end cuts the window);
  on a circular contig a row reached across the origin (flag bit 1) is first moved by one contig length to the
  anchor's side (+n downstream in coordinate order, -n upstream), so the span may run below 0 or past the contig's
  length;
- a base window (`--window 20kb`, W bases): `[anchor start - W, anchor end + W)`.

An array (its `feature` row's `[start, end)`; on a circular contig also moved by -n or +n, whichever puts the most
of its bases in the span, ties to the smaller gap to the anchor, then to no move) is **within the window** when it
shares at least one base with the span, and **partial** (`partial=true`) when part of it lies outside the span. A
gene in a base window is shown when its gap to the anchor is at most W, so a gene that starts exactly W bases after
the anchor's end is shown while an array that starts there shares no base with the span and is not: the one-base
difference of a half-open span. With `--circular off` every contig is linear here too.

Per array: `kind` (`crispr_array`), `feature_uid` (the array's stable key), `producer` (the `minced-...` producer
that called it; a release may hold arrays of two settings, card E-13), `start`, `end`, `strand` (as stored: 0-based
half-open on the contig, strand +1, the strand MinCED reads the repeats on; it gives no orientation), `n_repeats`,
`n_spacers` (from `l1_crispr`), `repeat` (the commonest repeat sequence of the array, ties to the one first in the
array; MinCED's own consensus is not stored), `delta_bp` (the gap to the anchor in bases, max(starts) - min(ends) as
for genes; negative when the array overlaps the anchor), `side_rel` (+1 downstream of the anchor along its strand,
-1 upstream, 0 overlapping it), `partial`, `crosses_origin` (true when the array is reached across the origin of a
circular contig, as a gene's flag bit 1).

Reads (AGENTS §4.4): one query per table for all the loci at once, never one per locus. The loci's contigs go into a
temp table and `feature` (the `minced-` producers' files), `l1_crispr` and `contig_tapes` (the tape version's: the
contig's length and topology) are filtered by `contig_uid IN (SELECT ...)`; each file is sorted by `contig_uid`, so
each scan reads the row groups that hold those contigs (`catalog.set_join_key_filter`), whatever the number of
loci. `near(..., report=d)` puts each query's `catalog.profile_query` scans into `d["scans"]`.

A release without arrays (no `l1_crispr` file) is never read: `near` returns None and every output is as before
card E-14.
"""

from __future__ import annotations

import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

import pyarrow as pa

from biobox import ids

KIND = "crispr_array"
KIND_CODE = ids.KIND[KIND]
PRODUCER_PREFIX = "minced-"
COLUMNS = ["kind", "feature_uid", "producer", "start", "end", "strand", "n_repeats", "n_spacers", "repeat", "delta_bp",
           "side_rel", "partial", "crosses_origin"]
MARKER = "# crispr_arrays"  # the line that opens a locus table's array rows (export.locus_tables)
SUMMARY_COLUMNS = ["kind", "producer", "n_loci", "n_loci_with_array", "share_of_loci", "n_loci_partial_only",
                   "n_arrays"]
THREADS = 2
MEMORY_LIMIT = "1GB"
FLAG_CROSSES_ORIGIN = 2


def feature_files(files) -> list[tuple[Path, dict]]:
    """The release's `feature` files of CRISPR array producers (`minced-...`)."""
    return [(p, parts) for p, parts in files.entries("feature")
            if str(parts.get("producer", "")).startswith(PRODUCER_PREFIX)]


def has_arrays(files) -> bool:
    """Whether the release holds arrays: an `l1_crispr` file in its file list (manifest or release_def stages)."""
    return bool(files.entries("l1_crispr"))


def _q(s) -> str:
    return "'" + str(s).replace("'", "''") + "'"


def _connect(tmp, threads: int | None):
    import duckdb

    scratch = Path(tempfile.mkdtemp(prefix="biobox-arrays-", dir=tmp))
    con = duckdb.connect(config={"threads": threads or THREADS, "memory_limit": MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    return con, scratch


def read(files, contigs: set[bytes], tape_version: str | None, *, tmp=None, threads: int | None = None,
         report: dict | None = None) -> tuple[dict[bytes, list[dict]], dict[bytes, tuple[int, bool]]]:
    """({contig_uid: its arrays}, {contig_uid: (contig_len, is_circular)}) for `contigs`, three queries over the
    release's files (module docstring). An array: feature_uid, producer, start, end, strand, n_repeats, n_spacers,
    repeat. The contig's length and topology are those of its tape in `tape_version`."""
    from biobox import catalog
    from biobox.search.export import _view

    con, scratch = _connect(tmp, threads)
    scans: list[dict] = []
    try:
        feats = feature_files(files)
        els = files.entries("l1_crispr")
        tapes = [] if tape_version is None else files.entries("contig_tapes", tape_version=tape_version)
        _view(con, "afeat", "feature", feats, ("feature_uid", "contig_uid", "start", "end", "strand", "kind"),
              ("producer",))
        _view(con, "ael", "l1_crispr", els, ("feature_uid", "contig_uid", "element", "ordinal", "seq"), ("producer",))
        _view(con, "atape", "contig_tapes", tapes, ("contig_uid", "contig_len", "is_circular"), ())
        catalog.set_join_key_filter(con, max(catalog.row_groups_of(p for p, _ in fs) for fs in (feats, els, tapes)))
        con.register("_c_in", pa.table({"contig_uid": pa.array(sorted(contigs), pa.binary(16))}))
        con.execute("CREATE TEMP TABLE _c AS SELECT contig_uid FROM _c_in")
        con.unregister("_c_in")
        f = catalog.profile_query(con, f"""SELECT contig_uid, feature_uid, producer, start, "end", strand FROM afeat
                                           WHERE contig_uid IN (SELECT contig_uid FROM _c) AND kind = {KIND_CODE}""",
                                  profile_dir=scratch)
        e = catalog.profile_query(con, """
            WITH el AS (SELECT feature_uid, producer, element, ordinal, seq FROM ael
                        WHERE contig_uid IN (SELECT contig_uid FROM _c)),
                 n AS (SELECT feature_uid, producer, count(*) FILTER (WHERE element = 'repeat') AS n_repeats,
                              count(*) FILTER (WHERE element = 'spacer') AS n_spacers
                       FROM el GROUP BY feature_uid, producer),
                 r AS (SELECT feature_uid, producer, seq, count(*) AS k, min(ordinal) AS first_at
                       FROM el WHERE element = 'repeat' GROUP BY feature_uid, producer, seq),
                 m AS (SELECT feature_uid, producer, first(seq ORDER BY k DESC, first_at) AS repeat
                       FROM r GROUP BY feature_uid, producer)
            SELECT n.feature_uid, n.producer, n.n_repeats, n.n_spacers, m.repeat
            FROM n LEFT JOIN m USING (feature_uid, producer)""", profile_dir=scratch)
        t = catalog.profile_query(con, "SELECT contig_uid, max(contig_len), bool_or(is_circular) FROM atape "
                                       "WHERE contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid",
                                  profile_dir=scratch)
        scans = [{"table": name, **s} for name, prof in (("feature", f), ("l1_crispr", e), ("contig_tapes", t))
                 for s in prof["scans"] if s["files"]]
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
    counts = {(fu, p): (nr, ns, rep) for fu, p, nr, ns, rep in e["result"]}
    out: dict[bytes, list[dict]] = defaultdict(list)
    for cu, fu, p, s, en, strand in f["result"]:
        nr, ns, rep = counts.get((fu, p), (0, 0, None))
        out[cu].append({"feature_uid": fu, "producer": p, "start": s, "end": en, "strand": strand,
                        "n_repeats": nr, "n_spacers": ns, "repeat": rep})
    if report is not None:
        report["scans"] = scans
    return dict(out), {cu: (int(n), bool(c)) for cu, n, c in t["result"]}


def _unwrap(r: dict, n: int | None) -> tuple[int, int]:
    """A row's span in its anchor's frame: moved by one contig length when reached across the origin."""
    s, e = r["start"], r["end"]
    if r["flags"] & FLAG_CROSSES_ORIGIN and r["offset"] != 0:
        if n is None:
            raise ValueError(f"contig {r['contig_uid'].hex()}: a row across the origin and no tape length")
        k = n if r["offset"] > 0 else -n
        s, e = s + k, e + k
    return s, e


def window(rows: list[dict], policy, n: int | None) -> tuple[int, int]:
    """The base span [lo, hi) of one locus's rows (module docstring)."""
    anchor = next(r for r in rows if r["offset"] == 0)
    if policy.window_bp is not None:
        return anchor["start"] - policy.window_bp, anchor["end"] + policy.window_bp
    spans = [_unwrap(r, n) for r in rows]
    return min(s for s, _ in spans), max(e for _, e in spans)


def place(arr: dict, anchor: dict, lo: int, hi: int, n: int | None, circular: bool) -> dict | None:
    """The array's row for one locus, or None when it shares no base with [lo, hi) (module docstring)."""
    s0, e0, a_s, a_e = arr["start"], arr["end"], anchor["start"], anchor["end"]
    best = None
    for k in ((0, -n, n) if circular and n else (0,)):
        s, e = s0 + k, e0 + k
        overlap = min(e, hi) - max(s, lo)
        if overlap <= 0:
            continue
        gap = max(s, a_s) - min(e, a_e)
        rank = (-overlap, abs(gap), k != 0)
        if best is None or rank < best[0]:
            best = (rank, k, s, e, gap)
    if best is None:
        return None
    _, k, s, e, gap = best
    side = 1 if s >= a_e else (-1 if e <= a_s else 0)
    return {"kind": KIND, "feature_uid": arr["feature_uid"], "producer": arr["producer"], "start": s0, "end": e0,
            "strand": arr["strand"], "n_repeats": arr["n_repeats"], "n_spacers": arr["n_spacers"],
            "repeat": arr["repeat"], "delta_bp": gap, "side_rel": side * anchor["strand"],
            "partial": s < lo or e > hi, "crosses_origin": k != 0, "_at": s}


def near(rows: list[dict], policy, *, files, tmp=None, threads: int | None = None,
         report: dict | None = None) -> dict[tuple[bytes, bytes], dict] | None:
    """The arrays within the window of each locus of `rows` (`nbr.fetch`'s, made with `policy`), keyed by
    (contig_uid, anchor_feature_uid): {"window": (lo, hi), "arrays": [rows of COLUMNS, by start in the anchor's
    frame, then feature_uid and producer]}. None when the release holds no arrays (`has_arrays`), without reading
    anything else."""
    if not has_arrays(files):
        return None
    loci: dict[tuple[bytes, bytes], list[dict]] = defaultdict(list)
    for r in rows:
        loci[(r["contig_uid"], r["anchor_feature_uid"])].append(r)
    tvs = {r["tape_version"] for r in rows}
    if len(tvs) > 1:
        raise ValueError(f"rows of several tape versions {sorted(tvs)}: one fetch reads one tape version")
    by_contig, tapes = read(files, {cu for cu, _ in loci}, next(iter(tvs), None), tmp=tmp, threads=threads,
                            report=report)
    out = {}
    for key in sorted(loci):
        cu = key[0]
        n, circ = tapes.get(cu, (None, False))
        circ = circ and policy.circular == "auto"
        lo, hi = window(loci[key], policy, n)
        anchor = next(r for r in loci[key] if r["offset"] == 0)
        got = [p for p in (place(a, anchor, lo, hi, n, circ) for a in by_contig.get(cu, [])) if p is not None]
        got.sort(key=lambda x: (x["_at"], x["feature_uid"], x["producer"]))
        for g in got:
            del g["_at"]
        out[key] = {"window": (lo, hi), "arrays": got}
    return out


def producers(files) -> list[str]:
    """The CRISPR array producers of the release's `l1_crispr` files."""
    return sorted({str(parts["producer"]) for _, parts in files.entries("l1_crispr")})


def summary(found: dict[tuple[bytes, bytes], dict], producers_: list[str]) -> list[dict]:
    """`nbr summarize`'s array table (SUMMARY_COLUMNS): one row per producer of the release (`producers`), with the
    anchor loci, those with at least one of its arrays within the window, their share of the loci, those whose only
    such arrays are partial, and the distinct arrays (feature_uid) within some window."""
    n_loci = len(found)
    rows = []
    for p in producers_:
        loci, whole, arrays_ = set(), set(), set()
        for key, v in found.items():
            mine = [a for a in v["arrays"] if a["producer"] == p]
            if mine:
                loci.add(key)
                arrays_ |= {a["feature_uid"] for a in mine}
                if any(not a["partial"] for a in mine):
                    whole.add(key)
        rows.append({"kind": KIND, "producer": p, "n_loci": n_loci, "n_loci_with_array": len(loci),
                     "share_of_loci": round(len(loci) / n_loci, 4) if n_loci else None,
                     "n_loci_partial_only": len(loci - whole), "n_arrays": len(arrays_)})
    return rows

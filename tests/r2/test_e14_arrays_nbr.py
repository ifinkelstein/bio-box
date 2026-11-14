"""Ring R2, card E-14: CRISPR arrays beside a locus in `bio-box nbr --format tsv`, `search hmm`'s locus tables and
`nbr summarize` (`biobox.search.arrays`).

The fixture is a copy of the micro-S0 search world (`test_micro_s0_search.world`: lambda and T4, linear; M.
genitalium, circular) with CRISPR arrays planted as card E-12's producer writes them: `feature` rows of kind
`crispr_array` under a `minced-` producer and `l1_crispr` rows, one per repeat and spacer, at coordinates chosen from
the tapes (after every ninth gene; across the circular origin; just after it; near a linear end; inside a gene; one
array of 40 repeats over several genes; one array called by a second setting, `-nr2`). Small row groups, so the
tables have many. MinCED is not needed: the arrays are constructed truth (the producer is tested against MinCED in
`tests/r2/test_crispr_producer.py`).

Oracle (AGENTS §5.2): each locus's window from the S-09 coordinate oracle (`tests/r1` `_oracle_nbr`) over the tape
rows, then for every planted array a loop over its bases and the circle's shifts that counts the bases inside the
window; the array's counts and its commonest repeat from the planted lists. The gene rows of every table, and every
table of the summary, must be byte-identical to the output of the same command on the world without arrays.

Byte identity on a release without arrays: the world itself holds no `l1_crispr` file. On it the array reader is
never called (it is replaced by one that fails), and the outputs equal those of the code with the array path cut
off (`arrays.near` returning None, as before the card): this is the "unflagged call" of the card, both run in this
process. The existing exact-oracle tests of these outputs (`test_micro_s0_search.py`, W-05, W-06, W-07) run
unchanged on the same world.
"""

from __future__ import annotations

import random
import shutil
import sys
from collections import Counter
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, cli, ids, io
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l3 import tapes
from biobox.registry.schemas import conformance_problems, table
from biobox.search import arrays, export, nbr

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_search as ms  # noqa: E402
from test_micro_s0_search import world  # noqa: E402,F401  (the module-scoped search world)

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-14", "S-09", "W-05", "W-06")]

REL = ms.REL
P = "minced-0.4.2-w8-nr3-r23x47-s26x50"
P2 = "minced-0.4.2-w8-nr2-r23x47-s26x50"  # a second setting (card E-13): its arrays are rows of their own
REP, SP = 30, 35  # repeat and spacer length of the planted arrays
FEATURE_CAP, ELEMENT_CAP = 4, 8  # rows per row group of the planted files
CRISPR = ids.KIND["crispr_array"]
TV = tapes.tape_version([V])


def _seq(rng: random.Random, k: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(k))


def _array(rng, cu: bytes, n: int, circular: bool, start: int, nrep: int, producer: str = P) -> dict:
    """One planted array: its span, repeats (the first one mutated, so the commonest is not the first; four repeats
    as A B B A, a tie the earlier one wins) and spacers, and its elements as the producer stores them (an element
    wholly past a circular origin at its position on the contig)."""
    a, b = _seq(rng, REP), _seq(rng, REP)
    reps = [a, b, b, a] if nrep == 4 else [("T" if a[0] != "T" else "G") + a[1:]] + [a] * (nrep - 1)
    sps = [_seq(rng, SP) for _ in range(nrep - 1)]
    end = start + nrep * REP + (nrep - 1) * SP
    assert end - start < n and (circular or end <= n)
    fu = ids.feature_uid(cu, start, end, 1, CRISPR)
    els, at = [], start
    for k in range(nrep):
        for kind, seq in (("repeat", reps[k]),) + ((("spacer", sps[k]),) if k < nrep - 1 else ()):
            s, e = at, at + len(seq)
            if circular and s >= n:
                s, e = s - n, e - n
            els.append({"feature_uid": fu, "contig_uid": cu, "region_uid": ids.region_uid(cu, s, e, 1),
                        "element": kind, "ordinal": len(els), "start": s, "end": e, "strand": 1, "seq": seq,
                        "contig_len": n})
            at += len(seq)
    return {"contig_uid": cu, "feature_uid": fu, "producer": producer, "start": start, "end": end, "strand": 1,
            "repeats": reps, "spacers": sps, "elements": els}


def plan(tps: dict[bytes, dict]) -> list[dict]:
    rng = random.Random(20261114)
    out = []
    for cu, t in sorted(tps.items()):
        n, circ = t["contig_len"], bool(t["is_circular"])
        sl = list(zip(t["start_bps"], t["end_bps"]))
        for k in range(3, len(sl) - 2, 9):  # after every ninth gene
            out.append(_array(rng, cu, n, circ, sl[k][1] + rng.randint(5, 120), rng.randint(3, 9)))
        out.append(_array(rng, cu, n, circ, sl[20][0] + 10, 3))  # inside a gene (overlaps it)
        out.append(_array(rng, cu, n, circ, sl[30][0], 40))  # 2,565 bases over several genes
        if circ:
            out.append(_array(rng, cu, n, circ, n - 200, 5))  # across the origin: end = n + 90
            out.append(_array(rng, cu, n, circ, 40, 3))  # just after it
        else:
            out.append(_array(rng, cu, n, circ, n - 400, 4))  # near the linear end
    first = out[0]
    out.append({**_array(random.Random(7), first["contig_uid"], tps[first["contig_uid"]]["contig_len"],
                         bool(tps[first["contig_uid"]]["is_circular"]), first["start"], 2, P2)})
    return out


def plant(root: Path, planted: list[dict]) -> None:
    """The producer's files for stage S0, one part per producer and table, sorted by the registry's sort keys."""
    for p in sorted({a["producer"] for a in planted}):
        mine = [a for a in planted if a["producer"] == p]
        feats = sorted(({"feature_uid": a["feature_uid"], "contig_uid": a["contig_uid"], "start": a["start"],
                         "end": a["end"], "strand": 1, "kind": CRISPR} for a in mine),
                       key=lambda r: (r["contig_uid"], r["start"], r["end"], r["strand"]))
        els = sorted((e for a in mine for e in a["elements"]),
                     key=lambda r: (r["contig_uid"], r["feature_uid"], r["ordinal"]))
        for name, rows, cap, rel in (("feature", feats, FEATURE_CAP, f"l1/feature/producer={p}/stage=S0"),
                                     ("l1_crispr", els, ELEMENT_CAP, f"l1/crispr/producer={p}/stage=S0")):
            t = pa.Table.from_pylist(rows, schema=table(name).schema)
            assert conformance_problems(t, table(name)) == []
            io.write_table([t], table(name), root / rel / "part-0.parquet", row_group_rows=cap)


@pytest.fixture(scope="module")
def arr(world, tmp_path_factory):
    root = tmp_path_factory.mktemp("e14") / "a"
    shutil.copytree(world[0], root)
    shutil.copytree(ms._packs(world), ms._packs_of(root))
    tps = ms._tape_rows(root, TV)
    planted = plan(tps)
    plant(root, planted)
    return {"root": root, "tapes": tps, "planted": planted}


def anchors(tps: dict[bytes, dict], root: Path) -> list[str]:
    """The search tests' anchors (near the circular origin, at linear ends, beside domain genes) and the genes next
    to the special arrays: slots 3, 4, 19, 20, 30 and 31 of every contig, and the ends of every tape."""
    out = set(ms._w06_anchors(root))
    for t in tps.values():
        p = t["protein_uids"]
        out |= {p[k] for k in (0, 1, 3, 4, 19, 20, 30, 31, len(p) - 2, len(p) - 1) if p[k] is not None}
    return sorted(u.hex() for u in out)


# --- the oracle -----------------------------------------------------------------------------------------------------


def _gap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(a[0], b[0]) - min(a[1], b[1])


def oracle_window(t: dict, i: int, policy: nbr.Policy) -> tuple[int, int]:
    """[lo, hi) of anchor slot i: the span of the S-09 oracle's window (rows across the origin moved to the anchor's
    side), or the anchor plus and minus W for a window in bases."""
    slots = list(zip(t["start_bps"], t["end_bps"], t["strands"]))
    n, circ = t["contig_len"], bool(t["is_circular"])
    if policy.window_bp is not None:
        return slots[i][0] - policy.window_bp, slots[i][1] + policy.window_bp
    lo = hi = None
    for j, d, _, _, fl in ms._r1()._oracle_nbr(slots, n, circ, i, policy):
        s, e = slots[j][0], slots[j][1]
        if fl & 2:
            s, e = (s + n, e + n) if d > 0 else (s - n, e - n)
        lo, hi = (s, e) if lo is None else (min(lo, s), max(hi, e))
    return lo, hi


def oracle_rows(t: dict, i: int, lo: int, hi: int, planted: list[dict]) -> list[list[str]]:
    """The array rows of one locus, by counting each array's bases inside [lo, hi) at every shift of the circle."""
    n, circ = t["contig_len"], bool(t["is_circular"])
    a_s, a_e, strand = t["start_bps"][i], t["end_bps"][i], t["strands"][i]
    got = []
    for a in planted:
        if a["contig_uid"] != t["contig_uid"]:
            continue
        best = None
        for k in ((0, -n, n) if circ else (0,)):
            inside = sum(1 for x in range(a["start"] + k, a["end"] + k) if lo <= x < hi)
            if inside == 0:
                continue
            s, e = a["start"] + k, a["end"] + k
            rank = (-inside, abs(_gap((s, e), (a_s, a_e))), k != 0)
            if best is None or rank < best[0]:
                best = (rank, k, s, e, inside)
        if best is None:
            continue
        _, k, s, e, inside = best
        c = Counter(a["repeats"])
        top = max(c.values())
        mode = next(r for r in a["repeats"] if c[r] == top)
        side = 0 if (s < a_e and e > a_s) else (1 if s >= a_e else -1)
        got.append(((s, a["feature_uid"], a["producer"]),
                    ["crispr_array", a["feature_uid"].hex(), a["producer"], str(a["start"]), str(a["end"]), "+1",
                     str(len(a["repeats"])), str(len(a["spacers"])), mode, str(_gap((s, e), (a_s, a_e))),
                     str(side * strand), "true" if inside < e - s else "false", "true" if k else "false"]))
    return [r for _, r in sorted(got)]


def facts(table_lines: list[str]) -> dict[str, str]:
    return dict(kv.split("=", 1) for kv in table_lines[1][2:].split("\t"))


def oracle_table(arr: dict, base: list[str], policy: nbr.Policy) -> tuple[list[str], list[list[str]]]:
    """A locus table of the release with arrays: the table without arrays, then the array section."""
    f = facts(base)
    t = arr["tapes"][bytes.fromhex(f["contig_uid"])]
    i = int(f["anchor_ordinal"])
    assert t["feature_uids"][i].hex() == f["anchor_feature_uid"]
    lo, hi = oracle_window(t, i, policy)
    rows = oracle_rows(t, i, lo, hi, arr["planted"])
    return (base + [f"# crispr_arrays\twindow_start={lo}\twindow_end={hi}", "\t".join(arrays.COLUMNS)]
            + ["\t".join(r) for r in rows]), rows


def _nbr_tsv(root: Path, uids: list[str], window: str, capsys, *extra: str, err: list | None = None) -> str:
    assert cli.main(["nbr", ",".join(uids), "--root", str(root), "--release", REL, "--window", window,
                     "--format", "tsv", *extra]) == 0
    got = capsys.readouterr()
    if err is not None:
        err.append(got.err)
    return got.out


def _summary(root: Path, uids: list[str], window: str, capsys, fmt: list[str]) -> str:
    assert cli.main(["nbr", "summarize", ",".join(uids), "--root", str(root), "--release", REL, "--window", window,
                     "--by", "pfam,nr50", "--top", "0", *fmt]) == 0
    return capsys.readouterr().out


# --- locus tables ---------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("window", ["3g", "2kb"])
def test_locus_tables_gain_the_arrays_within_the_window_and_keep_their_gene_rows(world, arr, capsys, window):
    uids = anchors(arr["tapes"], arr["root"])
    policy = nbr.Policy.parse(window=window)
    base = _nbr_tsv(world[0], uids, window, capsys).rstrip("\n").split("\n\n")
    new = _nbr_tsv(arr["root"], uids, window, capsys).rstrip("\n").split("\n\n")
    assert len(base) == len(new) > 20
    seen = []
    for b, g in zip(base, new):
        want, rows = oracle_table(arr, b.split("\n"), policy)
        assert g.split("\n") == want
        seen += rows
    cols = {c: k for k, c in enumerate(arrays.COLUMNS)}
    # the fixture reaches every case the definition names
    assert any(r[cols["partial"]] == "true" for r in seen) and any(r[cols["partial"]] == "false" for r in seen)
    assert any(r[cols["crosses_origin"]] == "true" for r in seen)
    assert any(r[cols["side_rel"]] == "0" for r in seen) and {"1", "-1"} <= {r[cols["side_rel"]] for r in seen}
    assert any(int(r[cols["delta_bp"]]) < 0 for r in seen) and any(r[cols["producer"]] == P2 for r in seen)
    assert any(len(g.split("\n")) == len(b.split("\n")) + 2 for b, g in zip(base, new))  # a locus with none


def test_out_dir_writes_the_same_tables_and_counts_array_rows(arr, capsys, tmp_path):
    uids = anchors(arr["tapes"], arr["root"])
    tables = _nbr_tsv(arr["root"], uids, "3g", capsys).rstrip("\n").split("\n\n")
    d, err = tmp_path / "loci", []
    _nbr_tsv(arr["root"], uids, "3g", capsys, "--out", str(d), err=err)
    files = sorted(p for p in d.iterdir())
    assert len(files) == len(tables)
    assert sorted(p.read_text() for p in files) == sorted(t + "\n" for t in tables)
    # card W-23's count: gene rows and array rows, not the marker and column lines of the array section
    rows = sum(len(t.split("\n")) - 5 for t in tables)
    assert all(arrays.MARKER in t for t in tables)
    assert f"wrote {d}: {len(tables)} tables, {rows} rows" in err[0], err


@pytest.mark.parametrize("window", ["3g", "2kb"])
def test_summary_gains_an_array_table_and_keeps_every_other_line(world, arr, capsys, window):
    uids = anchors(arr["tapes"], arr["root"])
    policy = nbr.Policy.parse(window=window)
    base = _summary(world[0], uids, window, capsys, ["--format", "tsv"])
    new = _summary(arr["root"], uids, window, capsys, ["--format", "tsv"])
    loci = [oracle_table(arr, b.split("\n"), policy)[1]
            for b in _nbr_tsv(world[0], uids, window, capsys).rstrip("\n").split("\n\n")]
    lines = ["", "\t".join(arrays.SUMMARY_COLUMNS)]
    pcol, partcol, fucol = arrays.COLUMNS.index("producer"), arrays.COLUMNS.index("partial"), 1
    for p in sorted([P, P2]):
        with_ = [rs for rs in loci if any(r[pcol] == p for r in rs)]
        partial_only = [rs for rs in with_ if all(r[partcol] == "true" for r in rs if r[pcol] == p)]
        n_arrays = len({r[fucol] for rs in loci for r in rs if r[pcol] == p})
        lines.append("\t".join(["crispr_array", p, str(len(loci)), str(len(with_)), str(round(len(with_) / len(loci), 4)),
                                str(len(partial_only)), str(n_arrays)]))
    assert new == base + "\n".join(lines) + "\n"
    assert int(lines[-1].split("\t")[3]) > 0
    # the JSON records are those of the release without arrays (the locked record schemas have no array key)
    assert _summary(arr["root"], uids, window, capsys, ["--json"]) == _summary(world[0], uids, window, capsys, ["--json"])


def test_a_release_manifest_gives_the_same_tables_as_the_release_def(arr, capsys, tmp_path):
    """Readers take the manifest's files when there is one (AGENTS §3.9): the arrays are found through it too."""
    root = tmp_path / "m"
    shutil.copytree(arr["root"], root)
    man = catalog.build_manifest(root, REL, ["S0"])
    assert {f["table"] for f in man["files"]} >= {"l1_crispr", "feature"}
    catalog.write_manifest(man, catalog.manifest_path(root, REL))
    assert catalog.ReleaseFiles.of(root, REL).source == "manifest"
    uids = anchors(arr["tapes"], arr["root"])
    assert _nbr_tsv(root, uids, "3g", capsys) == _nbr_tsv(arr["root"], uids, "3g", capsys)


def test_nbr_json_and_counts_are_unchanged(world, arr, capsys):
    uids = anchors(arr["tapes"], arr["root"])
    for fmt in (["--json"], []):
        outs = []
        for root in (world[0], arr["root"]):
            assert cli.main(["nbr", "--anchors", ",".join(uids), "--root", str(root), "--release", REL,
                             "--window", "3g", *fmt]) == 0
            outs.append(capsys.readouterr().out)
        assert outs[0] == outs[1] and outs[0]


# --- reads by contig (AGENTS §4.4) ----------------------------------------------------------------------------------


def _holding(files: list[Path], contigs: set[bytes]) -> tuple[int, int]:
    """(row groups whose contig_uid range holds one of `contigs`, all row groups), from the footers."""
    hold = total = 0
    for f in files:
        md = pq.ParquetFile(f).metadata
        k = md.schema.to_arrow_schema().get_field_index("contig_uid")
        for g in range(md.num_row_groups):
            st = md.row_group(g).column(k).statistics
            total += 1
            hold += any(st.min <= c <= st.max for c in contigs)
    return hold, total


def test_arrays_are_read_by_contig_whatever_the_number_of_loci(arr, capsys):
    root, tps = arr["root"], arr["tapes"]
    files = catalog.ReleaseFiles.of(root, REL)
    lam = min(cu for cu, t in tps.items() if not t["is_circular"] and t["contig_len"] < 60_000)
    one = {tps[lam]["protein_uids"][10]}
    many = {u for u in tps[lam]["protein_uids"] if u is not None}
    assert len(many) > 30
    scans, contigs, n_loci = {}, {}, {}
    for name, uids in (("one", one), ("many", many), ("all", {u for t in tps.values() for u in t["protein_uids"] if u})):
        rows, _ = nbr.fetch(root, uids, nbr.Policy.parse(window="3g"), files=files)
        report = {}
        found = arrays.near(rows, nbr.Policy.parse(window="3g"), files=files, report=report)
        n_loci[name] = len(found)
        assert n_loci[name] == len({(r["contig_uid"], r["anchor_feature_uid"]) for r in rows})
        contigs[name] = {r["contig_uid"] for r in rows}
        per: dict[str, list[int]] = {}
        for s in report["scans"]:  # one scan per producer's files; summed per table
            got = per.setdefault(s["table"], [0, 0])
            got[0] += s["row_groups_scanned"]
            got[1] += s["row_groups_total"]
        scans[name] = {k: tuple(v) for k, v in per.items()}
    assert contigs["one"] == contigs["many"] == {lam} and n_loci["one"] == 1 and n_loci["many"] > 30
    for tname in ("feature", "l1_crispr"):
        paths = [p for p in files.table(tname) if "producer=minced" in str(p)]
        hold, total = _holding(paths, {lam})
        assert total >= 12, (tname, total)
        assert scans["one"][tname] == (hold, total), (tname, scans)
        assert scans["many"][tname] == scans["one"][tname]  # one read for all loci, not one per locus
        assert hold < total and scans["all"][tname] == (_holding(paths, contigs["all"])[0], total)
    print("e14 scans", scans, "loci", n_loci)  # shown with pytest -rA


# --- no arrays: byte identity and nothing read ----------------------------------------------------------------------


def _outputs(root: Path, uids: list[str], capsys, tmp_path: Path, tag: str) -> dict[str, str]:
    out = {}
    for window in ("3g", "2kb"):
        out[f"tsv {window}"] = _nbr_tsv(root, uids, window, capsys)
        d, err = tmp_path / f"{tag}-{window}", []
        _nbr_tsv(root, uids, window, capsys, "--out", str(d), err=err)
        out[f"err {window}"] = err[0].replace(str(d), "DIR")
        out[f"dir {window}"] = "".join(p.name + "\n" + p.read_text() for p in sorted(d.iterdir()))
        out[f"sum {window}"] = _summary(root, uids, window, capsys, ["--format", "tsv"])
        out[f"json {window}"] = _summary(root, uids, window, capsys, ["--json"])
    return out


def test_on_a_release_without_arrays_nothing_is_read_and_every_output_is_the_unflagged_one(world, capsys, tmp_path,
                                                                                          monkeypatch):
    root = world[0]
    files = catalog.ReleaseFiles.of(root, REL)
    assert files.entries("l1_crispr") == [] and not arrays.has_arrays(files)
    uids = anchors(ms._tape_rows(root, TV), root)

    def no_read(*a, **k):
        raise AssertionError("a release without arrays must not be read for arrays")

    monkeypatch.setattr(arrays, "read", no_read)
    got = _outputs(root, uids, capsys, tmp_path, "new")
    assert not any(arrays.MARKER in v for v in got.values()) and "crispr_array" not in got["sum 3g"]
    monkeypatch.setattr(arrays, "near", lambda *a, **k: None)  # the call cut off: the code before card E-14
    assert _outputs(root, uids, capsys, tmp_path, "old") == got


# --- search hmm -----------------------------------------------------------------------------------------------------


def test_search_hmm_writes_the_array_rows_and_says_so_in_its_readme(world, arr, capsys, tmp_path):
    root = arr["root"]
    hmm, pos = ms._w07_hmm(world, tmp_path)
    spec = ms._w07_spec(root, pos)
    out, runs, work = tmp_path / "out", tmp_path / "runs", tmp_path / "work"
    code = cli.main(["search", "hmm", "--root", str(root), "--release", REL, "--hmm", str(hmm), "--bits", "25",
                     "--neighbors", "5", "--positive", spec, "--out", str(out), "--cpus", "1",
                     "--pack-dir", str(ms._packs_of(root)), "--runs", str(runs), "--work", str(work)])
    got = capsys.readouterr()
    assert code == 0, got.err
    rid = got.out.splitlines()[0].split("\t")[0]
    assert (out / "README.txt").read_text() == export.README_TEXT + export.README_ARRAYS
    loci = sorted((out / "loci").glob("locus_*.tsv"))
    assert loci and all(arrays.MARKER in p.read_text() for p in loci)
    assert (out / "neighbors.tsv").read_text().rstrip("\n").split("\n\n")[-1].startswith("\t".join(arrays.SUMMARY_COLUMNS))
    # the separate command on the same run writes the same tables
    sep = tmp_path / "sep"
    assert cli.main(["nbr", rid, "--root", str(root), "--pack-dir", str(ms._packs_of(root)), "--runs", str(runs),
                     "--window", "5g", "--format", "tsv", "--out", str(sep)]) == 0
    capsys.readouterr()
    assert {p.name: p.read_text() for p in loci} == {p.name: p.read_text() for p in sep.glob("locus_*.tsv")}
    # and the oracle's arrays for each locus
    policy = nbr.Policy.parse(window="5g")
    for p in loci:
        lines = p.read_text().rstrip("\n").split("\n")
        k = next(i for i, line in enumerate(lines) if line.startswith(arrays.MARKER))
        want, _ = oracle_table(arr, lines[:k], policy)
        assert lines == want


def test_search_hmm_readme_on_a_release_without_arrays_is_the_old_one(world, capsys, tmp_path):
    root = world[0]
    hmm, pos = ms._w07_hmm(world, tmp_path)
    out = tmp_path / "out"
    code = cli.main(["search", "hmm", "--root", str(root), "--release", REL, "--hmm", str(hmm), "--bits", "25",
                     "--neighbors", "5", "--positive", ms._w07_spec(root, pos), "--out", str(out), "--cpus", "1",
                     "--pack-dir", str(ms._packs(world)), "--runs", str(tmp_path / "runs"),
                     "--work", str(tmp_path / "work")])
    assert code == 0, capsys.readouterr().err
    assert (out / "README.txt").read_text() == export.README_TEXT
    assert not any(arrays.MARKER in p.read_text() for p in (out / "loci").glob("*.tsv"))
    assert "crispr_array" not in (out / "neighbors.tsv").read_text()

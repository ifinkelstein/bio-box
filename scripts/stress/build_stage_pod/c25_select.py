#!/usr/bin/env python3
"""S2C selection for card C-25 (host python3, stdlib only): at most CAP complete genomes per GTDB species.

    python3 -I c25_select.py GTDB_DIR PHASE1_ALL_TSV OUT_DIR [--cap 20] [--phase1-not-served ACC,ACC,...]
                             [--candidates] [--ani PAIRS_TSV [--min-af 15] [--require-pairs]]

Rule without --ani, by quality order only:
GTDB R226 bac120 + ar53 metadata, read by r50_select.complete_genomes (phase 1's reader); a genome is complete when
ncbi_assembly_level = 'Complete Genome' (no other field is used); within each gtdb_taxonomy the genomes are ranked
by checkm2_completeness - 5 * checkm2_contamination descending, a value that does not parse last, every tie by
GTDB accession (with its RS_/GB_ prefix; this orders as the bare accession does, since GB_ < RS_ and GCA < GCF);
ranks 1 to CAP are selected. Rank 1 is phase 1's rule (r50_select.py), and the whole rule is the one that gave the
plan's 25,458 (scripts/stress/gtdb_counts.py, S2C_complete_le20_per_species).

PHASE1_ALL_TSV is phase 1's all.tsv (11,973 genomes); --phase1-not-served lists the phase-1 accessions NCBI did not
serve. Writes, sorted by accession, OUT_DIR/s2c_selected.tsv and OUT_DIR/remainder.tsv (selected minus phase 1),
with phase 1's columns then species_rank, species_complete_genomes, quality_score, checkm2_completeness,
checkm2_contamination, gtdb_representative; OUT_DIR/by_phylum.tsv; OUT_DIR/selection.json (counts, the phase-1
comparison, what an ANI step could change, file digests); prints selection.json.

--candidates (card C-27) also writes OUT_DIR/ani_candidates.tsv, every complete genome of the species with more
than CAP, and OUT_DIR/ani_fetch.tsv, those not in PHASE1_ALL_TSV plus the --phase1-not-served genomes (retried
once, DEC-31), with a fetch_reason column; both sorted by accession, in the columns above.

--ani PAIRS_TSV (card C-27, DEC-31): in each species with more than CAP complete genomes the CAP genomes are
chosen by greedy farthest-first from the best-quality genome (rank 1, so phase 1 stays a subset; asserted). The
distance of a pair is 100 - ANI. A pair the table does not report counts as the largest distance, 100: a row
with an empty ani (skani was asked and printed nothing, below its screen), a row whose larger alignment fraction
is under --min-af (default 15, skani's own default for printing a value: an ANI over less than 15% of both
genomes is not a genome-wide value), or a pair absent from the table (--require-pairs makes that an error, for
a table that should hold every pair the pick consults). After each pick every unpicked genome keeps its distance
to the nearest pick; the next pick is the genome farthest from the picks, ties by quality order (which ends in
the accession). Species with at most CAP keep all their genomes. PAIRS_TSV is the normalised table that
c27_skani.py writes: columns genome_a, genome_b, ani, af_a, af_b (bare accessions; af in percent; either order;
a pair given twice keeps its larger ANI). The selection then adds columns ani_pick (1..CAP in the species over
the cap) and ani_to_nearest_earlier_pick (the ANI to the nearest earlier pick when it was picked, 'unreported'
for the largest distance), and OUT_DIR/ani_species.tsv: per species over the cap, how many picks differ from the
quality order, and the minimum pairwise ANI among the CAP picks under each rule with the pairs it rests on.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I leaves the script's directory off the path
from r50_select import COLS, FILES, cell, complete_genomes, summary  # noqa: E402

EXTRA = ["species_rank", "species_complete_genomes", "quality_score", "checkm2_completeness",
         "checkm2_contamination", "gtdb_representative"]
BINS = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 20), (21, 50), (51, 100), (101, None)]
PHYLUM_COLS = ["phylum", "species", "selected_genomes", "selected_bp", "remainder_genomes", "remainder_bp",
               "remainder_proteins"]


def rank_species(items: list[tuple]) -> dict[str, list[dict]]:
    """(key, row) pairs -> species -> rows in quality order, each annotated with species_rank and the species size."""
    by_sp: dict[str, list[tuple]] = {}
    for key, row in items:
        by_sp.setdefault(row["gtdb_taxonomy"], []).append((key, row))
    out = {}
    for sp in sorted(by_sp):
        ranked = [row for _, row in sorted(by_sp[sp], key=lambda kr: kr[0])]
        for i, row in enumerate(ranked, 1):
            row["species_rank"], row["species_complete_genomes"] = i, len(ranked)
        out[sp] = ranked
    return out


def fmt(row: dict, col: str) -> str:
    v = row[col]
    if col == "gtdb_representative":
        return "t" if v else "f"
    return cell(v)


def tsv_text(rows: list[dict], cols: list[str]) -> str:
    return "\t".join(cols) + "\n" + "".join("\t".join(fmt(r, c) for c in cols) + "\n" for r in rows)


def put(path: Path, text: str) -> str:
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return hashlib.sha256(text.encode()).hexdigest()


def same_score(a: dict, b: dict) -> bool:
    return a["quality_score"] == b["quality_score"]  # None == None: both unparsed, ordered by accession alone


def bin_label(n: int) -> str:
    for lo, hi in BINS:
        if n >= lo and (hi is None or n <= hi):
            return str(lo) if lo == hi else (f"{lo}+" if hi is None else f"{lo}-{hi}")
    raise ValueError(n)


# --- card C-27: farthest-first on ANI ---------------------------------------------------------------------------

ANI_COLS = ["genome_a", "genome_b", "ani", "af_a", "af_b"]
EXTRA_ANI = ["ani_pick", "ani_to_nearest_earlier_pick"]
SPECIES_COLS = ["species", "complete_genomes", "differ"] + [
    f"{rule}_{c}" for rule in ("ani", "quality")
    for c in ("min_pairwise_ani", "pairs_reported", "pairs_unreported", "pairs_absent")]
UNREPORTED = 0.0  # the ANI of a pair that is not reported: distance 100 - 0 = 100, the largest


class TableError(SystemExit):
    pass


def _pct(s: str, what: str, where) -> float:
    try:
        v = float(s)
    except ValueError:
        v = math.nan
    if not (0.0 < v <= 100.0):
        raise TableError(f"{where}: {what} {s!r} is not in (0, 100]")
    return v


def ani_value(ani: str, af_a: str, af_b: str, min_af: float, where: str = "") -> float | None:
    """One row's ANI as the pick uses it: None (not reported) for an empty ani or a larger alignment fraction under
    min_af. Shared with c27_skani.py, so that the driver's rounds and this script read a row the same way."""
    if ani == "":
        return None
    v = _pct(ani, "ani", where)
    return None if max(float(af_a), float(af_b)) < min_af else v


def read_ani_table(path: Path, min_af: float) -> tuple[dict[tuple[str, str], float | None], dict]:
    """The normalised table -> {(a, b) with a < b: ANI, or None when not reported}, and counts.

    Not reported: an empty ani, or a larger alignment fraction under min_af. A pair given twice keeps its larger
    ANI (a reported value beats an unreported one), so the result does not depend on row order."""
    pairs: dict[tuple[str, str], float | None] = {}
    st = {"rows": 0, "rows_empty_ani": 0, "rows_below_min_af": 0, "pairs_given_twice": 0}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        miss = [c for c in ANI_COLS if c not in (rd.fieldnames or [])]
        if miss:
            raise TableError(f"{path}: missing columns {miss}")
        for r in rd:
            st["rows"] += 1
            a, b = r["genome_a"], r["genome_b"]
            if a == b:
                raise TableError(f"{path}: a genome paired with itself: {a}")
            ani = ani_value(r["ani"], r["af_a"], r["af_b"], min_af, str(path))
            if r["ani"] == "":
                st["rows_empty_ani"] += 1
            elif ani is None:
                st["rows_below_min_af"] += 1
            key = (a, b) if a < b else (b, a)
            if key in pairs:
                st["pairs_given_twice"] += 1
                old = pairs[key]
                ani = old if ani is None else (ani if old is None else max(old, ani))
            pairs[key] = ani
    return pairs, st


def farthest_first(ranked: list[str], cap: int, round_fn) -> list[tuple[str, float | None]]:
    """Greedy farthest-first pick of `cap` genomes from `ranked` (quality order), starting from ranked[0].

    round_fn(newest_pick, rest) -> {genome: ANI or None} for the genomes in `rest` (quality order); None counts as
    UNREPORTED. Each unpicked genome keeps the largest ANI to any pick so far (its distance 100 - ANI to the
    nearest pick); the next pick has the smallest such ANI (the largest distance), ties by quality order. Returns
    [(accession, ANI to the nearest earlier pick when it was picked; None for the first)]. Shared with
    c27_skani.py, whose round_fn runs skani for the newest pick against the rest."""
    idx = {g: i for i, g in enumerate(ranked)}
    picks: list[tuple[str, float | None]] = [(ranked[0], None)]
    near = dict.fromkeys(ranked[1:], -1.0)
    while len(picks) < min(cap, len(ranked)):
        rest = sorted(near, key=idx.__getitem__)
        sims = round_fn(picks[-1][0], rest)
        for g in rest:
            s = sims.get(g)
            s = UNREPORTED if s is None else s
            if s > near[g]:
                near[g] = s
        best = min(rest, key=lambda g: (near[g], idx[g]))
        picks.append((best, near.pop(best)))
    return picks


def fmt_ani(v: float | None) -> str:
    if v is None:
        return ""
    return "unreported" if v == UNREPORTED else repr(v)


def pairwise_min(accs: list[str], pairs: dict) -> tuple[float | None, int, int, int]:
    """(minimum reported ANI, pairs reported, pairs unreported, pairs absent) over all pairs of `accs`."""
    lo, rep, unrep, absent = None, 0, 0, 0
    for i, a in enumerate(accs):
        for b in accs[i + 1:]:
            key = (a, b) if a < b else (b, a)
            if key not in pairs:
                absent += 1
            elif pairs[key] is None:
                unrep += 1
            else:
                rep += 1
                lo = pairs[key] if lo is None else min(lo, pairs[key])
    return lo, rep, unrep, absent


def quantiles(vals: list[float]) -> dict:
    """Nearest-rank quantiles (p0, 10, 25, 50, 75, 90, 100) of a list; {} when empty."""
    v = sorted(vals)
    if not v:
        return {}
    return {f"p{p}": v[max(0, math.ceil(p / 100 * len(v)) - 1)] for p in (0, 10, 25, 50, 75, 90, 100)} | {"n": len(v)}


def ani_select(a, ranked: dict, by_acc: dict, cap: int, p1_set: set) -> tuple[dict, dict, list[dict]]:
    """The --ani selection: species -> selected rows, the selection.json block, the ani_species.tsv rows."""
    pairs, st = read_ani_table(a.ani, a.min_af)
    sp_of = {acc: r["gtdb_taxonomy"] for acc, r in by_acc.items()}
    unused = 0
    for x, y in sorted(pairs):
        for g in (x, y):
            if g not in sp_of:
                raise TableError(f"{a.ani}: {g} is not a complete genome in the metadata")
        if sp_of[x] != sp_of[y]:
            raise TableError(f"{a.ani}: {x} and {y} are in two species")
        unused += len(ranked[sp_of[x]]) <= cap
    n = {"consulted": 0, "unreported": 0, "absent": 0}
    absent_eg: list[tuple[str, str]] = []

    def round_fn(newest: str, rest: list[str]) -> dict:
        got = {}
        for g in rest:
            key = (newest, g) if newest < g else (g, newest)
            n["consulted"] += 1
            if key not in pairs:
                n["absent"] += 1
                if len(absent_eg) < 5:
                    absent_eg.append(key)
                got[g] = None
            else:
                got[g] = pairs[key]
                n["unreported"] += pairs[key] is None
        return got

    for rows in ranked.values():
        for r in rows:
            r["ani_pick"] = r["ani_to_nearest_earlier_pick"] = ""
    sel, species_rows = {}, []
    mins: dict[str, list[float]] = {"ani": [], "quality": []}
    incomplete = {"ani": 0, "quality": 0}
    for sp, rows in ranked.items():
        if len(rows) <= cap:
            sel[sp] = rows
            continue
        picks = farthest_first([r["accession"] for r in rows], cap, round_fn)
        for i, (acc, s) in enumerate(picks, 1):
            by_acc[acc]["ani_pick"], by_acc[acc]["ani_to_nearest_earlier_pick"] = str(i), fmt_ani(s)
        sel[sp] = [by_acc[acc] for acc, _ in picks]
        f = [acc for acc, _ in picks]
        q = [r["accession"] for r in rows[:cap]]
        row = {"species": sp, "complete_genomes": len(rows), "differ": len(set(f) - set(q))}
        for rule, accs in (("ani", f), ("quality", q)):
            lo, rep, unrep, ab = pairwise_min(accs, pairs)
            row |= {f"{rule}_min_pairwise_ani": "" if lo is None else repr(lo), f"{rule}_pairs_reported": rep,
                    f"{rule}_pairs_unreported": unrep, f"{rule}_pairs_absent": ab}
            if lo is not None and unrep == 0 and ab == 0:
                mins[rule].append(lo)
            else:
                incomplete[rule] += 1
        species_rows.append(row)
    if n["absent"] and a.require_pairs:
        raise TableError(f"{a.ani}: {n['absent']} pairs the pick consults are not in the ANI table "
                         f"(--require-pairs), for example {absent_eg}")
    a_sel = {r["accession"] for rows in sel.values() for r in rows}
    q_p1 = {r["accession"] for rows in ranked.values() for r in rows[:cap]} & p1_set
    assert all(rows[0]["accession"] in a_sel for rows in ranked.values()), "a rank-1 genome was not picked"
    assert q_p1 <= a_sel, f"phase-1 genomes dropped by the ANI rule: {sorted(q_p1 - a_sel)}"
    top = sorted(species_rows, key=lambda r: (-r["complete_genomes"], r["species"]))[:10]
    block = {
        "rule": f"species with more than cap: greedy farthest-first from rank 1 on distance 100 - ANI, ties by "
                f"quality order; a pair not reported (empty ani, larger alignment fraction under {a.min_af:g}, or "
                f"absent) counts as distance 100; species with at most cap: all",
        "table": a.ani.name, "table_sha256": hashlib.sha256(a.ani.read_bytes()).hexdigest(), "min_af": a.min_af,
        "table_counts": st, "table_pairs": len(pairs), "table_pairs_in_species_at_or_under_cap": unused,
        "species_over_cap": len(species_rows), "consulted_pairs": n["consulted"],
        "consulted_unreported": n["unreported"], "consulted_absent": n["absent"],
        "require_pairs": a.require_pairs,
        "selected_differ_from_quality": sum(r["differ"] for r in species_rows),
        "phase1_subset": True,
        "largest_species": [{"species": r["species"].split(";")[-1], "complete_genomes": r["complete_genomes"],
                             "differ": r["differ"]} for r in top],
        "min_pairwise_ani": {rule: {"species_all_pairs_reported": quantiles(mins[rule]),
                                    "species_with_a_pair_unreported_or_absent": incomplete[rule]}
                             for rule in ("ani", "quality")},
    }
    return sel, block, species_rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("gtdb", type=Path)
    ap.add_argument("phase1", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--cap", type=int, default=20)
    ap.add_argument("--phase1-not-served", default="")
    ap.add_argument("--candidates", action="store_true", help="also write ani_candidates.tsv and ani_fetch.tsv")
    ap.add_argument("--ani", type=Path, help="normalised pairwise ANI table: pick for ANI diversity (card C-27)")
    ap.add_argument("--min-af", type=float, default=15.0,
                    help="a pair whose larger alignment fraction (percent) is under this is not reported")
    ap.add_argument("--require-pairs", action="store_true",
                    help="refuse when a pair the pick consults is absent from the ANI table")
    a = ap.parse_args()
    cap = a.cap
    a.out.mkdir(parents=True, exist_ok=True)

    p1_bytes = a.phase1.read_bytes()
    with open(a.phase1, newline="") as fh:
        p1 = [r["accession"] for r in csv.DictReader(fh, delimiter="\t")]
    p1_set = set(p1)
    not_served = sorted(x for x in a.phase1_not_served.split(",") if x)

    sha: dict[str, str] = {}
    levels: dict[str, str] = {}
    items = list(complete_genomes(a.gtdb, sha, watch=p1_set | set(not_served), levels=levels))
    header_ani = []
    for name in FILES:
        with gzip.open(a.gtdb / name, "rt") as fh:
            header_ani += [c for c in fh.readline().rstrip("\n").split("\t") if "ani" in c.lower().split("_")]
    unparsed = sum(r["quality_score"] is None for _, r in items)
    nonfinite = sum(r["quality_score"] is not None and not math.isfinite(r["quality_score"]) for _, r in items)

    ranked = rank_species(items)
    by_acc = {r["accession"]: r for rows in ranked.values() for r in rows}
    assert len(by_acc) == len(items), "an accession appears twice among the complete genomes"
    q_selected = sorted((r for rows in ranked.values() for r in rows[:cap]), key=lambda r: r["accession"])
    sel_by_sp = {sp: rows[:cap] for sp, rows in ranked.items()}
    ani_out = None
    if a.ani is not None:
        sel_by_sp, ani_out, species_rows = ani_select(a, ranked, by_acc, cap, p1_set)
    selected = sorted((r for rows in sel_by_sp.values() for r in rows), key=lambda r: r["accession"])
    rank1 = sorted((rows[0] for rows in ranked.values()), key=lambda r: r["accession"])
    remainder = [r for r in selected if r["accession"] not in p1_set]

    # phase 1 against the set
    missing = []
    sel_set = {r["accession"] for r in selected}
    for acc in sorted(p1_set):
        r = by_acc.get(acc)
        if acc in sel_set:
            continue
        if acc not in levels:
            why = "not in the metadata"
        elif levels[acc] != "Complete Genome":
            why = f"not complete in the metadata (ncbi_assembly_level {levels[acc]})"
        elif r["species_rank"] <= cap:
            why = "not picked by the ANI rule"  # cannot happen: ani_select asserts it
        else:
            why = f"rank {r['species_rank']} of {r['species_complete_genomes']} complete genomes of its species, " \
                  "beyond the cap"
        missing.append({"accession": acc, "reason": why})
    ns = [{"accession": acc, "species_rank": by_acc[acc]["species_rank"] if acc in by_acc else None,
           "species_complete_genomes": by_acc[acc]["species_complete_genomes"] if acc in by_acc else None}
          for acc in not_served]

    # what an ANI-diversity step could change: only species with more than cap complete genomes
    over = {sp: rows for sp, rows in ranked.items() if len(rows) > cap}
    lo = hi = 0
    tie_cut = cut_top = decided_by_acc = 0
    for sp, rows in ranked.items():
        best, rest = rows[0]["genome_size"], sorted(r["genome_size"] for r in rows[1:])
        k = min(cap - 1, len(rest))
        lo += best + sum(rest[:k])
        hi += best + sum(rest[len(rest) - k:])
        if sp in over:
            if same_score(rows[cap - 1], rows[cap]):
                tie_cut += 1
                decided_by_acc += sum(same_score(r, rows[cap - 1]) for r in rows[:cap])
            cut_top += same_score(rows[0], rows[cap])
    p1_bp = sum(r["genome_size"] for r in rank1)
    rem_over = [r for r in remainder if r["species_complete_genomes"] > cap]

    bins: dict[str, dict] = {}
    for lo_, hi_ in BINS:
        bins[bin_label(lo_)] = {"bin": bin_label(lo_), "species": 0, "selected_genomes": 0, "selected_bp": 0,
                                "remainder_genomes": 0, "remainder_bp": 0}
    for sp, rows in ranked.items():
        b = bins[bin_label(len(rows))]
        b["species"] += 1
        for r in sel_by_sp[sp]:
            b["selected_genomes"] += 1
            b["selected_bp"] += r["genome_size"]
            if r["accession"] not in p1_set:
                b["remainder_genomes"] += 1
                b["remainder_bp"] += r["genome_size"]

    phyla: dict[str, dict] = {}
    for sp, rows in ranked.items():
        ph = sp.split(";")[1] if ";" in sp else sp
        p = phyla.setdefault(ph, dict.fromkeys(PHYLUM_COLS, 0) | {"phylum": ph})
        p["species"] += 1
        for r in sel_by_sp[sp]:
            p["selected_genomes"] += 1
            p["selected_bp"] += r["genome_size"]
            if r["accession"] not in p1_set:
                p["remainder_genomes"] += 1
                p["remainder_bp"] += r["genome_size"]
                p["remainder_proteins"] += r["protein_count"] or 0
    phylum_rows = sorted(phyla.values(), key=lambda p: (-p["remainder_genomes"], -p["selected_genomes"], p["phylum"]))

    cols = COLS + EXTRA + (EXTRA_ANI if ani_out is not None else [])
    digests = {"s2c_selected.tsv": put(a.out / "s2c_selected.tsv", tsv_text(selected, cols)),
               "remainder.tsv": put(a.out / "remainder.tsv", tsv_text(remainder, cols)),
               "by_phylum.tsv": put(a.out / "by_phylum.tsv", tsv_text(phylum_rows, PHYLUM_COLS))}
    rule = "ncbi_assembly_level = 'Complete Genome'; per gtdb_taxonomy ranked by checkm2_completeness - 5 x " \
           "checkm2_contamination desc (unparsed last), ties by GTDB accession; ranks 1..cap selected; no ANI step"
    if ani_out is not None:
        digests["ani_species.tsv"] = put(a.out / "ani_species.tsv", tsv_text(species_rows, SPECIES_COLS))
        rule = rule.replace("ranks 1..cap selected; no ANI step", "") + ani_out["rule"]
    cand_out = None
    if a.candidates:
        cand = sorted((r for rows in ranked.values() if len(rows) > cap for r in rows), key=lambda r: r["accession"])
        want = {r["accession"]: "ani_candidate" for r in cand if r["accession"] not in p1_set}
        want |= {acc: "phase1_not_served_retry" for acc in not_served if acc in by_acc}
        fetch = [dict(by_acc[acc], fetch_reason=want[acc]) for acc in sorted(want)]
        digests["ani_candidates.tsv"] = put(a.out / "ani_candidates.tsv", tsv_text(cand, COLS + EXTRA))
        digests["ani_fetch.tsv"] = put(a.out / "ani_fetch.tsv", tsv_text(fetch, COLS + EXTRA + ["fetch_reason"]))
        cand_out = {"species": sum(len(rows) > cap for rows in ranked.values()), "genomes": len(cand),
                    "bp": sum(r["genome_size"] for r in cand),
                    "held_by_phase1": sum(r["accession"] in p1_set for r in cand),
                    "fetch_genomes": len(fetch), "fetch_bp": sum(r["genome_size"] for r in fetch),
                    "fetch_ani_candidates": sum(v == "ani_candidate" for v in want.values()),
                    "fetch_phase1_retry": sorted(acc for acc, v in want.items() if v != "ani_candidate"),
                    "not_served_not_in_metadata": sorted(acc for acc in not_served if acc not in by_acc)}
    out = {
        "gtdb_release": "R226", "gtdb_files_sha256": sha, "cap": cap,
        "rule": rule,
        "metadata_columns_named_ani": header_ani,
        "complete_genomes": len(items), "species_with_complete_genome": len(ranked),
        "complete_genomes_checkm2_unparsed": unparsed, "complete_genomes_checkm2_nonfinite": nonfinite,
        "selected": summary(selected),
        "selected_gtdb_representative": sum(r["gtdb_representative"] for r in selected),
        "phase1": {"file": a.phase1.name, "sha256": hashlib.sha256(p1_bytes).hexdigest(), "genomes": len(p1),
                   "rank1_equals_file": tsv_text(rank1, COLS).encode() == p1_bytes,
                   "rank1_bp": p1_bp, "not_in_selection": missing, "not_served": ns},
        "remainder": summary(remainder),
        "ani_step": {
            "species_over_cap": len(over),
            "genomes_in_species_over_cap": sum(len(v) for v in over.values()),
            "selected_in_species_over_cap": cap * len(over),
            "max_selected_differ_keep_rank1": sum(min(cap - 1, len(v) - cap) for v in over.values()),
            "max_selected_differ_any": sum(min(cap, len(v) - cap) for v in over.values()),
            "remainder_in_species_over_cap": {"genomes": len(rem_over), "bp": sum(r["genome_size"] for r in rem_over)},
            "bp_bounds_keep_rank1": {"min": lo, "max": hi, "quality_order": sum(r["genome_size"] for r in q_selected),
                                     "remainder_min": lo - p1_bp, "remainder_max": hi - p1_bp},
            "species_score_tie_at_cut": tie_cut,
            "selected_decided_by_accession_at_cut": decided_by_acc,
            "species_rank1_score_equals_rank_cap_plus_1": cut_top},
        "by_species_size": list(bins.values()),
    }
    if ani_out is not None:
        out["ani"] = ani_out
    if cand_out is not None:
        out["candidates"] = cand_out
    out["outputs_sha256"] = digests
    text = json.dumps(out, indent=1) + "\n"
    put(a.out / "selection.json", text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())

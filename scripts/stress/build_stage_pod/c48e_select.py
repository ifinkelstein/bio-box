"""Card C-48e part 1: choose a tranche of SPIRE samples by habitat under two caps, as a pure function of the tables.

    python c48e_select.py run --tables-dir /tmp/if2453/s3/spire/tables \
        --pilot-stats /tmp/if2453/s3/spire/contigs/sample_contig_stats.tsv \
        --habitats aquatic,terrestrial --max-contigs 3000000 --max-bp 40000000000 \
        [--exclude TAKEN.tsv ...] [--seed 20261120] [--k 15] --out-dir DIR
    python c48e_select.py validate --tables-dir T --pilot-stats S [--ks 5 15 25 40] --out OUT.json

The order (card C-48e; the pilot's allocation, card C-42, made sequential so that any prefix is a tranche):
1. Universe: SPIRE v1 study samples (`c42_select.merge_sample_tables`) whose level-1 habitat is one of `habitats`.
   A sample without a published assembly size (`assembly_bp` absent or 0) is not ordered: it is listed apart
   (all seven such aquatic or terrestrial samples the pilot requested answered HTTP 404).
2. Within a level-2 stratum, studies are taken in turn: studies in the order of sha256("<seed>:study:<study>"),
   each study's samples in the order of sha256("<seed>:sample:<sample_id>"), then round by round one sample from
   each study that has one left. No order depends on the order rows are read in.
3. Strata are interleaved by Webster's divisor method on weights sqrt(N_h), N_h every sample of the stratum (as
   the pilot counted them): the next sample comes from the stratum with the largest sqrt(N_h) / (2 k_h + 1) among
   strata with samples left, k_h the samples it has given so far (ties: larger N_h, then the stratum's name). A
   divisor method never takes a sample back as the tranche grows, so every prefix is proportional to sqrt(N_h) in
   the divisor sense, and a second tranche continues the same sequence.
4. Each sample's contigs and bases at >= 5 kb: the pilot's measured values when the pilot fetched it, else the
   ratio of sums over the k pilot samples of the same habitats nearest in ln(N50) (ties by sample id), times the
   sample's published assembly size (*estimate*).
5. The tranche is the longest prefix of the order, skipping samples an earlier tranche took (`--exclude`), whose
   estimated contigs and bases at >= 5 kb both stay within the caps; the rest of the order is the remainder, in
   the order a next tranche would take it.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import sys
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

DEFAULT_SEED = 20261120
DEFAULT_K = 15


def _sibling(name: str):
    """Import a sibling script of this directory (the pilot's `c42_*`) without changing `sys.path`."""
    if name in sys.modules:
        return sys.modules[name]
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def rank_key(seed: int, kind: str, value: str) -> str:
    return hashlib.sha256(f"{seed}:{kind}:{value}".encode()).hexdigest()


def _sized(r: dict) -> bool:
    return r.get("assembly_bp") is not None and int(r["assembly_bp"]) > 0


def order_stratum(rows: Iterable[dict], seed: int) -> list[dict]:
    """One stratum's samples, studies taken in turn (module docstring, item 2)."""
    by_study: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_study[r["study"]].append(r)
    studies = sorted(by_study, key=lambda s: (rank_key(seed, "study", s), s))
    queues = [sorted(by_study[s], key=lambda r: (rank_key(seed, "sample", r["sample_id"]), r["sample_id"]))
              for s in studies]
    out: list[dict] = []
    depth = max((len(q) for q in queues), default=0)
    for i in range(depth):
        for q in queues:
            if i < len(q):
                out.append(q[i])
    return out


def interleave(orders: dict[str, list[dict]], counts: dict[str, int]) -> list[tuple[str, dict]]:
    """Webster's sequence over strata (module docstring, item 3): (stratum, row) in order."""
    weights = {h: math.sqrt(counts[h]) for h in orders}
    taken = {h: 0 for h in orders}
    out: list[tuple[str, dict]] = []
    total = sum(len(v) for v in orders.values())
    while len(out) < total:
        best = None
        for h in orders:
            if taken[h] >= len(orders[h]):
                continue
            key = (weights[h] / (2 * taken[h] + 1), counts[h], _neg(h))
            if best is None or key > best[0]:
                best = (key, h)
        h = best[1]
        out.append((h, orders[h][taken[h]]))
        taken[h] += 1
    return out


def _neg(s: str) -> tuple:
    """A key that sorts strings in reverse, so that `max` breaks a tie by the smaller name."""
    return tuple(-ord(c) for c in s) + (1,)


def tranche_order(samples: Iterable[dict], habitats: Iterable[str], seed: int = DEFAULT_SEED
                  ) -> tuple[list[dict], list[dict]]:
    """(ordered samples of `habitats` with a published size, samples of `habitats` without one), as a pure function
    of the rows and the seed. Ordered rows gain `order_rank`, `stratum`, `stratum_n` and `stratum_rank` (1-based)."""
    habs = set(habitats)
    rows = sorted((dict(r) for r in samples if r["habitat"] in habs), key=lambda r: r["sample_id"])
    ids = [r["sample_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate sample_id in the table")
    counts: dict[str, int] = defaultdict(int)
    for r in rows:
        counts[r["habitat2"]] += 1
    by_stratum: dict[str, list[dict]] = defaultdict(list)
    unsized = []
    for r in rows:
        (by_stratum[r["habitat2"]] if _sized(r) else unsized).append(r)
    orders = {h: order_stratum(by_stratum[h], seed) for h in sorted(by_stratum)}
    out = []
    rank_in: dict[str, int] = defaultdict(int)
    for i, (h, r) in enumerate(interleave(orders, counts), 1):
        rank_in[h] += 1
        out.append({**r, "order_rank": i, "stratum": h, "stratum_n": counts[h], "stratum_rank": rank_in[h]})
    unsized_out = [{**r, "stratum": r["habitat2"], "stratum_n": counts[r["habitat2"]]}
                   for r in sorted(unsized, key=lambda r: (r["habitat2"], r["sample_id"]))]
    return out, unsized_out


def estimate_5k(rows: Iterable[dict], train: list[dict], k: int = DEFAULT_K) -> list[dict]:
    """Rows with `est_n_ge_5000`, `est_bp_ge_5000` and `est_source` (module docstring, item 4)."""
    by_id = {t["sample_id"]: t for t in train}
    pts = [(math.log(t["assembly_n50"]), t["sample_id"], t) for t in train
           if t.get("assembly_n50") and int(t["assembly_n50"]) > 0 and int(t["total_bp"]) > 0]
    pooled_bp = sum(int(t["total_bp"]) for t in train)
    pooled = ((sum(int(t["bp_ge_5000"]) for t in train) / pooled_bp, sum(int(t["n_ge_5000"]) for t in train) / pooled_bp)
              if pooled_bp else (0.0, 0.0))
    out = []
    for r in rows:
        r = dict(r)
        sid = r["sample_id"]
        if sid in by_id:
            t = by_id[sid]
            r.update(est_bp_ge_5000=int(t["bp_ge_5000"]), est_n_ge_5000=int(t["n_ge_5000"]), est_source="measured_pilot")
        elif not _sized(r):
            r.update(est_bp_ge_5000=None, est_n_ge_5000=None, est_source="no_assembly_size")
        else:
            bp = int(r["assembly_bp"])
            n50 = r.get("assembly_n50")
            if n50 and int(n50) > 0 and pts:
                x = math.log(int(n50))
                near = sorted(pts, key=lambda p: (abs(p[0] - x), p[1]))[:k]
                tb = sum(int(p[2]["total_bp"]) for p in near)
                rb = sum(int(p[2]["bp_ge_5000"]) for p in near) / tb
                rn = sum(int(p[2]["n_ge_5000"]) for p in near) / tb
                src = "nearest_n50"
            else:
                rb, rn = pooled
                src = "pooled_ratio"
            r.update(est_bp_ge_5000=round(bp * rb), est_n_ge_5000=round(bp * rn), est_source=src)
        out.append(r)
    return out


def take(ordered: list[dict], max_contigs: int, max_bp: int, exclude: Iterable[str] = ()
         ) -> tuple[list[dict], list[dict]]:
    """(the tranche, the remainder) (module docstring, item 5). Excluded samples are in neither."""
    skip = set(exclude)
    taken, rest = [], []
    n = bp = 0
    full = False
    for r in ordered:
        if r["sample_id"] in skip:
            continue
        if not full:
            rn, rb = int(r["est_n_ge_5000"]), int(r["est_bp_ge_5000"])
            if n + rn <= max_contigs and bp + rb <= max_bp:
                n, bp = n + rn, bp + rb
                taken.append({**r, "cum_est_n_ge_5000": n, "cum_est_bp_ge_5000": bp})
                continue
            full = True
        rest.append(r)
    return taken, rest


# --- the command ------------------------------------------------------------------------------------------------


def _write_tsv(path: Path, rows: list[dict], cols: list[str], gz: bool = False) -> None:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})
    data = buf.getvalue().encode()
    tmp = path.with_name(f".{path.name}.tmp")
    if gz:
        with open(tmp, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as fh:
            fh.write(data)
    else:
        tmp.write_bytes(data)
    tmp.replace(path)


def _read_ids(path: Path) -> set[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        return {r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")}


SEL_COLS = ["order_rank", "sample_id", "study", "habitat", "habitat2", "stratum_n", "stratum_rank", "assembly_bp",
            "assembly_n50", "est_n_ge_5000", "est_bp_ge_5000", "est_source", "cum_est_n_ge_5000", "cum_est_bp_ge_5000"]
REM_COLS = ["order_rank", "sample_id", "study", "habitat", "habitat2", "stratum_n", "stratum_rank", "assembly_bp",
            "assembly_n50", "est_n_ge_5000", "est_bp_ge_5000", "est_source", "status"]


def cmd_run(a: argparse.Namespace) -> int:
    c42 = _sibling("c42_select")
    tdir = Path(a.tables_dir)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    habitats = tuple(sorted(h.strip() for h in a.habitats.split(",") if h.strip()))

    merged = c42.merge_sample_tables(c42._load_study_tsvs(tdir / "study_tsv", tdir / "studies.html"),
                                     c42._load_api_samples(tdir / "api_samples"),
                                     c42._load_mag_bp(tdir / "spire_v1_genome_metadata.tsv.gz"))
    spire = [r for r in merged if r["in_study_tsv"]]
    by_id = {r["sample_id"]: r for r in spire}

    stats: dict[str, dict] = {}
    with open(a.pilot_stats, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stats[r["sample_id"]] = r
    train = []
    for sid in sorted(stats):
        r = by_id.get(sid)
        if r is None or r["habitat"] not in habitats or int(stats[sid]["total_bp"]) <= 0:
            continue
        train.append({"sample_id": sid, "assembly_n50": r["assembly_n50"], "total_bp": int(stats[sid]["total_bp"]),
                      "bp_ge_5000": int(stats[sid]["bp_ge_5000"]), "n_ge_5000": int(stats[sid]["n_ge_5000"])})

    exclude: set[str] = set()
    for p in a.exclude or []:
        exclude |= _read_ids(Path(p))

    order, unsized = tranche_order(spire, habitats, a.seed)
    est = estimate_5k(order, train, a.k)
    taken, rest = take(est, a.max_contigs, a.max_bp, exclude)
    rem = [{**r, "status": "not_taken"} for r in rest] + [
        {**r, "est_source": "no_assembly_size", "status": "no_assembly_size"} for r in unsized if r["sample_id"] not in exclude]

    _write_tsv(out / "selection.tsv", taken, SEL_COLS)
    _write_tsv(out / "remainder.tsv.gz", rem, REM_COLS, gz=True)

    # by stratum: the order's sizes, what the tranche takes, what is left
    strata = sorted({r["habitat2"] for r in spire if r["habitat"] in habitats})
    tk = defaultdict(lambda: [0, 0, 0, 0])
    for r in taken:
        v = tk[r["habitat2"]]
        v[0] += 1; v[1] += int(r["assembly_bp"]); v[2] += int(r["est_n_ge_5000"]); v[3] += int(r["est_bp_ge_5000"])
    lf = defaultdict(lambda: [0, 0, 0, 0, 0])
    for r in rem:
        v = lf[r["habitat2"]]
        v[0] += 1
        if r["status"] == "not_taken":
            v[1] += int(r["assembly_bp"]); v[2] += int(r["est_n_ge_5000"]); v[3] += int(r["est_bp_ge_5000"])
        else:
            v[4] += 1
    strat_rows = []
    for h in strata:
        n_h = sum(1 for r in spire if r["habitat2"] == h)
        strat_rows.append({"habitat2": h, "stratum_n": n_h, "sqrt_weight": round(math.sqrt(n_h), 4),
                           "taken": tk[h][0], "taken_assembly_bp": tk[h][1], "taken_est_n_ge_5000": tk[h][2],
                           "taken_est_bp_ge_5000": tk[h][3], "left": lf[h][0], "left_assembly_bp": lf[h][1],
                           "left_est_n_ge_5000": lf[h][2], "left_est_bp_ge_5000": lf[h][3],
                           "left_no_assembly_size": lf[h][4]})
    _write_tsv(out / "by_stratum.tsv", strat_rows, list(strat_rows[0].keys()))

    # other habitats: counts and bases by level-1 habitat, 5 kb by the pilot's measured share in that habitat
    pil = defaultdict(lambda: [0, 0, 0])
    for sid, s in stats.items():
        r = by_id.get(sid)
        if r is not None:
            v = pil[r["habitat"]]
            v[0] += int(s["total_bp"]); v[1] += int(s["bp_ge_5000"]); v[2] += int(s["n_ge_5000"])
    oth = []
    for h1 in sorted({r["habitat"] for r in spire} - set(habitats)):
        rows = [r for r in spire if r["habitat"] == h1]
        known = sum(r["assembly_bp"] or 0 for r in rows)
        tb, b5, n5 = pil[h1]
        oth.append({"habitat": h1, "samples": len(rows), "samples_with_assembly_bp": sum(1 for r in rows if _sized(r)),
                    "assembly_bp": known, "pilot_share_bp_ge_5000": round(b5 / tb, 6) if tb else None,
                    "pilot_n_ge_5000_per_gbp": round(n5 / (tb / 1e9), 1) if tb else None,
                    "est_bp_ge_5000": round(known * b5 / tb) if tb else None,
                    "est_n_ge_5000": round(known * n5 / tb) if tb else None})
    _write_tsv(out / "other_habitats.tsv", oth, list(oth[0].keys()))

    summary = {
        "card": "C-48e part 1",
        "habitats": list(habitats),
        "seed": a.seed,
        "k_nearest_n50": a.k,
        "caps": {"max_contigs_ge_5000": a.max_contigs, "max_bp_ge_5000": a.max_bp},
        "rule": __doc__.split("The order", 1)[1].split("\n", 1)[1].strip(),
        "spire_study_samples": len(spire),
        "habitat_samples": sum(1 for r in spire if r["habitat"] in habitats),
        "ordered_samples": len(order),
        "no_assembly_size_samples": len(unsized),
        "excluded_by_earlier_tranches": len(exclude),
        "train_pilot_samples": len(train),
        "taken": len(taken),
        "taken_strata": len({r["habitat2"] for r in taken}),
        "taken_studies": len({r["study"] for r in taken}),
        "taken_measured_in_pilot": sum(1 for r in taken if r["est_source"] == "measured_pilot"),
        "taken_assembly_bp": sum(int(r["assembly_bp"]) for r in taken),
        "taken_est_n_ge_5000": sum(int(r["est_n_ge_5000"]) for r in taken),
        "taken_est_bp_ge_5000": sum(int(r["est_bp_ge_5000"]) for r in taken),
        "next_sample": ({k: rest[0][k] for k in ("order_rank", "sample_id", "habitat2", "est_n_ge_5000", "est_bp_ge_5000")}
                        if rest else None),
        "left_ordered": len(rest),
        "left_assembly_bp": sum(int(r["assembly_bp"]) for r in rest),
        "left_est_n_ge_5000": sum(int(r["est_n_ge_5000"]) for r in rest),
        "left_est_bp_ge_5000": sum(int(r["est_bp_ge_5000"]) for r in rest),
        "habitats_all_est_n_ge_5000": sum(int(r["est_n_ge_5000"]) for r in est),
        "habitats_all_est_bp_ge_5000": sum(int(r["est_bp_ge_5000"]) for r in est),
        "habitats_all_assembly_bp": sum(int(r["assembly_bp"]) for r in est),
    }
    (out / "select_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "rule"}))
    return 0


def validate(train: list[dict], studies: dict[str, str], k: int, draws: int = 2000, size: int = 300,
             seed: int = DEFAULT_SEED) -> dict:
    """Leave-one-study-out check of `estimate_5k` on the pilot's samples: each sample estimated from the pilot samples
    of other studies; the ratio of estimated to measured totals over all, and its 5th and 95th percentiles over
    `draws` resamples (with replacement) of `size` samples."""
    est = []
    for t in train:
        others = [u for u in train if studies[u["sample_id"]] != studies[t["sample_id"]]]
        row = {"sample_id": t["sample_id"], "assembly_bp": t["total_bp"], "assembly_n50": t["assembly_n50"]}
        e = estimate_5k([row], others, k)[0]
        est.append((e["est_n_ge_5000"], int(t["n_ge_5000"]), e["est_bp_ge_5000"], int(t["bp_ge_5000"])))
    import random

    rng = random.Random(seed)
    rn, rb = [], []
    for _ in range(draws):
        idx = [rng.randrange(len(est)) for _ in range(size)]
        rn.append(sum(est[i][0] for i in idx) / sum(est[i][1] for i in idx))
        rb.append(sum(est[i][2] for i in idx) / sum(est[i][3] for i in idx))
    rn.sort()
    rb.sort()
    q = lambda v, p: round(v[int(p * (len(v) - 1))], 4)  # noqa: E731
    return {"k": k, "samples": len(est), "studies": len({studies[t["sample_id"]] for t in train}),
            "ratio_n_all": round(sum(e[0] for e in est) / sum(e[1] for e in est), 4),
            "ratio_bp_all": round(sum(e[2] for e in est) / sum(e[3] for e in est), 4),
            f"ratio_n_p05_p95_of_{size}": [q(rn, 0.05), q(rn, 0.95)],
            f"ratio_bp_p05_p95_of_{size}": [q(rb, 0.05), q(rb, 0.95)], "draws": draws, "seed": seed}


def cmd_validate(a: argparse.Namespace) -> int:
    c42 = _sibling("c42_select")
    tdir = Path(a.tables_dir)
    habitats = tuple(sorted(h.strip() for h in a.habitats.split(",") if h.strip()))
    merged = c42.merge_sample_tables(c42._load_study_tsvs(tdir / "study_tsv", tdir / "studies.html"),
                                     c42._load_api_samples(tdir / "api_samples"),
                                     c42._load_mag_bp(tdir / "spire_v1_genome_metadata.tsv.gz"))
    by_id = {r["sample_id"]: r for r in merged if r["in_study_tsv"]}
    train = []
    with open(a.pilot_stats, newline="") as fh:
        for s in sorted(csv.DictReader(fh, delimiter="\t"), key=lambda r: r["sample_id"]):
            r = by_id.get(s["sample_id"])
            if r is None or r["habitat"] not in habitats or int(s["total_bp"]) <= 0:
                continue
            train.append({"sample_id": s["sample_id"], "assembly_n50": r["assembly_n50"], "total_bp": int(s["total_bp"]),
                          "bp_ge_5000": int(s["bp_ge_5000"]), "n_ge_5000": int(s["n_ge_5000"])})
    studies = {t["sample_id"]: by_id[t["sample_id"]]["study"] for t in train}
    res = {"habitats": list(habitats), "by_k": [validate(train, studies, k) for k in a.ks]}
    Path(a.out).write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--tables-dir", required=True)
    v.add_argument("--pilot-stats", required=True)
    v.add_argument("--habitats", default="aquatic,terrestrial")
    v.add_argument("--ks", type=int, nargs="+", default=[5, 15, 25, 40])
    v.add_argument("--out", required=True)
    p = sub.add_parser("run")
    p.add_argument("--tables-dir", required=True)
    p.add_argument("--pilot-stats", required=True)
    p.add_argument("--habitats", default="aquatic,terrestrial")
    p.add_argument("--max-contigs", type=int, required=True)
    p.add_argument("--max-bp", type=int, required=True)
    p.add_argument("--exclude", action="append")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--k", type=int, default=DEFAULT_K)
    p.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    return cmd_validate(a) if a.cmd == "validate" else cmd_run(a)


if __name__ == "__main__":
    sys.exit(main())

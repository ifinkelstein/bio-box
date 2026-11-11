"""Card 0-09b: what scanning cluster representatives only would save, and how often an inferred Pfam domain
call would be wrong. Reads only; writes one JSON file.

    python representatives.py <root> <release> <out.json> [--domains] [--library pfam-38.2] [--threads N]

Files read are the producer's own lists, never globs: the release manifest's `protein_seq`, `cluster_member`
and `cluster_member_aln` files (each checked against the manifest's sha256), and with `--domains` each
stage's `l1_domains` files named in its published run record `runs/stage=<S>/domains-<L>.json` (checked
against the record's sha256).

Saving (always): for each level (nr90, nr50) and each scope (each stage of the release, and all proteins),
the proteins and residues that are not their cluster's representative. A protein's length is its sequence
without '*', the sequence `domains.scan` scans. A protein in two stages counts in both stages' scopes and
once in "all".

Inference (`--domains`): for every member that is not a representative, its own domain set D(m) (distinct
model accessions with an `l1_domains` row) against its representative's D(r), and its architecture (the
accessions ordered by env_start, env_end, accession, repeats kept) against the representative's. Inference
reports D(r) for m: pairs right |D(m) & D(r)|, false |D(r) - D(m)|, missed |D(m) - D(r)|. Classes, mutually
exclusive: equal (D(m) == D(r), both empty included), subset (D(m) < D(r), D(m) empty included), superset
(D(m) > D(r), D(r) empty included), disjoint (both non-empty, no common model), overlap (the rest).

Coordinates: for members whose set equals the representative's and is not empty, each member domain is paired
with the representative's domain of the same model whose envelope, mapped through the stored member
alignment (`cluster_member_aln`: q* are the member's coordinates, t* the representative's, 0-based half-open;
CIGAR from MMseqs2 with query = representative: M both, I the representative only, D the member only; the
convention is checked on every row), lies closest; reported as the absolute start and end difference in
residues, the same without mapping (coordinates copied), and as fractions of each protein's length.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

LEVELS = ("nr90", "nr50")
PIDENT_BINS = [0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0]  # last bin: pident == 1.0
SIZE_BINS = [(2, 2), (3, 4), (5, 9), (10, 99), (100, 10**12)]
CIGAR = re.compile(rb"(\d+)([MID])")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def checked(root: Path, f: dict, read: list) -> Path:
    p = root / f["path"]
    got = sha256(p)
    if got != f["sha256"]:
        raise SystemExit(f"{p}: sha256 {got} is not the record's {f['sha256']}")
    read.append({"path": f["path"], "bytes": p.stat().st_size, "sha256": got})
    return p


def stage_of(path: str) -> str:
    return re.search(r"/stage=([^/]+)/", path).group(1)


def sql_list(paths: list[Path]) -> str:
    return "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for p in paths) + "]"


def load_proteins(con, root: Path, manifest: dict, stages: list[str], read: list) -> dict:
    files = [f for f in manifest["files"] if f["table"] == "protein_seq"]
    parts = []
    differ = internal = 0
    for bit, stage in enumerate(stages):
        for f in sorted((f for f in files if stage_of(f["path"]) == stage), key=lambda f: f["path"]):
            p = checked(root, f, read)
            for b in pq.ParquetFile(p).iter_batches(batch_size=65536, columns=["protein_uid", "length", "seq"]):
                n = pc.utf8_length(pc.replace_substring(b.column("seq"), "*", "")).cast(pa.int64())
                differ += int(pc.sum(pc.not_equal(n, b.column("length").cast(pa.int64()))).as_py() or 0)
                internal += int(pc.sum(pc.match_substring_regex(b.column("seq"), r"\*.")).as_py() or 0)
                parts.append(pa.table({"uid": b.column("protein_uid"), "len": n,
                                       "mask": pa.array(np.full(b.num_rows, 1 << bit, dtype=np.int64))}))
    raw = pa.concat_tables(parts)
    con.register("raw_prot", raw)
    con.execute("CREATE TABLE prot AS SELECT uid, max(len) AS len, bit_or(mask) AS mask, count(DISTINCT len) AS nlen, "
                "count(*) AS nrows FROM raw_prot GROUP BY uid")
    con.unregister("raw_prot")
    bad = con.execute("SELECT count(*) FROM prot WHERE nlen <> 1").fetchone()[0]
    if bad:
        raise SystemExit(f"{bad} proteins have two lengths")
    multi = con.execute("SELECT count(*) FROM prot WHERE nrows > 1").fetchone()[0]
    return {"rows": raw.num_rows, "distinct": con.execute("SELECT count(*) FROM prot").fetchone()[0],
            "in_more_than_one_stage": multi, "stripped_length_differs_from_length_column": differ,
            "rows_with_a_stop_before_the_last_residue": internal}


def scopes(stages: list[str]) -> list[tuple[str, int]]:
    return [(s, 1 << i) for i, s in enumerate(stages)] + [("all", (1 << len(stages)) - 1)]


def saving(con, stages: list[str]) -> dict:
    out = {}
    for level in LEVELS:
        for name, bits in scopes(stages):
            r = con.execute(
                "SELECT count(*), sum(p.len), count(*) FILTER (WHERE cm.m = cm.r), "
                "coalesce(sum(p.len) FILTER (WHERE cm.m = cm.r), 0), count(DISTINCT cm.c), "
                "count(*) FILTER (WHERE cm.m <> cm.r AND (rp.mask & $bits) = 0), "
                "count(*) FILTER (WHERE cs.n = 1), "
                "count(*) FILTER (WHERE cm.m <> cm.r AND abs(p.len - rp.len) > 0.2 * rp.len), "
                "coalesce(sum(p.len) FILTER (WHERE cm.m <> cm.r AND abs(p.len - rp.len) > 0.2 * rp.len), 0) "
                "FROM cm JOIN prot p ON p.uid = cm.m JOIN prot rp ON rp.uid = cm.r "
                "JOIN cs ON cs.level = cm.level AND cs.c = cm.c "
                "WHERE cm.level = $level AND (p.mask & $bits) <> 0", {"bits": bits, "level": level}).fetchone()
            n, res, reps, rep_res, clusters, rep_elsewhere, singletons, hyb_n, hyb_res = r
            out[f"{level}/{name}"] = {
                "proteins": n, "residues": int(res), "reps": reps, "rep_residues": int(rep_res),
                "clusters_touched": clusters, "members_whose_rep_is_outside_scope": rep_elsewhere,
                "proteins_in_singleton_clusters": singletons,
                "members_length_differs_gt20": hyb_n, "members_length_differs_gt20_residues": int(hyb_res),
                "hybrid_residue_fraction_scanned": (int(rep_res) + int(hyb_res)) / int(res),
                "nonrep_protein_fraction": 1 - reps / n, "nonrep_residue_fraction": 1 - rep_res / int(res)}
    return out


def load_clusters(con, root: Path, manifest: dict, read: list, n_prot: int) -> dict:
    mem = [checked(root, f, read) for f in manifest["files"] if f["table"] == "cluster_member"]
    aln = [checked(root, f, read) for f in manifest["files"] if f["table"] == "cluster_member_aln"]
    con.execute(f"CREATE TABLE cm AS SELECT member_uid AS m, level, cluster_uid AS c, rep_uid AS r "
                f"FROM read_parquet({sql_list(mem)})")
    con.execute("CREATE TABLE cs AS SELECT level, c, count(*) AS n FROM cm GROUP BY level, c")
    checks = {}
    for level in LEVELS:
        n, d = con.execute("SELECT count(*), count(DISTINCT m) FROM cm WHERE level = ?", [level]).fetchone()
        missing = con.execute("SELECT count(*) FROM prot WHERE uid NOT IN (SELECT m FROM cm WHERE level = ?)",
                              [level]).fetchone()[0]
        extra = con.execute("SELECT count(*) FROM cm WHERE level = ? AND m NOT IN (SELECT uid FROM prot)",
                            [level]).fetchone()[0]
        selfless = con.execute("SELECT count(DISTINCT r) FROM cm WHERE level = ? AND r NOT IN "
                               "(SELECT m FROM cm WHERE level = ? AND m = r)", [level, level]).fetchone()[0]
        if not (n == d == n_prot and missing == extra == selfless == 0):
            raise SystemExit(f"{level}: cluster_member does not cover the proteins once "
                             f"({n}, {d}, {n_prot}, {missing}, {extra}, {selfless})")
        checks[level] = {"rows": n, "clusters": con.execute("SELECT count(*) FROM cs WHERE level = ?",
                                                             [level]).fetchone()[0]}
    return {"member_files": [str(p) for p in mem], "aln_files": [str(p) for p in aln], "levels": checks}


def load_domains(con, root: Path, stages: list[str], library: str, read: list) -> dict:
    parts = []
    for stage in stages:
        rec = json.loads((root / f"runs/stage={stage}/domains-{library}.json").read_text())
        if rec.get("exit_code") != 0:
            raise SystemExit(f"domains-{library} of stage {stage}: exit {rec.get('exit_code')}")
        for o in rec["outputs"]:
            if o["path"].startswith("l1/domains/"):
                parts.append((stage, checked(root, o, read)))
    union = " UNION ALL ".join(f"SELECT '{s}' AS stage, protein_uid AS uid, model_acc AS acc, env_start AS s, "
                               f"env_end AS e FROM read_parquet('{p}')" for s, p in parts)
    con.execute(f"CREATE TABLE doms AS {union}")
    # a protein in two stages was scanned twice: the rows must agree
    disagree = con.execute(
        "WITH multi AS (SELECT uid FROM prot WHERE nrows > 1), "
        "a AS (SELECT stage, uid, acc, s, e FROM doms WHERE uid IN (SELECT uid FROM multi)), "
        "cnt AS (SELECT uid, acc, s, e, count(DISTINCT stage) k FROM a GROUP BY ALL) "
        "SELECT count(*) FROM cnt WHERE k <> (SELECT bit_count(mask) FROM prot WHERE prot.uid = cnt.uid)").fetchone()[0]
    con.execute("CREATE TABLE dom AS SELECT DISTINCT uid, acc, s, e FROM doms")
    outside = con.execute("SELECT count(*) FROM dom LEFT JOIN prot USING (uid) WHERE prot.uid IS NULL OR dom.e > prot.len "
                          "OR dom.s >= dom.e").fetchone()[0]
    if outside:
        raise SystemExit(f"{outside} domain rows are on no protein or past its end")
    rows, prots = con.execute("SELECT count(*), count(DISTINCT uid) FROM dom").fetchone()
    pairs = con.execute("SELECT count(*) FROM (SELECT DISTINCT uid, acc FROM dom)").fetchone()[0]
    return {"rows": rows, "proteins_with_domain": prots, "pairs": pairs,
            "rows_disagreeing_between_stages_for_shared_proteins": disagree}


def domain_index(con) -> dict:
    """uid -> (frozenset of accessions, architecture tuple, list of (acc, s, e) sorted by s, e)."""
    idx: dict[bytes, tuple] = {}
    cur = None
    rows: list = []
    reader = con.execute("SELECT uid, acc, s, e FROM dom ORDER BY uid, s, e, acc").to_arrow_reader(262144)
    for batch in reader:
        for u, a, s, e in zip(*(batch.column(i).to_pylist() for i in range(4))):
            if u != cur:
                if cur is not None:
                    idx[cur] = (frozenset(x[0] for x in rows), tuple(x[0] for x in rows), rows)
                cur, rows = u, []
            rows.append((a, s, e))
    if cur is not None:
        idx[cur] = (frozenset(x[0] for x in rows), tuple(x[0] for x in rows), rows)
    return idx


class Mapper:
    """Rep coordinate -> member coordinate through one alignment."""

    def __init__(self, cigar: bytes, qs: int, qe: int, ts: int, te: int, lm: int):
        self.qs, self.qe, self.ts, self.te, self.lm = qs, qe, ts, te, lm
        self.segs = []
        rp, mp = ts, qs
        nm = ni = nd = 0
        for n, op in CIGAR.findall(cigar):
            n = int(n)
            if op == b"M":
                self.segs.append((rp, rp + n, mp, True))
                rp += n
                mp += n
                nm += n
            elif op == b"I":
                self.segs.append((rp, rp + n, mp, False))
                rp += n
                ni += n
            else:
                mp += n
                nd += n
        self.ok = (rp == te and mp == qe)

    def _pos(self, p: int) -> tuple[int, bool]:
        if p < self.ts:
            return self.qs - (self.ts - p), True
        if p >= self.te:
            return self.qe + (p - self.te), True
        for r0, r1, m0, match in self.segs:
            if r0 <= p < r1:
                return (m0 + (p - r0), True) if match else (m0, False)
        raise AssertionError("position not in alignment")

    def start(self, s: int) -> int:
        m, _ = self._pos(s)
        return min(max(m, 0), self.lm)

    def end(self, e: int) -> int:
        m, exact = self._pos(e - 1)
        m = m + 1 if exact else m
        return min(max(m, 0), self.lm)


def infer(con, idx: dict, stages: list[str], level: str) -> dict:
    """Per non-representative member arrays, and the coordinate differences."""
    res = con.execute(
        "SELECT cm.m, cm.r, cm.c, pm.mask, pr.mask, pm.len, pr.len, cs.n, a.pident, a.qstart, a.qend, a.tstart, "
        "a.tend, a.cigar FROM cm JOIN prot pm ON pm.uid = cm.m JOIN prot pr ON pr.uid = cm.r "
        "JOIN cs ON cs.level = cm.level AND cs.c = cm.c "
        "LEFT JOIN aln a ON a.member_uid = cm.m AND a.level = cm.level AND a.rep_uid = cm.r "
        "WHERE cm.level = ? AND cm.m <> cm.r ORDER BY cm.m", [level])
    cols = {k: [] for k in ("n_m", "n_r", "n_i", "arch_eq", "mask", "rmask", "lm", "lr", "csize", "pident", "cidx")}
    coord = {k: [] for k in ("ds", "de", "cs_", "ce", "fs", "fe", "mask", "short", "arch_eq")}
    cl_ids: dict[bytes, int] = {}
    empty = (frozenset(), (), [])
    no_aln = bad_cigar = 0
    reader = res.to_arrow_reader(131072)
    for batch in reader:
        v = [batch.column(i).to_pylist() for i in range(batch.num_columns)]
        for m, r, c, mk, rmk, lm, lr, n, pid, qs, qe, ts, te, cig in zip(*v):
            dm = idx.get(m, empty)
            dr = idx.get(r, empty)
            ni = len(dm[0] & dr[0])
            ae = dm[1] == dr[1]
            cols["n_m"].append(len(dm[0]))
            cols["n_r"].append(len(dr[0]))
            cols["n_i"].append(ni)
            cols["arch_eq"].append(ae)
            cols["mask"].append(mk)
            cols["rmask"].append(rmk)
            cols["lm"].append(lm)
            cols["lr"].append(lr)
            cols["csize"].append(n)
            cols["pident"].append(np.nan if pid is None else pid)
            cols["cidx"].append(cl_ids.setdefault(c, len(cl_ids)))
            if cig is None:
                no_aln += 1
                continue
            mp = Mapper(cig, qs, qe, ts, te, lm)
            if not mp.ok:
                bad_cigar += 1
                continue
            if not dm[0] or dm[0] != dr[0]:
                continue
            short = lm < 0.8 * lr
            for acc, s, e in dm[2]:
                best = None
                for acc_r, sr, er in dr[2]:
                    if acc_r != acc:
                        continue
                    ms, me = mp.start(sr), mp.end(er)
                    d = abs(ms - s) + abs(me - e)
                    if best is None or d < best[0]:
                        best = (d, ms, me, sr, er)
                _, ms, me, sr, er = best
                coord["ds"].append(abs(ms - s))
                coord["de"].append(abs(me - e))
                coord["cs_"].append(abs(sr - s))
                coord["ce"].append(abs(er - e))
                coord["fs"].append(abs(s / lm - sr / lr))
                coord["fe"].append(abs(e / lm - er / lr))
                coord["mask"].append(mk)
                coord["short"].append(short)
                coord["arch_eq"].append(ae)
    a = {k: np.asarray(v) for k, v in cols.items()}
    cc = {k: np.asarray(v) for k, v in coord.items()}
    return {"a": a, "coord": cc, "no_alignment": no_aln, "cigar_span_mismatch": bad_cigar}


def check_alignments(con, seq_files: list[Path], level: str, n: int) -> dict:
    """Independent check of the CIGAR convention and coordinates: on the first `n` non-representative
    members by uid, walk the CIGAR over the two sequences and count identical residues in M columns; that
    count over `alnlen` must be the stored `pident` (MMseqs2's fident, 3 decimals), and the walk must end at
    the stored ends. Done both ways (I on the representative, as assumed; and I on the member)."""
    rows = con.execute(
        f"WITH s AS (SELECT a.member_uid, a.rep_uid, a.qstart, a.qend, a.tstart, a.tend, a.pident, a.alnlen, a.cigar "
        f"FROM aln a WHERE a.level = ? AND a.member_uid <> a.rep_uid ORDER BY a.member_uid LIMIT {int(n)}), "
        f"q AS (SELECT protein_uid, seq FROM read_parquet({sql_list(seq_files)}) "
        f"WHERE protein_uid IN (SELECT member_uid FROM s UNION SELECT rep_uid FROM s)), "
        f"qd AS (SELECT DISTINCT protein_uid, seq FROM q) "
        f"SELECT s.*, mq.seq, rq.seq FROM s JOIN qd mq ON mq.protein_uid = s.member_uid "
        f"JOIN qd rq ON rq.protein_uid = s.rep_uid", [level]).fetchall()
    out = {}
    for name, i_on_rep in (("I_on_representative", True), ("I_on_member", False)):
        agree = ends = 0
        for _m, _r, qs, qe, ts, te, pid, alnlen, cig, ms, rs in rows:
            mp, rp, same = qs, ts, 0
            for k, op in CIGAR.findall(cig):
                k = int(k)
                if op == b"M":
                    same += sum(1 for j in range(k) if mp + j < len(ms) and rp + j < len(rs) and ms[mp + j] == rs[rp + j])
                    mp += k
                    rp += k
                elif (op == b"I") == i_on_rep:
                    rp += k
                else:
                    mp += k
            ends += (mp == qe and rp == te)
            agree += abs(same / alnlen - pid) <= 0.0015
        out[name] = {"members": len(rows), "ends_match": ends, "identity_matches_pident": agree}
    return out


def frac(x: float, y: float) -> float | None:
    return None if not y else x / y


def summarize(a: dict, sel: np.ndarray) -> dict:
    k = int(sel.sum())
    if not k:
        return {"members": 0}
    nm, nr, ni = a["n_m"][sel], a["n_r"][sel], a["n_i"][sel]
    eq = (ni == nm) & (ni == nr)
    sub = (ni == nm) & (nm < nr)
    sup = (ni == nr) & (nr < nm)
    dis = (nm > 0) & (nr > 0) & (ni == 0)
    ovl = ~(eq | sub | sup | dis)
    anyd = (nm > 0) | (nr > 0)
    inferred, true, right = int(nr.sum()), int(nm.sum()), int(ni.sum())
    _, inv = np.unique(a["cidx"][sel], return_inverse=True)
    cnt = np.bincount(inv)
    inf_c, true_c, right_c = (np.bincount(inv, weights=x.astype(float)) for x in (nr, nm, ni))
    eq_c = np.bincount(inv, weights=eq.astype(float)) / cnt

    def per_cluster(num: np.ndarray, den: np.ndarray) -> float | None:
        ok = den > 0
        return float((num[ok] / den[ok]).mean()) if ok.any() else None

    ae = a["arch_eq"][sel]
    return {
        "members": k, "clusters": int(cnt.size),
        "equal": eq.mean(), "subset": sub.mean(), "superset": sup.mean(), "disjoint": dis.mean(),
        "overlap": ovl.mean(),
        "member_empty_rep_not": ((nm == 0) & (nr > 0)).mean(), "rep_empty_member_not": ((nr == 0) & (nm > 0)).mean(),
        "both_empty": (~anyd).mean(),
        "equal_among_members_with_any_domain": frac(float((eq & anyd).sum()), float(anyd.sum())),
        "architecture_equal": ae.mean(), "architecture_equal_given_set_equal": frac(float((ae & eq).sum()), float(eq.sum())),
        "pairs_inferred": inferred, "pairs_true": true, "pairs_right": right,
        "pairs_false": inferred - right, "pairs_missed": true - right,
        "false_fraction_of_inferred": frac(inferred - right, inferred),
        "missed_fraction_of_true": frac(true - right, true),
        "cluster_weighted": {
            "equal": float(eq_c.mean()),
            "false_fraction_of_inferred": per_cluster(inf_c - right_c, inf_c),
            "missed_fraction_of_true": per_cluster(true_c - right_c, true_c),
        },
    }


def quantiles(x: np.ndarray) -> dict:
    if not x.size:
        return {"n": 0}
    q = np.quantile(x, [0.5, 0.9, 0.95, 0.99])
    return {"n": int(x.size), "median": float(q[0]), "p90": float(q[1]), "p95": float(q[2]), "p99": float(q[3]),
            "mean": float(x.mean())}


def coord_summary(cc: dict, sel: np.ndarray) -> dict:
    out = {k: quantiles(cc[k][sel].astype(float)) for k in ("ds", "de", "cs_", "ce", "fs", "fe")}
    both = np.maximum(cc["ds"][sel], cc["de"][sel])
    out["both_ends_within_5_mapped"] = float((both <= 5).mean()) if both.size else None
    out["both_ends_within_10_mapped"] = float((both <= 10).mean()) if both.size else None
    return out


def length_class(a: dict) -> dict:
    lm, lr = a["lm"].astype(float), a["lr"].astype(float)
    return {"shorter_gt20": lm < 0.8 * lr, "longer_gt20": lm > 1.2 * lr,
            "within20": (lm >= 0.8 * lr) & (lm <= 1.2 * lr)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("release")
    ap.add_argument("out", type=Path)
    ap.add_argument("--domains", action="store_true")
    ap.add_argument("--library", default="pfam-38.2")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    t0 = time.monotonic()
    root = a.root
    manifest = json.loads((root / f"manifests/release={a.release}.json").read_text())
    stages = list(manifest["stages"])
    tmp = a.out.parent / f"duck-{a.release}"
    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(config={"threads": a.threads, "memory_limit": "20GB", "temp_directory": str(tmp)})
    read: list = []
    doc: dict = {"root": str(root), "release": a.release, "stages": stages, "library": a.library if a.domains else None}
    doc["proteins"] = load_proteins(con, root, manifest, stages, read)
    print("proteins", doc["proteins"], round(time.monotonic() - t0), flush=True)
    doc["clusters"] = load_clusters(con, root, manifest, read, doc["proteins"]["distinct"])
    doc["saving"] = saving(con, stages)
    print("saving done", round(time.monotonic() - t0), flush=True)
    if a.domains:
        con.execute(f"CREATE VIEW aln AS SELECT * FROM read_parquet({sql_list([Path(p) for p in doc['clusters']['aln_files']])})")
        doc["domains"] = load_domains(con, root, stages, a.library, read)
        print("domains", doc["domains"], round(time.monotonic() - t0), flush=True)
        idx = domain_index(con)
        # true pairs and representative pairs per scope (for the decision table)
        nd = {u: len(v[0]) for u, v in idx.items()}
        con.register("nd", pa.table({"uid": pa.array(list(nd), pa.binary(16)), "nd": pa.array(list(nd.values()), pa.int64())}))
        con.execute("CREATE TABLE ndt AS SELECT * FROM nd")
        con.unregister("nd")
        doc["inference"] = {}
        seq_files = [root / f["path"] for f in manifest["files"] if f["table"] == "protein_seq"]
        for level in LEVELS:
            doc.setdefault("alignment_check", {})[level] = check_alignments(con, seq_files, level, 5000)
            print(level, "alignment check", doc["alignment_check"][level], flush=True)
            r = infer(con, idx, stages, level)
            arr, cc = r["a"], r["coord"]
            lc = length_class(arr)
            pid = arr["pident"]
            L: dict = {"no_alignment": r["no_alignment"], "cigar_span_mismatch": r["cigar_span_mismatch"],
                       "members_longer_than_rep": int((arr["lm"] > arr["lr"]).sum()), "scopes": {}}
            for name, bits in scopes(stages):
                sel = (arr["mask"] & bits) != 0
                S: dict = {"all": summarize(arr, sel)}
                bins = {}
                for i in range(len(PIDENT_BINS) - 1):
                    lo, hi = PIDENT_BINS[i], PIDENT_BINS[i + 1]
                    m = (pid >= lo) & (pid < hi)
                    bins[f"[{lo:.2f},{hi:.2f})"] = summarize(arr, sel & m)
                bins["1.00"] = summarize(arr, sel & (pid == 1.0))
                bins["none"] = summarize(arr, sel & np.isnan(pid))
                S["by_pident"] = bins
                S["by_cluster_size"] = {f"{lo}-{hi}" if hi < 10**12 else f">={lo}":
                                        summarize(arr, sel & (arr["csize"] >= lo) & (arr["csize"] <= hi))
                                        for lo, hi in SIZE_BINS}
                S["by_length"] = {k: summarize(arr, sel & v) for k, v in lc.items()}
                csel = (cc["mask"] & bits) != 0
                S["coordinates"] = {"all": coord_summary(cc, csel),
                                    "architecture_equal": coord_summary(cc, csel & cc["arch_eq"]),
                                    "member_shorter_gt20": coord_summary(cc, csel & cc["short"]),
                                    "member_not_shorter_gt20": coord_summary(cc, csel & ~cc["short"])}
                # the decision: reps scanned, members inferred; and the hybrid (members whose length differs
                # from the rep's by more than 20% scanned too)
                tot = con.execute(
                    "SELECT count(*), sum(p.len), coalesce(sum(ndt.nd), 0), "
                    "coalesce(sum(p.len) FILTER (WHERE cm.m = cm.r), 0), "
                    "coalesce(sum(ndt.nd) FILTER (WHERE cm.m = cm.r), 0) "
                    "FROM cm JOIN prot p ON p.uid = cm.m LEFT JOIN ndt ON ndt.uid = cm.m "
                    "WHERE cm.level = ? AND (p.mask & ?) <> 0", [level, bits]).fetchone()
                n, res, true_pairs, rep_res, rep_pairs = (int(x) for x in tot)
                nm, nr, ni, lm = arr["n_m"][sel], arr["n_r"][sel], arr["n_i"][sel], arr["lm"][sel]
                hyb = (np.abs(arr["lm"] - arr["lr"]) > 0.2 * arr["lr"])[sel]
                dec = {}
                for mode, h in (("reps_only", np.zeros_like(hyb)), ("hybrid", hyb)):
                    scanned = rep_res + int(lm[h].sum())
                    right = rep_pairs + int(nm[h].sum()) + int(ni[~h].sum())
                    false = int((nr - ni)[~h].sum())
                    missed = int((nm - ni)[~h].sum())
                    assert right + missed == true_pairs, (right, missed, true_pairs)
                    dec[mode] = {"proteins": n, "residues": res, "true_pairs": true_pairs,
                                 "members_scanned_too": int(h.sum()), "residues_scanned": scanned,
                                 "residue_fraction_scanned": scanned / res,
                                 "pairs_right": right, "pairs_false": false, "pairs_missed": missed,
                                 "right_fraction_of_true": right / true_pairs,
                                 "missed_fraction_of_true": missed / true_pairs,
                                 "false_fraction_of_reported": frac(false, right + false),
                                 "false_per_true_pair": false / true_pairs,
                                 "residues_not_scanned_millions": (res - scanned) / 1e6,
                                 "cpu_min_saved_at_7.44": (res - scanned) / 1e6 * 7.44,
                                 "cpu_min_saved_at_7.96": (res - scanned) / 1e6 * 7.96}
                    if mode == "hybrid":
                        dec[mode]["remaining_members"] = summarize(arr, sel & ~(np.abs(arr["lm"] - arr["lr"]) > 0.2 * arr["lr"]))
                S["decision"] = dec
                L["scopes"][name] = S
            doc["inference"][level] = L
            print(level, "done", round(time.monotonic() - t0), flush=True)
    doc["files_read"] = read
    doc["elapsed_s"] = round(time.monotonic() - t0, 1)
    con.close()

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        return o

    tmpf = a.out.with_suffix(".json.tmp")
    tmpf.write_text(json.dumps(clean(doc), indent=1, sort_keys=True) + "\n")
    tmpf.replace(a.out)
    print("wrote", a.out, doc["elapsed_s"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

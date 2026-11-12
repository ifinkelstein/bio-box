"""Card 0-08b: tRNAs lost to Infernal's filters (M2) and the ncRNA producer on short phage contigs (M5).
Run ON THE POD, in the image (`ncrna_set_lane.sh` runs every step under `timed.py` in a capped unit).

M2, the tRNA (RF00005) and tRNA-Sec (RF01852) models alone on the 8 genomes of measurement 0-09c:

    python ncrna_set.py records   <wg.fna> <wg.map.json> <genomes,|all> <chunk> <overlap> <out.fna> <out.map.json>
    python ncrna_set.py collect   <tbl> <records.map.json> <wg.map.json> <out.jsonl>
    python ncrna_set.py m2compare <m09c dir> <trna.cm> <accs,> <default.jsonl> <out.json> <max.jsonl>...

`records` cuts the producer's whole-genome records of 0-09c (`fa/wg_all.fna`: linear contigs whole, circular ones
with a flank of 5,746 bases each side) into chunks of at most `chunk` bases overlapping by `overlap` bases (a chunk
larger than a record leaves it whole), so that each `cmsearch --max` process stays under the CPU budget. `collect`
maps hits back to the whole record and then to the contig as the producer does (flank copies dropped by the
leftmost-base rule, a hit longer than a circle rejected); a hit that touches an internal chunk edge is dropped (the
overlap, larger than any hit of these models, holds it whole in the next chunk), identical hits are kept once, and
overlapping hits of one model on one strand are resolved as cmsearch does within a record (best score kept).

M5, the producer's configuration (`biobox.l1.ncrna.run`, set prok-1) on a stride sample of a stage's contigs:

    python ncrna_set.py sample    <root> <stage> <target bases> <out.json>
    python ncrna_set.py m5run     <root> <stage> <sample.json> <part> <cm> <cpus> <out root> <tmp> <result.json>
    python ncrna_set.py m5regions <root> <stage> <sample.json> <igr_m200.map.json of 0-09c> <out.json>
    python ncrna_set.py m5report  <sample.json> <cmsearch record dir> <out.json> <out root>...

`sample`: the stage's `genome_input` contigs sorted by contig_uid, lengths from its `contig_index` (default
variant); with T the stage's bases and k = max(1, T // target), every contig whose 0-based rank i in that order has
i % k == 0. The sample is cut into two parts by cumulative bases in the same order (part 0 until half the
sample's bases is reached). `m5run` calls the producer's `run` with `genomes=` the part's contigs (one GenomeInput
per contig, its circular flag from genome_input) and `out=` a scratch root: it reads the stage's sequences and gene
calls from the root and writes nothing there. Coordinates are 0-based half-open, strand +1/-1.
"""

from __future__ import annotations

import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

TRNA = ("RF00005", "RF01852")


# --- M2 ------------------------------------------------------------------------------------------------------


def _read_fasta(path: Path):
    name, seq = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(seq)
                name, seq = line[1:].split()[0].decode(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, b"".join(seq)


def records(wg_fna: Path, wg_map: Path, genomes: str, chunk: int, overlap: int, out_fna: Path, out_map: Path) -> None:
    m = json.loads(wg_map.read_text())["map"]
    want = None if genomes == "all" else set(genomes.split(","))
    recmap, nrec, bases = {}, 0, 0
    tmp = out_fna.with_name(out_fna.name + ".tmp")
    with open(tmp, "wb") as fh:
        for name, seq in _read_fasta(wg_fna):
            if want is not None and m[name]["genome"] not in want:
                continue
            L = len(seq)
            starts = [0] if L <= chunk else list(range(0, L - overlap, chunk - overlap))
            for cs in starts:
                ce = min(L, cs + chunk)
                rn = f"{name}@{cs}"
                fh.write(b">" + rn.encode() + b"\n" + seq[cs:ce] + b"\n")
                recmap[rn] = {"rec": name, "cs": cs, "ce": ce, "rlen": L}
                nrec += 1
                bases += ce - cs
    tmp.replace(out_fna)
    doc = {"chunk": chunk, "overlap": overlap, "genomes": genomes, "records": nrec, "bases_searched": bases,
           "map": recmap}
    out_map.write_text(json.dumps(doc, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in doc.items() if k != "map"}, sort_keys=True))


def _parse_tblout(path: Path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split(None, 17)
            yield {"target": f[0], "query": f[2], "query_acc": f[3], "mdl_from": int(f[5]), "mdl_to": int(f[6]),
                   "seq_from": int(f[7]), "seq_to": int(f[8]), "strand": f[9], "trunc": f[10],
                   "score": float(f[14]), "evalue": float(f[15])}


def _greedy(hits: list[dict]) -> list[dict]:
    """Per (contig, model, strand), best score first, drop a hit that overlaps one kept (cmsearch's rule within a
    record)."""
    groups = defaultdict(list)
    for h in hits:
        groups[(h["contig"], h["acc"], h["strand"])].append(h)
    out = []
    for g in groups.values():
        kept = []
        for h in sorted(g, key=lambda x: (-x["score"], x["start"], x["end"])):
            if all(h["end"] <= k["start"] or k["end"] <= h["start"] for k in kept):
                kept.append(h)
        out += kept
    return sorted(out, key=lambda x: (x["contig"], x["start"], x["end"], x["strand"], x["acc"]))


def collect(tbl: Path, rec_map: Path, wg_map: Path, out: Path) -> None:
    cm = json.loads(rec_map.read_text())["map"]
    wm = json.loads(wg_map.read_text())["map"]
    seen, edge, flank_copies, rejected = {}, 0, 0, 0
    for h in _parse_tblout(tbl):
        c = cm[h["target"]]
        left, right = min(h["seq_from"], h["seq_to"]) - 1, max(h["seq_from"], h["seq_to"])
        if (left == 0 and c["cs"] > 0) or (right == c["ce"] - c["cs"] and c["ce"] < c["rlen"]):
            edge += 1
            continue
        left, right = left + c["cs"], right + c["cs"]
        r = wm[c["rec"]]
        n = r["len"]
        start, end = left - r["offset"], right - r["offset"]
        if r["offset"] and not 0 <= start < n:
            flank_copies += 1
            continue
        if r["offset"] and end - start > n:
            rejected += 1
            continue
        strand = 1 if h["strand"] == "+" else -1
        row = {"contig": r["contig"], "genome": r["genome"], "start": start, "end": end, "strand": strand,
               "acc": h["query_acc"], "name": h["query"], "score": h["score"], "evalue": h["evalue"],
               "trunc": h["trunc"], "mdl_from": h["mdl_from"], "mdl_to": h["mdl_to"]}
        key = (row["contig"], start, end, strand, row["acc"])
        if key not in seen or seen[key]["score"] < row["score"]:
            seen[key] = row
    rows = _greedy(list(seen.values()))
    with open(out, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
    print(json.dumps({"hits": len(rows), "distinct_before_overlap_rule": len(seen), "edge_dropped": edge,
                      "flank_copies": flank_copies, "exceeds_circle": rejected}))


def _hits(path: Path) -> list[dict]:
    return [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()]


def _ov(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def _cover(s: int, e: int, ivs: list[tuple[int, int]]) -> int:
    """Bases of [s, e) inside the union of `ivs` (linear coordinates; circular wrap handled by the caller)."""
    pts = sorted(ivs)
    union = []
    for a, b in pts:
        if union and a <= union[-1][1]:
            union[-1][1] = max(union[-1][1], b)
        else:
            union.append([a, b])
    return sum(_ov(s, e, a, b) for a, b in union)


def _wrap(ivs: list[tuple[int, int]], n: int) -> list[tuple[int, int]]:
    out = []
    for a, b in ivs:
        out.append((a, min(b, n)))
        if b > n:
            out.append((0, b - n))
    return out


def m2compare(m09c: Path, trna_cm: Path, accs: str, dpath: Path, out: Path, xpaths: list[Path]) -> None:
    """The default scan `dpath` against the --max scans `xpaths` (together covering the 8 genomes), both restricted to
    the models `accs` (comma list); the default scan is first checked against 0-09c's whole-genome truth."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import infernal_igr as ig
    from biobox.l1 import ncrna

    want = tuple(accs.split(","))
    models = ncrna.read_models(trna_cm)
    ga = {a: models[a].ga for a in want}
    clen = {a: models[a].clen for a in want}
    root = m09c / "root"
    meta = json.loads((root / "m09c_genomes.json").read_text())
    clen_contig = {c["uid"]: c["len"] for g in meta["genomes"].values() for c in g["contigs"]}
    genome_of = {c["uid"]: name for name, g in meta["genomes"].items() for c in g["contigs"]}
    D_all = _hits(dpath)
    truth = _hits(m09c / "out/wgA_default.hits.jsonl") + _hits(m09c / "out/wgB_default.hits.jsonl")
    truth_trna = sorted((t["contig"], t["start"], t["end"], t["strand"], t["acc"]) for t in truth if t["acc"] in TRNA)
    d_all_keys = sorted((h["contig"], h["start"], h["end"], h["strand"], h["acc"]) for h in D_all)
    default_equals_truth = d_all_keys == truth_trna
    D = [h for h in D_all if h["acc"] in want]
    X = [h for p in xpaths for h in _hits(p) if h["acc"] in want]
    d_keys = sorted((h["contig"], h["start"], h["end"], h["strand"], h["acc"]) for h in D)
    genes = ig.genes(root)
    igr_runs = ig.igrs(root)
    m200 = defaultdict(list)
    for r in json.loads((m09c / "fa/igr_m200.map.json").read_text())["map"].values():
        m200[r["contig"]].append((r["start"], r["end"]))

    def same_locus(h, pool):
        return [p for p in pool if p["contig"] == h["contig"] and p["strand"] == h["strand"]
                and _ov(h["start"], h["end"], p["start"], p["end"]) > 0]

    d_by_contig = defaultdict(list)
    for h in D:
        d_by_contig[h["contig"]].append(h)
    t_by_contig = defaultdict(list)
    for t in truth:
        t_by_contig[t["contig"]].append(t)
    dset = set(d_keys)
    x_exact = sum((h["contig"], h["start"], h["end"], h["strand"], h["acc"]) in dset for h in X)
    only_x = []
    per_genome = {g: {"default": 0, "max": 0, "max_only": 0, "default_only": 0} for g in sorted(meta["genomes"])}
    for h in D:
        per_genome[h["genome"]]["default"] += 1
        if not same_locus(h, [x for x in X if x["contig"] == h["contig"]]):
            per_genome[h["genome"]]["default_only"] += 1
    for h in X:
        g = h["genome"]
        per_genome[g]["max"] += 1
        if same_locus(h, d_by_contig[h["contig"]]):
            continue
        per_genome[g]["max_only"] += 1
        n = clen_contig[h["contig"]]
        L = h["end"] - h["start"]
        gcalls = [(a, b, s) for a, b, s in genes.get(h["contig"], [])]
        ovg = []
        for a, b, s in gcalls:
            o = _cover(h["start"], h["end"], _wrap([(a, b)], n))
            if o > 0:
                ovg.append({"start": a, "end": b, "strand": s, "len": b - a, "overlap": o,
                            "same_strand": s == h["strand"]})
        in_genes = _cover(h["start"], h["end"], _wrap([(a, b) for a, b, _ in gcalls], n))
        in_igr = _cover(h["start"], h["end"], _wrap(igr_runs.get(h["contig"], []), n))
        in_m200 = _cover(h["start"], h["end"], _wrap(m200.get(h["contig"], []), n))
        at_truth = [{"name": t["name"], "acc": t["acc"], "strand": t["strand"], "start": t["start"], "end": t["end"],
                     "score": t["score"]} for t in t_by_contig[h["contig"]]
                    if _ov(h["start"], h["end"], t["start"], t["end"]) > 0]
        mspan = abs(h["mdl_to"] - h["mdl_from"]) + 1
        only_x.append({**h, "len": L, "ga": ga[h["acc"]], "margin_over_ga": round(h["score"] - ga[h["acc"]], 1),
                       "model_cover": round(mspan / clen[h["acc"]], 3), "genes": ovg, "bases_in_genes": in_genes,
                       "bases_intergenic": in_igr, "bases_in_margin200_regions": in_m200,
                       "truth_hits_here_any_model_any_strand": at_truth})
    scores_d = sorted(h["score"] - ga[h["acc"]] for h in D)
    res = {"models": list(want), "ga": ga, "clen": clen, "default_hits": len(D), "max_hits": len(X),
           "default_scan_both_models_hits": len(D_all), "truth_trna_hits_both_models": len(truth_trna),
           "default_equals_0-09c_truth_restricted": default_equals_truth,
           "max_hits_at_default_coords": x_exact, "max_only_loci": len(only_x),
           "default_only_loci": sum(v["default_only"] for v in per_genome.values()),
           "default_margin_over_ga": {"min": scores_d[0] if scores_d else None,
                                      "below_10_bits": sum(s < 10 for s in scores_d),
                                      "median": statistics.median(scores_d) if scores_d else None},
           "per_genome": per_genome, "max_only": only_x, "genome_of_contigs": len(genome_of)}
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "max_only"}, sort_keys=True))


# --- M5 ------------------------------------------------------------------------------------------------------


def _part_files(root: Path, rel: str) -> list[Path]:
    from biobox import io

    return io.part_files(Path(root) / rel)


def sample(root: Path, stage: str, target: int, out: Path) -> None:
    import duckdb

    from biobox.l0seq.store import INDEX_DIR
    from biobox.registry.schemas import table

    gi = _part_files(root, f"{table('genome_input').location}/stage={stage}")
    idx = _part_files(root, f"{INDEX_DIR}/stage={stage}")
    con = duckdb.connect(config={"threads": 1, "memory_limit": "2GB"})
    files = lambda fs: "[" + ", ".join("'" + str(f) + "'" for f in fs) + "]"  # noqa: E731
    rows = con.sql(
        f"SELECT g.contig_uid, g.is_circular, g.kind, g.genome_uid, c.contig_len FROM (SELECT DISTINCT contig_uid, "
        f"is_circular, kind, genome_uid FROM read_parquet({files(gi)}, hive_partitioning = false)) g JOIN (SELECT "
        f"contig_uid, max(contig_len) AS contig_len FROM read_parquet({files(idx)}, hive_partitioning = false) WHERE "
        "is_default_variant GROUP BY contig_uid) c ON c.contig_uid = g.contig_uid ORDER BY g.contig_uid").fetchall()
    con.close()
    total = sum(r[4] for r in rows)
    k = max(1, total // target)
    picked = [r for i, r in enumerate(rows) if i % k == 0]
    sb = sum(r[4] for r in picked)
    cum, out_rows = 0, []
    for cu, circ, kind, guid, n in picked:
        part = 0 if cum < sb / 2 else 1
        cum += n
        out_rows.append({"contig_uid": bytes(cu).hex(), "circular": bool(circ), "kind": kind,
                         "genome_uid": bytes(guid).hex(), "len": int(n), "part": part})
    lens = [r[4] for r in rows]
    doc = {"stage": stage, "rule": f"contigs of genome_input sorted by contig_uid; rank i % {k} == 0",
           "stage_contigs": len(rows), "stage_bases": total, "stage_genomes": len({bytes(r[3]) for r in rows}),
           "stage_kinds": dict(Counter(r[2] for r in rows)), "stage_circular": sum(bool(r[1]) for r in rows),
           "stage_len_median": statistics.median(lens), "target": target, "stride": k,
           "contigs": len(out_rows), "bases": sb,
           "parts": {p: {"contigs": sum(r["part"] == p for r in out_rows),
                         "bases": sum(r["len"] for r in out_rows if r["part"] == p)} for p in (0, 1)},
           "sample": out_rows}
    out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k2: v for k2, v in doc.items() if k2 != "sample"}, sort_keys=True))


def m5run(root: Path, stage: str, sample_json: Path, part: int, cm: Path, cpus: int, out: Path, tmp: Path,
          result: Path) -> None:
    from biobox.l0seq.inputs import GenomeInput
    from biobox.l1 import ncrna

    doc = json.loads(sample_json.read_text())
    genomes = [GenomeInput((bytes.fromhex(r["contig_uid"]),), kind=r["kind"],
                           circular=(bytes.fromhex(r["contig_uid"]),) if r["circular"] else ())
               for r in doc["sample"] if r["part"] == part]
    res = ncrna.run(root, stage, cm, cpus, genomes, out=out, tmp=tmp, seq_stages=[stage])
    res = {**res, "part": part, "genomes_given": len(genomes)}
    result.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(json.dumps(res, sort_keys=True))


def _dist(xs: list[int]) -> dict:
    if not xs:
        return {"n": 0}
    s = sorted(xs)
    q = lambda f: s[min(len(s) - 1, int(f * len(s)))]  # noqa: E731
    bins = [0, 100, 300, 500, 1000, 2000, 5000, 10000, 10 ** 12]
    hist = {f"{a}-{b}" if b < 10 ** 12 else f">={a}": sum(a <= x < b for x in s) for a, b in zip(bins, bins[1:])}
    return {"n": len(s), "sum": sum(s), "mean": round(sum(s) / len(s), 1), "min": s[0], "p10": q(0.1),
            "p25": q(0.25), "median": q(0.5), "p75": q(0.75), "p90": q(0.9), "max": s[-1], "hist": hist}


def m5regions(root: Path, stage: str, sample_json: Path, m09c_map: Path, out: Path) -> None:
    from biobox.l1 import genecall, ncrna

    doc = json.loads(sample_json.read_text())
    want = {bytes.fromhex(r["contig_uid"]): r for r in doc["sample"]}
    genes = defaultdict(list)
    files = ncrna._gene_files(root, stage, genecall.GENE_CALL_VERSION)
    for cu, s, e in ncrna._range_rows(files, ["contig_uid", "start", "end"]):
        if cu in want:
            genes[cu].append((s, e))
    reg_lens, per_contig, no_region, whole_circle, ngenes = [], [], 0, 0, []
    for cu, r in sorted(want.items()):
        regs = ncrna.regions(genes.get(cu, []), r["len"], r["circular"])
        ngenes.append(len(genes.get(cu, [])))
        per_contig.append(len(regs))
        no_region += not regs
        whole_circle += any(r["circular"] and e - s == r["len"] for s, e in regs)
        reg_lens += [e - s for s, e in regs]
    sb = doc["bases"]
    m = json.loads(m09c_map.read_text())
    ref = [v["end"] - v["start"] for v in m["map"].values()]
    res = {"sample_contigs": len(want), "sample_bases": sb, "contig_len": _dist([r["len"] for r in want.values()]),
           "genes_per_contig": _dist(ngenes), "genes": sum(ngenes), "regions": _dist(reg_lens),
           "regions_per_contig": _dist(per_contig), "contigs_without_region": no_region,
           "contigs_one_whole_circle_region": whole_circle,
           "fraction_bases_in_regions": round(sum(reg_lens) / sb, 4),
           "records_per_mbp": round(len(reg_lens) / sb * 1e6, 1),
           "m09c_8_genomes": {"genome_bases": m["genome_bases"], "regions": _dist(ref),
                              "fraction_bases_in_regions": round(sum(ref) / m["genome_bases"], 4),
                              "records_per_mbp": round(len(ref) / m["genome_bases"] * 1e6, 1)}}
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(json.dumps(res, sort_keys=True))


def m5report(sample_json: Path, rec_dir: Path, out: Path, roots: list[Path]) -> None:
    import pyarrow.parquet as pq

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import infernal_igr as ig

    doc = json.loads(sample_json.read_text())
    rows = []
    for r in roots:
        for f in sorted(Path(r).glob("l1/ncrna/producer=*/stage=*/*.parquet")):
            rows += pq.read_table(f, columns=["contig_uid", "start", "end", "strand", "model_acc", "model_name",
                                              "bitscore", "clan", "clan_winner", "trunc_5p", "trunc_3p"],
                                  partitioning=None).to_pylist()
    fam = Counter(r["model_name"] for r in rows)
    fam_win = Counter(r["model_name"] for r in rows if r["clan_winner"])
    contigs_with_hit = len({r["contig_uid"] for r in rows})
    calls = []
    for p in sorted(Path(rec_dir).glob("*.json")):
        rec = json.loads(p.read_text())
        kind = "intergenic" if "--notrunc" in rec["argv"] else "whole"
        outp = Path(str(p)[:-5] + ".out")
        pm = None
        if outp.is_file():
            tmpj = Path(str(p)[:-5] + ".permodel.json")
            ig.permodel(outp, tmpj)
            pm = json.loads(tmpj.read_text())
        calls.append({"record": p.name, "kind": kind, "cpu_min": rec["cpu_min"], "wall_s": rec["wall_s"],
                      "peak_rss_mib": rec["peak_rss_mib"], "exit_code": rec["exit_code"],
                      "top_models": None if pm is None else sorted(
                          ({"name": m["name"], "acc": m["acc"], "cpu_s": m["cpu_s"]} for m in pm),
                          key=lambda m: -(m["cpu_s"] or 0))[:12],
                      "models": None if pm is None else len(pm),
                      "permodel_cpu_s": None if pm is None else round(sum(m["cpu_s"] or 0 for m in pm), 1)})
    by_kind = defaultdict(float)
    for c in calls:
        by_kind[c["kind"]] += c["cpu_min"]
    res = {"sample_bases": doc["bases"], "sample_contigs": doc["contigs"], "hits": len(rows),
           "winners": sum(bool(r["clan_winner"]) for r in rows), "losers": sum(not r["clan_winner"] for r in rows),
           "families": len(fam), "contigs_with_a_hit": contigs_with_hit,
           "truncated_hits": sum(bool(r["trunc_5p"] or r["trunc_3p"]) for r in rows),
           "hits_by_family": dict(sorted(fam.items(), key=lambda kv: (-kv[1], kv[0]))),
           "winners_by_family": dict(sorted(fam_win.items(), key=lambda kv: (-kv[1], kv[0]))),
           "cmsearch_calls": calls, "cmsearch_cpu_min_by_kind": {k: round(v, 3) for k, v in sorted(by_kind.items())}}
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("hits_by_family", "winners_by_family",
                                                                  "cmsearch_calls")}, sort_keys=True))


def main(argv: list[str]) -> int:
    cmd, a = argv[0], argv[1:]
    if cmd == "records":
        records(Path(a[0]), Path(a[1]), a[2], int(a[3]), int(a[4]), Path(a[5]), Path(a[6]))
    elif cmd == "collect":
        collect(Path(a[0]), Path(a[1]), Path(a[2]), Path(a[3]))
    elif cmd == "m2compare":
        m2compare(Path(a[0]), Path(a[1]), a[2], Path(a[3]), Path(a[4]), [Path(x) for x in a[5:]])
    elif cmd == "sample":
        sample(Path(a[0]), a[1], int(a[2]), Path(a[3]))
    elif cmd == "m5run":
        m5run(Path(a[0]), a[1], Path(a[2]), int(a[3]), Path(a[4]), int(a[5]), Path(a[6]), Path(a[7]), Path(a[8]))
    elif cmd == "m5regions":
        m5regions(Path(a[0]), a[1], Path(a[2]), Path(a[3]), Path(a[4]))
    elif cmd == "m5report":
        m5report(Path(a[0]), Path(a[1]), Path(a[2]), [Path(x) for x in a[3:]])
    else:
        print(__doc__, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

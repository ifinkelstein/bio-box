"""Card 0-09c: Infernal on intergenic regions only, against whole-genome scanning (run ON THE POD, in the image).

    python infernal_igr.py prepare  <repo> <golden dir> <root> <tmp>
    python infernal_igr.py fasta    <root> <cm> wg  <batch> <out.fna> <out.map.json>
    python infernal_igr.py fasta    <root> <cm> igr <batch> <out.fna> <out.map.json> <margin>
    python infernal_igr.py tiny     <out.fna>                     (the fixed-cost probe: 1,000 seeded random bases)
    python infernal_igr.py collect  <tbl> <map.json> <out.jsonl>
    python infernal_igr.py subset   <repo> <truth.jsonl...> -- <out.keys>
    python infernal_igr.py permodel <cmsearch -o file> <out.json>
    python infernal_igr.py compare  <root> <truth.jsonl> <variant.jsonl> <scanned map.json or -> <out.json>

`batch` is A, B or all (GENOMES below). Every number this prints is measured on its inputs; nothing is
estimated here. The cmsearch calls themselves are made by lane.sh (under timed.py), with the producer's
argv (`biobox.l1.ncrna._cmsearch_argv`: --cpu, -Z 1e5, --cut_ga, --noali, --tformat fasta, --tblout) plus
`-o <file>` instead of /dev/null, so the per-model CPU of every run is kept; a variant appends its flags.

Whole genome (`wg`) records are built and their hits mapped exactly as the producer does it
(`ncrna._extended`, the library's own flank, the leftmost-base rule, `exceeds_circle` rejection).
Intergenic (`igr`) records are the stage's IGR rows (`l0seq/igr`, every run of >= 1 nt, gene-call
producer only), each widened by `margin` bases into the flanking genes on both sides, clipped at the ends
of a linear contig, taken modulo the length on a circular one, and merged where they overlap or touch
(also across the origin). Each record is cut from the + strand; cmsearch searches both strands; a hit
maps back by the record's start (start >= contig_len is reduced by contig_len, so a hit crossing the
origin has end > contig_len, AGENTS §3.1).

Coordinates throughout are 0-based half-open, strand +1/-1.
"""

from __future__ import annotations

import importlib.util
import json
import os
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

GENOMES = {  # batch: golden ids; each batch is one cmsearch process of about 9 Mbp (under 300 CPU-minutes)
    "A": ["ecoli_mg1655", "bsub_168", "lambda", "t4"],
    "B": ["mtb_h37rv", "spyogenes_sf370", "sacidocaldarius_dsm639", "lak_a1"],
}
STAGE = "M9C"
MIN_OVERLAP = 0.5  # a truth hit is recovered by a hit of its own model on its strand covering half its length


def _genomes_json(root: Path) -> dict:
    return json.loads((Path(root) / "m09c_genomes.json").read_text())


# --- prepare ---------------------------------------------------------------------------------------------


def prepare(repo: Path, gdir: Path, root: Path, tmp: Path) -> None:
    import yaml

    from biobox import fasta, ids
    from biobox.l0seq import SourceObject, ingest
    from biobox.l0seq.inputs import GenomeInput, write_genome_input
    from biobox.l1 import genecall, igr

    golden = yaml.safe_load((repo / "golden.yaml").read_text())
    lock = yaml.safe_load((repo / "golden.lock.yaml").read_text())
    entries = {e["id"]: e for e in golden["genomes"]}
    objs, genomes, info = [], [], {}
    for batch, names in GENOMES.items():
        for name in names:
            e = entries[name]
            path = gdir / lock["files"][name]["path"]
            recs = fasta.read(path)
            uids = tuple(ids.contig_uid(r.seq) for r in recs)
            circ = tuple(u for u, spec in zip(uids, e["records"]) if spec["topology"] == "circular")
            objs.append(SourceObject(path, "golden", str(golden["version"])))
            genomes.append(GenomeInput(uids, e["genetic_code"], kind=e["kind"], circular=circ))
            gc = sum(r.seq.upper().count(b"G") + r.seq.upper().count(b"C") for r in recs)
            n = sum(len(r.seq) for r in recs)
            info[name] = {"batch": batch, "kind": e["kind"], "genetic_code": e["genetic_code"], "bases": n,
                          "gc": round(gc / n, 4),
                          "contigs": [{"uid": u.hex(), "len": len(r.seq), "circular": u in circ}
                                      for u, r in zip(uids, recs)]}
    ingest(objs, root, STAGE)
    write_genome_input(root, STAGE, genomes)
    g = genecall.run(root, STAGE, workers=2, tmp=tmp)
    i = igr.run(root, STAGE, tmp=tmp)
    out = {"stage": STAGE, "genomes": info, "genecall": {k: v for k, v in g.items() if isinstance(v, (int, str))},
           "igr": i, "gene_call_version": genecall.GENE_CALL_VERSION, "igr_version": igr.igr_version()}
    (root / "m09c_genomes.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(json.dumps(out, sort_keys=True))


# --- tables of the store ---------------------------------------------------------------------------------


def _read_rows(root: Path, rel: str, cols: list[str]) -> list[dict]:
    import pyarrow.parquet as pq

    from biobox import io

    rows = []
    for f in io.part_files(Path(root) / rel):
        rows += pq.read_table(f, columns=cols, partitioning=None).to_pylist()
    return rows


def genes(root: Path) -> dict[str, list[tuple[int, int, int]]]:
    meta = _genomes_json(root)
    rows = _read_rows(root, f"l1/feature/producer={meta['gene_call_version']}/stage={STAGE}",
                      ["contig_uid", "start", "end", "strand"])
    out = defaultdict(list)
    for r in rows:
        out[r["contig_uid"].hex()].append((r["start"], r["end"], r["strand"]))
    return {k: sorted(v) for k, v in out.items()}


def igrs(root: Path) -> dict[str, list[tuple[int, int]]]:
    meta = _genomes_json(root)
    rows = _read_rows(root, f"l0seq/igr/igr_version={meta['igr_version']}/stage={STAGE}",
                      ["contig_uid", "start", "end"])
    out = defaultdict(list)
    for r in rows:
        out[r["contig_uid"].hex()].append((r["start"], r["end"]))
    return {k: sorted(v) for k, v in out.items()}


# --- regions ---------------------------------------------------------------------------------------------


def merge_regions(ivs: list[tuple[int, int]], n: int, circular: bool, margin: int) -> list[tuple[int, int]]:
    """Widen each [s, e) by `margin` on both sides and merge overlapping or touching runs. Linear: clipped
    to [0, n). Circular: 0 <= s < n and s < e <= s + n; a run crossing the origin has e > n. A circle
    covered whole is one run (0, n)."""
    if not ivs:
        return []
    if not circular:
        w = sorted((max(0, s - margin), min(n, e + margin)) for s, e in ivs)
        out = [list(w[0])]
        for s, e in w[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        return [tuple(x) for x in out]
    w = []
    for s, e in ivs:
        s, e = s - margin, e + margin
        if e - s >= n:
            return [(0, n)]
        s0 = s % n
        w.append((s0, s0 + (e - s)))
    w.sort()
    out = [list(w[0])]
    for s, e in w[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    # across the origin: the last run may reach into the first ones (their coordinates + n)
    while len(out) > 1 and out[-1][1] >= out[0][0] + n:
        first = out.pop(0)
        out[-1][1] = max(out[-1][1], first[1] + n)
    if out[-1][1] - out[-1][0] >= n:
        return [(0, n)]
    return [tuple(x) for x in out]


def write_fasta(root: Path, cm: Path, mode: str, batch: str, out_fna: Path, out_map: Path, margin: int = 0) -> None:
    from biobox.l0seq.store import SeqStore
    from biobox.l1 import ncrna

    meta = _genomes_json(root)
    names = [g for g, v in meta["genomes"].items() if batch == "all" or v["batch"] == batch]
    flank = max(max(m.w, m.maxl) for m in ncrna.read_models(cm).values()) if mode == "wg" else 0
    regions = igrs(root) if mode == "igr" else None
    store = SeqStore(root, [STAGE])
    recmap, bases, nrec, genome_bases = {}, 0, 0, 0
    tmp = out_fna.with_name(out_fna.name + ".tmp")
    try:
        with open(tmp, "wb") as fh:
            for g in sorted(names):
                for c in meta["genomes"][g]["contigs"]:
                    cu, n, circ = c["uid"], c["len"], c["circular"]
                    seq = store.get(bytes.fromhex(cu)).upper()
                    assert len(seq) == n
                    genome_bases += n
                    if mode == "wg":
                        rec = ncrna._extended(seq, flank) if circ else seq
                        fh.write(b">" + cu.encode() + b"\n" + rec + b"\n")
                        recmap[cu] = {"contig": cu, "genome": g, "len": n, "circular": circ, "offset": flank if circ else 0}
                        bases += len(rec)
                        nrec += 1
                        continue
                    ring = seq + seq if circ else seq
                    for s, e in merge_regions(regions.get(cu, []), n, circ, margin):
                        name = f"{cu}_{s}"
                        fh.write(b">" + name.encode() + b"\n" + ring[s:e] + b"\n")
                        recmap[name] = {"contig": cu, "genome": g, "len": n, "circular": circ, "start": s, "end": e}
                        bases += e - s
                        nrec += 1
        os.replace(tmp, out_fna)
    finally:
        store.close()
    doc = {"mode": mode, "batch": batch, "margin": margin, "flank": flank, "records": nrec, "bases_searched": bases,
           "genome_bases": genome_bases, "genomes": sorted(names), "map": recmap}
    out_map.write_text(json.dumps(doc, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in doc.items() if k != "map"}, sort_keys=True))


def tiny(out_fna: Path) -> None:
    rng = random.Random(20261111)
    out_fna.write_bytes(b">tiny\n" + bytes(rng.choices(b"ACGT", k=1000)) + b"\n")


# --- hits ------------------------------------------------------------------------------------------------


def collect(tbl: Path, map_path: Path, out: Path) -> None:
    """The run's hits in genome coordinates, one JSON object per line, sorted; the producer's rules for
    whole-genome records (flank copies dropped, a hit longer than the circle rejected)."""
    from biobox.l1 import ncrna

    doc = json.loads(map_path.read_text())
    m = doc["map"]
    rows, flank_copies, rejected = set(), 0, 0
    for h in ncrna._parse_tblout(tbl):
        r = m[h["target"]]
        n = r["len"]
        left, right = min(h["seq_from"], h["seq_to"]) - 1, max(h["seq_from"], h["seq_to"])
        if doc["mode"] == "wg":
            start, end = left - r["offset"], right - r["offset"]
            if r["offset"] and not 0 <= start < n:
                flank_copies += 1
                continue
            if r["offset"] and end - start > n:
                rejected += 1
                continue
        else:
            start, end = r["start"] + left, r["start"] + right
            if r["circular"] and start >= n:
                start, end = start - n, end - n
        strand = 1 if h["strand"] == "+" else -1
        rows.add((r["contig"], start, end, strand, h["query_acc"], h["query"], h["score"], h["evalue"], h["trunc"],
                  r["genome"]))
    keys = ("contig", "start", "end", "strand", "acc", "name", "score", "evalue", "trunc", "genome")
    with open(out, "w") as fh:
        for row in sorted(rows):
            fh.write(json.dumps(dict(zip(keys, row)), sort_keys=True) + "\n")
    print(json.dumps({"hits": len(rows), "flank_copies": flank_copies, "exceeds_circle": rejected,
                      "families": len({r[4] for r in rows})}))


def read_hits(path: Path) -> list[dict]:
    return [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()]


def g04_module(repo: Path):
    """The E-02 spike-in test module itself (its PANEL and helpers are imported, not copied)."""
    spec = importlib.util.spec_from_file_location("g04_spikeins", repo / "tests/r3/test_g04_spikeins.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def subset(repo: Path, truth_paths: list[Path], out_keys: Path) -> None:
    accs = {h["acc"] for p in truth_paths for h in read_hits(p)}
    panel = set(g04_module(repo).PANEL)
    keys = sorted(accs | panel)
    out_keys.write_text("".join(k + "\n" for k in keys))
    print(json.dumps({"truth_families": len(accs), "panel": len(panel), "panel_not_in_truth": sorted(panel - accs),
                      "models": len(keys)}))


_Q = re.compile(r"^Query:\s+(\S+)\s+\[CLEN=(\d+)\]")
_ACC = re.compile(r"^Accession:\s+(\S+)")
_CPU = re.compile(r"^# CPU time: ([\d.]+)u ([\d.]+)s")


def permodel(out_path: Path, out_json: Path) -> None:
    """Per model of one cmsearch -o file: name, accession, CLEN, CPU seconds (user + system, the
    '# CPU time' line Infernal prints for each query) and whether the model ran in HMM-only mode."""
    models, cur = [], None
    with open(out_path, errors="replace") as fh:
        for line in fh:
            q = _Q.match(line)
            if q:
                cur = {"name": q.group(1), "clen": int(q.group(2)), "acc": None, "cpu_s": None, "hmmonly": False}
                models.append(cur)
                continue
            if cur is None:
                continue
            a = _ACC.match(line)
            if a:
                cur["acc"] = a.group(1)
            elif line.startswith("Description:"):
                cur["desc"] = line.split(":", 1)[1].strip()
            elif line.startswith("Internal HMM-only pipeline"):
                cur["hmmonly"] = True
            c = _CPU.match(line)
            if c:
                cur["cpu_s"] = float(c.group(1)) + float(c.group(2))
    out_json.write_text(json.dumps(models, sort_keys=True) + "\n")
    tot = sum(m["cpu_s"] or 0 for m in models)
    print(json.dumps({"models": len(models), "cpu_s_sum": round(tot, 2),
                      "hmmonly_models": sum(m["hmmonly"] for m in models)}))


# --- comparison ------------------------------------------------------------------------------------------


def _pieces(s: int, e: int, n: int) -> list[tuple[int, int]]:
    """[s, e) on a circle of length n (0 <= s < n, e <= s + n) as at most two linear pieces."""
    return [(s, min(e, n))] + ([(0, e - n)] if e > n else [])


def _cover(s: int, e: int, n: int, ivs: list[tuple[int, int]]) -> int:
    """Bases of [s, e) that lie in the union of the runs `ivs` (on the same circle or line)."""
    pts = []
    for a, b in ivs:
        pts += _pieces(a, b, n)
    pts.sort()
    union = []
    for a, b in pts:
        if union and a <= union[-1][1]:
            union[-1][1] = max(union[-1][1], b)
        else:
            union.append([a, b])
    tot = 0
    for a, b in _pieces(s, e, n):
        for u0, u1 in union:
            tot += max(0, min(b, u1) - max(a, u0))
    return tot


def _ov(a0: int, a1: int, b0: int, b1: int, n: int) -> int:
    return _cover(a0, a1, n, [(b0, b1)])


def compare(root: Path, truth_path: Path, var_path: Path, scanned_map: str, out: Path) -> None:
    meta = _genomes_json(root)
    clen = {c["uid"]: c["len"] for g in meta["genomes"].values() for c in g["contigs"]}
    truth, var = read_hits(truth_path), read_hits(var_path)
    by_key = defaultdict(list)
    for v in var:
        by_key[(v["contig"], v["acc"], v["strand"])].append(v)

    def match(t, vs):
        n = clen[t["contig"]]
        need = MIN_OVERLAP * (t["end"] - t["start"])
        return [v for v in vs if _ov(t["start"], t["end"], v["start"], v["end"], n) >= need]

    scanned = None
    if scanned_map != "-":
        d = json.loads(Path(scanned_map).read_text())
        scanned = defaultdict(list)
        for r in d["map"].values():
            scanned[r["contig"]].append((r["start"], r["end"]))
    gene_runs = genes(root)
    matched_var = set()
    lost, found = [], 0
    fam_truth, fam_found = defaultdict(int), defaultdict(int)
    exact = 0
    for t in truth:
        fam_truth[t["acc"]] += 1
        ms = match(t, by_key[(t["contig"], t["acc"], t["strand"])])
        if ms:
            found += 1
            fam_found[t["acc"]] += 1
            exact += any(v["start"] == t["start"] and v["end"] == t["end"] for v in ms)
            for v in ms:
                matched_var.add(id(v))
            continue
        n = clen[t["contig"]]
        L = t["end"] - t["start"]
        g = gene_runs.get(t["contig"], [])
        in_genes = _cover(t["start"], t["end"], n, [(a, b) for a, b, _ in g])
        ovg = [(a, b, s) for a, b, s in g if _ov(t["start"], t["end"], a, b, n) > 0]
        rec = {**t, "len": L, "gene_bases": in_genes,
               "genes": [{"start": a, "end": b, "strand": s, "len": b - a, "same_strand": s == t["strand"]} for a, b, s in ovg]}
        same = [v for v in by_key[(t["contig"], t["acc"], t["strand"])] if _ov(t["start"], t["end"], v["start"], v["end"], n) > 0]
        rec["partial_hit"] = [{"start": v["start"], "end": v["end"], "score": v["score"]} for v in same]
        if scanned is not None:
            cov = _cover(t["start"], t["end"], n, scanned.get(t["contig"], []))
            rec["scanned_bases"] = cov
            rec["why"] = "not scanned (inside genes)" if cov == 0 else ("cut at a region boundary" if cov < L
                                                                         else "scanned whole, not reported")
        else:
            rec["why"] = "not reported (same sequence, other settings)"
        lost.append(rec)
    extra = [v for v in var if id(v) not in matched_var]
    fam_var = {v["acc"] for v in var}
    fam_rec = {a for a in fam_truth if fam_found.get(a)}
    fam_var_true = {v["acc"] for v in var if id(v) in matched_var}
    res = {
        "truth_hits": len(truth), "variant_hits": len(var), "recovered": found, "recovered_exact_coords": exact,
        "recall_hits": round(found / len(truth), 4) if truth else None,
        "precision_hits": round((len(var) - len(extra)) / len(var), 4) if var else None,
        "truth_families": len(fam_truth), "families_recovered": len(fam_rec),
        "families_lost_entirely": sorted(a for a in fam_truth if not fam_found.get(a)),
        "recall_families": round(len(fam_rec) / len(fam_truth), 4) if fam_truth else None,
        "variant_families": len(fam_var), "precision_families": round(len(fam_var_true) / len(fam_var), 4) if fam_var else None,
        "per_family": {a: {"truth": fam_truth[a], "found": fam_found.get(a, 0),
                           "name": next(t["name"] for t in truth if t["acc"] == a)} for a in sorted(fam_truth)},
        "lost": lost,
        "extra": extra,
    }
    why = defaultdict(int)
    for r in lost:
        why[r["why"]] += 1
    res["lost_by_reason"] = dict(sorted(why.items()))
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("lost", "extra", "per_family")}, sort_keys=True))


def e01_hits(root: Path, e01_root: Path, out: Path) -> None:
    """The E-01 producer's own `l1_ncrna` rows (card E-01 run, 2026-11-10, read only) for this store's contigs,
    in `collect`'s format: an independent run of the producer to check this script's whole-genome runner."""
    import pyarrow.parquet as pq

    meta = _genomes_json(root)
    want = {c["uid"]: g for g, v in meta["genomes"].items() for c in v["contigs"]}
    rows = set()
    for f in sorted(Path(e01_root).glob("l1/ncrna/producer=*/stage=*/part-*.parquet")):
        for r in pq.read_table(f, partitioning=None).to_pylist():
            cu = r["contig_uid"].hex()
            if cu in want:
                tr = ("5'" if r["trunc_5p"] else "") + ("&" if r["trunc_5p"] and r["trunc_3p"] else "") + \
                     ("3'" if r["trunc_3p"] else "")
                rows.add((cu, r["start"], r["end"], r["strand"], r["model_acc"], r["model_name"], r["bitscore"],
                          r["evalue_zref"], tr or "no", want[cu]))
    keys = ("contig", "start", "end", "strand", "acc", "name", "score", "evalue", "trunc", "genome")
    with open(out, "w") as fh:
        for row in sorted(rows):
            fh.write(json.dumps(dict(zip(keys, row)), sort_keys=True) + "\n")
    print(json.dumps({"hits": len(rows), "contigs": len({r[0] for r in rows})}))


def main(argv: list[str]) -> int:
    cmd, a = argv[0], argv[1:]
    if cmd == "e01":
        e01_hits(Path(a[0]), Path(a[1]), Path(a[2]))
        return 0
    if cmd == "prepare":
        prepare(Path(a[0]), Path(a[1]), Path(a[2]), Path(a[3]))
    elif cmd == "fasta":
        write_fasta(Path(a[0]), Path(a[1]), a[2], a[3], Path(a[4]), Path(a[5]), int(a[6]) if len(a) > 6 else 0)
    elif cmd == "tiny":
        tiny(Path(a[0]))
    elif cmd == "collect":
        collect(Path(a[0]), Path(a[1]), Path(a[2]))
    elif cmd == "subset":
        i = a.index("--")
        subset(Path(a[0]), [Path(p) for p in a[1:i]], Path(a[i + 1]))
    elif cmd == "permodel":
        permodel(Path(a[0]), Path(a[1]))
    elif cmd == "compare":
        compare(Path(a[0]), Path(a[1]), Path(a[2]), a[3], Path(a[4]))
    else:
        print(__doc__, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

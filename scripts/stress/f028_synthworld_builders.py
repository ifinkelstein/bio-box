"""Card F-02 session 8, part B: the emitted `micro` synthworld through the normal builders, compared with its truth.

    python f028_synthworld_builders.py world OUT
    python f028_synthworld_builders.py inputs WORLD OUT_DIR
    python f028_synthworld_builders.py compare WORLD ROOT RELEASE OUT_DIR

world: builds S1-S11 of the `micro` preset at its default seed into OUT (synthworld, the code of this tree) and
prints one JSON line per stage with its wall time, then the total.
inputs: the `bio-box build-stage` inputs of two stages from WORLD/emit (its ingest.tsv and FASTA): stage F0 holds
the original sources (src_a: every host genome, its contaminants, and the phages), stage F1 the duplicate sources
(src_b exact, src_c near), so an exact duplicate's contigs are first stored by F0 and only named by F1 (DEC-26);
within a stage a genome whose contig set an earlier source listed is listed once and a contig in several genomes
is listed in each (DEC-27, make_inputs.py's Listing). Writes OUT_DIR/genomes-<S>.json and objects-<S>.json
(object paths absolute) and prints a JSON report.
compare: reads under ROOT the run records of both stages and of RELEASE and the Parquet files they list (never a
glob, AGENTS 3.9), and compares with WORLD/emit/truth: every emitted contig stored once with its length and
digest, exact duplicates collapsed to one contig_uid and one genome, near duplicates kept apart, F1's calls of
F0's contigs skipped (DEC-26); gene calls against truth_gene by the stop-anchored match of synthworld_spec 3.1
(contig_uid, strand, stop_bp), recall and start agreement by class, a gene piece cut at its 3' end by a contig
edge matched by a call in its frame that overlaps half of it, and precision of the calls by contig origin.
Writes OUT_DIR/recall_by_class.tsv, precision_by_origin.tsv, steps.tsv and prints a JSON report.
"""

from __future__ import annotations

import collections
import csv
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import pyarrow.parquet as pq

USAGE = __doc__.split("\n\n")[1]


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def world(out: str) -> int:
    from biobox.synthworld import assembly, backbone, emit, fusions, intergenic, mobile, sequences, sources

    out_p = Path(out)
    params = backbone.load_preset("micro")
    total = 0.0
    for name, fn in (("S1-S2 backbone", lambda: backbone.build(params, params["default_seed"], out_p)),
                     ("S3-S4 sequences", lambda: sequences.build(out_p)),
                     ("S5 replicons", lambda: intergenic.build(out_p)),
                     ("S6 mobile", lambda: mobile.build(out_p)),
                     ("S7 fusions", lambda: fusions.build(out_p)),
                     ("S8-S9 assembly", lambda: assembly.build(out_p)),
                     ("S10 sources", lambda: sources.build(out_p)),
                     ("S11 emit", lambda: emit.build(out_p))):
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        total += dt
        print(json.dumps({"stage": name, "wall_s": round(dt, 2)}), flush=True)
    man = json.loads((out_p / "manifest.json").read_text())
    print(json.dumps({"stage": "S1-S11", "wall_s": round(total, 2), "world_digest": man["world_digest"],
                      "ru_maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)}))
    return 0


def _listing(sources: list[tuple], fasta_dir: Path) -> tuple[list[dict], list[dict], dict]:
    """GenomeInput entries and SourceObjects of one stage's sources, DEC-27's Listing rule."""
    from biobox import ids

    genomes, objects, seen, rep = [], [], set(), collections.Counter()
    for (_lab, _g), rows in sources:
        uids = list(dict.fromkeys(r["contig_uid"] for r in rows))
        circ = [u for u in uids if any(r["contig_uid"] == u and r["topology"] == "circular" for r in rows)]
        r0 = rows[0]
        objects.append({"path": str((fasta_dir / r0["file"]).resolve()), "source_name": "synthworld",
                        "source_release": "micro", "license": "synthetic data, CC0", "redistributable": True})
        guid = ids.genome_uid([bytes.fromhex(u) for u in uids])
        if guid in seen:
            rep["genomes_listed_once"] += 1
            continue
        seen.add(guid)
        genomes.append({"contig_uids": uids, "genetic_code": int(r0["genetic_code"]), "kind": r0["kind"],
                        "circular": circ, "assembly_quality": r0["assembly_quality"]})
    rep["genomes"], rep["objects"] = len(genomes), len(objects)
    return genomes, objects, dict(rep)


def inputs(world_dir: str, out: str) -> int:
    ed = Path(world_dir) / "emit"
    out_p = Path(out)
    out_p.mkdir(parents=True, exist_ok=True)
    by_src = collections.defaultdict(list)
    for r in _tsv(ed / "ingest.tsv"):
        by_src[(r["source_label"], r["genome_id"])].append(r)
    stages = {"F0": [(k, v) for k, v in sorted(by_src.items()) if k[0] == "src_a"],
              "F1": [(k, v) for k, v in sorted(by_src.items()) if k[0] != "src_a"]}
    report = {}
    for s, srcs in stages.items():
        g, o, rep = _listing(srcs, ed)
        (out_p / f"genomes-{s}.json").write_text(json.dumps(g))
        (out_p / f"objects-{s}.json").write_text(json.dumps(o, indent=1))
        report[s] = rep
    print(json.dumps(report))
    return 0


# --- compare -------------------------------------------------------------------------------------------------------

def _record(root: Path, stage: str, step: str) -> dict:
    return json.loads((root / "runs" / f"stage={stage}" / f"{step}.json").read_text())


def _files(rec: dict, prefix: str) -> list[str]:
    return sorted(o["path"] for o in rec["outputs"] if o["path"].startswith(prefix) and o["path"].endswith(".parquet"))


def _read(root: Path, paths: list[str], cols: list[str]) -> list[dict]:
    out = []
    for p in paths:
        out += pq.read_table(root / p, columns=cols).to_pylist()
    return out


def compare(world_dir: str, root_s: str, release: str, out: str) -> int:
    wd, root, out_p = Path(world_dir), Path(root_s), Path(out)
    out_p.mkdir(parents=True, exist_ok=True)
    ed = wd / "emit"
    ingest = _tsv(ed / "ingest.tsv")
    truth = {n: _tsv(ed / "truth" / f"{n}.tsv") for n in ("truth_gene", "truth_callerr", "truth_contig", "truth_dup",
                                                           "truth_genome")}
    fa_seq = {}
    for rel in sorted({r["file"] for r in ingest}):
        name = None
        chunks: list[str] = []
        for line in (ed / rel).read_text().splitlines():
            if line.startswith(">"):
                if name:
                    fa_seq[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
        if name:
            fa_seq[name] = "".join(chunks)
    uid_of_acc = {r["contig_accession"]: r["contig_uid"] for r in ingest}
    length = {r["contig_uid"]: int(r["length"]) for r in ingest}
    origin = {}
    for r in truth["truth_contig"]:
        if r["is_dropped"] == "false":
            origin.setdefault(r["contig_uid"], r["origin"])
    report: dict = {}

    # --- contigs: stored once, with their length and digest
    stages = ("F0", "F1")
    index, occ = [], []
    for s in stages:
        rec = _record(root, s, "ingest")
        index += [dict(r, stage=s) for r in _read(root, _files(rec, f"l0seq/contig_index/stage={s}/"),
                                                  ["contig_uid", "raw_digest", "contig_len", "is_default_variant"])]
        occ += [dict(r, stage=s) for r in _read(root, _files(rec, f"l0seq/contig_occurrence/"),
                                                ["accession", "contig_uid"])]
    default = collections.defaultdict(list)
    for r in index:
        if r["is_default_variant"]:
            default[r["contig_uid"].hex()].append(r)
    digest_ok = len_ok = 0
    for acc, s in fa_seq.items():
        u = uid_of_acc[acc]
        d = default.get(u, [])
        if len(d) == 1 and d[0]["raw_digest"].hex() == hashlib.sha256(s.encode()).hexdigest():
            digest_ok += 1
        if len(d) == 1 and d[0]["contig_len"] == len(s):
            len_ok += 1
    occ_uid = {r["accession"]: r["contig_uid"].hex() for r in occ}
    report["contigs"] = {
        "emitted_occurrences": len(ingest), "distinct_uids": len(length),
        "stored_default_variants": len(default), "uids_with_one_default": sum(1 for v in default.values() if len(v) == 1),
        "missing": len(set(length) - set(default)), "extra": len(set(default) - set(length)),
        "length_ok": len_ok, "digest_ok": digest_ok,
        "occurrence_rows": len(occ), "occurrences_matching_uid": sum(1 for a, u in occ_uid.items() if uid_of_acc.get(a) == u),
        "default_stage": dict(collections.Counter(v[0]["stage"] for v in default.values()))}

    # --- duplicates, DEC-26 and DEC-27
    dup = truth["truth_dup"]
    exact = [d for d in dup if d["kind"] == "exact"]
    near = [d for d in dup if d["kind"] == "near"]
    from biobox import ledger

    skipped = []
    for s in stages:
        for r in ledger.read(root, s, "genecall"):
            if str(r["reason_code"]).startswith("owned_by:"):
                key = r["object_key"]
                skipped.append({"stage": s, "contig_uid": key.hex() if isinstance(key, bytes) else str(key),
                                "owner_stage": str(r["reason_code"])[len("owned_by:"):]})
    exact_uids = {d["member_contig_uid"] for d in exact}
    dim_rec = _record(root, release, "dim")
    cdim = _read(root, _files(dim_rec, f"dim/contig_dim/release={release}/"), ["contig_uid", "genome_uid", "n_genomes"])
    gdim = _read(root, _files(dim_rec, f"dim/genome_dim/release={release}/"),
                 ["genome_uid", "n_contigs", "kind", "genetic_code", "assembly_quality"])
    n_gen = {r["contig_uid"].hex(): r["n_genomes"] for r in cdim}
    tg_uids = {r["genome_uid"] for r in truth["truth_genome"]}
    gd = {r["genome_uid"].hex(): r for r in gdim}
    attr_ok = sum(1 for r in truth["truth_genome"] if r["genome_uid"] in gd
                  and gd[r["genome_uid"]]["genetic_code"] == int(r["genetic_code"])
                  and gd[r["genome_uid"]]["kind"] == r["kind"]
                  and gd[r["genome_uid"]]["assembly_quality"] == r["assembly_quality"])
    report["duplicates"] = {
        "exact_member_contigs": len(exact), "exact_same_uid": sum(1 for d in exact if d["rep_contig_uid"] == d["member_contig_uid"]),
        "exact_uids_one_default": sum(1 for u in exact_uids if len(default.get(u, [])) == 1),
        "exact_uids_default_in_F0": sum(1 for u in exact_uids if default.get(u) and default[u][0]["stage"] == "F0"),
        "exact_uids_n_genomes": dict(collections.Counter(str(n_gen.get(u)) for u in exact_uids)),
        "near_member_contigs": len(near), "near_distinct_uids": sum(1 for d in near if d["rep_contig_uid"] != d["member_contig_uid"]),
        "near_both_stored": sum(1 for d in near if d["rep_contig_uid"] in default and d["member_contig_uid"] in default),
        "genecall_skipped": dict(collections.Counter(f"{x['stage']} owned by {x['owner_stage']}" for x in skipped)),
        # F1's contigs that F0 stored first: every exact duplicate's, and a near duplicate's left unchanged
        "skipped_are_f1_contigs_stored_by_f0": sorted({x["contig_uid"] for x in skipped if x["stage"] == "F1"}) == sorted(
            {d["member_contig_uid"] for d in dup if d["member_contig_uid"] == d["rep_contig_uid"]}),
        "contigs_n_genomes": dict(collections.Counter(str(v) for v in n_gen.values())),
        "genome_dim_rows": len(gdim), "truth_distinct_genomes": len(tg_uids),
        "genome_dim_matches_truth": len(set(gd) & tg_uids), "genome_attrs_agree": attr_ok}

    # --- gene calls
    gcv = None
    calls = collections.defaultdict(list)            # contig_uid -> (start, end, strand)
    for s in stages:
        rec = _record(root, s, "genecall")
        paths = [p for p in _files(rec, "l1/feature/")]
        gcv = paths[0].split("producer=")[1].split("/")[0] if paths else gcv
        for r in _read(root, paths, ["contig_uid", "start", "end", "strand", "kind"]):
            if r["kind"] == 1:
                calls[r["contig_uid"].hex()].append((r["start"], r["end"], r["strand"]))
    stop_index = {}
    for cu, cs in calls.items():
        for c in cs:
            stop_index[(cu, c[2], c[1] if c[2] == 1 else c[0])] = c
    callerr_kind = {}
    for r in truth["truth_callerr"]:
        callerr_kind[r["callerr_id"]] = r["kind"]
    code = {r["genome_id"]: int(r["genetic_code"]) for r in truth["truth_genome"] if r["source_label"] == "src_a"}
    alt_code_genes = {r["true_gene_id"] for r in truth["truth_callerr"] if r["kind"] == "alt_code"}
    matched_calls = set()
    cls = collections.defaultdict(lambda: [0, 0, 0])   # class -> [n, found, start exact]

    def frame_match(cu, s, e, st, p5, p3):
        """A call on the strand, in the piece's frame (anchored at its 5' end) and covering half of it."""
        if p5 and p3:
            return None
        anchor = s if st == 1 else e
        for c in calls.get(cu, ()):
            if c[2] != st:
                continue
            ca = c[0] if st == 1 else c[1]
            if (ca - anchor) % 3:
                continue
            ov = min(e, c[1]) - max(s, c[0])
            if ov * 2 >= e - s:
                return c
        return None

    for r in truth["truth_gene"]:
        cu, s, e, st = r["contig_uid"], int(r["start"]), int(r["end"]), int(r["strand"])
        p5, p3 = r["partial_5"] == "true", r["partial_3"] == "true"
        if not p3:
            c = stop_index.get((cu, st, int(r["stop_bp"])))
        else:
            c = frame_match(cu, s, e, st, p5, p3)
        if c is not None:
            matched_calls.add((cu, c))
        names = []
        occ_ = r["occurrence"]
        if occ_ in ("phage", "near_dup", "contaminant"):
            names.append(f"{occ_} genes")
        elif p5 or p3:
            names.append("assembly: cut by a contig edge, 3' end kept" if not p3 else
                         "assembly: cut by a contig edge, 3' end cut")
        else:
            g = r["genome_id"]
            kind = callerr_kind.get(r["callerr_id"], "") if r["callerr_id"] else ""
            if kind:
                names.append(f"call error: {kind}")
            elif code.get(g) == 4:
                names.append("code 4 genome: gene with internal TGA" if r["gene_id"] in alt_code_genes else
                             "code 4 genome: gene without internal TGA")
            elif code.get(g) == 25:
                names.append("code 25 genome: gene with internal TGA" if r["gene_id"] in alt_code_genes else
                             "code 25 genome: gene without internal TGA")
            elif r["source"] == "host" and r["feature_kind"] == "CDS":
                names.append("untouched (host CDS, code 11, whole)")
            else:
                names.append(f"source {r['source']} (code 11, whole)")
            if r["origin_spanning"] == "true":
                names.append("across a circular origin")
        names.append("all truth genes")
        for nm in names:
            k = cls[nm]
            k[0] += 1
            if c is not None:
                k[1] += 1
                if not (p5 or p3) and c == (s, e, st):
                    k[2] += 1
    rows = []
    for nm in sorted(cls):
        n, f, x = cls[nm]
        rows.append({"class": nm, "n": n, "found": f, "recall": round(f / n, 4) if n else "",
                     "start_exact": x, "start_exact_of_found": round(x / f, 4) if f else ""})
    _write_tsv(out_p / "recall_by_class.tsv", rows)
    # the code-11 pieces of broken genes (truth_callerr expected_calls)
    exp = collections.defaultdict(lambda: [0, 0])
    for r in truth["truth_callerr"]:
        if r["status"] != "whole":
            continue
        for item in filter(None, r["expected_calls"].split(";")):
            span, cst = item.split(":")
            x, y = (int(v) for v in span.split("-"))
            st = int(cst)
            exp[r["kind"]][0] += 1
            exp[r["kind"]][1] += (r["contig_uid"], st, y if st == 1 else x) in stop_index
    report["expected_calls_found"] = {k: f"{v[1]}/{v[0]}" for k, v in sorted(exp.items())}
    # precision by contig origin
    prec = collections.defaultdict(lambda: [0, 0])
    unmatched_len = collections.defaultdict(list)
    for cu, cs in calls.items():
        o = origin.get(cu, "not emitted")
        for c in cs:
            prec[o][0] += 1
            if (cu, c) in matched_calls:
                prec[o][1] += 1
            else:
                unmatched_len[o].append((c[1] - c[0]) // 3)
    prows = []
    for o in sorted(prec):
        n, m = prec[o]
        ul = sorted(unmatched_len[o])
        prows.append({"contig_origin": o, "calls": n, "matching_a_truth_gene": m, "precision": round(m / n, 4),
                      "unmatched": n - m, "unmatched_median_codons": ul[len(ul) // 2] if ul else ""})
    _write_tsv(out_p / "precision_by_origin.tsv", prows)
    report["gene_call_version"] = gcv
    report["calls"] = sum(len(v) for v in calls.values())

    # --- the steps' cost
    srows = []
    for st_dir in sorted((root / "runs").iterdir()):
        for f in sorted(st_dir.glob("*.json")):
            if f.name.endswith(".publishing.json") or ".failed" in f.name:
                continue
            r = json.loads(f.read_text())
            if "wall_s" in r:
                srows.append({"stage": r["stage"], "step": r["step"], "wall_s": r["wall_s"], "cpu_min": r["cpu_min"],
                              "peak_rss_mb": round(r["peak_rss_bytes"] / 2**20, 1) if r["peak_rss_bytes"] else "",
                              "threads": r["threads"]})
    if srows:
        _write_tsv(out_p / "steps.tsv", srows)
    print(json.dumps(report, indent=1, default=str))
    return 0


def main(argv: list[str]) -> int:
    cmds = {"world": (world, 1), "inputs": (inputs, 2), "compare": (compare, 4)}
    if not argv or argv[0] not in cmds or len(argv) - 1 != cmds[argv[0]][1]:
        print(USAGE, file=sys.stderr)
        return 2
    return cmds[argv[0]][0](*argv[1:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

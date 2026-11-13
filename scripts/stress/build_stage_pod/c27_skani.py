#!/usr/bin/env python3
"""Card C-27 (4): skani ANI for S2C phase 2's diversity pick, per species (stdlib only; skani from the image).

    python3 -I c27_skani.py normalise SKANI_TSV... -o OUT [--expect-query ACC --expect-refs FILE]
    python3 -I c27_skani.py paths OUT_TSV LOCK_JSON...
    python3 -I c27_skani.py plan CANDIDATES_TSV PATHS_TSV WORK [--cap 20] [--min-af 15] [--threads 8]
                            [--budget-cpu-min 240] [--sketch-cpu-s 2] [--pair-ms 10] [--image SIF]
                            [--unavailable FILE] [--dry-run]
    python3 -I c27_skani.py launch WORK [--parallel 2] [--dry-run]
    python3 -I c27_skani.py run-shard WORK K [--skani skani] [--threads N]
    python3 -I c27_skani.py collect WORK OUT_TSV

Which pairs: the farthest-first rounds, not all pairs. In a species of n candidates the pick of CAP genomes from
rank 1 needs the distances from each of the first CAP-1 picks to the genomes not yet picked: sum over species of
(CAP-1) n - CAP (CAP-1) / 2, 499,624 pairs for the 243 species at CAP 20 (derived from the per-species counts in
planning/measurements/2026-11-12_s2c_ani/ani_needs.json), against 15,835,149 for all pairs (E. coli alone
8,026,021). The rounds are sequential: round r runs `skani dist` with the r-th pick as query against the sketches
of the genomes not yet picked, and the next pick is chosen from all rounds so far by c25_select.farthest_first
(the same code c25_select.py --ani runs on the finished table, with the same reading of a row,
c25_select.ani_value). One `skani triangle -E` of the CAP best-quality genomes adds the pairs among the quality
pick that the rounds did not cover (CAP (CAP-1) / 2 per species at most), so the report can compare the minimum
pairwise ANI of both picks. A pair skani was asked for and did not print (below its screen or its --min-af) is
written as a row with an empty ani: it counts as the largest distance.

skani 0.4.0 as checked in the image (`skani <cmd> --help` and a run on four real genomes, 2026-11-13):
`sketch --separate-sketches -l LIST -o DIR` writes DIR/<FASTA file name>.sketch per genome (without the flag 0.4.0
writes one consolidated database, sketches.db2, that dist cannot take as -q); DIR must not exist; .fna.gz is read
directly. `dist -q SKETCH --rl LIST -o OUT` and `triangle -E -l LIST -o OUT` print the columns Ref_file,
Query_file, ANI, Align_fraction_ref, Align_fraction_query, Ref_name, Query_name, with the FASTA path a sketch was
made from in the file columns and each fraction on its own genome; a pair is printed only when one genome's
aligned fraction is > --min-af (default 15) and it passes the screen (-s, about 80% identity): any other pair is
absent. The driver passes --min-af 0, so every pair past the screen is printed and c25_select.py --min-af applies
the 15% rule itself. -t is the number of threads.

normalise: skani's dist/triangle -E output -> the table c25_select.py --ani reads: genome_a < genome_b (bare
accessions, from the file names), ani, af_a, af_b as skani printed them; a pair printed twice keeps its larger
ANI; with --expect-query/--expect-refs every pair asked for and not printed becomes a row with an empty ani.

paths: accession -> FASTA path from fetch lock.json files (phase 1's and this card's), for plan.

plan: species with more than CAP candidates, in quality order (c25_select.py --candidates' species_rank), less
the genomes listed with --unavailable (not served; refused for a rank-1 genome; any other candidate without a
file is refused), packed
into shards by first-fit decreasing on an estimate of CPU-minutes, (genomes x sketch-cpu-s + pairs x pair-ms) /
60, so that each shard process stays under --budget-cpu-min (default 240, under the 300 CPU-minute rule; a
species over the budget is a shard alone, and plan refuses one over 300). The cost per genome and per pair are
estimates until step 4 measures one species. Writes WORK/plan.json and WORK/members.tsv; --dry-run writes
nothing and prints the plan and the unit commands.

launch: starts up to --parallel shards that are neither DONE nor active, each as unit c27-skani-<k>
(systemd-run --user, nice 19, MemoryMax 32G, CPUQuota threads x 100%), running run-shard inside the image.

run-shard: one species at a time (bounded memory: one species' sketches on disk, its pairs in memory), each in
WORK/species/<id>/ with pairs.tsv, species.json (times) and DONE; a species with DONE is skipped, so a shard
resumes. At the end WORK/shards/<k>/DONE, or FAILED with the reason.

collect: every species DONE -> one table, sorted, for c25_select.py --ani TABLE --require-pairs.
"""

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I leaves the script's directory off the path
from c25_select import ANI_COLS, ani_value, farthest_first  # noqa: E402

ACC = re.compile(r"^(GC[AF]_\d{9}\.\d+)")
REQUIRED = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"]
STOP = Path("/stor/scratch/Wilke/if2453/all-the-data/STOP")
PROCESS_CPU_MIN = 300  # AGENTS.md section 7.2
# skani prints a pair only when one genome's aligned fraction is > --min-af (default 15): print every pair past the
# screen, so that c25_select.py --min-af is the one place that rule applies and the table shows what it removed
PRINT_ALL = ["--min-af", "0"]


def put_text(path: Path, text: str) -> str:
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return hashlib.sha256(text.encode()).hexdigest()


def accession_of(path: str) -> str:
    m = ACC.match(Path(path).name)
    if not m:
        raise SystemExit(f"no assembly accession at the start of the file name {path!r}")
    return m.group(1)


def _rank(row: tuple[str, str, str]) -> tuple:
    ani, af_a, af_b = row
    if ani == "":
        return (-1.0, -1.0, row)
    return (float(ani), max(float(af_a), float(af_b)), row)


def merge(pairs: dict, key: tuple[str, str], row: tuple[str, str, str]) -> None:
    """Keep the row with the larger ANI (then larger alignment fraction, then the text), whatever the order."""
    if key not in pairs or _rank(row) > _rank(pairs[key]):
        pairs[key] = row


def parse_skani(text: str, source) -> dict[tuple[str, str], tuple[str, str, str]]:
    """skani dist / triangle -E output -> {(a, b) with a < b: (ani, af_a, af_b)} as printed."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return {}
    head = lines[0].split("\t")
    miss = [c for c in REQUIRED if c not in head]
    if miss:
        raise SystemExit(f"{source}: not skani output, missing columns {miss}")
    ix = {c: head.index(c) for c in REQUIRED}
    out: dict = {}
    for ln in lines[1:]:
        f = ln.split("\t")
        r, q = accession_of(f[ix["Ref_file"]]), accession_of(f[ix["Query_file"]])
        if r == q:
            continue
        ani, afr, afq = f[ix["ANI"]], f[ix["Align_fraction_ref"]], f[ix["Align_fraction_query"]]
        merge(out, (r, q) if r < q else (q, r), (ani, afr, afq) if r < q else (ani, afq, afr))
    return out


def fill_expected(pairs: dict, query: str, refs: list[str]) -> None:
    for g in refs:
        if g != query:
            key = (query, g) if query < g else (g, query)
            if key not in pairs:
                pairs[key] = ("", "", "")


def table_text(pairs: dict) -> str:
    return "\t".join(ANI_COLS) + "\n" + "".join(f"{a}\t{b}\t{v[0]}\t{v[1]}\t{v[2]}\n" for (a, b), v in sorted(pairs.items()))


def read_table(path: Path) -> dict:
    with open(path, newline="") as fh:
        return {(r["genome_a"], r["genome_b"]): (r["ani"], r["af_a"], r["af_b"])
                for r in csv.DictReader(fh, delimiter="\t")}


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- plan ------------------------------------------------------------------------------------------------------

def slug(species: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", species.split(";")[-1].removeprefix("s__")).strip("_")


def build_plan(a) -> tuple[dict, list[tuple]]:
    paths = {r["accession"]: r["path"] for r in read_tsv(a.paths)}
    unavailable = set(a.unavailable.read_text().split()) if a.unavailable else set()
    by_sp: dict[str, list[dict]] = {}
    for r in read_tsv(a.candidates):
        by_sp.setdefault(r["gtdb_taxonomy"], []).append(r)
    cap = a.cap
    species, members, missing, left_out = [], [], [], []
    for i, sp in enumerate(sorted(s for s, rows in by_sp.items() if len(rows) > cap)):
        rows = sorted(by_sp[sp], key=lambda r: int(r["species_rank"]))
        assert [int(r["species_rank"]) for r in rows] == list(range(1, len(rows) + 1)), f"ranks of {sp}"
        if rows[0]["accession"] in unavailable:
            raise SystemExit(f"{sp}: its rank 1 genome {rows[0]['accession']} is listed as unavailable")
        left_out += [r["accession"] for r in rows if r["accession"] in unavailable]
        rows = [r for r in rows if r["accession"] not in unavailable]
        n = len(rows)
        ff = sum(n - r for r in range(1, min(cap, n)))  # round r: the r-th pick against the n - r not yet picked
        qtop = min(cap, n) * (min(cap, n) - 1) // 2
        sid = f"sp{i:03d}_{slug(sp)}"
        for r in rows:
            p = paths.get(r["accession"], "")
            if not p or not os.path.exists(p):
                missing.append(r["accession"])
            members.append((sid, int(r["species_rank"]), r["accession"], p))
        species.append({"id": sid, "species": sp, "genomes": n, "bp": sum(int(r["genome_size"]) for r in rows),
                        "pairs_farthest_first": ff, "pairs_quality_top": qtop, "pairs_all": n * (n - 1) // 2,
                        "est_cpu_min": (n * a.sketch_cpu_s + (ff + qtop) * a.pair_ms / 1000) / 60})
    shards: list[dict] = []
    for s in sorted(species, key=lambda s: (-s["est_cpu_min"], s["id"])):
        if s["est_cpu_min"] >= PROCESS_CPU_MIN:
            raise SystemExit(f"{s['species']}: estimated {s['est_cpu_min']:.0f} CPU-minutes, over the "
                             f"{PROCESS_CPU_MIN}-minute process rule; split it first")
        for sh in shards:
            if sh["est_cpu_min"] + s["est_cpu_min"] <= a.budget_cpu_min:
                break
        else:
            sh = {"k": len(shards), "unit": f"c27-skani-{len(shards)}", "est_cpu_min": 0.0, "species": []}
            shards.append(sh)
        sh["species"].append(s)
        sh["est_cpu_min"] += s["est_cpu_min"]
    for sh in shards:
        sh["species"].sort(key=lambda s: s["id"])
        sh["est_cpu_min"] = sum(s["est_cpu_min"] for s in sh["species"])
    plan = {
        "cap": cap, "min_af": a.min_af, "threads": a.threads, "image": a.image, "work": str(a.work.resolve()),
        "script": str(Path(__file__).resolve()), "budget_cpu_min": a.budget_cpu_min,
        "sketch_cpu_s": a.sketch_cpu_s, "pair_ms": a.pair_ms,
        "species": len(species), "genomes": sum(s["genomes"] for s in species),
        "bp": sum(s["bp"] for s in species),
        "pairs_farthest_first": sum(s["pairs_farthest_first"] for s in species),
        "pairs_quality_top": sum(s["pairs_quality_top"] for s in species),
        "pairs_all": sum(s["pairs_all"] for s in species),
        "est_cpu_min_total": sum(s["est_cpu_min"] for s in species),
        "unavailable": sorted(left_out),
        "missing_paths": len(missing), "missing_paths_first": sorted(missing)[:10],
        "candidates_sha256": hashlib.sha256(a.candidates.read_bytes()).hexdigest(),
        "shards": shards,
    }
    return plan, members


def unit_command(plan: dict, sh: dict) -> str:
    w, k, t = plan["work"], sh["k"], plan["threads"]
    return (f"systemd-run --user --unit {sh['unit']} --collect --quiet -p MemoryMax=32G -p MemorySwapMax=0 "
            f"-p CPUQuota={t * 100}% -p StandardOutput=append:{w}/logs/{sh['unit']}.log "
            f"-p StandardError=append:{w}/logs/{sh['unit']}.log "
            f"nice -n 19 singularity exec --cleanenv --bind /tmp/if2453 {plan['image']} "
            f"python -I {plan['script']} run-shard {w} {k} --threads {t}")


def cmd_plan(a) -> int:
    plan, members = build_plan(a)
    text = json.dumps(plan, indent=1)
    cmds = "\n".join(unit_command(plan, sh) for sh in plan["shards"])
    if a.dry_run:
        print(text + "\n# unit commands (launch starts at most --parallel at a time):\n" + cmds)
        return 0
    if plan["missing_paths"]:
        raise SystemExit(f"{plan['missing_paths']} candidate genomes have no file, e.g. {plan['missing_paths_first']} "
                         "(genomes that were not served are listed with --unavailable)")
    if not plan["image"]:
        raise SystemExit("--image is required for a plan that will run")
    a.work.mkdir(parents=True, exist_ok=True)
    put_text(a.work / "members.tsv", "species_id\trank\taccession\tpath\n" +
             "".join(f"{s}\t{r}\t{acc}\t{p}\n" for s, r, acc, p in members))
    put_text(a.work / "plan.json", text + "\n")
    print(text + "\n# unit commands:\n" + cmds)
    return 0


def is_active(unit: str) -> bool:
    return subprocess.run(["systemctl", "--user", "is-active", "--quiet", unit]).returncode == 0


def cmd_launch(a) -> int:
    plan = json.loads((a.work / "plan.json").read_text())
    pending = [sh for sh in plan["shards"] if not (a.work / "shards" / str(sh["k"]) / "DONE").exists()]
    if not a.dry_run:
        if STOP.exists():
            raise SystemExit(f"STOP present at {STOP}")
        active = [sh for sh in pending if is_active(sh["unit"])]
        pending = [sh for sh in pending if sh not in active]
        room = max(0, a.parallel - len(active))
    else:
        room = a.parallel
    for sh in pending[:room]:
        cmd = unit_command(plan, sh)
        print(cmd)
        if not a.dry_run:
            (a.work / "logs").mkdir(parents=True, exist_ok=True)
            subprocess.run(["systemctl", "--user", "reset-failed", sh["unit"]], capture_output=True)
            subprocess.run(cmd.split(), check=True)
    print(f"# {len(pending)} shards pending, {min(room, len(pending))} started" + (" (dry run)" if a.dry_run else ""))
    return 0


# --- run -------------------------------------------------------------------------------------------------------

def run(cmd: list, log: list) -> float:
    t = time.time()
    subprocess.run([str(c) for c in cmd], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    dt = time.time() - t
    log.append({"cmd": [str(c) for c in cmd][:2], "wall_s": round(dt, 3)})
    return dt


def run_species(work: Path, sp: dict, members: list[tuple[str, str]], plan: dict, skani: str, threads: int) -> dict:
    d = work / "species" / sp["id"]
    if (d / "DONE").exists():
        return json.loads((d / "DONE").read_text())
    d.mkdir(parents=True, exist_ok=True)
    cap, min_af = plan["cap"], plan["min_af"]
    ranked = [acc for acc, _ in members]
    path = dict(members)
    t0 = time.time()
    log: list[dict] = []
    put_text(d / "genomes.txt", "".join(path[g] + "\n" for g in ranked))
    sk = d / "sketches"
    shutil.rmtree(sk, ignore_errors=True)  # skani refuses an existing output folder
    run([skani, "sketch", "--separate-sketches", "-t", threads, "-l", d / "genomes.txt", "-o", sk], log)
    sketch = {}
    for p in sorted(sk.iterdir()):
        m = ACC.match(p.name)
        if m and p.name.endswith(".sketch"):
            sketch[m.group(1)] = p
    if set(sketch) != set(ranked):
        raise RuntimeError(f"{sp['id']}: sketches missing for {sorted(set(ranked) - set(sketch))[:5]}")
    pairs: dict = {}

    def round_fn(newest: str, rest: list[str]) -> dict:
        r = sum(1 for e in log if e["cmd"][1] == "dist") + 1
        refs, out = d / f"round_{r:02d}.refs", d / f"round_{r:02d}.tsv"
        put_text(refs, "".join(f"{sketch[g]}\n" for g in rest))
        run([skani, "dist", "-t", threads, *PRINT_ALL, "-q", sketch[newest], "--rl", refs, "-o", out], log)
        got = parse_skani(out.read_text(), out)
        fill_expected(got, newest, rest)
        sims = {}
        for g in rest:
            key = (newest, g) if newest < g else (g, newest)
            merge(pairs, key, got[key])
            sims[g] = ani_value(*pairs[key], min_af, str(out))
        return sims

    picks = farthest_first(ranked, cap, round_fn)
    qtop = ranked[:cap]
    put_text(d / "quality_top.txt", "".join(f"{sketch[g]}\n" for g in qtop))
    run([skani, "triangle", "-E", "-t", threads, *PRINT_ALL, "-l", d / "quality_top.txt", "-o", d / "quality_top.tsv"],
        log)
    tri = parse_skani((d / "quality_top.tsv").read_text(), d / "quality_top.tsv")
    for i, g in enumerate(qtop):
        fill_expected(tri, g, qtop[i + 1:])
    added = 0
    for key, row in tri.items():
        if key[0] in qtop and key[1] in qtop and key not in pairs:  # the rounds' values stay as the pick read them
            pairs[key] = row
            added += 1
    sha = put_text(d / "pairs.tsv", table_text(pairs))
    res = {"id": sp["id"], "species": sp["species"], "genomes": len(ranked), "picks": [g for g, _ in picks],
           "pairs": len(pairs), "pairs_quality_top_added": added,
           "pairs_unreported": sum(v[0] == "" for v in pairs.values()), "pairs_sha256": sha,
           "wall_s": round(time.time() - t0, 3), "calls": log}
    put_text(d / "species.json", json.dumps(res, indent=1) + "\n")
    put_text(d / "DONE", json.dumps({k: res[k] for k in ("id", "genomes", "pairs", "pairs_sha256", "wall_s")}) + "\n")
    return res


def members_of(work: Path) -> dict[str, list[tuple[str, str]]]:
    out: dict[str, list[tuple[int, str, str]]] = {}
    for r in read_tsv(work / "members.tsv"):
        out.setdefault(r["species_id"], []).append((int(r["rank"]), r["accession"], r["path"]))
    return {s: [(acc, p) for _, acc, p in sorted(v)] for s, v in out.items()}


def cmd_run_shard(a) -> int:
    plan = json.loads((a.work / "plan.json").read_text())
    sh = plan["shards"][a.k]
    mark = a.work / "shards" / str(a.k)
    mark.mkdir(parents=True, exist_ok=True)
    for m in ("DONE", "FAILED"):
        (mark / m).unlink(missing_ok=True)
    members = members_of(a.work)
    threads = a.threads or plan["threads"]
    t0, done = time.time(), []
    try:
        for sp in sh["species"]:
            if STOP.exists():
                raise RuntimeError(f"STOP present at {STOP}")
            r = run_species(a.work, sp, members[sp["id"]], plan, a.skani, threads)
            done.append({"id": sp["id"], "wall_s": r["wall_s"], "pairs": r["pairs"]})
    except Exception as e:  # noqa: BLE001
        msg = f"{type(e).__name__}: {e}"
        if isinstance(e, subprocess.CalledProcessError) and e.stderr:
            msg += " | " + e.stderr.decode(errors="replace")[-500:]
        put_text(mark / "FAILED", json.dumps({"reason": msg, "done": done}) + "\n")
        print(msg, file=sys.stderr)
        return 1
    put_text(mark / "DONE", json.dumps({"species": len(done), "wall_s": round(time.time() - t0, 3),
                                        "per_species": done}) + "\n")
    return 0


def cmd_collect(a) -> int:
    plan = json.loads((a.work / "plan.json").read_text())
    pairs: dict = {}
    not_done = []
    for sh in plan["shards"]:
        for sp in sh["species"]:
            d = a.work / "species" / sp["id"]
            if not (d / "DONE").exists():
                not_done.append(sp["id"])
                continue
            for key, row in read_table(d / "pairs.tsv").items():
                merge(pairs, key, row)
    if not_done:
        raise SystemExit(f"{len(not_done)} species not DONE, e.g. {not_done[:5]}")
    sha = put_text(a.out, table_text(pairs))
    print(json.dumps({"species": plan["species"], "pairs": len(pairs),
                      "pairs_unreported": sum(v[0] == "" for v in pairs.values()), "sha256": sha}))
    return 0


def cmd_normalise(a) -> int:
    pairs: dict = {}
    for p in a.inputs:
        for key, row in parse_skani(p.read_text(), p).items():
            merge(pairs, key, row)
    if a.expect_query:
        fill_expected(pairs, a.expect_query, [x for x in a.expect_refs.read_text().split() if x])
    put_text(a.out, table_text(pairs))
    return 0


def cmd_paths(a) -> int:
    out: dict[str, str] = {}
    for lp in a.locks:
        for acc, v in json.loads(lp.read_text())["genomes"].items():
            p = v.get("file") or str(lp.parent / v["file_name"])
            if acc in out and out[acc] != p:
                raise SystemExit(f"{acc} has two files: {out[acc]} and {p}")
            out[acc] = p
    put_text(a.out, "accession\tpath\n" + "".join(f"{k}\t{v}\n" for k, v in sorted(out.items())))
    return 0


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("normalise")
    p.add_argument("inputs", type=Path, nargs="+")
    p.add_argument("-o", "--out", type=Path, required=True)
    p.add_argument("--expect-query")
    p.add_argument("--expect-refs", type=Path)
    p = sub.add_parser("paths")
    p.add_argument("out", type=Path)
    p.add_argument("locks", type=Path, nargs="+")
    p = sub.add_parser("plan")
    p.add_argument("candidates", type=Path)
    p.add_argument("paths", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--cap", type=int, default=20)
    p.add_argument("--min-af", type=float, default=15.0)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--budget-cpu-min", type=float, default=240.0)
    p.add_argument("--sketch-cpu-s", type=float, default=2.0, help="estimate: CPU-seconds to sketch one genome")
    p.add_argument("--pair-ms", type=float, default=10.0, help="estimate: CPU-milliseconds per sketched pair")
    p.add_argument("--image", default="")
    p.add_argument("--unavailable", type=Path, help="file of accessions with no sequence (not served): left out")
    p.add_argument("--dry-run", action="store_true")
    p = sub.add_parser("launch")
    p.add_argument("work", type=Path)
    p.add_argument("--parallel", type=int, default=2)
    p.add_argument("--dry-run", action="store_true")
    p = sub.add_parser("run-shard")
    p.add_argument("work", type=Path)
    p.add_argument("k", type=int)
    p.add_argument("--skani", default="skani")
    p.add_argument("--threads", type=int, default=0)
    p = sub.add_parser("collect")
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    a = ap.parse_args(argv)
    return {"normalise": cmd_normalise, "paths": cmd_paths, "plan": cmd_plan, "launch": cmd_launch,
            "run-shard": cmd_run_shard, "collect": cmd_collect}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

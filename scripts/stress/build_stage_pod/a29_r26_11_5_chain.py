#!/usr/bin/env python3
"""Card A-29, part 1: the `r26.11.5` chain, exact prefix ownership over every built and staged stage, the confs promoted.

    python a29_r26_11_5_chain.py ownership --root ROOT --tree TREE --out OUT.json [--inputs-home DIR] [--tmp DIR]
                                           [--workers N] [--meta STAGE=DIR ...]
    python a29_r26_11_5_chain.py write-confs --tree TREE [--ownership OUT.json] [--date YYYY-MM-DD]

`ownership` computes DEC-26 ownership by build order (card A-26's method, in bounded memory: contig uids in Parquet
parts and DuckDB with a memory limit): every contig belongs to the first stage of the order that lists it. The order
is the root's 13 stored stages (`ROOT_STAGES`, read from the root's `contig_index`, the other root partitions left
alone: the `r26.11.4` stages the wave has built so far are taken from their inputs, in their place in the order),
then the 19 stages of `r26.11.4` (their `genomes.json` under `--inputs-home`), then `R26_11_5_CHAIN`, each read from
its staging directory: `genomes.json` there, and the FASTA objects that its `objects_staging.json` names (on `/stor`
or under `/tmp/if2453/<card>/`), each record's uid and length by `biobox.fasta.iter_records` and
`biobox.ids.contig_uid`, as ingest computes them. For each stage of the chain it also checks that the sha256 of
`genomes.json`, `objects.json` and every object equals its template's `SUM` line, that `objects.json` names the
objects under `/tmp/if2453/s2/inputs/<STAGE>/`, and it compares the template's `EXPECT=held_by` lines (written by the
staging cards over the root and `r26.11.4`, some over earlier `r26.11.5` stages too) with the counts found: a
template line that differs, or an owner of the root or `r26.11.4` that holds contigs and has no line, is a
discrepancy. The report is written under a temporary name and renamed; the exit status is 0 only when it has no
problem and no discrepancy.

`write-confs` writes `stages/<STAGE>.conf` for every stage of `R26_11_5_CHAIN` from `stages/templates/<STAGE>.conf`:
a first comment line saying so, the `AFTER=` of the chain, this card's note as a `NOTE`, the comment corrections of
`COMMENT_FIXES` (each must be found once), and, with `--ownership`, one `EXPECT=held_by` line for each `r26.11.5`
stage before it that holds its contigs (the template's lines are kept as they are). It refuses an ownership report
with a problem or a discrepancy. `GV1.conf` is written from its template too (card R-23's re-staged files).

The member list is `R26_11_5_CHAIN`, one ordered tuple read by the ownership run, `write-confs` and the tests of
`tests/r1/test_stage_driver.py`; a stage added to `r26.11.5` later is a line appended there, its metadata directory in
`STAGING_DIRS` and its template.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:  # absent under `python -I`
    sys.path.append(str(_HERE))
import a26_r26_11_4_chain as a26  # noqa: E402  (card A-26's chain and helpers; not changed)

CARD_NOTE = "planning/measurements/2026-11-27_a29_r26_11_5_confs.md"
ROOT_STAGES: tuple[str, ...] = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
R26_11_4_CHAIN: tuple[str, ...] = tuple(s for s, _ in a26.R26_11_4_CHAIN)
R26_11_4_LAST = R26_11_4_CHAIN[-1]  # GART2
# DEC-79 (1)'s order of staging (plan 2026-11-26_r26_11_5_staging_plan.md §5), less ASG1 (no genome sequence in what
# was fetched; card C-82) and MGE4 (not mobile elements). DEC-85 appends RMD1 and GBRMGD after WET1 (card C-83).
R26_11_5_CHAIN: tuple[str, ...] = (
    "GRARE1", "GRARE2", "GRARE3",
    "LR2", "SAG1", "NX1",
    "GV1", "EGV1", "ECE2", "MIR1",
    "SAG2",
    "DPN1", "SYM2", "LRP1", "MIR2",
    "MV2", "MVEUKV",
    "KAVHIFI", "INPH2604",
    "GWV1", "SAG3", "WET1",
    "RMD1", "GBRMGD",
)
OUT_OF_R26_11_5: tuple[str, ...] = ("ASG1", "MGE4")
# Where each stage's staged metadata is on the pod (its genomes.json, objects.json, objects_staging.json); the notes
# of the staging cards, and for LR2, NX1 and GV1 the re-stagings of 2026-11-27 (cards R-19, R-20, R-23).
STAGING_DIRS: dict[str, str] = {
    "GRARE1": "/tmp/if2453/c75/meta/GRARE1",
    "GRARE2": "/tmp/if2453/c75/meta/GRARE2",
    "GRARE3": "/tmp/if2453/c75/meta/GRARE3",
    "LR2": "/tmp/if2453/r19/LR2",
    "SAG1": "/tmp/if2453/c58/SAG1",
    "NX1": "/tmp/if2453/r20/NX1",
    "GV1": "/tmp/if2453/r23/GV1",
    "EGV1": "/tmp/if2453/c69/EGV1",
    "ECE2": "/tmp/if2453/c70/meta/ECE2",
    "MIR1": "/tmp/if2453/c70/meta/MIR1",
    "SAG2": "/tmp/if2453/c77/meta/SAG2",
    "DPN1": "/tmp/if2453/c78/meta/DPN1",
    "SYM2": "/tmp/if2453/c78/meta/SYM2",
    "LRP1": "/tmp/if2453/c78/meta/LRP1",
    "MIR2": "/tmp/if2453/c78/meta/MIR2",
    "MV2": "/tmp/if2453/c79/MV2",
    "MVEUKV": "/tmp/if2453/c79/MVEUKV",
    "KAVHIFI": "/tmp/if2453/c80/KAVHIFI",
    "INPH2604": "/tmp/if2453/c80/INPH2604",
    "GWV1": "/tmp/if2453/c81/meta/GWV1",
    "SAG3": "/tmp/if2453/c81/meta/SAG3",
    "WET1": "/tmp/if2453/c81/meta/WET1",
    "RMD1": "/tmp/if2453/c83/meta/RMD1",
    "GBRMGD": "/tmp/if2453/c83/meta/GBRMGD",
}
# Comment corrections made in the promoted confs (old text found exactly once in the template, new text).
# LR2, NX1 and DPN1 templates carry the DEC-84 (3) deletion wording directly since card A-32.
COMMENT_FIXES: dict[str, tuple[tuple[str, str], ...]] = {
    "SAG2": (("# Build order (plan §5): after ASG1.",
              "# Build order (plan §5, card A-29): after MIR1; ASG1 is out of r26.11.5 (no genome sequence in what was\n"
              "# fetched; card C-82)."),),
    "MV2": (("# Build order (plan 2026-11-26_r26_11_5_staging_plan.md §5): after MGE4 (or MIR2 if MGE4 is not kept, "
             "card A-29).",
             "# Build order (plan 2026-11-26_r26_11_5_staging_plan.md §5, card A-29): after MIR2; MGE4 is out of "
             "r26.11.5\n# (not mobile elements)."),),
    "GV1": (("# Build order (DEC-71, DEC-72, cards A-24 and A-26): after the r26.11.4 chain, as stages/GV1.conf.",
             "# Build order (DEC-71, DEC-72, plan 2026-11-26_r26_11_5_staging_plan.md §5, card A-29): in the r26.11.5\n"
             "# chain, after NX1."),),
}
STAGE_MEASURED: dict[str, str] = {
    **{s: "measured, wilkcomp02, 2026-11-27, unit a29-own" for s in R26_11_5_CHAIN[:22]},
    "RMD1": "measured, wilkcomp02, 2026-11-27, unit a32-own",
    "GBRMGD": "measured, wilkcomp02, 2026-11-27, unit a32-own",
}
FLUSH = 200_000


def chain_after(stage: str, chain: tuple[str, ...] = R26_11_5_CHAIN) -> str:
    """The stage this one is built after: the one before it in the chain, GART2 for the first."""
    i = chain.index(stage)
    return chain[i - 1] if i else R26_11_4_LAST


def _sql(s) -> str:
    return "'" + str(s).replace("'", "''") + "'"


def _sha(path: Path) -> str:
    return a26.sha256_file(path)


def _fasta_parts(stage: str, objects: list[Path], out_dir: Path) -> dict:
    """One staged stage's FASTA: (contig_uid, length) of every record in Parquet parts, and each object's sha256."""
    import pyarrow as pa

    from biobox import fasta, ids

    sdir = out_dir / stage
    sdir.mkdir(parents=True, exist_ok=True)
    schema = pa.schema([("contig_uid", pa.binary(16)), ("n", pa.int64())])
    parts: list[str] = []
    rows: list[tuple] = []
    sums: dict[str, str] = {}
    records = 0
    bases = 0
    t0 = time.monotonic()
    for obj in objects:
        sums[obj.name] = _sha(obj)
        for rec in fasta.iter_records(obj):
            n = len(rec.seq)
            rows.append((ids.contig_uid(rec.seq), n))
            records += 1
            bases += n
            if len(rows) >= FLUSH:
                p = sdir / f"l-{len(parts):05d}.parquet"
                a26._write_rows(rows, schema, p)
                parts.append(str(p))
                rows = []
    p = sdir / f"l-{len(parts):05d}.parquet"
    a26._write_rows(rows, schema, p)
    parts.append(str(p))
    return {"stage": stage, "parts": parts, "sums": sums, "records": records, "bases": bases,
            "seconds": round(time.monotonic() - t0, 1)}


def _template(tree: Path, stage: str):
    import stage_est

    return stage_est.read_conf(tree / "scripts" / "stress" / "build_stage_pod" / "stages" / "templates"
                               / f"{stage}.conf", stage, tree=tree)


def staged_objects(meta: Path) -> list[Path]:
    return [Path(o["path"]) for o in json.loads((meta / "objects_staging.json").read_text(encoding="utf-8"))]


def check_files(stage: str, meta: Path, conf, object_sums: dict[str, str]) -> list[str]:
    """The template's SUM lines against the staged files; objects.json names INPUTS/<basename> of each object."""
    problems = []
    want = dict((n, h) for h, n in conf.sums)
    got = {"genomes.json": _sha(meta / "genomes.json"), "objects.json": _sha(meta / "objects.json"), **object_sums}
    for n in sorted(set(want) | set(got)):
        if want.get(n) != got.get(n):
            problems.append(f"{stage}: {n}: sha256 {got.get(n)} against the template's {want.get(n)}")
    final = [Path(o["path"]) for o in json.loads((meta / "objects.json").read_text(encoding="utf-8"))]
    staging = staged_objects(meta)
    if [p.name for p in final] != [p.name for p in staging]:
        problems.append(f"{stage}: objects.json and objects_staging.json name different files")
    for p in final:
        if p.parent != conf.inputs:
            problems.append(f"{stage}: objects.json names {p}, not a file of INPUTS {conf.inputs}")
    return problems


def compute(root: Path, tree: Path, inputs_home: Path, metas: dict[str, Path], tmp: Path, *, workers: int = 3,
            root_stages: tuple[str, ...] = ROOT_STAGES, r26_11_4: tuple[str, ...] = R26_11_4_CHAIN,
            chain: tuple[str, ...] = R26_11_5_CHAIN) -> dict:
    """Exact prefix ownership over root stages, then `r26_11_4` (genomes.json under inputs_home), then `chain`."""
    import duckdb

    from biobox import io
    from biobox.l0seq.store import INDEX_DIR

    tmp.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="a29-chain-", dir=tmp))
    doc: dict = {"root": str(root), "root_stages": list(root_stages), "r26_11_4": list(r26_11_4),
                 "chain": list(chain), "after": {s: chain_after(s, chain) for s in chain}, "problems": [],
                 "discrepancies": [], "stages": {}, "root_partitions_not_read": {}}
    try:
        confs = {s: _template(tree, s) for s in chain}
        jobs = sorted(chain, key=lambda s: -sum(p.stat().st_size for p in staged_objects(metas[s])))
        extracted: dict[str, dict] = {}
        if workers > 1:
            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as ex:
                futs = [ex.submit(_fasta_parts, s, staged_objects(metas[s]), work / "parts") for s in jobs]
                for f in concurrent.futures.as_completed(futs):
                    r = f.result()
                    extracted[r["stage"]] = r
        else:
            for s in jobs:
                extracted[s] = _fasta_parts(s, staged_objects(metas[s]), work / "parts")
        for s in chain:
            doc["problems"] += check_files(s, metas[s], confs[s], extracted[s]["sums"])

        (work / "duck").mkdir()
        con = duckdb.connect(config={"threads": 2, "memory_limit": "4GB", "temp_directory": str(work / "duck"),
                                     "preserve_insertion_order": False})
        try:
            base = root / INDEX_DIR
            present = sorted(d.name.removeprefix("stage=") for d in base.glob("stage=*") if d.is_dir())
            missing = [s for s in root_stages if s not in present]
            if missing:
                raise RuntimeError(f"root {root} has no contig_index for {missing}")
            sel = []
            for s in root_stages:
                files = io.part_files(base / f"stage={s}")
                sel.append(f"SELECT contig_uid, {_sql(s)} AS owner FROM read_parquet({a26._lit(files)}, "
                           "hive_partitioning = false) WHERE is_default_variant")
            con.execute("CREATE TEMP TABLE owners AS SELECT contig_uid, min(owner) AS owner, count(DISTINCT owner) "
                        "AS k FROM (" + " UNION ALL ".join(sel) + ") GROUP BY contig_uid")
            dup, = con.execute("SELECT count(*) FROM owners WHERE k > 1").fetchone()
            if dup:
                raise RuntimeError(f"root has {dup} contigs stored by more than one stage")
            con.execute("ALTER TABLE owners DROP COLUMN k")
            for s in present:
                if s in root_stages:
                    continue
                files = io.part_files(base / f"stage={s}")
                n, = con.execute(f"SELECT count(DISTINCT contig_uid) FROM read_parquet({a26._lit(files)}, "
                                 "hive_partitioning = false) WHERE is_default_variant").fetchone() if files else (0,)
                doc["root_partitions_not_read"][s] = {"files": len(files), "contigs": int(n)}

            def load_genomes(gj: Path) -> None:
                con.execute("CREATE OR REPLACE TEMP TABLE cur_g AS SELECT DISTINCT from_hex(u) AS contig_uid FROM "
                            f"(SELECT unnest(contig_uids) AS u FROM read_json({_sql(gj)}, format='array', "
                            "columns={'contig_uids': 'VARCHAR[]'}, maximum_object_size=1073741824))")

            def held(with_len: bool):
                if with_len:
                    rows = con.execute("SELECT o.owner, count(*)::BIGINT, sum(l.n)::BIGINT FROM cur_g g JOIN owners o "
                                       "USING (contig_uid) JOIN cur_l l USING (contig_uid) GROUP BY o.owner "
                                       "ORDER BY o.owner").fetchall()
                    return {r[0]: int(r[1]) for r in rows}, {r[0]: int(r[2]) for r in rows}
                rows = con.execute("SELECT o.owner, count(*)::BIGINT FROM cur_g g JOIN owners o USING (contig_uid) "
                                   "GROUP BY o.owner ORDER BY o.owner").fetchall()
                return {r[0]: int(r[1]) for r in rows}, {}

            def take(st: str) -> None:
                con.execute(f"INSERT INTO owners SELECT contig_uid, {_sql(st)} FROM cur_g ANTI JOIN owners "
                            "USING (contig_uid)")

            for st in r26_11_4:
                load_genomes(inputs_home / st / "genomes.json")
                tot, = con.execute("SELECT count(*) FROM cur_g").fetchone()
                by, _ = held(False)
                doc["stages"][st] = {"total_contigs": int(tot), "held_by": by, "held_contigs": sum(by.values()),
                                     "new_contigs": int(tot) - sum(by.values())}
                take(st)

            for st in chain:
                load_genomes(metas[st] / "genomes.json")
                con.execute(f"CREATE OR REPLACE TEMP TABLE cur_l AS SELECT contig_uid, max(n) AS n, "
                            f"count(DISTINCT n) AS k FROM read_parquet({a26._lit(extracted[st]['parts'])}) "
                            "GROUP BY contig_uid")
                two, = con.execute("SELECT count(*) FROM cur_l WHERE k > 1").fetchone()
                absent, = con.execute("SELECT count(*) FROM cur_g ANTI JOIN cur_l USING (contig_uid)").fetchone()
                unlisted, = con.execute("SELECT count(*) FROM cur_l ANTI JOIN cur_g USING (contig_uid)").fetchone()
                for n, what in ((two, "contigs with two lengths"), (absent, "listed contigs in no object"),
                                (unlisted, "records in no genome")):
                    if n:
                        doc["problems"].append(f"{st}: {n} {what}")
                tot, totb = con.execute("SELECT count(*)::BIGINT, coalesce(sum(l.n), 0)::BIGINT FROM cur_g g "
                                        "JOIN cur_l l USING (contig_uid)").fetchone()
                by, byb = held(True)
                hc, hb = sum(by.values()), sum(byb.values())
                ex = extracted[st]
                doc["stages"][st] = {"total_contigs": int(tot), "total_bases": int(totb), "held_by": by,
                                     "held_bases_by": byb, "held_contigs": hc, "held_bases": hb,
                                     "new_contigs": int(tot) - hc, "new_bases": int(totb) - hb,
                                     "records": ex["records"], "record_bases": ex["bases"],
                                     "read_seconds": ex["seconds"], "object_sha256": ex["sums"]}
                take(st)
                # the template's EXPECT lines against the counts, over the owners it could know
                earlier = {*root_stages, *r26_11_4, *chain[:chain.index(st)]}
                exp = {o: n for o, (n, _) in confs[st].expect.items()}
                for o in sorted(set(exp) | {o for o in by if o in root_stages or o in r26_11_4}):
                    if o not in earlier:
                        doc["discrepancies"].append({"stage": st, "owner": o, "template": exp[o], "found": None,
                                                     "why": "the template names a stage that is not before it"})
                    elif exp.get(o, 0) != by.get(o, 0):
                        doc["discrepancies"].append({"stage": st, "owner": o, "template": exp.get(o),
                                                     "found": by.get(o, 0)})
            n_own, = con.execute("SELECT count(*) FROM owners").fetchone()
            doc["owned_contigs"] = int(n_own)
        finally:
            con.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    doc["new_bases_r26_11_5"] = sum(doc["stages"][s]["new_bases"] for s in chain)
    doc["new_contigs_r26_11_5"] = sum(doc["stages"][s]["new_contigs"] for s in chain)
    return doc


def added_expect(ownership: dict, stage: str, template_expect: dict) -> list[tuple[str, int]]:
    """(owner, contigs) of each r26.11.5 stage before `stage` that holds its contigs and that the template lacks."""
    chain = ownership["chain"]
    by = ownership["stages"][stage]["held_by"]
    return [(o, by[o]) for o in chain[:chain.index(stage)] if by.get(o, 0) and o not in template_expect]


def promote(template_text: str, stage: str, *, after: str, note: str = CARD_NOTE,
            expect_lines: list[tuple[str, int]] | None = None, measured: str = "",
            fixes: tuple[tuple[str, str], ...] = ()) -> str:
    """The promoted conf's text from its template's: see the module's docstring."""
    text = template_text
    for old, new in fixes:
        if text.count(old) != 1:
            raise ValueError(f"{stage}: the comment to correct is not in the template once: {old!r}")
        text = text.replace(old, new)
    lines = text.splitlines()
    if sum(1 for ln in lines if ln.startswith("AFTER=")) != 1:
        raise ValueError(f"{stage}: the template has no single AFTER= line")
    out = [f"# Promoted by card A-29 from stages/templates/{stage}.conf: the AFTER= of the r26.11.5 chain, the chain's",
           f"# EXPECT lines and this card's note ({note}) added;",
           "# every other line is the template's, but for corrected comments."]
    last_note = max(i for i, ln in enumerate(lines) if ln.startswith("NOTE="))
    exp_idx = [i for i, ln in enumerate(lines) if ln.startswith("EXPECT=")]
    after_idx = next(i for i, ln in enumerate(lines) if ln.startswith("AFTER="))
    anchor = exp_idx[-1] if exp_idx else after_idx
    for i, ln in enumerate(lines):
        out.append(f"AFTER={after}" if i == after_idx else ln)
        if i == last_note and f"NOTE={note}" not in lines:
            out.append(f"NOTE={note}")
        if i == anchor and expect_lines is not None:
            if expect_lines:
                out.append(f"# Held by the r26.11.5 stages before it, by exact contig_uid over the chain (card A-29; "
                           f"{measured}):")
                out += [f"EXPECT=held_by {o} {n} {note}" for o, n in expect_lines]
            else:
                out.append(f"# No r26.11.5 stage before it holds a contig of it (card A-29; {measured}).")
    return "\n".join(out) + "\n"


def write_confs(tree: Path, ownership: dict | None, *, measured: str = "",
                chain: tuple[str, ...] = R26_11_5_CHAIN) -> list[Path]:
    import stage_est

    if ownership is not None:
        if ownership.get("problems") or ownership.get("discrepancies"):
            raise ValueError("the ownership report has problems or discrepancies: no conf is written from it")
        if tuple(ownership["chain"]) != tuple(chain):
            raise ValueError(f"the ownership report's chain {ownership['chain']} is not {list(chain)}")
    stages = tree / "scripts" / "stress" / "build_stage_pod" / "stages"
    written = []
    for st in chain:
        tpl = stages / "templates" / f"{st}.conf"
        conf = stage_est.read_conf(tpl, st)
        exp = None if ownership is None else added_expect(ownership, st, conf.expect)
        m_label = measured if measured else STAGE_MEASURED.get(st, "")
        text = promote(tpl.read_text(encoding="utf-8"), st, after=chain_after(st, chain), expect_lines=exp,
                       measured=m_label, fixes=COMMENT_FIXES.get(st, ()))
        dst = stages / f"{st}.conf"
        tmp = dst.with_name(f".{dst.name}.tmp")
        tmp.write_text(text, encoding="utf-8")
        os.replace(tmp, dst)
        stage_est.read_conf(dst, st, tree=tree)
        written.append(dst)
    return written


def _write_json(doc: dict, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f".{out.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("ownership")
    o.add_argument("--root", type=Path, required=True)
    o.add_argument("--tree", type=Path, required=True)
    o.add_argument("--out", type=Path, required=True)
    o.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    o.add_argument("--tmp", type=Path, default=Path("/tmp/if2453/a29/work"))
    o.add_argument("--workers", type=int, default=3)
    o.add_argument("--meta", action="append", default=[], help="STAGE=DIR in place of STAGING_DIRS[STAGE]")
    w = sub.add_parser("write-confs")
    w.add_argument("--tree", type=Path, required=True)
    w.add_argument("--ownership", type=Path, default=None)
    w.add_argument("--measured", default="")
    a = ap.parse_args(argv)
    sys.path.insert(0, str(a.tree / "scripts" / "stress" / "build_stage_pod"))
    if a.cmd == "ownership":
        metas = {s: Path(d) for s, d in STAGING_DIRS.items()}
        for kv in a.meta:
            s, _, d = kv.partition("=")
            metas[s] = Path(d)
        t0 = time.time()
        doc = compute(a.root, a.tree, a.inputs_home, metas, a.tmp, workers=a.workers)
        doc["argv"] = sys.argv[1:] if argv is None else list(argv)
        doc["wall_seconds"] = round(time.time() - t0, 1)
        doc["complete"] = True
        _write_json(doc, a.out)
        print(json.dumps({"out": str(a.out), "problems": len(doc["problems"]),
                          "discrepancies": len(doc["discrepancies"]),
                          "new_bases_r26_11_5": doc["new_bases_r26_11_5"]}, sort_keys=True))
        return 0 if not doc["problems"] and not doc["discrepancies"] else 1
    own = json.loads(a.ownership.read_text(encoding="utf-8")) if a.ownership else None
    for p in write_confs(a.tree, own, measured=a.measured):
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

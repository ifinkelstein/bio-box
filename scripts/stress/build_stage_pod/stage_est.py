"""Card A-12: the stage driver's conf, estimates before launch, live rate guard, completeness and rehearsal fixture.

One helper for every stage (`stage_driver.sh STAGE MODE`), made from card C-47's `c47_est.py` and card C-48c's
`c48c_v1_est.py`; what differs between stages is in `stages/<STAGE>.conf`, never in this file.

    python3 -I stage_est.py conf FILE --stage S [--sums OUT]            (host, standard library only)
    python stage_est.py before   --conf F --stage S --root R --profile P --inputs DIR --tmp DIR [--out-plan F]
    python stage_est.py rates    --conf F --stage S --root R --profile P --objects F --tmp DIR
    python stage_est.py complete --conf F --stage S --root R
    python stage_est.py fixture  --tree TREE --out DIR --stage S --earlier E [--min-contig-len N]

`conf` reads and checks a conf (`read_conf`) and prints its scalar values as `KEY=VALUE` lines, each value checked
against a strict character set, so the driver assigns them without evaluating anything; `--sums` writes the
`sha256sum -c` list of the inputs.

`before` (reads the root and the inputs; writes only --out-plan and its scratch under --tmp): the code's own
estimates before launch, from the inputs, bounded in memory (the records' lengths go through Parquet files and one
DuckDB connection): ingest (card B-52's `stage.ingest_cpu_min` over the FASTA bytes and the bases another stage
holds), gene calling (card B-48's `genecall.cut_plan` over the genomes as `genecall.stage_items` reads them from the
root after ingest: same rows, same `genecall._item`), the CRISPR producers (card B-53's `stage.crispr_cpu_min`, flanks
`2 x min(crispr.CIRC_FLANK, n)` per circular contig stored, as `stage.crispr_flank_bases`). Exit 5 when a stage of
the root that owns a contig of this stage has no gene-call record (a stage parked half built, as `C5` was when `V1`
was built: a release without it would miss calls); exit 7 when build-stage would refuse the stage (a source object
or a genome over the budget, card B-30, card B-48). The conf's JUDGE figures are printed beside the code's.

`rates`: DEC-43 (`MIN_JUDGED_CPU_MIN = 24`: a step estimated under it is reported, never judged) and DEC-53 (a
judged step more than 20% over its estimate stops the build; more than 20% under is recorded and does not stop)
over the stage's run records, every estimate made by the code from the root at the time, as build-stage makes it:
ingest as one group (every `ingest.*` record and `ingest`), each gene-call batch against the published plan and the
merge, gene calling in one process against `genecall.plan`, each scan shard against its residues, each CRISPR
producer. A conf's `JUDGE=<step> <figure>` judges that step's group against the note's figure as well (card C-48c
part 3: `GV1`'s gene calling against 236 CPU-minutes); a group still running is judged over on its partial sum, never
under. A conf's `ACCEPTED=<step> <cpu> <estimate> <DEC-n>` reports that one record, at exactly those values, under
`accepted_over_20pct` instead of stopping (as `c47_est.py` did for DEC-56): only a decision row lifts a guard, and
only the PI or an H card adds one. Exit 3 to stop.

`complete`: `stage.stage_incomplete` (the release train's own check): every step of the stage's depth, libraries
and CRISPR choice has a run record with exit 0 and its outputs at their sizes, and no `.publishing.json` is left.

`fixture`: the rehearsal's inputs from `tests/fixtures/micro_s0` (card C-24's fixture): stage E (lambda, M.
genitalium), and stage S: T4; lambda again (every contig held by E: DEC-26, not called); T4 with a 20 kb slice of
itself (T4 shared within the stage: DEC-27, the case card C-47's check first missed); a 2 kb slice of M. genitalium,
listed only when the least contig length allows it, else read and dropped by ingest (card C-38d).
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

TOL = 0.20
MIN_JUDGED_CPU_MIN = 24.0  # DEC-43
DEPTHS = ("called", "annotated")  # `stored` builds no gene calls: not a stage this driver builds or checks
LIBRARIES = {("pfam", "38.2")}  # the libraries whose HMM file the driver knows (stage_driver.sh PFAM)
CUTOFFS = ("gathering", "trusted", "noise", "none")  # stage.CUTOFFS
JUDGE_STEPS = ("ingest", "genecall")  # the steps the guard judges as one group of records
SCALARS = ("STAGE", "INPUTS", "DEPTH", "MIN_CONTIG_LEN", "LIBRARY", "CRISPR", "AFTER")
LISTS = ("SUM", "JUDGE", "ACCEPTED", "NOTE")

STAGE_RE = re.compile(r"[A-Z][A-Za-z0-9]{0,15}")
PATH_RE = re.compile(r"/[A-Za-z0-9_.+/-]+")
NOTE_RE = re.compile(r"planning/measurements/[A-Za-z0-9_.+/-]+\.md")
SUM_RE = re.compile(r"([0-9a-f]{64}) ([A-Za-z0-9_.+-]+)")
STEP_RE = r"[A-Za-z0-9][A-Za-z0-9_.-]*"
JUDGE_RE = re.compile(rf"({STEP_RE}) ([0-9]+(?:\.[0-9]+)?) (\S+)")
ACCEPTED_RE = re.compile(rf"({STEP_RE}) ([0-9]+\.[0-9]{{3}}) ([0-9]+\.[0-9]{{3}}) (DEC-[0-9]+)")


class ConfError(ValueError):
    """A conf that the driver must not run on."""


@dataclass(frozen=True)
class Conf:
    stage: str
    inputs: Path
    sums: tuple[tuple[str, str], ...]  # (sha256 hex, file name in INPUTS)
    depth: str
    crispr: bool
    min_contig_len: int | None = None
    library: tuple[str, str, str] | None = None  # (name, version, cutoff)
    after: tuple[str, ...] = ()
    judge: dict = field(default_factory=dict)  # step: (figure, source)
    accepted: dict = field(default_factory=dict)  # (step, cpu_min, estimate): decision row
    notes: tuple[str, ...] = ()


def read_conf(path: str | Path, stage: str | None = None) -> Conf:
    """The conf at `path`, checked line by line: `KEY=VALUE`, `#` comments and blank lines; the scalar keys once, the
    list keys (SUM, JUDGE, ACCEPTED, NOTE) any number of times; every value against a strict pattern (no space where
    none belongs, no shell metacharacter, no `..`). ConfError says what is wrong."""
    path = Path(path)
    vals: dict[str, str] = {}
    lists: dict[str, list[str]] = {k: [] for k in LISTS}
    for n, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        key, eq, val = line.partition("=")
        if not eq:
            raise ConfError(f"{path}:{n}: not KEY=VALUE: {raw!r}")
        if key in LISTS:
            lists[key].append(val)
        elif key in SCALARS:
            if key in vals:
                raise ConfError(f"{path}:{n}: {key} given twice")
            vals[key] = val
        else:
            raise ConfError(f"{path}:{n}: unknown key {key!r} (keys: {', '.join(SCALARS + LISTS)})")
    st = vals.get("STAGE", "")
    if not STAGE_RE.fullmatch(st):
        raise ConfError(f"{path}: STAGE {st!r}: expected a capital letter, then letters and digits")
    if stage is not None and st != stage:
        raise ConfError(f"{path}: STAGE {st} is not the stage asked for, {stage}")
    if path.suffix == ".conf" and path.stem != st:
        raise ConfError(f"{path}: the file is not named for stage {st}")
    inputs = vals.get("INPUTS", "")
    if not PATH_RE.fullmatch(inputs) or ".." in inputs.split("/"):
        raise ConfError(f"{path}: INPUTS {inputs!r}: expected an absolute path of letters, digits and . _ + - /")
    sums = []
    for v in lists["SUM"]:
        m = SUM_RE.fullmatch(v)
        if m is None or m[2] in (".", ".."):
            raise ConfError(f"{path}: SUM {v!r}: expected '<64 hex digits> <file name in INPUTS>'")
        sums.append((m[1], m[2]))
    names = [s[1] for s in sums]
    if len(set(names)) != len(names):
        raise ConfError(f"{path}: a file is summed twice")
    for need in ("genomes.json", "objects.json"):
        if need not in names:
            raise ConfError(f"{path}: no SUM of {need} (the note records it; the driver checks it before launch)")
    depth = vals.get("DEPTH", "")
    if depth not in DEPTHS:
        raise ConfError(f"{path}: DEPTH {depth!r}: expected one of {', '.join(DEPTHS)}")
    crispr = vals.get("CRISPR", "")
    if crispr not in ("yes", "no"):
        raise ConfError(f"{path}: CRISPR {crispr!r}: expected yes or no")
    minlen = vals.get("MIN_CONTIG_LEN", "")
    if minlen and not re.fullmatch(r"[1-9][0-9]{0,8}", minlen):
        raise ConfError(f"{path}: MIN_CONTIG_LEN {minlen!r}: expected a positive integer")
    lib = None
    if vals.get("LIBRARY"):
        parts = vals["LIBRARY"].split(" ")
        if len(parts) != 3 or (parts[0], parts[1]) not in LIBRARIES or parts[2] not in CUTOFFS:
            raise ConfError(f"{path}: LIBRARY {vals['LIBRARY']!r}: expected '<name> <version> <cutoff>' of a library "
                            f"the driver knows ({sorted(LIBRARIES)}) and a cutoff of {', '.join(CUTOFFS)}")
        lib = (parts[0], parts[1], parts[2])
    if depth == "annotated" and lib is None:
        raise ConfError(f"{path}: DEPTH annotated needs a LIBRARY (build-stage refuses it without one)")
    if depth == "called" and lib is not None:
        raise ConfError(f"{path}: DEPTH called takes no LIBRARY (build-stage refuses it)")
    after: tuple[str, ...] = ()
    if vals.get("AFTER"):
        after = tuple(vals["AFTER"].split(","))
        if not all(STAGE_RE.fullmatch(a) for a in after) or st in after or len(set(after)) != len(after):
            raise ConfError(f"{path}: AFTER {vals['AFTER']!r}: expected other stages' names, comma-separated")
    judge = {}
    for v in lists["JUDGE"]:
        m = JUDGE_RE.fullmatch(v)
        if m is None or m[1] not in JUDGE_STEPS or m[1] in judge:
            raise ConfError(f"{path}: JUDGE {v!r}: expected '<{'|'.join(JUDGE_STEPS)}> <CPU-minutes> <the note that "
                            "states the figure>', once per step")
        judge[m[1]] = (float(m[2]), m[3])
    accepted = {}
    for v in lists["ACCEPTED"]:
        m = ACCEPTED_RE.fullmatch(v)
        if m is None:
            raise ConfError(f"{path}: ACCEPTED {v!r}: expected '<step> <CPU-minutes of the record, 3 decimals> "
                            "<its estimate, 3 decimals> <DEC-n, the decision row that accepts it>'")
        accepted[(m[1], float(m[2]), float(m[3]))] = m[4]
    notes = []
    for v in lists["NOTE"]:
        if not NOTE_RE.fullmatch(v) or ".." in v.split("/"):
            raise ConfError(f"{path}: NOTE {v!r}: expected planning/measurements/<file>.md")
        notes.append(v)
    if not notes:
        raise ConfError(f"{path}: no NOTE: every value comes from a stage's note")
    return Conf(stage=st, inputs=Path(inputs), sums=tuple(sums), depth=depth, crispr=crispr == "yes",
                min_contig_len=int(minlen) if minlen else None, library=lib, after=after, judge=judge,
                accepted=accepted, notes=tuple(notes))


def conf_cmd(a) -> int:
    try:
        c = read_conf(a.file, a.stage)
    except (ConfError, OSError) as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    if a.sums:
        Path(a.sums).write_text("".join(f"{h}  {n}\n" for h, n in c.sums))
    lib = c.library or ("", "", "")
    out = {"STAGE": c.stage, "INPUTS": str(c.inputs), "DEPTH": c.depth,
           "MIN_CONTIG_LEN": "" if c.min_contig_len is None else str(c.min_contig_len), "LIBNAME": lib[0],
           "LIBVER": lib[1], "LIBCUT": lib[2], "CRISPR": "yes" if c.crispr else "no", "AFTER": " ".join(c.after),
           "NOTES": " ".join(c.notes)}
    for k, v in out.items():
        print(f"{k}={v}")
    return 0


# the rate guard (pure) -------------------------------------------------------------------------------------------


def judge_items(items: list[dict], accepted: dict) -> dict:
    """Each item ({name, cpu_min, estimate, basis, complete, records}) under DEC-43 and DEC-53: judged when its estimate
    is at least `MIN_JUDGED_CPU_MIN`; over when judged and above 1 + TOL of it (on a partial sum too: a sum only
    grows); under when judged, complete and below 1 - TOL (recorded, never a stop). An over item whose (name, CPU,
    estimate), at 3 decimals, is in `accepted` is reported under `accepted_over_20pct` with its decision row."""
    rows, over, under, acc = [], [], [], {}
    for it in items:
        cpu, est = float(it["cpu_min"]), float(it["estimate"])
        ratio = cpu / est if est > 0 else math.inf
        judged = est >= MIN_JUDGED_CPU_MIN
        is_over = judged and ratio > 1.0 + TOL
        is_under = judged and it.get("complete", True) and ratio < 1.0 - TOL
        row = {"step": it["name"], "cpu_min": round(cpu, 3), "estimate": round(est, 3),
               "ratio": round(ratio, 3) if math.isfinite(ratio) else None, "basis": it.get("basis", ""),
               "judged": judged, "complete": it.get("complete", True), "records": it.get("records", 1),
               "over_20pct": is_over, "under_20pct": is_under}
        dec = accepted.get((it["name"], round(cpu, 3), round(est, 3)))
        if is_over and dec is not None and it.get("complete", True):
            row["accepted"] = dec
            row["over_20pct"] = is_over = False
            acc[it["name"]] = dec
        rows.append(row)
        if is_over:
            over.append(it["name"])
        if is_under:
            under.append(it["name"])
    return {"checked": len(rows), "over_20pct": over, "under_20pct": under, "accepted_over_20pct": acc, "rows": rows}


def rates_exit(doc: dict) -> int:
    """The driver stops its unit on 3: a judged step more than 20% over its estimate (DEC-53)."""
    return 3 if doc["over_20pct"] else 0


def group_names(records: dict, step: str) -> list[str]:
    """The run records of `step` run in pieces: `<step>` itself and every `<step>.<piece>`."""
    return sorted(n for n in records if n == step or n.startswith(step + "."))


# reading the root, in the image --------------------------------------------------------------------------------


def _records(root: Path, stage: str) -> dict:
    d = root / "runs" / f"stage={stage}"
    if not d.is_dir():
        return {}
    return {p.name[:-len(".json")]: json.loads(p.read_text())
            for p in sorted(d.glob("*.json")) if not p.name.endswith(".publishing.json")}


def _lit(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def _duck(tmp: Path, mem: str = "4GB"):
    import duckdb

    tmp.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(config={"threads": 2, "memory_limit": mem, "temp_directory": str(tmp),
                                  "preserve_insertion_order": False})


def rates(a) -> int:
    from biobox import profile as bprofile
    from biobox import stage as bs
    from biobox.l0seq import store
    from biobox.l1 import crispr, genecall
    from biobox.l1.genecall import GENE_CALL_VERSION as GCV

    conf = read_conf(a.conf, a.stage)
    root, st, tmp = Path(a.root), conf.stage, Path(a.tmp)
    prof = bprofile.load(a.profile)
    recs = _records(root, st)
    items, covered = [], set()

    def cpu_of(names):
        return sum(float(recs[n]["cpu_min"]) for n in names)

    def single(name: str, est: float, basis: str) -> None:
        if name in recs:
            items.append({"name": name, "cpu_min": cpu_of([name]), "estimate": est, "basis": basis, "complete": True,
                          "records": 1})
            covered.add(name)

    def group(label: str, step: str, est: float, basis: str) -> None:
        names = group_names(recs, step)
        if names:
            items.append({"name": label, "cpu_min": cpu_of(names), "estimate": est, "basis": basis,
                          "complete": step in recs, "records": len(names)})
            covered.update(names)

    if group_names(recs, "ingest"):
        # the FASTA bytes `estimate` measured with the code's own `stage.object_sizes` (it decompresses every object:
        # minutes for a stage of gzip inputs, so not again at every record); measured here only without it
        doc = json.loads(Path(a.before).read_text()) if a.before and Path(a.before).is_file() else {}
        n_bytes = doc.get("fasta_bytes")
        if not isinstance(n_bytes, int):
            n_bytes = sum(f for _, _, f in bs.object_sizes(Path(a.objects)))
        held = bs.held_bases(root, st, tmp)
        group("ingest", "ingest", bs.ingest_cpu_min(n_bytes, held),
              f"stage.ingest_cpu_min({n_bytes}, held={held}) over every ingest record")
    gp_path = root / genecall.plan_record_rel(st)
    if gp_path.exists():
        gp = json.loads(gp_path.read_text())
        for b in gp["batches"]:
            single(f"genecall.{genecall.batch_name(b['batch'])}", b["estimate_cpu_min"],
                   f"plan: {b['genomes']} genomes, {b['bases']} bases called")
        single("genecall", gp["merge_estimate_cpu_min"], f"plan: merge of {gp['genes_expected']} genes expected")
    elif "genecall" in recs:
        seq = [st, *(s for s in store.stages_on_disk(root) if s != st)]
        p = genecall.plan(root, st, version=GCV, seq_stages=seq, batch_cpu_min=bs.SHARD_MARGIN *
                          prof.cpu_min_per_process, tmp=tmp)
        single("genecall", p["estimate_cpu_min"], "genecall.plan in one process (card B-48)")
    for step, (fig, src) in sorted(conf.judge.items()):
        group(f"{step}:note", step, fig, f"the note's figure ({src}) over every {step} record")
    if conf.library is not None:
        lib, ver = conf.library[0], conf.library[1]
        name = f"domains-{lib}-{ver}"
        plan_p = root / "l1" / "domain_scan" / f"library={lib}-{ver}" / f"stage={st}" / "plan.json"
        if plan_p.exists():
            for k, s in enumerate(json.loads(plan_p.read_text())["shards"]):
                single(f"{name}.shard-{k:05d}", bs.scan_cpu_min(s["residues"]),
                       f"{s['residues']} residues x {bs.SCAN_CPU_MIN_PER_M} / 1e6")
        elif name in recs:
            r = bs.stage_residues(root, st, GCV)
            single(name, bs.scan_cpu_min(r), f"{r} residues in one process x {bs.SCAN_CPU_MIN_PER_M} / 1e6")
    if "ingest" in recs and conf.crispr:
        b, fb = bs.stage_bases(root, st), bs.crispr_flank_bases(root, st)
        single(f"crispr-{crispr.producer()}", bs.crispr_cpu_min(b, fb, near_cas=False),
               f"({b}+{fb}) bases x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60")
        if conf.library is not None and conf.depth == "annotated" and conf.library[0] == crispr.CAS_LIBRARY:
            single(f"crispr-{crispr.near_cas_producer(conf.library[1])}", bs.crispr_cpu_min(b, fb, near_cas=True),
                   f"({b}+{fb}) bases x {bs.CRISPR_CPU_S_PER_GBP} x {bs.CRISPR_NEAR_CAS_COST} / 1e9 / 60")
    doc = {"stage": st, **judge_items(items, conf.accepted),
           "not_judged": sorted(n for n in recs if n not in covered),
           "guard_policy": "DEC-43 (estimate >= 24 CPU-min judged) + DEC-53 (stop above 1.20; record below 0.80)"}
    print(json.dumps(doc, sort_keys=True))
    return rates_exit(doc)


def complete(a) -> int:
    from biobox import stage as bs

    conf = read_conf(a.conf, a.stage)
    libs = [] if conf.library is None else [(conf.library[0], conf.library[1])]
    problems = bs.stage_incomplete(Path(a.root), conf.stage, conf.depth, libraries=libs, crispr=conf.crispr)
    print(json.dumps({"stage": conf.stage, "depth": conf.depth, "problems": problems}, sort_keys=True))
    return 1 if problems else 0


# the estimate before launch --------------------------------------------------------------------------------------

_FLUSH = 200_000


def _write_rows(rows: list[tuple], schema, path: Path) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    cols = list(zip(*rows)) if rows else [[] for _ in schema.names]
    pq.write_table(pa.table({n: list(c) for n, c in zip(schema.names, cols)}, schema=schema), path)


def before(a) -> int:
    import pyarrow as pa

    from biobox import fasta, ids, io
    from biobox import profile as bprofile
    from biobox import stage as bs
    from biobox.l0seq.store import INDEX_DIR
    from biobox.l1 import crispr, genecall

    conf = read_conf(a.conf, a.stage)
    root, inp, st = Path(a.root), Path(a.inputs), conf.stage
    prof = bprofile.load(a.profile)
    budget = prof.cpu_min_per_process
    Path(a.tmp).mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="stage-est-", dir=a.tmp))
    try:
        sizes = bs.object_sizes(inp / "objects.json")
        fasta_bytes, disk_bytes = sum(f for _, _, f in sizes), sum(d for _, d, _ in sizes)
        big = max(sizes, key=lambda t: (t[2], str(t[0])), default=None)
        # the genomes, read as the inputs step reads them (cli._genomes), one row per membership into Parquet
        gschema = pa.schema([("genome_uid", pa.binary(16)), ("contig_uid", pa.binary(16)), ("circular", pa.bool_()),
                             ("kind", pa.string())])
        gparts, rows, n_genomes = [], [], 0
        for g in json.loads((inp / "genomes.json").read_text()):
            cus = [bytes.fromhex(u) for u in g["contig_uids"]]
            circ = {bytes.fromhex(u) for u in g.get("circular", ())}
            guid = ids.genome_uid(cus)
            n_genomes += 1
            rows += [(guid, cu, cu in circ, g["kind"]) for cu in cus]
            if len(rows) >= _FLUSH:
                gparts.append(work / f"g-{len(gparts):05d}.parquet")
                _write_rows(rows, gschema, gparts[-1])
                rows = []
        gparts.append(work / f"g-{len(gparts):05d}.parquet")
        _write_rows(rows, gschema, gparts[-1])
        # the records' lengths, streamed: a record under the least length is counted and dropped (ingest's rule)
        lschema = pa.schema([("contig_uid", pa.binary(16)), ("n", pa.int64())])
        lparts, rows = [], []
        n_records = short = short_bases = 0
        for p, _, _ in sizes:
            for r in fasta.iter_records(p):
                n_records += 1
                n = len(r.seq)
                if conf.min_contig_len is not None and n < conf.min_contig_len:
                    short, short_bases = short + 1, short_bases + n
                    continue
                rows.append((ids.contig_uid(r.seq), n))
                if len(rows) >= _FLUSH:
                    lparts.append(work / f"l-{len(lparts):05d}.parquet")
                    _write_rows(rows, lschema, lparts[-1])
                    rows = []
        lparts.append(work / f"l-{len(lparts):05d}.parquet")
        _write_rows(rows, lschema, lparts[-1])
        con = _duck(work / "duck")
        con.execute(f"CREATE VIEW gin AS SELECT * FROM read_parquet({_lit(gparts)})")
        con.execute(f"CREATE TEMP TABLE lens AS SELECT contig_uid, max(n) AS n, count(DISTINCT n) AS k "
                    f"FROM read_parquet({_lit(lparts)}) GROUP BY contig_uid")
        # every result is fetched whole (one row): DuckDB refuses the streaming reader below after a half-read one
        two, = con.execute("SELECT count(*) FROM lens WHERE k > 1").fetchall()[0]
        missing, = con.execute("SELECT count(DISTINCT contig_uid) FROM gin ANTI JOIN lens USING (contig_uid)"
                               ).fetchall()[0]
        if two or missing:
            print(f"refused: {two} contigs with two lengths, {missing} listed contigs with no record "
                  f"{'at or over the least length' if conf.min_contig_len else ''} in the objects", file=sys.stderr)
            return 6
        others = []
        base = root / INDEX_DIR
        stages = sorted(d.name.removeprefix("stage=") for d in base.glob("stage=*") if d.is_dir()) \
            if base.is_dir() else []
        for s in stages:
            files = io.part_files(base / f"stage={s}")
            if s != st and files:
                lit = "'" + s.replace("'", "''") + "'"
                others.append(f"SELECT contig_uid, {lit} AS st FROM read_parquet({_lit(files)}, "
                              "hive_partitioning = false) WHERE is_default_variant")
        held_sql = " UNION ALL ".join(others) if others else "SELECT NULL::BLOB AS contig_uid, NULL AS st LIMIT 0"
        con.execute(f"CREATE TEMP TABLE held AS SELECT contig_uid, min(st) AS st, count(DISTINCT st) AS m "
                    f"FROM ({held_sql}) SEMI JOIN (SELECT DISTINCT contig_uid FROM gin) USING (contig_uid) "
                    "GROUP BY contig_uid")
        con.execute("CREATE TEMP TABLE dist AS SELECT contig_uid, bool_or(circular) AS circ FROM gin "
                    "GROUP BY contig_uid")
        held_bases, held_n = con.execute("SELECT coalesce(sum(n), 0)::BIGINT, count(*) FROM held JOIN lens "
                                         "USING (contig_uid)").fetchall()[0]
        twice, = con.execute("SELECT count(*) FROM held WHERE m > 1").fetchall()[0]
        owned_by = dict(con.execute("SELECT st, count(*) FROM held GROUP BY st ORDER BY st").fetchall())
        n_stored, b_stored, flank, n_circ = con.execute(
            f"SELECT count(*), coalesce(sum(n), 0)::BIGINT, "
            f"coalesce(sum(CASE WHEN circ THEN 2 * least({crispr.CIRC_FLANK}, n) ELSE 0 END), 0)::BIGINT, "
            "count(*) FILTER (WHERE circ) FROM dist JOIN lens USING (contig_uid) ANTI JOIN held USING (contig_uid)"
        ).fetchall()[0]
        items_sql = ("SELECT g.genome_uid, g.contig_uid, g.circular, g.kind, l.n, "
                     "(h.contig_uid IS NULL AND m.g = g.genome_uid) AS mine FROM gin g JOIN lens l USING (contig_uid) "
                     "LEFT JOIN held h USING (contig_uid) "
                     "JOIN (SELECT contig_uid, min(genome_uid) AS g FROM gin GROUP BY contig_uid) m USING (contig_uid) "
                     "ORDER BY g.genome_uid, g.contig_uid")
        seq = [st, *[s for s in stages if s != st]]

        def items():  # the rows `genecall.stage_items` streams after ingest, grouped by `genecall._item`
            reader = con.sql(items_sql).to_arrow_reader(io.ROW_GROUP_ROWS)
            guid, grp = None, []
            for batch in reader:
                for r in zip(*(batch.column(i).to_pylist() for i in range(6))):
                    if r[0] != guid:
                        if grp:
                            yield genecall._item(root, seq, guid, grp)
                        guid, grp = r[0], []
                    grp.append(r)
            if grp:
                yield genecall._item(root, seq, guid, grp)

        plan = genecall.cut_plan(items, bs.SHARD_MARGIN * budget)
        con.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    plan.update({"stage": st, "gene_call_version": genecall.GENE_CALL_VERSION, "seq_stages": seq})
    if a.out_plan:
        Path(a.out_plan).write_text(json.dumps(plan, indent=1, sort_keys=True, default=str))
    gc_refusal = bs.genecall_refusal(prof, plan)
    ingest_refusal = None
    if big is not None and bs.ingest_cpu_min(big[2]) > budget:
        ingest_refusal = (f"{big[0].name} holds {big[2]:,} bytes of FASTA, estimated at {bs.ingest_cpu_min(big[2]):,.1f} "
                          f"CPU-minutes in one ingest process, over cpu_min_per_process {budget}")
    parked = sorted(s for s in owned_by if not (root / "runs" / f"stage={s}" / "genecall.json").exists())
    doc = {
        "stage": st, "depth": conf.depth, "min_contig_len": conf.min_contig_len,
        "library": None if conf.library is None else " ".join(conf.library), "crispr": conf.crispr,
        "fasta_bytes": fasta_bytes, "disk_bytes": disk_bytes, "records": n_records, "short_records": short,
        "short_bases": short_bases, "genomes": n_genomes, "contigs_listed": n_stored + held_n,
        "held_contigs": held_n, "held_bases": held_bases, "owned_by_stage": owned_by,
        "contigs_held_by_two_stages": twice, "contigs_stored": n_stored, "bases_stored": b_stored,
        "circular_stored": n_circ, "flank_bases": flank,
        "ingest_cpu_min": round(bs.ingest_cpu_min(fasta_bytes, held_bases), 3),
        "ingest_formula": f"{bs.INGEST_FIXED_CPU_MIN} + max(0, {fasta_bytes} - {held_bases}) x {bs.INGEST_CPU_MIN} / "
                          f"{bs.INGEST_BP}",
        "ingest_in_batches": bs.ingest_cpu_min(fasta_bytes, held_bases) > bs.SHARD_MARGIN * budget,
        "genecall": {k: plan[k] for k in ("genomes", "genomes_not_called", "bases", "bases_called", "estimate_cpu_min",
                                          "merge_estimate_cpu_min", "genes_expected", "largest_genome")}
                    | {"batches": [{k: b[k] for k in ("batch", "genomes", "bases", "estimate_cpu_min")}
                                   for b in plan["batches"]]},
        "margin_cpu_min": bs.SHARD_MARGIN * budget,
        "judge": {k: {"figure": f, "source": s} for k, (f, s) in sorted(conf.judge.items())},
        "accepted": [{"step": k[0], "cpu_min": k[1], "estimate": k[2], "decision": v}
                     for k, v in sorted(conf.accepted.items())],
        "parked_owners": parked, "refusals": [x for x in (ingest_refusal, gc_refusal) if x],
        "guard_policy": "DEC-43 (estimate >= 24 CPU-min judged) + DEC-53 (stop above 1.20; record below 0.80)",
    }
    if conf.crispr:
        doc["crispr_default_cpu_min"] = round(bs.crispr_cpu_min(b_stored, flank, near_cas=False), 3)
        if conf.library is not None and conf.library[0] == crispr.CAS_LIBRARY:
            doc["crispr_near_cpu_min"] = round(bs.crispr_cpu_min(b_stored, flank, near_cas=True), 3)
    print(json.dumps(doc, indent=1, sort_keys=True, default=str))
    if parked:
        print(f"refused: contigs of {st} are owned by {', '.join(parked)}, which has no gene-call record (a stage "
              "parked half built): a release cut without it would miss their calls", file=sys.stderr)
        return 5
    if twice:
        print(f"refused: {twice} contigs of {st} are held by two stages (DEC-26 allows one)", file=sys.stderr)
        return 5
    if doc["refusals"]:
        print("refused: build-stage would refuse the stage: " + "; ".join(doc["refusals"]), file=sys.stderr)
        return 7
    return 0


# the rehearsal's fixture -----------------------------------------------------------------------------------------


def fixture(a) -> int:
    import yaml

    from biobox import fasta, ids

    tree, out = Path(a.tree), Path(a.out)
    fix = tree / "tests" / "fixtures" / "micro_s0"
    contract = yaml.safe_load((fix / "sources.yaml").read_text())
    by_acc = {o["accession"]: o for o in contract["objects"]}
    lam, t4, mg = by_acc["NC_001416.1"], by_acc["NC_000866.4"], by_acc["NC_000908.2"]
    seq = {k: fasta.read(fix / o["file"])[0] for k, o in (("lam", lam), ("t4", t4), ("mg", mg))}

    def obj(path: Path) -> dict:
        return {"path": str(path), "source_name": contract["source_name"],
                "source_release": str(contract["source_release"]), "license": contract["license"],
                "redistributable": bool(contract["redistributable"])}

    def genome(cus: list[bytes], o: dict, circular: list[bytes] = ()) -> dict:
        return {"contig_uids": [c.hex() for c in cus], "genetic_code": o["genetic_code"], "kind": o["kind"],
                "circular": [c.hex() for c in circular], "assembly_quality": "complete"}

    def cu(s: bytes) -> bytes:
        return ids.contig_uid(s)

    e, s = out / "E", out / "R"
    e.mkdir(parents=True, exist_ok=False)
    s.mkdir(parents=True, exist_ok=False)
    eg = [genome([cu(seq["lam"].seq)], lam), genome([cu(seq["mg"].seq)], mg, [cu(seq["mg"].seq)])]
    (e / "genomes.json").write_text(json.dumps(eg, indent=1))
    (e / "objects.json").write_text(json.dumps([obj(fix / lam["file"]), obj(fix / mg["file"])], indent=1))
    t4s, slice_, short = seq["t4"].seq, seq["t4"].seq[:20_000], seq["mg"].seq[:2_000]
    recs = [("t4", t4s), ("t4_slice_20000", slice_), ("lambda", seq["lam"].seq), ("mg_slice_2000", short)]
    with open(s / "stage.fna", "wb") as fh:
        for h, x in recs:
            fh.write(b"\n".join(fasta.format_record(h, x, 70)) + b"\n")
    keep = (lambda n: a.min_contig_len is None or n >= a.min_contig_len)
    sg = [genome([cu(t4s)], t4), genome([cu(seq["lam"].seq)], lam)]
    if keep(len(slice_)):
        sg.append(genome([cu(t4s), cu(slice_)], t4))
    if keep(len(short)):
        sg.append(genome([cu(short)], mg))
    (s / "genomes.json").write_text(json.dumps(sg, indent=1))
    (s / "objects.json").write_text(json.dumps([obj(s / "stage.fna")], indent=1))
    print(json.dumps({"earlier": a.earlier, "stage": a.stage, "earlier_genomes": len(eg), "stage_genomes": len(sg),
                      "stage_records": len(recs), "min_contig_len": a.min_contig_len}, sort_keys=True))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("conf")
    p.add_argument("file")
    p.add_argument("--stage", required=True)
    p.add_argument("--sums")
    for name in ("before", "rates", "complete"):
        p = sub.add_parser(name)
        p.add_argument("--conf", required=True)
        p.add_argument("--stage", required=True)
        p.add_argument("--root", required=True)
        if name != "complete":
            p.add_argument("--profile", required=True)
            p.add_argument("--tmp", required=True)
        if name == "before":
            p.add_argument("--inputs", required=True)
            p.add_argument("--out-plan")
        if name == "rates":
            p.add_argument("--objects", required=True)
            p.add_argument("--before", help="estimate/before.json: its fasta_bytes")
    p = sub.add_parser("fixture")
    p.add_argument("--tree", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--stage", required=True)
    p.add_argument("--earlier", required=True)
    p.add_argument("--min-contig-len", type=int, default=None)
    a = ap.parse_args()
    return {"conf": conf_cmd, "before": before, "rates": rates, "complete": complete, "fixture": fixture}[a.mode](a)


if __name__ == "__main__":
    sys.exit(main())

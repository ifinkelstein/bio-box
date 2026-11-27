"""Card C-75 (DEC-79): stage GlobDB r226's rare phyla, unclassified domain and remaining archaea as GRARE1..3.

    python c75_grare_inputs.py select NON_GTDB_STAGES_TSV TAX_PLUS_STATS_TSV_GZ REKEY_TSV OUT_DIR
        [--minor-max-genomes 1000] [--n-stages 3] [--max-stage-bp 12000000000]
    python c75_grare_inputs.py stream SELECT_DIR GENOME_FASTA_TAR_GZ META_OUT FASTA_OUT [--bgzip] [--cpu-guard-min 280]
    python c75_grare_inputs.py screen SELECT_DIR META_OUT FASTA_OUT [--workers 8] [--screen-max-bp 500000]
    python c75_grare_inputs.py finish SELECT_DIR META_OUT FASTA_OUT [--root ROOT] [--root-stages S1,...,SPX1]
        [--chain-home /tmp/if2453/s2/inputs] [--inputs-home /tmp/if2453/s2/inputs] [--work DIR]

Selection (`select`, parameters stated in the note): the GlobDB r226 non-GTDB species representatives that no stage
holds (not `G1`, not a phylum of `GART1`/`GART2`: `r232_index.unbuilt_ids`, card C-74) that are archaeal, of the
unclassified domain, in a phylum with no GTDB genome in GlobDB r226, or in a bacterial phylum with fewer than
`minor_max_genomes` non-GTDB genomes in GlobDB r226 (the wishlist's Z3 rule: 1,000 gives its upper count, 500 its lower).
Order: archaea, the unclassified domain, the phyla with no GTDB genome (by name), then the minor phyla by their count
and name; inside a phylum by (taxonomy, ID). The cut keeps a phylum whole: the phylum whose genomes start at `b` of the
`B` selected bases goes to stage `1 + n_stages * b // B` (card C-41's genus rule at phylum grain); a stage over
`max_stage_bp`, or empty, is refused. Every genome carries its r232 key from card C-74's `rekey` output
(`kept`, `replaced_by <id>`, `dropped`, with the reason); a genome is never left out for it (the orchestrating session's
choice, STATUS W8). Also written: the unbuilt genomes not selected (`remainder.tsv`) and every archaeon of r226, GTDB
and not (`archaea.tsv`), the references of cards C-70 and C-78.

`stream` reads the archive once in this one process (Python's zlib, as card C-41's pass did), takes members by name,
checks records and bases against the tables, writes one FASTA per stage and one of the archaeal references, and stops
cleanly (no file left) if its own CPU time passes `cpu_guard_min` (AGENTS 7.2). A selected member the archive lacks is
reported, not invented. `screen` runs card C-63's code screen (`c63_globdb_art_phyla.screen_genetic_code`: codes 11, 25,
4; longest contigs up to 500 kb first, the whole genome when triggered, as for C-63's GART1) in worker processes that are
replaced every `max_tasks_per_child` genomes. `finish` writes each stage's `genomes.json`, `objects.json`,
`objects_staging.json`, `records.tsv`, `genetic_codes.tsv`, `repeated_records.tsv`, `shared_records.tsv` (C-63's
forms, with the r232 key), resolves shared-contig code conflicts as C-63 does, takes exact `contig_uid` ownership from
the root's stages, the `r26.11.4` stages in a26's chain order and the earlier GRARE stages, and runs
`s2_check_inputs.py` on `objects_staging.json`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import multiprocessing
import re
import resource
import shutil
import subprocess
import sys
import tarfile
import time
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c63_globdb_art_phyla as c63  # noqa: E402
import make_inputs  # noqa: E402
import r232_index  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

STAGES = ("GRARE1", "GRARE2", "GRARE3")
STAGE_FASTA = {"GRARE1": "grare1.fna.gz", "GRARE2": "grare2.fna.gz", "GRARE3": "grare3.fna.gz"}
ARCHAEA_FASTA = "globdb_r226_archaea.fna.gz"
GROUPS = ("archaea", "unclassified_domain", "no_gtdb_phylum", "minor_phylum")
GTDB_ID = re.compile(r"GC[AF]_\d{9}")
ROOT_STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
RETRIEVED_AT = "2026-10-23T19:51:00Z"  # the archive on /stor: fetched 2026-10-23, mtime 19:51 (card C-63's value)

SEL_HEAD = ["order", "genome", "stage", "group", "phylum", "phylum_non_gtdb_genomes", "c41_stage", "resource",
            "genome_size", "contig_number", "checkm2_completeness", "checkm2_contamination", "assembly_quality",
            "r232_status", "r232_id", "r232_reason", "taxonomy"]
REM_HEAD = ["genome", "c41_stage", "resource", "phylum", "phylum_non_gtdb_genomes", "genome_size", "contig_number",
            "r232_status", "r232_id", "r232_reason", "taxonomy"]
ARCH_HEAD = ["genome", "gtdb", "genome_size", "contig_number", "taxonomy"]
STREAM_REC_HEAD = ["genome", "record_ordinal", "header", "contig_uid", "length"]
ARCH_LIST_HEAD = ["genome", "gtdb", "first_record", "records", "bp", "headers_start_with_genome", "taxonomy"]
SCREEN_KEYS = ["genetic_code", "reason", "tested_alt", "cd_11", "mean_orf_11", "n_genes_11", "cd_25", "mean_orf_25",
               "n_genes_25", "cd_4", "mean_orf_4", "n_genes_4", "tga_per_kb"]
R232_COLS = ["group", "r232_status", "r232_id", "r232_reason"]
RECORDS_HEAD = [*c63.RECORDS_HEAD, *R232_COLS]
GC_HEAD = [*c63.GC_HEAD, *R232_COLS]
META_FILES = ("genomes.json", "objects.json", "objects_staging.json", "records.tsv", "repeated_records.tsv",
              "shared_records.tsv", "genetic_codes.tsv")


def _read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_rows(path: Path, head: list[str], rows) -> None:
    with open(path, "w", encoding="utf-8", buffering=1 << 20) as fh:
        fh.write("\t".join(head) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v) for v in r) + "\n")


def _base_status(status: str) -> str:
    s = status.split(" ", 1)[0]
    return "replaced" if s == "replaced_by" else s


def _need_cols(reader: csv.DictReader, need: set[str], path: Path) -> None:
    missing = sorted(need - set(reader.fieldnames or ()))
    if missing:
        raise Refused(f"{path}: missing column {missing[0]!r}")


# ---- select --------------------------------------------------------------------------------------------------------

def _read_rekey(path: Path) -> dict[str, tuple[str, str, str]]:
    out: dict[str, tuple[str, str, str]] = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        _need_cols(rd, {"r226_id", "status", "r232_id", "reason"}, path)
        for r in rd:
            st, rid = r["status"], r["r232_id"] or ""
            if st == "replaced":
                if not rid:
                    raise Refused(f"{path}: {r['r226_id']} replaced without an r232 id")
                st = f"replaced_by {rid}"
            elif st not in ("kept", "dropped"):
                raise Refused(f"{path}: {r['r226_id']} has status {st!r} (not a GlobDB r226 genome?)")
            out[r["r226_id"]] = (st, rid, r["reason"])
    return out


def select(ng_tsv, stats_gz, rekey_tsv, out_dir, *, minor_max_genomes: int = 1000, n_stages: int = 3,
           max_stage_bp: int = 12_000_000_000) -> dict:
    """The selection, the stage cut, the remainder and the archaeal references (module docstring)."""
    ng_tsv, stats_gz, rekey_tsv, out_dir = Path(ng_tsv), Path(stats_gz), Path(rekey_tsv), Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir} exists")
    if n_stages < 1 or n_stages > len(STAGES):
        raise Refused(f"n_stages {n_stages}: expected 1 to {len(STAGES)}")
    gtdb_phyla: set[str] = set()
    archaea: list[tuple] = []
    with gzip.open(stats_gz, "rt", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        _need_cols(rd, {"ID", "domain", "phylum", "class", "order", "family", "genus", "species", "genome_size",
                        "contig_number"}, stats_gz)
        for r in rd:
            is_gtdb = bool(GTDB_ID.fullmatch(r["ID"]))
            if is_gtdb:
                gtdb_phyla.add(r["phylum"])
            if r["domain"] == "Archaea":
                lineage = ";".join(f"{k}__{r[c]}" for k, c in (("d", "domain"), ("p", "phylum"), ("c", "class"),
                                                                ("o", "order"), ("f", "family"), ("g", "genus"),
                                                                ("s", "species")))
                archaea.append((r["ID"], "true" if is_gtdb else "false", int(r["genome_size"]),
                                int(r["contig_number"]), lineage))
    archaea.sort(key=lambda a: a[0])

    rows: list[dict] = []
    count: Counter = Counter()
    with open(ng_tsv, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        _need_cols(rd, {"ID", "stage", "resource", "genome_size", "contig_number", "checkm2_completeness",
                        "checkm2_contamination", "taxonomy"}, ng_tsv)
        seen: set[str] = set()
        for r in rd:
            if r["ID"] in seen:
                raise Refused(f"{ng_tsv}: duplicate ID {r['ID']!r}")
            seen.add(r["ID"])
            ranks = c63.parse_tax_ranks(r["taxonomy"])
            r["_d"], r["_p"] = ranks.get("d", ""), ranks.get("p", "")
            count[(r["_d"], r["_p"])] += 1
            rows.append(r)
    unbuilt = set(r232_index.unbuilt_ids(ng_tsv))
    rekey = _read_rekey(rekey_tsv)

    picked: list[tuple] = []  # (block key, taxonomy, id, row, group)
    remainder: list[dict] = []
    for r in rows:
        if r["ID"] not in unbuilt:
            continue
        d, p = r["_d"], r["_p"]
        if d == "Archaea":
            g = 0
        elif d == "Unclassified":
            g = 1
        elif p not in gtdb_phyla:
            g = 2
        elif count[(d, p)] < minor_max_genomes:
            g = 3
        else:
            remainder.append(r)
            continue
        if r["ID"] not in rekey:
            raise Refused(f"{rekey_tsv}: no r232 key for selected genome {r['ID']}")
        block = (g, count[(d, p)] if g == 3 else 0, d, p)
        picked.append((block, r["taxonomy"], r["ID"], r, GROUPS[g]))
    if not picked:
        raise Refused("nothing selected")
    picked.sort(key=lambda t: (t[0], t[1], t[2]))
    total = sum(int(t[3]["genome_size"]) for t in picked)
    stage_of_block: dict[tuple, str] = {}
    before = 0
    for block, _, _, r, _ in picked:
        if block not in stage_of_block:
            stage_of_block[block] = STAGES[min(n_stages, 1 + (n_stages * before) // total) - 1]
        before += int(r["genome_size"])

    sel_rows = []
    per_stage: dict[str, dict] = {s: {"genomes": 0, "bp": 0, "contigs": 0, "phyla": set(), "by_group": Counter(),
                                      "r232_status": Counter()} for s in STAGES[:n_stages]}
    for i, (block, t, gid, r, group) in enumerate(picked):
        st = stage_of_block[block]
        status, rid, reason = rekey[gid]
        q = c63.mag_quality(float(r["checkm2_completeness"]), float(r["checkm2_contamination"]))
        sel_rows.append((i, gid, st, group, r["_p"], count[(r["_d"], r["_p"])], r["stage"], r["resource"],
                         r["genome_size"], r["contig_number"], r["checkm2_completeness"], r["checkm2_contamination"],
                         q, status, rid, reason, t))
        ps = per_stage[st]
        ps["genomes"] += 1
        ps["bp"] += int(r["genome_size"])
        ps["contigs"] += int(r["contig_number"])
        ps["phyla"].add((r["_d"], r["_p"]))
        ps["by_group"][group] += 1
        ps["r232_status"][_base_status(status)] += 1
    for st, ps in per_stage.items():
        if ps["genomes"] == 0:
            raise Refused(f"{st} is empty")
        if ps["bp"] > max_stage_bp:
            raise Refused(f"{st} holds {ps['bp']} bp, more than {max_stage_bp}")

    remainder.sort(key=lambda r: (r["taxonomy"], r["ID"]))
    rem_rows = []
    rem_status: Counter = Counter()
    for r in remainder:
        status, rid, reason = rekey.get(r["ID"], ("", "", "no_rekey_row"))
        rem_status[_base_status(status) or "no_rekey_row"] += 1
        rem_rows.append((r["ID"], r["stage"], r["resource"], r["_p"], count[(r["_d"], r["_p"])], r["genome_size"],
                         r["contig_number"], status, rid, reason, r["taxonomy"]))

    tmp = out_dir.with_name(f".{out_dir.name}.tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    _write_rows(tmp / "selection.tsv", SEL_HEAD, sel_rows)
    _write_rows(tmp / "remainder.tsv", REM_HEAD, rem_rows)
    _write_rows(tmp / "archaea.tsv", ARCH_HEAD, archaea)
    report = {
        "parameters": {"minor_max_genomes": minor_max_genomes, "n_stages": n_stages, "max_stage_bp": max_stage_bp,
                       "keep_whole": "phylum",
                       "rule": "unbuilt (not G1, not a GART1/GART2 phylum) and (Archaea, or domain Unclassified, or "
                               "a phylum with no GTDB genome in GlobDB r226, or a bacterial phylum with fewer than "
                               "minor_max_genomes non-GTDB genomes in GlobDB r226)"},
        "inputs": {p.name: ids.sha256_file(p).hex() for p in (ng_tsv, stats_gz, rekey_tsv)},
        "unbuilt": len(unbuilt),
        "selected": {"genomes": len(sel_rows), "bp": total,
                     "contigs": sum(int(t[3]["contig_number"]) for t in picked),
                     "phyla": len({(t[3]["_d"], t[3]["_p"]) for t in picked}),
                     "by_group": dict(sorted(Counter(t[4] for t in picked).items()))},
        "stages": {st: {**{k: v for k, v in ps.items() if k not in ("phyla", "by_group", "r232_status")},
                        "phyla": len(ps["phyla"]), "by_group": dict(sorted(ps["by_group"].items())),
                        "r232_status": dict(sorted(ps["r232_status"].items()))} for st, ps in per_stage.items()},
        "remainder": {"genomes": len(rem_rows), "bp": sum(int(r["genome_size"]) for r in remainder),
                      "contigs": sum(int(r["contig_number"]) for r in remainder),
                      "phyla": len({(r["_d"], r["_p"]) for r in remainder}),
                      "r232_status": dict(sorted(rem_status.items()))},
        "archaea_references": {"genomes": len(archaea), "gtdb": sum(1 for a in archaea if a[1] == "true"),
                               "bp": sum(a[2] for a in archaea), "contigs": sum(a[3] for a in archaea)},
    }
    report["sha256"] = {n: ids.sha256_file(tmp / n).hex() for n in ("selection.tsv", "remainder.tsv", "archaea.tsv")}
    (tmp / "selection.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    tmp.rename(out_dir)
    return report


# ---- stream --------------------------------------------------------------------------------------------------------

def _cpu_min_self() -> float:
    ru = resource.getrusage(resource.RUSAGE_SELF)
    return (ru.ru_utime + ru.ru_stime) / 60.0


def _cpu_min_children() -> float:
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    return (ru.ru_utime + ru.ru_stime) / 60.0


def _stages_in(sel: list[dict]) -> list[str]:
    return [s for s in STAGES if any(r["stage"] == s for r in sel)]


def stream(select_dir, tar_gz, meta_out, fasta_out, *, use_bgzip: bool = False, cpu_guard_min: float = 280.0,
           log_every_s: float = 60.0) -> dict:
    """One pass over the archive (module docstring); returns the report also written as stream/stream.json."""
    select_dir, tar_gz, meta_out, fasta_out = Path(select_dir), Path(tar_gz), Path(meta_out), Path(fasta_out)
    sel = _read_tsv(select_dir / "selection.tsv")
    arch = _read_tsv(select_dir / "archaea.tsv")
    stages = _stages_in(sel)
    want = {r["genome"]: (r["stage"], int(r["genome_size"]), int(r["contig_number"])) for r in sel}
    want_a = {r["genome"]: (r["gtdb"], int(r["genome_size"]), int(r["contig_number"]), r["taxonomy"]) for r in arch}
    outs = {**{st: STAGE_FASTA[st] for st in stages}, "archaea": ARCHAEA_FASTA}
    sdir, adir = meta_out / "stream", meta_out / "archaea"
    for p in (sdir, adir, *(fasta_out / f for f in outs.values())):
        if p.exists():
            raise Refused(f"{p} exists")
    meta_out.mkdir(parents=True, exist_ok=True)
    fasta_out.mkdir(parents=True, exist_ok=True)
    stmp, atmp = meta_out / ".stream.tmp", meta_out / ".archaea.tmp"
    for d in (stmp, atmp):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir()
    partial = {k: fasta_out / f".{f}.partial" for k, f in outs.items()}
    writers: dict[str, c63._GzipWriter] = {}
    recf = {}
    found: dict[str, list[str]] = {k: [] for k in outs}
    counts = {k: {"records": 0, "bp": 0, "fasta_bytes": 0} for k in outs}
    arch_rows: list[tuple] = []
    header_mismatch = 0
    members = other = 0
    t0 = time.monotonic()
    last = t0
    try:
        for k in outs:
            partial[k].unlink(missing_ok=True)
            writers[k] = c63._GzipWriter(partial[k], use_bgzip=use_bgzip)
        for st in stages:
            recf[st] = open(stmp / f"{st}.records.tsv", "w", encoding="utf-8", buffering=1 << 20)
            recf[st].write("\t".join(STREAM_REC_HEAD) + "\n")
        with tarfile.open(tar_gz, mode="r|gz", bufsize=1 << 20) as tf:
            for m in tf:
                if not m.isfile():
                    continue
                name = m.name.rsplit("/", 1)[-1]
                if not name.endswith(".fa.gz"):
                    other += 1
                    continue
                members += 1
                used = _cpu_min_self()
                if used > cpu_guard_min:
                    raise Refused(f"CPU guard: this process used {used:.1f} CPU-minutes, more than {cpu_guard_min} "
                                  f"(at member {members}, {m.name})")
                gid = name[:-6]
                s = want.get(gid)
                a = want_a.get(gid)
                if s is None and a is None:
                    continue
                if (s is not None and gid in found[s[0]]) or (a is not None and gid in found["archaea"]):
                    raise Refused(f"{gid}: two members of the archive")
                ef = tf.extractfile(m)
                assert ef is not None
                recs = c63._iter_fasta_bytes(gzip.decompress(ef.read()))
                n, bp = len(recs), sum(len(q) for _, q, _ in recs)
                for want_n, want_bp, what in ((s[2], s[1], "selection") if s else (None, None, ""),
                                              (a[2], a[1], "archaea table") if a else (None, None, "")):
                    if want_n is not None and (n, bp) != (want_n, want_bp):
                        raise Refused(f"{gid}: the archive has {n} records / {bp} bp, the {what} says "
                                      f"{want_n} / {want_bp}")
                chunks = []
                for h, q, lw in recs:
                    if "\t" in h:
                        raise Refused(f"{gid}: a header with a tab")
                    chunks.append(b"\n".join(fasta.format_record(h, q, lw)) + b"\n")
                data = b"".join(chunks)
                if s is not None:
                    st = s[0]
                    writers[st].write(data)
                    for i, (h, q, _) in enumerate(recs):
                        recf[st].write(f"{gid}\t{i}\t{h}\t{ids.contig_uid(q).hex()}\t{len(q)}\n")
                    found[st].append(gid)
                    counts[st]["records"] += n
                    counts[st]["bp"] += bp
                    counts[st]["fasta_bytes"] += len(data)
                if a is not None:
                    writers["archaea"].write(data)
                    ok = all(h.startswith(gid) for h, _, _ in recs)
                    header_mismatch += 0 if ok else 1
                    arch_rows.append((gid, a[0], counts["archaea"]["records"], n, bp, str(ok).lower(), a[3]))
                    found["archaea"].append(gid)
                    counts["archaea"]["records"] += n
                    counts["archaea"]["bp"] += bp
                    counts["archaea"]["fasta_bytes"] += len(data)
                now = time.monotonic()
                if now - last >= log_every_s:
                    print(f"[c75] members={members} " + " ".join(f"{k}={len(v)}" for k, v in found.items())
                          + f" cpu_min={_cpu_min_self():.1f} elapsed={now - t0:.0f}s", file=sys.stderr, flush=True)
                    last = now
        for f in recf.values():
            f.close()
        for w in writers.values():
            w.close()
    except BaseException:
        for f in recf.values():
            f.close()
        for k, w in writers.items():
            try:
                w.close()
            except Exception:
                pass
            partial[k].unlink(missing_ok=True)
        shutil.rmtree(stmp, ignore_errors=True)
        shutil.rmtree(atmp, ignore_errors=True)
        raise
    missing = {}
    for st in stages:
        miss = sorted({g for g, v in want.items() if v[0] == st} - set(found[st]))
        if miss:
            missing[st] = miss
    miss_a = sorted(set(want_a) - set(found["archaea"]))
    if miss_a:
        missing["archaea"] = miss_a
    _write_rows(stmp / "missing.tsv", ["set", "genome"], [(k, g) for k, v in missing.items() for g in v])
    _write_rows(atmp / "archaea_list.tsv", ARCH_LIST_HEAD, arch_rows)
    for k, f in outs.items():
        partial[k].rename(fasta_out / f)
    rep = {
        "streams": 1, "archive": str(tar_gz), "archive_bytes": tar_gz.stat().st_size, "members_seen": members,
        "other_files": other, "missing": missing,
        "sets": {k: {"genomes": len(found[k]), **counts[k], "file": str(fasta_out / outs[k]),
                     "bytes": (fasta_out / outs[k]).stat().st_size,
                     "sha256": ids.sha256_file(fasta_out / outs[k]).hex()} for k in outs},
        "archaea_headers_not_starting_with_genome": header_mismatch,
        "archaea_gtdb_genomes": sum(1 for r in arch_rows if r[1] == "true"),
        "cpu_min_self": round(_cpu_min_self(), 3), "cpu_min_children": round(_cpu_min_children(), 3),
        "wall_s": round(time.monotonic() - t0, 1), "cpu_guard_min": cpu_guard_min,
    }
    (stmp / "stream.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    stmp.rename(sdir)
    atmp.rename(adir)
    return rep


# ---- screen --------------------------------------------------------------------------------------------------------

def _genomes_of(st: str, meta_out: Path, fasta_out: Path):
    """(genome, [sequences]) in FASTA order, each record checked against the stream's record list."""
    rows = _read_tsv(meta_out / "stream" / f"{st}.records.tsv")
    it = iter(rows)
    cur, seqs = None, []
    for rec in fasta.iter_records(fasta_out / STAGE_FASTA[st]):
        r = next(it, None)
        if r is None or r["header"] != rec.header or int(r["length"]) != len(rec.seq):
            raise Refused(f"{st}: FASTA record {rec.header!r} does not match the stream's record list")
        if r["genome"] != cur:
            if cur is not None:
                yield cur, seqs
            cur, seqs = r["genome"], []
        seqs.append(rec.seq)
    if next(it, None) is not None:
        raise Refused(f"{st}: the FASTA has fewer records than the stream's record list")
    if cur is not None:
        yield cur, seqs


def screen(select_dir, meta_out, fasta_out, *, workers: int = 1, screen_max_bp: int = 500_000,
           max_tasks_per_child: int = 200) -> dict:
    """C-63's code screen of every streamed stage genome; writes screen/<stage>.tsv sorted by genome."""
    select_dir, meta_out, fasta_out = Path(select_dir), Path(meta_out), Path(fasta_out)
    sel = _read_tsv(select_dir / "selection.tsv")
    taxonomy = {r["genome"]: r["taxonomy"] for r in sel}
    out = meta_out / "screen"
    if out.exists():
        raise Refused(f"{out} exists")
    tmp = meta_out / ".screen.tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    t0 = time.monotonic()
    rep: dict = {"workers": workers, "screen_max_bp": screen_max_bp, "max_tasks_per_child": max_tasks_per_child,
                 "stages": {}}
    pool = None
    if workers > 1:
        pool = ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                                   max_tasks_per_child=max_tasks_per_child)
    try:
        for st in _stages_in(sel):
            res: dict[str, dict] = {}
            if pool is None:
                for gid, seqs in _genomes_of(st, meta_out, fasta_out):
                    res[gid] = c63.screen_genetic_code(seqs, taxonomy[gid], screen_max_bp=screen_max_bp)
            else:
                pending: set = set()
                for gid, seqs in _genomes_of(st, meta_out, fasta_out):
                    pending.add(pool.submit(c63._screen_worker, (gid, seqs, taxonomy[gid], screen_max_bp)))
                    if len(pending) >= workers * 4:
                        done, pending = wait(pending, return_when=FIRST_COMPLETED)
                        for f in done:
                            g, r = f.result()
                            res[g] = r
                for f in pending:
                    g, r = f.result()
                    res[g] = r
            _write_rows(tmp / f"{st}.tsv", ["genome", *SCREEN_KEYS],
                        [(g, *(("true" if res[g][k] else "false") if k == "tested_alt" else res[g][k]
                               for k in SCREEN_KEYS)) for g in sorted(res)])
            rep["stages"][st] = {"genomes": len(res), "tested_alt": sum(1 for r in res.values() if r["tested_alt"]),
                                 "codes": dict(sorted(Counter(str(r["genetic_code"]) for r in res.values()).items())),
                                 "reasons": dict(sorted(Counter(r["reason"] for r in res.values()).items()))}
    except BaseException:
        if pool is not None:
            pool.shutdown(wait=False, cancel_futures=True)
            pool = None
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True)
    rep.update(cpu_min_self=round(_cpu_min_self(), 3), cpu_min_children=round(_cpu_min_children(), 3),
               wall_s=round(time.monotonic() - t0, 1))
    (tmp / "screen.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    tmp.rename(out)
    return rep


# ---- finish --------------------------------------------------------------------------------------------------------

def _check(stage: str, d: Path) -> dict:
    p = subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "s2_check_inputs.py"), stage,
                        str(d / "genomes.json"), str(d / "objects_staging.json")],
                       capture_output=True, text=True, check=False)
    doc = json.loads(p.stdout) if p.stdout.strip().startswith("{") else {"stdout": p.stdout}
    return {"exit": p.returncode, "problems": doc.get("problems"), "genomes": doc.get("genomes"),
            "genome_input_rows": doc.get("genome_input_rows"), "codes": doc.get("codes"),
            "stderr": p.stderr[-2000:]}


def _fmt(v) -> str:
    return "" if v is None or v == "" else str(v)


def finish(select_dir, meta_out, fasta_out, *, inputs_home=Path("/tmp/if2453/s2/inputs"), root=None,
           root_stages=ROOT_STAGES, chain: dict[str, Path] | None = None, retrieved_at: str = RETRIEVED_AT,
           work=None, threads: int = 2, memory_limit: str = "4GB", check: bool = True) -> dict:
    """Each stage's files in /tmp form (module docstring); returns the report also written as finish.json."""
    select_dir, meta_out, fasta_out = Path(select_dir), Path(meta_out), Path(fasta_out)
    inputs_home = Path(inputs_home)
    sel = _read_tsv(select_dir / "selection.tsv")
    stages = _stages_in(sel)
    for st in stages:
        if (meta_out / st).exists():
            raise Refused(f"{meta_out / st} exists")
    owners_in = dict(chain or {})
    for name, d in owners_in.items():
        if not (Path(d) / "genomes.json").is_file():
            raise Refused(f"chain stage {name}: no genomes.json in {d}")
    stream_rep = json.loads((meta_out / "stream" / "stream.json").read_text())
    missing = stream_rep["missing"]
    work_dir = Path(work) if work is not None else meta_out / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    report: dict = {"owner_stages": [*(list(root_stages) if root is not None else []), *owners_in]}
    for st in stages:
        miss = set(missing.get(st, []))
        rows = [r for r in sel if r["stage"] == st and r["genome"] not in miss]
        contigs: dict[str, list[tuple[bytes, int]]] = defaultdict(list)
        for r in _read_tsv(meta_out / "stream" / f"{st}.records.tsv"):
            contigs[r["genome"]].append((bytes.fromhex(r["contig_uid"]), int(r["length"])))
        gc = {r["genome"]: r for r in _read_tsv(meta_out / "screen" / f"{st}.tsv")}
        for r in rows:
            if r["genome"] not in contigs or r["genome"] not in gc:
                raise Refused(f"{st}: {r['genome']} has no streamed records or no screen result")
        tmp = meta_out / f".{st}.tmp"
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True)

        # DEC-27: genomes sharing a contig agree on the code; the least genome_uid's code wins (C-63's rule)
        guid = {r["genome"]: ids.genome_uid(list(dict.fromkeys(u for u, _ in contigs[r["genome"]]))) for r in rows}
        code = {r["genome"]: int(gc[r["genome"]]["genetic_code"]) for r in rows}
        first_code: dict[bytes, tuple[int, str]] = {}
        conflicts = []
        for r in sorted(rows, key=lambda x: guid[x["genome"]]):
            g = r["genome"]
            target = None
            for u, _ in contigs[g]:
                prev = first_code.get(u)
                if prev is not None and prev[0] != code[g]:
                    target = prev[0]
                    conflicts.append({"contig_uid": u.hex(), "first_genome": prev[1], "first_code": prev[0],
                                      "second_genome": g, "second_code": code[g]})
            if target is not None:
                code[g] = target
            for u, _ in contigs[g]:
                first_code.setdefault(u, (code[g], g))

        ordinal = 0
        with open(tmp / "records.tsv", "w", encoding="utf-8", buffering=1 << 20) as rf, \
                open(tmp / "genetic_codes.tsv", "w", encoding="utf-8", buffering=1 << 20) as gf:
            rf.write("\t".join(RECORDS_HEAD) + "\n")
            gf.write("\t".join(GC_HEAD) + "\n")
            for r in rows:
                g = r["genome"]
                ranks = c63.parse_tax_ranks(r["taxonomy"])
                s = gc[g]
                key = [r["group"], r["r232_status"], r["r232_id"], r["r232_reason"]]
                gf.write("\t".join([g, st, r["phylum"], ranks.get("c", ""), ranks.get("o", ""), r["resource"],
                                    r["genome_size"], r["contig_number"], r["assembly_quality"], str(code[g]),
                                    s["reason"], s["tested_alt"], *(_fmt(s[k]) for k in SCREEN_KEYS[3:]), *key])
                         + "\n")
                for i, (u, n) in enumerate(contigs[g]):
                    rf.write("\t".join([str(ordinal), g, str(i), u.hex(), str(n), r["phylum"], r["resource"],
                                        r["assembly_quality"], str(code[g]), r["c41_stage"], *key]) + "\n")
                    ordinal += 1

        owners = c63._owners_from_records(tmp / "records.tsv", None if root is None else Path(root),
                                          list(root_stages), owners_in, work=work_dir, threads=threads,
                                          memory_limit=memory_limit)
        listing = make_inputs.Listing()
        genomes, repeated, shared = [], [], []
        first_of: dict[bytes, str] = {}
        contig_len: dict[bytes, int] = {}
        ordinal = 0
        for r in rows:
            g = r["genome"]
            uids = [u for u, _ in contigs[g]]
            for i, (u, n) in enumerate(contigs[g]):
                contig_len.setdefault(u, n)
                if u in owners:
                    shared.append((ordinal, g, i, u.hex(), n, owners[u]))
                ordinal += 1
            entry = listing.add(uids, code[g], "prokaryote", r["assembly_quality"])
            if entry is None:
                repeated.append((ordinal - len(uids), g, guid[g].hex(), first_of[guid[g]]))
                continue
            genomes.append(entry)
            first_of[guid[g]] = g
        make_inputs._write_tsv(tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"],
                               repeated)
        make_inputs._write_tsv(tmp / "shared_records.tsv",
                               ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"], shared)
        src = {**c63.GLOBDB_SOURCE_BASE, "retrieved_at": retrieved_at}
        (tmp / "objects.json").write_text(json.dumps([{"path": str(inputs_home / st / STAGE_FASTA[st]), **src}],
                                                     indent=1))
        (tmp / "objects_staging.json").write_text(json.dumps([{"path": str(fasta_out / STAGE_FASTA[st]), **src}],
                                                             indent=1))
        (tmp / "genomes.json").write_text(json.dumps(genomes))
        tmp.rename(meta_out / st)
        owners_in[st] = meta_out / st

        own_c: Counter = Counter()
        own_b: Counter = Counter()
        for u, n in contig_len.items():
            if u in owners:
                own_c[owners[u]] += 1
                own_b[owners[u]] += n
        new_c = sum(1 for u in contig_len if u not in owners)
        new_b = sum(n for u, n in contig_len.items() if u not in owners)
        files = {n: meta_out / st / n for n in META_FILES}
        files[STAGE_FASTA[st]] = fasta_out / STAGE_FASTA[st]
        rep = {
            "selected_genomes": len(rows) + len(miss), "missing_genomes": sorted(miss), "genomes": len(genomes),
            "repeated_genomes": len(repeated), "records": ordinal, "contigs": len(contig_len),
            "contig_bases": sum(contig_len.values()), "shared_records": len(shared),
            "owner_contigs": dict(sorted(own_c.items())), "owner_bases": dict(sorted(own_b.items())),
            "new_contigs": new_c, "new_bases": new_b,
            "codes": dict(sorted(Counter(str(e["genetic_code"]) for e in genomes).items())),
            "tested_alt_genomes": sum(1 for r in rows if gc[r["genome"]]["tested_alt"] == "true"),
            "code_conflicts": conflicts,
            "quality": dict(sorted(Counter(e["assembly_quality"] for e in genomes).items())),
            "by_group": dict(sorted(Counter(r["group"] for r in rows).items())),
            "phyla": len({r["phylum"] for r in rows}),
            "r232_status": dict(sorted(Counter(_base_status(r["r232_status"]) for r in rows).items())),
            "estimates": c63._estimate_called_crispr(stream_rep["sets"][st]["fasta_bytes"], new_c, new_b,
                                                     sum(own_b.values())),
            "sha256": {n: ids.sha256_file(p).hex() for n, p in files.items()},
            "bytes": {n: p.stat().st_size for n, p in files.items()},
            **listing.report(),
        }
        if check:
            rep["s2_check"] = _check(st, meta_out / st)
        report[st] = rep
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    (meta_out / "finish.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    return report


# ---- command line ----------------------------------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage GlobDB r226 GRARE1..3 (card C-75, DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("select")
    for n in ("non_gtdb_stages", "tax_plus_stats", "rekey", "out_dir"):
        p.add_argument(n, type=Path)
    p.add_argument("--minor-max-genomes", type=int, default=1000)
    p.add_argument("--n-stages", type=int, default=3)
    p.add_argument("--max-stage-bp", type=int, default=12_000_000_000)
    p = sub.add_parser("stream")
    for n in ("select_dir", "tar_gz", "meta_out", "fasta_out"):
        p.add_argument(n, type=Path)
    p.add_argument("--bgzip", action="store_true")
    p.add_argument("--cpu-guard-min", type=float, default=280.0)
    p = sub.add_parser("screen")
    for n in ("select_dir", "meta_out", "fasta_out"):
        p.add_argument(n, type=Path)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--screen-max-bp", type=int, default=500_000)
    p.add_argument("--max-tasks-per-child", type=int, default=200)
    p = sub.add_parser("finish")
    for n in ("select_dir", "meta_out", "fasta_out"):
        p.add_argument(n, type=Path)
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--root-stages", default=",".join(ROOT_STAGES))
    p.add_argument("--chain-home", type=Path, default=None,
                   help="the r26.11.4 stages' inputs, read in a26's chain order as owners")
    p.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    p.add_argument("--retrieved-at", default=RETRIEVED_AT)
    p.add_argument("--work", type=Path, default=None)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--memory-limit", default="4GB")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "select":
            rep = select(a.non_gtdb_stages, a.tax_plus_stats, a.rekey, a.out_dir,
                         minor_max_genomes=a.minor_max_genomes, n_stages=a.n_stages, max_stage_bp=a.max_stage_bp)
        elif a.cmd == "stream":
            rep = stream(a.select_dir, a.tar_gz, a.meta_out, a.fasta_out, use_bgzip=a.bgzip,
                         cpu_guard_min=a.cpu_guard_min)
        elif a.cmd == "screen":
            rep = screen(a.select_dir, a.meta_out, a.fasta_out, workers=a.workers, screen_max_bp=a.screen_max_bp,
                         max_tasks_per_child=a.max_tasks_per_child)
        else:
            chain = {s: a.chain_home / s for s, _ in a26.R26_11_4_CHAIN} if a.chain_home is not None else {}
            rep = finish(a.select_dir, a.meta_out, a.fasta_out, inputs_home=a.inputs_home, root=a.root,
                         root_stages=[s for s in a.root_stages.split(",") if s], chain=chain,
                         retrieved_at=a.retrieved_at, work=a.work, threads=a.threads, memory_limit=a.memory_limit)
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(rep, indent=1, sort_keys=True))
    if a.cmd == "finish":
        return 0 if all(rep[s]["s2_check"]["exit"] == 0 and rep[s]["s2_check"]["problems"] == []
                        for s in STAGES if s in rep) else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

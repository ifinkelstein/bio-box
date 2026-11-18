"""Card C-40: what stage PL2's inputs share with the stored stages and between their sources, read from the inputs and
the data root (read only).

    python c40_analyze.py --root ROOT --inputs IN --scope-plsdb FASTA --scope-dir DL_PLASMIDSCOPE --gut-table TSV
                          --out OUT

IN holds what `make_inputs.py plasmid-sources ... --root ROOT` wrote (records.tsv, skipped_records.tsv, genomes.json,
objects.json and the objects). FASTA is PlasmidScope's PLSDB archive as one FASTA (c40_concat.py; not a source of PL2).
DL_PLASMIDSCOPE holds PlasmidScope's tables (`<SRC>.plasmid_list.download.tsv`, `ALL...`). TSV is the gut table
decompressed. The identity everywhere is `ids.contig_uid`. Nothing in ROOT is written. Writes OUT/summary.json:

- per source: records read, taken and skipped; genomes it gives (the first record of a sequence); the records, contigs
  and bases whose contig each stored stage owns (records.tsv `owner_stage`: the stage whose contig_index holds the
  default variant, DEC-26); the records and contigs that stage PL1 lists in its genome_input (PLSDB already has them,
  as a genome of PL1, whoever owns the contig); the contigs and bases the stage would add;
- between sources: distinct contigs each pair shares (records taken);
- PlasmidScope's own PLSDB file against PL1: records, distinct contigs, those PL1 lists, and PL1's contigs it lacks;
- lengths of the genomes PL2 would store (first records of contigs no stage owns), by source: quantiles, counts and
  bases in length bins, the share under 1.5 and 5 kbp (card C-38 (4)), the share complete and circular;
- terminal repeats: for records whose table says direct terminal repeat and for the circular gut plasmids, the
  longest exact repeat of the contig's start at its end (searched in the last 5,000 bases), which says whether the
  stored bytes hold the repeat twice;
- characters outside ACGTN (any case) in the stored contigs;
- PlasmidScope's tables: columns, Host values filled, and the ALL table's Data_Source lists against PL2's sources;
- the bytes each stage's files hold in ROOT beside the bases its contig_index stores (storage per Gbp).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import duckdb

from biobox import fasta, ids, io
from biobox.l0seq import store

BINS = (1_500, 5_000, 10_000, 50_000, 100_000, 500_000)
csv.field_size_limit(1 << 30)


def tsv(path: Path):
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def lengths_summary(lens: list[int], circ: list[bool], comp: list[bool]) -> dict:
    if not lens:
        return {"n": 0}
    s = sorted(lens)
    q = lambda f: s[min(len(s) - 1, int(f * len(s)))]  # noqa: E731
    bins, edges = Counter(), Counter()
    for n in lens:
        k = next((i for i, b in enumerate(BINS) if n < b), len(BINS))
        bins[k] += 1
        edges[k] += n
    names = [f"<{BINS[0]}"] + [f"{BINS[i - 1]}-{BINS[i]}" for i in range(1, len(BINS))] + [f">={BINS[-1]}"]
    tot = sum(lens)
    return {"n": len(lens), "bases": tot, "min": s[0], "p10": q(0.10), "p25": q(0.25), "median": q(0.5),
            "p75": q(0.75), "p90": q(0.90), "max": s[-1], "mean": round(statistics.fmean(lens), 1),
            "bins": {names[k]: {"n": bins[k], "bases": edges[k]} for k in range(len(names))},
            "under_1500": {"n": sum(1 for n in lens if n < 1500), "bases": sum(n for n in lens if n < 1500)},
            "under_5000": {"n": sum(1 for n in lens if n < 5000), "bases": sum(n for n in lens if n < 5000)},
            "circular": sum(circ), "complete": sum(comp), "share_circular": round(sum(circ) / len(lens), 4),
            "share_complete": round(sum(comp) / len(lens), 4),
            "share_complete_bases": round(sum(n for n, c in zip(lens, comp) if c) / tot, 4)}


def terminal_repeat(seq: bytes, window: int = 5000, seed: int = 20) -> int:
    """The longest k (at least `seed`) with seq[:k] == seq[-k:], k at most `window`; 0 when none."""
    if len(seq) < 2 * seed:
        return 0
    s = seq.upper()
    tail = s[-min(window, len(s) // 2):]
    w = len(tail)
    best, p = 0, tail.find(s[:seed])
    while p != -1:
        if tail[p:] == s[:w - p]:
            best = max(best, w - p)
            break  # the first hit from the left is the longest
        p = tail.find(s[:seed], p + 1)
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    for a in ("--root", "--inputs", "--scope-plsdb", "--scope-dir", "--gut-table", "--out"):
        ap.add_argument(a, required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    out: dict = {}

    recs = list(tsv(a.inputs / "records.tsv"))
    skipped = list(tsv(a.inputs / "skipped_records.tsv"))
    sources = list(dict.fromkeys(r["source"] for r in recs))
    out["skipped"] = {s: dict(Counter(r["reason"] for r in skipped if r["source"] == s)) for s in sources}

    con = duckdb.connect(config={"threads": 2, "memory_limit": "6GB", "temp_directory": str(a.out / "duckdb")})
    gi_pl1 = io.part_files(a.root / "l0seq/genome_input/stage=PL1")
    pl1 = {r[0] for r in con.execute(f"SELECT DISTINCT contig_uid FROM read_parquet({sql_list(gi_pl1)}, "
                                     "hive_partitioning = false)").fetchall()}
    out["pl1_genome_input_contigs"] = len(pl1)

    per = {}
    by_source_contigs: dict[str, set[bytes]] = defaultdict(set)
    for s in sources:
        rs = [r for r in recs if r["source"] == s]
        cus = {bytes.fromhex(r["contig_uid"]) for r in rs}
        by_source_contigs[s] = cus
        firsts = [r for r in rs if not r["first_accession"]]
        own_rec, own_ctg, own_bases = Counter(), defaultdict(set), Counter()
        for r in rs:
            if r["owner_stage"]:
                own_rec[r["owner_stage"]] += 1
                if r["contig_uid"] not in own_ctg[r["owner_stage"]]:
                    own_ctg[r["owner_stage"]].add(r["contig_uid"])
                    own_bases[r["owner_stage"]] += int(r["length"])
        in_pl1 = [r for r in rs if bytes.fromhex(r["contig_uid"]) in pl1]
        new_first = [r for r in firsts if not r["owner_stage"]]
        per[s] = {"records_taken": len(rs), "distinct_contigs": len(cus), "genomes": len(firsts),
                  "bases_taken": sum(int(r["length"]) for r in rs),
                  "owned_by_stage": {st: {"records": own_rec[st], "contigs": len(own_ctg[st]), "bases": own_bases[st]}
                                     for st in sorted(own_rec)},
                  "plsdb_has": {"records": len(in_pl1), "contigs": len({r["contig_uid"] for r in in_pl1})},
                  "adds": {"contigs": len(new_first), "bases": sum(int(r["length"]) for r in new_first)},
                  "repeats_of_earlier_source_by_first_source": dict(Counter(
                      r["first_source"] for r in rs if r["first_source"] and r["first_source"] != s)),
                  "topology_values": dict(Counter(r["topology_value"] for r in rs)),
                  "completeness_values": dict(Counter(r["completeness_value"] for r in rs)),
                  "fragments": sum(r["fragment"] == "true" for r in rs)}
    out["sources"] = per
    out["pairs_shared_contigs"] = {f"{x}|{y}": len(by_source_contigs[x] & by_source_contigs[y])
                                   for i, x in enumerate(sources) for y in sources[i + 1:]
                                   if by_source_contigs[x] & by_source_contigs[y]}
    allc = set().union(*by_source_contigs.values())
    out["all"] = {"records_taken": len(recs), "distinct_contigs": len(allc),
                  "plsdb_has_contigs": len(allc & pl1),
                  "owned_contigs": len({r["contig_uid"] for r in recs if r["owner_stage"]}),
                  "owned_by_stage_contigs": dict(Counter(st for st in {r["contig_uid"]: r["owner_stage"]
                                                                          for r in recs if r["owner_stage"]}.values()))}

    # lengths of what the stage would store: first records of contigs no stage owns
    new = [r for r in recs if not r["first_accession"] and not r["owner_stage"]]
    ln = {}
    for s in sources + ["all"]:
        rs = [r for r in new if s == "all" or r["source"] == s]
        ln[s] = lengths_summary([int(r["length"]) for r in rs], [r["circular"] == "true" for r in rs],
                                [r["assembly_quality"] == "complete" for r in rs])
    out["lengths_new"] = ln
    allf = [r for r in recs if not r["first_accession"]]
    out["lengths_all_genomes"] = lengths_summary([int(r["length"]) for r in allf], [r["circular"] == "true" for r in allf],
                                                 [r["assembly_quality"] == "complete" for r in allf])

    # terminal repeats and characters, from the objects (records the stage would store)
    want_tr = {(r["source"], r["accession"]) for r in recs
               if r["topology_value"] == "direct terminal repeat" or (r["source"].startswith("gut") and
                                                                       r["circular"] == "true")}
    new_keys = {(r["source"], r["accession"]) for r in new}
    tr = defaultdict(Counter)
    odd = Counter()
    odd_records = Counter()
    for o in json.loads((a.inputs / "objects.json").read_text()):
        p = Path(o["path"])
        src = p.stem
        for r in fasta.iter_records(p):
            key = (src, r.accession)
            if key in want_tr:
                k = terminal_repeat(r.seq)
                lab = "none" if k == 0 else "20-99" if k < 100 else "100-999" if k < 1000 else ">=1000"
                tr[f"{src}|{'dtr' if not src.startswith('gut') else 'circular'}"][lab] += 1
            if key in new_keys:
                bad = Counter(r.seq.upper().translate(None, b"ACGTN"))
                if bad:
                    odd_records[src] += 1
                    for c, n in bad.items():
                        odd[chr(c)] += n
    out["terminal_repeats"] = {k: dict(v) for k, v in tr.items()}
    out["non_acgtn"] = {"records_by_source": dict(odd_records), "characters": dict(odd.most_common(20))}

    # PlasmidScope's PLSDB file against PL1
    ps = {}
    n = 0
    for r in fasta.iter_records(a.scope_plsdb):
        n += 1
        ps.setdefault(ids.contig_uid(r.seq), 0)
    psset = set(ps)
    out["plasmidscope_plsdb_vs_pl1"] = {"records": n, "distinct_contigs": len(psset), "in_pl1": len(psset & pl1),
                                        "not_in_pl1": len(psset - pl1), "pl1_not_in_it": len(pl1 - psset),
                                        "in_pl2_sources": {s: len(psset & by_source_contigs[s]) for s in sources}}

    # PlasmidScope's tables
    tabs = {}
    for f in sorted(a.scope_dir.glob("*.plasmid_list.download.tsv")):
        rows = list(tsv(f))
        cols = list(rows[0].keys()) if rows else []
        host = Counter("filled" if (r.get("Host") or "-") not in ("-", "") else "empty" for r in rows)
        tabs[f.name.split(".")[0]] = {"rows": len(rows), "columns": cols, "host": dict(host)}
        if f.name.startswith("ALL"):
            ds = Counter()
            n_multi = 0
            for r in rows:
                srcs = sorted(set(r["Data_Source"].split(",")))
                n_multi += len(srcs) > 1
                for x in srcs:
                    ds[x] += 1
            tabs["ALL"].update({"rows_with_several_sources": n_multi, "rows_by_source": dict(ds)})
    out["plasmidscope_tables"] = tabs

    # the gut table: what the 68,350 representatives are, against all predictions
    g = Counter()
    with open(a.gut_table) as fh:
        lines = (x for x in fh if not x.startswith("#"))
        for r in csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE):
            g[(r["in_dereplicated_plasmids"], r["cluster_type"], r["circular"])] += 1
    out["gut_table_derep_type_circular"] = {"|".join(k): v for k, v in sorted(g.items())}

    # bytes per stage in the root, beside the bases each stores
    size = Counter()
    for dp, _, fs in os.walk(a.root):
        for f in fs:
            p = os.path.join(dp, f)
            for part in Path(p).relative_to(a.root).parts:
                if part.startswith("stage="):
                    size[part.removeprefix("stage=")] += os.path.getsize(p)
                    break
    stored = {}
    for s in store.stages_on_disk(a.root):
        ci = io.part_files(a.root / store.INDEX_DIR / f"stage={s}")
        (b,), = con.execute(f"SELECT coalesce(sum(contig_len), 0) FROM read_parquet({sql_list(ci)}, "
                            "hive_partitioning = false) WHERE is_default_variant").fetchall()
        stored[s] = {"bytes": size[s], "stored_bases": int(b),
                     "bytes_per_base": round(size[s] / int(b), 4) if b else None}
    out["root_bytes_by_stage"] = stored
    (a.out / "summary.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()

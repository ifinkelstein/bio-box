"""Card 0-08a (M1 of planning/design/0-08_ncrna_models.md): full Rfam 15.1 on whole contigs of a corpus sample.

Run ON THE POD, inside the image, with this tree's `src` first on PYTHONPATH (the driver does it):

    python rfam_sample.py sample <root> <inputs dir> <gtdb dir> <M> [--seed S] [--n-prok N] [--n-phage N]
                                 [--shard-bases B]
    python rfam_sample.py fasta  <root> <M> <cm>
    python rfam_sample.py genes  <root> <M>

`sample` chooses the genomes and cuts the scan into processes (`M/sample/sample.json`, `M/sample/plan.json`,
`M/sample/shards.txt`); `fasta` exports the sample's contigs, one FASTA per process (`M/fa/s<k>.fna`, with the
`.map.json` that `infernal_igr.py collect` reads); `genes` exports the gene calls on the sample's contigs and the
Pfam domains of their proteins (`M/sample/genes.jsonl`), read through the current code from the store. The
report is `rfam_sample_report.py`. Nothing here writes into `<root>`; everything goes under `M`.

The sample (reproducible from the seed and the sorted accession lists):

- Prokaryotes: the genomes of stages P1 and P2 (their accession is the file name of the source object at the same
  index of the stage's `objects.json` as the genome's entry in `genomes.json`; checked here by hashing the sampled
  genomes' FASTA records). Each is placed in its GTDB R226 phylum (domain and phylum of `gtdb_taxonomy`), matched
  by accession: exactly (GCA_/GCF_ with version against GTDB's accession without its RS_/GB_ prefix), else against
  GTDB's `ncbi_genbank_assembly_accession`, else by the accession's digits without prefix and version. Genomes
  with no match are left out (counted). Allocation of `n_prok` genomes over the P phyla present: one each, and the
  remaining `n_prok - P` in proportion to the phylum's genome count (largest remainder; ties by phylum name; never
  more than a phylum has). Within a phylum, `random.Random("0-08a:<seed>:<phylum>").sample(<sorted accessions>,
  k)`.
- Phages: the genomes of stage S1 sorted by `genome_uid` (hex), every `len // n_phage`-th from the first.

The plan: genomes in the order (prokaryotes by accession, then phages by genome_uid), packed greedily into
processes of at most `shard_bases` searched bases (a genome is never split; a circular contig counts its flanks).
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import random
import sys
from collections import defaultdict
from pathlib import Path

SEED = 20261112
N_PROK = 100
N_PHAGE = 300
# 18 Mbp per process: 207 CPU-minutes at 1,450 bases per CPU-second (0-09c, measured) and 227 at E-01's slower
# 1,319 (measured), against the 300 CPU-minute limit (AGENTS 7.2): a margin of at least 24%
SHARD_BASES = 18_000_000
STAGES = ("S1", "P1", "P2")
PROK_STAGES = ("P1", "P2")


def _read_cols(files, cols):
    import pyarrow.parquet as pq

    for f in files:
        t = pq.read_table(f, columns=cols, partitioning=None)
        yield from zip(*(t.column(c).to_pylist() for c in cols))


def _stage_files(root: Path, loc: str, stage: str):
    from biobox import io

    return io.part_files(Path(root) / loc / f"stage={stage}")


def _gtdb(gtdb_dir: Path) -> tuple[dict, dict, dict]:
    exact, genbank, digits = {}, {}, defaultdict(list)
    for name in ("bac120_metadata_r226.tsv.gz", "ar53_metadata_r226.tsv.gz"):
        with gzip.open(gtdb_dir / name, "rt", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                acc = r["accession"]
                bare = acc[3:] if acc[:3] in ("RS_", "GB_") else acc
                tax = dict(x.split("__", 1) for x in r["gtdb_taxonomy"].split(";"))
                v = {"gtdb_accession": acc, "domain": tax["d"], "phylum": tax["p"], "taxonomy": r["gtdb_taxonomy"]}
                exact[bare] = v
                gb = r.get("ncbi_genbank_assembly_accession", "")
                if gb and gb not in ("none", "na"):
                    genbank[gb] = v
                digits[bare[4:].split(".")[0]].append(v)
    return exact, genbank, digits


def _match(acc: str, exact, genbank, digits):
    if acc in exact:
        return exact[acc], "exact"
    if acc in genbank:
        return genbank[acc], "genbank"
    d = digits.get(acc[4:].split(".")[0], [])
    if len(d) == 1:
        return d[0], "digits"
    return None, "none" if not d else "ambiguous"


def _allocate(counts: dict[str, int], n: int) -> dict[str, int]:
    phyla = sorted(counts)
    if n < len(phyla):
        raise SystemExit(f"{n} genomes cannot cover {len(phyla)} phyla")
    rest = n - len(phyla)
    total = sum(counts.values())
    quota = {p: counts[p] / total * rest for p in phyla}  # the rest, in proportion to the phylum's genomes
    alloc = {p: 1 + min(int(quota[p]), counts[p] - 1) for p in phyla}
    left = n - sum(alloc.values())
    # largest remainder first (ties by name); a phylum with no genome left is passed over
    order = sorted(phyla, key=lambda p: (-(quota[p] - int(quota[p])), p))
    while left > 0:
        room = [p for p in order if alloc[p] < counts[p]]
        if not room:
            raise SystemExit(f"only {sum(alloc.values())} genomes in all to choose {n} from")
        for p in room[:left]:
            alloc[p] += 1
        left = n - sum(alloc.values())
    return alloc


def sample(root: Path, inputs: Path, gtdb_dir: Path, M: Path, seed: int, n_prok: int, n_phage: int,
           shard_bases: int) -> None:
    from biobox import fasta, ids

    out = M / "sample"
    out.mkdir(parents=True, exist_ok=True)
    clen, gin = {}, {}
    for s in STAGES:
        for cu, n in _read_cols(_stage_files(root, "l0seq/contig_index", s), ["contig_uid", "contig_len"]):
            clen.setdefault(cu.hex(), n)
        for cu, gu, circ in _read_cols(_stage_files(root, "l0seq/genome_input", s),
                                       ["contig_uid", "genome_uid", "is_circular"]):
            gin[(s, cu.hex())] = (gu.hex(), bool(circ))

    def genome_rec(stage, g, acc):
        cus = g["contig_uids"]
        gus = {gin[(stage, c)][0] for c in cus if (stage, c) in gin}
        gu = gus.pop() if len(gus) == 1 else ids.genome_uid([bytes.fromhex(c) for c in cus]).hex()
        contigs = [{"uid": c, "len": clen[c], "circular": c in set(g.get("circular", [])) or
                    gin.get((stage, c), (None, False))[1]} for c in cus]
        return {"id": acc or gu, "accession": acc, "genome_uid": gu, "stage": stage, "kind": g["kind"],
                "genetic_code": g["genetic_code"], "contigs": contigs, "bases": sum(c["len"] for c in contigs)}

    exact, genbank, digits = _gtdb(gtdb_dir)
    prok, unmatched, how = [], [], defaultdict(int)
    for s in PROK_STAGES:
        gs = json.loads((inputs / s / "genomes.json").read_text())
        objs = json.loads((inputs / s / "objects.json").read_text())
        if len(gs) != len(objs):
            raise SystemExit(f"{s}: {len(gs)} genomes but {len(objs)} objects; the index pairing does not hold")
        for g, o in zip(gs, objs):
            acc = Path(o["path"]).name.removesuffix(".fna")
            rec = genome_rec(s, g, acc)
            rec["path"] = o["path"]
            v, m = _match(acc, exact, genbank, digits)
            how[m] += 1
            if v is None:
                unmatched.append(acc)
                continue
            rec.update(v)
            rec["gtdb_match"] = m
            prok.append(rec)
    by_phylum = defaultdict(list)
    for r in prok:
        by_phylum[f"{r['domain']}:{r['phylum']}"].append(r)
    alloc = _allocate({p: len(v) for p, v in by_phylum.items()}, n_prok)
    chosen = []
    for p in sorted(by_phylum):
        accs = sorted(r["accession"] for r in by_phylum[p])
        pick = set(random.Random(f"0-08a:{seed}:{p}").sample(accs, alloc[p]))
        chosen += [r for r in by_phylum[p] if r["accession"] in pick]
    chosen.sort(key=lambda r: r["accession"])
    # the accession of each chosen genome is the object at its index: its FASTA must hash to its contigs
    for r in chosen:
        got = {ids.contig_uid(x.seq).hex() for x in fasta.iter_records(Path(r["path"]))}
        if not set(c["uid"] for c in r["contigs"]) <= got:
            raise SystemExit(f"{r['accession']}: the genome's contigs are not the records of {r['path']}")

    phages = [genome_rec("S1", g, None) for g in json.loads((inputs / "S1" / "genomes.json").read_text())]
    phages.sort(key=lambda r: r["genome_uid"])
    stride = len(phages) // n_phage
    ph = [phages[i * stride] for i in range(n_phage)]
    for r in ph:
        r.update(domain="Viruses", phylum="phage")

    genomes = chosen + ph
    comp = []
    for p in sorted(by_phylum):
        sel = [r for r in chosen if f"{r['domain']}:{r['phylum']}" == p]
        comp.append({"stratum": p, "available": len(by_phylum[p]), "sampled": len(sel),
                     "bases": sum(r["bases"] for r in sel), "contigs": sum(len(r["contigs"]) for r in sel)})
    comp.append({"stratum": "Viruses:phage (S1)", "available": len(phages), "sampled": len(ph),
                 "bases": sum(r["bases"] for r in ph), "contigs": sum(len(r["contigs"]) for r in ph)})
    doc = {
        "seed": seed, "n_prok": n_prok, "n_phage": n_phage, "phage_stride": stride,
        "prokaryotes_in_stages": sum(how.values()), "gtdb_match": dict(sorted(how.items())),
        "unmatched": sorted(unmatched), "phyla_present": len(by_phylum),
        "archaea_phyla": sorted(p for p in by_phylum if p.startswith("Archaea:")),
        "composition": comp,
        "bases": {"prokaryotes": sum(r["bases"] for r in chosen), "phages": sum(r["bases"] for r in ph),
                  "all": sum(r["bases"] for r in genomes)},
        "contigs": sum(len(r["contigs"]) for r in genomes),
        "circular_contigs": sum(c["circular"] for r in genomes for c in r["contigs"]),
        "genomes": genomes,
    }
    (out / "sample.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")

    shards, cur = [], None
    for r in genomes:
        b = r["bases"]
        if cur is None or (cur["bases"] + b > shard_bases and cur["genomes"]):
            cur = {"k": len(shards), "genomes": [], "bases": 0, "contigs": 0}
            shards.append(cur)
        cur["genomes"].append(r["id"])
        cur["bases"] += b
        cur["contigs"] += len(r["contigs"])
    big = [s for s in shards if s["bases"] > shard_bases]
    plan = {"shard_bases": shard_bases, "shards": shards, "oversize": [s["k"] for s in big],
            "rate_bases_per_cpu_s": 1450, "predicted_cpu_min_max": round(max(s["bases"] for s in shards) / 1450 / 60, 1)}
    (out / "plan.json").write_text(json.dumps(plan, indent=1, sort_keys=True) + "\n")
    (out / "shards.txt").write_text("".join(f"{s['k']:02d} {s['bases']}\n" for s in shards))
    print(json.dumps({k: v for k, v in doc.items() if k not in ("genomes", "unmatched")}, sort_keys=True))
    print(json.dumps({"shards": len(shards), "bases_max": max(s["bases"] for s in shards),
                      "oversize": plan["oversize"], "predicted_cpu_min_max": plan["predicted_cpu_min_max"]}))


def write_fasta(root: Path, M: Path, cm: Path) -> None:
    from biobox.l0seq.store import SeqStore
    from biobox.l1 import ncrna

    doc = json.loads((M / "sample/sample.json").read_text())
    plan = json.loads((M / "sample/plan.json").read_text())
    by_id = {g["id"]: g for g in doc["genomes"]}
    flank = max(max(m.w, m.maxl) for m in ncrna.read_models(cm).values())
    fa = M / "fa"
    fa.mkdir(parents=True, exist_ok=True)
    store = SeqStore(root, list(STAGES))
    tot = 0
    try:
        for s in plan["shards"]:
            stem = fa / f"s{s['k']:02d}"
            if (stem.with_suffix(".map.json")).exists():
                continue
            recmap, bases = {}, 0
            tmp = stem.with_name(stem.name + ".fna.tmp")
            with open(tmp, "wb") as fh:
                for gid in s["genomes"]:
                    for c in by_id[gid]["contigs"]:
                        seq = store.get(bytes.fromhex(c["uid"])).upper()
                        if len(seq) != c["len"]:
                            raise SystemExit(f"{c['uid']}: {len(seq)} bases, contig_index says {c['len']}")
                        f = flank if c["circular"] else 0
                        rec = ncrna._extended(seq, f) if f else seq
                        fh.write(b">" + c["uid"].encode() + b"\n" + rec + b"\n")
                        recmap[c["uid"]] = {"contig": c["uid"], "genome": gid, "len": c["len"],
                                            "circular": c["circular"], "offset": f}
                        bases += len(rec)
            os.replace(tmp, stem.with_suffix(".fna"))
            stem.with_suffix(".map.json").write_text(json.dumps(
                {"mode": "wg", "flank": flank, "records": len(recmap), "bases_searched": bases, "map": recmap},
                sort_keys=True) + "\n")
            tot += bases
            print(json.dumps({"shard": s["k"], "records": len(recmap), "bases_searched": bases}), flush=True)
    finally:
        store.close()
    print(json.dumps({"bases_written": tot}))


def genes(root: Path, M: Path) -> None:
    import duckdb

    from biobox import io
    from biobox.l1.genecall import GENE_CALL_VERSION

    doc = json.loads((M / "sample/sample.json").read_text())
    want = sorted({c["uid"] for g in doc["genomes"] for c in g["contigs"]})
    con = duckdb.connect()
    con.execute("SET threads=1")
    con.execute("SET memory_limit='4GB'")
    (M / "tmp").mkdir(exist_ok=True)
    con.execute(f"SET temp_directory='{M / 'tmp'}'")
    con.execute("CREATE TABLE want(cu BLOB)")
    con.executemany("INSERT INTO want VALUES (?)", [(bytes.fromhex(c),) for c in want])

    def files(rel, stages=STAGES):
        fs = [str(f) for s in stages for f in io.part_files(Path(root) / rel / f"stage={s}")]
        return "[" + ", ".join("'" + f.replace("'", "''") + "'" for f in fs) + "]"

    con.execute(f"""CREATE TABLE g AS SELECT DISTINCT f.contig_uid, f.start, f."end", f.strand, f.feature_uid
                    FROM read_parquet({files(f'l1/feature/producer={GENE_CALL_VERSION}')}) f
                    JOIN want ON f.contig_uid = want.cu""")
    con.execute(f"""CREATE TABLE fp AS SELECT DISTINCT p.feature_uid, p.protein_uid
                    FROM read_parquet({files(f'l1/feature_protein/gene_call_version={GENE_CALL_VERSION}')}) p
                    JOIN want ON p.contig_uid = want.cu""")
    con.execute(f"""CREATE TABLE d AS SELECT d.protein_uid, string_agg(DISTINCT d.model_name, '+'
                    ORDER BY d.model_name) AS pfam
                    FROM read_parquet({files('l1/domains/library=pfam-38.2')}) d
                    WHERE d.protein_uid IN (SELECT protein_uid FROM fp) GROUP BY d.protein_uid""")
    rows = con.execute("""SELECT g.contig_uid, g.start, g."end", g.strand, fp.protein_uid, d.pfam FROM g
                          LEFT JOIN fp ON g.feature_uid = fp.feature_uid LEFT JOIN d ON fp.protein_uid = d.protein_uid
                          ORDER BY g.contig_uid, g.start, g."end", g.strand""").fetchall()
    tmp = M / "sample/genes.jsonl.tmp"
    with open(tmp, "w") as fh:
        for cu, s, e, st, pu, pf in rows:
            fh.write(json.dumps({"contig": cu.hex(), "start": s, "end": e, "strand": st,
                                 "protein_uid": pu.hex() if pu else None, "pfam": pf}, sort_keys=True) + "\n")
    os.replace(tmp, M / "sample/genes.jsonl")
    print(json.dumps({"genes": len(rows), "contigs_with_genes": len({r[0] for r in rows}),
                      "contigs": len(want), "with_protein": sum(r[4] is not None for r in rows),
                      "with_pfam": sum(r[5] is not None for r in rows), "gene_call_version": GENE_CALL_VERSION}))


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2
    cmd, a = argv[0], argv[1:]

    def opt(name, default):
        if name in a:
            i = a.index(name)
            v = a[i + 1]
            del a[i:i + 2]
            return int(v)
        return default

    if cmd == "sample":
        seed, n_prok, n_phage = opt("--seed", SEED), opt("--n-prok", N_PROK), opt("--n-phage", N_PHAGE)
        sb = opt("--shard-bases", SHARD_BASES)
        sample(Path(a[0]), Path(a[1]), Path(a[2]), Path(a[3]), seed, n_prok, n_phage, sb)
    elif cmd == "fasta":
        write_fasta(Path(a[0]), Path(a[1]), Path(a[2]))
    elif cmd == "genes":
        genes(Path(a[0]), Path(a[1]))
    else:
        print(__doc__, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

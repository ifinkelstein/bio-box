"""Card C-48c (part 2) step 6: the cost of building stages V2 and GV1 from the code's own estimate functions, and what
the rest of IMG/VR v4.0 would add, from its table alone.

    python c48c_v2_gv1_est.py v2 --root ROOT --inputs V2_INPUTS --out OUT [--every 40]
    python c48c_v2_gv1_est.py gv1 --root ROOT --gvmags GVMAGS_TAR --contigs GV1_CANDIDATE_CONTIGS_TSV
        --genomes GVMAGS_GENOMES_TSV --held HELD_TSV --out OUT [--every 50]
    python c48c_v2_gv1_est.py imgvr-full --full-zip FULL_ZIP --hc-reps V2_INPUTS --out OUT

Nothing in ROOT is written. Every estimate is the code's (`biobox.stage`, `biobox.l1.genecall`, read from the tree):
ingest (`stage.ingest_cpu_min` of the FASTA bytes, and with the bases another stage holds), gene calling in gv-meta
mode (`genecall.one_process_cpu_min` over the bases the caller reads, every genome with a contig of its own called whole,
a circular contig with its wrap flanks; batches at `stage.SHARD_MARGIN` of the 300 CPU-minute budget), the two CRISPR
producers (`stage.crispr_cpu_min` of the stored bases and the circular wrap flanks, `crispr.CIRC_FLANK`), the Pfam scan
(`stage.SCAN_CPU_MIN_PER_M` per million residues of the distinct proteins; shards of 300 x 0.8 x 10^6 /
SCAN_CPU_MIN_PER_M residues). The distinct proteins and residues are not known before calling: a sample of every
EVERY-th genome that has a contig of its own is called here, one thread, with the step's own call path
(`genecall._call`), and its genes, distinct proteins and residues per Mbp are scaled to the stage (an over-estimate of
distinct proteins: a sample repeats fewer proteins than the whole); how many of the sample's proteins the root's
protein_seq already holds is counted too.

v2: from V2_INPUTS (genomes.json, records.tsv, shared_records.tsv, objects.json). `stage` is the same for any staged
inputs of that form (GV1's, card C-48c part 3); the sample also gives the gene-calling step at the rate it measured,
beside the code's (`genecall_cpu_min_sample_rate`, `genecall_cpu_min_code_rate`).
gv1: GV1 has no inputs (its selection waits on a decision): both candidate sets (all 18,726 genomes of
all_18kgenomes_fna.tar.gz; the 8,508 that fna_files.tar.gz names) from the overlap step's candidate contigs, every
contig linear (the archive records no topology), HELD_TSV the contig uids an earlier stage holds (one per line, hex).
imgvr-full: IMGVR_all_Sequence_information.tsv of IMG/VR v4.0 (the full set, high-confidence and putative) copied out of
FULL_ZIP to OUT (md5 checked against sources.yaml) and read with DuckDB: its vOTUs that V2 does not take (V2's rule
applied to this table: a vOTU with a high-confidence, high-quality member; checked against V2's records.tsv when
--hc-reps names V2's inputs and they exist), the longest
member of each standing for its representative (an upper bound on its bases; v4.0 names no representative either), and
for --min-contig-len values of 1, 5, 10 and 20 kb the vOTUs and bases left and the cost at --depth called (ingest, gene
calling, CRISPR; no Pfam), every genome linear.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import resource
import tarfile
import zipfile
from collections import Counter
from pathlib import Path

import duckdb

from biobox import fasta, ids, stage
from biobox import io as bbio
from biobox.l0seq import store
from biobox.l0seq.inputs import GenomeInput
from biobox.l1 import crispr, genecall

import c48c_v2_gv1_look as look  # noqa: E402
import c48c_virus_inputs as vi  # noqa: E402

FULL_TSV = "Cus_VR/IMG_VR_2022-12-19_7/IMGVR_all_Sequence_information.tsv"
FULL_TSV_MD5 = "658ce076cd26729892552ed5876ad17c"  # sources.yaml, jgi_img_vr
BUDGET = 300


def cpu() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def code_costs(called: dict, stored_bases: int, flank_bases: int, fasta_bytes: int, held_bases: int,
               genes: float | None, residues: float | None) -> dict:
    """The code's estimates of each step (module docstring), CPU-minutes."""
    gc = genecall.one_process_cpu_min(called)
    margin = stage.SHARD_MARGIN * BUDGET
    shard_res = int(BUDGET * stage.SHARD_MARGIN * 1e6 / stage.SCAN_CPU_MIN_PER_M)
    out = {"ingest_cpu_min": round(stage.ingest_cpu_min(fasta_bytes), 2),
           "ingest_cpu_min_with_held": round(stage.ingest_cpu_min(fasta_bytes, held_bases), 2),
           "fasta_bytes": fasta_bytes, "held_bases": held_bases,
           "called_bases": {f"{m} {'circular' if c else 'linear'}": b for (m, c), b in sorted(called.items())},
           "genecall_cpu_min": round(gc, 1), "genecall_batches": math.ceil(gc / margin),
           "genes_code_rate": round(genecall.GENES_PER_MBP["gv-meta"] * stored_bases / 1e6),
           "stored_bases": stored_bases, "crispr_flank_bases": flank_bases,
           "crispr_default_cpu_min": round(stage.crispr_cpu_min(stored_bases, flank_bases), 2),
           "crispr_near_cas_cpu_min": round(stage.crispr_cpu_min(stored_bases, flank_bases, near_cas=True), 2),
           "scan_shard_residues": shard_res}
    if genes is not None:
        out["merge_cpu_min"] = round(genecall.merge_cpu_min(genes), 2)
    if residues is not None:
        out["pfam_residues"] = int(residues)
        out["pfam_cpu_min"] = round(residues / 1e6 * stage.SCAN_CPU_MIN_PER_M, 1)
        out["pfam_shards"] = math.ceil(residues / shard_res)
        out["annotated_cpu_min"] = round(out["ingest_cpu_min_with_held"] + gc + out["crispr_default_cpu_min"]
                                         + out["crispr_near_cas_cpu_min"] + out["pfam_cpu_min"], 1)
        out["called_cpu_min"] = round(out["ingest_cpu_min_with_held"] + gc + out["crispr_default_cpu_min"], 1)
        out["called_saves_cpu_min"] = round(out["annotated_cpu_min"] - out["called_cpu_min"], 1)
    return out


class Sample:
    """Genomes called one at a time in this process and thread; their genes, distinct proteins, residues."""

    def __init__(self):
        self.genomes = self.bases = self.genes = 0
        self.cpu_s = 0.0
        self.prot: dict[bytes, int] = {}

    def call(self, gi: GenomeInput, seqs: dict[bytes, str]) -> None:
        t0 = cpu()
        calls = genecall._call(gi, seqs)
        self.cpu_s += cpu() - t0
        self.genomes += 1
        self.bases += sum(len(s) for s in seqs.values())
        self.genes += len(calls)
        for c in calls:
            self.prot.setdefault(ids.protein_uid(c["aa"]), len(c["aa"]))

    def report(self, root: Path, out: Path) -> dict:
        mbp = self.bases / 1e6
        res = sum(self.prot.values())
        return {"genomes": self.genomes, "bases": self.bases, "genes": self.genes, "cpu_s": round(self.cpu_s, 1),
                "s_per_mbp": round(self.cpu_s / mbp, 4) if mbp else None,
                "genes_per_mbp": round(self.genes / mbp, 2) if mbp else None,
                "distinct_proteins": len(self.prot), "distinct_per_mbp": round(len(self.prot) / mbp, 2) if mbp else None,
                "residues": res, "residues_per_mbp": round(res / mbp, 1) if mbp else None,
                "distinct_in_root": known(root, sorted(self.prot), out)}


def known(root: Path, uids: list[bytes], out: Path) -> int:
    """How many of `uids` the root's protein_seq holds (any stage, its one gene-call version)."""
    gcv = sorted(os.listdir(root / "l0seq" / "protein_seq"))
    if len(gcv) != 1:
        raise SystemExit(f"gene call versions {gcv}: expected one")
    ps = [f for s in store.stages_on_disk(root) for f in bbio.part_files(
        root / "l0seq/protein_seq" / gcv[0] / f"stage={s}", "protein_prefix=*/part-*.parquet")]
    con = duckdb.connect(config={"threads": 1, "memory_limit": "3GB", "temp_directory": str(out / "duckdb")})
    con.execute("CREATE TABLE q (protein_uid BLOB)")
    con.executemany("INSERT INTO q VALUES (?)", [(u,) for u in uids])
    (n,), = con.execute(f"SELECT count(DISTINCT protein_uid) FROM read_parquet({sql_list(ps)}, "
                        "hive_partitioning = false) WHERE protein_uid IN (SELECT protein_uid FROM q)").fetchall()
    con.close()
    return n


def scaled(sample: dict, stored_bases: int) -> tuple[float, float, float]:
    mbp = stored_bases / 1e6
    return sample["genes_per_mbp"] * mbp, sample["distinct_per_mbp"] * mbp, sample["residues_per_mbp"] * mbp


def est_v2(root: Path, inp: Path, out: Path, every: int) -> dict:
    with open(inp / "records.tsv", newline="") as fh:
        recs = list(csv.DictReader(fh, delimiter="\t"))
    with open(inp / "shared_records.tsv", newline="") as fh:
        held = {r["contig_uid"] for r in csv.DictReader(fh, delimiter="\t")}
    genomes = json.loads((inp / "genomes.json").read_text())
    length = {r["contig_uid"]: int(r["length"]) for r in recs}
    called = Counter()
    sample_want: dict[str, list[int]] = {}  # contig -> the sampled genomes (indices) that list it
    k = 0
    for j, g in enumerate(genomes):
        if all(u in held for u in g["contig_uids"]):
            continue  # every contig another stage owns: not called (genecall item 5)
        circ = set(g["circular"])
        for u in g["contig_uids"]:
            called[("gv-meta", u in circ)] += genecall.called_bases(length[u], u in circ)
        k += 1
        if k % every == 0:
            for u in g["contig_uids"]:
                sample_want.setdefault(u, []).append(j)
    stored = {u: n for u, n in length.items() if u not in held}
    flank = sum(2 * min(crispr.CIRC_FLANK, length[u]) for g in genomes for u in g["circular"] if u not in held)
    objects = json.loads((inp / "objects.json").read_text())
    fasta_bytes = sum(Path(o["path"]).stat().st_size for o in objects)
    # the sample: contigs read back from the inputs' FASTA objects
    s = Sample()
    seqs: dict[int, dict[bytes, str]] = {}
    done: set[int] = set()
    for o in objects:
        for r in fasta.iter_records(Path(o["path"])):
            cu = ids.contig_uid(r.seq)
            for j in sample_want.get(cu.hex(), ()):
                if j in done:
                    continue
                g = genomes[j]
                seqs.setdefault(j, {})[cu] = r.seq.upper().decode("ascii")
                if len(seqs[j]) == len(g["contig_uids"]):
                    uids = tuple(bytes.fromhex(u) for u in g["contig_uids"])
                    gi = GenomeInput(uids, g["genetic_code"], kind=g["kind"],
                                     circular=tuple(bytes.fromhex(u) for u in g["circular"]),
                                     assembly_quality=g["assembly_quality"])
                    s.call(gi, seqs.pop(j))
                    done.add(j)
    sam = s.report(root, out)
    sam["sampled_genomes_not_completed"] = len({j for js in sample_want.values() for j in js} - done)
    gc_code = genecall.one_process_cpu_min(dict(called))
    # the same step at the rate this sample measured (calls alone, then the code's own overhead and process cost)
    gc_sample = (genecall.PROCESS_S + genecall.STEP_OVERHEAD * sam["s_per_mbp"] * sum(called.values()) / 1e6) / 60
    sam["genecall_cpu_min_code_rate"] = round(gc_code, 1)
    sam["genecall_cpu_min_sample_rate"] = round(gc_sample, 1)
    sam["sample_rate_over_code_rate"] = round(gc_sample / gc_code, 3) if gc_code else None
    genes, distinct, residues = scaled(sam, sum(stored.values()))
    held_bases = sum(n for u, n in length.items() if u in held)
    return {"sample": {"every": every, **sam},
            "derived": {"genes": round(genes), "distinct_proteins": round(distinct), "residues": round(residues)},
            "code": code_costs(dict(called), sum(stored.values()), flank, fasta_bytes, held_bases, genes, residues),
            "storage_gb": {"at_s1_0.7175": round(sum(stored.values()) / 1e9 * 0.7175, 2),
                           "at_h1_0.7316": round(sum(stored.values()) / 1e9 * 0.7316, 2)}}


def est_gv1(root: Path, tar: Path, contigs: Path, genomes_tsv: Path, held_tsv: Path, out: Path, every: int) -> dict:
    with open(genomes_tsv, newline="") as fh:
        fna_set = {r["genome"] for r in csv.DictReader(fh, delimiter="\t") if r["set"] == "fna_files"}
    held = set(held_tsv.read_text().split())
    per: dict[str, list[tuple[str, int]]] = {}
    with open(contigs, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            per.setdefault(r["genome"], []).append((r["contig_uid"], int(r["length"])))
    sets = {"all_18k": set(per), "fna_files_8508": fna_set}
    res: dict = {}
    for name, members in sets.items():
        called = 0
        stored: dict[str, int] = {}
        hb = 0
        for g in sorted(members):
            cs = per[g]
            if all(u in held for u, _ in cs):
                continue  # every contig an earlier stage holds: not called
            for u, n in cs:
                called += n
                if u in held:
                    hb += n
                else:
                    stored[u] = n
        res[name] = {"genomes": len(members), "stored_contigs": len(stored), "stored_bases": sum(stored.values()),
                     "called_bases": called, "held_bases": hb}
    sample_names = {g for i, g in enumerate(sorted(sets["all_18k"]), 1) if i % every == 0}
    s = Sample()
    member_bytes = Counter()
    with look.open_member(tar, "gvmags_v2/all_18kgenomes_fna.tar.gz") as fh, tarfile.open(fileobj=fh, mode="r|gz") as tf:
        for m in tf:
            if not m.isfile() or not m.name.endswith(".fna"):
                continue
            g = m.name.split("/")[-1][:-4]
            member_bytes["all_18k"] += m.size
            if g in fna_set:
                member_bytes["fna_files_8508"] += m.size
            if g not in sample_names:
                continue
            seqs = {}
            for r in vi.iter_fasta_stream(io.BytesIO(tf.extractfile(m).read())):
                seqs.setdefault(ids.contig_uid(r.seq), r.seq.upper().decode("ascii"))
            # kind virus is gv-meta mode; the quality value does not enter the calls
            gi = GenomeInput(tuple(seqs), 11, kind="virus", circular=(), assembly_quality="contig")
            s.call(gi, seqs)
    sam = s.report(root, out)
    for name, r in res.items():
        genes, distinct, residues = scaled(sam, r["stored_bases"])
        r["derived"] = {"genes": round(genes), "distinct_proteins": round(distinct), "residues": round(residues)}
        r["code"] = code_costs({("gv-meta", False): r["called_bases"]}, r["stored_bases"], 0, member_bytes[name],
                               r["held_bases"], genes, residues)
        r["storage_gb_at_0.7316"] = round(r["stored_bases"] / 1e9 * 0.7316, 2)
    res["sample"] = {"every": every, "of": "all_18k genomes in name order", **sam}
    return res


def est_imgvr_full(zip_path: Path, v2: Path | None, out: Path) -> dict:
    dst = out / "IMGVR_all_Sequence_information.tsv"
    if not dst.is_file():
        tmp = out / ".IMGVR_all_Sequence_information.tsv.tmp"
        md5 = hashlib.md5()
        with zipfile.ZipFile(zip_path) as z, z.open(FULL_TSV) as r, open(tmp, "wb") as w:
            for blk in iter(lambda: r.read(1 << 22), b""):
                md5.update(blk)
                w.write(blk)
        if md5.hexdigest() != FULL_TSV_MD5:
            tmp.unlink()
            raise SystemExit(f"md5 of {FULL_TSV} is {md5.hexdigest()}, sources.yaml says {FULL_TSV_MD5}")
        tmp.replace(dst)
    con = duckdb.connect(config={"threads": 2, "memory_limit": "6GB", "temp_directory": str(out / "duckdb")})
    cols = "{" + ", ".join(f"'{c}': 'VARCHAR'" for c in look.IMGVR_COLS) + "}"
    con.execute(f"CREATE VIEW f AS SELECT * FROM read_csv('{dst}', delim='\t', header=true, quote='', escape='', "
                f"columns={cols})")
    q = lambda sql: con.execute(sql).fetchall()  # noqa: E731
    # V2's vOTUs by V2's rule applied to this table (a high-confidence, high-quality member); compared with V2's
    # records.tsv when it is given
    con.execute("CREATE TABLE v2 AS SELECT DISTINCT votu FROM f WHERE confidence = 'High-confidence' AND "
                "miuvig = 'High-quality'")
    (rows, votus, bases), = q("SELECT count(*), count(DISTINCT votu), sum(length::BIGINT) FROM f")
    res = {"rows": rows, "votus": votus, "bases": int(bases),
           "by_confidence": {k: [n, int(b)] for k, n, b in q(
               "SELECT confidence, count(*), sum(length::BIGINT) FROM f GROUP BY 1 ORDER BY 2 DESC")},
           "by_miuvig": {k: [n, int(b)] for k, n, b in q(
               "SELECT miuvig, count(*), sum(length::BIGINT) FROM f GROUP BY 1 ORDER BY 2 DESC")},
           "votus_by_v2_rule": q("SELECT count(*) FROM v2")[0][0],
           "votus_with_hq_member_any_confidence": q(
               "SELECT count(DISTINCT votu) FROM f WHERE miuvig = 'High-quality'")[0][0]}
    rec = v2 / "records.tsv" if v2 is not None else None
    if rec is not None and rec.is_file():
        con.execute(f"CREATE TABLE staged AS SELECT DISTINCT votu FROM read_csv('{rec}', delim='\t', header=true, "
                    "quote='', all_varchar=true) WHERE source = 'imgvr'")
        res["staged_v2_votus"] = q("SELECT count(*) FROM staged")[0][0]
        res["staged_v2_votus_not_by_rule_here"] = q(
            "SELECT count(*) FROM staged WHERE votu NOT IN (SELECT votu FROM v2)")[0][0]
    con.execute("CREATE TABLE rest AS SELECT votu, max(length::BIGINT) AS n, count(*) AS members, "
                "bool_or(confidence = 'High-confidence') AS hc FROM f WHERE votu NOT IN (SELECT votu FROM v2) "
                "GROUP BY votu")
    res["rest"] = {}
    for L in (1, 1000, 5000, 10000, 20000):
        (nv, nb, nhc), = q(f"SELECT count(*), coalesce(sum(n), 0), count(*) FILTER (WHERE hc) FROM rest WHERE n >= {L}")
        nb = int(nb)
        called = {("gv-meta", False): nb}
        gc = genecall.one_process_cpu_min(called)
        res["rest"][f"min_contig_len_{L}"] = {
            "votus": nv, "votus_with_a_high_confidence_member": nhc, "bases": nb,
            "fasta_bytes_estimate": round(nb * 61 / 60),
            "ingest_cpu_min": round(stage.ingest_cpu_min(round(nb * 61 / 60)), 1),
            "genecall_cpu_min": round(gc, 1), "genecall_batches": math.ceil(gc / (stage.SHARD_MARGIN * BUDGET)),
            "crispr_cpu_min": round(stage.crispr_cpu_min(nb), 1),
            "genes_code_rate": round(genecall.GENES_PER_MBP["gv-meta"] * nb / 1e6),
            "called_cpu_min": round(stage.ingest_cpu_min(round(nb * 61 / 60)) + gc + stage.crispr_cpu_min(nb), 1),
            "storage_gb_at_0.71": round(nb / 1e9 * 0.71, 1)}
    con.close()
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["v2", "stage", "gv1", "imgvr-full"])  # stage: v2's estimate of any staged inputs
    ap.add_argument("--root", type=Path)
    ap.add_argument("--inputs", type=Path)
    ap.add_argument("--gvmags", type=Path)
    ap.add_argument("--contigs", type=Path)
    ap.add_argument("--genomes", type=Path)
    ap.add_argument("--held", type=Path)
    ap.add_argument("--full-zip", type=Path)
    ap.add_argument("--hc-reps", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--every", type=int, default=40)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if a.mode in ("v2", "stage"):
        res = est_v2(a.root, a.inputs, a.out, a.every)
    elif a.mode == "gv1":
        res = est_gv1(a.root, a.gvmags, a.contigs, a.genomes, a.held, a.out, a.every)
    else:
        res = est_imgvr_full(a.full_zip, a.hc_reps, a.out)
    (a.out / f"{a.mode}.json").write_text(json.dumps(res, indent=1, sort_keys=True, default=str))
    print(json.dumps(res, sort_keys=True, default=str))


if __name__ == "__main__":
    main()

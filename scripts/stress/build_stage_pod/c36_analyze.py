"""Card C-36 (3) to (5): what the six phase-2 stages share with the stored stages, with PL1's inputs, with each other
and within themselves, by content digest; and each stage's estimates from the code's own functions. Read only.

    python c36_analyze.py --root ROOT --pl1 PL1_INPUTS --stages C36_STAGES --inputs INPUTS --out OUT
                          [--rehearsal CHECK_JSON]

ROOT is the data root (read only; only the `contig_index` of the stored stages S1, P1 to P6, H1 is read, by name, so a
stage being built beside this, PL1, is never read). PL1_INPUTS is the plasmid stage's inputs directory (card C-26 part
1: genomes.json, shared_records.tsv). C36_STAGES holds `c36_stages.py assign`'s output with each stage's contigs.tsv
(`c36_verify.py`: every record's `ids.contig_uid` and length); INPUTS holds each stage's inputs (`C<k>/genomes.json`,
`make_inputs.py prok`). Writes OUT/summary.json (and prints it) and OUT/C<k>_contigs_elsewhere.tsv (each contig of
the stage that another stage stores or lists, with where).

Identity is the store's: a contig is `ids.contig_uid` of its uppercase bytes; a stored stage owns a contig when its
`contig_index` holds it as the default variant (DEC-26, as `make_inputs.stage_owners` reads it). The build order
assumed is the stored stages, then PL1 (which stores every contig of its inputs that no stored stage owns), then C1 to
C6: a contig is stored by the first of these that lists it.

Estimates (labels in the note): ingest from `stage.ingest_cpu_min` over the stage's bytes of FASTA (`fasta_bytes`, the
decompressed bytes, as `stage.object_sizes` counts them) and whether it exceeds `SHARD_MARGIN` of the 300 CPU-minute
budget (then `build-stage` runs it in batches); gene calling from `genecall.cut_plan` over the stage's genomes in
genome_uid order with the item rule of `genecall.stage_items` (a genome is called when one of its contigs is the
stage's own: stored by the stage and its least genome is this genome, DEC-26, DEC-27; its cost `genome_cost_us`,
its genes `genome_genes`), with `batch_cpu_min` = `SHARD_MARGIN` x 300 as `build-stage` passes it; modes by
`genecall.mode_of`. With --rehearsal (`c36_rh_check.py`'s output), distinct proteins and residues from the rehearsal's
measured ratios, and the full Pfam scan from `stage.scan_cpu_min` and `stage.shard_residues`.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow.parquet as pq

from biobox import ids, io, stage
from biobox.l0seq import store
from biobox.l1 import genecall

STORED = ["S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1"]
PHASE2 = [f"C{k}" for k in range(1, 7)]
BUDGET = 300  # cpu_min_per_process of the pod's profile (AGENTS §7.2)


class _Prof:
    cpu_min_per_process = BUDGET


def tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser()
    for a in ("--root", "--pl1", "--stages", "--inputs", "--out"):
        ap.add_argument(a, required=True, type=Path)
    ap.add_argument("--rehearsal", type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    out: dict = {"budget_cpu_min": BUDGET, "shard_margin": stage.SHARD_MARGIN}

    # owners in the stored stages, read by name
    owner: dict[bytes, str] = {}
    stored_len: dict[bytes, int] = {}
    for s in STORED:
        files = io.part_files(a.root / store.INDEX_DIR / f"stage={s}")
        if not files:
            raise SystemExit(f"stage {s}: no contig_index in {a.root}")
        for f in files:
            t = pq.read_table(f, columns=["contig_uid", "is_default_variant", "contig_len"])
            for cu, d, n in zip(*(t.column(c).to_pylist() for c in ("contig_uid", "is_default_variant", "contig_len"))):
                if d:
                    if owner.get(cu, s) != s:
                        raise SystemExit(f"contig {cu.hex()} is a default variant of {owner[cu]} and {s}")
                    owner[cu], stored_len[cu] = s, n
    out["stored"] = {s: sum(1 for v in owner.values() if v == s) for s in STORED}

    # PL1's inputs: every contig its genomes list; those no stored stage owns PL1 would store
    pl1 = {bytes.fromhex(u) for g in json.loads((a.pl1 / "genomes.json").read_text()) for u in g["contig_uids"]}
    pl1_owned = {c for c in pl1 if c in owner}
    pl1_shared_tsv = {bytes.fromhex(r["contig_uid"]) for r in tsv(a.pl1 / "shared_records.tsv")}
    out["pl1"] = {"contigs": len(pl1), "owned_by_stored": len(pl1_owned), "would_store": len(pl1) - len(pl1_owned),
                  "owned_by_stored_equals_its_shared_records": pl1_owned == pl1_shared_tsv,
                  "owned_by_stored_by_stage": dict(sorted(Counter(owner[c] for c in pl1_owned).items()))}

    # the phase-2 stages
    recs: dict[str, list[dict]] = {}
    length: dict[bytes, int] = {}
    genomes: dict[str, list[dict]] = {}
    contigs: dict[str, set[bytes]] = {}
    for k in PHASE2:
        recs[k] = tsv(a.stages / k / "contigs.tsv")
        for r in recs[k]:
            cu, n = bytes.fromhex(r["contig_uid"]), int(r["length"])
            if length.setdefault(cu, n) != n:
                raise SystemExit(f"contig {r['contig_uid']} has two lengths")
        genomes[k] = json.loads((a.inputs / k / "genomes.json").read_text())
        contigs[k] = {bytes.fromhex(r["contig_uid"]) for r in recs[k]}
        listed = {bytes.fromhex(u) for g in genomes[k] for u in g["contig_uids"]}
        if listed != contigs[k]:
            raise SystemExit(f"{k}: genomes.json lists {len(listed)} contigs, its records hold {len(contigs[k])}")
    bad_len = [c for c in length if c in stored_len and stored_len[c] != length[c]]
    if bad_len:
        raise SystemExit(f"{len(bad_len)} contigs whose stored length differs, e.g. {bad_len[0].hex()}")

    first_c: dict[bytes, str] = {}
    for k in PHASE2:
        for c in contigs[k]:
            first_c.setdefault(c, k)

    def stored_by(c: bytes) -> str:  # the build order: stored stages, PL1, C1..C6
        return owner.get(c) or ("PL1" if c in pl1 else first_c[c])

    per: dict[str, dict] = {}
    estimates: dict[str, dict] = {}
    for k in PHASE2:
        cs = contigs[k]
        bases = lambda xs: sum(length[c] for c in xs)  # noqa: E731
        members = Counter(bytes.fromhex(u) for g in genomes[k] for u in g["contig_uids"])
        rec_count = Counter(bytes.fromhex(r["contig_uid"]) for r in recs[k])
        accs = defaultdict(list)
        for r in recs[k]:
            accs[bytes.fromhex(r["contig_uid"])].append(r["record"])
        in_stored = {c for c in cs if c in owner}
        in_pl1 = {c for c in cs if c in pl1}
        others = {j: {c for c in cs if c in contigs[j]} for j in PHASE2 if j != k}
        repeated = {c for c, n in members.items() if n > 1}
        adds = {c for c in cs if stored_by(c) == k}
        lock = json.loads((a.stages / k / "lock.json").read_text())["genomes"]
        acc_set: dict[str, list[bytes]] = defaultdict(list)
        for r in recs[k]:
            acc_set[r["accession"]].append(bytes.fromhex(r["contig_uid"]))
        first_of: dict[bytes, str] = {}
        same = []  # genomes whose contig set is an earlier genome's (accession order, as make_inputs.py lists them)
        for acc in sorted(acc_set):
            gu = ids.genome_uid(list(dict.fromkeys(acc_set[acc])))
            if gu in first_of:
                same.append([acc, first_of[gu]])
            else:
                first_of[gu] = acc
        per_same = same
        per[k] = {
            "identical_genomes": per_same,
            "genomes_in_lock": len(lock), "genomes_listed": len(genomes[k]),
            "records": len(recs[k]), "record_bases": sum(int(r["length"]) for r in recs[k]),
            "distinct_contigs": len(cs), "distinct_bases": bases(cs),
            "stored_by_S1_P_H1": {"contigs": len(in_stored), "bases": bases(in_stored),
                                  "records": sum(rec_count[c] for c in in_stored),
                                  "by_stage": {s: sum(1 for c in in_stored if owner[c] == s)
                                               for s in STORED if any(owner[c] == s for c in in_stored)}},
            "in_pl1_inputs": {"contigs": len(in_pl1), "bases": bases(in_pl1),
                              "of_which_stored_by_S1_P_H1": len(in_pl1 & in_stored),
                              "pl1_would_store": len(in_pl1 - in_stored),
                              "pl1_would_store_bases": bases(in_pl1 - in_stored)},
            "shared_with_other_phase2": {j: {"contigs": len(v), "bases": bases(v)} for j, v in others.items() if v},
            "shared_with_any_other_phase2": len(set().union(*others.values())),
            "repeated_within": {"contigs_in_several_genomes": len(repeated), "bases": bases(repeated),
                                "memberships_beyond_the_first": sum(n - 1 for n in members.values()),
                                "records_beyond_the_first": sum(n - 1 for n in rec_count.values()),
                                "largest_genome_count": max(members.values(), default=0)},
            "adds_after_pl1_in_order": {"contigs": len(adds), "bases": bases(adds)},
            "stored_elsewhere_before_it": {"contigs": len(cs - adds), "bases": bases(cs - adds),
                                           "by_owner": dict(sorted(Counter(stored_by(c) for c in cs - adds).items()))},
        }
        rows = []
        for c in sorted(cs):
            where = [s for s in ([owner[c]] if c in owner else []) + (["PL1"] if c in pl1 else [])
                     + [j for j in PHASE2 if j != k and c in contigs[j]]]
            if where or members[c] > 1:
                rows.append((c.hex(), length[c], stored_by(c), ",".join(where), members[c], ",".join(sorted(accs[c]))))
        with open(a.out / f"{k}_contigs_elsewhere.tsv", "w") as w:
            w.write("contig_uid\tlength\tstored_by\tlisted_elsewhere\tgenomes_in_stage\trecords\n")
            for r in rows:
                w.write("\t".join(map(str, r)) + "\n")

        # gene calling: the items of genecall.stage_items for the build order above
        least: dict[bytes, bytes] = {}
        gl = []
        for g in genomes[k]:
            cus = [bytes.fromhex(u) for u in g["contig_uids"]]
            guid = ids.genome_uid(cus)
            gl.append((guid, cus, g))
            for c in cus:
                least[c] = min(least.get(c, guid), guid)
        gl.sort(key=lambda x: x[0])
        items = []
        modes = Counter()
        mode_bases = Counter()
        for guid, cus, g in gl:
            circ = {bytes.fromhex(u) for u in g["circular"]}
            cl = [(length[c], c in circ) for c in cus]
            b = sum(n for n, _ in cl)
            mine = any(stored_by(c) == k and least[c] == guid for c in cus)
            if mine:
                m = genecall.mode_of(g["kind"], b)
                modes[m] += 1
                mode_bases[m] += b
                items.append((guid, genecall.genome_cost_us(g["kind"], cl), b, genecall.genome_genes(g["kind"], b), True))
            else:
                items.append((guid, 0, b, 0.0, False))
        plan = genecall.cut_plan(lambda: iter(items), stage.SHARD_MARGIN * BUDGET)
        fb = sum(int(e["bytes"]) for e in lock.values())
        ing = stage.ingest_cpu_min(fb)
        est = {
            "ingest": {"fasta_bytes": fb, "gz_bytes": sum(int(e["gz_bytes"]) for e in lock.values()),
                       "cpu_min": round(ing, 1), "in_batches": ing > stage.SHARD_MARGIN * BUDGET},
            "genecall": {"genomes": plan["genomes"], "not_called": plan["genomes_not_called"],
                         "modes": dict(modes), "bases_by_mode": dict(mode_bases), "bases_called": plan["bases_called"],
                         "one_process_cpu_min": plan["estimate_cpu_min"],
                         "in_batches": plan["estimate_cpu_min"] > stage.SHARD_MARGIN * BUDGET,
                         "batches": len(plan["batches"]),
                         "batch_cpu_min": [b["estimate_cpu_min"] for b in plan["batches"]],
                         "merge_cpu_min": plan["merge_estimate_cpu_min"],
                         "largest_genome_cpu_min": (plan["largest_genome"] or {}).get("estimate_cpu_min"),
                         "genes_expected": plan["genes_expected"]},
        }
        estimates[k] = est
    out["phase2"] = per
    out["estimates"] = estimates

    # all six together: distinct contigs, and what the six add after PL1
    allc = set().union(*contigs.values())
    out["phase2_all"] = {"distinct_contigs": len(allc), "distinct_bases": sum(length[c] for c in allc),
                         "stored_by_S1_P_H1": sum(1 for c in allc if c in owner),
                         "in_pl1_inputs": sum(1 for c in allc if c in pl1),
                         "in_several_phase2_stages": sum(1 for c in allc if sum(c in contigs[k] for k in PHASE2) > 1),
                         "adds_after_pl1": sum(per[k]["adds_after_pl1_in_order"]["contigs"] for k in PHASE2),
                         "adds_after_pl1_bases": sum(per[k]["adds_after_pl1_in_order"]["bases"] for k in PHASE2)}

    if a.rehearsal is not None:
        rh = json.loads(a.rehearsal.read_text())["ratios"]
        sel = {k: json.loads((a.stages / "assign.json").read_text())["stages"][k] for k in PHASE2}
        for k in PHASE2:
            g = estimates[k]["genecall"]["genes_expected"]
            share_first = 1 - sel[k]["share_bp_not_first_of_species"]
            g1 = g * share_first
            d = rh["distinct_per_gene_first"] * g1 + rh["new_per_gene_later"] * (g - g1)
            res = d * rh["residues_per_distinct"]
            estimates[k]["proteins"] = {
                "genes": round(g), "genes_first_of_species": round(g1), "distinct": round(d),
                "distinct_if_one_per_species_ratio": round(g * rh["distinct_per_gene_first"]),
                "residues": round(res), "gtdb_protein_count": sel[k]["gtdb_protein_count"]}
            estimates[k]["pfam_full"] = {
                "cpu_min": round(stage.scan_cpu_min(res), 1),
                "shards": math.ceil(res / stage.shard_residues(_Prof())),
                "shard_residues": stage.shard_residues(_Prof())}
        out["rehearsal_ratios"] = rh
    (a.out / "summary.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()

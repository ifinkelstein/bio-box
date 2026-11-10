"""Step 4: what changes when one ncRNA feature is added to E. coli K-12 (measurements file section 4).

    python step4_keys.py ROOT SCRATCH RELEASE > keys.json

ROOT is a golden store built by `build_s0.py ... RELEASE --cluster N` (stage A, genome inputs,
clusters, tokens and postings of RELEASE for the gene-call producer alone). SCRATCH must not exist; the
script copies ROOT to SCRATCH/after (without its manifests and catalog, and without
`l3/token_dict/release=RELEASE`, which the second token build writes again) and never writes under ROOT.

In the copy it adds one kind-2 (ncRNA) feature as a second producer `b18-ncrna`, 377 nt on the minus
strand, 60 nt into the first IGR of at least 500 nt at or after 3.27 Mb (the same place as the
2026-11-09 measurement), and runs the repo's own builders for the tape version of both producers:
`igr.run`, `tapes.run`, `tokens.run` and `postings.run`. Then every table is compared before (tape
version of the gene caller alone) and after (both producers), twice: keyed as the registry keys it
(by `feature_uid`, card B-17), and keyed by the slot's ordinal, the key phase B removed. A row is
"rewritten" when its key is gone or any of its other columns changed.
"""

import json
import shutil
import sys
from pathlib import Path

if len(sys.argv) < 4 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from biobox import ids, io  # noqa: E402
from biobox.l1 import igr  # noqa: E402
from biobox.l1.genecall import GENE_CALL_VERSION as V  # noqa: E402
from biobox.l3 import postings, tapes, tokens  # noqa: E402
from biobox.registry.schemas import table  # noqa: E402

root, scratch, rel = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
if scratch.exists():
    sys.exit(f"{scratch}: exists; give a fresh directory")
NEW = "b18-ncrna"
P1, P2 = [V], [V, NEW]
TV1, TV2 = tapes.tape_version(P1), tapes.tape_version(P2)
after = scratch / "after"


def rows(base: Path, pattern: str, columns=None) -> list[dict]:
    return [r for f in io.part_files(base, pattern) for r in pq.read_table(f, columns=columns, partitioning=None).to_pylist()]


# --- the copy and the new feature ----------------------------------------------------------------
shutil.copytree(root, after, ignore=lambda d, names: [n for n in names if Path(d) == root and n in ("manifests", "catalog")])
shutil.rmtree(after / f"l3/token_dict/release={rel}")
idx = {r["contig_uid"]: r["contig_len"] for r in rows(root / "l0seq/contig_index", "stage=S0/part-*.parquet")}
ecoli = next(cu for cu, n in idx.items() if n == 4_641_652)
igr1 = rows(root / f"l0seq/igr/igr_version={igr.igr_version(P1)}", "stage=S0/part-*.parquet")
target = next(g for g in sorted(igr1, key=lambda g: g["start"])
              if g["contig_uid"] == ecoli and g["start"] >= 3_270_000 and g["end"] - g["start"] >= 500)
s, e, strand, kind = target["start"] + 60, target["start"] + 60 + 377, -1, ids.KIND["ncrna"]
fu = ids.feature_uid(ecoli, s, e, strand, kind)
spec = table("feature")
io.write_table([pa.RecordBatch.from_pylist([{"feature_uid": fu, "contig_uid": ecoli, "start": s, "end": e, "strand": strand,
                                              "kind": kind}], schema=spec.schema)],
               spec, after / f"l1/feature/producer={NEW}/stage=S0/part-0.parquet")
built = {"igr": igr.run(after, "S0", producers=P2), "tapes": tapes.run(after, "S0", producers=P2),
         "tokens": tokens.run(after, rel, producers=P2), "postings": postings.run(after, rel, producers=P2)}


# --- comparison ----------------------------------------------------------------------------------
def diff(a: dict, b: dict) -> dict:
    """a, b: key -> rest of the row. Rewritten: key gone or row changed."""
    same = sum(1 for k, v in a.items() if k in b and b[k] == v)
    return {"before": len(a), "after": len(b), "unchanged": same, "rewritten": len(a) - same,
            "new": sum(1 for k in b if k not in a)}


def tape_slots(tape_rows: list[dict], tok_rows: list[dict], by_feature: bool) -> tuple[dict, dict]:
    """Slot entries of the tapes (all list columns) and of the tape tokens, keyed by feature_uid or ordinal."""
    lists = [c for c in table("contig_tapes").list_columns if c != "feature_uids"]
    toks = {r["contig_uid"]: r for r in tok_rows}
    slots, tslots = {}, {}
    for t in tape_rows:
        k = toks[t["contig_uid"]]
        for i, f in enumerate(t["feature_uids"]):
            key = (t["contig_uid"], f) if by_feature else (t["contig_uid"], i)
            rest = tuple((t[c] or [None] * t["n_slots"])[i] for c in lists)
            slots[key] = rest + ((i,) if by_feature else (f,))
            tslots[key] = tuple(k[f"{lvl}_tokens"][i] for lvl in tokens.LEVELS) + ((i,) if by_feature else (f,))
    return slots, tslots


def snapshot(base: Path, producers: list[str], tv: str) -> dict:
    return {
        "feature": rows(base / "l1/feature", "producer=*/stage=S0/part-*.parquet"),
        "feature_protein": rows(base / f"l1/feature_protein/gene_call_version={V}", "stage=S0/part-*.parquet"),
        "protein_occurrence": rows(base / f"l1/protein_occurrence/gene_call_version={V}", "stage=S0/part-*.parquet"),
        "igr": rows(base / f"l0seq/igr/igr_version={igr.igr_version(producers)}", "stage=S0/part-*.parquet"),
        "tapes": rows(base / f"l3/tapes/tape_version={tv}", "stage=S0/part-*.parquet"),
        "tape_tokens": rows(base / f"l3/tape_tokens/tape_version={tv}", f"release={rel}/part-*.parquet"),
        "token_dict": rows(base / "l3/token_dict", f"release={rel}/part-*.parquet"),
        "token_postings": rows(base / "l3/token_postings", f"release={rel}/tape_version={tv}/token_bucket=*/part-*.parquet"),
    }


b0, b1 = snapshot(root, P1, TV1), snapshot(after, P2, TV2)
out = {"contig": "E. coli K-12 MG1655 chromosome", "contig_is_circular": next(t["is_circular"] for t in b0["tapes"]
                                                                               if t["contig_uid"] == ecoli),
       "tape_versions": [TV1, TV2], "ncRNA": {"start": s, "end": e, "strand": strand, "feature_uid": fu.hex()},
       "ecoli_slots_before": next(t["n_slots"] for t in b0["tapes"] if t["contig_uid"] == ecoli),
       "built": built, "tables": {}}
ordinal_new = next(t for t in b1["tapes"] if t["contig_uid"] == ecoli)["feature_uids"].index(fu)
out["ncRNA"]["ordinal_in_new_tape"] = ordinal_new
T = out["tables"]


def keyed(rs: list[dict], key: tuple[str, ...]) -> dict:
    return {tuple(r[k] for k in key): tuple(v for c, v in sorted(r.items()) if c not in key) for r in rs}


T["feature (feature_uid; all producers)"] = diff(keyed(b0["feature"], ("feature_uid", "contig_uid")),
                                                  keyed(b1["feature"], ("feature_uid", "contig_uid")))
for name in ("feature_protein", "protein_occurrence"):
    T[f"{name} (feature_uid)"] = diff(keyed(b0[name], ("feature_uid",)), keyed(b1[name], ("feature_uid",)))
T["igr (contig_uid, start, end, strand); flanks by feature_uid"] = diff(
    keyed(b0["igr"], ("contig_uid", "start", "end", "strand")), keyed(b1["igr"], ("contig_uid", "start", "end", "strand")))
T["contig_tapes rows (contig_uid)"] = diff(keyed(b0["tapes"], ("contig_uid",)), keyed(b1["tapes"], ("contig_uid",)))
T["tape_tokens rows (contig_uid)"] = diff(keyed(b0["tape_tokens"], ("contig_uid",)), keyed(b1["tape_tokens"], ("contig_uid",)))
T["token_dict (token_id)"] = diff(keyed(b0["token_dict"], ("token_id",)), keyed(b1["token_dict"], ("token_id",)))
for label, by_f in (("feature_uid key (registry)", True), ("ordinal key (pre-phase-B)", False)):
    s0, k0 = tape_slots(b0["tapes"], b0["tape_tokens"], by_f)
    s1, k1 = tape_slots(b1["tapes"], b1["tape_tokens"], by_f)
    # positional columns: delta and junction to the previous slot, distances to the contig ends, and the
    # ordinal itself when the key is the feature_uid
    T[f"tape slot entries, all list columns and the {'ordinal' if by_f else 'feature_uid'}, {label}"] = diff(s0, s1)
    if by_f:
        T[f"tape slot entries, all list columns without the ordinal, {label}"] = diff(
            {k: v[:-1] for k, v in s0.items()}, {k: v[:-1] for k, v in s1.items()})
    T[f"tape slot entries, own coordinates, kind and protein only, {label}"] = diff(
        {k: v[:4] + v[-2:-1] for k, v in s0.items()}, {k: v[:4] + v[-2:-1] for k, v in s1.items()})
    T[f"tape_tokens slot entries, {label}"] = diff({k: v[:-1] for k, v in k0.items()} if by_f else k0,
                                                  {k: v[:-1] for k, v in k1.items()} if by_f else k1)
pk = ("token_id", "contig_uid", "feature_uid", "seg_idx")
ok = ("token_id", "contig_uid", "ordinal", "seg_idx")
T["token_postings, primary key (token_id, contig_uid, feature_uid, seg_idx), all columns"] = diff(
    keyed(b0["token_postings"], pk), keyed(b1["token_postings"], pk))


def without(rs: list[dict], cols: tuple[str, ...]) -> list[dict]:
    return [{c: x for c, x in r.items() if c not in cols} for r in rs]


T["token_postings, primary key, without n_slots"] = diff(
    keyed(without(b0["token_postings"], ("n_slots",)), pk), keyed(without(b1["token_postings"], ("n_slots",)), pk))
T["token_postings, primary key, without ordinal and n_slots"] = diff(
    keyed(without(b0["token_postings"], ("ordinal", "n_slots")), pk), keyed(without(b1["token_postings"], ("ordinal", "n_slots")), pk))
T["token_postings, ordinal key (token_id, contig_uid, ordinal, seg_idx), all columns"] = diff(
    keyed(b0["token_postings"], ok), keyed(b1["token_postings"], ok))
T["token_postings, ordinal key, without n_slots"] = diff(
    keyed(without(b0["token_postings"], ("n_slots",)), ok), keyed(without(b1["token_postings"], ("n_slots",)), ok))
T["token_postings rows on the E. coli contig"] = sum(1 for r in b0["token_postings"] if r["contig_uid"] == ecoli)
print(json.dumps(out, indent=1, default=lambda x: x.hex() if isinstance(x, bytes) else str(x)))

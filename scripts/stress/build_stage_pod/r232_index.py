"""Card C-74 (DEC-59, DEC-79): an accession index of GTDB R232 and GlobDB r232, and the r226-to-r232 re-key.

    r232_index.py build GLOBDB_TABLES_DIR GTDB_META_DIR GLOBDB_R226_DIR OUT_DIR
    r232_index.py lookup INDEX_DIR ACCESSIONS_FILE        (one accession per line; TSV to stdout)
    r232_index.py rekey INDEX_DIR R226_IDS_FILE           (one GlobDB r226 id per line; TSV to stdout)
    r232_index.py unbuilt NON_GTDB_STAGES_TSV             (the r226 ids no stage holds yet, one per line)

Inputs (all fetched by card C-72 part 1 and relayed by card C-73):
  GLOBDB_TABLES_DIR  globdb_r232_taxonomy.tsv.gz (id, lineage; no header), globdb_r232_checkm2.tsv.gz,
                     globdb_r232_dictionaries.tar.gz, globdb_r232_dataset_list.tsv
  GTDB_META_DIR      ar53_metadata_r232.tsv.gz, bac120_metadata_r232.tsv.gz
  GLOBDB_R226_DIR    globdb_r226_taxonomy.tsv.gz, globdb_r226_dictionaries.tar.gz (card C-41's download)

Output (Parquet, each sorted by its key with min/max statistics, byte-identical when rebuilt: AGENTS 4.3):
  gtdb_r232.parquet         one row per GTDB R232 genome; key `accession`, the versionless GCA_/GCF_ accession
  globdb_r232.parquet       one row per GlobDB r232 genome; key `globdb_id`
  globdb_r226_keys.parquet  one row per GlobDB r226 genome; key `r226_id`, with the original identifier its
                            dictionary gives (the genome's name in the source catalogue)
  globdb_r232_dictionary.parquet  every row of r232's dictionaries (also the genomes dereplication left out)
  summary.json              the counts, and the resources whose count differs from the provider's dataset list

A GlobDB identifier is not stable across releases (r232 numbers some resources again), so the re-key goes through the
source catalogue's own identifier, `original_id`, never through the GlobDB id alone. Nothing is guessed: a genome that
no rule below can place is `dropped` with the reason, and an identifier r226 never had is `unknown`.
"""

from __future__ import annotations

import gzip
import io
import json
import re
import sys
import tarfile
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

ROW_GROUP = 122_880
ACC = re.compile(r"GC[AF]_\d{9}")
GART_PHYLA = ("Myxococcota", "Myxococcota_A", "Cyanobacteriota", "Omnitrophota", "Planctomycetota",
              "Patescibacteriota")
FILES = {
    "r232_tax": "globdb_r232_taxonomy.tsv.gz", "r232_checkm2": "globdb_r232_checkm2.tsv.gz",
    "r232_dict": "globdb_r232_dictionaries.tar.gz", "r232_list": "globdb_r232_dataset_list.tsv",
    "ar53": "ar53_metadata_r232.tsv.gz", "bac120": "bac120_metadata_r232.tsv.gz",
    "r226_tax": "globdb_r226_taxonomy.tsv.gz", "r226_dict": "globdb_r226_dictionaries.tar.gz",
}


class Refused(Exception):
    pass


def _need(d: Path, key: str) -> Path:
    p = Path(d) / FILES[key]
    if not p.is_file():
        raise Refused(f"missing input file {p}")
    return p


def nover(acc: str) -> str:
    """`RS_GCF_000000033.2`, `GCA_000000011.1` or `GCA_000000011` to `GCF_000000033`, `GCA_000000011`."""
    m = ACC.search(acc)
    return m.group(0) if m else acc


def normalize_original(s: str) -> str:
    for suf in (".fna", ".fa", ".fasta"):
        if s.endswith(suf):
            return s[: -len(suf)]
    return s


def _lineage_field(tax: str, prefix: str) -> str | None:
    for f in tax.split(";"):
        if f.startswith(prefix):
            return f[len(prefix):] or None
    return None


def _read_tax(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            a, _, b = line.rstrip("\n").partition("\t")
            if a:
                out[a] = b
    return out


def _read_dictionaries(path: Path) -> list[tuple[str, str, str | None, str | None]]:
    """(globdb_id, original_id (normalized), cluster, insdc) for every dictionary row, in member-name order."""
    rows: list[tuple[str, str, str | None, str | None]] = []
    with tarfile.open(path) as tf:
        members = sorted((m for m in tf.getmembers() if m.isfile() and m.name.endswith("_dictionary.tsv")),
                         key=lambda m: m.name.rsplit("/", 1)[-1])
        for m in members:
            lines = tf.extractfile(m).read().decode().splitlines()  # type: ignore[union-attr]
            id_first = bool(lines) and lines[0].startswith("GlobDB_ID\t")
            for line in lines[1:] if id_first else lines:
                p = line.split("\t")
                cluster = None
                if id_first:
                    gid = p[0]
                    if len(p) >= 3:
                        cluster, orig = p[1], p[2]
                    elif len(p) == 2 and " " in p[1]:
                        cluster, _, orig = p[1].partition(" ")
                    else:
                        continue
                else:
                    if len(p) < 2:
                        continue
                    orig, gid = p[0], p[1]
                if gid.endswith(".fa"):
                    gid = gid[:-3]
                orig = normalize_original(orig)
                m2 = ACC.search(orig)
                rows.append((gid, orig, cluster, m2.group(0) if m2 else None))
    return rows


def _prefixes(list_path: Path) -> tuple[list[tuple[str, str]], dict[str, int]]:
    pre: list[tuple[str, str]] = []
    listed: dict[str, int] = {}
    lines = Path(list_path).read_text().splitlines()
    for line in lines[1:]:
        p = line.split("\t")
        if len(p) < 6 or not p[0]:
            continue
        listed[p[0]] = int(p[5])
        for tok in p[1].split("&"):
            tok = tok.strip()
            pre.append((tok + "_" if tok in ("GCA", "GCF") else tok, p[0]))
    pre.sort(key=lambda t: (-len(t[0]), t[0]))
    return pre, listed


def _resource(gid: str, pre: list[tuple[str, str]]) -> str | None:
    for pfx, name in pre:
        if gid.startswith(pfx):
            return name
    return None


def _write(table: pa.Table, key: str, path: Path) -> None:
    table = table.take(pc.sort_indices(table, sort_keys=[(key, "ascending")]))
    vals = table.column(key).to_pylist()
    if len(set(vals)) != len(vals):
        raise Refused(f"{path.name}: duplicate keys in {key}")
    tmp = path.with_name("." + path.name + ".partial")
    pq.write_table(table, tmp, row_group_size=ROW_GROUP, compression="zstd", write_statistics=True)
    tmp.replace(path)


def _gtdb_rows(path: Path) -> list[tuple]:
    cols = ["accession", "gtdb_genome_representative", "gtdb_representative", "gtdb_taxonomy",
            "ncbi_genbank_assembly_accession"]
    rd = pacsv.open_csv(
        str(path), read_options=pacsv.ReadOptions(block_size=1 << 24),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(include_columns=cols, column_types={c: pa.string() for c in cols}))
    out = []
    for b in rd:
        lists = [b.column(i).to_pylist() for i in range(len(cols))]
        out.extend(zip(*lists))
    return out


def build(globdb_dir, gtdb_dir, r226_dir, out_dir) -> dict:
    gl, gt, r6, out = Path(globdb_dir), Path(gtdb_dir), Path(r226_dir), Path(out_dir)
    paths = {k: _need(d, k) for d, ks in ((gl, ("r232_tax", "r232_checkm2", "r232_dict", "r232_list")),
                                          (gt, ("ar53", "bac120")), (r6, ("r226_tax", "r226_dict"))) for k in ks}
    out.mkdir(parents=True, exist_ok=True)

    # GTDB R232
    rec = []
    for key in ("ar53", "bac120"):
        for acc, rep, isrep, tax, gb in _gtdb_rows(paths[key]):
            rec.append((nover(acc), acc, nover(gb) if gb and gb != "none" else nover(acc), acc[:2],
                        _lineage_field(tax, "d__"), tax, _lineage_field(tax, "s__"), isrep == "t", nover(rep)))
    schema = pa.schema([("accession", pa.string()), ("accession_versioned", pa.string()),
                        ("genbank_accession", pa.string()), ("db", pa.string()), ("domain", pa.string()),
                        ("taxonomy", pa.string()), ("species", pa.string()), ("is_representative", pa.bool_()),
                        ("representative", pa.string())])
    cols = list(zip(*rec)) if rec else [()] * 9
    _write(pa.table([pa.array(list(c), type=f.type) for c, f in zip(cols, schema)], schema=schema), "accession",
           out / "gtdb_r232.parquet")
    species = {r[6] for r in rec}
    reps = sum(1 for r in rec if r[7])

    # GlobDB r232
    tax = _read_tax(paths["r232_tax"])
    pre, listed = _prefixes(paths["r232_list"])
    ck: dict[str, tuple] = {}
    rd = pacsv.open_csv(str(paths["r232_checkm2"]), parse_options=pacsv.ParseOptions(delimiter="\t"),
                        convert_options=pacsv.ConvertOptions(
                            include_columns=["Name", "Completeness", "Contamination", "Genome_Size", "Total_Contigs"],
                            column_types={"Name": pa.string(), "Completeness": pa.float64(),
                                          "Contamination": pa.float64(), "Genome_Size": pa.int64(),
                                          "Total_Contigs": pa.int64()}))
    for b in rd:
        for n, c, k, s, t in zip(*(b.column(i).to_pylist() for i in range(5))):
            ck[n] = (c, k, s, t)
    d232: dict[str, tuple] = {}
    for gid, orig, cl, ins in _read_dictionaries(paths["r232_dict"]):
        d232[gid] = (orig, cl, ins)
    grec = []
    by_res: dict[str, int] = {}
    for gid in sorted(tax):
        res = _resource(gid, pre) or "unlisted"
        by_res[res] = by_res.get(res, 0) + 1
        is_gtdb = res == "GTDB"
        orig, cl, ins = d232.get(gid, (None, None, None))
        if is_gtdb:
            ins = nover(gid)
        c = ck.get(gid, (None, None, None, None))
        t = tax[gid]
        grec.append((gid, res, is_gtdb, _lineage_field(t, "d__"), _lineage_field(t, "p__"),
                     _lineage_field(t, "s__"), t, orig, cl, ins, *c))
    gs = pa.schema([("globdb_id", pa.string()), ("resource", pa.string()), ("is_gtdb", pa.bool_()),
                    ("domain", pa.string()), ("phylum", pa.string()), ("species", pa.string()),
                    ("taxonomy", pa.string()), ("original_id", pa.string()), ("cluster_id", pa.string()),
                    ("insdc_accession", pa.string()), ("completeness", pa.float64()),
                    ("contamination", pa.float64()), ("genome_size", pa.int64()), ("total_contigs", pa.int64())])
    gcols = list(zip(*grec)) if grec else [()] * 14
    _write(pa.table([pa.array(list(c), type=f.type) for c, f in zip(gcols, gs)], schema=gs), "globdb_id",
           out / "globdb_r232.parquet")
    ds = pa.schema([("globdb_id", pa.string()), ("original_id", pa.string()), ("cluster_id", pa.string())])
    dcols = list(zip(*[(g, v[0], v[1]) for g, v in d232.items()])) or [()] * 3
    _write(pa.table([pa.array(list(c), type=f.type) for c, f in zip(dcols, ds)], schema=ds), "globdb_id",
           out / "globdb_r232_dictionary.parquet")
    diffs = {}
    for name in sorted(set(listed) | set(by_res)):
        if name == "unlisted":
            continue
        if listed.get(name, 0) != by_res.get(name, 0):
            diffs[name] = {"listed": listed.get(name, 0), "found": by_res.get(name, 0)}
    if "unlisted" in by_res:
        diffs["unlisted"] = {"listed": 0, "found": by_res["unlisted"]}

    # GlobDB r226 keys
    t226 = _read_tax(paths["r226_tax"])
    d226 = {gid: (orig, cl, ins) for gid, orig, cl, ins in _read_dictionaries(paths["r226_dict"])}
    krec = []
    for gid in sorted(t226):
        orig, cl, ins = d226.get(gid, (None, None, None))
        krec.append((gid, orig, cl, ins))
    ks = pa.schema([("r226_id", pa.string()), ("original_id", pa.string()), ("cluster_id", pa.string()),
                    ("insdc_accession", pa.string())])
    kcols = list(zip(*krec)) if krec else [()] * 4
    _write(pa.table([pa.array(list(c), type=f.type) for c, f in zip(kcols, ks)], schema=ks), "r226_id",
           out / "globdb_r226_keys.parquet")

    summary = {
        "gtdb_genomes": len(rec), "gtdb_species": len(species), "gtdb_representatives": reps,
        "globdb_genomes": len(grec), "globdb_by_resource": dict(sorted(by_res.items())),
        "dataset_list_differences": diffs, "globdb_r226_genomes": len(krec),
        "resource_prefixes": [[p, n] for p, n in pre],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    return summary


# ---- queries ---------------------------------------------------------------------------------------------------

def _rows(index_dir, name: str) -> list[dict]:
    return pq.read_table(Path(index_dir) / name).to_pylist()


def _col(index_dir, name: str, col: str) -> list:
    return pq.read_table(Path(index_dir) / name, columns=[col]).column(col).to_pylist()


def lookup(index_dir, accessions: list[str]) -> list[dict]:
    """Per input accession (any form), in the input's order: membership of GTDB R232 and GlobDB r232."""
    gt = _rows(index_dir, "gtdb_r232.parquet")
    by_key: dict[str, dict] = {}
    for r in gt:
        by_key[r["accession"]] = r
    for r in gt:
        by_key.setdefault(r["genbank_accession"], r)
    gl = _rows(index_dir, "globdb_r232.parquet")
    by_id = {r["globdb_id"]: r["globdb_id"] for r in gl}
    by_ins: dict[str, str] = {}
    for r in gl:
        if r["insdc_accession"]:
            by_ins.setdefault(r["insdc_accession"], r["globdb_id"])
    out = []
    for a in accessions:
        k = nover(a.strip())
        g = by_key.get(k)
        keys = [k] + ([g["accession"], g["genbank_accession"]] if g else [])
        gid = next((by_id[x] for x in keys if x in by_id), None) or next(
            (by_ins[x] for x in keys if x in by_ins), None)
        out.append({"accession": a.strip(), "in_gtdb_r232": g is not None, "in_globdb_r232": gid is not None,
                    "is_representative": g["is_representative"] if g else None,
                    "representative": g["representative"] if g else None,
                    "species": g["species"] if g else None, "globdb_r232_id": gid})
    return out


def rekey(index_dir, r226_ids: list[str]) -> list[dict]:
    """For each r226 id: kept (its r232 id), replaced (the r232 id that now stands for it), dropped, or unknown."""
    gt = {r["accession"]: r for r in _rows(index_dir, "gtdb_r232.parquet")}
    gl = _rows(index_dir, "globdb_r232.parquet")
    ids232 = {r["globdb_id"] for r in gl}
    dict232 = set(_col(index_dir, "globdb_r232_dictionary.parquet", "original_id"))
    by_orig: dict[str, str] = {}
    by_cluster: dict[str, str] = {}
    for r in gl:
        if r["original_id"]:
            by_orig.setdefault(r["original_id"], r["globdb_id"])
        if r["cluster_id"]:
            by_cluster.setdefault(r["cluster_id"], r["globdb_id"])
    keys = {r["r226_id"]: r for r in _rows(index_dir, "globdb_r226_keys.parquet")}
    pre = json.loads((Path(index_dir) / "summary.json").read_text())["resource_prefixes"]
    prefixes = [p for p, _ in pre]
    out = []

    def gtdb_rep(acc: str | None) -> str | None:
        if not acc:
            return None
        g = gt.get(acc) or next((v for v in gt.values() if v["genbank_accession"] == acc), None)
        if g is None:
            return None
        rep = g["representative"]
        return rep if rep in ids232 else None

    for rid in r226_ids:
        k = keys.get(rid)
        res = {"r226_id": rid, "status": "unknown", "r232_id": None, "reason": "not_a_globdb_r226_genome"}
        if k is None:
            out.append(res)
            continue
        orig, cl, ins = k["original_id"], k["cluster_id"], k["insdc_accession"]
        gtdb_form = bool(ACC.fullmatch(rid))
        if gtdb_form and not orig:
            ins = rid
        new = None
        if not gtdb_form and not any(rid.startswith(p) for p in prefixes):
            res.update(status="dropped", reason="resource_not_in_r232")
        elif gtdb_form and rid in ids232 and not orig:
            res.update(status="kept", r232_id=rid, reason="same_id")
        elif orig and by_orig.get(orig) in ids232:
            res.update(status="kept", r232_id=by_orig[orig], reason="same_original_id")
        elif not orig and rid in ids232 and not gtdb_form:
            res.update(status="kept", r232_id=rid, reason="same_id_no_dictionary")
        elif (new := gtdb_rep(ins)) is not None:
            res.update(status="replaced", r232_id=new,
                       reason="accession_in_gtdb_r232" if gtdb_form else "insdc_accession_in_gtdb_r232")
        elif cl and cl in by_cluster:
            res.update(status="replaced", r232_id=by_cluster[cl], reason="same_cluster_id")
        elif orig:
            res.update(status="dropped", reason="dictionary_lists_it_but_not_in_r232" if orig in dict232
                       else "original_not_in_r232")
        else:
            res.update(status="dropped", reason="accession_not_in_r232" if gtdb_form else "not_in_r232")
        out.append(res)
    return out


def unbuilt_ids(non_gtdb_stages_tsv) -> list[str]:
    """The r226 non-GTDB ids that no stage holds: not in G1 and not in a GART1/GART2 phylum (card C-63)."""
    out = []
    with open(non_gtdb_stages_tsv) as f:
        head = f.readline().rstrip("\n").split("\t")
        i_id, i_stage, i_tax = head.index("ID"), head.index("stage"), head.index("taxonomy")
        for line in f:
            p = line.rstrip("\n").split("\t")
            phylum = _lineage_field(p[i_tax], "p__")
            if p[i_stage] != "G1" and phylum not in GART_PHYLA:
                out.append(p[i_id])
    return sorted(out)


def _tsv(rows: list[dict], cols: list[str]) -> str:
    buf = io.StringIO()
    buf.write("\t".join(cols) + "\n")
    for r in rows:
        buf.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
    return buf.getvalue()


def main(argv: list[str]) -> int:
    try:
        if len(argv) == 6 and argv[1] == "build":
            s = build(*argv[2:6])
            print(json.dumps({k: v for k, v in s.items() if k != "resource_prefixes"}, indent=1, sort_keys=True))
            return 0
        if len(argv) == 4 and argv[1] == "lookup":
            ids = [x for x in Path(argv[3]).read_text().split() if x]
            sys.stdout.write(_tsv(lookup(argv[2], ids), ["accession", "in_gtdb_r232", "in_globdb_r232",
                                                         "is_representative", "representative", "species",
                                                         "globdb_r232_id"]))
            return 0
        if len(argv) == 4 and argv[1] == "rekey":
            ids = [x for x in Path(argv[3]).read_text().split() if x]
            sys.stdout.write(_tsv(rekey(argv[2], ids), ["r226_id", "status", "r232_id", "reason"]))
            return 0
        if len(argv) == 3 and argv[1] == "unbuilt":
            sys.stdout.write("".join(i + "\n" for i in unbuilt_ids(argv[2])))
            return 0
        print(__doc__, file=sys.stderr)
        return 2
    except Refused as e:
        print(f"r232_index: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))

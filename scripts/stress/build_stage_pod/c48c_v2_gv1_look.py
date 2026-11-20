"""Card C-48c (part 2): read-only looks into the archives behind stages V2 and GV1, before anything is decided.

    python c48c_v2_gv1_look.py head ARCHIVE MEMBER [N]        first N lines (default 20) of a member
    python c48c_v2_gv1_look.py headers ARCHIVE MEMBER [N]     the first N FASTA headers of a member
    python c48c_v2_gv1_look.py members ARCHIVE MEMBER [N]     the first N members of a tar(.gz) or zip member
    python c48c_v2_gv1_look.py dirs ARCHIVE MEMBER            files and bytes per directory of a tar(.gz) member
    python c48c_v2_gv1_look.py gvmags GVMAGS_TAR OUT          per genome of both FASTA sets: records, bases, ...
    python c48c_v2_gv1_look.py inovirus META_ZIP OUT          per GenBank record of the Inovirus set
    python c48c_v2_gv1_look.py imgvr HC_ZIP OUT               the HC tables copied out (md5 checked) and summarised
    python c48c_v2_gv1_look.py gvsubset GVMAGS_TAR GV1_INPUTS OUT  fna_files.tar.gz against GV1's staged genomes

ARCHIVE is a .zip or a .tar; MEMBER a path inside it, optionally `OUTER::INNER` for a member of an inner .zip or
.tar(.gz). A member ending in .gz is decompressed on the fly. Nothing is unpacked but the two IMG/VR tables (`imgvr`
copies them to OUT so DuckDB can read them; 2.2 GB); everything else is streamed one record or one member at a time
(AGENTS §4.1). OUT gets a TSV per genome or record and summary.json.

gvmags: for each `.fna` member of all_18kgenomes_fna.tar.gz and fna_files.tar.gz, its set, name, FASTA records,
bases, N bases, the sha256 of its bytes, its genome uid (`ids.genome_uid` of its distinct contig uids), its first
header; any other member is printed (it is data: never run). Then which genomes are in both sets and whether their
bytes agree, and the genomes long_ortho_8000.tsv names.

inovirus: every GenBank record of Gb_files_inoviruses.zip (inside the metadata zip), with its directory (family and
subfamily), file, LOCUS name, topology, length, contig uid and definition.

imgvr: IMGVR_all_Sequence_information-high_confidence.tsv and IMGVR_all_Host_information-high_confidence.tsv copied
to OUT (temporary name, then renamed), their md5 compared with sources.yaml's (refused, exit 2, when one differs),
and counted with DuckDB (2 threads, 4 GB): rows, UViGs, vOTUs, by confidence, MIUViG quality, topology and origin;
vOTUs with a high-quality member; the high-quality members' completeness, contamination and lengths.
"""

from __future__ import annotations

import contextlib
import csv
import gzip
import hashlib
import io
import json
import shutil
import sys
import tarfile
import zipfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import c48c_virus_inputs as vi  # noqa: E402
from biobox import ids  # noqa: E402

HC = "Cus_VR/IMG_VR_2022-12-19_7.1/"
HC_TABLES = {"IMGVR_all_Sequence_information-high_confidence.tsv": "3c516db128082fa29dc2c2f60520da1b",
             "IMGVR_all_Host_information-high_confidence.tsv": "71b54d0f5c186d813f058bf0379dfd24"}  # sources.yaml
INOVIRUS_ZIP = "Cusrus/data/Gb_files_inoviruses.zip"
GV_SETS = ("gvmags_v2/all_18kgenomes_fna.tar.gz", "gvmags_v2/fna_files.tar.gz")
IMGVR_COLS = ["uvig", "taxon_oid", "scaffold_oid", "coords", "ecosystem", "votu", "length", "topology",
              "genomad_score", "confidence", "completeness", "contamination", "miuvig", "gene_content", "taxonomy",
              "tax_method", "host_tax", "host_method", "origin"]


@contextlib.contextmanager
def open_member(archive: Path, member: str):
    """A binary stream of `member` of `archive` (zip or tar), `OUTER::INNER` for a nested member, .gz decompressed."""
    outer, _, inner = member.partition("::")
    with contextlib.ExitStack() as st:
        if zipfile.is_zipfile(archive):
            zf = st.enter_context(zipfile.ZipFile(archive))
            fh = st.enter_context(zf.open(outer))
        else:
            tf = st.enter_context(tarfile.open(archive, mode="r:"))
            fh = tf.extractfile(outer)
            if fh is None:
                raise SystemExit(f"{outer}: not a file in {archive}")
        if inner:
            if outer.endswith(".zip"):
                zf2 = st.enter_context(zipfile.ZipFile(fh))
                fh = st.enter_context(zf2.open(inner))
            else:
                tf2 = st.enter_context(tarfile.open(fileobj=fh, mode="r|*"))
                for m in tf2:
                    if m.name == inner:
                        fh = tf2.extractfile(m)
                        break
                else:
                    raise SystemExit(f"{inner}: not in {outer}")
            name = inner
        else:
            name = outer
        if name.endswith(".gz") and not name.endswith(".tar.gz"):
            fh = st.enter_context(gzip.GzipFile(fileobj=fh))
        yield fh


def write_tsv(path: Path, head: list[str], rows) -> None:
    tmp = path.with_name(f".{path.name}.tmp")
    with open(tmp, "w") as w:
        w.write("\t".join(head) + "\n")
        for r in rows:
            w.write("\t".join(map(str, r)) + "\n")
    tmp.replace(path)


def gvmags(tar: Path, out: Path) -> dict:
    rows, other = [], []
    for setname in GV_SETS:
        with open_member(tar, setname) as fh, tarfile.open(fileobj=fh, mode="r|gz") as tf:
            for m in tf:
                if not m.isfile():
                    continue
                if not m.name.endswith(".fna"):
                    data = tf.extractfile(m).read(4000)
                    other.append({"set": setname, "member": m.name, "bytes": m.size,
                                  "head": data.decode("utf-8", "replace")})
                    continue
                raw = tf.extractfile(m).read()  # one genome (at most a few MB) at a time
                n = bases = nn = 0
                uids, first = [], ""
                lens = []
                for r in vi.iter_fasta_stream(io.BytesIO(raw)):
                    n += 1
                    bases += len(r.seq)
                    nn += r.seq.upper().count(b"N")
                    uids.append(ids.contig_uid(r.seq))
                    lens.append(len(r.seq))
                    first = first or r.header
                rows.append((setname.split("/")[-1].split(".")[0], m.name.split("/")[-1][:-4], m.name, n, bases, nn,
                             max(lens, default=0), hashlib.sha256(raw).hexdigest(),
                             ids.genome_uid(list(dict.fromkeys(uids))).hex() if uids else "", first[:200]))
    rows.sort(key=lambda r: (r[0], r[1]))
    write_tsv(out / "gvmags_genomes.tsv", ["set", "genome", "member", "records", "bases", "n_bases", "longest",
                                           "sha256", "genome_uid", "first_header"], rows)
    with open_member(tar, "gvmags_v2/long_ortho_8000.tsv") as fh:
        ortho = Counter(line.split(b"\t")[1].decode() for k, line in enumerate(fh) if k and line.strip())
    by = {s: {r[1]: r for r in rows if r[0] == s} for s in {r[0] for r in rows}}
    a, b = by.get("all_18kgenomes_fna", {}), by.get("fna_files", {})
    both = sorted(set(a) & set(b))
    rep = {
        "sets": {s: {"genomes": len(v), "records": sum(r[3] for r in v.values()),
                     "bases": sum(r[4] for r in v.values()), "n_bases": sum(r[5] for r in v.values()),
                     "single_record_genomes": sum(1 for r in v.values() if r[3] == 1),
                     "records_per_genome_median": sorted(r[3] for r in v.values())[len(v) // 2] if v else 0,
                     "bases_per_genome_median": sorted(r[4] for r in v.values())[len(v) // 2] if v else 0,
                     "genomes_lt_100kb": sum(1 for r in v.values() if r[4] < 100_000),
                     "distinct_genome_uids": len({r[8] for r in v.values()})}
                 for s, v in sorted(by.items())},
        "names_in_both": len(both), "same_bytes_in_both": sum(1 for g in both if a[g][7] == b[g][7]),
        "same_genome_uid_in_both": sum(1 for g in both if a[g][8] == b[g][8]),
        "fna_files_not_in_18k": sorted(set(b) - set(a))[:50], "fna_files_not_in_18k_n": len(set(b) - set(a)),
        "differing_in_both": {
            "n": sum(1 for g in both if a[g][8] != b[g][8]),
            "fewer_records_in_fna_files": sum(1 for g in both if a[g][8] != b[g][8] and b[g][3] < a[g][3]),
            "more_records_in_fna_files": sum(1 for g in both if a[g][8] != b[g][8] and b[g][3] > a[g][3]),
            "same_records_in_fna_files": sum(1 for g in both if a[g][8] != b[g][8] and b[g][3] == a[g][3]),
            "same_bases": sum(1 for g in both if a[g][8] != b[g][8] and b[g][4] == a[g][4]),
            "bases_18k": sum(a[g][4] for g in both if a[g][8] != b[g][8]),
            "bases_fna_files": sum(b[g][4] for g in both if a[g][8] != b[g][8]),
            "examples": [(g, a[g][3], a[g][4], b[g][3], b[g][4], a[g][9][:80], b[g][9][:80])
                         for g in both if a[g][8] != b[g][8]][:8]},
        "genome_uid_repeats_18k": {"genomes": len(a) - len({r[8] for r in a.values()}),
                                   "examples": [sorted(g for g in a if a[g][8] == u)[:4] for u, n in Counter(
                                       r[8] for r in a.values()).most_common(5) if n > 1]},
        "header_examples_18k": [a[g][9] for g in sorted(a)[:: max(1, len(a) // 12)]],
        "ortho_genomes": len(ortho), "ortho_rows": sum(ortho.values()),
        "ortho_genomes_in_18k": sum(1 for g in ortho if g in a), "ortho_genomes_in_fna_files": sum(1 for g in ortho if g in b),
        "other_members": other,
    }
    return rep


def gvsubset(tar: Path, gv1: Path, out: Path) -> dict:
    """What fna_files.tar.gz is relative to GV1's staged genomes (card C-48c part 3): per genome of fna_files, its
    contig uids against the same-named GV1 genome's (GV1/records.tsv): the same set, a subset (its contigs all in the
    GV1 genome, fewer of them), or other sequences (a contig the GV1 genome lacks; then whether GV1 has it anywhere)."""
    per: dict[str, set[str]] = {}
    lens: dict[str, int] = {}
    with open(gv1 / "records.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            per.setdefault(r["genome"], set()).add(r["contig_uid"])
            lens[r["contig_uid"]] = int(r["length"])
    rows, kinds = [], Counter()
    with open_member(tar, GV_SETS[1]) as fh, tarfile.open(fileobj=fh, mode="r|gz") as tf:
        for m in tf:
            if not m.isfile() or not m.name.endswith(".fna"):
                continue
            g = m.name.split("/")[-1][:-4]
            mine = [ids.contig_uid(r.seq).hex() for r in vi.iter_fasta_stream(io.BytesIO(tf.extractfile(m).read()))]
            s, gv = set(mine), per.get(g)
            if gv is None:
                kind = "not_in_gv1"
            elif s == gv:
                kind = "same_contigs"
            elif s < gv:
                kind = "subset_of_gv1_genome"
            else:
                kind = "other_sequences"
            kinds[kind] += 1
            missing = sorted(gv - s) if gv else []
            rows.append((g, kind, len(s), len(gv or ()), len(missing), sum(lens[u] for u in missing),
                         sum(1 for u in s if u not in lens)))
    rows.sort()
    write_tsv(out / "fna_files_vs_gv1.tsv", ["genome", "relation", "contigs", "gv1_contigs", "gv1_contigs_not_in_it",
                                             "gv1_bases_not_in_it", "contigs_gv1_lacks"], rows)
    dropped = [r for r in rows if r[4]]
    return {"genomes": len(rows), "relation": dict(sorted(kinds.items())),
            "contigs_gv1_lacks": sum(r[6] for r in rows),
            "genomes_with_gv1_contigs_dropped": len(dropped),
            "gv1_contigs_dropped": sum(r[4] for r in dropped), "gv1_bases_dropped": sum(r[5] for r in dropped),
            "dropped_contigs_per_genome": dict(sorted(Counter(r[4] for r in dropped).items())),
            "dropped_bases_median": sorted(r[5] for r in dropped)[len(dropped) // 2] if dropped else 0}


def inovirus(meta_zip: Path, out: Path) -> dict:
    rows = []
    with zipfile.ZipFile(meta_zip) as z0:
        inner = z0.read(INOVIRUS_ZIP)  # 40 MB: held, since seeking in a deflated member re-inflates it from the start
    with zipfile.ZipFile(io.BytesIO(inner)) as z:
        for info in sorted(z.infolist(), key=lambda i: i.filename):
            if info.is_dir():
                continue
            with z.open(info) as fh:
                for k, r in enumerate(vi.iter_genbank(fh)):
                    parts = info.filename.split("/")
                    rows.append((info.filename, k, parts[0], parts[1] if len(parts) > 2 else "", r["locus"],
                                 r["accession"], r["topology"], r["locus_length"], len(r["seq"]),
                                 ids.contig_uid(r["seq"]).hex(), r["definition"][:200].replace("\t", " ")))
    write_tsv(out / "inovirus_records.tsv", ["file", "record", "family", "subfamily", "locus", "accession",
                                             "topology", "locus_length", "length", "contig_uid", "definition"], rows)
    files = Counter(r[0] for r in rows)
    uids = Counter(r[9] for r in rows)
    return {"files": len(files), "records": len(rows), "records_per_file": dict(Counter(files.values())),
            "bases": sum(r[8] for r in rows), "distinct_contigs": len(uids),
            "distinct_bases": sum({r[9]: r[8] for r in rows}.values()),
            "length_disagrees_with_locus": sum(1 for r in rows if r[7] != r[8]),
            "topology": dict(Counter(r[6] for r in rows)), "families": dict(Counter(r[2] for r in rows)),
            "length_min": min((r[8] for r in rows), default=0), "length_max": max((r[8] for r in rows), default=0),
            "first_rows": rows[:3]}


def imgvr(zip_path: Path, out: Path) -> dict:
    import duckdb

    rep: dict = {"md5": {}}
    with zipfile.ZipFile(zip_path) as z:
        for name, want in HC_TABLES.items():
            dst = out / name
            tmp = out / f".{name}.tmp"
            md5 = hashlib.md5()
            with z.open(HC + name) as r, open(tmp, "wb") as w:
                for blk in iter(lambda: r.read(1 << 22), b""):
                    md5.update(blk)
                    w.write(blk)
            rep["md5"][name] = {"got": md5.hexdigest(), "sources_yaml": want, "equal": md5.hexdigest() == want}
            if md5.hexdigest() != want:
                tmp.unlink()
                raise SystemExit(2)
            tmp.replace(dst)
    con = duckdb.connect(config={"threads": 2, "memory_limit": "4GB", "temp_directory": str(out / "duckdb")})
    cols = "{" + ", ".join(f"'{c}': 'VARCHAR'" for c in IMGVR_COLS) + "}"
    con.execute(f"CREATE VIEW s AS SELECT * FROM read_csv('{out / 'IMGVR_all_Sequence_information-high_confidence.tsv'}',"
                f" delim='\t', header=true, quote='', escape='', columns={cols})")
    q = lambda sql: con.execute(sql).fetchall()  # noqa: E731
    (rows, uvigs, votus, bases), = q("SELECT count(*), count(DISTINCT uvig), count(DISTINCT votu), "
                                     "sum(length::BIGINT) FROM s")
    rep.update({"rows": rows, "uvigs": uvigs, "votus": votus, "bases": int(bases)})
    for col in ("confidence", "miuvig", "topology", "tax_method", "host_method"):
        rep[f"by_{col}"] = {str(k): [n, int(b)] for k, n, b in q(
            f"SELECT {col}, count(*), sum(length::BIGINT) FROM s GROUP BY 1 ORDER BY 2 DESC")}
    rep["by_origin"] = dict(q("SELECT origin, count(*) FROM s GROUP BY 1 ORDER BY 2 DESC LIMIT 40"))
    rep["by_ecosystem_l1"] = dict(q("SELECT split_part(ecosystem, ';', 1), count(*) FROM s GROUP BY 1 ORDER BY 2 DESC"))
    rep["votus_with_hq_member"], = q("SELECT count(DISTINCT votu) FROM s WHERE miuvig = 'High-quality'")[0]
    rep["votus_with_hq_or_reference_member"], = q(
        "SELECT count(DISTINCT votu) FROM s WHERE miuvig IN ('High-quality', 'Reference')")[0]
    rep["votus_with_reference_member"], = q("SELECT count(DISTINCT votu) FROM s WHERE miuvig = 'Reference'")[0]
    rep["hq_members"] = dict(zip(["n", "bases", "completeness_min", "completeness_max", "contamination_max",
                                  "length_min", "length_max", "at_100_complete"], q(
        "SELECT count(*), sum(length::BIGINT), min(completeness::DOUBLE), max(completeness::DOUBLE), "
        "max(contamination::DOUBLE), min(length::BIGINT), max(length::BIGINT), "
        "count(*) FILTER (WHERE completeness::DOUBLE >= 100) FROM s WHERE miuvig = 'High-quality'")[0]))
    rep["hq_members_by_topology"] = dict(q("SELECT topology, count(*) FROM s WHERE miuvig = 'High-quality' "
                                           "GROUP BY 1 ORDER BY 2 DESC"))
    rep["hq_members_per_votu"] = dict(q(
        "SELECT n, count(*) FROM (SELECT votu, count(*) n FROM s WHERE miuvig = 'High-quality' GROUP BY 1) "
        "GROUP BY 1 ORDER BY 1 LIMIT 30"))
    rep["members_per_hq_votu"] = dict(q(
        "SELECT CASE WHEN n >= 10 THEN 10 ELSE n END, count(*) FROM (SELECT votu, count(*) n FROM s WHERE votu IN "
        "(SELECT votu FROM s WHERE miuvig = 'High-quality') GROUP BY 1) GROUP BY 1 ORDER BY 1"))
    rep["gvmag_rows"] = q("SELECT uvig, scaffold_oid, coords, length, miuvig, completeness FROM s "
                          "WHERE topology = 'GVMAG' LIMIT 5")
    rep["completeness_unparseable"], = q("SELECT count(*) FROM s WHERE completeness IS NOT NULL AND "
                                         "TRY_CAST(completeness AS DOUBLE) IS NULL")[0]
    rep["completeness_values_unparseable"] = q("SELECT completeness, count(*) FROM s WHERE completeness IS NOT NULL "
                                               "AND TRY_CAST(completeness AS DOUBLE) IS NULL GROUP BY 1 LIMIT 10")
    rep["coords_whole"], = q("SELECT count(*) FROM s WHERE coords = 'whole'")[0]
    rep["uvig_repeats"], = q("SELECT count(*) - count(DISTINCT uvig) FROM s")[0]
    rep["host_rows"], = q(f"SELECT count(*) FROM read_csv('{out / 'IMGVR_all_Host_information-high_confidence.tsv'}', "
                          "delim='\t', header=true, quote='', escape='', all_varchar=true)")[0]
    rep["host_columns"] = [r[0] for r in q(
        f"DESCRIBE SELECT * FROM read_csv('{out / 'IMGVR_all_Host_information-high_confidence.tsv'}', "
        "delim='\t', header=true, quote='', escape='', all_varchar=true)")]
    shutil.rmtree(out / "duckdb", ignore_errors=True)
    return rep


def main(argv: list[str]) -> int:
    modes = ("head", "headers", "members", "dirs", "gvmags", "inovirus", "imgvr", "gvsubset")
    if len(argv) < 3 or argv[0] not in modes:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    mode = argv[0]
    if mode == "gvsubset" and len(argv) == 4:
        out = Path(argv[3])
        out.mkdir(parents=True, exist_ok=True)
        rep = gvsubset(Path(argv[1]), Path(argv[2]), out)
        (out / "gvsubset_summary.json").write_text(json.dumps(rep, indent=1, sort_keys=True))
        print(json.dumps(rep, sort_keys=True))
        return 0
    if mode in ("gvmags", "inovirus", "imgvr"):
        out = Path(argv[2])
        out.mkdir(parents=True, exist_ok=True)
        rep = {"gvmags": gvmags, "inovirus": inovirus, "imgvr": imgvr}[mode](Path(argv[1]), out)
        (out / f"{mode}_summary.json").write_text(json.dumps(rep, indent=1, sort_keys=True, default=str))
        print(json.dumps(rep, sort_keys=True, default=str)[:20000])
        return 0
    archive, member = Path(argv[1]), argv[2]
    n = int(argv[3]) if len(argv) > 3 else 20
    with open_member(archive, member) as fh:
        if mode == "members" and member.endswith(".zip"):
            with zipfile.ZipFile(fh) as zf:
                infos = zf.infolist()
                print(len(infos), "members,", sum(i.file_size for i in infos), "bytes uncompressed")
                for i in infos[:n]:
                    print(i.file_size, i.filename)
            return 0
        if mode == "dirs":  # files and bytes per directory (all but the last path component), and file suffixes
            per: dict[str, list[int]] = {}
            suffix: dict[str, int] = {}
            with tarfile.open(fileobj=fh, mode="r|*") as tf:
                for m in tf:
                    if not m.isfile():
                        continue
                    d = m.name.rsplit("/", 1)[0] if "/" in m.name else "."
                    v = per.setdefault(d, [0, 0])
                    v[0] += 1
                    v[1] += m.size
                    s = m.name.rsplit(".", 1)[-1] if "." in m.name.rsplit("/", 1)[-1] else ""
                    suffix[s] = suffix.get(s, 0) + 1
            for d in sorted(per):
                print(per[d][0], per[d][1], d)
            print("suffixes", sorted(suffix.items()))
            return 0
        if mode == "members":
            with tarfile.open(fileobj=fh, mode="r|*") as tf:
                for k, m in enumerate(tf):
                    if k >= n:
                        break
                    print(m.size, m.name)
            return 0
        k = 0
        for line in fh:
            if mode == "headers" and not line.startswith(b">"):
                continue
            sys.stdout.write(line.decode("utf-8", "replace")[:2000])
            k += 1
            if k >= n:
                break
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

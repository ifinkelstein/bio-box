#!/usr/bin/env python3
"""Card A-26: finalize /tmp stage inputs and compute exact linear r26.11.4 + GV1 prefix ownership.

Subcommands:
  finalize-inputs  Copy a stage's FASTA from /stor staging and metadata from /tmp into
                   /tmp/if2453/s2/inputs/<STAGE>/ with a /tmp free-space guard (>= 355 GiB),
                   atomic dot-temp rename, sha256 verification, and objects.json path rewrite.
  chain-ownership  Compute exact DEC-26 prefix contig_uid ownership across stored root stages
                   plus the 19-stage r26.11.4 linear chain and deferred GV1 in bounded memory
                   (streaming Parquet parts + DuckDB).
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

R26_11_4_CHAIN: tuple[tuple[str, int | None], ...] = (
    ("V2", None),
    ("PL2", None),
    ("PL3", None),
    ("CG1", None),
    ("MGE1", None),
    ("SAT1", None),
    ("ECE1", None),
    ("O1", None),
    ("O2", None),
    ("LPL1", None),
    ("SYM1", None),
    ("LRV1", None),
    ("LGRT1", None),
    ("MGE2", None),
    ("MV1", None),
    ("V3", None),
    ("LR1", 5000),
    ("GART1", None),
    ("GART2", None),
)
DEFERRED_STAGES: tuple[tuple[str, int | None], ...] = (("GV1", None),)

_FLUSH = 200_000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(4 * 1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def rewrite_objects_json(src: Path, dst: Path, inputs_dir: Path) -> str:
    """Copy objects.json to dst, rewriting each object path into inputs_dir/<basename>.

    If all paths in src already point to inputs_dir/<basename>, copy src byte-identically
    so its existing sha256 is preserved. Returns the sha256 hex digest of dst.
    """
    raw = src.read_bytes()
    objs = json.loads(raw.decode("utf-8"))
    changed = False
    for obj in objs:
        target = str(inputs_dir / Path(obj["path"]).name)
        if obj["path"] != target:
            obj["path"] = target
            changed = True
    tmp = dst.parent / f".{dst.name}.tmp.{os.getpid()}"
    if changed:
        tmp.write_text(json.dumps(objs, indent=1) + "\n", encoding="utf-8")
    else:
        tmp.write_bytes(raw)
    os.replace(tmp, dst)
    return sha256_file(dst)


def _free_gib(path: Path) -> int:
    probe = path if path.exists() else path.parent
    return shutil.disk_usage(probe).free // (1024 ** 3)


def _copy_atomic(src: Path, dst: Path) -> None:
    tmp = dst.parent / f".{dst.name}.tmp.{os.getpid()}"
    try:
        with src.open("rb") as fin, tmp.open("wb") as fout:
            shutil.copyfileobj(fin, fout, length=8 * 1024 * 1024)
            fout.flush()
            os.fsync(fout.fileno())
        os.replace(tmp, dst)
    except BaseException:
        if tmp.exists():
            tmp.unlink()
        raise


def finalize_stage_inputs(
    stage: str,
    meta_dir: Path,
    fasta_src: Path,
    fasta_sha256: str,
    dst_dir: Path,
    *,
    min_free_gib: int = 355,
    poll_seconds: int = 60,
    max_wait_seconds: int = 3000,
) -> dict:
    """Populate dst_dir from meta_dir + fasta_src while keeping /tmp >= min_free_gib free."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    fasta_dst = dst_dir / fasta_src.name
    need_copy = True
    if fasta_dst.is_file() and fasta_dst.stat().st_size == fasta_src.stat().st_size:
        if sha256_file(fasta_dst) == fasta_sha256:
            need_copy = False
    if need_copy:
        fasta_gib = (fasta_src.stat().st_size + (1024 ** 3) - 1) // (1024 ** 3)
        required_gib = min_free_gib + fasta_gib
        deadline = time.monotonic() + max_wait_seconds
        while True:
            avail = _free_gib(dst_dir)
            if avail >= required_gib:
                break
            if time.monotonic() >= deadline:
                raise RuntimeError(
                    f"{stage}: /tmp has {avail} GiB free, need >= {required_gib} GiB "
                    f"(min_free_gib={min_free_gib} + fasta_gib={fasta_gib}) after waiting {max_wait_seconds}s"
                )
            print(
                f"[a26] {stage}: waiting for /tmp headroom ({avail} GiB free < {required_gib} GiB required)",
                file=sys.stderr,
                flush=True,
            )
            time.sleep(poll_seconds)
        _copy_atomic(fasta_src, fasta_dst)
        got_sha = sha256_file(fasta_dst)
        if got_sha != fasta_sha256:
            fasta_dst.unlink(missing_ok=True)
            raise ValueError(f"{stage}: {fasta_dst.name} sha256 mismatch: got {got_sha}, expected {fasta_sha256}")

    for src in sorted(meta_dir.iterdir()):
        if not src.is_file():
            continue
        if src.name == "objects_staging.json":
            continue
        dst = dst_dir / src.name
        if src.name == "objects.json":
            rewrite_objects_json(src, dst, dst_dir)
        else:
            if not (dst.is_file() and dst.stat().st_size == src.stat().st_size and sha256_file(dst) == sha256_file(src)):
                _copy_atomic(src, dst)

    files = {}
    for p in sorted(dst_dir.iterdir()):
        if p.is_file() and not p.name.startswith("."):
            files[p.name] = {"bytes": p.stat().st_size, "sha256": sha256_file(p)}
    return {
        "stage": stage,
        "dst_dir": str(dst_dir),
        "free_gib_after": _free_gib(dst_dir),
        "files": files,
    }


def _lit(files: list[Path]) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def _write_rows(rows: list[tuple], schema, path: Path) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    cols = list(zip(*rows)) if rows else [[] for _ in schema.names]
    pq.write_table(pa.table({n: list(c) for n, c in zip(schema.names, cols)}, schema=schema), path)


def _extract_stage_parquet(stage: str, inputs_dir: Path, min_contig_len: int | None, out_dir: Path) -> tuple[str, list[Path], list[Path]]:
    """Stream one stage's genomes.json contig_uids and FASTA record lengths into Parquet parts."""
    import pyarrow as pa

    from biobox import fasta, ids

    sdir = out_dir / stage
    sdir.mkdir(parents=True, exist_ok=True)

    gschema = pa.schema([("contig_uid", pa.binary(16))])
    gparts: list[Path] = []
    seen_g: set[bytes] = set()
    rows: list[tuple] = []
    for g in json.loads((inputs_dir / "genomes.json").read_text(encoding="utf-8")):
        for h in g["contig_uids"]:
            u = bytes.fromhex(h)
            if u not in seen_g:
                seen_g.add(u)
                rows.append((u,))
                if len(rows) >= _FLUSH:
                    p = sdir / f"g-{len(gparts):05d}.parquet"
                    _write_rows(rows, gschema, p)
                    gparts.append(p)
                    rows = []
    p = sdir / f"g-{len(gparts):05d}.parquet"
    _write_rows(rows, gschema, p)
    gparts.append(p)
    del seen_g, rows

    lschema = pa.schema([("contig_uid", pa.binary(16)), ("n", pa.int64())])
    lparts: list[Path] = []
    lrows: list[tuple] = []
    objs = json.loads((inputs_dir / "objects.json").read_text(encoding="utf-8"))
    for obj in objs:
        for rec in fasta.iter_records(Path(obj["path"])):
            n = len(rec.seq)
            if min_contig_len is not None and n < min_contig_len:
                continue
            lrows.append((ids.contig_uid(rec.seq), n))
            if len(lrows) >= _FLUSH:
                lp = sdir / f"l-{len(lparts):05d}.parquet"
                _write_rows(lrows, lschema, lp)
                lparts.append(lp)
                lrows = []
    lp = sdir / f"l-{len(lparts):05d}.parquet"
    _write_rows(lrows, lschema, lp)
    lparts.append(lp)
    return stage, gparts, lparts


def compute_chain_ownership(
    root: Path,
    stages: list[tuple[str, Path, int | None]],
    tmp_dir: Path,
    *,
    workers: int = 1,
) -> dict[str, dict]:
    """Compute exact DEC-26 prefix ownership across root stages + ordered staged inputs."""
    import duckdb

    from biobox import io
    from biobox.l0seq.store import INDEX_DIR

    tmp_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="a26-chain-", dir=tmp_dir))
    try:
        extracted: dict[str, tuple[list[Path], list[Path]]] = {}
        if workers > 1 and len(stages) > 1:
            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as ex:
                futs = {
                    ex.submit(_extract_stage_parquet, st, inp, mlen, work / "parts"): st
                    for st, inp, mlen in stages
                }
                for fut in concurrent.futures.as_completed(futs):
                    st, gparts, lparts = fut.result()
                    extracted[st] = (gparts, lparts)
        else:
            for st, inp, mlen in stages:
                _, gparts, lparts = _extract_stage_parquet(st, inp, mlen, work / "parts")
                extracted[st] = (gparts, lparts)

        duck_tmp = work / "duck"
        duck_tmp.mkdir(parents=True, exist_ok=True)
        con = duckdb.connect(
            config={
                "threads": 2,
                "memory_limit": "4GB",
                "temp_directory": str(duck_tmp),
                "preserve_insertion_order": False,
            }
        )
        try:
            base = root / INDEX_DIR
            root_stages = (
                sorted(d.name.removeprefix("stage=") for d in base.glob("stage=*") if d.is_dir())
                if base.is_dir()
                else []
            )
            parts_sql = []
            for s in root_stages:
                files = io.part_files(base / f"stage={s}")
                if files:
                    lit = "'" + s.replace("'", "''") + "'"
                    parts_sql.append(
                        f"SELECT contig_uid, {lit} AS owner FROM read_parquet({_lit(files)}, "
                        "hive_partitioning = false) WHERE is_default_variant"
                    )
            if parts_sql:
                con.execute(
                    "CREATE TEMP TABLE owners AS SELECT contig_uid, min(owner) AS owner, "
                    "count(DISTINCT owner) AS k FROM (" + " UNION ALL ".join(parts_sql) + ") GROUP BY contig_uid"
                )
                dup_root, = con.execute("SELECT count(*) FROM owners WHERE k > 1").fetchall()[0]
                if dup_root:
                    raise RuntimeError(f"root has {dup_root} contigs owned by >1 stage")
                con.execute("ALTER TABLE owners DROP COLUMN k")
            else:
                con.execute("CREATE TEMP TABLE owners (contig_uid BLOB, owner VARCHAR)")

            results: dict[str, dict] = {}
            for st, _, _ in stages:
                gparts, lparts = extracted[st]
                con.execute(f"CREATE OR REPLACE TEMP TABLE cur_g AS SELECT DISTINCT contig_uid FROM read_parquet({_lit(gparts)})")
                con.execute(
                    f"CREATE OR REPLACE TEMP TABLE cur_l AS SELECT contig_uid, max(n) AS n, count(DISTINCT n) AS k "
                    f"FROM read_parquet({_lit(lparts)}) GROUP BY contig_uid"
                )
                two, = con.execute("SELECT count(*) FROM cur_l WHERE k > 1").fetchall()[0]
                missing, = con.execute("SELECT count(*) FROM cur_g ANTI JOIN cur_l USING (contig_uid)").fetchall()[0]
                if two or missing:
                    raise RuntimeError(f"{st}: {two} contigs with two lengths, {missing} listed contigs missing in FASTA")

                held_rows = con.execute(
                    "SELECT o.owner, count(*)::BIGINT, coalesce(sum(l.n), 0)::BIGINT "
                    "FROM cur_g g JOIN owners o USING (contig_uid) JOIN cur_l l USING (contig_uid) "
                    "GROUP BY o.owner ORDER BY o.owner"
                ).fetchall()
                held_by = {str(r[0]): int(r[1]) for r in held_rows}
                held_bases_by = {str(r[0]): int(r[2]) for r in held_rows}
                held_contigs = sum(held_by.values())
                held_bases = sum(held_bases_by.values())

                tot_contigs, tot_bases = con.execute(
                    "SELECT count(*)::BIGINT, coalesce(sum(l.n), 0)::BIGINT "
                    "FROM cur_g g JOIN cur_l l USING (contig_uid)"
                ).fetchall()[0]
                new_contigs, new_bases = con.execute(
                    "SELECT count(*)::BIGINT, coalesce(sum(l.n), 0)::BIGINT "
                    "FROM cur_g g ANTI JOIN owners o USING (contig_uid) JOIN cur_l l USING (contig_uid)"
                ).fetchall()[0]

                # Add newly owned contigs to prefix owners table
                st_lit = "'" + st.replace("'", "''") + "'"
                con.execute(
                    f"INSERT INTO owners SELECT g.contig_uid, {st_lit} FROM cur_g g ANTI JOIN owners o USING (contig_uid)"
                )

                results[st] = {
                    "total_contigs": int(tot_contigs),
                    "total_bases": int(tot_bases),
                    "held_contigs": int(held_contigs),
                    "held_bases": int(held_bases),
                    "held_by": held_by,
                    "held_bases_by": held_bases_by,
                    "new_contigs": int(new_contigs),
                    "new_bases": int(new_bases),
                }
            return results
        finally:
            con.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)


A26_FINALIZE_SPECS: tuple[tuple[str, str, str, str], ...] = (
    (
        "GART1",
        "/tmp/if2453/c63/GART1",
        "/stor/scratch/Wilke/if2453/bio-box-staging/GART/gart1_globdb.fna.gz",
        "c679c6316be2f02434b1bcb32c6e188e4a3f140e700f7b9ed9701efd92a07c3b",
    ),
    (
        "GART2",
        "/tmp/if2453/c63/GART2",
        "/stor/scratch/Wilke/if2453/bio-box-staging/GART/gart2_globdb_cpr.fna.gz",
        "d6f3675de4e6a97b475d213779f96f1b47c5cce8bf495e2ca575ad4bb31a80c5",
    ),
    (
        "LR1",
        "/tmp/if2453/c56/LR1",
        "/stor/scratch/Wilke/if2453/bio-box-staging/LR1/lr1_longread_ge5kb.fna",
        "7164e7ed6ff7179ea62df5c39bd9fd99c10664337536867e319d15e8536d2543",
    ),
    (
        "V3",
        "/tmp/if2453/c55/V3",
        "/stor/scratch/Wilke/if2453/bio-box-staging/V3/vire_v3_dna_votus.fna",
        "b6004f1e811491aca898cda95587073cb82b8e9f824d7f38367b9344d28a7d46",
    ),
    (
        "MV1",
        "/tmp/if2453/c62/MV1",
        "/stor/scratch/Wilke/if2453/bio-box-staging/MV1/metavr_mv1.fna",
        "1461321bb892aef9e5cfa03045dee22c603e9049957683b577e37655fef81394",
    ),
)


def _run_s2_check(stage: str, inputs_dir: Path) -> dict:
    import subprocess

    checker = Path(__file__).resolve().parent / "s2_check_inputs.py"
    proc = subprocess.run(
        [sys.executable, str(checker), stage, str(inputs_dir / "genomes.json"), str(inputs_dir / "objects.json")],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"s2_check_inputs.py {stage} failed (rc={proc.returncode}): {proc.stdout}\n{proc.stderr}")
    doc = json.loads(proc.stdout)
    if doc.get("problems"):
        raise RuntimeError(f"s2_check_inputs.py {stage} reported problems: {doc['problems']}")
    return doc


def _remove_empty_tree(path: Path) -> bool:
    if not path.is_dir():
        return False
    for p in sorted(path.rglob("*"), reverse=True):
        if p.is_file() or p.is_symlink():
            return False
    shutil.rmtree(path)
    return True


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card A-26 input finalization and r26.11.4 linear chain ownership")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_fin = sub.add_parser("finalize-inputs")
    p_fin.add_argument("--stage", required=True)
    p_fin.add_argument("--meta-dir", type=Path, required=True)
    p_fin.add_argument("--fasta-src", type=Path, required=True)
    p_fin.add_argument("--fasta-sha256", required=True)
    p_fin.add_argument("--dst-dir", type=Path, required=True)
    p_fin.add_argument("--min-free-gib", type=int, default=355)
    p_fin.add_argument("--poll-seconds", type=int, default=60)
    p_fin.add_argument("--max-wait-seconds", type=int, default=3000)

    p_own = sub.add_parser("chain-ownership")
    p_own.add_argument("--root", type=Path, required=True)
    p_own.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    p_own.add_argument("--tmp", type=Path, required=True)
    p_own.add_argument("--workers", type=int, default=2)
    p_own.add_argument("--out", type=Path, required=True)

    p_all = sub.add_parser("finalize-all")
    p_all.add_argument("--root", type=Path, default=Path("/tmp/if2453/r10/root"))
    p_all.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    p_all.add_argument("--out-dir", type=Path, default=Path("/tmp/if2453/a26"))
    p_all.add_argument("--workers", type=int, default=3)
    p_all.add_argument("--min-free-gib", type=int, default=355)

    p_ver = sub.add_parser("verify-confs")
    p_ver.add_argument("--tree", type=Path, required=True)
    p_ver.add_argument("--root", type=Path, default=Path("/tmp/if2453/r10/root"))
    p_ver.add_argument("--profile", type=Path, required=True)
    p_ver.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    p_ver.add_argument("--tmp", type=Path, default=Path("/tmp/if2453/a26/est_tmp"))
    p_ver.add_argument("--stages", nargs="+", default=["MGE2", "MV1", "V3", "LR1", "GART1", "GART2"])
    p_ver.add_argument("--out", type=Path, required=True)

    args = ap.parse_args(argv)
    if args.cmd == "finalize-inputs":
        doc = finalize_stage_inputs(
            args.stage,
            args.meta_dir,
            args.fasta_src,
            args.fasta_sha256,
            args.dst_dir,
            min_free_gib=args.min_free_gib,
            poll_seconds=args.poll_seconds,
            max_wait_seconds=args.max_wait_seconds,
        )
        print(json.dumps(doc, indent=1, sort_keys=True))
        return 0
    if args.cmd == "chain-ownership":
        stages = [
            (st, args.inputs_home / st, mlen)
            for st, mlen in (*R26_11_4_CHAIN, *DEFERRED_STAGES)
        ]
        res = compute_chain_ownership(args.root, stages, args.tmp, workers=args.workers)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({"stages": len(res), "out": str(args.out)}, sort_keys=True))
        return 0
    if args.cmd == "finalize-all":
        args.out_dir.mkdir(parents=True, exist_ok=True)
        report: dict = {"finalized": {}, "s2_check": {}, "removed_dup_dirs": [], "removed_work_sym1": False}
        # First check MGE2 (already in /tmp/if2453/s2/inputs/MGE2)
        mge2_dir = args.inputs_home / "MGE2"
        print("[a26] checking MGE2", flush=True)
        report["s2_check"]["MGE2"] = _run_s2_check("MGE2", mge2_dir)
        report["finalized"]["MGE2"] = {
            "stage": "MGE2",
            "dst_dir": str(mge2_dir),
            "files": {
                p.name: {"bytes": p.stat().st_size, "sha256": sha256_file(p)}
                for p in sorted(mge2_dir.iterdir())
                if p.is_file() and not p.name.startswith(".")
            },
        }
        for stage, meta_str, fasta_str, expected_sha in A26_FINALIZE_SPECS:
            meta_dir = Path(meta_str)
            dst_dir = args.inputs_home / stage
            print(f"[a26] finalizing {stage} (free_gib={_free_gib(args.inputs_home)})", flush=True)
            fin = finalize_stage_inputs(
                stage,
                meta_dir,
                Path(fasta_str),
                expected_sha,
                dst_dir,
                min_free_gib=args.min_free_gib,
            )
            report["finalized"][stage] = fin
            print(f"[a26] s2_check_inputs {stage}", flush=True)
            report["s2_check"][stage] = _run_s2_check(stage, dst_dir)
            if meta_dir.is_dir():
                shutil.rmtree(meta_dir)
                report["removed_dup_dirs"].append(str(meta_dir))
                print(f"[a26] removed duplicate metadata dir {meta_dir}", flush=True)
        report["removed_work_sym1"] = _remove_empty_tree(args.inputs_home / "_work_sym1")
        report["free_gib_after_finalize"] = _free_gib(args.inputs_home)
        (args.out_dir / "finalize_report.json").write_text(
            json.dumps(report, indent=1, sort_keys=True) + "\n", encoding="utf-8"
        )
        print("[a26] computing chain ownership across 19 r26.11.4 stages + GV1", flush=True)
        stages = [(st, args.inputs_home / st, mlen) for st, mlen in (*R26_11_4_CHAIN, *DEFERRED_STAGES)]
        own = compute_chain_ownership(args.root, stages, args.out_dir / "work", workers=args.workers)
        (args.out_dir / "chain_ownership.json").write_text(
            json.dumps(own, indent=1, sort_keys=True) + "\n", encoding="utf-8"
        )
        print("[a26] DONE", flush=True)
        return 0
    if args.cmd == "verify-confs":
        import subprocess

        stage_est = args.tree / "scripts" / "stress" / "build_stage_pod" / "stage_est.py"
        stages_dir = args.tree / "scripts" / "stress" / "build_stage_pod" / "stages"
        confs_ok = []
        for conf_path in sorted(stages_dir.glob("*.conf")):
            st = conf_path.stem
            p = subprocess.run(
                [sys.executable, "-I", str(stage_est), "conf", str(conf_path), "--stage", st, "--tree", str(args.tree)],
                capture_output=True,
                text=True,
                check=False,
            )
            if p.returncode != 0:
                raise RuntimeError(f"stage_est.py conf failed on {conf_path.name}: {p.stderr}")
            confs_ok.append(st)

        root_stages = {
            d.name.removeprefix("stage=")
            for d in (args.root / "l0seq" / "contig_index").glob("stage=*")
            if d.is_dir()
        }
        args.tmp.mkdir(parents=True, exist_ok=True)
        before_docs = {}
        for st in args.stages:
            orig = stages_dir / f"{st}.conf"
            filtered_lines = []
            for raw in orig.read_text(encoding="utf-8").splitlines():
                line = raw.strip()
                if line.startswith("EXPECT=held_by "):
                    owner = line.split()[1]
                    if owner not in root_stages:
                        continue
                filtered_lines.append(raw)
            scratch_conf = args.tmp / f"{st}.conf"
            scratch_conf.write_text("\n".join(filtered_lines) + "\n", encoding="utf-8")
            print(f"[a26] stage_est.py before {st}", flush=True)
            p = subprocess.run(
                [
                    sys.executable,
                    str(stage_est),
                    "before",
                    "--conf",
                    str(scratch_conf),
                    "--stage",
                    st,
                    "--root",
                    str(args.root),
                    "--profile",
                    str(args.profile),
                    "--inputs",
                    str(args.inputs_home / st),
                    "--tmp",
                    str(args.tmp / f"work_{st}"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            if p.returncode != 0:
                raise RuntimeError(f"stage_est.py before failed on {st} (rc={p.returncode}): {p.stderr}\n{p.stdout}")
            doc = json.loads(p.stdout)
            if doc.get("refusals"):
                raise RuntimeError(f"{st}: refusals={doc['refusals']}")
            before_docs[st] = doc
        out_doc = {"confs_checked": len(confs_ok), "confs": confs_ok, "before": before_docs}
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(out_doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        print(f"[a26] verify-confs DONE ({len(confs_ok)} confs, {len(before_docs)} before checks)", flush=True)
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

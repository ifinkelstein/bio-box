"""Library domain scans with pyhmmer (operation O-SCAN).

Every distinct protein of a gene-calling version is searched with every model of
one HMM library, using the library's own bit-score cutoffs (Pfam GA, NCBIfam TC)
and a fixed search space (Z = domZ = 1e8, SEARCH_PLAN SD-04, DEC-04), so stored
E-values do not change as the corpus grows. Only included domains are kept.

Writes `<root>/l1/domains/library=<name>-<version>/stage=<S>/part-<k>.parquet` through
`biobox.io.write_parts`, cut at `max_file_bytes` of estimated size (card B-19), and the step's
ledger `ledger/stage=<S>/step=domains-<name>-<version>`.

Two ways to scan a stage, with byte-identical output (card B-24):

- `run`: the stage's proteins in one process.
- In shards, each in its own process under the 300 CPU-minute limit (AGENTS §7.2):
  `plan` cuts the stage's proteins, sorted by protein_uid, into shards of at most `shard_residues`
  residues each (cumulative `length`; a protein longer than that is a shard of its own) and writes
  the plan, `l1/domain_scan/library=<L>/stage=<S>/plan.json`; `scan_shard` scans shard k into a
  sorted scratch run, `<tmp>/biobox-domain-scan/<plan>/shard-<k>.parquet` (`<plan>`: the first 16
  hex of the plan's sha256), and writes its record `.../shard-<k>.json` beside the plan (the run's
  sha256 and counts); `merge` checks every run against its record and writes the runs, in shard
  order, as the files and ledger `run` writes. If scratch is cleared between the shards and the
  merge (a reboot clears /tmp), `rescan_shard` (`python -m biobox.l1.domains rescan`, card C-12)
  scans a shard again into scratch without writing its record again: the new run must have the
  sha256 the record holds, so the merge is of the same runs as before.

At a fixed Z a protein's hits do not depend on the other proteins searched (MR-12), the shards are
consecutive ranges of protein_uid, and `scan` sorts each one's rows by the table's sort key, so the
runs in shard order are the rows of `run` in its order, and `io.write_parts` cuts the same parts
from them. The plan and the shard records are the `domain_scan` artifact: small by nature, never in
a release manifest (a release is the same whether its scans were sharded or not), and what lets
`bio-box run-step` publish each step (a step must publish a file the registry claims). The runs stay
in scratch on local disk, never in the data root. The plan names its `protein_seq` files and the
sha256 of the HMM file, which every shard checks. Thread counts are the caller's (`cpus`).

Copies of measured rows (card B-47b, DEC-35 (2)). `plan` given `seq_stages` that name a stage other than its own
(build-stage passes every stage on disk, as it does to gene calling) leaves out every protein already measured with
the same scan in one of them, and `merge` copies those proteins' rows instead. At a fixed Z a protein's hits do not
depend on the other proteins (MR-12), so a copied row is the row a scan would write and the stage's domain files and
ledger are those of a full scan, byte for byte, whichever stages were built first. A named stage is a source only
when (`sources`):

- its `domains-<lib>` step is complete: its run record `runs/stage=<S>/domains-<lib>.json` exists (and no
  `.publishing.json` beside it), and the record lists its domain files and its ledger (with sha256);
- its shard plan exists: the plan names the proteins it measured (its `files`: every protein of those
  `protein_seq` files was scanned there, or copied there from a source of its own) and the scan's parameters. A
  stage scanned in one process (`run`) has no plan, so the sha256 of its HMM file is recorded nowhere, and it is not
  a source: its proteins are scanned again;
- its recorded scan is this one (`COMPARED`): the plan's `hmm_sha256` and `cutoff`; its `z` and `dom_z` (a plan
  written before this card holds neither, and is taken at 1e8: `Z_PROT` and `DOM_Z` have had that value since the
  module was written, `git log -G` on them); and the `tool_versions_json` of every row of its ledger (pyhmmer and
  pyarrow), against this process's `provenance.tool_versions`.

"Measured" is decided from those records, never from the presence of rows: a protein scanned with no hit has no row,
and is measured (it is copied as no rows). When several sources measured a protein it is copied from the first in
the order `seq_stages` gives; its rows are equal in each (the rule above). The proteins are matched by a merge of
sorted streams (`_assign`): the stage's `protein_uid`s and each source's, in uid order, one row group at a time, never
a set of a stage's proteins (AGENTS §4.1). The plan then holds, beside today's keys, `z`, `dom_z`, `tools` and a
`copy` block: the stages named, the fields compared, each source (its plan and run record with sha256, the
`protein_seq` files it measured, the domain files copied from with sha256, and the proteins, residues, rows and
proteins with a hit copied from it), the stages passed over and why, and the totals. Its `proteins` and `residues`
stay the stage's (the ledger's `n_in`); its shards are cut from the residues actually scanned, and each shard leaves
out the copied proteins of its uid range (`_shard_proteins`). The merge checks each source file against the plan
(bytes; sha256 for the domain files), copies the rows of the copied proteins in uid order between the shards' rows,
checks every copied row's `cutoff`, `z` and `dom_z`, and refuses when the copied counts are not the plan's. The
scanned and copied counts live in the plan and shard records only (`domain_scan`): never in a manifest, not in the
ledger, whose rows are a full scan's. Without another stage named the plan, the shards and the merge are exactly
those of card B-24.

Inferred rows (card B-47c) are measured nowhere and are never copied. The two places that card must narrow are
`_measured_uids` (the proteins a source measured: today every protein of its plan's files) and `_measured_rows`
(the rows taken from a source's domain files: today all of them; there, `inferred_from IS NULL`).
"""

from __future__ import annotations

import argparse
import bisect
import contextlib
import gzip
import hashlib
import io
import json
import os
import shutil
import sys
import tarfile
import tempfile
from collections.abc import Iterator
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyhmmer

from biobox import ids, layout, ledger, provenance
from biobox import io as bio
from biobox.registry.schemas import ARTIFACTS, table

Z_PROT = 1e8
DOM_Z = 1e8
# estimate: zstd bytes per domain row; cuts parts only. Measured 53.8 B/row for Pfam and NCBIfam on the
# golden set (planning/measurements/2026-11-09_design_stress_test.md §3).
EST_ROW_BYTES = 54
PLAN_FORMAT = 1
SCAN_DIR = "biobox-domain-scan"  # under the scratch directory: the shards' runs, one directory per plan
_READ_ROWS = 16_384  # protein_seq rows per read while planning or reading a shard
# card B-47b: what a source stage's records must equal for its rows to be copied (module docstring)
COMPARED = ("hmm_sha256", "cutoff", "z", "dom_z", "tools")
Z_UNRECORDED = 1e8  # the Z of a plan that records none: Z_PROT and DOM_Z have been 1e8 since 74b3e05 (`git log -G`)
FLUSH_PIECES = 1024  # the merge of a copying plan combines its row slices at least this often (`_merge_copies`)
DUCK_MEMORY_LIMIT = "256MB"  # single-thread DuckDB cap when streaming shared-contig proteins (card C-38h, AGENTS §4.1)


def iter_hmms(path: Path):
    """Read HMMs from a .hmm, .hmm.gz, or a .tgz of .HMM files (NCBIfam), in file order."""
    path = Path(path)
    name = path.name
    if name.endswith((".tgz", ".tar.gz")):
        with tarfile.open(path, "r:gz") as tar:
            for member in sorted((m for m in tar.getmembers() if m.isfile() and m.name.upper().endswith(".HMM")),
                                 key=lambda m: m.name):
                data = tar.extractfile(member).read()
                with pyhmmer.plan7.HMMFile(io.BytesIO(data)) as hf:
                    yield from hf
        return
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rb") as fh, pyhmmer.plan7.HMMFile(fh) as hf:
        yield from hf


def _protein_files(root: Path, gene_call_version: str, stage: str) -> list[Path]:
    return sorted((Path(root) / "l0seq/protein_seq").glob(
        f"gene_call_version={gene_call_version}/stage={stage}/protein_prefix=*/part-*.parquet"))


def load_proteins(root: Path, gene_call_version: str, stage: str) -> dict[bytes, str]:
    out: dict[bytes, str] = {}
    for f in _protein_files(root, gene_call_version, stage):
        t = pq.read_table(f, columns=["protein_uid", "seq"])
        out.update(zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist()))
    return out


def _decode(x) -> str:
    return x.decode() if isinstance(x, bytes) else (x or "")


def scan(
    proteins: dict[bytes, str],
    hmm_path: Path,
    cutoff: str,
    cpus: int = 1,
) -> list[dict]:
    """Search `proteins` with every model in `hmm_path`; return included domains as rows."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    uids = sorted(proteins)
    block = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=u.hex().encode(), sequence=proteins[u].replace("*", "")) for u in uids
    ).digitize(alphabet)
    kwargs = {"cpus": cpus, "Z": Z_PROT, "domZ": DOM_Z}
    if cutoff != "none":
        kwargs["bit_cutoffs"] = cutoff
    rows = []
    for hits in pyhmmer.hmmsearch(iter_hmms(hmm_path), block, **kwargs):
        q = hits.query
        acc, name, mlen = _decode(q.accession) or _decode(q.name), _decode(q.name), q.M
        for hit in hits:
            if not hit.included:
                continue
            for dom in hit.domains:
                if not dom.included:
                    continue
                a = dom.alignment
                rows.append({
                    "protein_uid": bytes.fromhex(_decode(hit.name)),
                    "model_acc": acc,
                    "model_name": name,
                    "model_len": mlen,
                    "seq_bitscore": hit.score,
                    "seq_evalue_zref": hit.evalue,
                    "dom_bitscore": dom.score,
                    "dom_ievalue_zref": dom.i_evalue,
                    "env_start": dom.env_from - 1,
                    "env_end": dom.env_to,
                    "ali_start": a.target_from - 1,
                    "ali_end": a.target_to,
                    "hmm_from": a.hmm_from - 1,
                    "hmm_to": a.hmm_to,
                    "cutoff": cutoff,
                    "z": Z_PROT,
                    "dom_z": DOM_Z,
                })
    rows.sort(key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
    return rows


def _batches(rows: list[dict], spec) -> list[pa.RecordBatch]:
    return [pa.RecordBatch.from_pylist(rows[i : i + bio.ROW_GROUP_ROWS], schema=spec.schema)
            for i in range(0, len(rows), bio.ROW_GROUP_ROWS)]


def _names(stage: str, library: str, library_version: str) -> tuple[str, str]:
    """The domain files' directory and the step's name."""
    lib = f"{library}-{library_version}"
    return f"l1/domains/library={lib}/stage={stage}", f"domains-{lib}"


def _tools() -> dict[str, str]:
    """The versions of the tools that shape the domain rows, as the ledger records them."""
    return provenance.tool_versions("pyhmmer", "pyarrow")


def _write(out: Path, stage: str, lib_dir: str, step: str, batches, n_in: int, hit: int,
           max_file_bytes: int) -> None:
    """The stage's domain files from sorted `batches`, and the step's ledger: `n_in` proteins in, `hit`
    with an included domain (the one writer of `run` and `merge`, so their files are the same)."""
    spec = table("l1_domains")
    bio.write_parts(batches, spec, out / lib_dir, EST_ROW_BYTES, max_file_bytes)
    lib = lib_dir.split("/")[2].removeprefix("library=").encode()
    tools = _tools()
    ledger.write(out, stage, step,
                 sorted([ledger.row(lib, "library", n_in, hit, tools=tools),
                         ledger.row(lib, "library", n_in, hit, n_rejected=n_in - hit,
                                    reason="no_hit", tools=tools)],
                        key=lambda r: (r["object_key"], r["reason_code"])))


def run(root: Path, stage: str, gene_call_version: str, library: str, library_version: str,
        hmm_path: Path, cutoff: str, cpus: int = 1, max_file_bytes: int = layout.MAX_FILE_BYTES,
        out: Path | None = None) -> dict:
    """Scan the stage's distinct proteins (read from `root`) with one library; write its domains and
    ledger under `out` (default `root`; card C-08). Refuses, before reading a protein, when either
    already holds a file in `root` or `out`."""
    out = Path(root) if out is None else Path(out)
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, out, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    proteins = load_proteins(root, gene_call_version, stage)
    rows = scan(proteins, hmm_path, cutoff, cpus=cpus)
    hit = len({r["protein_uid"] for r in rows})
    _write(out, stage, lib_dir, step, _batches(rows, table("l1_domains")), len(proteins), hit, max_file_bytes)
    return {"proteins": len(proteins), "domains": len(rows),
            "proteins_with_hit": hit,
            "models_with_hit": len({r["model_acc"] for r in rows})}


# --- the scan in shards (card B-24) -----------------------------------------------------------------


def scan_dir(stage: str, library: str, library_version: str) -> str:
    """The root-relative directory of the `domain_scan` artifact for one stage and library."""
    return f"{ARTIFACTS['domain_scan'].location}/library={library}-{library_version}/stage={stage}"


def plan_stages(root: Path, library: str, library_version: str) -> list[str]:
    """Stages on disk that have a `domain_scan` directory for `library`-`library_version` (card C-38h: includes
    stages that own 0 contigs of their own and therefore have no `seq_locator` directory)."""
    base = Path(root) / ARTIFACTS["domain_scan"].location / f"library={library}-{library_version}"
    return sorted(p.name.removeprefix("stage=") for p in base.glob("stage=*") if p.is_dir()) if base.is_dir() else []


def _plan_rel(stage: str, library: str, library_version: str) -> str:
    return f"{scan_dir(stage, library, library_version)}/plan.json"


def _shard_rel(stage: str, library: str, library_version: str, k: int) -> str:
    return f"{scan_dir(stage, library, library_version)}/shard-{k:05d}.json"


def _json_bytes(doc: dict) -> bytes:
    return json.dumps(doc, sort_keys=True, indent=1).encode() + b"\n"


def _write_once(root: Path, out: Path, rel: str, data: bytes, what: str) -> None:
    for base in dict.fromkeys([Path(root).resolve(), Path(out).resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: {what} is written once; refusing to overwrite")
    path = Path(out) / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _find(root: Path, out: Path | None, rel: str) -> Path | None:
    """`rel` in the root, else in `out` (written earlier by a direct call into another output root)."""
    for base in [Path(root)] + ([] if out is None else [Path(out)]):
        if (base / rel).is_file():
            return base / rel
    return None


def _plan_bytes(root: Path, stage: str, library: str, library_version: str, out: Path | None = None) -> bytes:
    p = _find(root, out, _plan_rel(stage, library, library_version))
    if p is None:
        raise FileNotFoundError(f"{Path(root) / _plan_rel(stage, library, library_version)}: no shard plan for "
                                f"domains-{library}-{library_version} of stage {stage}; write it with domains.plan")
    return p.read_bytes()


def read_plan(root: Path, stage: str, library: str, library_version: str, out: Path | None = None) -> dict:
    """The stage's shard plan for one library (from the root, else from `out`)."""
    return json.loads(_plan_bytes(root, stage, library, library_version, out))


def run_path(tmp: Path | None, root: Path, stage: str, library: str, library_version: str, k: int,
             out: Path | None = None) -> Path:
    """Where `scan_shard` writes shard k's sorted run: under the scratch directory (`tmp`, else the
    system's), in a directory named by the plan's sha256, so runs of another plan never mix in."""
    digest = hashlib.sha256(_plan_bytes(root, stage, library, library_version, out)).hexdigest()[:16]
    return Path(tmp if tmp is not None else tempfile.gettempdir()) / SCAN_DIR / digest / f"shard-{k:05d}.parquet"


def _uid_lengths(files: list[Path]) -> Iterator[tuple[bytes, int]]:
    """(protein_uid, length) of the stage's `protein_seq` files in uid order: in bucket and part order the
    files are one sorted run (the table's sort key, bucketed by uid prefix); a uid that does not follow
    the one before is a RuntimeError."""
    last = None
    for f in files:
        for b in pq.ParquetFile(f).iter_batches(batch_size=_READ_ROWS, columns=["protein_uid", "length"]):
            for u, n in zip(b.column(0).to_pylist(), b.column(1).to_pylist()):
                if last is not None and u <= last:
                    raise RuntimeError(f"{f}: protein_seq rows are not in strictly increasing protein_uid order")
                last = u
                yield u, n


def _cut(rows: Iterator[tuple[bytes, int]], shard_residues: int) -> Iterator[dict]:
    """Consecutive shards of at most `shard_residues` residues; a protein that would not fit closes the
    shard before it (a protein above the limit is a shard of its own)."""
    cur = None
    for u, n in rows:
        if cur is not None and cur["residues"] + n > shard_residues:
            yield cur
            cur = None
        if cur is None:
            cur = {"first": u.hex(), "last": u.hex(), "proteins": 0, "residues": 0}
        cur["last"] = u.hex()
        cur["proteins"] += 1
        cur["residues"] += n
    if cur is not None:
        yield cur


def plan(root: Path, stage: str, gene_call_version: str, library: str, library_version: str, hmm_path: Path,
         cutoff: str, shard_residues: int, out: Path | None = None, seq_stages: list[str] | None = None) -> dict:
    """Write the stage's shard plan for one library under `out` (default `root`): its proteins sorted by
    protein_uid, cut into shards of at most `shard_residues` residues. Refuses when the stage's domains of
    this library are written already, or a plan exists, in `root` or `out`. Returns the shard count and
    the totals.

    `seq_stages` (card B-47b): when it names a stage other than `stage`, the proteins that the first source among
    them measured with this scan (`sources`) are left out of the shards and recorded for the merge to copy (module
    docstring); the plan's `proteins` and `residues` stay the stage's, and the result adds the copied and scanned
    counts. Without another stage named the plan is card B-24's, byte for byte."""
    root = Path(root)
    out = root if out is None else Path(out)
    if not isinstance(shard_residues, int) or shard_residues < 1:
        raise ValueError(f"shard_residues {shard_residues!r}: must be a positive integer")
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, out, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    rel = _plan_rel(stage, library, library_version)
    for base in dict.fromkeys([root.resolve(), out.resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: the shard plan is written once; refusing to overwrite")
    files = _protein_files(root, gene_call_version, stage)
    named = list(dict.fromkeys(s for s in (seq_stages or ()) if s != stage))
    owned = _owned_plan(root, stage, gene_call_version, named) if named else None
    if not files and owned is None and not (named and _genome_input_files(root, stage)):
        raise FileNotFoundError(f"stage {stage}: no protein_seq files for gene_call_version={gene_call_version}")
    hmm_sha256 = ids.sha256_file(Path(hmm_path)).hex()
    head = {"format": PLAN_FORMAT, "stage": stage, "gene_call_version": gene_call_version, "library": library,
            "library_version": library_version, "cutoff": cutoff, "hmm_sha256": hmm_sha256,
            "shard_residues": shard_residues,
            "files": [{"path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size} for f in files]}
    has_rounds = seq_stages is not None and bool(round_scan_ids(root, library, library_version))
    if not named and not has_rounds:
        shards = list(_cut(_uid_lengths(files), shard_residues))
        doc = {**head, "proteins": sum(s["proteins"] for s in shards), "residues": sum(s["residues"] for s in shards),
               "shards": shards}
        _write_once(root, out, rel, _json_bytes(doc), "the shard plan")
        return {"shards": len(shards), "proteins": doc["proteins"], "residues": doc["residues"]}
    srcs, passed = sources(root, stage, library, library_version, hmm_sha256, cutoff, named)
    total = {"proteins": 0, "residues": 0}
    counts = [{"proteins": 0, "residues": 0, "rows": 0, "proteins_with_hit": 0} for _ in srcs]

    def scanned() -> Iterator[tuple[bytes, int]]:
        for u, n, k, rows in _walk(root, files, srcs, ["protein_uid"], owned=owned):
            total["proteins"] += 1
            total["residues"] += n
            if k < 0:
                yield u, n
                continue
            c = counts[k]
            c["proteins"] += 1
            c["residues"] += n
            got = sum(t.num_rows for t in rows)
            c["rows"] += got
            c["proteins_with_hit"] += got > 0

    shards = list(_cut(scanned(), shard_residues))
    src_docs = [{**s, **c} for s, c in zip(srcs, counts)]
    copy = {"seq_stages": list(seq_stages or ()), "compared": list(COMPARED), "sources": src_docs,
            "passed_over": passed,
            **{k: sum(c[k] for c in counts) for k in ("proteins", "residues", "rows", "proteins_with_hit")}}
    doc = {**head, "proteins": total["proteins"], "residues": total["residues"], "shards": shards,
           "z": Z_PROT, "dom_z": DOM_Z, "tools": _tools(), "copy": copy}
    if owned is not None:
        doc["owned"] = owned
    _write_once(root, out, rel, _json_bytes(doc), "the shard plan")
    scanned_p, scanned_r = sum(s["proteins"] for s in shards), sum(s["residues"] for s in shards)
    return {"shards": len(shards), "proteins": doc["proteins"], "residues": doc["residues"],
            "scanned_proteins": scanned_p, "scanned_residues": scanned_r, "copied_proteins": copy["proteins"],
            "copied_residues": copy["residues"], "copied_rows": copy["rows"],
            "sources": [_src_id(s) for s in srcs], "passed_over": [_src_id(p) for p in passed]}


# --- shared contigs owned by earlier stages (card C-38h, DEC-26, DEC-54 Q8) ------------------------------------------


def _sql_files(files: list[Path]) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


@contextlib.contextmanager
def _duck() -> Iterator[duckdb.DuckDBPyConnection]:
    scratch = Path(tempfile.mkdtemp(prefix="biobox-domains-"))
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        yield con
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)


def _feature_protein_files(root: Path, gene_call_version: str, stage: str) -> list[Path]:
    return bio.part_files(Path(root) / table("feature_protein").location
                          / f"gene_call_version={gene_call_version}" / f"stage={stage}")


def _genome_input_files(root: Path, stage: str) -> list[Path]:
    return bio.part_files(Path(root) / table("genome_input").location / f"stage={stage}")


def _file_entries(root: Path, files: list[Path]) -> list[dict]:
    return [{"path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size} for f in files]


def _owned_plan(root: Path, stage: str, gene_call_version: str, named: list[str]) -> dict | None:
    """Other stages in `named` whose `feature_protein` holds proteins on contigs of `stage`'s `genome_input`
    (card C-38h, DEC-26, DEC-54 Q8): None when `stage` has no `genome_input` or no stage in `named` contributes
    any protein, else `{"genome_input": [...], "stages": [{"stage": s, "feature_protein": [...],
    "protein_files": [...]}, ...]}` sorted by stage name."""
    gi = _genome_input_files(root, stage)
    if not gi or not named:
        return None
    cands = []
    for s in sorted(dict.fromkeys(x for x in named if x != stage)):
        fp = _feature_protein_files(root, gene_call_version, s)
        ps = _protein_files(root, gene_call_version, s)
        if fp and ps:
            cands.append((s, fp, ps))
    if not cands:
        return None
    unions = " UNION ALL ".join(
        f"SELECT {i} AS si, contig_uid, protein_uid FROM read_parquet({_sql_files(fp)}, hive_partitioning = false)"
        for i, (_, fp, _) in enumerate(cands)
    )
    sql = (f"WITH gi AS (SELECT DISTINCT contig_uid FROM read_parquet({_sql_files(gi)}, hive_partitioning = false)) "
           f"SELECT DISTINCT fp.si FROM ({unions}) fp SEMI JOIN gi USING (contig_uid) ORDER BY fp.si")
    with _duck() as con:
        hit_idxs = [int(i) for (i,) in con.execute(sql).fetchall()]
    if not hit_idxs:
        return None
    return {
        "genome_input": _file_entries(root, gi),
        "stages": [
            {"stage": cands[i][0],
             "feature_protein": _file_entries(root, cands[i][1]),
             "protein_files": _file_entries(root, cands[i][2])}
            for i in hit_idxs
        ],
    }


def shared_owner_stages(root: Path, stage: str, gene_call_version: str, seq_stages: list[str] | None) -> list[str]:
    """Stages in `seq_stages` (other than `stage`) that own contigs of `stage`'s `genome_input` with called
    proteins (card C-38h): used by `stage._sharded` so a stage sharing contigs with a called stage always runs
    through `domains.plan --seq-stages`."""
    owned = _owned_plan(Path(root), stage, gene_call_version, list(seq_stages or ()))
    return [] if owned is None else [s["stage"] for s in owned["stages"]]


def _check_owned(root: Path, owned: dict, what: str) -> tuple[list[Path], list[tuple[list[Path], list[Path]]]]:
    """Verify byte sizes of every file listed in an `owned` block; return `(gi_paths, [(fp_paths, ps_paths), ...])`."""
    gi = _checked(root, owned["genome_input"], f"{what} (genome_input)")
    owners = []
    for o in owned["stages"]:
        fp = _checked(root, o["feature_protein"], f"{what} (feature_protein of {o['stage']})")
        ps = _checked(root, o["protein_files"], f"{what} (protein_seq of {o['stage']})")
        owners.append((fp, ps))
    return gi, owners


def _stage_uid_lengths(root: Path, files: list[Path], owned: dict | None) -> Iterator[tuple[bytes, int]]:
    """Distinct `(protein_uid, length)` of `files` union the shared-contig proteins of `owned`, in strictly
    increasing `protein_uid` order. Without `owned`, streams `_uid_lengths(files)` directly."""
    if owned is None:
        yield from _uid_lengths(files)
        return
    gi, owners = _check_owned(root, owned, "the shard plan")
    parts = []
    if files:
        parts.append(f"SELECT protein_uid, length FROM read_parquet({_sql_files(files)}, hive_partitioning = false)")
    for fp, ps in owners:
        parts.append(
            f"SELECT ps.protein_uid, ps.length "
            f"FROM read_parquet({_sql_files(ps)}, hive_partitioning = false) ps "
            f"SEMI JOIN ("
            f"  SELECT DISTINCT fp.protein_uid "
            f"  FROM read_parquet({_sql_files(fp)}, hive_partitioning = false) fp "
            f"  SEMI JOIN gi USING (contig_uid)"
            f") u USING (protein_uid)"
        )
    sql = (f"WITH gi AS (SELECT DISTINCT contig_uid FROM read_parquet({_sql_files(gi)}, hive_partitioning = false)) "
           f"SELECT protein_uid, min(length)::BIGINT AS length "
           f"FROM ({' UNION ALL '.join(parts)}) "
           f"GROUP BY protein_uid ORDER BY protein_uid")
    with _duck() as con:
        reader = con.sql(sql).to_arrow_reader(_READ_ROWS)
        try:
            for b in reader:
                for u, n in zip(b.column(0).to_pylist(), b.column(1).to_pylist()):
                    yield bytes(u), int(n)
        finally:
            reader.close()


# --- copies of measured rows (card B-47b) ----------------------------------------------------------------------------


def _record_rel(stage: str, step: str) -> str:
    return f"runs/stage={stage}/{step}.json"  # where run-step writes a step's run record (`runstep.run_step`)


def _checked(root: Path, files: list[dict], what: str) -> list[Path]:
    """The paths of `files` (each `path`, `bytes`, and optionally `sha256`) under `root`, each checked against its
    entry; a missing file or other bytes is a ValueError naming `what`."""
    out = []
    for f in files:
        path = Path(root) / f["path"]
        if not path.is_file():
            raise ValueError(f"{path}: missing ({what} lists it)")
        if path.stat().st_size != f["bytes"]:
            raise ValueError(f"{path}: {path.stat().st_size} bytes, {what} lists {f['bytes']}")
        if "sha256" in f and ids.sha256_file(path).hex() != f["sha256"]:
            raise ValueError(f"{path}: sha256 differs from the one {what} lists")
        out.append(path)
    return out


def _source(root: Path, s: str, library: str, library_version: str, want: dict) -> tuple[dict | None, str | None]:
    """(stage `s` as a source of copies, None), or (None, why it is not one): module docstring."""
    lib_dir, step = _names(s, library, library_version)
    rec_rel = _record_rel(s, step)
    rec_path = Path(root) / rec_rel
    if not rec_path.is_file():
        return None, (f"no run record {rec_rel}: its {step} step is not complete (or did not run under run-step), so "
                      f"what it measured is not recorded")
    if (Path(root) / f"runs/stage={s}/{step}.publishing.json").exists():
        return None, f"its {step} step is being published (runs/stage={s}/{step}.publishing.json)"
    plan_rel = _plan_rel(s, library, library_version)
    if not (Path(root) / plan_rel).is_file():
        return None, (f"no shard plan {plan_rel}: scanned in one process, which records no sha256 of its HMM file "
                      f"and no list of the proteins it measured")
    raw_plan, raw_rec = (Path(root) / plan_rel).read_bytes(), rec_path.read_bytes()
    p, rec = json.loads(raw_plan), json.loads(raw_rec)
    if (p.get("stage"), p.get("library"), p.get("library_version")) != (s, library, library_version):
        return None, f"its plan {plan_rel} is of another stage or library"
    have = {"hmm_sha256": p.get("hmm_sha256"), "cutoff": p.get("cutoff"), "z": p.get("z", Z_UNRECORDED),
            "dom_z": p.get("dom_z", Z_UNRECORDED)}
    for k in ("hmm_sha256", "cutoff", "z", "dom_z"):
        if have[k] != want[k]:
            return None, f"its plan records {k} {have[k]!r}, this scan's is {want[k]!r}"
    outs = rec.get("outputs", [])
    files = [o for o in outs if o["path"].startswith(lib_dir + "/")]
    led = [o for o in outs if o["path"].startswith(f"ledger/stage={s}/step={step}/")]
    if not files or not led:
        return None, f"its run record {rec_rel} lists no domain files or no ledger of {step}"
    try:
        _checked(root, [{k: o[k] for k in ("path", "bytes")} for o in files], f"its run record {rec_rel}")
        _checked(root, p["files"], f"its plan {plan_rel}")
        if "owned" in p:
            _check_owned(root, p["owned"], f"its plan {plan_rel}")
        led_paths = _checked(root, led, f"its run record {rec_rel}")
    except ValueError as e:
        return None, str(e)
    tools = {r for f in led_paths for r in pq.read_table(f, columns=["tool_versions_json"], partitioning=None)
             .column(0).to_pylist()}
    mine = json.dumps(want["tools"], sort_keys=True)
    if tools != {mine}:
        return None, f"its ledger records tool versions {sorted(map(str, tools))}, this scan's are {mine}"
    out = {"stage": s, "plan": {"path": plan_rel, "sha256": hashlib.sha256(raw_plan).hexdigest()},
           "record": {"path": rec_rel, "sha256": hashlib.sha256(raw_rec).hexdigest()},
           "protein_files": p["files"],
           "files": [{"path": o["path"], "bytes": o["bytes"], "sha256": o["sha256"]} for o in files]}
    if "owned" in p:
        out["owned"] = p["owned"]
    return out, None


def _src_id(src: dict) -> str:
    """Identifier of a copy source dict (stage name or round id)."""
    return src["round"] if "round" in src else src["stage"]


def _round_source(root: Path, rid: str, library: str, library_version: str, want: dict) -> tuple[dict | None, str | None]:
    """(round `rid` as a source of copies, None), or (None, why it is not one): card C-38f, §3.6."""
    lib = f"{library}-{library_version}"
    step = f"domains-{lib}"
    if (Path(root) / f"runs/stage={rid}/{step}.publishing.json").exists():
        return None, f"round {rid} {step} is being published (runs/stage={rid}/{step}.publishing.json)"
    plan_rel = _round_plan_rel(rid, library, library_version)
    if not (Path(root) / plan_rel).is_file():
        return None, f"no round plan {plan_rel}"
    raw_plan = (Path(root) / plan_rel).read_bytes()
    p = json.loads(raw_plan)
    if (p.get("round"), p.get("library"), p.get("library_version")) != (rid, library, library_version):
        return None, f"its plan {plan_rel} is of another round or library"
    have = {"hmm_sha256": p.get("hmm_sha256"), "cutoff": p.get("cutoff"), "z": p.get("z", Z_UNRECORDED),
            "dom_z": p.get("dom_z", Z_UNRECORDED)}
    for k in ("hmm_sha256", "cutoff", "z", "dom_z"):
        if have[k] != want[k]:
            return None, f"its plan records {k} {have[k]!r}, this scan's is {want[k]!r}"
    dom_dir = round_domains_dir(rid, library, library_version)
    cov_dir = round_coverage_dir(rid, library, library_version)
    led_dir = round_ledger_dir(rid, library, library_version)
    rec_rel = _record_rel(rid, step)
    rec_path = Path(root) / rec_rel
    rec_info = None
    if rec_path.is_file():
        raw_rec = rec_path.read_bytes()
        rec = json.loads(raw_rec)
        outs = rec.get("outputs", [])
        dom_entries = [o for o in outs if o["path"].startswith(dom_dir + "/")]
        cov_entries = [o for o in outs if o["path"].startswith(cov_dir + "/")]
        led_entries = [o for o in outs if o["path"].startswith(led_dir + "/")]
        if not dom_entries or not cov_entries or not led_entries:
            return None, f"its run record {rec_rel} lacks round_domains, round_coverage or ledger of {step}"
        try:
            _checked(root, [{k: o[k] for k in ("path", "bytes")} for o in dom_entries], f"its run record {rec_rel}")
            _checked(root, [{k: o[k] for k in ("path", "bytes")} for o in cov_entries], f"its run record {rec_rel}")
            led_paths = _checked(root, led_entries, f"its run record {rec_rel}")
        except ValueError as e:
            return None, str(e)
        dom_files = [{"path": o["path"], "bytes": o["bytes"], "sha256": o["sha256"]} for o in dom_entries]
        cov_files = [{"path": o["path"], "bytes": o["bytes"]} for o in cov_entries]
        rec_info = {"path": rec_rel, "sha256": hashlib.sha256(raw_rec).hexdigest()}
    else:
        dom_paths = bio.part_files(Path(root) / dom_dir)
        cov_paths = bio.part_files(Path(root) / cov_dir)
        led_paths = bio.part_files(Path(root) / led_dir)
        if not dom_paths or not cov_paths or not led_paths:
            return None, f"round {rid} is incomplete (missing round_domains, round_coverage or ledger)"
        dom_files = [{"path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size,
                      "sha256": ids.sha256_file(f).hex()} for f in dom_paths]
        cov_files = [{"path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size} for f in cov_paths]
    tools = {r for f in led_paths for r in pq.read_table(f, columns=["tool_versions_json"], partitioning=None)
             .column(0).to_pylist()}
    mine = json.dumps(want["tools"], sort_keys=True)
    if tools != {mine}:
        return None, f"its ledger records tool versions {sorted(map(str, tools))}, this scan's are {mine}"
    out = {"kind": "round", "round": rid, "plan": {"path": plan_rel, "sha256": hashlib.sha256(raw_plan).hexdigest()},
           "coverage_files": cov_files, "files": dom_files}
    if rec_info is not None:
        out["record"] = rec_info
    return out, None


def sources(root: Path, stage: str, library: str, library_version: str, hmm_sha256: str, cutoff: str,
            seq_stages: list[str] | None, *, include_rounds: bool = True,
            exclude_round: str | None = None) -> tuple[list[dict], list[dict]]:
    """(the sources of copies among `seq_stages` and completed rounds on disk, in order, each named once and never
    `stage` or `exclude_round` itself; the sources passed over, each with why), for a scan of `library`-`library_version`
    with the HMM file of sha256 `hmm_sha256` and `cutoff` in this process (cards B-47b, C-38f §3.6). Reads only small
    records: run records, plans and the sources' ledgers."""
    want = {"hmm_sha256": hmm_sha256, "cutoff": cutoff, "z": Z_PROT, "dom_z": DOM_Z, "tools": _tools()}
    srcs, passed = [], []
    for s in dict.fromkeys(x for x in (seq_stages or ()) if x != stage):
        src, why = _source(Path(root), s, library, library_version, want)
        if src is None:
            passed.append({"stage": s, "why": why})
        else:
            srcs.append(src)
    if include_rounds:
        for rid in round_scan_ids(Path(root), library, library_version):
            if rid == exclude_round:
                continue
            src, why = _round_source(Path(root), rid, library, library_version, want)
            if src is None:
                passed.append({"round": rid, "why": why})
            else:
                srcs.append(src)
    return srcs, passed


def _measured_uids(root: Path, src: dict, lo: bytes | None = None, hi: bytes | None = None) -> Iterator[bytes]:
    """The proteins source `src` measured, in strictly increasing uid order (within [lo, hi] when given, reading only
    the row groups that meet it): for a stage, every protein of the protein_seq files its plan lists plus any
    shared-contig proteins in its `owned` block (card C-38h); for a round, every protein of its `round_coverage` files
    (card C-38f §3.6). Card B-47c: inferred proteins are measured nowhere and must be left out here."""
    if "coverage_files" in src:
        last = None
        for path in _checked(root, src["coverage_files"], f"the plan of {_src_id(src)}"):
            pf = pq.ParquetFile(path)
            col = pf.schema_arrow.get_field_index("protein_uid")
            for g in range(pf.metadata.num_row_groups):
                if lo is not None:
                    st = pf.metadata.row_group(g).column(col).statistics
                    if st is not None and st.has_min_max and (bytes(st.max) < lo or bytes(st.min) > hi):
                        continue
                for u in pf.read_row_group(g, columns=["protein_uid"]).column(0).to_pylist():
                    if lo is not None and not lo <= u <= hi:
                        continue
                    if last is not None and u <= last:
                        raise RuntimeError(f"{path}: round_coverage rows are not in strictly increasing protein_uid order")
                    last = u
                    yield u
        return
    if "owned" in src:
        own_files = _checked(root, src["protein_files"], f"the plan of {_src_id(src)}")
        gi, owners = _check_owned(root, src["owned"], f"the plan of {_src_id(src)}")
        filt = "" if lo is None else " WHERE protein_uid >= $lo AND protein_uid <= $hi"
        ps_filt = "" if lo is None else " WHERE ps.protein_uid >= $lo AND ps.protein_uid <= $hi"
        parts = []
        if own_files:
            parts.append(f"SELECT protein_uid FROM read_parquet({_sql_files(own_files)}, hive_partitioning = false){filt}")
        for fp, ps in owners:
            parts.append(
                f"SELECT ps.protein_uid "
                f"FROM read_parquet({_sql_files(ps)}, hive_partitioning = false) ps "
                f"SEMI JOIN ("
                f"  SELECT DISTINCT fp.protein_uid "
                f"  FROM read_parquet({_sql_files(fp)}, hive_partitioning = false) fp "
                f"  SEMI JOIN gi USING (contig_uid)"
                f") u USING (protein_uid){ps_filt}"
            )
        sql = (f"WITH gi AS (SELECT DISTINCT contig_uid FROM read_parquet({_sql_files(gi)}, hive_partitioning = false)) "
               f"SELECT DISTINCT protein_uid FROM ({' UNION ALL '.join(parts)}) ORDER BY protein_uid")
        params = None if lo is None else {"lo": lo, "hi": hi}
        with _duck() as con:
            reader = con.sql(sql, params=params).to_arrow_reader(_READ_ROWS)
            try:
                for b in reader:
                    for u in b.column(0).to_pylist():
                        yield bytes(u)
            finally:
                reader.close()
        return
    last = None
    for path in _checked(root, src["protein_files"], f"the plan of {src['stage']}"):
        pf = pq.ParquetFile(path)
        col = pf.schema_arrow.get_field_index("protein_uid")
        for g in range(pf.metadata.num_row_groups):
            if lo is not None:
                st = pf.metadata.row_group(g).column(col).statistics
                if st is not None and st.has_min_max and (bytes(st.max) < lo or bytes(st.min) > hi):
                    continue
            for u in pf.read_row_group(g, columns=["protein_uid"]).column(0).to_pylist():
                if lo is not None and not lo <= u <= hi:
                    continue
                if last is not None and u <= last:
                    raise RuntimeError(f"{path}: protein_seq rows are not in strictly increasing protein_uid order")
                last = u
                yield u


def _measured_rows(t: pa.Table) -> pa.Table:
    """The measured rows of a row group of a source's domain file: all of them today. Card B-47c appends
    `inferred_from` (null for a measured row); inferred rows are never copied, so there this keeps the rows whose
    `inferred_from` is null (and the column is then read too)."""
    return t


def _assign(mine: Iterator[tuple[bytes, int]], measured: list[Iterator[bytes]]) -> Iterator[tuple[bytes, int, int]]:
    """(uid, length, k) for each (uid, length) of `mine`: k is the index of the first stream of `measured` that holds
    the uid, or -1. Every stream is in strictly increasing uid order; this is a merge of sorted streams that holds one
    uid of each (AGENTS §4.1), never a set."""
    heads = [next(it, None) for it in measured]
    for u, n in mine:
        found = -1
        for k, it in enumerate(measured):
            h = heads[k]
            while h is not None and h < u:
                h = next(it, None)
            heads[k] = h
            if found < 0 and h == u:
                found = k
        yield u, n, found


class _Rows:
    """The rows of sorted domain files (`l1_domains`, or a shard's run), one row group at a time, for uids asked in
    increasing order. `take(u)` gives the rows of u and passes over the rows before it, without reading a row group
    whose protein_uid maximum is below u; `until(u)` gives every row before u; `rest()` every row left. `columns`
    None reads the table's columns, conformed to its spec."""

    def __init__(self, paths: list[Path], columns: list[str] | None, measured_only: bool = False):
        def groups() -> Iterator[tuple[pq.ParquetFile, int]]:
            for path in paths:
                pf = pq.ParquetFile(path)
                for g in range(pf.metadata.num_row_groups):
                    yield pf, g

        self._groups = groups()
        self._columns, self._measured = columns, measured_only
        self._t: pa.Table | None = None
        self._uids: list[bytes] = []
        self._pos = 0

    def _load(self, skip_below: bytes | None) -> bool:
        """Load the next row group (passing over those whose maximum is below `skip_below`); False at the end."""
        spec = table("l1_domains")
        for pf, g in self._groups:
            if skip_below is not None:
                col = pf.schema_arrow.get_field_index("protein_uid")
                st = pf.metadata.row_group(g).column(col).statistics
                if st is not None and st.has_min_max and bytes(st.max) < skip_below:
                    continue
            t = pf.read_row_group(g, columns=self._columns or spec.schema.names)
            if self._measured:
                t = _measured_rows(t)
            if self._columns is None:
                t = pa.Table.from_batches([bio.conform(b, spec) for b in t.to_batches()], schema=spec.schema)
            self._t, self._uids, self._pos = t, t.column("protein_uid").to_pylist(), 0
            return True
        self._t, self._uids, self._pos = None, [], 0
        return False

    def head(self) -> bytes | None:
        """The uid of the next row, or None at the end (loading the next row group when this one is used up)."""
        while self._pos >= len(self._uids):
            if not self._load(None):
                return None
        return self._uids[self._pos]

    def take(self, u: bytes) -> list[pa.Table]:
        out = []
        while True:
            if self._pos >= len(self._uids):
                if not self._load(u):
                    return out
                continue
            i = bisect.bisect_left(self._uids, u, self._pos)
            j = bisect.bisect_right(self._uids, u, i)
            if j > i:
                out.append(self._t.slice(i, j - i))
            self._pos = j
            if j < len(self._uids):
                return out

    def until(self, u: bytes) -> list[pa.Table]:
        out = []
        while self.head() is not None:
            i = bisect.bisect_left(self._uids, u, self._pos)
            if i > self._pos:
                out.append(self._t.slice(self._pos, i - self._pos))
            self._pos = i
            if i < len(self._uids):
                return out
        return out

    def rest(self) -> list[pa.Table]:
        out = []
        while self.head() is not None:
            out.append(self._t.slice(self._pos))
            self._pos = len(self._uids)
        return out


def _walk(root: Path, files: list[Path], srcs: list[dict], columns: list[str] | None,
          owned: dict | None = None) -> Iterator[tuple[bytes, int, int, list[pa.Table]]]:
    """(uid, length, k, rows) for every protein of the stage's protein_seq `files` (and shared-contig proteins of
    `owned` when given, card C-38h), in uid order: k the index of its source in `srcs` (the first that measured it)
    or -1, and for k >= 0 that source's measured rows of the protein (`columns` of its domain files; empty when it
    had no hit)."""
    cursors = [_Rows([Path(root) / f["path"] for f in s["files"]], columns, measured_only=True) for s in srcs]
    measured = [_measured_uids(root, s) for s in srcs]
    for u, n, k in _assign(_stage_uid_lengths(root, files, owned), measured):
        yield u, n, k, (cursors[k].take(u) if k >= 0 else [])


def _shard_proteins(root: Path, p: dict, s: dict) -> dict[bytes, str]:
    """The proteins of shard `s` (uids `first` to `last`) from the plan's files (and `owned` shared-contig files,
    card C-38h), reading only the row groups whose protein_uid range meets the shard's."""
    lo, hi = bytes.fromhex(s["first"]), bytes.fromhex(s["last"])
    out: dict[bytes, str] = {}
    for f in p["files"]:
        path = Path(root) / f["path"]
        if path.stat().st_size != f["bytes"]:
            raise ValueError(f"{path}: {path.stat().st_size} bytes, the shard plan lists {f['bytes']}")
        pf = pq.ParquetFile(path)
        col = pf.schema_arrow.get_field_index("protein_uid")
        for g in range(pf.metadata.num_row_groups):
            st = pf.metadata.row_group(g).column(col).statistics
            if st is not None and st.has_min_max and (bytes(st.max) < lo or bytes(st.min) > hi):
                continue
            t = pf.read_row_group(g, columns=["protein_uid", "seq"])
            for u, seq in zip(t.column(0).to_pylist(), t.column(1).to_pylist()):
                if lo <= u <= hi:
                    out[u] = seq
    owned = p.get("owned")
    if owned:
        gi, owners = _check_owned(root, owned, "the shard plan")
        parts = [
            f"SELECT ps.protein_uid, ps.seq "
            f"FROM read_parquet({_sql_files(ps)}, hive_partitioning = false) ps "
            f"SEMI JOIN ("
            f"  SELECT DISTINCT fp.protein_uid "
            f"  FROM read_parquet({_sql_files(fp)}, hive_partitioning = false) fp "
            f"  SEMI JOIN gi USING (contig_uid)"
            f") u USING (protein_uid) "
            f"WHERE ps.protein_uid >= $lo AND ps.protein_uid <= $hi"
            for fp, ps in owners
        ]
        sql = (f"WITH gi AS (SELECT DISTINCT contig_uid FROM read_parquet({_sql_files(gi)}, hive_partitioning = false)) "
               f"SELECT protein_uid, min(seq) AS seq FROM ({' UNION ALL '.join(parts)}) GROUP BY protein_uid")
        with _duck() as con:
            reader = con.sql(sql, params={"lo": lo, "hi": hi}).to_arrow_reader(_READ_ROWS)
            try:
                for b in reader:
                    for u, seq in zip(b.column(0).to_pylist(), b.column(1).to_pylist()):
                        out[bytes(u)] = seq
            finally:
                reader.close()
    copy = p.get("copy")
    if copy and copy["sources"] and out:  # card B-47b: the proteins of the range a source measured are copied
        mine = sorted(out)
        measured = [_measured_uids(root, src, lo, hi) for src in copy["sources"]]
        for u, _, k in _assign(((u, 0) for u in mine), measured):
            if k >= 0:
                del out[u]
    return out


def _shard_inputs(root: Path, stage: str, library: str, library_version: str, hmm_path: Path, k: int,
                  out: Path | None = None) -> tuple[bytes, dict, dict, dict[bytes, str]]:
    """(the plan's bytes, the plan, its shard k, the shard's proteins), checked: k is one of the plan's shards,
    `hmm_path` is the plan's HMM file, and the proteins read are the plan's counts."""
    raw = _plan_bytes(root, stage, library, library_version, out)
    p = json.loads(raw)
    if not 0 <= k < len(p["shards"]):
        step = _names(stage, library, library_version)[1]
        raise ValueError(f"shard {k}: the plan of {step} for stage {stage} has {len(p['shards'])} shards")
    have = ids.sha256_file(Path(hmm_path)).hex()
    if have != p["hmm_sha256"]:
        raise ValueError(f"{hmm_path}: sha256 {have} is not the plan's HMM file ({p['hmm_sha256']})")
    s = p["shards"][k]
    proteins = _shard_proteins(root, p, s)
    if len(proteins) != s["proteins"] or sum(len(x) for x in proteins.values()) != s["residues"]:
        raise ValueError(f"shard {k}: {len(proteins)} proteins read, the plan lists {s['proteins']} "
                         f"({s['residues']} residues)")
    return raw, p, s, proteins


def scan_shard(root: Path, stage: str, library: str, library_version: str, hmm_path: Path, k: int, cpus: int,
               tmp: Path | None = None, out: Path | None = None) -> dict:
    """Scan shard `k` of the stage's plan with `cpus` threads: write its sorted run under scratch
    (`run_path`) and its record (the run's sha256 and counts) under `out` (default `root`). Refuses a
    shard whose record exists in `root` or `out`, an HMM file other than the plan's, and proteins that
    do not match the plan's counts."""
    root = Path(root)
    out = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out)
    if not isinstance(cpus, int) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give the thread count explicitly (at least 1)")
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, out, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    rel = _shard_rel(stage, library, library_version, k)
    for base in dict.fromkeys([root.resolve(), out.resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: shard {k} is scanned once; refusing to overwrite")
    raw, p, s, proteins = _shard_inputs(root, stage, library, library_version, hmm_path, k, out)
    rows = scan(proteins, hmm_path, p["cutoff"], cpus=cpus)
    path = run_path(tmp, root, stage, library, library_version, k, out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.unlink(missing_ok=True)  # a run left by an attempt that died before its record: scanned again
    spec = table("l1_domains")
    bio.write_table(_batches(rows, spec), spec, path)
    rec = {"format": PLAN_FORMAT, "plan_sha256": hashlib.sha256(raw).hexdigest(), "shard": k, "first": s["first"],
           "last": s["last"], "proteins": s["proteins"], "residues": s["residues"], "domains": len(rows),
           "proteins_with_hit": len({r["protein_uid"] for r in rows}), "run_sha256": ids.sha256_file(path).hex(),
           "run_bytes": path.stat().st_size}
    _write_once(root, out, rel, _json_bytes(rec), f"the record of shard {k}")
    return {"shard": k, "proteins": s["proteins"], "residues": s["residues"], "domains": len(rows),
            "proteins_with_hit": rec["proteins_with_hit"]}


# --- a shard whose scratch run was lost (card C-12) --------------------------------------------------------


def shard_record(root: Path, stage: str, library: str, library_version: str, k: int) -> dict | None:
    """The published record of shard k (in the root), or None when there is none."""
    p = _find(root, None, _shard_rel(stage, library, library_version, k))
    return None if p is None else json.loads(p.read_bytes())


def run_problem(root: Path, stage: str, library: str, library_version: str, k: int, tmp: Path | None) -> str | None:
    """Why the scratch run of shard k, whose record is published, cannot be merged as it is (`merge` would refuse
    it, or could not check it), or None: "is missing", "differs from the sha256 its record holds", or "cannot be
    checked: its record holds no run sha256" (a record written without one)."""
    rec = shard_record(root, stage, library, library_version, k)
    if rec is None:
        raise FileNotFoundError(f"shard {k} of domains-{library}-{library_version}, stage {stage}: no record")
    path = run_path(tmp, root, stage, library, library_version, k)
    if not path.is_file():
        return "is missing"
    if "run_sha256" not in rec:
        return "cannot be checked: its record holds no run sha256"
    if ids.sha256_file(path).hex() != rec["run_sha256"]:
        return "differs from the sha256 its record holds"
    return None


def rescan_shard(root: Path, stage: str, library: str, library_version: str, hmm_path: Path, k: int, cpus: int,
                 tmp: Path | None = None) -> dict:
    """Scan shard k again into scratch (`run_path`) for a shard whose record is published in `root` but whose run
    is gone or damaged, e.g. after a reboot cleared /tmp (card C-12). The record is not written again: the new run
    must reproduce it. Its sha256 must equal the record's `run_sha256`; a record that holds none (written without
    one) is checked by its domain and hit counts instead, and the result says `"verified": "counts"`. On a mismatch
    the new run is deleted, the run in scratch is left as it was, and a ValueError names the shard. Refuses, like
    `scan_shard`, when the stage's domains are written, the HMM file is not the plan's, or the proteins do not
    match the plan."""
    root = Path(root)
    ledger.check_scratch(tmp, root)
    if not isinstance(cpus, int) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give the thread count explicitly (at least 1)")
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, root, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    rec = shard_record(root, stage, library, library_version, k)
    if rec is None:
        raise FileNotFoundError(f"shard {k} of {step}, stage {stage}: no published record to scan it again for; "
                                "scan it with domains.scan_shard")
    raw, p, s, proteins = _shard_inputs(root, stage, library, library_version, hmm_path, k)
    if rec["plan_sha256"] != hashlib.sha256(raw).hexdigest() or rec["shard"] != k:
        raise ValueError(f"shard {k} of {step}, stage {stage}: its record is of another plan or shard")
    rows = scan(proteins, hmm_path, p["cutoff"], cpus=cpus)
    path = run_path(tmp, root, stage, library, library_version, k)
    path.parent.mkdir(parents=True, exist_ok=True)
    new = path.with_name(f".rescan-{os.getpid()}-{path.name}")
    spec = table("l1_domains")
    bio.write_table(_batches(rows, spec), spec, new)
    have, hits = ids.sha256_file(new).hex(), len({r["protein_uid"] for r in rows})
    want = rec.get("run_sha256")
    if want is not None:
        bad = None if have == want else f"its run has sha256 {have}, the record holds {want}"
    else:
        got, rec_counts = (len(rows), hits), (rec["domains"], rec["proteins_with_hit"])
        bad = None if got == rec_counts else (f"it gives {got[0]} domains in {got[1]} proteins, the record "
                                              f"{rec_counts[0]} in {rec_counts[1]}")
    if bad is not None:
        new.unlink()
        raise ValueError(f"shard {k} of {step}, stage {stage}: scanned again, {bad}; the published record is not "
                         f"the scan of these inputs, so the merge cannot run (the run in scratch is left as it was)")
    os.replace(new, path)
    return {"shard": k, "domains": len(rows), "proteins_with_hit": hits, "run_sha256": have,
            "verified": "sha256" if want is not None else "counts"}


def merge(root: Path, stage: str, library: str, library_version: str, tmp: Path | None = None,
          max_file_bytes: int = layout.MAX_FILE_BYTES, out: Path | None = None) -> dict:
    """Write the shards' runs, in shard order, as the stage's domain files and ledger under `out`
    (default `root`): the files `run` writes. Every shard needs its record (in `root` or `out`) and a run
    whose sha256 is the record's (for a record that holds no run sha256, whose row count is the record's
    domains; build-stage scans such a shard again first, `rescan_shard`); otherwise nothing is written."""
    root = Path(root)
    out_root = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out_root)
    lib_dir, step = _names(stage, library, library_version)
    ledger.refuse_written(root, out_root, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")
    raw = _plan_bytes(root, stage, library, library_version, out)
    p, digest = json.loads(raw), hashlib.sha256(raw).hexdigest()
    missing = [k for k in range(len(p["shards"])) if _find(root, out, _shard_rel(stage, library, library_version, k)) is None]
    if missing:
        raise FileNotFoundError(f"{step} of stage {stage}: no record of shard {', '.join(map(str, missing))} "
                                f"of {len(p['shards'])}; scan it with domains.scan_shard")
    recs, runs = [], []
    for k in range(len(p["shards"])):
        rec = json.loads(_find(root, out, _shard_rel(stage, library, library_version, k)).read_bytes())
        if rec["plan_sha256"] != digest or rec["shard"] != k:
            raise ValueError(f"shard {k}: its record is of another plan or shard")
        path = run_path(tmp, root, stage, library, library_version, k, out)
        if not path.is_file():
            raise FileNotFoundError(f"{path}: the run of shard {k} is not in the scratch directory; its record "
                                    "is published, so scan the shard again into this scratch directory")
        if "run_sha256" not in rec:  # a record written without one: its row count is checked (card C-12)
            if pq.ParquetFile(path).metadata.num_rows != rec["domains"]:
                raise ValueError(f"{path}: row count differs from the domains in the record of shard {k}")
        elif ids.sha256_file(path).hex() != rec["run_sha256"]:
            raise ValueError(f"{path}: sha256 differs from the record of shard {k}")
        recs.append(rec)
        runs.append(path)
    models: set[str] = set()
    if "copy" in p:  # card B-47b
        return _merge_copies(root, out_root, stage, lib_dir, step, p, recs, runs, max_file_bytes)

    def batches() -> Iterator[pa.RecordBatch]:
        for path in runs:
            for b in pq.ParquetFile(path).iter_batches(batch_size=bio.ROW_GROUP_ROWS):
                models.update(pc.unique(b.column("model_acc")).to_pylist())
                yield b

    hit = sum(r["proteins_with_hit"] for r in recs)
    _write(out_root, stage, lib_dir, step, batches(), p["proteins"], hit, max_file_bytes)
    return {"proteins": p["proteins"], "domains": sum(r["domains"] for r in recs), "proteins_with_hit": hit,
            "models_with_hit": len(models), "shards": len(recs)}


def _merge_copies(root: Path, out_root: Path, stage: str, lib_dir: str, step: str, p: dict, recs: list[dict],
                  runs: list[Path], max_file_bytes: int) -> dict:
    """`merge` of a plan that copies (card B-47b): the shards' rows and the copied rows in uid order, as one sorted
    stream, written as `run` writes a full scan. Before anything is written: this process's tool versions are the
    plan's, every file of the stage and of its sources has the bytes the plan lists, and every domain file copied from
    has the plan's sha256. While writing: every row (copied or scanned) has this scan's cutoff, Z and domZ, and no
    protein is both copied and scanned. At the end of the stream, before its last rows are given to the writer: the copied proteins,
    residues, rows and proteins with a hit of each source are the plan's (else a ValueError; under run-step nothing is
    then published)."""
    copy, spec = p["copy"], table("l1_domains")
    if _tools() != p["tools"]:
        raise ValueError(f"{step} of stage {stage}: this process has tool versions {_tools()}, the plan was made and "
                         f"its sources compared for {p['tools']}")
    files = _checked(root, p["files"], "the shard plan")
    if "owned" in p:
        _check_owned(root, p["owned"], "the shard plan")
    srcs = copy["sources"]
    for s in srcs:
        _checked(root, s["files"], f"the shard plan (domain files of {_src_id(s)})")
        if "coverage_files" in s:
            _checked(root, s["coverage_files"], f"the shard plan (round_coverage of {_src_id(s)})")
        else:
            _checked(root, s["protein_files"], f"the shard plan (proteins of {_src_id(s)})")
        if "owned" in s:
            _check_owned(root, s["owned"], f"the shard plan (owned proteins of {_src_id(s)})")
    counts = [{"proteins": 0, "residues": 0, "rows": 0, "proteins_with_hit": 0} for _ in srcs]
    models: set[str] = set()
    scanned = _Rows(runs, None)
    want = {"cutoff": p["cutoff"], "z": Z_PROT, "dom_z": DOM_Z}

    def flush(pend: list[pa.Table]) -> pa.Table:
        """The pieces as one table (one chunk per column), every row checked: copied rows carry their source's
        scan, which is this one, and scanned rows this plan's."""
        t = pa.concat_tables(pend).combine_chunks()
        for col, v in want.items():
            if not pc.all(pc.equal(t.column(col), v)).as_py():
                raise ValueError(f"{step} of stage {stage}: a row to write has another {col} than {v!r}")
        models.update(pc.unique(t.column("model_acc")).to_pylist())
        return t

    def batches() -> Iterator[pa.Table]:
        # each piece is a slice (a table object of its own, about 7 KB of Arrow metadata with 17 columns), so they
        # are combined every FLUSH_PIECES pieces as well as every row group's worth of rows (AGENTS §4.1: a slice
        # per copied protein, held until a row group fills, grew 546 MiB from 8,540 to 85,400 proteins: *measured*,
        # wilkcomp02, 2026-11-15, tests/r2/test_b47b_copy_memory.py before this rule)
        pend: list[pa.Table] = []
        held = 0
        for u, n, k, rows in _walk(root, files, srcs, None, owned=p.get("owned")):
            if k < 0:
                continue
            c = counts[k]
            c["proteins"] += 1
            c["residues"] += n
            if not rows:
                continue
            c["rows"] += sum(t.num_rows for t in rows)
            c["proteins_with_hit"] += 1
            before = scanned.until(u)
            if scanned.head() == u:
                raise ValueError(f"protein {u.hex()} is copied from {_src_id(srcs[k])} and was also scanned")
            pend += before + rows
            held += sum(t.num_rows for t in before) + sum(t.num_rows for t in rows)
            if held >= bio.ROW_GROUP_ROWS or len(pend) >= FLUSH_PIECES:
                yield flush(pend)
                pend, held = [], 0
        pend += scanned.rest()
        bad = [(_src_id(s), {k: s[k] for k in c}, c) for s, c in zip(srcs, counts) if {k: s[k] for k in c} != c]
        if bad:
            raise ValueError(f"{step} of stage {stage}: copied counts differ from the plan's: {bad}")
        if pend:
            yield flush(pend)

    hit = sum(r["proteins_with_hit"] for r in recs) + copy["proteins_with_hit"]
    _write(out_root, stage, lib_dir, step, batches(), p["proteins"], hit, max_file_bytes)
    return {"proteins": p["proteins"], "domains": sum(r["domains"] for r in recs) + copy["rows"],
            "proteins_with_hit": hit, "models_with_hit": len(models), "shards": len(recs),
            "copied_proteins": copy["proteins"], "copied_rows": copy["rows"],
            "sources": [_src_id(s) for s in srcs]}


# --- annotation rounds (card C-38f, DEC-54, planning/design/C-38_import_depths.md §3.2, §3.3, §3.6) ---------

EST_COV_ROW_BYTES = 24  # estimate: zstd bytes per round_coverage row (protein_uid + uint32 n_rows); cuts parts only
_REQ_SCHEMA = pa.schema([pa.field("protein_uid", pa.binary(16), nullable=False),
                         pa.field("length", pa.uint32(), nullable=False)])


def round_id(library: str, library_version: str, hmm_sha256: str, cutoff: str, protein_uids,
             *, z: float = Z_PROT, dom_z: float = DOM_Z) -> str:
    """Content-derived round id `r-<16 hex>` (card C-38f, §3.2): sha256 over library name and version, HMM
    sha256, cutoff, Z, domZ, and the sorted distinct requested protein uids. Accepts a sequence/set (sorted
    and deduplicated here) or a streaming iterator of strictly increasing 16-byte `protein_uid`s."""
    h = hashlib.sha256()
    head = json.dumps({"library": library, "library_version": library_version, "hmm_sha256": hmm_sha256,
                       "cutoff": cutoff, "z": z, "dom_z": dom_z}, sort_keys=True).encode()
    h.update(head)
    h.update(b"\n")
    if isinstance(protein_uids, (list, tuple, set, frozenset)):
        stream = sorted(set(protein_uids))
    else:
        stream = protein_uids
    last: bytes | None = None
    for u in stream:
        b = bytes(u)
        if len(b) != 16:
            raise ValueError(f"round_id: protein_uid must be 16 bytes, got {len(b)}")
        if last is not None and b <= last:
            raise ValueError("round_id: streaming protein_uids must be strictly increasing")
        last = b
        h.update(b)
    return f"r-{h.hexdigest()[:16]}"


def round_scan_dir(rid: str, library: str, library_version: str) -> str:
    """Root-relative directory of the `round_scan` artifact for one round and library."""
    return f"{ARTIFACTS['round_scan'].location}/round={rid}/library={library}-{library_version}"


def round_domains_dir(rid: str, library: str, library_version: str) -> str:
    return f"{table('round_domains').location}/round={rid}/library={library}-{library_version}"


def round_coverage_dir(rid: str, library: str, library_version: str) -> str:
    return f"{table('round_coverage').location}/round={rid}/library={library}-{library_version}"


def round_ledger_dir(rid: str, library: str, library_version: str) -> str:
    return f"ledger/round={rid}/step=domains-{library}-{library_version}"


def _round_plan_rel(rid: str, library: str, library_version: str) -> str:
    return f"{round_scan_dir(rid, library, library_version)}/plan.json"


def _round_shard_rel(rid: str, library: str, library_version: str, k: int) -> str:
    return f"{round_scan_dir(rid, library, library_version)}/shard-{k:05d}.json"


def round_scan_ids(root: Path, library: str, library_version: str) -> list[str]:
    """Sorted round ids on disk under `root` that have a `round_scan` plan for `library`-`library_version`."""
    base = Path(root) / ARTIFACTS["round_scan"].location
    if not base.is_dir():
        return []
    lib = f"{library}-{library_version}"
    out = []
    for p in base.glob(f"round=*/library={lib}/plan.json"):
        if p.is_file():
            rid = p.parent.parent.name.removeprefix("round=")
            if rid:
                out.append(rid)
    return sorted(out)


def round_complete(root: Path, rid: str, library: str, library_version: str) -> bool:
    """Whether round `rid` for `library`-`library_version` has complete outputs (`round_domains`, `round_coverage`,
    and its ledger) under `root` and is not mid-publish."""
    root = Path(root)
    step = f"domains-{library}-{library_version}"
    if (root / f"runs/stage={rid}/{step}.publishing.json").exists():
        return False
    return bool(
        bio.part_files(root / round_domains_dir(rid, library, library_version))
        and bio.part_files(root / round_coverage_dir(rid, library, library_version))
        and bio.part_files(root / round_ledger_dir(rid, library, library_version))
    )


def _round_plan_bytes(root: Path, rid: str, library: str, library_version: str, out: Path | None = None) -> bytes:
    p = _find(root, out, _round_plan_rel(rid, library, library_version))
    if p is None:
        raise FileNotFoundError(
            f"{Path(root) / _round_plan_rel(rid, library, library_version)}: no round plan for "
            f"round {rid} domains-{library}-{library_version}; write it with domains.plan_round"
        )
    return p.read_bytes()


def read_round_plan(root: Path, rid: str, library: str, library_version: str, out: Path | None = None) -> dict:
    """The round's shard plan for one library (from `root`, else `out`)."""
    return json.loads(_round_plan_bytes(root, rid, library, library_version, out))


def _round_scratch_dir(tmp: Path | None, raw_plan: bytes) -> Path:
    digest = hashlib.sha256(raw_plan).hexdigest()[:16]
    return Path(tmp if tmp is not None else tempfile.gettempdir()) / SCAN_DIR / digest


def round_run_path(tmp: Path | None, root: Path, rid: str, library: str, library_version: str, k: int,
                   out: Path | None = None) -> Path:
    """Where `scan_round_shard` writes shard k's sorted run under scratch."""
    raw = _round_plan_bytes(root, rid, library, library_version, out)
    return _round_scratch_dir(tmp, raw) / f"shard-{k:05d}.parquet"


def round_requested_path(tmp: Path | None, root: Path, rid: str, library: str, library_version: str,
                         out: Path | None = None) -> Path:
    """Where `plan_round` writes the sorted `(protein_uid, length)` requested table in scratch."""
    raw = _round_plan_bytes(root, rid, library, library_version, out)
    return _round_scratch_dir(tmp, raw) / "requested.parquet"


def _write_input_uids_parquet(proteins_source: Path | list[bytes] | set[bytes], dst: Path) -> int:
    """Stream protein uids from a text file (one 32-hex uid per line, blank/# skipped), a `.parquet` file
    (column `protein_uid`), or an iterable of 16-byte uids into a temporary Parquet file at `dst`. Returns
    the raw count written (before deduplication). Bounded memory (AGENTS §4.1)."""
    schema = pa.schema([pa.field("protein_uid", pa.binary(16), nullable=False)])
    total = 0
    with pq.ParquetWriter(dst, schema, compression="zstd") as writer:
        if isinstance(proteins_source, (str, Path)):
            p = Path(proteins_source)
            if p.suffix == ".parquet":
                pf = pq.ParquetFile(p)
                for b in pf.iter_batches(batch_size=_READ_ROWS, columns=["protein_uid"]):
                    if b.num_rows:
                        writer.write_batch(pa.RecordBatch.from_arrays([b.column(0).cast(pa.binary(16))], schema=schema))
                        total += b.num_rows
            else:
                buf: list[bytes] = []
                with open(p, encoding="utf-8") as fh:
                    for lineno, line in enumerate(fh, 1):
                        s = line.strip()
                        if not s or s.startswith("#"):
                            continue
                        if len(s) != 32:
                            raise ValueError(f"{p}:{lineno}: {s!r} is not a protein uid (32 hex digits)")
                        try:
                            buf.append(bytes.fromhex(s))
                        except ValueError as e:
                            raise ValueError(f"{p}:{lineno}: {s!r} is not a protein uid (32 hex digits)") from e
                        if len(buf) >= _READ_ROWS:
                            writer.write_batch(pa.RecordBatch.from_arrays([pa.array(buf, pa.binary(16))], schema=schema))
                            total += len(buf)
                            buf = []
                if buf:
                    writer.write_batch(pa.RecordBatch.from_arrays([pa.array(buf, pa.binary(16))], schema=schema))
                    total += len(buf)
        else:
            buf = []
            for u in proteins_source:
                b = bytes(u)
                if len(b) != 16:
                    raise ValueError(f"protein_uid must be 16 bytes, got {len(b)}")
                buf.append(b)
                if len(buf) >= _READ_ROWS:
                    writer.write_batch(pa.RecordBatch.from_arrays([pa.array(buf, pa.binary(16))], schema=schema))
                    total += len(buf)
                    buf = []
            if buf:
                writer.write_batch(pa.RecordBatch.from_arrays([pa.array(buf, pa.binary(16))], schema=schema))
                total += len(buf)
    return total


def _filter_unscanned_to_parquet(
    root: Path,
    release: str,
    gene_call_version: str,
    library: str,
    library_version: str,
    proteins_source: Path | list[bytes] | set[bytes],
    out_parquet: Path,
) -> tuple[int, int, list[Path]]:
    """Filter `proteins_source` to proteins whose `domain_state` in `release` for `library`-`library_version`
    is `not scanned` (card C-38f §3.3, §3.6), writing sorted `(protein_uid, length)` to `out_parquet` in
    bounded memory via single-thread DuckDB (`DUCK_MEMORY_LIMIT`).

    Returns `(n_input_distinct, n_unscanned, all_protein_seq_files)`. Raises ValueError if any input protein
    is not present in `release`'s `protein_seq`."""
    from biobox import catalog

    files = catalog.ReleaseFiles.of(root, release)
    all_ps = files.table("protein_seq", gene_call_version=gene_call_version)
    if not all_ps:
        raise FileNotFoundError(f"release {release}: no protein_seq files for gene_call_version={gene_call_version}")
    lib = f"{library}-{library_version}"
    dom_files = files.table("l1_domains", library=lib)
    ann_stages = {parts["stage"] for _, parts in files.entries("l1_domains", library=lib) if "stage" in parts}
    ann_ps = [p for p, parts in files.entries("protein_seq", gene_call_version=gene_call_version)
              if parts.get("stage") in ann_stages]
    rcov_files = files.table("round_coverage", library=lib)

    in_pq = out_parquet.with_name(f".in-{os.getpid()}.parquet")
    try:
        raw_count = _write_input_uids_parquet(proteins_source, in_pq)
        if raw_count == 0:
            pq.write_table(pa.Table.from_pylist([], schema=_REQ_SCHEMA), out_parquet, compression="zstd")
            return 0, 0, all_ps
        with _duck() as con:
            in_sql = _sql_files([in_pq])
            ps_sql = _sql_files(all_ps)
            (n_in_distinct,), = con.execute(
                f"SELECT count(DISTINCT protein_uid)::BIGINT FROM read_parquet({in_sql}, hive_partitioning = false)"
            ).fetchall()
            (n_missing,), = con.execute(
                f"SELECT count(*)::BIGINT FROM ("
                f"  SELECT DISTINCT protein_uid FROM read_parquet({in_sql}, hive_partitioning = false)"
                f") r ANTI JOIN ("
                f"  SELECT protein_uid FROM read_parquet({ps_sql}, hive_partitioning = false)"
                f") ps USING (protein_uid)"
            ).fetchall()
            if n_missing:
                raise ValueError(
                    f"release {release}: {n_missing} requested protein(s) are not in the release's protein_seq"
                )
            settled_parts = []
            if dom_files:
                settled_parts.append(
                    f"SELECT protein_uid FROM read_parquet({_sql_files(dom_files)}, hive_partitioning = false)"
                )
            if ann_ps:
                settled_parts.append(
                    f"SELECT protein_uid FROM read_parquet({_sql_files(ann_ps)}, hive_partitioning = false)"
                )
            if rcov_files:
                settled_parts.append(
                    f"SELECT protein_uid FROM read_parquet({_sql_files(rcov_files)}, hive_partitioning = false)"
                )
            if settled_parts:
                settled_cte = f"settled AS ({' UNION ALL '.join(settled_parts)}), "
                req_sql = (
                    f"req AS ("
                    f"  SELECT DISTINCT r.protein_uid FROM read_parquet({in_sql}, hive_partitioning = false) r "
                    f"  ANTI JOIN settled USING (protein_uid)"
                    f")"
                )
            else:
                settled_cte = ""
                req_sql = f"req AS (SELECT DISTINCT protein_uid FROM read_parquet({in_sql}, hive_partitioning = false))"
            sql = (
                f"WITH {settled_cte}{req_sql} "
                f"SELECT ps.protein_uid, min(ps.length)::UINTEGER AS length "
                f"FROM read_parquet({ps_sql}, hive_partitioning = false) ps "
                f"SEMI JOIN req USING (protein_uid) "
                f"GROUP BY ps.protein_uid ORDER BY ps.protein_uid"
            )
            n_unscanned = 0
            reader = con.sql(sql).to_arrow_reader(_READ_ROWS)
            try:
                with pq.ParquetWriter(out_parquet, _REQ_SCHEMA, compression="zstd") as writer:
                    for b in reader:
                        if b.num_rows:
                            tb = pa.RecordBatch.from_arrays(
                                [b.column(0).cast(pa.binary(16)), b.column(1).cast(pa.uint32())],
                                schema=_REQ_SCHEMA,
                            )
                            writer.write_batch(tb)
                            n_unscanned += tb.num_rows
            finally:
                reader.close()
        return int(n_in_distinct), int(n_unscanned), all_ps
    finally:
        in_pq.unlink(missing_ok=True)


def filter_unscanned(
    root: Path,
    release: str,
    gene_call_version: str,
    library: str,
    library_version: str,
    proteins_source: Path | list[bytes] | set[bytes],
    *,
    tmp: Path | None = None,
) -> tuple[list[bytes], int]:
    """Return `(sorted unscanned protein_uids, n_input_distinct)` for `proteins_source` in `release`."""
    d = Path(tempfile.mkdtemp(prefix="biobox-round-filter-", dir=tmp))
    try:
        req_pq = d / "filtered.parquet"
        n_in, _, _ = _filter_unscanned_to_parquet(
            Path(root), release, gene_call_version, library, library_version, proteins_source, req_pq
        )
        uids = pq.read_table(req_pq, columns=["protein_uid"]).column(0).to_pylist()
        return [bytes(u) for u in uids], n_in
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _req_uid_lengths(req_pq: Path) -> Iterator[tuple[bytes, int]]:
    """Stream `(protein_uid, length)` from a sorted `requested.parquet` file."""
    last = None
    for b in pq.ParquetFile(req_pq).iter_batches(batch_size=_READ_ROWS, columns=["protein_uid", "length"]):
        for u, n in zip(b.column(0).to_pylist(), b.column(1).to_pylist()):
            if last is not None and u <= last:
                raise RuntimeError(f"{req_pq}: requested rows are not in strictly increasing protein_uid order")
            last = u
            yield u, int(n)


def prepare_round(
    root: Path,
    release: str,
    gene_call_version: str,
    library: str,
    library_version: str,
    hmm_path: Path,
    cutoff: str,
    proteins_file: Path,
    *,
    tmp: Path | None = None,
) -> dict:
    """Compute a round's content-derived `round_id`, filtered protein/residue counts, and whether it is already
    complete on disk, without writing anything under `root`."""
    root = Path(root)
    hmm_sha256 = ids.sha256_file(Path(hmm_path)).hex()
    d = Path(tempfile.mkdtemp(prefix="biobox-round-prep-", dir=tmp))
    try:
        req_pq = d / "requested.parquet"
        n_in, n_unscanned, _ = _filter_unscanned_to_parquet(
            root, release, gene_call_version, library, library_version, Path(proteins_file), req_pq
        )
        if n_unscanned == 0:
            return {"round_id": None, "input_proteins": n_in, "proteins": 0, "residues": 0,
                    "settled_in_release": n_in, "already_complete": False}
        total_res = 0

        def uids() -> Iterator[bytes]:
            nonlocal total_res
            for u, n in _req_uid_lengths(req_pq):
                total_res += n
                yield u

        rid = round_id(library, library_version, hmm_sha256, cutoff, uids())
        return {
            "round_id": rid,
            "input_proteins": n_in,
            "proteins": n_unscanned,
            "residues": total_res,
            "settled_in_release": n_in - n_unscanned,
            "already_complete": round_complete(root, rid, library, library_version),
        }
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _walk_requested(root: Path, req_pq: Path, srcs: list[dict], columns: list[str] | None) -> Iterator[
    tuple[bytes, int, int, list[pa.Table]]
]:
    """Stream `(uid, length, k, rows)` for every protein of `req_pq` against `srcs` in uid order."""
    cursors = [_Rows([Path(root) / f["path"] for f in s["files"]], columns, measured_only=True) for s in srcs]
    measured = [_measured_uids(root, s) for s in srcs]
    for u, n, k in _assign(_req_uid_lengths(req_pq), measured):
        yield u, n, k, (cursors[k].take(u) if k >= 0 else [])


def plan_round(
    root: Path,
    release: str,
    gene_call_version: str,
    library: str,
    library_version: str,
    hmm_path: Path,
    cutoff: str,
    proteins_file: Path,
    shard_residues: int,
    *,
    out: Path | None = None,
    seq_stages: list[str] | None = None,
    tmp: Path | None = None,
) -> dict:
    """Plan an annotation round (card C-38f §3.3, §3.6): filter `proteins_file` down to proteins whose
    `domain_state` in `release` for `library`-`library_version` is `not scanned`, compute the content-derived
    `round_id`, leave out proteins already measured on disk by any stage (including stages built after `release`)
    or earlier round with equal `COMPARED` parameters, cut the remaining proteins into shards of at most
    `shard_residues` residues, write `requested.parquet` in scratch (`tmp`), and publish `plan.json` under
    `out` (default `root`)."""
    from biobox.l0seq import store

    root = Path(root)
    out_root = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out_root)
    if not isinstance(shard_residues, int) or shard_residues < 1:
        raise ValueError(f"shard_residues {shard_residues!r}: must be a positive integer")
    if cutoff not in ("gathering", "trusted", "noise", "none"):
        raise ValueError(f"cutoff {cutoff!r}: expected gathering, trusted, noise or none")
    hmm_sha256 = ids.sha256_file(Path(hmm_path)).hex()
    tmp_base = Path(tmp) if tmp is not None else Path(tempfile.gettempdir())
    tmp_base.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="biobox-round-plan-", dir=tmp_base))
    try:
        req_tmp = work / "requested.parquet"
        n_in, n_unscanned, all_ps = _filter_unscanned_to_parquet(
            root, release, gene_call_version, library, library_version, Path(proteins_file), req_tmp
        )
        if n_unscanned == 0:
            return {"round_id": None, "input_proteins": n_in, "proteins": 0, "residues": 0,
                    "settled_in_release": n_in, "shards": 0, "scanned_proteins": 0, "scanned_residues": 0,
                    "copied_proteins": 0, "copied_residues": 0, "copied_rows": 0, "sources": [], "passed_over": []}
        rid = round_id(library, library_version, hmm_sha256, cutoff, (u for u, _ in _req_uid_lengths(req_tmp)))
        lib = f"{library}-{library_version}"
        dom_dir = round_domains_dir(rid, library, library_version)
        cov_dir = round_coverage_dir(rid, library, library_version)
        led_dir = round_ledger_dir(rid, library, library_version)
        ledger.refuse_written(root, out_root, [dom_dir, cov_dir, led_dir], f"round {rid} domains-{lib}")
        rel = _round_plan_rel(rid, library, library_version)
        for base in dict.fromkeys([root.resolve(), out_root.resolve()]):
            if (base / rel).exists():
                raise FileExistsError(f"{base / rel}: the round plan is written once; refusing to overwrite")

        cand_stages = (
            list(dict.fromkeys(seq_stages))
            if seq_stages is not None
            else sorted(set(store.stages_on_disk(root)) | set(plan_stages(root, library, library_version)))
        )
        srcs, passed = sources(root, "", library, library_version, hmm_sha256, cutoff, cand_stages,
                               include_rounds=True, exclude_round=rid)
        total = {"proteins": 0, "residues": 0}
        counts = [{"proteins": 0, "residues": 0, "rows": 0, "proteins_with_hit": 0} for _ in srcs]

        def scanned() -> Iterator[tuple[bytes, int]]:
            for u, n, k, rows in _walk_requested(root, req_tmp, srcs, ["protein_uid"]):
                total["proteins"] += 1
                total["residues"] += n
                if k < 0:
                    yield u, n
                    continue
                c = counts[k]
                c["proteins"] += 1
                c["residues"] += n
                got = sum(t.num_rows for t in rows)
                c["rows"] += got
                c["proteins_with_hit"] += got > 0

        shards = list(_cut(scanned(), shard_residues))
        src_docs = [{**s, **c} for s, c in zip(srcs, counts)]
        copy = {
            "seq_stages": cand_stages,
            "compared": list(COMPARED),
            "sources": src_docs,
            "passed_over": passed,
            **{k: sum(c[k] for c in counts) for k in ("proteins", "residues", "rows", "proteins_with_hit")},
        }
        doc = {
            "format": PLAN_FORMAT,
            "round": rid,
            "release": release,
            "gene_call_version": gene_call_version,
            "library": library,
            "library_version": library_version,
            "cutoff": cutoff,
            "hmm_sha256": hmm_sha256,
            "shard_residues": shard_residues,
            "protein_files": _file_entries(root, all_ps),
            "requested_sha256": ids.sha256_file(req_tmp).hex(),
            "requested_bytes": req_tmp.stat().st_size,
            "input_proteins": n_in,
            "settled_in_release": n_in - n_unscanned,
            "proteins": total["proteins"],
            "residues": total["residues"],
            "shards": shards,
            "z": Z_PROT,
            "dom_z": DOM_Z,
            "tools": _tools(),
            "copy": copy,
        }
        raw = _json_bytes(doc)
        sdir = _round_scratch_dir(tmp_base, raw)
        sdir.mkdir(parents=True, exist_ok=True)
        os.replace(req_tmp, sdir / "requested.parquet")
        _write_once(root, out_root, rel, raw, "the round plan")
        scanned_p, scanned_r = sum(s["proteins"] for s in shards), sum(s["residues"] for s in shards)
        return {
            "round_id": rid,
            "shards": len(shards),
            "input_proteins": n_in,
            "settled_in_release": n_in - n_unscanned,
            "proteins": doc["proteins"],
            "residues": doc["residues"],
            "scanned_proteins": scanned_p,
            "scanned_residues": scanned_r,
            "copied_proteins": copy["proteins"],
            "copied_residues": copy["residues"],
            "copied_rows": copy["rows"],
            "sources": [_src_id(s) for s in srcs],
            "passed_over": [_src_id(p) for p in passed],
        }
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _round_shard_proteins(root: Path, p: dict, s: dict, req_path: Path) -> dict[bytes, str]:
    """The proteins to scan in round shard `s` (uids `first`..`last` from `req_path`, excluding copied uids),
    reading sequences from `p['protein_files']` only in row groups meeting `[lo, hi]`."""
    lo, hi = bytes.fromhex(s["first"]), bytes.fromhex(s["last"])
    need: set[bytes] = set()
    pf = pq.ParquetFile(req_path)
    col = pf.schema_arrow.get_field_index("protein_uid")
    for g in range(pf.metadata.num_row_groups):
        st = pf.metadata.row_group(g).column(col).statistics
        if st is not None and st.has_min_max and (bytes(st.max) < lo or bytes(st.min) > hi):
            continue
        for u in pf.read_row_group(g, columns=["protein_uid"]).column(0).to_pylist():
            if lo <= u <= hi:
                need.add(u)
    copy = p.get("copy")
    if copy and copy["sources"] and need:
        mine = sorted(need)
        measured = [_measured_uids(root, src, lo, hi) for src in copy["sources"]]
        for u, _, k in _assign(((u, 0) for u in mine), measured):
            if k >= 0:
                need.discard(u)
    if not need:
        return {}
    out: dict[bytes, str] = {}
    for path in _checked(root, p["protein_files"], "the round plan"):
        psf = pq.ParquetFile(path)
        pcol = psf.schema_arrow.get_field_index("protein_uid")
        for g in range(psf.metadata.num_row_groups):
            st = psf.metadata.row_group(g).column(pcol).statistics
            if st is not None and st.has_min_max and (bytes(st.max) < lo or bytes(st.min) > hi):
                continue
            t = psf.read_row_group(g, columns=["protein_uid", "seq"])
            for u, seq in zip(t.column(0).to_pylist(), t.column(1).to_pylist()):
                if u in need:
                    out[u] = seq
    return out


def _round_shard_inputs(root: Path, rid: str, library: str, library_version: str, hmm_path: Path, k: int,
                        tmp: Path | None = None, out: Path | None = None) -> tuple[bytes, dict, dict, dict[bytes, str]]:
    raw = _round_plan_bytes(root, rid, library, library_version, out)
    p = json.loads(raw)
    if not 0 <= k < len(p["shards"]):
        raise ValueError(f"shard {k}: the plan of round {rid} domains-{library}-{library_version} has "
                         f"{len(p['shards'])} shards")
    have = ids.sha256_file(Path(hmm_path)).hex()
    if have != p["hmm_sha256"]:
        raise ValueError(f"{hmm_path}: sha256 {have} is not the round plan's HMM file ({p['hmm_sha256']})")
    req_path = _round_scratch_dir(tmp, raw) / "requested.parquet"
    if not req_path.is_file():
        raise FileNotFoundError(f"{req_path}: requested protein list for round {rid} is missing from scratch")
    if ids.sha256_file(req_path).hex() != p["requested_sha256"]:
        raise ValueError(f"{req_path}: sha256 differs from the round plan's requested_sha256")
    s = p["shards"][k]
    proteins = _round_shard_proteins(root, p, s, req_path)
    if len(proteins) != s["proteins"] or sum(len(x) for x in proteins.values()) != s["residues"]:
        raise ValueError(f"shard {k}: {len(proteins)} proteins read, the round plan lists {s['proteins']} "
                         f"({s['residues']} residues)")
    return raw, p, s, proteins


def scan_round_shard(root: Path, rid: str, library: str, library_version: str, hmm_path: Path, k: int, cpus: int,
                     *, tmp: Path | None = None, out: Path | None = None) -> dict:
    """Scan shard `k` of round `rid` with `cpus` threads; write its sorted run under scratch and its record
    under `out` (default `root`)."""
    root = Path(root)
    out_root = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out_root)
    if not isinstance(cpus, int) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give the thread count explicitly (at least 1)")
    lib = f"{library}-{library_version}"
    ledger.refuse_written(
        root, out_root,
        [round_domains_dir(rid, library, library_version),
         round_coverage_dir(rid, library, library_version),
         round_ledger_dir(rid, library, library_version)],
        f"round {rid} domains-{lib}",
    )
    rel = _round_shard_rel(rid, library, library_version, k)
    for base in dict.fromkeys([root.resolve(), out_root.resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: round shard {k} is scanned once; refusing to overwrite")
    raw, p, s, proteins = _round_shard_inputs(root, rid, library, library_version, hmm_path, k, tmp=tmp, out=out)
    rows = scan(proteins, hmm_path, p["cutoff"], cpus=cpus)
    path = _round_scratch_dir(tmp, raw) / f"shard-{k:05d}.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.unlink(missing_ok=True)
    spec = table("round_domains")
    bio.write_table(_batches(rows, spec), spec, path)
    rec = {
        "format": PLAN_FORMAT, "round": rid, "plan_sha256": hashlib.sha256(raw).hexdigest(), "shard": k,
        "first": s["first"], "last": s["last"], "proteins": s["proteins"], "residues": s["residues"],
        "domains": len(rows), "proteins_with_hit": len({r["protein_uid"] for r in rows}),
        "run_sha256": ids.sha256_file(path).hex(), "run_bytes": path.stat().st_size,
    }
    _write_once(root, out_root, rel, _json_bytes(rec), f"the record of round shard {k}")
    return {"round_id": rid, "shard": k, "proteins": s["proteins"], "residues": s["residues"],
            "domains": len(rows), "proteins_with_hit": rec["proteins_with_hit"]}


def round_shard_record(root: Path, rid: str, library: str, library_version: str, k: int) -> dict | None:
    p = _find(root, None, _round_shard_rel(rid, library, library_version, k))
    return None if p is None else json.loads(p.read_bytes())


def round_run_problem(root: Path, rid: str, library: str, library_version: str, k: int, tmp: Path | None) -> str | None:
    rec = round_shard_record(root, rid, library, library_version, k)
    if rec is None:
        raise FileNotFoundError(f"shard {k} of round {rid} domains-{library}-{library_version}: no record")
    path = round_run_path(tmp, root, rid, library, library_version, k)
    if not path.is_file():
        return "is missing"
    if "run_sha256" not in rec:
        return "cannot be checked: its record holds no run sha256"
    if ids.sha256_file(path).hex() != rec["run_sha256"]:
        return "differs from the sha256 its record holds"
    return None


def rescan_round_shard(root: Path, rid: str, library: str, library_version: str, hmm_path: Path, k: int, cpus: int,
                       *, tmp: Path | None = None) -> dict:
    """Rescan a round shard whose scratch run was lost; must reproduce its published record's sha256."""
    root = Path(root)
    ledger.check_scratch(tmp, root)
    if not isinstance(cpus, int) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give the thread count explicitly (at least 1)")
    lib = f"{library}-{library_version}"
    ledger.refuse_written(
        root, root,
        [round_domains_dir(rid, library, library_version),
         round_coverage_dir(rid, library, library_version),
         round_ledger_dir(rid, library, library_version)],
        f"round {rid} domains-{lib}",
    )
    rec = round_shard_record(root, rid, library, library_version, k)
    if rec is None:
        raise FileNotFoundError(f"shard {k} of round {rid} domains-{lib}: no published record")
    raw, p, s, proteins = _round_shard_inputs(root, rid, library, library_version, hmm_path, k, tmp=tmp)
    if rec["plan_sha256"] != hashlib.sha256(raw).hexdigest() or rec["shard"] != k:
        raise ValueError(f"shard {k} of round {rid} domains-{lib}: its record is of another plan or shard")
    rows = scan(proteins, hmm_path, p["cutoff"], cpus=cpus)
    path = _round_scratch_dir(tmp, raw) / f"shard-{k:05d}.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    new = path.with_name(f".rescan-{os.getpid()}-{path.name}")
    spec = table("round_domains")
    bio.write_table(_batches(rows, spec), spec, new)
    have, hits = ids.sha256_file(new).hex(), len({r["protein_uid"] for r in rows})
    if have != rec["run_sha256"]:
        new.unlink()
        raise ValueError(f"shard {k} of round {rid} domains-{lib}: rescan sha256 {have} differs from record "
                         f"{rec['run_sha256']}")
    os.replace(new, path)
    return {"round_id": rid, "shard": k, "domains": len(rows), "proteins_with_hit": hits,
            "run_sha256": have, "verified": "sha256"}


def merge_round(root: Path, rid: str, library: str, library_version: str, *, tmp: Path | None = None,
                max_file_bytes: int = layout.MAX_FILE_BYTES, out: Path | None = None) -> dict:
    """Merge a round's scanned shards and copied source rows into `round_domains`, `round_coverage`, and
    `ledger/round=<rid>/step=domains-<lib>/` under `out` (default `root`). Bounded memory (AGENTS §4.1)."""
    root = Path(root)
    out_root = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out_root)
    lib = f"{library}-{library_version}"
    step = f"domains-{lib}"
    dom_dir = round_domains_dir(rid, library, library_version)
    cov_dir = round_coverage_dir(rid, library, library_version)
    led_dir = round_ledger_dir(rid, library, library_version)
    ledger.refuse_written(root, out_root, [dom_dir, cov_dir, led_dir], f"round {rid} {step}")

    raw = _round_plan_bytes(root, rid, library, library_version, out)
    p, digest = json.loads(raw), hashlib.sha256(raw).hexdigest()
    if _tools() != p["tools"]:
        raise ValueError(f"round {rid} {step}: this process has tool versions {_tools()}, the plan was made "
                         f"for {p['tools']}")
    req_path = _round_scratch_dir(tmp, raw) / "requested.parquet"
    if not req_path.is_file():
        raise FileNotFoundError(f"{req_path}: requested protein list for round {rid} is missing from scratch")
    if ids.sha256_file(req_path).hex() != p["requested_sha256"]:
        raise ValueError(f"{req_path}: sha256 differs from the round plan's requested_sha256")

    missing = [k for k in range(len(p["shards"])) if _find(root, out, _round_shard_rel(rid, library, library_version, k)) is None]
    if missing:
        raise FileNotFoundError(f"round {rid} {step}: no record of shard {', '.join(map(str, missing))} "
                                f"of {len(p['shards'])}")
    recs, runs = [], []
    for k in range(len(p["shards"])):
        rec = json.loads(_find(root, out, _round_shard_rel(rid, library, library_version, k)).read_bytes())
        if rec["plan_sha256"] != digest or rec["shard"] != k:
            raise ValueError(f"round shard {k}: its record is of another plan or shard")
        path = _round_scratch_dir(tmp, raw) / f"shard-{k:05d}.parquet"
        if not path.is_file():
            raise FileNotFoundError(f"{path}: the run of round shard {k} is not in scratch")
        if ids.sha256_file(path).hex() != rec["run_sha256"]:
            raise ValueError(f"{path}: sha256 differs from the record of round shard {k}")
        recs.append(rec)
        runs.append(path)

    copy = p["copy"]
    srcs = copy["sources"]
    for s in srcs:
        _checked(root, s["files"], f"the round plan (domain files of {_src_id(s)})")
        if "coverage_files" in s:
            _checked(root, s["coverage_files"], f"the round plan (round_coverage of {_src_id(s)})")
        else:
            _checked(root, s["protein_files"], f"the round plan (proteins of {_src_id(s)})")
        if "owned" in s:
            _check_owned(root, s["owned"], f"the round plan (owned proteins of {_src_id(s)})")

    counts = [{"proteins": 0, "residues": 0, "rows": 0, "proteins_with_hit": 0} for _ in srcs]
    models: set[str] = set()
    scanned = _Rows(runs, None)
    want = {"cutoff": p["cutoff"], "z": Z_PROT, "dom_z": DOM_Z}
    dom_spec = table("round_domains")

    def flush(pend: list[pa.Table]) -> pa.Table:
        t = pa.concat_tables(pend).combine_chunks()
        for col, v in want.items():
            if not pc.all(pc.equal(t.column(col), v)).as_py():
                raise ValueError(f"round {rid} {step}: a row to write has another {col} than {v!r}")
        models.update(pc.unique(t.column("model_acc")).to_pylist())
        return t

    def dom_batches() -> Iterator[pa.Table]:
        pend: list[pa.Table] = []
        held = 0
        for u, n, k, rows in _walk_requested(root, req_path, srcs, None):
            if k < 0:
                continue
            c = counts[k]
            c["proteins"] += 1
            c["residues"] += n
            if not rows:
                continue
            c["rows"] += sum(t.num_rows for t in rows)
            c["proteins_with_hit"] += 1
            before = scanned.until(u)
            if scanned.head() == u:
                raise ValueError(f"protein {u.hex()} is copied from {_src_id(srcs[k])} and was also scanned")
            pend += before + rows
            held += sum(t.num_rows for t in before) + sum(t.num_rows for t in rows)
            if held >= bio.ROW_GROUP_ROWS or len(pend) >= FLUSH_PIECES:
                yield flush(pend)
                pend, held = [], 0
        pend += scanned.rest()
        bad = [(_src_id(s), {k: s[k] for k in c}, c) for s, c in zip(srcs, counts) if {k: s[k] for k in c} != c]
        if bad:
            raise ValueError(f"round {rid} {step}: copied counts differ from the plan's: {bad}")
        if pend:
            yield flush(pend)

    bio.write_parts(dom_batches(), dom_spec, out_root / dom_dir, EST_ROW_BYTES, max_file_bytes)

    # Stream round_coverage from req_path merged with the just-written round_domains parts (O(1) memory)
    written_doms = bio.part_files(out_root / dom_dir)
    cov_spec = table("round_coverage")

    def cov_batches() -> Iterator[pa.RecordBatch]:
        dom_cur = _Rows(written_doms, ["protein_uid"])
        u_buf: list[bytes] = []
        n_buf: list[int] = []
        for u, _ in _req_uid_lengths(req_path):
            rows_u = dom_cur.take(u)
            u_buf.append(u)
            n_buf.append(sum(t.num_rows for t in rows_u))
            if len(u_buf) >= bio.ROW_GROUP_ROWS:
                yield pa.RecordBatch.from_arrays(
                    [pa.array(u_buf, pa.binary(16)), pa.array(n_buf, pa.uint32())],
                    schema=cov_spec.schema,
                )
                u_buf, n_buf = [], []
        if u_buf:
            yield pa.RecordBatch.from_arrays(
                [pa.array(u_buf, pa.binary(16)), pa.array(n_buf, pa.uint32())],
                schema=cov_spec.schema,
            )

    bio.write_parts(cov_batches(), cov_spec, out_root / cov_dir, EST_COV_ROW_BYTES, max_file_bytes)

    # Write round ledger: ledger/round=<rid>/step=domains-<lib>/part-0.parquet (§3.3)
    n_in = p["proteins"]
    copied = copy["proteins"]
    scanned_n = n_in - copied
    scanned_hit = sum(r["proteins_with_hit"] for r in recs)
    total_hit = scanned_hit + copy["proteins_with_hit"]
    lib_bytes = lib.encode()
    tools = _tools()
    led_rows = [
        ledger.row(lib_bytes, "library", n_in, scanned_hit, tools=tools),
        ledger.row(lib_bytes, "library", n_in, scanned_hit, n_rejected=scanned_n - scanned_hit,
                   reason="no_hit", tools=tools),
    ]
    if copied:
        led_rows.append(ledger.row(lib_bytes, "library", n_in, scanned_hit, n_rejected=copied,
                                   reason="settled_elsewhere", tools=tools))
    led_rows.sort(key=lambda r: (r["object_key"], r["reason_code"]))
    led_spec = table("ledger")
    bio.write_table(pa.RecordBatch.from_pylist(led_rows, schema=led_spec.schema),
                    led_spec, out_root / led_dir / "part-0.parquet")

    return {
        "round_id": rid,
        "proteins": n_in,
        "domains": sum(r["domains"] for r in recs) + copy["rows"],
        "proteins_with_hit": total_hit,
        "scanned_proteins": scanned_n,
        "scanned_proteins_with_hit": scanned_hit,
        "copied_proteins": copied,
        "copied_rows": copy["rows"],
        "models_with_hit": len(models),
        "shards": len(recs),
        "sources": [_src_id(s) for s in srcs],
    }


def main(argv: list[str] | None = None) -> int:
    """`python -m biobox.l1.domains rescan --root R (--stage S | --round RID) --library L --library-version V
    --hmm F --shard K --cpus N --tmp D`: `rescan_shard` (or `rescan_round_shard`) in its own process, as
    build-stage / run_round runs it (cards C-12, C-38f). It writes nothing into the root, so it is not a
    run-step. Prints one JSON line; exit 2 with the reason on a refusal."""
    ap = argparse.ArgumentParser(prog="python -m biobox.l1.domains")
    ap.add_argument("action", choices=["rescan"])
    for name in ("--root", "--library", "--library-version", "--hmm", "--tmp"):
        ap.add_argument(name, required=True)
    who = ap.add_mutually_exclusive_group(required=True)
    who.add_argument("--stage")
    who.add_argument("--round", dest="round_id")
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--cpus", type=int, required=True)
    a = ap.parse_args(argv)
    try:
        if a.round_id is not None:
            res = rescan_round_shard(Path(a.root), a.round_id, a.library, a.library_version, Path(a.hmm),
                                     a.shard, a.cpus, tmp=Path(a.tmp))
            tag = {"round": a.round_id}
        else:
            res = rescan_shard(Path(a.root), a.stage, a.library, a.library_version, Path(a.hmm),
                               a.shard, a.cpus, tmp=Path(a.tmp))
            tag = {"stage": a.stage}
    except (ValueError, FileNotFoundError, FileExistsError) as e:
        print(f"domains rescan: {e}", file=sys.stderr)
        return 2
    print(json.dumps({"action": a.action, "result": res, **tag}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

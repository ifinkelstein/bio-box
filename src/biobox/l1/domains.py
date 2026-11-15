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
import gzip
import hashlib
import io
import json
import os
import sys
import tarfile
import tempfile
from collections.abc import Iterator
from pathlib import Path

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
    if not files:
        raise FileNotFoundError(f"stage {stage}: no protein_seq files for gene_call_version={gene_call_version}")
    hmm_sha256 = ids.sha256_file(Path(hmm_path)).hex()
    head = {"format": PLAN_FORMAT, "stage": stage, "gene_call_version": gene_call_version, "library": library,
            "library_version": library_version, "cutoff": cutoff, "hmm_sha256": hmm_sha256,
            "shard_residues": shard_residues,
            "files": [{"path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size} for f in files]}
    named = list(dict.fromkeys(s for s in (seq_stages or ()) if s != stage))
    if not named:
        shards = list(_cut(_uid_lengths(files), shard_residues))
        doc = {**head, "proteins": sum(s["proteins"] for s in shards), "residues": sum(s["residues"] for s in shards),
               "shards": shards}
        _write_once(root, out, rel, _json_bytes(doc), "the shard plan")
        return {"shards": len(shards), "proteins": doc["proteins"], "residues": doc["residues"]}
    srcs, passed = sources(root, stage, library, library_version, hmm_sha256, cutoff, named)
    total = {"proteins": 0, "residues": 0}
    counts = [{"proteins": 0, "residues": 0, "rows": 0, "proteins_with_hit": 0} for _ in srcs]

    def scanned() -> Iterator[tuple[bytes, int]]:
        for u, n, k, rows in _walk(root, files, srcs, ["protein_uid"]):
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
    _write_once(root, out, rel, _json_bytes(doc), "the shard plan")
    scanned_p, scanned_r = sum(s["proteins"] for s in shards), sum(s["residues"] for s in shards)
    return {"shards": len(shards), "proteins": doc["proteins"], "residues": doc["residues"],
            "scanned_proteins": scanned_p, "scanned_residues": scanned_r, "copied_proteins": copy["proteins"],
            "copied_residues": copy["residues"], "copied_rows": copy["rows"],
            "sources": [s["stage"] for s in srcs], "passed_over": [p["stage"] for p in passed]}


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
        led_paths = _checked(root, led, f"its run record {rec_rel}")
    except ValueError as e:
        return None, str(e)
    tools = {r for f in led_paths for r in pq.read_table(f, columns=["tool_versions_json"], partitioning=None)
             .column(0).to_pylist()}
    mine = json.dumps(want["tools"], sort_keys=True)
    if tools != {mine}:
        return None, f"its ledger records tool versions {sorted(map(str, tools))}, this scan's are {mine}"
    return {"stage": s, "plan": {"path": plan_rel, "sha256": hashlib.sha256(raw_plan).hexdigest()},
            "record": {"path": rec_rel, "sha256": hashlib.sha256(raw_rec).hexdigest()},
            "protein_files": p["files"],
            "files": [{"path": o["path"], "bytes": o["bytes"], "sha256": o["sha256"]} for o in files]}, None


def sources(root: Path, stage: str, library: str, library_version: str, hmm_sha256: str, cutoff: str,
            seq_stages: list[str] | None) -> tuple[list[dict], list[dict]]:
    """(the sources of copies among `seq_stages`, in their order, each named once and never `stage` itself; the
    stages passed over, each with why), for a scan of `library`-`library_version` with the HMM file of sha256
    `hmm_sha256` and `cutoff` in this process (module docstring). Reads only small records: run records, plans and
    the sources' ledgers."""
    want = {"hmm_sha256": hmm_sha256, "cutoff": cutoff, "z": Z_PROT, "dom_z": DOM_Z, "tools": _tools()}
    srcs, passed = [], []
    for s in dict.fromkeys(x for x in (seq_stages or ()) if x != stage):
        src, why = _source(Path(root), s, library, library_version, want)
        if src is None:
            passed.append({"stage": s, "why": why})
        else:
            srcs.append(src)
    return srcs, passed


def _measured_uids(root: Path, src: dict, lo: bytes | None = None, hi: bytes | None = None) -> Iterator[bytes]:
    """The proteins source `src` measured, in strictly increasing uid order (within [lo, hi] when given, reading only
    the row groups that meet it): every protein of the protein_seq files its plan lists. Card B-47c: once a plan marks
    the proteins its stage inferred, they are measured nowhere, and this is where they must be left out."""
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


def _walk(root: Path, files: list[Path], srcs: list[dict], columns: list[str] | None
          ) -> Iterator[tuple[bytes, int, int, list[pa.Table]]]:
    """(uid, length, k, rows) for every protein of the stage's protein_seq `files`, in uid order: k the index of its
    source in `srcs` (the first that measured it) or -1, and for k >= 0 that source's measured rows of the protein
    (`columns` of its domain files; empty when it had no hit)."""
    cursors = [_Rows([Path(root) / f["path"] for f in s["files"]], columns, measured_only=True) for s in srcs]
    measured = [_measured_uids(root, s) for s in srcs]
    for u, n, k in _assign(_uid_lengths(files), measured):
        yield u, n, k, (cursors[k].take(u) if k >= 0 else [])


def _shard_proteins(root: Path, p: dict, s: dict) -> dict[bytes, str]:
    """The proteins of shard `s` (uids `first` to `last`) from the plan's files, reading only the row
    groups whose protein_uid range meets the shard's."""
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
    srcs = copy["sources"]
    for s in srcs:
        _checked(root, s["files"], f"the shard plan (domain files of {s['stage']})")
        _checked(root, s["protein_files"], f"the shard plan (proteins of {s['stage']})")
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
        for u, n, k, rows in _walk(root, files, srcs, None):
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
                raise ValueError(f"protein {u.hex()} is copied from {srcs[k]['stage']} and was also scanned")
            pend += before + rows
            held += sum(t.num_rows for t in before) + sum(t.num_rows for t in rows)
            if held >= bio.ROW_GROUP_ROWS or len(pend) >= FLUSH_PIECES:
                yield flush(pend)
                pend, held = [], 0
        pend += scanned.rest()
        bad = [(s["stage"], {k: s[k] for k in c}, c) for s, c in zip(srcs, counts) if {k: s[k] for k in c} != c]
        if bad:
            raise ValueError(f"{step} of stage {stage}: copied counts differ from the plan's: {bad}")
        if pend:
            yield flush(pend)

    hit = sum(r["proteins_with_hit"] for r in recs) + copy["proteins_with_hit"]
    _write(out_root, stage, lib_dir, step, batches(), p["proteins"], hit, max_file_bytes)
    return {"proteins": p["proteins"], "domains": sum(r["domains"] for r in recs) + copy["rows"],
            "proteins_with_hit": hit, "models_with_hit": len(models), "shards": len(recs),
            "copied_proteins": copy["proteins"], "copied_rows": copy["rows"],
            "sources": [s["stage"] for s in srcs]}


def main(argv: list[str] | None = None) -> int:
    """`python -m biobox.l1.domains rescan --root R --stage S --library L --library-version V --hmm F --shard K
    --cpus N --tmp D`: `rescan_shard` in its own process, as build-stage runs it (card C-12). It writes nothing
    into the root, so it is not a run-step. Prints one JSON line; exit 2 with the reason on a refusal."""
    ap = argparse.ArgumentParser(prog="python -m biobox.l1.domains")
    ap.add_argument("action", choices=["rescan"])
    for name in ("--root", "--stage", "--library", "--library-version", "--hmm", "--tmp"):
        ap.add_argument(name, required=True)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--cpus", type=int, required=True)
    a = ap.parse_args(argv)
    try:
        res = rescan_shard(Path(a.root), a.stage, a.library, a.library_version, Path(a.hmm), a.shard, a.cpus,
                           tmp=Path(a.tmp))
    except (ValueError, FileNotFoundError, FileExistsError) as e:
        print(f"domains rescan: {e}", file=sys.stderr)
        return 2
    print(json.dumps({"action": a.action, "result": res, "stage": a.stage}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

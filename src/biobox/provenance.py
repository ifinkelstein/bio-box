"""Tool versions and the C-08 provenance check (data_layer_and_tests.md §4.1).

Each step records the versions of the tools that shaped its output in its ledger rows
(`tool_versions_json`). C-08 checks that every versioned table sits under its version
partition, that the versions a `gene_call_version` or `igr_version` names, and those in
the ledgers, equal the container lock: `uv.lock` for Python packages, the conda pins in
`container/environment.yml` for native tools.

An `igr_version` is `<tape_version>-igr1` (card B-20) and is checked through its tape version, as below.
A `tape_version` names no tool versions (card B-07): it is resolved to its producers through the
manifest's `tape_version_def` files, which must hash to it (`biobox.l3.tapes.tape_version`), and
each producer must name its tool versions as a gene-call version does and equal the lock.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tomllib
from importlib import metadata
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
NATIVE = {"mmseqs2": ["mmseqs", "version"]}

# table -> partition key that names the version its rows were made under
VERSION_PARTITION = {
    "protein_seq": "gene_call_version",
    "feature": "producer",
    "feature_protein": "gene_call_version",
    "protein_occurrence": "gene_call_version",
    "igr": "igr_version",
    "l1_domains": "library",
    "cluster_member": "release",
    "cluster_member_aln": "release",
    "contig_tapes": "tape_version",
    "tape_version_def": "tape_version",
    "tape_tokens": "release",
    "token_dict": "release",
    "token_postings": "release",
}
# ledger steps whose output depends on a tool, and the tools they must record
STEP_TOOLS = {"ingest": ("pyarrow",), "tape-*": ("pyarrow",), "genecall": ("pyrodigal", "pyrodigal-gv", "pyarrow"),
              "domains-*": ("pyhmmer", "pyarrow"), "cluster": ("mmseqs2", "pyarrow")}
GENE_CALL_RE = re.compile(r"^pyrodigal-(?P<pyrodigal>[^-]+)-gv-(?P<gv>[^-]+)-p\d+$")


def tool_versions(*names: str) -> dict[str, str]:
    """Installed versions: Python distributions by name, native tools by running them."""
    out = {}
    for n in names:
        if n in NATIVE:
            exe = shutil.which(NATIVE[n][0])
            if exe is None:
                raise FileNotFoundError(f"{NATIVE[n][0]} not on PATH")
            out[n] = subprocess.run([exe, *NATIVE[n][1:]], check=True, capture_output=True, text=True).stdout.strip()
        else:
            out[n] = metadata.version(n)
    return out


def locked_versions(repo: Path = REPO) -> dict[str, str]:
    """The container lock: uv.lock packages plus the conda pins (`name=version`)."""
    lock = tomllib.loads((repo / "uv.lock").read_text())
    out = {p["name"]: p["version"] for p in lock["package"] if "version" in p}
    env = yaml.safe_load((repo / "container" / "environment.yml").read_text())
    for dep in env["dependencies"]:
        if isinstance(dep, str) and "=" in dep:
            name, version = dep.split("=", 1)
            out[name.strip()] = version.split()[0]
    return out


def _step_tools(step: str) -> tuple[str, ...] | None:
    if step in STEP_TOOLS:
        return STEP_TOOLS[step]
    for prefix in ("domains-", "tape-"):  # one ledger per library and per tape version
        if step.startswith(prefix):
            return STEP_TOOLS[prefix + "*"]
    return None


def check_provenance(manifest: dict, lock: dict[str, str] | None = None) -> list[str]:
    """C-08 over a release manifest."""
    import pyarrow.parquet as pq

    lock = locked_versions() if lock is None else lock
    root = Path(manifest["root"])
    problems = []
    for f in manifest["files"]:
        key = VERSION_PARTITION.get(f["table"])
        if key and not f["partitions"].get(key):
            problems.append(f"{f['path']}: {f['table']} without its {key} partition")
        v = f["partitions"].get("gene_call_version")
        if v is not None:
            problems += _gene_call_problems(f["path"], f"gene_call_version={v}", v, lock)
        iv = f["partitions"].get("igr_version")
        if iv is not None and _igr_tape_version(iv) is None:
            problems.append(f"{f['path']}: igr_version={iv} is not <tape_version>-igr1")
        if f["table"] == "ledger":
            step = f["partitions"]["step"]
            need = _step_tools(step)
            if need is None:
                continue
            # one column, a batch at a time: a ledger has one row per contig of the stage (card B-19)
            for raw in _column_values(pq, root / f["path"], "tool_versions_json"):
                tv = json.loads(raw) if raw else None
                if tv is None:
                    problems.append(f"{f['path']}: {step} ledger row without tool versions")
                    break
                missing = [t for t in need if t not in tv]
                wrong = {t: (v, lock.get(t)) for t, v in tv.items() if lock.get(t) != v}
                if missing or wrong:
                    problems.append(f"{f['path']}: {step} tool versions missing {missing}, differ {wrong}")
                    break
    return problems + _tape_version_problems(manifest, root, lock)


def _column_values(pq, path: Path, column: str):
    """The values of one column of a Parquet file, in row order, read a batch at a time."""
    pf = pq.ParquetFile(path)
    if column not in pf.schema_arrow.names:
        yield from (None for _ in range(pf.metadata.num_rows))
        return
    for b in pf.iter_batches(batch_size=65_536, columns=[column]):
        yield from b.column(column).to_pylist()


IGR_SUFFIX = "-igr1"
TAPE_VERSION_RE = re.compile(r"^t-[0-9a-f]{12}$")


def _igr_tape_version(igr_version: str) -> str | None:
    """The tape version an `igr_version` (`<tape_version>-igr1`, card B-20) names; None if malformed."""
    tv = igr_version.removesuffix(IGR_SUFFIX) if igr_version.endswith(IGR_SUFFIX) else None
    return tv if tv is not None and TAPE_VERSION_RE.match(tv) else None


def _gene_call_problems(where: str, what: str, v: str, lock: dict[str, str]) -> list[str]:
    """`v` must name pyrodigal and pyrodigal-gv versions (GENE_CALL_RE) equal to the lock."""
    m = GENE_CALL_RE.match(v)
    if m is None:
        return [f"{where}: {what} does not name its tool versions"]
    if (m["pyrodigal"], m["gv"]) != (lock.get("pyrodigal"), lock.get("pyrodigal-gv")):
        return [f"{where}: {what} differs from the lock "
                f"(pyrodigal {lock.get('pyrodigal')}, pyrodigal-gv {lock.get('pyrodigal-gv')})"]
    return []


def _tape_version_problems(manifest: dict, root: Path, lock: dict[str, str]) -> list[str]:
    """C-08 for tape versions: each is defined by a `tape_version_def` file of the manifest that
    hashes to it, and each of its producers names tool versions equal to the lock."""
    import pyarrow.parquet as pq

    from biobox.l3.tapes import tape_version  # biobox.l3.tapes imports this module

    problems, defined = [], set()
    for f in manifest["files"]:
        tv = f["partitions"].get("tape_version")
        if f["table"] != "tape_version_def" or tv is None:
            continue
        defined.add(tv)
        producers = sorted(pq.read_table(root / f["path"], columns=["producer"], partitioning=None)
                           .column("producer").to_pylist())
        try:
            ok = tape_version(producers) == tv
        except ValueError:
            ok = False
        if not ok:
            problems.append(f"{f['path']}: tape_version_def lists {producers}, which does not hash to tape_version={tv}")
            continue
        for p in producers:
            problems += _gene_call_problems(f["path"], f"tape_version={tv} producer {p}", p, lock)
    for f in manifest["files"]:
        tv = f["partitions"].get("tape_version")
        if tv is not None and f["table"] != "tape_version_def" and tv not in defined:
            problems.append(f"{f['path']}: tape_version={tv} has no tape_version_def in the manifest")
        iv = f["partitions"].get("igr_version")
        tv = None if iv is None else _igr_tape_version(iv)  # a malformed one is reported by the caller
        if tv is not None and tv not in defined:
            problems.append(f"{f['path']}: igr_version={iv} names tape_version={tv}, which has no tape_version_def in the manifest")
    return problems


SEARCH_TOOLS = ("bio-box", "pyhmmer", "pyarrow")


def _versions_problems(where: str, tv: dict | None, lock: dict[str, str]) -> list[str]:
    if not tv:
        return [f"{where}: no tool versions"]
    missing = [t for t in SEARCH_TOOLS if t not in tv]
    wrong = {t: (v, lock.get(t)) for t, v in tv.items() if lock.get(t) != v}
    return [f"{where}: tool versions missing {missing}, differ {wrong}"] if missing or wrong else []


def check_search(root: Path, lock: dict[str, str] | None = None) -> list[str]:
    """C-08 for what lives outside the release manifest: each search pack's builder and each
    closed run's tool versions equal the container lock."""
    import pyarrow.parquet as pq

    lock = locked_versions() if lock is None else lock
    root = Path(root)
    problems = []
    for m in sorted(root.glob("search_pack/release=*/pack_manifest.json")):
        problems += _versions_problems(str(m.relative_to(root)), json.loads(m.read_text()).get("builder"), lock)
    idx = root / "search/run_index.json"
    for r in json.loads(idx.read_text())["runs"] if idx.exists() else []:
        path = root / f"search/run={r['run_id']}/run.parquet"
        row = pq.read_table(path, columns=["tool_versions_json"]).to_pylist()[0]
        tv = json.loads(row["tool_versions_json"]) if row["tool_versions_json"] else None
        problems += _versions_problems(str(path.relative_to(root)), tv, lock)
    return problems


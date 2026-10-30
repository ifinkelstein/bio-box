"""Tool versions and the C-08 provenance check (data_layer_and_tests.md §4.1).

Each step records the versions of the tools that shaped its output in its ledger rows
(`tool_versions_json`). C-08 checks that every versioned table sits under its version
partition, that the versions a `gene_call_version` or `igr_version` names, and those in
the ledgers, equal the container lock: `uv.lock` for Python packages, the conda pins in
`container/environment.yml` for native tools.
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
    "slot_protein": "gene_call_version",
    "protein_occurrence": "gene_call_version",
    "igr": "igr_version",
    "l1_domains": "library",
    "cluster_member": "release",
    "cluster_member_aln": "release",
}
# ledger steps whose output depends on a tool, and the tools they must record
STEP_TOOLS = {"ingest": ("pyarrow",), "genecall": ("pyrodigal", "pyrodigal-gv", "pyarrow"),
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
    return STEP_TOOLS["domains-*"] if step.startswith("domains-") else None


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
        for vkey in ("gene_call_version", "igr_version"):
            v = f["partitions"].get(vkey)
            if v is None:
                continue
            m = GENE_CALL_RE.match(v.removesuffix("-igr1") if vkey == "igr_version" else v)
            if m is None:
                problems.append(f"{f['path']}: {vkey}={v} does not name its tool versions")
            elif (m["pyrodigal"], m["gv"]) != (lock.get("pyrodigal"), lock.get("pyrodigal-gv")):
                problems.append(f"{f['path']}: {vkey}={v} differs from the lock "
                                f"(pyrodigal {lock.get('pyrodigal')}, pyrodigal-gv {lock.get('pyrodigal-gv')})")
        if f["table"] == "ledger":
            step = f["partitions"]["step"]
            need = _step_tools(step)
            if need is None:
                continue
            for r in pq.read_table(root / f["path"]).to_pylist():
                tv = json.loads(r["tool_versions_json"]) if r.get("tool_versions_json") else None
                if tv is None:
                    problems.append(f"{f['path']}: {step} ledger row without tool versions")
                    break
                missing = [t for t in need if t not in tv]
                wrong = {t: (v, lock.get(t)) for t, v in tv.items() if lock.get(t) != v}
                if missing or wrong:
                    problems.append(f"{f['path']}: {step} tool versions missing {missing}, differ {wrong}")
                    break
    return problems

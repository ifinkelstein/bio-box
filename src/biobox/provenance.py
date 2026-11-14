"""Tool versions and the C-08 provenance check (data_layer_and_tests.md §4.1).

Each step records the versions of the tools that shaped its output in its ledger rows
(`tool_versions_json`). C-08 checks that every versioned table sits under its version
partition, that the versions a `gene_call_version` or `igr_version` names, and those in
the ledgers, equal the container lock: `uv.lock` for Python packages, the conda pins in
`container/environment.yml` for native tools.

An `igr_version` is `<tape_version>-igr1` (card B-20) and is checked through its tape version, as below.
A `tape_version` names no tool versions (card B-07): it is resolved to its producers through the
manifest's `tape_version_def` files, which must hash to it (`biobox.l3.tapes.tape_version`), and
each producer must name its tool versions as a gene-call version does and equal the lock. The ncRNA
producer (card E-01) is `infernal-<v>-rfam-<v>`, with an optional `-<set>-<region rule>` tail (card E-07) and an
optional `-<clan rule>` after it (card E-08): Infernal's version must be the conda pin, and the Rfam
library's the version `refdb.yaml` freezes (the lock includes the refdb libraries for this, `rfam_clanin` among
them). A ledger row's tool versions may also name the gene-call version the step read (`gene_call_version`, the
ncRNA producer's intergenic regions, card E-08); that is checked as a gene-call version is.

The CRISPR array producer (card E-12, `biobox.l1.crispr`; checked since card E-14) is
`minced-<version>-w<searchWL>-nr<minNR>-r<minRL>x<maxRL>-s<minSL>x<maxSL>` (`MINCED_RE`): MinCED's version must be
the conda pin, and the rest names MinCED's parameters. Its step `crispr-<producer>` records `minced`, `pyarrow`,
`java` and `minced_args`. `minced` and `pyarrow` are checked against the lock as every tool is; `java` (the JVM, which
the conda solve pulls in with MinCED and nothing pins) must be recorded and non-empty, and is not compared with the
lock; `minced_args` (MinCED's arguments, not a version) must be exactly the arguments the producer's name encodes,
in `crispr.MINCED_ARGS`'s order. The producer name is checked wherever it appears: the step's name, the `producer=`
partition of a file (`feature` files of a `minced-` producer, every `l1_crispr` file), and a tape version's
definition.

The second CRISPR producer (card E-15, DEC-40; `crispr.run_near_cas`) is the same name with a tail,
`-cas<k>kb-<family list tag>-<library>-<version>` (for example `minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-
pfam-38.2`): the library must be one of the lock (refdb.yaml) at its locked version, and the family list tag one of
the repository's frozen lists (`crispr.CAS_FAMILY_SETS`, the name without hyphens). Its step records the default
producer's fields and `gene_call_version` (checked as a gene-call version), `domain_library`, `cas_families`,
`cas_filter` and `overlap_producer`: each must be recorded and equal what the name encodes (the library, the frozen
list whose tag it names, the filter, and the default producer of the same MinCED version and parameters at
`DEFAULT_MIN_REPEATS`).
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
NATIVE = {"mmseqs2": ["mmseqs", "version"], "infernal": ["cmsearch", "-h"]}
# native tools whose version is one line of their output, not all of it: the line's first group (card E-01)
NATIVE_VERSION_RE = {"infernal": re.compile(r"^# INFERNAL (\S+) ", re.M)}

# table -> partition key that names the version its rows were made under
VERSION_PARTITION = {
    "protein_seq": "gene_call_version",
    "feature": "producer",
    "feature_protein": "gene_call_version",
    "protein_occurrence": "gene_call_version",
    "igr": "igr_version",
    "l1_domains": "library",
    "l1_ncrna": "producer",
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
              "domains-*": ("pyhmmer", "pyarrow"), "cluster": ("mmseqs2", "pyarrow"),
              "ncrna-*": ("infernal", "rfam", "pyarrow"),
              "crispr-*": ("minced", "java", "minced_args", "pyarrow")}  # card E-14
GENE_CALL_RE = re.compile(r"^pyrodigal-(?P<pyrodigal>[^-]+)-gv-(?P<gv>[^-]+)-p\d+$")
# the ncRNA producer (biobox.l1.ncrna): Infernal's version and the Rfam library's (card E-01), then optionally the
# model set's version and the region rule (card E-07, `-prok1-igr200nt`), then optionally the clan rule (card E-08,
# `-clan1`). Only the two tool versions are checked against the lock; the set and rules name what was scanned.
NCRNA_RE = re.compile(r"^infernal-(?P<infernal>[^-]+)-rfam-(?P<rfam>[^-]+)"
                      r"(?:-(?P<set>[0-9a-z]+)-(?P<region>igr\d+nt)(?:-(?P<clan>clan\d+))?)?$")
# a ledger's tool versions may name the gene-call version a step read (the ncRNA producer's intergenic regions, card
# E-08): it is checked as a gene-call version (its pyrodigal versions equal the lock), not as a tool of the lock
LEDGER_GENE_CALL_KEY = "gene_call_version"
# the CRISPR array producer (biobox.l1.crispr, card E-12): MinCED's version, then its parameters (card E-14); the
# filtered producer of card E-15 adds a tail: the filter (`cas<k>kb`), the Cas family list's tag and the domain library
# it read (`<name>-<version>`, its `library=` partition)
MINCED_RE = re.compile(r"^minced-(?P<minced>[^-]+)-(?P<params>w(?P<w>\d+)-nr(?P<nr>\d+)-r(?P<rl>\d+)x(?P<rh>\d+)"
                       r"-s(?P<sl>\d+)x(?P<sh>\d+))"
                       r"(?:-(?P<filter>cas\d+kb)-(?P<famset>[0-9a-z]+)-(?P<lib>[a-z0-9_]+)-(?P<libv>[A-Za-z0-9_.]+))?$")
CRISPR_STEP_PREFIX = "crispr-"
# ledger fields of the crispr step that are not versions of the lock: the JVM (unpinned; recorded, not compared)
# and MinCED's arguments (compared with the producer's name)
CRISPR_UNLOCKED = ("java", "minced_args")
# the filtered producer's step (card E-15) records these too: none is a version of the lock, each is compared with what
# the producer's name encodes (its gene-call version, also recorded, is checked as a gene-call version)
NEAR_CAS_FIELDS = ("gene_call_version", "domain_library", "cas_families", "cas_filter", "overlap_producer")
NEAR_CAS_UNLOCKED = ("domain_library", "cas_families", "cas_filter", "overlap_producer")
# the minimum of repeats of the default producer, whose arrays the filtered one may not overlap (DEC-40: 3 everywhere)
DEFAULT_MIN_REPEATS = 3


def minced_args_of(producer: str) -> str | None:
    """The MinCED arguments a `minced-` producer name encodes, as the producer passes and records them
    (`crispr.MINCED_ARGS`: -searchWL, -minNR, -minRL, -maxRL, -minSL, -maxSL); None for another name."""
    m = MINCED_RE.match(producer)
    if m is None:
        return None
    return f"-searchWL {m['w']} -minNR {m['nr']} -minRL {m['rl']} -maxRL {m['rh']} -minSL {m['sl']} -maxSL {m['sh']}"


def tool_versions(*names: str) -> dict[str, str]:
    """Installed versions: Python distributions by name, native tools by running them."""
    out = {}
    for n in names:
        if n in NATIVE:
            exe = shutil.which(NATIVE[n][0])
            if exe is None:
                raise FileNotFoundError(f"{NATIVE[n][0]} not on PATH")
            text = subprocess.run([exe, *NATIVE[n][1:]], check=True, capture_output=True, text=True).stdout.strip()
            if n in NATIVE_VERSION_RE:
                m = NATIVE_VERSION_RE[n].search(text)
                if m is None:
                    raise ValueError(f"{exe}: no {n} version in the output of {' '.join(NATIVE[n])}")
                text = m[1]
            out[n] = text
        else:
            out[n] = metadata.version(n)
    return out


def locked_versions(repo: Path = REPO) -> dict[str, str]:
    """The container lock: uv.lock packages plus the conda pins (`name=version`), and the reference
    libraries frozen in `refdb.yaml` (`rfam`: "15.1"; card E-01 records the library in its ledger). A
    package name never clashes with a library name; a clash keeps the package's version."""
    lock = tomllib.loads((repo / "uv.lock").read_text())
    out = {p["name"]: p["version"] for p in lock["package"] if "version" in p}
    env = yaml.safe_load((repo / "container" / "environment.yml").read_text())
    for dep in env["dependencies"]:
        if isinstance(dep, str) and "=" in dep:
            name, version = dep.split("=", 1)
            out[name.strip()] = version.split()[0]
    for db in refdb_libraries(repo):
        out.setdefault(db["name"], db["version"])
    return out


def refdb_libraries(repo: Path = REPO) -> list[dict]:
    """The reference libraries of `refdb.yaml` (DEC-12): name, version (as a string), file, sha256."""
    doc = yaml.safe_load((repo / "refdb.yaml").read_text())
    return [{**db, "version": str(db["version"])} for db in doc.get("databases", [])]


def _near_cas(producer: str) -> bool:
    """Whether a `minced-` producer name has the filtered producer's tail (card E-15)."""
    m = MINCED_RE.match(producer)
    return m is not None and m["filter"] is not None


def _step_tools(step: str) -> tuple[str, ...] | None:
    if step in STEP_TOOLS:
        return STEP_TOOLS[step]
    if step.startswith(CRISPR_STEP_PREFIX) and _near_cas(step.removeprefix(CRISPR_STEP_PREFIX)):
        return STEP_TOOLS[CRISPR_STEP_PREFIX + "*"] + NEAR_CAS_FIELDS  # card E-15
    # one ledger per library, tape version, ncRNA producer and CRISPR array producer
    for prefix in ("domains-", "tape-", "ncrna-", CRISPR_STEP_PREFIX):
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
        prod = f["partitions"].get("producer")
        if prod is not None and (f["table"] == "l1_crispr" or (f["table"] == "feature" and prod.startswith("minced-"))):
            problems += _minced_problems(f["path"], f"{f['table']} producer={prod}", prod, lock)  # card E-14
        if f["table"] == "ledger":
            step = f["partitions"]["step"]
            need = _step_tools(step)
            if need is None:
                continue
            crispr = step.startswith(CRISPR_STEP_PREFIX)
            if crispr:  # card E-14: the producer the step is named after
                problems += _minced_problems(f["path"], f"step {step}", step.removeprefix(CRISPR_STEP_PREFIX), lock)
            near = crispr and _near_cas(step.removeprefix(CRISPR_STEP_PREFIX))
            unlocked = (LEDGER_GENE_CALL_KEY, *(CRISPR_UNLOCKED if crispr else ()), *(NEAR_CAS_UNLOCKED if near else ()))
            # one column, a batch at a time: a ledger has one row per contig of the stage (card B-19)
            for raw in _column_values(pq, root / f["path"], "tool_versions_json"):
                tv = json.loads(raw) if raw else None
                if tv is None:
                    problems.append(f"{f['path']}: {step} ledger row without tool versions")
                    break
                missing = [t for t in need if t not in tv]
                wrong = {t: (v, lock.get(t)) for t, v in tv.items() if t not in unlocked and lock.get(t) != v}
                gcv = tv.get(LEDGER_GENE_CALL_KEY)
                gcv_problems = [] if gcv is None else _gene_call_problems(
                    f["path"], f"{step} ledger {LEDGER_GENE_CALL_KEY}={gcv}", gcv, lock)
                field_problems = _crispr_field_problems(f["path"], step, tv) if crispr else []
                if missing or wrong or gcv_problems or field_problems:
                    problems += ([f"{f['path']}: {step} tool versions missing {missing}, differ {wrong}"]
                                 if missing or wrong else []) + gcv_problems + field_problems
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


def _minced_problems(where: str, what: str, p: str, lock: dict[str, str]) -> list[str]:
    """`p` must be a CRISPR array producer (`MINCED_RE`, card E-14) whose MinCED version equals the lock."""
    m = MINCED_RE.match(p)
    if m is None:
        return [f"{where}: {what} does not name its tool versions and parameters "
                "(minced-<version>-w<n>-nr<n>-r<min>x<max>-s<min>x<max>, then for the filtered producer of card E-15 "
                "-cas<k>kb-<family list>-<library>-<version>)"]
    if m["minced"] != lock.get("minced"):
        return [f"{where}: {what} differs from the lock (minced {lock.get('minced')})"]
    if m["filter"] is None:
        return []
    # card E-15: the domain library it read is a library of the lock (refdb.yaml) at the locked version, and its Cas
    # family list a frozen list of the repository
    if m["lib"] not in lock:
        return [f"{where}: {what} names the domain library {m['lib']}, which is not a library of the lock (refdb.yaml)"]
    if m["libv"] != lock[m["lib"]]:
        return [f"{where}: {what} differs from the lock ({m['lib']} {lock[m['lib']]})"]
    if m["famset"] not in _cas_family_tags():
        return [f"{where}: {what} names the Cas family list {m['famset']!r}, which is not a frozen family list of the "
                f"repository ({', '.join(sorted(_cas_family_tags().values()))})"]
    return []


def _cas_family_tags() -> dict[str, str]:
    """tag -> name of the frozen Cas family lists (`biobox.l1.crispr.CAS_FAMILY_SETS`; the tag is the name without
    hyphens, as the filtered producer's name carries it). Imported here: biobox.l1.crispr imports this module."""
    from biobox.l1.crispr import CAS_FAMILY_SETS

    return {name.replace("-", ""): name for name in CAS_FAMILY_SETS}


def _crispr_field_problems(where: str, step: str, tv: dict) -> list[str]:
    """The crispr step's fields that are not versions of the lock (card E-14): `java` recorded and non-empty, and
    `minced_args` exactly the arguments the step's producer encodes. A field that is absent is reported as missing
    by the caller."""
    out = []
    if "java" in tv and not (isinstance(tv["java"], str) and tv["java"].strip()):
        out.append(f"{where}: {step} ledger java is empty: the JVM's version must be recorded")
    want = minced_args_of(step.removeprefix(CRISPR_STEP_PREFIX))
    if "minced_args" in tv and want is not None and tv["minced_args"] != want:
        out.append(f"{where}: {step} ledger minced_args {tv['minced_args']!r} are not the arguments its producer "
                   f"names ({want!r})")
    m = MINCED_RE.match(step.removeprefix(CRISPR_STEP_PREFIX))
    if m is None or m["filter"] is None:
        return out
    # card E-15: the filter's fields, each what the name encodes (the gene-call version is checked by the caller)
    expect = {"domain_library": f"{m['lib']}-{m['libv']}", "cas_filter": m["filter"],
              "cas_families": _cas_family_tags().get(m["famset"]),
              "overlap_producer": f"minced-{m['minced']}-w{m['w']}-nr{DEFAULT_MIN_REPEATS}-r{m['rl']}x{m['rh']}"
                                  f"-s{m['sl']}x{m['sh']}"}
    for k, v in expect.items():
        if k in tv and tv[k] != v:
            out.append(f"{where}: {step} ledger {k} {tv[k]!r} is not what its producer names ({v!r})")
    return out


def _producer_problems(where: str, what: str, p: str, lock: dict[str, str]) -> list[str]:
    """A feature producer names its tool versions: the ncRNA producer (`NCRNA_RE`, card E-01) Infernal's and
    the Rfam library's, which must equal the lock; the CRISPR array producer (`MINCED_RE`, card E-14) MinCED's, which
    must equal the lock, and its parameters; any other producer is a gene-call version."""
    if p.startswith("minced-"):
        return _minced_problems(where, what, p, lock)
    m = NCRNA_RE.match(p)
    if m is None:
        return _gene_call_problems(where, what, p, lock)
    if (m["infernal"], m["rfam"]) != (lock.get("infernal"), lock.get("rfam")):
        return [f"{where}: {what} differs from the lock (infernal {lock.get('infernal')}, rfam {lock.get('rfam')})"]
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
            problems += _producer_problems(f["path"], f"tape_version={tv} producer {p}", p, lock)
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


def check_search(packs: Path | None = None, runs: Path | None = None,
                 lock: dict[str, str] | None = None) -> list[str]:
    """C-08 for what lives outside the release manifest, and outside the data root (card W-03): each search
    pack's builder in the pack directory `packs` (`release=<R>/<pack_id>/pack_manifest.json`) and each
    closed run's tool versions in the runs directory `runs` (`run=<run_id>/run.parquet`) equal the
    container lock."""
    import pyarrow.parquet as pq

    lock = locked_versions() if lock is None else lock
    problems = []
    if packs is not None:
        packs = Path(packs)
        for m in sorted(packs.glob("release=*/*/pack_manifest.json")):
            if m.parent.name.startswith("."):  # a pack being built, under its temporary name
                continue
            problems += _versions_problems(str(m.relative_to(packs)), json.loads(m.read_text()).get("builder"), lock)
    if runs is not None:
        runs = Path(runs)
        for path in sorted(runs.glob("run=*/run.parquet")):
            row = pq.read_table(path, columns=["tool_versions_json"]).to_pylist()[0]
            tv = json.loads(row["tool_versions_json"]) if row["tool_versions_json"] else None
            problems += _versions_problems(str(path.relative_to(runs)), tv, lock)
    return problems


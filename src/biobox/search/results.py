"""A run's results kept on durable storage: `bio-box results archive|flush|list|get` (card W-28, DEC-33).

Design (card W-28; DEC-33; AGENTS §3.9, §4.3, §7.3, §7.5):

- What a run's archive holds. `run/`: every regular file of the run's directory (`run.parquet`, the record, and
  the hits, profiles and neighbourhood tables beside it). `outputs/<name>/`: each output a user would want back,
  given with `--with` and named by its own base name: the `search hmm --out DIR` (`hits.tsv`, `hits.faa`,
  `neighbors.tsv`, `loci/`, `README.txt`, `recipe.yaml`, and anything written inside it, such as an
  `nbr flanks --out DIR/flanks`), an `nbr flanks` directory elsewhere, a `hits.tsv` written alone. They are named
  by hand because nothing on disk links a run to its outputs: the record's `argv` is the worker's command, not
  the user's (`runner.run` reuses the name), and the exports carry no run id. One check is possible and made: a
  directory holding a `recipe.yaml` that is not the run's recipe is another run's `search hmm` output, refused.
  Never the pack, the work directory or the pending record. The first member is `FILES.tsv`: each other member's
  path, size and sha256. A symbolic link or any other non-regular file is refused: the archive holds bytes, not
  where they came from.
- Deterministic (AGENTS §4.3). An uncompressed POSIX (pax) tar: compression is left out because zlib's bytes
  differ between library versions, and the same run must give the same bytes on any host. Members in sorted
  order of their names; mtime 0, uid and gid 0, no user or group name, mode 0644; regular files only, no
  directory entries. Names are relative to the run; no host name, absolute path or time of archiving is in the
  archive. (The run's own record holds its argv and host, as the runner wrote them; they are the run's, and are
  kept as they are.) The date recorded is the run's own, its record's `submitted_at` in UTC.
- The 64 MB rule (AGENTS §7.3). Nearly every run's archive is far under 64 MB, so a run is never written to
  `/stor` alone. `archive` writes the run's archive and its manifest line into a spool on the profile's local
  `tmp` (`<tmp>/results/spool/<run id>.tar` and `.tsv`), where small files are allowed. `flush` packs every
  spooled archive into one bundle (a tar whose first member, `MANIFEST.tsv`, holds one manifest line per run,
  then each run's archive, unchanged, as `<run id>.tar`, in run-id order), built on local disk first, and moves
  it to the results area only when the bundle is at least the profile's `small_file_floor_bytes`; below it,
  nothing is written there and the runs stay spooled. `--below-floor` writes a smaller bundle and says it breaks
  §7.3: only with the PI's agreement. So each run is still one archive (DEC-33), and what lands on `/stor` is
  one file of at least 64 MB or nothing. The spool is local disk: a reboot may clear it, as it clears the run.
- One archive per run, never overwritten, never deleted (AGENTS §7.5). Before a run is spooled or a bundle
  moved, the run ids are looked up in the spool's and the bundles' manifests: the same run with the same sha256
  is a no-op (said on stderr, exit 0); with another sha256 it is refused (exit 2), naming both. A bundle is copied
  to a dot-named temporary file in the results area, read back and checked, then renamed into place; an existing
  bundle of the same name is never replaced. A flush killed before the rename leaves only its dot-named
  temporary file, which `list` names and no reader takes; the tools never remove it (§7.5: ask the PI). The spool's
  entries of a bundle's runs are removed only after the bundle is in place and read back.
- Where. The results area is `<durable>/results`, the profile's durable storage (`durable`, which RUNBOOK §10
  calls the durable results), or a directory given with `--results` that must lie inside `durable`. That is
  decided on the resolved path (symbolic links followed) and on the file system: the area's nearest existing
  directory must be on the same device as `durable`, so a mount or a link inside `durable` that leads elsewhere
  is refused. The profile needs no new key. (DEC-33 says "under the data root's parent"; a profile names no data
  root, and a search runs against roots on `/tmp` as well as on `/stor`, so the durable storage the profile names
  stands in for it.)
- Reading. `list` reads only manifest lines: the spool's `.tsv` files and the first member of each bundle; it
  never reads a run's archive. The results area has no index file of its own (it would be a small file on
  `/stor`): its bundles are the files whose names match `BUNDLE_RE`. `get RUN --out DIR` copies the run's archive
  out of its bundle (or the spool) while hashing it, checks it against the manifest line, then restores every
  member into a new or empty DIR, checking each file's size and sha256 against `FILES.tsv`; DIR is filled under
  a temporary name beside it and renamed into place.
- Bounded memory (AGENTS §4.1): every file is streamed in 1 MiB blocks; nothing holds a whole archive.
"""

from __future__ import annotations

import fcntl
import hashlib
import io
import os
import re
import secrets
import shutil
import tarfile
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import timezone
from pathlib import Path

import pyarrow.parquet as pq

FILES = "FILES.tsv"            # first member of a run's archive
MANIFEST = "MANIFEST.tsv"      # first member of a bundle
ARCHIVE_HEADER = "# bio-box results archive 1"
MANIFEST_HEADER = "# bio-box results manifest 1"
FILE_COLUMNS = ("path", "size", "sha256")
COLUMNS = ("run_id", "release", "run_date", "status", "sha256", "size", "n_files", "outputs")
BUNDLE_RE = re.compile(r"^results-(\d{8})-(\d{8})-([0-9a-f]{16})\.tar$")
RUN_ID_RE = re.compile(r"^[0-9a-f]{32}$")
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]*$")
BLOCK = 1 << 20


class ResultsError(ValueError):
    """A refusal (exit 2): nothing was written."""


class Corrupt(ValueError):
    """An archive or a file whose bytes are not those its manifest line or FILES.tsv names (exit 1)."""


@dataclass(frozen=True)
class Places:
    spool: Path    # local: <profile tmp>/results/spool
    area: Path     # durable: <profile durable>/results, or --results inside durable
    durable: Path
    floor: int     # profile small_file_floor_bytes


# --- where ------------------------------------------------------------------------------------------------------


def check_area(area: Path, durable: Path) -> Path:
    """The results area, resolved, if it is on the profile's durable storage; else ResultsError. Decided on the
    resolved path and on the device of its nearest existing directory, never on the path's spelling."""
    d = Path(durable)
    if not d.is_dir():
        raise ResultsError(f"the profile's durable storage {d} does not exist (it is made by the PI, not by bio-box)")
    dr = d.resolve(strict=True)
    ar = Path(area).resolve()
    if ar != dr and dr not in ar.parents:
        raise ResultsError(f"the results area {area} (resolved: {ar}) is not inside the profile's durable storage "
                           f"{d} (resolved: {dr})")
    near = ar
    while not near.exists():
        near = near.parent
    if not near.is_dir():
        raise ResultsError(f"the results area {area}: {near} is not a directory")
    if os.stat(near).st_dev != os.stat(dr).st_dev:
        raise ResultsError(f"the results area {area} (resolved: {ar}) is on another file system than the profile's "
                           f"durable storage {dr}")
    return ar


def places(prof, results: Path | None = None) -> Places:
    """The spool, the checked results area and the floor of a host profile (`biobox.profile.Profile`)."""
    durable = Path(prof.durable)
    area = check_area(Path(results) if results is not None else durable / "results", durable)
    return Places(spool=Path(prof.tmp) / "results" / "spool", area=area, durable=durable,
                  floor=int(prof.small_file_floor_bytes))


@contextmanager
def _locked(spool: Path):
    """One `archive` or `flush` at a time on this host (the spool's lock file is on local disk, never on /stor)."""
    spool.mkdir(parents=True, exist_ok=True)
    with open(spool / ".lock", "a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


# --- manifest lines ---------------------------------------------------------------------------------------------


def manifest_text(lines: list[dict]) -> str:
    rows = [MANIFEST_HEADER, "\t".join(COLUMNS)]
    rows += ["\t".join(str(line[c]) for c in COLUMNS) for line in sorted(lines, key=lambda r: r["run_id"])]
    return "\n".join(rows) + "\n"


def parse_manifest(text: str, where: str) -> list[dict]:
    rows = text.splitlines()
    if len(rows) < 2 or rows[0] != MANIFEST_HEADER or rows[1].split("\t") != list(COLUMNS):
        raise Corrupt(f"{where}: not a bio-box results manifest")
    out = []
    for r in rows[2:]:
        v = r.split("\t")
        if len(v) != len(COLUMNS):
            raise Corrupt(f"{where}: a manifest line has {len(v)} fields, not {len(COLUMNS)}")
        line = dict(zip(COLUMNS, v))
        line["size"], line["n_files"] = int(line["size"]), int(line["n_files"])
        out.append(line)
    return out


# --- tar members ------------------------------------------------------------------------------------------------


def _info(name: str, size: int) -> tarfile.TarInfo:
    ti = tarfile.TarInfo(name)
    ti.size, ti.mtime, ti.mode, ti.type = size, 0, 0o644, tarfile.REGTYPE
    ti.uid = ti.gid = 0
    ti.uname = ti.gname = ""
    return ti


def _writer(path: Path) -> tarfile.TarFile:
    return tarfile.open(path, "w", format=tarfile.PAX_FORMAT, encoding="utf-8")


class _Hashing:
    """A reader that hashes what passes through it."""

    def __init__(self, fh):
        self.fh, self.h, self.n = fh, hashlib.sha256(), 0

    def read(self, n: int = -1) -> bytes:
        b = self.fh.read(n)
        self.h.update(b)
        self.n += len(b)
        return b


def sha256_file(path: Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as fh:
        while b := fh.read(BLOCK):
            h.update(b)
            n += len(b)
    return h.hexdigest(), n


def _add_bytes(tf: tarfile.TarFile, name: str, data: bytes) -> None:
    tf.addfile(_info(name, len(data)), io.BytesIO(data))


def _add_file(tf: tarfile.TarFile, name: str, path: Path, size: int, sha: str) -> None:
    with open(path, "rb") as fh:
        src = _Hashing(fh)
        tf.addfile(_info(name, size), src)
    if src.n != size or src.h.hexdigest() != sha:
        raise ResultsError(f"{path} changed while it was archived; archive the run again when nothing writes it")


def _fsync(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


# --- a run's archive --------------------------------------------------------------------------------------------


def read_record(run_dir: Path) -> dict:
    """The record of a closed run (ResultsError when the run is not closed)."""
    d = Path(run_dir)
    if not (d / "run.parquet").is_file():
        raise ResultsError(f"{d} holds no run.parquet: the run is not closed (or this is not a run directory); "
                           "archive it once `bio-box search status` says it has closed")
    return pq.read_table(d / "run.parquet").to_pylist()[0]


def run_date(record: dict) -> str:
    """The run's own date: its record's `submitted_at`, in UTC, to the second."""
    return record["submitted_at"].astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def check_output(record: dict, path: Path) -> None:
    """Refuse an output directory that is another run's: a `search hmm --out DIR` holds the recipe its run ran
    (`recipe.yaml`, the record's `recipe_yaml` byte for byte, card W-07), so a DIR whose recipe is not this run's
    belongs to another run. Other outputs carry nothing that names their run, and are taken as given."""
    recipe = Path(path) / "recipe.yaml"
    if Path(path).is_dir() and recipe.is_file() and recipe.read_text(encoding="utf-8") != record["recipe_yaml"]:
        raise ResultsError(f"{path} holds a recipe.yaml that is not the recipe run {record['run_id'].hex()} ran: it is "
                           "another run's output (`search hmm` writes the run's recipe there)")


def _members(base: str, path: Path) -> list[tuple[str, Path]]:
    """(name in the archive, file) of a file, or of every file under a directory; anything but regular files and
    directories is refused."""
    path = Path(path)
    if path.is_symlink():
        raise ResultsError(f"{path} is a symbolic link: give the file or directory it points to")
    if path.is_file():
        return [(base, path)]
    if not path.is_dir():
        raise ResultsError(f"{path} is not a regular file or directory")
    out = []
    for dirpath, dirnames, filenames in os.walk(path):
        for name in dirnames + filenames:
            p = Path(dirpath) / name
            if p.is_symlink() or not (p.is_dir() or p.is_file()):
                raise ResultsError(f"{p} is a symbolic link or not a regular file: the archive holds regular files "
                                   "only")
        for name in filenames:
            p = Path(dirpath) / name
            out.append((f"{base}/{p.relative_to(path).as_posix()}", p))
    return out


def plan(run_dir: Path, outputs: list[Path]) -> tuple[dict, list[tuple[str, Path]], list[str]]:
    """The run's record, its archive's members (sorted by name) and the outputs' names."""
    record = read_record(run_dir)
    members = _members("run", Path(run_dir))
    names: dict[str, Path] = {}
    for o in outputs:
        o = Path(o)
        name = o.name
        if not NAME_RE.match(name):
            raise ResultsError(f"output {o}: its name {name!r} must be letters, digits and . _ + - (not first)")
        if name in names and names[name].resolve() != o.resolve():
            raise ResultsError(f"two outputs are named {name!r} ({names[name]} and {o}): an output is archived "
                               "under its own name, so give one of them another name")
        rd, od = Path(run_dir).resolve(), o.resolve()
        if od == rd or rd in od.parents or od in rd.parents:
            raise ResultsError(f"output {o} is the run's directory, inside it or holds it: the run is archived as "
                               "run/ already")
        check_output(record, o)
        for other in names.values():
            a, b = other.resolve(), od
            if a in b.parents or b in a.parents:
                inner, outer = (o, other) if a in b.parents else (other, o)
                raise ResultsError(f"output {inner} is inside the output {outer}: it would be kept twice")
        if name not in names:
            names[name] = o
            members += _members(f"outputs/{name}", o)
    members.sort(key=lambda m: m[0])
    return record, members, sorted(names)


def write_archive(run_dir: Path, outputs: list[Path], dest: Path) -> dict:
    """Write the run's archive to `dest` (a new file) and return its manifest line."""
    record, members, names = plan(run_dir, outputs)
    listed = []
    for name, p in members:
        sha, size = sha256_file(p)
        listed.append((name, size, sha))
    files = "\n".join([ARCHIVE_HEADER, "\t".join(FILE_COLUMNS)] + [f"{n}\t{s}\t{h}" for n, s, h in listed]) + "\n"
    with _writer(dest) as tf:
        _add_bytes(tf, FILES, files.encode("utf-8"))
        for (name, p), (_, size, sha) in zip(members, listed):
            _add_file(tf, name, p, size, sha)
    sha, size = sha256_file(dest)
    return {"run_id": record["run_id"].hex(), "release": record["ran_against"], "run_date": run_date(record),
            "status": record["status"], "sha256": sha, "size": size, "n_files": len(listed),
            "outputs": ",".join(names) or "-"}


# --- what is kept -----------------------------------------------------------------------------------------------


def _bundle_manifest(path: Path) -> list[dict]:
    """The manifest lines of a bundle: its first member only; the runs' archives are not read."""
    try:
        with tarfile.open(path, "r:") as tf:
            first = tf.next()
            if first is None or first.name != MANIFEST or not first.isfile():
                raise Corrupt(f"{path}: its first member is not {MANIFEST}")
            return parse_manifest(tf.extractfile(first).read().decode("utf-8"), str(path))
    except tarfile.TarError as e:
        raise Corrupt(f"{path} is named as a bundle but is not one ({e})") from None


def bundles(area: Path) -> list[Path]:
    if not area.is_dir():
        return []
    return sorted(area / n for n in os.listdir(area) if BUNDLE_RE.match(n))


def partials(area: Path) -> list[Path]:
    """Temporary files of flushes that did not finish: never read, never removed by the tools."""
    if not area.is_dir():
        return []
    return sorted(area / n for n in os.listdir(area) if n.startswith(".results-") and ".part-" in n)


def spooled(spool: Path) -> list[dict]:
    """The spool's complete entries (an archive and its manifest line), each with `where` (the archive's path)."""
    if not spool.is_dir():
        return []
    out = []
    for n in sorted(os.listdir(spool)):
        if n.endswith(".tsv") and RUN_ID_RE.match(n[:-4]) and (spool / f"{n[:-4]}.tar").is_file():
            for line in parse_manifest((spool / n).read_text(encoding="utf-8"), str(spool / n)):
                out.append({**line, "where": str(spool / f"{n[:-4]}.tar"), "kept": "spool"})
    return out


def flushed(area: Path) -> list[dict]:
    out = []
    for b in bundles(area):
        out += [{**line, "where": str(b), "kept": "bundle"} for line in _bundle_manifest(b)]
    return out


def listing(where: Places) -> list[dict]:
    """Every kept run, flushed and spooled, from manifest lines only; sorted by run id, bundles first."""
    return sorted(flushed(where.area) + spooled(where.spool), key=lambda r: (r["run_id"], r["kept"] != "bundle"))


def _known(where: Places, run_id: str) -> list[dict]:
    return [r for r in listing(where) if r["run_id"] == run_id]


# --- archive ----------------------------------------------------------------------------------------------------


def archive(run_dir: Path, outputs: list[Path], where: Places, *, publish=os.replace) -> tuple[dict, str]:
    """Spool a run's archive. Returns (manifest line, state): `spooled`, or `same` when that run is already kept with
    the same sha256 (nothing written). ResultsError when it is kept with another sha256. `publish` is the rename
    of the archive into the spool (a test replaces it to simulate a crash before it)."""
    with _locked(where.spool):
        tmp = where.spool / f".archive.tmp-{os.getpid()}"
        try:
            line = write_archive(run_dir, outputs, tmp)
            rid = line["run_id"]
            for k in _known(where, rid):
                if k["sha256"] == line["sha256"]:
                    return {**line, "where": k["where"]}, "same"
                raise ResultsError(f"run {rid} is already archived in {k['where']} with sha256 {k['sha256']}; this "
                                   f"archive would be sha256 {line['sha256']} (the run's files or outputs differ from "
                                   "those archived). An archive is never overwritten")
            final = where.spool / f"{rid}.tar"
            publish(tmp, final)  # local spool: an orphan archive without its manifest line is replaced
            text = where.spool / f".{rid}.tsv.tmp-{os.getpid()}"
            text.write_text(manifest_text([line]), encoding="utf-8")
            os.replace(text, where.spool / f"{rid}.tsv")
            return {**line, "where": str(final)}, "spooled"
        finally:
            tmp.unlink(missing_ok=True)


# --- flush ------------------------------------------------------------------------------------------------------


def bundle_name(lines: list[dict], manifest: str) -> str:
    days = sorted(r["run_date"][:10].replace("-", "") for r in lines)
    return f"results-{days[0]}-{days[-1]}-{hashlib.sha256(manifest.encode('utf-8')).hexdigest()[:16]}.tar"


def _write_bundle(lines: list[dict], spool: Path, dest: Path) -> str:
    manifest = manifest_text(lines)
    with _writer(dest) as tf:
        _add_bytes(tf, MANIFEST, manifest.encode("utf-8"))
        for line in sorted(lines, key=lambda r: r["run_id"]):
            _add_file(tf, f"{line['run_id']}.tar", spool / f"{line['run_id']}.tar", line["size"], line["sha256"])
    return manifest


def _drop(spool: Path, run_id: str) -> None:
    """Remove a run's spool entry (local disk), its manifest line first: an archive without one is never read."""
    (spool / f"{run_id}.tsv").unlink(missing_ok=True)
    (spool / f"{run_id}.tar").unlink(missing_ok=True)


def flush(where: Places, *, below_floor: bool = False, publish=os.rename) -> dict:
    """Pack the spool into one bundle and move it to the results area when it reaches the floor (or `below_floor`).
    Returns what happened: `state` is `empty`, `below_floor` (nothing written), `written` or `same`. `publish` is the
    final rename (a test replaces it to simulate a crash before it)."""
    with _locked(where.spool):
        kept = {r["run_id"]: r for r in flushed(where.area)}
        todo, done = [], []
        for line in spooled(where.spool):
            k = kept.get(line["run_id"])
            if k is None:
                todo.append(line)
            elif k["sha256"] == line["sha256"]:
                done.append(line)
            else:
                raise ResultsError(f"run {line['run_id']} is spooled with sha256 {line['sha256']} and already in "
                                   f"{k['where']} with sha256 {k['sha256']}; nothing written")
        for line in done:  # a flush that stopped after its rename: these runs are in a bundle already
            _drop(where.spool, line["run_id"])
        if not todo:
            return {"state": "empty", "runs": 0, "size": 0, "dropped": len(done)}
        local = where.spool / f".bundle.tmp-{os.getpid()}"
        try:
            for line in todo:
                sha, size = sha256_file(where.spool / f"{line['run_id']}.tar")
                if (sha, size) != (line["sha256"], line["size"]):
                    raise Corrupt(f"the spooled archive of run {line['run_id']} is not the one its manifest line "
                                  f"names (sha256 {sha}, not {line['sha256']}); archive the run again")
            manifest = _write_bundle(todo, where.spool, local)
            sha, size = sha256_file(local)
            name = bundle_name(todo, manifest)
            out = {"runs": len(todo), "size": size, "sha256": sha, "dropped": len(done),
                   "bundle": str(where.area / name)}
            if size < where.floor and not below_floor:
                return {**out, "state": "below_floor"}
            final = where.area / name
            if final.exists():
                have, _ = sha256_file(final)
                if have != sha:
                    raise ResultsError(f"{final} exists with sha256 {have}; this bundle would be sha256 {sha}. A "
                                       "bundle is never overwritten")
                state = "same"
            else:
                where.area.mkdir(parents=True, exist_ok=True)
                part = where.area / f".{name}.part-{os.getpid()}-{secrets.token_hex(4)}"  # never an earlier one's
                with open(local, "rb") as src, open(part, "xb") as dst:
                    shutil.copyfileobj(src, dst, BLOCK)
                    dst.flush()
                    os.fsync(dst.fileno())
                back, _ = sha256_file(part)
                if back != sha:
                    raise Corrupt(f"{part} read back as sha256 {back}, not {sha}; it was not renamed into place")
                if final.exists():  # under the lock only another host could have; never replace it
                    raise ResultsError(f"{final} appeared while this flush ran; {part} is left for the PI")
                publish(part, final)
                _fsync(where.area)
                state = "written"
            for line in todo:
                _drop(where.spool, line["run_id"])
            return {**out, "state": state}
        finally:
            local.unlink(missing_ok=True)


# --- get --------------------------------------------------------------------------------------------------------


def _copy_archive(line: dict, dest: Path) -> None:
    """Copy the run's archive named by `line` to `dest`, hashing it; Corrupt unless it is the line's bytes."""
    src_path = Path(line["where"])
    with open(dest, "xb") as out:
        if line["kept"] == "bundle":
            with tarfile.open(src_path, "r:") as tf:
                m = tf.getmember(f"{line['run_id']}.tar")
                src = _Hashing(tf.extractfile(m))
                shutil.copyfileobj(src, out, BLOCK)
        else:
            with open(src_path, "rb") as fh:
                src = _Hashing(fh)
                shutil.copyfileobj(src, out, BLOCK)
    if (src.h.hexdigest(), src.n) != (line["sha256"], line["size"]):
        raise Corrupt(f"the archive of run {line['run_id']} in {src_path} has sha256 {src.h.hexdigest()} and {src.n} "
                      f"bytes, not {line['sha256']} and {line['size']} as its manifest line says")


def _safe(name: str) -> bool:
    parts = name.split("/")
    return bool(name) and not name.startswith("/") and all(p not in ("", ".", "..") for p in parts)


def _restore(archive: Path, into: Path) -> int:
    with tarfile.open(archive, "r:") as tf:
        first = tf.next()
        if first is None or first.name != FILES or not first.isfile():
            raise Corrupt(f"{archive}: its first member is not {FILES}")
        raw = tf.extractfile(first).read()
        rows = raw.decode("utf-8").splitlines()
        if rows[:2] != [ARCHIVE_HEADER, "\t".join(FILE_COLUMNS)]:
            raise Corrupt(f"{archive}: {FILES} is not a bio-box file list")
        want = {}
        for r in rows[2:]:
            n, s, h = r.split("\t")
            want[n] = (int(s), h)
        (into / FILES).write_bytes(raw)
        seen = set()
        while (m := tf.next()) is not None:
            if m.name not in want or m.name in seen or not m.isfile() or not _safe(m.name):
                raise Corrupt(f"{archive}: member {m.name!r} is not one {FILES} lists once as a file")
            seen.add(m.name)
            target = into / m.name
            target.parent.mkdir(parents=True, exist_ok=True)
            with open(target, "xb") as out:
                src = _Hashing(tf.extractfile(m))
                shutil.copyfileobj(src, out, BLOCK)
            if (src.n, src.h.hexdigest()) != want[m.name]:
                raise Corrupt(f"{archive}: {m.name} has {src.n} bytes and sha256 {src.h.hexdigest()}, not "
                              f"{want[m.name][0]} and {want[m.name][1]} as {FILES} says")
        if seen != set(want):
            raise Corrupt(f"{archive}: {FILES} lists files the archive does not hold: {sorted(set(want) - seen)[:5]}")
    return len(seen)


def get(where: Places, run_id: str, out: Path) -> tuple[dict, int]:
    """Restore run `run_id` into the new or empty directory `out`, every sha256 checked. Returns (its manifest line,
    files restored)."""
    if not RUN_ID_RE.match(run_id):
        raise ResultsError(f"{run_id!r} is not a run id (32 lowercase hex digits)")
    found = _known(where, run_id)
    if not found:
        raise ResultsError(f"run {run_id} is in no bundle of {where.area} and not in the spool {where.spool}")
    line = found[0]  # a bundle before the spool (listing's order)
    out = Path(out)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ResultsError(f"{out} is not empty: results get restores into a new or empty directory")
    tmp = out.parent / f".{out.name}.tmp-{os.getpid()}"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    try:
        arch = tmp.parent / f".{out.name}.archive-{os.getpid()}.tar"
        try:
            _copy_archive(line, arch)
            n = _restore(arch, tmp)
        finally:
            arch.unlink(missing_ok=True)
        os.replace(tmp, out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return line, n

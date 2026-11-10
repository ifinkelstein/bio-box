"""Per-release search pack: build and verify (SEARCH_PLAN SD-06, SD-07; data_layer_and_tests.md §2B.2; C-12).

    <pack dir>/release=<R>/<pack_id>/
      pack_manifest.json         shards, meta files, counts, sources, rules, builder versions, and the
                                 tape_version the region flanks (feature uids) refer to
      unique/part-<k>.faa.zst    every exact-unique protein, named by protein_uid hex
      igr/part-<k>.fna.zst       IGRs >= 30 nt (+ strand), named by region_uid hex
      igr/igr_meta.parquet       search_region rows for the IGRs
      dark_orf/part-<k>.faa.zst  dark ORFs (biobox.search.orfs), named by region_uid hex
      dark_orf/dark_orf_meta.parquet
      models/<library>.hmm       the release's HMM libraries, re-written as plain HMMER3 text

S0 targets are `unique`, `igr`, `dark_orf` and `models` (ST0: no reduction). `spacers` waits
for CRISPR calls (O-ARRAY) and is listed under `missing_targets`; `nr50`, `nr90`, `segments`
and `igr90` come from S1, and the DIAMOND index of `unique` when DIAMOND enters the image.

A shard is cut once it holds `shard_text_bytes` of uncompressed FASTA (2 GiB: about 1 GiB of
protein and 0.6 GiB of nucleotide zstd), so every shard but the last of a target is over
256 MB. `pack_id` is the BLAKE3-128 of `pack_manifest.json`; the manifest holds no
timestamps, so the same inputs give the same pack_id (D-01).

The pack lives outside the data root (card W-03), in a pack directory (`bio-box pack build
--pack-dir DIR`; by default `<profile tmp>/packs` with a host profile, else beside the search work
directory, `runner.default_pack_base`); a pack directory inside the root is refused. It is a derived,
rebuildable cache of a release, not a release file: the same inputs give the same bytes and the same
`pack_id`, its own manifest records what it was built from, and `verify` (C-12) checks it. It is kept
off the data root because on shared storage every release's pack has files under the 64 MiB floor
(LY-04; planning/measurements/2026-11-10_search_pack_size.md). `build` writes the pack under a
temporary name in `release=<R>/` and renames it to its `pack_id`; it never overwrites a pack, and it is a
no-op, before any work, when a pack built from the same inputs (`_fingerprint`) is there and verifies.

A pack holds the release's files only (card B-28): every input is listed by `catalog.ReleaseFiles`,
from the release manifest when there is one, else from the stages of the release's `release_def`, so
a stage on disk the release does not list never enters its pack.
"""

from __future__ import annotations

import json
import os
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import ids, layout, provenance
from biobox.catalog import ReleaseFiles
from biobox.l0seq.store import SeqStore
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l1.igr import igr_version
from biobox.l3.geometry import junction
from biobox.l3.tapes import tape_version
from biobox.registry.schemas import table
from biobox.search import fastz, orfs

PACK_FORMAT = "biobox-pack-1"
SHARD_TEXT_BYTES = 2 * 2**30
IGR_MIN_NT = 30
Z = {"z_prot": domains.Z_PROT, "dom_z": domains.DOM_Z, "z_nt_mb": 1e5}  # SD-04, DEC-04
TARGET_KIND = {"unique": "protein", "igr": "igr", "dark_orf": "dark_orf"}
MISSING_AT_S0 = {"spacers": "no CRISPR spacer calls in this release (O-ARRAY not built)"}


MANIFEST = "pack_manifest.json"
# the manifest keys that, with each target's `sources`, say what a pack was built from (`_fingerprint`)
FINGERPRINT_KEYS = ("pack_format", "release_id", "gene_call_version", "igr_version", "tape_version", "z", "rules",
                    "builder", "missing_targets")


class NoPack(FileNotFoundError):
    """No pack of the release in the pack directory; the message names the command that builds one."""


def release_dir(packs: Path, release: str) -> Path:
    """Where the packs of a release live in the pack directory `packs`."""
    return Path(packs) / f"release={release}"


def pack_dir(packs: Path, release: str, pack_id: bytes | str) -> Path:
    """The pack of `release` with `pack_id` in the pack directory `packs`: `<packs>/release=<R>/<pack_id>`."""
    return release_dir(packs, release) / (pack_id if isinstance(pack_id, str) else pack_id.hex())


def outside_root(root: Path, path: Path, what: str) -> Path:
    """`path`, after refusing (ValueError) a path that is the data root or lies under it (card W-03: nothing in
    the search path writes under the root)."""
    r, p = Path(root).resolve(), Path(path).resolve()
    if p == r or p.is_relative_to(r):
        raise ValueError(f"{what} {path} lies under the data root {root}; the search path writes nothing under the "
                         "root: give a directory outside it")
    return Path(path)


def build_command(root: Path, release: str, packs: Path) -> str:
    return (f"bio-box pack build --root {root} --release {release} --pack-dir {packs} "
            "--model <library>-<version>=<HMM file> ...")


def find(packs: Path, release: str, pack_id: bytes | str | None = None, *, root: Path | None = None) -> Path:
    """The pack directory of `release` in `packs` (of `pack_id`, when given). NoPack, naming the command that
    builds it, when there is none; ValueError when several packs of the release are there and none is named."""
    if pack_id is not None:
        d = pack_dir(packs, release, pack_id)
        found = [d] if (d / MANIFEST).is_file() else []
    else:
        base = release_dir(packs, release)
        found = sorted(p for p in base.iterdir() if not p.name.startswith(".") and (p / MANIFEST).is_file()) \
            if base.is_dir() else []
    if not found:
        what = "search pack " if pack_id is None else \
            f"search pack {pack_id if isinstance(pack_id, str) else pack_id.hex()} "
        raise NoPack(f"no {what}of release {release} in the pack directory {packs}; build it with: "
                     f"{build_command(root or '<data root>', release, packs)}")
    if len(found) > 1:
        raise ValueError(f"{len(found)} search packs of release {release} in {release_dir(packs, release)} "
                         f"({', '.join(p.name for p in found)}); keep the one to search and remove the others")
    return found[0]


def _rel(root: Path, p: Path) -> str:
    return p.relative_to(root).as_posix()


def _source_hash(root: Path, paths: list[Path]) -> dict:
    h = blake3()
    for p in sorted(paths):
        h.update(f"{_rel(root, p)}:{ids.sha256_file(p).hex()}\n".encode())
    return {"files": len(paths), "blake3": h.hexdigest()}


def _protein_files(rf: ReleaseFiles, gcv: str) -> list[Path]:
    return rf.table("protein_seq", gene_call_version=gcv)


def _occurrence_files(rf: ReleaseFiles, gcv: str) -> list[Path]:
    return rf.table("protein_occurrence", gene_call_version=gcv)


def _igr_files(rf: ReleaseFiles, iv: str) -> list[Path]:
    """The IGR files of igr version iv (`<tape_version>-igr1`, card B-20)."""
    return rf.table("igr", igr_version=iv)


def _tape_files(rf: ReleaseFiles, tv: str) -> list[Path]:
    return rf.table("contig_tapes", tape_version=tv)


def manifest_tape_version(manifest: dict) -> str:
    """The tape version a pack's flanks refer to, and so the one its runs read (card B-16); a manifest
    from before B-16 names none and was built from the gene-call producer's."""
    return manifest.get("tape_version") or tape_version([manifest["gene_call_version"]])


def _contigs(rf: ReleaseFiles, tv: str) -> tuple[dict[bytes, int], set[bytes], dict[bytes, int]]:
    """Contig lengths, the circular contigs and the strand of every slot's feature, from the tapes of tv."""
    lengths, circular, strands = {}, set(), {}
    for f in _tape_files(rf, tv):
        for r in pq.read_table(f, columns=["contig_uid", "contig_len", "is_circular", "feature_uids",
                                           "strands"]).to_pylist():
            lengths[r["contig_uid"]] = r["contig_len"]
            if r["is_circular"]:
                circular.add(r["contig_uid"])
            strands.update(zip(r["feature_uids"], r["strands"]))
    return lengths, circular, strands


def _slots(rf: ReleaseFiles, gcv: str) -> dict[bytes, list[dict]]:
    """The genes of each contig in coordinate order (start, end, strand)."""
    by_contig: dict[bytes, list[dict]] = defaultdict(list)
    for f in _occurrence_files(rf, gcv):
        for r in pq.read_table(f, columns=["contig_uid", "start", "end", "strand"]).to_pylist():
            by_contig[r["contig_uid"]].append(r)
    for v in by_contig.values():
        v.sort(key=lambda r: (r["start"], r["end"], r["strand"]))
    return by_contig


def _codes(rf: ReleaseFiles, gcv: str) -> dict[bytes, int]:
    """Genetic code per contig: the commonest code among its CDS features (11 for a contig with none)."""
    counts: dict[bytes, Counter] = defaultdict(Counter)
    for f in rf.table("feature_protein", gene_call_version=gcv):
        for r in pq.read_table(f, columns=["contig_uid", "genetic_code"]).to_pylist():
            counts[r["contig_uid"]][r["genetic_code"]] += 1
    return {cu: c.most_common(1)[0][0] for cu, c in counts.items()}


def _check_codes(rf: ReleaseFiles, gcv: str) -> None:
    """Fail fast, before anything is written, on a genetic code the dark-ORF stage cannot translate."""
    for cu, code in sorted(_codes(rf, gcv).items()):
        try:
            orfs.check_code(code)
        except ValueError as e:
            raise ValueError(f"{e}; first contig using it: {cu.hex()}") from None


def _gc(seq: bytes) -> float:
    s = seq.upper()
    acgt = sum(s.count(b) for b in (b"A", b"C", b"G", b"T"))
    return (s.count(b"G") + s.count(b"C")) / acgt if acgt else 0.0


def _region(kind: str, cu: bytes, start: int, end: int, strand: int, n_residues: int, **kw) -> dict:
    row = {"region_uid": ids.region_uid(cu, start, end, strand), "target_kind": kind, "contig_uid": cu,
           "start": start, "end": end, "strand": strand, "n_residues": n_residues,
           "orientation": None, "gc": None, "dist_to_contig_end_bp": None,
           "partial": False, "at_contig_end": False, "gene_overlap_frac": None, "genetic_code": None,
           "left_feature_uid": None, "right_feature_uid": None}
    row.update(kw)
    return row


def _write_meta(rows: list[dict], path: Path) -> None:
    rows.sort(key=lambda r: r["region_uid"])
    pq.write_table(pa.Table.from_pylist(rows, schema=table("search_region").schema), path, compression="zstd",
                   row_group_size=122_880)


def _build_unique(rf: ReleaseFiles, out: Path, gcv: str, shard_bytes: int) -> list[dict]:
    prots = {}
    for f in _protein_files(rf, gcv):
        t = pq.read_table(f, columns=["protein_uid", "seq"])
        prots.update(zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist()))
    return fastz.write_records(out / "unique", ".faa.zst",
                               ((u.hex(), prots[u].replace("*", "")) for u in sorted(prots)), shard_bytes)


def _build_igr(rf: ReleaseFiles, out: Path, iv: str, store: SeqStore, lengths, circular, strands, shard_bytes):
    meta, records = [], []
    for f in _igr_files(rf, iv):
        for r in pq.read_table(f, partitioning=None).to_pylist():
            cu, s, e = r["contig_uid"], r["start"], r["end"]
            if e - s < IGR_MIN_NT:
                continue
            L, circ = lengths[cu], cu in circular
            seq = store.get(cu, s, e, 1, circular=circ)
            lf, rf = r["left_feature_uid"], r["right_feature_uid"]
            # flanks of an IGR never nest, so only their strands matter
            orient = junction((0, 1, strands[lf]), (2, 3, strands[rf])) if lf and rf else None
            row = _region("igr", cu, s, e, 1, e - s, left_feature_uid=lf, right_feature_uid=rf,
                          orientation=orient, gc=_gc(seq), dist_to_contig_end_bp=None if circ else min(s, L - e),
                          partial=not circ and (s == 0 or e == L))
            meta.append(row)
            records.append((row["region_uid"].hex(), seq.decode("ascii")))
    records.sort()
    shards = fastz.write_records(out / "igr", ".fna.zst", records, shard_bytes)
    _write_meta(meta, out / "igr/igr_meta.parquet")
    return shards


def _build_dark_orf(rf: ReleaseFiles, out: Path, gcv: str, store: SeqStore, lengths, circular, slots, shard_bytes):
    codes = _codes(rf, gcv)
    meta, records = [], []
    for cu in sorted(lengths):
        circ = cu in circular
        seq = store.get(cu)
        genes = [(r["start"], r["end"], r["strand"]) for r in slots.get(cu, [])]
        code = codes.get(cu, 11)
        for o, frac, at_end in orfs.dark_orfs(seq, circ, genes, code):
            row = _region("dark_orf", cu, o.start, o.end, o.strand, len(o.aa), partial=o.partial,
                          at_contig_end=at_end, gene_overlap_frac=frac, genetic_code=code)
            meta.append(row)
            records.append((row["region_uid"].hex(), o.aa))
    records.sort()
    shards = fastz.write_records(out / "dark_orf", ".faa.zst", records, shard_bytes)
    _write_meta(meta, out / "dark_orf/dark_orf_meta.parquet")
    return shards


def _build_models(out: Path, models: dict[str, Path]) -> list[dict]:
    shards = []
    (out / "models").mkdir(parents=True, exist_ok=True)
    for lib in sorted(models):
        path = out / f"models/{lib}.hmm"
        n, m = 0, 0
        with open(path, "wb") as fh:
            for hmm in domains.iter_hmms(models[lib]):
                hmm.write(fh)
                n += 1
                m += hmm.M
        shards.append({"path": path, "n_seq": n, "n_residues": m})
    return shards


def _fingerprint(manifest: dict) -> dict:
    """What a pack was built from: the build's rules, versions and tools, and each target's sources (input
    files by path and sha256, library files by sha256). The pack's bytes are a function of these (D-01)."""
    return {**{k: manifest.get(k) for k in FINGERPRINT_KEYS},
            "sources": {t: info.get("sources") for t, info in sorted(manifest.get("targets", {}).items())}}


def build(root: Path, release: str, gene_call_version: str = GENE_CALL_VERSION,
          models: dict[str, Path] | None = None, shard_text_bytes: int = SHARD_TEXT_BYTES,
          producers: list[str] | None = None, release_manifest: dict | None = None, *, packs: Path) -> dict:
    """Build the release's pack from the data root into the pack directory `packs` (card W-03), at
    `<packs>/release=<R>/<pack_id>/`; returns the manifest (with `pack_id`).

    `packs` may not lie under the root (ValueError, before anything is read). When a pack of the release
    built from the same inputs is there and verifies, nothing is built or written and its manifest is
    returned. Otherwise the pack is written under a temporary name in `release=<R>/` and renamed to its
    pack_id; if a pack of that id is there already (and does not verify), FileExistsError: a pack is never
    overwritten.

    `producers` (default: the gene-call producer alone) name the feature set: the IGRs are those of
    its igr version and the region flanks are features of its tape version, which the manifest
    records as `tape_version` (card B-16). Every input file is the release's (card B-28): from
    `release_manifest` (default: `manifests/release=<R>.json` when it exists), else from the stages of
    the release's `release_def` (FileNotFoundError when it has neither)."""
    root = Path(root).resolve()
    packs = outside_root(root, packs, "the pack directory")
    rf = ReleaseFiles.of(root, release, release_manifest)
    _check_codes(rf, gene_call_version)
    gcv = gene_call_version
    producers = sorted(producers) if producers is not None else [gcv]
    tv, iv = tape_version(producers), igr_version(producers)
    seq_files = rf.artifact("seq_shards")
    sources = {
        "unique": {"protein_seq": _source_hash(root, _protein_files(rf, gcv))},
        "igr": {"igr": _source_hash(root, _igr_files(rf, iv)), "seq_shards": _source_hash(root, seq_files)},
        "dark_orf": {"protein_occurrence": _source_hash(root, _occurrence_files(rf, gcv)),
                     "contig_tapes": _source_hash(root, _tape_files(rf, tv)),
                     "seq_shards": _source_hash(root, seq_files)},
        "models": {lib: {"files": 1, "sha256": ids.sha256_file(p).hex()} for lib, p in sorted((models or {}).items())},
    }
    head = {
        "pack_format": PACK_FORMAT,
        "release_id": release,
        "gene_call_version": gcv,
        "igr_version": iv,
        "tape_version": tv,  # the tape version the region flanks refer to, and the one a run reads
        "z": Z,
        "rules": {"igr_min_nt": IGR_MIN_NT, "dark_orf_min_aa": orfs.MIN_AA, "dark_orf_end_bp": orfs.END_BP,
                  "dark_orf_max_gene_overlap": orfs.MAX_GENE_OVERLAP, "shard_text_bytes": shard_text_bytes},
        "builder": provenance.tool_versions("bio-box", "pyhmmer", "pyarrow"),
        "missing_targets": MISSING_AT_S0,
    }
    want = _fingerprint({**head, "targets": {t: {"sources": s} for t, s in sources.items()}})
    base = release_dir(packs, release)
    for d in sorted(p for p in base.iterdir() if not p.name.startswith(".")) if base.is_dir() else []:
        if (d / MANIFEST).is_file():
            have, pid = read_manifest(d)
            if _fingerprint(have) == want and not verify_dir(root, release, d, release_manifest=release_manifest):
                return {**have, "pack_id": pid.hex()}  # built from the same inputs: a no-op
    base.mkdir(parents=True, exist_ok=True)
    out = base / f".building-{os.getpid()}-{os.urandom(4).hex()}"
    out.mkdir()
    try:
        manifest, data = _build_into(out, root, rf, head, sources, gcv, iv, tv, models or {}, shard_text_bytes)
        final = pack_dir(packs, release, pack_id_of(data))
        if final.exists():
            if not verify_dir(root, release, final, release_manifest=release_manifest):
                shutil.rmtree(out)
                return {**manifest, "pack_id": final.name}
            raise FileExistsError(f"{final}: a pack with this pack_id exists and does not verify; a pack is never "
                                  "overwritten: remove it, then build again")
        os.rename(out, final)
    except BaseException:
        shutil.rmtree(out, ignore_errors=True)
        raise
    return {**manifest, "pack_id": final.name}


def _build_into(out: Path, root: Path, rf: ReleaseFiles, head: dict, sources: dict, gcv: str, iv: str, tv: str,
                models: dict[str, Path], shard_text_bytes: int) -> tuple[dict, bytes]:
    """Write the pack's files and its manifest into the empty directory `out`; (manifest, manifest bytes)."""
    lengths, circular, strands = _contigs(rf, tv)
    slots = _slots(rf, gcv)
    # a protein-only root (a synthetic world) has no sequence store: its nucleotide targets are empty
    store = SeqStore(root, rf) if rf.table("seq_locator") else None
    if store is None:
        lengths = {}
    try:
        built = {
            "unique": _build_unique(rf, out, gcv, shard_text_bytes),
            "igr": _build_igr(rf, out, iv, store, lengths, circular, strands, shard_text_bytes),
            "dark_orf": _build_dark_orf(rf, out, gcv, store, lengths, circular, slots, shard_text_bytes),
            "models": _build_models(out, models),
        }
    finally:
        if store is not None:
            store.close()
    manifest = {
        **head,
        "targets": {t: {"kind": TARGET_KIND.get(t, t), "n_seq": sum(s["n_seq"] for s in sh),
                        "n_residues": sum(s["n_residues"] for s in sh), "sources": sources[t]}
                    for t, sh in built.items()},
        "shards": [{"target": t, "path": _rel(out, s["path"]), "sha256": ids.sha256_file(s["path"]).hex(),
                    "bytes": s["path"].stat().st_size, "n_seq": s["n_seq"], "n_residues": s["n_residues"]}
                   for t, sh in built.items() for s in sh],
        "meta": [{"target": t, "path": f"{t}/{t}_meta.parquet",
                  "sha256": ids.sha256_file(out / f"{t}/{t}_meta.parquet").hex(),
                  "bytes": (out / f"{t}/{t}_meta.parquet").stat().st_size} for t in ("igr", "dark_orf")],
    }
    data = json.dumps(manifest, sort_keys=True, indent=1).encode()
    (out / MANIFEST).write_bytes(data)
    return manifest, data


def pack_id_of(manifest_bytes: bytes) -> bytes:
    return blake3(manifest_bytes).digest(length=16)


def read_manifest(pdir: Path) -> tuple[dict, bytes]:
    """(manifest, pack_id) of the pack in the directory `pdir` (`find`, `pack_dir`)."""
    data = (Path(pdir) / MANIFEST).read_bytes()
    return json.loads(data), pack_id_of(data)


def _count_source(rf: ReleaseFiles, manifest: dict, target: str, out: Path) -> int | None:
    gcv = manifest["gene_call_version"]
    if target == "unique":
        return len({u for f in _protein_files(rf, gcv)
                    for u in pq.read_table(f, columns=["protein_uid"]).column("protein_uid").to_pylist()})
    if target == "igr":
        return sum(1 for f in _igr_files(rf, manifest.get("igr_version") or igr_version([gcv]))
                   for r in pq.read_table(f, columns=["start", "end"]).to_pylist()
                   if r["end"] - r["start"] >= manifest["rules"]["igr_min_nt"])
    if target == "dark_orf":
        return pq.read_metadata(out / "dark_orf/dark_orf_meta.parquet").num_rows
    return None


def verify(root: Path, release: str, pack_id: bytes | None = None,
           min_shard_bytes: int = layout.TARGET_FILE_BYTES, release_manifest: dict | None = None, *,
           packs: Path) -> list[str]:
    """C-12: problems with the release's pack in the pack directory `packs` (the one of `pack_id`, when
    given; `find`); empty means it verifies. Sequence counts are checked against the release's files under
    `root` (`release_manifest`, as for `build`)."""
    try:
        out = find(packs, release, pack_id, root=root)
    except (FileNotFoundError, ValueError) as e:
        return [str(e)]
    return verify_dir(root, release, out, pack_id, min_shard_bytes, release_manifest)


def verify_dir(root: Path, release: str, out: Path, pack_id: bytes | None = None,
               min_shard_bytes: int = layout.TARGET_FILE_BYTES, release_manifest: dict | None = None) -> list[str]:
    """`verify` of the pack in the directory `out`: also, the directory is named by the manifest's hash."""
    root, out = Path(root).resolve(), Path(out)
    if not (out / MANIFEST).exists():
        return [f"{out}: no {MANIFEST}"]
    try:
        rf = ReleaseFiles.of(root, release, release_manifest)
    except FileNotFoundError as e:
        return [f"release {release}: {e}"]
    manifest, pid = read_manifest(out)
    problems = []
    if pack_id is not None and pack_id != pid:
        problems.append(f"pack_id {pack_id.hex()} differs from the manifest hash {pid.hex()}")
    if out.name != pid.hex():
        problems.append(f"{out}: the directory is not named by the manifest hash {pid.hex()}")
    if manifest.get("release_id") != release:
        problems.append(f"{out}: a pack of release {manifest.get('release_id')}, not {release}")
    listed = {MANIFEST}
    by_target: dict[str, list[dict]] = defaultdict(list)
    for entry in manifest["shards"] + manifest["meta"]:
        listed.add(entry["path"])
        p = out / entry["path"]
        if not p.exists():
            problems.append(f"{entry['path']}: listed but missing")
            continue
        if p.stat().st_size != entry["bytes"]:
            problems.append(f"{entry['path']}: {p.stat().st_size} bytes, manifest says {entry['bytes']}")
        if ids.sha256_file(p).hex() != entry["sha256"]:
            problems.append(f"{entry['path']}: sha256 differs from the manifest")
    for entry in manifest["shards"]:
        by_target[entry["target"]].append(entry)
    for p in sorted(x for x in out.rglob("*") if x.is_file()):
        if _rel(out, p) not in listed:
            problems.append(f"{_rel(out, p)}: on disk but not in the manifest")
    for target, shards in sorted(by_target.items()):
        info = manifest["targets"][target]
        if sum(s["n_seq"] for s in shards) != info["n_seq"]:
            problems.append(f"{target}: shard counts do not sum to the target's n_seq")
        # models are one file per library (small by nature, SEARCH_PLAN §3.2), not size-cut shards
        for s in (sorted(shards, key=lambda s: s["path"])[:-1] if target != "models" else []):
            if s["bytes"] < min_shard_bytes:
                problems.append(f"{s['path']}: {s['bytes']} bytes, under {min_shard_bytes} and not the last shard")
        if target == "models":
            for s in shards:
                p = out / s["path"]
                if p.exists() and sum(1 for _ in domains.iter_hmms(p)) != s["n_seq"]:
                    problems.append(f"{s['path']}: model count differs from the manifest")
            continue
        n_records = 0
        for s in shards:
            p = out / s["path"]
            if p.exists():
                n = sum(1 for _ in fastz.read_records(p))
                n_records += n
                if n != s["n_seq"]:
                    problems.append(f"{s['path']}: {n} records, manifest says {s['n_seq']}")
        expect = _count_source(rf, manifest, target, out)
        if expect is not None and expect != info["n_seq"]:
            problems.append(f"{target}: {info['n_seq']} sequences, its source has {expect}")
    return problems

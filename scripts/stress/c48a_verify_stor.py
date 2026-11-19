#!/usr/bin/env python3
"""Card C-48a: verify /stor/work JGI/NERSC/Zenodo archive inventory and sources.yaml locks.

Checks >= 64 MiB floor, absence of .cookie.cfg / temp files, archive member headers in
central-directory / streaming mode (no extraction), sha256 digest format, and sources.yaml
consistency.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
import zipfile
from pathlib import Path
from typing import Any, Sequence

import yaml

MIN_STOR_BYTES = 67_108_864  # 64 MiB (AGENTS.md §7.3)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

JGI_LOCK_SOURCES = (
    "jgi_gold_downloads",
    "jgi_img_m_tables",
    "jgi_img_vr",
    "jgi_img_pr",
    "jgi_img_umag",
    "jgi_iuig",
    "jgi_inovirus",
    "jgi_gems_nersc",
    "gvmags_v2",
)


def validate_sha256(digest: object) -> bool:
    """True iff digest is a 64-character lowercase hexadecimal SHA-256 string."""
    return isinstance(digest, str) and bool(SHA256_RE.match(digest))


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def inspect_archive_members(path: Path) -> list[dict[str, object]]:
    """List archive members (name, uncompressed size) without extracting payloads.

    Raises ValueError if any archive member looks like a cookie/secret file.
    """
    members: list[dict[str, object]] = []
    suffix = path.suffix.lower()
    if suffix == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            for info in zf.infolist():
                if info.is_dir():
                    continue
                if "cookie" in info.filename.lower():
                    raise ValueError(f"{path}: forbidden cookie member {info.filename!r}")
                members.append({"name": info.filename, "size": info.file_size})
    elif suffix == ".tar":
        with tarfile.open(path, "r:") as tf:
            for mem in tf.getmembers():
                if not mem.isfile():
                    continue
                if "cookie" in mem.name.lower():
                    raise ValueError(f"{path}: forbidden cookie member {mem.name!r}")
                members.append({"name": mem.name, "size": mem.size})
    else:
        raise ValueError(f"{path}: unsupported archive extension {suffix!r}")
    return members


def check_stor_dir(stor_dir: Path, *, min_bytes: int = MIN_STOR_BYTES) -> list[str]:
    """Audit a /stor/work archive directory for >= 64 MiB floor and forbidden files."""
    issues: list[str] = []
    if not stor_dir.is_dir():
        return [f"missing directory: {stor_dir}"]
    for p in sorted(stor_dir.iterdir()):
        if not p.is_file():
            issues.append(f"unexpected non-file entry: {p.name}")
            continue
        lower = p.name.lower()
        if "cookie" in lower:
            issues.append(f"forbidden cookie file present: {p.name}")
        if lower.startswith(".") or lower.endswith(".part") or lower.endswith(".tmp"):
            issues.append(f"temp/hidden file present: {p.name}")
        sz = p.stat().st_size
        if sz < min_bytes:
            issues.append(f"{p.name}: {sz} bytes is below 64 MiB ({min_bytes} bytes) floor")
    return issues


def load_stor_sha256sums(tsv_path: Path) -> dict[str, dict[str, object]]:
    """Parse STOR_SHA256SUMS.tsv (timestamp, filename, bytes, sha256)."""
    out: dict[str, dict[str, object]] = {}
    for line in tsv_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        ts, name, sz_s, sha = line.split("\t")
        out[name] = {"timestamp": ts, "bytes": int(sz_s), "sha256": sha}
    return out


def verify_sources_yaml(
    doc: dict[str, Any],
    *,
    required_names: Sequence[str] = JGI_LOCK_SOURCES,
    min_bytes: int = MIN_STOR_BYTES,
    check_files_on_disk: bool = True,
) -> list[str]:
    """Verify that required sources.yaml entries lock /stor/work archives with bytes & sha256."""
    issues: list[str] = []
    by_name = {entry.get("name"): entry for entry in doc.get("sources", []) if isinstance(entry, dict)}
    for name in required_names:
        entry = by_name.get(name)
        if entry is None:
            issues.append(f"sources.yaml missing entry: {name}")
            continue
        stored = entry.get("stored")
        if not isinstance(stored, list) or not stored:
            issues.append(f"{name}: missing non-empty 'stored' list")
            continue
        for item in stored:
            if not isinstance(item, dict):
                issues.append(f"{name}: malformed stored item {item!r}")
                continue
            p_str = item.get("path")
            sz = item.get("bytes")
            sha = item.get("sha256")
            if not isinstance(p_str, str) or not p_str:
                issues.append(f"{name}: stored item missing path")
                continue
            if not isinstance(sz, int) or sz < min_bytes:
                issues.append(f"{name}: {p_str} bytes={sz!r} below floor {min_bytes}")
            if not validate_sha256(sha):
                issues.append(f"{name}: {p_str} invalid sha256={sha!r}")
            if check_files_on_disk:
                p = Path(p_str)
                if not p.is_file():
                    issues.append(f"{name}: stored file missing on disk: {p}")
                elif isinstance(sz, int) and p.stat().st_size != sz:
                    issues.append(f"{name}: {p} size mismatch (yaml={sz}, disk={p.stat().st_size})")
    return issues


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stor-dir", type=Path, required=True)
    ap.add_argument("--sums-tsv", type=Path, default=None)
    ap.add_argument("--sources-yaml", type=Path, default=None)
    ap.add_argument("--hash-max-bytes", type=int, default=0, help="Recompute sha256 for files <= this size")
    args = ap.parse_args(argv)

    issues = check_stor_dir(args.stor_dir)
    sums = load_stor_sha256sums(args.sums_tsv) if args.sums_tsv and args.sums_tsv.is_file() else {}
    inventory = []

    for p in sorted(args.stor_dir.iterdir()):
        if not p.is_file():
            continue
        sz = p.stat().st_size
        members = inspect_archive_members(p)
        rec: dict[str, object] = {
            "name": p.name,
            "path": str(p),
            "bytes": sz,
            "member_count": len(members),
            "members_total_bytes": sum(int(m["size"]) for m in members),
            "members": members if len(members) <= 30 else members[:10],
        }
        if p.name in sums:
            rec["sha256_tsv"] = sums[p.name]["sha256"]
            rec["timestamp_tsv"] = sums[p.name]["timestamp"]
            if sums[p.name]["bytes"] != sz:
                issues.append(f"{p.name}: STOR_SHA256SUMS bytes={sums[p.name]['bytes']} != disk={sz}")
        if args.hash_max_bytes > 0 and sz <= args.hash_max_bytes:
            computed = compute_sha256(p)
            rec["sha256_computed"] = computed
            if p.name in sums and sums[p.name]["sha256"] != computed:
                issues.append(f"{p.name}: sha256 mismatch (tsv={sums[p.name]['sha256']} != computed={computed})")
        inventory.append(rec)

    if args.sources_yaml:
        doc = yaml.safe_load(args.sources_yaml.read_text(encoding="utf-8"))
        issues.extend(verify_sources_yaml(doc))

    print(json.dumps({"stor_dir": str(args.stor_dir), "archives": inventory, "issues": issues}, indent=2))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())

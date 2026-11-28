#!/usr/bin/env python3
"""Card R-26: synchronize and verify `sources.yaml`, `planning/CATALOG_GOAL.md`, and `STAGING_WORKER_RULES.md`.

Provides:
  - Read-only summarization of `tower` `SOURCE_URLS.tsv` + `SHA256SUMS` trees (`sources_incoming`, `sources_incoming_s2`,
    `sources`), including nested parent+child trees and active `.partial` transfer detection.
  - Surgical reconciliation of `sources.yaml` `fetched:` and `stored:` fields that preserves every `url*:` line and
    existing `approved:` string verbatim.
  - Repository ledger verification (`verify` subcommand):
      1. `sources.yaml` parses cleanly as YAML, has unique entry names, valid `redistributable` flags, valid URL host
         segments (no `{var}` in host), zero URL-form regressions against a baseline if given, and zero stale
         `"nothing yet"` / `"will be tower:"` fields on entries with landed `tower` files when an inventory is given.
      2. Every `*.conf` in `scripts/stress/build_stage_pod/stages/` and `stages/templates/` has a matching stage row in
         `planning/CATALOG_GOAL.md` §3 (extra pre-merge rows such as `PROMGE1` and `BBSAG1` are allowed), and every
         stage's `left_out.tsv` / deferred remainder is registered in `planning/CATALOG_GOAL.md` §4.1.
      3. `planning/tasks/STAGING_WORKER_RULES.md` requires merge-time / fetch-completion synchronization checked with
         `python scripts/orchestrator/r26_sync_sources_ledger.py verify`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCES = REPO_ROOT / "sources.yaml"
DEFAULT_CATALOG_GOAL = REPO_ROOT / "planning" / "CATALOG_GOAL.md"
DEFAULT_STAGES_DIR = REPO_ROOT / "scripts" / "stress" / "build_stage_pod" / "stages"
DEFAULT_WORKER_RULES = REPO_ROOT / "planning" / "tasks" / "STAGING_WORKER_RULES.md"

DEFAULT_REQUIRED_REMAINDERS: tuple[str, ...] = (
    "ET1",
    "O4",
    "GV2",
    "HROM1",
    "JGIUG1",
    "PROMGE1",
    "BBSAG1",
    "GlobDB",
    "GTDB",
    "ASG1",
    "DSV70",
    "MGnify",
)

_LEGACY_PROSE_URL_ENTRIES: frozenset[str] = frozenset({"inphared", "gtdb_metadata"})

_NAME_LINE_RE = re.compile(r"^  - name:[ \t]*(.+?)[ \t]*$")
_URL_LINE_RE = re.compile(r"^    (url[A-Za-z0-9_]*):[ \t]+(.+?)[ \t]*$")
_HTTP_FTP_RE = re.compile(r"^(https?|ftp)://")
_SCHEME_RE = re.compile(r"^([a-zA-Z][a-zA-Z0-9+.-]*)://([^/]+)(/.*)?$")
_TOWER_PATH_RE = re.compile(
    r"tower:/mnt/user/research-bulk/bio-box/(sources_incoming(?:_s2)?)/((?:[^{}\s,;\"')]|\{[^{}\s]+\})+)"
)
_ELLIPSIS_PATH_RE = re.compile(r"(?:^|[\s,;(\"'])\.\.\./([^\s,;\"')]+)")
_BACKTICK_STAGE_RE = re.compile(r"`([A-Z][A-Z0-9_]*)`")
_STAGE_RANGE_RE = re.compile(r"`([A-Z]+)(\d+)`\s*(?:to|\.\.)\s*`\1(\d+)`")


class LedgerError(ValueError):
    """Raised when `sources.yaml`, `CATALOG_GOAL.md`, or `STAGING_WORKER_RULES.md` fails a ledger check."""


Refused = LedgerError


def _strip_quotes(val: str) -> str:
    s = val.strip()
    if len(s) >= 2 and ((s[0] == '"' and s[-1] == '"') or (s[0] == "'" and s[-1] == "'")):
        return s[1:-1]
    return s


def _format_date_range(date_min: str | None, date_max: str | None) -> str:
    if not date_min and not date_max:
        return "unknown-date"
    d1 = date_min or date_max or "unknown-date"
    d2 = date_max or date_min or "unknown-date"
    return d1 if d1 == d2 else f"{d1}..{d2}"


def _parse_source_urls_file(tsv_path: Path) -> tuple[dict[str, dict[str, Any]], int]:
    by_file: dict[str, dict[str, Any]] = {}
    raw_rows = 0
    with tsv_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.rstrip("\r\n")
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if parts[0] == "filename" and len(parts) >= 2 and parts[1] == "url":
                continue
            if len(parts) < 6:
                continue
            raw_rows += 1
            fname, url, b_str, sha, lic, dt_utc = parts[:6]
            try:
                nbytes = int(b_str)
            except ValueError:
                nbytes = 0
            by_file[fname] = {
                "filename": fname,
                "url": url,
                "bytes": nbytes,
                "sha256": sha,
                "license": lic,
                "date_utc": dt_utc,
            }
    return by_file, raw_rows


def _parse_sha256sums_file(sums_path: Path) -> tuple[str, int, dict[str, str]]:
    raw = sums_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    entries: dict[str, str] = {}
    for line in raw.decode("utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split(None, 1)
        if len(parts) == 2:
            entries[parts[1].lstrip("* ")] = parts[0]
    return digest, len(raw), entries


def _summarize_single_dir(
    dir_path: Path,
    *,
    root_label: str,
    rel_dir: str,
) -> dict[str, Any]:
    """Summarize only the direct `SOURCE_URLS.tsv`, `SHA256SUMS`, and `.partial` files of `dir_path`."""
    rel = rel_dir.strip("/")
    filenames = sorted(p.name for p in dir_path.iterdir() if p.is_file())
    has_urls = "SOURCE_URLS.tsv" in filenames
    has_sums = "SHA256SUMS" in filenames
    has_partial = any(f.endswith(".partial") for f in filenames)

    by_f: dict[str, dict[str, Any]] = {}
    raw_n = 0
    if has_urls:
        by_f, raw_n = _parse_source_urls_file(dir_path / "SOURCE_URLS.tsv")

    sums_sha = None
    sums_bytes = 0
    sums_cnt = 0
    if has_sums:
        sums_sha, sums_bytes, sums_entries = _parse_sha256sums_file(dir_path / "SHA256SUMS")
        sums_cnt = len(sums_entries)

    dates = sorted(v["date_utc"] for v in by_f.values() if v["date_utc"])
    ts_min = dates[0] if dates else None
    ts_max = dates[-1] if dates else None
    dmin = ts_min[:10] if ts_min else None
    dmax = ts_max[:10] if ts_max else None

    single_file_sha256 = None
    single_filename = None
    if len(by_f) == 1:
        only = next(iter(by_f.values()))
        single_file_sha256 = only["sha256"]
        single_filename = only["filename"]

    return {
        "root_label": root_label,
        "rel_dir": rel,
        "stored_path": f"tower:/mnt/user/research-bulk/bio-box/{root_label}/{rel}/",
        "files": len(by_f),
        "raw_rows": raw_n,
        "bytes": sum(v["bytes"] for v in by_f.values()),
        "ts_min": ts_min,
        "ts_max": ts_max,
        "date_min": dmin,
        "date_max": dmax,
        "date_range": _format_date_range(dmin, dmax),
        "sha256sums_sha256": sums_sha,
        "sha256sums_bytes": sums_bytes,
        "sha256sums_count": sums_cnt,
        "top_sha256sums_sha256": sums_sha,
        "single_filename": single_filename,
        "single_file_sha256": single_file_sha256,
        "streaming": has_partial,
        "subdirs": [],
    }


def summarize_tower_dir(
    dir_path: Path,
    *,
    root_label: str = "sources_incoming",
    rel_dir: str | None = None,
) -> dict[str, Any]:
    """Summarize `SOURCE_URLS.tsv`, `SHA256SUMS`, and active `.partial` transfers recursively under `dir_path`."""
    if not dir_path.is_dir():
        raise LedgerError(f"directory does not exist: {dir_path}")

    rel = (rel_dir or dir_path.name).strip("/")
    merged_files: dict[str, dict[str, Any]] = {}
    raw_rows_total = 0
    streaming = False
    subdirs: list[dict[str, Any]] = []

    for cur_root, dirnames, filenames in os.walk(dir_path):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith(".fetch_to_tower"))
        p = Path(cur_root)
        sub_rel = p.relative_to(dir_path).as_posix()
        has_urls = "SOURCE_URLS.tsv" in filenames
        has_sums = "SHA256SUMS" in filenames
        has_partial = any(f.endswith(".partial") for f in filenames)
        if has_partial:
            streaming = True

        if sub_rel == ".":
            if has_urls:
                by_f, raw_n = _parse_source_urls_file(p / "SOURCE_URLS.tsv")
                raw_rows_total += raw_n
                for fn, info in by_f.items():
                    merged_files[fn] = info
        else:
            if has_urls or has_sums or has_partial:
                sub_files: dict[str, dict[str, Any]] = {}
                sub_raw = 0
                if has_urls:
                    sub_files, sub_raw = _parse_source_urls_file(p / "SOURCE_URLS.tsv")
                    raw_rows_total += sub_raw
                    for fn, info in sub_files.items():
                        merged_files[f"{sub_rel}/{fn}"] = info
                sub_sums_sha = None
                sub_sums_cnt = 0
                if has_sums:
                    sub_sums_sha, _, sub_entries = _parse_sha256sums_file(p / "SHA256SUMS")
                    sub_sums_cnt = len(sub_entries)
                sub_dates = sorted(v["date_utc"] for v in sub_files.values() if v["date_utc"])
                sub_dmin = sub_dates[0][:10] if sub_dates else None
                sub_dmax = sub_dates[-1][:10] if sub_dates else None
                subdirs.append(
                    {
                        "rel_dir": f"{rel}/{sub_rel}",
                        "files": len(sub_files),
                        "raw_rows": sub_raw,
                        "bytes": sum(v["bytes"] for v in sub_files.values()),
                        "date_min": sub_dmin,
                        "date_max": sub_dmax,
                        "date_range": _format_date_range(sub_dmin, sub_dmax),
                        "sha256sums_sha256": sub_sums_sha,
                        "sha256sums_count": sub_sums_cnt,
                        "streaming": has_partial,
                    }
                )

    top_sums_sha = None
    top_sums_bytes = 0
    top_sums_count = 0
    if (dir_path / "SHA256SUMS").is_file():
        top_sums_sha, top_sums_bytes, top_entries = _parse_sha256sums_file(dir_path / "SHA256SUMS")
        top_sums_count = len(top_entries)

    dates = sorted(v["date_utc"] for v in merged_files.values() if v["date_utc"])
    ts_min = dates[0] if dates else None
    ts_max = dates[-1] if dates else None
    date_min = ts_min[:10] if ts_min else None
    date_max = ts_max[:10] if ts_max else None

    single_file_sha256 = None
    single_filename = None
    if len(merged_files) == 1:
        only = next(iter(merged_files.values()))
        single_file_sha256 = only["sha256"]
        single_filename = only["filename"]

    return {
        "root_label": root_label,
        "rel_dir": rel,
        "stored_path": f"tower:/mnt/user/research-bulk/bio-box/{root_label}/{rel}/",
        "files": len(merged_files),
        "raw_rows": raw_rows_total,
        "bytes": sum(v["bytes"] for v in merged_files.values()),
        "ts_min": ts_min,
        "ts_max": ts_max,
        "date_min": date_min,
        "date_max": date_max,
        "date_range": _format_date_range(date_min, date_max),
        "sha256sums_sha256": top_sums_sha if not subdirs else None,
        "sha256sums_bytes": top_sums_bytes if not subdirs else 0,
        "sha256sums_count": top_sums_count + sum(c.get("sha256sums_count", 0) for c in subdirs),
        "top_sha256sums_sha256": top_sums_sha,
        "single_filename": single_filename,
        "single_file_sha256": single_file_sha256,
        "streaming": streaming,
        "subdirs": subdirs,
    }


def _rollup_parent_dirs(dirs_map: dict[str, dict[str, Any]]) -> None:
    """Roll up disjoint single-directory records into ancestor directory summaries (including parent+child trees)."""
    orig_dirs = dict(dirs_map)
    by_parent: dict[str, list[dict[str, Any]]] = {}
    for _key, info in orig_dirs.items():
        root_label = info["root_label"]
        rel = info["rel_dir"]
        parts = rel.split("/")
        for depth in range(1, len(parts)):
            prefix = "/".join(parts[:depth])
            pkey = f"{root_label}/{prefix}"
            by_parent.setdefault(pkey, []).append(info)

    for pkey, children in sorted(by_parent.items()):
        root_label, prefix = pkey.split("/", 1)
        self_info = orig_dirs.get(pkey)
        members = ([self_info] if self_info is not None else []) + children
        files_tot = sum(m["files"] for m in members)
        bytes_tot = sum(m["bytes"] for m in members)
        raw_tot = sum(m.get("raw_rows", m["files"]) for m in members)
        ts_mins = [m["ts_min"] for m in members if m.get("ts_min")]
        ts_maxs = [m["ts_max"] for m in members if m.get("ts_max")]
        ts_min = min(ts_mins) if ts_mins else None
        ts_max = max(ts_maxs) if ts_maxs else None
        dmin = ts_min[:10] if ts_min else None
        dmax = ts_max[:10] if ts_max else None

        active_children = sorted(
            [c for c in children if c["files"] > 0 or c.get("streaming") or c.get("sha256sums_count", 0) > 0],
            key=lambda x: x["rel_dir"],
        )
        top_sums_sha = self_info.get("sha256sums_sha256") if self_info is not None else None
        # Only expose single-directory `sha256sums_sha256` when no child directory carries files of its own
        only_top_has_files = (
            self_info is not None
            and self_info["files"] > 0
            and sum(c["files"] for c in children) == 0
        )

        dirs_map[pkey] = {
            "root_label": root_label,
            "rel_dir": prefix,
            "stored_path": f"tower:/mnt/user/research-bulk/bio-box/{root_label}/{prefix}/",
            "files": files_tot,
            "raw_rows": raw_tot,
            "bytes": bytes_tot,
            "ts_min": ts_min,
            "ts_max": ts_max,
            "date_min": dmin,
            "date_max": dmax,
            "date_range": _format_date_range(dmin, dmax),
            "sha256sums_sha256": top_sums_sha if only_top_has_files else None,
            "sha256sums_bytes": self_info["sha256sums_bytes"] if only_top_has_files and self_info else 0,
            "sha256sums_count": sum(m.get("sha256sums_count", 0) for m in members),
            "top_sha256sums_sha256": top_sums_sha,
            "single_filename": self_info["single_filename"] if only_top_has_files and self_info else None,
            "single_file_sha256": self_info["single_file_sha256"] if only_top_has_files and self_info else None,
            "streaming": any(m.get("streaming", False) for m in members),
            "subdirs": active_children,
        }


def build_inventory_from_roots(
    *,
    sources_incoming: Path | None = None,
    sources_incoming_s2: Path | None = None,
    sources_staged: Path | None = None,
) -> dict[str, Any]:
    """Scan local or synthetic `sources_incoming*` directory trees into an inventory dictionary."""
    dirs_map: dict[str, dict[str, Any]] = {}
    for root_label, root_path in (
        ("sources_incoming", sources_incoming),
        ("sources_incoming_s2", sources_incoming_s2),
        ("sources", sources_staged),
    ):
        if root_path is None or not root_path.is_dir():
            continue
        for cur_root, dirnames, filenames in os.walk(root_path):
            dirnames[:] = sorted(d for d in dirnames if not d.startswith(".fetch_to_tower"))
            p = Path(cur_root)
            if p == root_path:
                continue
            has_urls = "SOURCE_URLS.tsv" in filenames
            has_sums = "SHA256SUMS" in filenames
            has_partial = any(f.endswith(".partial") for f in filenames)
            if has_urls or has_sums or has_partial:
                rel = p.relative_to(root_path).as_posix()
                dirs_map[f"{root_label}/{rel}"] = _summarize_single_dir(p, root_label=root_label, rel_dir=rel)

    _rollup_parent_dirs(dirs_map)
    return {"collected_at_utc": None, "dirs": dirs_map}


def load_inventory_tsv(tsv_path: Path) -> dict[str, Any]:
    """Load a TSV inventory generated by `r26_collect_tower.sh`."""
    collected_at = None
    dirs_map: dict[str, dict[str, Any]] = {}
    with tsv_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            s = line.rstrip("\r\n")
            if not s:
                continue
            if s.startswith("# collected_at_utc="):
                collected_at = s.split("=", 1)[1].strip()
                continue
            if s.startswith("#") or s.startswith("tree\t"):
                continue
            parts = s.split("\t")
            if len(parts) < 14:
                continue
            (
                tree,
                rel,
                _lock_present,
                has_part,
                u_files,
                u_raw,
                u_bytes,
                ts_min,
                ts_max,
                sm_cnt,
                sm_bytes,
                sm_sha,
                s_file,
                s_sha,
            ) = parts[:14]
            ts_min_v = None if ts_min == "-" else ts_min
            ts_max_v = None if ts_max == "-" else ts_max
            dmin = ts_min_v[:10] if ts_min_v else None
            dmax = ts_max_v[:10] if ts_max_v else None
            sums_sha_v = None if sm_sha == "-" else sm_sha
            dirs_map[f"{tree}/{rel}"] = {
                "root_label": tree,
                "rel_dir": rel,
                "stored_path": f"tower:/mnt/user/research-bulk/bio-box/{tree}/{rel}/",
                "files": int(u_files),
                "raw_rows": int(u_raw),
                "bytes": int(u_bytes),
                "ts_min": ts_min_v,
                "ts_max": ts_max_v,
                "date_min": dmin,
                "date_max": dmax,
                "date_range": _format_date_range(dmin, dmax),
                "sha256sums_sha256": sums_sha_v,
                "sha256sums_bytes": int(sm_bytes),
                "sha256sums_count": int(sm_cnt),
                "top_sha256sums_sha256": sums_sha_v,
                "single_filename": None if s_file == "-" else s_file,
                "single_file_sha256": None if s_sha == "-" else s_sha,
                "streaming": has_part == "yes",
                "subdirs": [],
            }
    _rollup_parent_dirs(dirs_map)
    return {"collected_at_utc": collected_at, "dirs": dirs_map}


def extract_url_forms(sources_text: str) -> list[tuple[str, str, str]]:
    """Mirror `fetch_to_tower.sh`'s exact URL-form extraction (`(entry_name, key, token)`)."""
    forms: list[tuple[str, str, str]] = []
    cur_name = ""
    for line in sources_text.splitlines():
        m_name = _NAME_LINE_RE.match(line)
        if m_name:
            cur_name = _strip_quotes(m_name.group(1))
            continue
        m_url = _URL_LINE_RE.match(line)
        if m_url and cur_name:
            key = m_url.group(1)
            first_tok = m_url.group(2).split()[0]
            cleaned = _strip_quotes(first_tok)
            if _HTTP_FTP_RE.match(cleaned):
                forms.append((cur_name, key, cleaned))
    return forms


def extract_all_url_lines(sources_text: str) -> list[tuple[str, str]]:
    """Return `(entry_name, exact_line)` for every 4-space-indented `url*:` line in `sources_text`."""
    out: list[tuple[str, str]] = []
    cur_name = ""
    for line in sources_text.splitlines():
        m_name = _NAME_LINE_RE.match(line)
        if m_name:
            cur_name = _strip_quotes(m_name.group(1))
            continue
        if _URL_LINE_RE.match(line) and cur_name:
            out.append((cur_name, line))
    return out


def _expand_tower_targets(stored_val: str) -> list[tuple[str, str, bool]]:
    """Extract `(declared_root, rel_dir, is_dir_hint)` targets referenced in a `stored:` or `fetched:` string."""
    targets: list[tuple[str, str, bool]] = []
    seen: set[tuple[str, str]] = set()
    last_root: str | None = None
    last_top_seg: str | None = None

    for m in _TOWER_PATH_RE.finditer(stored_val):
        root_lbl = m.group(1)
        raw_token = m.group(2)
        is_dir_hint = raw_token.endswith("/")
        raw_rel = raw_token.strip("/")
        if not raw_rel:
            continue
        bmatch = re.match(r"^(.*)\{([^{}]+)\}(.*)$", raw_rel)
        if bmatch:
            pre, inner, post = bmatch.groups()
            rels = [f"{pre}{item.strip()}{post}".strip("/") for item in inner.split(",") if item.strip()]
        else:
            rels = [raw_rel]
        for r in rels:
            if "<" in r or ">" in r:
                continue
            last_root = root_lbl
            if "/" in r:
                last_top_seg = r.split("/", 1)[0]
            pair = (root_lbl, r)
            if pair not in seen:
                seen.add(pair)
                targets.append((root_lbl, r, is_dir_hint))

    if last_root is not None and last_top_seg is not None:
        for em in _ELLIPSIS_PATH_RE.finditer(stored_val):
            raw_sub = em.group(1)
            is_dir_hint = raw_sub.endswith("/")
            sub_rel = raw_sub.strip("/")
            if not sub_rel or "<" in sub_rel or ">" in sub_rel:
                continue
            full_rel = f"{last_top_seg}/{sub_rel}"
            pair = (last_root, full_rel)
            if pair not in seen:
                seen.add(pair)
                targets.append((last_root, full_rel, is_dir_hint))

    return targets


def _resolve_target_in_inventory(
    declared_root: str,
    rel_dir: str,
    dirs_map: dict[str, dict[str, Any]],
    *,
    is_dir_hint: bool = True,
) -> dict[str, Any] | None:
    """Resolve `(declared_root, rel_dir)` across `sources_incoming` and `sources_incoming_s2`."""
    candidates = [
        f"{declared_root}/{rel_dir}",
        f"sources_incoming_s2/{rel_dir}",
        f"sources_incoming/{rel_dir}",
    ]
    for ckey in candidates:
        info = dirs_map.get(ckey)
        if info and (info.get("files", 0) > 0 or info.get("streaming")):
            return info
    # Only fall back to parent directory when the reference clearly named a file (not a trailing-`/` directory)
    last_seg = rel_dir.rsplit("/", 1)[-1]
    if not is_dir_hint and "/" in rel_dir and "." in last_seg:
        parent_rel = rel_dir.rsplit("/", 1)[0]
        for root_try in (declared_root, "sources_incoming_s2", "sources_incoming"):
            info = dirs_map.get(f"{root_try}/{parent_rel}")
            if info and (info.get("files", 0) > 0 or info.get("streaming")):
                return info
    return None


def _replace_top_scalar_field(block_lines: list[str], field: str, new_value: str) -> list[str]:
    """Replace a 4-space-indented scalar field (`    <field>: ...` plus >=6-space continuations) in an entry block."""
    prefix = f"    {field}:"
    out: list[str] = []
    i = 0
    replaced = False
    safe_val = new_value
    if ": " in safe_val or safe_val.startswith(("{", "[", '"', "'", "*", "&", "!", "|", ">", "%", "@", "`")):
        escaped = safe_val.replace("\\", "\\\\").replace('"', '\\"')
        safe_val = f'"{escaped}"'

    while i < len(block_lines):
        line = block_lines[i]
        if not replaced and (line.startswith(prefix + " ") or line == prefix):
            out.append(f"    {field}: {safe_val}")
            replaced = True
            i += 1
            while i < len(block_lines):
                nxt = block_lines[i]
                if nxt.startswith("      ") and not nxt.lstrip().startswith("#"):
                    i += 1
                else:
                    break
            continue
        out.append(line)
        i += 1
    if not replaced:
        insert_at = len(out)
        while insert_at > 0 and out[insert_at - 1].strip() == "":
            insert_at -= 1
        out.insert(insert_at, f"    {field}: {safe_val}")
    return out


def _describe_sha256sums(info: dict[str, Any]) -> str:
    """Format SHA256SUMS digest provenance for a single or rolled-up directory."""
    n = info["files"]
    if n == 1 and info.get("single_file_sha256") and info.get("sha256sums_sha256"):
        return (
            f"file sha256 {info['single_file_sha256']}, "
            f"SHA256SUMS sha256 {info['sha256sums_sha256']}"
        )
    if info.get("sha256sums_sha256"):
        return f"SHA256SUMS sha256 {info['sha256sums_sha256']}"

    subdirs = info.get("subdirs", [])
    top_sha = info.get("top_sha256sums_sha256")
    if len(subdirs) <= 4 and (top_sha or any(s.get("sha256sums_sha256") for s in subdirs)):
        parts: list[str] = []
        if top_sha:
            parts.append(f"top SHA256SUMS sha256 {top_sha}")
        base_prefix = info["rel_dir"] + "/"
        for s in subdirs:
            if s.get("sha256sums_sha256"):
                short_rel = s["rel_dir"][len(base_prefix) :] if s["rel_dir"].startswith(base_prefix) else s["rel_dir"]
                parts.append(f"{short_rel} SHA256SUMS sha256 {s['sha256sums_sha256']}")
        if parts:
            return "; ".join(parts)

    return f"{len(subdirs)} subdirs SHA256SUMS"


def _format_reconciled_fields(
    entry_name: str,
    resolved: list[dict[str, Any]],
    *,
    as_of: str,
    expected_count: int | None,
) -> tuple[str, str]:
    if len(resolved) == 1:
        info = resolved[0]
        n = info["files"]
        b = info["bytes"]
        dr = info["date_range"]
        fword = "file" if n == 1 else "files"
        is_streaming = bool(info.get("streaming")) or (expected_count is not None and n < expected_count)
        stored = info["stored_path"]
        if is_streaming:
            cnt_txt = f"{n:,} of {expected_count:,} {fword}" if expected_count is not None else f"{n:,} {fword}"
            fetched = (
                f"in progress as of {as_of} ({cnt_txt}, {b:,} B landed {dr}, "
                f"*measured* SOURCE_URLS.tsv on tower; active streaming under DEC-81 (5) / DEC-89)"
            )
        else:
            sha_txt = _describe_sha256sums(info)
            fetched = (
                f"{dr} ({n:,} {fword}, {b:,} B, {sha_txt}, "
                f"*measured* SOURCE_URLS.tsv + SHA256SUMS on tower as of {as_of})"
            )
        return fetched, stored

    # Multi-directory entry (e.g. `tower_r26_11_6_ncbi_fna` or `tower_batch6b_ncbi_o4_org`)
    n_tot = sum(x["files"] for x in resolved)
    b_tot = sum(x["bytes"] for x in resolved)
    dmins = [x["date_min"] for x in resolved if x.get("date_min")]
    dmaxs = [x["date_max"] for x in resolved if x.get("date_max")]
    dr_tot = _format_date_range(min(dmins) if dmins else None, max(dmaxs) if dmaxs else None)
    fword = "file" if n_tot == 1 else "files"
    is_streaming = any(x.get("streaming") for x in resolved) or (
        expected_count is not None and n_tot < expected_count
    )
    parts_desc = []
    for x in resolved:
        sub_sha = _describe_sha256sums(x)
        parts_desc.append(f"{x['rel_dir']} {x['files']:,} files / {x['bytes']:,} B ({sub_sha})")
    joined_desc = "; ".join(parts_desc)
    stored = "; ".join(x["stored_path"] for x in resolved)
    if is_streaming:
        cnt_txt = f"{n_tot:,} of {expected_count:,} {fword}" if expected_count is not None else f"{n_tot:,} {fword}"
        fetched = (
            f"in progress as of {as_of} ({cnt_txt}, {b_tot:,} B landed {dr_tot} across {joined_desc}, "
            f"*measured* SOURCE_URLS.tsv on tower; active streaming under DEC-81 (5) / DEC-89)"
        )
    else:
        fetched = (
            f"{dr_tot} ({n_tot:,} {fword}, {b_tot:,} B across {joined_desc}, "
            f"*measured* SOURCE_URLS.tsv + SHA256SUMS on tower as of {as_of})"
        )
    return fetched, stored


def reconcile_sources_yaml(
    sources_text: str,
    inventory: dict[str, Any],
    *,
    as_of: str = "2026-11-28",
    expected_counts: dict[str, int] | None = None,
    target_overrides: dict[str, list[tuple[str, str]]] | None = None,
    refresh_names: set[str] | None = None,
) -> tuple[str, list[str]]:
    """Surgically update stale `fetched:` and `stored:` fields in `sources_text` from `inventory`."""
    dirs_map = inventory.get("dirs", {})
    lines = sources_text.splitlines()
    trailing_nl = sources_text.endswith("\n")

    entry_starts = [idx for idx, ln in enumerate(lines) if _NAME_LINE_RE.match(ln)]
    if not entry_starts:
        raise LedgerError("no `  - name:` entries found in sources.yaml")

    preamble = lines[: entry_starts[0]]
    blocks: list[list[str]] = []
    for pos_idx, start in enumerate(entry_starts):
        end = entry_starts[pos_idx + 1] if pos_idx + 1 < len(entry_starts) else len(lines)
        blocks.append(lines[start:end])

    orig_data = yaml.safe_load(sources_text)
    if not isinstance(orig_data, dict) or not isinstance(orig_data.get("sources"), list):
        raise LedgerError("sources.yaml must parse to a mapping with a `sources` list")
    orig_by_name = {e["name"]: e for e in orig_data["sources"] if isinstance(e, dict) and "name" in e}

    updated_names: list[str] = []
    new_blocks: list[list[str]] = []

    for blk in blocks:
        m_name = _NAME_LINE_RE.match(blk[0])
        assert m_name is not None
        name = _strip_quotes(m_name.group(1))
        entry = orig_by_name.get(name, {})
        fetched_val = entry.get("fetched")
        stored_val = entry.get("stored")
        fetched_str = fetched_val if isinstance(fetched_val, str) else ""
        stored_str = stored_val if isinstance(stored_val, str) else ""

        is_stale = (
            fetched_str.strip().lower().startswith("nothing yet")
            or stored_str.strip().lower().startswith("will be tower:")
            or (refresh_names is not None and name in refresh_names)
        )
        if not is_stale:
            new_blocks.append(blk)
            continue

        if target_overrides and name in target_overrides:
            targets = [(r, d, True) for r, d in target_overrides[name]]
        else:
            targets = _expand_tower_targets(stored_str) or _expand_tower_targets(fetched_str)
        if not targets:
            new_blocks.append(blk)
            continue

        resolved: list[dict[str, Any]] = []
        seen_paths: set[str] = set()
        for decl_root, rel_dir, is_dir_hint in targets:
            info = _resolve_target_in_inventory(decl_root, rel_dir, dirs_map, is_dir_hint=is_dir_hint)
            if info is not None and info.get("files", 0) > 0 and info["stored_path"] not in seen_paths:
                seen_paths.add(info["stored_path"])
                resolved.append(info)

        if not resolved:
            new_blocks.append(blk)
            continue

        new_fetched, new_stored = _format_reconciled_fields(
            name,
            resolved,
            as_of=as_of,
            expected_count=(expected_counts or {}).get(name),
        )
        blk2 = _replace_top_scalar_field(blk, "fetched", new_fetched)
        should_update_stored = (
            stored_val is None
            or (isinstance(stored_val, str) and not stored_val.strip())
            or stored_str.strip().lower().startswith("will be ")
        )
        if should_update_stored:
            blk2 = _replace_top_scalar_field(blk2, "stored", new_stored)
        new_blocks.append(blk2)
        updated_names.append(name)

    rebuilt_lines = list(preamble)
    for blk in new_blocks:
        rebuilt_lines.extend(blk)
    updated_text = "\n".join(rebuilt_lines) + ("\n" if trailing_nl else "")

    if extract_all_url_lines(updated_text) != extract_all_url_lines(sources_text):
        raise LedgerError("internal error: reconcile_sources_yaml altered a `url*:` line")
    new_data = yaml.safe_load(updated_text)
    if not isinstance(new_data, dict) or not isinstance(new_data.get("sources"), list):
        raise LedgerError("reconciled sources.yaml did not parse cleanly")
    new_by_name = {e["name"]: e for e in new_data["sources"]}
    for name, orig_e in orig_by_name.items():
        if "approved" in orig_e and new_by_name[name].get("approved") != orig_e["approved"]:
            raise LedgerError(f"internal error: reconcile_sources_yaml altered `approved:` on {name}")

    return updated_text, updated_names


def _is_valid_redistributable_value(val: Any) -> bool:
    if isinstance(val, bool):
        return True
    if isinstance(val, str):
        s = val.strip().lower()
        return s.startswith("true ") or s.startswith("false ") or s.startswith("unknown (")
    return False


def verify_sources_yaml(
    sources_path: Path,
    *,
    baseline_sources_path: Path | None = None,
    inventory: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate `sources.yaml` schema, URL host invariants, baseline URL preservation, and tower freshness."""
    if not sources_path.is_file():
        raise LedgerError(f"missing sources.yaml: {sources_path}")
    text = sources_path.read_text(encoding="utf-8")
    try:
        parsed = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise LedgerError(f"invalid YAML in {sources_path}: {exc}") from exc

    if not isinstance(parsed, dict) or parsed.get("version") != 1:
        raise LedgerError(f"{sources_path}: top-level `version: 1` mapping required")
    entries = parsed.get("sources")
    if not isinstance(entries, list) or not entries:
        raise LedgerError(f"{sources_path}: non-empty `sources` list required")

    seen_names: set[str] = set()
    for idx, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise LedgerError(f"{sources_path}: entry #{idx} is not a mapping")
        name = entry.get("name")
        if not isinstance(name, str) or not name.strip():
            raise LedgerError(f"{sources_path}: entry #{idx} lacks non-empty string `name`")
        if name in seen_names:
            raise LedgerError(f"{sources_path}: duplicate entry name `{name}`")
        seen_names.add(name)

        if "redistributable" in entry:
            redist = entry["redistributable"]
            if not _is_valid_redistributable_value(redist):
                raise LedgerError(
                    f"{sources_path}: `{name}` has invalid `redistributable` value (got {redist!r})"
                )
        if name.startswith(("biosino_", "jgi_portal_")):
            if not isinstance(entry.get("redistributable"), bool):
                raise LedgerError(
                    f"{sources_path}: portal entry `{name}` must have boolean `redistributable`"
                )

        url_keys = [k for k in entry if k.startswith("url")]
        if not url_keys:
            raise LedgerError(f"{sources_path}: `{name}` lacks any `url*` field")
        for uk in url_keys:
            uval = entry[uk]
            if not isinstance(uval, str) or not uval.strip():
                raise LedgerError(f"{sources_path}: `{name}.{uk}` must be a non-empty string")
            if name in _LEGACY_PROSE_URL_ENTRIES and uval.strip().lower().startswith("not recorded"):
                continue
            tok = _strip_quotes(uval.split()[0])
            m_sch = _SCHEME_RE.match(tok)
            if not m_sch:
                raise LedgerError(f"{sources_path}: `{name}.{uk}` has invalid URL scheme/host: {tok!r}")
            host = m_sch.group(2)
            if "{" in host or "}" in host:
                raise LedgerError(f"{sources_path}: `{name}.{uk}` has placeholder in host segment: {tok!r}")

        lic = entry.get("license")
        if not isinstance(lic, str) or not lic.strip():
            raise LedgerError(f"{sources_path}: `{name}` lacks non-empty `license`")

        has_prov = any(
            (isinstance(entry.get(k), str) and bool(entry[k].strip()))
            or (isinstance(entry.get(k), list) and len(entry[k]) > 0)
            for k in ("fetched", "stored", "files")
        )
        if not has_prov:
            raise LedgerError(f"{sources_path}: `{name}` lacks `fetched`, `stored`, or `files`")

        if name.startswith(("tower_", "biosino_", "jgi_portal_")):
            for req_field in ("fetched", "stored"):
                fval = entry.get(req_field)
                if not isinstance(fval, str) or not fval.strip():
                    raise LedgerError(f"{sources_path}: `{name}` lacks non-empty `{req_field}`")

    cur_forms = set(extract_url_forms(text))
    if baseline_sources_path is not None:
        if not baseline_sources_path.is_file():
            raise LedgerError(f"missing baseline sources.yaml: {baseline_sources_path}")
        base_forms = set(extract_url_forms(baseline_sources_path.read_text(encoding="utf-8")))
        missing_forms = sorted(base_forms - cur_forms)
        if missing_forms:
            sample = ", ".join(f"{n}.{k}={u}" for n, k, u in missing_forms[:5])
            raise LedgerError(
                f"{sources_path}: {len(missing_forms)} baseline URL form(s) missing or altered ({sample})"
            )

    if inventory is not None:
        dirs_map = inventory.get("dirs", {})
        stale_hits: list[str] = []
        for entry in entries:
            name = entry["name"]
            fetched_val = entry.get("fetched")
            stored_val = entry.get("stored")
            fetched_str = fetched_val.strip().lower() if isinstance(fetched_val, str) else ""
            stored_str = stored_val.strip() if isinstance(stored_val, str) else ""
            if not (fetched_str.startswith("nothing yet") or stored_str.lower().startswith("will be tower:")):
                continue
            targets = _expand_tower_targets(stored_str) or _expand_tower_targets(
                fetched_val if isinstance(fetched_val, str) else ""
            )
            for decl_root, rel_dir, is_dir_hint in targets:
                info = _resolve_target_in_inventory(decl_root, rel_dir, dirs_map, is_dir_hint=is_dir_hint)
                if info is not None and info.get("files", 0) > 0:
                    stale_hits.append(f"{name} ({info['root_label']}/{info['rel_dir']}: {info['files']} files)")
                    break
        if stale_hits:
            raise LedgerError(
                f"{sources_path}: stale `nothing yet` / `will be tower:` on landed entries: {', '.join(stale_hits)}"
            )

    return {
        "entry_count": len(entries),
        "url_form_count": len(cur_forms),
    }


def _split_markdown_sections(md_text: str) -> dict[str, str]:
    sections: dict[str, list[str]] = {}
    cur_h2 = ""
    cur_h3 = ""
    for line in md_text.splitlines():
        if line.startswith("## "):
            cur_h2 = line.strip()
            cur_h3 = ""
            sections.setdefault(cur_h2, [])
            continue
        if line.startswith("### "):
            cur_h3 = line.strip()
            sections.setdefault(cur_h3, [])
        if cur_h2:
            sections[cur_h2].append(line)
        if cur_h3:
            sections[cur_h3].append(line)
    return {k: "\n".join(v) for k, v in sections.items()}


def _find_section(sections: dict[str, str], prefix: str) -> str:
    for heading, body in sections.items():
        if heading.startswith(prefix):
            return body
    raise LedgerError(f"missing Markdown section starting with `{prefix}`")


def _discover_conf_stages(stages_dir: Path) -> set[str]:
    if not stages_dir.is_dir():
        raise LedgerError(f"missing stages directory: {stages_dir}")
    found: set[str] = set()
    search_dirs = [stages_dir]
    tpl = stages_dir / "templates"
    if tpl.is_dir():
        search_dirs.append(tpl)
    for d in search_dirs:
        for conf in sorted(d.glob("*.conf")):
            stem = conf.stem
            found.add(stem)
            for ln in conf.read_text(encoding="utf-8").splitlines():
                if ln.startswith("STAGE="):
                    declared = ln.split("=", 1)[1].strip()
                    if declared != stem:
                        raise LedgerError(f"{conf}: STAGE={declared!r} does not match filename stem {stem!r}")
    if not found:
        raise LedgerError(f"no .conf files found under {stages_dir}")
    return found


def _extract_sec3_stage_tokens(sec3: str) -> set[str]:
    """Extract stage tokens from the `Stage` column of §3's Markdown table (including `` `C1` to `C6` `` ranges)."""
    table_lines = [ln.strip() for ln in sec3.splitlines() if ln.lstrip().startswith("|")]
    if not table_lines:
        return set()

    stage_col_idx: int | None = None
    stages: set[str] = set()
    for ln in table_lines:
        cells = [c.strip() for c in ln.strip("|").split("|")]
        if all(set(c) <= {"-", ":", " "} for c in cells if c):
            continue
        if stage_col_idx is None:
            for idx, c in enumerate(cells):
                if c.strip().lower() == "stage":
                    stage_col_idx = idx
                    break
            if stage_col_idx is not None:
                continue
            stage_col_idx = 0

        target_cell = cells[stage_col_idx] if stage_col_idx < len(cells) else ln
        for m_rng in _STAGE_RANGE_RE.finditer(target_cell):
            prefix, start_s, end_s = m_rng.groups()
            lo, hi = int(start_s), int(end_s)
            if lo <= hi and (hi - lo) <= 200:
                for k in range(lo, hi + 1):
                    stages.add(f"{prefix}{k}")
        for m_tok in _BACKTICK_STAGE_RE.finditer(target_cell):
            stages.add(m_tok.group(1))
    return stages


def _discover_staging_remainders(staging_root: Path | None) -> set[str]:
    out: set[str] = set()
    if staging_root is None or not staging_root.is_dir():
        return out
    for child in sorted(staging_root.iterdir()):
        if not child.is_dir():
            continue
        for f in sorted(child.iterdir()):
            if not f.is_file():
                continue
            low = f.name.lower()
            if ("left_out" in low or "left_for_later" in low) and not low.endswith(".gz"):
                non_empty = [ln for ln in f.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
                if len(non_empty) >= 2:
                    out.add(child.name)
            elif ("left_out" in low or "left_for_later" in low) and low.endswith(".gz") and f.stat().st_size > 64:
                out.add(child.name)
    return out


def verify_catalog_goal(
    catalog_goal_path: Path,
    stages_dir: Path,
    *,
    staging_root: Path | None = None,
    required_remainders: list[str] | tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """Verify `CATALOG_GOAL.md` §3 covers every `.conf` stage and §4.1 registers every required remainder."""
    if not catalog_goal_path.is_file():
        raise LedgerError(f"missing CATALOG_GOAL.md: {catalog_goal_path}")
    md_text = catalog_goal_path.read_text(encoding="utf-8")
    sections = _split_markdown_sections(md_text)

    sec2 = _find_section(sections, "## 2.")
    sec3 = _find_section(sections, "## 3.")
    sec41 = _find_section(sections, "### 4.1")
    _find_section(sections, "### 4.3")

    for dec_tag in ("DEC-88", "DEC-89"):
        if dec_tag not in sec2 and dec_tag not in md_text:
            raise LedgerError(f"{catalog_goal_path}: §2 lacks `{dec_tag}`")

    conf_stages = _discover_conf_stages(stages_dir)
    sec3_stages = _extract_sec3_stage_tokens(sec3)

    missing_stages = sorted(conf_stages - sec3_stages)
    if missing_stages:
        raise LedgerError(
            f"{catalog_goal_path}: §3 is missing stage row(s) for .conf file(s): {', '.join(missing_stages)}"
        )

    req_rem = set(DEFAULT_REQUIRED_REMAINDERS if required_remainders is None else required_remainders)
    req_rem |= _discover_staging_remainders(staging_root)
    missing_rem = sorted(tok for tok in req_rem if tok not in sec41)
    if missing_rem:
        raise LedgerError(
            f"{catalog_goal_path}: §4.1 is missing deferred/left-out remainder registration for: "
            f"{', '.join(missing_rem)}"
        )

    return {
        "conf_stage_count": len(conf_stages),
        "sec3_stage_count": len(sec3_stages),
        "checked_remainders": sorted(req_rem),
    }


def verify_worker_rules(worker_rules_path: Path) -> None:
    """Verify `STAGING_WORKER_RULES.md` requires ledger synchronization checked via `r26_sync_sources_ledger.py verify`."""
    if not worker_rules_path.is_file():
        raise LedgerError(f"missing STAGING_WORKER_RULES.md: {worker_rules_path}")
    txt = worker_rules_path.read_text(encoding="utf-8")
    missing = [
        needle
        for needle in ("sources.yaml", "CATALOG_GOAL.md", "r26_sync_sources_ledger.py verify")
        if needle not in txt
    ]
    if missing:
        raise LedgerError(
            f"{worker_rules_path}: missing required synchronization rule token(s): {', '.join(missing)}"
        )


def _parse_expected_counts(raw_items: list[str] | None) -> dict[str, int]:
    out: dict[str, int] = {}
    for item in raw_items or []:
        if "=" not in item:
            raise LedgerError(f"--expected-count must have form NAME=INT (got {item!r})")
        k, v = item.split("=", 1)
        out[k.strip()] = int(v.strip())
    return out


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Synchronize and verify bio-box source & catalog ledgers (Card R-26).")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_ver = sub.add_parser("verify", help="Verify sources.yaml, CATALOG_GOAL.md, stage .confs, and worker rules.")
    p_ver.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    p_ver.add_argument("--baseline-sources", type=Path, default=None)
    p_ver.add_argument("--inventory-tsv", type=Path, default=None)
    p_ver.add_argument("--catalog-goal", type=Path, default=DEFAULT_CATALOG_GOAL)
    p_ver.add_argument("--stages-dir", type=Path, default=DEFAULT_STAGES_DIR)
    p_ver.add_argument("--staging-root", type=Path, default=None)
    p_ver.add_argument("--worker-rules", type=Path, default=DEFAULT_WORKER_RULES)
    p_ver.add_argument(
        "--require-remainder",
        action="append",
        dest="require_remainders",
        default=None,
        help="Stage/dataset token that must appear in CATALOG_GOAL.md §4.1 (may be repeated).",
    )

    p_rec = sub.add_parser("reconcile", help="Reconcile stale tower entries in sources.yaml from an inventory TSV.")
    p_rec.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    p_rec.add_argument("--inventory-tsv", type=Path, required=True)
    p_rec.add_argument("--as-of", default="2026-11-28")
    p_rec.add_argument("--output", type=Path, default=None, help="Output path (defaults to --sources in place).")
    p_rec.add_argument(
        "--expected-count",
        action="append",
        dest="expected_counts",
        default=None,
        help="Expected job count for an actively streaming entry as NAME=INT (may be repeated).",
    )
    p_rec.add_argument(
        "--refresh-name",
        action="append",
        dest="refresh_names",
        default=None,
        help="Entry name to refresh even if previously partially fetched (may be repeated).",
    )
    p_rec.add_argument(
        "--report-json",
        type=Path,
        default=None,
        help="Optional path to write a JSON reconciliation report.",
    )

    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        if args.cmd == "verify":
            inv = load_inventory_tsv(args.inventory_tsv) if args.inventory_tsv else None
            src_stats = verify_sources_yaml(
                args.sources,
                baseline_sources_path=args.baseline_sources,
                inventory=inv,
            )
            cg_stats = verify_catalog_goal(
                args.catalog_goal,
                args.stages_dir,
                staging_root=args.staging_root,
                required_remainders=args.require_remainders,
            )
            verify_worker_rules(args.worker_rules)
            print(
                f"OK: sources={src_stats['entry_count']} url_forms={src_stats['url_form_count']} "
                f"conf_stages={cg_stats['conf_stage_count']} remainders_checked={len(cg_stats['checked_remainders'])}"
            )
            return 0

        if args.cmd == "reconcile":
            inv = load_inventory_tsv(args.inventory_tsv)
            orig = args.sources.read_text(encoding="utf-8")
            exp_map = _parse_expected_counts(args.expected_counts)
            ref_set = set(args.refresh_names) if args.refresh_names else None
            updated, touched = reconcile_sources_yaml(
                orig,
                inv,
                as_of=args.as_of,
                expected_counts=exp_map,
                refresh_names=ref_set,
            )
            out_path = args.output or args.sources
            out_path.write_text(updated, encoding="utf-8")
            if args.report_json is not None:
                args.report_json.write_text(
                    json.dumps(
                        {
                            "collected_at_utc": inv.get("collected_at_utc"),
                            "as_of": args.as_of,
                            "updated_count": len(touched),
                            "updated_names": touched,
                        },
                        indent=2,
                    )
                    + "\n",
                    encoding="utf-8",
                )
            print(f"reconciled {len(touched)} entries into {out_path}")
            return 0
    except LedgerError as exc:
        print(f"r26_sync_sources_ledger: FAIL: {exc}", file=sys.stderr)
        return 1
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

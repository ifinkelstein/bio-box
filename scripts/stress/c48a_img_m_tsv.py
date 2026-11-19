#!/usr/bin/env python3
"""Card C-48a: convert on-disk IMG/M jqdt JSON exports into deterministic UTF-8 .tsv.gz
tables (zero PII/email columns, mtime=0) and bundle GOLD + IMG/M metadata into one
deterministic uncompressed .tar (>= 64 MiB, excluding .cookie* and temp files).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import html
import json
import re
import tarfile
from html.parser import HTMLParser
from pathlib import Path
from typing import Iterable, Sequence

DEFAULT_TAR_MTIME = 1794960000  # 2026-11-18T00:00:00Z
MIN_STOR_BYTES = 67_108_864  # 64 MiB (AGENTS.md §7.3)

_PII_PATTERNS = (
    re.compile(r"\be[-_ ]?mail\b", re.I),
    re.compile(r"\bcontact\b", re.I),
    re.compile(r"\bsubmitter\b", re.I),
    re.compile(r"\bprincipal\s+investigator\b", re.I),
    re.compile(r"\bpi[-_ ]?name\b", re.I),
    re.compile(r"\bpi[-_ ]?email\b", re.I),
    re.compile(r"\b(?:tele)?phone\b", re.I),
    re.compile(r"\bfax\b", re.I),
    re.compile(r"\bstreet\b", re.I),
    re.compile(r"\bpostal\b", re.I),
)

_MICROBIOME_ONLY_COLS = ("Predicted Viruses", "Predicted Plasmids")
_WS_RE = re.compile(r"\s+")


class _CellHTMLParser(HTMLParser):
    """Extract clean text and first anchor href from an HTML table cell."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.parts: list[str] = []
        self.href: str | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() == "a" and self.href is None:
            for k, v in attrs:
                if k.lower() == "href" and v:
                    self.href = v.strip()
                    break

    def handle_data(self, data: str) -> None:
        self.parts.append(data)


def parse_cell(val: object) -> tuple[str, str | None]:
    """Parse one IMG/M jqdt table cell into (clean_display_text, optional_href_or_extra)."""
    if val is None:
        return "", None
    if isinstance(val, bool):
        return ("true" if val else "false"), None
    if isinstance(val, (int, float)):
        return str(val), None

    raw = str(val)
    extra: str | None = None
    if "#==#" in raw:
        parts = raw.split("#==#")
        display = parts[0]
        # IMG/M encodes either "display#==#sort_or_oid" or "display#==##==#link_target"
        tail = [p.strip() for p in parts[1:] if p.strip()]
        if tail:
            extra = tail[-1]
    else:
        display = raw

    if "<" in display and ">" in display:
        parser = _CellHTMLParser()
        parser.feed(display)
        parser.close()
        text = "".join(parser.parts)
        if parser.href and extra is None:
            extra = parser.href
    else:
        text = html.unescape(display)

    cleaned = _WS_RE.sub(" ", text).strip()
    return cleaned, extra


def clean_cell(val: object) -> str:
    """Return clean single-line TSV cell text with HTML and #==# markup stripped."""
    return parse_cell(val)[0]


def assert_no_pii_columns(headers: Sequence[str]) -> None:
    """Raise ValueError if any column header looks like PII or contact/email metadata."""
    bad = [h for h in headers if any(pat.search(h) for pat in _PII_PATTERNS)]
    if bad:
        raise ValueError(f"Refusing PII/contact columns (AGENTS §8): {bad}")


def resolve_headers(base_headers: Sequence[str], raw_ncols: int) -> list[str]:
    """Align IMG/M header list to either 72-col (73 raw) or 70-col (71 raw) exports."""
    if len(base_headers) == raw_ncols:
        return list(base_headers)
    if len(base_headers) - len(_MICROBIOME_ONLY_COLS) == raw_ncols:
        reduced = [h for h in base_headers if h not in _MICROBIOME_ONLY_COLS]
        if len(reduced) == raw_ncols:
            return reduced
    raise ValueError(f"Cannot align {len(base_headers)} base headers to {raw_ncols} data columns")


def is_excluded_bundle_path(path: Path) -> bool:
    """True for cookies, secrets, temp/partial files, or logs that must never be bundled."""
    parts = [p.lower() for p in path.parts]
    name = path.name.lower()
    if any("cookie" in p for p in parts):
        return True
    if name.startswith(".tmp.") or name.endswith(".part") or name.endswith(".tmp"):
        return True
    if name.endswith(".headers.txt") or name.endswith(".sh") or name == "dl.log":
        return True
    return False


def _sort_key_for_row(row: list[str], id_idx: int | None) -> tuple[int, int, str, tuple[str, ...]]:
    if id_idx is not None and id_idx < len(row):
        val = row[id_idx]
        if val.isdigit():
            return (0, int(val), val, tuple(row))
        return (1, 0, val, tuple(row))
    return (2, 0, "", tuple(row))


def convert_img_m_json(
    json_path: Path,
    out_tsv_gz: Path,
    base_headers: Sequence[str],
) -> dict[str, object]:
    """Convert one IMG/M jqdt JSON file into a deterministic UTF-8 .tsv.gz table."""
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    if "data" in payload:
        raw_rows = payload["data"]
    elif "aaData" in payload:
        raw_rows = payload["aaData"]
    else:
        raise ValueError(f"{json_path}: expected 'data' or 'aaData' key, got {sorted(payload)}")

    if not raw_rows:
        headers = list(base_headers)
    else:
        headers = resolve_headers(base_headers, len(raw_rows[0]))

    drop_first = bool(headers) and headers[0] == ""
    out_headers = headers[1:] if drop_first else headers
    assert_no_pii_columns(out_headers)

    id_idx = out_headers.index("IMG Genome ID") if "IMG Genome ID" in out_headers else None

    cleaned_rows: list[list[str]] = []
    expected_len = len(headers)
    for idx, raw_row in enumerate(raw_rows):
        if len(raw_row) != expected_len:
            raise ValueError(f"{json_path}: row {idx} has {len(raw_row)} cols, expected {expected_len}")
        cells = raw_row[1:] if drop_first else raw_row
        cleaned_rows.append([clean_cell(c) for c in cells])

    cleaned_rows.sort(key=lambda r: _sort_key_for_row(r, id_idx))

    out_tsv_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp_out = out_tsv_gz.with_name(f".{out_tsv_gz.name}.tmp")
    with open(tmp_out, "wb") as raw_fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0, compresslevel=6) as gz:
            gz.write(("\t".join(out_headers) + "\n").encode("utf-8"))
            for row in cleaned_rows:
                gz.write(("\t".join(row) + "\n").encode("utf-8"))
    tmp_out.replace(out_tsv_gz)

    h = hashlib.sha256()
    with open(out_tsv_gz, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)

    unique_oids = (
        len({r[id_idx] for r in cleaned_rows if r[id_idx]}) if id_idx is not None else len(cleaned_rows)
    )
    return {
        "input": str(json_path),
        "output": str(out_tsv_gz),
        "rows": len(cleaned_rows),
        "cols": len(out_headers),
        "unique_oids": unique_oids,
        "bytes": out_tsv_gz.stat().st_size,
        "sha256": h.hexdigest(),
    }


def build_metadata_tar(
    src_dirs: Iterable[Path],
    out_tar: Path,
    *,
    mtime: int = DEFAULT_TAR_MTIME,
    min_bytes: int = MIN_STOR_BYTES,
) -> dict[str, object]:
    """Bundle metadata directories into one deterministic uncompressed .tar."""
    entries: list[tuple[str, Path]] = []
    for d in sorted((Path(p) for p in src_dirs), key=lambda p: p.name):
        if not d.is_dir():
            raise FileNotFoundError(f"missing source directory: {d}")
        for p in sorted(d.rglob("*")):
            if not p.is_file():
                continue
            rel = p.relative_to(d)
            if is_excluded_bundle_path(rel):
                continue
            arcname = f"{d.name}/{rel.as_posix()}"
            entries.append((arcname, p))

    entries.sort(key=lambda item: item[0])
    out_tar.parent.mkdir(parents=True, exist_ok=True)
    tmp_tar = out_tar.with_name(f".{out_tar.name}.tmp")

    with tarfile.open(tmp_tar, "w:", format=tarfile.GNU_FORMAT) as tf:
        for arcname, fpath in entries:
            ti = tf.gettarinfo(str(fpath), arcname=arcname)
            ti.uid = 0
            ti.gid = 0
            ti.uname = ""
            ti.gname = ""
            ti.mtime = mtime
            ti.mode = 0o644
            with open(fpath, "rb") as fh:
                tf.addfile(ti, fh)

    size = tmp_tar.stat().st_size
    if size < min_bytes:
        tmp_tar.unlink(missing_ok=True)
        raise ValueError(f"{out_tar} is {size} bytes (< {min_bytes} bytes floor)")

    tmp_tar.replace(out_tar)
    h = hashlib.sha256()
    with open(out_tar, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)

    return {
        "output": str(out_tar),
        "members": len(entries),
        "bytes": size,
        "sha256": h.hexdigest(),
    }


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_conv = sub.add_parser("convert", help="Convert all IMG/M *_raw.json files in a directory to .tsv.gz")
    p_conv.add_argument("--img-m-dir", type=Path, required=True)
    p_conv.add_argument("--cols-json", type=Path, required=True)

    p_tar = sub.add_parser("bundle", help="Bundle gold/ and img_m/ into deterministic uncompressed .tar")
    p_tar.add_argument("--src-dir", type=Path, action="append", required=True)
    p_tar.add_argument("--out-tar", type=Path, required=True)
    p_tar.add_argument("--min-bytes", type=int, default=MIN_STOR_BYTES)

    args = ap.parse_args(argv)
    if args.cmd == "convert":
        base_headers = json.loads(args.cols_json.read_text(encoding="utf-8"))
        results = []
        for raw_json in sorted(args.img_m_dir.glob("*_raw.json")):
            stem = raw_json.name.replace("_raw.json", "")
            out_path = args.img_m_dir / f"{stem}.tsv.gz"
            results.append(convert_img_m_json(raw_json, out_path, base_headers))
        print(json.dumps(results, indent=2))
        return 0

    if args.cmd == "bundle":
        res = build_metadata_tar(args.src_dir, args.out_tar, min_bytes=args.min_bytes)
        print(json.dumps(res, indent=2))
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())

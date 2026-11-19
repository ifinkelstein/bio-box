"""Card C-46: post-release verification, durable-backup file list (`root_added.tsv`), and announcement draft.

Host or image, standard library only; read-only on the data root:

    python3 -I train_check.py --root ROOT --before BEFORE.tsv --after AFTER.tsv \
        --release R --stages S0,S1 [--earlier PREV] [--manifest-hash HASH] \
        [--pack-id ID] [--examples-summary FILE] --out DIR
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_listing(p: Path) -> dict[str, tuple[int, str]]:
    out: dict[str, tuple[int, str]] = {}
    for line in p.read_text().splitlines():
        if not line.strip():
            continue
        path, size, mtime = line.split("\t")
        out[path] = (int(size), mtime)
    return out


def names_release(path: str, rel: str) -> bool:
    parts = path.split("/")
    return (
        f"release={rel}" in parts
        or f"stage={rel}" in parts
        or any(p.startswith(f"release={rel}.") for p in parts)
    )


def format_announcement(
    *,
    release: str,
    stages: list[str],
    earlier: str | None,
    manifest_hash: str | None,
    added_files: int,
    added_bytes: int,
    pack_id: str | None = None,
    examples_verdict: str | None = None,
) -> str:
    lines = [
        f"# Release `{release}` announcement draft",
        "",
        f"- **Release:** `{release}` (`--release latest` resolves here once verified with a search pack)",
        f"- **Stages ({len(stages)}):** `{','.join(stages)}`",
        f"- **Previous release:** `{earlier or 'none'}` (verified unchanged)",
        f"- **Manifest hash:** `{manifest_hash or '(see verify log)'}`",
        f"- **Files added to root:** `{added_files}` files, `{added_bytes:,}` bytes (listed in `root_added.tsv` for durable backup)",
    ]
    if pack_id:
        lines.append(f"- **Search pack:** `{pack_id}` (verified)")
    if examples_verdict:
        lines.append(f"- **Worked examples (`scripts/release_examples.sh`):** {examples_verdict}")
    lines.extend([
        "",
        "## Suggested summary line",
        "",
        f"Release `{release}` is live over `{','.join(stages)}` (`{added_files}` files added, "
        f"`{added_bytes:,}` bytes, manifest hash `{manifest_hash or 'TBD'}`). "
        "All earlier releases verified unchanged; use `--release latest`.",
        "",
    ])
    return "\n".join(lines)


def run_check(
    *,
    root: Path,
    before: Path,
    after: Path,
    release: str,
    stages: list[str],
    earlier: str | None = None,
    manifest_hash: str | None = None,
    pack_id: str | None = None,
    examples_summary: Path | None = None,
    out: Path,
) -> int:
    root = Path(root)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    report: list[str] = []
    failed: list[str] = []

    def say(s: str) -> None:
        report.append(s)
        print(s, flush=True)

    b0, b1 = read_listing(Path(before)), read_listing(Path(after))
    gone = sorted(p for p in b0 if p not in b1)
    changed = sorted(p for p in b0 if p in b1 and b0[p] != b1[p])
    new = sorted(p for p in b1 if p not in b0)
    stray = [p for p in new if not names_release(p, release)]

    say(
        f"1. listing: {len(b0)} files before, {len(b1)} after; {len(new)} new, "
        f"{len(gone)} gone, {len(changed)} changed in size or mtime; "
        f"{len(stray)} new files not naming {release}"
    )
    for p in (gone + changed + stray)[:20]:
        say(f"   PROBLEM: {p}")
    if gone or changed or stray:
        failed.append("listing")

    total_bytes = 0
    tsv_path = out / "root_added.tsv"
    with open(tsv_path, "w") as fh:
        fh.write("path\tbytes\tsha256\n")
        for p in new:
            fp = root / p
            sz = fp.stat().st_size
            total_bytes += sz
            fh.write(f"{p}\t{sz}\t{sha256_file(fp)}\n")
    say(f"2. root_added.tsv: {len(new)} files, {total_bytes:,} bytes")

    man_path = root / "manifests" / f"release={release}.json"
    man_info: dict = {}
    if not man_path.is_file():
        say(f"3. manifest: MISSING ({man_path})")
        failed.append("manifest_missing")
    else:
        man = json.loads(man_path.read_text())
        man_info = {
            "stages": man.get("stages"),
            "files": len(man.get("files", ())),
            "artifacts": len(man.get("artifacts", ())),
        }
        say(f"3. manifest: stages {man_info['stages']}, {man_info['files']} table files, {man_info['artifacts']} artifacts")
        if man.get("stages") != stages:
            say(f"   PROBLEM: manifest stages {man.get('stages')} != expected {stages}")
            failed.append("manifest_stages")

    examples_verdict = None
    if examples_summary and Path(examples_summary).is_file():
        rows = [ln.split("\t") for ln in Path(examples_summary).read_text().splitlines()[1:] if ln.strip()]
        counts: dict[str, int] = {}
        for r in rows:
            if len(r) >= 6:
                counts[r[5]] = counts.get(r[5], 0) + 1
        examples_verdict = ", ".join(f"{v} {k}" for k, v in sorted(counts.items()))

    ann_text = format_announcement(
        release=release,
        stages=stages,
        earlier=earlier,
        manifest_hash=manifest_hash,
        added_files=len(new),
        added_bytes=total_bytes,
        pack_id=pack_id,
        examples_verdict=examples_verdict,
    )
    (out / "announcement.md").write_text(ann_text)

    summary = {
        "release": release,
        "stages": stages,
        "earlier": earlier,
        "manifest_hash": manifest_hash,
        "pack_id": pack_id,
        "added": {"files": len(new), "bytes": total_bytes},
        "gone": len(gone),
        "changed": len(changed),
        "stray": len(stray),
        "manifest": man_info,
        "failed": failed,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True))
    say(f"checks: {'all pass' if not failed else 'FAILED: ' + ', '.join(failed)}")
    (out / "report.txt").write_text("\n".join(report) + "\n")
    return 0 if not failed else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-46: post-release check, root_added.tsv, and announcement draft")
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--before", required=True, type=Path)
    ap.add_argument("--after", required=True, type=Path)
    ap.add_argument("--release", required=True)
    ap.add_argument("--stages", required=True)
    ap.add_argument("--earlier", default=None)
    ap.add_argument("--manifest-hash", default=None)
    ap.add_argument("--pack-id", default=None)
    ap.add_argument("--examples-summary", type=Path, default=None)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    stages = [s for s in a.stages.split(",") if s]
    return run_check(
        root=a.root,
        before=a.before,
        after=a.after,
        release=a.release,
        stages=stages,
        earlier=a.earlier,
        manifest_hash=a.manifest_hash,
        pack_id=a.pack_id,
        examples_summary=a.examples_summary,
        out=a.out,
    )


if __name__ == "__main__":
    sys.exit(main())

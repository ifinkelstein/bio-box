"""Card C-40: one FASTA from a PlasmidScope archive, which holds one FASTA file per plasmid.

    python -I c40_concat.py SRC ARCHIVE OUT_DIR

ARCHIVE is `<SRC>.fasta.tar.gz` as fetched (untrusted data): it is read as a stream (`tarfile` mode `r|gz`) and never
unpacked, so no path in it is created on disk. Only regular members named `[<TOP>/]fasta/<stem>.<fa|fasta|fna>` are
taken, with one TOP for the whole archive (it is SRC, except EMBL for the ENA archive and none for mMGE's), in the
archive's order; their bytes are appended to OUT_DIR/<SRC>.fasta unchanged, with a newline added when a member
does not end with one. OUT_DIR must be new or empty. Every other member is counted and named in the report. For each
member taken: its name, size, the records in it (lines starting with '>') and its first header's first token, written
to OUT_DIR/members.tsv; a member with other than one record, or whose first token is not the file's stem, is counted.
Prints a JSON report and writes it to OUT_DIR/concat.json. Nothing in the archive is executed.
"""

from __future__ import annotations

import json
import re
import sys
import tarfile
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    src, archive, out = argv[0], Path(argv[1]), Path(argv[2])
    if out.exists() and any(out.iterdir()):
        print(f"c40_concat: {out} is not empty", file=sys.stderr)
        return 2
    out.mkdir(parents=True, exist_ok=True)
    name_re = re.compile(r"^(?:([^/]+)/)?fasta/([^/]+)\.(fa|fasta|fna)$")
    tops: set[str] = set()
    rep = {"source": src, "archive": str(archive), "members": 0, "taken": 0, "bytes": 0, "records": 0,
           "not_taken": [], "not_one_record": 0, "stem_not_first_token": 0, "newline_added": 0}
    tmp = out / f".{src}.fasta.tmp"
    with tarfile.open(archive, mode="r|gz") as tf, open(tmp, "wb") as w, open(out / "members.tsv", "w") as mt:
        mt.write("member\tbytes\trecords\tfirst_token\n")
        for m in tf:
            rep["members"] += 1
            hit = name_re.match(m.name) if m.isfile() else None
            if hit is None:
                if not m.isdir():
                    rep["not_taken"].append(m.name)
                continue
            tops.add(hit.group(1) or "")
            data = tf.extractfile(m).read()
            n_rec = sum(1 for ln in data.splitlines() if ln.startswith(b">"))
            head = data.split(b"\n", 1)[0]
            tok = head[1:].split(None, 1)[0].decode(errors="replace") if head.startswith(b">") and head[1:].split() else ""
            if not data.endswith(b"\n"):
                data += b"\n"
                rep["newline_added"] += 1
            w.write(data)
            rep["taken"] += 1
            rep["bytes"] += len(data)
            rep["records"] += n_rec
            rep["not_one_record"] += n_rec != 1
            rep["stem_not_first_token"] += tok != hit.group(2)
            mt.write(f"{m.name}\t{m.size}\t{n_rec}\t{tok}\n")
    rep["top_directories"] = sorted(tops)
    if len(tops) > 1:
        print(f"c40_concat: members under several top directories {sorted(tops)}", file=sys.stderr)
        return 1
    tmp.rename(out / f"{src}.fasta")
    rep["not_taken"] = rep["not_taken"][:50] + ([f"... {len(rep['not_taken']) - 50} more"]
                                               if len(rep["not_taken"]) > 50 else [])
    (out / "concat.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

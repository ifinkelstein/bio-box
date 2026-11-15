"""Card C-34: scratch roots of symlinks over the data root, read-only on it (the pattern of card C-33's c33_view.py).

    python3 -I c34_view.py rehearsal SRC DST S1,P1,P2
    python3 -I c34_view.py verify SRC REH DST
    python3 -I c34_view.py seven SRC DST

rehearsal  What `bio-box stage labels` reads, for each stage named: a link to SRC's
           l0seq/contig_occurrence/stage=<S> and, inside a real directory ledger/stage=<S>, a link to its step=ingest.
           The step's own files (l0seq/contig_label/stage=<S>, ledger/stage=<S>/step=labels) and run records (runs/)
           land in DST's real directories, so nothing is written through a link.
verify     A view of the whole root with the rehearsal's label files in place, for `release verify` by the code before
           and after the card: DST, DST/l0seq and DST/l0seq/contig_label are real directories; every other top-level
           entry of SRC, every other entry of SRC/l0seq and every stage directory of SRC/l0seq/contig_label is a link;
           each label stage directory of REH (the rehearsal root) is a real directory of links to its files. Python's
           rglob does not descend into a linked directory, so `verify_manifest` sees as unlisted files exactly the
           rehearsal's label files (and nothing under the links), which is the question asked of it.
seven      A view for a release defined afterwards: every top-level entry of SRC but manifests, catalog and release_def
           is a link; release_def is a real, empty directory, so `bio-box release define` writes there and `bio-box
           releases` lists that release alone.
Each refuses when DST exists.
"""
import os
import sys

mode = sys.argv[1]


def link(src: str, dst: str) -> None:
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    os.symlink(src, dst)


if mode == "rehearsal":
    src, dst, stages = sys.argv[2], sys.argv[3], sys.argv[4].split(",")
    os.makedirs(dst)
    for s in stages:
        occ = os.path.join(src, "l0seq", "contig_occurrence", f"stage={s}")
        ing = os.path.join(src, "ledger", f"stage={s}", "step=ingest")
        lab = os.path.join(src, "l0seq", "contig_label", f"stage={s}")
        if not os.path.isdir(occ) or not os.path.isdir(ing):
            sys.exit(f"stage {s}: {occ} or {ing} is missing")
        if os.path.exists(lab):
            sys.exit(f"stage {s}: {lab} exists: the stage has label files")
        link(occ, os.path.join(dst, "l0seq", "contig_occurrence", f"stage={s}"))
        link(ing, os.path.join(dst, "ledger", f"stage={s}", "step=ingest"))
    for d in ("runs", os.path.join("l0seq", "contig_label")):
        os.makedirs(os.path.join(dst, d), exist_ok=True)
    print(f"rehearsal view {dst}: contig_occurrence and the ingest ledger of {stages} linked from {src}")
elif mode == "verify":
    src, reh, dst = sys.argv[2], sys.argv[3], sys.argv[4]
    os.makedirs(os.path.join(dst, "l0seq", "contig_label"))
    n = 0
    for e in sorted(os.listdir(src)):
        if e != "l0seq":
            link(os.path.join(src, e), os.path.join(dst, e))
            n += 1
    for e in sorted(os.listdir(os.path.join(src, "l0seq"))):
        if e != "contig_label":
            link(os.path.join(src, "l0seq", e), os.path.join(dst, "l0seq", e))
            n += 1
    for e in sorted(os.listdir(os.path.join(src, "l0seq", "contig_label"))):
        link(os.path.join(src, "l0seq", "contig_label", e), os.path.join(dst, "l0seq", "contig_label", e))
        n += 1
    added = []
    base = os.path.join(reh, "l0seq", "contig_label")
    for e in sorted(os.listdir(base)):
        if os.path.exists(os.path.join(dst, "l0seq", "contig_label", e)):
            sys.exit(f"{e}: the root has label files of that stage already")
        os.makedirs(os.path.join(dst, "l0seq", "contig_label", e))
        for f in sorted(os.listdir(os.path.join(base, e))):
            link(os.path.join(base, e, f), os.path.join(dst, "l0seq", "contig_label", e, f))
            added.append(f"l0seq/contig_label/{e}/{f}")
    print(f"verify view {dst}: {n} links into {src}; the rehearsal's label files in place: {added}")
elif mode == "seven":
    src, dst = sys.argv[2], sys.argv[3]
    os.makedirs(os.path.join(dst, "release_def"))
    for e in sorted(os.listdir(src)):
        if e not in ("manifests", "catalog", "release_def"):
            link(os.path.join(src, e), os.path.join(dst, e))
    print(f"seven-stage view {dst}: links into {src}, its own empty release_def")
else:
    sys.exit("usage: c34_view.py rehearsal|verify|seven ...")

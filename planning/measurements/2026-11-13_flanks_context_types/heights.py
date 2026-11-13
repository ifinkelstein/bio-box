"""Card W-31: the merge distances of a `nbr flanks` tree and where the largest gaps lie.

    python -I heights.py TREE.nwk [TREE.nwk ...]

Node heights in tree.nwk are merge distance / 2, so a merge's distance is twice the height of its node (the sum of the
branch lengths from it down to any leaf). Prints, per tree: the number of merges, how many are at 0, below 0.5, in
[0.5, 1) and at 1, and the three largest gaps between consecutive distinct distances below 1 and 1 itself.
"""

import sys


def merges(text: str) -> list[float]:
    out, stack = [], []
    i = 0
    while i < len(text):
        c = text[i]
        if c == "(":
            stack.append([])
            i += 1
        elif c == ")":
            kids = stack.pop()
            j = i + 1
            bl = 0.0
            if j < len(text) and text[j] == ":":
                k = j + 1
                while k < len(text) and text[k] not in ",);":
                    k += 1
                bl = float(text[j + 1:k])
                i = k
            else:
                i = j
            h = max(kids)  # height of this node (all children equal in an ultrametric tree)
            out.append(2 * h)
            if stack:
                stack[-1].append(h + bl)
        elif c == ",":
            i += 1
        elif c == ";":
            break
        else:  # a leaf name and its branch length
            k = i
            while text[k] != ":":
                k += 1
            m = k + 1
            while text[m] not in ",);":
                m += 1
            stack[-1].append(float(text[k + 1:m]))
            i = m
    return out


for path in sys.argv[1:]:
    hs = sorted(merges(open(path).read().strip()))
    below = [h for h in hs if h < 1 - 1e-6]
    levels = []
    for h in below:
        if not levels or h > levels[-1] + 1e-6:
            levels.append(h)
    levels.append(1.0)
    gaps = sorted(((levels[i + 1] - levels[i], levels[i], levels[i + 1]) for i in range(len(levels) - 1)),
                  reverse=True)[:3]
    print(f"{path}\tmerges {len(hs)}\tat_0 {sum(h < 1e-6 for h in hs)}\tbelow_0.5 {sum(h < 0.5 for h in hs)}\t"
          f"0.5_to_1 {sum(0.5 <= h < 1 - 1e-6 for h in hs)}\tat_1 {sum(h >= 1 - 1e-6 for h in hs)}\tgaps "
          + " ".join(f"{g:.3f}[{a:.3f},{b:.3f}]" for g, a, b in gaps))

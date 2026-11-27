#!/usr/bin/env python3
"""The release-diversity figure: its data, its checks and its rendering (card A-31).

One file, four commands:

    plot_release_diversity.py root   --root ROOT --out AUDIT.json [--releases ...] [--built ...]
        Read-only, on the pod in the image: for every release manifest named, the stored bases of each stage
        (`contig_index` default variants, from the manifest's own file list), the `cluster_member` rows, and,
        where the release's `cluster_member` files are still in the root, the NR90 and NR50 clusters
        (rows whose member is its cluster's representative). `--built` adds stages built and not yet in a
        release (their partition directory is listed: no manifest names them yet). Writes the audit JSON.
    plot_release_diversity.py check  [--data JSON] [--sources sources.yaml] [--ownership A29.json]
        Checks the data file against the tree (every number has a label and a source; quotes are found in the
        files they cite; derived numbers equal their formula; sums by class equal the totals; each stage is in
        one class and one release) and prints the table the figure draws.
    plot_release_diversity.py tower  [--data JSON] [--sources sources.yaml]
        Prints what `sources.yaml` records as landed on `tower` and assigned to no release.
    plot_release_diversity.py render --out PREFIX [--data JSON] [--sources sources.yaml] [--ownership A29.json]
        Writes PREFIX.png (300 DPI) and PREFIX.svg (no date, fixed hash salt).

The figure's numbers are not typed into this file: every stage's bases, every release's proteins and clusters and
every benchmark are in `planning/measurements/release_diversity.json`, each with its label and source. When card
A-29's ownership JSON (`a29_r26_11_5_chain.py ownership --out`) is in the tree, set `r26_11_5_ownership` in the data
file to its path (or pass `--ownership`): its `stages.<STAGE>.new_bases` replace the staging notes' new bases of the
`r26.11.5` stages, and a stage it holds that was drawn as an estimate (card C-83's) is drawn as staged.
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import operator
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA = REPO_ROOT / "planning" / "measurements" / "release_diversity.json"
DEFAULT_SOURCES = REPO_ROOT / "sources.yaml"

LABELS = ("measured", "derived", "estimate")
STATES = ("released", "built", "staged", "planned")
STATE_TEXT = {
    "released": "In a release (measured)",
    "built": "Built, not yet in a release (measured)",
    "staged": "Staged, not built (measured bases)",
    "planned": "Planned (estimate)",
}


# ---------------------------------------------------------------------------------------------------------------
# Numbers: each is {"value": ..., "label": ..., "source": ...} with an optional quote (a string found in the file
# at `path`) or, for a derived number, a formula over named inputs.
# ---------------------------------------------------------------------------------------------------------------

_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}


def eval_formula(formula: str, inputs: dict[str, float]) -> float:
    """Evaluate `formula` (names, numbers, + - * / and parentheses only) over `inputs`."""

    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.BinOp) and type(node.op) in _OPS:
            return _OPS[type(node.op)](ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -ev(node.operand)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return node.value
        if isinstance(node, ast.Name):
            return inputs[node.id]
        raise ValueError(f"formula {formula!r}: {ast.dump(node)} is not allowed")

    return ev(ast.parse(formula, mode="eval"))


def quote_number(quote: str) -> float:
    """The number a quote states: digits with thousands commas, optional decimals (`1,664,109,938`, `41.24`)."""
    m = re.search(r"-?\d[\d,]*(?:\.\d+)?", quote)
    if not m:
        raise ValueError(f"no number in quote {quote!r}")
    return float(m.group(0).replace(",", ""))


def number_problems(num: object, where: str, root: Path = REPO_ROOT) -> list[str]:
    """What is wrong with one number record (empty when nothing is)."""
    if not isinstance(num, dict):
        return [f"{where}: not a number record"]
    out = []
    if num.get("label") not in LABELS:
        out.append(f"{where}: label {num.get('label')!r} is not one of {LABELS}")
    if not str(num.get("source", "")).strip():
        out.append(f"{where}: no source")
    if not isinstance(num.get("value"), (int, float)) or isinstance(num.get("value"), bool):
        out.append(f"{where}: no numeric value")
        return out
    value = float(num["value"])
    if "quote" in num:
        path = root / num["path"] if num.get("path") else None
        if path is None or not path.is_file():
            out.append(f"{where}: quote without a file in the tree ({num.get('path')!r})")
        elif num["quote"] not in path.read_text(encoding="utf-8"):
            out.append(f"{where}: {num['quote']!r} is not in {num['path']}")
        else:
            q = quote_number(num["quote"]) * float(num.get("quote_scale", 1))
            if not math.isclose(q, value, rel_tol=float(num.get("quote_rel_tol", 0)), abs_tol=1e-9):
                out.append(f"{where}: value {value} differs from the quote {num['quote']!r} (x{num.get('quote_scale', 1)})")
    if "formula" in num:
        ins = {}
        for k, sub in num.get("inputs", {}).items():
            out += number_problems(sub, f"{where}.inputs.{k}", root)
            if isinstance(sub, dict) and isinstance(sub.get("value"), (int, float)):
                ins[k] = float(sub["value"])
        for key, fkey in (("value", "formula"), ("low", "low_formula"), ("high", "high_formula")):
            if fkey not in num:
                continue
            try:
                got = _rounded(eval_formula(num[fkey], ins), num.get("round", 0))
                if not math.isclose(got, float(num.get(key, math.nan)), rel_tol=1e-12, abs_tol=1e-9):
                    out.append(f"{where}: {fkey} {num[fkey]} gives {got}, not {num.get(key)}")
            except (KeyError, ValueError, ZeroDivisionError) as e:
                out.append(f"{where}: {fkey} {num.get(fkey)!r} cannot be evaluated ({e})")
    if num.get("label") == "derived" and not ({"formula", "quote"} & set(num)):
        out.append(f"{where}: derived without a formula or the quote of a note that derives it")
    if num.get("label") == "measured" and not ({"quote", "host", "citation"} & set(num)):
        out.append(f"{where}: measured without a quote from a file, the host, date and command, or a citation")
    if "citation" in num and num.get("re_verified") is not False:
        out.append(f"{where}: an outside figure must say it was not re-verified (re_verified: false)")
    for k in ("low", "high"):
        if k in num and not isinstance(num[k], (int, float)):
            out.append(f"{where}: {k} is not a number")
    if "low" in num and "high" in num and not (num["low"] <= value <= num["high"]):
        out.append(f"{where}: value {value} outside its range [{num['low']}, {num['high']}]")
    return out


def _rounded(v: float, digits: int) -> float:
    return float(int(round(v))) if int(digits) == 0 else round(v, int(digits))


def fill(obj: object) -> int:
    """Set every formula's value (and low, high) from its inputs, inputs first; returns how many were set."""
    n = 0
    if isinstance(obj, dict):
        for v in obj.values():
            n += fill(v)
        if "formula" in obj and "value" in obj:
            ins = {k: float(s["value"]) for k, s in obj.get("inputs", {}).items()}
            for key, fkey in (("value", "formula"), ("low", "low_formula"), ("high", "high_formula")):
                if fkey in obj:
                    r = _rounded(eval_formula(obj[fkey], ins), obj.get("round", 0))
                    obj[key] = int(r) if int(obj.get("round", 0)) == 0 else r
                    n += 1
    elif isinstance(obj, list):
        for v in obj:
            n += fill(v)
    return n


def walk_numbers(obj: object, where: str = "") -> list[tuple[str, dict]]:
    """Every number record in a JSON tree: a dict with `value` and `label`."""
    found = []
    if isinstance(obj, dict):
        if "value" in obj and "label" in obj:
            found.append((where, obj))
            return found
        for k, v in obj.items():
            found += walk_numbers(v, f"{where}.{k}" if where else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            key = v.get("stage") or v.get("id") if isinstance(v, dict) else None
            found += walk_numbers(v, f"{where}[{key or i}]")
    return found


# ---------------------------------------------------------------------------------------------------------------
# The table the figure draws.
# ---------------------------------------------------------------------------------------------------------------

def load(path: Path = DEFAULT_DATA) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def apply_ownership(data: dict, ownership: dict | None, path: str = "") -> dict:
    """Card A-29's ownership JSON replaces the new bases of the `r26.11.5` stages it holds (a copy is returned)."""
    data = json.loads(json.dumps(data))
    if not ownership:
        return data
    held = ownership.get("stages", {})
    for st in data["stages"]:
        if st["release"] != data.get("ownership_release", "r26.11.5") or st["stage"] not in held:
            continue
        nb = int(held[st["stage"]]["new_bases"])
        st["bases"] = {"value": nb, "label": "measured", "path": path, "quote": f'"new_bases": {nb}',
                       "source": f"card A-29's exact ownership over the r26.11.5 chain ({path})"}
        if st["state"] == "planned":
            st["state"] = "staged"
    data["r26_11_5_upper_bound"] = False
    return data


def release_ids(data: dict) -> list[str]:
    return [r["id"] for r in data["releases"]]


def table(data: dict) -> dict:
    """Per release column: Gbp by class and state (cumulative), totals, low and high of the planned part."""
    rels = release_ids(data)
    classes = [c["id"] for c in data["classes"]]
    cols = {}
    for i, rid in enumerate(rels):
        upto = set(rels[: i + 1])
        cell = {c: {s: 0.0 for s in STATES} for c in classes}
        low = high = 0.0
        for st in data["stages"]:
            if st["release"] not in upto:
                continue
            b = st["bases"]
            gbp = float(b["value"]) / 1e9
            cell[st["class"]][st["state"]] += gbp
            low += float(b.get("low", b["value"])) / 1e9
            high += float(b.get("high", b["value"])) / 1e9
        total = sum(sum(v.values()) for v in cell.values())
        cols[rid] = {"by_class_state": cell, "total": total, "low": low, "high": high,
                     "by_class": {c: sum(cell[c].values()) for c in classes},
                     "by_state": {s: sum(cell[c][s] for c in classes) for s in STATES}}
    return cols


def effective_classes(gbp: list[float]) -> float:
    """Hill number of order 1, exp(Shannon entropy) of the class shares of stored Gbp."""
    total = sum(gbp)
    if total <= 0:
        return 0.0
    return math.exp(-sum((v / total) * math.log(v / total) for v in gbp if v > 0))


def data_problems(data: dict, root: Path = REPO_ROOT) -> list[str]:
    """Everything `check` asserts about the data file."""
    out = []
    nums = walk_numbers(data)
    for where, num in nums:
        out += number_problems(num, where, root)
    classes = {c["id"] for c in data["classes"]}
    rels = release_ids(data)
    seen = {}
    for st in data["stages"]:
        name = st["stage"]
        if name in seen:
            out.append(f"stage {name}: listed twice (releases {seen[name]} and {st['release']})")
        seen[name] = st["release"]
        if st.get("class") not in classes:
            out.append(f"stage {name}: class {st.get('class')!r} is not a class of the figure")
        if st.get("release") not in rels:
            out.append(f"stage {name}: release {st.get('release')!r} is not a column of the figure")
        if st.get("state") not in STATES:
            out.append(f"stage {name}: state {st.get('state')!r}")
        rel = next((r for r in data["releases"] if r["id"] == st.get("release")), None)
        if rel and rel["status"] == "released" and st.get("state") != "released":
            out.append(f"stage {name}: in released {rel['id']} but drawn as {st.get('state')}")
        if st.get("state") in ("released", "built") and st["bases"].get("label") != "measured":
            out.append(f"stage {name}: {st['state']} but its bases are {st['bases'].get('label')}")
        if st.get("state") == "planned" and st["bases"].get("label") == "measured":
            out.append(f"stage {name}: planned but its bases are measured")
    cols = table(data)
    for r in data["releases"]:
        col = cols[r["id"]]
        if abs(sum(col["by_class"].values()) - col["total"]) > 1e-9:
            out.append(f"{r['id']}: sums by class differ from the total")
        if r["status"] == "released":
            stored = r.get("stored_bases")
            if not stored:
                out.append(f"{r['id']}: released without stored_bases")
            elif abs(col["total"] * 1e9 - float(stored["value"])) > 0.5:
                out.append(f"{r['id']}: stages sum to {col['total'] * 1e9:.0f} bp, the release stores {stored['value']}")
    audit = data.get("audit", {})
    for r in data["releases"]:
        a = audit.get("releases", {}).get(r["id"])
        if r["status"] != "released":
            continue
        if not a:
            out.append(f"{r['id']}: released but not in the audit")
            continue
        if r["stored_bases"]["value"] != a["stored_bases"]:
            out.append(f"{r['id']}: stored_bases {r['stored_bases']['value']} != audit {a['stored_bases']}")
        if r["proteins"]["value"] != a["proteins"]:
            out.append(f"{r['id']}: proteins {r['proteins']['value']} != audit {a['proteins']}")
        for lvl in ("nr50", "nr90"):
            if a.get(lvl) is not None and r[lvl]["value"] != a[lvl]:
                out.append(f"{r['id']}: {lvl} {r[lvl]['value']} != audit {a[lvl]}")
        if sorted(a["stages"]) != sorted(s["stage"] for s in data["stages"] if s["release"] in rels[: rels.index(r["id"]) + 1]):
            out.append(f"{r['id']}: the audit's stages differ from the stages drawn")
    names = {s["stage"] for s in data["stages"]}
    for entry, to in data.get("tower_remainder", {}).get("assigned", {}).items():
        if to not in names and not str(to).startswith("("):
            out.append(f"tower: {entry} is assigned to {to!r}, which is no stage of the figure")
    for st in data["stages"]:
        a = audit.get("stages", {}).get(st["stage"])
        if st["state"] in ("released", "built"):
            if not a:
                out.append(f"stage {st['stage']}: {st['state']} but not in the audit")
            elif st["bases"]["value"] != a["default_bases"]:
                out.append(f"stage {st['stage']}: bases {st['bases']['value']} != audit {a['default_bases']}")
    return out


# ---------------------------------------------------------------------------------------------------------------
# `tower`: what sources.yaml records as landed on tower and assigned to no release.
# ---------------------------------------------------------------------------------------------------------------

_SEQUENCE_NAME = re.compile(r"\.(fna|fa|fasta|fas|fsa|ffn|seq)(\.(gz|bz2|xz|zst|bgz))?$|\.(tar(\.(gz|bz2|xz|zst))?|tgz|zip)$",
                            re.IGNORECASE)


def tower_remainder(data: dict, sources_path: Path = DEFAULT_SOURCES) -> dict:
    """Bytes of the files `sources.yaml` lists with a size for each `tower_*` entry that no release takes."""
    import yaml

    doc = yaml.safe_load(Path(sources_path).read_text(encoding="utf-8"))
    rule = data["tower_remainder"]
    assigned = rule["assigned"]
    rows = []
    for e in doc["sources"]:
        name = e.get("name", "")
        if not name.startswith("tower_"):
            continue
        files = [f for f in (e.get("files") or []) if isinstance(f, dict) and isinstance(f.get("bytes"), int)]
        seq = sum(f["bytes"] for f in files if _SEQUENCE_NAME.search(str(f.get("file", ""))))
        other = sum(f["bytes"] for f in files) - seq
        rows.append({"entry": name, "assigned": assigned.get(name), "files": len(files), "sequence_bytes": seq,
                     "other_bytes": other})
    free = [r for r in rows if r["assigned"] is None and r["files"]]
    return {
        "entries": len(rows),
        "unknown_assigned": sorted(set(assigned) - {r["entry"] for r in rows}),
        "landed_unassigned_entries": len(free),
        "sequence_bytes": sum(r["sequence_bytes"] for r in free),
        "other_bytes": sum(r["other_bytes"] for r in free),
        "not_landed_unassigned_entries": sum(1 for r in rows if r["assigned"] is None and not r["files"]),
        "rows": rows,
    }


# ---------------------------------------------------------------------------------------------------------------
# `root`: the audit's read of the data root (pod, image, read-only).
# ---------------------------------------------------------------------------------------------------------------

def collect_root(root: Path, releases: list[str], built: list[str], threads: int = 2, memory: str = "4GB") -> dict:
    import duckdb
    import time

    con = duckdb.connect()
    con.execute(f"SET threads={int(threads)}")
    con.execute(f"SET memory_limit='{memory}'")
    out: dict = {"root": str(root), "releases": {}, "stages": {}, "seconds": {}}

    def stage_sums(name: str, paths: list[Path]) -> None:
        t = time.monotonic()
        files = [str(p) for p in paths]
        rows, n_def, b_def, b_all = con.execute(
            "SELECT count(*), count(*) FILTER (WHERE is_default_variant), "
            "coalesce(sum(contig_len) FILTER (WHERE is_default_variant), 0), coalesce(sum(contig_len), 0) "
            "FROM read_parquet(?)", [files]).fetchone()
        out["stages"][name] = {"files": sorted(str(Path(p).relative_to(root)) for p in paths), "rows": int(rows),
                               "default_contigs": int(n_def), "default_bases": int(b_def), "all_bases": int(b_all)}
        out["seconds"][f"stage {name}"] = round(time.monotonic() - t, 1)

    for rel in releases:
        m = json.loads((root / "manifests" / f"release={rel}.json").read_text(encoding="utf-8"))
        files = m["files"]
        stages = list(m["stages"])
        ci: dict[str, list[Path]] = {}
        for f in files:
            if f["table"] == "contig_index":
                ci.setdefault(f["partitions"]["stage"], []).append(root / f["path"])
        for st in stages:
            if st not in out["stages"]:
                stage_sums(st, ci[st])
        cm = [f for f in files if f["table"] == "cluster_member"]
        cm_rows = sum(int(f["rows"]) for f in cm)
        rec = {"stages": stages, "stored_bases": sum(out["stages"][s]["default_bases"] for s in stages),
               "default_contigs": sum(out["stages"][s]["default_contigs"] for s in stages),
               "cluster_member_rows": cm_rows, "cluster_member_files": len(cm)}
        present = [root / f["path"] for f in cm if (root / f["path"]).is_file()]
        if cm and len(present) == len(cm):
            t = time.monotonic()
            by = con.execute(
                "SELECT level, count(*), count(*) FILTER (WHERE member_uid = rep_uid) FROM read_parquet(?) "
                "GROUP BY level ORDER BY level", [[str(p) for p in present]]).fetchall()
            out["seconds"][f"clusters {rel}"] = round(time.monotonic() - t, 1)
            rec["levels"] = {lvl: {"members": int(n), "clusters": int(k)} for lvl, n, k in by}
            members = {v["members"] for v in rec["levels"].values()}
            rec["proteins"] = members.pop() if len(members) == 1 else None
            rec["nr50"] = rec["levels"].get("nr50", {}).get("clusters")
            rec["nr90"] = rec["levels"].get("nr90", {}).get("clusters")
        else:
            rec["levels"] = None
            rec["proteins"] = cm_rows // 2 if cm_rows % 2 == 0 else None
            rec["proteins_rule"] = "cluster_member rows of the manifest / 2 (levels nr90 and nr50); files gone"
            rec["nr50"] = rec["nr90"] = None
        out["releases"][rel] = rec
    for st in built:
        part = root / "l0seq" / "contig_index" / f"stage={st}"
        stage_sums(st, sorted(part.glob("*.parquet")))
        out["stages"][st]["listing"] = "partition directory listed (no manifest names a built stage before its release)"
    con.close()
    return out


# ---------------------------------------------------------------------------------------------------------------
# `render`.
# ---------------------------------------------------------------------------------------------------------------

def _style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Source Sans 3", "Source Sans Pro", "DejaVu Sans", "Arial"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#333333",
        "axes.linewidth": 0.8,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "hatch.linewidth": 0.7,
        "svg.fonttype": "none",
        "svg.hashsalt": "bio-box-a31",
        "savefig.dpi": 300,
    })


def _state_style(state: str, color: str) -> dict:
    """Fill or outline for a state; the colour is always the class's."""
    if state == "released":
        return dict(facecolor=color, edgecolor="white", linewidth=0.6, hatch=None)
    if state == "built":
        return dict(facecolor=color, edgecolor="white", linewidth=0.6, hatch="xx", alpha=0.95)
    if state == "staged":
        return dict(facecolor=color, edgecolor="white", linewidth=0.6, hatch="////", alpha=0.55)
    return dict(facecolor="white", edgecolor=color, linewidth=1.1, hatch=None, linestyle=(0, (2, 1.5)))


def build_figure(data: dict, tower: dict | None = None):
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import numpy as np

    classes = data["classes"]
    cids = [c["id"] for c in classes]
    color = {c["id"]: c["color"] for c in classes}
    rels = data["releases"]
    rids = [r["id"] for r in rels]
    cols = table(data)
    n = len(rels)
    x = np.arange(n)
    bar_w = 0.66
    released_idx = [i for i, r in enumerate(rels) if r["status"] == "released"]
    last_rel = released_idx[-1]

    def label(r):
        return f"{r['id']}\n({r['date_label']})"

    xt = [label(r) for r in rels]

    fig = plt.figure(figsize=(16.5, 11.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.0], hspace=0.55, wspace=0.22,
                          left=0.06, right=0.985, top=0.86, bottom=0.13)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # a: share of stored Gbp by class; inside each class its states, released at the bottom.
    for i, rid in enumerate(rids):
        col = cols[rid]
        bottom = 0.0
        for c in cids:
            for s in STATES:
                v = col["by_class_state"][c][s]
                if v <= 0:
                    continue
                pct = 100.0 * v / col["total"]
                ax_a.bar(x[i], pct, bar_w, bottom=bottom, **_state_style(s, color[c]))
                bottom += pct
        bottom = 0.0
        for c in cids:
            pct = 100.0 * col["by_class"][c] / col["total"]
            if pct >= 9.0:
                ax_a.text(x[i], bottom + pct / 2, f"{pct:.0f}%", ha="center", va="center", fontsize=7.6,
                          fontweight="bold", color="#111111",
                          bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.75))
            bottom += pct
        d = effective_classes([col["by_class"][c] for c in cids])
        ax_a.text(x[i], 102.0, f"{d:.1f}", ha="center", va="bottom", fontsize=8.0, color="#333333",
                  fontweight="bold" if i <= last_rel else "normal")
    ax_a.text(-0.45, 110.5, f"D = exp(H) over the {len(cids)} source classes: one class gives 1.0, "
              f"{len(cids)} equal classes give {len(cids)}.0",
              ha="left", va="bottom", fontsize=7.4, color="#555555", style="italic")
    ax_a.set_ylim(0, 117)
    ax_a.set_yticks([0, 25, 50, 75, 100])
    ax_a.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_a.set_ylabel("Share of stored sequence (% of Gbp)", fontsize=9.2)
    ax_a.set_title("a   Source composition by release", loc="left", fontsize=10, fontweight="bold", pad=10)

    # b: cumulative Gbp, stacked by state (released at the bottom), classes inside each state.
    tops = np.zeros(n)
    for i, rid in enumerate(rids):
        col = cols[rid]
        bottom = 0.0
        for s in STATES:
            for c in cids:
                v = col["by_class_state"][c][s]
                if v <= 0:
                    continue
                ax_b.bar(x[i], v, bar_w, bottom=bottom, **_state_style(s, color[c]))
                bottom += v
            if col["by_state"][s] > 0 and bottom < col["total"] - 1e-9:
                ax_b.plot([x[i] - bar_w / 2, x[i] + bar_w / 2], [bottom, bottom], color="#222222", lw=0.8, zorder=4)
        tops[i] = bottom
        if col["high"] - col["low"] > 1e-6:
            ax_b.errorbar(x[i], col["total"], yerr=[[col["total"] - col["low"]], [col["high"] - col["total"]]],
                          fmt="none", ecolor="#333333", elinewidth=0.9, capsize=3, zorder=5)
    refs = data["benchmarks"]["gbp_refs"]
    ymax_b = max(max(cols[r]["high"] for r in rids), *(float(r["value"]) for r in refs)) * 1.14
    for ref in refs:
        rv = float(ref["value"])
        ax_b.axhline(rv, color="#666666", linestyle=":", linewidth=1.1, zorder=2)
        ax_b.text(-0.45, rv + ymax_b * 0.012, ref["figure_label"], ha="left", va="bottom", fontsize=7.6,
                  color="#444444", style="italic")
    for i, rid in enumerate(rids):
        col = cols[rid]
        txt = f"{col['total']:.1f}"
        if rels[i]["status"] != "released" and col["high"] - col["low"] > 1e-6:
            txt = f"{col['total']:.0f}\n({col['low']:.0f} to {col['high']:.0f})"
        y = col["high"] + ymax_b * 0.015
        ax_b.text(x[i], y, txt, ha="center", va="bottom", fontsize=7.6, color="#222222",
                  fontweight="bold" if i == last_rel else "normal", linespacing=1.0, zorder=6,
                  bbox=dict(boxstyle="round,pad=0.1", facecolor="white", edgecolor="none", alpha=0.85))
    ax_b.set_ylim(0, ymax_b)
    ax_b.set_ylabel("Stored sequence (Gbp)", fontsize=9.2)
    ax_b.set_title("b   Cumulative stored sequence by state and source class", loc="left", fontsize=10,
                   fontweight="bold", pad=10)
    if tower:
        gb = tower["sequence_bytes"] / 1e9
        ob = tower["other_bytes"] / 1e9
        ax_b.text(0.015, 0.80, data["tower_remainder"]["figure_text"].format(
            entries=tower["landed_unassigned_entries"], seq_gb=gb, other_gb=ob,
            waiting=tower["not_landed_unassigned_entries"]),
            transform=ax_b.transAxes, ha="left", va="top", fontsize=7.4, color="#333333", linespacing=1.25,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#F6F6F6", edgecolor="#BBBBBB", lw=0.7))

    for ax in (ax_a, ax_b):
        ax.set_xticks(x)
        ax.set_xticklabels(xt)
        ax.tick_params(labelsize=7.8)
        ax.set_xlim(-0.55, n - 0.45)
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.axvline(last_rel + 0.5, color="#555555", linestyle="--", linewidth=0.9)

    # c: proteins and NR50 families of the released releases; projections only where a formula is stated.
    sky, verm = "#56B4E9", "#D55E00"
    prot = np.array([float(r["proteins"]["value"]) / 1e6 if r.get("proteins") else np.nan for r in rels])
    nr50 = np.array([float(r["nr50"]["value"]) / 1e6 if r.get("nr50") else np.nan for r in rels])
    ri = np.array(released_idx)
    prefs = data["benchmarks"]["protein_refs"]
    proj = [(i, r["proteins_projection"]) for i, r in enumerate(rels) if r.get("proteins_projection")]
    ymax_c = max(float(np.nanmax(prot)), *(float(r["value"]) for r in prefs),
                 *(float(p["high"]) / 1e6 for _, p in proj)) * 1.17
    for ref in prefs:
        rv = float(ref["value"])
        ax_c.axhline(rv, color=ref["color"], linestyle=":", linewidth=1.0, alpha=0.8, zorder=2)
        ax_c.text(n - 0.5, rv + ymax_c * 0.008, ref["figure_label"], ha="right", va="bottom", fontsize=7.2,
                  color=ref["color"], style="italic")
    ax_c.plot(x[ri], prot[ri], color=sky, lw=2.2, marker="o", ms=6, zorder=4)
    ax_c.plot(x[ri], nr50[ri], color=verm, lw=2.2, marker="s", ms=5.6, zorder=4)
    for i in ri:
        ax_c.text(x[i], prot[i] + ymax_c * 0.03, f"{prot[i]:.1f}M", ha="center", va="bottom", fontsize=7.6,
                  color="#006699", fontweight="bold" if i == last_rel else "normal",
                  bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.8))
        ax_c.text(x[i] + 0.12, nr50[i] - ymax_c * 0.045, f"{nr50[i]:.1f}M", ha="left", va="top", fontsize=7.6,
                  color=verm, fontweight="bold" if i == last_rel else "normal")
    for i, p in proj:
        lo, hi, mid = float(p["low"]) / 1e6, float(p["high"]) / 1e6, float(p["value"]) / 1e6
        ax_c.plot([x[last_rel], x[i]], [prot[last_rel], mid], color=sky, lw=1.6, ls="--", zorder=3)
        ax_c.errorbar(x[i], mid, yerr=[[mid - lo], [hi - mid]], fmt="o", ms=6.5, mfc="white", mec=sky, mew=1.8,
                      ecolor=sky, elinewidth=1.2, capsize=3, zorder=5)
        ax_c.text(x[i] + 0.13, mid, f"{lo:.0f} to {hi:.0f}M\n(estimate)", ha="left", va="center", fontsize=7.2,
                  color="#006699")
    cf = data["benchmarks"].get("counterfactual")
    if cf:
        a = rids.index(cf["anchor_release"])
        end = rids.index(cf["through_release"])
        rate = float(cf["rate_nr50_per_gbp"]["value"]) / 1e6
        xs = x[a: end + 1]
        ys = [nr50[a] + (cols[rids[j]]["total"] - cols[rids[a]]["total"]) * rate for j in range(a, end + 1)]
        ax_c.plot(xs, ys, color="#666666", lw=1.3, ls=":", zorder=3)
        ax_c.plot(xs[-1], ys[-1], marker="d", ms=4.5, color="#666666", zorder=4)
        ax_c.text(xs[-1] + 0.12, ys[-1] - ymax_c * 0.01, cf["figure_label"].format(val=ys[-1]), ha="left",
                  va="top", fontsize=7.0, color="#555555", style="italic")
    for i, r in enumerate(rels):
        if r["status"] != "released" and not r.get("proteins_projection"):
            ax_c.text(x[i], ymax_c * 0.04, "not\nprojected", ha="center", va="bottom", fontsize=6.8,
                      color="#888888", style="italic")
    ax_c.text(x[0] + 0.1, prot[ri].max() * 0.62 + 8, "Distinct proteins", color="#006699", fontsize=8.6,
              fontweight="bold", ha="left")
    ax_c.text(x[last_rel] - 0.15, nr50[last_rel] + ymax_c * 0.075, "NR50 families", color=verm, fontsize=8.6,
              fontweight="bold", ha="right")
    co = data.get("callout")
    if co:
        ci = rids.index(co["release_id"])
        ax_c.annotate(co["figure_text"].format(growth=int(co["growth"]["value"])), xy=(x[ci] - 0.05, nr50[ci] + ymax_c * 0.01),
                      xytext=(max(0.2, ci - 3.2), ymax_c * 0.42), fontsize=7.6, color="#222222",
                      arrowprops=dict(arrowstyle="->", color="#666666", lw=0.9, connectionstyle="arc3,rad=-0.15"),
                      bbox=dict(boxstyle="round,pad=0.3", facecolor="#F8F8F8", edgecolor="#CCCCCC", lw=0.7), zorder=6)
    ax_c.axvline(last_rel + 0.5, color="#555555", linestyle="--", linewidth=0.9)
    ax_c.set_xlim(-0.55, n - 0.45)
    ax_c.set_ylim(0, ymax_c)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(xt)
    ax_c.tick_params(labelsize=7.8)
    ax_c.set_ylabel("Count (millions)", fontsize=9.2)
    ax_c.set_title("c   Distinct proteins and NR50 families against UniRef and ClusteredNR", loc="left",
                   fontsize=10, fontweight="bold", pad=10)
    ax_c.grid(axis="y", color="#E5E5E5", linewidth=0.6)
    ax_c.set_axisbelow(True)

    # d: new NR50 families per stored Gbp by stage.
    rows = sorted(data["stage_novelty"], key=lambda r: -float(r["per_gbp"]["value"]))
    base = next(r for r in rows if r["stage"] == data["novelty_baseline_stage"])
    bv = float(base["per_gbp"]["value"])
    yp = np.arange(len(rows))[::-1]
    vals = [float(r["per_gbp"]["value"]) / 1e3 for r in rows]
    stage_class = {s["stage"]: s["class"] for s in data["stages"]}
    ax_d.barh(yp, vals, height=0.64, color=[color[stage_class[r["stage"]]] for r in rows], edgecolor="white")
    vmax = max(vals)
    for y, r, v in zip(yp, rows, vals):
        mult = float(r["per_gbp"]["value"]) / bv
        note = f"{v:,.1f}k" if r is base else f"{v:,.1f}k  ({mult:.0f}x {base['stage']})"
        ax_d.text(v + 0.02 * vmax, y, note, va="center", ha="left", fontsize=8.0, color="#222222")
    ax_d.set_yticks(yp)
    ax_d.set_yticklabels([r["figure_label"] for r in rows], fontsize=7.8, linespacing=1.15)
    ax_d.set_xlim(0, vmax * 1.42)
    ax_d.tick_params(axis="x", labelsize=7.8)
    ax_d.set_xlabel("New NR50 families per stored Gbp (thousands)", fontsize=9.2)
    ax_d.set_title("d   New NR50 families per stored Gbp, by stage (c35_check.py)", loc="left", fontsize=10,
                   fontweight="bold", pad=10)
    ax_d.text(0.99, 0.02, data["novelty_caption"], transform=ax_d.transAxes, ha="right", va="bottom",
              fontsize=6.9, color="#555555", style="italic")
    ax_d.grid(axis="x", color="#E5E5E5", linewidth=0.6)
    ax_d.set_axisbelow(True)

    # Legends: classes (colour) and states (fill or outline).
    h_cls = [mpatches.Patch(facecolor=c["color"], edgecolor="white", label=c["label"]) for c in classes]
    h_st = [mpatches.Patch(label=STATE_TEXT[s], **_state_style(s, "#7A7A7A")) for s in STATES]
    leg1 = fig.legend(handles=h_cls, loc="upper left", bbox_to_anchor=(0.055, 0.945), ncol=3, frameon=False,
                      fontsize=7.9, handlelength=1.3, columnspacing=1.3, title="Source class (colour)",
                      title_fontsize=8.2, alignment="left")
    fig.add_artist(leg1)
    fig.legend(handles=h_st, loc="upper left", bbox_to_anchor=(0.70, 0.945), ncol=1, frameon=False, fontsize=7.9,
               handlelength=1.6, title="State (fill or outline)", title_fontsize=8.2, alignment="left")
    fig.suptitle(data["figure_title"], x=0.055, y=0.985, ha="left", fontsize=12.5, fontweight="bold", color="#111111")
    bound = data["r26_11_5_bound_text"]["upper" if data.get("r26_11_5_upper_bound", True) else "exact"]
    fig.text(0.055, 0.012, data["footnote"].format(bound=bound), fontsize=7.0, color="#555555", ha="left", va="bottom", linespacing=1.3)
    return fig


def render(data: dict, out: Path, tower: dict | None) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _style()
    fig = build_figure(data, tower)
    out.parent.mkdir(parents=True, exist_ok=True)
    png, svg = Path(f"{out}.png"), Path(f"{out}.svg")
    fig.savefig(png, dpi=300, metadata={"Software": None})
    fig.savefig(svg, metadata={"Date": None, "Creator": None})
    plt.close(fig)
    return png, svg


# ---------------------------------------------------------------------------------------------------------------

def _data_with_overlay(args) -> dict:
    data = load(args.data)
    own_path = args.ownership or data.get("r26_11_5_ownership")
    if own_path:
        p = Path(own_path)
        p = p if p.is_absolute() else REPO_ROOT / p
        data = apply_ownership(data, json.loads(p.read_text(encoding="utf-8")), str(own_path))
    return data


def _print_table(data: dict) -> None:
    cols = table(data)
    cids = [c["id"] for c in data["classes"]]
    print("release\t" + "\t".join(cids) + "\t" + "\t".join(STATES) + "\ttotal\tlow\thigh\tD")
    for rid, col in cols.items():
        d = effective_classes([col["by_class"][c] for c in cids])
        print(rid + "\t" + "\t".join(f"{col['by_class'][c]:.3f}" for c in cids) + "\t"
              + "\t".join(f"{col['by_state'][s]:.3f}" for s in STATES)
              + f"\t{col['total']:.3f}\t{col['low']:.3f}\t{col['high']:.3f}\t{d:.2f}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("check", "tower", "render", "fill"):
        p = sub.add_parser(name)
        p.add_argument("--data", type=Path, default=DEFAULT_DATA)
        p.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
        if name in ("check", "render"):
            p.add_argument("--ownership", type=Path, default=None, help="card A-29's ownership JSON")
        if name == "render":
            p.add_argument("--out", type=Path, required=True, help="output prefix (writes .png and .svg)")
    r = sub.add_parser("root")
    r.add_argument("--root", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--releases", default="rN.3,rM.2,r26.11.1,r26.11.2,r26.11.3")
    r.add_argument("--built", default="")
    r.add_argument("--threads", type=int, default=2)
    r.add_argument("--memory", default="4GB")
    args = ap.parse_args(argv)

    if args.cmd == "root":
        doc = collect_root(args.root, [s for s in args.releases.split(",") if s],
                           [s for s in args.built.split(",") if s], args.threads, args.memory)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.out.with_name(f".{args.out.name}.tmp")
        tmp.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        tmp.replace(args.out)
        print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("stages", "levels")}
                          for k, v in doc["releases"].items()}, sort_keys=True))
        return 0
    if args.cmd == "fill":
        data = load(args.data)
        n = fill(data)
        tmp = args.data.with_name(f".{args.data.name}.tmp")
        tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        tmp.replace(args.data)
        print(f"filled {n} values from their formulas in {args.data}")
        return 0
    if args.cmd == "tower":
        t = tower_remainder(load(args.data), args.sources)
        for row in t["rows"]:
            print(f"{row['entry']}\t{row['assigned'] or '-'}\t{row['files']}\t{row['sequence_bytes']}\t{row['other_bytes']}")
        print(json.dumps({k: v for k, v in t.items() if k != "rows"}, sort_keys=True))
        return 0
    data = _data_with_overlay(args)
    if args.cmd == "check":
        problems = data_problems(data)
        t = tower_remainder(data, args.sources)
        if t["unknown_assigned"]:
            problems.append(f"tower: assigned entries not in sources.yaml: {t['unknown_assigned']}")
        _print_table(data)
        for p in problems:
            print(f"problem: {p}")
        print(f"checks: {'all pass' if not problems else f'{len(problems)} problems'}")
        return 0 if not problems else 1
    png, svg = render(data, args.out, tower_remainder(data, args.sources))
    print(f"wrote {png} and {svg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

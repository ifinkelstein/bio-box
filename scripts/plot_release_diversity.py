#!/usr/bin/env python3
"""Render the four-panel bio-box release diversity figure (Claus Wilke / Okabe-Ito).

Reads canonical measured/queued release numbers and pinned external benchmarks
from ``planning/measurements/release_diversity.json`` and writes both a
300-DPI PNG and an editable SVG.

Standard release-closeout invocation (with dataviz-python / dvstyle):
    uv run --project /Users/ilya/projects/dataviz-python python3 \\
        scripts/plot_release_diversity.py \\
        --out planning/measurements/2026-11-24_release_r26.11.3/source_diversity_r26.11.3
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA = REPO_ROOT / "planning" / "measurements" / "release_diversity.json"
DEFAULT_DATAVIZ_DIR = Path("/Users/ilya/projects/dataviz-python/python")

OKABE_ITO = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
}


def _apply_style(dataviz_dir: Path | None) -> tuple[object | None, dict[str, str]]:
    """Apply dvstyle('paper') if available, else an equivalent Wilke matplotlib theme."""
    dvstyle = None
    oi = dict(OKABE_ITO)
    if dataviz_dir is not None and dataviz_dir.is_dir():
        sys.path.insert(0, str(dataviz_dir))
        try:
            import dvstyle as _dvstyle  # type: ignore[import-not-found]
            import palettes as _palettes  # type: ignore[import-not-found]

            _dvstyle.set_style("paper")
            dvstyle = _dvstyle
            oi = dict(_palettes.OKABE_ITO)
            return dvstyle, oi
        except ImportError:
            pass

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Source Sans 3", "Source Sans Pro", "DejaVu Sans", "Arial"],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#333333",
            "axes.linewidth": 0.8,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "svg.fonttype": "none",
            "savefig.dpi": 300,
        }
    )
    return dvstyle, oi


def _effective_classes(gbp_vec: list[float]) -> float:
    """Hill number q=1: exp(Shannon entropy) across biological source classes."""
    total = sum(gbp_vec)
    if total <= 0:
        return 0.0
    h = 0.0
    for v in gbp_vec:
        if v > 0:
            p = v / total
            h -= p * math.log(p)
    return math.exp(h)


def build_figure(
    data: dict,
    *,
    include_queued: bool = True,
    include_benchmarks: bool = True,
    max_novelty_bars: int = 10,
    oi: dict[str, str] = OKABE_ITO,
) -> plt.Figure:
    categories = data["categories"]
    releases = [
        r for r in data["releases"] if include_queued or r.get("status") != "queued"
    ]
    if not releases:
        raise ValueError("No releases selected to plot.")

    n_rel = len(releases)
    measured_idx = [i for i, r in enumerate(releases) if r.get("status") != "queued"]
    queued_idx = [i for i, r in enumerate(releases) if r.get("status") == "queued"]
    last_measured = measured_idx[-1] if measured_idx else -1

    rel_labels = [
        f"{r['id']}{'*' if r.get('status') == 'queued' else ''}\n({r['date_label']})"
        for r in releases
    ]

    gbp_by_cat = {
        cat["id"]: np.array(
            [float(r["gbp_by_category"].get(cat["id"], 0.0)) for r in releases],
            dtype=float,
        )
        for cat in categories
    }
    total_gbp = sum(gbp_by_cat.values())
    prot_m = np.array([float(r["proteins_m"]) for r in releases], dtype=float)
    nr50_m = np.array([float(r["nr50_m"]) for r in releases], dtype=float)
    benchmarks = data.get("benchmarks", {}) if include_benchmarks else {}

    # Scale figure width gracefully as more releases accumulate
    fig_w = max(12.6, 10.8 + 0.30 * n_rel)
    fig = plt.figure(figsize=(fig_w, 8.7))
    gs_outer = gridspec.GridSpec(2, 1, height_ratios=[1.08, 1.04], hspace=0.43)
    gs_top = gs_outer[0].subgridspec(1, 2, width_ratios=[1.02, 1.0], wspace=0.19)
    gs_bot = gs_outer[1].subgridspec(1, 2, width_ratios=[1.08, 1.0], wspace=0.43)

    ax_comp = fig.add_subplot(gs_top[0, 0])
    ax_gbp = fig.add_subplot(gs_top[0, 1])
    ax_nr50 = fig.add_subplot(gs_bot[0, 0])
    ax_nov = fig.add_subplot(gs_bot[0, 1])

    x = np.arange(n_rel)
    bar_w = min(0.62, 3.8 / max(n_rel, 1))

    # --- Panel a: Proportional source composition (% of stored Gbp) + exp(H) ---
    bottom_pct = np.zeros(n_rel)
    for cat in categories:
        key = cat["id"]
        color = cat["color"]
        pct = 100.0 * gbp_by_cat[key] / total_gbp
        bars = ax_comp.bar(
            x,
            pct,
            width=bar_w,
            bottom=bottom_pct,
            color=color,
            edgecolor="white",
            linewidth=0.8,
            label=cat["label"],
        )
        for qi in queued_idx:
            bars[qi].set_hatch("///")
            bars[qi].set_edgecolor("white")
        for i, val in enumerate(pct):
            if val >= 8.0:
                txt_col = "#111111" if cat.get("dark_text", False) else "white"
                bbox_kw = (
                    dict(boxstyle="round,pad=0.16", facecolor=color, edgecolor="none", alpha=0.92)
                    if i in queued_idx
                    else None
                )
                ax_comp.text(
                    x[i],
                    bottom_pct[i] + val / 2.0,
                    f"{val:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=8.0,
                    fontweight="bold",
                    color=txt_col,
                    bbox=bbox_kw,
                )
        bottom_pct += pct

    if include_benchmarks:
        for i, r in enumerate(releases):
            eff = _effective_classes(
                [float(r["gbp_by_category"].get(c["id"], 0.0)) for c in categories]
            )
            lbl = f"D={eff:.1f}" if i == 0 else f"{eff:.1f}"
            ax_comp.text(
                x[i],
                102.2,
                lbl,
                ha="center",
                va="bottom",
                fontsize=7.9,
                fontweight="bold" if i >= last_measured else "normal",
                color="#333333",
            )
        ax_comp.text(
            n_rel - 0.55,
            112.0,
            "D = exp(H) effective classes (single-domain DB = 1.0, max = 6.0)",
            ha="right",
            va="bottom",
            fontsize=7.3,
            color="#555555",
            style="italic",
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.88),
            zorder=5,
        )

    if queued_idx and last_measured >= 0:
        ax_comp.axvline(last_measured + 0.5, color="#555555", linestyle="--", linewidth=0.9, alpha=0.85)
    ax_comp.set_ylim(0, 118 if include_benchmarks else 105)
    ax_comp.set_yticks([0, 25, 50, 75, 100])
    ax_comp.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_comp.set_xticks(x)
    ax_comp.set_xticklabels(rel_labels)
    ax_comp.tick_params(labelsize=8.3)
    ax_comp.set_ylabel("Share of stored sequence (% of Gbp)", fontsize=9.2)
    ax_comp.set_title(
        "a   Source composition shifts from isolates to environmental & viral diversity",
        loc="left",
        fontsize=9.8,
        fontweight="bold",
        pad=10,
    )
    ax_comp.grid(axis="y", color="#E5E5E5", linewidth=0.6)
    ax_comp.set_axisbelow(True)

    # --- Panel b: Absolute stored sequence (Gbp) + OpenGenome 1 benchmark ---
    bottom_gbp = np.zeros(n_rel)
    for cat in categories:
        vals = gbp_by_cat[cat["id"]]
        bars = ax_gbp.bar(
            x,
            vals,
            width=bar_w,
            bottom=bottom_gbp,
            color=cat["color"],
            edgecolor="white",
            linewidth=0.8,
        )
        for qi in queued_idx:
            bars[qi].set_hatch("///")
            bars[qi].set_edgecolor("white")
        bottom_gbp += vals

    gbp_refs = benchmarks.get("gbp_refs", [])
    max_gbp = max([float(np.max(total_gbp)), *(float(r["value_gbp"]) for r in gbp_refs)])
    gbp_pad = max(6.5, 0.021 * max_gbp)

    for ref in gbp_refs:
        rv = float(ref["value_gbp"])
        ax_gbp.axhline(rv, color="#666666", linestyle=":", linewidth=1.1, alpha=0.85, zorder=2)
        ax_gbp.text(
            -0.38,
            rv + gbp_pad * 0.65,
            ref["label"],
            ha="left",
            va="bottom",
            fontsize=7.8,
            color="#444444",
            style="italic",
        )

    for i, tot in enumerate(total_gbp):
        suffix = " Gbp" if i == n_rel - 1 else ""
        # Keep bar value crisp even if close to a horizontal benchmark line
        y_txt = tot + gbp_pad
        bbox_b = None
        if any(abs(y_txt - float(r["value_gbp"])) < 14.0 for r in gbp_refs):
            bbox_b = dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="none", alpha=0.88)
        ax_gbp.text(
            x[i],
            y_txt,
            f"{tot:.1f}{suffix}",
            ha="center",
            va="bottom",
            fontsize=8.2,
            fontweight="bold" if i >= last_measured else "normal",
            color="#222222",
            bbox=bbox_b,
            zorder=5,
        )

    if queued_idx and last_measured >= 0:
        ax_gbp.axvline(last_measured + 0.5, color="#555555", linestyle="--", linewidth=0.9, alpha=0.85)
    ax_gbp.set_ylim(0, max_gbp * 1.15)
    ax_gbp.set_xticks(x)
    ax_gbp.set_xticklabels(rel_labels)
    ax_gbp.tick_params(labelsize=8.3)
    ax_gbp.set_ylabel("Stored sequence (Gbp)", fontsize=9.2)
    ax_gbp.set_title(
        f"b   Cumulative corpus size by source class ({total_gbp[0]:.1f} → {total_gbp[-1]:.1f} Gbp)",
        loc="left",
        fontsize=9.8,
        fontweight="bold",
        pad=10,
    )
    ax_gbp.grid(axis="y", color="#E5E5E5", linewidth=0.6)
    ax_gbp.set_axisbelow(True)

    # --- Panel c: Protein & NR50 growth + external & counterfactual benchmarks ---
    sky_blue = oi["sky_blue"]
    vermillion = oi["vermillion"]
    prot_refs = benchmarks.get("protein_refs", [])
    y_max_c = max([float(np.max(prot_m)) * 1.24, *(float(r["value_m"]) * 1.16 for r in prot_refs)])

    # Horizontal external benchmarks drawn softly behind data
    for ref in prot_refs:
        rv = float(ref["value_m"])
        rcol = ref.get("color", "#666666")
        ax_nr50.axhline(rv, color=rcol, linestyle=":", linewidth=1.0, alpha=0.72, zorder=2)
        if rv < 60.0:
            # UniRef50 (41.2M): two-line label in the completely open top-left pocket above rN.3 (10.6M)
            u_lbl = ref["label"].replace(", 50%", ",\n50%")
            ax_nr50.text(
                -0.30,
                rv + 4.5,
                u_lbl,
                ha="left",
                va="bottom",
                fontsize=7.1,
                linespacing=1.12,
                color=rcol,
                style="italic",
                zorder=3,
            )
        else:
            ax_nr50.text(
                -0.30,
                rv + 4.5,
                ref["label"],
                ha="left",
                va="bottom",
                fontsize=7.4,
                color=rcol,
                style="italic",
                zorder=3,
            )

    # Isolate-only counterfactual (C5 rate projected from r26.11.2 onward)
    cf = benchmarks.get("counterfactual")
    if cf:
        anchor_candidates = [i for i, r in enumerate(releases) if r["id"] == cf["anchor_release"]]
        if anchor_candidates:
            a_idx = anchor_candidates[0]
            rate_m = float(cf["rate_nr50_k_per_gbp"]) / 1000.0
            cf_x = x[a_idx:]
            cf_y = nr50_m[a_idx] + (total_gbp[a_idx:] - total_gbp[a_idx]) * rate_m
            ax_nr50.plot(cf_x, cf_y, color="#666666", linewidth=1.4, linestyle=":", zorder=3)
            ax_nr50.plot(cf_x[-1], cf_y[-1], marker="d", markersize=4.5, color="#666666", zorder=4)
            ax_nr50.text(
                cf_x[-1] + 0.28,
                cf_y[-1] + 4.0,
                f"C5 isolate rate (~{cf_y[-1]:.1f}M)",
                ha="right",
                va="bottom",
                fontsize=7.2,
                color="#555555",
                style="italic",
                bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.82),
            )

    if measured_idx:
        ax_nr50.plot(
            x[measured_idx], prot_m[measured_idx],
            color=sky_blue, linewidth=2.2, marker="o", markersize=6.0, zorder=4,
        )
        ax_nr50.plot(
            x[measured_idx], nr50_m[measured_idx],
            color=vermillion, linewidth=2.2, marker="s", markersize=5.8, zorder=4,
        )
    if queued_idx and last_measured >= 0:
        seg = [last_measured, *queued_idx]
        ax_nr50.plot(x[seg], prot_m[seg], color=sky_blue, linewidth=1.8, linestyle="--", zorder=3)
        ax_nr50.plot(x[seg], nr50_m[seg], color=vermillion, linewidth=1.8, linestyle="--", zorder=3)
        for qi in queued_idx:
            ax_nr50.plot(
                x[qi], prot_m[qi], marker="o", markersize=6.8,
                markerfacecolor="white", markeredgecolor=sky_blue, markeredgewidth=2.0, zorder=5,
            )
            ax_nr50.plot(
                x[qi], nr50_m[qi], marker="s", markersize=6.5,
                markerfacecolor="white", markeredgecolor=vermillion, markeredgewidth=2.0, zorder=5,
            )

    callout_rel = data.get("callout", {}).get("release_id")
    callout_idx = next((i for i, r in enumerate(releases) if r["id"] == callout_rel), None)

    for i, r in enumerate(releases):
        is_q = i in queued_idx
        p_lbl = f"~{prot_m[i]:.0f}M" if is_q else f"{prot_m[i]:.1f}M"
        n_lbl = f"~{nr50_m[i]:.0f}M" if is_q else f"{nr50_m[i]:.1f}M"

        # Distinct-protein label (offset left for final queued dot if near ClusteredNR 286M line)
        if is_q and any(abs(prot_m[i] + 15 - float(ref["value_m"])) < 18.0 for ref in prot_refs):
            ax_nr50.text(
                x[i] - 0.16,
                prot_m[i] - 2.0,
                p_lbl,
                ha="right",
                va="center",
                fontsize=8.0,
                color="#006699",
                fontweight="bold",
            )
        else:
            ax_nr50.text(
                x[i],
                prot_m[i] + (13.0 if i == 0 else 11.0),
                p_lbl,
                ha="center",
                va="bottom",
                fontsize=7.8,
                color="#006699",
                fontweight="bold" if i >= last_measured else "normal",
                bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.75)
                if i in (1, 2) and include_benchmarks
                else None,
            )

        # NR50 label
        if i == 0 and (prot_m[0] - nr50_m[0]) < 20.0:
            ax_nr50.text(
                x[i] + 0.36,
                nr50_m[i] + 4.5,
                n_lbl,
                ha="left",
                va="bottom",
                fontsize=7.8,
                color=vermillion,
            )
        elif i == callout_idx:
            ax_nr50.text(
                x[i] + 0.13,
                nr50_m[i] + 9.5,
                n_lbl,
                ha="left",
                va="bottom",
                fontsize=7.8,
                color=vermillion,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.82),
            )
        else:
            ax_nr50.text(
                x[i],
                nr50_m[i] + 9.5,
                n_lbl,
                ha="center",
                va="bottom",
                fontsize=7.8,
                color=vermillion,
                fontweight="bold" if i >= last_measured else "normal",
            )

    # Direct series labels in open whitespace (clear of UniRef100 227M and 22.2M)
    ax_nr50.text(
        max(0.5, n_rel - 1.62),
        202.0 if include_benchmarks else y_max_c * 0.775,
        "Distinct proteins",
        color="#006699",
        fontsize=8.6,
        fontweight="bold",
        ha="right",
    )
    ax_nr50.text(
        n_rel - 0.68,
        nr50_m[-1] + 32.0,
        "NR50 families",
        color=vermillion,
        fontsize=8.6,
        fontweight="bold",
        ha="right",
    )

    if callout_idx is not None and data.get("callout", {}).get("text"):
        ax_nr50.annotate(
            data["callout"]["text"],
            xy=(x[callout_idx] - 0.08, nr50_m[callout_idx] + 3.0),
            xytext=(max(0.8, callout_idx - 2.65), 148.0),
            fontsize=7.9,
            color="#222222",
            arrowprops=dict(arrowstyle="->", color="#666666", lw=0.9, connectionstyle="arc3,rad=-0.16"),
            bbox=dict(boxstyle="round,pad=0.30", facecolor="#F8F8F8", edgecolor="#CCCCCC", lw=0.7),
            zorder=6,
        )

    if queued_idx and last_measured >= 0:
        ax_nr50.axvline(last_measured + 0.5, color="#555555", linestyle="--", linewidth=0.9, alpha=0.85)
    ax_nr50.set_xlim(-0.38, n_rel - 0.60)
    ax_nr50.set_ylim(0, y_max_c)
    ax_nr50.set_xticks(x)
    ax_nr50.set_xticklabels(rel_labels)
    ax_nr50.tick_params(labelsize=8.3)
    ax_nr50.set_ylabel("Count (millions)", fontsize=9.2)
    ax_nr50.set_title(
        "c   Protein & NR50 diversity surge vs UniRef, ClusteredNR & isolate baseline",
        loc="left",
        fontsize=9.6,
        fontweight="bold",
        pad=10,
    )
    ax_nr50.grid(axis="y", color="#E5E5E5", linewidth=0.6)
    ax_nr50.set_axisbelow(True)

    # --- Panel d: Measured novelty yield per Gbp by added stage ---
    cat_color = {c["id"]: c["color"] for c in categories}
    baseline_stage = data["novelty_baseline_stage"]["stage"]
    baseline_val = float(data["novelty_baseline_stage"]["exact_value_k"])

    nov_rows = sorted(
        data["stage_novelty_nr50_per_gbp_k"],
        key=lambda row: float(row["value_k"]),
        reverse=True,
    )[:max_novelty_bars]
    y_pos = np.arange(len(nov_rows))[::-1]
    labels_d = [row["label"] for row in nov_rows]
    vals_d = [float(row["value_k"]) for row in nov_rows]
    cols_d = [cat_color[row["category"]] for row in nov_rows]

    ax_nov.barh(y_pos, vals_d, height=0.64, color=cols_d, edgecolor="white", linewidth=0.6)
    max_nov = max(vals_d)
    for yp, row in zip(y_pos, nov_rows):
        val = float(row["value_k"])
        mult = val / baseline_val
        note = (
            f"{val:,.1f}k"
            if row["stage"] == baseline_stage
            else f"{val:,.1f}k  ({mult:.0f}× {baseline_stage})"
        )
        ax_nov.text(val + 0.022 * max_nov, yp, note, va="center", ha="left", fontsize=8.2, color="#222222")

    ax_nov.set_yticks(y_pos)
    ax_nov.set_yticklabels(labels_d, fontsize=8.0, linespacing=1.15)
    ax_nov.set_xlim(0, max_nov * 1.41)
    ax_nov.tick_params(axis="x", labelsize=8.3)
    ax_nov.set_xlabel("New NR50 families per stored Gbp (thousands)", fontsize=9.2)
    ax_nov.set_title(
        "d   Measured novelty yield per Gbp by source (c35_check.py)",
        loc="left",
        fontsize=9.8,
        fontweight="bold",
        pad=10,
    )
    ax_nov.grid(axis="x", color="#E5E5E5", linewidth=0.6)
    ax_nov.set_axisbelow(True)

    # Shared horizontal legend centered between title and top panels
    legend_handles = [
        mpatches.Patch(facecolor=c["color"], edgecolor="white", label=c["label"])
        for c in categories
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.51, 0.942),
        ncol=3,
        frameon=False,
        fontsize=8.0,
        handlelength=1.25,
        columnspacing=1.35,
    )

    first_id = releases[0]["id"]
    last_meas_id = releases[last_measured]["id"] if last_measured >= 0 else releases[-1]["id"]
    queued_part = (
        f"; {', '.join(r['id'] + '*' for r in releases if r.get('status') == 'queued')} staged/queued"
        if queued_idx
        else ""
    )
    fig.suptitle(
        f"bio-box corpus diversity across releases ({first_id} → {last_meas_id} measured{queued_part})",
        x=0.055,
        y=0.982,
        ha="left",
        fontsize=12.2,
        fontweight="bold",
        color="#111111",
    )
    if data.get("footnote"):
        fig.text(0.055, 0.012, data["footnote"], fontsize=7.2, color="#555555", ha="left")

    fig.subplots_adjust(left=0.065, right=0.975, top=0.835, bottom=0.108)
    return fig


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Render bio-box release source diversity figure.")
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA, help="Path to release_diversity.json")
    ap.add_argument("--out", type=Path, required=True, help="Output path prefix (writes .png and .svg)")
    ap.add_argument("--dataviz-dir", type=Path, default=DEFAULT_DATAVIZ_DIR, help="Path to dataviz-python/python")
    ap.add_argument("--no-queued", action="store_true", help="Exclude status='queued' releases")
    ap.add_argument("--no-benchmarks", action="store_true", help="Omit external benchmarks, exp(H), and counterfactual")
    ap.add_argument("--max-novelty-bars", type=int, default=10, help="Max stage bars in Panel d")
    args = ap.parse_args(argv)

    dvstyle, oi = _apply_style(args.dataviz_dir)
    data = json.loads(args.data.read_text(encoding="utf-8"))
    fig = build_figure(
        data,
        include_queued=not args.no_queued,
        include_benchmarks=not args.no_benchmarks,
        max_novelty_bars=args.max_novelty_bars,
        oi=oi,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    png_path = Path(f"{args.out}.png")
    svg_path = Path(f"{args.out}.svg")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png_path} and {svg_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

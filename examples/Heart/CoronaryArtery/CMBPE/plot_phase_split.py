#!/usr/bin/env python3
"""Coronary closure comparison: where in the cycle the constant closures differ
from R_mu, and how dissipation and wall shear stress follow the flow change.

Reads <coronary>/coronary_summary.csv (compare_coronary_closures.py).
  (a) change of mean |Q_in| relative to R_mu, split into systole (dark) and
      diastole (light, hatched), for R_0 (operating-point constant) and R_inf
      (high-shear constant), per blood state;
  (b) change of cycle dissipation and of area-mean WSS against the change of
      mean |Q_in|, with the lines y = x and y = 2x.

  python3 plot_phase_split.py ~/Results/cmbpe/coronary/coronary_summary.csv
"""
import argparse
import csv
import os
import sys

GROUPS = [("pra700/healthy", "Healthy"), ("pra700/diabetic", "Type 2\ndiabetes"),
          ("pra700/hyperviscous", "Hyper-\nviscous"),
          ("pra2400/hyperviscous", "Hyperviscous\n$P_{RA}$=2400 Pa")]
CLOSURES = [("rest", r"$\mathcal{R}_0$ (resting point)", "#c4581f"),
            ("plateau", r"$\mathcal{R}_\infty$ (high-shear)", "#7a4fa0")]


def tint(hex_color, f):
    """Mix with white: f = 0 keeps the colour, f = 1 gives white."""
    c = [int(hex_color[i:i + 2], 16) for i in (1, 3, 5)]
    return "#%02x%02x%02x" % tuple(round(x + (255 - x) * f) for x in c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("summary")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rows = {(r["group"], r["case"]): r for r in csv.DictReader(open(a.summary))}

    groups = [(g, lab) for g, lab in GROUPS
              if any((g, c) in rows and rows[(g, c)].get("dQ_abs_%") for c, _, _ in CLOSURES)]
    if not groups:
        sys.exit("no rest/plateau rows with differences in " + a.summary)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
        "axes.labelsize": 8.5, "legend.fontsize": 8, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "hatch.linewidth": 0.5,
        "xtick.direction": "out", "ytick.direction": "out",
    })
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.0),
                                 gridspec_kw={"width_ratios": [1.35, 1]})

    # (a) systolic / diastolic split of the |Q| change
    width, gap = 0.34, 0.04
    for i, (g, _) in enumerate(groups):
        for k, (c, _, col) in enumerate(CLOSURES):
            r = rows.get((g, c))
            if not r or not r.get("dQ_abs_%"):
                continue
            s, d = float(r["dQ_abs_sys_%"]), float(r["dQ_abs_dia_%"])
            x = i + (k - 0.5) * (width + gap)
            ax.bar(x, s, width, color=col, edgecolor="white", linewidth=1.0, zorder=3)
            ax.bar(x, d, width, bottom=s, color=tint(col, 0.6), edgecolor="white",
                   linewidth=1.0, hatch="////", zorder=3)
            tot = s + d
            share = 100 * (d if d >= s else s) / tot
            ax.text(x, tot + 0.15, "%.1f%%\n%.0f%% %s" % (tot, share, "dia" if d >= s else "sys"),
                    ha="center", va="bottom", fontsize=7, color="#303030", linespacing=1.1)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([lab for _, lab in groups])
    ax.set_ylabel(r"Change of mean $|Q_{\rm in}|$ vs $\mathcal{R}_\mu$ (%)")
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_ylim(0, 1.3 * ax.get_ylim()[1])
    ax.grid(True, axis="y", lw=0.4, color="#d9d9d9", zorder=0)
    ax.axhline(0, color="#404040", lw=0.6, zorder=2)
    handles = [Patch(color=col, label=lab) for _, lab, col in CLOSURES]
    handles += [Patch(facecolor="#7d7d7d", edgecolor="white", label="systole"),
                Patch(facecolor=tint("#7d7d7d", 0.6), edgecolor="white", hatch="////",
                      label="diastole")]
    ax.legend(handles=handles, loc="upper left", frameon=False, ncol=2,
              handlelength=1.4, columnspacing=1.0)
    ax.text(-0.13, 1.0, "(a)", transform=ax.transAxes, fontsize=9, va="top")

    # (b) dissipation and WSS against the flow change
    markers = ["o", "s", "^", "D"]
    xmax = 0.0
    for i, (g, _) in enumerate(groups):
        for c, _, col in CLOSURES:
            r = rows.get((g, c))
            if not r or not r.get("dQ_abs_%"):
                continue
            q = float(r["dQ_abs_%"])
            xmax = max(xmax, q)
            bx.plot(q, float(r["dE_diss_%"]), markers[i], color=col, ms=6, mew=1.0, zorder=3)
            bx.plot(q, float(r["dWSS_mean_%"]), markers[i], mfc="white", mec=col, ms=6,
                    mew=1.1, zorder=3)
    xm = 1.15 * xmax
    bx.plot([0, xm], [0, xm], color="#909090", lw=0.8, ls="--", zorder=1)
    bx.plot([0, xm], [0, 2 * xm], color="#909090", lw=0.8, ls=":", zorder=1)
    bx.text(xm, xm, " $y=x$", fontsize=7, color="#606060", va="center")
    bx.text(xm / 2, xm, "$y=2x$ ", fontsize=7, color="#606060", ha="right", va="center")
    bx.set_xlim(0, xm * 1.12)
    bx.set_ylim(0, 2.1 * xm)
    bx.set_xlabel(r"Change of mean $|Q_{\rm in}|$ (%)")
    bx.set_ylabel("Change vs $\\mathcal{R}_\\mu$ (%)")
    bx.grid(True, lw=0.4, color="#e3e3e3", zorder=0)
    h2 = [Line2D([], [], marker="o", ls="none", color="#505050", ms=5.5, label="dissipation"),
          Line2D([], [], marker="o", ls="none", mfc="white", mec="#505050", ms=5.5,
                 label="mean WSS")]
    h2 += [Line2D([], [], marker=markers[i], ls="none", color="#9a9a9a", ms=5,
                  label=lab.replace("\n", " ").replace("- ", ""))
           for i, (_, lab) in enumerate(groups)]
    bx.legend(handles=h2, loc="upper left", frameon=False, handletextpad=0.3)
    bx.text(-0.2, 1.0, "(b)", transform=bx.transAxes, fontsize=9, va="top")

    for axx in (ax, bx):
        for side in ("top", "right"):
            axx.spines[side].set_visible(False)
    fig.tight_layout(w_pad=2.0)
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(a.summary)), "fig_phase_split.pdf")
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(os.path.splitext(out)[0] + ".png", dpi=300, bbox_inches="tight")
    print("written", out)


if __name__ == "__main__":
    sys.exit(main())

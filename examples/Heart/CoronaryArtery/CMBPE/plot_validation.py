#!/usr/bin/env python3
"""Figure of the truncation benchmark: error of the reduced closures against
the resolved solution, as a function of the resting operating point.

Reads <validation>/main/summary.csv (written by compare_rmu_axi.py) and, for
the Womersley numbers, the calibration lines of the full runs' log.txt.

Layout: one column per inlet radius (Womersley regime), two rows (waveform
error and dissipation error). Colour = closure, marker/line = blood.

  python3 plot_validation.py ~/Results/cmbpe/validation
  python3 plot_validation.py ~/Results/cmbpe/validation --with-plateau --out fig.pdf
"""
import argparse
import csv
import glob
import os
import re
import sys

CLOSURES = {
    "rest": (r"constant $\mathcal{R}$", "#d95f02", "-"),
    "single": (r"$\mathcal{R}_\mu$", "#1b4f9c", "-"),
    "tree": (r"$\mathcal{R}_\mu$ (all generations)", "#1b4f9c", ":"),
    "plateau": (r"constant $\mathcal{R}$ ($\mu_\infty$)", "#7f7f7f", "--"),
}
BLOODS = {
    "healthy": ("Healthy", "o"),
    "diabetic": ("Type 2 diabetes", "s"),
    "hyperviscous": ("Hyperviscous", "^"),
}
ROWS = [
    ("dQ_rms", r"Flow-waveform error (%)"),
    ("dE", r"Dissipation error (%)"),
]
FLOOR = 1e-2  # per cent; values below are drawn at the floor (log axis)


def womersley(validation, scale):
    """min-max Womersley of the narrowest tube over the full runs of a radius."""
    wo = []
    for log in glob.glob(os.path.join(validation, "main", "s%s_g*" % scale, "*", "full", "log.txt")):
        txt = open(log, errors="ignore").read()
        found = re.findall(r"tube R=([0-9.eE+-]+) L=[0-9.eE+-]+ .*?Wo=([0-9.eE+-]+)", txt)
        if found:
            rmin = min(float(r) for r, _ in found)
            wo += [float(w) for r, w in found if float(r) == rmin]
    return (min(wo), max(wo)) if wo else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("validation", help="the validation output folder")
    ap.add_argument("--with-plateau", action="store_true", help="add the mu_inf constant")
    ap.add_argument("--with-tree", action="store_true", help="add the exact multi-tube R_mu")
    ap.add_argument("--out", default=None, help="figure file (default <validation>/fig_validation.pdf)")
    a = ap.parse_args()

    path = os.path.join(a.validation, "main", "summary.csv")
    if not os.path.isfile(path):
        sys.exit("%s not found: run compare_rmu_axi.py first" % path)
    rows = list(csv.DictReader(open(path)))

    closures = ["rest", "single"] + (["tree"] if a.with_tree else []) + \
               (["plateau"] if a.with_plateau else [])

    data = {}  # (scale, closure, state) -> {gamma: row}
    scales = set()
    for r in rows:
        m = re.match(r"s([0-9.eE+-]+)_g([0-9.eE+-]+)$", r["group"])
        if not m:
            continue
        scale, gamma = m.group(1), float(m.group(2))
        scales.add(scale)
        data.setdefault((scale, r["closure"], r["state"]), {})[gamma] = r
    scales = sorted(scales, key=float)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullFormatter

    plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
                         "font.family": "serif", "mathtext.fontset": "cm"})
    fig, axes = plt.subplots(len(ROWS), len(scales), figsize=(2.6 * len(scales) + 0.6, 5.0),
                             sharex=True, sharey="row", squeeze=False)

    for j, scale in enumerate(scales):
        wo = womersley(a.validation, scale)
        title = r"$R_1 = %s$ mm" % ("%g" % (1e3 * float(scale)))
        if wo:
            title += "\n" + (r"$\mathrm{Wo} \approx %.2g$" % wo[0] if abs(wo[1] / wo[0] - 1) < 0.15
                             else r"$\mathrm{Wo} = %.2g$–$%.2g$" % wo)
        axes[0][j].set_title(title)
        for i, (key, label) in enumerate(ROWS):
            ax = axes[i][j]
            for c in closures:
                name, colour, style = CLOSURES[c]
                for s, (_, marker) in BLOODS.items():
                    d = data.get((scale, c, s))
                    if not d:
                        continue
                    g = sorted(d)
                    y = [max(abs(float(d[x][key])), FLOOR) for x in g]
                    ax.plot(g, y, color=colour, ls=style, marker=marker, ms=4, lw=1.2,
                            mfc="white" if c == "rest" else colour)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(True, which="both", lw=0.3, alpha=0.5)
            if j == 0:
                ax.set_ylabel(label)
            if i == len(ROWS) - 1:
                ax.set_xlabel(r"Resting shear rate $\dot\gamma_0$ (s$^{-1}$)")
            gammas = sorted({x for (sc, _, _), d in data.items() if sc == scale for x in d})
            ax.set_xticks(gammas)
            ax.set_xticklabels(["%g" % x for x in gammas])
            ax.xaxis.set_minor_formatter(NullFormatter())

    handles = [Line2D([], [], color=CLOSURES[c][1], ls=CLOSURES[c][2], lw=1.5, label=CLOSURES[c][0])
               for c in closures]
    handles += [Line2D([], [], color="k", ls="none", marker=m, mfc="white", label=n)
                for n, m in BLOODS.values()]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    out = a.out or os.path.join(a.validation, "fig_validation.pdf")
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(os.path.splitext(out)[0] + ".png", dpi=200, bbox_inches="tight")
    print("written", out, "and", os.path.splitext(out)[0] + ".png")

    # The numbers behind the figure, for the text.
    print("\n| R_1 | gamma0 | blood | " + " | ".join("%s rms / dE" % c for c in closures) + " |")
    print("|" + "---|" * (3 + len(closures)))
    for scale in scales:
        gammas = sorted({x for (sc, _, _), d in data.items() if sc == scale for x in d})
        for gm in gammas:
            for s in BLOODS:
                cells = []
                for c in closures:
                    r = data.get((scale, c, s), {}).get(gm)
                    cells.append("%.2g / %.2g" % (float(r["dQ_rms"]), float(r["dE"])) if r else "-")
                print("| %s | %g | %s | %s |" % (scale, gm, s, " | ".join(cells)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

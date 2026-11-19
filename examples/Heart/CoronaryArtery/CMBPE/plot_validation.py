#!/usr/bin/env python3
"""Figure of the truncation benchmark: error of the reduced closures against
the resolved solution, as a function of the resting operating point.

Reads <validation>/main/summary.csv (written by compare_rmu_axi.py) and, for
the Womersley numbers, the calibration lines of the full runs' log.txt.

Layout: one column per inlet radius (Womersley regime). Top row: rms error of
the inflow waveform (log axis). Bottom row: signed error of the dissipation in
the common region (linear axis; negative = the closure under-predicts it).
Hue = closure, shade (light to dark) and marker = blood; the shaded band spans
the three bloods of each closure.

  python3 plot_validation.py ~/Results/cmbpe/validation
  python3 plot_validation.py ~/Results/cmbpe/validation --scales 1e-4,1e-3
  python3 plot_validation.py ~/Results/cmbpe/validation --with-tree --with-plateau
"""
import argparse
import csv
import glob
import os
import re
import sys

BLOODS = [("healthy", "Healthy", "o"),
          ("diabetic", "Type 2 diabetes", "s"),
          ("hyperviscous", "Hyperviscous", "^")]

# closure -> (legend label, line style, shades for the three bloods, light to dark)
CLOSURES = {
    "rest": (r"constant $\mathcal{R}$", "--", ["#f2a65a", "#d1612a", "#8a3214"]),
    "single": (r"$\mathcal{R}_\mu$", "-", ["#7fb2d9", "#2f6fae", "#14375e"]),
    "tree": (r"$\mathcal{R}_\mu$, all generations", ":", ["#86c7a5", "#2e8b62", "#14482f"]),
    "plateau": (r"constant $\mathcal{R}$ at $\mu_\infty$", "-.", ["#bdbdbd", "#7d7d7d", "#3f3f3f"]),
}
NEUTRAL = ["#b0b0b0", "#707070", "#303030"]  # legend shades of the bloods
WAVEFORM_FLOOR = 1e-2                        # per cent, for the log axis


def womersley(validation, scale):
    """Range of the Womersley number of the narrowest tube for one radius."""
    wo = []
    for log in glob.glob(os.path.join(validation, "main", "s%s_g*" % scale, "*", "full", "log.txt")):
        found = re.findall(r"tube R=([0-9.eE+-]+) L=[0-9.eE+-]+ .*?Wo=([0-9.eE+-]+)",
                           open(log, errors="ignore").read())
        if found:
            rmin = min(float(r) for r, _ in found)
            wo += [float(w) for r, w in found if float(r) == rmin]
    return (min(wo), max(wo)) if wo else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("validation", help="the validation output folder")
    ap.add_argument("--scales", default="", help="comma-separated inlet radii to show (default all)")
    ap.add_argument("--with-plateau", action="store_true", help="add the mu_inf constant")
    ap.add_argument("--with-tree", action="store_true", help="add the exact multi-tube R_mu")
    ap.add_argument("--no-band", action="store_true", help="no shaded band across the bloods")
    ap.add_argument("--width", type=float, default=7.2, help="figure width (in)")
    ap.add_argument("--out", default=None, help="figure file (default <validation>/fig_validation.pdf)")
    a = ap.parse_args()

    path = os.path.join(a.validation, "main", "summary.csv")
    if not os.path.isfile(path):
        sys.exit("%s not found: run compare_rmu_axi.py first" % path)
    rows = list(csv.DictReader(open(path)))

    closures = ["rest", "single"] + (["tree"] if a.with_tree else []) + \
               (["plateau"] if a.with_plateau else [])

    data = {}  # (scale, closure, blood) -> {gamma: row}
    scales = set()
    for r in rows:
        m = re.match(r"s([0-9.eE+-]+)_g([0-9.eE+-]+)$", r["group"])
        if not m:
            continue
        scales.add(m.group(1))
        data.setdefault((m.group(1), r["closure"], r["state"]), {})[float(m.group(2))] = r
    if a.scales:
        wanted = [s.strip() for s in a.scales.split(",") if s.strip()]
        scales = [s for s in scales if any(float(s) == float(w) for w in wanted)]
    scales = sorted(scales, key=float)
    if not scales:
        sys.exit("no data for the requested radii")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullFormatter

    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
        "axes.labelsize": 8.5, "axes.titlesize": 8.5, "legend.fontsize": 8,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.direction": "out", "ytick.direction": "out",
    })

    ncol = len(scales)
    fig, axes = plt.subplots(2, ncol, figsize=(a.width, 4.7), sharex=True, sharey="row",
                             squeeze=False)

    def series(scale, closure, blood, key):
        d = data.get((scale, closure, blood))
        if not d:
            return [], []
        g = sorted(d)
        return g, [float(d[x][key]) for x in g]

    for j, scale in enumerate(scales):
        top, bottom = axes[0][j], axes[1][j]
        wo = womersley(a.validation, scale)
        title = r"$R_1 = %g$ mm" % (1e3 * float(scale))
        if wo:
            title += (r",  $\mathrm{Wo} \approx %.2g$" % wo[0] if abs(wo[1] / wo[0] - 1) < 0.15
                      else r",  $\mathrm{Wo} = %.2g$–$%.2g$" % wo)
        top.set_title(title, pad=6)

        for c in closures:
            label, style, shades = CLOSURES[c]
            for ax, key, signed in ((top, "dQ_rms", False), (bottom, "dE", True)):
                curves = []
                for k, (blood, _, marker) in enumerate(BLOODS):
                    g, y = series(scale, c, blood, key)
                    if not g:
                        continue
                    y = y if signed else [max(abs(v), WAVEFORM_FLOOR) for v in y]
                    curves.append((g, y))
                    ax.plot(g, y, ls=style, color=shades[k], lw=1.4, marker=marker, ms=4.2,
                            mfc=shades[k] if c != "rest" else "white", mec=shades[k], mew=1.0,
                            zorder=3, clip_on=False)
                if curves and not a.no_band and len(curves) > 1:
                    common = sorted(set.intersection(*[set(g) for g, _ in curves]))
                    lo = [min(y[g.index(x)] for g, y in curves) for x in common]
                    hi = [max(y[g.index(x)] for g, y in curves) for x in common]
                    ax.fill_between(common, lo, hi, color=shades[1], alpha=0.10, lw=0, zorder=1)

        top.set_yscale("log")
        bottom.axhline(0.0, color="#404040", lw=0.6, zorder=2)

        gammas = sorted({x for (sc, _, _), d in data.items() if sc == scale for x in d})
        for ax in (top, bottom):
            ax.set_xscale("log")
            ax.set_xticks(gammas)
            ax.set_xticklabels(["%g" % x for x in gammas])
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.tick_params(axis="x", which="minor", length=0)
            ax.grid(True, axis="y", which="major", lw=0.4, color="#d9d9d9", zorder=0)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
        bottom.set_xlabel(r"Resting shear rate $\dot\gamma_0$ (s$^{-1}$)")
        if j == 0:
            top.set_ylabel("Inflow waveform error (%)")
            bottom.set_ylabel("Dissipation error (%)")

    # Symmetric limits on the signed row make over- and under-prediction comparable.
    lim = max(abs(v) for ax in axes[1] for line in ax.get_lines() for v in line.get_ydata())
    for ax in axes[1]:
        ax.set_ylim(-1.08 * lim, 1.08 * lim)

    handles = [Line2D([], [], color=CLOSURES[c][2][1], ls=CLOSURES[c][1], lw=1.6,
                      label=CLOSURES[c][0]) for c in closures]
    handles += [Line2D([], [], color=NEUTRAL[k], ls="none", marker=m, ms=4.5, label=name)
                for k, (_, name, m) in enumerate(BLOODS)]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               bbox_to_anchor=(0.5, 0.0), handlelength=2.4, columnspacing=1.6)
    fig.tight_layout(rect=(0, 0.06, 1, 1), h_pad=1.2, w_pad=1.0)

    out = a.out or os.path.join(a.validation, "fig_validation.pdf")
    fig.savefig(out, bbox_inches="tight")
    png = os.path.splitext(out)[0] + ".png"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    print("written", out, "and", png)

    # The numbers behind the figure, for the text.
    print("\n| R_1 | gamma0 | blood | " + " | ".join("%s rms / dE" % c for c in closures) + " |")
    print("|" + "---|" * (3 + len(closures)))
    for scale in scales:
        gammas = sorted({x for (sc, _, _), d in data.items() if sc == scale for x in d})
        for gm in gammas:
            for blood, _, _ in BLOODS:
                cells = []
                for c in closures:
                    r = data.get((scale, c, blood), {}).get(gm)
                    cells.append("%.2g / %+.2g" % (float(r["dQ_rms"]), float(r["dE"])) if r else "-")
                print("| %s | %g | %s | %s |" % (scale, gm, blood, " | ".join(cells)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

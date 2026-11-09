#!/usr/bin/env python3
"""
Residence of blood in the aneurysm sac from the Pathlines3D output of every
A200 case of the campaign: statistics, probabilities, and publication figures
(ParaView renders of the pathlines composed with matplotlib).

    python3 pathlines_figures.py <root, e.g. ~/Results/arteries_v2/3D> [--out DIR]
        [--pvpython "pvpython --force-offscreen-rendering"] [--no-render]
        [--stride 2] [--boot 400]

Reads <root>/<artery>/<regime>/A200/<N|CY|R>/pathlines.csv (one row per
trajectory) and the period from pathlines_lesion.csv (or run.log).

Probabilities and statistics (times in cycles T, so that N, with its own T,
compares with CY and R):
  P_in        probability that a trajectory released at the inlet enters the
              sac (r > R inside the lesion slab), over all release phases;
  P_in(phi)   the same against the release phase of the cycle;
  S(tau)      washout survival: the probability that a trajectory which
              entered the sac has not left it for good tau after its first
              entry; Kaplan-Meier estimator (the trajectories stopped while
              still able to return are right-censored), Greenwood 95 % band;
  t_1/2       median of S (the half-life of the sac);
  S(1), S(2)  probability of staying more than one / two cycles;
  p(t_in)     probability density of the time spent inside the sac;
  t_in,50/90  its median and 90th percentile.
Each statistic gets a bootstrap 95 % interval (trajectories resampled), and
each difference CY - N, R - CY, R - N is marked significant when its
interval excludes zero. The intervals measure the sampling of trajectories,
not the cycle-to-cycle variability of the flow (one cycle of snapshots).

Outputs in --out (default <root>/../pathlines_figures):
  sac_residence.csv / .tex          the statistics per case
  sac_residence_diff.csv / .tex     the differences between fluids
  fig_washout.pdf/.png              S(tau), one panel per artery
  fig_sac_time_pdf.pdf/.png         p(t_in), one panel per artery
  fig_entry_phase.pdf/.png          P_in(phi), one panel per artery
  fig_pathlines_<artery>.pdf/.png   N, CY and R pathlines in the sac,
                                    coloured by the time inside (ParaView)

The render runs this same file under pvpython (--render-case); everything
else needs only numpy and matplotlib.
"""
import argparse
import csv
import glob
import json
import os
import re
import shlex
import subprocess
import sys

FLUIDS = ("N", "CY", "R")
PAIRS = (("CY", "N"), ("R", "CY"), ("R", "N"))
ARTERY_LABEL = {"aorta": "abdominal aorta", "carotid": "carotid", "coronary": "coronary",
                "femoral": "femoral", "small": "radial"}
ARTERY_ORDER = ("coronary", "carotid", "femoral", "aorta", "small")


# =========================================================================================
# ParaView render of one case (run under pvpython)
# =========================================================================================
def render_case(spec):
    """spec: case, out, vmin, vmax, stride, size. Half wall behind, pathlines
    that entered the sac coloured by the time inside it (log scale)."""
    import numpy as np
    from paraview import servermanager as sm
    from paraview.simple import (OpenDataFile, MergeBlocks, CleantoGrid, ExtractSurface, Calculator, Threshold,
                                 PythonCalculator, CreateView, Show, ColorBy, GetColorTransferFunction,
                                 SaveScreenshot, Render)
    from vtk.numpy_interface import dataset_adapter as dsa

    def threshold(src, array, lo, hi, assoc="POINTS"):
        th = Threshold(Input=src)
        th.Scalars = [assoc, array]
        try:                                    # ParaView >= 5.10
            th.LowerThreshold, th.UpperThreshold = lo, hi
            th.ThresholdMethod = "Between"
        except AttributeError:                  # older
            th.ThresholdRange = [lo, hi]
        return th

    case = spec["case"]
    # -- the wall: surface of the mesh, the lesion found from the radius
    # The rank pieces share their interface vertices: merge them, or the
    # partition interfaces would appear as surfaces inside the vessel.
    vol = CleantoGrid(Input=MergeBlocks(Input=OpenDataFile(os.path.join(case, "case.xdmf"))))
    surf = ExtractSurface(Input=vol)
    X = np.asarray(dsa.WrapDataObject(sm.Fetch(surf)).Points)
    ext = X.max(axis=0) - X.min(axis=0)
    a = int(np.argmax(ext))                     # the axis
    b, c = [k for k in range(3) if k != a]      # view along +b, up c
    yc, zc = np.median(X[:, b]), np.median(X[:, c])
    r = np.hypot(X[:, b] - yc, X[:, c] - zc)
    wall = r > 0.9 * np.median(r)               # drop the cap interiors
    R = np.median(r[wall])
    les = wall & (np.abs(r - R) > 0.05 * R)
    xs = X[les, a] if les.any() else X[wall, a]
    xc = 0.5 * (xs.min() + xs.max())
    Rmax = r[wall].max()
    centre = [0.0, 0.0, 0.0]
    centre[a], centre[b], centre[c] = xc, yc, zc
    coord = ("coordsX", "coordsY", "coordsZ")[b]
    side = Calculator(Input=surf)
    side.ResultArrayName = "side"
    side.Function = coord
    back = threshold(side, "side", yc, X[:, b].max() + 1.0, "POINTS")   # the far half

    # -- the pathlines that entered the sac, every stride-th seed
    pl = MergeBlocks(Input=OpenDataFile(os.path.join(case, "pathlines.xdmf")))
    sac = threshold(pl, "sacTime_T", 1.0e-6, 1.0e9)
    keep = PythonCalculator(Input=sac)
    keep.ArrayName = "keep"
    keep.Expression = "numpy.mod(inputs[0].PointData['seed'], %d)" % max(1, int(spec["stride"]))
    sel = threshold(keep, "keep", -0.5, 0.5)

    W, H = spec["size"]
    view = CreateView("RenderView")
    view.ViewSize = [W, H]
    for prop, val in (("UseColorPaletteForBackground", 0), ("Background", [1.0, 1.0, 1.0]),
                      ("OrientationAxesVisibility", 0), ("UseFXAA", 1), ("EnableRenderOnInteraction", 0)):
        try:
            setattr(view, prop, val)
        except AttributeError:
            pass
    wd = Show(back, view)
    ColorBy(wd, None)
    wd.DiffuseColor = [0.82, 0.82, 0.84]
    wd.AmbientColor = [0.82, 0.82, 0.84]
    wd.Specular = 0.15
    ld = Show(sel, view)
    ColorBy(ld, ("POINTS", "sacTime_T"))
    lut = GetColorTransferFunction("sacTime_T")
    lut.ApplyPreset("Viridis (matplotlib)", True)
    try:
        lut.AutomaticRescaleRangeMode = "Never"
    except AttributeError:
        pass
    lut.RescaleTransferFunction(spec["vmin"], spec["vmax"])
    lut.UseLogScale = 1
    ld.SetScalarBarVisibility(view, False)
    ld.LineWidth = 1.6
    try:
        ld.RenderLinesAsTubes = 1
    except AttributeError:
        pass
    # -- camera: orthographic, normal to the axis, framing the sac
    pos = list(centre)
    pos[b] -= 40.0 * R
    up = [0.0, 0.0, 0.0]
    up[c] = 1.0
    view.CameraPosition = pos
    view.CameraFocalPoint = centre
    view.CameraViewUp = up
    view.CameraParallelProjection = 1
    view.CameraParallelScale = 1.12 * Rmax
    Render(view)
    SaveScreenshot(spec["out"], view, ImageResolution=[W, H])
    print("rendered %s" % spec["out"])


# =========================================================================================
# statistics
# =========================================================================================
def load_case(folder):
    import numpy as np
    period = None
    les = os.path.join(folder, "pathlines_lesion.csv")
    if os.path.exists(les):
        with open(les) as f:
            for row in csv.DictReader(f):
                period = float(row["period"])
                break
    if period is None:
        with open(os.path.join(folder, "run.log"), errors="replace") as f:
            m = re.search(r"\[inflow\].*\bT=([\d.e+-]+) s", f.read())
        period = float(m.group(1))
    cols = {}
    with open(os.path.join(folder, "pathlines.csv")) as f:
        rd = csv.DictReader(f)
        for row in rd:
            for k, v in row.items():
                cols.setdefault(k, []).append(v)
    if "sacTime" not in cols:
        raise RuntimeError("%s/pathlines.csv has no sacTime: rerun the new Pathlines3D" % folder)
    exitc = np.array(cols["exit"])
    ok = exitc != "wall"                        # a tracing limit, not a fate
    d = {
        "T": period,
        "release": np.array(cols["release"], int)[ok],
        "t_in": np.array(cols["sacTime"], float)[ok] / period,
        "span": np.array(cols["sacSpan"], float)[ok] / period,
        "event": (exitc == "outlet")[ok],       # left for good (else censored)
    }
    d["entered"] = d["t_in"] > 0.0
    d["releases"] = int(d["release"].max()) + 1 if len(d["release"]) else 1
    return d


def kaplan_meier(t, event):
    """Survival S(t) right after each sorted time, and Greenwood variance."""
    import numpy as np
    o = np.argsort(t, kind="stable")
    t, e = t[o], event[o].astype(float)
    n = len(t)
    risk = n - np.arange(n, dtype=float)
    S = np.cumprod(1.0 - e / risk)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(risk > 1.0, e / (risk * (risk - e)), 0.0)
    var = S ** 2 * np.cumsum(g)
    return t, S, var


def surv_at(t, S, tau):
    import numpy as np
    k = np.searchsorted(t, tau, side="right") - 1
    return 1.0 if k < 0 else float(S[k])


def km_median(t, S):
    import numpy as np
    k = np.nonzero(S <= 0.5)[0]
    return float(t[k[0]]) if len(k) else float("nan")


def stats(d, idx=None):
    """The statistics of one case (on the trajectories idx, for the bootstrap)."""
    import numpy as np
    if idx is None:
        idx = np.arange(len(d["t_in"]))
    ent = d["entered"][idx]
    out = {"P_in": float(ent.mean()) if len(ent) else float("nan")}
    sub = idx[ent]
    if len(sub) < 5:
        return out
    t, S, _ = kaplan_meier(d["span"][sub], d["event"][sub])
    ti = np.sort(d["t_in"][sub])
    out.update({
        "t_half": km_median(t, S), "S1": surv_at(t, S, 1.0), "S2": surv_at(t, S, 2.0),
        "t_in50": float(ti[len(ti) // 2]), "t_in90": float(ti[min(len(ti) - 1, int(0.9 * len(ti)))]),
        "n_in": int(len(sub)), "censored": float(1.0 - d["event"][sub].mean()),
    })
    return out


STAT_KEYS = ("P_in", "t_half", "S1", "S2", "t_in50", "t_in90")


def bootstrap(d, B, rng):
    import numpy as np
    n = len(d["t_in"])
    reps = {k: [] for k in STAT_KEYS}
    for _ in range(B):
        s = stats(d, rng.integers(0, n, n))
        for k in STAT_KEYS:
            reps[k].append(s.get(k, np.nan))
    return {k: np.array(v) for k, v in reps.items()}


# =========================================================================================
# figures (matplotlib, JFM style)
# =========================================================================================
STYLE = {"N": dict(color="0.55", ls=(0, (4, 2)), lw=1.1, label="N"),
         "CY": dict(color="#1f5fa8", ls="-", lw=1.1, label="CY"),
         "R": dict(color="#b2182b", ls="-", lw=1.4, label="R")}


def jfm_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.linewidth": 0.6, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "legend.frameon": False,
        "savefig.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    })


def panel_grid(n):
    import matplotlib.pyplot as plt
    ncol = min(n, 2)
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2, 2.0 * nrow + 0.2), squeeze=False)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    return fig, list(axes.flat[:n])


def shared_legend(fig, ax):
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=len(l), bbox_to_anchor=(0.5, 1.0), handlelength=2.4)


def label(ax, k, text):
    ax.text(0.03, 0.95, "(%s) %s" % ("abcdefgh"[k], text), transform=ax.transAxes, va="top", ha="left")


def save(fig, out, name):
    import matplotlib.pyplot as plt
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, "%s.%s" % (name, ext)))
    plt.close(fig)


def fig_washout(data, out):
    import numpy as np
    fig, axes = panel_grid(len(data))
    for k, (art, cases) in enumerate(data):
        ax = axes[k]
        for f in FLUIDS:
            if f not in cases:
                continue
            d = cases[f]
            m = d["entered"]
            t, S, var = kaplan_meier(d["span"][m], d["event"][m])
            tt = np.concatenate([[0.0], t])
            SS = np.concatenate([[1.0], S])
            sd = np.sqrt(np.concatenate([[0.0], var]))
            ax.fill_between(tt, np.clip(SS - 1.96 * sd, 1e-4, 1), np.clip(SS + 1.96 * sd, 1e-4, 1),
                            step="post", color=STYLE[f]["color"], alpha=0.15, lw=0)
            ax.step(tt, SS, where="post", **STYLE[f])
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 1.05)
        tmax = np.percentile(np.concatenate([c["span"][c["entered"]] for c in cases.values()]), 99.5)
        ax.set_xlim(0, max(tmax, 1.0))
        ax.axhline(0.5, color="0.8", lw=0.5, zorder=0)
        ax.set_xlabel(r"$\tau/T$")
        ax.set_ylabel(r"$S(\tau)$")
        label(ax, k, ARTERY_LABEL.get(art, art))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    shared_legend(fig, axes[0])
    save(fig, out, "fig_washout")


def fig_time_pdf(data, out):
    import numpy as np
    fig, axes = panel_grid(len(data))
    for k, (art, cases) in enumerate(data):
        ax = axes[k]
        allt = np.concatenate([c["t_in"][c["entered"]] for c in cases.values()])
        if not len(allt):
            continue
        edges = np.logspace(np.log10(max(allt.min(), 1e-3)), np.log10(allt.max() * 1.01), 30)
        for f in FLUIDS:
            if f not in cases:
                continue
            ti = cases[f]["t_in"][cases[f]["entered"]]
            h, _ = np.histogram(ti, bins=edges, density=True)
            xc = np.sqrt(edges[1:] * edges[:-1])
            ax.plot(xc, h * xc, **STYLE[f])          # t p(t): area-true on a log axis
        ax.set_xscale("log")
        ax.set_xlabel(r"$t_{\mathrm{in}}/T$")
        ax.set_ylabel(r"$t_{\mathrm{in}}\,p(t_{\mathrm{in}})$")
        label(ax, k, ARTERY_LABEL.get(art, art))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    shared_legend(fig, axes[0])
    save(fig, out, "fig_sac_time_pdf")


def fig_entry_phase(data, out):
    import numpy as np
    fig, axes = panel_grid(len(data))
    for k, (art, cases) in enumerate(data):
        ax = axes[k]
        for f in FLUIDS:
            if f not in cases:
                continue
            d = cases[f]
            M = d["releases"]
            ph = np.arange(M) / M
            p = np.array([d["entered"][d["release"] == m].mean() if (d["release"] == m).any() else np.nan
                          for m in range(M)])
            n = np.array([(d["release"] == m).sum() for m in range(M)])
            err = 1.96 * np.sqrt(p * (1 - p) / np.maximum(n, 1))
            ph = np.append(ph, 1.0); p = np.append(p, p[0]); err = np.append(err, err[0])
            ax.errorbar(ph, p, yerr=err, marker="o", ms=2.5, capsize=1.5, elinewidth=0.6, **STYLE[f])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, None)
        ax.set_xlabel(r"release phase $t_0/T$")
        ax.set_ylabel(r"$P_{\mathrm{in}}$")
        label(ax, k, ARTERY_LABEL.get(art, art))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    shared_legend(fig, axes[0])
    save(fig, out, "fig_entry_phase")


def crop_white(img, pad=6):
    import numpy as np
    rgb = img[..., :3]
    mask = (rgb < 0.98).any(axis=2)
    if not mask.any():
        return img
    ys, xs = np.nonzero(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, img.shape[1])
    return img[y0:y1, x0:x1]


def fig_pathlines(art, pngs, vmin, vmax, out):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable
    imgs = [(f, crop_white(plt.imread(p))) for f, p in pngs]
    h = sum(im.shape[0] / im.shape[1] for _, im in imgs)
    fig = plt.figure(figsize=(5.2, 5.2 * h + 0.45))
    gs = fig.add_gridspec(len(imgs) + 1, 1, height_ratios=[im.shape[0] / im.shape[1] for _, im in imgs] + [0.05],
                          hspace=0.16)
    for k, (f, im) in enumerate(imgs):
        ax = fig.add_subplot(gs[k])
        ax.imshow(im, interpolation="lanczos")
        ax.set_axis_off()
        ax.text(0.0, 1.01, "(%s) %s" % ("abc"[k], f), transform=ax.transAxes, va="bottom", ha="left")
    cax = fig.add_subplot(gs[-1])
    cb = fig.colorbar(ScalarMappable(LogNorm(vmin, vmax), "viridis"), cax=cax, orientation="horizontal")
    cb.set_label(r"time inside the sac, $t_{\mathrm{in}}/T$")
    cb.outline.set_linewidth(0.5)
    save(fig, out, "fig_pathlines_%s" % art)


# =========================================================================================
# tables
# =========================================================================================
def fmt_ci(v, lo, hi, nd=2):
    import numpy as np
    if not np.isfinite(v):
        return "--"
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return ("%%.%df" % nd) % v
    return ("%%.%df\\,[%%.%df, %%.%df]" % (nd, nd, nd)) % (v, lo, hi)


def tables(rows, diffs, out):
    keys = ["artery", "fluid", "n_in", "censored"] + [x for k in STAT_KEYS for x in (k, k + "_lo", k + "_hi")]
    with open(os.path.join(out, "sac_residence.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    dk = ["artery", "pair"] + [x for k in STAT_KEYS for x in (k, k + "_lo", k + "_hi", k + "_sig")]
    with open(os.path.join(out, "sac_residence_diff.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=dk, extrasaction="ignore")
        w.writeheader()
        w.writerows(diffs)
    head = (r"$P_{\mathrm{in}}$", r"$t_{1/2}/T$", r"$S(T)$", r"$S(2T)$", r"$t_{\mathrm{in},50}/T$", r"$t_{\mathrm{in},90}/T$")
    with open(os.path.join(out, "sac_residence.tex"), "w") as f:
        f.write("\\begin{tabular}{ll" + "c" * len(STAT_KEYS) + "}\n\\hline\nartery & fluid & " + " & ".join(head) + "\\\\\n\\hline\n")
        for r in rows:
            f.write("%s & %s & %s\\\\\n" % (ARTERY_LABEL.get(r["artery"], r["artery"]), r["fluid"], " & ".join(
                fmt_ci(r.get(k, float("nan")), r.get(k + "_lo", float("nan")), r.get(k + "_hi", float("nan")))
                for k in STAT_KEYS)))
        f.write("\\hline\n\\end{tabular}\n% value [bootstrap 95 % interval over trajectories]\n")
    with open(os.path.join(out, "sac_residence_diff.tex"), "w") as f:
        f.write("\\begin{tabular}{ll" + "c" * len(STAT_KEYS) + "}\n\\hline\nartery & pair & " + " & ".join(
            r"$\Delta$" + h for h in head) + "\\\\\n\\hline\n")
        for r in diffs:
            cells = []
            for k in STAT_KEYS:
                c = fmt_ci(r.get(k, float("nan")), r.get(k + "_lo", float("nan")), r.get(k + "_hi", float("nan")))
                cells.append(c + ("$^{*}$" if r.get(k + "_sig") else ""))
            f.write("%s & %s & %s\\\\\n" % (ARTERY_LABEL.get(r["artery"], r["artery"]), r["pair"].replace("-", "$-$"),
                                            " & ".join(cells)))
        f.write("\\hline\n\\end{tabular}\n% $^{*}$: the bootstrap 95 % interval excludes zero\n")


# =========================================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", nargs="?")
    ap.add_argument("--out")
    ap.add_argument("--pvpython", default=os.environ.get("PVPYTHON", "pvpython --force-offscreen-rendering"))
    ap.add_argument("--no-render", action="store_true")
    ap.add_argument("--stride", type=int, default=2, help="draw every stride-th seed")
    ap.add_argument("--boot", type=int, default=400)
    ap.add_argument("--size", type=int, nargs=2, default=[2400, 1000])
    ap.add_argument("--render-case", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.render_case:
        render_case(json.loads(args.render_case))
        return
    if not args.root:
        ap.error("root is required")
    import numpy as np
    jfm_style()
    root = os.path.abspath(os.path.expanduser(args.root))
    out = os.path.abspath(args.out or os.path.join(root, "..", "pathlines_figures"))
    os.makedirs(out, exist_ok=True)
    found = {}
    for p in sorted(glob.glob(os.path.join(root, "*", "*", "A200", "*", "pathlines.csv"))):
        folder = os.path.dirname(p)
        fluid = os.path.basename(folder)
        art = os.path.relpath(folder, root).split(os.sep)[0]
        if fluid in FLUIDS:
            found.setdefault(art, {})[fluid] = folder
    if not found:
        sys.exit("no A200 case with pathlines.csv under %s" % root)
    rng = np.random.default_rng(1)
    data, rows, diffs = [], [], []
    for art in sorted(found, key=lambda a: ARTERY_ORDER.index(a) if a in ARTERY_ORDER else 99):
        cases, boots = {}, {}
        for f in FLUIDS:
            if f not in found[art]:
                continue
            try:
                cases[f] = load_case(found[art][f])
            except Exception as e:
                print("  %s/%s: skipped (%s)" % (art, f, e))
                continue
            s = stats(cases[f])
            boots[f] = bootstrap(cases[f], args.boot, rng)
            row = {"artery": art, "fluid": f, "n_in": s.get("n_in", 0), "censored": s.get("censored", np.nan)}
            for k in STAT_KEYS:
                row[k] = s.get(k, np.nan)
                row[k + "_lo"], row[k + "_hi"] = np.nanpercentile(boots[f][k], [2.5, 97.5]) \
                    if np.isfinite(boots[f][k]).any() else (np.nan, np.nan)
            rows.append(row)
            print("  %-9s %-2s P_in %.3f  t1/2 %.3f T  S(1T) %.3f  S(2T) %.3f  t_in50 %.3f T  t_in90 %.3f T  (n_in %d, censored %.1f %%)"
                  % (art, f, row["P_in"], row["t_half"], row["S1"], row["S2"], row["t_in50"], row["t_in90"],
                     row["n_in"], 100 * row["censored"]))
        for a, b in PAIRS:
            if a in cases and b in cases:
                r = {"artery": art, "pair": "%s-%s" % (a, b)}
                for k in STAT_KEYS:
                    sa, sb = stats(cases[a]).get(k, np.nan), stats(cases[b]).get(k, np.nan)
                    dd = boots[a][k] - boots[b][k]
                    lo, hi = np.nanpercentile(dd, [2.5, 97.5]) if np.isfinite(dd).any() else (np.nan, np.nan)
                    r.update({k: sa - sb, k + "_lo": lo, k + "_hi": hi, k + "_sig": bool(lo > 0 or hi < 0)})
                diffs.append(r)
        if cases:
            data.append((art, cases))
    tables(rows, diffs, out)
    fig_washout(data, out)
    fig_time_pdf(data, out)
    fig_entry_phase(data, out)
    print("statistics and figures in %s" % out)
    if args.no_render:
        return
    # ParaView renders, one per case, with a common colour range per artery.
    pv = shlex.split(args.pvpython)
    for art, cases in data:
        allt = np.concatenate([c["t_in"][c["entered"]] for c in cases.values()])
        if not len(allt):
            continue
        vmin = max(float(np.percentile(allt, 5)), 1e-3)
        vmax = float(np.percentile(allt, 99))
        pngs = []
        for f in FLUIDS:
            if f not in cases:
                continue
            png = os.path.join(out, "render_%s_%s.png" % (art, f))
            spec = {"case": found[art][f], "out": png, "vmin": vmin, "vmax": vmax,
                    "stride": args.stride, "size": args.size}
            cmd = pv + [os.path.abspath(__file__), "--render-case", json.dumps(spec)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0 or not os.path.exists(png):
                print("  render %s/%s failed:\n%s" % (art, f, (res.stderr or res.stdout)[-2000:]))
                continue
            pngs.append((f, png))
        if pngs:
            fig_pathlines(art, pngs, vmin, vmax, out)
            print("  fig_pathlines_%s" % art)


if __name__ == "__main__":
    main()

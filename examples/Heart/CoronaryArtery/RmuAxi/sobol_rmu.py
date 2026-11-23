#!/usr/bin/env python3
"""Sobol sensitivity of the truncation error of an outlet resistance.

The quasi-steady counterpart of the RmuAxi benchmark (RmuAxi.cpp): the stepped
pipe of tubes R:L = 1:10, 0.6:30, 0.35:30 (units of the inlet radius) is cut at
x_T = 13 and driven by p_in(t) = Pbar [1 + A sin(2 pi t / T)]. Every segment is
a WRMS tube element, Delta p_k = R_mu,k(Q) Q, so that the 0D model of the full
pipe is the series of all segments, and the truncated models replace the
segments downstream of x_T by

  rest     the constant R_0 = R_mu(Q_0) of the removed segments (Closure::Rest);
  single   one WRMS compartment of radius 0.35, its length matched at Q_0
           (Closure::Single, the R_mu of the paper).

The calibration is that of Benchmark::calibrate(): Q_0 = pi R_min^3 gd_0 / 4 and
Pbar = Q_0 times the series resistance of the full pipe at Q_0. In the
quasi-steady regime (Wo << 1) the series model reproduces the resolved
benchmark to within 0.6 per cent, so its errors are those of the closures.

Pressure, flow and dissipation are scaled by Pbar, Q_0 and Pbar Q_0. The
constant model is then Delta p / Pbar = Q / Q_0, and every WRMS element is
Q / Q_0 times f(Q / Q_0) with f(s) = mu_eff(gd_0 s) / mu_eff(gd_0): the inlet
radius, the lengths and the multiplicity drop out, and since mu_eff of the
Carreau-Yasuda law is mu_inf times a function of lambda gd_a, f depends only on
the Carreau number Cu = lambda gd_0, the viscosity ratio mu_0 / mu_inf, n and a.
The errors are those of compare_rmu_axi.py,

  dQ_rms  100 rms(Q_c - Q_f) / mean |Q_f|                  inflow waveform
  dE      100 (E_c - E_f) / E_f,  E = mean(Q Delta p_ret)   dissipation, x <= x_T

and their Sobol indices (first order S1, total ST; Saltelli sampling, SALib)
are computed over these five dimensionless groups and the pulsatility A
(BOUNDS). A deterministic map of dQ_rms over (gd_0, A) for the three fitted
bloods locates where the constant closure fails, together with resting shear
rates and pulsatilities measured in microvessels (LANDMARKS).

  pip install SALib
  python3 sobol_rmu.py --out sobol                      # indices, map and figures
  python3 sobol_rmu.py --out sobol --check summary.csv  # 0D vs RmuAxi runs
"""
import argparse
import csv
import math
import os
import sys

import numpy as np

Pi = math.pi

# Ranges of the Sobol analysis over the dimensionless groups that the errors
# depend on: (name, lower, upper, log-uniform). They cover the three fits of
# Table 1 at resting shear rates of 30-800 1/s (Cu = 280-14000,
# mu_0/mu_inf = 17-65).
BOUNDS = [
    ("carreau", 50.0, 2.0e4, True),     # Cu = lambda gd_0
    ("ratio", 10.0, 100.0, True),       # mu_0 / mu_inf
    ("n", 0.20, 0.40, False),
    ("a", 0.60, 2.00, False),
    ("amplitude", 0.20, 1.20, False),   # A
]
LABELS = {"carreau": r"$\lambda\dot\gamma_0$", "ratio": r"$\mu_0/\mu_\infty$",
          "n": r"$n$", "a": r"$a$", "amplitude": r"$A$"}
# The errors are invariant under (mu_0, mu_inf) -> k (mu_0, mu_inf) and
# (lambda, gd_0) -> (lambda / c, c gd_0); these representatives are used.
MuInfRef, GammaRef = 0.004, 100.0

# Table 1 of the paper (SI): mu0, muInf, lambda, n, a.
STATES = {
    "healthy": (0.0712616, 0.0041, 9.2780, 0.3572, 1.8773),
    "diabetic": (0.0973528, 0.0049622, 15.7825, 0.2670, 1.0829),
    "hyperviscous": (0.3467269, 0.0053264, 17.4172, 0.2200, 0.6957),
}


class Tube:
    """One tube of the series: radius and length, in units of the inlet radius."""

    def __init__(self, radius, length):
        self.radius = radius
        self.length = length


def splitTubes(tubes, truncation):
    """Retained (x <= x_T) and removed tubes, as run_rmu_axi.sh cuts them."""
    kept, removed, x = [], [], 0.0
    for t in tubes:
        if x < truncation - 1e-12:
            kept.append(Tube(t.radius, min(t.length, truncation - x)))
        if x + t.length > truncation + 1e-12:
            removed.append(Tube(t.radius, t.length - max(0.0, truncation - x)))
        x += t.length
    return kept, removed


def rowInterp(x, xp, fp):
    """Row-wise linear interpolation, clamped to the end values.

    x has shape (S, m); xp, fp have shape (S, n) with xp increasing in every
    row. One searchsorted on the rows laid end to end, each shifted by an
    offset larger than the range of the data.
    """
    S, n = xp.shape
    lo = min(xp.min(), x.min())
    span = max(xp.max(), x.max()) - lo + 1.0
    off = (np.arange(S) * span)[:, None]
    j = np.searchsorted((xp - lo + off).ravel(), (x - lo + off).ravel())
    i = np.clip(j.reshape(x.shape) - (np.arange(S) * n)[:, None], 1, n - 1)
    r = np.arange(S)[:, None]
    x0, x1 = xp[r, i - 1], xp[r, i]
    w = np.clip((x - x0) / np.where(x1 > x0, x1 - x0, 1.0), 0.0, 1.0)
    return fp[r, i - 1] + w * (fp[r, i] - fp[r, i - 1])


class HydraulicCurve:
    """Batched counterpart of HydraulicCurve in RmuAxi.cpp.

    One effective (hydraulic) viscosity curve mu_eff(gd_a) per parameter set,
    with mu_eff = tau_w^4 / (4 I) and gd_a = tau_w / mu_eff. The C++ class
    tabulates I on nodes in tau_w and inverts tau = mu(gd) gd by safeguarded
    Newton; here the nodes are placed in the wall shear rate gd_w, which needs
    no inversion since tau(gd) is increasing, and
    I = int_0^{gd_w} gd^3 mu^2 (mu + gd mu') d gd is the same integral (its
    integrand is that of the C++ Simpson rule), accumulated by the trapezoidal
    rule in ln gd. Outside the table the plateaus are the continuation, hence
    clamping, as in operator().
    """

    TableNodes = 2001          # quadrature error of mu_eff below 0.03 per cent
    GammaMin = 1.0e-6
    GammaMax = 1.0e6

    def __init__(self, mu0, muInf, lam, n, a):
        self.mu0, self.muInf, self.lam, self.n, self.a = (
            np.asarray(v, dtype=float)[:, None] for v in (mu0, muInf, lam, n, a))
        lg = np.linspace(math.log(self.GammaMin), math.log(self.GammaMax), self.TableNodes)
        g = np.exp(lg)[None, :]
        mu = self.viscosity(g)
        f = g ** 4 * mu ** 2 * (mu + g * self.derivative(g))      # integrand x gd
        h = lg[1] - lg[0]
        I0 = self.mu0 ** 3 * self.GammaMin ** 4 / 4.0              # mu = mu0 below gd_min
        I = I0 + np.concatenate([np.zeros_like(I0),
                                 np.cumsum(0.5 * h * (f[:, 1:] + f[:, :-1]), axis=1)], axis=1)
        tau = mu * g
        muAp = tau ** 4 / (4.0 * I)
        self.m_logGamma = np.maximum.accumulate(np.log(tau / muAp), axis=1)
        self.m_logMu = np.log(muAp)

    def viscosity(self, g):
        """Local Carreau-Yasuda viscosity mu(gd)."""
        return self.muInf + (self.mu0 - self.muInf) * \
            (1.0 + (self.lam * np.maximum(g, 0.0)) ** self.a) ** ((self.n - 1.0) / self.a)

    def derivative(self, g):
        """d mu / d gd of the local law."""
        base = 1.0 + (self.lam * g) ** self.a
        return (self.mu0 - self.muInf) * (self.n - 1.0) * \
            base ** ((self.n - 1.0 - self.a) / self.a) * self.lam ** self.a * g ** (self.a - 1.0)

    def __call__(self, gA):
        """Hydraulic (effective) viscosity at the apparent shear rate, one row per set."""
        return np.exp(rowInterp(np.log(np.maximum(gA, 1e-300)), self.m_logGamma, self.m_logMu))

    def resistance(self, R, L, Q):
        """Resistance of one tube (radius R, length L) at flow Q, one row per set."""
        gA = 4.0 * np.abs(Q) / (Pi * R ** 3)
        return 8.0 * self(gA) * L / (Pi * R ** 4)


class Benchmark:
    """The quasi-steady RmuAxi benchmark for a batch of parameter sets."""

    FullTubes = [Tube(1.0, 10.0), Tube(0.6, 30.0), Tube(0.35, 30.0)]
    Truncation = 13.0
    TimeSamples = 400
    FlowNodes = 801            # nodes of Q / Q_0 on which every model is tabulated

    def __init__(self, curve, gammaRef, amplitude):
        self.m_curve = curve
        self.m_amplitude = np.asarray(amplitude, dtype=float)[:, None]
        self.m_kept, self.m_removed = splitTubes(self.FullTubes, self.Truncation)
        rMin = min(t.radius for t in self.FullTubes)
        self.m_singleRadius = self.m_removed[-1].radius
        S = len(self.m_amplitude)
        self.m_q0 = (Pi * rMin ** 3 * np.asarray(gammaRef, dtype=float) / 4.0)[:, None]
        self.m_pbar = self.m_q0 * self.seriesResistance(self.FullTubes, self.m_q0)
        self.m_restResistance = self.seriesResistance(self.m_removed, self.m_q0)
        Rs = self.m_singleRadius
        muS = curve(4.0 * self.m_q0 / (Pi * Rs ** 3))
        self.m_singleLength = self.m_restResistance * Pi * Rs ** 4 / (8.0 * muS)
        self.m_s = np.concatenate([[0.0], np.logspace(-4, 3, self.FlowNodes - 1)])[None, :]
        self.m_s = np.repeat(self.m_s, S, axis=0)

    def seriesResistance(self, tubes, Q):
        return sum(self.m_curve.resistance(t.radius, t.length, Q) for t in tubes)

    def pressures(self, closure):
        """Delta p / Pbar of the retained part and of the whole model on m_s."""
        Q = self.m_s * self.m_q0
        kept = self.seriesResistance(self.m_kept, Q) * Q / self.m_pbar
        if closure == "full":
            down = self.seriesResistance(self.m_removed, Q) * Q
        elif closure == "rest":
            down = self.m_restResistance * Q
        elif closure == "single":
            down = self.m_curve.resistance(self.m_singleRadius, self.m_singleLength, Q) * Q
        else:
            raise ValueError(closure)
        return kept, kept + down / self.m_pbar

    def solve(self, closure):
        """Q(t)/Q_0 over one cycle and the cycle-mean dissipation of x <= x_T."""
        kept, total = self.pressures(closure)
        t = (np.arange(self.TimeSamples) + 0.5) / self.TimeSamples
        p = 1.0 + self.m_amplitude * np.sin(2.0 * Pi * t)[None, :]
        # every model is odd and strictly increasing in Q (Lemma 1)
        s = np.sign(p) * rowInterp(np.abs(p), total, self.m_s)
        pk = np.sign(s) * rowInterp(np.abs(s), self.m_s, kept)
        return s, (s * pk).mean(axis=1)

    def errors(self, closure, reference=None):
        sf, ef = reference or self.solve("full")
        sc, ec = self.solve(closure)
        dQ = 100.0 * np.sqrt(((sc - sf) ** 2).mean(axis=1)) / np.abs(sf).mean(axis=1)
        dE = 100.0 * (ec - ef) / ef
        return dQ, dE


def evaluate(X):
    """dQ_rms and dE of the rest and single closures.

    Rows of X: mu0, muInf, lambda, n, a, gammaRef, amplitude (SI).
    """
    mu0, muInf, lam, n, a, gammaRef, amplitude = X.T
    curve = HydraulicCurve(mu0, muInf, lam, n, a)
    bench = Benchmark(curve, gammaRef, amplitude)
    ref = bench.solve("full")
    out = {}
    for c in ("rest", "single"):
        dQ, dE = bench.errors(c, ref)
        out[c + "_dQ_rms"], out[c + "_dE"] = dQ, dE
    return out


def physical(Z):
    """Rows of dimensionless groups (BOUNDS order) -> rows for evaluate()."""
    cu, ratio, n, a, amp = Z.T
    one = np.ones(len(Z))
    return np.column_stack([ratio * MuInfRef, MuInfRef * one, cu / GammaRef, n, a,
                            GammaRef * one, amp])


def evaluateChunked(X, chunk=2048):
    parts = [evaluate(X[i:i + chunk]) for i in range(0, len(X), chunk)]
    return {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}


def check(summaryPath, scale="1e-4"):
    """0D model against the RmuAxi runs of one radius (compare_rmu_axi.py output)."""
    rows = list(csv.DictReader(open(summaryPath)))
    print("\n0D quasi-steady model vs RmuAxi, R_1 = %s m (A = 1.2):" % scale)
    print("%-13s %6s %-7s %17s %17s" % ("blood", "gd_0", "closure", "dQ_rms 0D / CFD", "dE 0D / CFD"))
    for state, (mu0, mi, lam, n, a) in STATES.items():
        for g in (30.0, 100.0, 200.0, 400.0):
            X = np.array([[mu0, mi, lam, n, a, g, 1.2]])
            e = evaluate(X)
            for c in ("rest", "single"):
                cfd = [r for r in rows if r["group"] == "s%s_g%g" % (scale, g)
                       and r["state"] == state and r["closure"] == c]
                cq = "%.2f" % float(cfd[0]["dQ_rms"]) if cfd else "-"
                ce = "%+.2f" % float(cfd[0]["dE"]) if cfd else "-"
                print("%-13s %6g %-7s %8.2f / %-6s %+8.2f / %-6s" % (
                    state, g, c, e[c + "_dQ_rms"][0], cq, e[c + "_dE"][0], ce))


# Resting wall shear rates (1/s) and pulsatilities measured in microvessels,
# and the operating points of the coronary model of this work. Each entry:
# label, shear rate, blood state(s) of the panel it is drawn on, pulsatility A
# or None (not reported: drawn as a vertical line), marker style. A is the
# relative amplitude of a sinusoidal flow, A = RI / (2 - RI) for a resistivity
# index RI = (v_sys - v_dia) / v_sys.
LANDMARKS = [
    # Nagaoka & Yoshida 2006, IOVS: human retina, first-order branches
    # (108 and 147 um), laser Doppler velocimetry; RI ~ 0.7 in retinal arteries
    # (Puyo et al. 2021, laser Doppler holography) gives A ~ 0.54.
    ("retinal arterioles", 1539.0, ("healthy",), 0.7 / 1.3, "s"),
    ("retinal venules", 632.0, ("healthy",), None, "s"),
    # Khansari et al. 2017, Sci Rep: human conjunctiva, 18-21 um.
    ("conjunctival arterioles", 280.0, ("healthy",), None, "c"),
    ("conjunctival venules", 183.0, ("healthy",), None, "c"),
    ("conjunctival arterioles", 193.0, ("diabetic",), None, "c"),
    ("conjunctival venules", 154.0, ("diabetic",), None, "c"),
    # Bands (style "b"): ranges of resting shear rate; pulsatility not reported.
    # Koutsiaris et al. 2007, Biorheology: human conjunctival capillaries and
    # post-capillary venules (4-24 um), WSS 0.28-9.55 Pa, i.e. ~80-2700 1/s
    # at 3.5 mPa s; the low end is the larger post-capillary venules.
    ("post-capillary venules", (80.0, 440.0), ("healthy",), None, "b"),
    # Cheung et al. 2001, Blood: conjunctival "trickle flow" in sickle cell
    # disease, 0.05-0.2 mm/s (normal 1.1-2.4 mm/s); 8 (0.9 v) / D with
    # D = 20 um assumed (Khansari et al. 2017 diameters).
    ("sickle-cell trickle flow", (18.0, 72.0), ("hyperviscous",), None, "b"),
    # Edul et al. 2012, Crit Care Med: sublingual RBC velocity 830 um/s in
    # septic shock (1332 um/s in controls); D = 10-20 um assumed.
    ("septic shock, sublingual", (300.0, 600.0), ("healthy",), None, "b"),
    # This work (calibration of the coronary outlets, flow reversal: A ~ 1).
    ("coronary arterioles", 800.0, ("healthy", "diabetic", "hyperviscous"), 1.0, "o"),
    ("coronary venules", 400.0, ("healthy", "diabetic", "hyperviscous"), 1.0, "o"),
]


def errorMap(out, nGamma=61, nAmp=41):
    """dQ_rms of the constant closure and of R_mu over (gd_0, A), per blood state.

    Prints, for A = 0.25, 0.5 and 1, the resting shear rate below which the
    constant closure exceeds 5 and 10 per cent, and draws the map with the
    landmarks.
    """
    gam = np.logspace(1, math.log10(2500.0), nGamma)
    amp = np.linspace(0.1, 1.4, nAmp)
    G, Aa = np.meshgrid(gam, amp)
    maps = {}
    print("\nresting shear rate (1/s) below which the constant closure errs by more than 5 / 10 %:")
    for state, (mu0, mi, lam, n, a) in STATES.items():
        X = np.column_stack([np.full(G.size, mu0), np.full(G.size, mi), np.full(G.size, lam),
                             np.full(G.size, n), np.full(G.size, a), G.ravel(), Aa.ravel()])
        e = evaluateChunked(X)
        rest = e["rest_dQ_rms"].reshape(G.shape)
        single = e["single_dQ_rms"].reshape(G.shape)
        maps[state] = (rest, single)
        for A in (0.25, 0.5, 1.0):
            i = int(np.argmin(np.abs(amp - A)))
            row = []
            for lim in (5.0, 10.0):
                above = gam[rest[i] > lim]
                row.append("%6.0f" % above.max() if above.size else "     -")
            print("  %-13s A = %.2f: %s / %s" % (state, A, *row))
    with open(os.path.join(out, "error_map.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["state", "gammaRef", "amplitude", "rest_dQ_rms", "single_dQ_rms"])
        for state, (rest, single) in maps.items():
            for i in range(nAmp):
                for j in range(nGamma):
                    w.writerow([state, gam[j], amp[i], rest[i, j], single[i, j]])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
                         "axes.labelsize": 8.5, "legend.fontsize": 7.5, "axes.linewidth": 0.6,
                         "xtick.direction": "out", "ytick.direction": "out"})
    levels = [0, 2, 5, 10, 20, 40]
    cmap = LinearSegmentedColormap.from_list("copper", ["#fdf1e8", "#f2b98f", "#d9773f",
                                                        "#a84a1a", "#5e2508"], N=len(levels) - 1)
    norm = BoundaryNorm(levels, cmap.N)
    names = {"healthy": "Healthy", "diabetic": "Type 2 diabetes", "hyperviscous": "Hyperviscous"}
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7), sharey=True)
    bandIndex = {}
    for ax, (state, (rest, single)) in zip(axes, maps.items()):
        cf = ax.contourf(G, Aa, np.clip(rest, 0, 39.9), levels=levels, cmap=cmap, norm=norm)
        cs = ax.contour(G, Aa, rest, levels=[5, 10, 20], colors="#3a1a08", linewidths=0.6)
        ax.clabel(cs, fmt="%g%%", fontsize=6.5, inline_spacing=2)
        cm = ax.contour(G, Aa, single, levels=[1, 2], colors="#2f6fae", linewidths=0.9,
                        linestyles=["--", "-"])
        ax.clabel(cm, fmt="%g%%", fontsize=6.5, inline_spacing=2)
        for label, g, states, A, style in LANDMARKS:
            if state not in states:
                continue
            if style == "b":
                lo, hi = g
                k = bandIndex.get(state, 0)
                bandIndex[state] = k + 1
                y = 1.30 - 0.13 * k
                ax.plot([lo, hi], [y, y], color="#6a1b9a", lw=2.2, solid_capstyle="butt",
                        zorder=6)
                ax.plot([lo, lo], [y - 0.03, y + 0.03], color="#6a1b9a", lw=0.8, zorder=6)
                ax.plot([hi, hi], [y - 0.03, y + 0.03], color="#6a1b9a", lw=0.8, zorder=6)
                ax.text(math.sqrt(lo * hi), y + 0.035, label, fontsize=6, ha="center",
                        va="bottom", color="#4a148c", zorder=6,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.4))
                continue
            short = "art." if "arteriole" in label else "ven."
            col = {"s": "#1b5e20", "c": "#202020", "o": "#202020"}[style]
            if A is None:
                ls = ":" if style == "c" else "-."
                ax.axvline(g, color=col, lw=0.8, ls=ls, zorder=4)
                left = short == "ven."
                ax.text(g * (0.96 if left else 1.04), 0.13, short, fontsize=6.5, rotation=90,
                        ha="right" if left else "left", va="bottom", color=col)
            else:
                mk = {"o": "o", "s": "D"}[style]
                ax.plot(g, A, mk, ms=5, mfc="white", mec=col, mew=1.0, zorder=5)
                ax.text(g, A - 0.06, short, fontsize=6.5, ha="center", va="top", color=col)
        ax.set_xscale("log")
        ax.set_title(names[state], fontsize=8.5, pad=12)
        ax.set_xlabel(r"Resting shear rate $\dot\gamma_0$ (s$^{-1}$)")
    axes[0].set_ylabel("Pulsatility $A$")
    cb = fig.colorbar(cf, ax=axes, fraction=0.025, pad=0.02, ticks=levels)
    cb.set_label(r"error of constant $\mathcal{R}_0$: waveform (%)")
    fig.legend(handles=[Line2D([], [], color="#3a1a08", lw=0.6, label=r"constant $\mathcal{R}_0$"),
                        Line2D([], [], color="#2f6fae", lw=0.9, ls="--",
                               label=r"$\mathcal{R}_\mu$: 1%"),
                        Line2D([], [], color="#2f6fae", lw=0.9, label=r"$\mathcal{R}_\mu$: 2%"),
                        Line2D([], [], color="#202020", lw=0.8, ls=":",
                               label="conjunctiva (Khansari et al. 2017)"),
                        Line2D([], [], color="#1b5e20", lw=0.8, ls="-.", marker="D", mfc="white",
                               label="retina (Nagaoka & Yoshida 2006; Puyo et al.)"),
                        Line2D([], [], marker="o", ls="none", mfc="white", mec="#202020",
                               label="coronary, this work"),
                        Line2D([], [], color="#6a1b9a", lw=2.2,
                               label="range estimated from velocity (see text)")],
               loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.45, 1.0),
               columnspacing=1.2, handlelength=1.8)
    path = os.path.join(out, "fig_error_map.pdf")
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(os.path.splitext(path)[0] + ".png", dpi=300, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="sobol")
    ap.add_argument("--samples", type=int, default=2048, help="Saltelli base sample (power of 2)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--check", default="", help="summary.csv of compare_rmu_axi.py")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    if a.check:
        check(a.check)
    errorMap(a.out)

    from SALib.sample import sobol as sobolSample
    from SALib.analyze import sobol as sobolAnalyze

    problem = {"num_vars": len(BOUNDS), "names": [b[0] for b in BOUNDS],
               "bounds": [[math.log10(lo), math.log10(hi)] if lg else [lo, hi]
                          for _, lo, hi, lg in BOUNDS]}
    U = sobolSample.sample(problem, a.samples, calc_second_order=False, seed=a.seed)
    X = U.copy()
    for k, (_, _, _, lg) in enumerate(BOUNDS):
        if lg:
            X[:, k] = 10.0 ** U[:, k]
    print("\n%d model evaluations ..." % len(X), flush=True)
    Y = evaluateChunked(physical(X))

    with open(os.path.join(a.out, "samples.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(problem["names"] + list(Y))
        for i in range(len(X)):
            w.writerow(list(X[i]) + [Y[k][i] for k in Y])

    results = {}
    with open(os.path.join(a.out, "sobol_indices.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["output", "parameter", "S1", "S1_conf", "ST", "ST_conf"])
        for k, y in Y.items():
            r = sobolAnalyze.analyze(problem, y, calc_second_order=False, seed=a.seed)
            results[k] = r
            for j, name in enumerate(problem["names"]):
                w.writerow([k, name, r["S1"][j], r["S1_conf"][j], r["ST"][j], r["ST_conf"][j]])
            m = np.abs(y)
            print("\n%s: |error| median %.2f, 95th percentile %.2f, above 5%%: %.0f%% of the "
                  "samples, range [%.2f, %.2f]" % (k, np.median(m), np.percentile(m, 95),
                                                   100 * (m > 5).mean(), y.min(), y.max()))
            for j in np.argsort(-r["ST"]):
                print("  %-10s S1 = %5.3f +- %.3f   ST = %5.3f +- %.3f" % (
                    problem["names"][j], r["S1"][j], r["S1_conf"][j], r["ST"][j], r["ST_conf"][j]))

    figure(results, problem["names"], os.path.join(a.out, "fig_sobol.pdf"))
    print("\nwritten to", a.out)


def figure(results, names, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
                         "axes.labelsize": 8.5, "legend.fontsize": 8, "axes.linewidth": 0.6,
                         "xtick.direction": "out", "ytick.direction": "out"})
    panels = [("rest_dQ_rms", "error of constant $\\mathcal{R}_0$: waveform", "#c4581f"),
              ("rest_dE", "error of constant $\\mathcal{R}_0$: dissipation", "#c4581f"),
              ("single_dQ_rms", "error of $\\mathcal{R}_\\mu$: waveform", "#2f6fae"),
              ("single_dE", "error of $\\mathcal{R}_\\mu$: dissipation", "#2f6fae")]
    order = np.argsort(-results["rest_dQ_rms"]["ST"])
    fig, axes = plt.subplots(1, len(panels), figsize=(9.0, 2.5), sharey=True)
    x = np.arange(len(names))
    for ax, (key, title, col) in zip(axes, panels):
        r = results[key]
        ax.bar(x - 0.19, r["S1"][order], 0.36, yerr=r["S1_conf"][order], color=col,
               edgecolor="white", lw=0.8, error_kw=dict(lw=0.6, capsize=1.5), zorder=3)
        ax.bar(x + 0.19, r["ST"][order], 0.36, yerr=r["ST_conf"][order], facecolor="white",
               edgecolor=col, lw=1.0, error_kw=dict(lw=0.6, capsize=1.5), zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([LABELS[names[j]] for j in order])
        ax.set_title(title, fontsize=8)
        ax.set_ylim(0, 0.8)
        ax.grid(True, axis="y", lw=0.4, color="#dddddd", zorder=0)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("Sobol index")
    fig.legend(handles=[Patch(color="#707070", label="first order $S_1$"),
                        Patch(facecolor="white", edgecolor="#707070", label="total $S_T$")],
               loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(os.path.splitext(path)[0] + ".png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Supplementary results of the coronary campaign (third cycle by default).

For the R_mu runs of pra700/{healthy,diabetic,hyperviscous} and
pra2400/hyperviscous, from CoronaryArtery.csv and calibration.csv:

  fig_supp_cycle.pdf      per state: intramyocardial and transmural pressure;
                          proximal and distal flow over Q_0; mu_eff / mu_inf of
                          both elements (systole shaded)
  fig_supp_residence.pdf  fraction of the cycle with mu_eff > (1 + eps) mu_inf
                          against eps, proximal and distal
  supplement.md/.tex      mass balance of the 3D domain, cycle-to-cycle change,
                          cycle-mean reduced wall stresses, residence above
                          1.2 mu_inf and 5 mu_inf, proximal reversal fraction,
                          minimum transmural pressure

mu_eff of each element is the logged factor Phi (ArteriolarViscosityRatio,
VenularViscosityRatio: mean over outlets) times the healthy calibration
viscosity at the element's reference shear rate, computed here from the
Carreau-Yasuda law as in HydraulicCurve (RmuAxi/sobol_rmu.py). The reduced
wall stress is mu_eff gd_a with gd_a = gd_0 |Q| / Q_0 for the bed totals.

  python3 supplement.py ~/Results/cmbpe/coronary [--cycle 3]
"""
import argparse
import csv
import math
import os
import sys

import numpy as np

T = 0.85
STATES = {  # mu0, muInf, lambda, n, a (Table 1)
    "healthy": (0.0712616, 0.0041, 9.2780, 0.3572, 1.8773),
    "diabetic": (0.0973528, 0.0049622, 15.7825, 0.2670, 1.0829),
    "hyperviscous": (0.3467269, 0.0053264, 17.4172, 0.2200, 0.6957),
}
CASES = [("pra700/healthy", "healthy", "Healthy, 700 Pa", "#2f6fae"),
         ("pra700/diabetic", "diabetic", "Diabetic, 700 Pa", "#c4581f"),
         ("pra700/hyperviscous", "hyperviscous", "Hyperviscous, 700 Pa", "#2e8b57"),
         ("pra2400/hyperviscous", "hyperviscous", "Hyperviscous, 2400 Pa", "#7a4fa0")]


def muEff(gA, mu0, muInf, lam, n, a):
    """Hydraulic viscosity at apparent shear rate gA (same table as sobol_rmu.py)."""
    lg = np.linspace(math.log(1e-6), math.log(1e6), 2001)
    g = np.exp(lg)
    base = 1.0 + (lam * g) ** a
    mu = muInf + (mu0 - muInf) * base ** ((n - 1.0) / a)
    dmu = (mu0 - muInf) * (n - 1.0) * base ** ((n - 1.0 - a) / a) * lam ** a * g ** (a - 1.0)
    f = g ** 4 * mu ** 2 * (mu + g * dmu)
    I = mu0 ** 3 * 1e-24 / 4.0 + np.concatenate([[0.0], np.cumsum(0.5 * (lg[1] - lg[0]) * (f[1:] + f[:-1]))])
    tau = mu * g
    muAp = tau ** 4 / (4.0 * I)
    return float(np.exp(np.interp(math.log(gA), np.log(tau / muAp), np.log(muAp))))


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) if r[k] not in ("", None) else 0.0 for r in rows])
            for k in rows[0]} if rows else None


def window(d, c):
    t = d["t"]
    m = (t >= (c - 1) * T - 1e-12) & (t < c * T - 1e-12)
    return {k: v[m] for k, v in d.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--cycle", type=int, default=3)
    a = ap.parse_args()
    out = os.path.abspath(os.path.expanduser(a.out))
    rep = os.path.join(out, "report", "supplement")
    os.makedirs(rep, exist_ok=True)

    hp = STATES["healthy"]
    data, rows = {}, []
    for g, state, label, col in CASES:
        d = os.path.join(out, g, "rmu")
        if not os.path.isfile(os.path.join(d, "CoronaryArtery.csv")):
            print("missing", d)
            continue
        full = read(os.path.join(d, "CoronaryArtery.csv"))
        cal = list(csv.DictReader(open(os.path.join(d, "calibration.csv"))))
        q0 = sum(float(r["q0"]) for r in cal)
        gA0, gV0 = float(cal[0]["gammaA"]), float(cal[0]["gammaV"])
        muCalA, muCalV = muEff(gA0, *hp), muEff(gV0, *hp)
        muInf = STATES[state][1]
        ncyc = int((full["t"][-1] + 2e-3) // T)
        c = min(a.cycle, ncyc)
        w = window(full, c)
        qin = -w["CoronaryInletFlux"]
        qp, qd = w["CoronaryOutletFluxTotal"], w["CoronaryDistalFluxTotal"]
        muP = w["ArteriolarViscosityRatio"] * muCalA
        muD = w["VenularViscosityRatio"] * muCalV
        tauP = muP * gA0 * np.abs(qp) / q0
        tauD = muD * gV0 * np.abs(qd) / q0
        pim = w["IntramyoPressure"]
        sysm = pim >= pim.min() + 0.5 * (pim.max() - pim.min())
        means = [np.abs(-window(full, k)["CoronaryInletFlux"]).mean() * 1e6 for k in range(1, ncyc + 1)]
        row = dict(case=label, cycle=c,
                   balance=100 * np.abs(w["FlowBalance"]).max() / np.abs(qin).mean(),
                   periodicity=100 * (means[-1] - means[-2]) / means[-2] if len(means) > 1 else float("nan"),
                   tauP=tauP.mean(), tauD=tauD.mean(),
                   tauP0=muCalA * gA0, tauD0=muCalV * gV0,
                   resP12=100 * (muP > 1.2 * muInf).mean(), resD12=100 * (muD > 1.2 * muInf).mean(),
                   resP5=100 * (muP > 5 * muInf).mean(), resD5=100 * (muD > 5 * muInf).mean(),
                   reversal=100 * (qp < 0).mean(), ptmMin=w["TransmuralPressure"].min())
        rows.append(row)
        data[label] = dict(t=(w["t"] - (c - 1) * T) / T, pim=pim, ptm=w["TransmuralPressure"],
                           qp=qp / q0, qd=qd / q0, rP=muP / muInf, rD=muD / muInf, sys=sysm, col=col)

    if not rows:
        sys.exit("no R_mu runs found under " + out)

    md = ["| case | cycle | mass balance max (% of mean inflow) | last-cycle change (%) | "
          "tau_w prox / dist (Pa) | calibration (Pa) | time > 1.2 mu_inf prox / dist (%) | "
          "time > 5 mu_inf prox / dist (%) | proximal reversal (%) | min p_tm (Pa) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    tex = [r"\begin{tabular}{@{}lccccccc@{}}", r"\toprule",
           r"Case & Balance (\%) & Cycle change (\%) & $\overline{\tau}_{w,p}$, $\overline{\tau}_{w,d}$ (Pa) & "
           r"$\mu_{\rm eff}>1.2\mu_\infty$ (\%) & $\mu_{\rm eff}>5\mu_\infty$ (\%) & Reversal (\%) & "
           r"$\min p_{\rm tm}$ (Pa)\\", r"\midrule"]
    for r in rows:
        md.append("| %s | %d | %.2e | %+.2f | %.2f / %.2f | %.2f / %.2f | %.1f / %.1f | %.2f / %.2f | %.1f | %.0f |" % (
            r["case"], r["cycle"], r["balance"], r["periodicity"], r["tauP"], r["tauD"], r["tauP0"],
            r["tauD0"], r["resP12"], r["resD12"], r["resP5"], r["resD5"], r["reversal"], r["ptmMin"]))
        tex.append("%s & %.1e & $%+.1f$ & %.2f, %.2f & %.1f, %.1f & %.2f, %.2f & %.1f & %.0f\\\\" % (
            r["case"].replace(" Pa", r"\,Pa"), r["balance"], r["periodicity"], r["tauP"], r["tauD"],
            r["resP12"], r["resD12"], r["resP5"], r["resD5"], r["reversal"], r["ptmMin"]))
    tex += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(rep, "supplement.md"), "w").write("\n".join(md) + "\n")
    open(os.path.join(rep, "supplement.tex"), "w").write("\n".join(tex) + "\n")
    print("\n".join(md))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 8,
                         "axes.labelsize": 8, "legend.fontsize": 7.5, "axes.linewidth": 0.6})
    labels = list(data)
    fig, ax = plt.subplots(3, len(labels), figsize=(1.85 * len(labels) + 0.6, 5.6), sharex=True,
                           sharey="row", squeeze=False)
    for j, lab in enumerate(labels):
        d = data[lab]
        s = d["t"][d["sys"]]
        for i in range(3):
            if s.size:
                ax[i][j].axvspan(s.min(), s.max(), color="#eeeeee", lw=0, zorder=0)
            for side in ("top", "right"):
                ax[i][j].spines[side].set_visible(False)
        ax[0][j].plot(d["t"], d["pim"] / 1e3, color="#606060", lw=1.0, ls="--")
        ax[0][j].plot(d["t"], d["ptm"] / 1e3, color=d["col"], lw=1.2)
        ax[1][j].plot(d["t"], d["qp"], color=d["col"], lw=1.2)
        ax[1][j].plot(d["t"], d["qd"], color=d["col"], lw=1.0, ls=":")
        ax[1][j].axhline(0, color="#404040", lw=0.5)
        ax[2][j].semilogy(d["t"], d["rP"], color=d["col"], lw=1.2)
        ax[2][j].semilogy(d["t"], d["rD"], color=d["col"], lw=1.0, ls=":")
        ax[2][j].axhline(1.2, color="#909090", lw=0.6, ls="--")
        ax[0][j].set_title(lab, fontsize=7.5)
        ax[2][j].set_xlabel(r"$(t-t_0)/T$")
    ax[0][0].set_ylabel("Pressure (kPa)")
    ax[1][0].set_ylabel(r"$Q/Q_0$")
    ax[2][0].set_ylabel(r"$\mu_{\rm eff}/\mu_\infty$")
    fig.legend(handles=[Line2D([], [], color="#606060", ls="--", label=r"$p_{\rm im}$"),
                        Line2D([], [], color="#303030", label=r"$p_{\rm tm}$; proximal element"),
                        Line2D([], [], color="#303030", ls=":", label="distal element"),
                        Line2D([], [], color="#909090", ls="--", lw=0.6, label=r"$1.2\,\mu_\infty$")],
               loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(os.path.join(rep, "fig_supp_cycle.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(rep, "fig_supp_cycle.png"), dpi=200, bbox_inches="tight")

    fig, axr = plt.subplots(figsize=(3.6, 2.6))
    eps = np.logspace(-2, 1, 200)
    for lab in labels:
        d = data[lab]
        axr.plot(eps, [100 * (d["rP"] > 1 + e).mean() for e in eps], color=d["col"], lw=1.2, label=lab)
        axr.plot(eps, [100 * (d["rD"] > 1 + e).mean() for e in eps], color=d["col"], lw=1.0, ls=":")
    axr.set_xscale("log")
    axr.set_xlabel(r"$\varepsilon$")
    axr.set_ylabel(r"cycle with $\mu_{\rm eff}>(1+\varepsilon)\mu_\infty$ (%)")
    axr.axvline(0.2, color="#909090", lw=0.6, ls="--")
    for side in ("top", "right"):
        axr.spines[side].set_visible(False)
    axr.legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(rep, "fig_supp_residence.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(rep, "fig_supp_residence.png"), dpi=200, bbox_inches="tight")
    print("\nwritten to", rep)


if __name__ == "__main__":
    main()

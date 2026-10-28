#!/usr/bin/env python3
"""
Fair Newtonian vs non-Newtonian (Carreau-Yasuda) comparison for the coronary
outlet flow law used in CoupledLV0DCoronary3D.

Purpose
-------
Answer the question "is the non-Newtonian outlet model worth it, or does it
behave like a Newtonian model here?".  The honest answer is regime-dependent:

  * At HIGH wall shear (systole / high flow) the Carreau-Yasuda viscosity is
    essentially fully thinned, mu -> muInf, so it matches a Newtonian model
    run at muInf.  If your simulation spends most of the cycle at high shear,
    you will NOT see a difference versus Newtonian-at-muInf.
  * At LOW wall shear (diastole, low flow, small vessels) mu climbs toward mu0
    and the shear-thinning is large.  This is where the non-Newtonian model
    earns its keep -- and in the wall-shear-stress field, not the gross flux.

This script quantifies exactly where the two diverge for YOUR parameters, so
you can decide whether the operating regime makes the non-Newtonian model
meaningful.  It compares against BOTH Newtonian limits (mu0 and muInf), which
is the fair comparison -- comparing only against muInf hides the whole effect.

Edit the PARAMETERS block (or wire in your own values), then run:

    python3 nonnewtonian_comparison.py            # writes nonnewtonian_comparison.png
    python3 nonnewtonian_comparison.py outlet.csv # also overlays real operating points

If a CSV is given it is scanned for a distal-flow column (e.g.
CoronaryOutlet*DistalFlux) to mark the actual operating flow band on the plots.
"""

import math
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# PARAMETERS  (defaults mirror CoupledLV0DCoronary3D::CarreauYasuda and the
#             distal surrogate geometry -- change these to match your run)
# ---------------------------------------------------------------------------
MU0 = 0.353        # low-shear viscosity  [Pa.s in the code's unit system]
MU_INF = 0.004181  # infinite-shear viscosity
LAMBDA = 15.6821   # relaxation time
N_CY = 0.2050      # power-law index
A_YASUDA = 0.6497  # Yasuda transition exponent

# Representative distal branch operating point (used only to build Q vs dp).
TARGET_Q = 3.0e-6 / 6.0   # per-branch flow ~ (total 3 mL/s over ~6 outlets)
RD_TARGET = 1.0e9         # target distal resistance  [Pa.s/m^3]
L_DISTAL = 0.2            # distal surrogate length (Config::distalLength)


# ---------------------------------------------------------------------------
# Carreau-Yasuda viscosity and the WRMS (Rabinowitsch-Mooney) tube flow law.
# Pure-python, mirrors outletFlow() in the C++.
# ---------------------------------------------------------------------------
def mu_cy(g, mu0=MU0, muInf=MU_INF, lam=LAMBDA, n=N_CY, a=A_YASUDA):
    return muInf + (mu0 - muInf) * (1.0 + (lam * np.asarray(g, float)) ** a) ** (
        (n - 1.0) / a
    )


def _mu_dmu(mu0, muInf, lam, n, a):
    delta = mu0 - muInf
    mu = lambda g: muInf + delta * (1.0 + (lam * g) ** a) ** ((n - 1.0) / a)
    dmu = lambda g: (
        delta * (n - 1.0) * (1.0 + (lam * g) ** a) ** ((n - 1.0 - a) / a)
        * lam ** a * g ** (a - 1.0)
    )
    return mu, dmu


def flow_through_tube(dp, L, R, mu0, muInf, lam, n, a, nquad=2000):
    """Volumetric flow through a tube of radius R, length L, pressure drop dp,
    for a Carreau-Yasuda fluid, via the WRMS integral (same as the C++)."""
    if dp <= 0.0:
        return 0.0
    mu, dmu = _mu_dmu(mu0, muInf, lam, n, a)
    tauW = R * dp / (2.0 * L)

    # wall shear rate: solve g*mu(g) = tauW by bisection
    f = lambda g: g * mu(g) - tauW
    gHi = max(tauW / muInf, 1e-8)
    while f(gHi) < 0.0:
        gHi *= 2.0
    lo, hi, flo = 1e-14, gHi, f(1e-14)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0.0:
            hi = mid
        else:
            lo, flo = mid, fm
    gW = 0.5 * (lo + hi)

    # Q = pi R^3 / tauW^3 * int_0^gW g^3 mu^2 (mu + g mu') dg  (Simpson)
    gs = np.linspace(0.0, gW, nquad + 1)
    m = mu(gs)
    dm = np.zeros_like(gs)
    dm[1:] = dmu(gs[1:])
    integ = gs ** 3 * m ** 2 * (m + gs * dm)
    integ[0] = 0.0
    I = (gW / nquad) / 3.0 * (
        integ[0] + integ[-1] + 4.0 * integ[1:-1:2].sum() + 2.0 * integ[2:-1:2].sum()
    )
    return math.pi * R ** 3 * I / tauW ** 3


def flow_newtonian(dp, L, R, mu):
    """Poiseuille flow for constant viscosity mu."""
    if dp <= 0.0:
        return 0.0
    return math.pi * R ** 4 * dp / (8.0 * mu * L)


def radius_for(mu, L, Rd):
    """Newtonian radius giving resistance Rd at viscosity mu."""
    return (8.0 * mu * L / (math.pi * Rd)) ** 0.25


# ---------------------------------------------------------------------------
# Optional: pull the real operating flow band from a CSV.
# ---------------------------------------------------------------------------
def operating_flow_band(csv_path):
    try:
        import pandas as pd

        df = pd.read_csv(csv_path)
    except Exception as exc:  # noqa: BLE001
        print(f"  (could not read {csv_path}: {exc})")
        return None
    cols = [c for c in df.columns if "DistalFlux" in c or c.endswith("Flux")]
    cols = [c for c in cols if "Total" not in c and "Fraction" not in c]
    if not cols:
        return None
    vals = np.abs(df[cols].to_numpy(float)).ravel()
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        return None
    return np.percentile(vals, 5), np.percentile(vals, 95)


# ---------------------------------------------------------------------------
# Build the figure.
# ---------------------------------------------------------------------------
def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else None
    band = operating_flow_band(csv_path) if csv_path else None

    # Radius calibrated at mu0 (as the model's Newtonian sizing would do).
    R = radius_for(MU0, L_DISTAL, RD_TARGET)

    # Sweep pressure drop; compute flow under the three models.
    dps = np.logspace(math.log10(RD_TARGET * TARGET_Q) - 2.5,
                      math.log10(RD_TARGET * TARGET_Q) + 1.5, 60)
    q_cy = np.array([flow_through_tube(dp, L_DISTAL, R, MU0, MU_INF,
                                       LAMBDA, N_CY, A_YASUDA) for dp in dps])
    q_mu0 = np.array([flow_newtonian(dp, L_DISTAL, R, MU0) for dp in dps])
    q_inf = np.array([flow_newtonian(dp, L_DISTAL, R, MU_INF) for dp in dps])

    # Effective viscosity implied by CY flow (mu_eff = Poiseuille mu matching Q).
    with np.errstate(divide="ignore", invalid="ignore"):
        mu_eff = math.pi * R ** 4 * dps / (8.0 * L_DISTAL * q_cy)

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("Carreau-Yasuda vs Newtonian coronary outlet law", fontsize=14)

    # (1) viscosity vs shear rate
    g = np.logspace(-1, 5, 400)
    ax[0, 0].loglog(g, mu_cy(g), "C0", lw=2, label="Carreau-Yasuda")
    ax[0, 0].axhline(MU0, color="C3", ls="--", label=f"mu0 = {MU0:g}")
    ax[0, 0].axhline(MU_INF, color="C2", ls="--", label=f"muInf = {MU_INF:g}")
    ax[0, 0].set(xlabel="shear rate  [1/s]", ylabel="viscosity  [Pa.s]",
                 title="Blood viscosity vs shear rate")
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(True, which="both", alpha=0.3)

    # (2) Q vs dp for the three models
    ax[0, 1].loglog(dps, q_cy * 1e6, "C0", lw=2, label="Carreau-Yasuda")
    ax[0, 1].loglog(dps, q_mu0 * 1e6, "C3", ls="--", label="Newtonian @ mu0")
    ax[0, 1].loglog(dps, q_inf * 1e6, "C2", ls="--", label="Newtonian @ muInf")
    ax[0, 1].set(xlabel="pressure drop  [Pa]", ylabel="flow  [mL/s]",
                 title="Distal-branch flow vs pressure drop")
    ax[0, 1].legend(fontsize=8)
    ax[0, 1].grid(True, which="both", alpha=0.3)

    # (3) relative difference of CY vs each Newtonian limit
    rel_mu0 = 100.0 * (q_cy - q_mu0) / q_mu0
    rel_inf = 100.0 * (q_cy - q_inf) / q_inf
    ax[1, 0].semilogx(q_cy * 1e6, rel_mu0, "C3", label="vs Newtonian @ mu0")
    ax[1, 0].semilogx(q_cy * 1e6, rel_inf, "C2", label="vs Newtonian @ muInf")
    ax[1, 0].axhline(0, color="k", lw=0.8)
    ax[1, 0].set(xlabel="Carreau-Yasuda flow  [mL/s]",
                 ylabel="flow difference  [%]",
                 title="Where non-Newtonian matters")
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(True, which="both", alpha=0.3)

    # (4) effective viscosity across the flow range
    ax[1, 1].semilogx(q_cy * 1e6, mu_eff, "C0", lw=2, label="mu_eff (CY)")
    ax[1, 1].axhline(MU0, color="C3", ls="--", label="mu0")
    ax[1, 1].axhline(MU_INF, color="C2", ls="--", label="muInf")
    ax[1, 1].set(xlabel="flow  [mL/s]", ylabel="effective viscosity  [Pa.s]",
                 title="Effective viscosity vs flow")
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].grid(True, which="both", alpha=0.3)

    # overlay the real operating band if we have it
    if band is not None:
        for a_ in (ax[0, 1], ax[1, 0], ax[1, 1]):
            a_.axvspan(band[0] * 1e6, band[1] * 1e6, color="gold", alpha=0.2,
                       label="operating band")
        print(f"  operating flow band (5-95 pct): "
              f"{band[0]*1e6:.3f} - {band[1]*1e6:.3f} mL/s")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = "nonnewtonian_comparison.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")

    # console summary
    print("\nSummary at representative flows:")
    for qtarget in (0.1 * TARGET_Q, TARGET_Q, 5.0 * TARGET_Q):
        i = int(np.argmin(np.abs(q_cy - qtarget)))
        print(f"  Q ~ {q_cy[i]*1e6:6.3f} mL/s :  mu_eff={mu_eff[i]:.4f} Pa.s"
              f"   CY vs mu0 = {rel_mu0[i]:+6.1f}%"
              f"   CY vs muInf = {rel_inf[i]:+6.1f}%")


if __name__ == "__main__":
    main()

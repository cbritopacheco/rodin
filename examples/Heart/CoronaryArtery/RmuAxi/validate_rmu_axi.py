#!/usr/bin/env python3
"""Report of validate_rmu_axi.sh.

Steady runs: the computed inflow and the wall shear stress of the stub wall (10 < x/R_1 < 20, developed flow) against the
exact solution for fully developed Carreau-Yasuda flow in a tube of radius R
and length L under the pressure drop Pbar (read from each log),

  tau_w = Pbar R / (2 L),   Q = pi R^3 I(tau_w) / tau_w^3,
  I = int_0^{gd_w} gd^3 mu^2 (mu + gd mu') d gd,   mu(gd_w) gd_w = tau_w,

evaluated by adaptive quadrature, independently of the C++ table. Pulsatile
runs: the waveform error e_Q of each outlet against the full tube over the
last cycle, as in the paper.

  python3 validate_rmu_axi.py rmu_axi_validation
"""
import csv
import math
import os
import re
import sys

from scipy.integrate import quad
from scipy.optimize import brentq

MU0, MUINF, LAM, N, A = 0.0712616, 0.0041, 9.2780, 0.3572, 1.8773
R, L = 1.0e-4, 40.0e-4


def mu(g):
    return MUINF + (MU0 - MUINF) * (1.0 + (LAM * g) ** A) ** ((N - 1.0) / A)


def dmu(g):
    return (MU0 - MUINF) * (N - 1.0) * (1.0 + (LAM * g) ** A) ** ((N - 1.0 - A) / A) \
        * LAM ** A * g ** (A - 1.0)


def exactFlow(pbar):
    tau = pbar * R / (2.0 * L)
    gw = brentq(lambda g: mu(g) * g - tau, 0.0, 1.0e9, xtol=1e-14, rtol=1e-14)
    f = lambda g: g ** 3 * mu(g) ** 2 * (mu(g) + g * dmu(g))
    pts = [x for x in (0.1 / LAM, 1.0 / LAM, 10.0 / LAM) if x < gw]
    I, _ = quad(f, 0.0, gw, points=pts or None, epsabs=0.0, epsrel=1e-12, limit=500)
    return math.pi * R ** 3 * I / tau ** 3, tau


def read(path):
    rows = list(csv.DictReader(open(path)))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def pbarOf(log):
    m = re.search(r"Pbar=([0-9.eE+-]+)", open(log).read())
    return float(m.group(1))


def main(out):
    print("Steady straight tube, healthy blood: error against the exact WRMS solution")
    print("%6s %4s %12s %12s %10s %10s" % ("gd_0", "nr", "Q (m^3/s)", "Q exact", "eQ (%)", "eTau (%)"))
    rows = []
    for g in sorted(os.listdir(os.path.join(out, "steady")), key=lambda s: float(s[1:])):
        for nr in ("nr8", "nr16", "nr32"):
            d = os.path.join(out, "steady", g, nr)
            if not os.path.exists(os.path.join(d, "RmuAxi.csv")):
                continue
            data = read(os.path.join(d, "RmuAxi.csv"))
            q, w = data["qIn"][-1], abs(data["wssStub"][-1])
            qe, te = exactFlow(pbarOf(os.path.join(d, "log.txt")))
            eq, et = 100.0 * (q - qe) / qe, 100.0 * (w - te) / te
            rows.append((g[1:], nr[2:], eq, et))
            print("%6s %4s %12.5e %12.5e %+10.4f %+10.4f" % (g[1:], nr[2:], q, qe, eq, et))
    print("\nObserved order of |eQ| between successive meshes:")
    for g in sorted({r[0] for r in rows}, key=float):
        e = [abs(r[2]) for r in rows if r[0] == g]
        if len(e) == 3 and min(e) > 0:
            print("  gd_0 = %s: %.2f, %.2f" % (g, math.log2(e[0] / e[1]), math.log2(e[1] / e[2])))

    print("\nPulsatile straight tube truncated at 20 R_1 (A = 1.2): waveform error e_Q (%)")
    for g in sorted(os.listdir(os.path.join(out, "pulsatile")), key=lambda s: float(s[1:])):
        base = os.path.join(out, "pulsatile", g)
        try:
            full = read(os.path.join(base, "none", "RmuAxi.csv"))
        except FileNotFoundError:
            continue
        last = max(full["cycle"])
        idx = [i for i, c in enumerate(full["cycle"]) if c == last]
        qf = [full["qIn"][i] for i in idx]
        line = "  gd_0 = %s:" % g[1:]
        for c, name in (("single", "R_mu"), ("rest", "R_0")):
            p = os.path.join(base, c, "RmuAxi.csv")
            if not os.path.exists(p):
                continue
            d = read(p)
            qc = [d["qIn"][i] for i in idx]
            err = 100.0 * math.sqrt(sum((a - b) ** 2 for a, b in zip(qc, qf)) / len(qf)) \
                / (sum(abs(b) for b in qf) / len(qf))
            line += "  %s %.3f" % (name, err)
        print(line)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "rmu_axi_validation"))

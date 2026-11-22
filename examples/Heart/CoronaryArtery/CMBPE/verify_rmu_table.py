#!/usr/bin/env python3
"""Check of the tabulated WRMS closure used by RmuAxi and CoronaryArtery.

For each Carreau-Yasuda state:
  1. reference mu_eff(tau_w) = tau_w^4 / (4 I(tau_w)), I = int_0^tau t^2 gdot(t) dt,
     with gdot(t) from inverting tau = mu(gdot) gdot (bisection, 1e-14) and
     I by adaptive Gauss-Kronrod (scipy.quad, rel 1e-12) in log variable;
  2. the solver's table: 241 log-spaced nodes in tau_w in [1e-6, 1e4] Pa,
     log-log linear interpolation in the apparent shear rate gdot_a = tau_w/mu_eff;
  3. max relative error of mu_eff(gdot_a) from the table against the reference
     at the log-midpoints between nodes (worst case for linear interpolation),
     over gdot_a in [1e-2, 1e4] 1/s;
  4. analytical tangent dQ/dp of the tube law against central finite differences.
"""
import math
import numpy as np
from scipy.integrate import quad

STATES = {  # mu0, muinf (Pa s), lambda (s), a, n  -- Table rheology
    "healthy": (71.2616e-3, 4.1000e-3, 9.2780, 1.8773, 0.3572),
    "diabetic": (97.3528e-3, 4.9622e-3, 15.7825, 1.0829, 0.2670),
    "hyperviscous": (346.7269e-3, 5.3264e-3, 17.4172, 0.6957, 0.2200),
}

def mu(g, p):
    m0, mi, lam, a, n = p
    return mi + (m0 - mi) * (1 + (lam * g) ** a) ** ((n - 1) / a)

def gdot(tau, p):           # tau = mu(g) g is strictly increasing
    lo, hi = 0.0, tau / p[1]
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mu(mid, p) * mid < tau: lo = mid
        else: hi = mid
        if hi - lo <= 1e-14 * hi: break
    return 0.5 * (lo + hi)

def I(tau, p):              # int_0^tau t^2 gdot(t) dt, s = ln t
    lo = math.log(tau) - 60
    f = lambda s: math.exp(3 * s) * gdot(math.exp(s), p)
    v, _ = quad(f, lo, math.log(tau), epsrel=1e-12, epsabs=0, limit=400)
    return v

def mu_eff(tau, p): return tau ** 4 / (4 * I(tau, p))

def run():
    for name, p in STATES.items():
        taus = np.logspace(-6, 4, 241)
        me = np.array([mu_eff(t, p) for t in taus])
        ga = taus / me
        lga, lme = np.log(ga), np.log(me)
        # off-node checks: log-midpoints between nodes, restricted to gdot_a in [1e-2, 1e4]
        errs = []
        for k in range(len(taus) - 1):
            tau = math.sqrt(taus[k] * taus[k + 1])
            ref = mu_eff(tau, p)
            g = tau / ref
            if not (1e-2 <= g <= 1e4):
                continue
            tab = math.exp(np.interp(math.log(g), lga, lme))
            errs.append(abs(tab / ref - 1))
        errs = np.array(errs)
        # tangent of Q(dp) for a tube R, L: Q = pi R^3 I(tau)/tau^3, tau = dp R / 2L
        R, L = 25e-6, 1e-3
        worst = 0.0
        for tau in np.logspace(-2, 2, 25):
            Q = lambda t: math.pi * R ** 3 * I(t, p) / t ** 3
            # analytic: dQ/dtau = pi R^3 (gdot(tau)/tau - 3 I/tau^4)
            dq = math.pi * R ** 3 * (gdot(tau, p) / tau - 3 * I(tau, p) / tau ** 4)
            h = 1e-5 * tau
            fd = (Q(tau + h) - Q(tau - h)) / (2 * h)
            worst = max(worst, abs(dq / fd - 1))
        print("%-12s table max rel err %.3e (median %.1e) on gdot_a in [1e-2,1e4]; "
              "tangent vs FD max rel err %.1e; mu_eff(800)=%.4f mPa s, mu_eff(400)=%.4f mPa s"
              % (name, errs.max(), np.median(errs), worst,
                 1e3 * math.exp(np.interp(math.log(800), lga, lme)),
                 1e3 * math.exp(np.interp(math.log(400), lga, lme))))

if __name__ == "__main__":
    run()

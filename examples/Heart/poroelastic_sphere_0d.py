"""
0D poroelastic sphere (perfused left ventricle) — Python port of Rodin's
``Rodin::Heart::PoroelasticSphereT`` with the constitutive laws of
Ruz et al., "A patient-specific perfused non-Newtonian FSI model of the LCA":

* skeletal energy: Usyk exponential law reduced to (lambda_r, lambda_theta)
  (Sec. 4.3, Eqs. 88-91; P_{M,theta} = (1/2) dPsi/dlambda_theta, i.e. the
  c0/4 factor consistent with the energy),
* pore storage: Bruinsma vascular compliance (Sec. 4.4, Eqs. 100, 104),
  or the linear K_Phi law,
* active Bestel-Clement-Sorine / Hill-Maxwell unit, mitral/aortic valves,
  2-element proximal + distal Windkessel (Newtonian branches), two-resistor
  coronary source, exactly as in the C++ model.

Unknowns per step: x = [y, phi, pv, par, pd, ec, kc, tauc, w].  Time stepping
is BDF1 on the first step and BDF2 afterwards, as in Rodin.  The Jacobian is
built by forward finite differences (9 residual evaluations), which keeps the
code short and makes any parameter a valid sensitivity input — just change a
field of :class:`Params` and call :func:`simulate`.

Only the Newtonian Windkessel branch law is ported (``Rp``, ``Rd``); the
Cross/Carreau-Yasuda WRMS branches of the C++ model are not.

Example::

    from poroelastic_sphere_0d import Params, simulate, summary
    p = Params()                       # defaults = examples/Heart/PoroelasticSphere_0D.cpp
    out = simulate(p, n_cycles=10)     # dict of numpy arrays, one entry per step
    print(summary(out, p))             # scalar QoIs of the last cycle
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Callable, Dict

import numpy as np

# ----------------------------------------------------------------------------
# Default time functions (identical to PoroelasticSphere_0D.cpp)
# ----------------------------------------------------------------------------
T_CYCLE = 0.85


def periodic_activation(t: float) -> float:
    tau = t - T_CYCLE * np.floor(t / T_CYCLE)
    if tau < 0.13:
        return 0.0
    if tau < 0.141:
        return 35.0 * ((tau - 0.13) / 0.011)
    if tau < 0.281:
        return 35.0
    if tau < 0.361:
        return 35.0 - 55.0 * ((tau - 0.281) / 0.08)
    if tau < 0.45:
        return -20.0
    return 0.0


def relaxation_m0(ec: float) -> float:
    if ec <= 0.0:
        return 1.6
    if ec >= 2.0:
        return 1.0
    s = ec / 2.0
    return (1.0 - s) * 1.6 + s * 1.0


def atrial_pressure(t: float) -> float:
    tau = t - T_CYCLE * np.floor(t / T_CYCLE)
    pmin, pmax, p2 = 500.0, 1000.0, 1250.0
    t1, t2, t3, t4, t5, t6 = 0.02, 0.15, 0.17, 0.56, 0.62, 0.85
    if tau < t1:
        a = -(tau - t1) / t1
        return a * pmin + (1 - a) * pmax
    if tau < t2:
        return pmax
    if tau < t3:
        a = -(tau - t3) / (t3 - t2)
        return a * pmax + (1 - a) * pmin
    if tau < t4:
        a = -(tau - t4) / (t4 - t3)
        return a * pmin + (1 - a) * p2
    if tau < t5:
        return p2
    if tau < t6:
        a = -(tau - t6) / (t6 - t5)
        return a * p2 + (1 - a) * pmin
    return pmin


# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
@dataclass
class Params:
    # Geometry / reference porosity
    R0: float = 2.4e-2          # endocardial reference radius R_in [m]
    d0: float = 1.45e-2         # reference wall thickness [m]
    phi0: float = 0.1
    # Usyk skeletal law
    c0: float = 0.88e3          # [Pa]
    brr: float = 3.58
    btt: float = 18.5
    brt: float = 2.8
    # Pore storage: "bruinsma" (p0, dPhiMax) or "linear" (KPhi)
    pore_storage: str = "bruinsma"
    p0: float = 1.0e4           # [Pa]
    dPhiMax: float = 0.05
    pore_cutoff: float = 1e-3   # C1 linear continuation beyond (1-cutoff)*dPhiMax
    KPhi: float = 2.0e5         # [Pa], linear law only
    # Active law (Bestel-Clement-Sorine / Hill-Maxwell)
    Es: float = 3.0e6
    mu: float = 70.0
    eta: float = 70.0
    alpha: float = 1.5
    alphaR: float = 0.12
    k0: float = 1.0e5
    sigma0: float = 1.25e5
    # Perfusion source
    gammaAr: float = 7.0e-10
    gammaVen: float = 7.0e-10
    # Windkessel (Newtonian branches)
    Rp: float = 5.0e7
    Cp: float = 6e-9
    Rd: float = 1.0e8
    Cd: float = 1.0e-9
    # Valves
    Kat: float = 6.0e-7
    Kp: float = 5.0e-11
    Kar: float = 1.0e-7
    cavityCapacity: float = 5.0e-12
    # Numerics
    absRegularization: float = 1e-14
    nq: int = 8                 # Gauss-Legendre points across the wall
    bdf2: bool = True
    # Time functions
    u: Callable[[float], float] = field(default=periodic_activation)
    pAt: Callable[[float], float] = field(default=atrial_pressure)
    pSv: Callable[[float], float] = field(default=lambda t: 1.0e3)
    m0: Callable[[float], float] = field(default=relaxation_m0)
    qArterialExternal: Callable[[float], float] = field(default=lambda t: 0.0)
    qPerfusionExternal: Callable[[float], float] = field(default=lambda t: 0.0)

    def with_(self, **kw) -> "Params":
        """Copy with some fields replaced (handy for sensitivity sweeps)."""
        return replace(self, **kw)


# Unknown indices
Y, PHI, PV, PAR, PD, EC, KC, TAUC, W = range(9)
NVAR = 9


# ----------------------------------------------------------------------------
# Constitutive laws
# ----------------------------------------------------------------------------
def usyk_energy(lr, lt, p: Params):
    A, B = lr * lr - 1.0, lt * lt - 1.0
    Q = 0.25 * (p.brr * A * A + 2.0 * p.btt * B * B + 2.0 * p.brt * A * B)
    return 0.5 * p.c0 * (np.exp(Q) - 1.0)


def usyk_stress(lr, lt, p: Params):
    """Radial Piola components (P_{M,R}, P_{M,theta}) of the reduced Usyk law."""
    A, B = lr * lr - 1.0, lt * lt - 1.0
    Q = 0.25 * (p.brr * A * A + 2.0 * p.btt * B * B + 2.0 * p.brt * A * B)
    e = 0.5 * p.c0 * np.exp(Q)
    PRR = e * lr * (p.brr * A + p.brt * B)              # dPsi/dlr
    Ptt = 0.5 * e * lt * (2.0 * p.btt * B + p.brt * A)  # (1/2) dPsi/dlt
    return PRR, Ptt


def pore_pressure(phi: float, p: Params):
    """Pore storage pressure p_P(Phi) and its derivative."""
    x = phi - p.phi0
    if p.pore_storage != "bruinsma":
        return p.KPhi * x, p.KPhi
    D, xmax = p.dPhiMax, (1.0 - p.pore_cutoff) * p.dPhiMax
    if x < xmax:
        return -p.p0 * np.log(1.0 - x / D), p.p0 / (D - x)
    slope = p.p0 / (D - xmax)
    return -p.p0 * np.log(1.0 - xmax / D) + slope * (x - xmax), slope


def recruitment(ec: float) -> float:
    """Piecewise-linear n_0(e_c) of CCMLC2014."""
    xs = [-0.4, 0.3, 0.73, 1.0, 1.3, 2.4]
    ys = [0.0, 0.38, 0.74, 1.0, 1.0, 0.0]
    if ec < xs[1]:
        v = (ys[1] - ys[0]) / (xs[1] - xs[0]) * (ec - xs[1]) + ys[1]
    elif ec < xs[2]:
        v = (ys[2] - ys[1]) / (xs[2] - xs[1]) * (ec - xs[2]) + ys[2]
    elif ec < xs[3]:
        v = (ys[3] - ys[2]) / (xs[3] - xs[2]) * (ec - xs[3]) + ys[3]
    elif ec < xs[4]:
        v = ys[3]
    elif ec < xs[5]:
        v = (ys[5] - ys[4]) / (xs[5] - xs[4]) * (ec - xs[5]) + ys[5]
    else:
        v = 0.0
    return max(v, 0.0)


# ----------------------------------------------------------------------------
# Wall: exact radial kinematics + Gauss-Legendre quadrature across the wall
# ----------------------------------------------------------------------------
class Wall:
    def __init__(self, p: Params):
        self.p = p
        xg, wg = np.polynomial.legendre.leggauss(p.nq)
        self.x = 0.5 * (xg + 1.0)
        self.w = 0.5 * wg
        self.Rin, self.Rout = p.R0, p.R0 + p.d0
        self.R = self.Rin + (self.Rout - self.Rin) * self.x
        self.wq = (self.Rout - self.Rin) * self.w * self.R ** 2
        self.Vw0 = 4.0 / 3.0 * np.pi * (self.Rout ** 3 - self.Rin ** 3)
        self.Rmid = 0.5 * (self.Rin + self.Rout)

    def kin(self, R, y, phi, yDot, phiDot):
        p = self.p
        rin = p.R0 + y
        J = 1.0 - p.phi0 + phi
        shell = R ** 3 - self.Rin ** 3
        r = np.cbrt(rin ** 3 + J * shell)
        v1 = rin * rin / (r * r)
        v2 = shell / (3.0 * r * r)
        lt = r / R
        lr = J / (lt * lt)
        dv1 = -2.0 * lr * v1 / r
        dv2 = R * R / (r * r) - 2.0 * lr * v2 / r
        rDot = v1 * yDot + v2 * phiDot
        ltDot = rDot / R
        lrDot = phiDot / (lt * lt) - 2.0 * J * ltDot / lt ** 3
        return r, lt, lr, v1, v2, dv1, dv2, ltDot, lrDot

    def evaluate(self, y, phi, yDot, phiDot, ec):
        """Returns wall pressure P_v(y, Phi), averaged multiplier lambdaBar,
        mid-wall hoop strain e_1D and active stress sigma."""
        p = self.p
        # Mid-wall contractile unit
        _, ltm, *_ = self.kin(self.Rmid, y, phi, yDot, phiDot)
        strain1D = 0.5 * (ltm * ltm - 1.0)
        h = 1.0 + 2.0 * ec
        sigma = p.Es / (h * h) * (strain1D - ec)

        r, lt, lr, v1, v2, dv1, dv2, ltDot, lrDot = self.kin(self.R, y, phi, yDot, phiDot)
        PRR_pas, Ptt_pas = usyk_stress(lr, lt, p)
        PRR = PRR_pas + p.eta * lr * lr * lrDot
        Ptt = Ptt_pas + p.eta * lt * lt * ltDot + 0.5 * sigma * lt
        I1 = np.sum(self.wq * (PRR * dv1 + 2.0 * Ptt * v1 / self.R))
        I2 = np.sum(self.wq * (PRR * dv2 + 2.0 * Ptt * v2 / self.R))
        rin = p.R0 + y
        wallPressure = I1 / (rin * rin)
        lambdaBar = -3.0 / (self.Rout ** 3 - self.Rin ** 3) * I2
        return wallPressure, lambdaBar, strain1D, sigma


# ----------------------------------------------------------------------------
# Coupled residual
# ----------------------------------------------------------------------------
class Model:
    def __init__(self, p: Params):
        self.p = p
        self.wall = Wall(p)

    def time_coeffs(self, dt, bdf2):
        if bdf2:
            return 1.5 / dt, -2.0 / dt, 0.5 / dt
        return 1.0 / dt, -1.0 / dt, 0.0

    def residual(self, x, xn, xnm1, t, dt, bdf2, diag=None):
        p, w = self.p, self.wall
        a0, a1, a2 = self.time_coeffs(dt, bdf2)
        xdot = a0 * x + a1 * xn + a2 * xnm1
        y, phi, pv, par, pd, ec, kc, tauc, ww = x
        yDot, phiDot, pvDot, parDot, pdDot, ecDot, kcDot, taucDot, wDot = xdot

        pAt, pSv = p.pAt(t), p.pSv(t)
        wallPressure, lambdaBar, strain1D, _ = w.evaluate(y, phi, yDot, phiDot, ec)
        pP, _ = pore_pressure(phi, p)
        pf = pP - lambdaBar
        Qin = p.gammaAr * (par - pf)
        Qout = p.gammaVen * (pf - pSv)
        rin = p.R0 + y

        # Valves (same branches as C++ evaluateValveDiagnostics)
        if pv <= pAt:
            cavityFlux = p.Kat * (pv - pAt)
        elif pv <= par:
            cavityFlux = p.Kp * (pv - pAt)
        else:
            cavityFlux = p.Kar * (pv - par) + p.Kp * (par - pAt)
        aortic = max(0.0, p.Kar * (pv - par))
        qp = (par - pd) / p.Rp
        qd = (pSv - pd) / p.Rd

        # Active unit
        h = 1.0 + 2.0 * ec
        uu = p.u(t)
        up, un = max(uu, 0.0), max(-uu, 0.0)
        n0 = recruitment(ec)
        rate = up + ww * un + p.alpha * np.sqrt(ecDot ** 2 + p.absRegularization ** 2)

        R = np.empty(NVAR)
        R[Y] = wallPressure - pv
        R[PHI] = phiDot - (Qin + p.qPerfusionExternal(t) - Qout) / w.Vw0
        R[PV] = p.cavityCapacity * pvDot + cavityFlux + 4.0 * np.pi * rin * rin * yDot
        R[PAR] = p.Cp * parDot + (qp - aortic) + Qin + p.qArterialExternal(t)
        R[PD] = p.Cd * pdDot + (-qp - qd)
        R[EC] = (tauc + p.mu * ecDot) * h ** 3 - p.Es * (strain1D - ec) * (1.0 + 2.0 * strain1D)
        R[KC] = kcDot + rate * kc - n0 * p.k0 * up
        R[TAUC] = taucDot - kc * ecDot - n0 * p.sigma0 * up + rate * tauc
        if p.alphaR <= np.finfo(float).eps:
            R[W] = ww - p.m0(ec)
        else:
            R[W] = wDot - (p.m0(ec) - ww) / p.alphaR

        if diag is not None:
            diag.update(lambdaBar=lambdaBar, pf=pf, Qcor=Qin, pat=pAt)
        return R

    def jacobian_fd(self, x, xn, xnm1, t, dt, bdf2, R0, eps=1e-7):
        J = np.empty((NVAR, NVAR))
        for j in range(NVAR):
            h = eps * max(1.0, abs(x[j]))
            xp = x.copy()
            xp[j] += h
            J[:, j] = (self.residual(xp, xn, xnm1, t, dt, bdf2) - R0) / h
        return J

    def newton(self, x0, xn, xnm1, t, dt, bdf2, atol=1e-8, rtol=1e-8, stol=1e-10, maxit=50):
        x = x0.copy()
        R = self.residual(x, xn, xnm1, t, dt, bdf2)
        n0 = np.linalg.norm(R)
        for it in range(maxit):
            if np.linalg.norm(R) < max(atol, rtol * n0):
                return x, True, it
            J = self.jacobian_fd(x, xn, xnm1, t, dt, bdf2, R)
            dx = np.linalg.solve(J, -R)
            x += dx
            R = self.residual(x, xn, xnm1, t, dt, bdf2)
            if np.linalg.norm(dx) < stol:
                return x, True, it + 1
        return x, np.linalg.norm(R) < max(atol, rtol * n0), maxit


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def initial_state(p: Params) -> np.ndarray:
    x = np.zeros(NVAR)
    x[PHI] = p.phi0
    x[PV] = p.pAt(0.0) - 100.0
    x[PAR] = 11000.0
    x[PD] = 10000.0
    x[W] = p.m0(0.0)
    return x


def simulate(p: Params, n_cycles: float = 10, dt: float = 1e-3, x0=None,
             verbose: bool = False) -> Dict[str, np.ndarray]:
    """Integrate ``n_cycles`` cardiac cycles; returns arrays keyed like the
    C++ CSV (t, y, phi, pv, par, pd, ec, kc, tauc, w, V, Q, pat, lambdaBar,
    pf, Qcor) plus ``converged`` and ``iterations``."""
    m = Model(p)
    x = initial_state(p) if x0 is None else np.asarray(x0, float).copy()
    xn, xnm1 = x.copy(), x.copy()
    t, nsteps = 0.0, int(round(n_cycles * T_CYCLE / dt))
    keys = ["t", "y", "phi", "pv", "par", "pd", "ec", "kc", "tauc", "w", "V", "Q",
            "pat", "lambdaBar", "pf", "Qcor", "converged", "iterations"]
    out = {k: np.empty(nsteps) for k in keys}
    Vprev = 4.0 / 3.0 * np.pi * p.R0 ** 3
    tn, tnm1 = 0.0, 0.0
    for i in range(nsteps):
        bdf2 = p.bdf2 and (tn > tnm1 + 0.5 * dt)
        tnp1 = t + dt
        xnew, ok, its = m.newton(xn, xn, xnm1, tnp1, dt, bdf2)
        if not ok:
            if verbose:
                print(f"Newton failed at step {i}, t={tnp1}")
            for k in keys:
                out[k] = out[k][:i]
            break
        diag = {}
        m.residual(xnew, xn, xnm1, tnp1, dt, bdf2, diag)
        xnm1, xn, tnm1, tn, t = xn, xnew, tn, tnp1, tnp1
        V = 4.0 / 3.0 * np.pi * (p.R0 + xn[Y]) ** 3
        vals = [t, *xn, V, (V - Vprev) / dt, diag["pat"], diag["lambdaBar"],
                diag["pf"], diag["Qcor"], float(ok), its]
        Vprev = V
        for k, v in zip(keys, vals):
            out[k][i] = v
    return out


def summary(out: Dict[str, np.ndarray], p: Params) -> Dict[str, float]:
    """Scalar quantities of interest over the last cycle (for sensitivity
    analysis)."""
    t = out["t"]
    m = t > t[-1] - T_CYCLE + 1e-9
    V, pv, pf, Q = out["V"][m], out["pv"][m], out["pf"][m], out["Qcor"][m]
    return dict(
        pv_max=pv.max(), pv_min=pv.min(),
        EDV=V.max(), ESV=V.min(), EF=(V.max() - V.min()) / V.max(),
        SV=V.max() - V.min(),
        pf_max=pf.max(), pf_min=pf.min(),
        phi_max=out["phi"][m].max(), phi_min=out["phi"][m].min(),
        Qcor_mean=Q.mean(), Qcor_max=Q.max(), Qcor_min=Q.min(),
        lambdaBar_min=out["lambdaBar"][m].min(),
    )


def check_usyk_consistency(p: Params | None = None, h: float = 1e-6) -> float:
    """Max relative error of (P_{M,R}, P_{M,theta}) vs finite differences of
    the energy. Returns ~1e-8; the manuscript's Eq. (92) with c0/2 would give
    ~1 for P_{M,theta}."""
    p = p or Params()
    err = 0.0
    for lr, lt in [(0.83, 1.17), (1.05, 0.97), (0.7, 1.3)]:
        PRR, Ptt = usyk_stress(lr, lt, p)
        dWr = (usyk_energy(lr + h, lt, p) - usyk_energy(lr - h, lt, p)) / (2 * h)
        dWt = (usyk_energy(lr, lt + h, p) - usyk_energy(lr, lt - h, p)) / (2 * h)
        s = abs(PRR) + abs(Ptt)
        err = max(err, abs(PRR - dWr) / s, abs(Ptt - 0.5 * dWt) / s)
    return err


if __name__ == "__main__":
    import sys
    import time

    p = Params()
    print("Usyk stress/energy consistency:", check_usyk_consistency(p))
    n = float(sys.argv[1]) if len(sys.argv) > 1 else 3
    t0 = time.time()
    out = simulate(p, n_cycles=n, verbose=True)
    print(f"{len(out['t'])} steps in {time.time() - t0:.1f} s, "
          f"max Newton its {int(out['iterations'].max())}")
    for k, v in summary(out, p).items():
        print(f"  {k:14s} {v:.6g}")
    np.savetxt("poroelastic_sphere_0d_py.csv",
               np.column_stack([out[k] for k in list(out)[:16]]), delimiter=",",
               header=",".join(list(out)[:16]), comments="")

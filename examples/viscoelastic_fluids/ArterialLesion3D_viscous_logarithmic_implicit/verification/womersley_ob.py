#!/usr/bin/env python3
"""Pulsatile, fully developed pipe flow of an Oldroyd-B fluid: exact solution.

For u = u(r, t) e_x the upper-convected terms of the shear stress vanish
(sigma_rr = 0), so

    sigma_xr + lambda d_t sigma_xr = eta_p d_r u,
    rho d_t u = G(t) + (1/r) d_r [r (eta_s d_r u + sigma_xr)],

is linear for any amplitude, and each Fourier mode e^{i n omega t} is the
Womersley flow with the complex viscosity eta*(n omega) = eta_s + eta_p/(1 + i
n omega lambda). With the flow rate Q(t) = Ubar (pi R^2) q(t),
q(t) = 1 + A sin(omega t) = 1 + 2 Re[c_1 e^{i omega t}], c_1 = -i A/2:

    u/Ubar      = 2 (1 - r^2/R^2) + 2 Re[c_1 phi_1(r) e^{i omega t}],
    phi_1(r)    = [1 - I_0(k r)/I_0(k R)] / [1 - 2 I_1(k R)/(k R I_0(k R))],
    k^2         = i omega rho / eta*(omega),
    sigma_xr    = eta_p Ubar (-4 r/R^2) + 2 Re[s_1 e^{i omega t}],
    s_1         = eta_p g_1/(1 + i omega lambda),  g_1 = Ubar c_1 phi_1'(r),

the inflow of ArterialLesion{3D,Axi} (PulsatilePipeInflow), which is
therefore exact for Oldroyd-B. The normal stress solves the linear ODE

    sigma_xx + lambda (d_t sigma_xx - 2 gamma sigma_xr) = 0,

whose periodic solution has the harmonics 0, 1, 2 of the product gamma sigma_xr:
with gamma = g_0 + 2 Re[g_1 e^{i w t}], sigma_xr = s_0 + 2 Re[s_1 e^{i w t}],

    sigma_xx = 2 lambda { P_0 + 2 Re[P_1 e^{i w t}/(1 + i w lambda)]
                              + 2 Re[P_2 e^{2 i w t}/(1 + 2 i w lambda)] },
    P_0 = g_0 s_0 + 2 Re(g_1 conj(s_1)),  P_1 = g_0 s_1 + g_1 s_0,  P_2 = g_1 s_1.

sigma_rr = sigma_tt = 0. Run as a script for a self-check (momentum and
constitutive residuals by finite differences in r and t).
"""
import math

import numpy as np

try:
    from scipy.special import iv as _iv
except Exception:  # scipy is optional: power series, fine for |z| < ~15
    _iv = None


def bessel_i(nu, z):
    z = np.asarray(z, dtype=complex)
    if _iv is not None:
        return _iv(nu, z)
    if np.max(np.abs(z)) > 15.0:
        raise RuntimeError("|k R| > 15: install scipy for the Bessel functions")
    out = np.zeros_like(z)
    term = (z / 2.0) ** nu / math.factorial(nu)
    for k in range(200):
        out = out + term
        term = term * (z / 2.0) ** 2 / ((k + 1) * (k + 1 + nu))
    return out


class WomersleyOB:
    def __init__(self, radius, rho, eta_s, eta_p, lam, omega, ubar, amplitude):
        self.R = float(radius)
        self.rho, self.eta_s, self.eta_p = float(rho), float(eta_s), float(eta_p)
        self.lam, self.w, self.U, self.A = float(lam), float(omega), float(ubar), float(amplitude)
        eta_star = self.eta_s + self.eta_p / (1.0 + 1j * self.w * self.lam)
        self.k = np.sqrt(1j * self.w * self.rho / eta_star)
        kR = self.k * self.R
        self.i0R = bessel_i(0, kR)
        self.C = 1.0 / (1.0 - 2.0 * bessel_i(1, kR) / (kR * self.i0R))
        self.c1 = -0.5j * self.A

    # complex amplitudes of the first harmonic
    def _phi1(self, r):
        return self.C * (1.0 - bessel_i(0, self.k * r) / self.i0R)

    def _dphi1(self, r):
        return -self.C * self.k * bessel_i(1, self.k * r) / self.i0R

    def _g(self, r):
        r = np.asarray(r, dtype=float)
        g0 = self.U * (-4.0 * r / self.R ** 2)
        g1 = self.U * self.c1 * self._dphi1(r)
        return g0, g1

    def _s(self, r):
        g0, g1 = self._g(r)
        return self.eta_p * g0, self.eta_p * g1 / (1.0 + 1j * self.w * self.lam)

    def flow_rate_factor(self, t):
        return 1.0 + self.A * math.sin(self.w * t)

    def u(self, r, t):
        r = np.asarray(r, dtype=float)
        e = np.exp(1j * self.w * t)
        return self.U * (2.0 * (1.0 - r ** 2 / self.R ** 2)
                         + 2.0 * np.real(self.c1 * self._phi1(r) * e))

    def shear_rate(self, r, t):
        g0, g1 = self._g(r)
        return g0 + 2.0 * np.real(g1 * np.exp(1j * self.w * t))

    def sigma_xr(self, r, t):
        s0, s1 = self._s(r)
        return s0 + 2.0 * np.real(s1 * np.exp(1j * self.w * t))

    def sigma_xx(self, r, t):
        g0, g1 = self._g(r)
        s0, s1 = self._s(r)
        wl = self.w * self.lam
        P0 = g0 * s0 + 2.0 * np.real(g1 * np.conj(s1))
        P1 = g0 * s1 + g1 * s0
        P2 = g1 * s1
        e = np.exp(1j * self.w * t)
        return 2.0 * self.lam * (P0 + 2.0 * np.real(P1 * e / (1.0 + 1j * wl))
                                 + 2.0 * np.real(P2 * e * e / (1.0 + 2j * wl)))

    def pressure_gradient(self, t, dr=1e-7):
        """G(t) = -dp/dx from the momentum balance at the wall-adjacent point."""
        # rho d_t u = G + (1/r) d_r(r tau), tau = eta_s gamma + sigma_xr; on the axis
        # (1/r) d_r(r tau) = 2 d_r tau.
        r = np.array([dr])
        tau = lambda rr: self.eta_s * self.shear_rate(rr, t) + self.sigma_xr(rr, t)
        dtau = (tau(2 * r) - tau(r)) / dr
        dt = 1e-6 / self.w
        dudt = (self.u(0.0, t + dt) - self.u(0.0, t - dt)) / (2 * dt)
        return float(self.rho * dudt - 2.0 * dtau[0])


def _self_check():
    R, rho, es, ep = 2e-3, 1060.0, 1.795e-3, 1.795e-3
    eta0 = es + ep
    D = 2 * R
    U = 10 * eta0 / (rho * D)
    w = 4 * 9 * eta0 / (rho * D * D)
    lam = 1.0 / w
    s = WomersleyOB(R, rho, es, ep, lam, w, U, 0.5)
    r = np.linspace(0.05 * R, 0.95 * R, 7)
    t = 0.37 * 2 * math.pi / w
    h, dt = 1e-6 * R, 1e-5 / w
    # flow rate
    rr = np.linspace(0, R, 4001)
    Q = np.trapezoid(2 * np.pi * rr * s.u(rr, t), rr) if hasattr(np, "trapezoid") else np.trapz(2 * np.pi * rr * s.u(rr, t), rr)
    print("flow rate  rel err %.2e" % (Q / (U * np.pi * R * R * s.flow_rate_factor(t)) - 1))
    # constitutive residuals
    sxr = s.sigma_xr(r, t)
    dsxr = (s.sigma_xr(r, t + dt) - s.sigma_xr(r, t - dt)) / (2 * dt)
    gam = (s.u(r + h, t) - s.u(r - h, t)) / (2 * h)
    print("shear rate vs d_r u  rel err %.2e" % (np.max(np.abs(gam - s.shear_rate(r, t))) / np.max(np.abs(gam))))
    res1 = sxr + lam * dsxr - ep * gam
    sxx = s.sigma_xx(r, t)
    dsxx = (s.sigma_xx(r, t + dt) - s.sigma_xx(r, t - dt)) / (2 * dt)
    res2 = sxx + lam * (dsxx - 2 * gam * sxr)
    print("Maxwell xr residual %.2e Pa,  xx residual %.2e Pa (scales %.2e, %.2e)"
          % (np.max(np.abs(res1)), np.max(np.abs(res2)), np.max(np.abs(sxr)), np.max(np.abs(sxx))))
    # momentum: rho d_t u - (1/r) d_r(r tau) must not depend on r
    tau = lambda q: es * s.shear_rate(q, t) + s.sigma_xr(q, t)
    dudt = (s.u(r, t + dt) - s.u(r, t - dt)) / (2 * dt)
    div = ((r + h) * tau(r + h) - (r - h) * tau(r - h)) / (2 * h * r)
    G = rho * dudt - div
    print("momentum: G(r) spread %.2e of |G| %.2e" % (np.ptp(G), np.max(np.abs(G))))


if __name__ == "__main__":
    _self_check()

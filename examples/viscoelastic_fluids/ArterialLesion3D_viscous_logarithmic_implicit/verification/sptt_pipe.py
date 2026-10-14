#!/usr/bin/env python3
"""Fully developed pipe flow of a linear simplified Phan-Thien--Tanner (sPTT)
fluid with a Newtonian solvent: the reference solution of the verification test.

Model (the one of ArterialLesion{3D,Axi}_viscous_logarithmic_implicit with
lambda0 = lambda):
    T = 2 eta_s D + sigma,
    f(tr sigma) sigma + lambda sigma^nabla = 2 eta_p D,
    f = 1 + (epsilon lambda / eta_p) tr sigma.
For u = u(r) e_x and gamma = du/dr, fully developed and steady, the only
non-zero polymer stresses are
    sigma_xr = s,   sigma_xx = 2 lambda s^2 / eta_p,   sigma_rr = sigma_tt = 0,
    f = 1 + a s^2,  a = 2 epsilon lambda^2 / eta_p^2,   eta_p gamma = f s,
and the axial momentum balance gives the total shear stress
    eta_s gamma + s = tau(r) = -G r / 2,   G = -dp/dx > 0.
Eliminating gamma, s solves the cubic (kappa = eta_s / eta_p)
    kappa a s^3 + (1 + kappa) s - tau = 0,
which is strictly increasing in s and has exactly one real root. Then
    u(r) = int_r^R (-gamma) dr',   Q = pi int_0^R r^2 (-gamma) dr.
Limits: epsilon = 0 is Oldroyd-B (Poiseuille with eta_0 = eta_s + eta_p);
eta_s = 0 is the closed form of Oliveira & Pinho (1999),
    u = G (R^2 - r^2)/(4 eta_p) + epsilon lambda^2 G^3 (R^4 - r^4)/(16 eta_p^3).

Usage as a module:  sol = SPTTPipe(R, eta_s, eta_p, lam, eps); sol.set_flow_rate(Q)
then sol.u(r), sol.sigma_xr(r), sol.sigma_xx(r), sol.G, sol.wall_shear_stress.
Run as a script for a self-check of the limits.
"""
import math

import numpy as np


class SPTTPipe:
    def __init__(self, radius, eta_s, eta_p, lam, eps, n_quad=4001):
        self.R = float(radius)
        self.eta_s = float(eta_s)
        self.eta_p = float(eta_p)
        self.lam = float(lam)
        self.eps = float(eps)
        self.a = 2.0 * self.eps * self.lam ** 2 / self.eta_p ** 2
        self.kappa = self.eta_s / self.eta_p
        self.G = None
        # Gauss-Legendre nodes on [0, 1] for the quadratures.
        x, w = np.polynomial.legendre.leggauss(64)
        self._gx = 0.5 * (x + 1.0)
        self._gw = 0.5 * w

    # ---- local relations ----------------------------------------------------
    def shear_polymer_stress(self, tau):
        """s = sigma_xr from the cubic kappa a s^3 + (1 + kappa) s - tau = 0."""
        tau = np.asarray(tau, dtype=float)
        k, a = self.kappa, self.a
        if a == 0.0:
            return tau / (1.0 + k)
        if k == 0.0:
            return tau.copy()
        # Newton from the Oldroyd-B value tau/(1+k): the cubic is odd, increasing
        # and convex for s > 0 (concave for s < 0), so the iteration converges
        # monotonically from that starting point.
        s = tau / (1.0 + k)
        for _ in range(60):
            g = k * a * s ** 3 + (1.0 + k) * s - tau
            dg = 3.0 * k * a * s ** 2 + (1.0 + k)
            ds = g / dg
            s = s - ds
            if np.all(np.abs(ds) <= 1e-15 * (np.abs(s) + 1e-300)):
                break
        return s

    def shear_rate(self, r):
        """gamma = du/dr (negative for G > 0)."""
        tau = -0.5 * self.G * np.asarray(r, dtype=float)
        s = self.shear_polymer_stress(tau)
        return (1.0 + self.a * s ** 2) * s / self.eta_p

    # ---- profiles -----------------------------------------------------------
    def u(self, r):
        """u(r) = int_r^R (-gamma) dr' by 64-point Gauss-Legendre."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        out = np.empty_like(r)
        for i, ri in enumerate(r.flat):
            ri = min(max(ri, 0.0), self.R)
            rr = ri + (self.R - ri) * self._gx
            out.flat[i] = -(self.R - ri) * np.dot(self._gw, self.shear_rate(rr))
        return out

    def flow_rate_for(self, G):
        old = self.G
        self.G = G
        rr = self.R * self._gx
        Q = -math.pi * self.R * np.dot(self._gw, rr ** 2 * self.shear_rate(rr))
        self.G = old
        return Q

    def set_flow_rate(self, Q):
        """G such that the flow rate is Q (bisection on a monotone map)."""
        lo, hi = 0.0, 1.0
        while self.flow_rate_for(hi) < Q:
            hi *= 2.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self.flow_rate_for(mid) < Q:
                lo = mid
            else:
                hi = mid
            if hi - lo <= 1e-15 * hi:
                break
        self.G = 0.5 * (lo + hi)
        return self.G

    def sigma_xr(self, r):
        return self.shear_polymer_stress(-0.5 * self.G * np.asarray(r, dtype=float))

    def sigma_xx(self, r):
        s = self.sigma_xr(r)
        return 2.0 * self.lam * s ** 2 / self.eta_p

    @property
    def wall_shear_stress(self):
        """|tau_w| = G R / 2 (total, solvent + polymer)."""
        return 0.5 * self.G * self.R


def _self_check():
    R, ep = 2e-3, 1.8e-3
    # Oldroyd-B: Poiseuille with eta_0.
    s = SPTTPipe(R, 1.8e-3, ep, 0.05, 0.0)
    Q = 1e-7
    s.set_flow_rate(Q)
    G_ex = 8.0 * (3.6e-3) * Q / (math.pi * R ** 4)
    r = np.linspace(0, R, 7)
    u_ex = 2.0 * Q / (math.pi * R ** 2) * (1 - (r / R) ** 2)
    print("Oldroyd-B  G rel err %.2e   u max rel err %.2e"
          % (s.G / G_ex - 1, np.max(np.abs(s.u(r) - u_ex)) / u_ex[0]))
    # eta_s = 0: Oliveira & Pinho closed form.
    lam, eps = 0.05, 0.25
    s = SPTTPipe(R, 0.0, ep, lam, eps)
    s.G = 300.0
    u_op = (s.G * (R ** 2 - r ** 2) / (4 * ep)
            + eps * lam ** 2 * s.G ** 3 * (R ** 4 - r ** 4) / (16 * ep ** 3))
    Q_op = math.pi * (s.G * R ** 4 / (8 * ep) + eps * lam ** 2 * s.G ** 3 * R ** 6 / (24 * ep ** 3))
    print("Oliveira-Pinho  u max rel err %.2e   Q rel err %.2e"
          % (np.max(np.abs(s.u(r) - u_op)) / u_op[0], s.flow_rate_for(s.G) / Q_op - 1))
    # General case: cubic residual and momentum balance.
    s = SPTTPipe(R, 1.8e-3, ep, lam, eps)
    s.set_flow_rate(Q)
    rr = np.linspace(0, R, 11)
    tau = -0.5 * s.G * rr
    res = s.eta_s * s.shear_rate(rr) + s.sigma_xr(rr) - tau
    print("sPTT + solvent  momentum residual %.2e Pa,  Q rel err %.2e,  G/G_OB = %.4f"
          % (np.max(np.abs(res)), s.flow_rate_for(s.G) / Q - 1, s.G / G_ex))


if __name__ == "__main__":
    _self_check()

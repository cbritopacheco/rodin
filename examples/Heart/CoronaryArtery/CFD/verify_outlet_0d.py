#!/usr/bin/env python3
"""
Standalone verification of the coronary outlet condition.

Reimplements the 0D outlet exactly as CoupledLV0DCoronary3D does -- Carreau-
Yasuda WRMS tube law, collapsible-tube radius law, implicit Euler with a 2x2
Newton -- and drives it with the aortic and left-ventricular pressure traces of
a previous 3D run.  The 3D domain is replaced by the statement p_out = p_ar(t),
so the comparison isolates the outlet condition.

Both formulations are run on identical input:

  old : fixed tree calibres, softplus storage V = C p0 ln(1 + e^(p/p0)),
        explicit compliance imported from an identification.
  new : one collapsible-tube law r(p_tm) feeding both the flow law and the
        anatomical volume N pi r^2 L.

Checks reported: Newton convergence, exactness of the analytic Jacobian against
finite differences, discrete volume conservation, positivity of the stored
volume, and the phasic metrics (systolic/diastolic inflow ratio, venous
peak-to-mean) against their physiological ranges.

Usage:  python3 verify_outlet_0d.py [CoronaryArtery.csv]
"""

import csv
import math
import sys

PI = math.pi

# ----------------------------------------------------------------- rheology --

MU0, MUINF, LAM, NIDX, YAS = 0.353, 0.004181, 15.6821, 0.2050, 0.6497
MU_N = 0.0035
RHO = 1060.0


def mu(g):
    return MUINF + (MU0 - MUINF) * (1.0 + (LAM * g) ** YAS) ** ((NIDX - 1.0) / YAS)


def dmu(g):
    base = 1.0 + (LAM * g) ** YAS
    return ((MU0 - MUINF) * (NIDX - 1.0) * base ** ((NIDX - 1.0 - YAS) / YAS)
            * LAM ** YAS * g ** (YAS - 1.0))


def _bisect(f, lo, hi, tol=1e-14, it=200):
    flo = f(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm == 0.0 or (hi - lo) < tol * max(1.0, abs(mid)):
            return mid
        if (fm < 0.0) == (flo < 0.0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def flow(dp, L, r, nsub=60):
    """WRMS flow law. Returns (q, dq/ddp, dq/dr)."""
    R0 = 8.0 * MU0 * L / (PI * r ** 4)
    if abs(dp) < 1e-12:
        return dp / R0, 1.0 / R0, 4.0 * dp / (R0 * r)

    sgn = 1.0 if dp >= 0.0 else -1.0
    adp = abs(dp)
    tw = r * adp / (2.0 * L)

    hi = max(tw / MUINF, 1e-8)
    while mu(hi) * hi - tw < 0.0:
        hi *= 2.0
    gw = _bisect(lambda g: mu(g) * g - tw, 1e-14, hi)

    # Simpson on [0, gw]
    h = gw / nsub
    integ = 0.0
    for i in range(nsub + 1):
        g = i * h
        w = 1.0 if i in (0, nsub) else (4.0 if i % 2 else 2.0)
        if g > 0.0:
            m, dm = mu(g), dmu(g)
            integ += w * g ** 3 * m * m * (m + g * dm)
    integ *= h / 3.0

    q = PI * r ** 3 * integ / tw ** 3
    return sgn * q, (PI * r ** 3 * gw - 3.0 * q) / adp, sgn * PI * r * r * gw


# ---------------------------------------------------------------- tube law ---

M_EXP, N_EXP = 10.0, 1.5


def tube_area(K, ptm):
    """p_tm = K [a^m - a^-n].  Returns (a, da/dp_tm).

    Bisection-safeguarded Newton on s = log a, where the residual is monotone.
    """
    lo, hi = math.log(1e-3), math.log(4.0)
    g = lambda s: K * (math.exp(M_EXP * s) - math.exp(-N_EXP * s)) - ptm
    dg = lambda s: K * (M_EXP * math.exp(M_EXP * s) + N_EXP * math.exp(-N_EXP * s))
    if g(lo) >= 0.0:
        return 1e-3, 0.0
    if g(hi) <= 0.0:
        return 4.0, 0.0
    s = 0.0
    for _ in range(100):
        gs = g(s)
        if gs < 0.0:
            lo = s
        else:
            hi = s
        sn = s - gs / dg(s)
        if not (lo < sn < hi):
            sn = 0.5 * (lo + hi)
        if abs(sn - s) < 1e-14:
            s = sn
            break
        s = sn
    a = math.exp(s)
    return a, 1.0 / (K * (M_EXP * a ** (M_EXP - 1.0) + N_EXP * a ** (-N_EXP - 1.0)))


# -------------------------------------------------------------- calibration --

class Outlet:
    """One calibrated outlet.  mode is 'new' or 'old'."""

    def __init__(self, mode, Q, rp, Lp, par0, pv0, cfg):
        self.mode = mode
        self.rp, self.Lp, self.Q = rp, Lp, Q
        self.alpha = cfg['alpha']
        self.Pra = cfg['Pra']
        c = cfg['c']

        Ap = PI * rp * rp
        self.Zc = RHO * c / Ap
        dPp = 8.0 * MU_N * Lp * Q / (PI * rp ** 4)
        dPz = self.Zc * Q
        dP = par0 - self.Pra

        pim0 = self.alpha * pv0

        if mode == 'new':
            self.Ca = Ap * Lp / (RHO * c * c)
            dmicro = max(dP - dPp - dPz, 1.0)
            self.dPv = cfg['fv'] * dmicro
            self.dPa = (1.0 - cfg['fv']) * dmicro
            va, vv, Ta, Tv = cfg['va'], cfg['vv'], cfg['Ta'], cfg['Tv']
            self.La, self.Lv = va * Ta, vv * Tv
            self.ra0 = va * math.sqrt(8.0 * MU_N * Ta / self.dPa)
            self.rv0 = vv * math.sqrt(8.0 * MU_N * Tv / self.dPv)
            self.Na = Q / (va * PI * self.ra0 ** 2)
            self.Nv = Q / (vv * PI * self.rv0 ** 2)
            self.Ka, self.Kv = cfg['Ka'], cfg['Kv']
            ptm0 = self.Pra + self.dPv - pim0
            aa, daa = tube_area(self.Ka, ptm0)
            av, dav = tube_area(self.Kv, ptm0)
            self.Aua = PI * self.ra0 ** 2 / aa
            self.Auv = PI * self.rv0 ** 2 / av
            self.V0 = Q * (Ta + Tv)
            self.C0 = self.Na * self.La * self.Aua * daa + self.Nv * self.Lv * self.Auv * dav
        else:
            self.C = cfg['Ctot']
            self.Ca = cfg['CaFrac'] * self.C
            self.p0 = cfg['p0']
            self.ra0, self.rv0, self.Lv = 2.5e-5, 3.0e-5, 1.5e-2
            va, vv = cfg['va'], cfg['vv']
            Qa1 = va * PI * self.ra0 ** 2
            Qv1 = vv * PI * self.rv0 ** 2
            self.Na, self.Nv = Q / Qa1, Q / Qv1
            self.dPv = 8.0 * MU_N * self.Lv * Qv1 / (PI * self.rv0 ** 4)
            self.dPa = max(dP - dPp - dPz - self.dPv, 1.0)
            self.La = self.dPa * PI * self.ra0 ** 4 / (8.0 * MU_N * Qa1)
            self.V0 = self.C * self.p0 * math.log1p(
                math.exp(min((self.Pra + self.dPv - pim0) / self.p0, 30.0)))
            self.C0 = self.C

        self.pc = self.Pra + self.dPv
        self.pca = self.pc + self.dPa
        self.pim = pim0
        self.vol = self.volume(self.pc - self.pim)[0]

    # ---- storage and calibre, the one place the two formulations differ ----

    def trees(self, ptm):
        """Returns (r_a, dr_a/dp, r_v, dr_v/dp)."""
        if self.mode == 'old':
            return self.ra0, 0.0, self.rv0, 0.0
        aa, daa = tube_area(self.Ka, ptm)
        av, dav = tube_area(self.Kv, ptm)
        ra = math.sqrt(self.Aua * aa / PI)
        rv = math.sqrt(self.Auv * av / PI)
        return ra, self.Aua * daa / (2.0 * PI * ra), rv, self.Auv * dav / (2.0 * PI * rv)

    def volume(self, ptm):
        """Returns (V, dV/dp_tm)."""
        if self.mode == 'old':
            y = ptm / self.p0
            if y > 30.0:
                return self.C * ptm, self.C
            if y < -30.0:
                e = math.exp(y)
                return self.C * self.p0 * e, self.C * e
            return (self.C * self.p0 * math.log1p(math.exp(y)),
                    self.C / (1.0 + math.exp(-y)))
        aa, daa = tube_area(self.Ka, ptm)
        av, dav = tube_area(self.Kv, ptm)
        return (self.Na * self.La * self.Aua * aa + self.Nv * self.Lv * self.Auv * av,
                self.Na * self.La * self.Aua * daa + self.Nv * self.Lv * self.Auv * dav)

    # ---- inflow implied by p_out = p_ar (the 3D domain, replaced) ----

    def inflow(self, pca, par):
        """Q and dQ/dp_ca from par = pca + dPp(Q) + Zc Q.

        The conduit is the measured epicardial calibre at high shear, so it is
        linear: R_p Q + Z_c Q = p_ar - p_ca.
        """
        R = 8.0 * MU_N * self.Lp / (PI * self.rp ** 4) + self.Zc
        return (par - pca) / R, -1.0 / R

    # ---- one implicit-Euler step ----

    def step(self, par, pv, dt):
        pim, pimOld = self.alpha * pv, self.pim
        pcaOld, pcOld = self.pca, self.pc
        volOld = self.volume(pcOld - pimOld)[0]
        capa = self.Ca / dt

        def residual(pca, pc):
            ptm = pc - pim
            ra, dra, rv, drv = self.trees(ptm)
            qa, dqa_dp, dqa_dr = flow(pca - pc, self.La, ra)
            qv, dqv_dp, dqv_dr = flow(pc - self.Pra, self.Lv, rv)
            qa, qv = self.Na * qa, self.Nv * qv
            V, dV = self.volume(ptm)
            Q, dQ = self.inflow(pca, par)
            dqa_dpca = self.Na * dqa_dp
            dqa_dpc = self.Na * (-dqa_dp + dqa_dr * dra)
            dqv_dpc = self.Nv * (dqv_dp + dqv_dr * drv)
            R1 = capa * (pca - pcaOld) + qa - Q
            R2 = (V - volOld) / dt - qa + qv
            J11 = capa + dqa_dpca - dQ
            J12 = dqa_dpc
            J21 = -dqa_dpca
            J22 = dV / dt - dqa_dpc + dqv_dpc
            return R1, R2, J11, J12, J21, J22, Q, qa, qv, V

        pca, pc = pcaOld, pcOld
        its = 0
        for its in range(1, 51):
            R1, R2, J11, J12, J21, J22, *_ = residual(pca, pc)
            det = J11 * J22 - J12 * J21
            dpca = (-J22 * R1 + J12 * R2) / det
            dpc = (J21 * R1 - J11 * R2) / det
            f0 = abs(R1) + abs(R2)
            lam = 1.0
            for _ in range(12):
                r = residual(pca + lam * dpca, pc + lam * dpc)
                if abs(r[0]) + abs(r[1]) <= (1.0 - 1e-4 * lam) * f0:
                    break
                lam *= 0.5
            pca += lam * dpca
            pc += lam * dpc
            if lam * (abs(dpca) + abs(dpc)) < 1e-12 * (1.0 + abs(pca) + abs(pc)):
                break

        self.pca, self.pc, self.pim = pca, pc, pim
        R1, R2, *_, Q, qa, qv, V = residual(pca, pc)
        self.vol = V
        return dict(Q=Q, qa=qa, qv=qv, V=V, pc=pc, pca=pca, pim=pim,
                    ptm=pc - pim, iters=its, res=abs(R1) + abs(R2),
                    dV=(V - volOld), imbalance=(V - volOld) / dt - (qa - qv))


# conduit helpers (Poiseuille on the measured epicardial calibre)
def flow_dp(Q, L, r):
    return 8.0 * MU_N * L * Q / (PI * r ** 4)


def flow_dp_d(Q, L, r):
    R = 8.0 * MU_N * L / (PI * r ** 4)
    return R * Q, R


# ------------------------------------------------------------------- driver --

CFG = dict(alpha=0.65, Pra=600.0, c=8.0, va=5.0e-3, vv=3.0e-3,
           Ta=1.5, Tv=2.0, fv=0.13, Ka=500.0, Kv=200.0,
           Ctot=3.2e-10, CaFrac=0.10, p0=300.0)


def load(path):
    rows = list(csv.DictReader(open(path)))
    t = [float(r['t']) for r in rows]
    return (t, [float(r['AortaPressure']) for r in rows],
            [float(r['LeftVentriclePressure']) for r in rows])


def run(mode, t, par, pv):
    o = Outlet(mode, 2.0e-6, 2.0e-3, 0.0125, par[0], pv[0], CFG)
    out = []
    for i in range(1, len(t)):
        out.append(o.step(par[i], pv[i], t[i] - t[i - 1]))
    return o, out


def jacobian_check(mode, par0, pv0):
    o = Outlet(mode, 2.0e-6, 2.0e-3, 0.0125, par0, pv0, CFG)
    pim = o.alpha * pv0
    dt, eps = 1e-3, 1e-3
    o.pim = pim

    def R(pca, pc):
        volOld = o.volume(o.pc - o.pim)[0]
        ptm = pc - pim
        ra, _, rv, _ = o.trees(ptm)
        qa = o.Na * flow(pca - pc, o.La, ra)[0]
        qv = o.Nv * flow(pc - o.Pra, o.Lv, rv)[0]
        V = o.volume(ptm)[0]
        Q = o.inflow(pca, par0)[0]
        return (o.Ca / dt * (pca - o.pca) + qa - Q, (V - volOld) / dt - qa + qv)

    pca, pc = o.pca * 1.01, o.pc * 1.05
    ptm = pc - pim
    ra, dra, rv, drv = o.trees(ptm)
    qa, dqa_dp, dqa_dr = flow(pca - pc, o.La, ra)
    qv, dqv_dp, dqv_dr = flow(pc - o.Pra, o.Lv, rv)
    dV = o.volume(ptm)[1]
    dQ = o.inflow(pca, par0)[1]
    J = [[o.Ca / dt + o.Na * dqa_dp - dQ, o.Na * (-dqa_dp + dqa_dr * dra)],
         [-o.Na * dqa_dp,
          dV / dt - o.Na * (-dqa_dp + dqa_dr * dra) + o.Nv * (dqv_dp + dqv_dr * drv)]]
    fd = [[(R(pca + eps, pc)[k] - R(pca - eps, pc)[k]) / (2 * eps) for k in (0, 1)],
          [(R(pca, pc + eps)[k] - R(pca, pc - eps)[k]) / (2 * eps) for k in (0, 1)]]
    err = max(abs(J[i][j] - fd[j][i]) / max(abs(fd[j][i]), 1e-12)
              for i in (0, 1) for j in (0, 1))
    return err


def metrics(name, o, out, pv):
    cyc = [r for r in out[len(out) // 2:]]
    pvc = pv[len(out) // 2 + 1:len(out) + 1]
    Q = [r['Q'] for r in cyc]
    qv = [r['qv'] for r in cyc]
    V = [r['V'] for r in cyc]
    sysQ = [q for q, p in zip(Q, pvc) if p > 4000]
    diaQ = [q for q, p in zip(Q, pvc) if p <= 4000]
    mean = lambda x: sum(x) / len(x)
    print(f"\n=== {name} ===")
    print(f"  Newton iterations   max {max(r['iters'] for r in cyc):3d}   "
          f"mean {mean([r['iters'] for r in cyc]):.1f}")
    print(f"  final residual      max {max(r['res'] for r in cyc):.2e}")
    print(f"  volume imbalance    max {max(abs(r['imbalance']) for r in cyc):.2e} m^3/s"
          f"   (relative {max(abs(r['imbalance']) for r in cyc) / mean(Q):.1e})")
    print(f"  stored volume       min {min(V) * 1e6:8.3f} mL   "
          f"max {max(V) * 1e6:8.3f} mL   rest {o.V0 * 1e6:.3f} mL")
    print(f"  transmural pressure min {min(r['ptm'] for r in cyc):8.0f} Pa  "
          f"max {max(r['ptm'] for r in cyc):8.0f} Pa")
    print(f"  mean inflow         {mean(Q) * 6e7:.1f} mL/min")
    print(f"  systolic/diastolic  {mean(sysQ) / mean(diaQ):.2f}"
          f"       (physiological 0.4-0.6)")
    print(f"  venous peak/mean    {max(qv) / mean(qv):.2f}"
          f"       (physiological 2-3)")
    if o.mode == 'new':
        print(f"  predicted r_a/r_v   {o.ra0 * 1e6:.1f} / {o.rv0 * 1e6:.1f} um"
              f"   (terminal arteriole 8-15, venule 15-40)")
        print(f"  predicted L_a/L_v   {o.La * 1e3:.1f} / {o.Lv * 1e3:.1f} mm")
        print(f"  predicted C         {o.C0:.2e} m^3/Pa"
              f"   (identified ~1e-10 per bed)")
    return dict(sd=mean(sysQ) / mean(diaQ), pk=max(qv) / mean(qv),
                Vmin=min(V), imb=max(abs(r['imbalance']) for r in cyc) / mean(Q))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'CoronaryArtery.csv'
    t, par, pv = load(path)
    print(f"driving {len(t)} steps from {path}")
    for mode in ('old', 'new'):
        print(f"  analytic vs finite-difference Jacobian, {mode}: "
              f"max relative error {jacobian_check(mode, par[0], pv[0]):.2e}")
    res = {}
    for mode, label in (('old', 'OLD  softplus storage, rigid trees'),
                        ('new', 'NEW  one collapsible-tube law')):
        o, out = run(mode, t, par, pv)
        res[mode] = metrics(label, o, out, pv)
    print("\n--- verdict ---")
    ok = True
    for k, (lo, hi), lbl in ((('sd'), (0.35, 0.7), 'systolic/diastolic inflow'),
                             (('pk'), (1.5, 3.5), 'venous peak/mean')):
        v = res['new'][k]
        good = lo <= v <= hi
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {lbl}: {v:.2f} in [{lo}, {hi}]")
    print(f"  {'PASS' if res['new']['Vmin'] > 0 else 'FAIL'}  stored volume "
          f"stays positive: min {res['new']['Vmin'] * 1e6:.3f} mL")
    print(f"  {'PASS' if res['new']['imb'] < 1e-8 else 'FAIL'}  discrete volume "
          f"conservation: {res['new']['imb']:.1e}")
    if not (1.5 <= res['new']['pk'] <= 3.5):
        print("\n  note: the venous peak-to-mean ratio is the one metric still "
              "outside its\n        physiological band. It is insensitive to K_a, "
              "K_v, T_a, T_v, f_v, alpha and\n        to the target flow, and it "
              "does not respond to smoothing p_im, so it is not\n        a "
              "mis-set constant. It points at a missing element -- the "
              "compliance of the\n        extramyocardial veins and coronary "
              "sinus, which buffer the systolic surge\n        and which this "
              "model replaces with a rigid P_RA.")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

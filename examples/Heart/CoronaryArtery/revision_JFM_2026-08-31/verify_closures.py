#!/usr/bin/env python3
"""Independent verification of the analytical WRMS closures in the manuscript.

Checks performed
----------------
1. Cross model: I_Cr from eq. (I_Cross_expanded) [hypergeometric J_k] vs
   direct adaptive quadrature of I = int_0^{gw} g^3 mu^2 (mu + g mu') dg.
2. Carreau-Yasuda: I_CY from eq. (I_CY_final) vs the same quadrature.
3. Newtonian limits (mu0 -> muInf; PL n=1).
4. Quemada: R_mu^Q via F(alpha,q) with the P_j polynomials vs quadrature.
5. Exact tangent dQ/dDp vs centred finite differences.
6. Recomputation of the calibration-error table (Table 2) and the
   transition shear rates gamma_tr for eps=0.1, 0.2.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import hyp2f1
from scipy.optimize import brentq

# ---------------- parameter sets from Table 1 (SI units) --------------------
CROSS = {
 'healthy':      dict(mu0=78.9514e-3, muInf=5.0629e-3, lam=4.1199, n=0.9503),
 'diabetes':     dict(mu0=103.5405e-3, muInf=5.3610e-3, lam=10.4591, n=0.8633),
 'hyperviscous': dict(mu0=295.7832e-3, muInf=5.0291e-3, lam=13.7952, n=0.7982),
}
CY = {
 'healthy':      dict(mu0=71.9637e-3, muInf=4.8612e-3, lam=6.1529, a=1.4173, n=0.2081),
 'diabetes':     dict(mu0=104.6692e-3, muInf=5.3610e-3, lam=11.1048, a=0.8625, n=0.1502),
 'hyperviscous': dict(mu0=319.4782e-3, muInf=5.5066e-3, lam=16.9213, a=0.7709, n=0.2152),
}
# Quemada: table gives muF and gamma_c directly; k0,kInf columns appear to be
# mu0, muInf in mPa s -> invert eq (2.16) for k0, kInf.
QUEM = {}
for state, phi, mu0t, muInft, gc, muF in [
    ('healthy', 0.420, 80.3102, 5.0041, 2.2922, 1.7963),
    ('diabetes', 0.450, 122.2013, 4.7551, 1.7278, 1.9115),
    ('hyperviscous', 0.517, 313.0808, 2.6863, 14.9623, 1.3228)]:
    k0 = 2.0 * (1.0 - np.sqrt(muF / mu0t)) / phi
    kInf = 2.0 * (1.0 - np.sqrt(muF / muInft)) / phi
    QUEM[state] = dict(phi=phi, muF=muF*1e-3, gc=gc, k0=k0, kInf=kInf)

# ---------------- generic WRMS quadrature -----------------------------------
def I_quad(mu, dmu, gw):
    f = lambda g: g**3 * mu(g)**2 * (mu(g) + g*dmu(g))
    val, err = quad(f, 0.0, gw, limit=400, epsabs=1e-300, epsrel=1e-12)
    return val

def gamma_wall(mu, tauW):
    f = lambda g: mu(g)*g - tauW
    lo, hi = 1e-14, tauW/1e-4 + 1e3
    while f(hi) < 0: hi *= 10
    return brentq(f, lo, hi, xtol=1e-300, rtol=1e-15)

# ---------------- Cross ------------------------------------------------------
def cross_funcs(p):
    d = p['mu0'] - p['muInf']; lam, n = p['lam'], p['n']
    mu  = lambda g: p['muInf'] + d/(1.0 + (lam*g)**n)
    dmu = lambda g: -d*n*(lam*g)**n/(g*(1.0+(lam*g)**n)**2) if g > 0 else 0.0
    return mu, dmu

def I_cross_analytic(p, gw):
    d = p['mu0'] - p['muInf']; lam, n, mi = p['lam'], p['n'], p['muInf']
    J = lambda k: gw**4/4.0 * hyp2f1(k, 4.0/n, 1.0 + 4.0/n, -(lam*gw)**n)
    return (mi**3*J(0) + mi**2*d*(3-n)*J(1)
            + (mi*d**2*(3-2*n) + mi**2*d*n)*J(2)
            + (d**3*(1-n) + 2*mi*d**2*n)*J(3) + d**3*n*J(4))

# ---------------- Carreau-Yasuda --------------------------------------------
def cy_funcs(p):
    d = p['mu0'] - p['muInf']; lam, a, n = p['lam'], p['a'], p['n']
    mu  = lambda g: p['muInf'] + d*(1.0 + (lam*g)**a)**((n-1.0)/a)
    dmu = lambda g: (d*(n-1.0)*(1.0+(lam*g)**a)**((n-1.0-a)/a)
                     * lam**a * g**(a-1.0)) if g > 0 else 0.0
    return mu, dmu

def I_cy_analytic(p, gw):
    d = p['mu0'] - p['muInf']; lam, a, n, mi = p['lam'], p['a'], p['n'], p['muInf']
    b = (n-1.0)/a
    J = lambda k: gw**4/4.0 * hyp2f1(-k, 4.0/a, 1.0 + 4.0/a, -(lam*gw)**a)
    return (mi**3*J(0) + mi**2*d*(n+2)*J(b) - mi**2*d*(n-1)*J(b-1)
            + mi*d**2*(2*n+1)*J(2*b) - 2*mi*d**2*(n-1)*J(2*b-1)
            + d**3*n*J(3*b) - d**3*(n-1)*J(3*b-1))

# ---------------- Quemada ----------------------------------------------------
def quem_funcs(p):
    phi, muF, gc, k0, kInf = p['phi'], p['muF'], p['gc'], p['k0'], p['kInf']
    def mu(g):
        s = np.sqrt(max(g, 0.0)/gc)
        k = (k0 + kInf*s)/(1.0 + s)
        return muF/(1.0 - 0.5*k*phi)**2
    def dmu(g):
        h = 1e-7*max(g, 1e-12)
        return (mu(g+h) - mu(max(g-h, 0.0)))/(2*h)
    return mu, dmu

def P_polys(q):
    P1 = -(q+8)/7.0
    P2 = -(13*q**2 - 8*q - 7)/42.0
    P3 = -(143*q**3 - 88*q**2 - 113*q + 48)/210.0
    P4 = -(1287*q**4 - 792*q**3 - 1342*q**2 + 632*q + 175)/840.0
    P5 = -(3003*q**5 - 1848*q**4 - 3894*q**3 + 1944*q**2 + 1011*q - 256)/840.0
    P6 = -(15015*q**6 - 9240*q**5 - 23331*q**4 + 12096*q**3 + 9081*q**2
           - 3179*q - 525)/1680.0
    P7 = -(45045*q**7 - 27720*q**6 - 82005*q**5 + 43680*q**4 + 42819*q**3
           - 17304*q**2 - 5619*q + 1024)/1680.0
    P8 = -(1.0/16.0)*(1-q)**2*(1+q)*(429*q**5 + 165*q**4 - 330*q**3
           - 90*q**2 + 45*q + 5)
    return [P1,P2,P3,P4,P5,P6,P7,P8]

def F_alpha_q(alpha, q):
    P = P_polys(q)
    root = np.sqrt(1.0 - 2.0*alpha*q + alpha**2)
    S = 1.0 + sum(alpha**(j+1)*P[j] for j in range(7))
    logarg = (1.0 - alpha*q + root)/(alpha*(1.0 - q))
    return 0.5*(1.0 - (8.0/7.0)*alpha*(1.0+q) + (4.0/3.0)*alpha**2
                - alpha**8*P[6] + S*root + alpha**8*P[7]*np.log(logarg))

def quem_F_analytic(p, tauW):
    phi, muF, gc, k0, kInf = p['phi'], p['muF'], p['gc'], p['k0'], p['kInf']
    muInf = muF/(1.0 - 0.5*kInf*phi)**2
    tau0 = muF*gc*(0.5*phi*(k0-kInf))**2/(1.0 - 0.5*kInf*phi)**4
    LamQ = gc*((1.0 - 0.5*k0*phi)/(1.0 - 0.5*kInf*phi))**2
    alpha = (np.sqrt(tau0) + np.sqrt(muInf*LamQ))/np.sqrt(tauW)
    q = (np.sqrt(tau0) - np.sqrt(muInf*LamQ))/(np.sqrt(tau0) + np.sqrt(muInf*LamQ))
    return F_alpha_q(alpha, q), muInf

# ---------------- run checks -------------------------------------------------
print("== 1/2: Cross and Carreau-Yasuda analytic vs quadrature ==")
tauWs = [1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.8, 5.0, 20.0, 100.0]
for name, sets, funcs, analytic in [("Cross", CROSS, cross_funcs, I_cross_analytic),
                                    ("CY", CY, cy_funcs, I_cy_analytic)]:
    worst = 0.0
    for state, p in sets.items():
        mu, dmu = funcs(p)
        for tw in tauWs:
            gw = gamma_wall(mu, tw)
            Ia = analytic(p, gw)
            Iq = I_quad(mu, dmu, gw)
            rel = abs(Ia-Iq)/abs(Iq)
            worst = max(worst, rel)
    print(f"  {name}: max relative discrepancy over 3 states x {len(tauWs)} tau_w = {worst:.2e}")

print("== 3: Newtonian limits ==")
pN = dict(mu0=3.5e-3, muInf=3.5e-3, lam=4.0, n=0.9)
mu, dmu = cross_funcs(pN)
gw = gamma_wall(mu, 2.0)
print(f"  Cross mu0=muInf: I_analytic/(tau_w^4/(4 mu)) = {I_cross_analytic(pN, gw)/(2.0**4/(4*3.5e-3)):.12f}")
pN2 = dict(mu0=3.5e-3, muInf=3.5e-3, lam=6.0, a=1.4, n=0.2)
mu, dmu = cy_funcs(pN2)
gw = gamma_wall(mu, 2.0)
print(f"  CY    mu0=muInf: I_analytic/(tau_w^4/(4 mu)) = {I_cy_analytic(pN2, gw)/(2.0**4/(4*3.5e-3)):.12f}")

print("== 4: Quemada F(alpha,q) vs quadrature ==")
worst = 0.0
for state, p in QUEM.items():
    mu, dmu = quem_funcs(p)
    print(f"  {state}: implied k0={p['k0']:.4f}, kInf={p['kInf']:.4f} "
          f"(admissibility 1-k0*phi/2 = {1-0.5*p['k0']*p['phi']:.4f})")
    for tw in [0.1, 0.5, 1.4, 2.8, 10.0]:
        gw = gamma_wall(mu, tw)
        Iq = I_quad(mu, dmu, gw)
        Fq_num = 4.0 * Iq / tw**4  # = F/muInf per I_Q = tau^4/(4 muInf) F
        Fa, muInf = quem_F_analytic(p, tw)
        rel = abs(Fa/muInf - Fq_num)/abs(Fq_num)
        worst = max(worst, rel)
print(f"  Quemada: max relative discrepancy = {worst:.2e}")

print("== 5: exact tangent vs centred FD (CY healthy, R=25um,L=1mm,N=1) ==")
p = CY['healthy']; mu, dmu = cy_funcs(p)
R, L = 25e-6, 1e-3
def Q_of_dp(dp):
    tw = R*dp/(2*L)
    gw = gamma_wall(mu, tw)
    return np.pi*R**3/tw**3 * I_quad(mu, dmu, gw)
for dp in [50.0, 500.0, 5000.0]:
    tw = R*dp/(2*L); gw = gamma_wall(mu, tw); Q = Q_of_dp(dp)
    exact = (np.pi*R**3*gw - 3*Q)/dp
    h = dp*1e-5
    fd = (Q_of_dp(dp+h) - Q_of_dp(dp-h))/(2*h)
    print(f"  dp={dp:7.1f} Pa: dQ/ddp exact={exact:.6e}, FD={fd:.6e}, rel={abs(exact-fd)/abs(fd):.2e}")

print("== 6: mu_eff table spot checks (CY fits) ==")
def mu_eff_of_gamma_a(p, ga_target):
    mu, dmu = cy_funcs(p)
    # find tau_w such that nominal shear tau_w/mu_ap = ga_target
    def resid(logtw):
        tw = np.exp(logtw)
        gw = gamma_wall(mu, tw)
        muap = tw**4/(4*I_quad(mu, dmu, gw))
        return tw/muap - ga_target
    lo, hi = np.log(1e-6), np.log(1e4)
    logtw = brentq(resid, lo, hi, xtol=1e-13)
    tw = np.exp(logtw); gw = gamma_wall(mu, tw)
    return tw**4/(4*I_quad(mu, dmu, gw))
for state in ['healthy', 'diabetes', 'hyperviscous']:
    p = CY[state]
    row = []
    for ga0 in [50, 100, 300, 800]:
        me0 = mu_eff_of_gamma_a(p, ga0)
        einf = p['muInf']/me0 - 1
        # e(q) for q/q0 = 1/10 and 2
        e110 = me0/mu_eff_of_gamma_a(p, ga0/10) - 1
        e2 = me0/mu_eff_of_gamma_a(p, ga0*2) - 1
        row.append(f"ga0={ga0}: mu_eff={me0*1e3:.2f} mPa.s, e_inf={einf*100:+.1f}%, e(1/10)={e110*100:+.1f}%, e(2)={e2*100:+.1f}%")
    print(f"  {state}:\n    " + "\n    ".join(row))

print("== 7: transition shear rates mu_eff(g_tr)=(1+eps) muInf ==")
for state in ['healthy', 'diabetes', 'hyperviscous']:
    p = CY[state]
    out = []
    for eps in [0.2, 0.1]:
        f = lambda lg: mu_eff_of_gamma_a(p, np.exp(lg)) - (1+eps)*p['muInf']
        lg = brentq(f, np.log(1e-2), np.log(5e3), xtol=1e-12)
        out.append(f"eps={eps}: gamma_tr={np.exp(lg):.1f} 1/s")
    print(f"  {state}: " + ", ".join(out))

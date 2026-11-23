"""Constitutive parameters and hydraulic curves for the outlet replay."""

import numpy as np
from scipy.integrate import cumulative_simpson


CY = {
    "normal": (71.2616e-3, 4.1000e-3, 9.2780, 1.8773, .3572),
    "diabetic": (97.3528e-3, 4.9622e-3, 15.7825, 1.0829, .2670),
    "hyperviscous": (346.7269e-3, 5.3264e-3, 17.4172, .6957, .2200),
}
T0, T1, PERIOD = 1.70, 2.55, .85
MU_REF, Q_REF = .0035, 1.5e-6
trap = getattr(np, "trapezoid", None) or np.trapz


def viscosity(g, state="normal"):
    m0, mi, lam, a, n = CY[state]
    return mi + (m0-mi) * (1+(lam*g)**a)**((n-1)/a)


def law_and_slope(g, law="CY", state="normal"):
    if law == "CY":
        m0, mi, lam, a, n = CY[state]
        z = (lam*g)**a
        m = mi + (m0-mi)*(1+z)**((n-1)/a)
        slope = m + (m0-mi)*(n-1)*z*(1+z)**((n-1)/a-1)
    elif law == "Cross":
        m0, mi, lam, n = .0977252, .0041000, 8.3756, .7209
        z = (lam*g)**n
        m = mi+(m0-mi)/(1+z)
        slope = m-(m0-mi)*n*z/(1+z)**2
    elif law == "Quemada":
        phi, mf, k0, ki, gc = .425, .0013, 4.2250, 2.0560, 3.4861
        s = np.sqrt(g/gc)
        den = 1-.5*phi*(k0+ki*s)/(1+s)
        m = mf/den**2
        slope = m+mf*.5*phi*(ki-k0)*s/(den**3*(1+s)**2)
    elif law == "Power law":
        m = .0213051*g**(.6667-1)
        slope = .6667*m
    else:
        raise ValueError(law)
    return m, slope


def hydraulic_curve(law="CY", state="normal"):
    """Integrate on a logarithmic shear grid, with the analytic low-end tail.

    The right-hand abscissa is the apparent rate, not the wall rate. The
    logarithmic grid resolves the low-stress plateau missed by old scripts.
    """
    g = np.geomspace(1e-10, 1e9, 45000)
    m, slope = law_and_slope(g, law, state)
    tau = m*g
    if law == "Power law":
        n, consistency = .6667, .0213051
        integral = consistency**3*n/(3*n+1)*g**(3*n+1)
    else:
        integral = cumulative_simpson(g**3*m**2*slope, x=g, initial=0)
        integral += m[0]**3*g[0]**4/4
    me = tau**4/(4*integral)
    ga = 4*integral/tau**3
    if not (np.all(np.isfinite(me)) and np.all(me > 0) and np.all(np.diff(ga)>0)):
        raise ValueError(f"Invalid hydraulic curve for {law}, {state}")
    return tau, ga, me

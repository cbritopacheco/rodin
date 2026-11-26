"""Shared settings of the PTAGX validation tests.

The atrial model of LeftAtrium2DPTAGX (P, T, A, X, G and immobile F) has the
same nodal kinetics as PTAGP with platelets off and theta_0 = 1, factor Xa on,
immobile fibrin without monomer or cap, and no thrombomodulin. The tests on
tissue-factor patches (E2-E5) run PTAGXBench: a copy of PTAGP with the
surface amplification layer of PTAGX (-pg_sigma_s, -pg_delta_s, -pg_patch),
which provides strips and channels with a reactive patch. The atrial benchmark
(Q1) runs LeftAtrium2DPTAGX itself.

Parameter sets:
  MODEL    the proposed model, on which every verdict is given: k_a from
           Thrombodynamics (E2), no surface amplification;
  SURFACE  MODEL plus the surface amplification sigma calibrated on Shen 2008
           (E5, calibrate_sigma.py); a sensitivity case (it breaks E2-E4);
  LEGACY   the PTAGX defaults (k_a = 0.10, no surface amplification).

Runner helpers (parser, run_many, read_csv, ...) are those of ../Test.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_GUESS = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BENCH = os.path.join(REPO_GUESS, "build", "examples", "Heart", "Thrombosis", "PTAGXBench",
                     "PTAGXBench")
os.environ.setdefault("PTAGP_EXE", BENCH)          # --exe of ../Test/ptagp_common
sys.path.insert(0, os.path.join(HERE, "..", "Test"))
import ptagp_common as pc  # noqa: E402,F401  (re-exported)

REPO = pc.REPO
EXE_X = os.environ.get(
    "PTAGX_EXE",
    os.path.join(REPO, "build", "examples", "Heart", "Thrombosis", "LeftAtrium2DPTAGX",
                 "LeftAtrium2DPTAGX"))

# Constants of LeftAtrium2DPTAGX.h (SI). Keep in sync.
K = dict(
    D_P=4.6e-11, D_G=2.0e-11,
    P0=1.4e-3, A0=2.8e-3, G0=7.35e-3, Pref=1.4e-3,
    k_a=0.10, K_A=1.0e-5, k_on=7.1, k_F=42.0, K_G=3.6e-3,
    kappa=3.0, k_dX=1.0e-2, K_F=2.0e-3,
)

# Proposed model. k_a: Thrombodynamics (E2, all four criteria at 0.035).
MODEL = dict(k_a=0.035, sigma_s=0.0, delta_s=10e-6, phi_i=0.0,
             kappa=10.0, eps_x=0.1, k_tm=1.0e-5, monomer=False)
# PTAGX before the changes (first round): no cofactor requirement, no
# thrombomodulin. The bench has no endothelium, so k_tm only acts in the atrium.
BASE = dict(k_a=0.035, sigma_s=0.0, delta_s=10e-6, phi_i=0.0,
            kappa=3.0, eps_x=1.0, k_tm=0.0, monomer=False)
# Sensitivity: surface amplification (option 2). sigma_s: Shen 2008 at 5 1/s
# (E5/calibrate_sigma.py): initiation switches on between 1.91e-3 (never) and
# 1.94e-3 m (120 s); 2.0e-3 is 4 % above it. delta_s: thickness of the layer
# that carries the surface membrane in the bench (in the atrium the endocardial
# layer delta = 0.5 mm carries it; sigma is the same constant).
SURFACE = dict(k_a=0.035, sigma_s=2.0e-3, delta_s=10e-6, phi_i=0.0,
               kappa=3.0, eps_x=1.0)
LEGACY = dict(k_a=0.10, sigma_s=0.0, delta_s=10e-6, phi_i=0.0,
              kappa=3.0, eps_x=1.0)
PLASMA = MODEL          # former name of MODEL

# PTAGP/PTAGXBench options that reduce it to the PTAGX kinetics.
PTAGX = dict(
    pg_platelets=False, pg_theta0=1.0, pg_xa=True,
    pg_immobile_fibrin=True, pg_fibrin_monomer=False, pg_F_cap=0,
    pg_tm="none", pg_contact=False, pg_RF_linear=False,
)


def consts(model):
    """0D constants (K) of a parameter set."""
    return dict(K, k_a=model["k_a"], kappa=model.get("kappa", 3.0),
                eps_x=model.get("eps_x", 1.0))


def kinetics(model, patch=None):
    """PTAGXBench options of a parameter set (MODEL, SURFACE or LEGACY)."""
    o = dict(PTAGX, pg_k_a=model["k_a"], pg_sigma_s=model["sigma_s"],
             pg_delta_s=model["delta_s"], pg_phi_i=model.get("phi_i", 0.0),
             pg_kappa_x=model.get("kappa", 3.0), pg_eps_x=model.get("eps_x", 1.0))
    if patch is not None:
        o["pg_patch"] = f"{patch[0]:.9g},{patch[1]:.9g}"
    return o


def verdict(out, test, checks, notes=None):
    """Print and store the pass/fail of each criterion (on MODEL).

    checks: list of (criterion, passed, detail); notes: informative lines
    (sensitivity sets), not part of the verdict. Writes <out>/verdict.md and
    returns True if all pass.
    """
    ok_all = all(c[1] for c in checks)
    lines = [f"# {test}: {'PASA' if ok_all else 'FALLA'} (MODEL)", ""]
    for crit, ok, detail in checks:
        lines.append(f"- [{'PASA' if ok else 'FALLA'}] {crit}" + (f" -- {detail}" if detail else ""))
    if notes:
        lines += ["", "Sensibilidad (no entra en el veredicto):"] + [f"- {n}" for n in notes]
    text = "\n".join(lines) + "\n"
    with open(os.path.join(out, "verdict.md"), "w") as f:
        f.write(text)
    print("\n" + text, flush=True)
    return ok_all


def advance(s, dt, c=K, sX=0.0, theta=1.0):
    """Nodal step of LeftAtrium2DPTAGX (ConservingKinetics::advance), 0D.

    s = (P, T, A, G, X, F); sX is a volumetric Xa source (mol m^-3 s^-1);
    theta multiplies k_a (1 + surface term).
    """
    import math
    P, T, A, G, X, F = (max(v, 0.0) for v in s)
    Xn = (X + dt * sX) / (1.0 + dt * c["k_dX"])
    eps = c.get("eps_x", 1.0)
    thetaX = eps + (1.0 - eps) * T / (c["K_A"] + T)
    Pn = P / (1.0 + dt * (c["k_a"] * theta * T * T / ((c["K_A"] + T) * c["Pref"])
                          + c["kappa"] * Xn * thetaX / c["Pref"]))
    dP = P - Pn
    a = dt * c["k_on"]
    b = 1.0 + a * (A + T + dP)
    cc = a * A * (T + dP)
    dI = 2.0 * cc / (b + math.sqrt(max(b * b - 4.0 * a * cc, 0.0)))
    Tn, An = T + dP - dI, A - dI
    Gn = G / (1.0 + dt * c["k_F"] * Tn / (c["K_G"] + G))
    return (Pn, Tn, An, Gn, Xn, F + (G - Gn))

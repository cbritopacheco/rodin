#!/usr/bin/env python3
"""Cycle measures of the coronary closure comparison (run_coronary_closures.sh).

For every directory <out>/.../<case>/ holding CoronaryArtery.csv (and
CoronaryArtery_fields.csv) the last cycle is reduced to:

  Q_signed, Q_abs        mean of Q_in and of |Q_in| (mL/s)
  Q_sys, Q_dia           |Q_in|-transport of the compressive and relaxed phases
                         (mL per cycle); the compressive phase is the part of
                         the cycle where p_im >= p_im,min + 0.5 (p_im,max -
                         p_im,min), or --systole t0,t1 in units of T
  phiA_flow              flow-weighted arteriolar viscosity factor,
                         int Phi_a |Q| dt / int |Q| dt
  E_diss                 dissipated energy per cycle in the 3D domain (mJ)
  WSS_mean               cycle mean of the wall-area mean |tau_w| (Pa)
  dP_drive               cycle mean of (p_in - p_out), averaged over outlets (Pa)

Every case is then compared with the rmu case of the same group (same blood,
drainage pressure, time step and mesh): differences in per cent, and the
systolic and diastolic shares of the signed and |Q| transport differences.
Writes <out>/coronary_summary.csv.
"""
import argparse
import csv
import os
import sys


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    return {k: [float(r[k]) if r[k] not in ("", None) else 0.0 for r in rows] for k in rows[0]}


def window(data, t0, t1):
    idx = [i for i, t in enumerate(data["t"]) if t0 - 1e-12 <= t < t1 - 1e-12]
    return {k: [v[i] for i in idx] for k, v in data.items()}


def reduce_case(d, f, T, cycles, systole):
    t0 = (cycles - 1) * T
    d = window(d, t0, t0 + T)
    n = len(d["t"])
    if n < 2:
        return None
    dt = (d["t"][-1] - d["t"][0]) / (n - 1)
    # the 3D flux is outward, so the inflow is its negative (mL/s)
    q = [-x * 1e6 for x in d["CoronaryInletFlux"]]
    qa = [abs(x) for x in q]
    if systole:
        a, b = systole
        sys_mask = [a <= (t - t0) / T < b for t in d["t"]]
    else:
        pim = d["IntramyoPressure"]
        lo, hi = min(pim), max(pim)
        sys_mask = [p >= lo + 0.5 * (hi - lo) for p in pim]
    phi = d.get("ArteriolarViscosityRatio", [0.0] * n)
    out = {
        "q": q,
        "sys": sys_mask,
        "dt": dt,
        "Q_signed": sum(q) / n,
        "Q_abs": sum(qa) / n,
        "Q_sys": dt * sum(x for x, s in zip(qa, sys_mask) if s),
        "Q_dia": dt * sum(x for x, s in zip(qa, sys_mask) if not s),
        "sys_fraction": sum(sys_mask) / n,
        "phiA_flow": sum(p * x for p, x in zip(phi, qa)) / (sum(qa) or 1e-300),
    }
    if f:
        f = window(f, t0, t0 + T)
        m = len(f["t"])
        if m >= 2:
            fdt = (f["t"][-1] - f["t"][0]) / (m - 1)
            out["E_diss"] = 1e3 * fdt * sum(f["dissipation"])
            out["WSS_mean"] = sum(f["meanAbsWSS"]) / m
            outs = [k for k in f if k.startswith("outlet") and k.endswith("Pressure")]
            drops = [sum(f["inletPressure"][i] - f[k][i] for k in outs) / len(outs) for i in range(m)]
            out["dP_drive"] = sum(drops) / m
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--period", type=float, default=0.85)
    ap.add_argument("--cycles", type=int, default=3)
    ap.add_argument("--systole", default="", help="t0,t1 in units of T (default: p_im criterion)")
    a = ap.parse_args()
    systole = tuple(map(float, a.systole.split(","))) if a.systole else None

    cases = {}
    for root, _, files in os.walk(a.out, followlinks=True):
        if "CoronaryArtery.csv" not in files:
            continue
        d = read(os.path.join(root, "CoronaryArtery.csv"))
        fpath = os.path.join(root, "CoronaryArtery_fields.csv")
        f = read(fpath) if os.path.isfile(fpath) else None
        if not d:
            continue
        r = reduce_case(d, f, a.period, a.cycles, systole)
        if r:
            cases[os.path.relpath(root, a.out)] = r

    if not cases:
        print("no runs found under", a.out)
        return 1

    scalars = ["Q_signed", "Q_abs", "Q_sys", "Q_dia", "sys_fraction", "phiA_flow",
               "E_diss", "WSS_mean", "dP_drive"]
    rows = []
    for name in sorted(cases):
        c = cases[name]
        group, case = os.path.split(name)
        row = {"group": group, "case": case}
        row.update({k: c.get(k, float("nan")) for k in scalars})
        ref = cases.get(os.path.join(group, "rmu"))
        if ref and case != "rmu":
            n = min(len(c["q"]), len(ref["q"]))
            dt = c["dt"]
            dsig = [c["q"][i] - ref["q"][i] for i in range(n)]
            dabs = [abs(c["q"][i]) - abs(ref["q"][i]) for i in range(n)]
            msk = ref["sys"][:n]
            tot_abs = dt * sum(abs(x) for x in ref["q"][:n]) or 1e-300
            row["dQ_signed_%"] = 100.0 * (c["Q_signed"] - ref["Q_signed"]) / (abs(ref["Q_signed"]) or 1e-300)
            row["dQ_abs_%"] = 100.0 * (c["Q_abs"] - ref["Q_abs"]) / (ref["Q_abs"] or 1e-300)
            row["dQ_signed_sys_%"] = 100.0 * dt * sum(x for x, s in zip(dsig, msk) if s) / tot_abs
            row["dQ_signed_dia_%"] = 100.0 * dt * sum(x for x, s in zip(dsig, msk) if not s) / tot_abs
            row["dQ_abs_sys_%"] = 100.0 * dt * sum(x for x, s in zip(dabs, msk) if s) / tot_abs
            row["dQ_abs_dia_%"] = 100.0 * dt * sum(x for x, s in zip(dabs, msk) if not s) / tot_abs
            for k in ("phiA_flow", "E_diss", "WSS_mean", "dP_drive"):
                if k in c and k in ref:
                    row["d%s_%%" % k] = 100.0 * (c[k] - ref[k]) / (abs(ref[k]) or 1e-300)
        rows.append(row)

    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(a.out, "coronary_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    print("| " + " | ".join(keys) + " |")
    print("|" + "---|" * len(keys))
    for r in rows:
        print("| " + " | ".join(("%.4g" % r[k]) if isinstance(r.get(k), float) else str(r.get(k, ""))
                                for k in keys) + " |")
    return 0


if __name__ == "__main__":
    sys.exit(main())

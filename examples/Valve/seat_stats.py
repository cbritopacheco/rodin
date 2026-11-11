#!/usr/bin/env python3
"""
Seat (seal) statistics of the KnifeGateValve runs: how much of the seat
pocket is swept by the flow, for the seat designs of make_meshes.py --seat
(run_valves_seat.sh). One PDF and one CSV.

Seat zone: the fluid of the lower half of the bore (y < 0) within |z| <=
Z_SEAT of the knife plane that lies outside the bore surface (the pocket, if
the design has one) or within LAYER of it (the near-wall layer, so that the
designs are compared on the same footing whatever the pocket). Time
averages over the saved fields after the transient (t_avg of the run).

Per case:
  u_mean, u_p10      volume-weighted mean and 10th percentile of the
                     time-averaged |u| in the seat zone [m/s]
  f_slow_0.1, 0.3    fraction of the zone volume with mean |u| below 0.1,
                     0.3 m/s
  tau_p10            10th percentile of the time-averaged shear stress
                     mu(gd) gd in the zone [Pa] (gd from the viscosity field
                     of the Bingham law, read from case.log)
  f_unyielded        fraction of the zone volume where the mean shear stress
                     is below the yield stress: material that sits there as
                     a plug
  pocket_*           the same restricted to the pocket proper (outside the
                     bore surface), when the design has one
and figures: metrics against the design; time-averaged |u|
and the shear rate on the plane x = 0 zoomed on the seat; mean |u| in the pocket against the
angle around the bore (bottom = -90 deg).

Usage:
  python3 seat_stats.py runs_seat [-o seat_report.pdf] [--zone-z 25] [--layer 1.0]
Needs numpy, matplotlib, h5py.
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import postprocess_valves as pv  # noqa: E402

BORE = {"FPC8": 101.6e-3, "FPC16": 191.0e-3}  # m, radius
DESIGN = {"": "CAD", "fillet": "A radios", "chamfer": "B chaflanes", "upramp": "C rampa",
          "rampup": "A+C arriba", "rampdn": "A+C abajo"}


def fluid_params(log):
    """tau_y, mu_p, m of the Bingham law from the [fluid] line of the log."""
    txt = open(log, errors="ignore").read()
    m = re.search(r"tau_y=([0-9.eE+-]+) Pa\s+mu_p=([0-9.eE+-]+) Pa s\s+m=([0-9.eE+-]+)", txt)
    if not m:
        return None
    return tuple(float(x) for x in m.groups())


def shear_rate(mu, tau_y, mu_p, m):
    """Inverse of mu(gd) = mu_p + tau_y (1 - exp(-m gd))/gd (decreasing)."""
    lo = np.zeros_like(mu)
    hi = np.full_like(mu, 1e6)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f = mu_p + tau_y * (-np.expm1(-m * mid)) / np.maximum(mid, 1e-30)
        up = f > mu  # viscosity at mid still above the target: gd is larger
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    return 0.5 * (lo + hi)


def case_fields(d, tag, t0):
    """Time-averaged |u| and viscosity per cell (deduplicated), cell centroids and volumes."""
    base, steps = pv.read_xdmf(os.path.join(d, tag + ".xdmf"))
    steps = [s for s in steps if s[0] > t0 + 1e-12] or steps[-1:]
    cache = {}
    meshes = pv.load_mesh(base, steps[0][1], cache)
    cen, vol, keyidx = [], [], []
    for V, T in meshes:
        P = V[T]
        c = P.mean(1)
        v = np.abs(np.einsum("ij,ij->i", np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), P[:, 3] - P[:, 0])) / 6.0
        cen.append(c); vol.append(v)
    cen = np.vstack(cen); vol = np.concatenate(vol)
    # ghost cells appear on several ranks: keep one copy
    key = np.round(cen * 1e6).astype(np.int64)
    _, uniq = np.unique(key, axis=0, return_index=True)
    umag = np.zeros(len(cen)); mu = np.zeros(len(cen)); n = 0
    for t, ranks in steps:
        off = 0
        for (V, T), r in zip(meshes, ranks):
            u = pv.h5(base, r["velocity"]).reshape(-1, 3)
            um = np.linalg.norm(u, axis=1)
            umag[off:off + len(T)] += um[T].mean(1)
            if "viscosity" in r:
                m_ = pv.h5(base, r["viscosity"]).reshape(-1)
                mu[off:off + len(T)] += m_[T].mean(1)
            off += len(T)
        n += 1
    umag /= n; mu /= n
    return cen[uniq], vol[uniq], umag[uniq], mu[uniq], len(steps), (steps[0][0], steps[-1][0]), meshes, steps, base


def wquantile(x, w, q):
    o = np.argsort(x)
    cw = np.cumsum(w[o]); cw /= cw[-1]
    return float(np.interp(q, cw, x[o]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs")
    ap.add_argument("-o", "--out", default="seat_report.pdf")
    ap.add_argument("--zone-z", type=float, default=25.0, help="half-width of the seat zone along the pipe, mm")
    ap.add_argument("--layer", type=float, default=1.0, help="near-wall layer of the bore included in the zone, mm")
    ap.add_argument("--skip", type=float, default=None, help="transient cut in D/U (default: t_avg of the run)")
    a = ap.parse_args()

    cases = sorted(d for d in glob.glob(os.path.join(a.runs, "*_open*")) if os.path.isdir(d))
    rows, planes, profiles = [], [], []
    for d in cases:
        tag = os.path.basename(d)
        valve, op = tag.split("_open"); op = int(op)
        base_valve, _, design = valve.partition("-")
        if not os.path.isfile(os.path.join(d, tag + ".xdmf")):
            print("[%s] no XDMF, skipped" % tag); continue
        R = BORE[base_valve]
        env = pv.read_env(d)
        tc = env.get("tc", 0.0)
        t0 = a.skip * tc if (a.skip is not None and tc > 0) else env.get("tavg", 0.0)
        par = fluid_params(os.path.join(d, "case.log"))
        cen, vol, umag, mu, nsteps, (ta, tb), meshes, steps, base = case_fields(d, tag, t0)
        r = np.hypot(cen[:, 0], cen[:, 1])
        zone = (cen[:, 1] < 0) & (np.abs(cen[:, 2]) <= a.zone_z * 1e-3) & (r > R - a.layer * 1e-3)
        pocket = zone & (r > R + 0.2e-3)
        if par:
            gd = shear_rate(mu, *par)
            tau = mu * gd
        else:
            tau = np.full_like(mu, np.nan)
        row = {"case": tag, "valve": base_valve, "design": design or "cad", "opening_pct": op,
               "snapshots": nsteps, "t_from_s": ta, "t_to_s": tb,
               "zone_volume_cm3": zone[zone].size and vol[zone].sum() * 1e6,
               "pocket_volume_cm3": vol[pocket].sum() * 1e6}
        for name, sel in (("", zone), ("pocket_", pocket)):
            if sel.sum() == 0:
                for k in ("u_mean", "u_p10", "f_slow_0.1", "f_slow_0.3", "tau_p10", "f_unyielded"):
                    row[name + k] = float("nan")
                continue
            w = vol[sel]; u = umag[sel]; tt = tau[sel]
            row[name + "u_mean"] = float((w * u).sum() / w.sum())
            row[name + "u_p10"] = wquantile(u, w, 0.1)
            row[name + "f_slow_0.1"] = float(w[u < 0.1].sum() / w.sum())
            row[name + "f_slow_0.3"] = float(w[u < 0.3].sum() / w.sum())
            if par:
                row[name + "tau_p10"] = wquantile(tt, w, 0.1)
                row[name + "f_unyielded"] = float(w[tt < par[0]].sum() / w.sum())
            else:
                row[name + "tau_p10"] = row[name + "f_unyielded"] = float("nan")
        rows.append(row)
        print("[%s] %d snapshots in [%.3f, %.3f] s: zone %.1f cm3, u_mean %.3f m/s, f_slow(0.1) %.2f, f_unyielded %.2f"
              % (tag, nsteps, ta, tb, row["zone_volume_cm3"], row["u_mean"], row["f_slow_0.1"], row["f_unyielded"]))
        # pocket mean |u| against the angle around the bore
        ang = np.degrees(np.arctan2(cen[:, 1], cen[:, 0]))
        bins = np.arange(-180, 1, 10)
        prof = []
        for b0, b1 in zip(bins[:-1], bins[1:]):
            s = pocket & (ang >= b0) & (ang < b1)
            prof.append((vol[s] * umag[s]).sum() / vol[s].sum() if s.sum() else np.nan)
        profiles.append((tag, 0.5 * (bins[:-1] + bins[1:]), np.array(prof)))
        # plane x = 0 of the time-averaged |u| and shear stress (nodal averages)
        V, T = meshes[0]
        allV = np.vstack([m[0] for m in meshes]); allT = np.vstack([m[1] + off for m, off in
                                                                  zip(meshes, np.cumsum([0] + [len(m[0]) for m in meshes[:-1]]))])
        un = np.zeros(len(allV)); mun = np.zeros(len(allV)); n = 0
        for t, ranks in steps:
            off = 0
            for (Vr, Tr), rk in zip(meshes, ranks):
                un[off:off + len(Vr)] += np.linalg.norm(pv.h5(base, rk["velocity"]).reshape(-1, 3), axis=1)
                if "viscosity" in rk:
                    mun[off:off + len(Vr)] += pv.h5(base, rk["viscosity"]).reshape(-1)
                off += len(Vr)
            n += 1
        un /= n; mun /= n
        gdn = shear_rate(mun, *par) if par else np.full_like(mun, np.nan)
        def cut(val):
            r_ = pv.slice_tets(allV, allT, val, 0)
            if r_ is None:
                return None
            X, Y, F, tri = r_
            keep = np.isfinite(X)[tri].all(1)
            return X, Y, F, tri[keep]
        planes.append((tag, cut(un), cut(gdn), R, par))

    if not rows:
        sys.exit("no case with output in " + a.runs)
    csvf = os.path.splitext(a.out)[0] + ".csv"
    with open(csvf, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader()
        for r_ in rows:
            wr.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r_.items()})

    pdf = PdfPages(a.out)
    # ---- table ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(11.7, 8.3)); ax.axis("off")
    cols = ["caso", "snapshots", "zona [cm³]", "ū [m/s]", "ū p10 [m/s]", "f(ū<0.1)", "f(ū<0.3)",
            "τ p10 [Pa]", "f(τ<τ_y)", "bolsa [cm³]", "bolsa f(ū<0.1)", "bolsa f(τ<τ_y)"]
    cells = [[r_["case"], r_["snapshots"], "%.1f" % r_["zone_volume_cm3"], "%.3f" % r_["u_mean"],
              "%.3f" % r_["u_p10"], "%.2f" % r_["f_slow_0.1"], "%.2f" % r_["f_slow_0.3"],
              "%.2g" % r_["tau_p10"], "%.2f" % r_["f_unyielded"], "%.1f" % r_["pocket_volume_cm3"],
              "%.2f" % r_["pocket_f_slow_0.1"], "%.2f" % r_["pocket_f_unyielded"]] for r_ in rows]
    tb_ = ax.table(cellText=cells, colLabels=cols, loc="upper center", cellLoc="center")
    tb_.auto_set_font_size(False); tb_.set_fontsize(7.5); tb_.scale(1.0, 1.5)
    ax.set_title("Asiento: fracción de la zona del sello barrida por el flujo (promedios tras el transitorio)")
    ax.text(0, 0.02, "Zona del sello: fluido con y < 0, |z| ≤ %.0f mm, fuera del diámetro o a menos de %.1f mm de él.\n"
            "ū: módulo de la velocidad promediado en el tiempo; f(ū<v): fracción del volumen de la zona con ū < v;\n"
            "τ = μ(γ̇) γ̇ (Bingham regularizado); f(τ<τ_y): fracción del volumen que no fluye (tapón). Bolsa: sólo fuera del diámetro.\n"
            "Datos: %s" % (a.zone_z, a.layer, os.path.basename(csvf)), transform=ax.transAxes, fontsize=8, va="bottom")
    pdf.savefig(fig); plt.close(fig)
    # ---- metrics against the design ------------------------------------------
    openings = sorted({r_["opening_pct"] for r_ in rows})
    designs = [d_ for d_ in ("cad", "fillet", "chamfer", "upramp", "rampup", "rampdn") if any(r_["design"] == d_ for r_ in rows)]
    fig, axs = plt.subplots(1, 3, figsize=(11.7, 4.2))
    for ax, key, title in zip(axs, ("u_mean", "f_slow_0.1", "f_unyielded"),
                              ("ū en la zona del sello [m/s]", "fracción con ū < 0.1 m/s", "fracción sin fluir (τ < τ_y)")):
        x = np.arange(len(designs)); wd = 0.8 / max(len(openings), 1)
        for i, o in enumerate(openings):
            vals = [next((r_[key] for r_ in rows if r_["design"] == d_ and r_["opening_pct"] == o), np.nan) for d_ in designs]
            ax.bar(x + (i - 0.5 * (len(openings) - 1)) * wd, vals, wd, label="%d %%" % o)
        ax.set_xticks(x); ax.set_xticklabels([DESIGN.get(d_ if d_ != "cad" else "", d_) for d_ in designs])
        ax.set_title(title, fontsize=10); ax.grid(axis="y", alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("Diseños del asiento: FPC8, pulpa Bingham", fontsize=11); fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)
    # ---- pocket profile against the angle --------------------------------------
    fig, axs = plt.subplots(1, len(openings), figsize=(11.7, 4.2), squeeze=False)
    for ax, o in zip(axs[0], openings):
        for tag, angc, prof in profiles:
            if int(tag.split("_open")[1]) != o:
                continue
            ax.plot(angc, prof, marker="o", ms=3, label=tag.split("_open")[0])
        ax.set_xlabel("ángulo alrededor del diámetro [°] (−90: fondo, 0/−180: costados)")
        ax.set_ylabel("ū en la bolsa [m/s]"); ax.set_title("%d %% de apertura" % o, fontsize=10)
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("Velocidad media en la bolsa del asiento a lo largo del perímetro (sin bolsa: nada)", fontsize=11)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    # ---- planes ---------------------------------------------------------------
    for tag, pu, pt, R, par in planes:
        fig, axs = plt.subplots(1, 2, figsize=(11.7, 4.6))
        zoom = ((-0.04, 0.04), (-R - 0.006, -R + 0.05))
        for ax, pf, vmax, title, cmap in ((axs[0], pu, None, "ū [m/s]", "jet"),
                                          (axs[1], pt, None, "log10 γ̇ [1/s] (de la viscosidad promedio)", "viridis")):
            if pf is None:
                continue
            X, Y, F, tris = pf
            if title.startswith("log"):
                F = np.log10(np.maximum(F, 1e-2))
            sel = (X > zoom[0][0]) & (X < zoom[0][1]) & (Y > zoom[1][0]) & (Y < zoom[1][1])
            vmin = np.nanmin(F[sel]) if sel.any() else None
            vmax = np.nanmax(F[sel]) if sel.any() else None
            im = pv.draw_plane(ax, (X, Y, F, tris), vmin, vmax, cmap)
            ax.set_xlim(*zoom[0]); ax.set_ylim(*zoom[1]); ax.set_aspect("equal")
            ax.axhline(-R, color="w", lw=0.6, ls="--")
            if par and title.startswith("log"):
                ax.set_title("%s; tapón (τ < τ_y) para γ̇ < %.2g 1/s (log10 = %.1f)"
                             % (title, 1.0 / par[2], math.log10(1.0 / par[2])), fontsize=8)
            else:
                ax.set_title(title, fontsize=9)
            ax.set_xlabel("z [m]"); ax.set_ylabel("y [m]")
            fig.colorbar(im, ax=ax, shrink=0.8)
        fig.suptitle("%s — plano x = 0, zona del sello (fondo del diámetro, cuchillo en z = 0)" % tag, fontsize=11)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    pdf.close()
    print("PDF:", a.out); print("CSV:", csvf)


if __name__ == "__main__":
    main()

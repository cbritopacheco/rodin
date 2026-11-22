#!/usr/bin/env python3
"""Tables and figures of the coronary campaign from whatever has run so far.

Every comparison uses the last cycle that ALL runs being compared have
completed, so partial runs never mix a complete and an incomplete cycle; the
cycle used is written next to every number.

  python3 report.py ~/Results/cmbpe/coronary            # -> OUT/report/

Writes to OUT/report/:
  summary_latest.csv   one row per run, compare_coronary_closures format plus
                       the column `cycle` (plot_phase_split.py reads it)
  tables.md            closures, verification (dt, mesh) and periodicity
  tables.tex           the same as booktabs tables for the manuscript
  fig_waveforms.pdf    |Q_in| over the cycle and the difference of every
                       constant closure from R_mu, one column per group
  fig_phase_split.pdf  systolic/diastolic split (if plot_phase_split.py is
                       next to this script)
"""
import argparse
import csv
import os
import subprocess
import sys

GROUPS = [("pra700/healthy", "Healthy, 700 Pa"), ("pra700/diabetic", "Diabetic, 700 Pa"),
          ("pra700/hyperviscous", "Hyperviscous, 700 Pa"),
          ("pra2400/hyperviscous", "Hyperviscous, 2400 Pa")]
CLOSURES = [("rest", r"$\mathcal{R}_0$", "#c4581f", "-"),
            ("rest_matched", r"$\mathcal{R}_0^*$", "#c4581f", "--"),
            ("plateau", r"$\mathcal{R}_\infty$", "#7a4fa0", "-")]
RMU = ("rmu", r"$\mathcal{R}_\mu$", "#2f6fae", "-")
DT_REF, MESH_REF = "pra2400/hyperviscous", "pra700/healthy"


# ---------------------------------------------------------------- reading --
def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    return {k: [float(r[k]) if r[k] not in ("", None) else 0.0 for r in rows] for k in rows[0]}


def window(d, t0, t1):
    idx = [i for i, t in enumerate(d["t"]) if t0 - 1e-12 <= t < t1 - 1e-12]
    return {k: [v[i] for i in idx] for k, v in d.items()}


def complete_cycles(d, T):
    t = d["t"]
    dt = t[1] - t[0] if len(t) > 1 else 0.0
    return int((t[-1] + dt + 1e-9) // T) if t else 0


def reduce(run, T, c):
    """Cycle-c measures of one run (same definitions as compare_coronary_closures)."""
    t0 = (c - 1) * T
    d = window(run["d"], t0, t0 + T)
    n = len(d["t"])
    if n < 2:
        return None
    dt = (d["t"][-1] - d["t"][0]) / (n - 1)
    q = [-x * 1e6 for x in d["CoronaryInletFlux"]]
    qa = [abs(x) for x in q]
    pim = d.get("IntramyoPressure", [0.0] * n)
    lo, hi = min(pim), max(pim)
    sysm = [p >= lo + 0.5 * (hi - lo) for p in pim] if hi > lo else [False] * n
    phi = d.get("ArteriolarViscosityRatio", [0.0] * n)
    r = dict(t=[(x - t0) / T for x in d["t"]], q=q, sys=sysm, dt=dt,
             Q_signed=sum(q) / n, Q_abs=sum(qa) / n,
             Q_sys=dt * sum(x for x, s in zip(qa, sysm) if s),
             Q_dia=dt * sum(x for x, s in zip(qa, sysm) if not s),
             sys_fraction=sum(sysm) / n,
             phiA_flow=sum(p * x for p, x in zip(phi, qa)) / (sum(qa) or 1e-300))
    f = run["f"]
    if f:
        f = window(f, t0, t0 + T)
        m = len(f["t"])
        if m >= 2:
            fdt = (f["t"][-1] - f["t"][0]) / (m - 1)
            r["E_diss"] = 1e3 * fdt * sum(f["dissipation"])
            r["WSS_mean"] = sum(f["meanAbsWSS"]) / m
            outs = [k for k in f if k.startswith("outlet") and k.endswith("Pressure")]
            if outs:
                r["dP_drive"] = sum(sum(f["inletPressure"][i] - f[k][i] for k in outs) / len(outs)
                                    for i in range(m)) / m
    return r


def diff(c, ref):
    n = min(len(c["q"]), len(ref["q"]))
    dt, msk = c["dt"], ref["sys"][:n]
    tot = dt * sum(abs(x) for x in ref["q"][:n]) or 1e-300
    dabs = [abs(c["q"][i]) - abs(ref["q"][i]) for i in range(n)]
    dsig = [c["q"][i] - ref["q"][i] for i in range(n)]
    out = {"dQ_signed_%": 100 * (c["Q_signed"] - ref["Q_signed"]) / (abs(ref["Q_signed"]) or 1e-300),
           "dQ_abs_%": 100 * (c["Q_abs"] - ref["Q_abs"]) / ref["Q_abs"],
           "dQ_signed_sys_%": 100 * dt * sum(x for x, s in zip(dsig, msk) if s) / tot,
           "dQ_signed_dia_%": 100 * dt * sum(x for x, s in zip(dsig, msk) if not s) / tot,
           "dQ_abs_sys_%": 100 * dt * sum(x for x, s in zip(dabs, msk) if s) / tot,
           "dQ_abs_dia_%": 100 * dt * sum(x for x, s in zip(dabs, msk) if not s) / tot}
    for k in ("phiA_flow", "E_diss", "WSS_mean", "dP_drive"):
        if k in c and k in ref:
            out["d%s_%%" % k] = 100 * (c[k] - ref[k]) / (abs(ref[k]) or 1e-300)
    return out


def load(out):
    runs, seen = {}, {}
    for root, dirs, files in os.walk(out, followlinks=True):
        dirs[:] = sorted(x for x in dirs if not x.startswith(".") and x != "report")
        if "CoronaryArtery.csv" not in files:
            continue
        name = os.path.relpath(root, out)
        real = os.path.realpath(root)
        if real in seen:
            if seen[real].startswith("pra") or not name.startswith("pra"):
                continue                   # keep the first pra*/ name of a linked run
            runs.pop(seen[real], None)
        d = read(os.path.join(root, "CoronaryArtery.csv"))
        if not d or "CoronaryInletFlux" not in d:
            continue
        fp = os.path.join(root, "CoronaryArtery_fields.csv")
        runs[name] = dict(d=d, f=read(fp) if os.path.isfile(fp) else None)
        seen[real] = name
    return runs


# ------------------------------------------------------------------ main --
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--period", type=float, default=0.85)
    a = ap.parse_args()
    out = os.path.abspath(os.path.expanduser(a.out))
    T = a.period
    rep = os.path.join(out, "report")
    os.makedirs(rep, exist_ok=True)
    runs = load(out)
    if not runs:
        sys.exit("no runs under " + out)
    for r in runs.values():
        r["n"] = complete_cycles(r["d"], T)

    def run(name):
        r = runs.get(name)
        return r if r and r["n"] >= 1 else None

    md, tex = [], []

    # ---- closures per group ------------------------------------------------
    rows, wave = [], {}
    md += ["## Closures vs R_mu (last cycle completed by both the closure and R_mu)", "",
           "| group | cycle | |Q| R_mu (mL/s) | closure | d|Q| % | sys / dia share | dE % | dWSS % | dP % |",
           "|---|---|---|---|---|---|---|---|---|"]
    tex += [r"\begin{tabular}{@{}llccccccc@{}}", r"\toprule",
            r"State & Cycle & $\langle|Q_{\rm in}|\rangle_\mu$ (mL\,s$^{-1}$) & Closure & "
            r"$\Delta\langle|Q|\rangle$ (\%) & Systole / diastole (\%) & $\Delta E$ (\%) & "
            r"$\Delta\overline{\tau}_w$ (\%) & $\Delta(\Delta p)$ (\%)\\", r"\midrule"]
    for g, label in GROUPS:
        ref = run(g + "/rmu")
        present = [(c, run(g + "/" + c)) for c, *_ in CLOSURES if run(g + "/" + c)]
        if not ref:
            continue
        top = max([min(ref["n"], r["n"]) for _, r in present] or [ref["n"]])
        rtop = reduce(ref, T, top)
        wave[g] = (label, top, rtop, {})
        row = {"group": g, "case": "rmu", "cycle": top}
        row.update({k: rtop.get(k, "") for k in ("Q_signed", "Q_abs", "Q_sys", "Q_dia",
                                                  "sys_fraction", "phiA_flow", "E_diss",
                                                  "WSS_mean", "dP_drive")})
        rows.append(row)
        first = True
        for c, r in present:
            cyc = min(ref["n"], r["n"])      # last cycle completed by both
            rr = reduce(ref, T, cyc)
            cr = reduce(r, T, cyc)
            dd = diff(cr, rr)
            wave[g][3][c] = (cyc, cr, rr)
            row = {"group": g, "case": c, "cycle": cyc}
            row.update({k: cr.get(k, "") for k in ("Q_signed", "Q_abs", "Q_sys", "Q_dia",
                                                   "sys_fraction", "phiA_flow", "E_diss",
                                                   "WSS_mean", "dP_drive")})
            row.update(dd)
            rows.append(row)
            tot = dd["dQ_abs_%"] or 1e-300
            share = "%.0f / %.0f" % (100 * dd["dQ_abs_sys_%"] / tot, 100 * dd["dQ_abs_dia_%"] / tot)
            cl = dict((k, l) for k, l, *_ in CLOSURES)[c]
            md.append("| %s | %d | %.4f | %s | %+.2f | %s | %s | %s | %s |" % (
                g, cyc, rr["Q_abs"], c, dd["dQ_abs_%"], share,
                *("%+.2f" % dd[k] if k in dd else "-" for k in ("dE_diss_%", "dWSS_mean_%",
                                                                 "ddP_drive_%"))))
            tex.append("%s & %d & %.3f & %s & $%+.2f$ & %s & %s & %s & %s\\\\" % (
                label if first else "", cyc, rr["Q_abs"], cl, dd["dQ_abs_%"], share,
                *("$%+.1f$" % dd[k] if k in dd else "--" for k in ("dE_diss_%", "dWSS_mean_%",
                                                                   "ddP_drive_%"))))
            first = False
        if not present:
            md.append("| %s | %d | %.4f | (no constant closure yet) | | | | | |" % (g, top, rtop["Q_abs"]))
    tex += [r"\bottomrule", r"\end{tabular}", ""]

    keys = ["group", "case", "cycle", "Q_signed", "Q_abs", "Q_sys", "Q_dia", "sys_fraction",
            "phiA_flow", "E_diss", "WSS_mean", "dP_drive", "dQ_signed_%", "dQ_abs_%",
            "dQ_signed_sys_%", "dQ_signed_dia_%", "dQ_abs_sys_%", "dQ_abs_dia_%",
            "dphiA_flow_%", "dE_diss_%", "dWSS_mean_%", "ddP_drive_%"]
    with open(os.path.join(rep, "summary_latest.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})

    # ---- verification: time step and mesh ------------------------------------
    md += ["", "## Verification (common completed cycle of the runs compared)", "",
           "| study | level | closure | cycle | |Q| (mL/s) | change vs reference % | R_0 - R_mu (pp) |",
           "|---|---|---|---|---|---|---|"]
    tex += [r"\begin{tabular}{@{}llcccc@{}}", r"\toprule",
            r"Study & Level & Cycle & $\langle|Q_{\rm in}|\rangle_\mu$ (mL\,s$^{-1}$) & "
            r"Change (\%) & $\Delta\langle|Q|\rangle_{0-\mu}$ (\%)\\", r"\midrule"]
    studies = [("time step", DT_REF, "dt = 1e-3",
                sorted({(x.split("/")[1], "dt/" + x.split("/")[1]) for x in runs
                       if x.startswith("dt/") and x.count("/") == 2 and "dt1e-3" not in x})),
               ("mesh", MESH_REF, "medium", [("coarse", "mesh/coarse"), ("fine", "mesh/fine")])]
    for study, refg, reflab, levels in studies:
        levels = [(lab, g) for lab, g in levels if run(g + "/rmu")]
        if not run(refg + "/rmu") or not levels:
            continue
        for lab, g in levels:
            names = [refg + "/rmu", g + "/rmu"] + [x + "/rest" for x in (refg, g) if run(x + "/rest")]
            cyc = min(run(x)["n"] for x in names)
            a0, a1 = reduce(run(refg + "/rmu"), T, cyc), reduce(run(g + "/rmu"), T, cyc)
            d0 = d1 = None
            if run(refg + "/rest") and run(g + "/rest"):
                d0 = diff(reduce(run(refg + "/rest"), T, cyc), a0)["dQ_abs_%"]
                d1 = diff(reduce(run(g + "/rest"), T, cyc), a1)["dQ_abs_%"]
            ch = 100 * (a1["Q_abs"] - a0["Q_abs"]) / a0["Q_abs"]
            md.append("| %s | %s (ref) | rmu | %d | %.4f | | %s |" % (
                study, reflab, cyc, a0["Q_abs"], "%+.2f" % d0 if d0 is not None else "-"))
            md.append("| %s | %s | rmu | %d | %.4f | %+.2f | %s |" % (
                study, lab, cyc, a1["Q_abs"], ch, "%+.2f" % d1 if d1 is not None else "-"))
            tex.append("%s & %s / %s & %d & %.4f / %.4f & $%+.2f$ & %s\\\\" % (
                study.capitalize(), reflab, lab, cyc, a0["Q_abs"], a1["Q_abs"], ch,
                "$%+.2f$ / $%+.2f$" % (d0, d1) if d0 is not None else "--"))
    tex += [r"\bottomrule", r"\end{tabular}", ""]

    # ---- periodicity ---------------------------------------------------------
    md += ["", "## Periodicity of R_mu (mean |Q_in| per complete cycle, mL/s)", "",
           "| group | cycles | change of the last cycle % |", "|---|---|---|"]
    for g, _ in GROUPS + [("mesh/coarse", ""), ("mesh/fine", "")]:
        r = run(g + "/rmu")
        if not r:
            continue
        qs = [reduce(r, T, c)["Q_abs"] for c in range(1, r["n"] + 1)]
        md.append("| %s | %s | %s |" % (g, " ".join("%.4f" % x for x in qs),
                                        "%+.2f" % (100 * (qs[-1] - qs[-2]) / qs[-2]) if len(qs) > 1 else "-"))

    open(os.path.join(rep, "tables.md"), "w").write("\n".join(md) + "\n")
    open(os.path.join(rep, "tables.tex"), "w").write("\n".join(tex) + "\n")
    print("\n".join(md))

    # ---- waveform figure -----------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5,
                         "axes.labelsize": 8.5, "legend.fontsize": 8, "axes.linewidth": 0.6,
                         "xtick.direction": "out", "ytick.direction": "out"})
    gs = [g for g, _ in GROUPS if g in wave]
    if gs:
        fig, axes = plt.subplots(2, len(gs), figsize=(1.9 * len(gs) + 0.6, 4.2), sharex=True,
                                 sharey="row", squeeze=False)
        for j, g in enumerate(gs):
            label, cyc, rr, others = wave[g]
            top, bot = axes[0][j], axes[1][j]
            s = [x for x, m in zip(rr["t"], rr["sys"]) if m]
            for ax in (top, bot):
                if s:
                    ax.axvspan(min(s), max(s), color="#eeeeee", lw=0, zorder=0)
                for side in ("top", "right"):
                    ax.spines[side].set_visible(False)
            top.plot(rr["t"], [abs(x) for x in rr["q"]], color=RMU[2], lw=1.3, zorder=3)
            cycs = []
            for c, lab, col, ls in CLOSURES:
                if c not in others:
                    continue
                cc, cr, rc = others[c]
                if cc != cyc:
                    cycs.append("%s: %d" % (lab, cc))
                n = min(len(cr["q"]), len(rc["q"]))
                top.plot(cr["t"][:n], [abs(x) for x in cr["q"][:n]], color=col, ls=ls, lw=1.0, zorder=2)
                bot.plot(cr["t"][:n], [100 * (abs(cr["q"][i]) - abs(rc["q"][i])) / cr["Q_abs"]
                                       for i in range(n)], color=col, ls=ls, lw=1.1)
            bot.axhline(0, color="#404040", lw=0.6)
            top.set_title("%s\ncycle %d%s" % (label, cyc, " (" + ", ".join(cycs) + ")" if cycs else ""),
                          fontsize=7.5)
            bot.set_xlabel(r"$(t-t_0)/T$")
        axes[0][0].set_ylabel(r"$|Q_{\rm in}|$ (mL s$^{-1}$)")
        axes[1][0].set_ylabel(r"$100\,\delta_Q$ (%)")
        handles = [Line2D([], [], color=RMU[2], lw=1.4, label=RMU[1])]
        handles += [Line2D([], [], color=col, ls=ls, lw=1.2, label=lab) for _, lab, col, ls in CLOSURES]
        fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0))
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        fig.savefig(os.path.join(rep, "fig_waveforms.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(rep, "fig_waveforms.png"), dpi=200, bbox_inches="tight")

    ps = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_phase_split.py")
    if os.path.isfile(ps):
        subprocess.call([sys.executable, ps, os.path.join(rep, "summary_latest.csv"),
                         "--out", os.path.join(rep, "fig_phase_split.pdf")])
    print("\nwritten to", rep)


if __name__ == "__main__":
    main()

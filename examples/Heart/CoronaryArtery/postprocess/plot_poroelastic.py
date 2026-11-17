#!/usr/bin/env python3
"""Poroelastic LV and coronary 0D/3D coupling: cycle curves and metrics.

Reads <results>/<geometry>/<scenario>/CoronaryArtery_Explicit_PoroElastic.csv
for every case found, keeps the last complete cardiac cycle and writes

  fig_lv_poroelastic.pdf     LV mechanics and myocardial perfusion (2 x 3)
  fig_coronary_coupling.pdf  coronary inflow, outlet pressure, microcirculation
  fig_phase_loops.pdf        Phi - p~ and Q_ar - p~ loops (poroelastic hysteresis)
  fig_scenario_changes.pdf   relative change of the cycle metrics vs. healthy
  metrics_poroelastic.csv / .tex

Usage: plot_poroelastic.py RESULTS_DIR FIGURE_DIR [--period 0.85] [--metrics3d CSV]
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmame_style as st  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

CSV_NAME = "CoronaryArtery_Explicit_PoroElastic.csv"


def read_csv(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size < 2:
        raise ValueError(f"{path}: fewer than two rows")
    return {name: np.asarray(data[name], float) for name in data.dtype.names}


def discover(results):
    cases = []
    for path in sorted(glob.glob(os.path.join(results, "*", "*", CSV_NAME))):
        scen_dir = os.path.dirname(path)
        geom = os.path.basename(os.path.dirname(scen_dir))
        scen = os.path.basename(scen_dir)
        cases.append((geom, scen, path))

    def key(c):
        g = st.GEOMETRY_ORDER.index(c[0]) if c[0] in st.GEOMETRY_ORDER else 99
        s = st.SCENARIO_ORDER.index(c[1]) if c[1] in st.SCENARIO_ORDER else 99
        return (g, s, c[0], c[1])

    return sorted(cases, key=key)


def last_cycle(d, period):
    """Rows of the last complete cycle (or everything after the first cycle)."""
    t = d["t"]
    n = int(np.floor(t[-1] / period + 1e-9))
    if n >= 1:
        t0 = (n - 1) * period
        mask = (t >= t0 - 1e-12) & (t <= t0 + period + 1e-12)
    else:
        t0 = t[0]
        mask = np.ones_like(t, bool)
        print(f"  warning: less than one cycle simulated (t_end = {t[-1]:.3f} s)")
    out = {k: v[mask] for k, v in d.items()}
    out["tau"] = (out["t"] - t0) / period
    out["t0"] = t0
    return out


def systole_mask(c):
    """Ejection + isovolumic phases: p_v above 10% of its cycle range."""
    pv = c["lv_pv"]
    return pv > pv.min() + 0.1 * (pv.max() - pv.min())


def shade_systole(ax, c):
    m = systole_mask(c)
    if not m.any():
        return
    tau = c["tau"]
    edges = np.flatnonzero(np.diff(np.r_[0, m.astype(int), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        ax.axvspan(tau[a], tau[b - 1], color="0.92", lw=0, zorder=0)


def mean_over(tau, y, mask=None):
    if mask is None:
        mask = np.ones_like(tau, bool)
    if mask.sum() < 2:
        return float(np.mean(y[mask])) if mask.any() else np.nan
    integrate = getattr(np, "trapezoid", None) or np.trapz
    return float(integrate(y[mask], tau[mask]) / (tau[mask][-1] - tau[mask][0]))


# ---------------------------------------------------------------------------
# Derived quantities (plot units)
# ---------------------------------------------------------------------------
def q_in(c):
    return -c["q_in"] * st.M3S_TO_MLMIN  # inflow positive, mL/min


def p_out_mean(c):
    keys = [k for k in c if k.startswith("p_out_")]
    return np.mean([c[k] for k in keys], axis=0) / st.PA_PER_MMHG if keys else None


QUANTITIES = {
    # key: (function, axis label)
    "pv": (lambda c: c["lv_pv"] / st.PA_PER_MMHG, r"$p_v$ [mmHg]"),
    "par": (lambda c: c["lv_par"] / st.PA_PER_MMHG, r"$p_{ar}$ [mmHg]"),
    "h": (lambda c: 1e3 * c["lv_h"], r"Wall thickness $h$ [mm]"),
    "phi": (lambda c: 100 * c["lv_phi"], r"Perfusion fraction $\Phi$ [%]"),
    "pf": (lambda c: c["lv_pf"] / st.PA_PER_MMHG,
           r"$\tilde p = p_P(\Phi) - \bar\lambda$ [mmHg]"),
    "pP": (lambda c: c["lv_pP"] / st.PA_PER_MMHG, r"Storage pressure $p_P$ [mmHg]"),
    "qar": (lambda c: c["lv_q_ar"] * st.M3S_TO_MLMIN, r"$Q_{ar}$ [mL/min]"),
    "qven": (lambda c: c["lv_q_ven"] * st.M3S_TO_MLMIN, r"$Q_{ven}$ [mL/min]"),
    "gamma": (lambda c: c["lv_shear_rate"], r"Pore shear rate $\dot\gamma$ [s$^{-1}$]"),
    "nu": (lambda c: 1e3 * c["lv_nu_eff"], r"Pore viscosity $\nu_{CY}$ [mPa$\,$s]"),
    "qin": (q_in, r"LCA inflow $Q_{in}$ [mL/min]"),
    "pout": (p_out_mean, r"Mean outlet pressure [mmHg]"),
    "ptm": (lambda c: c["ptm_mean"] / st.PA_PER_MMHG,
            r"Transmural pressure $p_{tm}$ [mmHg]"),
    "phia": (lambda c: c["phi_a_mean"], r"Arteriolar viscosity ratio $\varphi_a$ [-]"),
    "phiv": (lambda c: c["phi_v_mean"], r"Venular viscosity ratio $\varphi_v$ [-]"),
    "vol": (lambda c: 1e6 * c["stored_volume"], r"Distal stored volume [mL]"),
}


def curve_figure(cases, keys, path, ncols=3, height_per_row=48 * st.MM):
    nrows = int(np.ceil(len(keys) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(st.DOUBLE, nrows * height_per_row),
                             sharex=True, squeeze=False)
    letters = "abcdefghijklmnop"
    for i, key in enumerate(keys):
        ax = axes.flat[i]
        fn, label = QUANTITIES[key]
        shade_systole(ax, cases[0][2])
        for geom, scen, c in cases:
            y = fn(c)
            if y is None:
                continue
            ax.plot(c["tau"], y, color=st.color(scen), ls=st.dash(geom))
        ax.set_ylabel(label)
        ax.set_xlim(0, 1)
        st.panel_label(ax, letters[i], x=-0.20)
        if i // ncols == nrows - 1:
            ax.set_xlabel(r"$t/T$ (last cycle)")
    for j in range(len(keys), nrows * ncols):
        axes.flat[j].axis("off")
    legend(fig, cases)
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=1.2, h_pad=0.8)
    st.save(fig, path)
    plt.close(fig)


def legend(fig, cases, systole=True):
    scens = list(dict.fromkeys(s for _, s, _ in cases))
    geoms = list(dict.fromkeys(g for g, _, _ in cases))
    handles = [Line2D([], [], color=st.color(s), lw=1.2, label=st.scenario_label(s))
               for s in scens]
    if len(geoms) > 1:
        handles += [Line2D([], [], color="0.35", ls=st.dash(g), lw=1.0,
                           label=st.geometry_label(g)) for g in geoms]
    if systole:
        handles.append(plt.Rectangle((0, 0), 1, 1, color="0.92", lw=0, label="Systole"))
    fig.legend(handles=handles, loc="upper center", ncol=len(handles),
               bbox_to_anchor=(0.5, 1.0), columnspacing=1.4)


def loop_figure(cases, path):
    fig, axes = plt.subplots(1, 2, figsize=(st.ONEHALF, 58 * st.MM))
    for geom, scen, c in cases:
        pf = QUANTITIES["pf"][0](c)
        axes[0].plot(pf, QUANTITIES["phi"][0](c), color=st.color(scen), ls=st.dash(geom))
        axes[1].plot(pf, QUANTITIES["qar"][0](c), color=st.color(scen), ls=st.dash(geom))
        # arrow showing the sense of the loop
        k = len(pf) // 4
        for ax, y in ((axes[0], QUANTITIES["phi"][0](c)), (axes[1], QUANTITIES["qar"][0](c))):
            ax.annotate("", xy=(pf[k + 1], y[k + 1]), xytext=(pf[k], y[k]),
                        arrowprops=dict(arrowstyle="-|>", color=st.color(scen), lw=0.8,
                                        mutation_scale=7))
    axes[0].set_xlabel(QUANTITIES["pf"][1])
    axes[0].set_ylabel(QUANTITIES["phi"][1])
    axes[1].set_xlabel(QUANTITIES["pf"][1])
    axes[1].set_ylabel(QUANTITIES["qar"][1])
    st.panel_label(axes[0], "a", x=-0.18)
    st.panel_label(axes[1], "b", x=-0.18)
    legend(fig, cases, systole=False)
    fig.tight_layout(rect=(0, 0, 1, 0.88), w_pad=2.0)
    st.save(fig, path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Cycle metrics
# ---------------------------------------------------------------------------
METRICS = [
    # name, latex label, unit
    ("Qin_mean", r"$\overline{Q}_{in}$", "mL/min"),
    ("Qin_dia_sys", r"$\overline{Q}_{in}^{dia}/\overline{Q}_{in}^{sys}$", "-"),
    ("Qar_mean", r"$\overline{Q}_{ar}$", "mL/min"),
    ("Qven_mean", r"$\overline{Q}_{ven}$", "mL/min"),
    ("pv_max", r"$\max p_v$", "mmHg"),
    ("pf_max", r"$\max\tilde p$", "mmHg"),
    ("pf_mean", r"$\overline{\tilde p}$", "mmHg"),
    ("phi_mean", r"$\overline{\Phi}$", "%"),
    ("phi_amp", r"$\Delta\Phi$", "%"),
    ("nu_mean", r"$\overline{\nu}_{CY}$", "mPa s"),
    ("pout_mean", r"$\overline{p}_{out}$", "mmHg"),
    ("ptm_mean", r"$\overline{p}_{tm}$", "mmHg"),
]


def metrics_of(c):
    tau = c["tau"]
    sys_m = systole_mask(c)
    dia_m = ~sys_m
    qin = q_in(c)
    phi = QUANTITIES["phi"][0](c)
    pf = QUANTITIES["pf"][0](c)
    m = {
        "Qin_mean": mean_over(tau, qin),
        "Qin_dia_sys": (mean_over(tau, qin, dia_m) / mean_over(tau, qin, sys_m)
                        if sys_m.sum() > 1 and dia_m.sum() > 1 else np.nan),
        "Qar_mean": mean_over(tau, QUANTITIES["qar"][0](c)),
        "Qven_mean": mean_over(tau, QUANTITIES["qven"][0](c)),
        "pv_max": float(QUANTITIES["pv"][0](c).max()),
        "pf_max": float(pf.max()),
        "pf_mean": mean_over(tau, pf),
        "phi_mean": mean_over(tau, phi),
        "phi_amp": float(phi.max() - phi.min()),
        "nu_mean": mean_over(tau, QUANTITIES["nu"][0](c)),
        "pout_mean": mean_over(tau, p_out_mean(c)) if p_out_mean(c) is not None else np.nan,
        "ptm_mean": mean_over(tau, QUANTITIES["ptm"][0](c)),
    }
    return m


# ---------------------------------------------------------------------------
# Diagnostics: periodicity, interface mass balance, solver effort
# ---------------------------------------------------------------------------
def cycle_rows(d, period, k_from_end):
    """Rows of the k-th last complete cycle (k=1: last)."""
    t = d["t"]
    n = int(np.floor(t[-1] / period + 1e-9))
    if n - k_from_end < 0:
        return None
    t0 = (n - k_from_end) * period
    m = (t >= t0 - 1e-12) & (t <= t0 + period + 1e-12)
    return m if m.sum() > 2 else None


def diagnostics(d, period):
    """Checks one should look at before trusting a run."""
    integrate = getattr(np, "trapezoid", None) or np.trapz
    t = d["t"]
    out = {}
    m1, m2 = cycle_rows(d, period, 1), cycle_rows(d, period, 2)
    m = m1 if m1 is not None else np.ones_like(t, bool)

    def avg(y, mask):
        return float(integrate(y[mask], t[mask]) / (t[mask][-1] - t[mask][0]))

    qar = d["lv_q_ar"]
    qven = d["lv_q_ven"]
    qin = -d["q_in"]
    out["n_cycles"] = float(np.floor(t[-1] / period + 1e-9))
    # 0D mass balance over the last cycle: stored volume must return -> Q_ar = Q_ven
    out["Qar_Qven_imbalance"] = 100 * abs(avg(qar, m) - avg(qven, m)) / max(abs(avg(qar, m)), 1e-30)
    # cycle-to-cycle change of the mean inflow and of the perfusion fraction
    if m1 is not None and m2 is not None:
        out["dQin_cycle"] = 100 * abs(avg(qin, m1) - avg(qin, m2)) / max(abs(avg(qin, m1)), 1e-30)
        out["dphi_cycle"] = 100 * abs(avg(d["lv_phi"], m1) - avg(d["lv_phi"], m2)) / max(
            abs(avg(d["lv_phi"], m1)), 1e-30)
    else:
        out["dQin_cycle"] = out["dphi_cycle"] = np.nan
    # 3D-0D consistency: LCA inflow vs 0D arterial inflow
    out["Qin_vs_Qar"] = 100 * (avg(qin, m) - avg(qar, m)) / max(abs(avg(qar, m)), 1e-30)
    # interface leakage int (u - w).n over the wet wall (exactly 0 for a perfect coupling)
    if "q_interface_slip" in d:
        qs = d["q_interface_slip"]
        out["leak_net"] = 100 * abs(avg(qs, m)) / max(abs(avg(qin, m)), 1e-30)
        out["leak_rms"] = 100 * np.sqrt(avg(qs ** 2, m)) / max(np.sqrt(avg(qin ** 2, m)), 1e-30)
    for k, name in (("coupling_its", "its_mean"), ("gmres_its", "gmres_mean")):
        if k in d:
            out[name] = float(np.mean(d[k][m]))
            out[name.replace("mean", "max")] = float(np.max(d[k]))
    if "coupling_rel" in d:
        out["res_max"] = float(np.max(d["coupling_rel"][m]))
    return out


DIAGNOSTICS = [
    # name, label, unit, warn-if-above
    ("n_cycles", "Cycles simulated", "-", None),
    ("dQin_cycle", "Cycle-to-cycle change of mean Q_in", "%", 2.0),
    ("dphi_cycle", "Cycle-to-cycle change of mean Phi", "%", 1.0),
    ("Qar_Qven_imbalance", "|<Q_ar> - <Q_ven>| / <Q_ar> (0D periodicity)", "%", 3.0),
    ("Qin_vs_Qar", "(<Q_in> - <Q_ar>) / <Q_ar> (3D vs 0D inflow)", "%", None),
    ("leak_net", "Net interface leak |<(u-w).n>| / <Q_in>", "%", 2.0),
    ("leak_rms", "RMS interface leak / RMS Q_in", "%", 5.0),
    ("its_mean", "Mean coupling (Newton) iterations per step", "-", None),
    ("its_max", "Max coupling (Newton) iterations per step", "-", None),
    ("gmres_mean", "Mean GMRES iterations per step", "-", None),
    ("gmres_max", "Max GMRES iterations per step", "-", None),
    ("res_max", "Max final coupling residual (last cycle)", "-", None),
]


def write_diagnostics(diag, out_dir):
    names = [n for n, _, _, _ in DIAGNOSTICS if any(n in dd for _, _, dd in diag)]
    with open(os.path.join(out_dir, "diagnostics.csv"), "w") as f:
        f.write("geometry,scenario," + ",".join(names) + "\n")
        for g, s, dd in diag:
            f.write(f"{g},{s}," + ",".join(f"{dd.get(n, np.nan):.6g}" for n in names) + "\n")
    w = max(len(l) for n, l, _, _ in DIAGNOSTICS if n in names) + 8
    head = "".join(f"{s[:16]:>18s}" for _, s, _ in diag)
    lines = ["", "  Run diagnostics (last cycle)", "  " + " " * w + head]
    for n, lab, unit, warn in DIAGNOSTICS:
        if n not in names:
            continue
        row = ""
        for _, _, dd in diag:
            v = dd.get(n, np.nan)
            flag = "!" if (warn is not None and np.isfinite(v) and abs(v) > warn) else " "
            row += f"{v:17.3g}{flag}"
        lines.append(f"  {lab + ' [' + unit + ']':<{w}s}{row}")
    lines.append("  '!' = above the usual tolerance: more cycles, smaller dt or tighter coupling.")
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(out_dir, "diagnostics.txt"), "w") as f:
        f.write(txt + "\n")


def numerics_figure(full, period, path):
    """Coupling effort and interface leakage over the whole run."""
    has_its = any("coupling_its" in d for _, _, d in full)
    has_gm = any("gmres_its" in d for _, _, d in full)
    panels = []
    if has_its:
        panels.append(("coupling_its", lambda d: d["coupling_its"],
                       "Coupling iterations\nper step", "steps"))
    if has_gm:
        panels.append(("gmres_its", lambda d: d["gmres_its"], "GMRES iterations\nper step",
                       "steps"))
    panels.append(("q_interface_slip",
                   lambda d: 100 * np.abs(d["q_interface_slip"]) / np.maximum(
                       np.abs(d["q_in"]), 1e-3 * np.max(np.abs(d["q_in"]))),
                   r"$|\int_\Gamma(\mathbf{u}-\mathbf{w})\cdot\mathbf{n}|/|Q_{in}|$ [%]",
                   "line"))
    if any("coupling_rel" in d for _, _, d in full):
        panels.append(("coupling_rel", lambda d: d["coupling_rel"],
                       "Final coupling residual", "log"))
    panels = [p for p in panels if any(p[0] in d for _, _, d in full)]
    if not panels:
        return
    ncols = 2
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(st.DOUBLE, nrows * 42 * st.MM),
                             sharex=True, squeeze=False)
    letters = "abcdefgh"
    for i, (key, fn, label, kind) in enumerate(panels):
        ax = axes.flat[i]
        for geom, scen, d in full:
            if key not in d:
                continue
            y = fn(d)
            if kind == "steps":
                # thousands of steps: rolling mean over ~period/40 with min-max band
                t = d["t"]
                dt = np.median(np.diff(t)) if len(t) > 1 else period
                w = max(1, int(round(period / 40 / max(dt, 1e-12))))
                if w > 1 and len(y) > w:
                    from numpy.lib.stride_tricks import sliding_window_view as swv
                    Y = swv(np.pad(y, (w // 2, w - 1 - w // 2), mode="edge"), w)
                    ax.fill_between(t, Y.min(1), Y.max(1), color=st.color(scen),
                                    alpha=0.15, lw=0)
                    y = Y.mean(1)
                ax.plot(t, y, color=st.color(scen), ls=st.dash(geom), lw=0.9)
            else:
                ax.plot(d["t"], y, color=st.color(scen), ls=st.dash(geom), lw=0.8)
        if kind == "log":
            ax.set_yscale("log")
        tmax = max(d["t"][-1] for _, _, d in full)
        for k in range(1, int(tmax / period) + 1):
            ax.axvline(k * period, color="0.8", lw=0.5, zorder=0)
        ax.set_ylabel(label)
        st.panel_label(ax, letters[i], x=-0.20)
        if i // ncols == nrows - 1 or i + ncols >= len(panels):
            ax.set_xlabel(r"$t$ [s]")
            ax.tick_params(labelbottom=True)
    for j in range(len(panels), nrows * ncols):
        axes.flat[j].axis("off")
    legend(fig, [(g, s, None) for g, s, _ in full], systole=False)
    fig.tight_layout(rect=(0, 0, 1, 0.92), w_pad=1.5, h_pad=0.8)
    st.save(fig, path)
    plt.close(fig)


def tex_escape(text):
    return text.replace("%", "\\%")


def write_metrics(rows, names, out_dir, extra_units=None):
    units = {n: u for n, _, u in METRICS}
    labels = {n: l for n, l, _ in METRICS}
    if extra_units:
        for n, (l, u) in extra_units.items():
            units[n], labels[n] = u, l
    csv_path = os.path.join(out_dir, "metrics_poroelastic.csv")
    with open(csv_path, "w") as f:
        f.write("geometry,scenario," + ",".join(names) + "\n")
        for geom, scen, m in rows:
            f.write(f"{geom},{scen}," + ",".join(f"{m.get(n, np.nan):.6g}" for n in names) + "\n")
    tex_path = os.path.join(out_dir, "metrics_poroelastic.tex")
    with open(tex_path, "w") as f:
        f.write("% Requires \\usepackage{booktabs}\n")
        f.write("\\begin{tabular}{ll" + "r" * len(rows) + "}\n\\toprule\n")
        heads = []
        for g, s, _ in rows:
            heads.append("\\shortstack{" + tex_escape(st.geometry_label(g)) + "\\\\"
                         + tex_escape(st.scenario_label(s)) + "}")
        f.write(" & & " + " & ".join(heads) + " \\\\\n\\midrule\n")
        for n in names:
            vals = " & ".join(f"{m.get(n, np.nan):.3g}" for _, _, m in rows)
            f.write(labels[n] + " & [" + tex_escape(units[n]) + "] & " + vals + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    print(f"  wrote {csv_path}\n  wrote {tex_path}")


def change_figure(rows, names, path, labels):
    """Relative change (%) of each metric vs. the healthy case of the same
    geometry; for the healthy cases, the stenosis effect vs. no stenosis."""
    by = {(g, s): m for g, s, m in rows}
    groups = []  # (label, colour, alpha, values)
    for g, s, m in rows:
        if s == "healthy":
            continue
        ref = by.get((g, "healthy"))
        if ref is None:
            continue
        vals = [100 * (m[n] - ref[n]) / abs(ref[n]) if ref[n] else np.nan for n in names]
        groups.append((f"{st.scenario_label(s)} ({st.geometry_label(g)})", st.color(s),
                       1.0 if g == "nostenosis" else 0.45, vals))
    if ("nostenosis", "healthy") in by and ("stenosis80", "healthy") in by:
        ref, m = by[("nostenosis", "healthy")], by[("stenosis80", "healthy")]
        vals = [100 * (m[n] - ref[n]) / abs(ref[n]) if ref[n] else np.nan for n in names]
        groups.append(("80% stenosis (healthy)", "0.45", 1.0, vals))
    if not groups:
        print("  (no healthy reference: relative-change figure skipped)")
        return
    n = len(names)
    fig, ax = plt.subplots(figsize=(st.SINGLE * 1.25, (20 + 7.5 * n) * st.MM))
    h = 0.8 / len(groups)
    ypos = np.arange(n)[::-1]
    for k, (lab, col, alpha, vals) in enumerate(groups):
        ax.barh(ypos + 0.4 - (k + 0.5) * h, vals, height=h * 0.92, color=col, alpha=alpha,
                edgecolor=col, lw=0.4, label=lab)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels([labels[x] for x in names])
    ax.tick_params(axis="y", which="both", length=0)
    ax.minorticks_off()
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.set_xlabel("Change relative to healthy [%]")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.12 - 0.05 * 7 / n), ncol=1,
              handlelength=1.4)
    fig.tight_layout()
    st.save(fig, path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results")
    ap.add_argument("figures")
    ap.add_argument("--period", type=float, default=0.85)
    ap.add_argument("--metrics3d", default=None,
                    help="metrics_3d.csv from render_3d.py, merged into the table")
    args = ap.parse_args()
    st.apply()
    os.makedirs(args.figures, exist_ok=True)

    found = discover(args.results)
    if not found:
        sys.exit(f"No {CSV_NAME} under {args.results}/<geometry>/<scenario>/")
    cases, rows, full, diag = [], [], [], []
    for geom, scen, path in found:
        print(f"  {geom}/{scen}: {path}")
        d = read_csv(path)
        c = last_cycle(d, args.period)
        cases.append((geom, scen, c))
        rows.append((geom, scen, metrics_of(c)))
        full.append((geom, scen, d))
        diag.append((geom, scen, diagnostics(d, args.period)))
    write_diagnostics(diag, args.figures)
    numerics_figure(full, args.period, os.path.join(args.figures, "fig_numerics"))

    curve_figure(cases, ["pv", "h", "phi", "pf", "qar", "nu"],
                 os.path.join(args.figures, "fig_lv_poroelastic"))
    curve_figure(cases, ["qin", "pout", "ptm", "phia", "phiv", "vol"],
                 os.path.join(args.figures, "fig_coronary_coupling"))
    loop_figure(cases, os.path.join(args.figures, "fig_phase_loops"))

    names = [n for n, _, _ in METRICS]
    labels = {n: l for n, l, _ in METRICS}
    extra = {}
    if args.metrics3d and os.path.exists(args.metrics3d):
        m3 = np.genfromtxt(args.metrics3d, delimiter=",", names=True, dtype=None,
                           encoding="utf-8")
        m3 = np.atleast_1d(m3)
        cols3 = [n for n in m3.dtype.names if n not in ("geometry", "scenario")]
        from render_3d import METRICS3D  # labels and units
        for g, s, m in rows:
            sel = m3[(m3["geometry"] == g) & (m3["scenario"] == s)]
            for n in cols3:
                m[n] = float(sel[n][0]) if len(sel) else np.nan
        for n in cols3:
            lab, unit = METRICS3D.get(n, (n, "-"))
            extra[n] = (lab, unit)
            labels[n] = lab
        names += cols3
    write_metrics(rows, names, args.figures, extra)
    key_metrics = ["Qin_mean", "Qin_dia_sys", "Qar_mean", "pf_max", "phi_mean", "phi_amp",
                   "nu_mean", "ptm_mean"]
    key_metrics += [n for n in ("TAWSS_mean", "low_TAWSS_area", "OSI_mean", "disp_amp_max")
                    if n in names]
    change_figure(rows, key_metrics, os.path.join(args.figures, "fig_scenario_changes"),
                  labels)
    print(f"  figures in {args.figures}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Publication figures for the NaturalConvectionPCM roof thickness sweep.

Reads every case OUT/e_XXmm (NaturalConvectionPCM.csv + .xdmf/.h5 written by
run_pcm_sweep.sh / run_santiago_pc9.sh), analyses one day of the periodic
regime (default: the last complete day) and writes vector (PDF) and raster
(PNG, 600 dpi) figures plus the metrics table to OUT/figures:

  fig01_weather          solar irradiance on the roof, air temperature, wind
  fig02_heat_fluxes      q_roof->PCM, q_PCM->attic, q_ceiling->dwelling vs solar time
  fig03_temperatures     roof surface, attic-side surface and attic air temperatures
  fig04_melting          mean liquid fraction and melted depth (Y e) vs solar time
  fig05_spacetime        T(eta, t) at mid-slope for every thickness + melt front
  fig06_profiles         T, liquid fraction and along-slope velocity across the layer
  fig07_slope_T_eXXmm    T over the whole sandwich along the slope at several hours (one per thickness)
  fig08_metrics          peak flux, lags, decrement factor, energies vs thickness
  metrics.csv / .md      the numbers behind fig08 (one row per thickness)
  figNN_*.csv            the data behind every figure, for ALL thicknesses
                         (fig08 -> metrics.csv)

Comparability across thicknesses. Fluxes and energies are per m^2 of roof.
Across the layer every field uses eta = (y - y_bottom)/e (0: attic side,
1: roof side), so all thicknesses share one axis. Mid-slope quantities are
averaged over the central part of the slope (--window, fraction of the
length) to remove the end effects and the cell-to-cell scatter of the
night-time convection rolls. The whole-slope maps (fig07, cases
--field-cases, hours --field-hours) share one colour scale across thicknesses
and hours. Time series show --lines thicknesses (default: five,
evenly spread) with a one-hue ramp, light = thin, dark = thick; fig05 and
fig08 show all of them.

Metrics (last day): peak and daily energy of each flux; peak lag (hour of
max q_attic minus hour of max q_roof); decrement factor f = amplitude of the
attic-side surface temperature / amplitude of the roof surface temperature,
and its lag; latent share = daily swing of latent heat stored in the PCM
(from the fields) / daily swing of the heat stored in the sandwich
(cumulative q_roof - q_attic).

    python3 figures_paper.py OUT [--day N] [--lines 5 15 25 40 50]
        [--profile-hours 9 13 17 21] [--field-hours 9 13 17 21] [--field-cases ...]

Needs numpy, matplotlib, h5py.
"""
import argparse
import csv
import glob
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

DAY = 86400.0
MM = 1.0 / 25.4  # inch per mm
W1, W2 = 90 * MM, 190 * MM  # single and double column widths

# Colours: one-hue ordinal ramp for thickness (validated: monotone lightness,
# visible steps, light end >= 2:1 on white), one accent for single series,
# neutral inks for text, grid and night shading.
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
ACCENT = "#2a78d6"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
GRID, NIGHT = "#e6e5e1", "#f1f0ec"
FRONT = "#9ec5f4"
TMAP = "inferno"


# ------------------------------------------------------------------ style ---
def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 7.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK2,
        "xtick.color": INK2, "ytick.color": INK2,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "grid.linestyle": "-",
        "axes.axisbelow": True, "lines.linewidth": 1.3, "lines.solid_capstyle": "round",
        "legend.frameon": False, "legend.handlelength": 1.6,
        "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic",
        "mathtext.bf": "sans:bold", "mathtext.default": "it", "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.dpi": 150, "savefig.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
    })


def ptitle(ax, letter, text="", pad=4):
    """Panel title: bold letter and a short description on one line, left-aligned with the axes."""
    t = "$\\mathbf{%s}$" % letter if letter else ""
    if text:
        t = (t + "    " if t else "") + text
    ax.set_title(t, loc="left", pad=pad, fontsize=7.5, color=INK)


def save(fig, out, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, "%s.%s" % (name, ext)))
    plt.close(fig)
    print("  %s.pdf/.png" % name)


def hour_axis(ax, label=True):
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 3))
    if label:
        ax.set_xlabel("Solar time (h)")


def shade_night(ax, h, I):
    """Light band where the roof gets no sun."""
    night = I <= 1.0
    start = None
    for k in range(len(h)):
        if night[k] and start is None:
            start = h[k]
        if (not night[k] or k == len(h) - 1) and start is not None:
            ax.axvspan(start, h[k], color=NIGHT, lw=0, zorder=0)
            start = None


def ramp(n):
    if n <= len(RAMP):
        idx = np.round(np.linspace(0, len(RAMP) - 1, n)).astype(int) if n > 1 else [2]
        return [RAMP[i] for i in idx]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("e", RAMP)
    return [cmap(i / (n - 1)) for i in range(n)]


# -------------------------------------------------------------------- data ---
def read_csv(path):
    with open(path) as f:
        rows = [r for r in csv.DictReader(f)]
    out = {}
    if len(rows) < 3:
        return out
    for k in rows[0]:
        try:
            out[k] = np.array([float(r[k]) for r in rows])
        except (TypeError, ValueError):
            pass
    return out


def case_info(d):
    info = {}
    p = os.path.join(d, "case.txt")
    if os.path.exists(p):
        for tok in open(p).read().split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                try:
                    info[k] = float(v)
                except ValueError:
                    info[k] = v
    if "thickness_mm" not in info:
        info["thickness_mm"] = float(re.findall(r"(\d+)mm", d)[-1])
    info.setdefault("al_mm", 0.0)
    return info


def integral(t, q):
    return float(np.sum(0.5 * (q[1:] + q[:-1]) * np.diff(t)))


def _dataitem(base, item):
    import h5py
    fname, path = item.text.strip().split(":", 1)
    with h5py.File(os.path.join(base, fname), "r") as f:
        return np.array(f[path])


def xdmf_steps(path):
    root = ET.parse(path).getroot()
    steps = []
    for step in root.find("Domain/Grid").findall("Grid"):
        steps.append((float(step.find("Time").get("Value")), step))
    return os.path.dirname(path), steps


def load_step(base, step, names=("temperature", "liquid_fraction", "velocity")):
    """Nodes (deduplicated across MPI ranks), triangles and nodal fields."""
    xy_all, tri_all, f_all, off = [], [], {n: [] for n in names}, 0
    for g in step.findall("Grid"):
        topo = g.find("Topology")
        conn = _dataitem(base, topo.find("DataItem")).astype(np.int64).ravel()
        if topo.get("TopologyType", "").lower() == "mixed":
            ne, tri, i = int(topo.get("NumberOfElements")), [], 0
            for _ in range(ne):
                typ = int(conn[i])
                nv = {4: 3, 5: 4, 1: 1}.get(typ, 3)
                if typ == 2:
                    nv = int(conn[i + 1]) + 1
                if typ == 4:
                    tri.append(conn[i + 1:i + 4])
                i += 1 + nv
            tri = np.array(tri)
        else:
            tri = conn.reshape(-1, 3)
        xy = _dataitem(base, g.find("Geometry/DataItem")).reshape(-1, 2)
        for a in g.findall("Attribute"):
            n = a.get("Name")
            if n in f_all:
                v = _dataitem(base, a.find("DataItem"))
                f_all[n].append(v.reshape(len(xy), -1) if n == "velocity" else v.ravel())
        xy_all.append(xy)
        tri_all.append(tri + off)
        off += len(xy)
    xy = np.vstack(xy_all)
    tri = np.vstack(tri_all)
    key = np.round(xy / 1e-9).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    fields = {n: np.concatenate(v)[first] for n, v in f_all.items() if v}
    return xy[first], inv[tri], fields


class Case:
    def __init__(self, d, day, props):
        self.dir = d
        self.info = case_info(d)
        self.e = self.info["thickness_mm"]
        self.al = self.info.get("al_mm", 0.0)
        self.done = os.path.exists(os.path.join(d, "DONE"))
        self.csv = read_csv(os.path.join(d, "NaturalConvectionPCM.csv"))
        t = self.csv.get("t", np.zeros(1))
        full = int(math.floor(t[-1] / DAY + 1e-6))
        self.nd = full if day is None else day
        self.ok = 1 <= self.nd <= full
        self.t0 = (self.nd - 1) * DAY
        self.props = props

    # time series of the analysis day
    def day(self, key):
        t = self.csv["t"]
        m = (t >= self.t0 - 1e-6) & (t <= self.t0 + DAY + 1e-6)
        return (t[m] - self.t0) / 3600.0, self.csv[key][m]

    def has(self, key):
        return key in self.csv

    def load_fields(self, window):
        """Hourly snapshots of the day: row averages over the mid-slope window."""
        xdmf = os.path.join(self.dir, "NaturalConvectionPCM.xdmf")
        self.snap = []
        if not os.path.exists(xdmf):
            return
        base, steps = xdmf_steps(xdmf)
        times = np.array([s[0] for s in steps])
        sel = np.where((times >= self.t0 - 1.0) & (times <= self.t0 + DAY + 1.0))[0]
        e, al = self.e * 1e-3, self.al * 1e-3
        rho_s, rho_l, lat = self.props
        for k in sel:
            xy, tri, f = load_step(base, steps[k][1])
            if "temperature" not in f:
                continue
            L = xy[:, 0].max() - xy[:, 0].min()
            xc = 0.5 * (xy[:, 0].max() + xy[:, 0].min())
            eta = (xy[:, 1] - al) / e
            pcm = (eta > -1e-6) & (eta < 1 + 1e-6)
            win = pcm & (np.abs(xy[:, 0] - xc) <= 0.5 * window * L + 1e-9)
            rows = np.round(eta[win], 6)
            ur, inv = np.unique(rows, return_inverse=True)
            cnt = np.bincount(inv)
            avg = lambda v: np.bincount(inv, weights=v[win]) / cnt  # noqa: E731
            T = f["temperature"] - 273.15
            Y = f.get("liquid_fraction", np.zeros_like(T))
            U = f.get("velocity", np.zeros((len(T), 2)))
            # latent heat per m^2 of roof, whole PCM: L (rho_s y + (rho_l - rho_s) y^2)
            tr = mtri.Triangulation(xy[:, 0], xy[:, 1], tri)
            yc = xy[tri, 1].mean(axis=1)
            cpcm = (yc > al) & (yc < al + e)
            d1, d2 = xy[tri[:, 1]] - xy[tri[:, 0]], xy[tri[:, 2]] - xy[tri[:, 0]]
            area = 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0])
            lat_node = lat * (rho_s * Y + (rho_l - rho_s) * Y * Y)
            Elat = float(np.sum(area[cpcm] * lat_node[tri[cpcm]].mean(axis=1))) / L
            self.snap.append(dict(
                h=(times[k] - self.t0) / 3600.0, eta=ur, T=avg(T), Y=avg(Y),
                us=avg(U[:, 0]), Elat=Elat, xy=xy, tri=tri, tr=tr, f=f, L=L, xc=xc))
        # keep the raw fields only where needed later (field windows)

    def snap_at(self, h):
        if not self.snap:
            return None
        k = int(np.argmin([abs(s["h"] - h) for s in self.snap]))
        return self.snap[k]


# ----------------------------------------------------------------- metrics ---
def metrics(c, ref_peak=None):
    h, qt = c.day("qTop")
    _, qb = c.day("qBottom")
    _, qc = c.day("qCeiling") if c.has("qCeiling") else (h, np.zeros_like(h))
    _, Tt = c.day("topT_C")
    _, Tb = c.day("bottomT_C")
    _, Ta = c.day("Tattic_C") if c.has("Tattic_C") else (h, np.full_like(h, np.nan))
    _, Y = c.day("liquidFraction")
    t = h * 3600.0
    it, ib = int(np.argmax(qt)), int(np.argmax(qb))
    jt, jb = int(np.argmax(Tt)), int(np.argmax(Tb))
    store = np.concatenate([[0.0], np.cumsum(0.5 * ((qt - qb)[1:] + (qt - qb)[:-1]) * np.diff(t))])
    r = dict(
        e_mm=c.e, al_mm=c.al,
        q_roof_peak=qt[it], h_q_roof_peak=h[it],
        q_attic_peak=qb[ib], h_q_attic_peak=h[ib],
        lag_q_h=(h[ib] - h[it]) % 24.0,
        E_attic_in=integral(t, np.maximum(qb, 0)) / 3600.0,
        E_attic_net=integral(t, qb) / 3600.0,
        E_dwelling_in=integral(t, np.maximum(qc, 0)) / 3600.0,
        q_dwelling_peak=float(qc.max()),
        T_roof_max=float(Tt.max()), T_attic_side_max=float(Tb.max()),
        T_attic_air_max=float(np.nanmax(Ta)),
        decrement=float((Tb.max() - Tb.min()) / max(Tt.max() - Tt.min(), 1e-9)),
        lag_T_h=(h[jb] - h[jt]) % 24.0,
        Y_max=float(Y.max()), Y_min=float(Y.min()),
        stored_swing=float(store.max() - store.min()) / 3600.0,
    )
    if c.snap:
        El = np.array([s["Elat"] for s in c.snap])
        r["latent_swing"] = float(El.max() - El.min()) / 3600.0
        r["latent_share"] = 100.0 * r["latent_swing"] / max(r["stored_swing"], 1e-9)
    else:
        r["latent_swing"] = r["latent_share"] = float("nan")
    # periodicity: daily net heat to the attic, this day vs the previous one
    t_all, qb_all = c.csv["t"], c.csv["qBottom"]
    if c.nd >= 2:
        m0 = (t_all >= c.t0 - DAY - 1e-6) & (t_all <= c.t0 + 1e-6)
        prev = integral(t_all[m0], qb_all[m0]) / 3600.0
        r["periodicity_pct"] = 100.0 * abs(r["E_attic_net"] - prev) / max(abs(r["E_attic_net"]), 1e-9)
    else:
        r["periodicity_pct"] = float("nan")
    if ref_peak:
        r["peak_reduction_pct"] = 100.0 * (1.0 - r["q_attic_peak"] / ref_peak)
    return r


COLUMNS = [
    ("e_mm", "e (mm)", "%.0f"),
    ("q_roof_peak", "peak q roof->PCM (W/m2)", "%.0f"),
    ("q_attic_peak", "peak q PCM->attic (W/m2)", "%.1f"),
    ("peak_reduction_pct", "peak reduction vs thinnest (%)", "%.0f"),
    ("h_q_attic_peak", "time of peak q attic (h)", "%.1f"),
    ("lag_q_h", "flux peak lag (h)", "%.1f"),
    ("decrement", "decrement factor (-)", "%.3f"),
    ("lag_T_h", "temperature lag (h)", "%.1f"),
    ("E_attic_in", "heat into attic (Wh/m2 day)", "%.0f"),
    ("E_dwelling_in", "heat into dwelling (Wh/m2 day)", "%.1f"),
    ("T_roof_max", "max roof surface T (C)", "%.1f"),
    ("T_attic_side_max", "max attic-side surface T (C)", "%.1f"),
    ("T_attic_air_max", "max attic air T (C)", "%.1f"),
    ("Y_max", "max liquid fraction (-)", "%.2f"),
    ("Y_min", "min liquid fraction (-)", "%.2f"),
    ("stored_swing", "daily stored heat swing (Wh/m2)", "%.0f"),
    ("latent_swing", "daily latent heat swing (Wh/m2)", "%.0f"),
    ("latent_share", "latent share (%)", "%.0f"),
    ("periodicity_pct", "day-to-day change (%)", "%.1f"),
]


def write_metrics(rows, out, site):
    with open(os.path.join(out, "metrics.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([c[0] for c in COLUMNS])
        for r in rows:
            w.writerow([("%.6g" % r[c[0]]) if isinstance(r.get(c[0]), float) else r.get(c[0]) for c in COLUMNS])
    lines = ["# PCM roof sandwich, thickness sweep: %s" % site, "",
             "| " + " | ".join(c[1] for c in COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
    for r in rows:
        lines.append("| " + " | ".join(c[2] % r.get(c[0], float("nan")) for c in COLUMNS) + " |")
    open(os.path.join(out, "metrics.md"), "w").write("\n".join(lines) + "\n")
    print("  metrics.csv/.md")


def write_table(out, name, header, rows):
    """Data behind a figure: OUT/figures/<name>.csv (one row per plotted point)."""
    with open(os.path.join(out, name + ".csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([("%.6g" % v) if isinstance(v, (float, np.floating)) else v for v in r])
    print("  %s.csv" % name)


def export_series(cases, out):
    """CSV for fig01-fig04, every thickness (the figures show --lines only)."""
    c0 = cases[0]
    h, I = c0.day("I")
    Ta = c0.day("Ta_C")[1] if c0.has("Ta_C") else np.full_like(h, np.nan)
    w = (c0.day("ha")[1] - 8.91) / 2.0 if c0.has("ha") else np.full_like(h, np.nan)
    write_table(out, "fig01_weather", ["hour", "I_roof_W_m2", "T_outdoor_C", "wind_m_s"],
                zip(h, I, Ta, w))
    get = lambda c, k: c.day(k)[1] if c.has(k) else np.full_like(c.day("t")[1], np.nan)  # noqa: E731
    rows2, rows3, rows4 = [], [], []
    for c in cases:
        h = c.day("t")[0]
        qt, qb, qc = get(c, "qTop"), get(c, "qBottom"), get(c, "qCeiling")
        Tt, Tb, Tat, Ta = get(c, "topT_C"), get(c, "bottomT_C"), get(c, "Tattic_C"), get(c, "Ta_C")
        Y = get(c, "liquidFraction")
        for k in range(len(h)):
            rows2.append((c.e, h[k], qt[k], qb[k], qc[k]))
            rows3.append((c.e, h[k], Tt[k], Tb[k], Tat[k], Ta[k]))
            rows4.append((c.e, h[k], Y[k], Y[k] * c.e))
    write_table(out, "fig02_heat_fluxes", ["e_mm", "hour", "q_roof_to_PCM_W_m2", "q_PCM_to_attic_W_m2",
                                           "q_ceiling_to_dwelling_W_m2"], rows2)
    write_table(out, "fig03_temperatures", ["e_mm", "hour", "T_roof_surface_C", "T_attic_side_surface_C",
                                            "T_attic_air_C", "T_outdoor_C"], rows3)
    write_table(out, "fig04_melting", ["e_mm", "hour", "liquid_fraction", "melted_depth_mm"], rows4)


# ----------------------------------------------------------------- figures ---
def legend_top(fig, cases, colors, y=1.0):
    """One legend for the whole figure, centred above the panels, outside the data."""
    hs = [plt.Line2D([], [], color=col, lw=1.6) for col in colors]
    fig.legend(hs, ["%g mm" % c.e for c in cases], loc="lower center", bbox_to_anchor=(0.5, y),
               ncol=len(cases), title="PCM thickness $e$", title_fontsize=7, columnspacing=1.6,
               handletextpad=0.5, handlelength=1.8, borderaxespad=0.0)


def fig_weather(c, out):
    h, I = c.day("I")
    series = [(I, "$I$ (W m$^{-2}$)", "Solar irradiance on the roof plane")]
    if c.has("Ta_C"):
        series.append((c.day("Ta_C")[1], "$T_a$ (°C)", "Outdoor air temperature"))
    if c.has("ha"):
        series.append(((c.day("ha")[1] - 8.91) / 2.0, "$w$ (m s$^{-1}$)", "Wind speed"))
    fig, ax = plt.subplots(len(series), 1, figsize=(W1, 31 * MM * len(series)), sharex=True)
    for k, (a, (v, ylab, title)) in enumerate(zip(ax, series)):
        shade_night(a, h, I)
        a.plot(h, v, color=ACCENT)
        a.set_ylabel(ylab)
        ptitle(a, "abc"[k], title)
        hour_axis(a, label=(k == len(series) - 1))
    ax[0].set_ylim(bottom=0)
    fig.align_ylabels(ax)
    fig.tight_layout(h_pad=0.8)
    save(fig, out, "fig01_weather")


def fig_fluxes(cases, colors, out):
    keys = [("qTop", "Roof surface \u2192 PCM"),
            ("qBottom", "PCM \u2192 attic"),
            ("qCeiling", "Ceiling \u2192 dwelling")]
    keys = [k for k in keys if cases[0].has(k[0])]
    fig, ax = plt.subplots(len(keys), 1, figsize=(W1, 38 * MM * len(keys)), sharex=True)
    h0, I0 = cases[0].day("I")
    for k, (a, (key, lab)) in enumerate(zip(ax, keys)):
        shade_night(a, h0, I0)
        a.axhline(0, color=INK3, lw=0.6, zorder=1)
        for c, col in zip(cases, colors):
            h, q = c.day(key)
            a.plot(h, q, color=col, zorder=3)
        a.set_ylabel("$q$ (W m$^{-2}$)")
        ptitle(a, "abc"[k], lab)
        hour_axis(a, label=(k == len(keys) - 1))
    fig.align_ylabels(ax)
    fig.tight_layout(h_pad=0.8)
    legend_top(fig, cases, colors)
    save(fig, out, "fig02_heat_fluxes")


def fig_temperatures(cases, colors, out):
    keys = [("topT_C", "Roof surface"), ("bottomT_C", "Attic-side surface of the sandwich"),
            ("Tattic_C", "Attic air")]
    keys = [k for k in keys if cases[0].has(k[0])]
    fig, ax = plt.subplots(len(keys), 1, figsize=(W1, 38 * MM * len(keys)), sharex=True)
    h0, I0 = cases[0].day("I")
    Tm = None
    for k, (a, (key, lab)) in enumerate(zip(ax, keys)):
        shade_night(a, h0, I0)
        for c, col in zip(cases, colors):
            h, T = c.day(key)
            a.plot(h, T, color=col, zorder=3)
        a.set_ylabel("$T$ (°C)")
        ptitle(a, "abc"[k], lab)
        hour_axis(a, label=(k == len(keys) - 1))
    if cases[0].has("Ta_C"):
        h, Ta = cases[0].day("Ta_C")
        ax[-1].plot(h, Ta, color=INK3, lw=0.9, zorder=2)
        k = int(np.argmax(Ta))
        ax[-1].annotate("outdoor air", xy=(h[k], Ta[k]), xytext=(0, -5), textcoords="offset points",
                        fontsize=6.5, color=INK2, ha="center", va="top")
    fig.align_ylabels(ax)
    fig.tight_layout(h_pad=0.8)
    legend_top(fig, cases, colors)
    del Tm
    save(fig, out, "fig03_temperatures")


def fig_melting(cases, colors, out):
    fig, ax = plt.subplots(2, 1, figsize=(W1, 76 * MM), sharex=True)
    h0, I0 = cases[0].day("I")
    for a in ax:
        shade_night(a, h0, I0)
    for c, col in zip(cases, colors):
        h, Y = c.day("liquidFraction")
        ax[0].plot(h, Y, color=col, zorder=3)
        ax[1].plot(h, Y * c.e, color=col, zorder=3)
    ax[0].set_ylabel("$\\bar{Y}$ (–)")
    ax[0].set_ylim(-0.02, 1.02)
    ax[1].set_ylabel("$\\bar{Y}e$ (mm)")
    ax[1].set_ylim(bottom=0)
    hour_axis(ax[0], label=False)
    hour_axis(ax[1])
    ptitle(ax[0], "a", "Mean liquid fraction of the PCM")
    ptitle(ax[1], "b", "Melted depth")
    fig.align_ylabels(ax)
    fig.tight_layout(h_pad=0.8)
    legend_top(fig, cases, colors)
    save(fig, out, "fig04_melting")


def eta_axis(a, left=True):
    a.set_ylim(0, 1)
    a.set_yticks([0, 0.5, 1])
    if left:
        a.set_ylabel("$\\eta$ (–)")
    else:
        a.tick_params(labelleft=False)


def fig_spacetime(cases, out, ts):
    cs = [c for c in cases if len(c.snap) >= 3]
    if not cs:
        return
    n = len(cs)
    ncol = min(5, n)
    nrow = int(math.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(W2, (34 + 30 * nrow) * MM), sharex=True, sharey=True,
                            squeeze=False)
    Tall = np.concatenate([np.concatenate([s["T"] for s in c.snap]) for c in cs])
    lo, hi = np.floor(np.percentile(Tall, 0.5)), np.ceil(np.percentile(Tall, 99.5))
    norm = Normalize(lo, hi)
    im = None
    for k, c in enumerate(cs):
        a = axs[k // ncol, k % ncol]
        a.grid(False)
        eta = c.snap[0]["eta"]
        hh = np.array([s["h"] for s in c.snap])
        T = np.array([np.interp(eta, s["eta"], s["T"]) for s in c.snap]).T
        Y = np.array([np.interp(eta, s["eta"], s["Y"]) for s in c.snap]).T
        im = a.pcolormesh(hh, eta, np.clip(T, lo, hi), cmap=TMAP, norm=norm, shading="gouraud",
                          rasterized=True)
        if Y.min() < 0.5 < Y.max():
            a.contour(hh, eta, Y, levels=[0.5], colors=[FRONT], linewidths=0.9)
        ptitle(a, "abcdefghijklmnop"[k] if k < 16 else "", "$e$ = %g mm" % c.e, pad=3)
        a.set_xlim(0, 24)
        a.set_xticks([0, 6, 12, 18, 24])
        eta_axis(a, left=(k % ncol == 0))
        if k // ncol == nrow - 1:
            a.set_xlabel("Solar time (h)")
    for k in range(n, nrow * ncol):
        axs[k // ncol, k % ncol].axis("off")
    write_table(out, "fig05_spacetime", ["e_mm", "hour", "eta", "T_C", "liquid_fraction", "u_s_mm_s"],
                [(c.e, s["h"], s["eta"][i], s["T"][i], s["Y"][i], 1e3 * s["us"][i])
                 for c in cs for s in c.snap for i in range(len(s["eta"]))])
    cb = fig.colorbar(im, ax=axs.ravel().tolist(), fraction=0.025, pad=0.015, aspect=30)
    cb.set_label("$T$ (°C), mid-slope average")
    step = 10 if hi - lo > 30 else (5 if hi - lo > 12 else 2)
    cb.set_ticks(np.arange(math.ceil(lo / step) * step, hi + 0.1, step))
    cb.outline.set_linewidth(0.5)
    save(fig, out, "fig05_spacetime")


def fig_profiles(cases, colors, hours, out, all_cases=None):
    cs = [(c, col) for c, col in zip(cases, colors) if c.snap]
    if not cs:
        return
    nh = len(hours)
    fig, axs = plt.subplots(3, nh, figsize=(W2, 128 * MM), sharey=True, squeeze=False)
    quantities = [("T", "$T$ (°C)", 1.0), ("Y", "$Y$ (–)", 1.0),
                  ("us", "$u_s$ (mm s$^{-1}$)", 1e3)]
    for j, hh in enumerate(hours):
        for i, (key, lab, scale) in enumerate(quantities):
            a = axs[i, j]
            if key == "us":
                a.axvline(0, color=INK3, lw=0.6)
            for c, col in cs:
                s = c.snap_at(hh)
                a.plot(s[key] * scale, s["eta"], color=col)
            eta_axis(a, left=(j == 0))
            a.set_xlabel(lab, labelpad=1.5)
            if key == "Y":
                a.set_xlim(-0.05, 1.05)
            s = cs[0][0].snap_at(hh)
            letter = "abcdefghijkl"[i * nh + j] if nh <= 4 else ""
            ptitle(a, letter, "%02d:%02d" % (int(s["h"]), int(round((s["h"] % 1) * 60)) % 60) if i == 0 else "")
    # same velocity range in a row so that magnitudes compare across hours
    vmax = max(abs(v) for a in axs[2] for v in a.get_xlim())
    for a in axs[2]:
        a.set_xlim(-vmax, vmax)
    tmin = min(a.get_xlim()[0] for a in axs[0])
    tmax = max(a.get_xlim()[1] for a in axs[0])
    for a in axs[0]:
        a.set_xlim(tmin, tmax)
    pass
    rows = []
    for c in (all_cases or [c for c, _ in cs]):
        if not c.snap:
            continue
        for hh in hours:
            s = c.snap_at(hh)
            rows += [(c.e, hh, s["h"], s["eta"][i], s["T"][i], s["Y"][i], 1e3 * s["us"][i])
                     for i in range(len(s["eta"]))]
    write_table(out, "fig06_profiles", ["e_mm", "hour_requested", "hour", "eta", "T_C", "liquid_fraction",
                                        "u_s_mm_s"], rows)
    fig.tight_layout(h_pad=1.2, w_pad=0.8)
    legend_top(fig, [c for c, _ in cs], [col for _, col in cs])
    save(fig, out, "fig06_profiles")


def fig_slope(cases, hours, out):
    """Temperature over the whole sandwich (plates + PCM) at several hours, one figure per thickness.

    Common colour scale for every thickness and hour; vertical scale exaggerated.
    """
    cs = [c for c in cases if c.snap]
    if not cs:
        return
    snaps = {(c.e, hh): c.snap_at(hh) for c in cs for hh in hours}
    Tall = np.concatenate([s["f"]["temperature"] - 273.15 for s in snaps.values()])
    lo, hi = np.floor(np.percentile(Tall, 0.2)), np.ceil(np.percentile(Tall, 99.8))
    levels = np.linspace(lo, hi, 65)
    step = 10 if hi - lo > 30 else (5 if hi - lo > 12 else 2)
    for c in cs:
        nh = len(hours)
        fig, axs = plt.subplots(nh, 1, figsize=(W2, nh * 19 * MM + 12 * MM), sharex=True, squeeze=False)
        axs = axs[:, 0]
        fig.subplots_adjust(hspace=0.45, left=0.07, right=0.9)
        rows, cf = [], None
        for k, (a, hh) in enumerate(zip(axs, hours)):
            s = snaps[(c.e, hh)]
            xy, tri, f = s["xy"], s["tri"], s["f"]
            T = f["temperature"] - 273.15
            y_mm = xy[:, 1] * 1e3
            tr = mtri.Triangulation(xy[:, 0], y_mm, tri)
            cf = a.tricontourf(tr, np.clip(T, lo, hi), levels=levels, cmap=TMAP)
            cf.set_rasterized(True)
            if "liquid_fraction" in f:
                Y = f["liquid_fraction"]
                yc = y_mm[tri].mean(axis=1)
                pcm = (yc > c.al) & (yc < c.al + c.e)
                trp = mtri.Triangulation(xy[:, 0], y_mm, tri, mask=~pcm)
                if Y.min() < 0.5 < Y.max():
                    a.tricontour(trp, Y, levels=[0.5], colors=[FRONT], linewidths=0.9)
            else:
                Y = np.full_like(T, np.nan)
            if c.al > 0:
                for yl in (c.al, c.al + c.e):
                    a.axhline(yl, color="white", lw=0.5, alpha=0.6)
            a.grid(False)
            for sp_ in ("top", "right"):
                a.spines[sp_].set_visible(False)
            a.set_xlim(xy[:, 0].min(), xy[:, 0].max())
            a.set_ylim(0, y_mm.max())
            a.set_yticks([0, round(y_mm.max())])
            if k < nh - 1:
                a.tick_params(axis="x", length=0)
            ptitle(a, "abcdefgh"[k] if k < 8 else "",
                   "%02d:%02d solar time" % (int(s["h"]), int(round((s["h"] % 1) * 60)) % 60))
            rows += [(c.e, s["h"], xy[n, 0], y_mm[n], T[n], Y[n]) for n in range(len(T))]
        axs[-1].set_xlabel("Distance along the slope (m)")
        fig.supylabel("$y$ (mm)", fontsize=7.5, x=0.015)
        cb = fig.colorbar(cf, ax=list(axs), fraction=0.02, pad=0.015, aspect=25)
        cb.set_label("$T$ (°C)")
        cb.set_ticks(np.arange(math.ceil(lo / step) * step, hi + 0.1, step))
        cb.outline.set_linewidth(0.5)
        name = "fig07_slope_T_e%02dmm" % int(round(c.e))
        write_table(out, name, ["e_mm", "hour", "x_m", "y_mm", "T_C", "liquid_fraction"], rows)
        save(fig, out, name)


def fig_metrics(rows, out):
    e = np.array([r["e_mm"] for r in rows])
    g = lambda k: np.array([r.get(k, np.nan) for r in rows], dtype=float)  # noqa: E731
    panels = [
        ("q_attic_peak", "Peak flux to the attic", "W m$^{-2}$"),
        ("lag_q_h", "Flux peak lag", "h"),
        ("decrement", "Decrement factor", "$f$ (–)"),
        ("lag_T_h", "Temperature peak lag", "h"),
        ("E_attic_in", "Daily heat to the attic", "Wh m$^{-2}$"),
        ("E_dwelling_in", "Daily heat to the dwelling", "Wh m$^{-2}$"),
        ("T_attic_air_max", "Peak attic air temperature", "°C"),
        ("latent_share", "Latent share of stored heat", "%"),
    ]
    panels = [p for p in panels if np.isfinite(g(p[0])).any()]
    nc = 4
    nr = int(math.ceil(len(panels) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(W2, 47 * MM * nr), squeeze=False, sharex=True)
    for k, (key, lab, unit) in enumerate(panels):
        a = axs[k // nc, k % nc]
        v = g(key)
        a.plot(e, v, color=ACCENT, lw=1.2, zorder=2)
        a.plot(e, v, "o", ms=4, mfc=ACCENT, mec="white", mew=0.8, zorder=3)
        ptitle(a, "abcdefgh"[k], lab)
        a.set_ylabel(unit)
        a.set_xlim(0, max(e) * 1.06)
        a.margins(y=0.15)
        if k // nc == nr - 1:
            a.set_xlabel("PCM thickness $e$ (mm)")
    for k in range(len(panels), nr * nc):
        axs[k // nc, k % nc].axis("off")
    fig.align_ylabels(axs[:, 0])
    fig.tight_layout(h_pad=1.6, w_pad=1.4)
    save(fig, out, "fig08_metrics")


def write_captions(out, site, lines, a, day):
    e_lines = ", ".join("%g" % c.e for c in lines)
    caps = [
        ("fig01_weather", "Boundary conditions of the analysed day (%s): (a) solar irradiance on the roof plane, "
         "(b) outdoor air temperature and (c) wind speed. Shaded: night." % site),
        ("fig02_heat_fluxes", "Heat fluxes per unit roof area over the periodic day (day %d) for e = %s mm: "
         "(a) from the outer aluminium surface into the sandwich, (b) from the sandwich into the attic and "
         "(c) through the ceiling into the dwelling. Positive values point inwards. Shaded: night." % (day, e_lines)),
        ("fig03_temperatures", "Temperatures over the periodic day for e = %s mm: (a) roof surface, (b) attic-side "
         "surface of the sandwich and (c) attic air (grey: outdoor air)." % e_lines),
        ("fig04_melting", "(a) Mean liquid fraction of the PCM and (b) melted depth (mean liquid fraction times "
         "thickness) over the periodic day for e = %s mm." % e_lines),
        ("fig05_spacetime", "Temperature across the PCM layer over the periodic day for every thickness, averaged over "
         "the central %.0f %% of the slope. The coordinate across the layer is normalised by the thickness "
         "(eta = 0: attic side, 1: roof side); the blue line is the melt front (Y = 0.5). Common colour scale."
         % (100 * a.window)),
        ("fig06_profiles", "Profiles across the PCM layer at %s solar time, averaged over the central %.0f %% of the "
         "slope: (a-d) temperature, (e-h) liquid fraction and (i-l) along-slope velocity (positive up-slope), for "
         "e = %s mm." % (", ".join("%02d:00" % h for h in a.profile_hours), 100 * a.window, e_lines)),
        ("fig07_slope_T_eXXmm", "Temperature in the roof sandwich (aluminium plates and PCM; thin white lines: "
         "plate/PCM interfaces; y = 0 is the attic side) along the whole slope at %s solar time, for a thickness of XX mm. Vertical scale "
         "exaggerated; blue line: melt front (Y = 0.5). Same colour scale for every thickness and hour."
         % ", ".join("%02d:00" % h for h in a.field_hours)),
        ("fig08_metrics", "Thermal performance versus PCM thickness over the periodic day: (a) peak heat flux into the "
         "attic, (b) lag of that peak behind the peak flux entering the roof, (c) decrement factor (amplitude of the "
         "attic-side surface temperature over that of the roof surface) and (d) its lag, (e) heat entering the attic "
         "and (f) the dwelling per day, (g) maximum attic air temperature and (h) latent share of the daily swing of "
         "stored heat."),
    ]
    with open(os.path.join(out, "captions.md"), "w") as f:
        f.write("# Draft captions\n\n")
        for n, c in caps:
            f.write("**%s.** %s\n\n" % (n, c))
    print("  captions.md")


# -------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", nargs="?", default="pcm_santiago")
    ap.add_argument("--day", type=int, default=None, help="day analysed (default: last complete day)")
    ap.add_argument("--lines", type=float, nargs="*", default=None,
                    help="thicknesses (mm) in the time-series and profile figures (default: 5 evenly spread)")
    ap.add_argument("--profile-hours", type=float, nargs="+", default=[9, 13, 17, 21])
    ap.add_argument("--field-hours", type=float, nargs="+", default=[9, 13, 17, 21],
                    help="solar hours of the whole-slope temperature maps (fig07)")
    ap.add_argument("--field-cases", type=float, nargs="*", default=None,
                    help="thicknesses (mm) in fig07 (default: same as --lines)")
    ap.add_argument("--window", type=float, default=0.2,
                    help="mid-slope averaging window, fraction of the slope length")
    ap.add_argument("--rho-s", type=float, default=880.0)
    ap.add_argument("--rho-l", type=float, default=760.0)
    ap.add_argument("--latent", type=float, default=165e3)
    ap.add_argument("--no-fields", action="store_true", help="CSV-based figures only (fast)")
    a = ap.parse_args()

    dirs = sorted(d for d in glob.glob(os.path.join(a.out, "e_*mm"))
                  if os.path.exists(os.path.join(d, "NaturalConvectionPCM.csv")))
    if not dirs:
        sys.exit("No cases with NaturalConvectionPCM.csv in %s" % a.out)
    site = "?"
    sw = os.path.join(a.out, "sweep.txt")
    if os.path.exists(sw):
        txt = open(sw).read()
        m = re.search(r"SITE=(\w+)", txt)
        site = m.group(1).capitalize() if m else "?"
        w = re.search(r"WEATHER=(\S+)", txt)
        if w and os.path.exists(w.group(1)):
            head = open(w.group(1)).readline()
            p = re.search(r"#\s*(.*?),\s*(.*?) month \d+ of PVGIS.*?real year (\w+)", head)
            if p:
                what = {"mean": "mean January day"}.get(p.group(2).split(" of")[0], "January " + p.group(2))
                site = "%s, %s (PVGIS TMY, %s)" % (p.group(1), what, p.group(3))
            else:
                site = head.lstrip("# ").split(":")[0].strip()
    figs = os.path.join(a.out, "figures")
    os.makedirs(figs, exist_ok=True)
    style()

    props = (a.rho_s, a.rho_l, a.latent)
    cases = []
    for d in dirs:
        c = Case(d, a.day, props)
        if not c.ok:
            print("  skip %s: day %s not complete (t = %.1f h)" % (d, c.nd, c.csv.get("t", np.zeros(1))[-1] / 3600))
            continue
        if not c.done:
            print("  note: %s has no DONE (still running?); using day %d" % (d, c.nd))
        cases.append(c)
    if not cases:
        sys.exit("No case has a complete day")
    cases.sort(key=lambda c: c.e)
    es = [c.e for c in cases]
    print("Cases: %s mm; day %s" % (", ".join("%g" % e for e in es), ", ".join(str(c.nd) for c in cases)))

    if a.lines:
        lines = [c for c in cases if any(abs(c.e - x) < 1e-6 for x in a.lines)]
    else:
        idx = sorted(set(np.round(np.linspace(0, len(cases) - 1, min(5, len(cases)))).astype(int)))
        lines = [cases[i] for i in idx]
    colors = ramp(len(lines))
    fcases = lines if a.field_cases is None else [c for c in cases if any(abs(c.e - x) < 1e-6 for x in a.field_cases)]

    if not a.no_fields:
        for c in cases:
            print("  reading fields e = %g mm ..." % c.e, flush=True)
            c.load_fields(a.window)
    else:
        for c in cases:
            c.snap = []

    rows = []
    for c in cases:
        rows.append(metrics(c))
    ref = rows[0]["q_attic_peak"]
    for r in rows:
        r["peak_reduction_pct"] = 100.0 * (1.0 - r["q_attic_peak"] / ref) if ref else float("nan")

    print("Writing %s" % figs)
    write_metrics(rows, figs, site)
    export_series(cases, figs)
    fig_weather(cases[0], figs)
    fig_fluxes(lines, colors, figs)
    fig_temperatures(lines, colors, figs)
    fig_melting(lines, colors, figs)
    fig_metrics(rows, figs)
    if not a.no_fields:
        fig_spacetime(cases, figs, None)
        fig_profiles(lines, colors, a.profile_hours, figs, all_cases=cases)
        fig_slope(fcases, a.field_hours, figs)
    write_captions(figs, site, lines, a, cases[0].nd)


if __name__ == "__main__":
    main()

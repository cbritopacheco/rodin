#!/usr/bin/env python3
"""
atrium_post.py -- figuras para el artículo (JFM) de las corridas del atrio.

    python3 atrium_post.py [--dir ~/runs/fa] [--runs N:prod_N,CY:prod_CY,VE:prod_VE,TVE:prod_TVE]
                           [--what all] [--phases auto] [--skip 2] [--out figs]

Lee la salida de KokkosAtriumRheologyPETSc (XDMF + HDF5, varias partes si la
corrida se retomó) y escribe en --out, en PDF (vectorial) y PNG:

  material   funciones del material: eta_CY, lambda_st y t_s contra la tasa de
             corte, con el periodo marcado (no necesita corridas);
  waveforms  caudal mitral, presiones y energía cinética en el ciclo: media de
             fase y banda de +-1 desviación entre ciclos, por fluido;
  slices     corte por la orejuela y el plano mitral: |u| medio de fase con
             líneas de corriente del campo en el plano, y la tasa de corte;
  elastic    en el mismo corte, para VE y TVE: esfuerzo elástico / esfuerzo
             inercial, Weissenberg local y tamaño de agregado (TVE);
  laa        acercamiento a la orejuela: tasa de corte y fracción del ciclo con
             gamma < 1 1/s (estasis), por fluido;
  stasis     mapa de estasis en el corte: fracción del ciclo con gamma < 1 1/s;
  q3d        superficies del criterio Q (mismo umbral para todos los fluidos),
             coloreadas por |u|, dentro de la pared translúcida;
  energy     balance de energía por ciclo: cinética, elástica, disipación
             viscosa y potencia elástica (energy.csv);
  wall       TAWSS, OSI y ECAP = OSI/TAWSS en la pared, con el esfuerzo elástico
             incluido en VE y TVE (exp(psi)); mapas en la orejuela y
             wall_metrics.csv con medias en la orejuela y en toda la pared;
  residence  lavado de la orejuela: partículas sembradas en ella al inicio de
             los ciclos analizados, advectadas con los campos guardados
             (RK2, --substeps por cuadro); fracción que sigue en la orejuela y
             en el atrio, vida media (residence.csv);
  drift      <n> de TVE (y VE) en la orejuela cuadro a cuadro: si sigue
             derivando, la orejuela no está en régimen;
  laa        incluye además la viscosidad efectiva eta_eff/eta_inf (para VE y
             TVE, eta_pl + G lambda).

Todo lo de campos es media de fase sobre los ciclos usados (los mismos que
stats_models.py: ciclos completos comunes a todas las corridas, menos --skip),
salvo --instant, que usa el último ciclo.

Requisitos: pip install pyvista matplotlib h5py scipy  (y ffmpeg solo para el video).
"""
import argparse
import csv
import os
import re
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
try:
    import h5py
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors
    import pyvista as pv
except ImportError as e:
    sys.exit(f"Falta un paquete ({e.name}): pip install pyvista matplotlib h5py")

RHO = 1060.0
# Sangre de referencia (los mismos valores que el driver).
ETA0, ETAINF, LAMCY, NCY, ACY = 0.072, 0.0049, 6.1529, 0.2081, 1.4173
ETAPL, LAMH = 1.2e-3, 5e-3
G = (ETAINF - ETAPL) / LAMH
A0, AC, GDC, GDMAX = 1.0, 2.0, 5.78, 900.0
FLUIDS = ("N", "CY", "VE", "TVE")
COLORS = {"N": "#1b1b1b", "CY": "#0072b2", "VE": "#009e73", "TVE": "#d55e00"}

# ------------------------------------------------------------------ estilo JFM
JFM_WIDTH = 5.3          # ancho de texto de JFM en pulgadas


def style(usetex):
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "text.usetex": usetex,
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.linewidth": 0.6, "lines.linewidth": 1.0,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "legend.frameon": False, "savefig.dpi": 400, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02, "image.interpolation": "bilinear",
    })


def panel(ax, k, color="k", x=0.0, y=1.01):
    """Rótulo (a), (b), ... sobre la esquina superior izquierda, fuera de los datos."""
    ax.text(x, y, f"({chr(97 + k)})", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9, color="k")


def save(fig, out, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  {name}.pdf")


# ------------------------------------------------------------------ material
def eta_cy(g):
    return ETAINF + (ETA0 - ETAINF) * (1 + (LAMCY * g) ** ACY) ** ((NCY - 1) / ACY)


def smooth(s):
    return s * s * (3 - 2 * s)


def agg_rate(g):
    g = np.asarray(g, float)
    return np.where(g <= GDC, A0 + (AC - A0) * smooth(np.clip(g / GDC, 0, 1)),
                    AC * (1 - smooth(np.clip((g - GDC) / (GDMAX - GDC), 0, 1))))


def n_st(g):
    r = np.maximum((eta_cy(g) - ETAPL) / (ETAINF - ETAPL), 1)
    return r * (1 + 1.5 * agg_rate(g) * LAMH)


def t_struct(g):
    ns = n_st(g)
    b = agg_rate(g) / np.maximum(ns * (ns - 1), 1e-300)
    rate = 0.5 * b * (2 * ns - 1)
    return np.where(rate > 0, 1 / np.maximum(rate, 1e-300), np.inf)


# ------------------------------------------------------------------ corridas
def parts(p):
    if os.path.exists(p + ".parts"):
        return [x.strip() for x in open(p + ".parts") if x.strip()]
    return [p]


class Run:
    """Una corrida (todas sus partes): serie temporal y campos por (ciclo, fase)."""

    def __init__(self, label, prefix, period):
        self.label, self.prefix, self.T = label, prefix, period
        byt = {}
        for q in parts(prefix):                       # la parte posterior manda
            if not os.path.exists(q + ".csv"):
                sys.exit(f"No encuentro {q}.csv")
            for r in csv.DictReader(open(q + ".csv")):
                byt[round(float(r["t"]), 9)] = r
        rows = [byt[t] for t in sorted(byt)]
        self.S = {k: np.array([float(r.get(k) or 0) for r in rows]) for k in rows[0]}
        self.frames = {}
        for q in parts(prefix):
            times = [float(t) for t in re.findall(r'Time Value="([^"]+)"', open(q + ".xdmf").read())]
            for k, t in enumerate(times):
                if t > 0 and os.path.exists(f"{q}.velocity.{k:06d}.h5"):
                    c = int(np.ceil(t / period - 1e-6))
                    self.frames[(c, round(t - (c - 1) * period, 6))] = (q, k)
        self.mesh_file = parts(prefix)[0] + ".mesh.h5"

    def field(self, name, c, ph):
        q, k = self.frames[(c, ph)]
        f = f"{q}.{name}.{k:06d}.h5"
        if not os.path.exists(f):
            return None
        return h5py.File(f, "r")["GridFunction/Values/Data"][:]

    def has(self, name):
        q, k = next(iter(self.frames.values()))
        return os.path.exists(f"{q}.{name}.{k:06d}.h5")


def common_cycles(runs, period, skip):
    steps = None
    full = []
    for r in runs.values():
        cyc, cnt = np.unique(r.S["cycle"], return_counts=True)
        steps = steps or int(cnt.max())
        full.append({int(c) for c, n in zip(cyc, cnt) if n == steps})
    common = sorted(set.intersection(*full))
    cycles = common[skip:]
    if not cycles:
        sys.exit(f"Ciclos completos comunes {common}; con skip={skip} no queda ninguno.")
    phases = sorted({ph for (c, ph) in next(iter(runs.values())).frames if c in cycles})
    phases = [ph for ph in phases if all((c, ph) in r.frames for r in runs.values() for c in cycles)]
    return cycles, phases


def read_mesh(run):
    m = h5py.File(run.mesh_file, "r")
    X = m["Mesh/Geometry/Vertices"][:]
    unit = 1e-3 if np.ptp(X, axis=0).max() > 1.0 else 1.0
    T = m["Mesh/XDMF/Topology"][:].reshape(-1, 5)[:, 1:].astype(np.int64)
    cells = np.hstack([np.full((len(T), 1), 4), T]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(T), pv.CellType.TETRA, np.uint8), X * unit)
    grid.field_data["T"] = T.ravel()
    return grid


def tets(mesh):
    return np.asarray(mesh.field_data["T"]).reshape(-1, 4).astype(np.int64)


def lumped_volume(mesh):
    T, X = tets(mesh), np.asarray(mesh.points)
    E = X[T[:, 1:]] - X[T[:, [0]]]
    vol = np.abs(np.linalg.det(E)) / 6
    lump = np.zeros(len(X))
    np.add.at(lump, T.ravel(), np.repeat(vol / 4, 4))
    return lump


def conformation(run, c, ph):
    """c = exp(psi) en los vértices, (n, 3, 3)."""
    d, o = run.field("psi_diagonal", c, ph), run.field("psi_off_diagonal", c, ph)
    if d is None or o is None:
        return None
    P = np.zeros((len(d), 3, 3))
    P[:, [0, 1, 2], [0, 1, 2]] = d
    P[:, 0, 1] = P[:, 1, 0] = o[:, 0]
    P[:, 0, 2] = P[:, 2, 0] = o[:, 1]
    P[:, 1, 2] = P[:, 2, 1] = o[:, 2]
    w, R = np.linalg.eigh(P)
    return np.einsum("nij,nj,nkj->nik", R, np.exp(w), R)


def read_medit(path):
    """Vértices (m) y triángulos con etiqueta de un .mesh de medit, sin leer los tetraedros."""
    V, F = None, None
    with open(path) as f:
        it = iter(f)
        for line in it:
            key = line.strip()
            if key in ("Vertices", "Triangles", "Tetrahedra", "Edges"):
                n = int(next(it).split()[0])
                if key == "Vertices":
                    V = np.loadtxt([next(it) for _ in range(n)])[:, :3]
                elif key == "Triangles":
                    F = np.loadtxt([next(it) for _ in range(n)], dtype=np.int64)
                else:
                    for _ in range(n):
                        next(it)
            if V is not None and F is not None:
                break
    unit = 1e-3 if np.ptp(V, axis=0).max() > 1.0 else 1.0
    return V * unit, F


# ------------------------------------------------------------------ geometría
class Setting:
    def __init__(self, mesh_path, laa):
        V, F = read_medit(mesh_path)
        self.V, self.F = V, F
        lab = F[:, 3]
        self.openings = {}
        for l in np.unique(lab):
            if l == 4:
                continue
            ids = np.unique(F[lab == l, :3] - 1)
            self.openings[int(l)] = V[ids].mean(axis=0)
        self.mitral = self.openings.get(6)
        tri = F[:, :3] - 1
        a = 0.5 * np.linalg.norm(np.cross(V[tri[:, 1]] - V[tri[:, 0]], V[tri[:, 2]] - V[tri[:, 0]]), axis=1)
        self.mitral_area = a[lab == 6].sum()
        self.D_m = 2 * np.sqrt(self.mitral_area / np.pi)
        faces = np.hstack([np.full((len(tri), 1), 3), tri]).ravel()
        surf = pv.PolyData(V, faces)
        surf.cell_data["label"] = lab
        self.wall = surf.extract_cells(lab == 4).extract_surface().clean()
        self.laa_c = np.array(laa[:3]) * 1e-3
        self.laa_r = laa[3] * 1e-3
        veins = [self.openings[l] for l in (1, 2, 3, 5) if l in self.openings]
        p0, p1, p2 = self.laa_c, self.mitral, np.mean(veins, axis=0)
        # Plano por la orejuela, el centro mitral y el centro de las venas.
        e1 = p1 - p0
        e1 /= np.linalg.norm(e1)
        nrm = np.cross(e1, p2 - p0)
        nrm /= np.linalg.norm(nrm)
        self.e1, self.e2, self.normal, self.origin = e1, np.cross(nrm, e1), nrm, p0

    def plane_grid(self, mesh, h=2.5e-4, window=None):
        """Malla regular en el plano (coordenadas s, r en mm desde la orejuela)."""
        P = mesh.points - self.origin
        s, r = P @ self.e1, P @ self.e2
        near = np.abs(P @ self.normal) < 3e-3
        if window is None:
            s0, s1, r0, r1 = s[near].min(), s[near].max(), r[near].min(), r[near].max()
        else:
            s0, s1, r0, r1 = window
        S, R = np.meshgrid(np.arange(s0, s1 + h, h), np.arange(r0, r1 + h, h))
        pts = self.origin + S[..., None] * self.e1 + R[..., None] * self.e2
        probe = pv.PolyData(pts.reshape(-1, 3)).sample(mesh)
        mask = probe["vtkValidPointMask"].reshape(S.shape).astype(bool)
        return S * 1e3, R * 1e3, probe, mask

    def project(self, probe, name, shape, mask):
        a = np.asarray(probe[name])
        if a.ndim == 2:
            out = np.stack([np.ma.array((a @ e).reshape(shape), mask=~mask) for e in (self.e1, self.e2)])
            return out
        return np.ma.array(a.reshape(shape), mask=~mask)


# ------------------------------------------------------------------ figuras
def fig_material(out, period):
    g = np.logspace(-2, np.log10(500), 400)
    lam = (eta_cy(g) - ETAPL) / G
    ts = t_struct(g)
    fig, ax = plt.subplots(1, 2, figsize=(JFM_WIDTH, 2.1), constrained_layout=True)
    ax[0].loglog(g, 1e3 * eta_cy(g), color=COLORS["CY"], label=r"$\eta_{CY}$")
    ax[0].axhline(1e3 * ETAINF, color=COLORS["N"], ls="--", label=r"$\eta_\infty$ (N)")
    ax[0].axhline(1e3 * ETAPL, color="0.5", ls=":", label=r"$\eta_{pl}$")
    ax[0].set_xlabel(r"$\dot\gamma$ (s$^{-1}$)")
    ax[0].set_ylabel(r"$\eta$ (mPa s)")
    ax[0].legend(loc="upper right")
    ax[1].loglog(g, lam, color=COLORS["VE"], label=r"$\lambda_{\rm st}$")
    ax[1].loglog(g, ts, color=COLORS["TVE"], label=r"$t_s$")
    ax[1].axhline(period, color="k", lw=0.6)
    ax[1].text(1.5e-2, period * 1.25, r"$T$", fontsize=8)
    ax[1].fill_between(g, 1e-3, 1e3, where=ts > period, color=COLORS["TVE"], alpha=0.08, lw=0)
    ax[1].set_ylim(2e-3, 1e2)
    ax[1].set_xlabel(r"$\dot\gamma$ (s$^{-1}$)")
    ax[1].set_ylabel("time (s)")
    ax[1].legend(loc="upper right")
    for k, a in enumerate(ax):
        panel(a, k)
    save(fig, out, "material")


def fig_waveforms(out, runs, cycles, period, phases_mark):
    fig, ax = plt.subplots(3, 1, figsize=(JFM_WIDTH * 0.62, 3.9), sharex=True, constrained_layout=True)
    for lab, r in runs.items():
        S = r.S
        sel = np.isin(S["cycle"].astype(int), cycles)
        tau = (S["t"][sel] - (S["cycle"][sel] - 1) * period) / period
        bins = np.round(tau * 1e4).astype(int)
        for i, (key, scale) in enumerate((("qOut", 1e6), ("pOut", 1e-3), ("kineticEnergy", 1e3))):
            v = S[key][sel] * scale
            ub = np.unique(bins)
            m = np.array([v[bins == b].mean() for b in ub])
            s = np.array([v[bins == b].std() for b in ub])
            x = ub / 1e4
            ax[i].plot(x, m, color=COLORS.get(lab, None), label=lab, lw=0.9)
            ax[i].fill_between(x, m - s, m + s, color=COLORS.get(lab, None), alpha=0.18, lw=0)
            if i == 1 and lab == list(runs)[0]:
                pv_ = S["pVenous"][sel] * 1e-3
                ax[1].plot(x, [pv_[bins == b].mean() for b in ub], color="0.45", ls="--", lw=0.8,
                           label=r"$p_{\rm pv}$")
    ax[0].set_ylabel(r"$Q_m$ (mL s$^{-1}$)")
    ax[1].set_ylabel(r"$p$ (kPa)")
    ax[2].set_ylabel(r"$E_k$ (mJ)")
    ax[2].set_xlabel(r"$t/T$")
    ax[2].set_xlim(0, 1)
    for a in ax:
        for ph in phases_mark:
            a.axvline(ph / period, color="0.6", lw=0.5, ls=":")
    ax[0].legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.42), handlelength=1.4)
    from matplotlib.lines import Line2D
    ax[1].legend([Line2D([], [], color="0.3"), Line2D([], [], color="0.45", ls="--")],
                 [r"$p_v$", r"$p_{\rm pv}$"], loc="upper right", handlelength=1.6)
    for k, a in enumerate(ax):
        panel(a, k)
    save(fig, out, "waveforms")


def phase_mean(run, names, cycles, ph):
    acc = {}
    for c in cycles:
        for nm in names:
            a = run.field(nm, c, ph)
            if a is None:
                continue
            acc[nm] = acc.get(nm, 0) + a / len(cycles)
    return acc


def stasis_fraction(run, cycles, phases, threshold=1.0):
    acc = 0
    for c in cycles:
        for ph in phases:
            acc = acc + (run.field("shear_rate", c, ph) < threshold)
    return acc / (len(cycles) * len(phases))


def attach(mesh, fields):
    m = mesh.copy(deep=False)
    for k, v in fields.items():
        m.point_data[k] = v
    return m


def outline(ax, mask, S, R, lw=0.5):
    ax.contour(S, R, mask.astype(float), levels=[0.5], colors="k", linewidths=lw)


def imshow(ax, S, R, Z, **kw):
    return ax.pcolormesh(S, R, Z, shading="gouraud", rasterized=True, **kw)


def fig_slices(out, runs, mesh, st, cycles, ph, tag, instant):
    use = [cycles[-1]] if instant else cycles
    S, R, probe0, mask = st.plane_grid(mesh)
    data = {}
    for lab, r in runs.items():
        f = phase_mean(r, ["velocity", "shear_rate"], use, ph)
        pr = pv.PolyData(probe0.points).sample(attach(mesh, f))
        uv = st.project(pr, "velocity", S.shape, mask)
        sp = np.ma.array(np.linalg.norm(pr["velocity"], axis=1).reshape(S.shape), mask=~mask)
        gd = st.project(pr, "shear_rate", S.shape, mask)
        data[lab] = (uv, sp, gd)
    vmax = max(float(np.percentile(d[1].compressed(), 99.5)) for d in data.values())
    ncol = len(runs)
    aspect = (R.max() - R.min()) / (S.max() - S.min())
    colw = (JFM_WIDTH - 0.8) / ncol
    fig, ax = plt.subplots(2, ncol, figsize=(JFM_WIDTH, 2 * colw * aspect + 0.45), constrained_layout=True,
                           squeeze=False)
    for j, (lab, (uv, sp, gd)) in enumerate(data.items()):
        a = ax[0, j]
        im0 = imshow(a, S, R, sp, cmap="inferno", vmin=0, vmax=vmax)
        a.streamplot(S, R, uv[0].filled(np.nan), uv[1].filled(np.nan), density=1.1, color="w",
                     linewidth=0.35, arrowsize=0.45, broken_streamlines=True)
        outline(a, mask, S, R)
        a.set_title(lab)
        b = ax[1, j]
        im1 = imshow(b, S, R, np.ma.maximum(gd, 1e-2), cmap="cividis",
                     norm=mcolors.LogNorm(vmin=0.1, vmax=300))
        outline(b, mask, S, R)
        for x in (a, b):
            x.set_aspect("equal")
            x.set_xticks([])
            x.set_yticks([])
            for sp_ in x.spines.values():
                sp_.set_visible(False)
            x.add_patch(plt.Circle((0, 0), st.laa_r * 1e3, fill=False, ls=":", lw=0.5, color="w"))
    for k, x in enumerate(ax.ravel()):
        panel(x, k)
    c0 = fig.colorbar(im0, ax=ax[0, :], shrink=0.8, pad=0.01)
    c0.set_label(r"$|\bar{\boldsymbol{u}}|$ (m s$^{-1}$)" if not plt.rcParams["text.usetex"] else r"$|\bar u|$ (m s$^{-1}$)")
    c1 = fig.colorbar(im1, ax=ax[1, :], shrink=0.8, pad=0.01)
    c1.set_label(r"$\dot\gamma$ (s$^{-1}$)")
    scalebar(ax[1, 0], S, R)
    save(fig, out, f"slices_{tag}")


def scalebar(ax, S, R, length=10):
    x0, y0 = S.min() + 2, R.min() + 2
    ax.plot([x0, x0 + length], [y0, y0], color="k", lw=1.2, solid_capstyle="butt")
    ax.text(x0 + length / 2, y0 + 1.2, f"{length} mm", ha="center", va="bottom", fontsize=7)


def fig_elastic(out, runs, mesh, st, cycles, ph, tag, instant):
    el = {k: r for k, r in runs.items() if r.has("elastic_stress")}
    if not el:
        print("  elastic: no hay corridas elásticas")
        return
    use = [cycles[-1]] if instant else cycles
    S, R, probe0, mask = st.plane_grid(mesh)
    rows = []
    for lab, r in el.items():
        f = phase_mean(r, ["velocity", "elastic_stress", "relaxation_time", "shear_rate", "aggregate_size"], use, ph)
        pr = pv.PolyData(probe0.points).sample(attach(mesh, f))
        u2 = np.sum(np.asarray(pr["velocity"]) ** 2, axis=1)
        ratio = np.asarray(pr["elastic_stress"]) / (RHO * u2 + 1e-12)
        wi = np.asarray(pr["relaxation_time"]) * np.asarray(pr["shear_rate"])
        row = [("ratio", np.ma.array(ratio.reshape(S.shape), mask=~mask)),
               ("wi", np.ma.array(wi.reshape(S.shape), mask=~mask))]
        if "aggregate_size" in f and lab == "TVE":
            nn = np.asarray(pr["aggregate_size"]) / n_st(0.0)
            row.append(("n", np.ma.array(nn.reshape(S.shape), mask=~mask)))
        rows.append((lab, row))
    ncol = max(len(r[1]) for r in rows)
    aspect = (R.max() - R.min()) / (S.max() - S.min())
    fig, ax = plt.subplots(len(rows), ncol, figsize=(JFM_WIDTH, (JFM_WIDTH - 0.3) / ncol * aspect * len(rows) + 0.75),
                           constrained_layout=True, squeeze=False)
    cmaps = {"ratio": ("magma", mcolors.LogNorm(1e-4, 1), r"$|\boldsymbol{\tau}_p|/\rho|\bar{\boldsymbol{u}}|^2$"),
             "wi": ("viridis", mcolors.LogNorm(1e-3, 3), r"$\lambda\dot\gamma$"),
             "n": ("YlOrBr", mcolors.Normalize(0, 1), r"$n/n_{\rm st}(0)$")}
    ims = {}
    for i, (lab, row) in enumerate(rows):
        for j in range(ncol):
            a = ax[i, j]
            a.set_aspect("equal")
            a.set_xticks([])
            a.set_yticks([])
            for sp_ in a.spines.values():
                sp_.set_visible(False)
            if j >= len(row):
                continue
            key, Z = row[j]
            cm, nrm, _ = cmaps[key]
            Zc = np.ma.maximum(Z, nrm.vmin) if isinstance(nrm, mcolors.LogNorm) else Z
            ims[key] = imshow(a, S, R, Zc, cmap=cm, norm=nrm)
            outline(a, mask, S, R)
            a.add_patch(plt.Circle((0, 0), st.laa_r * 1e3, fill=False, ls=":", lw=0.5, color="k"))
            if j == 0:
                a.text(-0.04, 0.5, lab, transform=a.transAxes, rotation=90, ha="right", va="center")
    for key, im in ims.items():
        cb = fig.colorbar(im, ax=ax[:, list(cmaps).index(key)], orientation="horizontal", shrink=0.85, pad=0.01)
        cb.set_label(cmaps[key][2])
    for k, x in enumerate(a for a in ax.ravel() if a.collections):
        panel(x, k)
    save(fig, out, f"elastic_{tag}")


def fig_laa(out, runs, mesh, st, cycles, phases, ph, tag):
    w = 18e-3
    S, R, probe0, mask = st.plane_grid(mesh, h=1.2e-4, window=(-w, w, -w, w))
    ncol = len(runs)
    colw = (JFM_WIDTH - 0.8) / ncol
    fig, ax = plt.subplots(3, ncol, figsize=(JFM_WIDTH, 3 * colw + 0.5),
                           constrained_layout=True, squeeze=False)
    for j, (lab, r) in enumerate(runs.items()):
        names = ["velocity", "shear_rate", "viscosity"] + (["relaxation_time"] if r.has("relaxation_time") else [])
        f = phase_mean(r, names, cycles, ph)
        # Viscosidad efectiva: la del fluido; para VE y TVE, eta_pl + G lambda (la de corte
        # estacionario con la estructura del instante), no la viscosidad numérica del esquema.
        if "relaxation_time" in f:
            f["eta"] = ETAPL + G * f["relaxation_time"]
        elif lab == "N" or "viscosity" not in f:
            f["eta"] = np.full(mesh.n_points, ETAINF)
        else:
            f["eta"] = f["viscosity"]
        f["eta"] = f["eta"] / ETAINF
        f["stasis"] = stasis_fraction(r, cycles, phases)
        pr = pv.PolyData(probe0.points).sample(attach(mesh, f))
        gd = st.project(pr, "shear_rate", S.shape, mask)
        uv = st.project(pr, "velocity", S.shape, mask)
        sf = st.project(pr, "stasis", S.shape, mask)
        et = st.project(pr, "eta", S.shape, mask)
        a, b, e = ax[0, j], ax[1, j], ax[2, j]
        im0 = imshow(a, S, R, np.ma.maximum(gd, 1e-2), cmap="cividis", norm=mcolors.LogNorm(0.1, 100))
        a.streamplot(S, R, uv[0].filled(np.nan), uv[1].filled(np.nan), density=0.9, color="w",
                     linewidth=0.35, arrowsize=0.45)
        im1 = imshow(b, S, R, sf, cmap="Blues", vmin=0, vmax=1)
        im2 = imshow(e, S, R, np.ma.maximum(et, 1), cmap="OrRd",
                     norm=mcolors.LogNorm(1, ETA0 / ETAINF))
        a.set_title(lab)
        for x in (a, b, e):
            outline(x, mask, S, R)
            x.set_aspect("equal")
            x.set_xlim(S.min(), S.max())
            x.set_ylim(R.min(), R.max())
            x.set_xticks([])
            x.set_yticks([])
            for sp_ in x.spines.values():
                sp_.set_visible(False)
    for k, x in enumerate(ax.ravel()):
        panel(x, k)
    fig.colorbar(im0, ax=ax[0, :], shrink=0.8, pad=0.01).set_label(r"$\dot\gamma$ (s$^{-1}$)")
    fig.colorbar(im1, ax=ax[1, :], shrink=0.8, pad=0.01).set_label(r"fraction of $T$ with $\dot\gamma<1$ s$^{-1}$")
    fig.colorbar(im2, ax=ax[2, :], shrink=0.8, pad=0.01).set_label(r"$\eta_{\rm eff}/\eta_\infty$")
    scalebar(ax[2, 0], S, R, length=5)
    save(fig, out, f"laa_{tag}")


def fig_stasis(out, runs, mesh, st, cycles, phases):
    S, R, probe0, mask = st.plane_grid(mesh)
    ncol = len(runs)
    aspect = (R.max() - R.min()) / (S.max() - S.min())
    fig, ax = plt.subplots(1, ncol, figsize=(JFM_WIDTH, (JFM_WIDTH - 0.8) / ncol * aspect + 0.4),
                           constrained_layout=True, squeeze=False)
    for j, (lab, r) in enumerate(runs.items()):
        pr = pv.PolyData(probe0.points).sample(attach(mesh, {"stasis": stasis_fraction(r, cycles, phases)}))
        sf = st.project(pr, "stasis", S.shape, mask)
        a = ax[0, j]
        im = imshow(a, S, R, sf, cmap="Blues", vmin=0, vmax=1)
        outline(a, mask, S, R)
        a.set_title(lab)
        a.set_aspect("equal")
        a.set_xticks([])
        a.set_yticks([])
        for sp_ in a.spines.values():
            sp_.set_visible(False)
        a.add_patch(plt.Circle((0, 0), st.laa_r * 1e3, fill=False, ls=":", lw=0.5, color="k"))
        panel(a, j)
    fig.colorbar(im, ax=ax[0, :], shrink=0.8, pad=0.01).set_label(r"fraction of $T$ with $\dot\gamma<1$ s$^{-1}$")
    save(fig, out, "stasis")


def render_q(mesh, st, u, qlevel, umax, size, camera=None):
    m = attach(mesh, {"velocity": u})
    d = m.compute_derivative(scalars="velocity", qcriterion=True, gradient=False)
    iso = d.contour([qlevel], scalars="qcriterion")
    pl = pv.Plotter(off_screen=True, window_size=size, lighting="three lights")
    pl.set_background("white")
    pl.add_mesh(st.wall, color="#c8c8c8", opacity=0.13, smooth_shading=True, specular=0.2)
    pl.add_mesh(st.wall.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False),
                color="#555555", line_width=1.2)
    if iso.n_points:
        iso.point_data["speed"] = np.linalg.norm(iso.point_data["velocity"], axis=1)
        pl.add_mesh(iso, scalars="speed", cmap="inferno", clim=(0, umax), smooth_shading=True,
                    specular=0.35, specular_power=25, show_scalar_bar=False)
    if camera is None:
        c = st.wall.center
        dist = 2.6 * st.wall.length / 2
        pl.camera.position = tuple(c + dist * st.normal)
        pl.camera.focal_point = tuple(c)
        pl.camera.up = tuple(-st.e1)
        pl.reset_camera()
        pl.camera.zoom(1.15)
    else:
        pl.camera = camera
    try:
        pl.enable_anti_aliasing("ssaa")
    except Exception:
        pass
    img = pl.screenshot(return_img=True)
    cam = pl.camera.copy()
    pl.close()
    return img, cam


def fig_q3d(out, runs, mesh, st, cycles, ph, tag, instant, qplus, u_m):
    use = [cycles[-1]] if instant else cycles
    U = {lab: phase_mean(r, ["velocity"], use, ph)["velocity"] for lab, r in runs.items()}
    scale = (u_m / st.D_m) ** 2
    if qplus is None:
        d = attach(mesh, {"velocity": U[list(U)[0]]}).compute_derivative(scalars="velocity", qcriterion=True,
                                                                          gradient=False)
        q = np.asarray(d["qcriterion"])
        qlevel = float(np.percentile(q[q > 0], 97))
        print(f"  umbral Q = {qlevel:.3g} s^-2 = {qlevel / scale:.3g} (U_m/D_m)^2")
    else:
        qlevel = qplus * scale
    umax = max(float(np.percentile(np.linalg.norm(u, axis=1), 99.5)) for u in U.values())
    imgs, cam = [], None
    for lab in U:
        img, cam = render_q(mesh, st, U[lab], qlevel, umax, (1400, 1400), cam)
        imgs.append(img)
    ncol = len(imgs)
    fig, ax = plt.subplots(1, ncol, figsize=(JFM_WIDTH, (JFM_WIDTH - 0.8) / ncol + 0.45), constrained_layout=True,
                           squeeze=False)
    for j, (lab, img) in enumerate(zip(U, imgs)):
        a = ax[0, j]
        a.imshow(img, interpolation="lanczos")
        a.set_title(lab)
        a.axis("off")
        panel(a, j)
    sm = plt.cm.ScalarMappable(cmap="inferno", norm=mcolors.Normalize(0, umax))
    fig.colorbar(sm, ax=ax[0, :], shrink=0.7, pad=0.01).set_label(r"$|\bar{\boldsymbol{u}}|$ (m s$^{-1}$)")
    fig.text(0.01, 0.01, rf"$Q={qlevel / scale:.2g}\,(U_m/D_m)^2$", fontsize=7, color="0.3")
    save(fig, out, f"q3d_{tag}")


# ------------------------------------------------------------------ tiempos del análisis
def frames_in(run, cycles, period):
    """Cuadros (t, c, ph) de los ciclos dados, en orden temporal, más el del inicio."""
    out = []
    for (c, ph) in run.frames:
        if c in cycles or (c == cycles[0] - 1 and abs(ph - period) < 1e-6):
            out.append(((c - 1) * period + ph, c, ph))
    return sorted(out)


# ------------------------------------------------------------------ esfuerzo en la pared
def wall_setup(mesh, st):
    """Triángulos de pared, su tetraedro vecino y los gradientes de base de ese tetraedro."""
    T, X = tets(mesh), np.asarray(mesh.points)
    F = st.F[st.F[:, 3] == 4, :3] - 1
    nv = len(X)
    if nv >= 2 ** 21:
        sys.exit("Malla demasiado grande para la codificación de caras.")
    loc = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    keys = np.concatenate([np.sort(T[:, l], axis=1) for l in loc])
    k = (keys[:, 0] << 42) | (keys[:, 1] << 21) | keys[:, 2]
    owner = np.tile(np.arange(len(T)), 4)
    order = np.argsort(k)
    k, owner = k[order], owner[order]
    fs = np.sort(F, axis=1)
    fk = (fs[:, 0] << 42) | (fs[:, 1] << 21) | fs[:, 2]
    pos = np.searchsorted(k, fk)
    ok = (pos < len(k)) & (k[np.minimum(pos, len(k) - 1)] == fk)
    F, pos = F[ok], pos[ok]
    tet = owner[pos]
    Tt = T[tet]
    E = X[Tt[:, 1:]] - X[Tt[:, [0]]]
    Einv = np.linalg.inv(E)
    Gb = np.concatenate([-Einv.sum(axis=2, keepdims=True), Einv], axis=2)  # (nt, 3, 4)
    P = X[F]
    nrm = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    area = 0.5 * np.linalg.norm(nrm, axis=1)
    nrm /= (2 * area)[:, None]
    cen = P.mean(axis=1)
    inl = np.linalg.norm(cen - st.laa_c, axis=1) < st.laa_r
    return dict(F=F, tet=T[tet], G=Gb, n=nrm, area=area, laa=inl)


def wall_shear(W, u, eta_fn, cmat=None):
    """Vector del esfuerzo de corte en la pared, por triángulo."""
    J = np.einsum("tin,tnk->tki", W["G"], u[W["tet"]])          # J[k, i] = d u_k / d x_i
    D = 0.5 * (J + J.transpose(0, 2, 1))
    gd = np.sqrt(2 * np.einsum("tij,tij->t", D, D))
    t = 2 * eta_fn(gd)[:, None] * np.einsum("tij,tj->ti", D, W["n"])
    if cmat is not None:                                         # esfuerzo elástico, media de los vértices
        tau = G * (cmat[W["F"]].mean(axis=1) - np.eye(3))
        t = t + np.einsum("tij,tj->ti", tau, W["n"])
    return t - np.einsum("ti,ti->t", t, W["n"])[:, None] * W["n"]


def fig_wall(out, runs, mesh, st, cycles, period):
    W = wall_setup(mesh, st)
    res = {}
    for lab, r in runs.items():
        elastic = r.has("psi_diagonal")
        if lab == "N":
            eta = lambda g: np.full_like(g, ETAINF)
        elif elastic:
            eta = lambda g: np.full_like(g, ETAPL)
        else:
            eta = eta_cy
        frs = [f for f in frames_in(r, cycles, period) if f[1] in cycles]
        sv, sm = 0, 0
        for (_, c, ph) in frs:
            tw = wall_shear(W, r.field("velocity", c, ph), eta, conformation(r, c, ph) if elastic else None)
            sv = sv + tw
            sm = sm + np.linalg.norm(tw, axis=1)
        tawss = sm / len(frs)
        osi = 0.5 * (1 - np.linalg.norm(sv / len(frs), axis=1) / np.maximum(tawss, 1e-12))
        ecap = osi / np.maximum(tawss, 1e-6)
        res[lab] = dict(tawss=tawss, osi=osi, ecap=ecap)
    A = W["area"]
    rows = []
    for lab, d in res.items():
        for zone, m in (("LAA", W["laa"]), ("wall", np.ones(len(A), bool))):
            a = A[m]
            rows.append([lab, zone] + [float(np.sum(a * d[k][m]) / a.sum()) for k in ("tawss", "osi", "ecap")]
                        + [float(a[d["tawss"][m] < 0.1].sum() / a.sum()), float(a[d["ecap"][m] > 1.0].sum() / a.sum())])
    with open(os.path.join(out, "wall_metrics.csv"), "w") as f:
        f.write("fluid,zone,TAWSS_Pa,OSI,ECAP_1/Pa,area_TAWSS<0.1Pa,area_ECAP>1/Pa\n")
        for r_ in rows:
            f.write(",".join([r_[0], r_[1]] + [f"{x:.5g}" for x in r_[2:]]) + "\n")
    print("  wall_metrics.csv")
    for r_ in rows:
        if r_[1] == "LAA":
            print(f"    {r_[0]:>4s} orejuela: TAWSS {r_[2]:.3g} Pa, OSI {r_[3]:.3f}, ECAP {r_[4]:.3g} 1/Pa, "
                  f"área con TAWSS<0.1 Pa {100 * r_[5]:.0f}%")
    # Mapas en la pared de la orejuela.
    near = np.linalg.norm(st.V[W["F"]].mean(axis=1) - st.laa_c, axis=1) < 1.7 * st.laa_r
    tri = W["F"][near]
    faces = np.hstack([np.full((len(tri), 1), 3), tri]).ravel()
    surf = pv.PolyData(np.asarray(mesh.points), faces).clean()
    view = st.laa_c - np.asarray(st.wall.center)
    view /= np.linalg.norm(view)
    imgs = []
    tmax = max(np.percentile(d["tawss"][near], 98) for d in res.values())
    for lab, d in res.items():
        row = []
        for key, cmap, clim, log in (("tawss", "viridis", (1e-3, tmax), True), ("ecap", "magma", (1e-2, 1e2), True)):
            sf = pv.PolyData(np.asarray(mesh.points), faces)
            sf.cell_data[key] = np.maximum(d[key][near], clim[0])
            pl = pv.Plotter(off_screen=True, window_size=(900, 900))
            pl.set_background("white")
            pl.add_mesh(sf, scalars=key, cmap=cmap, clim=clim, log_scale=log, smooth_shading=True,
                        ambient=0.45, diffuse=0.6, specular=0.05, show_scalar_bar=False)
            pl.camera.position = tuple(st.laa_c + 4 * st.laa_r * view)
            pl.camera.focal_point = tuple(st.laa_c)
            pl.camera.up = tuple(-st.e1)
            pl.reset_camera()
            pl.camera.zoom(1.3)
            try:
                pl.enable_anti_aliasing("ssaa")
            except Exception:
                pass
            row.append(pl.screenshot(return_img=True))
            pl.close()
        imgs.append((lab, row))
    ncol = len(imgs)
    colw = (JFM_WIDTH - 0.8) / ncol
    fig, ax = plt.subplots(2, ncol, figsize=(JFM_WIDTH, 2 * colw + 0.45), constrained_layout=True, squeeze=False)
    for j, (lab, row) in enumerate(imgs):
        for i in range(2):
            ax[i, j].imshow(row[i], interpolation="lanczos")
            ax[i, j].axis("off")
        ax[0, j].set_title(lab)
    for k, x in enumerate(ax.ravel()):
        panel(x, k)
    fig.colorbar(plt.cm.ScalarMappable(cmap="viridis", norm=mcolors.LogNorm(1e-3, tmax)), ax=ax[0, :],
                 shrink=0.8, pad=0.01).set_label("TAWSS (Pa)")
    fig.colorbar(plt.cm.ScalarMappable(cmap="magma", norm=mcolors.LogNorm(1e-2, 1e2)), ax=ax[1, :],
                 shrink=0.8, pad=0.01).set_label(r"ECAP (Pa$^{-1}$)")
    save(fig, out, "wall_laa")


# ------------------------------------------------------------------ residencia en la orejuela
class Locator:
    """Ubica puntos en los tetraedros (árbol de centroides + coordenadas baricéntricas)."""

    def __init__(self, mesh):
        from scipy.spatial import cKDTree
        self.T, self.X = tets(mesh), np.asarray(mesh.points)
        P = self.X[self.T]
        self.x0 = P[:, 0]
        self.inv = np.linalg.inv((P[:, 1:] - P[:, [0]]).transpose(0, 2, 1))
        self.tree = cKDTree(P.mean(axis=1))

    def find(self, p, k=16):
        _, cand = self.tree.query(p, k=k)
        lam = np.einsum("pkij,pkj->pki", self.inv[cand], p[:, None, :] - self.x0[cand])
        bary = np.concatenate([1 - lam.sum(axis=2, keepdims=True), lam], axis=2)
        ok = (bary > -1e-9).all(axis=2)
        hit = ok.any(axis=1)
        first = ok.argmax(axis=1)
        idx = np.arange(len(p))
        cell = np.where(hit, cand[idx, first], -1)
        return cell, bary[idx, first]

    def interp(self, field, cell, bary):
        return np.einsum("pa,pac->pc", bary, field[self.T[cell]])


def fig_residence(out, runs, mesh, st, cycles, period, spacing, sub):
    try:
        loc = Locator(mesh)
    except ImportError:
        print("  residence: falta scipy (pip install scipy)")
        return
    # Semillas: malla regular dentro de la esfera de la orejuela y del atrio.
    h = spacing * 1e-3
    g = np.arange(-st.laa_r, st.laa_r + h / 2, h)
    P = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    P = P[np.linalg.norm(P, axis=1) < st.laa_r] + st.laa_c
    cell, _ = loc.find(P)
    seeds = P[cell >= 0]
    print(f"  residence: {len(seeds)} partículas en la orejuela (cada {spacing} mm), "
          f"{sub} subpasos por cuadro")
    op = np.array([st.openings[l] for l in st.openings if l != 4])
    curves = {}
    for lab, r in runs.items():
        frs = frames_in(r, cycles, period)
        t0 = frs[0][0]
        x = seeds.copy()
        alive = np.ones(len(x), bool)                  # sigue en el atrio
        ts, inside, inatrium = [0.0], [1.0], [1.0]
        U0 = r.field("velocity", frs[0][1], frs[0][2])
        for j in range(len(frs) - 1):
            ta, tb = frs[j][0], frs[j + 1][0]
            U1 = r.field("velocity", frs[j + 1][1], frs[j + 1][2])
            dt = (tb - ta) / sub
            for s_ in range(sub):
                th0, th1 = s_ / sub, (s_ + 0.5) / sub

                def vel(p, th, m):
                    c_, b_ = loc.find(p[m])
                    v = np.zeros((m.sum(), 3))
                    ok_ = c_ >= 0
                    v[ok_] = (1 - th) * loc.interp(U0, c_[ok_], b_[ok_]) + th * loc.interp(U1, c_[ok_], b_[ok_])
                    return v, ok_

                m = alive.copy()
                v1, _ = vel(x, th0, m)
                xm = x.copy()
                xm[m] = x[m] + 0.5 * dt * v1
                v2, ok2 = vel(xm, th1, m)
                xn = x.copy()
                xn[m] = x[m] + dt * v2
                c_, _ = loc.find(xn[m])
                lost = c_ < 0
                if lost.any():
                    ids = np.where(m)[0][lost]
                    # Fuera del dominio: salió si está cerca de una abertura; si no, queda en su lugar (pared).
                    d_op = np.min(np.linalg.norm(xn[ids][:, None, :] - op[None], axis=2), axis=1)
                    gone = d_op < 0.6 * st.D_m
                    alive[ids[gone]] = False
                    xn[ids[~gone]] = x[ids[~gone]]
                x = xn
            U0 = U1
            ts.append(tb - t0)
            inside.append(float(np.mean(alive & (np.linalg.norm(x - st.laa_c, axis=1) < st.laa_r))))
            inatrium.append(float(np.mean(alive)))
        curves[lab] = (np.array(ts), np.array(inside), np.array(inatrium))
        half = np.interp(0.5, curves[lab][1][::-1], curves[lab][0][::-1]) if curves[lab][1].min() < 0.5 else np.nan
        print(f"    {lab:>4s}: en la orejuela tras 1, 2, 4 ciclos (los que alcancen): "
              + ", ".join(f"{100 * np.interp(k * period, ts, inside):.0f}%" for k in (1, 2, 4) if k * period <= ts[-1] + 1e-9)
              + (f"; vida media {half / period:.2f} T" if np.isfinite(half) else "; vida media > integración"))
    with open(os.path.join(out, "residence.csv"), "w") as f:
        f.write("fluid,t_over_T,fraction_in_LAA,fraction_in_atrium\n")
        for lab, (ts, ins, ina) in curves.items():
            for a_, b_, c_ in zip(ts, ins, ina):
                f.write(f"{lab},{a_ / period:.4f},{b_:.5f},{c_:.5f}\n")
    fig, ax = plt.subplots(1, 2, figsize=(JFM_WIDTH, 2.1), constrained_layout=True)
    for lab, (ts, ins, ina) in curves.items():
        ax[0].plot(ts / period, ins, color=COLORS.get(lab), label=lab)
        ax[1].plot(ts / period, ina, color=COLORS.get(lab), label=lab)
    ax[0].set_ylabel(r"fraction in $\Omega_a$")
    ax[1].set_ylabel("fraction in the atrium")
    for a_ in ax:
        a_.set_xlabel(r"$t/T$ since release")
        a_.set_ylim(0, 1.02)
        a_.set_xlim(0, None)
    ax[0].axhline(0.5, color="0.6", lw=0.5, ls=":")
    ax[0].legend(loc="lower left")
    for k, a_ in enumerate(ax):
        panel(a_, k)
    save(fig, out, "residence")


# ------------------------------------------------------------------ energía
def fig_energy(out, runs, cycles):
    keys = [("kineticEnergy", r"$E_k$", 1e3), ("elasticEnergy", r"$E_p$", 1e3),
            ("viscousDissipation", r"$\mathcal{D}_v$", 1e3), ("polymerPower", r"$\mathcal{P}_p$", 1e3)]
    tab = {}
    for lab, r in runs.items():
        S = r.S
        per = []
        for c in cycles:
            m = S["cycle"].astype(int) == c
            per.append([S[k][m].mean() for k, _, _ in keys])
        per = np.array(per)
        tab[lab] = per
    with open(os.path.join(out, "energy.csv"), "w") as f:
        f.write("fluid," + ",".join(f"{k}_mean,{k}_std" for k, _, _ in keys)
                + ",polymer_share_of_power,elastic_over_kinetic_energy\n")
        for lab, per in tab.items():
            m, sd = per.mean(axis=0), per.std(axis=0, ddof=1) if len(per) > 1 else np.zeros(len(keys))
            frac = per[:, 3] / (per[:, 2] + per[:, 3])
            ratio = per[:, 1] / per[:, 0]
            f.write(lab + "," + ",".join(f"{a:.6g},{b:.3g}" for a, b in zip(m, sd))
                    + f",{frac.mean():.4g},{ratio.mean():.4g}\n")
            # La potencia de los agregados es la parte viscosa que cargan ellos (eta_CY - eta_pl);
            # la memoria se mide por la energía elástica almacenada frente a la cinética.
            print(f"    {lab:>4s}: E_k {1e3 * m[0]:.3f} mJ, E_p {1e3 * m[1]:.3g} mJ (E_p/E_k = {ratio.mean():.2e}), "
                  f"D_v + P_p = {1e3 * (m[2] + m[3]):.3f} mW, de ella {100 * frac.mean():.0f}% por los agregados")
    print("  energy.csv")
    fig, ax = plt.subplots(1, 2, figsize=(JFM_WIDTH, 2.1), constrained_layout=True)
    labs = list(tab)
    x = np.arange(len(labs))
    ekw = dict(lw=0.6, capsize=1.5)
    # (a) Energía cinética y elástica, en escala logarítmica: la elástica es órdenes menor.
    for i, (name, col) in enumerate(((r"$E_k$", "0.25"), (r"$E_p$", COLORS["TVE"]))):
        v = np.array([1e3 * tab[l][:, i].mean() for l in labs])
        e = np.array([1e3 * tab[l][:, i].std() for l in labs])
        keep = v > 0
        ax[0].bar(x[keep] + (i - 0.5) * 0.38, v[keep], 0.36, yerr=e[keep], color=col, label=name, error_kw=ekw)
    ax[0].set_yscale("log")
    ax[0].set_ylabel("energy (mJ)")
    # (b) Potencia disipada: viscosa del solvente (o del fluido) y la de los agregados, apiladas.
    dv = np.array([1e3 * tab[l][:, 2].mean() for l in labs])
    pp = np.array([1e3 * tab[l][:, 3].mean() for l in labs])
    tot_e = np.array([1e3 * (tab[l][:, 2] + tab[l][:, 3]).std() for l in labs])
    ax[1].bar(x, dv, 0.55, color="0.25", label=r"$\mathcal{D}_v$ (solvent)")
    ax[1].bar(x, pp, 0.55, bottom=dv, color=COLORS["VE"], label=r"$\mathcal{P}_p$ (aggregates)",
              yerr=tot_e, error_kw=ekw)
    ax[1].set_ylabel("power (mW)")
    for k, a_ in enumerate(ax):
        a_.set_xticks(x)
        a_.set_xticklabels(labs)
        a_.tick_params(axis="x", which="minor", bottom=False, top=False)
        a_.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, handlelength=1.2)
        panel(a_, k, x=-0.2, y=1.0)
    save(fig, out, "energy")


# ------------------------------------------------------------------ deriva de la estructura
def fig_drift(out, runs, mesh, st, cycles, period):
    el = {k: r for k, r in runs.items() if r.has("aggregate_size")}
    if not el:
        print("  drift: no hay corridas con estructura")
        return
    lump = lumped_volume(mesh)
    inside = np.linalg.norm(np.asarray(mesh.points) - st.laa_c, axis=1) < st.laa_r
    w = lump * inside
    ns0 = float(n_st(0.0))
    fig, ax = plt.subplots(figsize=(JFM_WIDTH * 0.62, 2.0), constrained_layout=True)
    for lab, r in el.items():
        pts = sorted(((c - 1) * period + ph, c, ph) for (c, ph) in r.frames)
        t = np.array([p[0] for p in pts])
        nm = np.array([np.sum(w * r.field("aggregate_size", c, ph)) / w.sum() for (_, c, ph) in pts]) / ns0
        ax.plot(t / period, nm, color=COLORS.get(lab), label=lab, lw=0.9)
        sel = np.isin([p[1] for p in pts], cycles)
        per = [nm[np.array([p[1] for p in pts]) == c].mean() for c in cycles]
        slope = np.polyfit(np.arange(len(per)), per, 1)[0] if len(per) > 1 else np.nan
        print(f"    {lab:>4s}: <n>/n_st(0) en la orejuela {np.mean(per):.3f}, deriva {100 * slope / np.mean(per):+.2f}% por ciclo"
              f" en los ciclos analizados")
    ax.axvspan(cycles[0] - 1, cycles[-1], color="0.9", lw=0, zorder=0)
    ax.set_xlabel(r"$t/T$")
    ax.set_ylabel(r"$\langle n\rangle_{\Omega_a}/n_{\rm st}(0)$")
    ax.legend(loc="best")
    save(fig, out, "drift_laa")


# ------------------------------------------------------------------ principal
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.expanduser("~/runs/fa"))
    ap.add_argument("--runs", default="N:prod_N,CY:prod_CY,VE:prod_VE,TVE:prod_TVE")
    ap.add_argument("--mesh", default="Atrium.mesh", help="malla medit con etiquetas de borde")
    ap.add_argument("--laa", default="44,188.4,-197.3,13", help="esfera de la orejuela: x,y,z,r en mm")
    ap.add_argument("--what", default="all",
                    help="material,waveforms,energy,slices,elastic,laa,stasis,q3d,wall,residence,drift o all")
    ap.add_argument("--seed-spacing", type=float, default=0.7, help="separación de las partículas en mm")
    ap.add_argument("--substeps", type=int, default=10, help="subpasos de integración por cuadro")
    ap.add_argument("--phases", default="auto",
                    help="fases en s dentro del ciclo; auto = pico mitral y válvula cerrada")
    ap.add_argument("--skip", type=int, default=2)
    ap.add_argument("--period", type=float, default=0.85)
    ap.add_argument("--qplus", type=float, default=None, help="umbral Q en unidades de (U_m/D_m)^2")
    ap.add_argument("--instant", action="store_true", help="último ciclo en vez de media de fase")
    ap.add_argument("--usetex", action="store_true", help="texto con LaTeX (más lento, idéntico a JFM)")
    ap.add_argument("--out", default="figs")
    a = ap.parse_args()

    style(a.usetex)
    what = set(("material,waveforms,energy,slices,elastic,laa,stasis,q3d,wall,residence,drift"
                if a.what == "all" else a.what).split(","))
    os.chdir(os.path.expanduser(a.dir))
    os.makedirs(a.out, exist_ok=True)
    if "material" in what:
        fig_material(a.out, a.period)
        what.discard("material")
    if not what:
        return
    labels = dict(x.split(":", 1) for x in a.runs.split(","))
    runs = {k: Run(k, v, a.period) for k, v in labels.items() if os.path.exists(parts(v)[0] + ".csv")}
    missing = set(labels) - set(runs)
    if missing:
        print("  sin datos (se omiten): " + ", ".join(sorted(missing)))
    if not runs:
        sys.exit("No hay corridas.")
    cycles, phases = common_cycles(runs, a.period, a.skip)
    print(f"Ciclos {cycles[0]}-{cycles[-1]} (n = {len(cycles)}), {len(phases)} fases por ciclo")

    if a.phases == "auto":
        r0 = runs[list(runs)[0]]
        S = r0.S
        sel = np.isin(S["cycle"].astype(int), cycles)
        tau = S["t"][sel] - (S["cycle"][sel] - 1) * a.period
        q = S["qOut"][sel]
        t_peak = tau[np.argmax(q)]
        closed = tau[np.abs(q) < 0.02 * q.max()]
        t_closed = np.median(closed) if len(closed) else tau[np.argmin(q)]
        wanted = [t_peak, t_closed]
    else:
        wanted = [float(x) for x in a.phases.split(",")]
    chosen = []
    for w in wanted:
        ph = min(phases, key=lambda p: abs(p - w))
        if ph not in chosen:
            chosen.append(ph)
    print("Fases: " + ", ".join(f"{p:.3f} s" for p in chosen))

    if "waveforms" in what:
        fig_waveforms(a.out, runs, cycles, a.period, chosen)
    if "energy" in what:
        fig_energy(a.out, runs, cycles)
    need3d = what & {"slices", "elastic", "laa", "stasis", "q3d", "wall", "residence", "drift"}
    if not need3d:
        return
    mesh = read_mesh(runs[list(runs)[0]])
    st = Setting(a.mesh, [float(x) for x in a.laa.split(",")])
    q0 = runs[list(runs)[0]].S["qOut"]
    u_m = q0.max() / st.mitral_area
    print(f"Malla: {mesh.n_points} vértices, {mesh.n_cells} celdas; D_m = {1e3 * st.D_m:.1f} mm, "
          f"U_m = {u_m:.3f} m/s")
    for ph in chosen:
        tag = f"{ph / a.period:.2f}T"
        if "slices" in what:
            fig_slices(a.out, runs, mesh, st, cycles, ph, tag, a.instant)
        if "elastic" in what:
            fig_elastic(a.out, runs, mesh, st, cycles, ph, tag, a.instant)
        if "laa" in what:
            fig_laa(a.out, runs, mesh, st, cycles, phases, ph, tag)
        if "q3d" in what:
            fig_q3d(a.out, runs, mesh, st, cycles, ph, tag, a.instant, a.qplus, u_m)
    if "stasis" in what:
        fig_stasis(a.out, runs, mesh, st, cycles, phases)
    if "drift" in what:
        fig_drift(a.out, runs, mesh, st, cycles, a.period)
    if "wall" in what:
        fig_wall(a.out, runs, mesh, st, cycles, a.period)
    if "residence" in what:
        fig_residence(a.out, runs, mesh, st, cycles, a.period, a.seed_spacing, a.substeps)


if __name__ == "__main__":
    main()

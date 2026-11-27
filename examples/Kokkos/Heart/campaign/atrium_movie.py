#!/usr/bin/env python3
"""
atrium_movie.py -- video de un ciclo de TVE en el atrio.

    python3 atrium_movie.py [--dir ~/runs/fa] [--run movie_TVE] [--out tve_cycle.mp4]
                            [--fps 30] [--loops 2] [--orbit 50] [--jobs 4]

Cada cuadro (1920x1080, fondo oscuro):
  izquierda   estructuras de vórtice (criterio Q, umbral fijo en todo el video)
              coloreadas por |u|, dentro de la pared translúcida; la cámara gira
              lentamente (--orbit grados en todo el video);
  derecha     corte por la orejuela y el plano mitral: |u| con líneas de
              corriente, y el tamaño de agregado n/n_st(0), la memoria de la
              sangre; abajo, el caudal mitral del ciclo con el instante marcado.

Conviene una corrida con salida densa (output_every=5, 170 cuadros por ciclo):
ver la ayuda de make_movie_data.sh. Con la salida de producción (17 cuadros por
ciclo) el video sale entrecortado; --interp k agrega k-1 cuadros interpolados
linealmente entre cada par (útil para probar, no para publicar).

Requisitos: pip install pyvista matplotlib h5py, y ffmpeg.
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import warnings
from multiprocessing import Pool

import numpy as np

warnings.filterwarnings("ignore")
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import pyvista as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atrium_post import Setting, n_st, parts  # noqa: E402

BG = "#0d0f14"
FG = "#e8e8ec"
ACCENT = "#ff8a3d"

# ------------------------------------------------------------------ datos
G = {}          # estado por proceso


def frames_of(prefix):
    out = []
    for q in parts(prefix):
        times = [float(t) for t in re.findall(r'Time Value="([^"]+)"', open(q + ".xdmf").read())]
        for k, t in enumerate(times):
            if os.path.exists(f"{q}.velocity.{k:06d}.h5"):
                out.append((t, q, k))
    out.sort()
    keep = {}
    for t, q, k in out:                       # la parte posterior manda
        keep[round(t, 9)] = (t, q, k)
    return [keep[t] for t in sorted(keep)]


def load(q, k, name):
    f = f"{q}.{name}.{k:06d}.h5"
    return h5py.File(f, "r")["GridFunction/Values/Data"][:] if os.path.exists(f) else None


def init(args):
    os.chdir(os.path.expanduser(args.dir))
    m = h5py.File(parts(args.run)[0] + ".mesh.h5", "r")
    X = m["Mesh/Geometry/Vertices"][:]
    unit = 1e-3 if np.ptp(X, axis=0).max() > 1.0 else 1.0
    T = m["Mesh/XDMF/Topology"][:].reshape(-1, 5)[:, 1:].astype(np.int64)
    cells = np.hstack([np.full((len(T), 1), 4), T]).ravel()
    G["mesh"] = pv.UnstructuredGrid(cells, np.full(len(T), pv.CellType.TETRA, np.uint8), X * unit)
    G["st"] = Setting(args.mesh, [float(x) for x in args.laa.split(",")])
    st = G["st"]
    S, R, probe, mask = st.plane_grid(G["mesh"], h=args.h * 1e-3)
    G["plane"] = (S, R, probe.points.copy(), mask)
    G["args"] = args


def fields_at(frame, frac=None, other=None):
    """Campos del cuadro, o interpolados entre dos cuadros con peso frac."""
    t, q, k = frame
    names = ("velocity", "aggregate_size")
    f = {n: load(q, k, n) for n in names}
    if other is not None and frac:
        t2, q2, k2 = other
        for n in names:
            b = load(q2, k2, n)
            if f[n] is not None and b is not None:
                f[n] = (1 - frac) * f[n] + frac * b
        t = (1 - frac) * t + frac * t2
    return t, f


# ------------------------------------------------------------------ render
def render_3d(u, angle, qlevel, umax, size):
    st, mesh = G["st"], G["mesh"]
    m = mesh.copy(deep=False)
    m.point_data["velocity"] = u
    d = m.compute_derivative(scalars="velocity", qcriterion=True, gradient=False)
    iso = d.contour([qlevel], scalars="qcriterion")
    pl = pv.Plotter(off_screen=True, window_size=size, lighting="three lights")
    pl.set_background(BG)
    pl.add_mesh(st.wall, color="#9aa3b5", opacity=0.10, smooth_shading=True, specular=0.4)
    pl.add_mesh(st.wall.extract_feature_edges(boundary_edges=True, feature_edges=False,
                                              manifold_edges=False), color="#c9cfdb", line_width=1.5)
    if iso.n_points:
        iso.point_data["speed"] = np.linalg.norm(iso.point_data["velocity"], axis=1)
        pl.add_mesh(iso, scalars="speed", cmap="inferno", clim=(0, umax), smooth_shading=True,
                    specular=0.5, specular_power=30, ambient=0.15, show_scalar_bar=False)
    c = np.array(st.wall.center)
    dist = 2.5 * st.wall.length / 2
    a = np.deg2rad(angle)
    view = np.cos(a) * st.normal + np.sin(a) * st.e2      # giro alrededor de e1
    pl.camera.position = tuple(c + dist * view)
    pl.camera.focal_point = tuple(c)
    pl.camera.up = tuple(-st.e1)
    pl.reset_camera()
    pl.camera.zoom(1.15)
    try:
        pl.enable_anti_aliasing("ssaa")
    except Exception:
        pass
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


def frame_job(job):
    idx, frame, frac, other, angle = job
    a = G["args"]
    out = os.path.join(a.tmp, f"frame_{idx:05d}.png")
    if os.path.exists(out):
        return out
    t, f = fields_at(frame, frac, other)
    st = G["st"]
    S, R, pts, mask = G["plane"]
    img3d = render_3d(f["velocity"], angle, a.qlevel, a.umax, (1080, 1080))

    m = G["mesh"].copy(deep=False)
    m.point_data["velocity"] = f["velocity"]
    if f["aggregate_size"] is not None:
        m.point_data["n"] = f["aggregate_size"] / n_st(0.0)
    pr = pv.PolyData(pts).sample(m)
    uv = [np.ma.array((np.asarray(pr["velocity"]) @ e).reshape(S.shape), mask=~mask) for e in (st.e1, st.e2)]
    sp = np.ma.array(np.linalg.norm(pr["velocity"], axis=1).reshape(S.shape), mask=~mask)

    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 13,
                         "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG,
                         "axes.edgecolor": FG})
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor=BG)
    ax3 = fig.add_axes([0.0, 0.0, 0.5625, 1.0])
    ax3.imshow(img3d)
    ax3.axis("off")
    # Paneles del corte.
    hgt = 0.34
    w = min(0.34, hgt * 10.8 / 19.2 * (S.max() - S.min()) / (R.max() - R.min()))
    x0 = 0.585
    axu = fig.add_axes([x0, 0.585, w, hgt])
    axn = fig.add_axes([x0, 0.585 - hgt - 0.045, w, hgt])
    for ax in (axu, axn):
        ax.set_facecolor(BG)
        ax.set_aspect("equal")
        ax.axis("off")
    axu.pcolormesh(S, R, sp, shading="gouraud", cmap="inferno", vmin=0, vmax=a.umax, rasterized=True)
    axu.streamplot(S, R, uv[0].filled(np.nan), uv[1].filled(np.nan), density=1.4, color=(1, 1, 1, 0.75),
                   linewidth=0.6, arrowsize=0.6)
    axu.contour(S, R, mask.astype(float), [0.5], colors=[FG], linewidths=0.8)
    axu.text(0.0, 1.02, r"$|\boldsymbol{u}|$ and streamlines", transform=axu.transAxes, fontsize=14)
    if "n" in pr.point_data:
        nn = np.ma.array(np.asarray(pr["n"]).reshape(S.shape), mask=~mask)
        axn.pcolormesh(S, R, nn, shading="gouraud", cmap="YlOrBr_r" if a.dark_n else "YlOrBr",
                       vmin=0, vmax=1, rasterized=True)
        axn.text(0.0, 1.02, r"aggregate size $n/n_{\rm st}(0)$: the memory of blood",
                 transform=axn.transAxes, fontsize=14)
    axn.contour(S, R, mask.astype(float), [0.5], colors=[FG], linewidths=0.8)
    for ax in (axu, axn):
        ax.add_patch(plt.Circle((0, 0), st.laa_r * 1e3, fill=False, ls=":", lw=0.9, color=FG))
    axu.text(0, -st.laa_r * 1e3 - 2, "LAA", ha="center", va="top", fontsize=11, color=FG)
    # Barras de color.
    cbu = fig.add_axes([x0 + w + 0.012, 0.585 + 0.04, 0.007, hgt - 0.08])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap="inferno", norm=mcolors.Normalize(0, a.umax)), cax=cbu)
    cb.set_label(r"m s$^{-1}$", color=FG)
    cb.outline.set_edgecolor(FG)
    if "n" in pr.point_data:
        cbn = fig.add_axes([x0 + w + 0.012, 0.585 - hgt - 0.045 + 0.04, 0.007, hgt - 0.08])
        cb2 = fig.colorbar(plt.cm.ScalarMappable(cmap="YlOrBr_r" if a.dark_n else "YlOrBr",
                                                 norm=mcolors.Normalize(0, 1)), cax=cbn)
        cb2.outline.set_edgecolor(FG)
    # Caudal mitral con el instante.
    axq = fig.add_axes([0.61, 0.075, 0.33, 0.12], facecolor=BG)
    tc, qc = a.wave
    axq.plot(tc, qc, color=FG, lw=1.3)
    tau = ((t - a.t0) % a.period) / a.period
    axq.fill_between(tc, 0, qc, where=tc <= tau, color=ACCENT, alpha=0.18, lw=0, interpolate=True)
    axq.plot([tau], [np.interp(tau, tc, qc)], "o", color=ACCENT, ms=8)
    axq.set_xlim(0, 1)
    axq.set_xlabel(r"$t/T$")
    axq.set_ylabel(r"$Q_m$ (mL s$^{-1}$)")
    for sp_ in ("top", "right"):
        axq.spines[sp_].set_visible(False)
    # Títulos.
    fig.text(0.02, 0.955, "Left atrium, thixotropic-viscoelastic blood (Owens)", fontsize=20, color=FG)
    fig.text(0.02, 0.92, a.subtitle, fontsize=13, color="#a8adbb")
    fig.text(0.02, 0.04, rf"$t/T = {tau:.2f}$", fontsize=18, color=FG)
    fig.text(0.02, 0.012, rf"vortices: $Q = {a.qplus:.2g}\,(U_m/D_m)^2$, coloured by $|\boldsymbol{{u}}|$",
             fontsize=11, color="#a8adbb")
    fig.savefig(out, facecolor=BG, dpi=100)
    plt.close(fig)
    return out


# ------------------------------------------------------------------ principal
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.expanduser("~/runs/fa"))
    ap.add_argument("--run", default="movie_TVE")
    ap.add_argument("--mesh", default="Atrium.mesh")
    ap.add_argument("--laa", default="44,188.4,-197.3,13")
    ap.add_argument("--out", default="tve_cycle.mp4")
    ap.add_argument("--period", type=float, default=0.85)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--loops", type=int, default=2, help="veces que se repite el ciclo")
    ap.add_argument("--orbit", type=float, default=50.0, help="grados de giro de la cámara en todo el video")
    ap.add_argument("--interp", type=int, default=1, help="cuadros por intervalo de salida (1 = sin interpolar)")
    ap.add_argument("--qplus", type=float, default=None, help="umbral Q en (U_m/D_m)^2; por defecto automático")
    ap.add_argument("--h", type=float, default=0.25, help="resolución del corte en mm")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--max-frames", type=int, default=0, help="solo para probar")
    ap.add_argument("--crf", type=int, default=18)
    ap.add_argument("--preset", default="veryfast", help="slow para la versión final, más fina")
    ap.add_argument("--subtitle", default="reduced venous load, 2.0 M elements, one cardiac cycle")
    ap.add_argument("--dark-n", action="store_true", help="n con la escala invertida (agregado = oscuro)")
    ap.add_argument("--keep", action="store_true", help="no borrar los cuadros PNG")
    a = ap.parse_args()
    if not shutil.which("ffmpeg"):
        sys.exit("Falta ffmpeg: sudo apt install ffmpeg")
    a.out = os.path.abspath(a.out)
    os.chdir(os.path.expanduser(a.dir))

    frames = frames_of(a.run)
    frames = [f for f in frames if f[0] > 0]
    if len(frames) < 2:
        sys.exit(f"{a.run}: no hay cuadros.")
    t0 = frames[0][0]
    one = [f for f in frames if f[0] < t0 + a.period - 1e-9]       # un ciclo
    a.t0 = (np.ceil(t0 / a.period - 1e-6) - 1) * a.period
    dtout = np.median(np.diff([f[0] for f in one]))
    print(f"{a.run}: {len(one)} cuadros en el ciclo (salida cada {1e3 * dtout:.0f} ms)")
    if len(one) < 60 and a.interp == 1:
        print("  aviso: pocos cuadros; el video saldrá entrecortado (ver --interp o make_movie_data.sh)")

    # Caudal mitral del ciclo (del csv).
    byt = {}
    for q in parts(a.run):
        for r in csv.DictReader(open(q + ".csv")):
            byt[round(float(r["t"]), 9)] = r
    ts = np.array(sorted(byt))
    sel = (ts >= one[0][0] - 1e-9) & (ts <= one[-1][0] + 1e-9)
    tcyc = ((ts[sel] - a.t0) % a.period) / a.period
    order = np.argsort(tcyc)
    a.wave = (tcyc[order], np.array([float(byt[t]["qOut"]) for t in ts[sel]])[order] * 1e6)

    init(a)
    st, mesh = G["st"], G["mesh"]
    q0 = max(float(r["qOut"]) for r in byt.values())
    u_m = q0 / st.mitral_area
    scale = (u_m / st.D_m) ** 2
    # Umbral Q y escala de velocidad fijos: de cuadros repartidos en el ciclo.
    samp = one[:: max(1, len(one) // 6)]
    qs, us = [], []
    for fr in samp:
        u = load(fr[1], fr[2], "velocity")
        m = mesh.copy(deep=False)
        m.point_data["velocity"] = u
        d = m.compute_derivative(scalars="velocity", qcriterion=True, gradient=False)
        qq = np.asarray(d["qcriterion"])
        qs.append(np.percentile(qq[qq > 0], 97))
        us.append(np.percentile(np.linalg.norm(u, axis=1), 99.5))
    a.qlevel = a.qplus * scale if a.qplus else float(np.median(qs))
    a.qplus = a.qlevel / scale
    a.umax = float(np.max(us))
    print(f"  Q = {a.qlevel:.3g} s^-2 = {a.qplus:.3g} (U_m/D_m)^2;  |u| hasta {a.umax:.2f} m/s")

    # Lista de cuadros: ciclo (interpolado si se pide), repetido --loops veces.
    seq = []
    cyc = one + [one[0]]
    for i in range(len(one)):
        for j in range(a.interp):
            seq.append((cyc[i], j / a.interp, cyc[i + 1] if j else None))
    seq = seq * a.loops
    if a.max_frames:
        seq = seq[: a.max_frames]
    a.tmp = os.path.abspath(a.out + ".frames")
    os.makedirs(a.tmp, exist_ok=True)
    jobs = [(i, fr, frac, other, a.orbit * i / max(len(seq) - 1, 1)) for i, (fr, frac, other) in enumerate(seq)]
    print(f"  {len(jobs)} cuadros, {len(jobs) / a.fps:.1f} s a {a.fps} cuadros/s; {a.jobs} procesos")
    if a.jobs > 1:
        with Pool(a.jobs, initializer=init, initargs=(a,)) as pool:
            for k, _ in enumerate(pool.imap(frame_job, jobs)):
                if k % 20 == 0:
                    print(f"    cuadro {k}/{len(jobs)}", flush=True)
    else:
        for k, j in enumerate(jobs):
            frame_job(j)
            if k % 20 == 0:
                print(f"    cuadro {k}/{len(jobs)}", flush=True)
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps),
           "-i", os.path.join(a.tmp, "frame_%05d.png"), "-c:v", "libx264", "-preset", a.preset,
           "-crf", str(a.crf), "-pix_fmt", "yuv420p", "-movflags", "+faststart", a.out]
    subprocess.run(cmd, check=True)
    if not a.keep:
        shutil.rmtree(a.tmp)
    print(f"Video: {a.out}")


if __name__ == "__main__":
    main()

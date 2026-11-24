#!/usr/bin/env python3
"""
compare_models.py -- N, CY, VE y TVE del atrio en el último ciclo, con las
mismas métricas para los cuatro, para juzgar si las diferencias entre
reologías son de un tamaño que sostenga un artículo.

    python3 compare_models.py [N=pref] [CY=pref] [VE=pref] [TVE=pref] [out=prefijo]

pref es el prefijo de salida de KokkosAtriumRheologyPETSc (output=...): se leen
pref.csv, pref.xdmf, pref.mesh.h5 y pref.velocity.NNNNNN.h5. Por defecto,
en el directorio actual: coarse3_N (o coarse_N si no está), coarse3_CY,
coarse3_VE, coarse3_TVE.

Escribe <out>_summary.txt, <out>_frames.csv y, si hay matplotlib, <out>.png.
La tasa de corte se calcula aquí, desde la velocidad de cada celda, para que
sea la misma en los cuatro modelos (no depende del campo shear_rate del driver).
"""
import csv
import os
import re
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("Falta h5py: pip install --user h5py")

MODELS = ["N", "CY", "VE", "TVE"]
PAIRS = [("CY", "N", "adelgazamiento"), ("VE", "CY", "elasticidad"),
         ("TVE", "VE", "tixotropía"), ("TVE", "N", "todo")]
LOW_SHEAR = (1.0, 10.0)        # s^-1: umbrales de estasis
SLOW = 0.05                    # m/s

args = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
pref = {}
for m in MODELS:
    if m in args:
        pref[m] = args[m]
    else:
        for p in (f"coarse3_{m}", f"coarse_{m}"):
            if os.path.exists(p + ".csv"):
                pref[m] = p
                break
OUT = args.get("out", "compare")
missing = [m for m in MODELS if m not in pref or not os.path.exists(pref[m] + ".csv")]
if missing:
    sys.exit("No encuentro las corridas de: " + ", ".join(missing))
lines = []


def say(s=""):
    print(s)
    lines.append(s)


# ---------------------------------------------------------------- series de tiempo
def series(p):
    rows = list(csv.DictReader(open(p + ".csv")))
    col = lambda k: np.array([float(r.get(k) or 0) for r in rows])
    return {k: col(k) for k in rows[0].keys()}


S = {m: series(pref[m]) for m in MODELS}
say("Corridas: " + ", ".join(f"{m}={pref[m]}" for m in MODELS))
ncyc = {m: int(S[m]["cycle"].max()) for m in MODELS}
say("Ciclos: " + ", ".join(f"{m} {ncyc[m]}" for m in MODELS))
last = min(ncyc.values())
if len(set(ncyc.values())) > 1:
    say(f"AVISO: no todas tienen los mismos ciclos; se compara el ciclo {last}.")


def cyc(m, c):
    s = S[m]
    k = s["cycle"] == c
    return {key: v[k] for key, v in s.items()}


L = {m: cyc(m, last) for m in MODELS}
dt = {m: np.median(np.diff(L[m]["t"])) for m in MODELS}
if max(dt.values()) - min(dt.values()) > 1e-9:
    say("AVISO: los pasos de tiempo difieren: " + str(dt))


def integral(m, key):
    return float(np.sum(L[m][key]) * dt[m])


period = {m: L[m]["t"][-1] - L[m]["t"][0] + dt[m] for m in MODELS}
rows = [
    ("Q mitral máximo (mL/s)", lambda m: 1e6 * L[m]["qOut"].max()),
    ("Volumen por la mitral (mL)", lambda m: 1e6 * integral(m, "qOut")),
    ("|u_i| máximo (m/s)", lambda m: L[m]["maxSpeedComponent"].max()),
    ("Energía cinética media (mJ)", lambda m: 1e3 * integral(m, "kineticEnergy") / period[m]),
    ("Energía cinética máxima (mJ)", lambda m: 1e3 * L[m]["kineticEnergy"].max()),
    # En VE y TVE la viscosa es la del plasma; la polimérica media en un ciclo
    # periódico es la disipación del polímero. La total es la comparable.
    ("Disipación total media (mW)", lambda m: 1e3 * (integral(m, "viscousDissipation")
                                                    + integral(m, "polymerPower")) / period[m]),
    ("  viscosa (plasma en VE/TVE)", lambda m: 1e3 * integral(m, "viscousDissipation") / period[m]),
    ("  polimérica", lambda m: 1e3 * integral(m, "polymerPower") / period[m]),
    ("Energía elástica máxima (uJ)", lambda m: 1e6 * L[m]["elasticEnergy"].max()),
    ("Wi máximo", lambda m: L[m]["maxWeissenberg"].max()),
    ("tr(c) máximo", lambda m: L[m]["maxTraceConformation"].max()),
    ("Tamaño medio de agregado", lambda m: L[m]["meanAggregateSize"].mean()),
    ("Presión ventricular máx. (Pa)", lambda m: L[m]["pOut"].max()),
]
say(f"\n=== Escalares del ciclo {last}")
say(f"{'':32s}" + "".join(f"{m:>11s}" for m in MODELS) + "   CY/N   VE/CY  TVE/VE")
ELASTIC_ONLY = ("polimérica", "elástica", "Wi", "tr(c)", "agregado")
for name, f in rows:
    v = {m: f(m) for m in MODELS}
    if any(k in name for k in ELASTIC_ONLY):
        v["N"] = v["CY"] = float("nan")
    rel = lambda a, b: (f"{100 * (v[a] - v[b]) / abs(v[b]):+6.1f}%"
                        if np.isfinite(v[b]) and abs(v[b]) > 1e-30 else "      -")
    say(f"{name:32s}" + "".join(f"{v[m]:11.4g}" for m in MODELS)
        .replace("        nan", "          -") + f"  {rel('CY', 'N')} {rel('VE', 'CY')} {rel('TVE', 'VE')}")


def waveform(a, b, key):
    x, y = L[a][key], L[b][key]
    n = min(len(x), len(y))
    return np.linalg.norm(x[:n] - y[:n]) / max(np.linalg.norm(y[:n]), 1e-300)


say(f"\n=== Diferencia de forma de onda en el ciclo {last} (L2 relativa en el tiempo)")
for a, b, what in PAIRS:
    say(f"  {a:>3s} vs {b:<3s} ({what:14s}): Q mitral {100 * waveform(a, b, 'qOut'):5.1f}%"
        f"   energía cinética {100 * waveform(a, b, 'kineticEnergy'):5.1f}%")

say("\n=== Periodicidad: ciclo último vs anterior (L2 relativa)")
for m in MODELS:
    if ncyc[m] < 2:
        say(f"  {m}: un solo ciclo")
        continue
    a, b = cyc(m, last), cyc(m, last - 1)
    n = min(len(a["qOut"]), len(b["qOut"]))
    d = lambda k: np.linalg.norm(a[k][:n] - b[k][:n]) / max(np.linalg.norm(a[k][:n]), 1e-300)
    extra = f"   n medio {100 * d('meanAggregateSize'):.2f}%" if m == "TVE" else ""
    say(f"  {m:>3s}: Q mitral {100 * d('qOut'):.2f}%   energía cinética {100 * d('kineticEnergy'):.2f}%{extra}")

# ---------------------------------------------------------------- campos
def frames(p):
    xd = open(p + ".xdmf").read()
    times = [float(t) for t in re.findall(r'Time Value="([^"]+)"', xd)]
    return {round(t, 9): k for k, t in enumerate(times)
            if os.path.exists(f"{p}.velocity.{k:06d}.h5")}


F = {m: frames(pref[m]) for m in MODELS}
t0, t1 = max(L[m]["t"][0] for m in MODELS) - 1e-9, min(L[m]["t"][-1] for m in MODELS) + 1e-9
common = sorted(t for t in set.intersection(*(set(F[m]) for m in MODELS)) if t0 <= t <= t1)
say(f"\n=== Campos: {len(common)} instantes comunes en el ciclo {last}")
if not common:
    say("No hay instantes comunes (¿output_every distinto?): solo escalares.")
else:
    mesh = h5py.File(pref["N"] + ".mesh.h5", "r")
    X = mesh["Mesh/Geometry/Vertices"][:]
    T = mesh["Mesh/XDMF/Topology"][:].reshape(-1, 5)[:, 1:].astype(np.int64)
    for m in MODELS[1:]:
        Xm = h5py.File(pref[m] + ".mesh.h5", "r")["Mesh/Geometry/Vertices"]
        if Xm.shape != X.shape or np.abs(Xm[:] - X).max() > 1e-9 * np.abs(X).max():
            sys.exit(f"La malla de {m} no es la de N: no se pueden comparar campos.")
    if np.ptp(X, axis=0).max() > 1.0:     # mm -> m
        X = X * 1e-3
    E = X[T[:, 1:]] - X[T[:, [0]]]                        # aristas de cada celda
    det = np.linalg.det(E)
    vol = np.abs(det) / 6
    Einv = np.linalg.inv(E)                               # gradientes baricéntricos
    G = np.concatenate([-Einv.sum(axis=2, keepdims=True), Einv], axis=2)   # (c, 3, 4)
    lumped = np.zeros(len(X))
    np.add.at(lumped, T.ravel(), np.repeat(vol / 4, 4))
    say(f"  malla: {len(X)} vértices, {len(T)} celdas, volumen {1e6 * vol.sum():.1f} mL")

    def velocity(m, t):
        return h5py.File(f"{pref[m]}.velocity.{F[m][t]:06d}.h5", "r")["GridFunction/Values/Data"][:]

    def shear(u):                                         # sqrt(2 D:D) por celda
        J = np.einsum("cin,cnk->cki", G, u[T])           # J[c, k, i] = d u_k / d x_i
        D = 0.5 * (J + J.transpose(0, 2, 1))
        return np.sqrt(2 * np.einsum("cij,cij->c", D, D))

    norm = lambda w: np.sqrt(np.sum(lumped[:, None] * w * w))
    rec, mean = [], {m: np.zeros_like(X) for m in MODELS}
    for t in common:
        U = {m: velocity(m, t) for m in MODELS}
        r = {"t": t}
        for a, b, _ in PAIRS:
            r[f"du_{a}_{b}"] = norm(U[a] - U[b]) / max(norm(U[b]), 1e-300)
        # Piso: el mismo modelo, mismo instante del ciclo anterior. Si el flujo
        # no se repite de un ciclo al otro, una diferencia entre modelos de ese
        # tamaño no dice nada de la reología.
        for m in MODELS:
            tp = round(t - period[m], 9)
            r[f"cyc_{m}"] = (norm(U[m] - velocity(m, tp)) / max(norm(U[m]), 1e-300)
                             if tp in F[m] else float("nan"))
        for m in MODELS:
            g = shear(U[m])
            for s in LOW_SHEAR:
                r[f"low{s:g}_{m}"] = vol[g < s].sum() / vol.sum()
            r[f"slow_{m}"] = lumped[np.linalg.norm(U[m], axis=1) < SLOW].sum() / lumped.sum()
            mean[m] += U[m] / len(common)
        rec.append(r)
    with open(OUT + "_frames.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rec[0].keys()))
        w.writeheader()
        w.writerows(rec)
    col = lambda k: np.array([r[k] for r in rec])
    say("  Diferencia de velocidad entre modelos, ||u_a - u_b|| / ||u_b|| en el volumen:")
    for a, b, what in PAIRS:
        d = col(f"du_{a}_{b}")
        dm = norm(mean[a] - mean[b]) / max(norm(mean[b]), 1e-300)
        say(f"    {a:>3s} vs {b:<3s} ({what:14s}): media en el ciclo {100 * d.mean():5.1f}%, "
            f"máx {100 * d.max():5.1f}%; campo medio del ciclo {100 * dm:5.1f}%")
    floor = {m: np.nanmean(col(f"cyc_{m}")) if np.isfinite(col(f"cyc_{m}")).any() else float("nan")
             for m in MODELS}
    say("  Piso (mismo modelo, ciclo anterior, mismo instante): "
        + ", ".join(f"{m} {100 * floor[m]:.1f}%" for m in MODELS))
    say("  Una diferencia entre modelos que no supera claramente el piso es variabilidad")
    say("  ciclo a ciclo del flujo, no efecto de la reología.")
    say("  Fracción de volumen de baja tasa de corte y baja velocidad, media en el ciclo:")
    say(f"    {'':22s}" + "".join(f"{m:>9s}" for m in MODELS))
    for s in LOW_SHEAR:
        say(f"    {'gamma < ' + format(s, 'g') + ' 1/s':22s}"
            + "".join(f"{100 * col(f'low{s:g}_{m}').mean():8.2f}%" for m in MODELS))
    say(f"    {'|u| < ' + format(SLOW, 'g') + ' m/s':22s}"
        + "".join(f"{100 * col(f'slow_{m}').mean():8.2f}%" for m in MODELS))

open(OUT + "_summary.txt", "w").write("\n".join(lines) + "\n")
print(f"\nEscrito: {OUT}_summary.txt" + (f", {OUT}_frames.csv" if common else ""))

# ---------------------------------------------------------------- figura
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)
fig, ax = plt.subplots(1, 3 if common else 2, figsize=(15 if common else 10, 3.6))
for m in MODELS:
    tt = L[m]["t"] - L[m]["t"][0]
    ax[0].plot(tt, 1e6 * L[m]["qOut"], label=m)
    ax[1].plot(tt, 1e3 * L[m]["kineticEnergy"], label=m)
ax[0].set(xlabel="t en el ciclo (s)", ylabel="Q mitral (mL/s)")
ax[1].set(xlabel="t en el ciclo (s)", ylabel="energía cinética (mJ)")
ax[0].legend()
if common:
    tc = np.array(common) - L["N"]["t"][0]
    for a, b, what in PAIRS:
        ax[2].plot(tc, 100 * np.array([r[f"du_{a}_{b}"] for r in rec]), "o-", label=f"{a}-{b} ({what})")
    for m in MODELS:
        y = np.array([r[f"cyc_{m}"] for r in rec])
        if np.isfinite(y).any():
            ax[2].plot(tc, 100 * y, ":", color="gray", lw=1)
    ax[2].plot([], [], ":", color="gray", label="piso: ciclo anterior")
    ax[2].set(xlabel="t en el ciclo (s)", ylabel="||u_a-u_b|| / ||u_b|| (%)")
    ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT + ".png", dpi=150)
print(f"Figura: {OUT}.png")

#!/usr/bin/env python3
"""
stats_models.py -- estadística entre ciclos de varias corridas en la misma malla.

    python3 stats_models.py [runs=N:prod_N,CY:prod_CY,VE:prod_VE,TVE:prod_TVE]
                            [pairs=CY/N,VE/CY,TVE/VE,TVE/N] [skip=2] [detrend=1]
                            [region=x,y,z,r] [out=stats]

runs     etiqueta:prefijo de salida de KokkosAtriumRheologyPETSc (output=...).
         Por defecto N, CY, VE y TVE con prefijo prod_ (o fine_ si no existe).
pairs    comparaciones a/b.
skip     ciclos que se descartan al comienzo de los ciclos comunes a todas las
         corridas (tras un reinicio, lo que tarda cada reología en adaptarse).
detrend  1: en cada fase se quita la tendencia lineal entre ciclos antes de
         medir la variabilidad (la deriva lenta del sistema no cuenta como
         variabilidad entre latidos). Necesita 4 ciclos o más.
region   esfera (centro y radio en mm, en las coordenadas de la malla) donde
         además se miden la estasis y la variabilidad, p. ej. la orejuela.

En cada fase del ciclo (cada instante escrito):
  media de fase   ubar = promedio de u sobre los ciclos;
  fluctuación     u' = u - ubar (- tendencia): lo que cambia de un latido a otro;
  intensidad      I = sqrt(<||u'||^2>) / ||ubar||.
Entre corridas:
  diferencia de las medias de fase ||ubar_a - ubar_b|| / ||ubar_b|| contra su
  ruido con n ciclos; señal/ruido > 2 indica efecto, no muestreo.
Escalares por ciclo: comparación pareada ciclo a ciclo (prueba t pareada): las
corridas que parten del mismo estado comparten la deriva, que así se cancela.

Escribe <out>_summary.txt, <out>_phase.csv y, si hay matplotlib, <out>.png.
"""
import csv
import os
import re
import sys
import warnings

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("Falta h5py: sudo apt install python3-h5py")
try:
    from scipy.stats import ttest_rel
except ImportError:
    ttest_rel = None

warnings.filterwarnings("ignore")
RHO = 1060.0
args = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
if "runs" in args:
    RUNS = dict(x.split(":", 1) for x in args["runs"].split(","))
else:
    base = "prod" if os.path.exists("prod_N.csv") else "fine"
    RUNS = {m: f"{base}_{m}" for m in ("N", "CY", "VE", "TVE")}
MODELS = list(RUNS)
WHAT = {"CY/N": "adelgazamiento", "VE/CY": "elasticidad", "TVE/VE": "tixotropía", "TVE/N": "todo"}
if "pairs" in args:
    PAIRS = [tuple(p.split("/")) for p in args["pairs"].split(",")]
else:
    PAIRS = [tuple(p.split("/")) for p in WHAT if all(m in RUNS for m in p.split("/"))]
SKIP = int(args.get("skip", 2))
OUT = args.get("out", "stats")
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def label(a, b):
    return f"{a:>7s} vs {b:<7s}" + (f" ({WHAT[f'{a}/{b}']})" if f"{a}/{b}" in WHAT else "")


for m in MODELS:
    if not os.path.exists(RUNS[m] + ".csv"):
        sys.exit(f"No encuentro {RUNS[m]}.csv")

# ---------------------------------------------------------------- series por ciclo
S = {}
for m in MODELS:
    rows = list(csv.DictReader(open(RUNS[m] + ".csv")))
    S[m] = {k: np.array([float(r.get(k) or 0) for r in rows]) for k in rows[0]}
dt = np.median(np.diff(S[MODELS[0]]["t"]))
steps = int(max(np.sum(S[MODELS[0]]["cycle"] == c) for c in np.unique(S[MODELS[0]]["cycle"])))
period = steps * dt
full = [set(int(c) for c in np.unique(S[m]["cycle"]) if np.sum(S[m]["cycle"] == c) == steps) for m in MODELS]
common = sorted(set.intersection(*full))
cycles = common[SKIP:]
n = len(cycles)
if n < 2:
    sys.exit(f"Ciclos completos comunes {common}; con skip={SKIP} quedan {n}: hacen falta 2.")
DETREND = int(args.get("detrend", 1)) == 1 and n >= 4
say("Corridas: " + ", ".join(f"{m}={RUNS[m]}" for m in MODELS))
say(f"Ciclos comunes {common[0]}-{common[-1]}; se usan {cycles[0]}-{cycles[-1]} (n = {n}); "
    f"periodo {period:.3f} s; tendencia {'quitada' if DETREND else 'no quitada'}")


def per_cycle(m, f):
    s = S[m]
    return np.array([f({k: v[s["cycle"] == c] for k, v in s.items()}) for c in cycles])


scalars = {
    "Energía cinética media (mJ)": lambda s: 1e3 * s["kineticEnergy"].mean(),
    "Disipación total media (mW)": lambda s: 1e3 * (s["viscousDissipation"] + s["polymerPower"]).mean(),
    "|u_i| máximo (m/s)": lambda s: s["maxSpeedComponent"].max(),
    "Q mitral máximo (mL/s)": lambda s: 1e6 * s["qOut"].max(),
    "Presión ventricular máx. (kPa)": lambda s: 1e-3 * s["pOut"].max(),
}
SC = {name: {m: per_cycle(m, f) for m in MODELS} for name, f in scalars.items()}


# ---------------------------------------------------------------- campos por fase
def frame_index(p):
    times = [float(t) for t in re.findall(r'Time Value="([^"]+)"', open(p + ".xdmf").read())]
    out = {}
    for k, t in enumerate(times):
        if t <= 0 or not os.path.exists(f"{p}.velocity.{k:06d}.h5"):
            continue
        c = int(np.ceil(t / period - 1e-6))
        out[(c, round(t - (c - 1) * period, 6))] = k
    return out


FI = {m: frame_index(RUNS[m]) for m in MODELS}
phases = sorted({ph for (c, ph) in FI[MODELS[0]] if c in cycles})
phases = [ph for ph in phases if all((c, ph) in FI[m] for m in MODELS for c in cycles)]
say(f"Fases con todos los ciclos y corridas: {len(phases)}")
if not phases:
    sys.exit("No hay instantes comunes.")

mesh = h5py.File(RUNS[MODELS[0]] + ".mesh.h5", "r")
X = mesh["Mesh/Geometry/Vertices"][:]
T = mesh["Mesh/XDMF/Topology"][:].reshape(-1, 5)[:, 1:].astype(np.int64)
unit = 1e-3 if np.ptp(X, axis=0).max() > 1.0 else 1.0       # -> m
X = X * unit
E = X[T[:, 1:]] - X[T[:, [0]]]
vol = np.abs(np.linalg.det(E)) / 6
Einv = np.linalg.inv(E)
G = np.concatenate([-Einv.sum(axis=2, keepdims=True), Einv], axis=2)
del E, Einv
lump = np.zeros(len(X))
np.add.at(lump, T.ravel(), np.repeat(vol / 4, 4))
say(f"Malla: {len(X)} vértices, {len(T)} celdas, {1e6 * vol.sum():.1f} mL")
REGION = None
if "region" in args:
    x, y, z, r = (float(v) * 1e-3 for v in args["region"].split(","))
    inside_v = np.linalg.norm(X - [x, y, z], axis=1) < r
    inside_c = np.linalg.norm(X[T].mean(axis=1) - [x, y, z], axis=1) < r
    REGION = (inside_v, inside_c)
    say(f"Región: esfera de {1e3 * r:g} mm, {1e6 * vol[inside_c].sum():.2f} mL, {inside_v.sum()} vértices")
    if inside_c.sum() == 0:
        sys.exit("La región no contiene celdas.")


def load(m, c, ph):
    return h5py.File(f"{RUNS[m]}.velocity.{FI[m][(c, ph)]:06d}.h5", "r")["GridFunction/Values/Data"][:]


def sq(w, mask=None):
    wl = lump if mask is None else lump * mask
    return float(np.sum(wl[:, None] * w * w))


def shear(u):
    J = np.einsum("cin,cnk->cki", G, u[T])
    D = 0.5 * (J + J.transpose(0, 2, 1))
    return np.sqrt(2 * np.einsum("cij,cij->c", D, D))


kc = np.arange(n) - (n - 1) / 2                      # índice de ciclo centrado
dof = n - 2 if DETREND else n - 1
rec = []
fluct = {m: np.zeros(n) for m in MODELS}
low = {(m, s, z): np.zeros(n) for m in MODELS for s in (1.0, 10.0) for z in ("total", "región")}
for ph in phases:
    U = {m: np.stack([load(m, c, ph) for c in cycles]) for m in MODELS}
    r = {"fase": ph}
    mean, var = {}, {}
    for m in MODELS:
        mean[m] = U[m].mean(axis=0)
        dev = U[m] - mean[m]
        if DETREND:
            slope = np.tensordot(kc, dev, axes=1) / np.sum(kc * kc)
            dev = dev - kc[:, None, None] * slope
        e = np.array([sq(dev[i]) for i in range(n)]) * n / dof
        fluct[m] += 0.5 * RHO * e / len(phases)
        var[m] = e.mean()
        r[f"I_{m}"] = np.sqrt(var[m] / sq(mean[m]))
        if REGION is not None:
            vr = np.mean([sq(dev[i], REGION[0]) for i in range(n)]) * n / dof
            r[f"Ireg_{m}"] = np.sqrt(vr / sq(mean[m], REGION[0]))
        for i in range(n):
            g = shear(U[m][i])
            for s in (1.0, 10.0):
                low[(m, s, "total")][i] += vol[g < s].sum() / vol.sum() / len(phases)
                if REGION is not None:
                    c = REGION[1]
                    low[(m, s, "región")][i] += vol[c & (g < s)].sum() / vol[c].sum() / len(phases)
    for a, b in PAIRS:
        d = np.sqrt(sq(mean[a] - mean[b]) / sq(mean[b]))
        noise = np.sqrt((var[a] + var[b]) / n / sq(mean[b]))
        r[f"D_{a}_{b}"], r[f"ruido_{a}_{b}"] = d, noise
    rec.append(r)
    print(f"  fase {ph:.3f} s lista", flush=True)

SC["Energía de fluctuación media (uJ)"] = {m: 1e6 * fluct[m] for m in MODELS}
for s in (1.0, 10.0):
    SC[f"Volumen con gamma < {s:g} 1/s (%)"] = {m: 100 * low[(m, s, "total")] for m in MODELS}
    if REGION is not None:
        SC[f"Región: volumen con gamma < {s:g} 1/s (%)"] = {m: 100 * low[(m, s, "región")] for m in MODELS}

# ---------------------------------------------------------------- resumen
col = lambda key: np.array([r[key] for r in rec])
say("\n=== Variabilidad entre latidos: I = sqrt(<||u'||^2>) / ||ubar||, media en el ciclo")
say("    " + "  ".join(f"{m}: {100 * col(f'I_{m}').mean():5.1f}%" for m in MODELS))
if REGION is not None:
    say("    región: " + "  ".join(f"{m}: {100 * col(f'Ireg_{m}').mean():5.1f}%" for m in MODELS))
for a, b in PAIRS:
    ia, ib = col(f"I_{a}"), col(f"I_{b}")
    say(f"    {label(a, b)}: {100 * (ia.mean() / ib.mean() - 1):+6.1f}%, {a} menor en "
        f"{np.sum(ia < ib)} de {len(ia)} fases")

say("\n=== Medias de fase: diferencia entre corridas contra el ruido de n ciclos")
for a, b in PAIRS:
    d, z = col(f"D_{a}_{b}"), col(f"ruido_{a}_{b}")
    say(f"    {label(a, b)}: diferencia {100 * d.mean():5.1f}%, ruido {100 * z.mean():5.1f}%, "
        f"señal/ruido {np.mean(d / z):4.2f} (> 2 en {np.sum(d > 2 * z)} de {len(d)} fases)")

say("\n=== Escalares por ciclo: media ± desviación entre ciclos")
say(f"    {'':42s}" + "".join(f"{m:>18s}" for m in MODELS))
for name, v in SC.items():
    say(f"    {name:42s}" + "".join(f"{v[m].mean():10.4g} ±{v[m].std(ddof=1):7.2g}" for m in MODELS))

say("\n=== Comparación pareada ciclo a ciclo: diferencia media ± error estándar (prueba t pareada)")
for name, v in SC.items():
    parts = []
    for a, b in PAIRS:
        if np.any(v[b] == 0):
            parts.append(f"{a}/{b}      -")
            continue
        rel = 100 * (v[a] - v[b]) / np.abs(v[b])
        se = rel.std(ddof=1) / np.sqrt(n)
        p = ttest_rel(v[a], v[b]).pvalue if ttest_rel is not None else float("nan")
        parts.append(f"{a}/{b} {rel.mean():+6.1f}±{se:4.1f}% (p={p:.2g})")
    say(f"    {name:42s} " + "   ".join(parts))

say("\n=== Valores por ciclo (deriva)")
for name in ("Energía cinética media (mJ)", "Q mitral máximo (mL/s)", "Energía de fluctuación media (uJ)"):
    say(f"    {name}, ciclos {cycles[0]}-{cycles[-1]}")
    for m in MODELS:
        v = SC[name][m]
        drift = 100 * (v[-1] / v[0] - 1) / max(n - 1, 1)
        say(f"      {m:>7s}: " + "  ".join(f"{x:8.4g}" for x in v) + f"   ({drift:+.1f}%/ciclo)")

with open(OUT + "_phase.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rec[0].keys()))
    w.writeheader()
    w.writerows(rec)
open(OUT + "_summary.txt", "w").write("\n".join(lines) + "\n")
print(f"\nEscrito: {OUT}_summary.txt, {OUT}_phase.csv")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)
ph = np.array(phases)
fig, ax = plt.subplots(1, 3, figsize=(15, 3.8))
for m in MODELS:
    ax[0].plot(ph, 100 * col(f"I_{m}"), "o-", label=m)
ax[0].set(xlabel="fase del ciclo (s)", ylabel="intensidad entre latidos I (%)")
ax[0].legend()
for a, b in PAIRS:
    l, = ax[1].plot(ph, 100 * col(f"D_{a}_{b}"), "o-", label=f"{a}-{b}")
    ax[1].plot(ph, 200 * col(f"ruido_{a}_{b}"), ":", color=l.get_color())
ax[1].set(xlabel="fase del ciclo (s)", ylabel="||ubar_a - ubar_b|| / ||ubar_b|| (%)",
          title="punteado: 2 x ruido de n ciclos")
ax[1].legend(fontsize=8)
name = "Energía de fluctuación media (uJ)"
for i, m in enumerate(MODELS):
    y = SC[name][m]
    ax[2].errorbar([i], [y.mean()], yerr=[y.std(ddof=1)], fmt="o", capsize=5)
    ax[2].plot(np.full(len(y), i) + 0.15, y, ".", color="gray")
ax[2].set_xticks(range(len(MODELS)), MODELS)
ax[2].set(ylabel=name, title="media ± desviación entre ciclos")
fig.tight_layout()
fig.savefig(OUT + ".png", dpi=150)
print(f"Figura: {OUT}.png")

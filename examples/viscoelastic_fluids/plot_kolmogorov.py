#!/usr/bin/env python3
"""Plot test_kolmogorov_viscoelastic output against Berti & Boffetta (2010).

Usage:
    python3 plot_kolmogorov.py kf_steady_Wi8.csv kf_transition_Wi16.csv [...]
Writes kf_kolmogorov.png with four panels:
  (a) K(t) and (b) Sigma(t) = <tr c>(t)                      -- their Fig. 3
  (c) <Sigma> and <K> against the a-posteriori Wi = tau U/L,
      with the laminar predictions 2 + Wi^2/4 and U^2/4         -- their Fig. 1
  (d) frequency spectrum S_K(omega) of K(t), with omega_A = 2 pi U/L0
      and omega_tau = 2 pi/tau marked (they find omega_A <~ omega* < omega_tau)
                                                                 -- their Fig. 4
Averages and spectra use the second half of each run. U0 = 4, L = 1, L0 = 2 pi
and L = 1/4 as in their Table I and Figs. 5, 7; Wi0 is read from the file name (kf_<case>_Wi<Wi0>.csv).
"""
import csv
import math
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

U0, L, L0 = 4.0, 0.25, 2.0 * math.pi   # L = 1/4: four periods in the box

files = sys.argv[1:]
if not files:
    sys.exit(__doc__)

runs = []
for f in files:
    rows = list(csv.DictReader(open(f)))
    col = lambda k: np.array([float(r[k]) for r in rows])
    m = re.search(r"Wi([0-9.]+?)\.csv$", f)
    wi0 = float(m.group(1)) if m else float("nan")
    t, K, S, U = col("t"), col("K"), col("Sigma"), col("U")
    half = t >= 0.5 * t[-1]
    tau = wi0 * L / U0
    Umean = U[half].mean()
    runs.append(dict(name=f, t=t, K=K, S=S, half=half, tau=tau, wi0=wi0,
                     Wi=tau * Umean / L, U=Umean, Kmean=K[half].mean(), Kstd=K[half].std(),
                     Smean=S[half].mean(), Sstd=S[half].std()))

fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
for i, r in enumerate(runs):
    c = f"C{i}"
    label = f"Wi0={r['wi0']:g} (Wi={r['Wi']:.1f})"
    ax[0, 0].plot(r["t"], r["K"], color=c, lw=1, label=label)
    ax[0, 1].plot(r["t"], r["S"], color=c, lw=1, label=label)
ax[0, 0].set_xlabel("$t$"); ax[0, 0].set_ylabel(r"$K = \langle|\mathbf{u}|^2\rangle/2$")
ax[0, 1].set_xlabel("$t$"); ax[0, 1].set_ylabel(r"$\Sigma = \langle \mathrm{tr}\,\mathbf{c}\rangle$")
ax[0, 0].set_title("(a) kinetic energy", fontsize=9)
ax[0, 1].set_title("(b) mean square elongation", fontsize=9)
for a in ax[0]:
    a.grid(alpha=0.3); a.legend(fontsize=7)

# (c) means against the a-posteriori Wi, with the laminar predictions
wi = np.linspace(4, 25, 100)
a = ax[1, 0]
a.semilogy(wi, 2 + wi ** 2 / 4, "k--", lw=0.8, label=r"laminar $2 + \mathrm{Wi}^2/4$")
for i, r in enumerate(runs):
    a.errorbar(r["Wi"], r["Smean"], yerr=r["Sstd"], fmt="o", color=f"C{i}",
               label=rf"$\langle\Sigma\rangle$, Wi0={r['wi0']:g}")
a.axvline(10, color="0.6", lw=0.8, ls=":")
a.text(10.2, 3, r"Wi$_c\approx 10$", fontsize=7, color="0.4")
a.set_xlabel(r"Wi $= \tau U/L$ (a posteriori)"); a.set_ylabel(r"$\langle\Sigma\rangle$")
a.set_title(r"(c) $\langle\Sigma\rangle$ vs Wi; inset text: $\langle K\rangle$ vs laminar $U^2/4$", fontsize=9)
for i, r in enumerate(runs):
    a.text(0.03, 0.95 - 0.07 * i,
           rf"Wi0={r['wi0']:g}: $\langle K\rangle$={r['Kmean']:.3f} $\pm$ {r['Kstd']:.1e}"
           rf" (laminar $U_0^2/4$={U0**2/4:g})",
           transform=a.transAxes, fontsize=7, va="top", color=f"C{i}")
a.grid(alpha=0.3, which="both"); a.legend(fontsize=7, loc="lower right")

# (d) spectrum of K(t) on the second half
a = ax[1, 1]
for i, r in enumerate(runs):
    t, K = r["t"][r["half"]], r["K"][r["half"]]
    if len(t) < 16 or K.std() < 1e-10 * max(1.0, abs(K.mean())):
        continue   # steady run: nothing to transform
    tt = np.linspace(t[0], t[-1], len(t))
    Ku = np.interp(tt, t, K) - K.mean()
    dt = tt[1] - tt[0]
    Tmax = tt[-1] - tt[0]
    w = 2 * math.pi * np.fft.rfftfreq(len(tt), dt)
    S = np.abs(np.fft.rfft(Ku) * dt) ** 2 / Tmax
    a.plot(w, S, color=f"C{i}", lw=1, label=f"Wi0={r['wi0']:g}")
    a.axvline(2 * math.pi * r["U"] / L0, color=f"C{i}", ls="--", lw=0.8)
    a.axvline(2 * math.pi / r["tau"], color=f"C{i}", ls=":", lw=0.8)
a.set_xlim(0, 20)
a.set_xlabel(r"$\omega$"); a.set_ylabel(r"$S_K(\omega)$")
a.set_title(r"(d) spectrum of $K(t)$; dashed $\omega_A = 2\pi U/L_0$, dotted $\omega_\tau = 2\pi/\tau$",
            fontsize=9)
a.grid(alpha=0.3)
if a.lines:
    a.legend(fontsize=7)

fig.tight_layout()
fig.savefig("kf_kolmogorov.png", dpi=150)
for r in runs:
    print(f"{r['name']}: Wi0={r['wi0']:g}  U={r['U']:.3f}  Wi={r['Wi']:.2f}  "
          f"<K>={r['Kmean']:.4f}+-{r['Kstd']:.1e}  <Sigma>={r['Smean']:.3f}+-{r['Sstd']:.1e}  "
          f"laminar Sigma(Wi)={2 + r['Wi']**2/4:.3f}")
print("wrote kf_kolmogorov.png")

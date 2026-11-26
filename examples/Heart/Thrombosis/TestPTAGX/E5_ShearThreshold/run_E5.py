#!/usr/bin/env python3
"""E5 - Initiation on a TF patch under flow: shear-rate threshold (Shen 2008).

Experiment (Shen, Kastrup, Liu, Ismagilov, ATVB 28 (2008) 2035; supplement):
normal pooled plasma (corn trypsin inhibitor, < 1e4 platelets/uL), silica
capillaries, a 200 um TF patch (lipid 79.5 % DLPC / 20 % PS, the Kastrup 2006
composition), steady flow:
  - 5 1/s: clotting on the patch at ~210 s;
  - 40 1/s: no initiation in 300 s (lowering the flow to 5 1/s then gives
    clotting at ~100 s);
  - 15 1/s (700 um capillary) initiates, 70 1/s (150 um, same velocity) does
    not: the control variable is the shear rate.

Model: PTAGX kinetics (PTAGXBench), planar channel, exact Poiseuille profile
(H = 450 um, the capillary diameter; the 0.25 mm next to the wall is meshed),
a 200 um patch 1 mm downstream of the inlet. Xa dose from the TF density of
Kastrup 2006 and kappa_TF of E1: j_X = 4e-13 mol m^-2 s^-1 (not fitted).
Initiation = thrombin above 100 nM.

Sets: MODEL (verdict) and SURFACE (MODEL plus sigma, sensitivity). sigma of
SURFACE is calibrated on the 5 1/s case (calibrate_sigma.py), so for SURFACE
the first criterion is a calibration check; for MODEL nothing is fitted here.

Criteria (MODEL):
  1. initiates at 5 1/s within 300 s (measured ~210 s) [calibration];
  2. does not initiate at 40 1/s within 300 s [prediction];
  3. initiates at 15 1/s in a 700 um capillary [prediction];
  4. does not initiate at 70 1/s in a 150 um capillary [prediction].
(the patch length of the last two is taken as 200 um; the meshed near-wall
band, 0.25 mm, is wider than half of the 150 um channel, but the exact
Poiseuille profile is what enters.)
Cost: 6 runs of ~2 min.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

pc = px.pc
L, H = 200e-6, 450e-6
X0 = 1.0e-3
J_X = 4.0e-13
SHEARS = [5.0, 20.0, 40.0, 80.0]
BURST = 1.0e-4          # 100 nM
T_END = 600.0
CITES = [
    "Shen F, Kastrup CJ, Liu Y, Ismagilov RF. Threshold response of initiation of blood coagulation "
    "by tissue factor in patterned microfluidic capillaries is controlled by shear rate. "
    "Arterioscler Thromb Vasc Biol 28 (2008) 2035-2041, and its Supplemental Information",
]


def channel_mesh(out):
    return pc.make_mesh(os.path.join(out, "channel.mesh"), Lx=X0 + L + 1.5e-3, Ly=0.25e-3,
                        x0=X0, x1=X0 + L, margin=50e-6, hx_patch=8e-6, hx_far=60e-6,
                        rx=1.12, hy0=2e-6, ry=1.18, hy_max=30e-6)


def case_options(mesh, gamma, j, model, height=H):
    return dict(px.kinetics(model, patch=(X0, X0 + L)), pg_mesh=mesh, pg_inlet=1,
                pg_reactive=5, pg_H=height, pg_gamma=gamma, pg_crosswind=0.7, pg_dt=2.0,
                pg_t_end=T_END, pg_csv_interval=6, pg_output_interval=0, pg_log_every=0,
                pg_j_x=j, pg_stop_thrombin=BURST)


def main():
    ap = pc.parser(__doc__)
    args = ap.parse_args()
    out = pc.results_dir(args, __file__)
    mesh = channel_mesh(out)
    cases = {}
    for sn, model in (("model", px.MODEL),):
        for g in SHEARS:
            cases[f"{sn}_g{g:g}"] = case_options(mesh, g, J_X, model)
        # the two capillaries of the supplement with the same mean velocity
        cases[f"{sn}_g15_H700um"] = case_options(mesh, 15.0, J_X, model, height=700e-6)
        cases[f"{sn}_g70_H150um"] = case_options(mesh, 70.0, J_X, model, height=150e-6)

    csvs = pc.run_many(args, out, cases)

    res = {}
    lines = ["| set | shear (1/s) | initiates | time (s) | max T (nM) |", "|---|---|---|---|---|"]
    for n in cases:
        d = pc.read_csv(csvs[n])
        t0 = pc.first_time(d["t"], d["maxT"] >= BURST)
        res[n] = t0
        s, g = n.split("_g", 1)
        lines.append(f"| {s} | {g} | {'si' if np.isfinite(t0) else 'no'} | {t0:.0f} | "
                     f"{d['maxT'].max()*1e6:.1f} |")
    lines += ["", f"MODEL: k_a = {px.MODEL['k_a']}, sigma = 0; SURFACE: sigma = "
              f"{px.SURFACE['sigma_s']:.2e} m, delta_s = {px.SURFACE['delta_s']*1e6:.0f} um; "
              "LEGACY: k_a = 0.10, sigma = 0.",
              "", "References:"] + [f"- {c}" for c in CITES]

    plt = pc.setup_matplotlib()
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    for n in cases:
        d = pc.read_csv(csvs[n])
        ax.plot(d["t"], d["maxT"] * 1e6, "-" if n.startswith("model") else (":" if n.startswith("surface") else "--"), label=n)
    ax.axhline(BURST * 1e6, color="k", lw=0.8, ls=":")
    ax.axvline(300, color="k", lw=0.8, ls=":")
    ax.set(xlabel="t (s)", ylabel="max thrombin (nM)", yscale="log", title="E5, 200 um patch")
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "E5_shear.png"))
    pc.write_summary(out, "E5 - Shear threshold of initiation, plasma (Shen 2008)", lines)

    def fmt(t):
        return f"{t:.0f} s" if np.isfinite(t) else "no"

    t5, t40 = res["model_g5"], res["model_g40"]
    t15, t70 = res["model_g15_H700um"], res["model_g70_H150um"]
    px.verdict(out, "E5 Shen 2008", [
        ("modelo inicia a 5 1/s en <= 300 s (medido ~210 s)",
         bool(np.isfinite(t5) and t5 <= 300), f"{t5:.0f} s"),
        ("modelo no inicia a 40 1/s en 300 s (prediccion)",
         bool(not np.isfinite(t40) or t40 > 300), f"{t40:.0f} s" if np.isfinite(t40) else "no inicia"),
        ("modelo inicia a 15 1/s en capilar de 700 um (prediccion)",
         bool(np.isfinite(t15)), f"{t15:.0f} s" if np.isfinite(t15) else "no inicia"),
        ("modelo no inicia a 70 1/s en capilar de 150 um (prediccion)",
         bool(not np.isfinite(t70) or t70 > 300), f"{t70:.0f} s" if np.isfinite(t70) else "no inicia"),
    ], notes=["ronda anterior con SURFACE (sigma 2e-3 calibrado a 5 1/s): "
              "g5 62 s, g10/g20/g40/g80 y los dos capilares no inician"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""E4 - Threshold patch size for clotting on tissue factor, no flow (Kastrup 2006).

Experiment (Kastrup et al., PNAS 103 (2006) 15747): normal pooled plasma
(corn trypsin inhibitor), quiescent chamber, lipid patches with tissue factor
(TF:lipid 2.5e-7, about 0.08 fmol/cm^2 = 8e-13 mol/m^2) in an inert lipid
background.
  - patches >= 100 um initiate clotting in < 3 min (40/44);
  - patches <= 50 um do not (28/28; background clotting at 32-75 min);
  - the clot then spreads as a front at 25-35 um/min.

Model: PTAGX kinetics, axisymmetric chamber (radius 1.25 mm, height 1.2 mm,
-pg_axisym), the p x p patch as a disc of equal area (radius p/sqrt(pi)) on
the axis. A planar 2D model would make the patch an infinite strip, which
feeds far more thrombin to its centre than a square does, and lowers the
threshold artificially. Xa dose from the measured TF density with the TF
turnover of E1 (kappa_TF = 0.5 1/s): j_X = 4e-13 mol m^-2 s^-1, not fitted.
Plasma held at the top. Parameter sets: MODEL (verdict) and SURFACE (MODEL plus
the sigma of E5, sensitivity). Nothing is fitted on this test.

Clotting = thrombin burst (max T > 100 nM, the excited state of the model)
within 25 min. A gel next to the patch (F >= K_F) can form without a burst, from
the patch's own thrombin; it is reported but is not the clotting of the assay.
The front speed is the slope of the gel height after the burst.
Criteria (MODEL): patches <= 50 um do not initiate within the run; patches
>= 100 um initiate within 3 min.
Cost: 10 runs of ~8 min.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

pc = px.pc
GAMMA_TF = 8.0e-13       # mol/m^2 (0.08 fmol/cm^2)
KAPPA_TF = 0.5           # 1/s, as in E1
J_X = GAMMA_TF * KAPPA_TF
WIDTHS = [25e-6, 50e-6, 100e-6, 200e-6, 400e-6]
SETS = {"model": px.MODEL, "surface": px.SURFACE}
CITES = [
    "Kastrup CJ, Runyon MK, Shen F, Ismagilov RF. Modular chemical mechanism predicts spatiotemporal "
    "dynamics of initiation in the complex network of hemostasis. PNAS 103 (2006) 15747-15752",
    "Kastrup CJ, Shen F, Runyon MK, Ismagilov RF. Characterization of the threshold response of "
    "initiation of blood clotting to stimulus patch size. Biophys J 93 (2007) 2969-2977",
]


def main():
    ap = pc.parser(__doc__)
    args = ap.parse_args()
    out = pc.results_dir(args, __file__)
    R, Ly = 1.25e-3, 1.2e-3
    t_end = 900 if args.quick else 1500
    cases = {}
    for p in WIDTHS:
        a = p / np.sqrt(np.pi)          # disc of the area of the p x p patch
        h = min(8e-6, a / 5)
        mesh = pc.make_mesh(os.path.join(out, f"axi_{p*1e6:.0f}um.mesh"), Lx=R, Ly=Ly,
                            x0=0.0, x1=a, margin=120e-6, hx_patch=h, hx_far=80e-6, rx=1.12,
                            hy0=h, ry=1.1, hy_max=80e-6, left=3, right=3, bottom=3, top=1)
        for sn, model in SETS.items():
            cases[f"p{p*1e6:.0f}_{sn}"] = dict(
                px.kinetics(model, patch=(0.0, a)), pg_axisym=True, pg_mesh=mesh, pg_inlet=1,
                pg_reactive=5, pg_gamma=0, pg_crosswind=0, pg_dt=1.0, pg_t_end=t_end,
                pg_csv_interval=5, pg_output_interval=600, pg_log_every=0, pg_j_x=J_X)
    csvs = pc.run_many(args, out, cases)

    BURST = 1.0e-4
    lines = ["| set | patch (um) | burst at (min) | gel at (min) | max T (nM) | front speed (um/min) | "
             "measured |", "|---|---|---|---|---|---|---|"]
    res = {}
    for ka in SETS:
        for p in WIDTHS:
            n = f"p{p*1e6:.0f}_{ka}"
            d = pc.read_csv(csvs[n])
            t = d["t"] / 60
            tc = pc.first_time(t, d["maxT"] >= BURST)
            tg = pc.first_time(t, d["gelArea"] > 0)
            v = pc.slope(t, d["gelHeight"] * 1e6, tc + 3, t[-1]) if np.isfinite(tc) else np.nan
            res[n] = (tc, d["maxT"].max() * 1e6, v)
            meas = "clots < 3 min" if p >= 100e-6 else "no clot"
            lines.append(f"| {ka} | {p*1e6:.0f} | {tc:.1f} | {tg:.1f} | {d['maxT'].max()*1e6:.1f} | "
                         f"{v:.1f} | {meas} |")
    lines.append("")
    for ka in SETS:
        clot = [p for p in WIDTHS if np.isfinite(res[f"p{p*1e6:.0f}_{ka}"][0])
                and res[f"p{p*1e6:.0f}_{ka}"][0] <= 3.0]
        none = [p for p in WIDTHS if not np.isfinite(res[f"p{p*1e6:.0f}_{ka}"][0])]
        thr = f"{max(none)*1e6:.0f}-{min(clot)*1e6:.0f} um" if clot and none else (
            "all clot" if clot and not none else ("none clots" if not clot else "mixed"))
        ok = (all(p >= 100e-6 for p in clot) and all(p <= 50e-6 for p in none)
              and len(clot) + len(none) == len(WIDTHS))
        lines.append(f"- {ka}: threshold {thr} (measured 50-100 um) -> {'PASA' if ok else 'FALLA'}")
    lines += [f"- j_X = {J_X:.1e} mol m^-2 s^-1 from Gamma_TF = {GAMMA_TF:.0e} mol/m^2 and kappa_TF = "
              f"{KAPPA_TF} 1/s (E1), not fitted.", "", "References:"] + [f"- {c}" for c in CITES]

    plt = pc.setup_matplotlib()
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    for j, ka in enumerate(SETS):
        for p in WIDTHS:
            d = pc.read_csv(csvs[f"p{p*1e6:.0f}_{ka}"])
            ax[j].plot(d["t"] / 60, d["maxT"] * 1e6, label=f"{p*1e6:.0f} um")
        ax[j].axvline(3, color="k", ls=":", lw=0.8)
        ax[j].set(xlabel="t (min)", ylabel="max thrombin (nM)", yscale="log", title=ka)
        ax[j].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "E4_patch.png"))
    pc.write_summary(out, "E4 - Patch-size threshold, no flow (Kastrup 2006)", lines)
    def fmt(t):
        return f"{t:.1f} min" if np.isfinite(t) else "no inicia"

    m = {p: res[f"p{p*1e6:.0f}_model"][0] for p in WIDTHS}
    sf = {p: res[f"p{p*1e6:.0f}_surface"][0] for p in WIDTHS}
    px.verdict(out, "E4 Kastrup 2006", [
        ("parches <= 50 um no inician", all(not np.isfinite(m[p]) for p in WIDTHS if p <= 50e-6),
         ", ".join(f"{p*1e6:.0f} um: {fmt(m[p])}" for p in WIDTHS if p <= 50e-6)),
        ("parches >= 100 um inician en < 3 min",
         all(np.isfinite(m[p]) and m[p] < 3.0 for p in WIDTHS if p >= 100e-6),
         ", ".join(f"{p*1e6:.0f} um: {fmt(m[p])}" for p in WIDTHS if p >= 100e-6)),
        ("frente de coagulacion a 25-35 um/min (no ajustado)",
         bool(np.isfinite(res["p400_model"][2]) and 25 <= res["p400_model"][2] <= 35),
         f"{res['p400_model'][2]:.1f} um/min (400 um)"),
    ], notes=["SURFACE, inicio: " + ", ".join(f"{p*1e6:.0f} um {fmt(sf[p])}" for p in WIDTHS)])


if __name__ == "__main__":
    main()

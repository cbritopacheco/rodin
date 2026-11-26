#!/usr/bin/env python3
"""E2 - Thrombodynamics: clot growth from a TF surface into quiescent plasma.

Experiment (HemaCore Thrombodynamics): platelet-free plasma, no stirring, a
TF-coated insert; the fibrin front grows away from the surface for 30 min.
  Tlag 0.6-1.5 min, Vi 38-56 um/min, Vst 20-29 um/min, CS30 800-1200 um
  (manufacturer ranges, cited by Gebekova et al. 2022); pooled normal plasma
  Vst 34.8 um/min (Shibeko et al. 2022). Spatial thrombin is a travelling
  wave that needs factor XI and VIII (Dashkevich et al. 2012).

Model: PTAGX kinetics (PTAGXBench), a 5 mm strip, Xa emitted by the bottom
wall at j_X = 1e-13 mol m^-2 s^-1, plasma held at the far end, no flow.
Cases: MODEL (k_a of ptagx_common), SURFACE (MODEL plus the surface
amplification of E5), LEGACY (PTAGX defaults), and a k_a sweep (sensitivity).

Observables as in the assay: Tlag (gel at the activator), Vi (front slope over
[Tlag+2, Tlag+6] min), Vst (slope over [15, 30] min), CS30 (front at 30 min).
NOTE: k_a of MODEL was chosen with this test (all four criteria at 0.035), so
the four ranges are a calibration check; Vi > Vst is not fitted.

Criteria (MODEL): Tlag, Vi, Vst and CS30 inside the reference ranges.
Cost: 7 runs of ~3 min.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

pc = px.pc
REF = {"Tlag": (0.6, 1.5), "Vi": (38.0, 56.0), "Vst": (20.0, 29.0), "CS30": (800.0, 1200.0)}
J_X = 1.0e-13
KAS = [0.05, 0.04, 0.03, 0.025, 0.02]
PHIS = []
CITES = [
    "Gebekova ZA, et al. Thrombodynamics test in assessing the risk of thrombus formation in "
    "patients with atrial fibrillation taking direct oral anticoagulants. Ration Pharmacother "
    "Cardiol 18 (2022) 544-552. doi:10.20996/1819-6446-2022-09-07",
    "Shibeko AM, et al. Systems biology approach for personalized hemostasis correction. "
    "J Pers Med 12 (2022) 1903. doi:10.3390/jpm12111903",
    "Dashkevich NM, et al. Thrombin activity propagates in space during blood coagulation as an "
    "excitation wave. Biophys J 103 (2012) 2233-2240",
]


def strip_mesh(out, quick):
    Ly = 3e-3 if quick else 5e-3
    hy = 10e-6 if quick else 5e-6
    return pc.make_mesh(os.path.join(out, "strip.mesh"), Lx=20e-6, Ly=Ly, nx=2, hy0=hy,
                        ry=1.0, hy_max=hy, patch_all=True, left=3, right=3, top=1)


def base_options(mesh, model=None):
    o = dict(px.PTAGX, pg_mesh=mesh, pg_inlet=1, pg_reactive=5, pg_gamma=0,
             pg_crosswind=0, pg_dt=1.0, pg_t_end=1800, pg_csv_interval=5,
             pg_output_interval=600, pg_log_every=0, pg_j_x=J_X)
    if model is not None:
        o.update(px.kinetics(model))
    return o


def observables(csv):
    d = pc.read_csv(csv)
    t = d["t"] / 60.0
    xf = d["gelHeight"] * 1e6
    lag = pc.first_time(t, d["gelArea"] > 0)
    return dict(Tlag=lag, Vi=pc.slope(t, xf, lag + 2, lag + 6), Vst=pc.slope(t, xf, 15, 30),
                CS30=np.interp(30.0, t, xf)), (t, xf)


def main():
    ap = pc.parser(__doc__)
    args = ap.parse_args()
    out = pc.results_dir(args, __file__)
    mesh = strip_mesh(out, args.quick)
    cases = {"model": base_options(mesh, px.MODEL), "surface": base_options(mesh, px.SURFACE),
             "legacy": base_options(mesh, px.LEGACY)}
    # gel without the soluble monomer (thrombin's cut gels where it happens),
    # and the feedback partly routed through Xa (phi, delayed by 1/k_dX)
    cases["nomono"] = dict(base_options(mesh, px.MODEL), pg_fibrin_monomer=False)
    # the Xa channel without the cofactor requirement (epsilon = 1, kappa = 3)
    cases["nocof"] = base_options(mesh, dict(px.MODEL, kappa=3.0, eps_x=1.0))
    for phi in PHIS:
        cases[f"phi{phi:g}"] = base_options(mesh, dict(px.MODEL, phi_i=phi))
    for ka in KAS:
        cases[f"ka{ka:g}"] = base_options(mesh, dict(px.MODEL, k_a=ka))
    csvs = pc.run_many(args, out, cases)

    obs = {n: observables(csvs[n]) for n in cases}
    plt = pc.setup_matplotlib()
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    for n, (o, (t, xf)) in obs.items():
        ax.plot(t, xf, "-" if n in ("model", "surface", "legacy", "nomono", "nocof") else ":", lw=2 if n == "model" else 1,
                label=n.replace("ka", "k_a="))
    ax.fill_between([29, 31], *REF["CS30"], color="k", alpha=0.15, lw=0, label="CS30 ref.")
    for v in REF["Vst"]:
        ax.plot([1, 30], [0, v * 29], "--", color="k", lw=0.6)
    ax.set(xlabel="t (min)", ylabel="clot front (um)", xlim=(0, 30), ylim=(0, 3000),
           title="Thrombodynamics (dashed: 20 and 29 um/min)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "E2_thrombodynamics.png"))

    lines = ["| case | Tlag (min) | Vi (um/min) | Vst (um/min) | CS30 (um) | in ref. (Tlag/Vi/Vst/CS30) |",
             "|---|---|---|---|---|---|"]
    for n, (o, _) in obs.items():
        ok = "/".join("si" if REF[k][0] <= o[k] <= REF[k][1] else "NO" for k in REF)
        lines.append(f"| {n} | {o['Tlag']:.2f} | {o['Vi']:.1f} | {o['Vst']:.1f} | {o['CS30']:.0f} | {ok} |")
    lines += ["", "References: " + "; ".join(f"{k} {v[0]}-{v[1]}" for k, v in REF.items()), ""]
    lines += [f"- {c}" for c in CITES]
    pc.write_summary(out, "E2 - Thrombodynamics, PTAGX kinetics", lines)

    m = obs["model"][0]
    px.verdict(out, "E2 Thrombodynamics", [
        (f"{k} en {REF[k][0]}-{REF[k][1]}", bool(REF[k][0] <= m[k] <= REF[k][1]), f"{m[k]:.2f}")
        for k in REF] + [("Vi > Vst (el frente decelera, como en el ensayo)", bool(m["Vi"] > m["Vst"]),
                          f"Vi {m['Vi']:.1f}, Vst {m['Vst']:.1f}")],
        notes=[f"{n}: " + ", ".join(f"{k} {obs[n][0][k]:.2f}" for k in REF)
               for n in ["nomono", "nocof", "surface", "legacy"]])


if __name__ == "__main__":
    main()

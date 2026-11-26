#!/usr/bin/env python3
"""E3 - Thrombodynamics under DOACs: peak vs trough in AF patients (ADAGE).

Experiment (Thrombodynamics-4D, fibrinography; AF patients >= 80 y on
reduced-dose DOACs; medians at peak (Tmax) and trough (Tmin), with measured
drug levels):

  drug          level peak/trough (ng/mL)  Tlag (min)   Vi (um/min)   CS30 (um)
  rivaroxaban   248 / 42                   2.0 / 1.2    60.6 / 70.4   1236 / 1691
  apixaban      170 / 74                   1.3 / 1.1    64.9 / 68.3   1362 / 1655
  dabigatran    171 / 39                   1.3 / 1.0    60.9 / 70.3   1436 / 1707

The comparison uses the peak/trough ratios, which remove the absolute offset
of an elderly, hyperfibrinogenaemic cohort (CS30 above the adult range at
trough).

Model: E2 set-up with the proposed parameters (ptagx_common.MODEL).
  anti-Xa (rivaroxaban, apixaban): prothrombinase is Xa.Va, so the drug scales
     both kappa and k_a by f (variant "kappa only": initiation alone);
  anti-IIa (dabigatran): faster effective thrombin removal, k_on A_0 x g
     (A_0 scaled, since only k_on A_0 enters the inhibition).
For each drug, f (or g) is chosen to reproduce the Tlag ratio; the Vi and
CS30 ratios are then predictions, compared with the measured ones.
Criteria: for each drug, at the dose that reproduces its Tlag ratio, the Vi and
CS30 ratios within +-0.05 of the measured ones (anti-Xa: the better of the two
variants).
Cost: 13 runs of ~3 min.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "E2_Thrombodynamics"))
import run_E2 as e2  # noqa: E402

pc = px.pc
ADAGE = {  # (Tlag, Vi, CS30) at peak and trough, medians
    "rivaroxaban": dict(peak=(2.0, 60.6, 1236), trough=(1.2, 70.4, 1691), cls="Xa"),
    "apixaban": dict(peak=(1.3, 64.9, 1362), trough=(1.1, 68.3, 1655), cls="Xa"),
    "dabigatran": dict(peak=(1.3, 60.9, 1436), trough=(1.0, 70.3, 1707), cls="IIa"),
}
F_BOTH = [0.97, 0.94, 0.9, 0.85, 0.8]
F_KAPPA = [0.9, 0.75, 0.6, 0.5, 0.25]
G_DTI = [1.05, 1.1, 1.2, 1.35]
CITES = [
    "ADAGE substudy: Fibrinography and thrombography (Thrombodynamics-4D) in atrial fibrillation: "
    "assessment of direct oral anticoagulants in geriatric patients aged 80 years and older "
    "(PMC12351382), Tables 1-3",
    "Gebekova ZA, et al. Ration Pharmacother Cardiol 18 (2022) 544-552",
]


def main():
    ap = pc.parser(__doc__)
    args = ap.parse_args()
    out = pc.results_dir(args, __file__)
    mesh = e2.strip_mesh(out, args.quick)
    base = e2.base_options(mesh, px.MODEL)
    KA = px.MODEL["k_a"]
    KAP = px.MODEL.get("kappa", px.K["kappa"])   # MODEL's kappa, not the default
    K = px.K
    cases = {"control": dict(base)}
    for f in F_BOTH:
        cases[f"Xa_both_{f:g}"] = dict(base, pg_kappa_x=f * KAP, pg_k_a=f * KA)
    for f in F_KAPPA:
        cases[f"Xa_kappa_{f:g}"] = dict(base, pg_kappa_x=f * KAP)
    for g in G_DTI:
        cases[f"IIa_{g:g}"] = dict(base, pg_A0=g * K["A0"])
    csvs = pc.run_many(args, out, cases)
    obs = {n: e2.observables(csvs[n])[0] for n in cases}
    c0 = obs["control"]

    def ratios(n):
        o = obs[n]
        return (o["Tlag"] / c0["Tlag"], o["Vi"] / c0["Vi"], o["CS30"] / c0["CS30"])

    lines = ["| case | Tlag (min) | Vi (um/min) | Vst (um/min) | CS30 (um) | Tlag/ctrl | Vi/ctrl | CS30/ctrl |",
             "|---|---|---|---|---|---|---|---|"]
    for n, o in obs.items():
        r = ratios(n)
        lines.append(f"| {n} | {o['Tlag']:.2f} | {o['Vi']:.1f} | {o['Vst']:.1f} | {o['CS30']:.0f} | "
                     f"{r[0]:.2f} | {r[1]:.2f} | {r[2]:.2f} |")

    lines += ["", "Peak/trough ratios: measured vs model at the dose that matches the Tlag ratio.", "",
              "| drug | family | measured Tlag/Vi/CS30 | model dose | model Tlag/Vi/CS30 | Vi, CS30 within 0.05 |",
              "|---|---|---|---|---|---|"]
    fams = {"Xa_both": [f"Xa_both_{f:g}" for f in F_BOTH],
            "Xa_kappa": [f"Xa_kappa_{f:g}" for f in F_KAPPA],
            "IIa": [f"IIa_{g:g}" for g in G_DTI]}
    passed = {}
    for drug, d in ADAGE.items():
        meas = tuple(p / t for p, t in zip(d["peak"], d["trough"]))
        for fam in (["Xa_both", "Xa_kappa"] if d["cls"] == "Xa" else ["IIa"]):
            names = ["control"] + fams[fam]
            R = np.array([ratios(n) for n in names])
            # interpolate along the dose sequence on the Tlag ratio (monotone if the model is sane)
            tl = R[:, 0]
            if np.all(np.diff(tl) > 0) and tl[0] <= meas[0] <= tl[-1]:
                k = np.searchsorted(tl, meas[0])
                w = (meas[0] - tl[k - 1]) / (tl[k] - tl[k - 1])
                mod = R[k - 1] * (1 - w) + R[k] * w
                dose = f"between {names[k-1]} and {names[k]}"
                ok = abs(mod[1] - meas[1]) <= 0.05 and abs(mod[2] - meas[2]) <= 0.05
                lines.append(f"| {drug} | {fam} | {meas[0]:.2f}/{meas[1]:.2f}/{meas[2]:.2f} | {dose} | "
                             f"{mod[0]:.2f}/{mod[1]:.2f}/{mod[2]:.2f} | {'si' if ok else 'NO'} |")
                detail = (f"medido Tlag/Vi/CS30 {meas[0]:.2f}/{meas[1]:.2f}/{meas[2]:.2f}, "
                          f"modelo ({fam}) {mod[0]:.2f}/{mod[1]:.2f}/{mod[2]:.2f}")
            else:
                ok = False
                detail = f"{fam}: el modelo no alcanza el cociente de Tlag {meas[0]:.2f} ({tl.min():.2f}-{tl.max():.2f})"
                lines.append(f"| {drug} | {fam} | {meas[0]:.2f}/{meas[1]:.2f}/{meas[2]:.2f} | "
                             f"Tlag ratio not reached (model {tl.min():.2f}-{tl.max():.2f}) | - | NO |")
            # report the family that passes, else the one that reaches the measured Tlag ratio
            reached = "no alcanza" not in detail
            if ok or drug not in passed or (reached and not passed[drug][2] and not passed[drug][0]):
                passed[drug] = (ok, detail, reached)

    plt = pc.setup_matplotlib()
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
    for fam, names, mk in [("anti-Xa (kappa and k_a)", fams["Xa_both"], "o-"),
                           ("anti-Xa (kappa only)", fams["Xa_kappa"], "s--"),
                           ("anti-IIa (k_on A_0)", fams["IIa"], "^-")]:
        R = np.array([ratios(n) for n in ["control"] + names])
        for j in range(3):
            ax[j].plot(R[:, 0] if j else np.arange(len(R)), R[:, j], mk, label=fam)
    for drug, d in ADAGE.items():
        meas = [p / t for p, t in zip(d["peak"], d["trough"])]
        for j in (1, 2):
            ax[j].plot(meas[0], meas[j], "k*", ms=9)
            ax[j].annotate(drug[:4], (meas[0], meas[j]), xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax[0].set(xlabel="dose index", ylabel="Tlag / control", title="Lag")
    ax[1].set(xlabel="Tlag ratio", ylabel="Vi ratio", title="Vi vs Tlag (stars: ADAGE peak/trough)")
    ax[2].set(xlabel="Tlag ratio", ylabel="CS30 ratio", title="CS30 vs Tlag")
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "E3_doac.png"))
    lines += ["", "References:"] + [f"- {c}" for c in CITES]
    pc.write_summary(out, "E3 - Thrombodynamics under DOACs, PTAGX kinetics (MODEL)", lines)
    px.verdict(out, "E3 DOAC (ADAGE)", [
        (f"{drug}: Vi y CS30 pico/valle dentro de 0.05 al Tlag medido", ok, det)
        for drug, (ok, det, _) in passed.items()])


if __name__ == "__main__":
    main()

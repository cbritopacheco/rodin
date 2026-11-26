#!/usr/bin/env python3
"""E1 - Calibrated automated thrombogram (CAT), 0D, PTAGX kinetics.

Experiment: platelet-poor plasma, 5 pM tissue factor, 4 uM phospholipid,
well mixed, no flow. Reference intervals of healthy adults (ELSA-Brasil
2022): lag 1.4-3.7 min, peak 256-479 nM, ETP 1414-2658 nM min.

Model: the nodal step of LeftAtrium2DPTAGX (ptagx_common.advance), no
transport. Tissue factor is a volumetric Xa source s_X = kappa_TF [TF].
kappa_TF (Xa made per TF-VIIa per second at plasma factor X) belongs to the
assay, not to the model; it sets the lag only and is taken as 0.5 1/s
(kcat ~ 2-5 1/s times [X]/(Km + [X]) ~ 0.1-0.4).

Checks:
  E1a  5 pM: lag, peak, ETP inside all three intervals.
  E1b  1 pM: longer lag and lower peak than at 5 pM (Butenas 1999; Hemker 2003).
  E1c  k_a sweep: which k_a the CAT admits (the PTAGX default is 0.10).
Cost: seconds, no executable.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

pc = px.pc
KAPPA_TF = 0.5
PM = 1.0e-9
REF = dict(lag=(1.4, 3.7), peak=(256.0, 479.0), etp=(1414.0, 2658.0))
CITES = [
    "ELSA-Brasil. Performance and reference intervals of thrombin generation test. "
    "Sao Paulo Med J (2022)",
    "Hemker HC, et al. Calibrated automated thrombin generation measurement in clotting plasma. "
    "Pathophysiol Haemost Thromb 33 (2003) 4-15",
    "Butenas S, van 't Veer C, Mann KG. 'Normal' thrombin generation. Blood 94 (1999) 2169-2178",
]


def thrombogram(tf_pM, dt=0.1, t_end=3600.0, **over):
    c = dict(px.MODEL and px.consts(px.MODEL))
    c.update(over)
    s = (c["P0"], 0.0, c["A0"], c["G0"], 0.0, 0.0)
    n = int(round(t_end / dt))
    t = np.arange(n + 1) * dt
    T = np.zeros(n + 1)
    for i in range(1, n + 1):
        s = px.advance(s, dt, c, sX=KAPPA_TF * tf_pM * PM)
        T[i] = s[1]
    return t, T


def metrics(t, T):
    Tn = T * 1e6
    peak = Tn.max()
    lag = pc.first_time(t, Tn >= peak / 6.0) / 60.0
    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return lag, peak, t[int(Tn.argmax())] / 60.0, trap(Tn, t) / 60.0


def main():
    ap = pc.parser(__doc__)
    args = ap.parse_args()
    out = pc.results_dir(args, __file__)
    plt = pc.setup_matplotlib()
    rows = ["| case | lag (min) | peak (nM) | ttPeak (min) | ETP (nM min) | in ref. |",
            "|---|---|---|---|---|---|"]
    curves, res = {}, {}

    def add(name, tf, **kw):
        t, T = thrombogram(tf, **kw)
        lag, peak, ttp, etp = metrics(t, T)
        ok = (REF["lag"][0] <= lag <= REF["lag"][1], REF["peak"][0] <= peak <= REF["peak"][1],
              REF["etp"][0] <= etp <= REF["etp"][1])
        curves[name] = (t / 60, T * 1e6)
        res[name] = (lag, peak, etp)
        rows.append(f"| {name} | {lag:.2f} | {peak:.0f} | {ttp:.2f} | {etp:.0f} | "
                    f"{'/'.join('si' if o else 'NO' for o in ok)} |")

    add("5 pM (MODEL)", 5.0)
    add("1 pM (MODEL)", 1.0)
    add("5 pM, LEGACY", 5.0, **px.consts(px.LEGACY))
    # without the cofactor requirement (epsilon = 1): the old kinetics
    add("5 pM, eps 1, kappa 3", 5.0, eps_x=1.0, kappa=3.0)
    for eps in (0.3, 0.03):
        add(f"5 pM, eps {eps:g}", 5.0, eps_x=eps)
    for ka in (0.05, 0.04, 0.03, 0.025):
        add(f"5 pM, k_a {ka:g}", 5.0, k_a=ka)

    l5, p5, _ = res["5 pM (MODEL)"]
    l1, p1, _ = res["1 pM (MODEL)"]
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    for n, (t, T) in curves.items():
        ax.plot(t, T, label=n)
    ax.axhspan(*REF["peak"], color="k", alpha=0.08, lw=0, label="peak ref.")
    ax.axvspan(*REF["lag"], color="C1", alpha=0.08, lw=0, label="lag ref.")
    ax.set(xlabel="t (min)", ylabel="thrombin (nM)", xlim=(0, 30), title="CAT, PTAGX kinetics")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "E1_cat.png"))

    rows += ["",
             "E1a (5 pM in all three intervals, MODEL): "
             + ("PASA" if all(REF[k][0] <= v <= REF[k][1] for k, v in
                              zip(("lag", "peak", "etp"), res["5 pM (MODEL)"])) else "FALLA"),
             "E1c cases at 5 pM inside all three intervals: "
             + (", ".join(n for n, v in res.items() if n.startswith("5 pM")
                          and all(REF[k][0] <= x <= REF[k][1] for k, x in zip(("lag", "peak", "etp"), v)))
                or "none"),
             f"E1b (1 pM: lag longer, peak lower): lag {l1:.2f} vs {l5:.2f} min, "
             f"peak {p1:.0f} vs {p5:.0f} nM -> {'PASA' if (l1 > l5 and p1 < p5) else 'FALLA'}"
             + (" (peak independent of TF: known limit of one amplification term)"
                if abs(p1 - p5) / p5 < 0.05 else ""),
             "", "References:"] + [f"- {c}" for c in CITES]
    pc.write_summary(out, "E1 - CAT thrombogram (0D), PTAGX kinetics", rows)
    ka = px.MODEL["k_a"]
    t, T = thrombogram(5.0)
    lag, peak, _, etp = metrics(t, T)
    px.verdict(out, "E1 CAT", [
        (f"lag en {REF['lag'][0]}-{REF['lag'][1]} min (k_a {ka:g})", REF["lag"][0] <= lag <= REF["lag"][1], f"{lag:.2f}"),
        (f"pico en {REF['peak'][0]:.0f}-{REF['peak'][1]:.0f} nM", REF["peak"][0] <= peak <= REF["peak"][1], f"{peak:.0f}"),
        (f"ETP en {REF['etp'][0]:.0f}-{REF['etp'][1]:.0f} nM min", REF["etp"][0] <= etp <= REF["etp"][1], f"{etp:.0f}"),
    ], notes=[f"pico a 1 pM frente a 5 pM: {p1:.0f} vs {p5:.0f} nM "
              f"({'depende del FT' if abs(p1 - p5) / p5 >= 0.05 else 'no depende del FT'})"])


if __name__ == "__main__":
    main()

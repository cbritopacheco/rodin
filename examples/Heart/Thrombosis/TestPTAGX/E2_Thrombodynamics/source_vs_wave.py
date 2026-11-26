#!/usr/bin/env python3
"""Is the clot front an autonomous wave, or is the activator still pushing it?

In the assay the thrombin front is an excitation wave (Dashkevich et al. 2012):
once it has started, its speed is a property of the plasma, not of the
activator. The test here is direct: the same kinetics with three activator
doses j_X. If Vi, Vst and CS30 move with j_X, the front is being driven by the
source and the model has no autonomous wave; if they do not, it has one.

This is the diagnostic behind the E3 result: a front that the source still
drives is a front that every perturbation of the source slows down, which is
why the model loses far more clot-growth speed than the ADAGE patients do for
the same delay in the lag.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_E2 as e2  # noqa: E402

pc = px.pc
FACTORS = [0.5, 1.0, 2.0, 4.0]


def main():
    ap = pc.parser(__doc__)
    args = ap.parse_args()
    out = os.path.join(pc.results_dir(args, __file__), "source")
    os.makedirs(out, exist_ok=True)
    mesh = e2.strip_mesh(out, args.quick)
    cases = {f"jx{f:g}": dict(e2.base_options(mesh, px.MODEL), pg_j_x=f * e2.J_X)
             for f in FACTORS}
    csvs = pc.run_many(args, out, cases)

    lines = ["| j_X / j_X0 | Tlag (min) | Vi (um/min) | Vst (um/min) | CS30 (um) |",
             "|---|---|---|---|---|"]
    ref = None
    for f in FACTORS:
        d = pc.read_csv(csvs[f"jx{f:g}"])
        t, xf = d["t"] / 60, d["gelHeight"] * 1e6
        lag = pc.first_time(t, d["gelArea"] > 0)
        o = (lag, pc.slope(t, xf, lag + 2, lag + 6), pc.slope(t, xf, 15, 30),
             float(np.interp(30.0, t, xf)))
        if f == 1.0:
            ref = o
        lines.append(f"| {f:g} | {o[0]:.2f} | {o[1]:.1f} | {o[2]:.1f} | {o[3]:.0f} |")
    spread = max(abs(float(np.interp(30.0, pc.read_csv(csvs[f'jx{f:g}'])["t"] / 60,
                                     pc.read_csv(csvs[f'jx{f:g}'])["gelHeight"] * 1e6)) / ref[3] - 1.0)
                 for f in FACTORS)
    lines += ["", f"CS30 changes by up to {100*spread:.0f} % over an eightfold change of the "
              "activator dose, with the same kinetics.",
              "An autonomous wave would not move; a source-driven front does.",
              "", f"k_a/(k_on A_0) = {px.MODEL['k_a'] / (px.K['k_on'] * px.K['A0']):.2f}: "
              "the volume loop is only just above the inhibition threshold, which is why."]
    pc.write_summary(out, "Autonomous wave or source-driven front?", lines)
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()

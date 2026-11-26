#!/usr/bin/env python3
"""Shen et al. 2008 as a calibration of the tissue-factor dose j_X.

The paper does not report the tissue-factor density of its patches, so what
E5 constrains is not the kinetics but the product of the dose and the
kinetics. Calibrating sigma (surface amplification) on it, as SURFACE does,
spends a free parameter of the kinetics on an unmeasured property of the
experiment; calibrating j_X spends it on the unmeasured property itself.

Scan of j_X with MODEL (no surface amplification) on the 200 um patch:
  5 1/s  must initiate (thrombin > 100 nM) within 300 s (measured ~210 s),
  40 1/s must not initiate within 300 s.
If some j_X does both, the model reproduces the shear threshold with one dose
and the question becomes how that dose compares with the only measured one,
Kastrup et al. 2006: j_X = Gamma_TF kappa_TF = 4e-13 mol m^-2 s^-1 (E4).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_E5 as e5  # noqa: E402

pc = px.pc
DOSES = [4e-13, 1.2e-12, 4e-12, 6e-12, 8e-12, 1.0e-11, 1.2e-11, 4e-11, 1.2e-10]
SHEARS = [5.0, 40.0]


def main():
    ap = pc.parser(__doc__)
    args = ap.parse_args()
    out = os.path.join(pc.results_dir(args, __file__), "dose")
    os.makedirs(out, exist_ok=True)
    mesh = e5.channel_mesh(out)

    cases = {f"jx{j:.1e}_g{g:g}": e5.case_options(mesh, g, j, px.MODEL)
             for j in DOSES for g in SHEARS}
    csvs = pc.run_many(args, out, cases)

    t = {n: pc.first_time(d["t"], d["maxT"] >= e5.BURST)
         for n, d in ((n, pc.read_csv(c)) for n, c in csvs.items())}
    lines = ["Measured: initiation at 5 1/s at ~210 s. A model that only switches "
             "from 'never' to 'at once' as the dose rises has no dose that reproduces "
             "that number, however well it reproduces the threshold itself.", "",
             "| j_X (mol m^-2 s^-1) | x Kastrup | 5 1/s (s) | 40 1/s (s) | both criteria |",
             "|---|---|---|---|---|"]
    ok = []
    for j in DOSES:
        t5, t40 = t[f"jx{j:.1e}_g5"], t[f"jx{j:.1e}_g40"]
        good = bool(np.isfinite(t5) and t5 <= 300 and (not np.isfinite(t40) or t40 > 300))
        if good:
            ok.append(j)
        lines.append(f"| {j:.1e} | {j/4e-13:.0f} | {t5:.0f} | {t40:.0f} | "
                     f"{'si' if good else 'no'} |")
    if ok:
        lines += ["", f"window: j_X in [{min(ok):.1e}, {max(ok):.1e}] mol m^-2 s^-1, "
                  f"{min(ok)/4e-13:.0f}-{max(ok)/4e-13:.0f} times the measured dose of "
                  "Kastrup 2006 (E4)"]
    else:
        lines += ["", "no dose in the scan satisfies both criteria"]
    pc.write_summary(out, "E5 - Tissue-factor dose that reproduces the shear threshold", lines)
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()

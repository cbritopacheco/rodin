#!/usr/bin/env python3
"""Largest sigma that keeps the 50 um patch quiet without flow (E4), to compare
with the smallest sigma that lets the 200 um patch initiate at 5 1/s (E5:
1.92e-3 m). If the first is below the second, no sigma passes both tests.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

pc = px.pc
SIGMAS = [3e-5, 1e-4, 3e-4, 1e-3]
P = 50e-6
LX, LY = 2.5e-3, 1.2e-3
J_X = 4.0e-13


def main():
    ap = pc.parser(__doc__)
    args = ap.parse_args()
    out = os.path.join(pc.results_dir(args, __file__), "sigma_scan")
    os.makedirs(out, exist_ok=True)
    x0, x1 = 0.5 * (LX - P), 0.5 * (LX + P)
    mesh = pc.make_mesh(os.path.join(out, "chamber_50um.mesh"), Lx=LX, Ly=LY, x0=x0, x1=x1,
                        margin=120e-6, hx_patch=8e-6, hx_far=80e-6, rx=1.12, hy0=8e-6, ry=1.1,
                        hy_max=80e-6, left=3, right=3, bottom=3, top=1)
    cases = {f"sigma{s:.0e}": dict(
        px.kinetics(dict(px.MODEL, sigma_s=s), patch=(x0, x1)), pg_mesh=mesh, pg_inlet=1,
        pg_reactive=5, pg_gamma=0, pg_crosswind=0, pg_dt=1.0, pg_t_end=1500, pg_csv_interval=5,
        pg_output_interval=0, pg_log_every=0, pg_j_x=J_X, pg_stop_thrombin=1e-4) for s in SIGMAS}
    csvs = pc.run_many(args, out, cases)
    lines = ["| sigma (m) | 50 um patch, no flow: initiation (min) |", "|---|---|"]
    quiet = []
    for s in SIGMAS:
        d = pc.read_csv(csvs[f"sigma{s:.0e}"])
        t = pc.first_time(d["t"] / 60, d["maxT"] >= 1e-4)
        lines.append(f"| {s:.0e} | {t:.1f} |")
        if not np.isfinite(t):
            quiet.append(s)
    smax = max(quiet) if quiet else 0.0
    lines += ["", f"largest sigma with no initiation (of those tried): {smax:.0e} m; "
              f"E5 needs >= 1.92e-3 m -> {'compatible' if smax >= 1.92e-3 else 'INCOMPATIBLE'}"]
    pc.write_summary(out, "E4 - sigma that the 50 um patch tolerates", lines)


if __name__ == "__main__":
    main()

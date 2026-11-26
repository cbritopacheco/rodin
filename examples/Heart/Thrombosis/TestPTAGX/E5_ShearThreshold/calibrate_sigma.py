#!/usr/bin/env python3
"""Calibration of sigma (surface amplification) on Shen et al. 2008.

Target: on the 200 um TF patch at 5 1/s, initiation (thrombin > 100 nM) at
~210 s. Bisection on log(sigma) at k_a = MODEL k_a. The result goes into
ptagx_common.SURFACE["sigma_s"]; E4 (no flow) and the 40 1/s case of E5 then
test it without further tuning.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_E5 as e5  # noqa: E402

pc = px.pc
TARGET = 210.0


def main():
    ap = pc.parser(__doc__)
    ap.add_argument("--lo", type=float, default=1e-5)
    ap.add_argument("--hi", type=float, default=1e-2)
    ap.add_argument("--iterations", type=int, default=8)
    args = ap.parse_args()
    out = pc.results_dir(args, __file__)
    out = os.path.join(out, "calibration")
    os.makedirs(out, exist_ok=True)
    mesh = e5.channel_mesh(out)

    def t_init(sigma, gamma=5.0):
        name = f"sigma{sigma:.4e}_g{gamma:g}"
        opts = e5.case_options(mesh, gamma, e5.J_X, dict(px.MODEL, sigma_s=sigma))
        d = pc.read_csv(pc.run_case(args, out, name, opts))
        return pc.first_time(d["t"], d["maxT"] >= e5.BURST)

    lo, hi = math.log(args.lo), math.log(args.hi)
    hist = []
    t_hi = t_init(math.exp(hi))
    hist.append((math.exp(hi), t_hi))
    if not math.isfinite(t_hi):
        print(f"no initiation even at sigma = {args.hi:g}")
        return
    for _ in range(args.iterations):
        mid = 0.5 * (lo + hi)
        t = t_init(math.exp(mid))
        hist.append((math.exp(mid), t))
        print(f"  sigma {math.exp(mid):.3e} m -> {t:.0f} s", flush=True)
        # larger sigma initiates earlier
        if math.isfinite(t) and t <= TARGET:
            hi = mid
        else:
            lo = mid
    sigma = math.exp(hi)
    lines = ["| sigma (m) | initiation at 5 1/s (s) |", "|---|---|"]
    lines += [f"| {s:.3e} | {t:.0f} |" for s, t in sorted(hist)]
    lines += ["", f"sigma = {sigma:.3e} m (first value with initiation <= {TARGET:.0f} s)",
              f"theta at the wall with delta_s = {px.MODEL['delta_s']*1e6:.0f} um: "
              f"{sigma / px.MODEL['delta_s']:.1f}"]
    pc.write_summary(out, "Calibration of sigma on Shen 2008 (5 1/s, ~210 s)", lines)


if __name__ == "__main__":
    main()

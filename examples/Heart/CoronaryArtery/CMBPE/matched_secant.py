#!/usr/bin/env python3
"""Constant outlet resistance matched to the mean inflow of R_mu.

The rest closure (constant at the resting operating point) is scaled by s
until its mean |Q_in| over the last cycle equals that of the R_mu run of the
same group. Secant iteration on s, starting from the s = 1 run of the group;
each iterate is a full coupled run in <group>/rest_matched_it<k>, and the
converged one is linked as <group>/rest_matched (with scale.txt).

    matched_secant.py --group DIR --run "MPIRUN -n 10 EXE <common args>" \\
        --period 0.85 --cycles 3
"""
import argparse
import csv
import os
import subprocess
import sys


def mean_inflow(run_dir, period, cycles):
    rows = list(csv.DictReader(open(os.path.join(run_dir, "CoronaryArtery.csv"))))
    t0 = (cycles - 1) * period
    q = [abs(float(r["CoronaryInletFlux"])) for r in rows if float(r["t"]) >= t0 - 1e-12]
    if not q:
        raise RuntimeError("no last-cycle samples in %s" % run_dir)
    return sum(q) / len(q)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    ap.add_argument("--run", required=True, help="command up to (excluding) closure and prefix")
    ap.add_argument("--period", type=float, default=0.85)
    ap.add_argument("--cycles", type=int, default=3)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--maxit", type=int, default=4)
    a = ap.parse_args()

    target = mean_inflow(os.path.join(a.group, "rmu"), a.period, a.cycles)
    s0, m0 = 1.0, mean_inflow(os.path.join(a.group, "rest"), a.period, a.cycles)
    print("target mean|Q| = %.6g, rest(s=1) = %.6g" % (target, m0), flush=True)
    if abs(m0 - target) <= a.tol * target:
        s1, best = 1.0, os.path.join(a.group, "rest")
    else:
        # mean|Q| falls roughly like 1/s: first guess from that, then secant
        s1 = s0 * m0 / target
        best = None
        for it in range(1, a.maxit + 1):
            d = os.path.join(a.group, "rest_matched_it%d" % it)
            os.makedirs(d, exist_ok=True)
            done = os.path.join(d, ".done")
            if not os.path.exists(done):
                cmd = "%s -coronary_outlet_closure rest -coronary_outlet_constant_scale %.12g " \
                      "-coronary_output_prefix %s" % (a.run, s1, d)
                with open(os.path.join(d, "log.txt"), "w") as log:
                    rc = subprocess.call(["bash", "-c", cmd], stdout=log, stderr=subprocess.STDOUT)
                if rc != 0:
                    sys.exit("run failed in %s" % d)
                open(done, "w").close()
            m1 = mean_inflow(d, a.period, a.cycles)
            print("it %d: s = %.8g  mean|Q| = %.6g  (%.3e relative)"
                  % (it, s1, m1, (m1 - target) / target), flush=True)
            best = d
            if abs(m1 - target) <= a.tol * target:
                break
            s2 = s1 - (m1 - target) * (s1 - s0) / (m1 - m0)
            s0, m0, s1 = s1, m1, s2

    link = os.path.join(a.group, "rest_matched")
    if os.path.islink(link):
        os.unlink(link)
    if not os.path.exists(link):
        os.symlink(os.path.relpath(best, a.group), link)
    with open(os.path.join(best, "scale.txt"), "w") as f:
        f.write("%.12g\n" % s1)
    print("matched: %s (s = %.8g)" % (best, s1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Run the PTAGX validation battery and print each verdict as soon as a test ends.

  python3 run_all.py                       # E1-E5 and Q1 (all of it)
  python3 run_all.py --only E1,E2          # some tests
  python3 run_all.py --skip-q1             # the experimental tests only (~1 h on 2 cores)
  python3 run_all.py --only Q1 --q1-kinetics qureshi,ptagx,ptagx_legacy

Each test writes <test>/results/verdict.md (Q1: one per kinetics in
Q1_Qureshi2021/results/<kinetics>/verdict.md); this script prints it right
after the test and a table of all verdicts at the end. The model parameters
are in ptagx_common.py (MODEL = proposed, SURFACE = MODEL + sigma, LEGACY =
PTAGX defaults); verdicts are on MODEL.

Environment: PTAGP_EXE (PTAGXBench executable), PTAGX_EXE (LeftAtrium2DPTAGX),
PTAGP_MPIRUN (e.g. "mpirun --allow-run-as-root -n 1" in a container).
"""
import argparse
import glob
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = {
    "E1": ["E1_CAT/run_E1.py"],
    "E2": ["E2_Thrombodynamics/run_E2.py"],
    "E3": ["E3_DOAC/run_E3.py"],
    "E4": ["E4_PatchThreshold/run_E4.py"],
    "E5": ["E5_ShearThreshold/run_E5.py"],
    "Q1": ["Q1_Qureshi2021/run_Q1.py --stage calib", "Q1_Qureshi2021/run_Q1.py --stage coag",
           "Q1_Qureshi2021/run_Q1.py --stage report"],
}


def verdict_files(test):
    if test == "Q1":
        return sorted(glob.glob(os.path.join(HERE, "Q1_Qureshi2021", "results", "*", "verdict.md")))
    d = os.path.dirname(TESTS[test][0])
    return [os.path.join(HERE, d, "results", "verdict.md")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="E1,E2,E3,E4,E5,Q1")
    ap.add_argument("--skip-q1", action="store_true")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--q1-kinetics", default="qureshi,ptagx")
    ap.add_argument("--q1-shapes", default="CW,WS,BR,CA")
    ap.add_argument("--force", action="store_true", help="rerun cases that already exist")
    a = ap.parse_args()
    tests = [t for t in a.only.split(",") if not (a.skip_q1 and t == "Q1")]
    summary = []
    for t in tests:
        t0 = time.time()
        print(f"\n===== {t} =====", flush=True)
        for step in TESTS[t]:
            cmd = [sys.executable] + [os.path.join(HERE, step.split()[0])] + step.split()[1:]
            cmd += ["--jobs", str(a.jobs)]
            if t == "Q1":
                cmd += ["--kinetics", a.q1_kinetics, "--shapes", a.q1_shapes]
            if a.force:
                cmd.append("--force")
            log = os.path.join(HERE, os.path.dirname(step.split()[0]), "run.log")
            with open(log, "a") as f:
                r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
            if r.returncode != 0:
                print(f"  {step}: ERROR (see {log})", flush=True)
        for vf in verdict_files(t):
            if os.path.exists(vf):
                text = open(vf).read()
                print(text, flush=True)
                summary.append((t, text.splitlines()[0].lstrip("# ")))
            else:
                print(f"  (no verdict: {vf})", flush=True)
                summary.append((t, "sin veredicto"))
        print(f"  [{t}: {(time.time() - t0) / 60:.1f} min]", flush=True)
    print("\n===== Resumen =====")
    for t, line in summary:
        print(f"  {line}")


if __name__ == "__main__":
    main()

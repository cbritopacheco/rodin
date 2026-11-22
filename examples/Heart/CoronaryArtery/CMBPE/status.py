#!/usr/bin/env python3
"""State of every coronary run, finished or not.

For each run directory under OUT (one holding log.txt):
  state      done / running / STALLED (no log update for 15 min) / ERROR
  progress   last step and simulated time, per cent of the run
  s/step     mean wall time of the last 20 steps, and the remaining hours
  |Q| cycles mean |Q_in| (mL/s) of every COMPLETE cycle, and the change
             between the last two (periodicity)
  vs rmu     change of mean |Q_in| against the rmu run of the same group,
             on the last cycle that both have completed

  python3 status.py ~/Results/cmbpe/coronary [--period 0.85]
"""
import argparse
import csv
import os
import re
import time


def tail(path, nbytes=400000):
    with open(path, "rb") as f:
        f.seek(0, 2)
        f.seek(max(0, f.tell() - nbytes))
        return f.read().decode(errors="ignore")


def cycles(path, T):
    """Mean |Q_in| (mL/s) of each complete cycle in CoronaryArtery.csv."""
    if not os.path.isfile(path):
        return []
    t, q = [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                t.append(float(r["t"]))
                q.append(abs(float(r["CoronaryInletFlux"])) * 1e6)
            except (KeyError, ValueError):
                pass
    if not t:
        return []
    dt = t[1] - t[0] if len(t) > 1 else 0.0
    out, k = [], 0
    while (k + 1) * T <= t[-1] + dt + 1e-9:
        v = [x for s, x in zip(t, q) if k * T - 1e-12 <= s < (k + 1) * T - 1e-12]
        if len(v) > 1:
            out.append(sum(v) / len(v))
        k += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--period", type=float, default=0.85)
    a = ap.parse_args()
    out = os.path.abspath(os.path.expanduser(a.out))

    runs, seen = {}, set()
    for root, dirs, files in os.walk(out, followlinks=True):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        if "log.txt" not in files:
            continue
        real = os.path.realpath(root)
        if real in seen:
            continue
        seen.add(real)
        name = os.path.relpath(root, out)
        log = os.path.join(root, "log.txt")
        txt = tail(log)
        steps = re.findall(r"Step (\d+)\s+t = ([0-9.eE+-]+) s / ([0-9.eE+-]+) s", txt)
        times = [float(x) for x in re.findall(r"step = \d+\s+total = ([0-9.eE+-]+)", txt)][-20:]
        age = time.time() - os.path.getmtime(log)
        if os.path.exists(os.path.join(root, ".done")):
            state = "done"
        elif re.search(r"\b(Fatal|FATAL|Abort|Segmentation|terminate called|MPI_ABORT)", txt):
            state = "ERROR"
        elif age > 900:
            state = "STALLED"
        else:
            state = "running"
        r = dict(state=state, step="", pct="", sps="", eta="",
                 cyc=cycles(os.path.join(root, "CoronaryArtery.csv"), a.period))
        if steps:
            n, t, tf = int(steps[-1][0]), float(steps[-1][1]), float(steps[-1][2])
            r["step"] = "%d (t=%.3f)" % (n, t)
            r["pct"] = 100.0 * t / tf if tf else 0.0
            if times:
                sps = sum(times) / len(times)
                dt = t / max(n - 1, 1)
                r["sps"] = sps
                r["eta"] = max(tf - t, 0.0) / dt * sps / 3600.0 if dt > 0 else ""
        runs[name] = r

    def ref(name):
        g, c = os.path.split(name)
        if c == "rmu":
            return None
        return runs.get(os.path.join(g, "rmu")) or \
            runs.get(os.path.join(os.path.dirname(g), "rmu"))  # rest_matched_itK

    order = sorted(runs, key=lambda n: (n.split("/")[0], n))
    print("%-38s %-8s %-15s %5s %7s %6s  %-34s %8s %9s" %
          ("run", "state", "step", "%", "s/step", "eta h", "|Q| per complete cycle (mL/s)",
           "dQ cyc%", "vs rmu %"))
    for name in order:
        r = runs[name]
        cyc = r["cyc"]
        per = "%+.2f" % (100 * (cyc[-1] - cyc[-2]) / cyc[-2]) if len(cyc) > 1 else ""
        vs = ""
        rr = ref(name)
        if rr and cyc and rr["cyc"]:
            k = min(len(cyc), len(rr["cyc"]))
            vs = "%+.2f (c%d)" % (100 * (cyc[k - 1] - rr["cyc"][k - 1]) / rr["cyc"][k - 1], k)
        print("%-38s %-8s %-15s %5s %7s %6s  %-34s %8s %9s" % (
            name[:38], r["state"], r["step"],
            "%.0f" % r["pct"] if r["pct"] != "" else "",
            "%.1f" % r["sps"] if r["sps"] != "" else "",
            "%.1f" % r["eta"] if r["eta"] not in ("", None) else "",
            " ".join("%.4f" % x for x in cyc)[:34], per, vs))
    bad = [n for n in order if runs[n]["state"] in ("ERROR", "STALLED")]
    if bad:
        print("\ncheck the log of:", ", ".join(bad))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Errors of the reduced closures against the full solution (RmuAxi).

Reads <out>/s<scale>_g<gamma>/<state>/<case>/RmuAxi.csv, keeps the last cycle
and compares every truncated case with the full one of the same scale, shear
rate and blood. All errors are relative to the cycle mean of |Q_in| of the
full solution (or to the corresponding full quantity), in per cent:

  dQ_signed    mean(Q_c - Q_f)                      net transport
  dQ_abs       mean(|Q_c|) - mean(|Q_f|)            transport based on |Q|
  dQ_fwd/rev   part of the signed difference accumulated while the full
               flow is forward / reversed (they add up to dQ_signed)
  dQ_rms       rms(Q_c - Q_f)                        waveform
  dQ_max       max |Q_c - Q_f| / max |Q_f|
  dE           cycle-integrated dissipation in the common region x <= x_T
  dWSS_parent  cycle mean of |tau_w| on the parent wall
  dWSS_stub    cycle mean of |tau_w| on the stub wall

  periodicity  change of mean |Q_in| between the last two cycles (own run;
               full_periodicity for the reference)

Writes <out>/summary.csv and, if matplotlib is available, one figure of
Q_in(t) per (scale, gamma, state).
"""
import argparse
import csv
import math
import os
import sys


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]} if rows else None


def last_cycle(data, cycles):
    cyc = data["cycle"]
    last = cycles - 1 if cycles else max(cyc)
    idx = [i for i, c in enumerate(cyc) if int(round(c)) == last]
    return {k: [v[i] for i in idx] for k, v in data.items()}


def mean(x):
    return sum(x) / len(x)


def periodicity(data):
    """|mean|Q_in| of the last cycle - that of the previous one| / the last, %."""
    cyc = [int(round(c)) for c in data["cycle"]]
    last = max(cyc)
    if last < 1:
        return float("nan")
    a = [abs(q) for q, c in zip(data["qIn"], cyc) if c == last]
    b = [abs(q) for q, c in zip(data["qIn"], cyc) if c == last - 1]
    return 100.0 * abs(mean(a) - mean(b)) / (mean(a) or 1e-300)


def integrate(t, y):
    # the CSV is uniform in t; rectangle rule over the steps of the cycle
    dt = (t[-1] - t[0]) / max(1, len(t) - 1)
    return dt * sum(y)


def metrics(full, case):
    qf, qc = full["qIn"], case["qIn"]
    n = min(len(qf), len(qc))
    qf, qc = qf[:n], qc[:n]
    t = full["t"][:n]
    ref = mean([abs(q) for q in qf]) or 1e-300
    diff = [c - f for c, f in zip(qc, qf)]
    fwd = sum(d for d, f in zip(diff, qf) if f >= 0.0) / n
    rev = sum(d for d, f in zip(diff, qf) if f < 0.0) / n
    ef = integrate(t, full["dissipation"][:n])
    ec = integrate(t, case["dissipation"][:n])
    wpf = mean([abs(w) for w in full["wssParent"][:n]])
    wpc = mean([abs(w) for w in case["wssParent"][:n]])
    wsf = mean([abs(w) for w in full["wssStub"][:n]])
    wsc = mean([abs(w) for w in case["wssStub"][:n]])
    pct = lambda a, b: 100.0 * (a - b) / (abs(b) if b else 1e-300)
    return {
        "Qf_mean_abs": ref,
        "Qf_mean_signed": mean(qf),
        "reversal_fraction": sum(1 for f in qf if f < 0.0) / n,
        "dQ_signed": 100.0 * mean(diff) / ref,
        "dQ_abs": 100.0 * (mean([abs(q) for q in qc]) - ref) / ref,
        "dQ_fwd": 100.0 * fwd / ref,
        "dQ_rev": 100.0 * rev / ref,
        "dQ_rms": 100.0 * math.sqrt(mean([d * d for d in diff])) / ref,
        "dQ_max": 100.0 * max(abs(d) for d in diff) / (max(abs(q) for q in qf) or 1e-300),
        "dE": pct(ec, ef),
        "dWSS_parent": pct(wpc, wpf),
        "dWSS_stub": pct(wsc, wsf),
        "picard_mean": mean(case["picardIts"]),
        "picard_max": max(case["picardIts"]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--cycles", type=int, default=0, help="cycles run (0: last in the CSV)")
    a = ap.parse_args()

    rows = []
    figures = []
    for group in sorted(os.listdir(a.out)):
        gdir = os.path.join(a.out, group)
        if not (group.startswith("s") and os.path.isdir(gdir)):
            continue
        for state in sorted(os.listdir(gdir)):
            sdir = os.path.join(gdir, state)
            fpath = os.path.join(sdir, "full", "RmuAxi.csv")
            if not os.path.isfile(fpath):
                continue
            full = read(fpath)
            if not full:
                continue
            full_period = periodicity(full)
            full = last_cycle(full, a.cycles)
            cases = {}
            for case in sorted(os.listdir(sdir)):
                cpath = os.path.join(sdir, case, "RmuAxi.csv")
                if case == "full" or not os.path.isfile(cpath):
                    continue
                data = read(cpath)
                if not data:
                    continue
                period = periodicity(data)
                data = last_cycle(data, a.cycles)
                cases[case] = data
                m = metrics(full, data)
                rows.append(dict(group=group, state=state, closure=case, **m,
                                 periodicity=period, full_periodicity=full_period))
            figures.append((group, state, full, cases))

    if not rows:
        print("no complete case pairs found under", a.out)
        return 1

    keys = list(rows[0].keys())
    with open(os.path.join(a.out, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    show = ["group", "state", "closure", "reversal_fraction", "dQ_signed", "dQ_abs",
            "dQ_fwd", "dQ_rev", "dQ_rms", "dQ_max", "dE", "dWSS_parent", "dWSS_stub",
            "picard_max", "periodicity", "full_periodicity"]
    print("| " + " | ".join(show) + " |")
    print("|" + "---|" * len(show))
    for r in rows:
        cells = []
        for k in show:
            v = r[k]
            cells.append("%.3g" % v if isinstance(v, float) else str(v))
        print("| " + " | ".join(cells) + " |")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return 0

    for group, state, full, cases in figures:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        t0 = full["t"][0]
        tt = [(t - t0) for t in full["t"]]
        ax[0].plot(tt, [q * 1e9 for q in full["qIn"]], "k", lw=2, label="full")
        for name, d in cases.items():
            ax[0].plot(tt[:len(d["qIn"])], [q * 1e9 for q in d["qIn"]], lw=1, label=name)
            ax[1].plot(tt[:len(d["qIn"])],
                       [100.0 * (c - f) / (mean([abs(q) for q in full["qIn"]]) or 1e-300)
                        for c, f in zip(d["qIn"], full["qIn"])], lw=1, label=name)
        ax[0].set_xlabel("t - t0 (s)")
        ax[0].set_ylabel("Q_in (mm^3/s)")
        ax[1].set_xlabel("t - t0 (s)")
        ax[1].set_ylabel("(Q_c - Q_full) / mean|Q_full| (%)")
        ax[0].legend()
        ax[1].axhline(0.0, color="k", lw=0.5)
        fig.suptitle("%s  %s" % (group, state))
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "%s_%s_qin.png" % (group, state)), dpi=150)
        plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Transfer of an outlet resistance away from its calibration.

Quasi-steady 0D experiments built on the WRMS elements of sobol_rmu.py. Every
model is calibrated once, at one operating point and with one blood, and is
then used, without recalibration, where that calibration no longer holds:

  pressure   the stepped pipe of the RmuAxi benchmark at steady inlet pressures
             phi Pbar, phi in [0.1, 4] (operating-point transfer);
  split      a Y whose daughters narrow beyond the cut, calibrated at the
             reference head and driven at phi times that head, phi in
             [0.05, 1]: change in the narrow daughter's share of the flow;
  blood      the stepped pipe calibrated with healthy blood and run with the
             diabetic and hyperviscous laws under the same inlet pressure;
  matched    the constant rescaled so that its cycle-mean |Q| equals that of
             the full pipe (the best constant, which needs the answer).

Closures, as in the paper: rmu (one WRMS element of the radius of the last
removed tube, its length matched at the reference flow), rest (the constant
R_0 = R_mu(Q_0) of the removed tubes) and, for the blood transfer, rest
rescaled by mu_inf (the conventional plasma/high-shear correction) and rest
recalibrated with the new blood at Q_0.

  python3 transfer_rmu.py --out transfer
"""
import argparse
import csv
import math
import os
import sys

import numpy as np

import sobol_rmu as S

Pi = math.pi
Healthy = "healthy"
States = list(S.STATES)


def curveOf(state):
    mu0, mi, lam, n, a = S.STATES[state]
    return S.HydraulicCurve([mu0], [mi], [lam], [n], [a])


def muInfOf(state):
    return S.STATES[state][1]


class SeriesBranch:
    """Tubes in series, Delta p(Q) and its inverse on a log grid of flows."""

    FlowNodes = 1601

    def __init__(self, curve, tubes, outlet=None):
        self.m_curve = curve
        self.m_tubes = tubes
        self.m_outlet = outlet              # callable Q -> Delta p of the outlet

    def pressure(self, Q):
        Q = np.atleast_2d(Q)
        p = sum(self.m_curve.resistance(t.radius, t.length, Q) for t in self.m_tubes) * Q
        if self.m_outlet is not None:
            p = p + self.m_outlet(Q)
        return p

    def flow(self, p, qScale):
        """Q such that pressure(Q) = p (p >= 0), by tabulation and interpolation."""
        q = qScale * np.concatenate([[0.0], np.logspace(-6, 4, self.FlowNodes - 1)])[None, :]
        dp = self.pressure(q)
        return S.rowInterp(np.atleast_2d(p), dp, q)


class SteppedPipe:
    """The stepped pipe of the benchmark in SI units, R_1 = 0.1 mm."""

    Scale = 1.0e-4

    def __init__(self, calibrationState, gammaRef):
        self.m_full = [S.Tube(t.radius * self.Scale, t.length * self.Scale)
                       for t in S.Benchmark.FullTubes]
        self.m_kept, self.m_removed = S.splitTubes(self.m_full, S.Benchmark.Truncation * self.Scale)
        rMin = min(t.radius for t in self.m_full)
        self.m_q0 = Pi * rMin ** 3 * gammaRef / 4.0
        cal = curveOf(calibrationState)
        self.m_calibrationState = calibrationState
        self.m_pbar = float(SeriesBranch(cal, self.m_full).pressure(self.m_q0)[0, 0])
        self.m_rest = float(SeriesBranch(cal, self.m_removed).pressure(self.m_q0)[0, 0]) / self.m_q0
        rs = self.m_removed[-1].radius
        muS = float(cal(np.array([[4.0 * self.m_q0 / (Pi * rs ** 3)]]))[0, 0])
        self.m_single = S.Tube(rs, self.m_rest * Pi * rs ** 4 / (8.0 * muS))

    def branch(self, state, closure):
        curve = curveOf(state)
        if closure == "full":
            return SeriesBranch(curve, self.m_full)
        if closure == "rmu":
            t = self.m_single
            return SeriesBranch(curve, self.m_kept,
                                lambda Q: curve.resistance(t.radius, t.length, Q) * Q)
        if closure == "rest":
            return SeriesBranch(curve, self.m_kept, lambda Q: self.m_rest * Q)
        if closure == "restInf":
            k = muInfOf(state) / muInfOf(self.m_calibrationState)
            return SeriesBranch(curve, self.m_kept, lambda Q: k * self.m_rest * Q)
        if closure == "restNew":
            r = float(SeriesBranch(curve, self.m_removed).pressure(self.m_q0)[0, 0]) / self.m_q0
            return SeriesBranch(curve, self.m_kept, lambda Q: r * Q)
        raise ValueError(closure)

    def waveform(self, state, closure, phi=1.0, amplitude=0.0, samples=400, scale=1.0):
        t = (np.arange(samples) + 0.5) / samples
        p = phi * self.m_pbar * (1.0 + amplitude * np.sin(2.0 * Pi * t))
        b = self.branch(state, closure)
        if scale != 1.0:
            inner = b.m_outlet
            b.m_outlet = lambda Q: scale * inner(Q)
        return np.sign(p) * b.flow(np.abs(p), self.m_q0)[0]


def errors(qc, qf):
    mean = 100.0 * (np.abs(qc).mean() - np.abs(qf).mean()) / np.abs(qf).mean()
    wave = 100.0 * np.sqrt(((qc - qf) ** 2).mean()) / np.abs(qf).mean()
    return mean, wave


class SteppedY:
    """A Y whose daughters narrow beyond the cut (radii in m).

    Parent 100 um x 450 um; daughters 72 and 46 um, retained over 300 um. The
    removed part of each daughter is the extension of the 3D benchmark (360
    and 520 um), whose distal half narrows to 0.6 of the daughter radius.
    """

    Parent = S.Tube(100e-6, 450e-6)
    Daughters = [(72e-6, 300e-6, 360e-6), (46e-6, 300e-6, 520e-6)]   # wide, narrow
    Narrowing = 0.6
    PressureNodes = 1201

    def __init__(self, state, gammaRef):
        """gammaRef: apparent shear rate of the narrow daughter at calibration."""
        self.m_curve = curveOf(state)
        self.m_kept, self.m_removed = [], []
        for r, lk, lr in self.Daughters:
            self.m_kept.append([S.Tube(r, lk)])
            self.m_removed.append([S.Tube(r, 0.5 * lr), S.Tube(self.Narrowing * r, 0.5 * lr)])
        self.m_qScale = Pi * self.Parent.radius ** 3 * 1000.0 / 4.0
        full = self.closures("full", None)
        # head at which the narrow daughter carries gammaRef
        rn = self.Daughters[1][0]
        qn = Pi * rn ** 3 * gammaRef / 4.0
        pj = float(full[1].pressure(qn)[0, 0])
        qw = float(full[0].flow(pj, self.m_qScale)[0, 0])
        self.m_head = pj + float(SeriesBranch(self.m_curve, [self.Parent]).pressure(qn + qw)[0, 0])
        self.m_q0 = [qw, qn]

    def closures(self, closure, q0):
        out = []
        for k in range(2):
            kept, removed = self.m_kept[k], self.m_removed[k]
            if closure == "full":
                out.append(SeriesBranch(self.m_curve, kept + removed))
                continue
            rest = float(SeriesBranch(self.m_curve, removed).pressure(q0[k])[0, 0]) / q0[k]
            if closure == "rest":
                out.append(SeriesBranch(self.m_curve, kept, lambda Q, r=rest: r * Q))
            elif closure == "rmu":
                rs = removed[-1].radius
                mu = float(self.m_curve(np.array([[4.0 * q0[k] / (Pi * rs ** 3)]]))[0, 0])
                t = S.Tube(rs, rest * Pi * rs ** 4 / (8.0 * mu))
                out.append(SeriesBranch(self.m_curve, kept,
                                        lambda Q, t=t: self.m_curve.resistance(t.radius, t.length, Q) * Q))
            else:
                raise ValueError(closure)
        return out

    def narrowShare(self, closure, phi):
        """Narrow daughter's percentage of the outflow at inlet heads phi * head."""
        wide, narrow = self.closures(closure, self.m_q0)
        pj = self.m_head * np.logspace(-6, 0.5, self.PressureNodes)[None, :]
        qw = wide.flow(pj, self.m_qScale)
        qn = narrow.flow(pj, self.m_qScale)
        pin = pj + SeriesBranch(self.m_curve, [self.Parent]).pressure(qw + qn)
        heads = np.atleast_2d(np.asarray(phi) * self.m_head)
        qwI = S.rowInterp(heads, pin, qw)
        qnI = S.rowInterp(heads, pin, qn)
        return (100.0 * qnI / (qwI + qnI))[0]


def bestScale(pipe, state, amplitude):
    """Factor on R_0 for which the cycle-mean |Q| equals that of the full pipe."""
    qf = pipe.waveform(state, "full", amplitude=amplitude)
    target = np.abs(qf).mean()
    lo, hi = 0.05, 20.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        qm = pipe.waveform(state, "rest", amplitude=amplitude, scale=mid)
        if np.abs(qm).mean() > target:
            lo = mid
        else:
            hi = mid
    k = math.sqrt(lo * hi)
    return k, errors(pipe.waveform(state, "rest", amplitude=amplitude, scale=k), qf)


def run(out):
    os.makedirs(out, exist_ok=True)
    gammas = (30.0, 100.0, 200.0, 400.0)

    # 1. operating-point transfer of the stepped pipe (steady)
    phis = np.logspace(-1, math.log10(4.0), 41)
    pressureRows = []
    for state in States:
        for g in (100.0, 400.0):
            pipe = SteppedPipe(state, g)
            for phi in phis:
                qf = pipe.waveform(state, "full", phi)[0]
                row = [state, g, phi]
                for c in ("rest", "rmu"):
                    row.append(100.0 * (pipe.waveform(state, c, phi)[0] - qf) / qf)
                pressureRows.append(row)
    write(out, "pressure.csv", ["state", "gamma0", "phi", "rest_dQ", "rmu_dQ"], pressureRows)

    # 2. flow split of the stepped Y below and at the calibration head
    phisY = np.linspace(0.05, 1.0, 39)
    splitRows = []
    for state in States:
        for g in (100.0, 400.0):
            y = SteppedY(state, g)
            ref = {c: y.narrowShare(c, 1.0)[0] for c in ("full", "rest", "rmu")}
            shares = {c: y.narrowShare(c, phisY) for c in ("full", "rest", "rmu")}
            for i, phi in enumerate(phisY):
                splitRows.append([state, g, phi] + [shares[c][i] - ref[c] for c in ("full", "rest", "rmu")]
                                 + [shares[c][i] for c in ("full", "rest", "rmu")])
    write(out, "split.csv", ["state", "gamma0", "phi", "full_dShare", "rest_dShare", "rmu_dShare",
                             "full_share", "rest_share", "rmu_share"], splitRows)

    # 3. blood-state transfer: calibrated with healthy blood, run with each law
    bloodRows = []
    for g in gammas:
        pipe = SteppedPipe(Healthy, g)
        for A in (0.5, 1.0):
            for state in States:
                qf = pipe.waveform(state, "full", amplitude=A)
                for c in ("rmu", "rest", "restInf", "restNew"):
                    m, w = errors(pipe.waveform(state, c, amplitude=A), qf)
                    bloodRows.append([g, A, state, c, m, w])
    write(out, "blood.csv", ["gamma0", "A", "state", "closure", "mean_dQ", "wave_dQ"], bloodRows)

    # 4. constant matched to the cycle-mean |Q| of the full pipe
    matchedRows = []
    for state in States:
        for g in gammas:
            pipe = SteppedPipe(state, g)
            for A in (0.5, 1.0):
                k, (m, w) = bestScale(pipe, state, A)
                qf = pipe.waveform(state, "full", amplitude=A)
                _, w0 = errors(pipe.waveform(state, "rest", amplitude=A), qf)
                _, wm = errors(pipe.waveform(state, "rmu", amplitude=A), qf)
                matchedRows.append([state, g, A, k, w0, w, wm])
    write(out, "matched.csv", ["state", "gamma0", "A", "scale", "rest_wave", "matched_wave", "rmu_wave"],
          matchedRows)
    report(out)


def write(out, name, header, rows):
    with open(os.path.join(out, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def report(out):
    read = lambda n: list(csv.DictReader(open(os.path.join(out, n))))
    print("\nOperating-point transfer (steady), |error in Q| over phi in [0.2, 2]:")
    for r0 in ("100.0", "400.0"):
        for state in States:
            rows = [r for r in read("pressure.csv") if r["state"] == state and r["gamma0"] == r0
                    and 0.2 - 1e-9 <= float(r["phi"]) <= 2.0 + 1e-9]
            print("  %-13s gd0=%5s  R_0 up to %5.1f%%   R_mu up to %4.2f%%" % (
                state, r0, max(abs(float(r["rest_dQ"])) for r in rows),
                max(abs(float(r["rmu_dQ"])) for r in rows)))
    print("\nStepped Y: change in narrow share (pp) at phi = 0.2:")
    for r in read("split.csv"):
        if abs(float(r["phi"]) - 0.2) < 0.02:
            print("  %-13s gd0=%5s full %+5.2f  R_0 %+5.2f  R_mu %+5.2f  (share %.2f)" % (
                r["state"], r["gamma0"], float(r["full_dShare"]), float(r["rest_dShare"]),
                float(r["rmu_dShare"]), float(r["full_share"])))
    print("\nBlood transfer (healthy calibration), mean / waveform error %:")
    for r in read("blood.csv"):
        print("  gd0=%5s A=%s %-13s %-8s %+6.1f / %5.1f" % (r["gamma0"], r["A"], r["state"], r["closure"],
                                                         float(r["mean_dQ"]), float(r["wave_dQ"])))
    print("\nConstant matched to the mean |Q|: waveform error % (R_0 -> matched; R_mu):")
    for r in read("matched.csv"):
        print("  %-13s gd0=%5s A=%s  k=%.3f  %5.1f -> %5.1f ; %4.2f" % (
            r["state"], r["gamma0"], r["A"], float(r["scale"]), float(r["rest_wave"]),
            float(r["matched_wave"]), float(r["rmu_wave"])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="transfer")
    a = ap.parse_args()
    run(a.out)


if __name__ == "__main__":
    sys.exit(main())

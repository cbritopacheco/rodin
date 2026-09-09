#!/usr/bin/env python3
"""Paired R_mu vs constant-resistance comparison at fixed state.

For each scenario, takes the CSV of the rheology-dependent run and the CSV of
its constant-outlet twin (identical 3D rheology, identical calibration, only
the reduced closure differs) and reports, over the last complete cycle, the
quantities of section 7.5 of the manuscript:

  d<|Q_in|>   = <|Q_in|>^{const} / <|Q_in|>^{R_mu} - 1
                the change in cycle-averaged coronary inflow caused by
                discarding the shear dependence (paper: +3.4 / +2.0 / +7.4 %
                for healthy / diabetic / hyperviscous at rest);
  d<mu_eff>   = <Phi>^{const} / <Phi>^{R_mu} - 1, per resistive element
                how much the constant closure understates the cycle-averaged
                effective hydraulic viscosity (paper: -8.6 / -5.2 / -24.6 %);
  E_Q(t)      = (|Q_in^{R_mu}| - |Q_in^{const}|) / |Q_in^{const}|
                reported as cycle mean, and as the extreme value with the time
                (in reduced phase t/T) at which it occurs;
  max ratio   = max_t Phi^{R_mu}/Phi^{const} per element, the instantaneous
                viscosity understatement (paper: exceeds an order of magnitude
                near proximal flow closure).

Phi is read from the CSV columns ArteriolarViscosityRatio /
VenularViscosityRatio (Phi = mu_ap / mu_N), so Phi^{R_mu}/Phi^{const} is
exactly mu_eff^{R_mu}/mu_inf. The two runs are aligned by cycle phase and the
constant run is linearly interpolated onto the phases of the R_mu run, so a
different time-step history (adaptivity) does not bias the comparison.

Usage:
  python3 compare_closures.py [--period T] [--csv-out FILE] \\
      name=rmu.csv,const.csv [name2=... ...]
"""
import csv
import sys

COLS = {
    't': 't',
    'qin': 'CoronaryInletFlux',
    'qdist': 'CoronaryDistalFluxTotal',
    'mua': 'ArteriolarViscosityRatio',
    'muv': 'VenularViscosityRatio',
}


def read_csv(path):
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        missing = [c for c in COLS.values() if c not in reader.fieldnames]
        if missing:
            sys.exit(f'{path}: missing columns {missing}')
        data = {k: [] for k in COLS}
        for row in reader:
            for k, c in COLS.items():
                data[k].append(float(row[c]))
    return data


def last_cycle(data, T):
    t = data['t']
    t0 = t[-1] - T
    idx = [i for i, ti in enumerate(t) if t0 - 1e-12 <= ti <= t[-1] + 1e-12]
    if len(idx) < 10:
        sys.exit(f'last cycle has only {len(idx)} samples; check --period')
    return {k: [v[i] for i in idx] for k, v in data.items()}


def trapz(y, t):
    return sum(0.5 * (y[i] + y[i + 1]) * (t[i + 1] - t[i]) for i in range(len(t) - 1))


def mean(y, t):
    return trapz(y, t) / (t[-1] - t[0])


def phase(t, T):
    return [ti % T for ti in t]


def interp_by_phase(src_t, src_y, T, targets):
    """Linear interpolation of src_y(phase) at the given target phases,
    treating phase as periodic on [0, T)."""
    pairs = sorted(zip(phase(src_t, T), src_y))
    ph = [p for p, _ in pairs]
    yy = [y for _, y in pairs]
    # wrap guards
    ph = [ph[-1] - T] + ph + [ph[0] + T]
    yy = [yy[-1]] + yy + [yy[0]]
    out = []
    for q in targets:
        lo, hi = 0, len(ph) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if ph[mid] <= q:
                lo = mid
            else:
                hi = mid
        span = ph[hi] - ph[lo]
        w = 0.0 if span <= 0 else (q - ph[lo]) / span
        out.append(yy[lo] + w * (yy[hi] - yy[lo]))
    return out


def compare(rmu_path, const_path, T):
    r = last_cycle(read_csv(rmu_path), T)
    c = last_cycle(read_csv(const_path), T)
    t = r['t']
    ph_r = phase(t, T)

    ci = interp_by_phase(c['t'], c['qin'], T, ph_r)
    ca = interp_by_phase(c['t'], c['mua'], T, ph_r)
    cv = interp_by_phase(c['t'], c['muv'], T, ph_r)

    q_r = mean([abs(q) for q in r['qin']], t)
    q_c = mean([abs(q) for q in ci], t)

    eq = []
    for a, b in zip(r['qin'], ci):
        den = abs(b)
        eq.append((abs(a) - den) / den if den > 1e-18 else 0.0)
    eq_ext = max(eq, key=abs)
    eq_ext_phase = ph_r[eq.index(eq_ext)] / T

    out = {
        'mean |Q_in| R_mu [uL/s]': q_r * 1e9,
        'mean |Q_in| const [uL/s]': q_c * 1e9,
        'd<|Q_in>| const vs R_mu [%]': 100.0 * (q_c / q_r - 1.0) if q_r else float('nan'),
        'mean E_Q [%]': 100.0 * mean(eq, t),
        'extreme E_Q [%]': 100.0 * eq_ext,
        'extreme E_Q at t/T': eq_ext_phase,
    }
    for lbl, yr, yc in (('proximal', r['mua'], ca), ('distal', r['muv'], cv)):
        mr, mc = mean(yr, t), mean(yc, t)
        out[f'd<mu_eff> {lbl} [%]'] = 100.0 * (mc / mr - 1.0) if mr else float('nan')
        ratios = [a / b if b > 1e-30 else float('nan') for a, b in zip(yr, yc)]
        out[f'max mu_eff/mu_inf {lbl}'] = max(ratios)
    return out


def main(argv):
    T = 0.85
    csv_out = None
    runs = []
    it = iter(argv)
    for a in it:
        if a == '--period':
            T = float(next(it))
        elif a == '--csv-out':
            csv_out = next(it)
        elif a in ('-h', '--help'):
            print(__doc__)
            return 0
        else:
            name, _, paths = a.rpartition('=')
            parts = paths.split(',')
            if len(parts) != 2:
                sys.exit(f'expected name=rmu.csv,const.csv but got "{a}"')
            runs.append((name or parts[0], parts[0], parts[1]))
    if not runs:
        print(__doc__)
        return 1

    rows = {n: compare(a, b, T) for n, a, b in runs}
    keys = list(next(iter(rows.values())).keys())
    w = max(len(k) for k in keys) + 2
    print(' ' * w + ''.join(f'{n:>16s}' for n in rows))
    for k in keys:
        line = f'{k:<{w}s}'
        for n in rows:
            line += f'{rows[n][k]:16.3f}'
        print(line)

    if csv_out:
        with open(csv_out, 'w', newline='') as f:
            wcsv = csv.writer(f)
            wcsv.writerow(['metric'] + list(rows))
            for k in keys:
                wcsv.writerow([k] + [rows[n][k] for n in rows])
        print(f'\nwritten: {csv_out}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

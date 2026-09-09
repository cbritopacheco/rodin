#!/usr/bin/env python3
"""Cycle metrics for the coupled LV-0D/coronary-3D runs.

Reads one or more CoronaryArtery.csv files (as written by
CoupledLV0DCoronary3D::writeCSVRow) and reports, for the last complete
cardiac cycle of each run:

  - mean coronary inflow <|Q_in|> and forward/retrograde volumes per beat,
    with the retrograde fraction V_retro / V_forward  (normal hearts: a few
    per cent; grows strongly in LVH / aortic stenosis);
  - the reversal phase: total duration with Q_in < 0, number of episodes,
    peak retrograde flow;
  - systolic/diastolic inflow ratio (systole defined as p_LV above a
    fraction of its cycle maximum; LAD physiology ~0.3-0.5);
  - distal nominal shear-rate scale gamma_d(t) = gamma0 * |Q_d(t)| / Q0
    (Q0 taken from --baseline when given, so scenarios share the reference),
    its minimum, cycle mean, and the fraction of the cycle spent below the
    low-shear thresholds (defaults: 47 s^-1 = healthy eps=0.2 and
    108 s^-1 = healthy eps=0.1 / hyperviscous eps=0.2);
  - minimum transmural pressure and the viscosity-ratio extrema of the
    reduced elements.

Usage:
  python3 retrograde_metrics.py [opts] name=path.csv [name2=path2.csv ...]
  python3 retrograde_metrics.py CoronaryArtery.csv

Options:
  --period T          cardiac period in s (default 0.85)
  --gamma0-distal G   distal reference shear rate in 1/s (default 400)
  --gamma0-proximal G proximal reference shear rate in 1/s (default 800)
  --thresholds A,B    low-shear thresholds in 1/s (default 47,108)
  --baseline PATH     baseline CSV whose cycle-mean distal/proximal flows
                      define Q0 for every scenario (recommended)
  --systole-frac F    p_LV fraction defining systole (default 0.25)
  --csv-out PATH      also write the metrics table as CSV
"""
import csv
import math
import sys

COLS = {
    't': 't',
    'qin': 'CoronaryInletFlux',
    'qdist': 'CoronaryDistalFluxTotal',
    'qout': 'CoronaryOutletFluxTotal',
    'plv': 'LeftVentriclePressure',
    'pim': 'IntramyoPressure',
    'ptm': 'TransmuralPressure',
    'muv': 'VenularViscosityRatio',
    'mua': 'ArteriolarViscosityRatio',
}


def read_csv(path):
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        missing = [c for c in COLS.values() if c not in reader.fieldnames]
        if missing:
            sys.exit(f'{path}: missing columns {missing}; found {reader.fieldnames}')
        data = {k: [] for k in COLS}
        for row in reader:
            for k, c in COLS.items():
                data[k].append(float(row[c]))
    return data


def last_cycle(data, T):
    t = data['t']
    t1 = t[-1]
    t0 = t1 - T
    idx = [i for i, ti in enumerate(t) if t0 - 1e-12 <= ti <= t1 + 1e-12]
    if len(idx) < 10:
        sys.exit(f'last cycle has only {len(idx)} samples; check --period')
    return {k: [v[i] for i in idx] for k, v in data.items()}


def trapz(y, t):
    return sum(0.5 * (y[i] + y[i + 1]) * (t[i + 1] - t[i]) for i in range(len(t) - 1))


def metrics(cyc, T, g0d, g0p, thresholds, q0d, q0p, sys_frac):
    t, qin = cyc['t'], cyc['qin']
    dur = t[-1] - t[0]

    v_fwd = trapz([max(q, 0.0) for q in qin], t)
    v_ret = trapz([max(-q, 0.0) for q in qin], t)
    mean_abs_qin = trapz([abs(q) for q in qin], t) / dur

    # reversal phase
    rev_time = trapz([1.0 if q < 0 else 0.0 for q in qin], t)
    episodes = sum(1 for i in range(1, len(qin)) if qin[i] < 0 <= qin[i - 1])
    peak_ret = min(qin)

    # systole from p_LV
    plv = cyc['plv']
    pmax, pmin = max(plv), min(plv)
    thr = pmin + sys_frac * (pmax - pmin)
    in_sys = [p > thr for p in plv]
    q_sys = [abs(q) for q, s in zip(qin, in_sys) if s]
    q_dia = [abs(q) for q, s in zip(qin, in_sys) if not s]
    sysdia = (sum(q_sys) / len(q_sys)) / (sum(q_dia) / len(q_dia)) \
        if q_sys and q_dia else float('nan')
    sys_time = trapz([1.0 if s else 0.0 for s in in_sys], t)

    # distal / proximal nominal shear
    gd = [g0d * abs(q) / q0d for q in cyc['qdist']]
    gp = [g0p * abs(q) / q0p for q in cyc['qout']]
    out = {
        'mean |Q_in| [uL/s]': mean_abs_qin * 1e9,
        'V_forward [uL/beat]': v_fwd * 1e9,
        'V_retro [uL/beat]': v_ret * 1e9,
        'retro fraction [%]': 100.0 * v_ret / v_fwd if v_fwd > 0 else float('nan'),
        'reversal time [ms]': rev_time * 1e3,
        'reversal episodes': episodes,
        'peak retro Q [uL/s]': -peak_ret * 1e9 if peak_ret < 0 else 0.0,
        'sys/dia inflow ratio': sysdia,
        'systole span [ms]': sys_time * 1e3,
        'gamma_d min [1/s]': min(gd),
        'gamma_d mean [1/s]': trapz(gd, t) / dur,
        'gamma_p mean [1/s]': trapz(gp, t) / dur,
        'min p_tm [Pa]': min(cyc['ptm']),
        'max Phi_v': max(cyc['muv']),
        'mean Phi_a': trapz(cyc['mua'], t) / dur,
    }
    for th in thresholds:
        frac = trapz([1.0 if g < th else 0.0 for g in gd], t) / dur
        out[f'gamma_d < {th:g}/s [% cycle]'] = 100.0 * frac
    return out


def main(argv):
    T, g0d, g0p = 0.85, 400.0, 800.0
    thresholds = [47.0, 108.0]
    baseline = None
    sys_frac = 0.25
    csv_out = None
    runs = []

    it = iter(argv)
    for a in it:
        if a == '--period':
            T = float(next(it))
        elif a == '--gamma0-distal':
            g0d = float(next(it))
        elif a == '--gamma0-proximal':
            g0p = float(next(it))
        elif a == '--thresholds':
            thresholds = [float(x) for x in next(it).split(',')]
        elif a == '--baseline':
            baseline = next(it)
        elif a == '--systole-frac':
            sys_frac = float(next(it))
        elif a == '--csv-out':
            csv_out = next(it)
        elif a in ('-h', '--help'):
            print(__doc__)
            return 0
        else:
            name, _, path = a.rpartition('=')
            runs.append((name or path, path))
    if not runs:
        print(__doc__)
        return 1

    # Q0 reference from the baseline cycle means (so scenarios share it)
    ref = baseline or runs[0][1]
    ref_cyc = last_cycle(read_csv(ref), T)
    dur = ref_cyc['t'][-1] - ref_cyc['t'][0]
    q0d = trapz([abs(q) for q in ref_cyc['qdist']], ref_cyc['t']) / dur
    q0p = trapz([abs(q) for q in ref_cyc['qout']], ref_cyc['t']) / dur
    print(f'# reference (Q0) from {ref}: '
          f'distal {q0d*1e9:.3f} uL/s, proximal {q0p*1e9:.3f} uL/s\n')

    rows = {}
    for name, path in runs:
        cyc = last_cycle(read_csv(path), T)
        rows[name] = metrics(cyc, T, g0d, g0p, thresholds, q0d, q0p, sys_frac)

    keys = list(next(iter(rows.values())).keys())
    w = max(len(k) for k in keys) + 2
    header = ' ' * w + ''.join(f'{n:>16s}' for n in rows)
    print(header)
    for k in keys:
        line = f'{k:<{w}s}'
        for n in rows:
            v = rows[n][k]
            line += f'{v:16.3f}' if isinstance(v, float) else f'{v:16d}'
        print(line)

    if csv_out:
        with open(csv_out, 'w', newline='') as f:
            wcsv = csv.writer(f)
            wcsv.writerow(['metric'] + [n for n in rows])
            for k in keys:
                wcsv.writerow([k] + [rows[n][k] for n in rows])
        print(f'\nwritten: {csv_out}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

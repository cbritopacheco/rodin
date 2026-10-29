#!/usr/bin/env python3
# Copyright Carlos BRITO PACHECO 2026.
# Distributed under the Boost Software License, Version 1.0.
"""Pressure-driven fluid checks on an affine tetrahedral tube, without a solid solve."""
import argparse
import csv
import json
import math
import os
import re
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET


def write_mesh(path):
    """Write a rotated octagonal tube in mm with inlet 1, outlet 6, and wall 4."""
    points = []
    for layer in range(3):
        for sector in range(9):
            angle = 2 * math.pi * sector / 8
            x, y = (8 * math.cos(angle), 8 * math.sin(angle)) if sector < 8 else (0, 0)
            z = 4 * layer
            points.append((math.cos(.6) * x + math.sin(.6) * z, y,
                           -math.sin(.6) * x + math.cos(.6) * z))
    cells = []
    for layer in range(2):
        for sector in range(8):
            a, b, c = sorted((9 * layer + sector, 9 * layer + (sector + 1) % 8,
                              9 * layer + 8))
            A, B, C = a + 9, b + 9, c + 9
            for cell in [(a, b, c, A), (b, c, A, B), (c, A, B, C)]:
                q = [points[i] for i in cell]
                u, v, w = [tuple(q[j][k] - q[0][k] for k in range(3)) for j in range(1, 4)]
                det = sum(u[k] * (v[(k+1)%3] * w[(k+2)%3] - v[(k+2)%3] * w[(k+1)%3])
                          for k in range(3))
                assert det != 0
                cells.append(cell if det > 0 else (cell[1], cell[0], cell[2], cell[3]))
    faces = {}
    for cell in cells:
        for omit in range(4):
            key = tuple(sorted(cell[j] for j in range(4) if j != omit))
            faces[key] = faces.get(key, 0) + 1
    with path.open('w') as out:
        out.write(f'MeshVersionFormatted 2\nDimension 3\nVertices\n{len(points)}\n')
        for p in points:
            out.write(' '.join(format(x, '.17g') for x in p) + ' 0\n')
        out.write(f'Tetrahedra\n{len(cells)}\n')
        for cell in cells:
            out.write(' '.join(str(i+1) for i in cell) + ' 1\n')
        boundary = [key for key, count in faces.items() if count == 1]
        out.write(f'Triangles\n{len(boundary)}\n')
        for face in boundary:
            layers = {i // 9 for i in face}
            label = 1 if layers == {0} else 6 if layers == {2} else 4
            out.write(' '.join(str(i+1) for i in face) + f' {label}\n')
        out.write('End\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[4]
    parser.add_argument('--executable', type=Path,
                        default=root/'build/examples/Heart/Atrium/Atrium')
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = (args.output or Path(tempfile.mkdtemp(prefix='rodin-atrium-fluid-'))).resolve()
    output.mkdir(parents=True, exist_ok=True)
    mesh = output/'Tube.mesh'
    write_mesh(mesh)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')

    def run(name, ranks, pressure, outlet=0, dt=.001, steps=3, amplitude=0,
            kinetics=False, lag=1, max_it=100000):
        folder = output/name
        folder.mkdir(exist_ok=True)
        cmd = [args.mpiexec, '-n', str(ranks), str(args.executable.resolve()),
               '-atrium_mode', 'fluid', '-atrium_mesh', str(mesh),
               '-atrium_output', str(folder), '-atrium_steps', str(steps),
               '-atrium_output_every', '0', '-atrium_dt', str(dt),
               '-atrium_fluid_pin', str(pressure), '-atrium_fluid_pin_amplitude', str(amplitude),
               '-atrium_fluid_outlet_pressure', str(outlet), '-atrium_fluid_lv_model', 'false',
               '-atrium_fluid_kinetics', str(kinetics).lower(),
               '-atrium_fluid_flow_cycles', '0', '-atrium_fluid_lag_preconditioner', str(lag),
               '-ksp_type', 'gmres' if lag > 1 else 'preonly',
               '-ksp_max_it', str(max_it), '-ksp_rtol', '1e-12', '-pc_type', 'lu',
               '-pc_factor_mat_solver_type', 'mumps', '-mat_mumps_icntl_20', '0',
               '-ksp_diagonal_scale', 'false', '-ksp_diagonal_scale_fix', 'false']
        if lag > 1:
            cmd += ['-ksp_norm_type', 'unpreconditioned',
                    '-ksp_reuse_preconditioner', 'true', '-log_view', f':{folder}/PETSc.log']
        with (folder/'run.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env,
                           cwd=output, check=True)
        with (folder/'Diagnostics.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        fields = {a.attrib['Name'] for a in
                  ET.parse(folder/'Fluid.xdmf').iter('Attribute')}
        expected = {'Attribute', 'Velocity', 'Pressure', 'WallShearStress',
                    'TimeAveragedWallShearStress', 'OscillatoryShearIndex'}
        if kinetics:
            expected |= {'Thrombin', 'Fibrinogen', 'Fibrin', 'ThrombosisActivation'}
        assert fields == expected, (name, fields)
        assert len(rows) == steps + 1, (name, len(rows))
        assert abs(float(rows[-1]['time_s']) - dt*steps) < 1e-12
        for row in rows[1:]:
            for key in ('fluid_max_velocity_m_per_s', 'fluid_inlet_flux_m3_per_s',
                        'fluid_outlet_flux_m3_per_s'):
                assert math.isfinite(float(row[key])), (name, key)
        return rows[1:]

    rest = run('rest', 2, 0)
    run('kinetics', 2, 1, steps=1, kinetics=True)
    assert all(float(row['fluid_max_velocity_m_per_s']) == 0 for row in rest)
    forward = run('forward', 2, 1)
    serial = run('serial', 1, 1)
    reverse = run('reverse', 2, -1, steps=1)
    gauge = run('gauge', 2, 1001, outlet=1000, steps=1)
    speed = 'fluid_max_velocity_m_per_s'
    inlet, outlet = 'fluid_inlet_flux_m3_per_s', 'fluid_outlet_flux_m3_per_s'
    refreshed = run('refresh', 2, 1, steps=12, lag=5)
    fresh = run('refresh-fresh', 2, 1, steps=12)
    assert [int(r['fluid_pc_refreshed']) for r in refreshed] == [1,0,0,0,0,1,0,0,0,0,1,0]
    assert [int(r['fluid_pc_age_timesteps']) for r in refreshed] == [0,1,2,3,4,0,1,2,3,4,0,1]
    assert int(refreshed[-1]['fluid_pc_refreshes']) == 3
    log = (output/'refresh/PETSc.log').read_text()
    for event, expected in [('MatLUFactorSym', 1), ('MatLUFactorNum', 3)]:
        assert int(re.search(r'^'+event+r'\s+(\d+)', log, re.M)[1]) == expected
    for a, b in zip(refreshed, fresh):
        for key in (speed, inlet, outlet):
            assert abs(float(a[key])-float(b[key])) <= 1e-7 * max(abs(float(b[key])), 1e-12)
    retried = run('refresh-retry', 2, 1, steps=3, lag=5, max_it=1)
    assert [int(r['fluid_pc_retried']) for r in retried] == [0,1,1]
    assert int(retried[-1]['fluid_pc_refreshes']) == 3
    for a, b in zip(forward, serial):
        for key in (speed, inlet, outlet):
            assert abs(float(a[key])-float(b[key])) <= 1e-9 * max(abs(float(a[key])), 1e-12)
    assert float(forward[0][inlet]) < 0 < float(forward[0][outlet])
    for key in (inlet, outlet):
        assert abs(float(forward[0][key])+float(reverse[0][key])) < 1e-15
    assert abs(float(forward[0][speed])-float(gauge[0][speed])) < 1e-9
    imbalance = max(abs(float(row[inlet])+float(row[outlet])) /
                    max(abs(float(row[inlet])), abs(float(row[outlet]))) for row in forward)
    assert imbalance < 1e-6
    # A smooth pressure source starts from zero and avoids the abrupt-load startup.
    q = []
    for dt in (.001, .0005, .00025, .000125):
        row = run(f'smooth-{dt}', 2, 0, dt=dt, steps=round(.008/dt), amplitude=1)[-1]
        q.append(float(row[outlet]))
    differences = [abs(a-b) for a, b in zip(q, q[1:])]
    assert all(0 < b < a for a, b in zip(differences, differences[1:])), q
    orders = [math.log2(a/b) for a, b in zip(differences, differences[1:])]
    assert all(.75 < order < 1.25 for order in orders), orders
    report = dict(rest=True, mpi_consistency=True, reversal=True, gauge=True,
                  maximum_relative_flux_imbalance=imbalance,
                  smooth_source_outlet_flux=q, smooth_source_observed_orders=orders)
    (output/'checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print('Output:', output)


if __name__ == '__main__':
    main()

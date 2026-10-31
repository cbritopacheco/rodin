#!/usr/bin/env python3
# Copyright Carlos BRITO PACHECO 2026.
# Distributed under the Boost Software License, Version 1.0.
"""Check analytic interface Newton, state reuse, and MPI consistency."""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
from FSI import write_common_mesh


def main():
    root = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path,
                        default=root/'build/examples/Heart/Atrium/Atrium')
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = (args.output or Path(tempfile.mkdtemp(prefix='rodin-atrium-newton-'))).resolve()
    output.mkdir(parents=True, exist_ok=True)
    mesh = output/'TubeWall.mesh'
    write_common_mesh(mesh)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')
    apex = ','.join(map(str, [math.cos(.6)*8, 0, -math.sin(.6)*8]))

    def run(name, ranks, taylor=False, guard_recovery=False):
        folder = output/name
        folder.mkdir()
        cmd = [args.mpiexec, '-n', str(ranks), str(args.executable.resolve()),
               '-atrium_mode', 'fsi', '-atrium_mesh', str(mesh),
               '-atrium_output', str(folder), '-atrium_steps', '2',
               '-atrium_dt', '.02', '-atrium_fsi_newton',
               '-atrium_fsi_atol', '1e-11', '-atrium_fsi_solid_atol', '1e-11',
               '-atrium_solid_apex', apex, '-atrium_solid_contractility', '100',
               '-atrium_solid_cache_constitutive', '-atrium_output_every', '1',
               '-ksp_type', 'preonly', '-pc_type', 'lu',
               '-pc_factor_mat_solver_type', 'mumps', '-mat_mumps_icntl_20', '0',
               '-ksp_diagonal_scale', 'false', '-ksp_diagonal_scale_fix', 'false']
        if taylor:
            cmd.append('-atrium_fsi_check_tangent')
        if guard_recovery:
            # Deliberate overshoot: valid solutions stay below the guard, while
            # trial states cross it and must be rejected without committing.
            cmd += ['-atrium_fluid_max_velocity', '.0015',
                    '-atrium_fsi_interface_snes_linesearch_damping', '50',
                    '-atrium_fsi_interface_snes_max_it', '25']
        with (folder/'run.log').open('w') as log:
            subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT,
                           check=True, cwd=root, timeout=180)
        log = (folder/'run.log').read_text()
        results = re.findall(r'converged: iterations=(\d+) residual=(\S+) evaluations=(\d+)', log)
        assert len(results) == 2, log
        assert all(float(r) < 1e-11 for _, r, _ in results)
        if not taylor and not guard_recovery:
            assert all(int(e) == int(n)+1 for n, _, e in results), results
        elif taylor:
            remainders = [float(r) for r in re.findall(r'^TAYLOR h=\S+ remainder=(\S+)', log, re.M)]
            assert len(remainders) == 8, remainders
            for offset in (0, 4):
                ratios = [a/b for a,b in zip(remainders[offset:offset+3], remainders[offset+1:offset+4])]
                assert all(3.8 < ratio < 4.2 for ratio in ratios), ratios
        if guard_recovery:
            rejected = [int(n) for n in re.findall(r'rejected_trials=(\d+)', log)]
            assert len(rejected) == 2 and all(n > 0 for n in rejected), log
            assert '[FSI/Newton trial]: rejected:' in log
        assert 'ANALYTIC_J_ZERO' not in log
        rows = list(csv.DictReader((folder/'Diagnostics.csv').open()))
        assert len(rows) == 3
        for row in rows[1:]:
            assert int(row['fsi_converged']) == 1
            assert int(row['fsi_budget_exhausted']) == 0
            assert float(row['fsi_min_mesh_jacobian']) > 0
            assert float(row['fsi_velocity_mismatch_m_per_s']) < 1e-8
            assert abs(float(row['fsi_cavity_volume_difference_ml'])) < 1e-9
        return rows

    serial = run('serial', 1)
    parallel = run('parallel', 2, taylor=True)
    recovered = run('guard-recovery', 2, guard_recovery=True)
    for a,b,c in zip(serial[1:], parallel[1:], recovered[1:]):
        assert float(c['fluid_max_velocity_m_per_s']) <= .0015
        for key in ('max_displacement_mm', 'solid_cavity_volume_ml', 'fluid_max_velocity_m_per_s'):
            assert math.isclose(float(a[key]), float(b[key]), rel_tol=1e-8, abs_tol=1e-12), (key,a[key],b[key])
            assert math.isclose(float(a[key]), float(c[key]), rel_tol=1e-8, abs_tol=1e-12), (key,a[key],c[key])
    print('Analytic tangents, state reuse, guard backtracking, shared geometry and MPI consistency: passed')
    print('Output:', output)


if __name__ == '__main__':
    main()

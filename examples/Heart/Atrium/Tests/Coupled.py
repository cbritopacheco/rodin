#!/usr/bin/env python3
# Copyright Oscar RUZ NUNES 2026.
# Distributed under the Boost Software License, Version 1.0.
"""Check pressure-driven FSI with split subscales, the lumped ventricle, and species.

The conforming tube is advanced with the mixed wall, term-by-term orthogonal
subscales, the 0D ventricle behind the mitral diode, and P-T-A-G transport on
the moving mesh. The checks cover the analytic interface tangent, transport of
a uniform state, nonnegative species, and MPI consistency.
"""
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
    output = (args.output or Path(tempfile.mkdtemp(prefix='rodin-atrium-coupled-'))).resolve()
    output.mkdir(parents=True, exist_ok=True)
    mesh = output/'TubeWall.mesh'
    write_common_mesh(mesh)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')
    apex = ','.join(map(str, [math.cos(.6)*8, 0, -math.sin(.6)*8]))
    fibrinogen = 7.35e-3

    def run(name, ranks, steps, extra=()):
        folder = output/name
        folder.mkdir()
        cmd = [args.mpiexec, '-n', str(ranks), str(args.executable.resolve()),
               '-atrium_mode', 'fsi', '-atrium_mesh', str(mesh),
               '-atrium_output', str(folder), '-atrium_steps', str(steps),
               '-atrium_dt', '.02', '-atrium_fsi_newton',
               '-atrium_fsi_atol', '1e-11', '-atrium_fsi_solid_atol', '1e-11',
               '-atrium_solid_apex', apex, '-atrium_solid_contractility', '100',
               '-atrium_solid_cache_constitutive', '-atrium_output_every', '1',
               '-atrium_fluid_vms', '-atrium_fluid_lv_model',
               '-atrium_fluid_pressure_ramp', '.1',
               '-atrium_fluid_kinetics', '-atrium_fluid_period', '.04',
               '-atrium_fluid_flow_cycles', '1',
               '-ksp_type', 'preonly', '-pc_type', 'lu',
               '-pc_factor_mat_solver_type', 'mumps', '-mat_mumps_icntl_20', '0',
               '-ksp_diagonal_scale', 'false', '-ksp_diagonal_scale_fix', 'false',
               *extra]
        with (folder/'run.log').open('w') as log:
            subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT,
                           check=True, cwd=root, timeout=600)
        log = (folder/'run.log').read_text()
        results = re.findall(r'converged: iterations=(\d+) residual=(\S+) evaluations=(\d+)', log)
        assert len(results) == steps, log
        assert all(float(r) < 1e-11 for _, r, _ in results)
        rows = list(csv.DictReader((folder/'Diagnostics.csv').open()))
        assert len(rows) == steps+1
        for row in rows[1:]:
            assert int(row['fsi_converged']) == 1
            assert float(row['fsi_min_mesh_jacobian']) > 0
            assert float(row['fsi_velocity_mismatch_m_per_s']) < 1e-8
        return rows, log

    # Tangent of the pressure-driven split-subscale problem, including a
    # step with nonzero lagged projections.
    _, log = run('tangent', 2, 2, ['-atrium_fsi_check_tangent'])
    remainders = [float(r) for r in re.findall(r'^TAYLOR h=\S+ remainder=(\S+)', log, re.M)]
    assert len(remainders) == 8, remainders
    for offset in (0, 4):
        ratios = [a/b for a,b in zip(remainders[offset:offset+3], remainders[offset+1:offset+4])]
        assert all(3.5 < ratio < 4.5 for ratio in ratios), ratios

    # Species start after the first two-step cycle. Without wall initiation
    # the uniform state is transported exactly on the moving mesh.
    rows, _ = run('uniform', 2, 6, ['-atrium_fluid_ptag_wall_flux', '0'])
    for row in rows[3:]:
        assert abs(float(row['fluid_min_fibrinogen']) - fibrinogen) < 1e-9*fibrinogen, row
        assert abs(float(row['fluid_max_thrombin'])) < 1e-15, row

    serial, _ = run('serial', 1, 6)
    parallel, _ = run('parallel', 2, 6)
    for a,b in zip(serial[1:], parallel[1:]):
        for key in ('max_displacement_mm', 'solid_cavity_volume_ml', 'fluid_max_velocity_m_per_s'):
            assert math.isclose(float(a[key]), float(b[key]), rel_tol=1e-7, abs_tol=1e-12), (key,a[key],b[key])
    for row in parallel[3:]:
        assert float(row['fluid_max_thrombin']) >= 0
        assert float(row['fluid_min_fibrinogen']) <= fibrinogen*(1+1e-9)
        assert float(row['fluid_max_fibrin']) >= 0
    print('Split-subscale tangent, ventricular driver, moving-mesh species, and MPI consistency: passed')
    print('Output:', output)


if __name__ == '__main__':
    main()

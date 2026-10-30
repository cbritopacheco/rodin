#!/usr/bin/env python3
# Copyright Carlos BRITO PACHECO 2026.
# Distributed under the Boost Software License, Version 1.0.
"""Check shared run options, state cadence, and component outputs on two MPI ranks."""
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from FSI import write_common_mesh


def main():
    root = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path,
                        default=root/'build/examples/Heart/Atrium/Atrium')
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = (args.output or Path(tempfile.mkdtemp(prefix='rodin-atrium-parameters-'))).resolve()
    output.mkdir(parents=True, exist_ok=True)
    mesh = output/'TubeWall.mesh'
    write_common_mesh(mesh)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')
    base = [args.mpiexec, '-n', '2', str(args.executable.resolve())]
    help_run = subprocess.run(base + ['-help', '-atrium_mesh', '/nonexistent.mesh'],
                              env=env, text=True, capture_output=True, timeout=60)
    assert help_run.returncode == 0, help_run.stdout + help_run.stderr
    for flag in ('-atrium_mode', '-atrium_dt', '-atrium_solid_contractility',
                 '-atrium_fluid_pin', '-atrium_fluid_max_velocity', '-atrium_fsi_rtol', '-atrium_fsi_newton',
                 '-atrium_fsi_atol', '-atrium_fsi_solid_atol'):
        assert flag in help_run.stdout, flag
    solver = ['-ksp_type', 'gmres', '-ksp_rtol', '1e-12', '-pc_type', 'lu',
              '-pc_factor_mat_solver_type', 'mumps', '-mat_mumps_icntl_20', '0',
              '-ksp_diagonal_scale', 'false', '-ksp_diagonal_scale_fix', 'false']
    apex = ','.join(map(str, [math.cos(.6)*8, 0, -math.sin(.6)*8]))
    for mode in ('solid', 'fluid', 'fsi'):
        folder = output/mode
        cmd = base + ['-atrium_mode', mode, '-atrium_mesh', str(mesh),
                      '-atrium_output', str(folder), '-atrium_dt', '.001',
                      '-atrium_steps', '3', '-atrium_output_every', '2',
                      '-atrium_solid_apex', apex, '-atrium_solid_cache_constitutive'] + solver
        folder.mkdir(exist_ok=True)
        with (folder/'run.log').open('w') as log:
            subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=180)
        rows = list(csv.DictReader((folder/'Diagnostics.csv').open()))
        assert len(rows) == 4
        assert [float(r['time_s']) for r in rows] == [0, .001, .002, .003]
        assert list(folder.glob('*.csv')) == [folder/'Diagnostics.csv']
        for component in ('Solid', 'Fluid'):
            path = folder/(component+'.xdmf')
            enabled = mode == 'fsi' or mode == component.lower()
            assert path.exists() == enabled, (mode, component)
            if enabled:
                times = [float(t.attrib['Value']) for t in ET.parse(path).iter('Time')]
                assert times == [0, .002, .003], (mode, component, times)
        print(mode, 'shared time, cadence, outputs: passed', flush=True)
    for options in (['-atrium_mode', 'invalid'], ['-atrium_dt', '0'],
                    ['-atrium_steps', '-1'], ['--dt', '.001'],
                    ['-atrium_solid_apex', '1,2'], ['-atrium_pin', '1'],
                    ['-atrium_fluid_max_velocity', '0'],
                    ['-atrium_fluid_max_velocity', '-1'],
                    ['-atrium_fluid_max_velocity', 'nan'],
                    ['-atrium_typo', '1']):
        failed = subprocess.run(base + options, env=env, text=True,
                                capture_output=True, timeout=60)
        assert failed.returncode != 0, options
    print('Help and invalid configuration checks: passed', flush=True)


if __name__ == '__main__':
    main()

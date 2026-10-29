#!/usr/bin/env python3
# Copyright Carlos BRITO PACHECO 2026.
# Distributed under the Boost Software License, Version 1.0.
"""Contraction-driven MPI FSI on a conforming tetrahedron/wedge tube."""
import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from Fluid import write_mesh


def write_common_mesh(path):
    write_mesh(path)
    lines = path.read_text().splitlines()
    def section(name, cast):
        i = lines.index(name)
        return [list(map(cast, line.split())) for line in
                lines[i+2:i+2+int(lines[i+1])]]
    vertices = section('Vertices', float)
    tetrahedra = section('Tetrahedra', int)
    triangles = section('Triangles', int)
    faces = {}
    for tet in tetrahedra:
        for i in range(4):
            faces[tuple(sorted(tet[j] for j in range(4) if j != i))] = tet[:4]
    wall = [t for t in triangles if t[3] == 4]
    inner = sorted({v for t in wall for v in t[:3]})
    outer = {}
    for v in inner:
        x, y, z = vertices[v-1][:3]
        radial_x = math.cos(.6)*x - math.sin(.6)*z
        r = math.hypot(radial_x, y)
        nx, ny = radial_x/r, y/r
        outer[v] = len(vertices)+1
        vertices.append([x+math.cos(.6)*nx, y+ny, z-math.sin(.6)*nx, 4])
    wedges = []
    for triangle in wall:
        a, b, c = triangle[:3]
        p, q, r = [vertices[v-1][:3] for v in (a,b,c)]
        e = [q[i]-p[i] for i in range(3)]
        f = [r[i]-p[i] for i in range(3)]
        normal = [e[1]*f[2]-e[2]*f[1], e[2]*f[0]-e[0]*f[2], e[0]*f[1]-e[1]*f[0]]
        tet = faces[tuple(sorted((a,b,c)))]
        inside = [sum(vertices[v-1][i] for v in tet)/4 for i in range(3)]
        center = [(p[i]+q[i]+r[i])/3 for i in range(3)]
        if sum(normal[i]*(center[i]-inside[i]) for i in range(3)) < 0:
            b, c = c, b
        wedges.append([a,b,c,outer[a],outer[b],outer[c],4])
        triangles.append([outer[a],outer[b],outer[c],7])
    quads = []
    for label in (1,6):
        counts = Counter(tuple(sorted((t[i],t[(i+1)%3])))
                         for t in triangles if t[3] == label for i in range(3))
        for (a,b), count in counts.items():
            if count == 1:
                quads.append([a,b,outer[b],outer[a],label])
    with path.open('w') as out:
        out.write('MeshVersionFormatted 2\nDimension 3\nVertices\n'+str(len(vertices))+'\n')
        for row in vertices:
            out.write(' '.join(format(v,'.17g') for v in row[:3])+' 0\n')
        for name, rows in [('Tetrahedra',tetrahedra), ('Wedges',wedges),
                           ('Triangles',triangles), ('Quadrilaterals',quads)]:
            out.write(name+'\n'+str(len(rows))+'\n')
            for row in rows:
                out.write(' '.join(map(str,row))+'\n')
        out.write('End\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[4]
    parser.add_argument('--executable', type=Path,
                        default=root/'build/examples/Heart/Atrium/Atrium')
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = (args.output or Path(tempfile.mkdtemp(prefix='rodin-atrium-fsi-'))).resolve()
    output.mkdir(parents=True, exist_ok=True)
    mesh = output/'TubeWall.mesh'
    write_common_mesh(mesh)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')
    apex = ','.join(map(str, [math.cos(.6)*8, 0, -math.sin(.6)*8]))

    def run(name, ranks=2, contractility=20, budget=5, tolerance=1e-6,
            steps=3, dt=.001, lag=5, mode='fsi'):
        folder = output/name
        folder.mkdir(exist_ok=True)
        cmd = [args.mpiexec, '-n', str(ranks), str(args.executable.resolve()),
               '-atrium_mode', mode, '-atrium_mesh', str(mesh),
               '-atrium_output', str(folder), '-atrium_steps', str(steps),
               '-atrium_dt', str(dt), '-atrium_fsi_iterations', str(budget),
               '-atrium_fsi_rtol', str(tolerance), '-atrium_solid_apex', apex,
               '-atrium_solid_contractility', str(contractility),
               '-atrium_solid_cache_constitutive', '-atrium_output_every', '1',
               '-atrium_fluid_lag_preconditioner', str(lag),
               '-ksp_type', 'gmres', '-ksp_rtol', '1e-12', '-ksp_max_it', '500',
               '-pc_type', 'lu', '-pc_factor_mat_solver_type', 'mumps',
               '-mat_mumps_icntl_20', '0', '-ksp_diagonal_scale', 'false',
               '-ksp_diagonal_scale_fix', 'false']
        with (folder/'run.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                           env=env, check=True, cwd=root)
        rows = list(csv.DictReader((folder/'Diagnostics.csv').open()))
        assert len(rows) == steps+1
        assert math.isclose(float(rows[-1]['time_s']), steps*dt)
        return rows

    rest = run('rest', contractility=0)
    assert all(float(r['max_displacement_mm']) == 0 for r in rest)
    assert all(float(r['fluid_max_velocity_m_per_s'] or 0) == 0 for r in rest)
    assert all(int(r['fsi_converged']) == 1 for r in rest[1:])
    fixed = run('fixed-budget', budget=3, tolerance=0)
    assert all(int(r['fsi_iterations']) == 3 and int(r['fsi_budget_exhausted']) == 1 for r in fixed[1:])
    one = run('one-sweep', budget=1)
    assert all(int(r['fsi_iterations']) == 1 and int(r['fsi_budget_exhausted']) == 1
               and int(r['fsi_converged']) == 0 for r in one[1:])
    assert 'last sweep accepted' in (output/'one-sweep/run.log').read_text()
    active = run('active')
    serial = run('active-serial', ranks=1)
    four = run('active-four', ranks=4)
    fresh = run('active-fresh', lag=1)
    tight = run('tight', budget=60, tolerance=1e-8)
    loose = run('loose', budget=5, tolerance=1)
    finite = run('finite-motion', contractility=100, dt=.02, steps=10)
    solid = run('solid', mode='solid')
    for rows in (active, serial, four, fresh, tight, loose, finite):
        for row in rows[1:]:
            iterations = int(row['fsi_iterations'])
            converged = bool(int(row['fsi_converged']))
            assert converged != bool(int(row['fsi_budget_exhausted']))
            if converged:
                assert iterations > 1
                assert float(row['fsi_displacement_change']) < float(row['fsi_displacement_tolerance'])
            assert float(row['fsi_min_mesh_jacobian']) > 0
            assert float(row['fsi_displacement_mismatch_m']) < 1e-12
            assert abs(float(row['fsi_cavity_volume_difference_ml'])) < 1e-9
            fluxes = [float(row[k]) for k in ('fluid_inlet_flux_m3_per_s',
                      'fluid_outlet_flux_m3_per_s', 'fluid_wall_flux_m3_per_s')]
            assert abs(sum(fluxes)) < 1e-7*max(map(abs,fluxes)) + 1e-15
    for comparison in (serial, four, fresh):
        for a,b in zip(active[1:],comparison[1:]):
            for key in ('max_displacement_mm','solid_cavity_volume_ml','fluid_max_velocity_m_per_s'):
                assert abs(float(a[key])-float(b[key])) < 1e-6*max(abs(float(b[key])),1e-9), (key,a[key],b[key])
    assert all(int(r['fsi_converged']) for r in tight[1:])
    assert all(int(r['fsi_iterations']) <= 2 for r in loose[1:])
    assert max(float(r['fluid_max_velocity_m_per_s']) for r in active[1:]) > 0
    assert abs(float(active[-1]['max_displacement_mm'])-float(solid[-1]['max_displacement_mm'])) > 1e-10
    comparison = {key:dict(budgeted=float(active[-1][key]), tight=float(tight[-1][key]))
                  for key in ('max_displacement_mm','solid_cavity_volume_ml',
                              'fluid_max_velocity_m_per_s','fsi_velocity_mismatch_m_per_s')}
    report = dict(rest=True, budget_exhaustion_commits=True, early_displacement_convergence=True,
                  mpi_consistency=True, shared_interface_volume=True, global_flux_balance=True,
                  finite_motion_geometry=True, fluid_feedback=True, budget_vs_tight=comparison)
    (output/'checks.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    print('Output:',output)


if __name__ == '__main__':
    main()

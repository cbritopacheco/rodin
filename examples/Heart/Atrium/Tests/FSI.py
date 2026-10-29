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
    apex = [math.cos(.6)*8, 0, -math.sin(.6)*8]
    def run(name, ranks, contractility, mode='fsi', lag=5, steps=5, atol=1e-4, dt=.001):
        folder = output/name
        folder.mkdir(exist_ok=True)
        cmd = [args.mpiexec, '-n', str(ranks), str(args.executable.resolve()),
               '-atrium_mode', mode, '-atrium_mesh', str(mesh),
               '-atrium_output', str(folder/'Flow'), '--output', str(folder/'Solid'),
               '-atrium_steps', str(steps), '-atrium_dt', str(dt), '--dt', str(dt),
               '--steps', str(steps), '-atrium_fsi_atol', str(atol), '--apex', *map(str,apex), '--apex-radius', '3',
               '--contractility', str(contractility), '--cache-constitutive',
               '-atrium_output_every', '1', '-atrium_fluid_lag_preconditioner', str(lag),
               '-ksp_type', 'gmres', '-ksp_rtol', '1e-12', '-ksp_max_it', '500',
               '-pc_type', 'lu', '-pc_factor_mat_solver_type', 'mumps',
               '-mat_mumps_icntl_20', '0', '-ksp_diagonal_scale', 'false',
               '-ksp_diagonal_scale_fix', 'false']
        with (folder/'run.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                           env=env, check=True, cwd=root)
        with (folder/'Solid.csv').open() as stream:
            return list(csv.DictReader(stream))
    rest = run('rest', 2, 0)
    assert all(float(r['max_displacement_mm']) == 0 for r in rest)
    assert all(float(r['fluid_max_velocity_m_per_s'] or 0) == 0 for r in rest)
    parallel = run('active', 2, 20)
    serial = run('active-serial', 1, 20)
    four = run('active-four-ranks', 4, 20)
    fresh = run('active-fresh', 2, 20, lag=1)
    tight = run('active-tight', 2, 20, atol=1e-6)
    probe = run('robin-probe', 2, 20, steps=1, atol=1)
    assert float(probe[-1]['fsi_velocity_mismatch_m_per_s']) > 1e-12
    finite = run('finite-motion', 2, 100, steps=10, dt=.02)
    for r in finite[1:]:
        assert float(r['fsi_min_mesh_jacobian'])>0
        assert abs(float(r['fsi_cavity_volume_difference_ml']))<1e-9
    solid = run('solid', 2, 20, mode='solid')
    for r in parallel[1:]:
        assert float(r['fluid_inlet_pressure_pa']) == 0
        assert float(r['fluid_outlet_pressure_pa']) == 0
        assert not r['fluid_arterial_pressure_pa']
        assert float(r['fsi_min_mesh_jacobian']) > 0
        assert float(r['fsi_robin_residual_pa']) <= float(r['fsi_robin_tolerance_pa'])
        tolerance = math.sqrt(3)*float(r['fsi_robin_tolerance_pa'])/float(r['fsi_robin_coefficient_pa_s_per_m'])
        assert float(r['fsi_velocity_mismatch_m_per_s']) <= tolerance*(1+1e-6)
        assert float(r['fsi_displacement_mismatch_m']) < 1e-12
        assert abs(float(r['fsi_cavity_volume_difference_ml'])) < 1e-9
        flux = sum(float(r[k]) for k in
                   ('fluid_inlet_flux_m3_per_s','fluid_outlet_flux_m3_per_s','fluid_wall_flux_m3_per_s'))
        scale = max(abs(float(r[k])) for k in
                    ('fluid_inlet_flux_m3_per_s','fluid_outlet_flux_m3_per_s','fluid_wall_flux_m3_per_s'))
        assert abs(flux) < 1e-7*scale + 1e-15
    assert max(float(r['max_displacement_mm']) for r in parallel) > 1e-10
    assert max(float(r['fluid_max_velocity_m_per_s'] or 0) for r in parallel) > 1e-9
    for comparison in (serial, four):
        for a,b in zip(parallel[1:],comparison[1:]):
            for key in ('max_displacement_mm','solid_cavity_volume_ml','fluid_max_velocity_m_per_s'):
                assert abs(float(a[key])-float(b[key])) < 1e-6*max(abs(float(b[key])),1e-9)
    for a,b in zip(parallel[1:],fresh[1:]):
        for key in ('max_displacement_mm','fluid_max_velocity_m_per_s'):
            assert abs(float(a[key])-float(b[key])) < 1e-7*max(abs(float(b[key])),1e-9)
    assert float(tight[-1]['fsi_velocity_mismatch_m_per_s']) < float(parallel[-1]['fsi_velocity_mismatch_m_per_s'])
    assert abs(float(tight[-1]['max_displacement_mm'])-float(parallel[-1]['max_displacement_mm'])) < .01*float(tight[-1]['max_displacement_mm'])
    assert float(parallel[1]['max_displacement_mm']) < float(solid[1]['max_displacement_mm'])
    assert abs(float(parallel[-1]['max_displacement_mm'])-
               float(solid[-1]['max_displacement_mm'])) > 1e-10
    report = dict(rest=True, motion_driven_flow=True, mpi_consistency=True,
                  shared_interface_volume=True, global_flux_balance=True,
                  interface_continuity_to_robin_tolerance=True, retained_pc_matches_fresh=True,
                  robin_probe_has_finite_slip=True, tighter_tolerance_reduces_slip=True, finite_motion_geometry=True,
                  fluid_feedback_changes_solid=True, initial_fluid_load_resists_contraction=True)
    (output/'checks.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    print('Output:',output)


if __name__ == '__main__':
    main()

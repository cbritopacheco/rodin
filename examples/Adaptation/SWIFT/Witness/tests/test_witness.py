# Copyright (c) 2026 Rodin contributors.
# Distributed under the Boost Software License, Version 1.0.
"""C++/JSON integration, independent numerical checks and plotting tests."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot import Plot
from reference import Reference


class WitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[5]
        default = root/'build-p1-3d-clang19/examples/Adaptation/SWIFT/Witness/SWIFT_Witness'
        cls.binary = Path(os.environ.get('SWIFT_WITNESS_BINARY', default)).resolve()
        if not cls.binary.is_file():
            raise unittest.SkipTest('Set SWIFT_WITNESS_BINARY to the built C++ executable.')

    def execute(self, *arguments, supplied=None):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'result.json'
            command = [str(self.binary), *map(str, arguments), '--output', str(output)]
            if supplied is not None:
                input_path = Path(directory)/'input.json'
                input_path.write_text(json.dumps(supplied))
                command += ['--evaluate', str(input_path)]
            process = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=60)
            self.assertEqual(process.returncode, 0, process.stderr)
            return json.loads(output.read_text())

    def test_disabled_tool_does_not_require_cgal(self):
        cmake = shutil.which('cmake')
        if cmake is None:
            self.skipTest('CMake is unavailable.')
        source = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'CMakeLists.txt').write_text(
                'cmake_minimum_required(VERSION 3.16)\n'
                'project(WitnessOptional NONE)\n'
                f'add_subdirectory("{source.as_posix()}" witness)\n')
            process = subprocess.run(
                [cmake, '-S', str(root), '-B', str(root/'build'),
                 '-DRODIN_BUILD_SWIFT_WITNESS=OFF', '-DCMAKE_DISABLE_FIND_PACKAGE_CGAL=TRUE'],
                capture_output=True, text=True, timeout=60)
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertNotIn('CGAL_DIR:', (root/'build/CMakeCache.txt').read_text())

    def test_reference_geometry_and_wireframe(self):
        records = self.execute('1', '--geometry', 'all', '--iterations', '0')
        self.assertEqual(len(records), 8)
        edges = {'point': 0, 'segment': 1, 'triangle': 3, 'quadrilateral': 4,
                 'tetrahedron': 6, 'pyramid': 8, 'hexahedron': 12, 'wedge': 9}
        for record in records:
            name = record['geometry']
            np.testing.assert_array_equal(record['reference_vertices'], Reference.VERTICES[name])
            self.assertEqual(len(record['reference_edges']), edges[name])
            self.assertEqual(record['count'], len(Reference.VERTICES[name]))
            self.assertIsInstance(record['covering_radius'], (int, float))

    def test_analytic_interval(self):
        for subdivisions in (1, 3, 16):
            count = subdivisions+1
            record = self.execute(subdivisions, '--geometry', 'segment')
            np.testing.assert_allclose(record['points'], ((np.arange(count)+.5)/count)[:, None])
            self.assertAlmostEqual(record['covering_radius'], 1/(2*count), places=12)

    def test_uniform_lattice_initialization(self):
        fixtures = [('segment', 3, 2), ('triangle', 6, 2),
                    ('quadrilateral', 9, 2), ('tetrahedron', 10, 2),
                    ('pyramid', 14, 2), ('hexahedron', 27, 2), ('wedge', 18, 2)]
        for geometry, count, resolution in fixtures:
            with self.subTest(geometry=geometry):
                record = self.execute(resolution, '--geometry', geometry, '--iterations', 0)
                initial = np.asarray(record['initial_points'])
                self.assertEqual(record['initialization'], {
                    'method': 'complete_uniform_lattice', 'subdivisions': resolution,
                    'lattice_size': count})
                self.assertEqual(len(np.unique(initial, axis=0)), count)
                np.testing.assert_allclose(initial*resolution, np.round(initial*resolution))
                vertices = Reference.VERTICES[geometry]
                self.assertTrue(all(any(np.array_equal(v, p) for p in initial) for v in vertices))

    def test_lattice_count_and_sampling_independence(self):
        first = self.execute(4, '--iterations', 0, '--resolution', 4)
        second = self.execute('--subdivisions', 4, '--iterations', 0, '--resolution', 12)
        np.testing.assert_array_equal(first['initial_points'], second['initial_points'])
        self.assertEqual(first['initialization'], {
            'method': 'complete_uniform_lattice', 'subdivisions': 4, 'lattice_size': 15})
        self.assertEqual(len(np.unique(first['initial_points'], axis=0)), 15)
        for record in self.execute(1, '--geometry', 'all', '--iterations', 0):
            self.assertEqual(sorted(record['initial_points']), sorted(record['reference_vertices']))

    def test_triangle_six_equispaced(self):
        record = self.execute(supplied={'geometry': 'triangle', 'points':
            [[0, 0], [1, 0], [0, 1], [.5, 0], [0, .5], [.5, .5]]})
        self.assertAlmostEqual(record['covering_radius'], np.sqrt(2)/4, places=12)
        holes = sorted(record['worst_locations'])
        np.testing.assert_allclose(holes, [[.25, .25], [.25, .75], [.75, .25]], atol=1e-12)

    def test_independent_qhull_oracle_all_geometries(self):
        rng = np.random.default_rng(13)
        for name in Reference.VERTICES:
            reference = Reference(name)
            fixtures = [reference.vertices]
            if reference.dimension:
                for count in (1, 6, 16):
                    fixtures.append(rng.dirichlet(np.ones(len(reference.vertices)), count) @ reference.vertices)
            for points in fixtures:
                with self.subTest(geometry=name, count=len(points)):
                    record = self.execute(supplied={'geometry': name, 'points': points.tolist()})
                    self.assertAlmostEqual(record['covering_radius'], reference.radius(points), places=9)
                    self.assertTrue(all(record['cells']))

    def test_unrestricted_search_determinism_and_budget(self):
        options = ('4', '--iterations', '5')
        first, second = self.execute(*options), self.execute(*options)
        np.testing.assert_array_equal(first['points'], second['points'])
        self.assertLessEqual(first['covering_radius'], first['initial_radius'])
        self.assertEqual(first['optimality'], 'not_certified')
        limited = self.execute('4', '--max-evaluations', '3')
        self.assertEqual(limited['search']['evaluations'], 3)
        self.assertEqual(limited['search']['termination'], 'evaluation_limit')
        self.assertEqual(limited['search']['implementation'], 'checked_enclosing_ball')
        self.assertEqual(limited['search']['enclosing_ball_backend'], 'CGAL::Min_sphere_of_spheres_d')
        radii = [row['radius'] for row in limited['history']]
        self.assertTrue(all(b <= a+1e-10 for a, b in zip(radii, radii[1:])))

    def test_canonical_enclosing_ball_all_geometries(self):
        records = self.execute(2, '--geometry', 'all')
        for record in records:
            with self.subTest(geometry=record['geometry']):
                radius = Reference(record['geometry']).radius(np.array(record['points']))
                self.assertAlmostEqual(radius, record['covering_radius'], places=9)
                self.assertLessEqual(radius, record['initial_radius']+1e-10)
                if record['geometry'] == 'quadrilateral':
                    self.assertAlmostEqual(radius, np.sqrt(2)/6, places=5)

    def test_before_after_plot_all_geometries(self):
        records = self.execute('3', '--geometry', 'all', '--iterations', '0')
        with tempfile.TemporaryDirectory() as directory:
            for record in records:
                output = Path(directory)/f'{record["geometry"]}-single.svg'
                Plot(record).save(output)
                root = ET.parse(output).getroot()
                ids = [node.attrib['id'] for node in root.iter() if 'id' in node.attrib]
                self.assertEqual(len(ids), len(set(ids)))
                if record['dimension'] == 3:
                    root = ET.parse(output).getroot()
                    labels = [node.text or '' for node in root.iter()]
                    for view in ('XYZ', 'XY', 'YZ', 'XZ'):
                        self.assertTrue(any(label.endswith(f': {view}') for label in labels))
                    self.assertTrue(any(node.get('class') == 'uncovered-ball' for node in root.iter()))
                    points = np.array([[.1, .2, .3]])
                    np.testing.assert_array_equal(Plot(record).project(points, 'XYZ'),
                                                  Plot(record).project(points))
                    for view, axes in [('XY', [0, 1]), ('YZ', [1, 2]), ('XZ', [0, 2])]:
                        np.testing.assert_array_equal(Plot(record).project(points, view), points[:, axes])

    def test_single_geometry_plot_command(self):
        scripts = Path(__file__).resolve().parents[1]
        for geometry in ('triangle', 'tetrahedron'):
            with self.subTest(geometry=geometry), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                record = self.execute(4, '--geometry', geometry)
                source = root/'single.json'
                source.write_text(json.dumps(record))
                output = root/'before-after.svg'
                process = subprocess.run(
                    [sys.executable, str(scripts/'plot.py'), str(source), '--output', str(output)],
                    capture_output=True, text=True, timeout=60)
                self.assertEqual(process.returncode, 0, process.stderr)
                svg = ET.parse(output).getroot()
                labels = [node.text or '' for node in svg.iter()]
                self.assertTrue(any('initialization' in label for label in labels))
                self.assertTrue(any('Numerical covering candidate' in label for label in labels))
                if geometry == 'triangle':
                    disks = [node for node in svg.iter()
                             if node.tag.endswith('circle') and node.get('clip-path')]
                    self.assertGreaterEqual(len(disks), 2)

    def test_invalid_input(self):
        for arguments in (('0',), ('--resolution', '0'), ('2', '--subdivisions', '3'),
                          ('--geometry', 'unknown'), ('--objective', 'unknown'),
                          ('--step-tolerance', 'nan'), ('--subdivisions', '4x'), ('--n', '4'),
                          ('16', '--lobatto-order', '2')):
            with tempfile.TemporaryDirectory() as directory:
                process = subprocess.run([str(self.binary), *arguments], cwd=directory,
                    capture_output=True, text=True, timeout=10)
                self.assertNotEqual(process.returncode, 0)

    def test_near_coincident_triangle_corner_regression(self):
        points = [
            [.3324177715466269, .27437811140650886],
            [.9999999700466811, 2.995331891719523e-8],
            [.9999999993624066, 6.375934034053183e-10],
            [.9999999700466811, 2.97097808167532e-8],
            [.9999999700466811, 1.7472769379504406e-8],
            [7/12, 0], [1/3, 2/3], [2/3, 1/3], [1/12, 1/4],
            [1/4, 0], [1/6, 5/6], [.5, .5], [5/6, 1/6],
            [1/6, 5/12], [5/12, 1/12], [.5, .25]]
        record = self.execute(supplied={'geometry': 'triangle', 'points': points})
        self.assertEqual(len(record['cells']), 16)
        self.assertTrue(all(record['cells']))
        # A dense lattice gives an independent lower bound even when Qhull's
        # handling of nearly coincident sites depends on the installed version.
        lattice = np.array([(i/40, j/40) for i in range(41) for j in range(41-i)])
        sampled = np.linalg.norm(lattice[:, None, :]-np.array(points)[None, :, :], axis=2).min(axis=1).max()
        self.assertLessEqual(sampled, record['covering_radius']+1e-10)


if __name__ == '__main__':
    unittest.main()

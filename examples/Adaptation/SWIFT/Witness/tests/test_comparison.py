# Copyright (c) 2026 Rodin contributors.
# Distributed under the Boost Software License, Version 1.0.
"""Independent audit of saved experimental comparison results, when supplied."""

import json
import os
from pathlib import Path
import unittest

import numpy as np

from reference import Reference


class ComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.environ.get('SWIFT_WITNESS_COMPARISON_RESULTS')
        if not path:
            raise unittest.SkipTest('Set SWIFT_WITNESS_COMPARISON_RESULTS to audit a completed comparison.')
        cls.records = json.loads(Path(path).read_text())

    def test_complete_grid_and_common_initialization(self):
        self.assertEqual(len(self.records), 24)
        keys = {(r['geometry'], r['method']) for r in self.records}
        self.assertEqual(len(keys), 24)
        for name in Reference.VERTICES:
            records = [r for r in self.records if r['geometry'] == name]
            self.assertEqual({r['method'] for r in records},
                             {'miniball', 'smooth_newton', 'current_search'})
            for record in records:
                self.assertEqual(record['count'], 1 if name == 'point' else 16)
                np.testing.assert_array_equal(record['initial_points'], records[0]['initial_points'])

    def test_independent_continuous_coverage(self):
        for record in self.records:
            with self.subTest(geometry=record['geometry'], method=record['method']):
                reference = Reference(record['geometry'])
                radius = reference.radius(np.array(record['points']))
                self.assertAlmostEqual(record['covering_radius'], radius, places=9)
                self.assertLessEqual(radius, record['initial_radius']+1e-10)
                self.assertLessEqual(record['validation_lattice_lower_bound'], radius+1e-10)

    def test_refinement_history_monotonicity(self):
        for record in self.records:
            if record['method'] == 'current_search':
                continue
            radii = [entry['radius'] for entry in record['history']]
            self.assertTrue(all(b <= a+1e-10 for a, b in zip(radii, radii[1:])))
            self.assertAlmostEqual(radii[-1], record['covering_radius'], places=12)
            self.assertLessEqual(record['search']['sweeps'], 100)

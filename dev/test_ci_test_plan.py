"""Preserve configuration coverage while changing test scheduling groups."""

import importlib.util
import itertools
from pathlib import Path
import re
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent))
from ci_test_plan import plan, targeted_plan, native_convergence_groups


class TestPlanTest(unittest.TestCase):
    def test_incremental_publications_fit_the_service_limit(self):
        import yaml
        from ci_tests import ROOT
        import ci_build_graph as graph
        import ci_schedule
        ordinary, targeted = ci_schedule.requests()
        by_id = {r['id']: f for f, rows in ordinary.items() for r in rows}
        names = []
        for entry in graph.prepare(ordinary, targeted):
            if entry['definition']['cmake']['RODIN_CODE_COVERAGE'] != 'ON':
                continue
            families = {by_id[r['matrix']['id']] if r['job'] == 'Ordinary' else r['job'] for r in entry['requests']}
            names += ['ci-batch-' + f + '-' + entry['id'] for f in families]
        self.assertGreater(len(names), 17)
        self.assertLess(len(names), 150)
        self.assertEqual(len(names), len(set(names)))
        self.assertTrue(all(len(name.encode('utf-8')) <= 100 for name in names))
        action = yaml.safe_load((ROOT / '.github/actions/coverage/action.yml').read_text())
        self.assertFalse(any(s.get('uses', '').startswith('codecov/') for s in action['runs']['steps']))
        ordinary = yaml.safe_load((ROOT / '.github/workflows/TestGroup.yml').read_text())
        self.assertFalse(any(s.get('uses', '').startswith('codecov/') for s in ordinary['jobs']['Run']['steps']))

    def test_configuration_matrix_and_job_limits(self):
        matrices = plan()
        rows = [row for entries in matrices.values() for row in entries]
        self.assertEqual(len({row["id"] for row in rows}), len(rows))
        for entries in matrices.values():
            self.assertGreater(len(entries), 0)
            self.assertLessEqual(len(entries), 256)
        required = set(itertools.product(["Debug", "Release"], ["ON", "OFF"], ["ON", "OFF"], ["ON", "OFF"]))
        for profile in ["unit", "manufactured"]:
            selected = [row for row in rows if row["profile"] == profile and not row["baseline"]]
            for group in {row["group"] for row in selected}:
                actual = {(r["build_type"], r["multithreaded"], r["asan"], r["ubsan"])
                          for r in selected if r["group"] == group and r["backend"] != "baseline"}
                self.assertEqual(actual, required)
        baseline = [r for r in rows if r["baseline"] and r["group"] == "src/Rodin"]
        self.assertEqual(len(baseline), 10)
        self.assertNotIn("coverage", matrices)
        self.assertFalse(any(r["profile"] == "coverage" and not r["baseline"] for r in rows))
        self.assertTrue(any(r["group"] == "cmake/ThirdPartyIsolation" for r in rows))

    def test_old_coverage_configuration_moves_to_ordinary_lanes(self):
        rows = [r for entries in plan().values() for r in entries]
        for profile in ["unit", "manufactured"]:
            normal = [r for r in rows if r["profile"] == profile and not r["baseline"]]
            for group in {r["group"] for r in normal if not r["group"].startswith("cmake/")}:
                baseline = [r for r in normal if r["group"] == group and r["backend"] == "baseline"]
                self.assertEqual({(r["build_type"], r["multithreaded"], r["asan"], r["ubsan"])
                                  for r in baseline}, {("Release", mt, "OFF", "OFF") for mt in ["ON", "OFF"]})

    def test_every_targeted_execution_has_independent_collection(self):
        import yaml
        from ci_tests import ROOT
        jobs = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]
        native = [{"group": group, "multithreaded": mt} for group in native_convergence_groups()
                  for mt in ["ON", "OFF"]]
        definitions = targeted_plan(native)
        self.assertNotIn("Coverage", jobs)
        for name, job in jobs.items():
            steps = job.get("steps", [])
            if not any("cmake -S . -B build" in s.get("run", "") for s in steps):
                continue
            self.assertTrue(any(d["job"] == name for d in definitions))
            for step in steps:
                command = step.get("run", "")
                if "cmake -S . -B build" in command:
                    self.assertIn("-DRODIN_CODE_COVERAGE=ON", command)
                if "apt-get install" in command:
                    self.assertIn("lcov", command)
            publisher = next(s for s in steps if s.get("id") == "coverage_record")
            self.assertEqual(publisher["uses"], "./.github/actions/coverage")
        action = yaml.safe_load((ROOT / ".github/actions/coverage/action.yml").read_text())
        ids = [s.get("id") for s in action["runs"]["steps"]]
        self.assertLess(ids.index("coverage_capture"), ids.index("coverage_receipt"))
        for name, job in jobs.items():
            if not any(d["job"] == name for d in definitions):
                continue
            self.assertEqual(sum(s.get("uses") == "./.github/actions/test-build" for s in job["steps"]), 1)
            import ci_schedule
            dag = ci_schedule.configuration_workflow()['jobs']
            _, executions = ci_schedule.layout()
            members = [k for k, e in executions.items() if e['job'] == name]
            self.assertTrue(members)
            self.assertTrue(set(members) <= set(dag['Publish_' + name]['needs']))

    def test_native_restriction_rank_partitions_preserve_every_invocation(self):
        import re
        import yaml
        from ci_tests import ROOT
        job = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]["MPIH1Restriction"]
        rows = job["strategy"]["matrix"]["include"]
        command = next(s["run"] for s in job["steps"] if s.get("name") == "Run native MPI restriction rank partition")
        self.assertIn("-j 1", command)
        self.assertEqual(job["timeout-minutes"],
                         "${{ matrix.partition == 'pyramid' && matrix.np == '8' && 180 || 90 }}")
        self.assertEqual(job["env"]["OMP_NUM_THREADS"], 2)
        self.assertEqual(job["env"]["OMP_WAIT_POLICY"], "PASSIVE")
        self.assertEqual(job["env"]["OMPI_MCA_hwloc_base_binding_policy"], "none")
        self.assertEqual(job["env"]["OMPI_MCA_mpi_yield_when_idle"], "1")
        self.assertTrue(any(s.get("name") == "Record native MPI CPU placement" for s in job["steps"]))
        geometries = ["Segment", "Triangle", "Quadrilateral", "Tetrahedron", "Pyramid", "Hexahedron", "Wedge"]
        expected = {f"RodinMPIH1SubMeshRestrictionTest_{geometry}_K{degree}_{mapping}_{value_range}_np{rank}"
                    for geometry, degree, mapping, value_range, rank in itertools.product(
                        geometries, range(1, 7), ["Affine", "Curved"],
                        ["RealScalar", "ComplexScalar", "RealVector", "ComplexVector", "RealMatrix", "ComplexMatrix"],
                        [1, 2, 3, 4, 8])}
        for mt in ["ON", "OFF"]:
            selected = []
            partitions = [r for r in rows if r["multithreaded"] == mt]
            self.assertEqual({r["np"] for r in partitions}, {"1", "2", "3", "4", "8"})
            for row in partitions:
                pattern = command.split("-R '", 1)[1].split("'", 1)[0].replace("${{ matrix.np }}", row["np"]).replace("${{ matrix.geometry }}", row["geometry"])
                selected.extend(name for name in expected if re.fullmatch(pattern, name))
            self.assertEqual(len(selected), 2520)
            self.assertEqual(set(selected), expected)

    def test_curved_nonlinear_registration_covers_each_case_once(self):
        import re
        from ci_tests import ROOT
        directory = ROOT / "tests/convergence/isoparametric/NonlinearPoisson"
        source = (directory / "NonlinearPoisson.cpp").read_text()
        expected = {f"{suite}.{case}" for suite, case in
                    re.findall(r"TEST_P\(Local(Test|Q[13]Test), (\w+)\)", source)}
        registration = (directory / "CMakeLists.txt").read_text()
        block = registration.split("set(_RODIN_CURVED_NONLINEAR_POISSON_CASES", 1)[1].split(")", 1)[0]
        cases = block.split()
        self.assertEqual(len(cases), 24)
        self.assertEqual(len(cases), len(set(cases)))
        self.assertEqual(set(cases), expected)
        self.assertIn("${suite}${entry}/${geometry}", registration)

    def test_stokes_split_preserves_both_orders_at_every_hierarchy_point(self):
        from ci_tests import ROOT
        for directory, suites in [("Stokes", ["StokesStabilityTest"]),
                                  ("PETScStokes", ["LocalStabilityTest", "MPIStabilityTest"])]:
            source = (ROOT / f"tests/convergence/h/{directory}/Stability.cpp").read_text()
            for mapping in ["Affine", "Curved"]:
                for resolution in [3, 4, 5]:
                    for suite in suites:
                        self.assertEqual(source.count(f"TEST_P({suite}, P4P3PressureSpectrum{mapping}N{resolution})"), 1)
                    self.assertEqual(source.count(f"checkHigherOrderHierarchy<4>({resolution}, {str(mapping == 'Curved').lower()});"), len(suites))
            self.assertIn("i == 0 ? AssemblyOrder : RefinedOrder", source)

    def test_degree_nonlinear_registration_covers_every_independent_case(self):
        from ci_tests import ROOT
        directory = ROOT / "tests/convergence/p/NonlinearPoisson"
        expected = set(re.findall(r"TEST_P\(NonlinearPoissonPTest, (\w+)\)",
                                  (directory / "NonlinearPoisson.cpp").read_text()))
        registration = (directory / "CMakeLists.txt").read_text()
        cases = registration.split("foreach(case", 1)[1].split(")", 1)[0].split()
        self.assertEqual(len(cases), len(set(cases)))
        self.assertEqual(set(cases), expected)

    def test_petsc_pressure_partitions_select_every_split_case_once(self):
        import yaml
        from ci_tests import ROOT
        workflow = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())
        job = workflow["jobs"]["PETScConvergence"]
        command = next(step["run"] for step in job["steps"]
                       if "geometry_pattern" in step.get("run", ""))
        patterns = re.findall(r'pattern="([^"\n]+StokesStability_[^"\n]+)"', command)
        self.assertEqual(len(patterns), 2)
        geometries = ["Triangle", "Quadrilateral", "Tetrahedron", "Pyramid", "Hexahedron", "Wedge"]
        points = [f"Rates_{mapping}N{n}" for mapping in ["Affine", "Curved"] for n in [3, 4, 5]] + ["Controls"]
        expected = set()
        for rank in [None, 1, 2, 3, 4]:
            prefix = "RodinConvergenceHPETSc" + ("MPI" if rank else "") + "StokesStability_"
            middle = f"np{rank}_" if rank else ""
            expected.update(prefix + middle + geometry for geometry in geometries)
            expected.update(prefix + "P4P3_" + middle + point + "_" + geometry
                            for point in points for geometry in geometries)
            expected.add(prefix + "Utility" + (f"_np{rank}" if rank else ""))
        rows = job["strategy"]["matrix"]["include"]
        for mt in ["ON", "OFF"]:
            selected = []
            for row in rows:
                if row["multithreaded"] != mt or not row["suite"].startswith("pressure_"):
                    continue
                local = row["context"] == "local"
                geometry = {"pressure_light": "(Triangle|Quadrilateral|Hexahedron|Wedge)",
                            "pressure_tetrahedron": "Tetrahedron", "pressure_pyramid": "Pyramid"}[row["suite"]]
                rank = row["context"].removeprefix("np")
                utility = ("|Utility" + ("" if local else "_np" + rank)) if row["suite"] == "pressure_light" else ""
                pattern = patterns[0 if local else 1].replace("${geometry_pattern}", geometry).replace("${utility_pattern}", utility).replace("${rank}", rank)
                selected.extend(name for name in expected if re.fullmatch(pattern, name))
            self.assertEqual(len(selected), len(expected))
            self.assertEqual(set(selected), expected)

    def test_curved_boundary_selector_keeps_degree_cases_and_rank_partitions(self):
        import os
        import subprocess
        import tempfile
        from collections import Counter
        import yaml
        from ci_tests import ROOT
        job = yaml.safe_load((ROOT / '.github/ci/targeted-lanes.yml').read_text())['jobs']['PETScCurvedBoundaryConvergence']
        command = next(s['run'] for s in job['steps'] if s.get('name') == 'Run curved boundary convergence partition')
        geometries = ['Segment', 'Triangle', 'Quadrilateral', 'Tetrahedron', 'Pyramid', 'Hexahedron', 'Wedge']
        degrees = {'LinearElasticity': ['MixedTractionP' + str(k) + 'Rates' for k in [1, 2, 3]] + ['Controls'],
                   'NonlinearPoisson': [b + 'P' + str(k) + 'Rates' for b in ['MixedNeumann', 'Robin', 'PureNeumann']
                                        for k in [1, 2, 3]] + ['Controls']}
        with tempfile.TemporaryDirectory() as temporary:
            ctest = Path(temporary) / 'ctest'
            ctest.write_text('#!/bin/sh\nwhile [ "$1" != "-R" ]; do shift; done\nprintf "%s\\n" "$2"\n')
            ctest.chmod(0o755)
            counts = {mt: Counter() for mt in ['ON', 'OFF']}
            for row in job['strategy']['matrix']['include']:
                script = command
                for name, value in row.items():
                    script = script.replace('${{ matrix.' + name + ' }}', value)
                pattern = subprocess.check_output(['bash', '-eu', '-c', script], text=True,
                                                  env=dict(os.environ, PATH=temporary + ':' + os.environ['PATH'])).strip()
                shapes = [g for g in geometries if g not in ['Tetrahedron', 'Pyramid']] if row['geometry'] == 'light' else [row['geometry'].capitalize()]
                # Names reflect PETScRefinement.cmake's RATE_CASES and Controls
                # registrations. Local selection must not capture an MPI rank.
                ranks = [None, 1, 2, 3, 4] if row['selection'] == 'all' else [
                    None if row['selection'] == 'local' else int(row['selection'].removeprefix('np'))]
                prefix = 'RodinConvergenceIsoparametricPETSc' + row['physics'] + 'Boundary_'
                names = {prefix + ('' if rank is None else 'MPI_np' + str(rank) + '_') +
                         (case + '_' if case else '') + geometry
                         for rank in [None, 1, 2, 3, 4] for geometry in geometries
                         for case in degrees.get(row['physics'], [''])}
                expected = {name for name in names if any(name.endswith('_' + shape) for shape in shapes)
                            and any(name.startswith(prefix + ('MPI_np' + str(rank) + '_' if rank else ''))
                                    and (rank is not None or not name.startswith(prefix + 'MPI_')) for rank in ranks)}
                if row['study'] != 'all':
                    expected = {name for name in expected if '_' + row['study'] + '_' in name}
                actual = {name for name in names if re.fullmatch(pattern, name)}
                self.assertTrue(expected)
                self.assertEqual(actual, expected, row)
                counts[row['multithreaded']].update(actual)
            self.assertTrue(all(count == 1 for c in counts.values() for count in c.values()))


if __name__ == "__main__":
    unittest.main()

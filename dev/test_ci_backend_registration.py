"""Check optional-backend headers and concrete CMake test registrations."""

from collections import Counter
import importlib.util
import itertools
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("ci_tests", ROOT / "dev/ci_tests.py")
ci = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ci)
GEOMETRIES = ["Segment", "Triangle", "Quadrilateral", "Tetrahedron", "Pyramid", "Hexahedron", "Wedge"]


class BackendRegistrationTest(unittest.TestCase):
    def configure(self, directory, petsc=True, complex=False, coverage=False):
        temporary = tempfile.TemporaryDirectory(prefix="rodin-ci-registration-")
        self.addCleanup(temporary.cleanup)
        source = Path(temporary.name)
        targets = ["GTest::gtest", "GTest::gtest_main", "Rodin::Geometry", "Rodin::Variational",
                   "Rodin::MPI", "Rodin::PETSc", "Rodin::Solver", "Rodin::Assembly",
                   "RodinConvergence", "RodinConvergenceH"]
        directories = [directory] if isinstance(directory, str) else directory
        folders = ["cases"] if len(directories) == 1 else [f"cases{i}" for i in range(len(directories))]
        (source / "CMakeLists.txt").write_text(
            'cmake_minimum_required(VERSION 3.16)\nproject(Registration LANGUAGES CXX)\nenable_testing()\n'
            + ''.join(f'add_library({target} INTERFACE IMPORTED)\n' for target in targets)
            + 'function(rodin_suppress_external_mpi_lsan_for_test name)\nendfunction()\n'
            + f'set(RODIN_USE_PETSC {"ON" if petsc else "OFF"})\nset(RODIN_USE_MPI ON)\n'
            + f'set(RODIN_CONVERGENCE_PETSC_COMPLEX {"ON" if complex else "OFF"})\n'
            + f'set(RODIN_CODE_COVERAGE {"ON" if coverage else "OFF"})\n'
            + 'set(MPIEXEC_EXECUTABLE /usr/bin/true)\nset(MPIEXEC_NUMPROC_FLAG -n)\n'
            + f'include("{ROOT / "tests/convergence/PETScRefinement.cmake"}")\n'
            + ''.join(f'add_subdirectory("{ROOT / entry}" {folder})\n' for entry, folder in zip(directories, folders)))
        build = source / "build"
        subprocess.run(["cmake", "-S", str(source), "-B", str(build)],
                       check=True, capture_output=True, text=True)
        tests = {}
        for folder in folders:
            commands = ci.commands((build / folder / "CTestTestfile.cmake").read_text())
            tests.update({args[0]: args[1:] for command, args, *_ in commands if command == "add_test"})
        return tests, build

    def test_native_hp_stokes_cases_have_independent_time_budgets(self):
        tests, build = self.configure("tests/convergence/hp/Stokes")
        source = (ROOT / "tests/convergence/hp/Stokes/Stokes.cpp").read_text()
        cases = re.findall(r'TEST_P\(StokesHPTest,\s*(\w+)\)', source)
        expected = Counter(f"AllGeometries/StokesHPTest.{case}/{geometry}"
                           for case, geometry in itertools.product(cases, GEOMETRIES[1:]))
        actual = Counter()
        for command in tests.values():
            self.assertTrue(command[0].endswith("RodinConvergenceHPStokes"))
            filters = [arg.removeprefix("--gtest_filter=") for arg in command
                       if arg.startswith("--gtest_filter=")]
            self.assertEqual(len(filters), 1)
            actual[filters[0]] += 1
        self.assertEqual(actual, expected)
        properties = {}
        for command, args, *_ in ci.commands((build / "cases/CTestTestfile.cmake").read_text()):
            if command == "set_tests_properties":
                boundary = args.index("PROPERTIES")
                row = dict(zip(args[boundary + 1::2], args[boundary + 2::2]))
                for name in args[:boundary]:
                    properties[name] = row
        self.assertEqual(set(properties), set(tests))
        for row in properties.values():
            self.assertEqual(row["TIMEOUT"], "600")
            self.assertEqual(set(row["LABELS"].split(";")), {"convergence", "slow"})

    def test_matrix_range_target_and_all_ranks_are_registered(self):
        for petsc in [False, True]:
            tests, build = self.configure("tests/unit/Rodin/MPI/Variational", petsc)
            for rank in [1, 2, 3, 4]:
                command = tests[f"RodinMPIMatrixRangeTest_np{rank}"]
                self.assertIn(str(rank), command)
                self.assertTrue(any(arg.endswith("RodinMPIMatrixRangeTest") for arg in command))
            flags = (build / "cases/CMakeFiles/RodinMPIMatrixRangeTest.dir/flags.make").read_text()
            self.assertEqual("-DRODIN_TEST_WITH_PETSC" in flags, petsc)

    def test_mpi_assembly_header_is_inert_when_backend_is_disabled(self):
        compiler = shutil.which("c++")
        if not compiler:
            self.skipTest("C++ compiler unavailable")
        with tempfile.TemporaryDirectory(prefix="rodin-no-mpi-header-") as directory:
            source = Path(directory)
            (source / "Rodin").mkdir()
            (source / "Rodin/Configure.h").write_text('// MPI disabled\n')
            (source / "main.cpp").write_text(
                '#include "Rodin/PETSc/Assembly/MPI.h"\nint main() { return 0; }\n')
            # No MPI, PETSc or Boost include/library paths: a disabled backend
            # must not instantiate archive registration or require their symbols.
            subprocess.run([compiler, "-std=c++20", "-I" + str(source),
                            "-I" + str(ROOT / "src"), str(source / "main.cpp"),
                            "-o", str(source / "header")], check=True, capture_output=True, text=True)

    def test_curved_elasticity_boundary_enables_curved_workload(self):
        _, build = self.configure("tests/convergence/isoparametric/LinearElasticity")
        target = "RodinConvergenceIsoparametricPETScLinearElasticityBoundary"
        flags = (build / f"cases/CMakeFiles/{target}.dir/flags.make").read_text()
        self.assertIn("-DRODIN_ELASTICITY_BOUNDARY_CURVED", flags)

    def test_curved_traction_budgets_preserve_studies_and_ranks(self):
        tests, build = self.configure("tests/convergence/isoparametric/LinearElasticity")
        prefix = "RodinConvergenceIsoparametricPETScLinearElasticityBoundary"
        properties = {}
        for command, args, *_ in ci.commands((build / "cases/CTestTestfile.cmake").read_text()):
            if command == "set_tests_properties":
                boundary = args.index("PROPERTIES")
                values = dict(zip(args[boundary + 1::2], args[boundary + 2::2]))
                for name in args[:boundary]:
                    properties.setdefault(name, {}).update(values)
        for rank, degree, geometry in itertools.product(range(5), range(1, 4), GEOMETRIES):
            suffix = f"_MPI_np{rank}" if rank else ""
            name = f"{prefix}{suffix}_MixedTractionP{degree}Rates_{geometry}"
            command = tests[name]
            self.assertIn("--gtest_fail_if_no_test_selected", command)
            suite = "MPITest" if rank else "LocalTest"
            self.assertIn(f"--gtest_filter=AllGeometries/{suite}.MixedTractionP{degree}Rates/{geometry}", command)
            expected = 14400 if degree == 1 else 7200
            self.assertEqual(properties[name]["TIMEOUT"], str(expected if geometry == "Tetrahedron" else 1800))
            if rank:
                self.assertIn(str(rank), command)
                self.assertEqual(properties[name]["PROCESSORS"], str(rank))
        job = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]["PETScCurvedBoundaryConvergence"]
        self.assertIn("(matrix.physics == 'LinearElasticity' && matrix.geometry == 'tetrahedron') && 360", job["timeout-minutes"])

    def test_curved_elasticity_studies_have_exactly_one_process_each(self):
        tests, _ = self.configure("tests/convergence/isoparametric/LinearElasticity")
        text = (ROOT / "tests/convergence/isoparametric/LinearElasticity/LinearElasticity.cpp").read_text()
        cases = re.findall(r'TEST_P\((Local(?:Q[13])?Test|MPI(?:Q[13])?Test),\s*(\w+)\)', text)
        expected = Counter()
        for suite, case in cases:
            ranks = [1, 2, 3, 4] if suite.startswith("MPI") else [0, 0]
            for geometry, rank in itertools.product(GEOMETRIES, ranks):
                expected[(f"AllGeometries/{suite}.{case}/{geometry}", rank)] += 1
        actual = Counter()
        for name, command in tests.items():
            if "Boundary" in name:
                continue
            filter_ = next(a.removeprefix("--gtest_filter=") for a in command if a.startswith("--gtest_filter="))
            rank = int(command[command.index("-n") + 1]) if command[0] == "/usr/bin/true" else 0
            actual[(filter_, rank)] += 1
        self.assertEqual(actual, expected)

    def test_curved_job_partitions_select_every_study_once(self):
        job = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]["PETScCurvedConvergence"]
        step = next(s["run"] for s in job["steps"] if s.get("name") == "Run curved PETSc convergence suites")
        tests, _ = self.configure("tests/convergence/isoparametric/LinearElasticity")
        names = [n for n in tests if "PETSc" in n and "Boundary" not in n]
        source = (ROOT / "tests/convergence/isoparametric/NonlinearPoisson/CMakeLists.txt").read_text()
        entries = re.search(r'set\(_RODIN_CURVED_NONLINEAR_POISSON_CASES\s+(.*?)\)', source, re.S)[1].split()
        for geometry, rank in itertools.product(GEOMETRIES, [0, 1, 2, 3, 4]):
            prefix = "RodinConvergenceIsoparametricPETSc" + ("MPI" if rank else "")
            context = f"np{rank}_" if rank else ""
            names.extend(f"{prefix}NonlinearPoisson_{context}{case.replace('.', '_')}_{geometry}" for case in entries)
        counts = Counter()
        with tempfile.TemporaryDirectory(prefix="rodin-ci-selectors-") as directory:
            ctest = Path(directory) / "ctest"
            ctest.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
            ctest.chmod(0o755)
            environment = dict(os.environ, PATH=directory + os.pathsep + os.environ["PATH"])
            for row in job["strategy"]["matrix"]["include"]:
                if row["multithreaded"] != "ON":
                    continue
                command = re.sub(r'\$\{\{\s*matrix\.(\w+)\s*\}\}', lambda m: str(row.get(m[1], "")), step)
                args = subprocess.run(["bash", "-ec", command], env=environment, check=True,
                                      capture_output=True, text=True).stdout.splitlines()
                include = re.compile(args[args.index("-R") + 1])
                exclude = re.compile(args[args.index("-E") + 1])
                for name in names:
                    if include.search(name) and not exclude.search(name):
                        counts[name] += 1
        self.assertEqual(counts, Counter({name: 1 for name in names}))

    def test_degree_job_partitions_select_every_registered_study_once(self):
        physics = ["Poisson", "Conductivity", "NonlinearPoisson", "ReactionDiffusion", "LinearElasticity", "Stokes"]
        tests, _ = self.configure([f"tests/convergence/{family}/PETSc{p}" for family, p in itertools.product(["p", "hp"], physics)])
        job = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]["PETScDegreeConvergence"]
        step = next(s["run"] for s in job["steps"] if s.get("name") == "Run PETSc degree-refinement suites")
        counts = Counter()
        with tempfile.TemporaryDirectory(prefix="rodin-ci-degree-selectors-") as directory:
            ctest = Path(directory) / "ctest"
            ctest.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
            ctest.chmod(0o755)
            environment = dict(os.environ, PATH=directory + os.pathsep + os.environ["PATH"])
            for row in job["strategy"]["matrix"]["include"]:
                if row["multithreaded"] != "ON":
                    continue
                command = re.sub(r'\$\{\{\s*matrix\.(\w+)\s*\}\}', lambda m: str(row.get(m[1], "")), step)
                args = subprocess.run(["bash", "-ec", command], env=environment, check=True,
                                      capture_output=True, text=True).stdout.splitlines()
                include = re.compile(args[args.index("-R") + 1])
                for name in tests:
                    if include.search(name):
                        counts[name] += 1
        self.assertEqual(counts, Counter({name: 1 for name in tests}))

    def test_boundary_job_partitions_select_every_registered_study_once(self):
        modules = ["Poisson", "Conductivity", "LinearElasticity", "Stokes", "ReactionDiffusion", "NonlinearPoisson"]
        tests, _ = self.configure([f"tests/convergence/h/PETSc{module}" for module in modules])
        names = [name for name in tests if re.match(r"^RodinConvergenceHPETSc(?:" + "|".join(modules) + r")Boundary_", name)]
        job = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]["PETScConvergence"]
        step = next(s["run"] for s in job["steps"] if s.get("name") == "Run PETSc boundary convergence suites")
        counts = Counter()
        with tempfile.TemporaryDirectory(prefix="rodin-ci-traction-selectors-") as directory:
            ctest = Path(directory) / "ctest"
            ctest.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
            ctest.chmod(0o755)
            environment = dict(os.environ, PATH=directory + os.pathsep + os.environ["PATH"])
            for row in job["strategy"]["matrix"]["include"]:
                if row["multithreaded"] != "ON" or not row.get("physics"):
                    continue
                command = re.sub(r'\$\{\{\s*matrix\.(\w+)\s*\}\}', lambda m: str(row.get(m[1], "")), step)
                args = subprocess.run(["bash", "-ec", command], env=environment, check=True,
                                      capture_output=True, text=True).stdout.splitlines()
                include = re.compile(args[args.index("-R") + 1])
                for name in names:
                    if include.search(name):
                        counts[name] += 1
        self.assertEqual(counts, Counter({name: 1 for name in names}))

    def test_aggregate_petsc_lanes_partition_registered_tests_without_overlap(self):
        jobs = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]
        modules = ["Poisson", "MPIPoisson", "Conductivity", "LinearElasticity",
                   "ReactionDiffusion", "NonlinearPoisson", "Stokes"]
        real, _ = self.configure([f"tests/convergence/h/PETSc{m}" for m in modules])
        real = {n for n in real if not re.search(r"Boundary_|StokesStability_", n)}
        complex_tests, _ = self.configure("tests/convergence/isoparametric/Helmholtz", complex=True)
        complex_tests = {n for n in complex_tests if n.startswith("RodinConvergenceIsoparametricPETScHelmholtzBoundary_")}
        for family, suite, step_name, names in [
            ("PETScConvergence", "bulk", "Run PETSc convergence suite", real),
            ("PETScComplexConvergence", "curved_boundary", "Run curved complex PETSc Helmholtz boundary suite", complex_tests),
        ]:
            self.assertTrue(names)
            job = jobs[family]
            step = next(s["run"] for s in job["steps"] if s.get("name") == step_name)
            with tempfile.TemporaryDirectory(prefix="rodin-ci-context-selectors-") as directory:
                ctest = Path(directory) / "ctest"
                ctest.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
                ctest.chmod(0o755)
                for mt in ["ON", "OFF"]:
                    counts = Counter()
                    for row in job["strategy"]["matrix"]["include"]:
                        if row["multithreaded"] != mt or row["suite"] != suite:
                            continue
                        command = re.sub(r'\$\{\{\s*matrix\.(\w+)\s*\}\}', lambda m: str(row.get(m[1], "")), step)
                        args = subprocess.check_output(["bash", "-ec", command], text=True,
                            env=dict(os.environ, PATH=directory + os.pathsep + os.environ["PATH"])).splitlines()
                        include = re.compile(args[args.index("-R") + 1])
                        exclude = re.compile(args[args.index("-E") + 1]) if "-E" in args else None
                        selected = {n for n in names if include.search(n) and not (exclude and exclude.search(n))}
                        self.assertTrue(selected, row)
                        counts.update(selected)
                    self.assertEqual(counts, Counter({n: 1 for n in names}))

    def test_curved_boundary_lanes_cover_configured_names_including_controls(self):
        modules = ["Diffusion", "ReactionDiffusion", "LinearElasticity", "NonlinearPoisson", "Stokes"]
        configured, _ = self.configure(["tests/convergence/isoparametric/" + m for m in modules])
        names = {n for n in configured if re.match(r"^RodinConvergenceIsoparametricPETSc[A-Za-z]+Boundary_", n)}
        job = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]["PETScCurvedBoundaryConvergence"]
        step = next(s["run"] for s in job["steps"] if s.get("name") == "Run curved boundary convergence partition")
        # The dedicated tetrahedron family owns nonlinear tetrahedra.
        names = {n for n in names if not ("NonlinearPoissonBoundary_" in n and n.endswith("_Tetrahedron"))}
        with tempfile.TemporaryDirectory(prefix="rodin-ci-curved-boundary-selectors-") as directory:
            ctest = Path(directory) / "ctest"
            ctest.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
            ctest.chmod(0o755)
            for mt in ["ON", "OFF"]:
                counts = Counter()
                for row in job["strategy"]["matrix"]["include"]:
                    if row["multithreaded"] != mt:
                        continue
                    command = re.sub(r'\$\{\{\s*matrix\.(\w+)\s*\}\}', lambda m: str(row.get(m[1], "")), step)
                    args = subprocess.check_output(["bash", "-ec", command], text=True,
                        env=dict(os.environ, PATH=directory + os.pathsep + os.environ["PATH"])).splitlines()
                    include = re.compile(args[args.index("-R") + 1])
                    selected = {n for n in names if include.search(n)}
                    self.assertTrue(selected, row)
                    counts.update(selected)
                self.assertEqual(counts, Counter({n: 1 for n in names}))

    def test_petsc_timeout_overrides_preserve_commands_ranks_and_other_budgets(self):
        modules = ["p/PETScNonlinearPoisson", "hp/PETScNonlinearPoisson",
                   "h/PETScNonlinearPoisson", "h/PETScStokes",
                   "isoparametric/NonlinearPoisson", "isoparametric/ReactionDiffusion"]
        tests, build = self.configure(["tests/convergence/" + m for m in modules])
        properties = {}
        for path in build.rglob("CTestTestfile.cmake"):
            for command, args, *_ in ci.commands(path.read_text()):
                if command != "set_tests_properties":
                    continue
                boundary = args.index("PROPERTIES")
                values = dict(zip(args[boundary + 1::2], args[boundary + 2::2]))
                for name in args[:boundary]:
                    properties.setdefault(name, {}).update(values)
        checked = 0
        for name, args in tests.items():
            expected = None
            if re.match(r"RodinConvergence(P|HP)PETScNonlinearPoisson_", name):
                expected = 10800 if name.endswith("_Pyramid") else 1800
            elif re.match(r"RodinConvergenceHPETSc(MPI)?NonlinearPoisson(_np[1-4])?_", name):
                expected = 3600 if name.endswith("_Pyramid") else 600
            elif name.startswith("RodinConvergenceHPETScStokesBoundary_"):
                expected = 1800
                if "_P3P2Rates_" in name and name.rsplit("_", 1)[-1] in ["Tetrahedron", "Pyramid", "Hexahedron", "Wedge"]:
                    expected = 7200 if name.endswith("_Pyramid") else 3600
            elif name.startswith("RodinConvergenceIsoparametricPETScReactionDiffusionBoundary_"):
                expected = 14400 if name.endswith("_Pyramid") else 7200 if name.endswith(("_Tetrahedron", "_Hexahedron")) else 1800
            elif name.startswith("RodinConvergenceIsoparametricPETScNonlinearPoissonBoundary_") and name.endswith("_Pyramid"):
                expected = 10800 if "Rates_" in name else 1800
            if expected is None:
                continue
            checked += 1
            self.assertEqual(int(properties[name]["TIMEOUT"]), expected, name)
            rank = re.search(r"(?:MPI_np|_np)([1-4])", name)
            if rank:
                self.assertIn(rank[1], args)
                self.assertEqual(properties[name]["PROCESSORS"], rank[1])
            filters = [a for a in args if a.startswith("--gtest_filter=")]
            self.assertEqual(len(filters), 1, name)
            self.assertIn("/" + name.rsplit("_", 1)[-1], filters[0])
        self.assertGreater(checked, 200)

    def test_curved_stokes_coverage_watchdog_preserves_every_case_and_rank(self):
        configurations = [self.configure("tests/convergence/isoparametric/Stokes", coverage=value)
                          for value in [False, True]]
        commands = [{name: [arg.replace(str(build.parent), "<FIXTURE>") for arg in args]
                     for name, args in tests.items()} for tests, build in configurations]
        self.assertEqual(commands[0], commands[1])
        properties = []
        for _, build in configurations:
            rows = {}
            for path in build.rglob("CTestTestfile.cmake"):
                for command, args, *_ in ci.commands(path.read_text()):
                    if command != "set_tests_properties":
                        continue
                    boundary = args.index("PROPERTIES")
                    values = {key: value.replace(str(build.parent), "<FIXTURE>")
                              for key, value in zip(args[boundary + 1::2], args[boundary + 2::2])}
                    for name in args[:boundary]:
                        rows.setdefault(name, {}).update(values)
            properties.append(rows)
        selected = "RodinConvergenceIsoparametricPETScStokesBoundary_P3P2Rates_Pyramid"
        self.assertEqual(properties[0][selected]["TIMEOUT"], "1800")
        self.assertEqual(properties[1][selected]["TIMEOUT"], "3600")
        properties[1][selected]["TIMEOUT"] = "1800"
        self.assertEqual(properties[0], properties[1])
        for rank in range(1, 5):
            name = f"RodinConvergenceIsoparametricPETScStokesBoundary_MPI_np{rank}_P3P2Rates_Pyramid"
            self.assertIn(str(rank), configurations[1][0][name])
            self.assertEqual(properties[1][name]["PROCESSORS"], str(rank))

    def test_expanded_jobs_fit_github_matrix_limit(self):
        jobs = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]
        for job in jobs.values():
            matrix = job.get("strategy", {}).get("matrix", {})
            if isinstance(matrix, dict):
                self.assertLessEqual(len(matrix.get("include", [])), 256)


if __name__ == "__main__":
    unittest.main()

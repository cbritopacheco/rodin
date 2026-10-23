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
    def configure(self, directory, petsc=True):
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
            + 'set(MPIEXEC_EXECUTABLE /usr/bin/true)\nset(MPIEXEC_NUMPROC_FLAG -n)\n'
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

    def test_expanded_jobs_fit_github_matrix_limit(self):
        jobs = yaml.safe_load((ROOT / ".github/ci/targeted-lanes.yml").read_text())["jobs"]
        for job in jobs.values():
            matrix = job.get("strategy", {}).get("matrix", {})
            if isinstance(matrix, dict):
                self.assertLessEqual(len(matrix.get("include", [])), 256)


if __name__ == "__main__":
    unittest.main()

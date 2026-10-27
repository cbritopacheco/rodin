"""Exercise scoped discovery with real CMake/CTest and unbuilt sibling targets."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

SPEC = importlib.util.spec_from_file_location("ci_tests", Path(__file__).with_name("ci_tests.py"))
ci = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ci)


class CMakeSyntaxTest(unittest.TestCase):
    def test_generated_arguments(self):
        text = '# comment (\nadd_test([=[a)b]=] "x\\\"y" [==[semi;colon\nline]==])'
        parsed = list(ci.commands(text))
        self.assertEqual(parsed[0][1], ['a)b', 'x"y', 'semi;colon\nline'])
        self.assertIn('==', ci.quote('x]=]y'))

    def test_unterminated_input(self):
        with self.assertRaises(ValueError):
            list(ci.commands('add_test("broken)'))


class DiscoveryCoverageTest(unittest.TestCase):
    def test_ctest_discovery_failure_preserves_diagnostic(self):
        failure = subprocess.CalledProcessError(8, ['ctest'], output='', stderr='discovery executable crashed')
        with patch.object(ci.subprocess, 'run', side_effect=failure):
            with self.assertRaisesRegex(RuntimeError, 'discovery executable crashed'):
                ci.discover(Path('/tmp/stage'), Path('/tmp/build'))

    def test_discovery_counters_are_outside_capture_tree(self):
        build = Path("/tmp/rodin-discovery-test/build")
        response = SimpleNamespace(stdout='{"tests": []}')
        with patch.object(ci.subprocess, "run", return_value=response) as run:
            ci.discover(build / "ci-tests", build)
        prefix = Path(run.call_args.kwargs["env"]["GCOV_PREFIX"])
        self.assertFalse(prefix.is_relative_to(build))


class ExecutionProfilesTest(unittest.TestCase):
    def test_dedicated_convergence_targets_are_not_required_in_baseline_plan(self):
        with tempfile.TemporaryDirectory() as directory:
            build = Path(directory)
            row = {'target': 'RodinConvergenceIsoparametricLinearElasticity', 'kind': 'test',
                   'owner': 'tests/convergence/isoparametric/LinearElasticity',
                   'executable': str(build / 'elasticity'), 'profiles': 'convergence-baseline'}
            ci.configured_profiles.cache_clear()
            self.addCleanup(ci.configured_profiles.cache_clear)
            with patch.object(ci, 'targets', return_value=[row]):
                result = ci.configured_profiles(build)
            self.assertIn('convergence-baseline', result[row['target']])
            ci.configured_profiles.cache_clear()
            with patch.object(ci, 'targets', return_value=[dict(row, target='UnscheduledConvergenceTarget')]):
                with self.assertRaisesRegex(ValueError, 'missing from the execution plan'):
                    ci.configured_profiles(build)

    def test_preserve_old_selections(self):
        def test(name, *labels):
            return {"name": name, "properties": [{"name": "LABELS", "value": list(labels)}]}
        self.assertTrue(ci.selected("unit", test("mpi", "unit", "distributed")))
        self.assertFalse(ci.selected("unit", test("slow", "unit", "slow")))
        self.assertTrue(ci.selected("unit", test("RodinMPIAssemblyGeometrySlowTest_np8_Wedge", "slow")))
        self.assertFalse(ci.selected("slow", test("RodinMPIH1SubMeshRestrictionTest_scalar_np8", "slow")))
        self.assertTrue(ci.selected("coverage", test("serial", "unit", "slow")))
        self.assertFalse(ci.selected("coverage", test("mpi", "unit", "slow", "distributed")))
        self.assertTrue(ci.selected("coverage", test("manufactured-mpi", "manufactured", "distributed")))
        self.assertFalse(ci.selected("manufactured", test("slow", "manufactured", "slow")))

    def test_configuration_coverage_eligibility_and_backend_contract(self):
        import itertools
        for profile in ["unit", "manufactured", "slow"]:
            for build_type, mt, asan, ubsan in itertools.product(
                    ["Debug", "Release"], ["ON", "OFF"], ["ON", "OFF"], ["ON", "OFF"]):
                with patch.dict(os.environ, {"SUITESPARSE_PREFIX": "/tmp/suitesparse"}), \
                        patch.object(ci.subprocess, "run") as run:
                    ci.configure(Path("build"), "unit/Rodin/Math", profile, build_type, mt, asan, ubsan)
                    options = dict(a[2:].split("=", 1) for a in run.call_args.args[0] if a.startswith("-D"))
                self.assertEqual(options["RODIN_CODE_COVERAGE"], "ON" if asan == ubsan == "OFF" else "OFF")
                self.assertEqual(options["CMAKE_BUILD_TYPE"], build_type)
                self.assertEqual(options["RODIN_MULTITHREADED"], mt)
                if profile == "unit":
                    self.assertEqual(options["RODIN_USE_PARU"], "ON")
                    self.assertEqual(options["RODIN_USE_MUMPS"], "OFF" if asan == "ON" else "ON")
                if profile == "slow":
                    self.assertEqual(options["RODIN_USE_SCOTCH"], "ON")
        for profile in ["unit", "manufactured"]:
            with patch.object(ci.subprocess, "run") as run:
                ci.configure(Path("build"), "unit/Rodin/Math", profile, "Release", "ON", "OFF", "OFF", "baseline")
                options = dict(a[2:].split("=", 1) for a in run.call_args.args[0] if a.startswith("-D"))
            for option in ["RODIN_USE_SCOTCH", "RODIN_USE_UMFPACK", "RODIN_USE_CHOLMOD", "RODIN_USE_SPQR"]:
                self.assertEqual(options[option], "ON")
            self.assertNotIn("RODIN_USE_PARU", options)


class ScopedCTestTest(unittest.TestCase):
    def test_unbuilt_discovery_and_invocation_preservation(self):
        with tempfile.TemporaryDirectory(prefix="rodin-ci-test-") as directory:
            root = Path(directory)
            source = root / "source"
            build = root / "build"
            owner = source / "tests/unit/Rodin/Math"
            owner.mkdir(parents=True)
            (owner / "main.cpp").write_text('int main() { return 0; }\n')
            (source / "CMakeLists.txt").write_text('''cmake_minimum_required(VERSION 3.16)
project(ScopedTests LANGUAGES CXX)
enable_testing()
add_subdirectory(tests/unit/Rodin/Math)
''')
            (owner / "CMakeLists.txt").write_text('''add_executable(Fast main.cpp)
add_executable(Unbuilt main.cpp)
add_test(NAME explicit COMMAND Fast "--rank=3" "semi;colon")
set_tests_properties(explicit PROPERTIES LABELS "unit;distributed" TIMEOUT 17
  ENVIRONMENT "RODIN_FIXTURE=ok" PROCESSORS 3 RESOURCE_LOCK mesh)
file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/RodinTestTargets-.tsv" CONTENT
"target\\towner\\texecutable\\tkind\\nFast\\ttests/unit/Rodin/Math\\t$<TARGET_FILE:Fast>\\ttest\\nUnbuilt\\ttests/unit/Rodin/Other\\t$<TARGET_FILE:Unbuilt>\\ttest\\n")
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/Fast[1]_include.cmake" CONTENT
"add_test(dynamic [==[$<TARGET_FILE:Fast>]==] --gtest_filter=Fixture.Dynamic)\\nset_tests_properties(dynamic PROPERTIES LABELS unit WORKING_DIRECTORY [==[${CMAKE_CURRENT_BINARY_DIR}]==] TIMEOUT 19)\\n")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/Unbuilt[1]_include.cmake"
  "message(FATAL_ERROR \\"Unbuilt discovery must not run\\")\\n")
set_property(DIRECTORY APPEND PROPERTY TEST_INCLUDE_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/Fast[1]_include.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/Unbuilt[1]_include.cmake")
''')
            subprocess.run(["cmake", "-S", str(source), "-B", str(build)], check=True, capture_output=True)
            subprocess.run(["cmake", "--build", str(build), "--target", "Fast", "--parallel", "1"],
                           check=True, capture_output=True)
            original = subprocess.run(["ctest", "--test-dir", str(build), "--show-only=json-v1"],
                                      capture_output=True, text=True)
            self.assertNotEqual(original.returncode, 0)
            self.assertIn("Unbuilt discovery must not run", original.stderr)
            with patch.object(ci, "ROOT", source):
                stage, report = ci.prepare(build, "unit/Rodin/Math")
                tests = {test["name"]: test for test in report["tests"]}
                self.assertEqual(set(tests), {"explicit", "dynamic"})
                direct = tests["explicit"]
                self.assertEqual(direct["command"][1:], ["--rank=3", "semi;colon"])
                props = ci.properties(direct)
                self.assertEqual(props["TIMEOUT"], 17)
                self.assertEqual(props["PROCESSORS"], 3)
                self.assertEqual(props["ENVIRONMENT"], ["RODIN_FIXTURE=ok"])
                self.assertEqual(props["RESOURCE_LOCK"], ["mesh"])
                self.assertEqual(props["WORKING_DIRECTORY"], str(build / "tests/unit/Rodin/Math"))
                self.assertIn("Unit/Rodin/Math", props["LABELS"])
                self.assertFalse((build / "tests/unit/Rodin/Math/Unbuilt").exists())
                self.assertEqual(ci.run_group(build, "unit/Rodin", "unit"), 0)
                grouped = json.loads((build / "ci-invocations.json").read_text())
                for test in grouped["tests"]:
                    self.assertIn("Unit/Rodin/Math", ci.properties(test)["LABELS"])
                self.assertEqual(len(json.loads((build / "ci-selection.json").read_text())["selected"]), 2)


class GoogleTestForwardingTest(unittest.TestCase):
    def test_semicolon_properties_and_empty_arguments_are_preserved(self):
        with tempfile.TemporaryDirectory(prefix="rodin-ci-forward-") as directory:
            source = Path(directory) / "source"
            source.mkdir()
            build = Path(directory) / "build"
            registration = (Path(__file__).resolve().parents[1] / "tests/CMakeLists.txt").read_text()
            begin = registration.index("function(rodin_discover_tests")
            end = registration.index("endfunction()", begin) + len("endfunction()")
            (source / "Registration.cmake").write_text(registration[begin:end])
            (source / "fake.cpp").write_text('#include <iostream>\nint main() { std::cout << "Suite.\\n  One\\n"; }\n')
            (source / "CMakeLists.txt").write_text('''cmake_minimum_required(VERSION 3.16)
project(Forwarding LANGUAGES CXX)
if(POLICY CMP0178)
  cmake_policy(SET CMP0178 NEW)
endif()
enable_testing()
include(GoogleTest)
include(Registration.cmake)
add_executable(Original fake.cpp)
add_executable(Wrapped fake.cpp)
gtest_discover_tests(Original TEST_PREFIX Original.
  EXTRA_ARGS "semi;colon" "" "--rank=3"
  PROPERTIES LABELS "unit;distributed" TIMEOUT 17 ENVIRONMENT "A=b;C=d")
rodin_discover_tests(Wrapped TEST_PREFIX Wrapped.
  EXTRA_ARGS "semi;colon" "" "--rank=3"
  PROPERTIES LABELS "unit;distributed" TIMEOUT 17 ENVIRONMENT "A=b;C=d")
get_property(_profiles TARGET Wrapped PROPERTY RODIN_TEST_PROFILES)
file(WRITE "${CMAKE_BINARY_DIR}/profiles.txt" "${_profiles}")
''')
            subprocess.run(["cmake", "-S", str(source), "-B", str(build)], check=True, capture_output=True)
            for target in ["Original", "Wrapped"]:
                subprocess.run(["cmake", "--build", str(build), "--target", target, "--parallel", "1"],
                               check=True, capture_output=True)
            result = json.loads(subprocess.check_output(
                ["ctest", "--test-dir", str(build), "--show-only=json-v1"], text=True))
            self.assertEqual(set((build / "profiles.txt").read_text().split(";")), {"unit", "coverage"})
            original, wrapped = result["tests"]
            self.assertEqual(original["command"][1:], wrapped["command"][1:])
            self.assertEqual(ci.properties(original), ci.properties(wrapped))


if __name__ == "__main__":
    unittest.main()

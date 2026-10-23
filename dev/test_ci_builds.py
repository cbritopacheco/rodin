"""Regression contracts for immutable shared builds and isolated test execution."""
from pathlib import Path
import copy
import hashlib
import json
import re
import sys
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import yaml
import ci_builds as ci


class SharedBuildPlanTest(unittest.TestCase):
    def setUp(self):
        self.workflow = yaml.safe_load((ci.ROOT / '.github/ci/targeted-lanes.yml').read_text())

    def test_original_compile_contracts_are_preserved(self):
        # Independently derived by executing the pre-rework configure/build
        # selectors with a recording CMake. Execution-only coordinates are absent.
        specs = {json.dumps(ci.definition(job, row), sort_keys=True)
                 for job in ci.FAMILIES
                 for row in self.workflow['jobs'][job]['strategy']['matrix']['include']}
        self.assertEqual(len(specs), 84)
        self.assertEqual(hashlib.sha256('\n'.join(sorted(specs)).encode()).hexdigest(),
                         '817bd43886398c45fbbbd6c30b005e597402761cd76ba8448a4418ed810cdee4')

    def test_every_execution_lane_has_exactly_one_compatible_producer(self):
        builds = ci.plan(self.workflow)
        for job, family in ci.FAMILIES.items():
            definition = self.workflow['jobs'][job]
            import ci_schedule
            _, executions = ci_schedule.layout()
            self.assertTrue(any(e['job'] == job and e['dependencies'] for e in executions.values()))
            self.assertFalse(any(step.get('name', '').startswith(('Build ', 'Configure '))
                                 for step in definition['steps']))
            self.assertEqual(sum(step.get('uses') == './.github/actions/test-build'
                                 for step in definition['steps']), 1)
            for row in definition['strategy']['matrix']['include']:
                expected = ci.definition(job, row)
                self.assertEqual(expected['runner'], definition['runs-on'])
                setup = next(step for step in definition['steps'] if step.get('uses') == 'aminya/setup-cpp@v1')
                self.assertEqual(expected['compiler'], setup['with']['compiler'])
                matches = [entry for entry in builds[family] if entry['id'] == ci.identity(expected)]
                self.assertEqual(matches, [{'id': ci.identity(expected), 'definition': expected}])
        self.assertEqual(sum(map(len, builds.values())), 84)
        self.assertEqual(sum(len(self.workflow['jobs'][job]['strategy']['matrix']['include'])
                             for job in ci.FAMILIES), 742)
        top = yaml.safe_load((ci.ROOT / '.github/workflows/Tests.yml').read_text())['jobs']
        self.assertTrue(all(name in top['CoverageFinalize']['needs'] for name in top if name.startswith('Configurations_')))

    def test_rank_geometry_and_study_share_the_compiled_executable(self):
        row = dict(multithreaded='ON', suite='nonlinear_boundary', physics='NonlinearPoisson',
                   context='local', geometry='light', study='all')
        expected = ci.definition('PETScConvergence', row)
        for context in ['local', 'np1', 'np2', 'np3', 'np4']:
            for geometry in ['light', 'tetrahedron', 'pyramid']:
                for study in ['all', 'P1Rates', 'P2Rates', 'P3Rates', 'Controls']:
                    self.assertEqual(ci.definition('PETScConvergence', dict(row, context=context,
                                                                          geometry=geometry, study=study)), expected)
        for key, value in [('compiler', 'gcc-15'), ('runner', 'ubuntu-26.04')]:
            self.assertNotEqual(ci.identity(dict(expected, **{key: value})), ci.identity(expected))
        changed = copy.deepcopy(expected)
        changed['cmake']['RODIN_MULTITHREADED'] = 'OFF'
        self.assertNotEqual(ci.identity(changed), ci.identity(expected))
        changed = copy.deepcopy(expected)
        changed['targets'].append('AnotherExecutable')
        self.assertNotEqual(ci.identity(changed), ci.identity(expected))

    def test_shared_curved_boundary_target_is_produced_once(self):
        ordinary = ci.definition('PETScCurvedBoundaryConvergence',
                                 dict(multithreaded='OFF', physics='NonlinearPoisson'))
        tetrahedron = ci.definition('PETScCurvedNonlinearTetrahedronConvergence',
                                   dict(multithreaded='OFF', selection='np4', partition='remaining'))
        self.assertEqual(ordinary, tetrahedron)
        entries = ci.plan(self.workflow)['petsc_boundary']
        self.assertEqual(sum(entry['id'] == ci.identity(ordinary) for entry in entries), 1)

    def test_unknown_build_coordinates_fail_before_artifact_lookup(self):
        with self.assertRaises(KeyError):
            ci.definition('UnknownFamily', {})
        for row in [dict(multithreaded='maybe'), dict(multithreaded='ON', suite='unknown'),
                    dict(multithreaded='ON', suite='bulk', physics='UnknownPhysics')]:
            with self.assertRaises(ValueError):
                ci.definition('PETScConvergence', row)


class SharedBuildArchiveTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='rodin-shared-build-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.build = self.root / 'build'
        self.build.mkdir()
        self.artifact = self.root / 'artifact'
        self.spec = ci.definition('MPIH1Restriction', dict(multithreaded='ON'))
        self.library = self.root / 'libdependency.so'
        self.library.write_bytes(b'compatible runtime')
        for function, value in [('commit', 'current-commit'), ('compiler', '14.2.0'),
                                ('runtime_libraries', {str(self.library): ci.digest(self.library)})]:
            patcher = patch.object(ci, function, return_value=value)
            patcher.start()
            self.addCleanup(patcher.stop)

    def fixture(self):
        (self.build / 'notes.gcno').write_bytes(b'zero-hit template notes')
        (self.build / 'notes.gcda').write_bytes(b'producer discovery counters')
        (self.build / 'notes.o').write_bytes(b'unneeded compile object')
        (self.build / 'empty-working-directory').mkdir()
        executable = self.build / 'test'
        executable.write_bytes(b'test executable')
        executable.chmod(0o755)
        (self.build / 'CTestTestfile.cmake').write_text('add_test(first test)\n')
        (self.build / 'RodinTestTargets-Release.tsv').write_text('target\texecutable\n')

    def package(self):
        self.fixture()
        ci.pack(self.spec, self.build, self.artifact)
        self.build.rename(self.root / 'producer-build')

    def test_two_restores_have_independent_clean_counters_and_complete_notes(self):
        self.package()
        for lane in ['first', 'second']:
            ci.restore(self.spec, self.build, self.artifact)
            self.assertFalse((self.build / 'notes.gcda').exists())
            self.assertFalse((self.build / 'notes.o').exists())
            self.assertTrue((self.build / 'empty-working-directory').is_dir())
            self.assertEqual((self.build / 'notes.gcno').read_bytes(), b'zero-hit template notes')
            self.assertEqual((self.build / 'test').stat().st_mode & 0o777, 0o755)
            self.assertTrue((self.build / 'RodinTestTargets-Release.tsv').is_file())
            (self.build / 'notes.gcda').write_text(lane)
            self.build.rename(self.root / lane)
        self.assertEqual((self.root / 'first/notes.gcda').read_text(), 'first')
        self.assertEqual((self.root / 'second/notes.gcda').read_text(), 'second')

    def test_stale_commit_workspace_compiler_and_configuration_are_rejected(self):
        self.package()
        manifest = self.artifact / 'manifest.json'
        original = json.loads(manifest.read_text())
        for key in ['commit', 'workspace', 'build', 'compiler_version', 'architecture', 'identity']:
            with self.subTest(key=key):
                manifest.write_text(json.dumps(dict(original, **{key: 'incompatible'})))
                with self.assertRaisesRegex(ValueError, key):
                    ci.restore(self.spec, self.build, self.artifact)
                self.assertFalse(self.build.exists())
        manifest.write_text(json.dumps(original))
        changed = copy.deepcopy(self.spec)
        changed['cmake']['RODIN_MULTITHREADED'] = 'OFF'
        with self.assertRaisesRegex(ValueError, 'definition'):
            ci.restore(changed, self.build, self.artifact)

    def test_corrupt_artifact_and_changed_runtime_library_are_rejected(self):
        self.package()
        self.library.write_bytes(b'incompatible runtime')
        with self.assertRaisesRegex(ValueError, 'runtime library'):
            ci.restore(self.spec, self.build, self.artifact)
        shutil.rmtree(self.build)
        archive = self.artifact / 'build.tar.gz'
        with archive.open('ab') as stream:
            stream.write(b'corruption')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            ci.restore(self.spec, self.build, self.artifact)

    def test_missing_coverage_notes_and_reused_build_location_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'coverage notes'):
            ci.pack(self.spec, self.build, self.artifact)
        self.package()
        self.build.mkdir()
        with self.assertRaisesRegex(ValueError, 'empty build'):
            ci.restore(self.spec, self.build, self.artifact)

    def test_real_cmake_build_runs_two_ctest_selections_without_recompilation(self):
        if not shutil.which('cmake') or not shutil.which('c++'):
            self.skipTest('CMake and C++ compiler required')
        source = self.root / 'source'
        source.mkdir()
        (source / 'src/Rodin').mkdir(parents=True)
        (source / 'CMakeLists.txt').write_text(
            'cmake_minimum_required(VERSION 3.16)\nproject(SharedTest LANGUAGES CXX)\n'
            'set(CMAKE_SUPPRESS_REGENERATION ON)\nenable_testing()\nadd_executable(shared src/Rodin/Operator.cpp)\n'
            'target_compile_options(shared PRIVATE --coverage)\n'
            'target_link_options(shared PRIVATE --coverage)\n'
            'add_test(NAME first COMMAND shared a)\nadd_test(NAME second COMMAND shared b)\n')
        (source / 'src/Rodin/Operator.cpp').write_text(
            'int main(int argc, char** argv) {\n'
            'if (argc != 2) return 1;\n'
            "if (argv[1][0] == 'a')\n  return 0;\n"
            "return argv[1][0] == 'b' ? 0 : 1;\n}\n")
        def run(command):
            subprocess.run(command, check=True, capture_output=True, text=True)
        run(['cmake', '-S', str(source), '-B', str(self.build)])
        run(['cmake', '--build', str(self.build), '-j', '1'])
        self.assertTrue(list(self.build.rglob('*.o')))
        notes = {str(p.relative_to(self.build)): ci.digest(p) for p in self.build.rglob('*.gcno')}
        self.assertTrue(notes)
        ci.pack(self.spec, self.build, self.artifact)
        cache = (self.build / 'CMakeCache.txt').read_text()
        compiler = re.search(r'^CMAKE_CXX_COMPILER:[^=]+=(.*)$', cache, re.M)[1]
        version = subprocess.check_output([compiler, '-dumpversion'], text=True).strip().split('.')[0]
        gcov = shutil.which('gcov-' + version)
        collect = sys.platform.startswith('linux') and shutil.which('lcov') and gcov
        reports = []
        self.build.rename(self.root / 'producer-build')
        for lane in ['first', 'second']:
            ci.restore(self.spec, self.build, self.artifact)
            self.assertFalse(list(self.build.rglob('*.o')))
            self.assertFalse(list(self.build.rglob('*.gcda')))
            run(['ctest', '--test-dir', str(self.build), '-R', '^' + lane + '$',
                 '--no-tests=error', '--output-on-failure', '-j', '1'])
            self.assertTrue(list(self.build.rglob('*.gcda')))
            if collect:
                for mode in ['initial', 'runtime']:
                    command = ['lcov', '--capture', '--directory', str(self.build),
                               '--base-directory', str(source), '--gcov-tool', gcov,
                               '--output-file', str(self.build / (mode + '.info'))]
                    if mode == 'initial':
                        command.append('--initial')
                    run(command)
                    reports.append(self.root / lane / (mode + '.info'))
            self.assertEqual({str(p.relative_to(self.build)): ci.digest(p)
                              for p in self.build.rglob('*.gcno')}, notes)
            self.assertFalse(list(self.build.rglob('*.o')))
            self.build.rename(self.root / lane)
        if collect:
            from ci_coverage import production_lines
            first = production_lines(reports[:2])
            second = production_lines(reports[2:])
            pooled = production_lines(reports)
            self.assertEqual(len(pooled), len(first))
            self.assertEqual(len(pooled), len(second))
            self.assertGreater(sum(pooled.values()), max(sum(first.values()), sum(second.values())))
            self.assertEqual(production_lines(list(reversed(reports))), pooled)


if __name__ == '__main__':
    unittest.main()

"""Uniform sharing must preserve execution contracts and avoid recompilation."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import yaml

import ci_build_graph as graph
import ci_builds as archive
import ci_test_plan as planning
import ci_tests as tests


class BuildGraphPlanTest(unittest.TestCase):
    def test_hierarchy_cuts_preserve_category_and_module_boundaries(self):
        for owner, expected in [
                ('unit/Rodin', 'unit/Rodin'),
                ('unit/Rodin/MPI/Geometry', 'unit/Rodin/MPI'),
                ('unit/Rodin/PETSc/IO', 'unit/Rodin/PETSc'),
                ('unit/Rodin/Models/Solid', 'unit/Rodin/Models'),
                ('unit/Rodin/Variational/H1', 'unit/Rodin/Variational'),
                ('manufactured/Rodin/Solver/NewtonSolver', 'manufactured/Rodin/Solver'),
                ('manufactured/Rodin/MPI/H1/Poisson', 'manufactured/Rodin/MPI'),
                ('convergence', 'convergence'),
                ('convergence/h/PETScPoisson', 'convergence/h'),
                ('convergence/p/Poisson', 'convergence/p'),
                ('convergence/hp/LinearElasticity', 'convergence/hp'),
                ('convergence/isoparametric/Diffusion', 'convergence/isoparametric')]:
            self.assertEqual(graph.build_group(owner), expected)
            self.assertEqual(graph.build_group('tests/' + owner), expected)
        with self.assertRaisesRegex(ValueError, 'Unknown test build hierarchy'):
            graph.build_group('benchmarks/Rodin/Math')

    def test_siblings_and_large_target_sets_share_one_producer(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            targets = []
            for i in range(15):
                owner = 'unit/Rodin/MPI/' + ('Geometry' if i % 2 else 'Variational')
                folder = root / 'tests' / owner
                folder.mkdir(parents=True, exist_ok=True)
                (folder / 'CMakeLists.txt').touch()
                targets.append({'target': 'Target' + str(i), 'owner': 'tests/' + owner, 'kind': 'test'})
            # Formerly isolated components remain execution partitions, not
            # additional compilation producers.
            targets.append(dict(targets[-1], target='RodinMPIMatrixRangeTest'))
            first = planning.entry('unit', 'unit/Rodin/MPI', 'Debug', 'OFF', backend='default')
            slow = planning.entry('slow', 'unit/Rodin/MPI', 'Debug', 'OFF', backend='default')
            records = [{'request': {'job': 'Ordinary', 'matrix': first}, 'selected': targets, 'configured': targets},
                       {'request': {'job': 'Ordinary', 'matrix': slow}, 'selected': targets[:1], 'configured': targets}]
            frozen = json.dumps(records, sort_keys=True)
            entry = {'id': 'common', 'definition': graph.common_specification(graph.specification('Ordinary', first))}
            with patch.object(graph, 'ROOT', root), patch.object(archive, 'commit', return_value='current'):
                document = graph.directory_plan(entry, records)
            self.assertEqual(list(document['directories']), ['unit/Rodin/MPI'])
            producer = document['directories']['unit/Rodin/MPI']
            self.assertEqual(set(producer['targets']), {t['target'] for t in targets})
            self.assertNotIn('part', producer)
            regular = document['requests'][graph.request_key('Ordinary', first)]
            slow_record = document['requests'][graph.request_key('Ordinary', slow)]
            self.assertEqual(regular['builds'], slow_record['builds'])
            self.assertEqual(regular['selected'], targets)
            self.assertEqual(slow_record['selected'], targets[:1])
            self.assertEqual(json.dumps(records, sort_keys=True), frozen)

    def test_oversized_optional_sanitizer_build_is_split_without_losing_targets(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root / 'tests/unit/Rodin/Variational'
            folder.mkdir(parents=True)
            (folder / 'CMakeLists.txt').touch()
            targets = [{'target': f'Target{i:03}', 'owner': 'tests/unit/Rodin/Variational',
                        'kind': 'test'} for i in range(133)]
            first = planning.entry('unit', 'unit/Rodin/Variational', 'Release', 'ON',
                                   asan='ON', ubsan='ON', backend='optional')
            request = {'job': 'Ordinary', 'matrix': first}
            entry = {'id': 'common', 'definition': graph.common_specification(graph.specification(**request))}
            records = [{'request': request, 'selected': targets, 'configured': targets}]
            with patch.object(graph, 'ROOT', root), patch.object(archive, 'commit', return_value='current'):
                document = graph.directory_plan(entry, records)
            parts = list(document['directories'].values())
            self.assertEqual(len(parts), 6)
            self.assertEqual([p['part'] for p in parts], list(range(1, 7)))
            self.assertTrue(all(p['common'] == 'common' and len(p['targets']) <= 24 for p in parts))
            self.assertEqual([t for part in parts for t in part['targets']],
                             sorted(t['target'] for t in targets))
            self.assertEqual(set(document['requests'][graph.request_key(**request)]['builds']),
                             {archive.identity(part) for part in parts})

    def requests(self):
        native = [{'multithreaded': mt, 'suite': 'baseline', 'group': group, 'path': tests.display_group(group)}
                  for group in planning.native_convergence_groups() for mt in ['ON', 'OFF']]
        return planning.plan(), planning.targeted_plan(native)

    def test_all_execution_requests_have_one_common_configuration(self):
        ordinary, targeted = self.requests()
        # Reviewed execution partition snapshot: numerical registrations and MPI
        # ranks are independently checked by test_ci_backend_registration.py.
        for value, digest in [(ordinary, '009a23d4c5fc380bcca406302cfb7062208a2b14f3526f2ce6078cde95005c54'),
                              (targeted, 'b7f67cd76ea3c97eb5f7ea1a6525d094b1293c3f95e0e5050ed7c5f7166dc90f')]:
            self.assertEqual(hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest(), digest)
        rows = graph.prepare(ordinary, targeted)
        self.assertEqual(len(rows), 46)
        requests = [r for row in rows for r in row['requests']]
        self.assertEqual(len(requests), 1354)
        self.assertEqual(len({graph.request_key(**r) for r in requests}), 1354)
        for row in rows:
            for request in row['requests']:
                original = graph.specification(**request)
                common = row['definition']
                for key, value in original['cmake'].items():
                    if not key.startswith('RODIN_BUILD_') and key != 'CMAKE_SUPPRESS_REGENERATION':
                        self.assertEqual(common['cmake'][key], value, request)
                self.assertEqual(common['compiler'], original['compiler'])
                self.assertEqual(common['runner'], original['runner'])

    def test_slow_and_regular_share_without_combining_configurations(self):
        ordinary, targeted = self.requests()
        rows = graph.prepare(ordinary, targeted)
        lookup = {graph.request_key(**r): row['id'] for row in rows for r in row['requests']}
        for slow in ordinary['slow']:
            regular = next(d for entries in ordinary.values() for d in entries
                           if not d['baseline'] and d['profile'] != 'slow' and
                           d['group'] == slow['group'] and d['backend'] == slow['backend'] and
                           d['build_type'] == slow['build_type'] and d['multithreaded'] == slow['multithreaded'] and
                           d['asan'] == slow['asan'] and d['ubsan'] == slow['ubsan'])
            self.assertEqual(lookup[graph.request_key('Ordinary', slow)],
                             lookup[graph.request_key('Ordinary', regular)])

    def test_runtime_jobs_restore_and_publish_without_compiling(self):
        jobs = yaml.safe_load((tests.ROOT / '.github/ci/targeted-lanes.yml').read_text())['jobs']
        ordinary, targeted = self.requests()
        for name in {d['job'] for d in targeted}:
            job = jobs[name]
            self.assertEqual(sum(s.get('uses') == './.github/actions/test-build' for s in job['steps']), 1)
            self.assertFalse(any('cmake --build' in s.get('run', '') or 'cmake -S' in s.get('run', '')
                                 for s in job['steps']))
        group = yaml.safe_load((tests.ROOT / '.github/workflows/TestGroup.yml').read_text())
        build = next(s for s in group['jobs']['Run']['steps'] if s.get('id') == 'build')
        self.assertEqual(build['uses'], './.github/actions/test-build')
        publisher = yaml.safe_load((tests.ROOT / '.github/workflows/CoveragePublish.yml').read_text())
        steps = publisher['jobs']['Publish']['steps']
        self.assertLess(next(i for i, s in enumerate(steps) if 'collect --kind' in s.get('run', '')),
                        next(i for i, s in enumerate(steps) if s.get('id') == 'upload'))
        top = yaml.safe_load((tests.ROOT / '.github/workflows/Tests.yml').read_text())['jobs']
        gate = next(s for s in top['CoverageFinalize']['steps'] if s.get('env', {}).get('CI_JOB_RESULTS'))
        self.assertEqual(gate['env']['CI_JOB_RESULTS'], '${{ toJSON(needs.*.result) }}')

    def test_every_ordinary_coordinate_and_publisher_is_wired(self):
        import ci_schedule
        jobs = ci_schedule.configuration_workflow()['jobs']
        _, executions = ci_schedule.layout()
        ordinary, _ = self.requests()
        seen = set()
        for name, job in jobs.items():
            if job.get('uses') != './.github/workflows/TestGroup.yml':
                continue
            family = executions[name]['family']
            self.assertIn(family, ordinary)
            seen.add(family)
            for coordinate in ['profile', 'backend', 'group', 'build_type', 'multithreaded',
                               'asan', 'ubsan', 'id', 'name', 'baseline']:
                self.assertEqual(job['with'][coordinate], '${{ matrix.' + coordinate + ' }}')
            publisher = next(p for p in jobs.values() if p.get('uses') == './.github/workflows/CoveragePublish.yml'
                             and p['with']['kind'] == 'ordinary' and p['with']['family'] == family)
            self.assertIn(name, publisher['needs'])
        self.assertEqual(seen, ordinary.keys())
        for name in executions:
            self.assertEqual(set(jobs[name]['needs']), {'Common', *executions[name]['dependencies']})

    def test_reviewed_targeted_commands_and_watchdogs_match_snapshot(self):
        jobs = yaml.safe_load((tests.ROOT / '.github/ci/targeted-lanes.yml').read_text())['jobs']
        snapshot = {}
        for name, job in jobs.items():
            if not any(s.get('id') == 'coverage_record' for s in job.get('steps', [])):
                continue
            snapshot[name] = {k: job.get(k) for k in ['strategy', 'env', 'timeout-minutes', 'runs-on', 'name']}
            snapshot[name]['execution'] = [s for s in job['steps'] if s.get('name', '').startswith(
                ('Run ', 'Discover ', 'Record native MPI CPU', 'Hydrate '))]
        digest = hashlib.sha256(json.dumps(snapshot, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        self.assertEqual(digest, '3fce779995376f36c30a94b0e5359e33a2612a3a8db047f86ca725ac0cb9606c')


class BuildGraphExecutionTest(unittest.TestCase):
    def test_common_objects_are_reused_by_two_directory_builds_and_isolated_runs(self):
        if not shutil.which('cmake') or not shutil.which('c++'):
            self.skipTest('CMake and C++ compiler required')
        with tempfile.TemporaryDirectory(prefix='rodin-build-graph-') as temporary:
            root = Path(temporary)
            build = root / 'build'
            geometry = root / 'tests/unit/Rodin/Geometry'
            variational = root / 'tests/unit/Rodin/Variational'
            nested = geometry / 'Nested'
            for path in [geometry, nested, variational, root / 'src/Rodin']:
                path.mkdir(parents=True)
            (root / 'CMakeLists.txt').write_text(
                'cmake_minimum_required(VERSION 3.20)\nproject(Graph LANGUAGES CXX)\n'
                'enable_testing()\nadd_compile_options(--coverage)\nadd_link_options(--coverage)\n'
                'add_library(Common STATIC src/Rodin/Common.cpp)\n'
                'add_library(Unused STATIC src/Rodin/Unused.cpp)\n'
                'add_subdirectory(tests/unit/Rodin/Geometry)\n'
                'add_subdirectory(tests/unit/Rodin/Variational)\n'
                'file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/RodinTestTargets-Debug.tsv" CONTENT '
                '"target\\towner\\texecutable\\tkind\\tprofiles\\n'
                'Common\\tsrc/Rodin\\t\\tlibrary\\t\\n'
                'First\\ttests/unit/Rodin/Geometry\\t$<TARGET_FILE:First>\\ttest\\tunit;slow\\n'
                'Third\\ttests/unit/Rodin/Geometry/Nested\\t$<TARGET_FILE:Third>\\ttest\\tunit\\n'
                'Second\\ttests/unit/Rodin/Variational\\t$<TARGET_FILE:Second>\\ttest\\tunit\\n")\n')
            (root / 'src/Rodin/Common.cpp').write_text('int common(int n) { return n == 2 ? 0 : 1; }\n')
            (root / 'src/Rodin/Unused.cpp').write_text('#error An unrelated library must not be compiled\n')
            for owner, name in [(geometry, 'First'), (variational, 'Second'), (nested, 'Third')]:
                (owner / (name + '.cpp')).write_text('int common(int); int main(int n, char**) { return common(n); }\n')
                body = f'add_executable({name} {name}.cpp)\ntarget_link_libraries({name} Common)\n'
                body += f'add_test(NAME {name} COMMAND {name} regular)\nset_tests_properties({name} PROPERTIES LABELS unit)\n'
                if name == 'First':
                    body += 'add_test(NAME Slow COMMAND First slow)\nset_tests_properties(Slow PROPERTIES LABELS slow)\n'
                    body += 'add_subdirectory(Nested)\n'
                (owner / 'CMakeLists.txt').write_text(body)
            inputs = [*root.rglob('CMakeLists.txt'), *root.rglob('*.cpp')]
            def stamp():
                for path in inputs:
                    os.utime(path, (946684800, 946684800))
            def download(name, destination):
                source = (root / 'ci-common-artifact') if name.startswith('ci-common-') else root / 'artifacts' / name.removeprefix('ci-build-')
                if destination.exists():
                    shutil.rmtree(destination)
                shutil.copytree(source, destination)
            first = planning.entry('unit', 'unit/Rodin/Geometry', 'Debug', 'OFF', backend='default')
            slow = planning.entry('slow', 'unit/Rodin/Geometry', 'Debug', 'OFF', backend='default')
            second = planning.entry('unit', 'unit/Rodin/Variational', 'Debug', 'OFF', backend='default')
            original_specification = graph.specification
            def specification(job, matrix):
                if job == 'Together':
                    return original_specification('Ordinary', first) | {'targets': ['First', 'Second', 'Third']}
                return original_specification(job, matrix)
            patches = [patch.object(graph, 'ROOT', root), patch.object(tests, 'ROOT', root),
                       patch.object(archive, 'ROOT', root), patch.object(archive, 'commit', return_value='current'),
                       patch.object(archive, 'compiler', return_value='fixture'),
                       patch.object(archive, 'runtime_libraries', return_value={}),
                       patch.object(graph, 'stamp_sources', side_effect=stamp),
                       patch.object(graph, 'download', side_effect=download),
                       patch.object(graph, 'DIRECTORY_BANDS', 2),
                       patch.object(graph, 'specification', side_effect=specification)]
            for p in patches:
                p.start()
                self.addCleanup(p.stop)
            tests.configured_profiles.cache_clear()
            common = graph.prepare({'unit': [first, second], 'slow': [slow]}, [{'job': 'Together', 'matrix': {}}])
            self.assertEqual(len(common), 1)
            (root / 'ci-common-plan.json').write_text(json.dumps(common))
            graph.compile_common(common[0], build, root / 'ci-common-records' / common[0]['id'])
            index, bands = graph.index_graph(root / 'ci-common-records')
            self.assertEqual(len(index['builds']), 2)
            self.assertEqual(next(r['definition']['targets'] for r in index['builds'].values()
                                  if r['definition']['directory'] == 'unit/Rodin/Geometry'), ['First', 'Third'])
            self.assertEqual(index['requests'][graph.request_key('Ordinary', first)]['builds'],
                             index['requests'][graph.request_key('Ordinary', slow)]['builds'])
            original_objects = {p.relative_to(build): archive.digest(p) for p in build.rglob('*.o')}
            self.assertTrue(original_objects)
            build.rename(root / 'common-producer-build')
            (root / 'artifacts').mkdir()
            for key, row in index['builds'].items():
                graph.compile_directory(row['definition'], build)
                self.assertTrue(all(archive.digest(build / p) == digest for p, digest in original_objects.items()))
                (root / 'ci-build-artifact').rename(root / 'artifacts' / key)
                (root / 'ci-built-record').rename(root / ('receipt-' + key))
                build.rename(root / ('producer-' + key))
            graph.verify_builds(index, root, 'current')
            # Combining two directory archives must preserve both executables.
            graph.restore_request(index, 'Together', {}, build)
            self.assertFalse(list(build.rglob('*.o')))
            self.assertFalse(list(build.rglob('*.gcda')))
            subprocess.run(['ctest', '--test-dir', str(build), '--no-tests=error', '-j', '1', '--output-on-failure'], check=True)
            self.assertTrue(list(build.rglob('*.gcda')))
            build.rename(root / 'first-execution')
            graph.restore_request(index, 'Ordinary', first, build)
            self.assertFalse(list(build.rglob('*.gcda')))
            self.assertEqual(tests.run_group(build, first['group'], 'unit', 'default'), 0)
            self.assertFalse((build / 'tests/unit/Rodin/Variational/Second').exists())
            self.assertTrue((build / 'tests/unit/Rodin/Geometry/Nested/Third').is_file())
            build.rename(root / 'regular-execution')
            graph.restore_request(index, 'Ordinary', slow, build)
            self.assertFalse(list(build.rglob('*.gcda')))
            self.assertEqual(tests.run_group(build, slow['group'], 'slow', 'default'), 0)


class BuildArtifactDownloadTest(unittest.TestCase):
    def test_requires_preloaded_artifacts_without_rest_fallback(self):
        with patch.dict(os.environ, CI_RESTORED_ARTIFACTS=''), patch.object(graph.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'not preloaded'):
                graph.download('ci-build-test', Path('unused'))
            run.assert_not_called()

    def test_missing_or_partial_archive_remains_failure_even_with_read_token(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'ci-build-test'
            source.mkdir()
            env = {'CI_RESTORED_ARTIFACTS': directory, 'GH_TOKEN': 'test-token'}
            with patch.dict(os.environ, env), patch.object(graph.subprocess, 'run') as run:
                for present in [None, 'manifest.json', 'build.tar.gz']:
                    for path in source.iterdir():
                        path.unlink()
                    if present:
                        (source / present).write_bytes(b'incomplete')
                    with self.assertRaisesRegex(ValueError, 'not downloaded'):
                        graph.download('ci-build-test', Path(directory) / 'unused')
                run.assert_not_called()

    def test_complete_preloaded_archive_is_returned_without_network(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'ci-build-test'
            source.mkdir()
            (source / 'manifest.json').write_text('{}')
            (source / 'build.tar.gz').write_bytes(b'archive')
            with patch.dict(os.environ, CI_RESTORED_ARTIFACTS=directory), patch.object(graph.subprocess, 'run') as run:
                self.assertEqual(graph.download('ci-build-test', Path(directory) / 'unused'), source)
                run.assert_not_called()


if __name__ == '__main__':
    unittest.main()

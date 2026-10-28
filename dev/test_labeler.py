"""Exercise the pinned labeler against a local GitHub API, never a real PR."""
from __future__ import annotations

import base64
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import threading
import unittest
from urllib.parse import parse_qs, unquote, urlsplit

import yaml

ROOT = Path(__file__).resolve().parents[1]
MODULES = {
    'Adaptation': 'adaptation', 'Advection': 'advection', 'Alert': 'alerts',
    'Assembly': 'assembly', 'Context': 'context', 'Distance': 'distance',
    'Eikonal': 'eikonal', 'FormLanguage': 'form-language', 'Geometry': 'geometry',
    'Heart': 'heart', 'Hilbert': 'hilbert', 'IO': 'io', 'Location': 'location',
    'Math': 'math', 'MFEM': 'mfem', 'MMG': 'mmg', 'MPI': 'mpi', 'PETSc': 'petsc',
    'QF': 'qf', 'Scotch': 'scotch', 'Serialization': 'serialization', 'Solid': 'solid',
    'Solver': 'solver', 'Test': 'test', 'Utility': 'utility', 'Variational': 'variational',
}


class LabelCatalogTest(unittest.TestCase):
    def test_every_labeler_rule_has_documented_catalog_metadata(self):
        rules = yaml.safe_load((ROOT / '.github/labeler.yml').read_text())
        catalog = yaml.safe_load((ROOT / '.github/labels.yml').read_text())
        self.assertLessEqual(set(rules), set(catalog))
        for name, metadata in catalog.items():
            with self.subTest(label=name):
                self.assertEqual(set(metadata), {'color', 'description'})
                self.assertRegex(metadata['color'], r'^[0-9A-Fa-f]{6}$')
                self.assertTrue(metadata['description'].strip())
                self.assertLessEqual(len(metadata['description']), 100)
                self.assertLessEqual(len(name), 50)
        for name in ['convergence', 'installation', 'nafems']:
            self.assertIn(name, rules)
            self.assertNotIn(name + '-tests', rules)
            self.assertNotIn(name + '-tests', catalog)

    def test_related_labels_share_their_palette(self):
        catalog = yaml.safe_load((ROOT / '.github/labels.yml').read_text())
        for group in [
                ['core', 'geometry', 'location', 'math', 'variational', 'assembly', 'solver'],
                ['adaptation', 'advection', 'distance', 'eikonal', 'heart', 'hilbert', 'solid'],
                ['mfem', 'mmg', 'mpi', 'petsc', 'scotch'],
                ['unit-tests', 'manufactured-tests', 'convergence', 'installation', 'nafems', 'benchmarks'],
                ['ci', 'cmake', 'developer-tooling', 'dependencies', 'python', 'resources'],
                ['bug', 'hotfix']]:
            with self.subTest(group=group):
                self.assertEqual(len({catalog[name]['color'] for name in group}), 1)


class LabelerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        action = os.environ.get('LABELER_ACTION_PATH')
        if not action:
            raise RuntimeError('Set LABELER_ACTION_PATH to the pinned actions/labeler checkout; see dev/labeler.md')
        cls.action = Path(action).resolve()
        cls.node = os.environ.get('NODE_BINARY') or shutil.which('node')
        if not cls.node or not (cls.action / 'dist/index.js').is_file():
            raise RuntimeError('Node.js and the pinned labeler dist/index.js are required')
        version = subprocess.check_output([cls.node, '--version'], text=True)
        if int(version.lstrip('v').split('.')[0]) < 24:
            raise RuntimeError('The pinned labeler requires Node.js 24 or later')

    def label_files(self, files, existing=(), late=()):
        """Run the real action with paginated files and inspect its API writes."""
        config = (ROOT / '.github/labeler.yml').read_bytes()
        calls, writes, pages, unexpected = [], [], [], []
        get_count = 0
        prefix = '/repos/owner/rodin'

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *_):
                pass

            def reply(self, payload, status=200, headers=None):
                body = json.dumps(payload).encode()
                self.send_response(status)
                self.send_header('Content-Type', 'application/json')
                self.send_header('Content-Length', str(len(body)))
                for key, value in (headers or {}).items():
                    self.send_header(key, value)
                self.end_headers()
                self.wfile.write(body)

            def do_GET(self):
                nonlocal get_count
                request = urlsplit(self.path)
                path = unquote(request.path)
                query = parse_qs(request.query)
                calls.append(('GET', path, query))
                if path == prefix + '/pulls/368':
                    labels = list(existing) + (list(late) if get_count else [])
                    get_count += 1
                    self.reply({'number': 368, 'labels': [{'name': label} for label in labels],
                                'head': {'ref': 'topic'}, 'base': {'ref': 'develop'}})
                elif path == prefix + '/pulls/368/files':
                    page = int(query.get('page', ['1'])[0])
                    size = int(query.get('per_page', ['30'])[0])
                    pages.append(page)
                    headers = {}
                    if page * size < len(files):
                        headers['Link'] = (f'<http://127.0.0.1:{self.server.server_port}' +
                                           prefix + f'/pulls/368/files?page={page + 1}&per_page={size}>; rel="next"')
                    self.reply([{'filename': f} for f in files[(page - 1) * size:page * size]], headers=headers)
                elif path == prefix + '/contents/.github/labeler.yml':
                    self.reply({'content': base64.b64encode(config).decode(), 'encoding': 'base64'})
                else:
                    unexpected.append(('GET', self.path))
                    self.reply({'message': 'Unexpected mock API request'}, status=400)

            def do_PUT(self):
                calls.append(('PUT', self.path, {}))
                payload = json.loads(self.rfile.read(int(self.headers['Content-Length'])))
                if self.path != prefix + '/issues/368/labels':
                    unexpected.append(('PUT', self.path))
                    self.reply({'message': 'Unexpected mock API write'}, status=400)
                    return
                writes.append(payload['labels'])
                self.reply([{'name': label} for label in payload['labels']])

        server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        thread = threading.Thread(target=lambda: server.serve_forever(poll_interval=0.01), daemon=True)
        thread.start()
        try:
            with tempfile.TemporaryDirectory() as folder:
                work = Path(folder)
                (work / 'event.json').write_text(json.dumps({'pull_request': {'number': 368}}))
                (work / 'outputs').touch()
                env = {key: value for key, value in os.environ.items() if not key.startswith(
                    ('GITHUB_', 'INPUT_', 'ACTIONS_')) and key.upper() not in
                    {'HTTP_PROXY', 'HTTPS_PROXY', 'ALL_PROXY', 'NODE_ENV'}}
                env.update({
                    'GITHUB_REPOSITORY': 'owner/rodin', 'GITHUB_EVENT_NAME': 'pull_request_target',
                    'GITHUB_EVENT_PATH': str(work / 'event.json'), 'GITHUB_SHA': 'b' * 40,
                    'GITHUB_WORKSPACE': folder, 'GITHUB_OUTPUT': str(work / 'outputs'),
                    'GITHUB_API_URL': f'http://127.0.0.1:{server.server_port}',
                    'INPUT_REPO-TOKEN': 'local-test-token', 'INPUT_CONFIGURATION-PATH': '.github/labeler.yml',
                    'INPUT_SYNC-LABELS': 'false', 'INPUT_DOT': 'true', 'NO_PROXY': '127.0.0.1,localhost',
                })
                result = subprocess.run([self.node, str(self.action / 'dist/index.js')], cwd=folder,
                                        env=env, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stdout[-6000:] + result.stderr[-2000:] + str(unexpected))
        finally:
            server.shutdown()
            server.server_close()
            thread.join(timeout=5)
        self.assertFalse(unexpected, unexpected)
        self.assertTrue(pages, 'The action did not enumerate any changed files')
        fetched = [query for method, path, query in calls if '/contents/' in path]
        self.assertEqual(fetched, [{'ref': ['b' * 40]}], 'Rules must be fetched from the trusted base SHA')
        self.assertLessEqual(len(writes), 1)
        self.pages = pages
        return set(writes[-1] if writes else existing)

    def test_source_directories_and_umbrella_headers_match_independently(self):
        self.assertEqual(set(MODULES), {p.name for p in (ROOT / 'src/Rodin').iterdir() if p.is_dir()})
        for module, label in MODULES.items():
            for path in ['src/Rodin/' + module + '/Only.h', 'src/Rodin/' + module + '.h']:
                with self.subTest(path=path):
                    self.assertEqual(self.label_files([path]), {label})

    def test_test_only_changes_receive_module_and_category_labels(self):
        for module, label in MODULES.items():
            with self.subTest(module=module):
                self.assertEqual(self.label_files(['tests/unit/Rodin/' + module + '/OnlyTest.cpp']),
                                 {'unit-tests', label})

    def test_nested_backend_tests_receive_both_module_labels(self):
        cases = {
            'tests/unit/Rodin/MPI/Geometry/OnlyTest.cpp': {'unit-tests', 'mpi', 'geometry'},
            'tests/unit/Rodin/PETSc/IO/OnlyTest.cpp': {'unit-tests', 'petsc', 'io'},
            'tests/unit/Rodin/PETSc/MPIGridFunctionTest.cpp': {'unit-tests', 'petsc', 'mpi'},
            'tests/manufactured/Rodin/MPI/H1/PoissonTest.cpp': {'manufactured-tests', 'mpi', 'variational'},
            'tests/unit/Rodin/Models/Solid/OnlyTest.cpp': {'unit-tests', 'solid'},
            'tests/manufactured/Rodin/Models/Advection/OnlyTest.cpp': {'manufactured-tests', 'advection'},
            'tests/manufactured/Rodin/P1/LinearElasticity3D.cpp':
                {'manufactured-tests', 'solid', 'linear-elasticity', 'variational'},
            'tests/benchmarks/Rodin/Math/OnlyBenchmark.cpp': {'benchmarks', 'math'},
        }
        for path, expected in cases.items():
            with self.subTest(path=path):
                self.assertEqual(self.label_files([path]), expected)

    def test_convergence_physics_and_backend_paths(self):
        cases = {
            'tests/convergence/h/PETScMPIPoisson/Only.cpp': {'convergence', 'petsc', 'mpi', 'variational'},
            'tests/convergence/hp/PETScLinearElasticity/Only.cpp':
                {'convergence', 'petsc', 'solid', 'linear-elasticity', 'variational'},
            'tests/convergence/isoparametric/GeometryApproximation/Only.cpp':
                {'convergence', 'geometry', 'variational'},
        }
        for path, expected in cases.items():
            with self.subTest(path=path):
                self.assertEqual(self.label_files([path]), expected)

    def test_core_headers_and_threads_tests(self):
        for path in sorted((ROOT / 'src/Rodin').glob('*')):
            if path.is_file() and path.suffix in {'.h', '.cpp', '.in'} and path.stem not in MODULES:
                with self.subTest(path=path):
                    self.assertEqual(self.label_files([path.relative_to(ROOT).as_posix()]), {'core'})
        self.assertEqual(self.label_files(['tests/unit/Rodin/Threads/OnlyTest.cpp']), {'unit-tests', 'core'})

    def test_documentation_matches_root_nested_and_dot_directories(self):
        for path, expected in {
            'README.md': {'docs'}, '.github/workflows/README.md': {'docs', 'ci'},
            'tests/convergence/README.md': {'docs', 'convergence', 'variational'},
            'examples/Geometry/README.md': {'docs', 'examples'},
        }.items():
            with self.subTest(path=path):
                self.assertEqual(self.label_files([path]), expected)

    def test_infrastructure_categories(self):
        cases = {
            '.github/workflows/Tests.yml': {'ci'}, 'codecov.yml': {'ci'},
            '.github/labels.yml': {'ci'},
            'dev/check_something.py': {'developer-tooling'},
            'src/Rodin/Math/CMakeLists.txt': {'cmake', 'math'},
            'tests/unit/Rodin/Geometry/CMakeLists.txt': {'cmake', 'geometry', 'unit-tests'},
            'tests/cmake/ThirdPartyIsolation/Only.cmake': {'cmake'},
            'tests/installation/Only.cpp': {'installation'},
            'tests/nafems/cases/Only.cpp': {'nafems'},
            'py/Rodin/Only.cpp': {'python'}, 'resources/Only.mesh': {'resources'},
            'third-party/mmg/Only.h': {'dependencies', 'mmg'},
        }
        for path, expected in cases.items():
            with self.subTest(path=path):
                self.assertEqual(self.label_files([path]), expected)

    def test_unrelated_paths_and_similar_module_names_do_not_match(self):
        for path, expected in {
            'unrelated/file.txt': set(), 'unrelated/Math/file.cpp': set(),
            'tests/unit/Rodin/Mathematics/Only.cpp': {'unit-tests'},
        }.items():
            with self.subTest(path=path):
                self.assertEqual(self.label_files([path]), expected)

    def test_large_changes_include_later_file_pages(self):
        files = ['unrelated/file' + str(i) + '.txt' for i in range(101)]
        files.append('src/Rodin/Math/Only.h')
        self.assertEqual(self.label_files(files), {'math'})
        self.assertGreater(len(self.pages), 1)

    def test_manual_and_concurrently_added_labels_are_preserved(self):
        self.assertEqual(self.label_files(['src/Rodin/Math/Only.h'], existing=['bug', 'geometry'],
                                         late=['help wanted']), {'bug', 'geometry', 'help wanted', 'math'})

    def test_workflows_use_same_pinned_action_and_separate_write_permissions(self):
        workflow = yaml.safe_load((ROOT / '.github/workflows/Labeler.yml').read_text())
        job = workflow['jobs']['label']
        self.assertEqual(set(workflow['on']), {'pull_request_target'})
        self.assertFalse(any('actions/checkout@' in s.get('uses', '') for s in job['steps']))
        action = next(s for s in job['steps'] if s.get('uses', '').startswith('actions/labeler@'))
        revision = action['uses'].split('@')[1]
        self.assertRegex(revision, r'^[a-f0-9]{40}$')
        self.assertIs(action['with']['dot'], True)
        self.assertIs(action['with']['sync-labels'], False)
        self.assertEqual(job['permissions'], {'contents': 'read', 'pull-requests': 'write', 'issues': 'write'})
        checks = yaml.safe_load((ROOT / '.github/workflows/LabelerRules.yml').read_text())
        self.assertEqual(set(checks['on']), {'pull_request', 'push'})
        self.assertNotIn('paths', checks['on']['pull_request'] or {})
        self.assertEqual(checks['permissions'], {'contents': 'read'})
        checkout = next(s for s in checks['jobs']['Rules']['steps'] if
                        s.get('with', {}).get('repository') == 'actions/labeler')
        self.assertEqual(checkout['with']['ref'], revision)
        self.assertFalse(checkout['with']['persist-credentials'])
        self.assertEqual(json.loads((self.action / 'package.json').read_text())['version'], '6.2.0')


if __name__ == '__main__':
    unittest.main()

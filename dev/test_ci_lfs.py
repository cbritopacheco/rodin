"""Bound LFS retries without bypassing hydration or permanent failures."""
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import ci_lfs


class LFSRetryTest(unittest.TestCase):
    def response(self, output='', code=0):
        return subprocess.CompletedProcess([], code, stdout='', stderr=output)

    def test_transient_recovery_preserves_include_and_bounds(self):
        include = 'resources/a.mesh,resources/tests/**'
        with patch.object(ci_lfs.subprocess, 'run', side_effect=[
                self.response('batch response: Fatal error: Server error from HTTP 502', 2),
                self.response('connection reset by peer', 2), self.response('HTTP 502', 2),
                self.response('HTTP 503', 2), self.response()]) as run, \
                patch.object(ci_lfs.time, 'sleep') as sleep:
            ci_lfs.pull(include)
        self.assertEqual(run.call_count, 5)
        for call in run.call_args_list:
            self.assertEqual(call.args[0], ['git', 'lfs', 'pull', '--include=' + include])
            self.assertEqual(call.kwargs['timeout'], 300)
        self.assertEqual([c.args[0] for c in sleep.call_args_list], [15, 30, 60, 120])

    def test_exhaustion_fails_closed(self):
        for failure in [self.response('HTTP 503', 2), subprocess.TimeoutExpired('git', 300)]:
            with self.subTest(failure=failure), \
                    patch.object(ci_lfs.subprocess, 'run', side_effect=[failure] * 5) as run, \
                    patch.object(ci_lfs.time, 'sleep') as sleep:
                with self.assertRaises((RuntimeError, subprocess.CalledProcessError)):
                    ci_lfs.pull('resources/**')
                self.assertEqual(run.call_count, 5)
                self.assertEqual(sleep.call_count, 4)

    def test_permanent_and_unknown_failures_are_not_retried(self):
        for output in ['HTTP 403', 'HTTP 404', 'bandwidth limit exceeded', 'invalid pointer',
                       'unknown failure', 'HTTP 502; authentication failed']:
            with self.subTest(output=output), \
                    patch.object(ci_lfs.subprocess, 'run', return_value=self.response(output, 2)) as run, \
                    patch.object(ci_lfs.time, 'sleep') as sleep:
                with self.assertRaises(subprocess.CalledProcessError):
                    ci_lfs.pull('resources/**')
                run.assert_called_once()
                sleep.assert_not_called()

    def test_cli_recovers_using_the_same_resource_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / 'git'
            state = root / 'calls.json'
            executable.write_text('''#!/usr/bin/env python3
import json, os, sys
from pathlib import Path
state = Path(os.environ['LFS_CALLS'])
calls = json.loads(state.read_text()) if state.exists() else []
calls.append(sys.argv[1:])
state.write_text(json.dumps(calls))
if len(calls) <= 3:
    print('batch response: Fatal error: Server error from HTTP 502', file=sys.stderr)
    sys.exit(2)
''')
            executable.chmod(0o755)
            include = 'resources/one.mesh,resources/tests/**'
            env = dict(os.environ, PATH=str(root) + os.pathsep + os.environ['PATH'],
                       LFS_CALLS=str(state))
            # Exercise the CLI and subprocesses without real backoff sleeps.
            launch = ('import runpy, sys, time; time.sleep=lambda seconds: None; '
                      'sys.argv=sys.argv[1:]; runpy.run_path(sys.argv[0], run_name="__main__")')
            result = subprocess.run(['python3', '-c', launch, str(Path(ci_lfs.__file__)), '--include=' + include],
                                    env=env, text=True, capture_output=True, timeout=15)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(json.loads(state.read_text()),
                             [['lfs', 'pull', '--include=' + include]] * 4)
            self.assertIn('retry 4/5', result.stdout)

    def test_uniform_grid_elasticity_partition_has_no_lfs_dependency(self):
        import yaml
        root = Path(__file__).resolve().parents[1]
        workflow = yaml.safe_load((root / '.github/workflows/TestGroup.yml').read_text())
        hydrate = next(s for s in workflow['jobs']['Run']['steps'] if s.get('name') == 'Hydrate test resources')
        self.assertEqual(hydrate['if'], "inputs.group != 'manufactured/Rodin/P1/LinearElasticity3D'")
        source = (root / 'tests/manufactured/Rodin/P1/LinearElasticity3D.cpp').read_text()
        self.assertIn('UniformGrid(', source)
        self.assertNotIn('resources/', source)

    def test_all_explicit_workflow_pulls_use_bounded_helper(self):
        root = Path(__file__).resolve().parents[1]
        for name in ['TestGroup', 'TestDirectoryBuild', 'BenchmarkGroup']:
            source = (root / '.github/workflows' / (name + '.yml')).read_text()
            self.assertNotIn('run: git lfs pull', source)
            self.assertIn('run: python3 dev/ci_lfs.py --include=', source)


if __name__ == '__main__':
    unittest.main()

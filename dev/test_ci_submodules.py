"""Retry transient submodule transfers without changing pinned gitlink commits."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import tempfile
import unittest
from unittest.mock import patch

import yaml

import ci_submodules

ROOT = Path(__file__).resolve().parents[1]


class SubmoduleRetryTest(unittest.TestCase):
    def response(self, output='', code=0):
        return subprocess.CompletedProcess([], code, stdout='', stderr=output)

    def test_transient_errors_and_timeouts_recover_with_bounded_attempts(self):
        for error in [self.response('fatal: server returned error: 503', 1),
                      self.response('Could not resolve host: github.com', 1),
                      subprocess.TimeoutExpired('git', 300)]:
            with self.subTest(error=error), \
                    patch.object(ci_submodules, 'run', side_effect=[
                        self.response(), error, self.response()]) as run, \
                    patch.object(ci_submodules.time, 'sleep') as sleep:
                ci_submodules.hydrate()
            self.assertEqual(run.call_args_list[0].args[0], ['git', 'submodule', 'sync', '--recursive'])
            self.assertEqual(run.call_args_list[0].kwargs['timeout'], 60)
            for call in run.call_args_list[1:]:
                self.assertEqual(call.args[0], ['git', 'submodule', 'update', '--init', '--recursive'])
                self.assertEqual(call.kwargs['timeout'], 300)
            sleep.assert_called_once_with(15)

    def test_exhaustion_and_permanent_failures_fail_closed(self):
        for error, attempts in [('HTTP 503', 5), ('HTTP 403', 1),
                                ('HTTP 503; authentication failed', 1),
                                ('not our ref 012345', 1), ('unknown git failure', 1)]:
            with self.subTest(error=error), \
                    patch.object(ci_submodules, 'run', side_effect=[
                        self.response()] + [self.response(error, 1)] * attempts) as run, \
                    patch.object(ci_submodules.time, 'sleep') as sleep:
                with self.assertRaises(subprocess.CalledProcessError):
                    ci_submodules.hydrate()
            self.assertEqual(run.call_count, attempts + 1)
            self.assertEqual([call.args[0] for call in sleep.call_args_list],
                             list(ci_submodules.RETRY_DELAYS[:attempts - 1]))

    def test_watchdog_stops_clone_children_before_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory) / 'abandoned-child'
            child = ('import time; from pathlib import Path; time.sleep(0.6); '
                     f'Path({str(marker)!r}).write_text("still running")')
            parent = ('import subprocess, sys, time; '
                      f'subprocess.Popen([sys.executable, "-c", {child!r}]); '
                      'print("started", flush=True); time.sleep(10)')
            with self.assertRaises(subprocess.TimeoutExpired) as failure:
                ci_submodules.run([sys.executable, '-c', parent], timeout=0.3)
            output = failure.exception.output
            self.assertIn('started', output.decode() if isinstance(output, bytes) else output)
            time.sleep(0.8)
            self.assertFalse(marker.exists(), 'Timed-out clone child survived into the next attempt')

    def test_cli_recovers_and_checks_out_recursive_pins_instead_of_branch_tips(self):
        git = shutil.which('git')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            env = dict(os.environ, GIT_ALLOW_PROTOCOL='file', GIT_CONFIG_NOSYSTEM='1',
                       GIT_CONFIG_GLOBAL=os.devnull, GIT_AUTHOR_NAME='CI test',
                       GIT_AUTHOR_EMAIL='ci@example.invalid', GIT_COMMITTER_NAME='CI test',
                       GIT_COMMITTER_EMAIL='ci@example.invalid')
            def command(cwd, *args):
                return subprocess.check_output([git, *args], cwd=cwd, env=env,
                                               text=True, stderr=subprocess.STDOUT).strip()
            def repository(name):
                path = root / name
                path.mkdir()
                command(path, 'init', '-q')
                (path / 'data').write_text('first\n')
                command(path, 'add', '.')
                command(path, 'commit', '-qm', 'first')
                return path
            leaf = repository('leaf')
            leaf_pin = command(leaf, 'rev-parse', 'HEAD')
            middle = repository('middle')
            command(middle, 'submodule', 'add', str(leaf), 'nested')
            command(middle, 'commit', '-qam', 'pin nested')
            middle_pin = command(middle, 'rev-parse', 'HEAD')
            parent = repository('parent')
            command(parent, 'submodule', 'add', str(middle), 'dependency')
            command(parent, 'commit', '-qam', 'pin dependency')
            for path in [leaf, middle]:
                (path / 'data').write_text('new branch tip\n')
                command(path, 'commit', '-qam', 'advance tip')
            checkout = root / 'checkout'
            command(root, 'clone', str(parent), str(checkout))
            tools = root / 'bin'
            tools.mkdir()
            state = root / 'calls.json'
            shim = tools / 'git'
            shim.write_text('''#!/usr/bin/env python3
import json, os, sys
from pathlib import Path
state = Path(os.environ['SUBMODULE_CALLS'])
calls = json.loads(state.read_text()) if state.exists() else []
calls.append(sys.argv[1:])
state.write_text(json.dumps(calls))
if len(calls) == 2:
    print('fatal: requested URL returned error: 503', file=sys.stderr)
    sys.exit(1)
os.execv(os.environ['REAL_GIT'], [os.environ['REAL_GIT'], *sys.argv[1:]])
''')
            shim.chmod(0o755)
            env.update(PATH=str(tools) + os.pathsep + os.environ['PATH'],
                       SUBMODULE_CALLS=str(state), REAL_GIT=git)
            launch = ('import runpy, sys, time; time.sleep=lambda seconds: None; '
                      'runpy.run_path(sys.argv[1], run_name="__main__")')
            result = subprocess.run(['python3', '-c', launch, str(Path(ci_submodules.__file__))],
                                    cwd=checkout, env=env, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(command(checkout / 'dependency', 'rev-parse', 'HEAD'), middle_pin)
            self.assertEqual(command(checkout / 'dependency/nested', 'rev-parse', 'HEAD'), leaf_pin)
            self.assertEqual(json.loads(state.read_text()), [
                ['submodule', 'sync', '--recursive'],
                ['submodule', 'update', '--init', '--recursive'],
                ['submodule', 'update', '--init', '--recursive']])

    def test_test_and_benchmark_checkouts_hydrate_before_using_local_actions(self):
        paths = list((ROOT / '.github/workflows').glob('Test*.yml'))
        paths += list((ROOT / '.github/workflows').glob('Benchmark*.yml'))
        paths += [ROOT / '.github/ci/targeted-lanes.yml']
        count = 0
        for path in paths:
            document = yaml.safe_load(path.read_text())
            for job in document.get('jobs', {}).values():
                steps = job.get('steps', [])
                for index, step in enumerate(steps):
                    if not step.get('uses', '').startswith('actions/checkout@'):
                        continue
                    if 'submodules' not in step.get('with', {}):
                        continue  # Planning/publication needs only the root checkout.
                    count += 1
                    self.assertIs(step.get('with', {}).get('submodules'), False, str(path))
                    following = index + 1
                    if steps[following].get('name') == 'Start resource monitor':
                        self.assertTrue(steps[following]['continue-on-error'])
                        following += 1
                    self.assertEqual(steps[following].get('run'), 'python3 dev/ci_submodules.py', str(path))
        self.assertGreater(count, 30)


if __name__ == '__main__':
    unittest.main()

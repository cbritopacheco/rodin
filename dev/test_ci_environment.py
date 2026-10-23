"""Provisioning failures must not become stale packages or skipped numerical work."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import yaml

import ci_environment as environment


ROOT = Path(__file__).resolve().parents[1]


class EnvironmentTest(unittest.TestCase):
    def setUp(self):
        read_text = Path.read_text
        self.release = patch.object(Path, 'read_text', autospec=True,
                                   side_effect=lambda path, *args, **kwargs:
                                   'ID=ubuntu\nVERSION_ID="24.04"\n' if str(path) == '/etc/os-release'
                                   else read_text(path, *args, **kwargs))
        self.release.start()
        self.addCleanup(self.release.stop)

    def test_package_input_is_canonical_and_rejects_shell_or_apt_options(self):
        self.assertEqual(environment.packages('libx-dev lcov libx-dev'), ['lcov', 'libx-dev'])
        self.assertEqual(environment.packages('libx:amd64=1:2.3-4~ubuntu1'),
                         ['libx:amd64=1:2.3-4~ubuntu1'])
        for value in ['', '-y lcov', 'lcov; true', '$(id)', 'foo==1', 'foo >= 1']:
            with self.subTest(value=value), self.assertRaises(ValueError):
                environment.packages(value)

    def test_ample_space_never_removes_tools_or_requires_hosted_runner(self):
        with patch.object(environment.shutil, 'disk_usage', return_value=SimpleNamespace(free=20 * environment.GIB)), \
                patch.object(environment.subprocess, 'run') as run, patch.dict(os.environ, {}, clear=True):
            environment.disk_space(Path.cwd())
            run.assert_not_called()

    def test_cleanup_stops_as_soon_as_headroom_is_sufficient(self):
        with patch.object(environment.shutil, 'disk_usage', side_effect=[
                SimpleNamespace(free=3 * environment.GIB), SimpleNamespace(free=9 * environment.GIB)]), \
                patch.object(environment.Path, 'exists', return_value=True), \
                patch.object(environment.subprocess, 'run') as run, \
                patch.dict(os.environ, GITHUB_ACTIONS='true', RUNNER_ENVIRONMENT='github-hosted'):
            environment.disk_space(Path.cwd())
            run.assert_called_once_with(['sudo', 'rm', '-rf', '--', environment.CLEANUP_PATHS[0]], check=True)

    def test_low_space_on_local_or_self_hosted_runner_cannot_remove_system_tools(self):
        for runner in ['', 'self-hosted']:
            with self.subTest(runner=runner), \
                    patch.object(environment.shutil, 'disk_usage', return_value=SimpleNamespace(free=0)), \
                    patch.object(environment.subprocess, 'run') as run, \
                    patch.dict(os.environ, GITHUB_ACTIONS='true', RUNNER_ENVIRONMENT=runner):
                with self.assertRaisesRegex(RuntimeError, 'restricted'):
                    environment.disk_space(Path.cwd())
                run.assert_not_called()

    def test_cleanup_cannot_silently_accept_insufficient_space(self):
        with patch.object(environment.shutil, 'disk_usage', return_value=SimpleNamespace(free=0)), \
                patch.object(environment.Path, 'exists', return_value=True), \
                patch.object(environment.subprocess, 'run') as run, \
                patch.dict(os.environ, GITHUB_ACTIONS='true', RUNNER_ENVIRONMENT='github-hosted'):
            with self.assertRaisesRegex(RuntimeError, 'after allowlisted'):
                environment.disk_space(Path.cwd())
            self.assertEqual(run.call_count, len(environment.CLEANUP_PATHS))
            self.assertFalse(any('clean' in c.args[0] for c in run.call_args_list))

    def test_mirror_copy_preserves_signatures_suites_components_and_third_party_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); system = root / 'etc'; parts = system / 'sources.list.d'
            parts.mkdir(parents=True)
            deb822 = ('Types: deb\nURIs: http://azure.archive.ubuntu.com/ubuntu/\n'
                      'Suites: noble noble-updates\nComponents: main universe\n'
                      'Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg\n')
            (parts / 'ubuntu.sources').write_text(deb822)
            third_party = 'deb [signed-by=/keys/vendor.gpg] https://vendor.example/apt stable main\n'
            (parts / 'vendor.list').write_text(third_party)
            (system / 'sources.list').write_text('deb http://security.ubuntu.com/ubuntu noble-security main\n')
            original = {str(p): p.read_bytes() for p in system.rglob('*') if p.is_file()}
            destination = root / 'fallback'
            options = environment.sources(destination, environment.MIRRORS[1], system)
            self.assertIn('Dir::Etc::sourceparts=' + str(destination / 'sources.list.d'), options)
            self.assertEqual((destination / 'sources.list.d/ubuntu.sources').read_text(),
                             deb822.replace('http://azure.archive.ubuntu.com/ubuntu/',
                                            'https://archive.ubuntu.com/ubuntu/'))
            self.assertEqual((destination / 'sources.list.d/vendor.list').read_text(), third_party)
            self.assertIn('noble-security main', (destination / 'sources.list').read_text())
            self.assertEqual(original, {str(p): p.read_bytes() for p in system.rglob('*') if p.is_file()})

    def test_apt_bounds_entire_privileged_process_group_without_weakening_validation(self):
        with patch.object(environment.subprocess, 'run', return_value=SimpleNamespace(stdout='ok')) as run:
            self.assertEqual(environment.apt(['update'], ['-o', 'source=copy'], Path('/cache'), 90, capture=True), 'ok')
        command = run.call_args.args[0]
        self.assertEqual(command[:10], ['sudo', 'env', 'DEBIAN_FRONTEND=noninteractive', 'LC_ALL=C',
                                       'timeout', '--signal=TERM', '--kill-after=5', '90', 'apt-get', '-o'])
        for option in ['Acquire::Retries=1', 'Acquire::https::Timeout=20',
                       'APT::Update::Error-Mode=any', 'APT::Keep-Downloaded-Packages=true',
                       'Dir::Cache::archives=/cache']:
            self.assertIn(option, command)
        self.assertTrue(run.call_args.kwargs['check'])
        self.assertFalse(any('AllowUnauthenticated' in x or 'allow-unauthenticated' in x for x in command))

    def test_prepare_pins_transitive_versions_and_cache_key_changes_with_plan(self):
        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(environment, 'sources', return_value=[]), \
                patch.object(environment, 'candidates', side_effect=[{'lcov': '2.0-1'}, {'lcov': '2.0-1'}]), \
                patch.object(environment.subprocess, 'check_output', return_value='amd64\n'), \
                patch.object(environment, 'apt', side_effect=[None, 'Inst libdep (3.1-2 Ubuntu [amd64])\n',
                                                             None, 'Inst libdep (3.1-3 Ubuntu [amd64])\n']):
            root = Path(tmp); output = root / 'output'
            with patch.dict(os.environ, GITHUB_OUTPUT=str(output), ImageOS='ubuntu24', ImageVersion='same'):
                plan = environment.prepare(['lcov'], root / 'cache', root / 'state')
                environment.prepare(['lcov'], root / 'cache', root / 'state')
            self.assertEqual(plan['versions'], {'lcov': '2.0-1', 'libdep': '3.1-2'})
            keys = [x for x in output.read_text().splitlines() if x.startswith('key=')]
            prefixes = [x for x in output.read_text().splitlines() if x.startswith('restore-prefix=')]
            self.assertNotEqual(keys[0], keys[1])
            self.assertEqual(prefixes[0], prefixes[1])

    def test_index_failure_retries_another_signed_mirror(self):
        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(environment, 'sources', return_value=[]), \
                patch.object(environment, 'candidates', return_value={'lcov': '2.0-1'}), \
                patch.object(environment.subprocess, 'check_output', return_value='amd64\n'), \
                patch.object(environment, 'apt', side_effect=[subprocess.CalledProcessError(100, 'apt'), None, '']) as apt:
            root = Path(tmp)
            plan = environment.prepare(['lcov'], root / 'cache', root / 'state')
            self.assertEqual(plan['source_index'], 1)
            self.assertEqual([c.args[0][0] for c in apt.call_args_list], ['update', 'update', '--simulate'])

    def test_download_fallback_keeps_the_same_exact_versions_and_installs_offline(self):
        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(environment, 'sources', return_value=[]), \
                patch.object(environment, 'apt', side_effect=[subprocess.CalledProcessError(124, 'timeout'),
                                                             None, None, None]) as apt:
            root = Path(tmp)
            (root / 'plan.json').write_text(json.dumps({'source_index': 0, 'versions': {'libdep': '1:2.0-1'}}))
            environment.install(root / 'cache', root)
            commands = [c.args[0] for c in apt.call_args_list]
            self.assertEqual(commands[0], commands[2])
            self.assertEqual(commands[1], ['update'])
            self.assertIn('libdep=1:2.0-1', commands[0])
            self.assertIn('--download-only', commands[0])
            self.assertIn('--no-download', commands[3])
            self.assertTrue(all('--no-remove' in command for command in commands if command != ['update']))

    def test_configuration_failure_is_not_retried_as_a_mirror_problem(self):
        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(environment, 'sources', return_value=[]), \
                patch.object(environment, 'apt', side_effect=[None, subprocess.CalledProcessError(100, 'dpkg')]) as apt:
            root = Path(tmp)
            (root / 'plan.json').write_text(json.dumps({'source_index': 0, 'versions': {'lcov': '1'}}))
            with self.assertRaises(subprocess.CalledProcessError):
                environment.install(root / 'cache', root)
            self.assertEqual(apt.call_count, 2)

    def test_all_download_sources_failing_fails_the_job(self):
        with tempfile.TemporaryDirectory() as tmp, \
                patch.object(environment, 'sources', return_value=[]), \
                patch.object(environment, 'apt', side_effect=subprocess.CalledProcessError(100, 'apt')):
            root = Path(tmp)
            (root / 'plan.json').write_text(json.dumps({'source_index': 0, 'versions': {'lcov': '1'}}))
            with self.assertRaisesRegex(RuntimeError, 'All APT archive sources failed'):
                environment.install(root / 'cache', root)

    def test_cache_contains_only_archives_and_uses_resolved_version_key(self):
        action = yaml.safe_load((ROOT / '.github/actions/apt-dependencies/action.yml').read_text())
        steps = action['runs']['steps']
        restore = next(s for s in steps if s.get('uses') == 'actions/cache/restore@v4')
        save = next(s for s in steps if s.get('uses') == 'actions/cache/save@v4')
        self.assertEqual(restore['with']['path'], '${{ runner.temp }}/rodin-apt-archives/*.deb')
        self.assertEqual(restore['with']['key'], '${{ steps.plan.outputs.key }}')
        self.assertEqual(restore['with'], {**save['with'], 'restore-keys': '${{ steps.plan.outputs.restore-prefix }}'})
        self.assertLess(next(i for i, s in enumerate(steps) if s.get('id') == 'plan'), steps.index(restore))
        self.assertLess(next(i for i, s in enumerate(steps) if s.get('name') == 'Install resolved dependencies'),
                        steps.index(save))

    def test_hierarchy_and_execution_no_longer_delete_all_tools_unconditionally(self):
        for path in ['TestDirectoryBuild.yml', 'TestGroup.yml']:
            text = (ROOT / '.github/workflows' / path).read_text()
            self.assertNotIn('sudo rm -rf', text)
            self.assertNotIn('apt-get clean', text)
            self.assertIn('python3 dev/ci_environment.py disk', text)


if __name__ == '__main__':
    unittest.main()

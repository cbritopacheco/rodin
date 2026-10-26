"""A branch push or a stale successful Pages build is not deployment success."""
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml
import deploy_pages

ROOT = Path(__file__).resolve().parents[1]
COMMIT = 'a' * 40
SITE = {'build_type': 'legacy', 'source': {'branch': 'gh-pages', 'path': '/'},
        'html_url': 'https://example.test/rodin/'}


class PagesTest(unittest.TestCase):
    def test_waits_through_stale_build_and_verifies_exact_commit(self):
        with patch('deploy_pages.api', side_effect=[SITE, {'object': {'sha': COMMIT}}, {},
                {'commit': 'old', 'status': 'built'}, {'commit': COMMIT, 'status': 'building'},
                {'commit': COMMIT, 'status': 'built'}]) as api, \
                patch('deploy_pages.time.sleep'):
            self.assertEqual(deploy_pages.deploy('owner/repo'), COMMIT)
            self.assertIn(unittest.mock.call('owner/repo', 'pages/builds', 'POST'), api.call_args_list)

    def test_failed_current_build_is_not_success(self):
        with patch('deploy_pages.api', side_effect=[SITE, {'object': {'sha': COMMIT}}, {},
                {'commit': COMMIT, 'status': 'errored', 'error': {'message': 'bad build'}}]):
            with self.assertRaisesRegex(RuntimeError, 'bad build'):
                deploy_pages.deploy('owner/repo')

    def test_no_confirmed_build_times_out(self):
        with patch('deploy_pages.api', side_effect=[SITE, {'object': {'sha': COMMIT}}, {}]), \
                patch('deploy_pages.time.monotonic', side_effect=[0, 1]):
            with self.assertRaises(TimeoutError):
                deploy_pages.deploy('owner/repo', timeout=1)

    def test_changed_branch_rejected_before_build_request(self):
        with patch('deploy_pages.api', side_effect=[SITE, {'object': {'sha': COMMIT}}]) as api, \
                patch('deploy_pages.subprocess.check_output', return_value='different\n'):
            with self.assertRaisesRegex(ValueError, 'changed'):
                deploy_pages.deploy('owner/repo', Path('checkout'))
            self.assertEqual(api.call_count, 2)

    def test_unexpected_source_rejected(self):
        with patch('deploy_pages.api', return_value=dict(SITE, build_type='workflow')):
            with self.assertRaisesRegex(ValueError, 'legacy'):
                deploy_pages.deploy('owner/repo')

    def test_every_pages_writer_has_same_lock_permissions_and_explicit_deployment(self):
        for filename, job in [('Benchmarks.yml', 'Publish'), ('StaticAnalysis.yml', 'IncludeDependencyGraph')]:
            workflow = yaml.safe_load((ROOT / '.github/workflows' / filename).read_text())
            definition = workflow['jobs'][job]
            self.assertEqual(definition['concurrency']['group'], 'rodin-gh-pages-publication')
            permissions = definition.get('permissions', workflow.get('permissions', {}))
            self.assertEqual(permissions['pages'], 'write')
            self.assertTrue(any('dev/deploy_pages.py' in step.get('run', '') for step in definition['steps']))

    def test_shared_apt_guard_precedes_version_resolution(self):
        action = yaml.safe_load((ROOT / '.github/actions/apt-dependencies/action.yml').read_text())
        steps = action['runs']['steps']
        guard = next(i for i, step in enumerate(steps) if step.get('uses') == './.github/actions/prepare-apt')
        prepare = next(i for i, step in enumerate(steps) if 'prepare --packages' in step.get('run', ''))
        self.assertLess(guard, prepare)


if __name__ == '__main__':
    unittest.main()

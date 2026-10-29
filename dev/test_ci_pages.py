"""A branch push or a stale successful Pages build is not deployment success."""
import unittest
import json
import tempfile
from pathlib import Path
from unittest.mock import patch
from urllib.error import URLError

import yaml
import deploy_pages

ROOT = Path(__file__).resolve().parents[1]
COMMIT = 'a' * 40
SITE = {'build_type': 'legacy', 'source': {'branch': 'gh-pages', 'path': '/'},
        'html_url': 'https://example.test/rodin/'}


class PagesTest(unittest.TestCase):
    def marker(self):
        directory = tempfile.TemporaryDirectory(prefix='rodin-pages-verification-')
        self.addCleanup(directory.cleanup)
        checkout = Path(directory.name)
        marker = {'commit': COMMIT, 'repository': 'owner/repo'}
        (checkout / deploy_pages.MARKER).write_text(json.dumps(marker))
        return checkout, marker

    def test_prepared_artifact_identifies_exact_published_commit(self):
        checkout, marker = self.marker()
        with patch('deploy_pages.api', side_effect=[dict(SITE, build_type='workflow'),
                {'object': {'sha': COMMIT}}]), \
                patch('deploy_pages.subprocess.check_output', return_value=COMMIT + '\n'):
            self.assertEqual(deploy_pages.prepare('owner/repo', checkout), marker)
            self.assertEqual(json.loads((checkout / deploy_pages.MARKER).read_text()), marker)

    def test_waits_through_stale_content_and_verifies_exact_commit(self):
        checkout, marker = self.marker()
        with patch('deploy_pages.fetch_marker', side_effect=[dict(marker, commit='old'), marker]) as fetch, \
                patch('deploy_pages.time.sleep'):
            self.assertEqual(deploy_pages.verify(checkout, SITE['html_url']), COMMIT)
            self.assertEqual(fetch.call_count, 2)

    def test_unavailable_or_stale_content_is_not_success(self):
        checkout, marker = self.marker()
        for content in [dict(marker, commit='old'), ValueError('not JSON'), URLError('not yet available')]:
            with self.subTest(content=content), \
                    patch('deploy_pages.fetch_marker', side_effect=[content]), \
                    patch('deploy_pages.time.monotonic', side_effect=[0, 0, 2]), \
                    patch('deploy_pages.time.sleep'):
                with self.assertRaises(TimeoutError):
                    deploy_pages.verify(checkout, SITE['html_url'], timeout=1)

    def test_verification_rejects_insecure_site_url(self):
        checkout, _ = self.marker()
        with self.assertRaisesRegex(ValueError, 'HTTPS'):
            deploy_pages.verify(checkout, 'http://example.test/rodin/')

    def test_changed_branch_rejected_before_marker_written(self):
        checkout, _ = self.marker()
        (checkout / deploy_pages.MARKER).unlink()
        with patch('deploy_pages.api', side_effect=[dict(SITE, build_type='workflow'),
                {'object': {'sha': COMMIT}}]), \
                patch('deploy_pages.subprocess.check_output', return_value='different\n'):
            with self.assertRaisesRegex(ValueError, 'changed'):
                deploy_pages.prepare('owner/repo', checkout)
            self.assertFalse((checkout / deploy_pages.MARKER).exists())

    def test_legacy_source_requires_explicit_settings_change(self):
        with patch('deploy_pages.api', return_value=SITE):
            with self.assertRaisesRegex(ValueError, 'GitHub Actions'):
                deploy_pages.prepare('owner/repo', Path('checkout'))

    def test_every_pages_writer_deploys_its_artifact_under_same_lock(self):
        for filename, job in [('Documentation.yml', 'PublishDocumentation'),
                              ('Benchmarks.yml', 'Publish')]:
            workflow = yaml.safe_load((ROOT / '.github/workflows' / filename).read_text())
            definition = workflow['jobs'][job]
            self.assertEqual(definition['concurrency']['group'], 'rodin-gh-pages-publication')
            permissions = definition.get('permissions', workflow.get('permissions', {}))
            self.assertEqual(permissions['pages'], 'write')
            self.assertEqual(permissions['id-token'], 'write')
            self.assertEqual(definition['environment']['name'], 'Documentation and benchmarks')
            self.assertTrue(any(step.get('uses') == './.github/actions/publish-pages'
                                for step in definition['steps']))
        action = yaml.safe_load((ROOT / '.github/actions/publish-pages/action.yml').read_text())
        steps = action['runs']['steps']
        self.assertEqual(steps[1]['uses'], 'actions/upload-pages-artifact@v4')
        self.assertEqual(steps[2]['uses'], 'actions/deploy-pages@v4')
        self.assertIn(' prepare ', steps[0]['run'])
        self.assertIn(' verify ', steps[3]['run'])

    def test_documentation_is_independent_and_pull_requests_do_not_deploy(self):
        docs = yaml.safe_load((ROOT / '.github/workflows/Documentation.yml').read_text())
        benchmarks = yaml.safe_load((ROOT / '.github/workflows/Benchmarks.yml').read_text())
        events = docs.get('on', docs.get(True))
        self.assertEqual(events['push']['branches'], ['master', 'develop'])
        self.assertIn('pull_request', events)
        publication = docs['jobs']['PublishDocumentation']
        self.assertEqual(publication['needs'], 'Build-And-Deploy')
        self.assertEqual(publication['if'], "github.event_name == 'push'")
        self.assertNotIn('Documentation', benchmarks['jobs'])
        self.assertNotIn('PublishDocumentation', benchmarks['jobs'])
        self.assertEqual(benchmarks['jobs']['Publish']['needs'], 'Pool')

    def test_shared_apt_guard_precedes_version_resolution(self):
        action = yaml.safe_load((ROOT / '.github/actions/apt-dependencies/action.yml').read_text())
        steps = action['runs']['steps']
        guard = next(i for i, step in enumerate(steps) if step.get('uses') == './.github/actions/prepare-apt')
        prepare = next(i for i, step in enumerate(steps) if 'prepare --packages' in step.get('run', ''))
        self.assertLess(guard, prepare)


if __name__ == '__main__':
    unittest.main()

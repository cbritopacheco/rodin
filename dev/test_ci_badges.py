"""API conclusions, current-commit isolation and badge publication contracts."""
import base64
import json
import unittest

from ci_badges import Publisher, badge, state


class BadgeTest(unittest.TestCase):
    def test_status_colours(self):
        cases = [('queued', None, 'queued', '2874c7'), ('in_progress', None, 'running', 'd49b16'),
                 ('completed', 'success', 'success', '2e8540'), ('completed', 'failure', 'failure', 'c62828'),
                 ('completed', 'timed_out', 'timed_out', 'c62828'),
                 ('completed', 'cancelled', 'cancelled', 'e87924'),
                 ('completed', 'skipped', 'skipped', 'e87924'),
                 ('completed', None, 'unknown', '808080')]
        for status, conclusion, message, color in cases:
            with self.subTest(status=status, conclusion=conclusion):
                record = badge('Tests', state(status, conclusion), 'abc')
                self.assertEqual(record['message'], message.replace('_', ' '))
                self.assertEqual(record['color'], color)

    def fixture(self, source='abc', documentation=None):
        calls = []
        def api(path, method='GET', data=None):
            calls.append((path, method, data))
            if '/git/ref/heads/' in path:
                return {'object': {'sha': 'abc'}}
            if '/actions/workflows/' in path:
                return {'workflow_runs': [dict(id=10, head_sha=source, status='completed', conclusion='cancelled')]}
            if '/jobs?' in path:
                return {'jobs': [] if documentation is None else [documentation]}
            if '/contents/' in path and method == 'GET':
                return {'sha': 'old', 'content': base64.b64encode(b'{}').decode()}
            if method == 'PUT':
                return {}
            raise AssertionError(path)
        return Publisher('owner/repo', 'develop', api), calls

    def test_cancelled_runs_do_not_inherit_github_badge_failure(self):
        publisher, _ = self.fixture()
        records = publisher.collect('abc')
        self.assertEqual(records['Build']['message'], 'cancelled')
        self.assertEqual(records['Tests']['color'], 'e87924')
        self.assertEqual(records['Documentation']['message'], 'unknown')

    def test_previous_commit_success_is_unknown(self):
        publisher, _ = self.fixture(source='old')
        self.assertTrue(all(x['message'] == 'unknown' for x in publisher.collect('abc').values()))

    def test_documentation_uses_publication_job_independently(self):
        publisher, _ = self.fixture(documentation=dict(name='Publish documentation', status='completed', conclusion='success'))
        self.assertEqual(publisher.collect('abc')['Documentation']['message'], 'success')
        self.assertEqual(publisher.collect('abc')['Benchmarks']['message'], 'cancelled')

    def test_verified_deployment_override_is_tied_to_source_run_and_commit(self):
        publisher, _ = self.fixture()
        self.assertEqual(publisher.collect('abc', 'success', '10', 'abc')['Documentation']['message'], 'success')
        self.assertEqual(publisher.collect('abc', 'success', '11', 'abc')['Documentation']['message'], 'unknown')
        self.assertEqual(publisher.collect('abc', 'success', '10', 'old')['Documentation']['message'], 'unknown')

    def test_only_badge_branch_files_are_written(self):
        publisher, calls = self.fixture()
        publisher.publish()
        writes = [(path, data) for path, method, data in calls if method == 'PUT']
        self.assertEqual(len(writes), 4)
        for path, data in writes:
            self.assertTrue(path.startswith('repos/owner/repo/contents/develop/'))
            self.assertEqual(data['branch'], 'ci-badges')
            self.assertEqual(data['sha'], 'old')
            self.assertEqual(json.loads(base64.b64decode(data['content']))['commit'], 'abc')

    def test_branch_advance_prevents_writing_old_results(self):
        publisher, calls = self.fixture()
        heads = iter(['abc', 'new'])
        publisher.head = lambda: next(heads)
        publisher.publish()
        self.assertFalse(any(method == 'PUT' for _, method, _ in calls))

    def test_non_branch_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            Publisher('owner/repo', '../other')


if __name__ == '__main__':
    unittest.main()

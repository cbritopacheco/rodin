"""Compact reports distinguish complete evidence from partial progress."""
import unittest
import subprocess
from unittest.mock import patch
from ci_workflow_summary import render, workflow_jobs


class WorkflowSummaryTest(unittest.TestCase):
    def setUp(self):
        self.run = dict(name='Tests', id=7, head_sha='abc', head_branch='develop', run_attempt=2, conclusion='failure')
        self.jobs = [dict(id=9, name='Geometry <Compile>|Debug', status='completed', conclusion='failure',
                          steps=[dict(name='Compile', conclusion='failure')])]

    def test_rest_inventory_keeps_every_page(self):
        with patch('ci_workflow_summary.api', return_value=[{'jobs': self.jobs}, {'jobs': [dict(id=10)]}]):
            self.assertEqual(len(workflow_jobs(self.run, 'owner/repo')), 2)

    def test_large_matrix_fallback_preserves_latest_ids_and_start_times(self):
        self.run['check_suite_node_id'] = 'suite'
        pages = iter([
            {'nodes': [dict(databaseId=8, name='Compile', status='COMPLETED', conclusion='FAILURE', startedAt='old')],
             'pageInfo': {'hasNextPage': True, 'endCursor': 'next'}},
            {'nodes': [dict(databaseId=9, name='Compile', status='COMPLETED', conclusion='SUCCESS', startedAt='new'),
                       dict(databaseId=10, name='Coverage / Finalize', status='COMPLETED', conclusion='FAILURE', startedAt='new')],
             'pageInfo': {'hasNextPage': False}}])
        cursors = []
        def request(endpoint, method, data):
            cursors.append(data['variables']['cursor'])
            return {'data': {'node': {'checkRuns': next(pages)}}}
        error = subprocess.CalledProcessError(1, ['gh'])
        with patch('ci_workflow_summary.api', side_effect=error), patch('ci_badges.Publisher.github', side_effect=request):
            jobs = workflow_jobs(self.run, 'owner/repo')
        self.assertEqual(cursors, [None, 'next'])
        self.assertEqual([j['id'] for j in jobs], [9, 10])
        self.assertEqual(jobs[0]['started_at'], 'new')
        text = render(self.run, jobs, 'owner/repo')
        self.assertIn('| 2 | 2 | 1 | 1 |', text)
        self.assertIn('/job/10', text)
        self.assertIn('**Coverage / Finalize:** failure', text)

    def test_unavailable_inventories_do_not_report_zero_jobs(self):
        error = subprocess.CalledProcessError(1, ['gh'])
        with patch('ci_workflow_summary.api', side_effect=error), patch('ci_badges.Publisher.checks', return_value=None):
            with self.assertRaises(subprocess.CalledProcessError):
                workflow_jobs(self.run, 'owner/repo')

    def test_failure_has_log_link_and_does_not_invent_coverage(self):
        text = render(self.run, self.jobs, 'owner/repo')
        self.assertIn('/actions/runs/7/job/9', text)
        self.assertIn('&lt;Compile&gt;&#124;Debug', text)
        self.assertIn('Compile', text)
        self.assertIn('report unavailable', text)
        self.assertNotIn('0.00%', text)

    def test_complete_report_keeps_gate_status_explicit(self):
        self.jobs.append(dict(id=10, name='Coverage / Finalize', status='completed', conclusion='success'))
        text = render(self.run, self.jobs, 'owner/repo', dict(hit=96, lines=100))
        self.assertIn('96/100 (96.00%)', text)
        self.assertIn('**Coverage / Finalize:** success', text)

    def test_benchmarks_separate_registered_lanes_from_completeness(self):
        self.run['name'] = 'Benchmarks'
        jobs = [dict(id=10, name='Local / Solver/NP1/0', status='completed', conclusion='success'),
                dict(id=11, name='Validate complete campaign', status='completed', conclusion='failure')]
        text = render(self.run, jobs, 'owner/repo')
        self.assertIn('1/1 registered lanes passed', text)
        self.assertIn('**Complete campaign validation / pooling:** failure', text)

    def test_site_link_requires_successful_publication(self):
        self.run['name'] = 'Documentation'
        job = dict(id=12, name='Publish documentation', status='completed', conclusion='failure')
        self.assertNotIn('[Published site]', render(self.run, [job], 'owner/repo'))
        job['conclusion'] = 'success'
        self.assertIn('/develop/index.html', render(self.run, [job], 'owner/repo'))

    def test_large_failure_lists_are_bounded(self):
        self.assertIn('10 additional failures', render(self.run, self.jobs * 40, 'owner/repo'))


if __name__ == '__main__':
    unittest.main()

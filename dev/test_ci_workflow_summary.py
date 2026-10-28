"""Compact reports distinguish complete evidence from partial progress."""
import unittest
from ci_workflow_summary import render


class WorkflowSummaryTest(unittest.TestCase):
    def setUp(self):
        self.run = dict(name='Tests', id=7, head_sha='abc', head_branch='develop', run_attempt=2, conclusion='failure')
        self.jobs = [dict(id=9, name='Geometry <Compile>|Debug', status='completed', conclusion='failure',
                          steps=[dict(name='Compile', conclusion='failure')])]

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
        self.run['name'] = 'Publish documentation manually'
        job = dict(id=12, name='Publish documentation', status='completed', conclusion='failure')
        self.assertNotIn('[Published site]', render(self.run, [job], 'owner/repo'))
        job['conclusion'] = 'success'
        self.assertIn('/develop/index.html', render(self.run, [job], 'owner/repo'))

    def test_large_failure_lists_are_bounded(self):
        self.assertIn('10 additional failures', render(self.run, self.jobs * 40, 'owner/repo'))


if __name__ == '__main__':
    unittest.main()

"""Verify automatic CI routing, draft skips and post-failure finalization."""
import ast
import re
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / '.github/workflows'
MAIN_BRANCHES = ['master', 'develop']
PR_EVENTS = {'opened', 'synchronize', 'reopened', 'ready_for_review', 'converted_to_draft'}


def workflow(name):
    return yaml.safe_load((WORKFLOWS / name).read_text())


def events(document):
    return document.get('on', document.get(True, {}))


def permits(job, event, draft=False, ref='refs/heads/develop', cancelled=False,
            upstream_event='push'):
    """Evaluate the boolean subset used by these top-level job conditions."""
    expression = str(job.get('if', 'true')).removeprefix('${{').removesuffix('}}').strip()
    values = {'github.event.pull_request.draft': draft,
              'github.event.workflow_run.event': upstream_event,
              'github.event_name': event, 'github.ref': ref}
    for token, value in values.items():
        expression = expression.replace(token, repr(value))
    expression = expression.replace('always()', 'True').replace('success()', 'True')
    expression = expression.replace('cancelled()', repr(cancelled))
    expression = expression.replace('&&', ' and ').replace('||', ' or ')
    expression = re.sub(r'!(?!=)', 'not ', expression)
    expression = re.sub(r'\btrue\b', 'True', expression)
    expression = re.sub(r'\bfalse\b', 'False', expression)
    def evaluate(node):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.BoolOp):
            values = [bool(evaluate(value)) for value in node.values]
            return all(values) if isinstance(node.op, ast.And) else any(values)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
            return not evaluate(node.operand)
        if isinstance(node, ast.Compare) and len(node.ops) == 1:
            left, right = evaluate(node.left), evaluate(node.comparators[0])
            if isinstance(node.ops[0], ast.Eq):
                return left == right
            if isinstance(node.ops[0], ast.NotEq):
                return left != right
        raise AssertionError('Unsupported trigger expression: ' + ast.dump(node))

    return bool(evaluate(ast.parse(expression, mode='eval').body))


class TriggerPolicyTest(unittest.TestCase):
    def automatic_workflows(self):
        for path in sorted(WORKFLOWS.glob('*.yml')):
            document = workflow(path.name)
            if set(events(document)) & {'push', 'pull_request', 'pull_request_target'}:
                yield path.name, document

    def test_pushes_are_restricted_to_main_branches(self):
        for name, document in self.automatic_workflows():
            if 'push' in events(document):
                with self.subTest(workflow=name):
                    self.assertEqual(events(document)['push']['branches'], MAIN_BRANCHES)
                    self.assertNotIn('branches-ignore', events(document)['push'])

    def test_pr_events_cover_ready_updates_and_return_to_draft(self):
        for name, document in self.automatic_workflows():
            for event in ['pull_request', 'pull_request_target']:
                if event not in events(document):
                    continue
                with self.subTest(workflow=name, event=event):
                    config = events(document)[event]
                    self.assertEqual(config['branches'], MAIN_BRANCHES)
                    self.assertEqual(set(config['types']), PR_EVENTS)

    def test_every_draft_job_is_skipped_including_finalizers_and_matrices(self):
        for name, document in self.automatic_workflows():
            event = 'pull_request_target' if 'pull_request_target' in events(document) else 'pull_request'
            for job_name, job in document['jobs'].items():
                with self.subTest(workflow=name, job=job_name):
                    self.assertFalse(permits(job, event, draft=True, ref='refs/pull/400/merge'))

    def test_ready_prs_run_all_validation_jobs(self):
        deployment_jobs = {('Documentation.yml', 'PublishDocumentation'), ('Benchmarks.yml', 'Publish')}
        for name, document in self.automatic_workflows():
            event = 'pull_request_target' if 'pull_request_target' in events(document) else 'pull_request'
            for job_name, job in document['jobs'].items():
                with self.subTest(workflow=name, job=job_name):
                    self.assertEqual(permits(job, event, ref='refs/pull/400/merge'),
                                     (name, job_name) not in deployment_jobs)

    def test_main_branch_pushes_run_validation_and_publication(self):
        for name, document in self.automatic_workflows():
            if 'push' not in events(document):
                continue
            for branch in MAIN_BRANCHES:
                for job_name, job in document['jobs'].items():
                    with self.subTest(workflow=name, job=job_name, branch=branch):
                        self.assertTrue(permits(job, 'push', ref='refs/heads/' + branch))

    def test_main_validation_workflows_have_unfiltered_pr_updates(self):
        for name in ['Build.yml', 'Tests.yml', 'Style.yml', 'Documentation.yml', 'Benchmarks.yml', 'Preflight.yml']:
            with self.subTest(workflow=name):
                self.assertIn('pull_request', events(workflow(name)))
                self.assertNotIn('paths', events(workflow(name))['pull_request'])

    def test_returning_to_draft_can_supersede_the_existing_pr_campaign(self):
        for name, document in self.automatic_workflows():
            with self.subTest(workflow=name):
                self.assertTrue(document['concurrency']['cancel-in-progress'])
                self.assertIn('github.', document['concurrency']['group'])

    def test_failed_ready_campaign_still_runs_mandatory_finalization(self):
        coverage = workflow('Tests.yml')['jobs']['CoverageFinalize']
        pooling = workflow('Benchmarks.yml')['jobs']['Pool']
        for job in [coverage, pooling]:
            self.assertIn('always()', job['if'])
            self.assertTrue(permits(job, 'pull_request'))
            self.assertFalse(permits(job, 'pull_request', draft=True))
        for name, job in workflow('Tests.yml')['jobs'].items():
            if name.startswith(('Configurations_', 'Installation')):
                self.assertTrue(permits(job, 'pull_request'))
                self.assertFalse(permits(job, 'pull_request', cancelled=True))

    def test_benchmark_publication_is_limited_to_main_branches(self):
        publish = workflow('Benchmarks.yml')['jobs']['Publish']
        self.assertFalse(permits(publish, 'pull_request', ref='refs/pull/400/merge'))
        self.assertFalse(permits(publish, 'workflow_dispatch', ref='refs/heads/topic'))
        self.assertTrue(permits(publish, 'workflow_dispatch'))

    def test_reports_and_badges_ignore_pr_workflow_runs(self):
        recovery = workflow('RetryInterruptedJobs.yml')['jobs']['Recover']
        badges = workflow('CIBadges.yml')['jobs']['Publish']
        for job in [recovery, badges]:
            self.assertFalse(permits(job, 'workflow_run', upstream_event='pull_request'))
            self.assertTrue(permits(job, 'workflow_run'))
        self.assertTrue(permits(badges, 'workflow_dispatch', upstream_event='pull_request'))


if __name__ == '__main__':
    unittest.main()

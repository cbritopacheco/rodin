"""Failure summaries retain diagnostics without publishing partial successes."""

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import ci_summary


class SummaryTest(unittest.TestCase):
    def test_failure_before_execution_has_metadata_and_safe_links(self):
        with tempfile.TemporaryDirectory() as folder:
            previous = Path.cwd()
            os.chdir(folder)
            try:
                environment = {
                    'BENCHMARK_LANE': json.dumps({'id': 'lane', 'label': '<script>|lane',
                                                 'entry': {'target': 'Solver'}, 'ranks': 4,
                                                 'configuration': 'distributed', 'cases': ['Stokes/3']}),
                    'BENCHMARK_STEPS': json.dumps({'environment': {'outcome': 'failure'}}),
                    'GITHUB_SHA': 'abc', 'GITHUB_RUN_ATTEMPT': '2', 'GITHUB_RUN_ID': '123',
                    'GITHUB_REPOSITORY': 'owner/repo', 'GITHUB_STEP_SUMMARY': str(Path(folder) / 'summary.md'),
                }
                with patch.dict(os.environ, environment):
                    ci_summary.benchmark()
                summary = Path('summary.md').read_text()
                self.assertIn('&lt;script&gt;&#124;lane', summary)
                self.assertIn('environment', summary)
                self.assertIn('distributed / 4', summary)
                self.assertIn('/actions/runs/123/attempts/2', summary)
                self.assertIn('pending inspection', summary)
                context = json.loads(Path('benchmark-diagnostics/context.json').read_text())
                self.assertEqual(context['lane']['cases'], ['Stokes/3'])
                self.assertFalse(Path('benchmark-result/receipt.json').exists())
                environment['BENCHMARK_STEPS'] = json.dumps({'execution': {'outcome': 'success'}})
                with patch.dict(os.environ, environment):
                    ci_summary.benchmark()
                self.assertEqual(Path('summary.md').read_text(), summary)
            finally:
                os.chdir(previous)

    def test_diagnostic_steps_do_not_weaken_execution_or_pooling(self):
        import yaml
        root = Path(__file__).resolve().parents[1]
        data = yaml.safe_load((root / '.github/workflows/BenchmarkGroup.yml').read_text())
        steps = data['jobs']['Run']['steps']
        execution = next(s for s in steps if s.get('id') == 'execution')
        self.assertNotIn('continue-on-error', execution)
        self.assertIn('set -o pipefail', execution['run'])
        self.assertIn('tee benchmark-diagnostics/execution.log', execution['run'])
        self.assertIn('/usr/bin/time -v', execution['run'])
        diagnostics = next(s for s in steps if s.get('id') == 'diagnostics')
        self.assertIn('failure()', diagnostics['if'])
        self.assertIn('cancelled()', diagnostics['if'])
        self.assertIn('github.run_attempt', diagnostics['with']['name'])
        results = next(s for s in steps if s.get('with', {}).get('name', '').startswith('benchmark-result-'))
        self.assertNotIn('if', results)
        self.assertEqual(results['with']['if-no-files-found'], 'error')


if __name__ == '__main__':
    unittest.main()

"""Publication preserves component histories and rejects partial campaigns."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from benchmark_ci import dump, configuration
from benchmark_publish import publish, HISTORY_COMMITS
from ci_builds import digest

COMMIT = 'c' * 40


class BenchmarkPublicationTest(unittest.TestCase):
    def fixture(self, root):
        site, pages, results = root / 'site', root / 'pages', root / 'results'
        component = 'Rodin/Math/SpatialMatrix'
        page = site / 'benchmarks' / component / 'index.html'
        page.parent.mkdir(parents=True)
        page.write_text('<div data-benchmark-history></div>')
        (site / 'benchmarks/index.html').write_text('Benchmark navigation')
        record = {'schema': 1, 'commit': COMMIT, 'configuration': 'config-np1',
                  'component': component, 'format': 'googlecpp', 'ranks': 1,
                  'definition': configuration('local'),
                  'results': [{'name': 'Product/3', 'metric': 'cpu_time', 'unit': 'ns',
                               'value': 20, 'samples': [20]}]}
        path = results / component / 'config-np1/results.json'
        dump(path, record)
        dump(results / 'complete.json', {'schema': 1, 'commit': COMMIT,
             'configurations': ['distributed', 'local'], 'cases': 1,
             'components': {str(path.relative_to(results)): digest(path)}})
        old = pages / 'benchmarks/refs/heads/develop'
        old.mkdir(parents=True)
        (old / 'index.html').write_text('Retired aggregate dashboard')
        return site, pages, results

    def test_mcss_page_and_histories_publish_together_and_retire_aggregate(self):
        with tempfile.TemporaryDirectory() as folder, patch('benchmark_publish.subprocess.check_output', return_value=COMMIT):
            root = Path(folder)
            site, pages, results = self.fixture(root)
            self.assertEqual(publish(results, site, pages, 'refs/heads/develop', '123.1'), 1)
            page = pages / 'docs/refs/heads/develop/benchmarks/Rodin/Math/SpatialMatrix/index.html'
            markup = page.read_text()
            self.assertIn('data-history-index=', markup)
            history = next((pages / 'benchmark-data/refs/heads/develop/Rodin/Math/SpatialMatrix/config-np1/cases').glob('*.json'))
            self.assertEqual(json.loads(history.read_text())['entries'][0]['run'], '123.1')
            self.assertFalse((pages / 'benchmarks').exists())
            self.assertEqual(json.loads((history.parent.parent.parent / 'index.json').read_text())['configurations'][0]['id'], 'config-np1')

    def test_republishing_same_commit_is_idempotent_and_other_refs_survive(self):
        with tempfile.TemporaryDirectory() as folder, patch('benchmark_publish.subprocess.check_output', return_value=COMMIT):
            root = Path(folder)
            site, pages, results = self.fixture(root)
            master = pages / 'docs/refs/heads/master/index.html'
            master.parent.mkdir(parents=True)
            master.write_text('Other branch documentation')
            publish(results, site, pages, 'refs/heads/develop', '123.1')
            publish(results, site, pages, 'refs/heads/develop', '123.2')
            history = next((pages / 'benchmark-data').glob('**/cases/*.json'))
            entries = json.loads(history.read_text())['entries']
            self.assertEqual(len(entries), 1)
            self.assertEqual(entries[0]['run'], '123.2')
            self.assertEqual(master.read_text(), 'Other branch documentation')

    def test_case_histories_are_independent_and_have_a_bounded_window(self):
        with tempfile.TemporaryDirectory() as folder, patch('benchmark_publish.subprocess.check_output', return_value=COMMIT):
            root=Path(folder);site,pages,results=self.fixture(root)
            publish(results,site,pages,'refs/heads/develop','1')
            history=next((pages / 'benchmark-data').glob('**/cases/*.json'))
            data=json.loads(history.read_text())
            data['entries']=[dict(commit=str(i),results=[]) for i in range(HISTORY_COMMITS + 5)]
            dump(history,data)
            publish(results,site,pages,'refs/heads/develop','2')
            data=json.loads(history.read_text())
            self.assertEqual(len(data['entries']),HISTORY_COMMITS)
            self.assertEqual(data['entries'][-1]['commit'],COMMIT)
            index=json.loads((history.parent.parent / 'index.json').read_text())
            self.assertEqual(index['cases'][0]['name'],'Product/3')
            self.assertNotIn('entries',index)

    def test_changed_or_missing_component_is_rejected_before_mutation(self):
        for action in ('change', 'remove'):
            with self.subTest(action=action), tempfile.TemporaryDirectory() as folder, patch('benchmark_publish.subprocess.check_output', return_value=COMMIT):
                root = Path(folder)
                site, pages, results = self.fixture(root)
                path = next(results.glob('**/results.json'))
                if action == 'remove': path.unlink()
                else: path.write_text('{}')
                with self.assertRaises(ValueError):
                    publish(results, site, pages, 'refs/heads/develop', '123.1')
                self.assertTrue((pages / 'benchmarks/refs/heads/develop/index.html').is_file())
                self.assertFalse((pages / 'docs').exists())

    def test_unsafe_ref_and_stale_commit_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder, patch('benchmark_publish.subprocess.check_output', return_value=COMMIT):
            root = Path(folder)
            site, pages, results = self.fixture(root)
            with self.assertRaises(ValueError): publish(results, site, pages, '../escape', '1')
            complete = json.loads((results / 'complete.json').read_text())
            complete['commit'] = 'd' * 40
            dump(results / 'complete.json', complete)
            with self.assertRaises(ValueError): publish(results, site, pages, 'refs/heads/develop', '1')
            self.assertTrue((pages / 'benchmarks').exists())


if __name__ == '__main__':
    unittest.main()

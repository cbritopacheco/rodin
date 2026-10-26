"""Publication preserves component histories and rejects partial campaigns."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from benchmark_ci import dump, configuration
from benchmark_publish import publish, publish_documentation, HISTORY_COMMITS
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

    def test_documentation_deploys_without_a_campaign(self):
        with tempfile.TemporaryDirectory() as folder:
            site, pages, results = self.fixture(Path(folder))
            (site / 'index.html').write_text('New documentation')
            publish_documentation(site, pages, 'refs/heads/develop')
            destination = pages / 'docs/refs/heads/develop'
            self.assertEqual((destination / 'index.html').read_text(), 'New documentation')
            self.assertFalse((pages / 'benchmark-data').exists())
            self.assertNotIn('data-history-index', (destination / 'benchmarks/Rodin/Math/SpatialMatrix/index.html').read_text())

    def test_new_documentation_keeps_existing_complete_measurements(self):
        with tempfile.TemporaryDirectory() as folder, patch('benchmark_publish.subprocess.check_output', return_value=COMMIT):
            site, pages, results = self.fixture(Path(folder))
            publish(results, site, pages, 'refs/heads/develop', '1')
            data = pages / 'benchmark-data'
            before = {str(path.relative_to(data)): path.read_bytes() for path in data.rglob('*.json')}
            (site / 'index.html').write_text('Newer documentation')
            for _ in range(2):
                publish_documentation(site, pages, 'refs/heads/develop')
            self.assertEqual(before, {str(path.relative_to(data)): path.read_bytes() for path in data.rglob('*.json')})
            page = pages / 'docs/refs/heads/develop/benchmarks/Rodin/Math/SpatialMatrix/index.html'
            self.assertEqual(page.read_text().count('data-history-index='), 1)
            import re
            from html import unescape
            url = unescape(re.search(r'data-history-index="([^"]+)"', page.read_text()).group(1))
            self.assertTrue((page.parent / url).resolve().is_file())

    def test_campaign_updates_published_site_without_replacing_documentation(self):
        with tempfile.TemporaryDirectory() as folder, patch('benchmark_publish.subprocess.check_output', return_value=COMMIT):
            site, pages, results = self.fixture(Path(folder))
            (site / 'index.html').write_text('Documentation before measurements')
            publish_documentation(site, pages, 'refs/heads/develop')
            destination = pages / 'docs/refs/heads/develop'
            (destination / 'index.html').write_text('Newest documentation must survive')
            publish(results, destination, pages, 'refs/heads/develop', '1')
            self.assertEqual((destination / 'index.html').read_text(), 'Newest documentation must survive')
            self.assertIn('data-history-index=', (destination / 'benchmarks/Rodin/Math/SpatialMatrix/index.html').read_text())
            self.assertEqual(json.loads((pages / 'benchmark-data/refs/heads/develop/complete.json').read_text())['commit'], COMMIT)
            next(results.glob('**/results.json')).unlink()
            with self.assertRaises(ValueError):
                publish(results, destination, pages, 'refs/heads/develop', '2')
            self.assertEqual((destination / 'index.html').read_text(), 'Newest documentation must survive')

    def test_invalid_documentation_cannot_erase_published_site(self):
        with tempfile.TemporaryDirectory() as folder:
            site, pages, results = self.fixture(Path(folder))
            publish_documentation(site, pages, 'refs/heads/develop')
            destination = pages / 'docs/refs/heads/develop'
            with self.assertRaises(ValueError):
                publish_documentation(destination, pages, 'refs/heads/develop')
            with self.assertRaises(ValueError):
                publish_documentation(site, pages, '../escape')
            (site / 'benchmarks/index.html').unlink()
            with self.assertRaises(ValueError):
                publish_documentation(site, pages, 'refs/heads/develop')
            self.assertTrue((destination / 'benchmarks/index.html').is_file())

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

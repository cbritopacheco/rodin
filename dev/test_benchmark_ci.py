"""Complete campaigns require every discovered case and verified raw contribution."""
import csv
import io
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
from unittest.mock import patch

import benchmark_ci as ci

COMMIT = 'c' * 40


class BenchmarkCampaignTest(unittest.TestCase):
    def report(self, name='Kernel/3', **changes):
        return {'context': {'num_cpus': 2}, 'benchmarks': [dict(
            name=name, run_name=name, iterations=20, cpu_time=10, real_time=30,
            time_unit='ns', **changes)]}

    def fixture(self, root):
        commons, groups, contributions = [], [], []
        for mode in ('local', 'distributed'):
            entry = {'target': mode + 'Benchmark', 'path': 'Rodin/' + ('Math/Kernel' if mode == 'local' else 'MPI/Location/Kernel'),
                     'state': 'registered', 'format': 'googlecpp', 'metric': 'auto',
                     'profile': 'regular', 'ranks': '1' if mode == 'local' else '2;4', 'args': '', 'mpi_args': ''}
            plan = ci.make_plan([entry], ci.configuration(mode))
            common = root / ('common-' + mode)
            ci.dump(common / 'plan.json', plan)
            commons.append(common)
            group = plan['groups'][0]
            lanes, discoveries = [], []
            for ranks in map(int, entry['ranks'].split(';')):
                cases = ['Kernel/3']
                lane = {'entry': entry, 'ranks': ranks, 'cases': cases, 'configuration': mode, 'definition': group['definition']}
                lane['id'] = ci.archive.identity(lane)
                lane['label'] = mode
                lanes.append(lane)
                discoveries.append({'target': entry['target'], 'ranks': ranks, 'cases': cases})
                contribution = root / ('result-' + mode + str(ranks))
                report = contribution / 'report.json'
                ci.dump(report, self.report())
                values = ci.google_results(self.report(), cases, entry['metric'])
                ci.dump(contribution / 'receipt.json', {'schema': 1, 'commit': COMMIT, 'lane': lane,
                        'report': report.name, 'sha256': ci.archive.digest(report), 'results': values})
                contributions.append(contribution)
            built = root / ('group-' + mode)
            ci.dump(built / 'execution.json', {'schema': 1, 'commit': COMMIT, 'group': group,
                                             'lanes': lanes, 'discoveries': discoveries})
            groups.append(built)
        return commons, groups, contributions

    def test_exact_download_names_include_every_group_and_execution_lane(self):
        with patch.object(ci.archive, 'commit', return_value=COMMIT), tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            commons, groups, _ = self.fixture(root)
            for common in commons:
                plan = json.loads((common / 'plan.json').read_text())
                common.rename(root / ('benchmark-plan-common-' + plan['configuration']))
            expected_groups = []
            expected_results = []
            for group in groups:
                execution = json.loads((group / 'execution.json').read_text())
                name = 'benchmark-plan-group-' + execution['group']['id']
                expected_groups.append(name)
                expected_results += ['benchmark-result-' + lane['id'] for lane in execution['lanes']]
                group.rename(root / name)
            self.assertEqual(ci.artifact_names(root, 'groups'), sorted(expected_groups))
            self.assertEqual(ci.artifact_names(root, 'results'), sorted(expected_results))
            (root / expected_groups[0] / 'execution.json').unlink()
            with self.assertRaises(FileNotFoundError):
                ci.artifact_names(root, 'results')

    def test_complete_campaign_pools_without_mixing_ranks(self):
        with patch.object(ci.archive, 'commit', return_value=COMMIT), tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            inputs = self.fixture(root)
            ci.pool(*inputs, root / 'pooled')
            summary = json.loads((root / 'pooled/complete.json').read_text())
            self.assertEqual(summary['lanes'], 3)
            self.assertEqual(summary['cases'], 3)
            self.assertEqual(len(list((root / 'pooled').glob('**/results.json'))), 3)
            configurations = [json.loads(path.read_text()) for path in (root / 'pooled').glob('**/results.json')]
            self.assertEqual({row['ranks'] for row in configurations}, {1, 2, 4})
            self.assertTrue(all(row['results'][0]['value'] == 10 for row in configurations))

    def test_missing_configuration_group_or_receipt_fails(self):
        with patch.object(ci.archive, 'commit', return_value=COMMIT), tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            commons, groups, receipts = self.fixture(root)
            for modified in [(commons[:1], groups, receipts), (commons, groups[:1], receipts), (commons, groups, receipts[:1])]:
                with self.subTest(inputs=modified), self.assertRaises(ValueError):
                    ci.pool(*modified, root / 'pooled')
            self.assertFalse((root / 'pooled').exists())

    def test_omitting_a_lane_and_its_receipt_does_not_bypass_completeness(self):
        with patch.object(ci.archive, 'commit', return_value=COMMIT), tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            commons, groups, receipts = self.fixture(root)
            execution = json.loads((groups[1] / 'execution.json').read_text())
            execution['lanes'].pop()
            ci.dump(groups[1] / 'execution.json', execution)
            with self.assertRaisesRegex(ValueError, 'partitions'):
                ci.pool(commons, groups, receipts[:-1], root / 'pooled')

    def test_corrupt_raw_report_and_tampered_receipt_fail(self):
        for corrupt in ('report', 'receipt'):
            with self.subTest(corrupt=corrupt), patch.object(ci.archive, 'commit', return_value=COMMIT), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                inputs = self.fixture(root)
                path = inputs[2][0]
                if corrupt == 'report':
                    ci.dump(path / 'report.json', self.report(name='Other/3'))
                else:
                    receipt = json.loads((path / 'receipt.json').read_text())
                    receipt['results'][0]['value'] = 999
                    ci.dump(path / 'receipt.json', receipt)
                with self.assertRaises(ValueError):
                    ci.pool(*inputs, root / 'pooled')

    def test_stale_commit_and_duplicate_contributions_fail(self):
        with patch.object(ci.archive, 'commit', return_value=COMMIT), tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            commons, groups, receipts = self.fixture(root)
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                ci.pool(commons, groups, receipts + receipts[:1], root / 'pooled')
            plan = json.loads((commons[0] / 'plan.json').read_text())
            plan['commit'] = 'd' * 40
            ci.dump(commons[0] / 'plan.json', plan)
            with self.assertRaisesRegex(ValueError, 'stale'):
                ci.pool(commons, groups, receipts, root / 'pooled')

    def test_repetitions_use_median_and_ignore_aggregate_statistics(self):
        report = self.report()
        second = {**report['benchmarks'][0], 'cpu_time': 14, 'repetition_index': 1}
        report['benchmarks'] += [second, {'name': 'Kernel/3_cv', 'run_type': 'aggregate', 'time_unit': '%'}]
        result = ci.google_results(report, ['Kernel/3'], 'auto')[0]
        self.assertEqual(result['value'], 12)
        self.assertEqual(result['samples'], [10, 14])
        self.assertEqual(result['metric'], 'cpu_time')
        report['benchmarks'].append(second)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            ci.google_results(report, ['Kernel/3'], 'auto')

    def test_manual_and_wall_measurements_do_not_become_cpu_time(self):
        for suffix, metric in [('/manual_time', 'manual_time'), ('/real_time', 'real_time')]:
            with self.subTest(suffix=suffix):
                result = ci.google_results(self.report('Kernel/3' + suffix), ['Kernel/3' + suffix], 'auto')[0]
                self.assertEqual((result['value'], result['metric']), (30, metric))

    def test_numerical_errors_skips_missing_cases_and_nonfinite_times_fail(self):
        reports = [self.report(error_occurred=True), self.report(skipped=True), self.report()]
        reports[-1]['benchmarks'][0]['cpu_time'] = float('nan')
        for report in reports:
            with self.assertRaises(ValueError):
                ci.google_results(report, ['Kernel/3'], 'auto')
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            ci.google_results(self.report(), ['Kernel/3', 'Kernel/4'], 'auto')

    def test_missing_repetitions_fail_and_raw_counters_are_retained(self):
        report = self.report(repetitions=2, repetition_index=0, points=15)
        with self.assertRaisesRegex(ValueError, 'repetitions'):
            ci.google_results(report, ['Kernel/3'], 'auto')
        report['benchmarks'].append({**report['benchmarks'][0], 'repetition_index': 1})
        result=ci.google_results(report, ['Kernel/3'], 'auto')[0]
        self.assertEqual([row['points'] for row in result['records']], [15,15])

    def test_csv_workload_requires_all_operations(self):
        text = io.StringIO()
        writer = csv.DictWriter(text, fieldnames=['geometry','degree','resolution','distortion','warp_degree','query','pruning','samples','ns'])
        writer.writeheader()
        for query in ['interior', 'miss', 'build']:
            for pruning in ['0', '1']:
                writer.writerow(dict(geometry='Point', degree='1', resolution='1', distortion='0.000000', warp_degree='1', query=query, pruning=pruning, samples='1', ns='20'))
        self.assertEqual(len(ci.csv_results(text.getvalue(), ['Point/1/1/0/1'])), 6)
        with self.assertRaises(ValueError):
            ci.csv_results('\n'.join(text.getvalue().splitlines()[:-1]), ['Point/1/1/0/1'])

    def test_selector_partitions_cover_cases_and_respect_argument_budget(self):
        names = ['name/' + str(i) + 'a' * 400 for i in range(1000)]
        chunks = list(ci.partitions(names, 256))
        self.assertEqual([name for chunk in chunks for name in chunk], names)
        self.assertTrue(all(len(chunk) <= 256 and sum(len(name) + 1 for name in chunk) <= 60000 for chunk in chunks))


    def test_common_build_selects_only_the_required_library_closure(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);build=root / 'build'
            (root / 'unit.cpp').write_text('int value;')
            (root / 'main.cpp').write_text('int main() {}')
            (root / 'CMakeLists.txt').write_text('cmake_minimum_required(VERSION 3.16)\n'
                'project(Dependencies LANGUAGES CXX)\n'
                'add_library(Used STATIC unit.cpp)\nadd_library(Unused STATIC unit.cpp)\n'
                'add_library(Vendor STATIC unit.cpp)\n'
                'target_link_libraries(Used PRIVATE Vendor)\n'
                'add_executable(Benchmark main.cpp)\n'
                'target_link_libraries(Benchmark PRIVATE Used)\n')
            query=build / '.cmake/api/v1/query/codemodel-v2'
            query.parent.mkdir(parents=True);query.touch()
            subprocess.run(['cmake','-S',str(root),'-B',str(build)],check=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            self.assertEqual(ci.library_targets(build,['Benchmark']),['Used','Vendor'])

    def test_execution_discovers_exact_names_and_writes_a_validated_receipt(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(ci.archive, 'restore') as restore, patch.object(ci.archive, 'commit', return_value=COMMIT):
            root = Path(folder)
            driver = root / 'driver.py'
            driver.write_text("""import json,re,sys
names=['Kernel/3/repeats:2','Kernel/30','Other/3']
if '--benchmark_list_tests=true' in sys.argv:
 print('\\n'.join(names))
else:
 selector=next(arg.split('=',1)[1] for arg in sys.argv if arg.startswith('--benchmark_filter='))
 target=next(arg.split('=',1)[1] for arg in sys.argv if arg.startswith('--benchmark_out='))
 rows=[dict(name=name,iterations=1,cpu_time=2,real_time=3,time_unit='ns') for name in names if re.search(selector,name)]
 with open(target,'w') as stream:json.dump(dict(context={'driver':'fixture'},benchmarks=rows),stream)
""")
            entry = dict(target='Fixture', executable=sys.executable, args=str(driver), format='googlecpp', metric='auto')
            self.assertEqual(ci.discover(entry, 1), ['Kernel/3/repeats:2','Kernel/30','Other/3'])
            lane = dict(id='lane', definition={'kind':'benchmark'}, entry=entry, ranks=1, cases=['Kernel/3/repeats:2'])
            ci.dump(root / 'build/execution.json', {'lanes':[lane]})
            ci.run_lane(root / 'build', 'lane', root / 'output')
            restore.assert_called_once()
            receipt=json.loads((root / 'output/receipt.json').read_text())
            self.assertEqual([row['name'] for row in receipt['results']], ['Kernel/3/repeats:2'])
            self.assertEqual(receipt['sha256'], ci.archive.digest(root / 'output/report.json'))

    def test_source_mtimes_are_stable_and_dirty_checkouts_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            subprocess.run(['git','init','-q',str(root)],check=True)
            source=root / 'source.cpp';source.write_text('int value;')
            subprocess.run(['git','-C',str(root),'add','.'],check=True)
            subprocess.run(['git','-C',str(root),'-c','commit.gpgsign=false','-c','user.name=Fixture','-c','user.email=fixture@example.org','commit','-qm','fixture'],check=True)
            with patch.object(ci,'ROOT',root):
                ci.normalize_checkout()
                self.assertEqual(source.stat().st_mtime,946684800)
                source.write_text('int changed;')
                with self.assertRaisesRegex(ValueError,'clean committed'):
                    ci.normalize_checkout()


if __name__ == '__main__':
    unittest.main()

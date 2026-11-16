"""Resource summaries distinguish sampled pressure, kernel peaks and missing data."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import yaml

import ci_resource_monitor as monitor
import ci_summary

ROOT = Path(__file__).resolve().parents[1]


class ResourceStatsTest(unittest.TestCase):
    def sample(self, clock, total, idle, available, peak, oom, kills):
        return {'monotonic': clock, 'cpu': {'total': total, 'idle': idle, 'cores': 4},
                'host_memory': {'MemTotal': '100000 kB', 'MemAvailable': f'{available} kB'},
                'memory.peak': str(peak), 'memory.events': f'oom {oom}\noom_kill {kills}\n'}

    def test_counter_differences_and_peaks_have_explicit_scopes(self):
        result = monitor.statistics([self.sample(10, 100, 25, 50000, 100, 5, 2),
                                     self.sample(20, 300, 75, 10000, 200, 7, 3)])
        self.assertEqual(result['elapsed_seconds'], 10)
        self.assertEqual(result['cpu_percent'], 75)
        self.assertEqual(result['cpu_cores'], 4)
        self.assertEqual(result['host_total_bytes'], 100000 * 1024)
        self.assertEqual(result['peak_host_used_bytes'], 90000 * 1024)
        self.assertEqual(result['minimum_available_bytes'], 10000 * 1024)
        self.assertEqual(result['peak_cgroup_bytes'], 200)
        self.assertEqual(result['oom_events'], 2)
        self.assertEqual(result['oom_kills'], 1)

    def test_missing_samples_and_counter_resets_are_not_zero_usage(self):
        for samples in [[], [{'monotonic': 1, 'cpu': None, 'host_memory': None}],
                        [self.sample(1, 100, 25, 50000, 100, 5, 2),
                         self.sample(2, 50, 10, 50000, 100, 1, 0)]]:
            with self.subTest(samples=samples):
                result = monitor.statistics(samples)
                self.assertIsNone(result['cpu_percent'])
                self.assertIsNone(result['oom_events'])
                self.assertIsNone(result['oom_kills'])
        self.assertIsNone(monitor.statistics([])['peak_host_used_bytes'])

    def test_linux_snapshot_reads_cpu_without_double_counting_guest_time(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'proc').mkdir()
            (root / 'proc/stat').write_text('cpu 10 20 30 40 50 60 70 80 90 100\ncpu0 0\ncpu1 0\n')
            (root / 'proc/meminfo').write_text('MemTotal: 100 kB\nMemAvailable: 40 kB\n')
            group = root / 'sys/fs/cgroup'
            group.mkdir(parents=True)
            (group / 'memory.peak').write_text('12345\n')
            sample = monitor.snapshot(root)
        self.assertEqual(sample['cpu'], {'total': 360, 'idle': 90, 'cores': 2})
        self.assertEqual(sample['memory.peak'], '12345')
        self.assertIsNone(sample['memory.events'])

    def test_summary_labels_samples_and_unavailable_counters(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'summary.md'
            with patch.dict(os.environ, GITHUB_STEP_SUMMARY=str(path)):
                ci_summary.resources(monitor.statistics([]))
            text = path.read_text()
        self.assertIn('Unavailable', text)
        self.assertIn('Peak host memory used (sampled)', text)
        self.assertIn('100% means all CPUs busy', text)
        self.assertIn('before sampling started', text)
        self.assertNotIn('0.0%', text)
        self.assertNotIn('0.00 GiB', text)

    def test_start_finish_cli_saves_samples_and_preserves_status(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder) / 'resources'
            summary = Path(folder) / 'summary.md'
            env = dict(os.environ, GITHUB_STEP_SUMMARY=str(summary))
            command = [sys.executable, str(ROOT / 'dev/ci_resource_monitor.py')]
            try:
                subprocess.run([*command, 'start', str(directory)], env=env, check=True,
                               capture_output=True, timeout=10)
                time.sleep(0.3)
                duplicate = subprocess.run([*command, 'start', str(directory)], env=env,
                                           capture_output=True, timeout=10)
                self.assertNotEqual(duplicate.returncode, 0)
                subprocess.run([*command, 'finish', str(directory)], env=env, check=True,
                               capture_output=True, timeout=10)
                report = json.loads((directory / 'summary.json').read_text())
                self.assertGreaterEqual(report['samples'], 2)
                self.assertGreater(report['elapsed_seconds'], 0)
                self.assertIn('Runner resource usage', summary.read_text())
                time.sleep(0.4)
                finished = (directory / 'samples.jsonl').read_text()
                time.sleep(0.4)
                self.assertEqual((directory / 'samples.jsonl').read_text(), finished)
            finally:
                directory.mkdir(exist_ok=True)
                (directory / 'stop').touch()

    def test_build_test_and_benchmark_monitors_surround_environment_and_execution(self):
        selected = {
            'Build.yml': ['Kokkos'], 'TestCommon.yml': ['Build'],
            'TestDirectoryBuild.yml': ['Build'], 'TestGroup.yml': ['Run'],
            'BenchmarkConfiguration.yml': ['Common'], 'BenchmarkGroup.yml': ['Build', 'Run']}
        for filename, jobs in selected.items():
            data = yaml.safe_load((ROOT / '.github/workflows' / filename).read_text())
            for name in jobs:
                with self.subTest(workflow=filename, job=name):
                    steps = data['jobs'][name]['steps']
                    start = next(i for i,s in enumerate(steps) if s.get('name') == 'Start resource monitor')
                    end = next(i for i,s in enumerate(steps) if s.get('name') == 'Summarize resource usage')
                    self.assertEqual(start, 1)
                    self.assertEqual(end, len(steps) - 1)
                    self.assertEqual(steps[end]['if'], 'always()')
                    self.assertTrue(steps[start]['continue-on-error'])
                    self.assertTrue(steps[end]['continue-on-error'])
                    self.assertEqual(steps[start]['run'], 'python3 dev/ci_resource_monitor.py start')
                    self.assertEqual(steps[end]['run'], 'python3 dev/ci_resource_monitor.py finish')
        targeted = yaml.safe_load((ROOT / '.github/ci/targeted-lanes.yml').read_text())
        for job in targeted['jobs'].values():
            self.assertEqual(job['steps'][1]['name'], 'Start resource monitor')
            self.assertEqual(job['steps'][-1]['name'], 'Summarize resource usage')


if __name__ == '__main__':
    unittest.main()

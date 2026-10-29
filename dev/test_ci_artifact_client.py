"""Network retries must never accept an incomplete or unpinned client."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import ci_artifact_client as client
from ci_resource_monitor import snapshot


class ArtifactClientTest(unittest.TestCase):
    def fixture(self, root):
        source, destination = root / 'action', root / 'client'
        source.mkdir()
        (source / 'package-lock.json').write_text('{"pinned":true}\n')
        (source / 'package.json').write_text('{}\n')
        return source, destination

    def test_reset_retries_then_records_only_complete_install(self):
        with tempfile.TemporaryDirectory() as folder:
            source, destination = self.fixture(Path(folder))
            calls, waits = [], []
            def run(command, **kwargs):
                calls.append(command)
                self.assertFalse((destination / 'installed-package-lock.json').exists())
                if len(calls) == 1:
                    return subprocess.CompletedProcess(command, 152, '', 'npm error code ECONNRESET')
                package = destination / 'node_modules/@actions/artifact/package.json'
                package.parent.mkdir(parents=True);package.write_text('{}')
                return subprocess.CompletedProcess(command, 0, '', '')
            client.install(source, destination, run, waits.append)
            self.assertEqual(len(calls), 2)
            self.assertEqual(waits, [5])
            self.assertEqual((destination / 'installed-package-lock.json').read_bytes(), (source / 'package-lock.json').read_bytes())
            client.install(source, destination, lambda *a, **k: self.fail('Valid install should be reused'))

    def test_permanent_failures_do_not_retry_or_write_marker(self):
        for error in ('npm error code EINTEGRITY', 'npm error code EUSAGE', 'npm error code ERESOLVE ECONNRESET'):
            with self.subTest(error=error), tempfile.TemporaryDirectory() as folder:
                source, destination = self.fixture(Path(folder))
                calls = []
                def run(command, **kwargs):
                    calls.append(command)
                    return subprocess.CompletedProcess(command, 1, '', error)
                with self.assertRaises(subprocess.CalledProcessError):
                    client.install(source, destination, run, lambda _: self.fail('Unexpected sleep'))
                self.assertEqual(len(calls), 1)
                self.assertFalse((destination / 'installed-package-lock.json').exists())

    def test_network_failures_exhaust_exact_budget(self):
        with tempfile.TemporaryDirectory() as folder:
            source, destination = self.fixture(Path(folder));calls=[]
            def run(command, **kwargs):
                calls.append(command)
                return subprocess.CompletedProcess(command, 1, '', 'npm error network ETIMEDOUT')
            with self.assertRaises(subprocess.CalledProcessError):
                client.install(source, destination, run, lambda _: None)
            self.assertEqual(len(calls), 3)
            self.assertFalse((destination / 'installed-package-lock.json').exists())

    def test_false_success_does_not_record_install(self):
        with tempfile.TemporaryDirectory() as folder:
            source, destination = self.fixture(Path(folder))
            with self.assertRaisesRegex(RuntimeError, 'without installing'):
                client.install(source, destination, lambda command, **k: subprocess.CompletedProcess(command, 0, '', ''))
            self.assertFalse((destination / 'installed-package-lock.json').exists())

    def test_memory_snapshot_handles_missing_cgroups_without_inventing_values(self):
        with tempfile.TemporaryDirectory() as folder:
            result = snapshot(Path(folder))
            self.assertIsNone(result['memory.current'])
            self.assertIsNone(result['host_memory'])
            group = Path(folder) / 'sys/fs/cgroup';group.mkdir(parents=True)
            (group / 'memory.events').write_text('oom 2\noom_kill 1\n')
            self.assertEqual(snapshot(Path(folder))['memory.events'], 'oom 2\noom_kill 1')


if __name__ == '__main__':
    unittest.main()

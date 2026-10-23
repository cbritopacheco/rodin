"""Complete publication preserves independent lane evidence and zero-hit lines."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import ci_coverage as coverage
import ci_coverage_publish as publish


class PublicationTest(unittest.TestCase):
    def test_configuration_local_collection_requires_only_its_complete_partition_set(self):
        import ci_build_graph as graph
        import ci_builds as archive
        definitions = [{'job': 'Convergence', 'matrix': {'multithreaded': mt, 'suite': 'curved_projection'}}
                       for mt in ['ON', 'OFF']]
        common = archive.identity(graph.common_specification(graph.specification(**definitions[0])))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lane = root / 'lane'
            lane.mkdir()
            report = lane / 'RodinCoverage.info'
            report.write_text('SF:src/Rodin/Foo.h\nDA:1,1\nDA:2,0\nend_of_record\n')
            receipt = {'definition': definitions[0], 'commit': 'current',
                       'report_sha256': hashlib.sha256(report.read_bytes()).hexdigest()}
            (lane / 'target-receipt.json').write_text(json.dumps(receipt))
            publish.collect('targeted', 'Convergence', lane, root / 'out', 'current', {}, definitions, common)
            batch = json.loads((root / 'out/publication.json').read_text())
            self.assertEqual(batch['common'], common)
            self.assertEqual(set(batch['contributions']), {'target:' + coverage.target_identity(definitions[0])})
            receipt['definition'] = definitions[1]
            (lane / 'target-receipt.json').write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, 'Unexpected'):
                publish.collect('targeted', 'Convergence', lane, root / 'invalid', 'current', {}, definitions, common)

    def test_partition_artifact_names_do_not_depend_on_matrix_positions(self):
        import os
        definitions = [{'job': 'Convergence', 'matrix': {'multithreaded': 'ON', 'suite': suite}}
                       for suite in ['pressure_light', 'pressure_tetrahedron']]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build = root / 'build'
            build.mkdir()
            (build / 'RodinCoverage.info').write_text('SF:src/Rodin/Foo.h\nDA:1,1\nend_of_record\n')
            (build / 'coverage-audit.json').write_text('{}')
            output = root / 'outputs'
            for index, definition in enumerate(definitions):
                with patch.dict(os.environ, CI_TARGET_DEFINITION=json.dumps(definition),
                                GITHUB_SHA='current', GITHUB_OUTPUT=str(output)):
                    coverage.record_target(SimpleNamespace(build=build, output=root / str(index)))
            names = [line.split('=', 1)[1] for line in output.read_text().splitlines()]
            self.assertEqual(len(set(names)), 2)
            for name, definition in zip(names, definitions):
                self.assertTrue(name.endswith(coverage.target_identity(definition)[:24]))

    def fixture(self, root):
        ordinary = {'one': {'coverage_expected': True, 'report_sha256': ''}}
        definition = {'job': 'Convergence', 'matrix': {'multithreaded': 'ON'}}
        key = coverage.target_identity(definition)
        targeted = {key: {'definition': definition, 'report_sha256': ''}}
        for name, records, text in [
                ('one', ordinary, 'DA:1,1\nDA:2,0\nDA:3,0\n'),
                (key, targeted, 'DA:1,0\nDA:2,1\nDA:3,0\n')]:
            path = root / name
            path.mkdir()
            report = path / 'RodinCoverage.info'
            report.write_text('SF:src/Rodin/Math/Operator.h\n' + text + 'end_of_record\n')
            records[name]['report_sha256'] = hashlib.sha256(report.read_bytes()).hexdigest()
        batch = root / 'publication'
        batch.mkdir()
        lines = coverage.production_lines(list(root.rglob('RodinCoverage.info')))
        report = batch / 'CombinedCoverage.info'
        publish.write_lcov(lines, report)
        receipt = {'commit': 'current', 'uploaded': True,
                   'contributions': {'ordinary:one': ordinary['one']['report_sha256'],
                                     'target:' + key: targeted[key]['report_sha256']},
                   'report_sha256': hashlib.sha256(report.read_bytes()).hexdigest()}
        (batch / 'publication.json').write_text(json.dumps(receipt))
        return ordinary, targeted, batch, receipt

    def test_all_hits_and_zero_hit_lines_survive_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ordinary, targeted, batch, _ = self.fixture(root)
            publish.verify_publications(ordinary, targeted, root, 'current')
            lines = coverage.production_lines([batch / 'CombinedCoverage.info'])
            self.assertEqual(len(lines), 3)
            self.assertEqual(sum(lines.values()), 2)

    def test_missing_failed_stale_changed_and_duplicate_publications_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ordinary, targeted, batch, receipt = self.fixture(root)
            file = batch / 'publication.json'
            for change in [{'uploaded': False}, {'commit': 'old'},
                           {'contributions': {}}, {'contributions': {'ordinary:one': 'wrong'}}]:
                file.write_text(json.dumps(receipt | change))
                with self.assertRaises(ValueError):
                    publish.verify_publications(ordinary, targeted, root, 'current')
            file.write_text(json.dumps(receipt))
            duplicate = root / 'duplicate'
            duplicate.mkdir()
            (duplicate / 'publication.json').write_text(json.dumps(receipt))
            (duplicate / 'CombinedCoverage.info').write_bytes((batch / 'CombinedCoverage.info').read_bytes())
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                publish.verify_publications(ordinary, targeted, root, 'current')

    def test_successful_upload_cannot_conceal_a_different_union(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ordinary, targeted, batch, receipt = self.fixture(root)
            report = batch / 'CombinedCoverage.info'
            report.write_text(report.read_text().replace('DA:2,1', 'DA:2,0'))
            receipt['report_sha256'] = hashlib.sha256(report.read_bytes()).hexdigest()
            (batch / 'publication.json').write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, 'complete raw lane union'):
                publish.verify_publications(ordinary, targeted, root, 'current')

    def test_targeted_collection_requires_every_planned_partition(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lane = root / 'lane'
            lane.mkdir()
            definitions = [{'job': 'Convergence', 'matrix': {'np': n}} for n in [1, 2]]
            report = lane / 'RodinCoverage.info'
            report.write_text('SF:src/Rodin/Foo.h\nDA:1,0\nend_of_record\n')
            receipt = {'definition': definitions[0], 'commit': 'current',
                       'report_sha256': hashlib.sha256(report.read_bytes()).hexdigest()}
            (lane / 'target-receipt.json').write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, 'Missing targeted'):
                publish.collect('targeted', 'Convergence', lane, root / 'out', 'current', {}, definitions)


if __name__ == '__main__':
    unittest.main()

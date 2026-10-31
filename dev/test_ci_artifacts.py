"""Artifact lookup must preserve every planned build and coverage contribution."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import yaml

import ci_artifacts as artifacts
import ci_build_graph as graph
import ci_coverage as coverage
import ci_coverage_publish as publish
import ci_schedule as schedule
import ci_tests as tests


class ArtifactPlanTest(unittest.TestCase):
    def test_publication_names_match_recorded_target_artifacts_and_partition(self):
        definitions = [{'job': 'Convergence', 'matrix': {'multithreaded': mt, 'suite': 'curved_projection'}}
                       for mt in ['ON', 'OFF']]
        common = artifacts.common_id(**definitions[0])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build = root / 'build'
            build.mkdir()
            (build / 'RodinCoverage.info').write_text('SF:src/Rodin/Foo.h\nDA:1,1\nend_of_record\n')
            (build / 'coverage-audit.json').write_text('{}')
            output = root / 'outputs'
            from types import SimpleNamespace
            import os
            with patch.dict(os.environ, CI_TARGET_DEFINITION=json.dumps(definitions[0]),
                            GITHUB_SHA='current', GITHUB_OUTPUT=str(output)):
                coverage.record_target(SimpleNamespace(build=build, output=root / 'record'))
            recorded = output.read_text().strip().removeprefix('artifact=')
        self.assertEqual(artifacts.publication_names('targeted', 'Convergence', common, {}, definitions), [recorded])
        self.assertNotEqual(recorded, artifacts.target_name(definitions[1]))

    def test_every_configuration_contribution_appears_in_exactly_one_publisher_and_finalization(self):
        ordinary, targeted = schedule.requests()
        common = graph.prepare(ordinary, targeted)
        families = {row['id']: family for family, rows in ordinary.items() for row in rows}
        published = []
        expected_publications = []
        for entry in common:
            if entry['definition']['cmake'].get('RODIN_CODE_COVERAGE') != 'ON':
                continue
            selected = {families[r['matrix']['id']] if r['job'] == 'Ordinary' else r['job']
                        for r in entry['requests']}
            for family in selected:
                kind = 'ordinary' if family in ordinary else 'targeted'
                names = artifacts.publication_names(kind, family, entry['id'], ordinary, targeted)
                self.assertTrue(names)
                published += names
                expected_publications.append('ci-publication-' + family + '-' + entry['id'])
                part, targets = publish.publication_plan(kind, family, ordinary, targeted, entry['id'])
                self.assertEqual(len(names), sum(map(len, part.values())) + len(targets))
        self.assertEqual(len(published), len(set(published)))
        final = artifacts.final_names(ordinary, targeted, common, {'builds': {'build-a': {}, 'build-b': {}}})
        self.assertEqual(len(final), len(set(final)))
        self.assertEqual(set(final), set(artifacts.contribution_names(ordinary, targeted) +
                                        expected_publications + ['ci-built-build-a', 'ci-built-build-b']))
        self.assertLessEqual(set(published), set(final))
        # Sanitized receipts are still mandatory final evidence, despite no upload.
        self.assertGreater(len(artifacts.contribution_names(ordinary, targeted)), len(published))

    def test_ordinary_names_match_the_actual_upload_step(self):
        ordinary, _ = schedule.requests()
        workflow = yaml.safe_load((tests.ROOT / '.github/workflows/TestGroup.yml').read_text())
        upload = next(s for s in workflow['jobs']['Run']['steps']
                      if s.get('name') == 'Save invocation and coverage records')
        for family, rows in ordinary.items():
            for row in rows:
                expected = upload['with']['name'].replace('${{ inputs.family }}', family).replace(
                    '${{ steps.selection.outputs.common }}', artifacts.common_id('Ordinary', row)).replace(
                    '${{ inputs.id }}', row['id'])
                self.assertEqual(artifacts.ordinary_name(family, row), expected)

    def test_more_than_1000_partitions_are_not_truncated(self):
        definitions = [{'job': 'Convergence', 'matrix': {'multithreaded': 'ON',
                        'suite': 'curved_projection', 'partition': i}} for i in range(1005)]
        names = artifacts.publication_names('targeted', 'Convergence', artifacts.common_id(**definitions[0]),
                                            {}, definitions)
        self.assertEqual(len(names), 1005)
        self.assertIn(artifacts.target_name(definitions[-1]), names)

    def test_workflows_use_exact_names_and_retain_complete_coverage_gates(self):
        root = tests.ROOT
        documents = [(path, yaml.safe_load(path.read_text())) for path in (root / '.github/workflows').glob('*.yml')]
        documents += [(path, yaml.safe_load(path.read_text())) for path in (root / '.github/actions').glob('*/action.yml')]
        for path, document in documents:
            jobs = document.get('jobs', {})
            groups = [job.get('steps', []) for job in jobs.values()] + [document.get('runs', {}).get('steps', [])]
            for steps in groups:
                for step in steps:
                    if step.get('uses', '').startswith('actions/download-artifact@'):
                        self.assertIn('name', step['with'], str(path))
                        self.assertNotIn('pattern', step['with'], str(path))
        finalize = yaml.safe_load((root / '.github/workflows/Tests.yml').read_text())['jobs']['CoverageFinalize']
        commands = [step.get('run', '') for step in finalize['steps']]
        self.assertIn('python3 dev/ci_coverage.py finalize', commands)
        self.assertIn('python3 dev/ci_coverage.py threshold', commands)
        self.assertIn('python3 dev/ci_artifacts.py final --output ci-final-inputs.json', commands)

    def test_benchmark_restore_uses_retrying_exact_names_and_preserves_uploads(self):
        workflow = yaml.safe_load((tests.ROOT / '.github/workflows/BenchmarkGroup.yml').read_text())
        for job, artifact, output in [('Build', 'benchmark-common-${{ inputs.configuration }}', 'benchmark-group'),
                                      ('Run', 'benchmark-group-${{ inputs.group }}', 'benchmark-result')]:
            steps = workflow['jobs'][job]['steps']
            self.assertFalse(any(s.get('uses', '').startswith('actions/download-artifact@') for s in steps))
            restore = next(s for s in steps if s.get('uses') == './.github/actions/download-artifacts')
            self.assertEqual(json.loads(restore['with']['names']), [artifact])
            self.assertEqual(restore['with']['path'], 'benchmark-inputs')
            execute = next(s for s in steps if '--source "$SOURCE"' in s.get('run', ''))
            self.assertEqual(execute['env']['SOURCE'], 'benchmark-inputs/' + artifact)
            upload = next(s for s in steps if s.get('uses') == 'actions/upload-artifact@v4')
            self.assertEqual(upload['with']['path'], output)
            self.assertIn('name', upload['with'])

    def test_runtime_downloads_and_retries(self):
        subprocess.run(['node', '--test', '--test-concurrency=1', 'dev/test_ci_artifact_downloads.js'],
                       cwd=tests.ROOT, check=True)


if __name__ == '__main__':
    unittest.main()

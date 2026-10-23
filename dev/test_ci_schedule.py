"""Hierarchical readiness must not reintroduce backend or configuration barriers."""
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import yaml

import ci_build_graph as graph
import ci_lanes
import ci_schedule as schedule
import ci_tests as tests


class HierarchicalScheduleTest(unittest.TestCase):
    def test_variants_do_not_declare_unrelated_execution_or_family_branches(self):
        ordinary, targeted = schedule.requests()
        entries = graph.prepare(ordinary, targeted)
        families = {r['id']: f for f, rows in ordinary.items() for r in rows}
        expected = {}
        for entry in entries:
            variant = graph.configuration_variant(entry['definition'])
            expected.setdefault(variant, set()).update(schedule.execution_key(
                r, families[r['matrix']['id']] if r['job'] == 'Ordinary' else r['job']) for r in entry['requests'])
        declared = 0
        for variant, identities in expected.items():
            _, executions = schedule.layout(variant)
            self.assertEqual(set(executions), identities)
            jobs = schedule.configuration_workflow(variant)['jobs']
            declared += len(jobs) * sum(graph.configuration_variant(e['definition']) == variant for e in entries)
            if variant.endswith('_sanitized'):
                self.assertFalse(any(name.startswith('Publish_') for name in jobs))
        self.assertLess(declared, 46 * 235 // 4)
        self.assertEqual(len(schedule.configuration_workflow('mpi')['jobs']), 4)
        for family in ci_lanes.jobs():
            jobs = schedule.targeted_workflow(family)['jobs']
            self.assertEqual(set(jobs), {family})
            self.assertNotIn('if', jobs[family])

    def test_preloaded_artifacts_never_fall_back_to_network(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'ci-build-one'
            source.mkdir()
            (source / 'manifest.json').write_text('{}')
            with patch.dict(os.environ, CI_RESTORED_ARTIFACTS=str(root)), patch.object(graph.subprocess, 'run') as run:
                self.assertEqual(graph.download('ci-build-one', root / 'unused'), source)
                with self.assertRaisesRegex(ValueError, 'was not downloaded'):
                    graph.download('ci-build-missing', root / 'unused')
                run.assert_not_called()

    def test_only_required_hierarchies_block_each_execution(self):
        producers, executions = schedule.layout()
        workflow = schedule.configuration_workflow()['jobs']
        top = yaml.safe_load((tests.ROOT / '.github/workflows/Tests.yml').read_text())['jobs']
        configurations = [job for name, job in top.items() if name.startswith('Configurations_')]
        self.assertEqual(len(configurations), 11)
        for job in configurations:
            self.assertEqual(job['needs'], 'Plan')
            self.assertFalse(job['strategy']['fail-fast'])
        for identity in producers:
            self.assertEqual(workflow[identity]['needs'], ['Common'])
            self.assertFalse(workflow[identity]['strategy']['fail-fast'])
        for identity, execution in executions.items():
            self.assertEqual(set(workflow[identity]['needs']), {'Common', *execution['dependencies']})
            self.assertTrue(set(execution['dependencies']) <= producers.keys())
            if execution['job'] == 'ConvergenceBaseline':
                self.assertEqual(len(execution['dependencies']), 1)
        for family in ['Convergence', 'PETScCurvedConvergence', 'PETScCurvedBoundaryConvergence']:
            members = [e for e in executions.values() if e['job'] == family]
            self.assertGreater(len(members), 1)
            # Geometry/physics branches do not wait for every other branch.
            self.assertTrue(any(set(a['dependencies']).isdisjoint(b['dependencies'])
                                for a in members for b in members))

    def test_generated_runtime_preserves_every_numerical_step(self):
        generated = schedule.targeted_workflow()['jobs']
        for name, source in ci_lanes.jobs().items():
            runtime = generated[name]
            for field in ['env', 'timeout-minutes', 'runs-on', 'name']:
                self.assertEqual(runtime.get(field), source.get(field))
            self.assertEqual([s for s in runtime['steps'] if s.get('id') != 'coverage_record'],
                             [s for s in source['steps'] if s.get('id') != 'coverage_record'])
            self.assertEqual(runtime['if'], "${{ inputs.family == '" + name + "' }}")
            self.assertEqual(runtime['strategy']['matrix'], '${{ fromJSON(inputs.matrix) }}')
            self.assertFalse(runtime['strategy']['fail-fast'])

    def test_configured_binding_rejects_missing_and_incompatible_source_edges(self):
        ordinary, _ = schedule.requests()
        family, row = next((f, r) for f, rows in ordinary.items() for r in rows
                           if r['group'] == 'unit/Rodin/Geometry')
        request = {'job': 'Ordinary', 'matrix': row}
        entry = {'id': 'common', 'definition': {'mode': 'default', 'cmake': {'RODIN_CODE_COVERAGE': 'ON'}}, 'requests': [request]}
        producer = schedule.key('b', row['group'])
        execution = schedule.execution_key(request, family)
        index = {'builds': {'one': {'definition': {'directory': row['group'],
                   'targets': ['GeometryFixture'], 'runner': 'ubuntu-latest'}}},
                 'requests': {graph.request_key(**request): {'builds': ['one'], 'selected': ['fixture']}}}
        good = ({producer: {'path': row['group']}}, {execution: {'family': family, 'dependencies': [producer]}})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'ci-common-plan.json').write_text(json.dumps([entry]))
            output = root / 'outputs'
            with patch.object(schedule, 'ROOT', root), patch.object(graph, 'index_graph', return_value=(index, {})), \
                    patch.dict(os.environ, GITHUB_OUTPUT=str(output)), patch.object(schedule, 'layout', return_value=good):
                schedule.prepare('common')
            values = dict(line.split('=', 1) for line in output.read_text().splitlines())
            self.assertEqual(json.loads(values[producer])['include'][0]['id'], 'one')
            self.assertEqual(json.loads(values[execution])['include'], [dict(row, family=family)])
            self.assertEqual(values['p_' + family], 'true')
            broken = copy.deepcopy(good)
            broken[1][execution]['dependencies'] = []
            with patch.object(schedule, 'ROOT', root), patch.object(graph, 'index_graph', return_value=(index, {})), \
                    patch.object(schedule, 'layout', return_value=broken):
                with self.assertRaisesRegex(ValueError, 'Missing execution dependency edge'):
                    schedule.prepare('common')
            with patch.object(schedule, 'ROOT', root), patch.object(graph, 'index_graph', return_value=(index, {})), \
                    patch.object(schedule, 'layout', return_value=({}, good[1])):
                with self.assertRaisesRegex(ValueError, 'no unique declared hierarchy edge'):
                    schedule.prepare('common')

    def test_generated_workflows_are_current_and_bounded(self):
        for name, document in schedule.generated().items():
            text = (tests.ROOT / '.github/workflows' / name).read_text()
            self.assertEqual(yaml.safe_load(text), document)
            self.assertLess(len(text.encode()), 500_000)
            for job in document['jobs'].values():
                self.assertTrue(set(job.get('needs', [])) <= document['jobs'].keys())


if __name__ == '__main__':
    unittest.main()

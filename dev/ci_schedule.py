#!/usr/bin/env python3
"""Generate explicit hierarchical dependency edges, then bind configured builds.

GitHub cannot express a dependency on one matrix child. Each configuration
therefore owns a reusable DAG with distinct producer jobs per registration
directory/component. Source directories bound that DAG; configured manifests
alone choose executable targets. Missing edges fail before any directory builds.
"""
from __future__ import annotations

import argparse
import copy
import functools
import html
import json
import os
from pathlib import Path

import yaml

import ci_build_graph as graph
import ci_builds as archive
import ci_lanes
import ci_test_plan as planning
import ci_tests as tests

ROOT = tests.ROOT
# Source scheduling bounds for optional directories. Configured manifests still
# certify every selected target and reject an omitted producer.
DIRECTORY_REQUIREMENTS = {'unit/Rodin/Scotch': 'RODIN_USE_SCOTCH',
                          'unit/Rodin/Kokkos': 'RODIN_USE_KOKKOS'}


def key(prefix, value):
    return prefix + '_' + archive.identity(value)[:16]


def scope(target):
    owner = target['owner'].removeprefix('tests/')
    component = tests.COMPONENTS.get(target['target'])
    if target['target'] == 'RodinMPIMatrixRangeTest':
        component = 'MatrixRanges'
    return owner + '/' + component if component else owner


def source_scopes():
    """Declare directory edges, never executable eligibility or build lists."""
    scopes, declared = {}, {}
    for suite in ['unit', 'manufactured', 'convergence']:
        for path in sorted((ROOT / 'tests' / suite).rglob('CMakeLists.txt')):
            owner = str(path.parent.relative_to(ROOT / 'tests'))
            commands = list(tests.commands(path.read_text()))
            # Only executable declarations create build producers. add_test and
            # target_link_libraries references are not additional executable owners.
            declarations = {args[0] for name, args, _ in commands if args and
                            name.lower() in ['add_executable', 'rodin_add_petsc_refinement'] and
                            args[0].startswith('Rodin')}
            if any(name.lower() == 'add_executable' and args and not args[0].startswith('Rodin')
                   for name, args, _ in commands):
                # Loop-generated executable names retain their physical owner;
                # their configured target inventory remains authoritative.
                scopes[owner] = {'owner': owner, 'execution': tests.execution_group(owner)}
            for name in sorted(declarations):
                target = {'owner': 'tests/' + owner, 'target': name}
                branch = scope(target)
                if name in declared and declared[name] != branch:
                    raise ValueError('Ambiguous source dependency edge: ' + name)
                declared[name] = branch
                scopes[branch] = {'owner': owner, 'execution': tests.execution_group(branch)}
    # These two target names are assembled in a CMake loop in Diffusion.
    # This supplies scheduling edges only; prepare() verifies actual ownership.
    for physics in ['Poisson', 'Conductivity']:
        declared['RodinConvergenceIsoparametricPETSc' + physics + 'Boundary'] = 'convergence/isoparametric/Diffusion'
    return scopes, declared


@functools.lru_cache(maxsize=1)
def requests():
    ordinary = planning.plan()
    native = [{'multithreaded': mt, 'suite': 'baseline', 'group': group,
               'path': tests.display_group(group)} for group in planning.native_convergence_groups()
              for mt in ['ON', 'OFF']]
    targeted = planning.targeted_plan(native)
    return ordinary, targeted


def execution_key(request, family):
    job, matrix = request['job'], request['matrix']
    selection = (matrix['group'] if job in ['Ordinary', 'ConvergenceBaseline'] else
                 graph.specification(job, matrix)['targets'])
    return key('e', [job, family, selection])


def layout(variant=None):
    scopes, declared = source_scopes()
    ordinary, targeted = requests()
    executions = {}
    rows = [(family, {'job': 'Ordinary', 'matrix': row})
            for family, entries in ordinary.items() for row in entries]
    rows += [(r['job'], r) for r in targeted]
    for family, request in rows:
        job, matrix = request['job'], request['matrix']
        spec = graph.specification(job, matrix)
        if variant is not None and graph.configuration_variant(spec) != variant:
            continue
        identity = execution_key(request, family)
        if identity in executions:
            continue
        if job == 'Ordinary':
            group = matrix['group']
            dependencies = {s for s in scopes if spec['cmake'].get('RODIN_BUILD_' +
                            s.split('/')[0].upper() + '_TESTS') == 'ON'} if group == 'src/Rodin' else {
                s for s, info in scopes.items() if info['execution'] == group or s == group}
            label = matrix['name']
        elif job == 'ConvergenceBaseline':
            dependencies = {matrix['group']}
            label = matrix['path']
        else:
            required = graph.specification(job, matrix)['targets']
            missing = set(required) - declared.keys()
            if missing:
                raise ValueError('Dedicated targets need explicit source dependency edges: ' + str(missing))
            dependencies = {declared[t] for t in required}
            label = ' + '.join(tests.display_group(s) for s in sorted(dependencies))
        dependencies = {s for s in dependencies if all(
                        spec['cmake'].get(option) == 'ON'
                        for directory, option in DIRECTORY_REQUIREMENTS.items()
                        if s == directory or s.startswith(directory + '/'))}
        executions[identity] = {'job': job, 'family': family, 'name': label,
                                'dependencies': sorted({key('b', graph.build_group(s)) for s in dependencies})}
    required = {k for e in executions.values() for k in e['dependencies']}
    return {key('b', group): {'path': group} for group in sorted({graph.build_group(s) for s in scopes})
            if variant is None or key('b', group) in required}, executions


def prepare(common_id):
    """Validate every actual dependency before emitting the configuration DAG."""
    entries = json.loads((ROOT / 'ci-common-plan.json').read_text())
    entry = next(e for e in entries if e['id'] == common_id)
    index, _ = graph.index_graph(ROOT / 'ci-common-record', common_id=common_id)
    producers, executions = layout(graph.configuration_variant(entry['definition']))
    build_rows, execution_rows = {}, {}
    build_owner = {}
    for identity, record in index['builds'].items():
        spec = record['definition']
        owner = key('b', graph.build_group(spec['directory']))
        if owner not in producers or spec['directory'] != graph.build_group(spec['directory']):
            raise ValueError('Configured producer has no unique declared hierarchy edge')
        build_owner[identity] = owner
        label = tests.display_group(producers[owner]['path'])
        if owner in build_rows:
            raise ValueError('Hierarchy group has more than one build producer')
        build_rows.setdefault(owner, []).append({'id': identity, 'runner': spec['runner'],
                                                 'name': label})
    ordinary, _ = requests()
    families = {row['id']: family for family, rows in ordinary.items() for row in rows}
    for request in entry['requests']:
        family = families[request['matrix']['id']] if request['job'] == 'Ordinary' else request['job']
        identity = execution_key(request, family)
        actual = {build_owner[k] for k in index['requests'][graph.request_key(**request)]['builds']}
        if not actual <= set(executions[identity]['dependencies']):
            raise ValueError('Missing execution dependency edge: ' + identity + ': ' + str(actual))
        row = dict(request['matrix'])
        if request['job'] == 'Ordinary':
            row['family'] = family
        execution_rows.setdefault(identity, []).append(row)
    expected = {graph.request_key(**r) for r in entry['requests']}
    if expected != index['requests'].keys():
        raise ValueError('Configuration scheduling lost execution requests')
    if build_rows.keys() != producers.keys():
        missing = [producers[k]['path'] for k in producers.keys() - build_rows.keys()]
        raise ValueError('Declared producer has no configured executable: ' + str(sorted(missing)))
    if execution_rows.keys() != executions.keys():
        raise ValueError('Declared execution group has no required requests')
    (ROOT / 'ci-build-index.json').write_text(json.dumps(index, indent=2) + '\n')
    values = {**build_rows, **execution_rows}
    with open(os.environ['GITHUB_OUTPUT'], 'a') as output:
        for name, rows in values.items():
            if not 1 <= len(rows) <= 256:
                raise ValueError('Hierarchical job matrix exceeds 256 entries: ' + name)
            output.write(name + '=' + json.dumps({'include': rows}, separators=(',', ':')) + '\n')
        flags = entry['definition']['cmake']
        if flags.get('RODIN_CODE_COVERAGE') == 'ON':
            for family in sorted({executions[k]['family'] for k in execution_rows}):
                contributing = [r for r in entry['requests'] if
                                (families[r['matrix']['id']] if r['job'] == 'Ordinary' else r['job']) == family and
                                index['requests'][graph.request_key(**r)]['selected']]
                if not contributing:
                    raise ValueError('Declared coverage family has no contributing requests: ' + family)
                output.write('p_' + family + '=true\n')
    print('Independent hierarchical producers:', len(build_rows), '; execution groups:', len(execution_rows))


def expression(value):
    return '${{ ' + value + ' }}'


def configuration_workflow(variant=None):
    producers, executions = layout(variant)
    common_source = yaml.safe_load((ROOT / '.github/workflows/TestCommon.yml').read_text())
    common = copy.deepcopy(common_source['jobs']['Build'])
    common['name'] = 'Build / Common Libraries'
    common['outputs'] = {name: expression('steps.schedule.outputs.' + name) for name in [
        *producers, *executions, *('p_' + family for family in sorted({e['family'] for e in executions.values()}))]}
    common['steps'] += [
        {'name': 'Bind configured targets to independent hierarchy edges', 'id': 'schedule',
         'env': {'COMMON_ID': expression('inputs.id')},
         'run': 'python3 dev/ci_schedule.py --prepare "$COMMON_ID"'},
        {'uses': 'actions/upload-artifact@v4', 'with': {'overwrite': True,
         'name': 'ci-build-index-${{ inputs.id }}', 'path': 'ci-build-index.json', 'if-no-files-found': 'error'}}]
    jobs = {'Common': common}
    for identity, producer in sorted(producers.items(), key=lambda x: x[1]['path']):
        jobs[identity] = {'name': 'Build / ${{ matrix.name }}', 'needs': ['Common'],
                         'if': expression('!cancelled()'),
                         'strategy': {'fail-fast': False, 'matrix': expression('fromJSON(needs.Common.outputs.' + identity + ')')},
                         'uses': './.github/workflows/TestDirectoryBuild.yml',
                         'with': {'id': expression('matrix.id'), 'runner': expression('matrix.runner'),
                                  'common': expression('inputs.id')}}
    for identity, execution in sorted(executions.items()):
        dependencies = ['Common', *execution['dependencies']]
        # Every declaration is mandatory. Missing matrices or artifacts must
        # fail; dependency failures must not turn required work into skipped jobs.
        job = {'name': execution['name'], 'needs': dependencies, 'if': expression('!cancelled()')}
        if execution['job'] == 'Ordinary':
            job.update({'strategy': {'fail-fast': False, 'matrix': expression('fromJSON(needs.Common.outputs.' + identity + ')')},
                        'uses': './.github/workflows/TestGroup.yml',
                        'with': {k: expression('matrix.' + k) for k in ['family', 'backend', 'profile', 'group',
                                 'build_type', 'multithreaded', 'asan', 'ubsan', 'id', 'name', 'baseline']}})
        else:
            job.update({'uses': './.github/workflows/Test' + execution['job'] + '.yml',
                        'with': {'matrix': expression('needs.Common.outputs.' + identity)}})
        job['secrets'] = 'inherit'
        jobs[identity] = job
    for family in sorted({e['family'] for e in executions.values()}):
        if variant is not None and variant.endswith('_sanitized'):
            continue
        members = [k for k, e in executions.items() if e['family'] == family]
        jobs['Publish_' + family] = {
            'name': 'Coverage / Publish / ' + tests.display_group(family.replace('_', '/')),
            'needs': ['Common', *members],
            'if': expression('!cancelled()'),
            'uses': './.github/workflows/CoveragePublish.yml', 'secrets': 'inherit',
            'with': {'kind': 'ordinary' if family in requests()[0] else 'targeted', 'family': family,
                     'common': expression('inputs.id')}}
    return {'name': 'Tests / Configuration', 'on': {'workflow_call': {'inputs': {
        'definition': {'type': 'string', 'required': True}, 'id': {'type': 'string', 'required': True}},
        'secrets': {'CODECOV_TOKEN': {'required': False}}}},
        'permissions': {'contents': 'read', 'actions': 'read'}, 'jobs': jobs}


def targeted_workflow(family=None):
    jobs = {k: copy.deepcopy(v) for k, v in ci_lanes.jobs().items() if family is None or family == k}
    for name, job in jobs.items():
        job.pop('needs', None)
        if family is None:
            job['if'] = expression("inputs.family == '" + name + "'")
        job['strategy']['matrix'] = expression('fromJSON(inputs.matrix)')
        for step in job['steps']:
            if step.get('id') == 'coverage_record':
                step['with'].pop('artifact_name', None)
    inputs = {'matrix': {'type': 'string', 'required': True}}
    if family is None:
        inputs['family'] = {'type': 'string', 'required': True}
    return {'name': 'Tests / Dedicated Selection', 'on': {'workflow_call': {'inputs': inputs,
        'secrets': {'CODECOV_TOKEN': {'required': False}}}},
        'permissions': {'contents': 'read', 'actions': 'read'}, 'jobs': jobs}


def generated():
    ordinary, targeted = requests()
    variants = sorted({graph.configuration_variant(e['definition']) for e in graph.prepare(ordinary, targeted)})
    return {**{'TestConfiguration_' + v + '.yml': configuration_workflow(v) for v in variants},
            **{'Test' + f + '.yml': targeted_workflow(f) for f in ci_lanes.jobs()}}


def workflow_diagram(variant, document):
    """Render the generated dependency edges without expanding runtime matrices."""
    owners = {key('b', graph.build_group(path)): tests.display_group(graph.build_group(path))
              for path in source_scopes()[0]}
    rows = [f'# Configuration graph: `{variant}`', '',
            'Generated by `python3 dev/ci_schedule.py --write`; edit the scheduling sources.', '',
            'Each node is a declared stage or matrix group within one compatible configuration.',
            'Matrix children and live status appear in the GitHub Actions run graph.',
            'Arrows are the actual `needs` dependencies; this is a scheduling diagram, not a live report.', '',
            '[Graph index](README.md) · [Architecture](../../workflows/README.md)', '',
            '```mermaid', 'flowchart LR']
    for identity, job in document['jobs'].items():
        label = 'Build / ' + owners[identity] if identity in owners else job['name']
        style = 'build' if identity == 'Common' or identity in owners else (
            'coverage' if identity.startswith('Publish_') else 'execution')
        rows.append(f'    {identity}["{html.escape(label, quote=True)}"]:::{style}')
    for identity, job in document['jobs'].items():
        for dependency in job.get('needs', []):
            rows.append(f'    {dependency} --> {identity}')
    rows += ['    classDef build fill:#dbeafe,stroke:#2563eb,color:#111827',
             '    classDef execution fill:#dcfce7,stroke:#16a34a,color:#111827',
             '    classDef coverage fill:#f3e8ff,stroke:#9333ea,color:#111827', '```', '',
             'Blue: compilation. Green: execution or baseline inventory. Purple: coverage publication.', '']
    return '\n'.join(rows)


def graph_summary(variants):
    """Link the live run graph and commit-matched diagrams from the existing Plan job."""
    repository = os.environ['GITHUB_SERVER_URL'] + '/' + os.environ['GITHUB_REPOSITORY']
    run = repository + '/actions/runs/' + os.environ['GITHUB_RUN_ID']
    source = repository + '/blob/' + os.environ['GITHUB_SHA']
    rows = ['## Test job graphs', '', f'[Open the live workflow graph and job logs]({run})', '',
            'The native graph shows live job status and dependencies as matrices expand.',
            'The diagrams below show declared matrix groups, not individual test cases or live status.', '',
            '| Configuration variant | Dependency diagram | Workflow source |', '| --- | --- | --- |']
    for variant in variants:
        rows.append(f'| `{variant}` | [Graph]({source}/.github/ci/graphs/{variant}.md) | '
                    f'[YAML]({source}/.github/workflows/TestConfiguration_{variant}.yml) |')
    with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as output:
        output.write('\n'.join(rows) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--prepare')
    parser.add_argument('--summary', action='store_true', help='Add graph navigation to the CI job summary')
    args = parser.parse_args()
    if args.prepare:
        prepare(args.prepare)
        return
    documents = generated()
    variants = {}
    for name, document in documents.items():
        path = ROOT / '.github/workflows' / name
        text = '# Generated by dev/ci_schedule.py --write; edit registration or .github/ci/targeted-lanes.yml.\n'
        text += yaml.safe_dump(document, sort_keys=False, width=120)
        if len(text.encode()) > 500_000:
            raise ValueError('Generated workflow exceeds GitHub file limit: ' + name)
        if args.write:
            path.write_text(text)
        elif not path.is_file() or path.read_text() != text:
            raise ValueError('Stale hierarchical workflow; run python3 dev/ci_schedule.py --write: ' + name)
        if name.startswith('TestConfiguration_'):
            variant = name.removeprefix('TestConfiguration_').removesuffix('.yml')
            variants[variant] = workflow_diagram(variant, document)
        print(name, ':', len(document['jobs']), 'explicit dependency groups;', len(text.encode()), 'bytes')
    index = ['# Test job graphs', '',
             'Open **Actions → Tests → a workflow run** for the native live dependency graph.',
             'Click a job to inspect its logs. Matrix groups expand as their inputs become available.', '',
             'The existing CMake / CI job summary links the live run and these commit-matched diagrams.',
             'These generated diagrams show declared stage groups within each compatible configuration;',
             'they do not show expanded matrix children, completion status or coverage percentages.', '',
             '[Architecture and category overview](../../workflows/README.md#shared-builds-and-isolated-execution)', '',
             '## Configuration diagrams', '']
    index += [f'- [{variant}]({variant}.md)' for variant in variants]
    index += ['', 'Regenerate with `python3 dev/ci_schedule.py --write`; plain invocation checks consistency.', '']
    diagrams = {**variants, 'README': '\n'.join(index)}
    for name, text in diagrams.items():
        path = ROOT / '.github/ci/graphs' / (name + '.md')
        if args.write:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        elif not path.is_file() or path.read_text() != text:
            raise ValueError('Stale dependency diagram; run python3 dev/ci_schedule.py --write: ' + str(path))
    if args.summary:
        graph_summary(variants)


if __name__ == '__main__':
    main()

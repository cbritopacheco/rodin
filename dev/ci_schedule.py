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
            if not any(name.lower() in ['add_executable', 'add_test', 'rodin_add_petsc_refinement']
                       for name, _, _ in commands):
                continue
            # A helper can declare a target, e.g. rodin_add_petsc_refinement.
            # Literal first-argument references bound its source directory;
            # the configured manifest must confirm that bound before use.
            declarations = {args[0] for _, args, _ in commands if args and args[0].startswith('Rodin')}
            scopes[owner] = {'owner': owner, 'execution': tests.execution_group(owner)}
            for name in sorted(declarations):
                target = {'owner': 'tests/' + owner, 'target': name}
                branch = scope(target)
                if name in declared and declared[name] != branch:
                    raise ValueError('Ambiguous source dependency edge: ' + name)
                declared[name] = branch
                if branch != owner:
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


def layout():
    scopes, declared = source_scopes()
    ordinary, targeted = requests()
    executions = {}
    rows = [(family, {'job': 'Ordinary', 'matrix': row})
            for family, entries in ordinary.items() for row in entries]
    rows += [(r['job'], r) for r in targeted]
    for family, request in rows:
        job, matrix = request['job'], request['matrix']
        identity = execution_key(request, family)
        if identity in executions:
            continue
        if job == 'Ordinary':
            group = matrix['group']
            dependencies = set(scopes) if group == 'src/Rodin' else {
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
        executions[identity] = {'job': job, 'family': family, 'name': label,
                                'dependencies': sorted(key('b', s) for s in dependencies)}
    return {key('b', s): dict(info, path=s) for s, info in scopes.items()}, executions


def prepare(common_id):
    """Validate every actual dependency before emitting the configuration DAG."""
    entries = json.loads((ROOT / 'ci-common-plan.json').read_text())
    entry = next(e for e in entries if e['id'] == common_id)
    index, _ = graph.index_graph(ROOT / 'ci-common-record', common_id=common_id)
    producers, executions = layout()
    build_rows, execution_rows = {}, {}
    build_owner = {}
    for identity, record in index['builds'].items():
        spec = record['definition']
        targets = [{'owner': 'tests/' + spec['directory'], 'target': t} for t in spec['targets']]
        owners = {key('b', scope(t)) for t in targets}
        if len(owners) != 1 or not owners <= producers.keys():
            raise ValueError('Configured producer has no unique declared hierarchy edge')
        owner = owners.pop()
        build_owner[identity] = owner
        label = tests.display_group(producers[owner]['path'])
        if producers[owner]['path'] == spec['directory'] and 'part' in spec:
            label += ' / ' + str(spec['part'])
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
                if contributing:
                    output.write('p_' + family + '=true\n')
    print('Independent hierarchical producers:', len(build_rows), '; execution groups:', len(execution_rows))


def expression(value):
    return '${{ ' + value + ' }}'


def configuration_workflow():
    producers, executions = layout()
    common_source = yaml.safe_load((ROOT / '.github/workflows/TestCommon.yml').read_text())
    common = copy.deepcopy(common_source['jobs']['Build'])
    common['name'] = 'Build / Common libraries'
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
                         'if': expression("needs.Common.outputs." + identity + " != ''"),
                         'strategy': {'fail-fast': False, 'matrix': expression('fromJSON(needs.Common.outputs.' + identity + ')')},
                         'uses': './.github/workflows/TestDirectoryBuild.yml',
                         'with': {'id': expression('matrix.id'), 'runner': expression('matrix.runner'),
                                  'common': expression('inputs.id')}}
    for identity, execution in sorted(executions.items()):
        dependencies = ['Common', *execution['dependencies']]
        condition = ("!cancelled() && needs.Common.result == 'success' && needs.Common.outputs." + identity +
                     " != '' && !contains(needs.*.result, 'failure') && !contains(needs.*.result, 'cancelled')")
        job = {'name': execution['name'], 'needs': dependencies, 'if': expression(condition)}
        if execution['job'] == 'Ordinary':
            job.update({'strategy': {'fail-fast': False, 'matrix': expression('fromJSON(needs.Common.outputs.' + identity + ')')},
                        'uses': './.github/workflows/TestGroup.yml',
                        'with': {k: expression('matrix.' + k) for k in ['family', 'backend', 'profile', 'group',
                                 'build_type', 'multithreaded', 'asan', 'ubsan', 'id', 'name', 'baseline']}})
        else:
            job.update({'uses': './.github/workflows/TestTargeted.yml',
                        'with': {'family': execution['job'], 'matrix': expression('needs.Common.outputs.' + identity)}})
        job['secrets'] = 'inherit'
        jobs[identity] = job
    for family in sorted({e['family'] for e in executions.values()}):
        members = [k for k, e in executions.items() if e['family'] == family]
        jobs['Publish_' + family] = {
            'name': 'Coverage / Publish / ' + family,
            'needs': ['Common', *members],
            'if': expression("!cancelled() && needs.Common.outputs.p_" + family +
                             " == 'true' && !contains(needs.*.result, 'failure') && !contains(needs.*.result, 'cancelled')"),
            'uses': './.github/workflows/CoveragePublish.yml', 'secrets': 'inherit',
            'with': {'kind': 'ordinary' if family in requests()[0] else 'targeted', 'family': family,
                     'common': expression('inputs.id')}}
    return {'name': 'Tests / Configuration', 'on': {'workflow_call': {'inputs': {
        'definition': {'type': 'string', 'required': True}, 'id': {'type': 'string', 'required': True}},
        'secrets': {'CODECOV_TOKEN': {'required': False}}}},
        'permissions': {'contents': 'read', 'actions': 'read'}, 'jobs': jobs}


def targeted_workflow():
    jobs = copy.deepcopy(ci_lanes.jobs())
    for name, job in jobs.items():
        job.pop('needs', None)
        job['if'] = expression("inputs.family == '" + name + "'")
        job['strategy']['matrix'] = expression('fromJSON(inputs.matrix)')
        for step in job['steps']:
            if step.get('id') == 'coverage_record':
                step['with'].pop('artifact_name', None)
    return {'name': 'Tests / Dedicated selection', 'on': {'workflow_call': {'inputs': {
        'family': {'type': 'string', 'required': True}, 'matrix': {'type': 'string', 'required': True}},
        'secrets': {'CODECOV_TOKEN': {'required': False}}}},
        'permissions': {'contents': 'read', 'actions': 'read'}, 'jobs': jobs}


def generated():
    return {'TestConfiguration.yml': configuration_workflow(), 'TestTargeted.yml': targeted_workflow()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--prepare')
    args = parser.parse_args()
    if args.prepare:
        prepare(args.prepare)
        return
    for name, document in generated().items():
        path = ROOT / '.github/workflows' / name
        text = '# Generated by dev/ci_schedule.py --write; edit registration or .github/ci/targeted-lanes.yml.\n'
        text += yaml.safe_dump(document, sort_keys=False, width=120)
        if len(text.encode()) > 500_000:
            raise ValueError('Generated workflow exceeds GitHub file limit: ' + name)
        if args.write:
            path.write_text(text)
        elif not path.is_file() or path.read_text() != text:
            raise ValueError('Stale hierarchical workflow; run python3 dev/ci_schedule.py --write: ' + name)
        print(name, ':', len(document['jobs']), 'explicit dependency groups;', len(text.encode()), 'bytes')


if __name__ == '__main__':
    main()

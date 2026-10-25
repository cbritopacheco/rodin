#!/usr/bin/env python3
"""Build once, execute independently, and certify complete benchmark campaigns."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess

import ci_builds as archive
from benchmark_site import inventory, publication_path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / 'build-benchmarks'


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def emit(key, value):
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
            stream.write(key + '=' + json.dumps(value, separators=(',', ':')) + '\n')


def configuration(name):
    if name not in ('local', 'distributed'):
        raise ValueError('Unknown benchmark configuration')
    distributed = name == 'distributed'
    flags = {'CMAKE_BUILD_TYPE': 'Release', 'CMAKE_CXX_COMPILER': 'g++-14',
             'CMAKE_C_COMPILER': 'gcc-14', 'CMAKE_SUPPRESS_REGENERATION': 'ON',
             'RODIN_BUILD_SRC': 'ON', 'RODIN_BUILD_BENCHMARKS': 'ON',
             'RODIN_BUILD_UNIT_TESTS': 'OFF', 'RODIN_BUILD_MANUFACTURED_TESTS': 'OFF',
             'RODIN_BUILD_CONVERGENCE_TESTS': 'OFF', 'RODIN_BUILD_EXAMPLES': 'OFF',
             'RODIN_BUILD_DOC': 'OFF', 'RODIN_INSTALL_RESOURCES': 'OFF',
             'RODIN_CODE_COVERAGE': 'OFF', 'RODIN_MULTITHREADED': 'OFF',
             'RODIN_USE_ASAN': 'OFF', 'RODIN_USE_UBSAN': 'OFF',
             'RODIN_USE_MPI': 'ON' if distributed else 'OFF',
             'RODIN_USE_PETSC': 'ON' if distributed else 'OFF'}
    return {'kind': 'benchmark', 'mode': name, 'compiler': 'gcc-14',
            'runner': 'ubuntu-24.04', 'cmake': flags, 'targets': []}


def group_path(row):
    parts = publication_path(row['path']).parts
    if parts[0] == 'Rodin':
        return '/'.join(parts[:3] if parts[1] in ('MPI', 'PETSc') else parts[:2])
    return parts[0]


def make_plan(rows, spec):
    groups = defaultdict(list)
    for row in rows:
        if row['state'] != 'registered':
            continue
        distributed = row['path'].startswith(('Rodin/MPI/', 'Rodin/PETSc/'))
        if distributed != (spec['mode'] == 'distributed'):
            continue
        groups[group_path(row)].append(row)
    if not groups:
        raise ValueError('Empty benchmark configuration')
    result = []
    for path, selected in sorted(groups.items()):
        definition = {**spec, 'targets': sorted(row['target'] for row in selected), 'stage': 'components'}
        result.append({'id': archive.identity(definition), 'label': path.replace('Rodin/', 'Modules/'),
                       'definition': definition, 'entries': selected})
    return {'schema': 1, 'commit': archive.commit(), 'configuration': spec['mode'], 'groups': result}


def normalize_checkout():
    """Give identical tracked sources stable mtimes before reusing Make objects."""
    if subprocess.run(['git', 'diff', '--quiet', 'HEAD', '--'], cwd=ROOT).returncode:
        raise ValueError('Shared benchmark builds require a clean committed checkout')
    # An old fixed epoch is independent of clone time and cannot be a future mtime.
    source_epoch = 946684800
    names = subprocess.check_output(['git', 'ls-files', '-z', '--recurse-submodules'], cwd=ROOT, text=True)
    for name in names.split('\0'):
        if name:
            path = ROOT / name
            if path.exists() or path.is_symlink():
                os.utime(path, (source_epoch, source_epoch), follow_symlinks=False)


def library_targets(build, selected):
    """Use CMake's actual dependency graph, including linked third-party libraries."""
    reply = build / '.cmake/api/v1/reply'
    index = json.loads(sorted(reply.glob('index-*.json'))[-1].read_text())
    model = json.loads((reply / index['reply']['codemodel-v2']['jsonFile']).read_text())
    nodes, names = {}, {}
    for target in model['configurations'][0]['targets']:
        row = json.loads((reply / target['jsonFile']).read_text())
        nodes[row['id']] = row
        names[row['name']] = row['id']
    pending = [names[target] for target in selected]
    closure = set()
    while pending:
        key = pending.pop()
        if key in closure:
            continue
        closure.add(key)
        pending.extend(item['id'] for item in nodes[key].get('dependencies', []))
    return sorted(nodes[key]['name'] for key in closure if nodes[key]['type'] in
                  ('STATIC_LIBRARY', 'SHARED_LIBRARY', 'MODULE_LIBRARY', 'OBJECT_LIBRARY'))


def prepare(name, destination):
    spec = configuration(name)
    normalize_checkout()
    query = BUILD / '.cmake/api/v1/query/codemodel-v2'
    query.parent.mkdir(parents=True, exist_ok=True)
    query.touch()
    subprocess.run(['cmake', '-S', str(ROOT), '-B', str(BUILD)] +
                   ['-D' + key + '=' + value for key, value in spec['cmake'].items()], check=True)
    rows = inventory(BUILD / 'RodinBenchmarkTargets-Release.tsv')
    plan = make_plan(rows, spec)
    selected = [entry['target'] for group in plan['groups'] for entry in group['entries']]
    common = {**spec, 'targets': library_targets(BUILD, selected),
              'stage': 'common'}
    if not common['targets']:
        raise ValueError('Shared benchmark build has no production library targets')
    archive.build(common, BUILD)
    archive.pack(common, BUILD, destination)
    plan['common'] = common
    dump(destination / 'plan.json', plan)
    emit('groups', [{'id': group['id'], 'label': group['label']} for group in plan['groups']])


def command(entry, ranks):
    result = []
    if ranks > 1:
        cache = (BUILD / 'CMakeCache.txt').read_text()
        match = re.search(r'^MPIEXEC_EXECUTABLE:[^=]+=(.*)$', cache, re.M)
        if not match:
            raise ValueError('Distributed benchmark has no configured MPI launcher')
        result = [match[1], '-n', str(ranks)]
    return result + [entry['executable']] + [v for v in entry['args'].split(';') if v] + (
        [v for v in entry.get('mpi_args', '').split(';') if v] if ranks > 1 else [])


def csv_case(row):
    # Discovery and measurement use different stream precision for distortion.
    return '/'.join([row['geometry'], str(int(row['degree'])), str(int(row['resolution'])),
                     format(finite(float(row['distortion'])), '.12g'), str(int(row['warp_degree']))])


def discover(entry, ranks):
    argv = command(entry, ranks)
    if entry['format'] == 'csv':
        output = subprocess.check_output(argv + ['--list'], text=True)
        cases = list(csv.DictReader(io.StringIO(output)))
        if not cases or any(set(row) != {'geometry', 'degree', 'resolution', 'distortion', 'warp_degree'} for row in cases):
            raise ValueError('Malformed CSV workload discovery')
        names = [csv_case(row) for row in cases]
    else:
        output = subprocess.check_output(argv + ['--benchmark_list_tests=true'], text=True)
        names = [line.strip() for line in output.splitlines() if line.strip()]
        counts = Counter(names)
        if any(count not in (1, ranks) for count in counts.values()):
            raise ValueError('Ambiguous or inconsistent distributed benchmark discovery')
        names = sorted(counts)
    if not names or len(names) != len(set(names)):
        raise ValueError('Empty or duplicate benchmark case discovery')
    return names


def partitions(names, count):
    """Bound both case counts and the platform's single-argument regex length."""
    chunk, size = [], 0
    for name in names:
        length = len(re.escape(name).encode()) + 1
        if length > 60000:
            raise ValueError('A benchmark case name exceeds the selector budget')
        if chunk and (len(chunk) == count or size + length > 60000):
            yield chunk
            chunk, size = [], 0
        chunk.append(name)
        size += length
    if chunk:
        yield chunk


def build_group(source, group_id, destination):
    plan = json.loads((source / 'plan.json').read_text())
    group = next(item for item in plan['groups'] if item['id'] == group_id)
    normalize_checkout()
    archive.restore(plan['common'], BUILD, source)
    archive.build(group['definition'], BUILD)
    lanes = []
    discoveries = []
    for entry in group['entries']:
        for ranks in map(int, entry['ranks'].split(';')):
            names = discover(entry, ranks)
            discoveries.append({'target': entry['target'], 'ranks': ranks, 'cases': names})
            batch = int(entry.get('case_batch') or 256)
            if entry['format'] == 'csv':
                batch = min(batch, 8)
            for part, selected in enumerate(partitions(names, batch), 1):
                lane = {'entry': entry, 'ranks': ranks, 'cases': selected,
                        'configuration': plan['configuration'], 'definition': group['definition']}
                lane['id'] = archive.identity(lane)
                lane['label'] = entry['path'].replace('Rodin/', 'Modules/') + '/' + entry['profile'].title() + '/NP' + str(ranks) + '/' + str(part)
                lanes.append(lane)
    if not lanes or len(lanes) > 256:
        raise ValueError('Benchmark group needs a narrower hierarchy: lane matrix exceeds 256')
    archive.pack(group['definition'], BUILD, destination)
    dump(destination / 'execution.json', {'schema': 1, 'commit': plan['commit'], 'group': group, 'lanes': lanes, 'discoveries': discoveries})
    emit('lanes', [{'id': lane['id'], 'label': lane['label']} for lane in lanes])


def finite(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError('Invalid benchmark measurement')
    return value


def google_results(report, expected, metric):
    groups = defaultdict(list)
    repetitions = defaultdict(set)
    if not isinstance(report.get('benchmarks'), list):
        raise ValueError('Missing benchmark records')
    for row in report['benchmarks']:
        if row.get('error_occurred') or row.get('error_message') or row.get('skipped'):
            raise ValueError('Benchmark error or skipped case: ' + row.get('name', '?'))
        if row.get('run_type') == 'aggregate':
            continue
        name = row.get('run_name', row.get('name', ''))
        repetition = row.get('repetition_index', 0)
        if repetition in repetitions[name]:
            raise ValueError('Duplicate benchmark repetition: ' + name)
        repetitions[name].add(repetition)
        if row.get('iterations', 0) <= 0:
            raise ValueError('Benchmark completed no iterations')
        selected = metric
        if metric == 'auto':
            selected = 'real_time' if '/real_time' in name or '/manual_time' in name else 'cpu_time'
        value = finite(row[selected])
        unit = row['time_unit']
        if unit not in ('ns', 'us', 'ms', 's'):
            raise ValueError('Unknown benchmark time unit')
        label = 'manual_time' if '/manual_time' in name and selected == 'real_time' else selected
        groups[name].append((value, unit, label, row))
    if set(groups) != set(expected):
        raise ValueError('Incomplete benchmark cases: missing=' + str(sorted(set(expected) - set(groups))) + ' unexpected=' + str(sorted(set(groups) - set(expected))))
    results = []
    for name, samples in sorted(groups.items()):
        if len({(unit, label) for _, unit, label, _ in samples}) != 1:
            raise ValueError('Mixed measurement units or metrics')
        declared = {row.get('repetitions') for _, _, _, row in samples}
        if declared != {None}:
            count = next(iter(declared))
            if len(declared) != 1 or not isinstance(count, int) or isinstance(count, bool) or count < 1 or repetitions[name] != set(range(count)):
                raise ValueError('Incomplete benchmark repetitions: ' + name)
        values = [value for value, _, _, _ in samples]
        results.append({'name': name, 'value': statistics.median(values), 'unit': samples[0][1],
                        'metric': samples[0][2], 'samples': values, 'records': [row for _, _, _, row in samples],
                        'context': report.get('context', {})})
    return results


def csv_results(text, expected):
    rows = list(csv.DictReader(io.StringIO(text)))
    groups = defaultdict(set)
    results = []
    for row in rows:
        key = csv_case(row)
        operation = (row['query'], row['pruning'])
        if operation in groups[key] or row['pruning'] not in ('0', '1') or int(row['samples']) <= 0:
            raise ValueError('Duplicate or invalid CSV workload row')
        groups[key].add(operation)
        value = finite(float(row['ns']))
        results.append({'name': key + '/' + '/'.join(operation), 'value': value,
                        'unit': 'ns', 'metric': 'process_cpu_time', 'samples': [value], 'records': [row]})
    if set(groups) != set(expected):
        raise ValueError('Incomplete CSV workload configurations')
    for key, operations in groups.items():
        queries = ['interior', 'miss', 'build'] if key.startswith('Point/') else ['interior', 'shared', 'miss', 'build']
        if operations != {(query, pruning) for query in queries for pruning in ('0', '1')}:
            raise ValueError('Incomplete CSV workload operations: ' + key)
    return results


def run_lane(source, lane_id, destination):
    execution = json.loads((source / 'execution.json').read_text())
    lane = next(item for item in execution['lanes'] if item['id'] == lane_id)
    archive.restore(lane['definition'], BUILD, source)
    destination.mkdir(parents=True, exist_ok=False)
    entry = lane['entry']
    argv = command(entry, lane['ranks'])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if entry['format'] == 'googlecpp':
        report = destination / 'report.json'
        # Escape exact discovered case names; never use a broad partial-match selector.
        selector = '^(' + '|'.join(re.escape(name) for name in lane['cases']) + ')$'
        subprocess.run(argv + ['--benchmark_filter=' + selector, '--benchmark_min_time=0.2s',
                              '--benchmark_out_format=json', '--benchmark_out=' + str(report.resolve())],
                       check=True, env=env)
        values = google_results(json.loads(report.read_text()), lane['cases'], entry['metric'])
    else:
        report = destination / 'report.csv'
        buffer = io.StringIO();writer = None
        for case in lane['cases']:
            geometry, degree, resolution, distortion, warp = case.split('/')
            output = subprocess.check_output(argv + ['0.003', geometry, '/'.join([degree, distortion, resolution]), warp], text=True, env=env)
            reader = csv.DictReader(io.StringIO(output))
            if writer is None:
                writer = csv.DictWriter(buffer, fieldnames=reader.fieldnames);writer.writeheader()
            writer.writerows(reader)
        report.write_text(buffer.getvalue())
        values = csv_results(report.read_text(), lane['cases'])
    dump(destination / 'receipt.json', {'schema': 1, 'commit': archive.commit(), 'lane': lane,
                                       'report': report.name, 'sha256': archive.digest(report), 'results': values})


def pool(common_sources, group_sources, result_sources, output):
    plans = [json.loads((path / 'plan.json').read_text()) for path in common_sources]
    commits = {plan['commit'] for plan in plans}
    if commits != {archive.commit()}:
        raise ValueError('Mixed or stale benchmark commits')
    if {plan['configuration'] for plan in plans} != {'local', 'distributed'}:
        raise ValueError('A complete campaign requires local and distributed configuration plans')
    if len({plan['configuration'] for plan in plans}) != len(plans):
        raise ValueError('Duplicate benchmark configuration plans')
    expected_groups = {group['id']: group for plan in plans for group in plan['groups']}
    executions = [json.loads((path / 'execution.json').read_text()) for path in group_sources]
    if len(executions) != len(expected_groups) or {execution['group']['id'] for execution in executions} != set(expected_groups):
        raise ValueError('Missing or duplicate benchmark build groups')
    lanes = {}
    for execution in executions:
        if execution['commit'] not in commits or execution['group'] != expected_groups[execution['group']['id']]:
            raise ValueError('Stale or changed benchmark group')
        expected = {(entry['target'], int(rank)): entry for entry in execution['group']['entries']
                    for rank in entry['ranks'].split(';')}
        discovery = {(item['target'], item['ranks']): item['cases'] for item in execution['discoveries']}
        if set(discovery) != set(expected) or len(discovery) != len(execution['discoveries']):
            raise ValueError('Incomplete case discovery inventory')
        actual = defaultdict(list)
        for lane in execution['lanes']:
            pair = (lane['entry']['target'], lane['ranks'])
            if pair not in expected or lane['entry'] != expected[pair] or lane['definition'] != execution['group']['definition']:
                raise ValueError('Execution differs from configured benchmark plan')
            if archive.identity({key: value for key, value in lane.items() if key not in ('id', 'label')}) != lane['id']:
                raise ValueError('Changed benchmark lane identity')
            actual[pair].extend(lane['cases'])
            if lane['id'] in lanes:
                raise ValueError('Duplicate execution lane')
            lanes[lane['id']] = lane
        if set(actual) != set(discovery) or any(sorted(actual[key]) != sorted(cases) or len(set(cases)) != len(cases)
                                               for key, cases in discovery.items()):
            raise ValueError('Case partitions do not exactly cover discovery')
    receipts = [json.loads((path / 'receipt.json').read_text()) for path in result_sources]
    if len(receipts) != len(lanes) or {receipt['lane']['id'] for receipt in receipts} != set(lanes):
        raise ValueError('Missing or duplicate benchmark contributions')
    components = defaultdict(list)
    configurations = {}
    seen = set()
    for path, receipt in zip(result_sources, receipts):
        lane = lanes[receipt['lane']['id']]
        if receipt['commit'] not in commits or receipt['lane'] != lane:
            raise ValueError('Changed benchmark execution receipt')
        report = path / receipt['report']
        if report.parent.resolve() != path.resolve() or archive.digest(report) != receipt['sha256']:
            raise ValueError('Benchmark report checksum or path mismatch')
        entry = lane['entry']
        values = google_results(json.loads(report.read_text()), lane['cases'], entry['metric']) if entry['format'] == 'googlecpp' else csv_results(report.read_text(), lane['cases'])
        if values != receipt['results']:
            raise ValueError('Receipt does not match raw benchmark measurements')
        configuration_id = archive.identity({key: lane['definition'][key] for key in ['compiler', 'runner', 'cmake']}) + '-np' + str(lane['ranks'])
        component = (entry['path'], configuration_id)
        configurations[component] = {'definition': lane['definition'], 'ranks': lane['ranks']}
        for value in values:
            identity = (*component, value['name'], value['metric'])
            if identity in seen:
                raise ValueError('Benchmark case contributed more than once')
            seen.add(identity)
            components[component].append(value)
    output.mkdir(parents=True, exist_ok=False)
    for (path, configuration_id), values in components.items():
        dump(output / path / configuration_id / 'results.json',
             {'schema': 1, 'commit': archive.commit(), 'configuration': configuration_id,
              'component': path, 'format': 'csv' if any(v['metric'] == 'process_cpu_time' for v in values) else 'googlecpp', **configurations[(path, configuration_id)], 'results': sorted(values, key=lambda row: row['name'])})
    dump(output / 'complete.json', {'schema': 1, 'commit': archive.commit(), 'configurations': sorted(plan['configuration'] for plan in plans),
                                    'groups': len(executions), 'lanes': len(receipts), 'cases': len(seen),
                                    'components': {str(path.relative_to(output)): archive.digest(path) for path in sorted(output.glob('**/results.json'))}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['prepare', 'build', 'run', 'pool'])
    parser.add_argument('--configuration', choices=['local', 'distributed'])
    parser.add_argument('--source', type=Path)
    parser.add_argument('--id')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--artifacts', type=Path)
    args = parser.parse_args()
    if args.operation == 'prepare': prepare(args.configuration, args.output)
    elif args.operation == 'build': build_group(args.source, args.id, args.output)
    elif args.operation == 'run': run_lane(args.source, args.id, args.output)
    else:
        pool(sorted(args.artifacts.glob('benchmark-plan-common-*')), sorted(args.artifacts.glob('benchmark-plan-group-*')),
             sorted(args.artifacts.glob('benchmark-result-*')), args.output)


if __name__ == '__main__':
    main()

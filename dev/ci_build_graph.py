#!/usr/bin/env python3
"""Common objects -> directory executables -> isolated execution workspaces.

Only the configured CMake manifest determines target ownership and eligibility.
Execution requests retain the old selectors/configurations; rank, geometry and
Slow do not belong to a build identity. Common objects are retained for builders,
whereas execution archives contain notes and executables, never runtime counters.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess

import ci_builds as archive
import ci_tests as tests

ROOT = tests.ROOT
DIRECTORY_BANDS = 2
# Limit the executable work queued behind one producer; configuration variants
# stay separate. Very large individual translation units keep their own budget.
DIRECTORY_TARGET_LIMIT = 12
# These translation units already had separate execution jobs for compile/run
# cost. Keep their compiler memory budgets isolated inside the physical tree.
ISOLATED_DIRECTORY_TARGETS = set(tests.COMPONENTS) | {'RodinMPIMatrixRangeTest'}
BACKEND_NAMES = {'optional': 'Optional', 'default': 'Default', 'baseline': 'Baseline',
                 'native': 'Native', 'mpi': 'MPI', 'petsc': 'PETSc', 'complex': 'PETSc Complex',
                 'ranges_real': 'PETSc Real', 'ranges_complex': 'PETSc Complex'}


def request_key(job, matrix):
    if job == 'Ordinary':
        matrix = {k: v for k, v in matrix.items() if k != 'family'}
    return archive.identity({'job': job, 'matrix': matrix})


def specification(job, matrix):
    if job == 'Ordinary':
        flags = tests.configuration(matrix['group'], matrix['profile'], matrix['build_type'],
                                    matrix['multithreaded'], matrix['asan'], matrix['ubsan'], matrix['backend'])
        return {'mode': matrix['backend'], 'runner': 'ubuntu-latest', 'compiler': 'gcc-14',
                'cmake': flags, 'targets': []}
    if job == 'ConvergenceBaseline':
        return archive.definition('Convergence', dict(multithreaded=matrix['multithreaded'],
                                                      suite='curved_projection')) | {'targets': []}
    if job == 'PETScMatrixRanges':
        spec = archive.definition('PETScComplexConvergence', dict(multithreaded=matrix['multithreaded'], suite='bulk'))
        spec['mode'] = 'ranges_' + matrix['scalar']
        spec['runner'] = 'ubuntu-latest'
        spec['cmake']['RODIN_BUILD_CONVERGENCE_TESTS'] = 'OFF'
        # These are the original explicit build targets; execution regex/ranks
        # stay in the workflow and are independently guarded by the inventory.
        spec['targets'] = ['RodinMPIMatrixRangeTest', 'RodinPETScGridFunctionTest',
                           'RodinPETScMPIGridFunctionTest', 'RodinPETScMPIP0P0gTest',
                           'RodinPETScMPISubMeshGridFunctionTest', 'RodinPETScIOHDF5IOTest',
                           'RodinPETScIOMEDITIOTest', 'RodinPETScIOMFEMIOTest', 'RodinPETScIOMFEMMPIIOTest',
                           'RodinPETScSingleFieldConstraintTest', 'RodinPETScRestrictedNonlocalTest',
                           'RodinPETScMixedTargetAssemblyTest', 'RodinPETScMPIMixedTargetAssemblyTest',
                           'RodinPETScSingleFieldTargetAssemblyTest', 'RodinPETScSNESStateTest']
        return spec
    return archive.definition(job, matrix)


def common_specification(spec):
    flags = dict(spec['cmake'])
    # Registration switches do not affect the production ABI. A superset of
    # registrations lets regular/Slow and different suites share the same tree.
    for suite in ['UNIT', 'MANUFACTURED', 'CONVERGENCE']:
        flags['RODIN_BUILD_' + suite + '_TESTS'] = 'ON'
    flags['CMAKE_SUPPRESS_REGENERATION'] = 'ON'
    return {**spec, 'cmake': flags, 'targets': [], 'stage': 'common'}


def prepare(ordinary, targeted):
    requests = [{'job': 'Ordinary', 'matrix': row} for rows in ordinary.values() for row in rows] + targeted
    common = {}
    for request in requests:
        spec = common_specification(specification(**request))
        identity = archive.identity(spec)
        common.setdefault(identity, {'id': identity, 'definition': spec, 'requests': []})['requests'].append(request)
    if len(common) > 256:
        raise ValueError('Common configuration matrix exceeds GitHub limit')
    return list(common.values())


def stamp_sources():
    # Checkouts on a later runner otherwise make headers newer than restored
    # objects. Only tracked inputs (including submodules) get this stable time.
    paths = subprocess.check_output(['git', '-C', str(ROOT), 'ls-files', '--recurse-submodules', '-z']).split(b'\0')
    for raw in paths:
        if raw:
            path = ROOT / os.fsdecode(raw)
            if path.is_file() and not path.is_symlink():
                os.utime(path, (946684800, 946684800))


def original_inventory(build, request):
    spec = specification(**request)
    switches = {'unit': 'UNIT', 'manufactured': 'MANUFACTURED', 'convergence': 'CONVERGENCE'}
    return [target for target in tests.targets(build)
            if not target['owner'].startswith('tests/') or any(
                target['owner'].startswith('tests/' + suite + '/') and
                spec['cmake'].get('RODIN_BUILD_' + switch + '_TESTS') == 'ON'
                for suite, switch in switches.items())]


def selected(build, request):
    job, matrix = request['job'], request['matrix']
    configured = {t['target']: t for t in original_inventory(build, request)}
    if job == 'Ordinary':
        chosen = tests.selected_targets(build, matrix['group'], profile=matrix['profile'], baseline=matrix['baseline'])
        return [t for t in chosen if t['target'] in configured]
    if job == 'ConvergenceBaseline':
        return tests.selected_targets(build, matrix['group'], profile='convergence-baseline')
    names = specification(job, matrix)['targets']
    missing = set(names) - configured.keys()
    if missing:
        raise ValueError('Missing explicitly required targets: ' + str(sorted(missing)))
    return [configured[name] for name in names]


def directory(target):
    # Build ownership follows physical registration directories. Execution can
    # batch short siblings or split logical components independently of it.
    group = target['owner'].removeprefix('tests/')
    if not (ROOT / 'tests' / group / 'CMakeLists.txt').is_file():
        raise ValueError('Target has no source directory: ' + str(target))
    return group


def compile_common(entry, build, destination):
    spec = entry['definition']
    stamp_sources()
    query = build / '.cmake/api/v1/query/codemodel-v2'
    query.parent.mkdir(parents=True, exist_ok=True)
    query.touch()
    flags = dict(spec['cmake'])
    if spec['mode'] in ['complex', 'ranges_complex']:
        flags['PETSc_DIR'] = subprocess.check_output(['pkg-config', '--variable=prefix', 'PETSc'], text=True).strip()
    subprocess.run(['cmake', '-S', str(ROOT), '-B', str(build)] +
                   ['-D' + k + '=' + v for k, v in flags.items()], check=True)
    reply = build / '.cmake/api/v1/reply'
    index = json.loads(next(reply.glob('index-*.json')).read_text())
    model = json.loads((reply / index['reply']['codemodel-v2']['jsonFile']).read_text())
    nodes = {}
    by_name = {}
    for target in model['configurations'][0]['targets']:
        data = json.loads((reply / target['jsonFile']).read_text())
        nodes[data['id']] = data
        by_name[data['name']] = data['id']
    selections = [(request, selected(build, request)) for request in entry['requests']]
    pending = [by_name[t['target']] for _, chosen in selections for t in chosen]
    closure = set()
    while pending:
        key = pending.pop()
        if key in closure:
            continue
        closure.add(key)
        pending.extend(d['id'] for d in nodes[key].get('dependencies', []))
    libraries = [nodes[key]['name'] for key in closure if nodes[key]['type'] in
                 ['STATIC_LIBRARY', 'SHARED_LIBRARY', 'MODULE_LIBRARY', 'OBJECT_LIBRARY']]
    if not libraries:
        raise ValueError('Common build has no library targets')
    print('Common library dependency closure:', sorted(libraries), flush=True)
    parallel = '1' if (spec['mode'] == 'optional' and flags['CMAKE_BUILD_TYPE'] == 'Release' and
                       flags['RODIN_USE_ASAN'] == 'OFF' and flags['RODIN_USE_UBSAN'] == 'OFF') else '2'
    subprocess.run(['cmake', '--build', str(build), '--parallel', parallel, '--target'] + sorted(libraries),
                   check=True, env=dict(os.environ, GCOV_PREFIX=str(ROOT / 'coverage-discovery-excluded')))
    groups, requests = {}, {}
    for request, chosen in selections:
        owned = set()
        for target in chosen:
            if target['kind'] != 'test':
                continue
            group = directory(target)
            owned.add(group)
            groups.setdefault(group, {})[target['target']] = target
        requests[request_key(**request)] = {'request': request, 'selected': chosen,
                                           'configured': original_inventory(build, request),
                                           'directories': sorted(owned), 'common': entry['id']}
    specs = {}
    for group, targets in groups.items():
        names = sorted(set(targets) - ISOLATED_DIRECTORY_TARGETS)
        chunks = [names[i:i + DIRECTORY_TARGET_LIMIT] for i in range(0, len(names), DIRECTORY_TARGET_LIMIT)]
        chunks += [[name] for name in sorted(set(targets) & ISOLATED_DIRECTORY_TARGETS)]
        for part, chunk in enumerate(chunks, 1):
            row = {**spec, 'stage': 'directory', 'common': entry['id'], 'directory': group,
                   'targets': chunk}
            label = group
            if len(chunks) > 1:
                row['part'] = part
                label += '/' + str(part)
            specs[label] = row
    for record in requests.values():
        wanted = {t['target'] for t in record['selected']}
        record['builds'] = [archive.identity(row) for row in specs.values()
                            if wanted.intersection(row['targets'])]
        record.pop('directories')
    document = {'commit': archive.commit(), 'common': entry['id'], 'definition': spec,
                'requests': requests, 'directories': specs}
    destination.mkdir(parents=True)
    (destination / 'graph.json').write_text(json.dumps(document, indent=2) + '\n')
    archive.pack(spec, build, ROOT / 'ci-common-artifact')


def index_graph(folder, mode=None, common_id=None):
    requests, builds, common = {}, {}, {}
    for path in folder.rglob('graph.json'):
        graph = json.loads(path.read_text())
        if graph['commit'] != archive.commit() or graph['common'] in common:
            raise ValueError('Stale or duplicate common build inventory')
        common[graph['common']] = graph['definition']
        for key, value in graph['requests'].items():
            if key in requests:
                raise ValueError('Duplicate execution request')
            requests[key] = value
        for group, spec in graph['directories'].items():
            key = archive.identity(spec)
            if key in builds:
                raise ValueError('Duplicate directory build')
            builds[key] = {'id': key, 'name': tests.display_group(group), 'definition': spec}
    expected = [e for e in json.loads((ROOT / 'ci-common-plan.json').read_text())
                if (mode is None or e['definition']['mode'] == mode)
                and (common_id is None or e['id'] == common_id)]
    expected_requests = {request_key(**r) for e in expected for r in e['requests']}
    if set(common) != {e['id'] for e in expected} or set(requests) != expected_requests:
        raise ValueError('Missing common builds or execution requests')
    for entry in expected:
        if common[entry['id']] != entry['definition']:
            raise ValueError('Common build configuration changed')
        for request in entry['requests']:
            record = requests[request_key(**request)]
            if record['request'] != request or record['common'] != entry['id']:
                raise ValueError('Changed execution/build binding')
            selected_names = {t['target'] for t in record['selected'] if t['kind'] == 'test'}
            restored_names = set()
            for key in record['builds']:
                if key not in builds or builds[key]['definition']['common'] != entry['id']:
                    raise ValueError('Incompatible directory build binding')
                restored_names.update(builds[key]['definition']['targets'])
            if not selected_names <= restored_names:
                raise ValueError('Selected executable has no directory producer')
    owners = set()
    for row in builds.values():
        for target in row['definition']['targets']:
            owner = (row['definition']['common'], target)
            if owner in owners:
                raise ValueError('Executable compiled by more than one directory group')
            owners.add(owner)
    # Optional Unit spans sixteen configurations and hundreds of executables;
    # four matrix bands keep it below GitHub's 256-child limit. Other backend
    # branches need two, including the two-producer native MPI restriction set.
    band_count = min(len(builds), 4 if mode == 'optional' else DIRECTORY_BANDS)
    bands = {f'directory_{i + 1}': [] for i in range(band_count)}
    for i, row in enumerate(sorted(builds.values(), key=lambda r: (r['name'], r['id']))):
        spec = row['definition']
        flags = spec['cmake']
        label = (row['name'] + f" ({spec['compiler'].upper()}, {BACKEND_NAMES[spec['mode']]}, {flags['CMAKE_BUILD_TYPE']}, "
                 f"MT-{flags['RODIN_MULTITHREADED']}, ASan-{flags['RODIN_USE_ASAN']}, "
                 f"UBSan-{flags['RODIN_USE_UBSAN']})")
        # Large definitions travel as an artifact, not through GitHub's 1 MiB
        # job-output allowance (which counts matrix JSON as UTF-16).
        bands[f'directory_{i % band_count + 1}'].append(
            {'id': row['id'], 'name': label, 'runner': spec['runner'], 'mode': spec['mode']})
    if any(not rows or len(rows) > 256 for rows in bands.values()):
        raise ValueError('Directory build bands must contain 1–256 jobs')
    return {'commit': archive.commit(), 'requests': requests, 'common': common, 'builds': builds}, bands


def download(name, destination):
    if os.environ.get('CI_RESTORED_ARTIFACTS'):
        source = Path(os.environ['CI_RESTORED_ARTIFACTS']) / name
        if not (source / 'manifest.json').is_file():
            raise ValueError('Required shared build artifact was not downloaded: ' + name)
        return source
    subprocess.run(['gh', 'run', 'download', os.environ['GITHUB_RUN_ID'], '--name', name,
                    '--dir', str(destination)], check=True)
    return destination


def compile_directory(spec, build):
    stamp_sources()
    artifact = ROOT / 'ci-common-input'
    artifact = download('ci-common-' + spec['common'], artifact) or artifact
    common = {k: v for k, v in spec.items() if k not in ['common', 'directory', 'part']}
    common.update(stage='common', targets=[])
    archive.restore(common, build, artifact)
    objects = {p: (archive.digest(p), p.stat().st_mtime_ns) for p in build.rglob('*.o')}
    # Keep the existing single-compiler safeguards for memory-heavy MPI and
    # optimized optional-backend instantiations.
    flags = spec['cmake']
    parallel = '1' if ('MPI/Variational' in spec['directory'] or
                      spec['mode'] == 'optional' and flags['CMAKE_BUILD_TYPE'] == 'Release' and
                      flags['RODIN_USE_ASAN'] == 'OFF' and flags['RODIN_USE_UBSAN'] == 'OFF') else '2'
    subprocess.run(['cmake', '--build', str(build), '--parallel', parallel, '--target'] + spec['targets'],
                   check=True, env=dict(os.environ, GCOV_PREFIX=str(ROOT / 'coverage-discovery-excluded')))
    if any(not p.is_file() or (archive.digest(p), p.stat().st_mtime_ns) != state for p, state in objects.items()):
        raise ValueError('Directory build recompiled common objects')
    archive.pack(spec, build, ROOT / 'ci-build-artifact')
    destination = ROOT / 'ci-built-record'
    destination.mkdir()
    (destination / 'build-receipt.json').write_text(json.dumps(
        {'commit': archive.commit(), 'definition': spec, 'common_objects_unchanged': True,
         'archive_sha256': archive.digest(ROOT / 'ci-build-artifact/build.tar.gz')}, indent=2) + '\n')


def restore_request(index, job, matrix, build):
    if job == 'Ordinary':
        matrix = {k: v for k, v in matrix.items() if k != 'family'}
    if index['commit'] != archive.commit():
        raise ValueError('Stale execution/build index')
    record = index['requests'][request_key(job, matrix)]
    if record['request'] != {'job': job, 'matrix': matrix}:
        raise ValueError('Wrong execution/build request')
    common = record['common']
    if not record['builds']:
        folder = ROOT / 'ci-restore' / common
        folder = download('ci-common-' + common, folder) or folder
        archive.restore(index['common'][common], build, folder)
    for i, key in enumerate(record['builds']):
        folder = ROOT / 'ci-restore' / key
        folder = download('ci-build-' + key, folder) or folder
        archive.restore(index['builds'][key]['definition'], build, folder, overlay=i > 0)
    for path in build.rglob('*.o'):
        path.unlink()
    for target in record['selected']:
        if target['executable'] and not Path(target['executable']).is_file():
            raise ValueError('Requested executable missing after restore: ' + target['target'])
    (build / 'ci-built-targets.json').write_text(json.dumps(record['selected'], indent=2) + '\n')
    (build / 'ci-configured-targets.json').write_text(json.dumps(record['configured'], indent=2) + '\n')
    (build / 'ci-build-binding.json').write_text(json.dumps(record, indent=2) + '\n')


def verify_builds(index, folder, commit):
    if index['commit'] != commit:
        raise ValueError('Stale shared-build plan')
    received = set()
    for path in folder.rglob('build-receipt.json'):
        receipt = json.loads(path.read_text())
        key = archive.identity(receipt['definition'])
        if (key not in index['builds'] or key in received or receipt['commit'] != commit or
                receipt['definition'] != index['builds'][key]['definition'] or
                not receipt['common_objects_unchanged']):
            raise ValueError('Invalid shared directory build receipt')
        received.add(key)
    if received != index['builds'].keys():
        raise ValueError('Missing shared directory builds')


def load_indices(folder, commit, expected):
    combined = {'commit': commit, 'requests': {}, 'common': {}, 'builds': {}}
    for path in folder.rglob('ci-build-index.json'):
        index = json.loads(path.read_text())
        if index['commit'] != commit:
            raise ValueError('Stale backend build index')
        for category in ['requests', 'common', 'builds']:
            if combined[category].keys() & index[category].keys():
                raise ValueError('Duplicate backend build inventory')
            combined[category].update(index[category])
    requests = {request_key(**r) for e in expected for r in e['requests']}
    if combined['common'].keys() != {e['id'] for e in expected} or combined['requests'].keys() != requests:
        raise ValueError('Missing backend build inventory')
    return combined


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['prepare', 'common', 'index', 'directory', 'directory-definition', 'mode', 'identify', 'restore'])
    parser.add_argument('--id')
    parser.add_argument('--mode')
    parser.add_argument('--definition')
    parser.add_argument('--job')
    parser.add_argument('--matrix')
    parser.add_argument('--build', type=Path, default=ROOT / 'build')
    args = parser.parse_args()
    if args.operation == 'prepare':
        rows = prepare(json.loads((ROOT / 'ci-test-plan.json').read_text()),
                       json.loads((ROOT / 'ci-coverage-targets.json').read_text()))
        (ROOT / 'ci-common-plan.json').write_text(json.dumps(rows, indent=2) + '\n')
        print('Common configurations:', len(rows), '; execution requests:', sum(len(e['requests']) for e in rows))
        if os.environ.get('GITHUB_OUTPUT'):
            with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
                stream.write('configurations=' + json.dumps({'include': [
                    {'id': e['id'], 'definition': e['definition']} for e in rows]}, separators=(',', ':')) + '\n')
                for mode in sorted({e['definition']['mode'] for e in rows}):
                    stream.write('common_' + mode + '=' + json.dumps(
                        {'include': [{'id': e['id'], 'definition': e['definition']}
                                     for e in rows if e['definition']['mode'] == mode]}, separators=(',', ':')) + '\n')
    elif args.operation == 'common':
        entries = json.loads((ROOT / 'ci-common-plan.json').read_text())
        spec = json.loads(args.definition)
        entry = next(e for e in entries if e['id'] == archive.identity(spec))
        compile_common(entry, args.build, ROOT / 'ci-common-record')
    elif args.operation == 'directory':
        index = json.loads((ROOT / 'ci-build-index.json').read_text())
        compile_directory(index['builds'][args.id]['definition'], args.build)
    elif args.operation == 'index':
        index, bands = index_graph(ROOT / 'ci-common-records', args.mode)
        (ROOT / 'ci-build-index.json').write_text(json.dumps(index, indent=2) + '\n')
        (ROOT / 'ci-directory-plan.json').write_text(json.dumps(bands, indent=2) + '\n')
        print('Directory builds:', len(index['builds']), '; matrix bands:', {k: len(v) for k, v in bands.items()})
        if os.environ.get('GITHUB_OUTPUT'):
            with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
                for name, rows in bands.items():
                    stream.write(name + '=' + json.dumps({'include': rows}, separators=(',', ':')) + '\n')
    elif args.operation == 'mode':
        with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
            spec = specification(args.job, json.loads(args.matrix))
            stream.write('mode=' + spec['mode'] + '\n')
            stream.write('common=' + archive.identity(common_specification(spec)) + '\n')
    elif args.operation == 'directory-definition':
        index = json.loads((ROOT / 'ci-build-index.json').read_text())
        if index['commit'] != archive.commit():
            raise ValueError('Stale directory build index')
        with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
            stream.write('definition=' + json.dumps(index['builds'][args.id]['definition'], separators=(',', ':')) + '\n')
    else:
        index = json.loads((ROOT / 'ci-build-index.json').read_text())
        matrix = json.loads(args.matrix)
        if args.operation == 'restore':
            restore_request(index, args.job, matrix, args.build)
        else:
            spec = index['common'][index['requests'][request_key(args.job, matrix)]['common']]
            with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
                stream.write('definition=' + json.dumps(spec, separators=(',', ':')) + '\n')
                record = index['requests'][request_key(args.job, matrix)]
                names = ['ci-build-' + k for k in record['builds']] or ['ci-common-' + record['common']]
                pattern = names[0] if len(names) == 1 else '{' + ','.join(names) + '}'
                stream.write('artifacts=' + pattern + '\n')


if __name__ == '__main__':
    main()

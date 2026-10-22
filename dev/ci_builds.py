#!/usr/bin/env python3
"""Share compiled test families between independent execution partitions.

A build is identified by its complete CMake configuration and target set. Rank,
geometry and case selectors belong to execution lanes, not to the build key.
Each consumer restores an immutable archive into its own runner workspace.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = {
    'Convergence': 'native',
    'MPIH1Restriction': 'restriction',
    'PETScConvergence': 'petsc_h',
    'PETScDegreeConvergence': 'petsc_degree',
    'PETScCurvedConvergence': 'petsc_curved',
    'PETScCurvedBoundaryConvergence': 'petsc_boundary',
    'PETScCurvedNonlinearTetrahedronConvergence': 'petsc_boundary',
    'PETScComplexConvergence': 'complex',
}


def definition(job, row):
    family = FAMILIES[job]
    mt = row['multithreaded']
    if mt not in ['ON', 'OFF']:
        raise ValueError('Invalid multithreading configuration')
    suite = row.get('suite', '')
    physics = row.get('physics', '')
    allowed = ['Poisson', 'Conductivity', 'LinearElasticity', 'Stokes',
               'ReactionDiffusion', 'NonlinearPoisson', 'Diffusion', 'P0Projection', 'GeometryApproximation']
    if physics and physics not in allowed:
        raise ValueError('Unknown physics: ' + physics)
    if family == 'native':
        if suite.startswith('pressure_'):
            names = ['RodinConvergenceHStokesStability']
        elif suite == 'curved_projection':
            names = ['RodinConvergenceIsoparametricP0Projection',
                     'RodinConvergenceIsoparametricGeometryApproximation']
        else:
            choices = {'diffusion': 'Diffusion', 'reaction': 'ReactionDiffusion',
                       'nonlinear': 'NonlinearPoisson', 'elasticity': 'LinearElasticity',
                       'stokes': 'Stokes', 'helmholtz': 'Helmholtz'}
            name = next((value for key, value in choices.items()
                         if suite == 'curved_' + key or suite.startswith('curved_' + key + '_')), None)
            if name is None:
                raise ValueError('Unknown native build suite: ' + suite)
            names = ['RodinConvergenceIsoparametric' + name]
    elif family == 'restriction':
        names = ['RodinMPISubMeshGridFunctionTest']
    elif family == 'petsc_h':
        if suite == 'bulk':
            names = ['RodinConvergenceHPETSc' + p for p in
                     ['Poisson', 'MPIPoisson', 'Conductivity', 'LinearElasticity',
                      'ReactionDiffusion', 'NonlinearPoisson', 'Stokes']]
        elif suite.startswith('pressure_'):
            names = ['RodinConvergenceHPETScStokesStability']
        elif suite.endswith('_boundary') and physics:
            names = ['RodinConvergenceHPETSc' + physics + 'Boundary']
        else:
            raise ValueError('Unknown PETSc h build suite: ' + suite)
    elif family == 'petsc_degree':
        if not physics:
            raise ValueError('Degree refinement requires physics')
        names = ['RodinConvergence' + degree + 'PETSc' + physics for degree in ['P', 'HP']]
    elif family == 'petsc_curved':
        if not physics:
            raise ValueError('Curved refinement requires physics')
        names = ['RodinConvergenceIsoparametricPETSc' + physics]
    elif family == 'petsc_boundary':
        physics = physics or 'NonlinearPoisson'
        names = ['RodinConvergenceIsoparametricPETSc' + physics + 'Boundary']
    else:
        if suite == 'bulk':
            names = ['RodinConvergenceHPETScHelmholtz', 'RodinPETScSNESStateTest',
                     'RodinConvergencePPETScHelmholtz', 'RodinConvergenceHPPETScHelmholtz',
                     'RodinPETScScalarIntegralTest', 'RodinVariationalIntegralTest']
        elif suite.startswith('curved_helmholtz'):
            names = ['RodinConvergenceIsoparametricPETScHelmholtz']
        else:
            choices = {'curved_projection': 'RodinConvergenceIsoparametricPETScP0Projection',
                       'boundary': 'RodinConvergenceHPETScHelmholtzBoundary',
                       'submesh': 'RodinPETScMPIH1SubMeshRestrictionTest',
                       'curved_boundary': 'RodinConvergenceIsoparametricPETScHelmholtzBoundary'}
            if suite not in choices:
                raise ValueError('Unknown complex build suite: ' + suite)
            names = [choices[suite]]
    mode = 'native' if family == 'native' else 'mpi' if family == 'restriction' else 'complex' if family == 'complex' else 'petsc'
    flags = {
        'CMAKE_BUILD_TYPE': 'RelWithDebInfo', 'CMAKE_SUPPRESS_REGENERATION': 'ON',
        'RODIN_CODE_COVERAGE': 'ON',
        'RODIN_LCOV_FLAGS': '--ignore-errors gcov,gcov --ignore-errors mismatch --ignore-errors unused',
        'RODIN_BUILD_SRC': 'ON', 'RODIN_BUILD_UNIT_TESTS': 'ON' if mode in ['mpi', 'complex'] else 'OFF',
        'RODIN_BUILD_MANUFACTURED_TESTS': 'OFF', 'RODIN_BUILD_CONVERGENCE_TESTS': 'OFF' if mode == 'mpi' else 'ON',
        'RODIN_BUILD_EXAMPLES': 'OFF', 'RODIN_BUILD_BENCHMARKS': 'OFF', 'RODIN_INSTALL_RESOURCES': 'OFF',
        'RODIN_MULTITHREADED': mt, 'RODIN_USE_MPI': 'OFF' if mode == 'native' else 'ON',
        'RODIN_USE_PETSC': 'ON' if mode in ['petsc', 'complex'] else 'OFF',
        'RODIN_USE_ASAN': 'OFF', 'RODIN_USE_UBSAN': 'OFF',
    }
    return {'mode': mode, 'compiler': 'gcc-14',
            'runner': 'ubuntu-latest' if mode == 'native' else 'ubuntu-24.04',
            'cmake': flags, 'targets': sorted(names)}


def identity(spec):
    return hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:24]


def plan(workflow):
    families = {}
    for job, family in FAMILIES.items():
        for row in workflow['jobs'][job]['strategy']['matrix']['include']:
            spec = definition(job, row)
            key = identity(spec)
            families.setdefault(family, {})[key] = {'id': key, 'definition': spec}
    for family, rows in families.items():
        if len(rows) > 256:
            raise ValueError('Shared build family exceeds matrix limit: ' + family)
    return {family: list(rows.values()) for family, rows in families.items()}


def build(spec, build_dir):
    flags = dict(spec['cmake'])
    if spec['mode'] == 'complex':
        flags['PETSc_DIR'] = subprocess.check_output(
            ['pkg-config', '--variable=prefix', 'PETSc'], text=True).strip()
    subprocess.run(['cmake', '-S', str(ROOT), '-B', str(build_dir)] +
                   ['-D' + key + '=' + value for key, value in flags.items()], check=True)
    # The configured manifest remains authoritative for enabled target names.
    from ci_tests import targets
    enabled = {row['target'] for row in targets(build_dir)}
    missing = set(spec['targets']) - enabled
    if missing:
        raise ValueError('Unregistered build targets: ' + ', '.join(sorted(missing)))
    env = dict(os.environ, GCOV_PREFIX=str(ROOT / 'coverage-discovery-excluded'))
    subprocess.run(['cmake', '--build', str(build_dir), '-j', '2', '--target'] + spec['targets'],
                   env=env, check=True)


def commit():
    return subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()


def compiler(spec):
    return subprocess.check_output([spec['compiler'].replace('gcc-', 'g++-'), '-dumpfullversion', '-dumpversion'], text=True).strip()


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def runtime_libraries(spec, build_dir):
    from ci_tests import targets
    paths = set()
    for target in targets(build_dir):
        if target['target'] not in spec['targets']:
            continue
        output = subprocess.check_output(['ldd', target['executable']], text=True)
        if 'not found' in output:
            raise ValueError('Unresolved test runtime dependency')
        paths.update(re.findall(r'(/\S+)\s+\(', output))
    # CTest uses the configured launcher path. Detect a changed MPI alternative
    # even when both implementations' libraries are installed on the runner.
    if spec['cmake'].get('RODIN_USE_MPI') == 'ON':
        cache = (build_dir / 'CMakeCache.txt').read_text()
        launcher = re.search(r'^MPIEXEC_EXECUTABLE:[^=]+=(.*)$', cache, re.M)
        if launcher and Path(launcher[1]).is_file():
            paths.add(launcher[1])
    return {path: digest(Path(path)) for path in sorted(paths)}


def pack(spec, build_dir, destination):
    if not next(build_dir.rglob('*.gcno'), None):
        raise ValueError('Shared instrumented build has no coverage notes')
    destination.mkdir(parents=True, exist_ok=False)
    metadata = {'schema': 1, 'definition': spec, 'identity': identity(spec),
                'commit': commit(), 'workspace': str(ROOT), 'build': str(build_dir.resolve()),
                'architecture': platform.machine(), 'compiler_version': compiler(spec),
                'runtime_libraries': runtime_libraries(spec, build_dir)}
    archive = destination / 'build.tar.gz'
    with tarfile.open(archive, 'w:gz', compresslevel=1) as tar:
        for path in sorted(build_dir.rglob('*')):
            # Executables, CTest files, CMake metadata and zero-hit coverage
            # notes travel together. Objects and runtime counters do not.
            if path.is_dir() or (path.is_file() and path.suffix not in ['.o', '.gcda']):
                if path.is_symlink():
                    path.resolve().relative_to(build_dir.resolve())
                tar.add(path, arcname=str(path.relative_to(build_dir)), recursive=False)
    metadata['sha256'] = digest(archive)
    (destination / 'manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')


def restore(spec, build_dir, source):
    metadata = json.loads((source / 'manifest.json').read_text())
    expected = {'schema': 1, 'definition': spec, 'identity': identity(spec),
                'commit': commit(), 'workspace': str(ROOT), 'build': str(build_dir.resolve()),
                'architecture': platform.machine(), 'compiler_version': compiler(spec)}
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError('Incompatible shared test build: ' + key)
    archive = source / 'build.tar.gz'
    if digest(archive) != metadata['sha256']:
        raise ValueError('Shared test build checksum mismatch')
    if build_dir.exists():
        raise ValueError('Restore requires an empty build location')
    build_dir.mkdir(parents=True)
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            destination = build_dir / member.name
            destination.resolve().relative_to(build_dir.resolve())
            if member.name.endswith('.gcda'):
                raise ValueError('Shared test build contains runtime coverage counters')
            if member.issym() or member.islnk():
                target = (destination.parent if member.issym() else build_dir) / member.linkname
                target.resolve().relative_to(build_dir.resolve())
        options = {'filter': 'data'} if hasattr(tarfile, 'data_filter') else {}
        tar.extractall(build_dir, **options)
    for path, expected_hash in metadata['runtime_libraries'].items():
        if not Path(path).is_file() or digest(Path(path)) != expected_hash:
            raise ValueError('Incompatible test runtime library: ' + path)
    for path in build_dir.rglob('*.gcda'):
        raise ValueError('Shared test build contains runtime coverage counters: ' + str(path))
    print('Restored independent test build ' + identity(spec))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['key', 'build', 'pack', 'restore'])
    parser.add_argument('--job')
    parser.add_argument('--matrix')
    parser.add_argument('--definition')
    parser.add_argument('--build', type=Path, default=ROOT / 'build')
    parser.add_argument('--artifact', type=Path, default=ROOT / 'ci-build-artifact')
    args = parser.parse_args()
    spec = json.loads(args.definition) if args.definition else definition(args.job, json.loads(args.matrix))
    if args.operation == 'key':
        print(identity(spec))
        if os.environ.get('GITHUB_OUTPUT'):
            with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
                stream.write('id=' + identity(spec) + '\n')
    elif args.operation == 'build':
        build(spec, args.build)
    elif args.operation == 'pack':
        pack(spec, args.build, args.artifact)
    else:
        restore(spec, args.build, args.artifact)


if __name__ == '__main__':
    main()

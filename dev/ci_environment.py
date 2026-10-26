#!/usr/bin/env python3
"""Bounded Ubuntu dependency provisioning and conditional hosted-runner cleanup.

Package archives are an optimization only: APT still validates versions,
dependencies, signatures and archive hashes. No installed system tree or test
execution state is cached. Run this helper as the runner user, never with sudo.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from urllib.parse import urlsplit


GIB = 1024 ** 3
# Conservative initial headroom policy, not a measured peak for every module.
# Log remaining space so this can be refined using real hierarchy-build peaks.
MIN_FREE_GIB = 8
CLEANUP_PATHS = ('/usr/share/dotnet', '/usr/local/lib/android', '/opt/ghc',
                 '/opt/hostedtoolcache/CodeQL')
# A healthy install took about 40 seconds in the measured campaign. These
# budgets allow slower runners while bounding the observed 25–54 minute stalls.
INDEX_TIMEOUT = 90
DOWNLOAD_TIMEOUT = 300
INSTALL_TIMEOUT = 900
MIRRORS = (None, 'https://archive.ubuntu.com/ubuntu',
           'https://mirrors.edge.kernel.org/ubuntu')
UBUNTU_URI = re.compile(r'https?://(?:azure\.archive|archive|security)\.ubuntu\.com/ubuntu/?')
YAML_VERSION = '6.0.2'
# A wheel install normally takes seconds. Bound network stalls to three minutes
# per attempt and retry whole pip invocations, including malformed index JSON.
PYTHON_INSTALL_TIMEOUT = 180
PYTHON_INSTALL_ATTEMPTS = 3
PYTHON_RETRY_DELAY = 2


def install_python():
    """Reuse the exact working PyYAML version or install it with bounded retries."""
    probe = [sys.executable, '-c',
             f'import yaml; raise SystemExit(yaml.__version__ != {YAML_VERSION!r})']
    if subprocess.run(probe, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                      timeout=30).returncode == 0:
        return
    command = [sys.executable, '-m', 'pip', 'install', '--disable-pip-version-check',
               '--no-cache-dir', '--retries', '3', '--timeout', '30', 'PyYAML==' + YAML_VERSION]
    last_error = None
    for attempt in range(1, PYTHON_INSTALL_ATTEMPTS + 1):
        try:
            subprocess.run(command, check=True, timeout=PYTHON_INSTALL_TIMEOUT)
            subprocess.run(probe, check=True, timeout=30)
            return
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            last_error = error
            print(f'Pinned PyYAML installation attempt {attempt} failed.', flush=True)
            if attempt < PYTHON_INSTALL_ATTEMPTS:
                time.sleep(PYTHON_RETRY_DELAY * attempt)
    raise RuntimeError('Could not install and import PyYAML==' + YAML_VERSION) from last_error


def select_compiler(compiler):
    """Verify a matched GCC toolchain and export it for subsequent CMake steps."""
    match = re.fullmatch(r'gcc-([0-9]+)', compiler)
    if not match:
        raise ValueError('Expected an explicit GCC major version')
    cxx = compiler.replace('gcc-', 'g++-', 1)
    versions = [subprocess.check_output([tool, '-dumpfullversion', '-dumpversion'],
                                      text=True, timeout=30).strip() for tool in [compiler, cxx]]
    if any(version.split('.')[0] != match[1] for version in versions) or versions[0] != versions[1]:
        raise RuntimeError('Requested GCC and G++ versions do not match: ' + repr(versions))
    subprocess.run(['cmake', '--version'], check=True, timeout=30)
    output = os.environ.get('GITHUB_ENV')
    if not output:
        raise RuntimeError('Compiler selection requires GITHUB_ENV')
    with Path(output).open('a') as stream:
        stream.write(f'CC={compiler}\nCXX={cxx}\n')
    print(f'Selected {compiler} / {cxx}: {versions[0]}', flush=True)


def disk_space(path, minimum_gib=MIN_FREE_GIB):
    """Delete allowlisted unused tools only below the explicit free-space floor."""
    required = int(minimum_gib * GIB)
    free = shutil.disk_usage(path).free
    print(f'Disk: {free / GIB:.2f} GiB free; required headroom {minimum_gib:g} GiB', flush=True)
    if free >= required:
        return
    if os.environ.get('GITHUB_ACTIONS') != 'true' or os.environ.get('RUNNER_ENVIRONMENT') != 'github-hosted':
        raise RuntimeError('Insufficient disk space; tool cleanup is restricted to GitHub-hosted runners')
    for candidate in CLEANUP_PATHS:
        if Path(candidate).exists():
            start = time.monotonic()
            subprocess.run(['sudo', 'rm', '-rf', '--', candidate], check=True)
            free = shutil.disk_usage(path).free
            print(f'Cleanup {candidate}: {time.monotonic() - start:.1f}s; {free / GIB:.2f} GiB free', flush=True)
            if free >= required:
                return
    raise RuntimeError('Insufficient disk space after allowlisted tool cleanup')


def packages(value):
    result = sorted(set(shlex.split(value)))
    if not result or any(not re.fullmatch(r'[a-z0-9][a-z0-9+.:\-]*(?:=[A-Za-z0-9.+:~\-]+)?', p) for p in result):
        raise ValueError('Expected a nonempty list of Debian package names or exact versions')
    return result


def without_microsoft_ubuntu(text, deb822=False):
    """Exclude only Microsoft's unused Ubuntu feeds, preserving other trust rules."""
    def unrelated(uri):
        parsed = urlsplit(uri)
        return parsed.hostname == 'packages.microsoft.com' and parsed.path.startswith('/ubuntu/')
    if not deb822:
        return ''.join(line for line in text.splitlines(keepends=True)
                       if not (re.match(r'^\s*deb(?:-src)?\s', line) and
                               any(unrelated(uri) for uri in re.findall(r'https?://[^\s]+', line))))
    result = []
    changed = False
    for paragraph in re.split(r'\n[ \t]*\n', text):
        field = re.search(r'^URIs:[ \t]*([^\n]*(?:\n[ \t]+[^\n]*)*)', paragraph, re.M)
        if field:
            uris = field[1].split()
            retained = [uri for uri in uris if not unrelated(uri)]
            if len(retained) != len(uris):
                changed = True
                if not retained:
                    continue
                paragraph = paragraph[:field.start()] + 'URIs: ' + ' '.join(retained) + paragraph[field.end():]
        result.append(paragraph)
    return '\n\n'.join(result) if changed else text


def setup_sources(system=Path('/etc/apt')):
    """Prepare disposable hosted-runner feeds before setup-cpp runs apt-get."""
    if (os.environ.get('GITHUB_ACTIONS') != 'true' or
            os.environ.get('RUNNER_ENVIRONMENT') != 'github-hosted'):
        print('APT source preparation is restricted to GitHub-hosted runners', flush=True)
        return
    entries = [system / 'sources.list']
    entries += sorted((system / 'sources.list.d').glob('*.list'))
    entries += sorted((system / 'sources.list.d').glob('*.sources'))
    for source in entries:
        if not source.is_file():
            continue
        original = source.read_text()
        prepared = without_microsoft_ubuntu(original, source.suffix == '.sources')
        if prepared == original:
            continue
        backup = source.with_name(source.name + '.rodin-backup')
        if not backup.exists():
            subprocess.run(['sudo', 'cp', '--', str(source), str(backup)], check=True)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.apt') as temporary:
            temporary.write(prepared)
            temporary.flush()
            subprocess.run(['sudo', 'install', '-m', '644', temporary.name, str(source)], check=True)
        print('Excluded unused Microsoft Ubuntu feed from ' + str(source), flush=True)


def sources(directory, mirror, system=Path('/etc/apt')):
    """Copy source definitions; change Ubuntu transport only, retaining trust rules."""
    directory.mkdir(parents=True, exist_ok=True)
    parts = directory / 'sources.list.d'
    parts.mkdir(exist_ok=True)
    main = system / 'sources.list'
    entries = [main] if main.exists() else []
    entries += sorted((system / 'sources.list.d').glob('*.list'))
    entries += sorted((system / 'sources.list.d').glob('*.sources'))
    for source in entries:
        text = source.read_text()
        if mirror:
            # Hosted images can hide the selected archive behind APT's local
            # mirror-list transport. Replacing direct URLs alone leaves those
            # definitions downloading from the same stalled Azure endpoint.
            for reference in re.findall(r'mirror\+file:(/[^\s]+)', text):
                mirror_file = Path(reference)
                if not mirror_file.is_file():
                    continue
                urls = [line.split()[0] for line in mirror_file.read_text().splitlines()
                        if line.strip() and not line.lstrip().startswith('#')]
                if urls and all(UBUNTU_URI.fullmatch(url) for url in urls):
                    text = text.replace('mirror+file:' + reference, mirror.rstrip('/') + '/')
            text = UBUNTU_URI.sub(mirror.rstrip('/') + '/', text)
        target = directory / 'sources.list' if source == main else parts / source.name
        target.write_text(text)
    (directory / 'sources.list').touch(exist_ok=True)
    return ['-o', 'Dir::Etc::sourcelist=' + str(directory / 'sources.list'),
            '-o', 'Dir::Etc::sourceparts=' + str(parts)]


def apt(arguments, source_options, archives, timeout, capture=False):
    # GNU timeout runs with APT's privileges and terminates its process group,
    # including transport workers, before fallback tries to acquire the locks.
    command = ['sudo', 'env', 'DEBIAN_FRONTEND=noninteractive', 'LC_ALL=C',
               'timeout', '--signal=TERM', '--kill-after=5', str(timeout),
               'apt-get', '-o', 'Acquire::Retries=1',
               '-o', 'Acquire::http::Timeout=20', '-o', 'Acquire::https::Timeout=20',
               '-o', 'APT::Update::Error-Mode=any',
               '-o', 'APT::Keep-Downloaded-Packages=true',
               '-o', 'Dir::Cache::archives=' + str(archives), *source_options, *arguments]
    result = subprocess.run(command, check=True, text=True,
                            stdout=subprocess.PIPE if capture else None,
                            env=dict(os.environ, DEBIAN_FRONTEND='noninteractive', LC_ALL='C'))
    return result.stdout if capture else None


def candidates(requested, source_options):
    """Pin candidates from the selected source before resolving the full plan."""
    result = {}
    for package in requested:
        if '=' in package:
            name, version = package.split('=', 1)
        else:
            name = package
            policy = subprocess.check_output(['apt-cache', *source_options, 'policy', name], text=True,
                                             env=dict(os.environ, LC_ALL='C'))
            match = re.search(r'^\s*Candidate: (\S+)', policy, re.M)
            if not match or match[1] == '(none)':
                raise ValueError('No APT candidate for ' + name)
            version = match[1]
        result[name] = version
    return result


def prepare(requested, archives, state):
    archives.mkdir(parents=True, exist_ok=True)
    (archives / 'partial').mkdir(exist_ok=True)
    state.mkdir(parents=True, exist_ok=True)
    last_error = None
    for index, mirror in enumerate(MIRRORS):
        options = sources(state / str(index), mirror)
        try:
            apt(['update'], options, archives, INDEX_TIMEOUT)
            break
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            last_error = error
            print('APT index attempt failed:', mirror or 'runner sources', flush=True)
    else:
        raise RuntimeError('All signed APT index sources failed') from last_error
    versions = candidates(requested, options)
    pinned = [name + '=' + version for name, version in sorted(versions.items())]
    simulation = apt(['--simulate', '--no-remove', 'install', *pinned], options,
                     archives, INDEX_TIMEOUT, capture=True)
    for name, version in re.findall(r'^Inst (\S+)(?: \[[^\]]+\])? \((\S+)', simulation, re.M):
        if name in versions and versions[name] != version:
            raise ValueError('Conflicting resolved version for ' + name)
        versions[name] = version
    plan = {'requested': requested, 'versions': versions, 'source_index': index}
    (state / 'plan.json').write_text(json.dumps(plan, sort_keys=True) + '\n')
    system = Path('/etc/os-release').read_text()
    architecture = subprocess.check_output(['dpkg', '--print-architecture'], text=True).strip()
    # Image identity covers the preinstalled package baseline. Compiler flags,
    # MPI rank counts and test selections do not fragment dependency caches.
    baseline = [system, architecture, os.environ.get('ImageOS'), os.environ.get('ImageVersion'), requested]
    prefix = 'rodin-apt-v1-' + hashlib.sha256(json.dumps(baseline, sort_keys=True).encode()).hexdigest()[:24] + '-'
    key = prefix + hashlib.sha256(json.dumps(versions, sort_keys=True).encode()).hexdigest()[:24]
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'], 'a') as output:
            output.write(f'key={key}\nrestore-prefix={prefix}\n')
    print('APT archive cache key:', key, flush=True)
    return plan


def install(archives, state):
    plan = json.loads((state / 'plan.json').read_text())
    pinned = [name + '=' + version for name, version in sorted(plan['versions'].items())]
    last_error = None
    for index in range(plan['source_index'], len(MIRRORS)):
        options = sources(state / str(index), MIRRORS[index])
        try:
            if index != plan['source_index']:
                apt(['update'], options, archives, INDEX_TIMEOUT)
            # Bound the whole transfer, not only idle sockets: slow downloads
            # can otherwise run for almost an hour while still making progress.
            apt(['--yes', '--no-remove', '--download-only', 'install', *pinned],
                options, archives, DOWNLOAD_TIMEOUT)
            break
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            last_error = error
            print('APT download attempt failed:', MIRRORS[index] or 'runner sources', flush=True)
    else:
        raise RuntimeError('All APT archive sources failed for the pinned package plan') from last_error
    # Only download/index errors trigger mirror fallback. Do not retry a failed
    # package configuration or change its versions underneath shared artifacts.
    apt(['--yes', '--no-remove', '--no-download', 'install', *pinned],
        options, archives, INSTALL_TIMEOUT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['setup', 'disk', 'prepare', 'install', 'python', 'compiler'])
    parser.add_argument('--packages', default='')
    parser.add_argument('--compiler', default='')
    parser.add_argument('--minimum-gib', type=float, default=MIN_FREE_GIB)
    args = parser.parse_args()
    if args.minimum_gib <= 0:
        parser.error('Disk headroom must be positive')
    root = Path(os.environ.get('RUNNER_TEMP', tempfile.gettempdir()))
    if args.operation == 'setup':
        setup_sources()
    elif args.operation == 'disk':
        disk_space(Path.cwd(), args.minimum_gib)
    elif args.operation == 'python':
        install_python()
    elif args.operation == 'compiler':
        select_compiler(args.compiler)
    elif args.operation == 'prepare':
        prepare(packages(args.packages), root / 'rodin-apt-archives', root / 'rodin-apt-state')
    else:
        install(root / 'rodin-apt-archives', root / 'rodin-apt-state')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Bounded Ubuntu dependency provisioning and conditional hosted-runner cleanup.

APT resolves and downloads packages from signed sources. Shared test consumers
can install that exact package bundle with local APT sources only. No installed
system tree or execution state is cached. Run as the runner user, never with sudo.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import random
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
# Slow but progressing signed snapshot refreshes exceeded 90 seconds in CI.
# Allow three minutes per network refresh; local policy queries remain bounded
# separately and the privileged process group is still terminated on timeout.
INDEX_REFRESH_TIMEOUT = 180
DOWNLOAD_TIMEOUT = 300
# Cold MatrixRanges downloads of roughly 200 MB were still progressing at
# 100–170 kB/s when the five-minute watchdog killed each snapshot attempt.
# Completed archives survive between the three bounded attempts.
SNAPSHOT_DOWNLOAD_TIMEOUT = 900
INSTALL_TIMEOUT = 900
MIRRORS = (None, 'https://archive.ubuntu.com/ubuntu',
           'https://mirrors.edge.kernel.org/ubuntu')
# Snapshot normalization can send every mirror attempt to the same service.
# Spread retries beyond short 502/503 bursts and desynchronize matrix runners.
APT_RETRY_DELAYS = (30, 120)
APT_RETRY_JITTER = 30
# Hosted runners take their archive from the Azure mirror or from a regional
# EC2 mirror such as us-east-2.ec2.archive.ubuntu.com.
UBUNTU_URI = re.compile(r'https?://(?:(?:azure\.archive|[a-z0-9-]+\.ec2\.archive|archive|security)'
                       r'\.ubuntu\.com|mirrors\.edge\.kernel\.org)/ubuntu/?')
YAML_VERSION = '6.0.2'
# Shared build campaigns must not cross a live archive security update. Refresh
# this signed Ubuntu snapshot deliberately for producers and consumers together.
SHARED_APT_SNAPSHOT = '20261113T000000Z'
UBUNTU_KEYRING = Path('/usr/share/keyrings/ubuntu-archive-keyring.gpg')
# A wheel install normally takes seconds. Bound network stalls to three minutes
# per attempt and retry whole pip invocations, including malformed index JSON.
PYTHON_INSTALL_TIMEOUT = 180
PYTHON_INSTALL_ATTEMPTS = 3
PYTHON_RETRY_DELAY = 2
# The read-only version probe stalled once on a hosted runner. Retry timeouts
# within a 96-second total budget; unsuccessful CMake exits remain fatal.
CMAKE_PROBE_TIMEOUT = 30
CMAKE_PROBE_ATTEMPTS = 3
CMAKE_PROBE_RETRY_DELAY = 2


class MissingCandidate(ValueError):
    """A refreshed index did not provide every requested package."""


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
    for attempt in range(1, CMAKE_PROBE_ATTEMPTS + 1):
        try:
            subprocess.run(['cmake', '--version'], check=True, timeout=CMAKE_PROBE_TIMEOUT)
            break
        except subprocess.TimeoutExpired as error:
            print(f'CMake version probe timed out (attempt {attempt}/{CMAKE_PROBE_ATTEMPTS}).', flush=True)
            if attempt == CMAKE_PROBE_ATTEMPTS:
                raise RuntimeError('CMake version probe exhausted its bounded timeout retries') from error
            time.sleep(CMAKE_PROBE_RETRY_DELAY * attempt)

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


def sources(directory, mirror, system=Path('/etc/apt'), *, snapshot=None):
    """Copy source definitions; change Ubuntu transport only, retaining trust rules."""
    prefer_archive = snapshot is not None
    if snapshot is None:
        snapshot = os.environ.get('RODIN_APT_SNAPSHOT', '')
    if snapshot == 'shared':
        snapshot = SHARED_APT_SNAPSHOT
    if snapshot and not re.fullmatch(r'[0-9]{8}T[0-9]{6}Z', snapshot):
        raise ValueError('Invalid Ubuntu package snapshot: ' + snapshot)
    directory.mkdir(parents=True, exist_ok=True)
    parts = directory / 'sources.list.d'
    parts.mkdir(exist_ok=True)
    main = system / 'sources.list'
    entries = [main] if main.exists() else []
    entries += sorted((system / 'sources.list.d').glob('*.list'))
    entries += sorted((system / 'sources.list.d').glob('*.sources'))
    ubuntu = False
    hosted = (os.environ.get('GITHUB_ACTIONS') == 'true' and
              os.environ.get('RUNNER_ENVIRONMENT') == 'github-hosted')
    for source in entries:
        text = source.read_text()
        # Read hosted mirror lists for source detection. Preserve the runner's
        # transport on the first live attempt, and resolve it for fallback or
        # snapshot selection so those attempts do not reuse the same mirror.
        original = text
        for reference in re.findall(r'mirror\+file:(/[^\s]+)', text):
            mirror_file = Path(reference)
            if not mirror_file.is_file():
                continue
            urls = [line.split()[0] for line in mirror_file.read_text().splitlines()
                    if line.strip() and not line.lstrip().startswith('#')]
            if urls and all(UBUNTU_URI.fullmatch(url) for url in urls):
                text = text.replace('mirror+file:' + reference, (mirror or urls[0]).rstrip('/') + '/')
        if mirror:
            text = UBUNTU_URI.sub(mirror.rstrip('/') + '/', text)
        if snapshot or hosted:
            # APT::Snapshot alone schedules live indexes as well. Point the
            # Ubuntu definitions directly at the signed snapshot in that mode.
            snapshot_uri = 'https://snapshot.ubuntu.com/ubuntu/' + snapshot + '/'
            if source.suffix == '.sources':
                paragraphs = re.split(r'(\n[ \t]*\n)', text)
                for index, paragraph in enumerate(paragraphs):
                    field = re.search(r'^URIs:[ \t]*([^\n]*(?:\n[ \t]+[^\n]*)*)',
                                      paragraph, re.M | re.I)
                    if not field or not any(UBUNTU_URI.fullmatch(value) for value in field[1].split()):
                        continue
                    binary = re.search(r'^Types:[ \t]*([^\n]+)', paragraph, re.M | re.I)
                    disabled = re.search(r'^Enabled:[ \t]*no[ \t]*$', paragraph, re.M | re.I)
                    ubuntu = ubuntu or bool(binary and 'deb' in binary[1].split() and not disabled)
                    if not snapshot:
                        continue
                    paragraph = (paragraph[:field.start(1)] +
                                 re.sub(r'\S+', lambda value: snapshot_uri if UBUNTU_URI.fullmatch(value[0])
                                        else value[0], field[1]) + paragraph[field.end(1):])
                    option = re.compile(r'^Snapshot:[^\n]*(?:\n[ \t]+[^\n]*)*', re.M | re.I)
                    if option.search(paragraph):
                        paragraph = option.sub('Snapshot: no', paragraph)
                    else:
                        ending = '\n' if paragraph.endswith('\n') else ''
                        paragraph = paragraph.rstrip('\n') + '\nSnapshot: no' + ending
                    paragraphs[index] = paragraph
                text = ''.join(paragraphs)
            else:
                lines = []
                for line in text.splitlines(keepends=True):
                    entry = re.match(r'^(\s*deb(?:-src)?\s+)(\[[^]\n]*\]\s+)?(\S+)(.*)', line)
                    if entry and UBUNTU_URI.fullmatch(entry[3]):
                        ubuntu = ubuntu or bool(re.match(r'^\s*deb\s', line))
                        if snapshot:
                            flags = entry[2].strip()[1:-1].split() if entry[2] else []
                            flags = [value for value in flags if not value.lower().startswith('snapshot=')]
                            flags.append('snapshot=no')
                            line = (entry[1] + '[' + ' '.join(flags) + '] ' + snapshot_uri + entry[4] +
                                    line[len(entry[0]):])
                    lines.append(line)
                text = ''.join(lines)
        target = directory / 'sources.list' if source == main else parts / source.name
        target.write_text(text if mirror or snapshot else original)
    if (snapshot or hosted) and not ubuntu:
        # APT can refresh an empty binary source set successfully. Recover in
        # the disposable copy using the existing Ubuntu archive trust root.
        release = dict(re.findall(r'^([A-Z_]+)=(.*)$', Path('/etc/os-release').read_text(), re.M))
        codename = release.get('VERSION_CODENAME', '').strip('"')
        if release.get('ID', '').strip('"') != 'ubuntu' or not re.fullmatch(r'[a-z][a-z0-9]*', codename):
            raise RuntimeError('Cannot recover APT sources without an Ubuntu codename')
        if not UBUNTU_KEYRING.is_file():
            raise RuntimeError('Cannot recover APT sources without the Ubuntu archive keyring')
        uri = ('https://snapshot.ubuntu.com/ubuntu/' + snapshot if snapshot
               else mirror or MIRRORS[1])
        fallback = parts / 'rodin-ubuntu.sources'
        fallback.write_text(
            'Types: deb\nURIs: ' + uri.rstrip('/') + '/\n'
            'Suites: ' + ' '.join(codename + pocket for pocket in ['', '-updates', '-backports', '-security']) + '\n'
            'Components: main restricted universe multiverse\n'
            'Signed-By: ' + str(UBUNTU_KEYRING) + '\n' + ('Snapshot: no\n' if snapshot else ''))
        print('No active Ubuntu binary sources copied; using signed Ubuntu sources for ' + codename, flush=True)
    print(f'APT source copy: {len(entries)} source files; Ubuntu snapshot source available: {bool(snapshot)}', flush=True)
    (directory / 'sources.list').touch(exist_ok=True)
    options = ['-o', 'Dir::Etc::sourcelist=' + str(directory / 'sources.list'),
               '-o', 'Dir::Etc::sourceparts=' + str(parts),
               # apt-cache must see the same required binary indexes as
               # apt-get, even when a runner disables the default target.
               '-o', 'Acquire::IndexTargets::deb::Packages::DefaultEnabled=true']
    if snapshot or prefer_archive:
        # Choose a version from the signed archive even when a runner has
        # a newer build whose package file would be unavailable.
        preferences = directory / 'preferences'
        original = system / 'preferences'
        text = original.read_text() if original.is_file() else ''
        preferences.write_text(text + '\nPackage: *\nPin: release o=Ubuntu\nPin-Priority: 1001\n')
        options += ['-o', 'Dir::Etc::preferences=' + str(preferences)]
    if snapshot:
        options += ['-o', 'APT::Snapshot=' + snapshot]
    return options


def apt(arguments, source_options, archives, timeout, capture=False):
    # GNU timeout runs with APT's privileges and terminates its process group,
    # including transport workers, before fallback tries to acquire the locks.
    command = ['sudo', 'env', 'DEBIAN_FRONTEND=noninteractive', 'LC_ALL=C',
               'timeout', '--signal=TERM', '--kill-after=5', str(timeout),
               'apt-get', '-o', 'Acquire::Retries=1',
               '-o', 'Acquire::http::Timeout=20', '-o', 'Acquire::https::Timeout=20',
               '-o', 'APT::Update::Error-Mode=any',
               # Provisioning needs package indexes, not translated descriptions
               # or command-not-found suggestions, within the same time budget.
               '-o', 'Acquire::Languages=none',
               '-o', 'Acquire::IndexTargets::deb::CNF::DefaultEnabled=false',
               '-o', 'APT::Keep-Downloaded-Packages=true',
               '-o', 'Dir::Cache::archives=' + str(archives), *source_options]
    if 'install' in arguments and any(option.startswith('APT::Snapshot=') for option in source_options):
        # Only disposable CI snapshot environments permit exact-version rollback.
        command.append('--allow-downgrades')
    command += arguments
    result = subprocess.run(command, check=True, text=True,
                            stdout=subprocess.PIPE if capture else None,
                            env=dict(os.environ, DEBIAN_FRONTEND='noninteractive', LC_ALL='C'))
    return result.stdout if capture else None


def candidates(requested, source_options):
    """Pin candidates from the selected source before resolving the full plan."""
    result = dict(package.split('=', 1) for package in requested if '=' in package)
    names = sorted(set(package for package in requested if '=' not in package))
    if not names:
        return result
    # One cache read for the full closure, rather than spawning apt-cache once
    # per library on every producer and consumer.
    policy = subprocess.check_output(['apt-cache', *source_options, 'policy', *names],
                                     text=True, timeout=INDEX_TIMEOUT,
                                     env=dict(os.environ, LC_ALL='C'))
    records = dict(re.findall(r'^(\S+):\n((?:[ \t].*(?:\n|$))*)', policy, re.M))
    for name in names:
        match = re.search(r'^\s*Candidate: (\S+)', records.get(name, records.get(name.split(':', 1)[0], '')), re.M)
        if not match or match[1] == '(none)':
            raise MissingCandidate('No APT candidate for ' + name)
        result[name] = match[1]
    return result


def installed_dependencies(requested, source_options):
    """Include preinstalled transitive dependencies in explicit candidate selection.

    APT normally leaves a dependency untouched when its installed version meets
    the minimum constraint. Hosted-image rollouts can therefore produce different
    runtime bytes from identical top-level requests. Traverse only installed
    Depends/PreDepends edges; do not install unused alternatives or recommendations.
    Seed this traversal with every package selected by APT, including packages
    not installed yet: otherwise an uninstalled intermediate package hides its
    preinstalled runtime dependencies.
    """
    output = subprocess.check_output(
        ['apt-cache', *source_options, 'depends', '--recurse', '--installed',
         '--important', '--no-recommends', '--no-suggests', '--no-conflicts',
         '--no-breaks', '--no-replaces', '--no-enhances', *requested],
        text=True, timeout=INDEX_TIMEOUT, env=dict(os.environ, LC_ALL='C'))
    # Unindented names identify traversed real packages. Indented dependency
    # edges and <virtual-package> placeholders are not install requests.
    dependencies = set(re.findall(r'^([a-z0-9][a-z0-9+.\-]*(?::[a-z0-9]+)?)$', output, re.M))
    explicit = {package.split('=', 1)[0] for package in requested}
    return sorted(dependencies - explicit)


def installed_systemd_packages(architecture):
    """List installed binaries from the systemd source package."""
    output = subprocess.check_output(
        ['dpkg-query', '-W', '-f=${db:Status-Status}\t${source:Package}\t${Package}\t${Architecture}\n'],
        text=True, timeout=INDEX_TIMEOUT)
    packages = set()
    for line in output.splitlines():
        status, source, name, arch = line.split('\t')
        if status == 'installed' and source == 'systemd':
            packages.add(name if arch in {architecture, 'all'} else f'{name}:{arch}')
    return packages


def wait_for_apt_retry(attempt):
    delay = APT_RETRY_DELAYS[attempt] + random.uniform(0, APT_RETRY_JITTER)
    print(f'Waiting {delay:.1f}s before the next signed APT attempt.', flush=True)
    time.sleep(delay)


class PackageSourceUnavailable(RuntimeError):
    """The signed package indexes or archives could not be retrieved."""


def prepare(requested, archives, state, *, mirrors=MIRRORS, snapshot=None, for_bundle=False):
    archives.mkdir(parents=True, exist_ok=True)
    (archives / 'partial').mkdir(exist_ok=True)
    state.mkdir(parents=True, exist_ok=True)
    last_error = None
    for index, mirror in enumerate(mirrors):
        options = sources(state / str(index), mirror, snapshot=snapshot)
        try:
            apt(['update'], options, archives, INDEX_REFRESH_TIMEOUT)
            versions = candidates(requested, options)
            break
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, MissingCandidate) as error:
            last_error = error
            print('APT index attempt failed:', mirror or 'runner sources', str(error), flush=True)
            if index < len(mirrors) - 1:
                wait_for_apt_retry(index)
    else:
        raise PackageSourceUnavailable('All signed APT index sources failed') from last_error
    architecture = subprocess.check_output(['dpkg', '--print-architecture'], text=True).strip() if for_bundle else None
    systemd = installed_systemd_packages(architecture) if for_bundle else set()
    # APT omits satisfied installed dependencies from its Inst records, while
    # --installed traversal omits new intermediate dependencies. Combine both
    # views until the explicit plan stops growing; never remove pins or accept
    # a solver version that conflicts with an earlier signed candidate.
    while True:
        if systemd:
            selected = {name.removesuffix(':' + architecture) for name in versions}
            companions = systemd - selected
            if selected & systemd and companions:
                # Resolve lockstep dependencies before asking APT to simulate
                # an upgrade or downgrade of any single systemd library.
                versions.update(candidates(sorted(companions), options))
        pinned = [name + '=' + version for name, version in sorted(versions.items())]
        simulation = apt(['--simulate', '--no-remove', 'install', *pinned], options,
                         archives, INDEX_TIMEOUT, capture=True)
        resolved = dict(versions)
        for name, version in re.findall(r'^Inst (\S+)(?: \[[^\]]+\])? \((\S+)', simulation, re.M):
            if name in resolved and resolved[name] != version:
                raise ValueError('Conflicting resolved version for ' + name)
            resolved[name] = version
        roots = [name + '=' + version for name, version in sorted(resolved.items())]
        missing = sorted(set(installed_dependencies(roots, options)) - resolved.keys())
        if missing:
            resolved.update(candidates(missing, options))
        if resolved == versions:
            break
        versions = resolved
    print('Pinned APT package plan:', json.dumps(versions, sort_keys=True), flush=True)
    plan = {'requested': requested, 'versions': versions, 'source_index': index,
            'snapshot': next((option.partition('=')[2] for option in options
                              if option.startswith('APT::Snapshot=')), None)}
    (state / 'plan.json').write_text(json.dumps(plan, sort_keys=True) + '\n')
    system = Path('/etc/os-release').read_text()
    if architecture is None:
        architecture = subprocess.check_output(['dpkg', '--print-architecture'], text=True).strip()
    # Image identity covers the preinstalled package baseline. Compiler flags,
    # MPI rank counts and test selections do not fragment dependency caches.
    baseline = [system, architecture, os.environ.get('ImageOS'), os.environ.get('ImageVersion'), requested, os.environ.get('RODIN_APT_SNAPSHOT', '')]
    prefix = 'rodin-apt-v3-' + hashlib.sha256(json.dumps(baseline, sort_keys=True).encode()).hexdigest()[:24] + '-'
    key = prefix + hashlib.sha256(json.dumps(versions, sort_keys=True).encode()).hexdigest()[:24]
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'], 'a') as output:
            output.write(f'key={key}\nrestore-prefix={prefix}\n')
    print('APT archive cache key:', key, flush=True)
    return plan


def install(archives, state, bundle=None, *, mirrors=MIRRORS, snapshot=None):
    plan = json.loads((state / 'plan.json').read_text())
    pinned = [name + '=' + version for name, version in sorted(plan['versions'].items())]
    policy = ['--allow-downgrades'] if snapshot == '' else []
    last_error = None
    for index in range(plan['source_index'], len(mirrors)):
        options = sources(state / str(index), mirrors[index], snapshot=snapshot)
        timeout = (SNAPSHOT_DOWNLOAD_TIMEOUT if any(option.startswith('APT::Snapshot=') for option in options)
                   else DOWNLOAD_TIMEOUT)
        try:
            if index != plan['source_index']:
                apt(['update'], options, archives, INDEX_REFRESH_TIMEOUT)
            # Bound the whole transfer, not only idle sockets: slow downloads
            # can otherwise run for almost an hour while still making progress.
            download = ['--yes', '--no-remove', '--download-only', *policy]
            if bundle is not None:
                # Consumers may need packages already present on this runner.
                download.append('--reinstall')
            apt([*download, 'install', *pinned], options, archives, timeout)
            break
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            last_error = error
            print('APT download attempt failed:', mirrors[index] or 'runner sources', flush=True)
            if index < len(mirrors) - 1:
                wait_for_apt_retry(index - plan['source_index'])
    else:
        raise PackageSourceUnavailable('All APT archive sources failed for the pinned package plan') from last_error
    # Only download/index errors trigger mirror fallback. Do not retry a failed
    # package configuration or change its versions underneath shared artifacts.
    apt(['--yes', '--no-remove', '--no-download', *policy, 'install', *pinned],
        options, archives, INSTALL_TIMEOUT)
    if bundle is not None:
        export_bundle(bundle, archives, plan)


def provision_bundle(requested, root, bundle):
    archives = root / 'rodin-apt-archives'
    state = root / 'rodin-apt-state'
    # All configured mirrors map to the same endpoint in snapshot mode. On
    # failure start a new plan against signed live sources in separate APT state.
    try:
        prepare(requested, archives, state / 'snapshot', mirrors=MIRRORS[:1],
                snapshot=SHARED_APT_SNAPSHOT, for_bundle=True)
        install(archives, state / 'snapshot', bundle, mirrors=MIRRORS[:1], snapshot=SHARED_APT_SNAPSHOT)
    except PackageSourceUnavailable as error:
        print(f'::warning::{error}. Resolving a new exact plan from signed Ubuntu mirrors.', flush=True)
        prepare(requested, archives, state / 'live', snapshot='', for_bundle=True)
        install(archives, state / 'live', bundle, snapshot='')
        return
    if output := os.environ.get('GITHUB_OUTPUT'):
        with open(output, 'a') as stream:
            stream.write('snapshot=true\n')


def bundle_key(requested):
    architecture = subprocess.check_output(
        ['dpkg', '--print-architecture'], text=True, timeout=INDEX_TIMEOUT).strip()
    identity = [SHARED_APT_SNAPSHOT, requested, architecture, Path('/etc/os-release').read_text(),
                os.environ.get('ImageOS'), os.environ.get('ImageVersion'), package_digest(Path(__file__))]
    return 'rodin-apt-bundle-' + hashlib.sha256(json.dumps(identity).encode()).hexdigest()


def package_digest(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def export_bundle(bundle, archives, plan):
    """Publish the exact, signed-source package closure selected by APT."""
    architecture = subprocess.check_output(
        ['dpkg', '--print-architecture'], text=True, timeout=INDEX_TIMEOUT).strip()
    files = {}
    versions = plan['versions']
    for path in sorted(archives.glob('*.deb')):
        name, version, arch = subprocess.check_output(
            ['dpkg-deb', '--show', '--showformat=${Package}\t${Version}\t${Architecture}',
             str(path)], text=True, timeout=INDEX_TIMEOUT).split('\t')
        keys = [f'{name}:{arch}']
        if arch in {architecture, 'all'}:
            keys.append(name)
        for key in keys:
            if versions.get(key) == version:
                files[key] = path
    if missing := sorted(versions.keys() - files.keys()):
        raise RuntimeError('Missing exact package archives: ' + ', '.join(missing))
    bundle.mkdir(parents=True, exist_ok=False)
    checksums = {}
    for path in sorted(set(files.values())):
        shutil.copy2(path, bundle / path.name)
        checksums[path.name] = package_digest(path)
    (bundle / 'manifest.json').write_text(json.dumps({
        'architecture': architecture, 'snapshot': plan['snapshot'],
        'requested': plan['requested'], 'versions': versions,
        'files': {key: path.name for key, path in files.items()}, 'sha256': checksums,
    }, sort_keys=True) + '\n')


def install_bundle(bundle, requested, archives, *, snapshot_only=False):
    """Verify and install a common-build bundle without querying a remote index."""
    manifest = json.loads((bundle / 'manifest.json').read_text())
    architecture = subprocess.check_output(
        ['dpkg', '--print-architecture'], text=True, timeout=INDEX_TIMEOUT).strip()
    valid_snapshots = {SHARED_APT_SNAPSHOT} if snapshot_only else {SHARED_APT_SNAPSHOT, None}
    if (manifest['snapshot'] not in valid_snapshots or
            manifest['architecture'] != architecture or manifest['requested'] != requested):
        raise ValueError('Package bundle does not match the shared test environment')
    files = manifest['files']
    checksums = manifest['sha256']
    if (not files or files.keys() != manifest['versions'].keys() or
            set(files.values()) != checksums.keys()):
        raise ValueError('Incomplete shared test package bundle')
    paths = []
    for name, checksum in sorted(checksums.items()):
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9+._%~:-]*\.deb', name):
            raise ValueError('Invalid package filename in shared test bundle')
        path = bundle / name
        if not path.is_file() or path.is_symlink() or package_digest(path) != checksum:
            raise ValueError('Package bundle checksum or file mismatch: ' + name)
        paths.append(str(path.resolve()))
    archives.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='rodin-apt-offline-') as folder:
        empty = Path(folder)
        (empty / 'sources.list').touch()
        (empty / 'sources.list.d').mkdir()
        (empty / 'lists').mkdir()
        options = ['-o', 'Dir::Etc::sourcelist=' + str(empty / 'sources.list'),
                   '-o', 'Dir::Etc::sourceparts=' + str(empty / 'sources.list.d'),
                   '-o', 'Dir::State::lists=' + str(empty / 'lists')]
        apt(['--yes', '--no-remove', '--allow-downgrades', 'install', *paths], options, archives, INSTALL_TIMEOUT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=[
        'setup', 'disk', 'prepare', 'install', 'provision-bundle', 'install-bundle', 'bundle-key', 'python', 'compiler'])
    parser.add_argument('--packages', default='')
    parser.add_argument('--compiler', default='')
    parser.add_argument('--bundle', type=Path)
    parser.add_argument('--snapshot-only', action='store_true')
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
    elif args.operation == 'provision-bundle':
        if args.bundle is None:
            parser.error('provision-bundle requires --bundle')
        provision_bundle(packages(args.packages), root, args.bundle)
    elif args.operation == 'bundle-key':
        print(bundle_key(packages(args.packages)))
    elif args.operation == 'install-bundle':
        if args.bundle is None:
            parser.error('install-bundle requires --bundle')
        install_bundle(args.bundle, packages(args.packages), root / 'rodin-apt-archives',
                       snapshot_only=args.snapshot_only)
    else:
        install(root / 'rodin-apt-archives', root / 'rodin-apt-state', args.bundle)


if __name__ == '__main__':
    main()

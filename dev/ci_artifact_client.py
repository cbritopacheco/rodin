#!/usr/bin/env python3
"""Install the lockfile-pinned artifact client with bounded network retries."""

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import time

ATTEMPTS = 3
TRANSIENT = re.compile(r'ECONNRESET|EAI_AGAIN|ETIMEDOUT|ECONNREFUSED|ENETUNREACH|EHOSTUNREACH|ECONNABORTED|ERR_SOCKET_TIMEOUT|\bE50[234]\b')
PERMANENT = re.compile(r'\b(?:EUSAGE|EINTEGRITY|ERESOLVE|EBADENGINE|ENOENT|ELOCKVERIFY)\b')


def install(action, client, run=subprocess.run, sleep=time.sleep):
    lock = (action / 'package-lock.json').read_bytes()
    marker = client / 'installed-package-lock.json'
    package = client / 'node_modules/@actions/artifact/package.json'
    if marker.exists() and marker.read_bytes() == lock and package.exists():
        print('Pinned artifact client already installed.')
        return
    client.mkdir(parents=True, exist_ok=True)
    marker.unlink(missing_ok=True)
    if action.resolve() != client.resolve():
        for name in ('package.json', 'package-lock.json'):
            shutil.copyfile(action / name, client / name)
    command = ['npm', 'ci', '--prefix', str(client), '--ignore-scripts', '--no-audit', '--no-fund',
               '--fetch-retries=2', '--fetch-retry-mintimeout=1000',
               '--fetch-retry-maxtimeout=10000', '--fetch-timeout=60000']
    for attempt in range(1, ATTEMPTS + 1):
        result = run(command, text=True, capture_output=True)
        output = (result.stdout or '') + (result.stderr or '')
        print(output, end='', flush=True)
        if result.returncode == 0:
            if not package.exists():
                raise RuntimeError('npm succeeded without installing the pinned artifact client')
            marker.write_bytes(lock)
            return
        if attempt == ATTEMPTS or PERMANENT.search(output) or not TRANSIENT.search(output):
            raise subprocess.CalledProcessError(result.returncode, command)
        print(f'Transient npm download failure; retry {attempt + 1}/{ATTEMPTS}.', flush=True)
        sleep(5 * attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--action', type=Path, required=True)
    parser.add_argument('--client', type=Path, required=True)
    args = parser.parse_args()
    install(args.action, args.client)


if __name__ == '__main__':
    main()

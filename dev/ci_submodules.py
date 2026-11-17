#!/usr/bin/env python3
"""Hydrate pinned recursive submodules with bounded transient-failure recovery."""
import os
import re
import signal
import subprocess
import time


# Match the LFS network budget while retaining gitlink commits on every retry.
RETRY_DELAYS = (15, 30, 60, 120)
# An interrupted clone may have a git directory but no HEAD. Restore only
# that state, from the enclosing repository's exact gitlink, without --force.
RECOVER = ['git', 'submodule', 'foreach', '--recursive',
           'git rev-parse --verify HEAD >/dev/null 2>&1 || '
           '{ git fetch --no-tags origin "$sha1" && git checkout --detach "$sha1"; }']


def run(command, timeout):
    # git submodule invokes shells and git clone children. Terminate the whole
    # attempt on timeout so an abandoned clone cannot overlap the next retry.
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, start_new_session=True) as process:
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass  # The attempt exited as the watchdog fired.
            process.communicate()
            raise
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def hydrate():
    result = run(['git', 'submodule', 'sync', '--recursive'], timeout=60)
    print(result.stdout + result.stderr, end='', flush=True)
    result.check_returncode()
    command = ['git', 'submodule', 'update', '--init', '--recursive']
    attempts = len(RETRY_DELAYS) + 1
    interrupted = False
    for attempt in range(attempts):
        try:
            active = command
            if interrupted:
                active = RECOVER
                result = run(active, timeout=60)
                if result.returncode == 0:
                    print(result.stdout + result.stderr, end='', flush=True)
                    active = command
                    result = run(active, timeout=300)
            else:
                result = run(active, timeout=300)
        except subprocess.TimeoutExpired as error:
            interrupted = True
            print('Submodule attempt timed out; checking incomplete clones before retry.', flush=True)
            retryable = True
            failure = RuntimeError(f'Submodule command exceeded {error.timeout:g} seconds')
        else:
            output = result.stdout + result.stderr
            print(output, end='', flush=True)
            if result.returncode == 0:
                return
            permanent = re.search(r'\b(?:401|403|404|422)\b|authentication|'
                                  r'permission denied|repository not found|'
                                  r'not our ref|did not contain', output, re.I)
            retryable = not permanent and re.search(
                r'(?:HTTP\s+|returned error:\s*)(?:408|429|500|502|503|504)\b|'
                r'connection reset|connection refused|timed? out|timeout|'
                r'TLS handshake|temporary failure|could not resolve host|'
                r'lookup[^\n]*(?:server misbehaving|no such host)|unexpected EOF', output, re.I)
            failure = subprocess.CalledProcessError(result.returncode, active)
        if not retryable or attempt == attempts - 1:
            raise failure
        delay = RETRY_DELAYS[attempt]
        print(f'Submodule transient failure; retry {attempt + 2}/{attempts} in {delay}s.', flush=True)
        time.sleep(delay)


if __name__ == '__main__':
    hydrate()

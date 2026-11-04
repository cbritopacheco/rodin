#!/usr/bin/env python3
"""Hydrate the requested LFS paths with bounded transient-failure recovery."""
import argparse
import re
import subprocess
import time


def pull(include):
    command = ['git', 'lfs', 'pull', '--include=' + include]
    for attempt in range(3):
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=300)
        except subprocess.TimeoutExpired:
            retryable = True
            failure = RuntimeError('Git LFS hydration exceeded 300 seconds')
        else:
            output = result.stdout + result.stderr
            print(output, end='', flush=True)
            if result.returncode == 0:
                return
            # Permission, missing-object, quota and malformed-pointer errors
            # require intervention even if another line mentions a connection.
            permanent = re.search(r'\b(?:401|403|404|422)\b|quota|bandwidth limit|'
                                  r'authentication|permission denied|invalid pointer|'
                                  r'object does not exist', output, re.I)
            retryable = not permanent and re.search(
                r'\bHTTP\s+(?:408|429|500|502|503|504)\b|connection reset|'
                r'connection refused|timed? out|timeout|TLS handshake|'
                r'temporary failure|could not resolve host|unexpected EOF', output, re.I)
            failure = subprocess.CalledProcessError(result.returncode, command)
        if not retryable or attempt == 2:
            raise failure
        delay = 2 ** (attempt + 1)
        print(f'Git LFS transient failure; retry {attempt + 2}/3 in {delay}s.', flush=True)
        time.sleep(delay)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include', required=True)
    pull(parser.parse_args().include)

#!/usr/bin/env python3
"""Identify a Pages artifact and verify the exact snapshot served after deployment."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


VERIFY_TIMEOUT = 300  # Allow CDN propagation after deploy-pages reports success.
POLL_INTERVAL = 15
MARKER = 'rodin-publication.json'


def api(repository, resource):
    result = subprocess.run(
        ['gh', 'api', f'repos/{repository}/{resource}'],
        check=True, capture_output=True, text=True, timeout=60)
    return json.loads(result.stdout)


def prepare(repository, checkout):
    if not re.fullmatch(r'[\w.-]+/[\w.-]+', repository):
        raise ValueError('Invalid GitHub repository')
    site = api(repository, 'pages')
    if site.get('build_type') != 'workflow':
        raise ValueError('Set repository Settings > Pages > Source to GitHub Actions; '
                         'the Documentation and benchmarks environment must allow develop and master')
    commit = api(repository, 'git/ref/heads/gh-pages')['object']['sha']
    local = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    if local != commit:
        raise ValueError('Published branch changed before deployment')
    marker = {'commit': commit, 'repository': repository}
    (checkout / MARKER).write_text(json.dumps(marker) + '\n')
    output = os.environ.get('GITHUB_OUTPUT')
    if output:
        with Path(output).open('a') as stream:
            stream.write(f'commit={commit}\n')
    print(f'Prepared Pages snapshot {commit}')
    return marker


def fetch_marker(url):
    request = Request(url, headers={'Cache-Control': 'no-cache'})
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def verify(checkout, url, timeout=VERIFY_TIMEOUT):
    parsed = urlsplit(url)
    if parsed.scheme != 'https' or not parsed.netloc or parsed.query or parsed.fragment:
        raise ValueError('Expected an HTTPS Pages site URL without a query or fragment')
    expected = json.loads((checkout / MARKER).read_text())
    marker_url = url.rstrip('/') + '/' + MARKER + '?commit=' + expected['commit']
    deadline = time.monotonic() + timeout
    last_error = 'stale publication marker'
    while time.monotonic() < deadline:
        try:
            if fetch_marker(marker_url) == expected:
                print(f'Pages deployed {expected["commit"]}: {url}')
                summary = os.environ.get('GITHUB_STEP_SUMMARY')
                if summary:
                    with Path(summary).open('a') as stream:
                        stream.write(f'Pages deployed `{expected["commit"]}`: {url}\n')
                return expected['commit']
        except (HTTPError, URLError, TimeoutError, ValueError) as error:
            last_error = str(error)
        time.sleep(POLL_INTERVAL)
    raise TimeoutError(f'Pages did not serve {expected["commit"]} within {timeout} seconds: {last_error}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['prepare', 'verify'])
    parser.add_argument('--repository', default=os.environ.get('GITHUB_REPOSITORY'))
    parser.add_argument('--checkout', type=Path, required=True)
    parser.add_argument('--url')
    args = parser.parse_args()
    if args.operation == 'prepare':
        if not args.repository:
            parser.error('--repository or GITHUB_REPOSITORY is required')
        prepare(args.repository, args.checkout)
    else:
        if not args.url:
            parser.error('--url is required for verification')
        verify(args.checkout, args.url)


if __name__ == '__main__':
    main()

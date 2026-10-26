#!/usr/bin/env python3
"""Request a legacy Pages build and verify the exact published gh-pages commit."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time


BUILD_TIMEOUT = 3600  # Hosted Pages builds can queue behind another deployment.
POLL_INTERVAL = 15


def api(repository, resource, method='GET'):
    result = subprocess.run(
        ['gh', 'api', '--method', method, f'repos/{repository}/{resource}'],
        check=True, capture_output=True, text=True, timeout=60)
    return json.loads(result.stdout)


def deploy(repository, checkout=None, timeout=BUILD_TIMEOUT):
    if not re.fullmatch(r'[\w.-]+/[\w.-]+', repository):
        raise ValueError('Invalid GitHub repository')
    site = api(repository, 'pages')
    if site.get('build_type') != 'legacy' or site.get('source') != {'branch': 'gh-pages', 'path': '/'}:
        raise ValueError('Expected legacy Pages publishing from gh-pages /')
    commit = api(repository, 'git/ref/heads/gh-pages')['object']['sha']
    if checkout is not None:
        local = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
        if local != commit:
            raise ValueError('Published branch changed before deployment')
    # A GITHUB_TOKEN push alone does not trigger a legacy Pages build.
    api(repository, 'pages/builds', 'POST')
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        build = api(repository, 'pages/builds/latest')
        # The API can briefly return the previous build after accepting a request.
        if build.get('commit') == commit:
            if build.get('status') == 'built':
                print(f'Pages deployed {commit}: {site["html_url"]}')
                summary = os.environ.get('GITHUB_STEP_SUMMARY')
                if summary:
                    with Path(summary).open('a') as stream:
                        stream.write(f'Pages deployed `{commit}`: {site["html_url"]}\n')
                return commit
            if build.get('status') in ['errored', 'cancelled', 'canceled']:
                raise RuntimeError(f'Pages deployment failed for {commit}: {build.get("error")}')
        time.sleep(POLL_INTERVAL)
    raise TimeoutError(f'Pages did not confirm deployment of {commit} within {timeout} seconds')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', default=os.environ.get('GITHUB_REPOSITORY'))
    parser.add_argument('--checkout', type=Path)
    args = parser.parse_args()
    if not args.repository:
        parser.error('--repository or GITHUB_REPOSITORY is required')
    deploy(args.repository, args.checkout)


if __name__ == '__main__':
    main()

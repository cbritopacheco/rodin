#!/usr/bin/env python3
"""Publish Shields endpoint data from current-commit Actions results."""

import argparse
import base64
import json
import os
import subprocess
from urllib.parse import quote

COLORS = {'success': '2e8540', 'failure': 'c62828', 'timed_out': 'c62828',
          'queued': '2874c7', 'running': 'd49b16', 'cancelled': 'e87924',
          'skipped': 'e87924', 'unknown': '808080'}
WORKFLOWS = {'Build': 'Build.yml', 'Tests': 'Tests.yml', 'Benchmarks': 'Benchmarks.yml'}


def state(status, conclusion):
    if status != 'completed':
        return 'running' if status == 'in_progress' else 'queued' if status in {'queued', 'requested', 'pending', 'waiting'} else 'unknown'
    return conclusion if conclusion in COLORS else 'failure' if conclusion in {'startup_failure', 'action_required'} else 'unknown'


def badge(label, status, commit):
    return {'schemaVersion': 1, 'label': label, 'message': status.replace('_', ' '),
            'color': COLORS[status], 'labelColor': '30363d', 'style': 'flat', 'commit': commit}


class Publisher:
    """Read current source runs, then atomically update isolated badge files.

    The badge branch holds JSON only. Before every write, the source branch is
    checked again. One workflow-level lock serializes publishers; file SHA
    checks also prevent silent replacement of concurrent changes.
    """

    def __init__(self, repository, branch, request=None):
        if branch not in {'develop', 'master'}:
            raise ValueError('Expected develop or master')
        self.repository, self.branch = repository, branch
        self.prefix = 'repos/' + repository
        self.request = request or self.github

    @staticmethod
    def github(endpoint, method='GET', data=None):
        command = ['gh', 'api', '--method', method, endpoint]
        if data is not None:
            command += ['--input', '-']
        result = subprocess.run(command, input=json.dumps(data) if data is not None else None,
                                text=True, capture_output=True, check=True)
        return json.loads(result.stdout) if result.stdout.strip() else None

    def head(self):
        return self.request(f'{self.prefix}/git/ref/heads/{self.branch}')['object']['sha']

    def branch_ready(self):
        try:
            self.request(self.prefix + '/git/ref/heads/ci-badges')
        except subprocess.CalledProcessError as error:
            if '404' not in error.stderr:
                raise
            tree = self.request(self.prefix + '/git/trees', 'POST', {'tree': [
                {'path': 'README.md', 'mode': '100644', 'type': 'blob',
                 'content': 'Generated CI status badge data.\n'}]})
            commit = self.request(self.prefix + '/git/commits', 'POST',
                                  {'message': 'Initialize CI badge data', 'tree': tree['sha'], 'parents': []})
            self.request(self.prefix + '/git/refs', 'POST', {'ref': 'refs/heads/ci-badges', 'sha': commit['sha']})

    def collect(self, sha, documentation=None, source_run=None, source_sha=None):
        result = {}
        benchmark_run = None
        for label, file in WORKFLOWS.items():
            runs = self.request(f'{self.prefix}/actions/workflows/{file}/runs?branch={self.branch}&event=push&per_page=1')['workflow_runs']
            run = runs[0] if runs else None
            status = state(run['status'], run.get('conclusion')) if run and run['head_sha'] == sha else 'unknown'
            result[label] = badge(label, status, sha)
            if label == 'Benchmarks':
                benchmark_run = run if run and run['head_sha'] == sha else None
        doc_status = 'unknown'
        if benchmark_run:
            # Documentation publication is independent of the long benchmark campaign.
            page = 1
            while True:
                jobs = self.request(f"{self.prefix}/actions/runs/{benchmark_run['id']}/jobs?filter=latest&per_page=100&page={page}")['jobs']
                found = next((j for j in jobs if j['name'] == 'Publish documentation'), None)
                if found:
                    doc_status = state(found['status'], found.get('conclusion'))
                    break
                if len(jobs) < 100:
                    if benchmark_run['status'] != 'completed':
                        doc_status = 'queued'
                    break
                page += 1
            # A publication step can report deployment immediately, before its
            # enclosing job or benchmark campaign has completed.
            if (documentation in {'success', 'failure'}
                    and str(benchmark_run['id']) == source_run
                    and sha == source_sha):
                doc_status = documentation
        result['Documentation'] = badge('Documentation', doc_status, sha)
        return result

    def publish(self, documentation=None, source_run=None, source_sha=None):
        sha = self.head()
        records = self.collect(sha, documentation, source_run, source_sha)
        self.branch_ready()
        for label, record in records.items():
            if self.head() != sha:
                print('Source branch advanced; leaving newer events to refresh badges.')
                return
            path = f'{self.prefix}/contents/{self.branch}/{label.lower()}.json'
            previous = None
            try:
                previous = self.request(path + '?ref=ci-badges')
            except subprocess.CalledProcessError as error:
                if '404' not in error.stderr:
                    raise
            content = json.dumps(record, sort_keys=True) + '\n'
            if previous and base64.b64decode(previous['content']).decode() == content:
                continue
            payload = {'message': f'Update {self.branch} {label} badge', 'branch': 'ci-badges',
                       'content': base64.b64encode(content.encode()).decode()}
            if previous:
                payload['sha'] = previous['sha']
            self.request(path, 'PUT', payload)
            print(label + ': ' + record['message'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch', choices=['develop', 'master'], required=True)
    parser.add_argument('--documentation', choices=['success', 'failure'])
    parser.add_argument('--source-run')
    parser.add_argument('--source-sha')
    args = parser.parse_args()
    Publisher(os.environ['GITHUB_REPOSITORY'], args.branch).publish(args.documentation, args.source_run, args.source_sha)


if __name__ == '__main__':
    main()

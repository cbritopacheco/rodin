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
    return {'schemaVersion': 1, 'label': label, 'message': {'success': 'passing', 'failure': 'failing'}.get(status, status.replace('_', ' ')),
            'color': COLORS[status], 'labelColor': '555', 'style': 'flat', 'commit': commit}


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

    def jobs(self, run):
        key = run['id']
        if key not in self.job_cache:
            result, page = [], 1
            try:
                while True:
                    rows = self.request(f"{self.prefix}/actions/runs/{key}/jobs?filter=latest&per_page=100&page={page}")['jobs']
                    result.extend(rows)
                    if len(rows) < 100:
                        break
                    page += 1
                self.job_cache[key] = result
            except subprocess.CalledProcessError:
                self.job_cache[key] = self.checks(run)
                if self.job_cache[key] is None:
                    print('Job lookup unavailable; publishing unknown for affected unfinished statuses.')
        return self.job_cache[key]

    def checks(self, run):
        # Large matrix runs can fail the REST job-list endpoint. Check results
        # provide an independent read path for status and publication names.
        if not run.get('check_suite_node_id'):
            return None
        latest, cursor = {}, None
        query = ('query($suite:ID!,$cursor:String){node(id:$suite){... on CheckSuite{'
                 'checkRuns(first:100,after:$cursor){nodes{databaseId name status conclusion startedAt}'
                 'pageInfo{hasNextPage endCursor}}}}}')
        try:
            while True:
                response = self.request('graphql', 'POST', {'query': query, 'variables':
                                        {'suite': run['check_suite_node_id'], 'cursor': cursor}})
                page = response['data']['node']['checkRuns']
                for job in page['nodes']:
                    if job['name'] not in latest or job['databaseId'] > latest[job['name']]['databaseId']:
                        latest[job['name']] = job
                if not page['pageInfo']['hasNextPage']:
                    break
                cursor = page['pageInfo']['endCursor']
            return [dict(id=j['databaseId'], name=j['name'], started_at=j.get('startedAt'), status=j['status'].lower(),
                         conclusion=(j['conclusion'] or '').lower()) for j in latest.values()]
        except (subprocess.CalledProcessError, KeyError, TypeError):
            return None

    def collect(self, sha, documentation=None, source_run=None, source_sha=None):
        result = {}
        self.job_cache = {}
        for label, file in WORKFLOWS.items():
            runs = self.request(f'{self.prefix}/actions/workflows/{file}/runs?branch={self.branch}&event=push&per_page=1')['workflow_runs']
            run = runs[0] if runs else None
            status = state(run['status'], run.get('conclusion')) if run and run['head_sha'] == sha else 'unknown'
            if run and run['head_sha'] == sha and run['status'] != 'completed':
                jobs = self.jobs(run)
                if jobs is None:
                    status = 'unknown'
                elif any(j['status'] == 'in_progress' for j in jobs):
                    status = 'running'
                elif any(j['status'] in {'queued', 'waiting', 'pending', 'requested'} for j in jobs):
                    status = 'queued'
            result[label] = badge(label, status, sha)
        doc_status = 'unknown'
        doc_runs = self.request(f'{self.prefix}/actions/workflows/Documentation.yml/runs?branch={self.branch}&event=push&per_page=1')['workflow_runs']
        doc_run = doc_runs[0] if doc_runs and doc_runs[0]['head_sha'] == sha else None
        if doc_run:
            # Only the documentation workflow's deployment establishes publication.
            jobs = self.jobs(doc_run)
            if jobs is not None:
                found = next((j for j in jobs if j['name'] == 'Publish documentation'), None)
                if any(j.get('conclusion') in {'failure', 'timed_out', 'startup_failure', 'action_required'} for j in jobs):
                    doc_status = 'failure'
                elif found and found.get('conclusion') == 'success':
                    doc_status = 'success'
                elif any(j['status'] == 'in_progress' for j in jobs):
                    doc_status = 'running'
                elif doc_run['status'] == 'completed':
                    doc_status = state(doc_run['status'], doc_run.get('conclusion'))
                    if doc_status == 'success':
                        doc_status = 'unknown'
                else:
                    doc_status = 'queued'
            # A publication step can report deployment immediately, before its
            # enclosing publication job has completed.
            if (documentation in {'success', 'failure'}
                    and str(doc_run['id']) == source_run
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

#!/usr/bin/env python3
"""One compact, read-only report for a completed workflow and its latest jobs."""

from collections import Counter
import io
import json
import os
from pathlib import Path
import subprocess
import zipfile

from ci_summary import append, cell


def api(endpoint):
    result = subprocess.run(['gh', 'api', '--paginate', '--slurp', endpoint],
                            check=True, capture_output=True, text=True)
    return json.loads(result.stdout)


def render(run, jobs, repository, coverage=None):
    base = f"{os.environ.get('GITHUB_SERVER_URL', 'https://github.com')}/{repository}/actions/runs/{run['id']}"
    counts = Counter(j.get('conclusion') or j['status'] for j in jobs)
    ran = sum(bool(j.get('started_at')) and j.get('conclusion') != 'skipped' for j in jobs)
    rows = [f"## {cell(run['name'])}", '',
            f"Commit: {cell(run['head_sha'])} · Attempt {run['run_attempt']} · **{cell(run['conclusion'])}**", '',
            '### Policies', '',
            '- **Automatic retry:** After workflow completion, at most one interrupted job is selected per attempt. Each job can be retried once, only on a current develop or master push commit. GitHub also reruns its dependents.',
            '- **Retry exclusions:** Compilation, numerical, sanitizer, timeout, checksum and out-of-memory failures are not retried. Intentional workflow cancellations and superseded commits are excluded. Missing logs prevent recovery.',
            ('- **Completeness:** Every planned benchmark case and MPI rank contribution is required before pooling and publication.' if run['name'] == 'Benchmarks' else '- **Completeness:** Every mandatory test and coverage contribution is required. Coverage / Finalize must pass the complete production line coverage gate of at least 95%.' if run['name'] == 'Tests' else '- **Publication:** The site link is shown only after successful publication.'), '',
            '### Job results', '',
            '| Registered | Ran | Passed | Failed / timed out | Cancelled | Skipped | Other / pending |',
            '| --- | --- | --- | --- | --- | --- | --- |',
            f"| {len(jobs)} | {ran} | {counts['success']} | {counts['failure'] + counts['timed_out']} | {counts['cancelled']} | {counts['skipped']} | {len(jobs) - sum(counts[k] for k in ('success', 'failure', 'timed_out', 'cancelled', 'skipped'))} |", '',
            'Registered jobs include expanded matrices only. Ran counts jobs with a start time, excluding skipped jobs.', '',
            f'[Workflow, attempts and artifacts]({base})', '']
    failures = [j for j in jobs if j.get('conclusion') in {'failure', 'timed_out', 'cancelled'}]
    if failures:
        rows += ['| Failed job / logs | Failed step |', '| --- | --- |']
        for job in failures[:30]:
            steps = [s['name'] for s in job.get('steps', []) if s.get('conclusion') in {'failure', 'timed_out', 'cancelled'}]
            rows.append(f"| [{cell(job['name'])}]({base}/job/{job['id']}) | {cell('; '.join(steps) or 'See job log')} |")
        if len(failures) > 30:
            rows += ['', f'{len(failures) - 30} additional failures are listed in the workflow.']
    if run['name'] == 'Tests':
        gate = next((j for j in jobs if j['name'] == 'Coverage / Finalize'), None)
        status = gate.get('conclusion') if gate else 'unavailable'
        rows += ['', f'**Coverage / Finalize:** {cell(status)}.']
        if coverage:
            hit, lines = coverage['hit'], coverage['lines']
            rows += [f'Combined production lines: **{hit}/{lines} ({100 * hit / lines:.2f}%)**.']
        else:
            rows += ['Complete combined coverage report unavailable; no partial percentage is substituted.']
        rows += ['The finalization job must pass to confirm all contributions, publication receipts and the 95% gate.']
    if run['name'] == 'Benchmarks':
        lanes = [j for j in jobs if '/NP' in j['name']]
        passed = sum(j.get('conclusion') == 'success' for j in lanes)
        pool = next((j for j in jobs if j['name'] == 'Validate complete campaign'), None)
        rows += ['', f'Execution lanes: **{passed}/{len(lanes)} registered lanes passed**. '
                 'Unexpanded matrices are not counted as completed work.',
                 f"**Complete campaign validation / pooling:** {cell(pool.get('conclusion') if pool else 'unavailable')}."]
    publication = [j for j in jobs if j['name'] in {'Publish documentation', 'Publish documentation and benchmarks'}]
    for job in publication:
        rows += ['', f"**{cell(job['name'])}:** {cell(job.get('conclusion') or job['status'])}."]
        if job.get('conclusion') == 'success' and run['head_branch'] in {'develop', 'master'}:
            rows += [f"[Published site](https://cbritopacheco.github.io/rodin/docs/refs/heads/{run['head_branch']}/index.html)"]
    return '\n'.join(rows)


def main():
    event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text())
    repository = os.environ['GITHUB_REPOSITORY']
    run_id = int(event['workflow_run']['id'])
    prefix = f'repos/{repository}/actions/runs/{run_id}'
    run = api(prefix)[0]
    jobs = [j for page in api(prefix + '/jobs?filter=latest&per_page=100') for j in page['jobs']]
    coverage = None
    if run['name'] == 'Tests':
        artifacts = [a for page in api(prefix + '/artifacts?per_page=100') for a in page['artifacts']
                     if a['name'] == 'coverage-final-summary' and not a['expired']
                     and a['created_at'] >= run['run_started_at']]
        if len(artifacts) == 1:
            try:
                result = subprocess.run(['gh', 'api', f"repos/{repository}/actions/artifacts/{artifacts[0]['id']}/zip"],
                                        check=True, capture_output=True)
                with zipfile.ZipFile(io.BytesIO(result.stdout)) as archive:
                    info = archive.getinfo('coverage-summary.json')
                    if info.file_size > 1024 * 1024:
                        raise ValueError('Unexpected coverage summary size')
                    data = json.loads(archive.read(info))
                if (type(data['hit']) is not int or type(data['lines']) is not int
                        or not 0 <= data['hit'] <= data['lines'] or data['lines'] == 0):
                    raise ValueError('Invalid coverage counts')
                coverage = data
            except (subprocess.CalledProcessError, ValueError, KeyError, zipfile.BadZipFile):
                print('Coverage summary unavailable or invalid; consult finalization logs.')
    append(render(run, jobs, repository, coverage))


if __name__ == '__main__':
    main()

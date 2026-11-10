#!/usr/bin/env python3
"""Render bounded diagnostic Markdown for native GitHub Actions summaries."""

import html
import json
import os
from pathlib import Path


def cell(value):
    """Keep metadata literal inside a Markdown table cell."""
    return html.escape(str(value), quote=True).replace('|', '&#124;').replace('\n', ' ').replace('\r', ' ')


def append(markdown):
    if path := os.environ.get('GITHUB_STEP_SUMMARY'):
        with Path(path).open('a') as stream:
            stream.write(markdown + '\n')


def benchmark():
    """Record failed-lane context, including failures before execution begins."""
    lane = json.loads(os.environ['BENCHMARK_LANE'])
    steps = json.loads(os.environ['BENCHMARK_STEPS'])
    entry = lane.get('entry', {})
    diagnostics = Path('benchmark-diagnostics')
    diagnostics.mkdir(exist_ok=True)
    record = {'commit': os.environ['GITHUB_SHA'], 'attempt': os.environ['GITHUB_RUN_ATTEMPT'],
              'run_id': os.environ['GITHUB_RUN_ID'], 'lane': lane, 'steps': steps}
    (diagnostics / 'context.json').write_text(json.dumps(record, indent=2) + '\n')
    if os.environ.get('BENCHMARK_CONTEXT_ONLY') == '1':
        return
    # Successful lanes remain quiet; hundreds of success cards obscure failures.
    failed = [name for name, step in steps.items() if step.get('outcome') in {'failure', 'cancelled'}]
    if not failed:
        return
    run = (f"{os.environ.get('GITHUB_SERVER_URL', 'https://github.com')}/"
           f"{os.environ['GITHUB_REPOSITORY']}/actions/runs/{os.environ['GITHUB_RUN_ID']}")
    rows = ['## Benchmark lane failure', '', '| Field | Value |', '| --- | --- |',
            f"| Lane | {cell(lane.get('label', lane['id']))} |",
            f"| Commit / attempt | {cell(record['commit'])} / {cell(record['attempt'])} |",
            f"| Executable | {cell(entry.get('target', 'Unknown'))} |",
            f"| Configuration / MPI ranks | {cell(lane.get('configuration', 'Unknown'))} / {cell(lane.get('ranks', 'Unknown'))} |",
            f"| Selected cases | {len(lane.get('cases', []))} |",
            f"| Failed steps | {cell(', '.join(failed))} |", '',
            f'[Workflow and job logs]({run}/attempts/{record["attempt"]})', '',
            'Retry: pending inspection after workflow completion; no retry has been requested by this lane.', '',
            'The **Recover interrupted CI jobs** workflow records eligibility and retry decisions. '
            'An interruption alone does not establish its cause.']
    resource_file = diagnostics / 'resources.txt'
    if resource_file.exists():
        resource_rows = []
        for line in resource_file.read_text().splitlines():
            for label in ('Elapsed (wall clock) time', 'Maximum resident set size (kbytes)', 'Exit status'):
                if line.strip().startswith(label):
                    resource_rows.append(cell(line.strip()))
        if resource_rows:
            rows += ['', 'Lane process statistics (including build restore):', '',
                     *['- ' + line for line in resource_rows]]
    rows += ['', 'Partial reports are diagnostic only; a validated receipt remains mandatory for pooling.', '',
             '<details><summary>Selected benchmark cases</summary>', '',
             *[f'- {cell(case)}' for case in lane.get('cases', [])[:50]], '', '</details>']
    if len(lane.get('cases', [])) > 50:
        rows += ['', 'The complete case list is saved in context.json.']
    append('\n'.join(rows))


if __name__ == '__main__':
    benchmark()

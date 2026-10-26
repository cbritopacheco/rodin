#!/usr/bin/env python3
"""Publish documentation independently and update only complete benchmark campaigns."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from html import escape
import hashlib
import json
from pathlib import Path
import os
import re
import shutil
import subprocess

from benchmark_site import publication_path
from benchmark_ci import dump
from ci_builds import digest


# Bound Pages storage while retaining a useful recent trend. Complete raw reports
# remain workflow artifacts; retiring historical requirements permits this window.
HISTORY_COMMITS = 20


def link_histories(destination, history_root):
    """Bind generated component pages to the latest published complete histories."""
    for page in (destination / 'benchmarks').glob('**/index.html'):
        component = page.parent.relative_to(destination / 'benchmarks')
        index = history_root / component / 'index.json'
        if not index.is_file():
            continue
        url = os.path.relpath(index, page.parent).replace(os.sep, '/')
        markup = page.read_text()
        # Rebinding is idempotent, including pages already published by a campaign.
        markup = re.sub(r' data-history-index="[^"]*"', '', markup)
        page.write_text(markup.replace('data-benchmark-history',
                        'data-benchmark-history data-history-index="' + escape(url, quote=True) + '"', 1))


def publish_documentation(site, pages, ref):
    """Replace documentation independently, retaining all published benchmark data."""
    publication_path(ref)
    if not (site / 'benchmarks/index.html').is_file():
        raise ValueError('Documentation build has no benchmark section')
    destination = pages / 'docs' / ref
    destination.resolve().relative_to(pages.resolve())
    if (site.resolve().is_relative_to(destination.resolve()) or
        destination.resolve().is_relative_to(site.resolve())):
        raise ValueError('Documentation source and destination must not overlap')
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(site, destination)
    link_histories(destination, pages / 'benchmark-data' / ref)
    (pages / '.nojekyll').touch()


def publish(results, site, pages, ref, run):
    publication_path(ref)
    complete = json.loads((results / 'complete.json').read_text())
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if complete.get('schema') != 1 or complete.get('commit') != commit:
        raise ValueError('Only the complete current-commit campaign can be published')
    if complete.get('configurations') != ['distributed', 'local']:
        raise ValueError('Incomplete publication configuration set')
    history_root = pages / 'benchmark-data' / ref
    destination = pages / 'docs' / ref
    if destination.exists():
        site = destination
    if not (site / 'benchmarks/index.html').is_file():
        raise ValueError('Documentation build has no benchmark section')
    components = {}
    files = sorted(results.glob('**/results.json'))
    if not files:
        raise ValueError('Complete campaign has no component results')
    if {str(path.relative_to(results)): digest(path) for path in files} != complete.get('components'):
        raise ValueError('Incomplete or changed pooled result package')
    for path in files:
        data = json.loads(path.read_text())
        component = str(publication_path(data['component']))
        configuration = str(publication_path(data['configuration']))
        if '/' in configuration or data['commit'] != commit or not data['results']:
            raise ValueError('Invalid component publication')
        if not (site / 'benchmarks' / component / 'index.html').is_file():
            raise ValueError('Measured component is absent from documentation: ' + component)
        components.setdefault(component, {})[configuration] = data
    # Validate the entire campaign before changing public histories.
    updates = []
    now = datetime.now(timezone.utc).isoformat()
    for component, configurations in components.items():
        index = []
        for configuration, data in sorted(configurations.items()):
            cases = []
            for value in data['results']:
                case_id = hashlib.sha256((value['name'] + '\0' + value['metric']).encode()).hexdigest()
                target = history_root / component / configuration / 'cases' / (case_id + '.json')
                history = json.loads(target.read_text()) if target.exists() else {
                    'schema': 1, 'component': component, 'configuration': configuration,
                    'name': value['name'], 'metric': value['metric'], 'entries': []}
                if (history.get('schema') != 1 or history.get('component') != component or
                    history.get('configuration') != configuration or history.get('name') != value['name'] or
                    history.get('metric') != value['metric']):
                    raise ValueError('Incompatible case history')
                history['format'] = data['format']
                history['window'] = HISTORY_COMMITS
                history['entries'] = [entry for entry in history['entries'] if entry['commit'] != commit]
                history['entries'].append({'commit': commit, 'date': now, 'run': run, 'results': [value]})
                history['entries'] = history['entries'][-HISTORY_COMMITS:]
                cases.append({'id': case_id, 'name': value['name'], 'metric': value['metric'],
                              'url': 'cases/' + case_id + '.json'})
                updates.append((target, history))
            updates.append((history_root / component / configuration / 'index.json', {
                'schema': 1, 'definition': data['definition'], 'ranks': data['ranks'], 'cases': cases}))
            label = data['definition']['compiler'] + ' / ' + data['definition']['cmake']['CMAKE_BUILD_TYPE'] + ' / ' + data['definition']['mode'] + ' / NP' + str(data['ranks'])
            index.append({'id': configuration, 'label': label, 'url': configuration + '/index.json'})
        updates.append((history_root / component / 'index.json', {'schema': 1, 'configurations': index}))
    # Check every destination before modifying the checkout or replacing docs.
    for path in [destination, pages / 'benchmarks', pages / '.nojekyll', history_root / 'complete.json', *[path for path, _ in updates]]:
        path.resolve().relative_to(pages.resolve())
    if not destination.exists():
        publish_documentation(site, pages, ref)
    for path, value in updates:
        path.resolve().relative_to(pages.resolve())
        dump(path, value)
    link_histories(destination, history_root)
    dump(history_root / 'complete.json', complete)
    (pages / '.nojekyll').touch()
    # The user authorized retiring the generated historical dashboards, not source cases.
    historical = pages / 'benchmarks'
    historical.resolve().relative_to(pages.resolve())
    if historical.exists():
        shutil.rmtree(historical)
    return len(components)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--documentation-only', action='store_true')
    parser.add_argument('--results', type=Path)
    for name in ['site', 'pages']:
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('--run')
    args = parser.parse_args()
    if args.documentation_only:
        publish_documentation(args.site, args.pages, args.ref)
        print('Published documentation with the latest complete benchmark histories.')
    else:
        if args.results is None or args.run is None:
            parser.error('Benchmark publication requires --results and --run')
        print('Published', publish(args.results, args.site, args.pages, args.ref, args.run), 'benchmark components.')


if __name__ == '__main__':
    main()

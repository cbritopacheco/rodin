#!/usr/bin/env python3
"""Derive exact mandatory artifact names from the existing CI plans.

No workflow-wide listing or glob discovery is needed, including after a run
passes 1000 artifacts. Missing artifacts, receipts and reports remain failures.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import ci_build_graph as graph
import ci_builds as archive
import ci_coverage as coverage
from ci_coverage_publish import publication_plan


def common_id(job, matrix):
    return archive.identity(graph.common_specification(graph.specification(job, matrix)))


def ordinary_name(family, definition):
    return 'ci-result-' + family + '-' + common_id('Ordinary', definition) + '-' + definition['id']


def target_name(definition):
    return ('ci-target-' + definition['job'] + '-' + common_id(**definition) + '-' +
            coverage.target_identity(definition)[:24])


def contribution_names(ordinary, targeted):
    return ([ordinary_name(family, row) for family, rows in ordinary.items() for row in rows] +
            [target_name(row) for row in targeted])


def publication_names(kind, family, common, ordinary, targeted):
    ordinary, targeted = publication_plan(kind, family, ordinary, targeted, common)
    return contribution_names(ordinary, targeted)


def final_names(ordinary, targeted, common, index):
    names = contribution_names(ordinary, targeted)
    names += ['ci-built-' + key for key in index['builds']]
    families = {row['id']: family for family, rows in ordinary.items() for row in rows}
    for entry in common:
        if entry['definition']['cmake'].get('RODIN_CODE_COVERAGE') != 'ON':
            continue
        selected_families = {families[r['matrix']['id']] if r['job'] == 'Ordinary' else r['job']
                             for r in entry['requests']}
        names += ['ci-publication-' + family + '-' + entry['id'] for family in sorted(selected_families)]
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['publication', 'indices', 'final'])
    parser.add_argument('--kind', choices=['ordinary', 'targeted'])
    parser.add_argument('--family')
    parser.add_argument('--common')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    common = json.loads(Path('ci-common-plan.json').read_text())
    if args.operation == 'indices':
        names = ['ci-build-index-' + entry['id'] for entry in common]
    else:
        ordinary = json.loads(Path('ci-test-plan.json').read_text())
        targeted = json.loads(Path('ci-coverage-targets.json').read_text())
        if args.operation == 'publication':
            if not all([args.kind, args.family, args.common]):
                parser.error('Publication requires kind, family and common')
            names = publication_names(args.kind, args.family, args.common, ordinary, targeted)
        else:
            index = graph.load_indices(Path('ci-build-indices'), os.environ['GITHUB_SHA'], common)
            names = final_names(ordinary, targeted, common, index)
    if not names or len(names) != len(set(names)):
        raise ValueError('Expected nonempty, unique mandatory artifact names')
    args.output.write_text(json.dumps(sorted(names), indent=2) + '\n')
    print('Required exact artifacts:', len(names))


if __name__ == '__main__':
    main()

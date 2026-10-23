#!/usr/bin/env python3
"""Publish bounded incremental unions without dropping any lane contribution.

Lane receipts and hashes remain independent. A publication is accepted by the
final gate only after its uploader succeeded and every expected report belongs
to exactly one batch. The final gate still recomputes the full raw LCOV union.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import ci_coverage as coverage


def write_lcov(lines, path):
    sources = {}
    for (source, line), hit in lines.items():
        sources.setdefault(source, {})[line] = hit
    with path.open('w') as stream:
        for source, records in sorted(sources.items()):
            stream.write('TN:\nSF:' + source + '\n')
            for line, hit in sorted(records.items()):
                stream.write(f'DA:{line},{int(hit)}\n')
            stream.write(f'LF:{len(records)}\nLH:{sum(records.values())}\nend_of_record\n')


def collect(kind, family, folder, output, commit, ordinary_plan, target_plan):
    if kind == 'ordinary':
        receipts = coverage.verify_receipts({family: ordinary_plan[family]}, folder, commit, check_inventory=False)
        expected = {'ordinary:' + key: r['report_sha256'] for key, r in receipts.items() if r['coverage_expected']}
    else:
        definitions = [d for d in target_plan if d['job'] == family]
        if not definitions:
            raise ValueError('Unknown targeted publication family')
        receipts = coverage.verify_target_receipts(definitions, folder, commit)
        expected = {'target:' + key: r['report_sha256'] for key, r in receipts.items()}
    reports = list(folder.rglob('RodinCoverage.info'))
    if len(reports) != len(expected):
        raise ValueError('Unexpected or missing lane coverage files')
    output.mkdir(parents=True)
    report = output / 'CombinedCoverage.info'
    lines = coverage.production_lines(reports)
    write_lcov(lines, report)
    if coverage.production_lines([report]) != lines:
        raise ValueError('Publication changed production coverage')
    receipt = {'commit': commit, 'kind': kind, 'family': family, 'uploaded': False,
               'contributions': expected, 'report_sha256': hashlib.sha256(report.read_bytes()).hexdigest()}
    (output / 'publication.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(f'Validated {len(expected)} independent contributions; {sum(lines.values())}/{len(lines)} production lines')


def verify_publications(ordinary, targeted, folder, commit):
    expected = {'ordinary:' + key: r['report_sha256'] for key, r in ordinary.items() if r['coverage_expected']}
    expected.update({'target:' + key: r['report_sha256'] for key, r in targeted.items()})
    received = {}
    paths = list(folder.rglob('publication.json'))
    for path in paths:
        receipt = json.loads(path.read_text())
        report = path.parent / 'CombinedCoverage.info'
        if receipt['commit'] != commit or not receipt['uploaded']:
            raise ValueError('Stale or unsuccessful Codecov publication')
        if not report.is_file() or hashlib.sha256(report.read_bytes()).hexdigest() != receipt['report_sha256']:
            raise ValueError('Missing or changed published LCOV union')
        for key, digest in receipt['contributions'].items():
            if key in received or expected.get(key) != digest:
                raise ValueError('Unexpected, duplicate or changed published contribution: ' + key)
            received[key] = digest
    if received != expected:
        raise ValueError('Missing Codecov contributions: ' + str(len(expected.keys() - received.keys())))
    if coverage.production_lines([p.parent / 'CombinedCoverage.info' for p in paths]) != coverage.production_lines(
            list(folder.rglob('RodinCoverage.info'))):
        raise ValueError('Published coverage does not equal the complete raw lane union')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['collect', 'record'])
    parser.add_argument('--kind', choices=['ordinary', 'targeted'])
    parser.add_argument('--family')
    parser.add_argument('--folder', type=Path, default=Path('ci-contributions'))
    parser.add_argument('--output', type=Path, default=Path('ci-publication'))
    args = parser.parse_args()
    if args.operation == 'collect':
        collect(args.kind, args.family, args.folder, args.output, os.environ['GITHUB_SHA'],
                json.loads(Path('ci-test-plan.json').read_text()),
                json.loads(Path('ci-coverage-targets.json').read_text()))
    else:
        receipt = args.output / 'publication.json'
        data = json.loads(receipt.read_text())
        if os.environ.get('COVERAGE_UPLOADED') != 'success':
            raise ValueError('Codecov upload did not succeed')
        data['uploaded'] = True
        receipt.write_text(json.dumps(data, indent=2) + '\n')


if __name__ == '__main__':
    main()

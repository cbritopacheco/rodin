#!/usr/bin/env python3
"""Build benchmark navigation from configured targets and rendered API overviews."""
from __future__ import annotations

import argparse
import csv
from html import escape
import json
from pathlib import Path, PurePosixPath
import re
import shutil
from urllib.parse import urljoin, urlparse

from jinja2 import Environment, FileSystemLoader, select_autoescape

from benchmark_search import build_search


ROOT = Path(__file__).resolve().parents[1]


def publication_path(value):
    path = PurePosixPath(value)
    if not value or path.is_absolute() or any(p in ('', '.', '..') for p in value.split('/')):
        raise ValueError('Invalid benchmark publication path: ' + value)
    if any(not re.fullmatch(r'[A-Za-z0-9_-]+', p) for p in path.parts):
        raise ValueError('Invalid benchmark publication path: ' + value)
    return path


def inventory(manifest):
    with manifest.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        required = {'target', 'owner', 'executable', 'state', 'path', 'description', 'format', 'profile', 'metric', 'ranks', 'args'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Incomplete benchmark manifest header')
        records = list(reader)
    targets, paths = set(), set()
    for row in records:
        if not row.get('target') or row['target'] in targets:
            raise ValueError('Missing or duplicate benchmark target')
        targets.add(row['target'])
        if row['state'] not in ('registered', 'unregistered', 'library'):
            raise ValueError('Unknown benchmark registration state')
        if row['state'] != 'registered':
            continue
        path = publication_path(row['path'])
        if str(path) in paths:
            raise ValueError('Duplicate benchmark publication path: ' + str(path))
        paths.add(str(path))
        if row['metric'] not in ('real_time', 'cpu_time', 'auto'):
            raise ValueError('Unknown benchmark timing metric: ' + row['metric'])
        if row['format'] not in ('googlecpp', 'csv'):
            raise ValueError('Unknown benchmark report format')
        ranks = row['ranks'].split(';')
        if not ranks or any(not rank.isdigit() or int(rank) < 1 for rank in ranks):
            raise ValueError('Benchmark ranks must be positive integers')
        if not Path(row['description']).is_file():
            raise ValueError('Missing benchmark description: ' + row['description'])
    return records


def rebase_links(markup, documentation_url):
    """Preserve links from API overviews to the canonical API site."""
    if not documentation_url.startswith(('https://', 'http://', '/', '../', '__RODIN_DOCS__/')):
        raise ValueError('Documentation URL must be absolute')
    return re.sub(r'\bhref=([\'"])([^\'"]+)\1',
                  lambda m: 'href=' + m[1] + (m[2] if m[2].startswith('#') else
                                              urljoin(documentation_url, m[2])) + m[1], markup)


def benchmark_navigation(documentation, documentation_url):
    """Reuse the documentation navigation, replacing only its search glyphs."""
    header = re.search(r'<header\b[^>]*>.*?</header>', documentation, re.S)
    if not header:
        raise ValueError('Documentation navigation is missing')
    markup = rebase_links(header[0], documentation_url)
    icon = ('<svg class="rodin-benchmark-icon" viewBox="0 0 24 24" aria-hidden="true">'
            '<path d="M4 3v17h17 M7 15l4-5 4 3 6-7"></path></svg>')
    def search_control(match):
        anchor = match[0].replace('m-doc-search-icon', 'm-doc-search-icon benchmark-search-trigger', 1)
        anchor = re.sub(r'title="[^"]*"', 'title="Search benchmarks"', anchor)
        anchor = anchor.replace('<svg ', '<svg class="benchmark-magnifier-icon" ', 1)
        anchor = anchor.replace('>', ' aria-label="Search benchmarks">', 1)
        return anchor.replace('</a>', icon + '</a>')
    return re.sub(r'<a\b[^>]*class="[^"]*m-doc-search-icon[^"]*"[^>]*>.*?</a>',
                  search_control, markup, flags=re.S)


def build_site(manifest, documentation, output, configuration, documentation_url):
    records = inventory(manifest)
    catalogue = json.loads((documentation / 'rodin-module-catalogue.json').read_text())
    if catalogue.get('schema') != 1:
        raise ValueError('Unsupported documentation module catalogue')
    footer = re.search(r'<footer\b[^>]*>.*?</footer>',
                       (documentation / 'index.html').read_text(), re.S)
    if not footer:
        raise ValueError('Documentation footer is missing')
    footer = rebase_links(footer[0], documentation_url if not documentation_url.startswith('../')
                         else '__RODIN_DOCS__/')
    navigation = benchmark_navigation((documentation / 'index.html').read_text(),
                                     documentation_url if not documentation_url.startswith('../')
                                     else '__RODIN_DOCS__/')
    output.mkdir(parents=True, exist_ok=True)
    assets = output / 'assets'
    assets.mkdir(exist_ok=True)
    stylesheets = []
    for stylesheet in catalogue['stylesheets']:
        if urlparse(stylesheet).netloc:
            stylesheets.append(stylesheet)
        else:
            name = Path(stylesheet).name
            shutil.copyfile(documentation / name, assets / name)
            stylesheets.append('assets/' + name)
    for name in ['rodin.css', 'rodin-theme.css', 'rodin-theme.js', 'rodin-navigation.js']:
        shutil.copyfile(ROOT / 'doc' / name, assets / name)
    shutil.copyfile(documentation / 'Doxygen.svg', assets / 'Doxygen.svg')
    for name in ['wordmark-white.svg', 'favicon.svg', 'Manrope-Regular.ttf', 'Manrope-SemiBold.ttf',
                 'Manrope-Bold.ttf', 'OFL.txt']:
        shutil.copyfile(documentation / name, assets / name)
    for name in ['site.css', 'charts.js']:
        shutil.copyfile(ROOT / 'dev' / 'benchmarks' / name, assets / name)
    stylesheets += ['assets/rodin.css', 'assets/rodin-theme.css', 'assets/site.css']
    environment = Environment(loader=FileSystemLoader(ROOT / 'dev' / 'benchmarks'),
                              autoescape=select_autoescape(['html']))
    template = environment.get_template('page.html')
    pages = {'': {'children': set(), 'description': '', 'record': None}}
    modules = {}
    groups = sorted(catalogue['groups'], key=lambda group: group['name'].casefold())
    for group in groups:
        group['modules'].sort(key=lambda module: module['title'].casefold())
        for module in group['modules']:
            path = str(publication_path(module['path']))
            modules[path] = module
            for key in ('brief', 'representative', 'description'):
                module[key] = rebase_links(module[key], documentation_url if not documentation_url.startswith('../') else '__RODIN_DOCS__/')
            pages.setdefault(path, {'children': set(), 'description': '', 'record': None})
    for row in records:
        if row['state'] != 'registered':
            continue
        path = publication_path(row['path'])
        pages.setdefault(str(path), {'children': set(), 'description': '', 'record': None})
        pages[str(path)]['record'] = row
        pages[str(path)]['description'] = Path(row['description']).read_text()
    # Include ancestors once, without creating a separate physical-context taxonomy.
    for value in list(pages):
        if not value:
            continue
        parts = PurePosixPath(value).parts
        for depth in range(len(parts)):
            parent, child = '/'.join(parts[:depth]), '/'.join(parts[:depth + 1])
            pages.setdefault(parent, {'children': set(), 'description': '', 'record': None})
            pages[parent]['children'].add(child)
    # Resolve parent prose before rendering so child cards reuse the same source.
    for value, page in pages.items():
        description = ROOT / 'tests' / 'benchmarks' / value / 'index.html'
        if not page['description'] and description.is_file():
            page['description'] = description.read_text()
    # Retire pages removed by a hierarchy change in an existing output tree.
    expected = {output / path / 'index.html' for path in pages}
    for previous in output.rglob('index.html'):
        if previous not in expected:
            previous.unlink()
    # Count browsable registered benchmark pages, not timing cases or MPI ranks.
    counts = {}
    for path in sorted(pages, key=lambda value: len(PurePosixPath(value).parts), reverse=True):
        counts[path] = int(pages[path]['record'] is not None) + sum(counts[child] for child in pages[path]['children'])
    for path, module in modules.items():
        module['benchmark_count'] = counts[path]
    search = build_search(pages, modules, assets)
    contexts = []
    for path in sorted(child for child in pages['']['children'] if child != 'Rodin'):
        description = ROOT / 'tests' / 'benchmarks' / path / 'index.html'
        contexts.append({'path': path, 'title': PurePosixPath(path).name,
                         'description': description.read_text() if description.is_file() else '',
                         'benchmark_count': counts[path]})
    for value, page in pages.items():
        directory = output / value
        directory.mkdir(parents=True, exist_ok=True)
        prefix = '../' * len(PurePosixPath(value).parts) if value else ''
        titles = {path: module['title'] for path, module in modules.items()}
        titles['Rodin'] = 'Modules'
        title = titles.get(value, PurePosixPath(value).name if value else 'Benchmarks')
        breadcrumb = [('', 'Benchmarks')]
        if value:
            parts = PurePosixPath(value).parts
            breadcrumb += [('/'.join(parts[:i]), titles.get('/'.join(parts[:i]), parts[i - 1]))
                           for i in range(1, len(parts) + 1)]
        children = sorted(page['children'], key=lambda child: titles.get(child, PurePosixPath(child).name).casefold())
        cards = []
        for child in children:
            brief = re.search(r'<p(?:\s[^>]*)?>.*?</p>', pages[child]['description'], re.S)
            cards.append({'path': child, 'title': titles.get(child, PurePosixPath(child).name),
                          'brief': brief[0] if brief else '', 'benchmark_count': counts[child],
                          'record': pages[child]['record'] is not None})
        rendered = template.render(
            title=title, value=value, prefix=prefix, configuration=configuration,
            breadcrumbs=breadcrumb, page=page, module=modules.get(value), groups=groups,
            children=cards,
            contexts=contexts, footer=footer, navigation=navigation,
            stylesheets=stylesheets, escape=escape, search=search)
        (directory / 'index.html').write_text(rendered.replace('__RODIN_DOCS__/', prefix + documentation_url))
    # Keep unregistered executables visible to maintainers; never claim a complete campaign.
    (output / 'inventory.json').write_text(json.dumps(records, indent=2) + '\n')
    return len(pages)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--documentation', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--configuration', required=True)
    parser.add_argument('--documentation-url', required=True)
    args = parser.parse_args()
    count = build_site(args.manifest, args.documentation, args.output,
                       args.configuration, args.documentation_url)
    print('Rendered', count, 'benchmark navigation pages.')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Render the API site and its concise overview from the same module descriptions."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from html.parser import HTMLParser
import importlib.util
import logging
from pathlib import Path
import re
import runpy
import subprocess
import sys
import shutil
import xml.etree.ElementTree as ET

from documentation_hierarchy import prepare_hierarchy, prepare_member_definitions, symbol_index
from documentation_links import publish_enum_anchors, unlink_unpublished_references
from documentation_highlighting import cpp_symbols
from documentation_reference import prepare_reference_headers, render_references, render_admonition_notes
from documentation_search import assemble_search, prepare_class_search


class Representative(HTMLParser):
    """Locate the first complete displayed formula or highlighted usage block."""
    VOID = {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta',
            'param', 'source', 'track', 'wbr'}

    def __init__(self, source):
        super().__init__(convert_charrefs=False)
        self.source = source
        self.lines = [0] + [m.end() for m in re.finditer('\n', source)]
        self.start = None
        self.end = None
        self.depth = 0
        self.feed(source)

    def source_offset(self):
        line, column = self.getpos()
        return self.lines[line - 1] + column

    def handle_starttag(self, tag, attributes):
        if self.end is not None:
            return
        classes = dict(attributes).get('class', '').split()
        if self.start is None and ((tag == 'div' and 'm-math' in classes) or
                                   (tag == 'pre' and 'm-code' in classes)):
            self.start = self.source_offset()
        if self.start is not None and tag not in self.VOID:
            self.depth += 1

    def handle_endtag(self, tag):
        if tag in self.VOID:
            return
        if self.start is not None and self.end is None:
            self.depth -= 1
            if self.depth == 0:
                self.end = self.source.index('>', self.source_offset()) + 1

    def html(self):
        if self.start is None or self.end is None:
            raise ValueError('Module description needs a displayed formula or highlighted usage block')
        return self.source[self.start:self.end]


def scope_html(source, prefix, source_url):
    """Keep copied SVG/section IDs unique and preserve links to API-only anchors."""
    attribute = re.compile(r'(?<![\w:-])(id|href|xlink:href)\s*=\s*([\'\"])(.*?)\2')
    ids = {m[3]: prefix + m[3] for tag in re.findall(r'<[^!][^>]*>', source)
           for m in attribute.finditer(tag) if m[1] == 'id'}

    def tag_markup(match):
        def replace_attribute(item):
            name, quote, value = item.groups()
            if name == 'id':
                value = ids[value]
            elif value.startswith('#'):
                anchor = value[1:]
                value = '#' + ids[anchor] if anchor in ids else source_url + value
            return name + '=' + quote + value + quote
        tag = attribute.sub(replace_attribute, match[0])
        return re.sub(r'url\(#([^)]+)\)', lambda m: 'url(#' + ids.get(m[1], m[1]) + ')', tag)
    return re.sub(r'<[^!][^>]*>', tag_markup, source)


def overview(compound, target):
    if not compound.brief or not compound.description:
        raise ValueError('Module needs a brief and detailed description: ' + target)
    details = compound.description
    representative = Representative(details).html()
    paragraph = re.search(r'<p(?:\s[^>]*)?>.*?</p>', details, re.DOTALL)
    if paragraph is None:
        raise ValueError('Module needs an introductory description paragraph: ' + target)
    return {'description': scope_html(paragraph[0], target + '-description-', target + '.html'),
            'brief': scope_html(compound.brief, target + '-brief-', target + '.html'),
            'representative': scope_html(representative, target + '-representative-', target + '.html')}


def version_stylesheets(stylesheets, output):
    """Publish local CSS under content-based names so browsers refresh changes."""
    from urllib.parse import urlparse
    versioned = []
    for stylesheet in stylesheets:
        source = Path(stylesheet)
        if urlparse(stylesheet).netloc or not source.is_file():
            versioned.append(stylesheet)
            continue
        content = source.read_bytes()
        digest = hashlib.sha256(content).hexdigest()[:12]
        target = output / 'styles' / (source.stem + '-' + digest + source.suffix)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        versioned.append(str(target))
    return versioned


def assemble_stylesheets(sources, output):
    """Bundle the pinned upstream CSS without replacing its theme variables."""
    target = output / 'styles' / 'mcss.css'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text('\n'.join(Path(source).read_text() for source in sources))
    return str(target)


def prepare_inline_images(xml, image_root):
    """Stage authored inline images that Doxygen leaves outside its XML output."""
    from urllib.parse import urlparse
    for path in xml.glob('*.xml'):
        tree = ET.parse(path)
        changed = False
        for image in tree.findall('.//image[@type="html"][@inline="yes"]'):
            name = image.get('name')
            if urlparse(name).netloc or (xml / name).is_file():
                continue
            source = image_root / name
            if not source.is_file():
                continue
            target = xml / source.name
            shutil.copyfile(source, target)
            image.set('name', source.name)
            changed = True
        if changed:
            tree.write(path, encoding='utf-8', xml_declaration=True)


def clear_generated_pages(xml, html):
    """Discard generated pages so removed sources cannot survive a rebuild.

    Formula caches and copied assets are retained. m.css renders the complete
    site after Doxygen regenerates its XML.
    """
    for directory, pattern in [(xml, '*.xml'), (html, '*.html')]:
        for page in directory.glob(pattern):
            page.unlink()


def export_module_catalogue(groups, categories, overviews, namespaces, stylesheets, output, paths=None):
    """Expose the same rendered module descriptions to the benchmark site."""
    from urllib.parse import urlparse
    catalogue = {'schema': 1, 'groups': [], 'stylesheets': []}
    for stylesheet in stylesheets:
        catalogue['stylesheets'].append(stylesheet if urlparse(stylesheet).netloc else Path(stylesheet).name)
    for category, modules in groups:
        entries = []
        for title, target in modules:
            path = (paths or {}).get(target, namespaces[target].replace('::', '/'))
            if not path.startswith('Rodin/'):
                raise ValueError('Module needs a Rodin component path: ' + target)
            entries.append({'title': title, 'path': path, **overviews[target]})
        catalogue['groups'].append({'name': category, 'description': categories[category],
                                    'modules': entries})
    (output / 'rodin-module-catalogue.json').write_text(json.dumps(catalogue, indent=2) + '\n')
    return catalogue


def build(configuration, renderer, templates):
    # Use the repository-pinned m.css implementation for parsing and rendering;
    # no vendor edits or second copies of the scientific descriptions are needed.
    sys.path.insert(0, str(renderer.parent))
    spec = importlib.util.spec_from_file_location('rodin_mcss', renderer)
    mcss = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mcss
    spec.loader.exec_module(mcss)
    config = copy.deepcopy(mcss.default_config)
    config.update((key, value) for key, value in runpy.run_path(str(configuration)).items() if key.isupper())
    state = mcss.State(config)
    doxyfile = configuration.parent / config['DOXYFILE']
    mcss.parse_doxyfile(state, str(doxyfile))
    output = Path(state.doxyfile['OUTPUT_DIRECTORY'])
    if not output.is_absolute():
        output = doxyfile.parent / output
    output.mkdir(parents=True, exist_ok=True)
    dark_figures = config.get('RODIN_DARK_FIGURES', [])
    figure_manifest = output / 'rodin-figure-themes.json'
    figure_manifest.write_text(json.dumps({
        Path(path).name.replace('-dark.', '.'): Path(path).name
        for path in dark_figures
    }, sort_keys=True))
    state.config['EXTRA_FILES'] += dark_figures + [str(figure_manifest)]
    xml = output / state.doxyfile['XML_OUTPUT']
    clear_generated_pages(xml, output / state.doxyfile['HTML_OUTPUT'])
    subprocess.run(['doxygen', str(doxyfile)], cwd=doxyfile.parent, check=True)
    print('Resolved Doxygen member references:', prepare_member_definitions(xml))
    prepare_reference_headers(xml)
    print('Indexed forward-declared class families:', prepare_class_search(xml))
    prepare_inline_images(xml, Path(__file__).resolve().parents[1] / 'doc')
    print('Completed documentation hierarchy:', prepare_hierarchy(xml))
    stylesheets = state.config['STYLESHEETS']
    if config.get('RODIN_MCSS_STYLESHEETS'):
        stylesheets = [assemble_stylesheets(config['RODIN_MCSS_STYLESHEETS'], output)] + stylesheets
    state.config['STYLESHEETS'] = version_stylesheets(stylesheets, output)
    if 'RODIN_CPP_HIGHLIGHTER' in config:
        state.config['M_CODE_FILTERS_POST']['C++'] = config['RODIN_CPP_HIGHLIGHTER'](*cpp_symbols(xml))
    # Render the full API and search index once. Modules is produced below,
    # after the module descriptions are available; search data is not replaced.
    mcss.run(state, templates=str(templates),
             index_pages=[name for name in mcss.default_index_pages if name not in ['modules', 'annotated']])
    data = {}
    for _, modules in config['RODIN_MODULES']:
        for _, target in modules:
            parsed = mcss.parse_xml(state, str(xml / (target + '.xml')))
            data[target] = overview(parsed.compound, target)
    state.config['RODIN_MODULE_OVERVIEWS'] = data

    # This is the same Jinja environment and context used by m.css for index
    # pages. Rendering these indexes here preserves the complete API search data.
    from jinja2 import Environment, FileSystemLoader
    from urllib.parse import urljoin, urlparse
    environment = Environment(loader=FileSystemLoader([str(templates), mcss.default_templates]),
                              trim_blocks=True, lstrip_blocks=True, enable_async=True)
    environment.filters['basename_or_url'] = lambda p: p if urlparse(p).netloc else Path(p).name
    environment.filters['rtrim'] = str.rstrip
    environment.filters['urljoin'] = urljoin
    index = mcss.parse_index_xml(state, str(xml / 'index.xml'))
    index.index.symbols = symbol_index(state.compounds, xml / 'index.xml')
    for filename in ['modules.html', 'annotated.html']:
        rendered = environment.get_template(filename).render(
            index=index.index, DOXYGEN_VERSION=index.version, FILENAME=filename,
            SEARCHDATA_FORMAT_VERSION=mcss.searchdata_format_version, **state.doxyfile, **state.config)
        (output / state.doxyfile['HTML_OUTPUT'] / filename).write_text(rendered.rstrip() + '\n')
    print('Rendered reference layouts:', render_references(xml, output / state.doxyfile['HTML_OUTPUT']))
    print('Ranked general class search references:', assemble_search(
        xml, output / state.doxyfile['HTML_OUTPUT'], renderer.parent / 'search.js',
        Path(__file__).resolve().parents[1] / 'doc' / 'rodin-search.js', mcss.searchdata_format_version))
    namespaces = {target: ET.parse(xml / (target + '.xml')).findtext('./compounddef/compoundname')
                  for target in data}
    export_module_catalogue(config['RODIN_MODULES'], config['RODIN_MODULE_CATEGORIES'], data,
                            namespaces, state.config['STYLESHEETS'],
                            output / state.doxyfile['HTML_OUTPUT'], config.get('RODIN_MODULE_PATHS'))
    if config.get('RODIN_BENCHMARK_MANIFEST'):
        from benchmark_site import build_site
        html_output = output / state.doxyfile['HTML_OUTPUT']
        count = build_site(Path(config['RODIN_BENCHMARK_MANIFEST']), html_output,
                           html_output / 'benchmarks', 'Published configurations',
                           config.get('RODIN_BENCHMARK_DOCUMENTATION_URL', '../'))
        print('Rendered benchmark documentation pages:', count)
    print('Rendered m.css admonition notes:', render_admonition_notes(output / state.doxyfile['HTML_OUTPUT']))
    print('Published enum-value anchors:', publish_enum_anchors(xml, output / state.doxyfile['HTML_OUTPUT']))
    print('Unlinked unpublished scopes:', unlink_unpublished_references(xml, output / state.doxyfile['HTML_OUTPUT']))
    print('Rendered concise overview for', len(data), 'modules.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=Path)
    parser.add_argument('--renderer', type=Path, required=True)
    parser.add_argument('--templates', type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    build(args.configuration.resolve(), args.renderer.resolve(), args.templates.resolve())


if __name__ == '__main__':
    main()

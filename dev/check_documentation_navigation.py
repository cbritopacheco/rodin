#!/usr/bin/env python3
"""Verify rendered module navigation, configuration and specialization references."""
from __future__ import annotations

import argparse
import ast
from html.parser import HTMLParser
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from documentation_hierarchy import scope_prefixes


class Page(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.links = set()
        self.root_links = set()
        self.in_main = False
        self.list_depth = 0
        self.text = []
        self.identifiers = []
        self.sections = []
        self.module_roles = {}
        self.module_examples = set()
        self.categories = {}
        self.category_capture = None
        self.technical_entries = set()
        self.feed(path.read_text())

    def handle_starttag(self, tag, attributes):
        attributes = dict(attributes)
        if tag == 'main':
            self.in_main = True
        if self.in_main and tag == 'ul':
            self.list_depth += 1
        if 'id' in attributes:
            self.identifiers.append(attributes['id'])
        if tag == 'section':
            self.sections.append(attributes.get('data-module') or (self.sections[-1] if self.sections else None))
        module = self.sections[-1] if self.sections else None
        role = attributes.get('data-role')
        if module and role in ['representative', 'brief', 'description']:
            self.module_roles.setdefault(module, []).append(role)
        if module and tag in ['svg', 'pre'] and self.module_roles.get(module, [])[-1:] == ['representative']:
            self.module_examples.add(module)
        if tag == 'p' and role == 'category-description':
            self.category_capture = (attributes['data-category'], [])
        if tag == 'a' and 'href' in attributes:
            self.links.add(attributes['href'])
            if self.in_main and self.list_depth == 1:
                self.root_links.add(attributes['href'])
            if not module and attributes['href'].startswith('group__'):
                self.technical_entries.add(attributes['href'])

    def handle_endtag(self, tag):
        if self.in_main and tag == 'ul':
            self.list_depth -= 1
        if tag == 'main':
            self.in_main = False
        if tag == 'section' and self.sections:
            self.sections.pop()
        if tag == 'p' and self.category_capture is not None:
            name, text = self.category_capture
            self.categories[name] = ''.join(text).strip()
            self.category_capture = None

    def handle_data(self, text):
        self.text.append(text)
        if self.category_capture is not None:
            self.category_capture[1].append(text)


def configuration_value(configuration, name):
    for node in ast.parse(configuration.read_text()).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise ValueError('Missing documentation configuration: ' + name)


def module_catalogue(configuration):
    return configuration_value(configuration, 'RODIN_MODULES')


def specialization_targets(xml):
    compounds = ET.parse(xml / 'index.xml').findall('compound')
    classes = {c.findtext('name'): c.get('refid') for c in compounds if c.get('kind') in ('class', 'struct')}
    targets = []
    for group in [c for c in compounds if c.get('kind') == 'group']:
        definition = ET.parse(xml / (group.get('refid') + '.xml')).find('compounddef')
        title = definition.findtext('title') or ''
        if 'Template Specializations' not in title and group.findtext('name') != 'MeshTypes':
            continue
        name = title.split()[0]
        if name == 'Traits':
            candidates = ['Rodin::FormLanguage::Traits']
        else:
            candidates = [n for n, identity in classes.items() if n.endswith('::' + name) and '<' not in n
                          and ET.parse(xml / (identity + '.xml')).find('./compounddef/templateparamlist') is not None]
        # Identically named form-language traits and variational operators are
        # distinct families. Group members identify the owning namespace.
        members = {c.text.split('<', 1)[0].strip() for c in definition.findall('innerclass') if c.text}
        owned = [candidate for candidate in candidates if candidate in members]
        if len(owned) == 1:
            candidates = owned
        elif members:
            namespaces = {member.rsplit('::', 1)[0] for member in members}
            scoped = [candidate for candidate in candidates if candidate.rsplit('::', 1)[0] in namespaces]
            if len(scoped) == 1:
                candidates = scoped
        if len(candidates) != 1:
            raise ValueError('Cannot resolve primary template for ' + title)
        targets.append((classes[candidates[0]], group.get('refid'), title))
    if not targets:
        raise ValueError('No specialization references were verified')
    return targets


def verify_reference_indexes(site, xml):
    files = Page(site / 'files.html')
    classes = Page(site / 'annotated.html')
    counts = {'files': 0, 'types': 0}
    for item in ET.parse(xml / 'index.xml').getroot():
        kind = item.get('kind')
        if kind not in ['file', 'class', 'struct', 'union']:
            continue
        identity = item.get('refid')
        definition = ET.parse(xml / (identity + '.xml')).find('compounddef')
        if definition.get('language', 'C++') != 'C++' or definition.get('prot') == 'private':
            continue
        brief = definition.find('briefdescription')
        details = definition.find('detaileddescription')
        if not ((brief is not None and ''.join(brief.itertext()).strip()) or
                (details is not None and len(details))):
            continue
        if kind != 'file' and scope_prefixes(item.findtext('name')) and identity + '.html' in classes.root_links:
            raise ValueError('Type escaped its namespace hierarchy: ' + item.findtext('name'))
        links = files.links if kind == 'file' else classes.links
        if identity + '.html' not in links or not (site / (identity + '.html')).is_file():
            raise ValueError('Documented declaration missing from reference index: ' + item.findtext('name'))
        counts['files' if kind == 'file' else 'types'] += 1
    return counts


def verify(site, xml, configuration):
    index = Page(site / 'modules.html')
    # These are technical reference families, not Rodin components.
    if index.technical_entries:
        raise ValueError('Technical template/directive groups leaked into the Modules index')
    if 'configuration.html' not in index.links:
        raise ValueError('Modules index does not link to Configuration')
    categories = configuration_value(configuration, 'RODIN_MODULE_CATEGORIES')
    for category, _ in module_catalogue(configuration):
        if not categories.get(category) or index.categories.get(category) != categories[category]:
            raise ValueError('Missing category description: ' + category)
    modules = [row for _, rows in module_catalogue(configuration) for row in rows]
    if len(index.identifiers) != len(set(index.identifiers)):
        raise ValueError('Duplicate HTML or SVG IDs in the Modules overview')
    for link in index.links:
        if link.startswith('#') and link not in ['#', '#!'] and link[1:] not in index.identifiers:
            raise ValueError('Broken overview anchor: ' + link)
    for label, target in modules:
        if index.module_roles.get(target) != ['brief', 'representative', 'description']:
            raise ValueError('Wrong or missing module overview order: ' + label)
        if target not in index.module_examples:
            raise ValueError('Overview lacks a rendered formula or usage: ' + label)
        path = site / (target + '.html')
        if target + '.html' not in index.links or not path.is_file():
            raise ValueError('Missing module page/link: ' + label)
        page = Page(path)
        definition = ET.parse(xml / (target + '.xml')).find('compounddef')
        brief = definition.find('briefdescription')
        detail = definition.find('detaileddescription')
        if brief is None or not ''.join(brief.itertext()).strip() or detail is None:
            raise ValueError('Missing module description: ' + label)
        if detail.find('.//formula') is None and detail.find('.//programlisting') is None:
            raise ValueError('Module needs a representative formula or usage: ' + label)
        if not re.search(r'<span class="m-thin">module</span>', path.read_text()):
            raise ValueError('Public module is still labeled as a namespace: ' + label)
        if re.search(r'Directory\s+(?:src/|Rodin/)', ' '.join(page.text)):
            raise ValueError('Stray directory paragraph: ' + label)
    configuration_page = Page(site / 'configuration.html')
    if '_configure_8h_8in.html' not in configuration_page.links or not (site / '_configure_8h_8in.html').is_file():
        raise ValueError('Configuration does not link to the header reference')
    macros = ET.parse(xml / '_configure_8h_8in.xml').findall('./compounddef/sectiondef/memberdef[@kind="define"]')
    if not any(m.findtext('name') == 'RODIN_USE_MPI' for m in macros):
        raise ValueError('Configuration macros are missing from their file reference')
    targets = specialization_targets(xml)
    for primary, group, title in targets:
        definition = ET.parse(xml / (primary + '.xml')).find('compounddef')
        see = definition.findall('./detaileddescription//simplesect[@kind="see"]//ref')
        if not any(ref.get('refid') == group for ref in see):
            raise ValueError('Primary template has no See also link: ' + title)
        page = Page(site / (primary + '.html'))
        if group + '.html' not in page.links or not (site / (group + '.html')).is_file():
            raise ValueError('Specialization reference link is missing/broken: ' + title)
    if not (site / 'rodin-navigation.js').is_file():
        raise ValueError('Missing sidebar navigation asset')
    for page in ['modules.html', 'configuration.html'] + [target + '.html' for _, target in modules]:
        if '<script src="rodin-navigation.js" defer></script>' not in (site / page).read_text():
            raise ValueError('Page does not load sidebar navigation: ' + page)
    print('Reference index coverage:', verify_reference_indexes(site, xml))
    print('Verified', len(modules), 'module pages, Configuration and', len(targets), 'specialization references.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site', type=Path, default=Path('build/doc/html'))
    parser.add_argument('--xml', type=Path, default=Path('build/doc/xml'))
    parser.add_argument('--configuration', type=Path, default=Path('doc/conf.py.in'))
    args = parser.parse_args()
    verify(args.site, args.xml, args.configuration)


if __name__ == '__main__':
    main()

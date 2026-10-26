"""Assemble Rodin search ranking around the unchanged upstream search engine."""
import hashlib
import json
import re
from pathlib import Path
import xml.etree.ElementTree as ET

from documentation_reference import compact, scope_name, text


def prepare_class_search(xml):
    """Index family landing pages whose bare forward declarations lack prose."""
    definitions = {}
    for path in Path(xml).glob('*.xml'):
        tree = ET.parse(path)
        definition = tree.find('compounddef')
        if definition is not None and definition.get('kind') in ('class', 'struct', 'union'):
            definitions[compact(definition.findtext('compoundname'))] = (path, tree, definition)
    families = set()
    for _, _, definition in definitions.values():
        scope, name = scope_name(definition.findtext('compoundname'))
        if '<' in name:
            families.add(compact((scope + '::' if scope else '') + name.split('<', 1)[0]))
    count = 0
    for name in families:
        if name not in definitions:
            continue
        path, tree, definition = definitions[name]
        brief = definition.find('briefdescription')
        if brief is not None and not text(brief) and not text(definition.find('detaileddescription')):
            ET.SubElement(brief, 'para').text = 'Primary template for ' + scope_name(definition.findtext('compoundname'))[1] + '.'
            tree.write(path, encoding='utf-8', xml_declaration=True)
            count += 1
    return count


def primary_classes(xml, site):
    """Catalogue general class references; specializations retain normal lookup."""
    entries = []
    for path in sorted(Path(xml).glob('*.xml')):
        definition = ET.parse(path).find('compounddef')
        if definition is None or definition.get('kind') not in ('class', 'struct', 'union'):
            continue
        name = definition.findtext('compoundname') or ''
        if not name or '<' in name:
            continue
        url = definition.get('id') + '.html'
        if not (Path(site) / url).is_file():
            continue
        scopes = name.split('::')
        entries.append({'name': name, 'url': url,
                        'keys': ['::'.join(scopes[index:]).lower() for index in range(len(scopes))]})
    return entries


def assemble_search(xml, site, upstream, extension, version):
    """Publish one bundle, retaining the upstream implementation and licence."""
    catalogue = primary_classes(xml, site)
    script = Path(extension).read_text()
    if script.count('__RODIN_SEARCH_CLASSES__') != 1:
        raise ValueError('Search extension must have one catalogue placeholder')
    script = script.replace('__RODIN_SEARCH_CLASSES__', json.dumps(catalogue, ensure_ascii=True))
    target = Path(site) / f'search-v{version}.js'
    bundle = Path(upstream).read_text() + '\n' + script
    target.write_text(bundle)
    digest = hashlib.sha256(bundle.encode()).hexdigest()[:12]
    versioned = f'search-v{version}-{digest}.js'
    (Path(site) / versioned).write_text(bundle)
    for page in Path(site).glob('*.html'):
        source = page.read_text()
        rendered = re.sub(r'(?<=src=")search-v' + str(version) + r'(?:-[0-9a-f]{12})?\.js(?=")', versioned, source)
        if rendered != source:
            page.write_text(rendered)
    return len(catalogue)

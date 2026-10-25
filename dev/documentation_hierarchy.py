"""Complete navigation relationships in generated Doxygen XML before rendering."""
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def scope_prefixes(name):
    """Yield enclosing C++ scopes, ignoring :: within template arguments."""
    depth = 0
    positions = []
    for index, character in enumerate(name):
        if character == '<':
            depth += 1
        elif character == '>':
            depth -= 1
        elif character == ':' and name[index:index + 2] == '::' and depth == 0:
            positions.append(index)
    return [name[:index] for index in reversed(positions)]


def normalized_name(name):
    return re.sub(r'\s*(::|[<>,])\s*', r'\1', name.strip())


def prepare_hierarchy(xml):
    """Retain header directories and attach orphan types to known enclosing scopes.

    m.css requires documentation on directory nodes to retain their children.
    Doxygen also omits some innerclass links for nested specialization members.
    Only generated XML is changed; declarations and vendor sources are untouched.
    """
    xml = Path(xml)
    compounds = {}
    documented_files = set()
    names = {}
    owned = set()
    changed = set()
    counts = {'directories': 0, 'ownership': 0}
    for entry in ET.parse(xml / 'index.xml').getroot():
        if entry.get('kind') not in ['dir', 'file', 'namespace', 'class', 'struct', 'union']:
            continue
        identity = entry.get('refid')
        tree = ET.parse(xml / (identity + '.xml'))
        definition = tree.find('compounddef')
        if definition.get('kind') == 'file':
            if definition.get('language', 'C++') == 'C++' and (
                    ''.join(definition.find('briefdescription').itertext()).strip() or
                    len(definition.find('detaileddescription'))):
                documented_files.add(identity)
            continue
        compounds[identity] = (tree, definition)
        if definition.get('kind') != 'dir':
            name = normalized_name(definition.findtext('compoundname'))
            # Doxygen can emit anonymous namespace records whose compoundname
            # repeats a public namespace. Such aliases must not replace its ID.
            if definition.get('kind') == 'namespace' and name != normalized_name(entry.findtext('name')):
                continue
            names[name] = identity
            owned.update(child.get('refid') for child in definition.findall('innerclass'))

    directories = {}

    def has_documented_files(identity):
        if identity in directories:
            return directories[identity]
        definition = compounds[identity][1]
        found = any(child.get('refid') in documented_files for child in definition.findall('innerfile'))
        found = found or any(has_documented_files(child.get('refid')) for child in definition.findall('innerdir'))
        directories[identity] = found
        return found

    for identity, (_, definition) in compounds.items():
        if definition.get('kind') != 'dir' or not has_documented_files(identity):
            continue
        if not ''.join(definition.find('briefdescription').itertext()).strip() and not len(definition.find('detaileddescription')):
            ET.SubElement(definition.find('briefdescription'), 'para').text = (
                'Headers in ' + definition.findtext('compoundname') + '.')
            changed.add(identity)
            counts['directories'] += 1

    for identity, (_, definition) in compounds.items():
        if definition.get('kind') not in ['class', 'struct', 'union'] or identity in owned:
            continue
        name = definition.findtext('compoundname')
        for prefix in scope_prefixes(normalized_name(name)):
            parent = names.get(prefix)
            if parent is None or parent == identity:
                continue
            owner = compounds[parent][1]
            child = ET.SubElement(owner, 'innerclass', refid=identity,
                                  prot=definition.get('prot', 'public'))
            child.text = name
            owned.add(identity)
            changed.add(parent)
            counts['ownership'] += 1
            break

    ET.register_namespace('mcss', 'http://mcss.mosra.cz/doxygen/')
    for identity in changed:
        compounds[identity][0].write(xml / (identity + '.xml'), encoding='utf-8', xml_declaration=True)
    return counts


def symbol_index(compounds, index_xml):
    """Keep documented types reachable through their nearest documented ancestor.

    An enclosing specialization or internal namespace may have no description of
    its own. Its documented descendants still belong under the enclosing public
    namespace, with the omitted scope retained in their displayed names.
    """
    from types import SimpleNamespace
    entries = {}
    kinds = {'namespace', 'class', 'struct', 'union'}
    for item in ET.parse(index_xml).getroot():
        identity = item.get('refid')
        compound = compounds.get(identity)
        if compound is None or compound.kind not in kinds or not compound.has_details:
            continue
        if compound.kind == 'namespace' and compound.name != item.findtext('name'):
            continue
        definition = ET.parse(Path(index_xml).parent / (identity + '.xml')).find('compounddef')
        if definition.get('prot') == 'private':
            continue
        entries[identity] = SimpleNamespace(id=identity, kind=compound.kind,
            name=compound.name, url=compound.url, brief=compound.brief, children=[],
            deprecated=compound.deprecated, since=compound.since,
            is_inline=getattr(compound, 'is_inline', False),
            is_final=getattr(compound, 'is_final', False), has_nestable_children=False)
    roots = []
    for identity, entry in entries.items():
        parent = compounds[identity].parent
        seen = {identity}
        while parent and parent not in entries:
            if parent in seen:
                raise ValueError('Cyclic documentation ownership: ' + identity)
            seen.add(parent)
            parent = compounds[parent].parent
        if parent:
            prefix = compounds[parent].name + '::'
            if not entry.name.startswith(prefix):
                raise ValueError('Inconsistent documentation scope: ' + entry.name)
            entry.name = entry.name[len(prefix):]
            entries[parent].children.append(entry)
            if entry.kind == 'namespace':
                entries[parent].has_nestable_children = True
        else:
            roots.append(entry)
    for children in [roots] + [entry.children for entry in entries.values()]:
        children.sort(key=lambda entry: entry.kind != 'namespace')
    root_order = {'Rodin': 0, 'Eigen': 1, 'std': 2}
    roots.sort(key=lambda entry: root_order.get(entry.name, 3))
    return roots

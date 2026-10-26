"""Semantic additions to Pygments highlighting for documented Rodin C++ APIs."""
import re
import xml.etree.ElementTree as ET


DEFAULT_TYPES = {'Mesh', 'Polytope', 'Type', 'P1', 'TrialFunction', 'TestFunction',
                 'Problem', 'Integral', 'Grad', 'DirichletBC', 'Zero', 'CG',
                 'RealFunction', 'XDMF'}


def cpp_symbols(xml):
    """Collect type names and enum-value scopes from the generated declarations."""
    types = set(DEFAULT_TYPES)
    enums = {}
    identifier = re.compile(r'[A-Za-z_][A-Za-z_0-9]*\Z')
    for path in xml.glob('*.xml'):
        definition = ET.parse(path).find('compounddef')
        if definition is None or definition.get('language', 'C++') != 'C++':
            continue
        if definition.get('kind') in ['class', 'struct', 'union']:
            depth = 0
            outer = ''
            for character in definition.findtext('compoundname') or '':
                if character == '<':
                    depth += 1
                elif character == '>':
                    depth -= 1
                elif depth == 0:
                    outer += character
            name = outer.rsplit('::', 1)[-1].strip()
            if identifier.fullmatch(name) and name[0].isupper():
                types.add(name)
        for member in definition.findall('./sectiondef/memberdef'):
            name = member.findtext('name') or ''
            if member.get('kind') == 'typedef' and identifier.fullmatch(name) and name[0].isupper():
                types.add(name)
            if member.get('kind') != 'enum':
                continue
            if identifier.fullmatch(name):
                types.add(name)
            owner = member.findtext('qualifiedname') or ''
            if member.get('strong') != 'yes':
                owner = owner.rpartition('::')[0]
            for value in member.findall('enumvalue'):
                name = value.findtext('name') or ''
                if identifier.fullmatch(name):
                    enums.setdefault(name, set()).add(owner)
    return types, enums


def cpp_highlighter(type_names=None, enum_values=None):
    """Classify identifiers while preserving the formatter's literals and comments."""
    types = DEFAULT_TYPES if type_names is None else type_names
    enums = {'Triangle': {'Rodin::Geometry::Polytope::Type'}} if enum_values is None else enum_values
    span = re.compile(r'<span class="([^"]+)">([^<]*)</span>')
    identifier = re.compile(r'[A-Za-z_][A-Za-z_0-9]*\Z')

    def highlight(markup):
        tokens = list(span.finditer(markup))
        replacements = {}
        for index, token in enumerate(tokens):
            kind, name = token.groups()
            if kind not in ['n', 'nf', 'nc', 'nn', 'no'] or not identifier.fullmatch(name):
                continue
            qualifiers = []
            previous = index - 1
            while previous >= 0:
                while previous >= 0 and tokens[previous][1] == 'w':
                    previous -= 1
                if previous < 0 or tokens[previous][2] != '::':
                    break
                previous -= 1
                while previous >= 0 and tokens[previous][1] == 'w':
                    previous -= 1
                if previous < 0 or not identifier.fullmatch(tokens[previous][2]):
                    break
                qualifiers.insert(0, tokens[previous][2])
                previous -= 1
            scope = '::'.join(qualifiers)
            owners = enums.get(name, set())
            is_enum = bool(owners) and (
                any(owner == scope or owner.endswith('::' + scope) for owner in owners)
                if scope else name not in types)
            css = 'no' if is_enum else 'nc' if name in types else 'fm' if name == 'main' else None
            if css:
                replacements[token.start()] = '<span class="' + css + '">' + name + '</span>'
        return span.sub(lambda match: replacements.get(match.start(), match[0]), markup)
    return highlight

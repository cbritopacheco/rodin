"""Present Doxygen reference declarations without repeating their enclosing scope.

The XML remains authoritative for template parameters, specialization arguments
and source locations. Only rendered reference markup is changed; search data,
member anchors and the upstream renderer are retained.
"""
from dataclasses import dataclass, field
from html import escape, unescape
from html.parser import HTMLParser
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import CppLexer

from documentation_highlighting import cpp_highlighter, cpp_symbols


def text(element):
    return ''.join(element.itertext()).strip() if element is not None else ''


def compact(value):
    """Ignore Doxygen's optional whitespace when comparing C++ spellings."""
    return re.sub(r'\s+', '', value)


def split_arguments(value):
    """Split comma-separated arguments, respecting nested C++ delimiters."""
    result, start, stack, quote = [], 0, [], None
    escaped = False
    for index, character in enumerate(value):
        if quote:
            if character == quote and not escaped:
                quote = None
            escaped = character == '\\' and not escaped
            continue
        if character in '\"\'':
            quote = character
        elif character in '<([{':
            stack.append(character)
        elif character in '>)]}' and stack:
            stack.pop()
        elif character == ',' and not stack:
            result.append(value[start:index].strip())
            start = index + 1
    if value[start:].strip():
        result.append(value[start:].strip())
    return result


def scope_name(value):
    """Find the last scope separator outside template arguments."""
    depth, last = 0, -1
    for index, character in enumerate(value):
        if character == '<':
            depth += 1
        elif character == '>':
            depth -= 1
        elif depth == 0 and value[index:index + 2] == '::':
            last = index
    return (value[:last], value[last + 2:]) if last >= 0 else ('', value)


def parameters(definition):
    result = []
    for index, param in enumerate(definition.findall('./templateparamlist/param')):
        kind = text(param.find('type'))
        name = param.findtext('declname') or param.findtext('defname')
        if not name:
            match = re.search(r'(?:class|typename)\s*(?:\.\.\.\s*)?([A-Za-z_]\w*)$', kind)
            name = match[1] if match else f'_{index + 1}'
        declaration = kind if re.search(r'\b' + re.escape(name) + r'$', kind) else kind + ' ' + name
        default = text(param.find('defval'))
        result.append((name, declaration, default, '...' in kind))
    return result


@dataclass
class Node:
    tag: str
    attrs: dict
    start: int
    end: int = 0
    parent: object = None
    children: list = field(default_factory=list)


class Markup(HTMLParser):
    """Locate replaceable elements while retaining untouched HTML byte-for-byte."""
    VOID = {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link',
            'meta', 'param', 'source', 'track', 'wbr'}

    def __init__(self, source):
        super().__init__(convert_charrefs=False)
        self.source, self.nodes, self.stack = source, [], []
        self.lines = [0] + [m.end() for m in re.finditer('\n', source)]
        self.feed(source)

    def source_offset(self):
        line, column = self.getpos()
        return self.lines[line - 1] + column

    def handle_starttag(self, tag, attrs):
        node = Node(tag, dict(attrs), self.source_offset(), parent=self.stack[-1] if self.stack else None)
        self.nodes.append(node)
        if node.parent:
            node.parent.children.append(node)
        if tag in self.VOID:
            node.end = self.source.index('>', node.start) + 1
        else:
            self.stack.append(node)

    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        if tag not in self.VOID:
            self.stack.pop().end = self.source_offset() + len(self.get_starttag_text())

    def handle_endtag(self, tag):
        if tag in self.VOID:
            return
        if self.stack and self.stack[-1].tag == tag:
            self.stack.pop().end = self.source.index('>', self.source_offset()) + 1

    def raw(self, node):
        return self.source[node.start:node.end]

    def content(self, node):
        return unescape(re.sub(r'<[^>]+>', '', self.raw(node))).strip()


def replace_ranges(source, edits):
    """Apply disjoint edits, refusing accidental deletion of overlapping content."""
    previous = len(source)
    for start, end, replacement in sorted(edits, reverse=True):
        if end > previous or start > end:
            raise ValueError('Overlapping reference layout edits')
        source = source[:start] + replacement + source[end:]
        previous = start
    return source


class ReferenceLayout:
    def __init__(self, xml, site):
        self.site = Path(site)
        self.definitions = {}
        for path in Path(xml).glob('*.xml'):
            definition = ET.parse(path).find('compounddef')
            if definition is not None and definition.get('kind') in ('class', 'struct', 'union', 'namespace'):
                self.definitions[definition.get('id')] = definition
        self.members = {member.get('id'): member for definition in self.definitions.values()
                        for member in definition.findall('./sectiondef/memberdef')}
        self.names = {compact(d.findtext('compoundname')): identity for identity, d in self.definitions.items()}
        self.inheritance_edges = {'basecompoundref': {}, 'derivedcompoundref': {}}
        for identity, definition in self.definitions.items():
            for direction in self.inheritance_edges:
                for relation in definition.findall(direction):
                    target = relation.get('refid')
                    if not target or relation.get('prot') == 'private':
                        continue
                    base, child = (target, identity) if direction == 'basecompoundref' else (identity, target)
                    self.inheritance_edges['basecompoundref'].setdefault(child, set()).add(base)
                    self.inheritance_edges['derivedcompoundref'].setdefault(base, set()).add(child)

        self.families, self.primaries = {}, {}
        for identity, definition in self.definitions.items():
            scope, name = scope_name(definition.findtext('compoundname'))
            if '<' not in name or not name.rstrip().endswith('>'):
                continue
            primary = self.names.get(compact((scope + '::' if scope else '') + name.split('<', 1)[0]))
            if primary and primary != identity and self.exists(identity) and self.exists(primary):
                self.families.setdefault(primary, []).append(identity)
                self.primaries[identity] = primary
        for members in self.families.values():
            members.sort(key=lambda member: compact(self.definitions[member].findtext('compoundname')))
        self.colorize = cpp_highlighter(*cpp_symbols(Path(xml)))

    def exists(self, identity):
        return (self.site / (identity + '.html')).is_file()

    def code(self, value, definition=None):
        markup = self.colorize(highlight(value, CppLexer(), HtmlFormatter(nowrap=True)))
        if definition is not None:
            candidates = {}
            for reference in definition.findall('.//ref'):
                identity = reference.get('refid', '')
                target = identity + '.html'
                if reference.get('kindref') == 'member':
                    match = re.fullmatch(r'(.*)_1([a-z]+[0-9a-f]{32})', identity)
                    if not match:
                        continue
                    target = match[1] + '.html#' + match[2]
                if not (self.site / target.split('#')[0]).is_file():
                    continue
                label = scope_name(text(reference))[1].split('<', 1)[0]
                candidates.setdefault(label, set()).add((text(reference), target))

            def linked(token):
                label = unescape(token[2])
                prefix = unescape(re.sub(r'<[^>]+>', '', markup[:token.start()])).rstrip()
                qualified = re.search(r'([A-Za-z_]\w*(?:::[A-Za-z_]\w*)*::)$', prefix)
                references = candidates.get(label, set())
                if prefix.endswith('::'):
                    # A dependent Traits<T>::Member is not an unrelated alias
                    # named Member in the current class.
                    spelling = qualified[1] + label if qualified else None
                    references = {(name, target) for name, target in references
                                  if spelling and (name == spelling or name.endswith('::' + spelling))}
                targets = {target for _, target in references}
                if len(targets) != 1:
                    return token[0]
                return '<a href="' + escape(next(iter(targets))) + '">' + token[0] + '</a>'
            markup = re.sub(r'<span class="(nc|n|no)">([^<]+)</span>', linked, markup)
        return '<pre class="m-code rodin-declaration">' + markup + '</pre>'

    def declaration(self, definition):
        _, name = scope_name(definition.findtext('compoundname'))
        params = parameters(definition)
        result = ''
        if definition.find('templateparamlist') is not None:
            declarations = [p[1] + (' = ' + p[2] if p[2] else '') for p in params]
            result = 'template <' + (('\n  ' + ',\n  '.join(declarations)) if declarations else '') + '>\n'
        if '<' in name and name.rstrip().endswith('>'):
            base, args = name.split('<', 1)
            args = split_arguments(args.rstrip()[:-1])
            name = base + '<\n  ' + ',\n  '.join(args) + '>'
        requirement = text(definition.find('requiresclause'))
        if requirement:
            result += 'requires ' + requirement + '\n'
        return result + definition.get('kind') + ' ' + name + (' final' if definition.get('final') == 'yes' else '') + ';'

    def bindings(self, identity):
        primary = self.primaries.get(identity)
        if not primary:
            return ''
        definition = self.definitions[identity]
        _, name = scope_name(definition.findtext('compoundname'))
        arguments = split_arguments(name.split('<', 1)[1].rstrip()[:-1])
        params, rows, cursor = parameters(self.definitions[primary]), [], 0
        for param, _, _, pack in params:
            if cursor >= len(arguments):
                break
            value = ', '.join(arguments[cursor:]) if pack else arguments[cursor]
            rows.append('<tr><td><code>' + escape(param + ('...' if pack else '')) + '</code></td><td><code>' + escape(value) + '</code></td></tr>')
            cursor = len(arguments) if pack else cursor + 1
        # Do not invent a mapping if metadata cannot identify every argument.
        if cursor != len(arguments):
            rows = ['<tr><td>Arguments</td><td><code>' + escape(', '.join(arguments)) + '</code></td></tr>']
        return ('<section id="rodin-specialization"><h2><a href="#rodin-specialization">Specialization arguments</a></h2>'
                '<p>Specialization of <a href="' + primary + '.html">' + escape(scope_name(self.definitions[primary].findtext('compoundname'))[1]) + '</a>.</p>'
                '<table class="m-table m-fullwidth m-flat rodin-bindings"><tbody>' + ''.join(rows) + '</tbody></table></section>')

    def overview_argument(self, value):
        # Overview cells reserve their horizontal space for the argument shape.
        # Complete constructor-independent C++ spellings remain in the disclosure.
        maximum_length = 40
        value = re.sub(r'\s*::\s*', '::', value)
        value = re.sub(r'<\s*', '<', value)
        value = re.sub(r'\s*>', '>', value)
        if len(value) <= maximum_length or '<' not in value or not value.endswith('>'):
            return value
        name, arguments = value.split('<', 1)
        args = split_arguments(arguments[:-1])
        if len(args) > 1 and len(args[0]) < maximum_length // 2:
            return name + '<' + args[0] + ', …>'
        return name + '<…>'

    def family(self, identity, members=None):
        rows = []
        for member in members if members is not None else self.families.get(identity, []):
            definition = self.definitions[member]
            _, name = scope_name(definition.findtext('compoundname'))
            arguments = name[name.index('<'):]
            display = ', '.join(self.overview_argument(a) for a in split_arguments(arguments[1:-1]))
            brief = text(definition.find('briefdescription'))
            rows.append('<tr><td><a href="' + member + '.html"><code>' + escape(display) + '</code></a>'
                        '<details><summary>Full specialization</summary>' + self.code(self.declaration(definition), definition) + '</details></td><td>' + escape(brief) + '</td></tr>')
        return ('<table class="m-table m-fullwidth rodin-specializations"><thead><tr><th>Specialization</th><th>Description</th></tr></thead><tbody>' + ''.join(rows) + '</tbody></table>')

    def class_header(self, identity, markup):
        definition = self.definitions[identity]
        header = next((n for n in markup.nodes if n.tag == 'h1'), None)
        if header is None:
            return []
        _, name = scope_name(definition.findtext('compoundname'))
        title = name.split('<', 1)[0].strip()
        # Retain the renderer's linked scope breadcrumb and status labels.
        links = re.findall(r'<span class="m-breadcrumb">(.*?)</span>', markup.raw(header), re.S)
        breadcrumb = '<div class="rodin-reference-breadcrumb">' + '<span>›</span>'.join(re.sub(r'::<wbr\s*/?>', '', link) for link in links) + '</div>' if links else ''
        labels = ''.join(re.findall(r'<span class="m-label[^>]*>.*?</span>', markup.raw(header), re.S))
        result = breadcrumb + '<h1>' + escape(title) + (' ' + labels if labels else '') + '</h1>'
        brief = next((n for n in header.parent.children if n.start >= header.end), None)
        end = header.end
        if brief is not None and brief.tag == 'p':
            result += markup.raw(brief).replace('<p>', '<p class="rodin-reference-brief">', 1)
            end = brief.end
        location = definition.find('location')
        if location is not None:
            source = location.get('file', '')
            includes = definition.find('includes')
            target = includes.get('refid') if includes is not None else None
            if target and self.exists(target):
                result += '<p class="rodin-reference-include">Defined in <a href="' + target + '.html"><code>' + escape(source) + '</code></a></p>'
        result += self.code(self.declaration(definition), definition) + self.bindings(identity)
        return [(header.start, end, result)]

    def member_signatures(self, definition, markup):
        owner = definition.findtext('compoundname')
        parameter_names = ','.join(p[0] for p in parameters(definition))
        accepted = {compact(owner + '::'), compact(owner + '<' + parameter_names + '>::')}
        edits = []
        for node in markup.nodes:
            if node.tag != 'h3' or not node.parent or not node.parent.parent or 'm-doc-details' not in node.parent.parent.attrs.get('class', '').split():
                continue
            source = markup.raw(node)
            # The first wrap bumper ends in the owning scope; its preceding HTML
            # can contain linked return types and must remain untouched.
            match = re.search(r'(<span class="m-doc-wrap-bumper">)(.*?)(</span>)', source, re.S)
            if match:
                body = match[2]
                # The generated owner prefix is unlinked, apart from wbr tags.
                for candidate in re.finditer(r'(?:[A-Za-z_]\w*|<wbr\s*/?>|::|&lt;|&gt;|[<>,.\s])+', body):
                    suffix = candidate[0]
                    for start in re.finditer(r'[A-Za-z_]\w*(?=::)', suffix):
                        raw = suffix[start.start():]
                        plain = unescape(re.sub(r'<[^>]+>', '', raw))
                        if compact(plain) in accepted:
                            offset = match.start(2) + candidate.start() + start.start()
                            source = source[:offset] + source[match.end(2):]
                            break
                    else:
                        continue
                    break
            # Typedef headings use a direct prefix instead of a wrap bumper.
            # Keep link targets and external qualifications intact.
            alias = re.search(r'<a[^>]+class="m-doc-self"', source)
            if alias and re.search(r'\b(?:using|typedef)\b', source[:alias.start()]):
                before = source[:alias.start()]
                for prefix in accepted:
                    tokens = re.findall(r'&[^;]+;|.', escape(prefix))
                    gap = r'(?:\s|<wbr\s*/?>)*'
                    pattern = r'(?<![\w:])' + gap.join(re.escape(token) for token in tokens) + gap + '$'
                    before = re.sub(pattern, '', before)
                source = before + source[alias.start():]
            # Remove only the enclosing class-template declaration. Function
            # templates, including their own parameters/defaults, remain visible.
            if definition.get('kind') != 'namespace' and parameters(definition):
                block = re.search(r'<div class="m-doc-template">(.*?)</div>', source, re.S)
                if block:
                    expected = 'template<' + ','.join(p[1] for p in parameters(definition)) + '>'
                    plain = unescape(re.sub(r'<[^>]+>', '', block[1]))
                    if compact(plain) == compact(expected):
                        source = source[:block.start()] + source[block.end():]
            if source != markup.raw(node):
                edits.append((node.start, node.end, source))
        return edits

    def relative_name(self, value, definition):
        scope = (definition.findtext('compoundname') if definition.get('kind') == 'namespace'
                 else scope_name(definition.findtext('compoundname'))[0])
        if scope and value.startswith(scope + '::'):
            return value[len(scope) + 2:]
        if scope.startswith('Rodin::') and value.startswith('Rodin::'):
            return value[len('Rodin::'):]
        return value

    def member_target(self, identity):
        match = re.fullmatch(r'(.*)_1([a-z]+[0-9a-f]{32})', identity)
        if match and self.exists(match[1]):
            return match[1], match[1] + '.html#' + match[2]
        return None

    def relation_name(self, relation, definition):
        value = text(relation)
        target = self.definitions.get(relation.get('refid', ''))
        # Doxygen can shorten an out-of-module relation to its bare name.
        # Recover its scope from the referenced declaration, retaining the
        # relation's actual template arguments.
        if target is not None and not scope_name(value)[0]:
            scope = scope_name(target.findtext('compoundname'))[0]
            if scope:
                value = scope + '::' + value
        return self.relative_name(value, definition)

    def relation_groups(self, relations):
        groups = {}
        for relation in relations:
            target = relation.get('refid', '')
            scope, name = scope_name(text(relation))
            key = self.primaries.get(target, target) or compact(scope + '::' + name.split('<', 1)[0])
            groups.setdefault(key, []).append(relation)
        return list(groups.items())

    def inheritance_diagram(self, identity, groups):
        definition = self.definitions[identity]
        maximum_nodes = 3
        rows = []
        for relations, target in groups:
            if not relations:
                rows.append('')
                continue
            nodes = []
            families = self.relation_groups(relations)
            for family, members in families[:maximum_nodes]:
                relation = members[0]
                full = self.relation_name(relation, definition)
                if family in self.definitions:
                    label = self.relative_name(self.definitions[family].findtext('compoundname'), definition)
                    if '<' in scope_name(label)[0]:
                        label = scope_name(label)[1]
                    label = label.split('<', 1)[0]
                else:
                    label = full.split('<', 1)[0]
                url = family + '.html'
                tooltip = '; '.join(self.relation_name(member, definition) for member in members)
                content = '<a href="' + escape(url) + '" title="' + escape(tooltip) + '">' + escape(label) + '</a>' if self.exists(family) else escape(label)
                nodes.append('<li>' + content + '</li>')
            if len(relations) > maximum_nodes:
                remaining = len(relations) - min(len(families), maximum_nodes)
                nodes.append('<li class="rodin-inheritance-more"><a href="#' + target + '">+' + str(remaining) + ' more</a></li>')
            rows.append('<ul class="rodin-inheritance-row">' + ''.join(nodes) + '</ul>')
        label = scope_name(definition.findtext('compoundname'))[1].split('<', 1)[0]
        return ('<figure class="rodin-inheritance"><figcaption>Direct inheritance by class family</figcaption>' + rows[0] +
                ('<div class="rodin-inheritance-arrow" aria-label="Inherited by">↓</div>' if rows[0] else '') +
                '<div class="rodin-inheritance-current">' + escape(label) + '</div>' +
                ('<div class="rodin-inheritance-arrow" aria-label="Inherited by">↓</div>' + rows[1] if rows[1] else '') + '</figure>')

    def inheritance(self, identity, markup):
        definition = self.definitions[identity]
        if definition.get('kind') == 'namespace':
            return []
        edits, groups, sections = [], [], []
        for tag, section_id, title, first_column in [
                ('basecompoundref', 'base-classes', 'Inherits from', 'Base class'),
                ('derivedcompoundref', 'derived-classes', 'Derived classes', 'Derived class')]:
            section = next((n for n in markup.nodes if n.tag == 'section' and n.attrs.get('id') == section_id), None)
            original = markup.raw(section) if section else ''
            relations = [ref for ref in definition.findall(tag) if ref.get('prot') != 'private' or ref.get('refid', '') + '.html' in original]
            groups.append((relations, section_id))
            if tag == 'derivedcompoundref':
                families = self.relation_groups(relations)
                relations = [members[0] for _, members in families] + [member for _, members in families for member in members[1:]]
            if not section or not relations:
                continue
            sections.append(section)
            rows = []
            # Preserve the renderer's linked and mathematically rendered briefs.
            briefs = {}
            for listing in [n for n in markup.nodes if n.tag == 'dl' and section.start < n.start < section.end]:
                for index, node in enumerate(listing.children[:-1]):
                    following = listing.children[index + 1]
                    if node.tag == 'dt' and following.tag == 'dd':
                        links = re.findall(r'href="([^"#]+)\.html"', markup.raw(node))
                        if links:
                            briefs.setdefault(links[-1], re.sub(r'^<dd[^>]*>|</dd>$', '', markup.raw(following)))
            for relation in relations:
                target = relation.get('refid', '')
                full = self.relation_name(relation, definition)
                label = self.overview_argument(full)
                link = '<a href="' + escape(target) + '.html"><code>' + escape(label) + '</code></a>' if self.exists(target) else '<code>' + escape(label) + '</code>'
                qualifiers = ', '.join(value for value in [relation.get('prot'), 'virtual' if relation.get('virt') == 'virtual' else None] if value)
                description = briefs.get(target, escape(text(self.definitions[target].find('briefdescription'))) if target in self.definitions else '')
                rows.append('<tr><td>' + link + '<span class="rodin-inheritance-qualifiers">' + escape(qualifiers) + '</span>' +
                            '<details><summary>Full type</summary>' + self.code(full, definition) + '</details></td><td>' + description + '</td></tr>')
            heading = ('<table class="m-table m-fullwidth rodin-relations"><thead><tr><th>' + first_column + '</th><th>' +
                       ('Responsibility' if tag == 'basecompoundref' else 'Description') + '</th></tr></thead><tbody>')
            visible_rows = 3 if tag == 'derivedcompoundref' else len(rows)
            table = ('<section id="' + section_id + '"><h2><a href="#' + section_id + '">' + title + '</a></h2>' +
                     heading + ''.join(rows[:visible_rows]) + '</tbody></table>')
            if len(rows) > visible_rows:
                table += ('<details class="rodin-derived-disclosure"><summary class="rodin-disclosure-button">'
                          '<span class="rodin-see-more">See more</span><span class="rodin-see-less">See less</span></summary>' +
                          heading + ''.join(rows[visible_rows:]) + '</tbody></table></details>')
            table += '</section>'
            edits.append((section.start, section.end, table))
        if not sections:
            return edits
        # Show one local, bounded diagram. The tables retain every direct edge.
        first = min(sections, key=lambda node: node.start)
        for index, (start, end, replacement) in enumerate(edits):
            if start == first.start:
                edits[index] = (start, end, self.inheritance_diagram(identity, groups) + replacement)
                break
        overrides, inherited = [], {}
        for member in definition.findall('./sectiondef/memberdef'):
            if member.get('prot') == 'private' or not member.get('id', '').startswith(identity + '_1'):
                continue
            targets = [self.member_target(ref.get('refid', '')) for ref in member.findall('reimplements')]
            targets = [target for target in targets if target]
            own = self.member_target(member.get('id', ''))
            if own and targets:
                links = ['<a href="' + target[1] + '">' + escape(scope_name(self.definitions[target[0]].findtext('compoundname'))[1].split('<', 1)[0]) + '::' + escape(member.findtext('name')) + '</a>' for target in targets]
                overrides.append('<tr><td><a href="' + own[1] + '"><code>' + escape(member.findtext('name')) + '()</code></a></td><td>' + ', '.join(links) + '</td></tr>')
        for member in definition.findall('./listofallmembers/member'):
            target = self.member_target(member.get('refid', ''))
            if not target or target[0] == identity or target[0] not in self.definitions or member.get('prot') == 'private':
                continue
            name = member.findtext('name') or ''
            base_name = scope_name(self.definitions[target[0]].findtext('compoundname'))[1].split('<', 1)[0]
            if name == base_name or name.startswith('~'):
                continue
            declared = self.members.get(member.get('refid'))
            if declared is not None and declared.get('kind') == 'friend':
                continue
            label = name + ('()' if declared is not None and declared.get('kind') == 'function' else '')
            inherited.setdefault(target[0], {})[target[1]] = '<li><a href="' + target[1] + '"><code>' + escape(label) + '</code></a>' + (' <span class="rodin-inheritance-qualifiers">protected</span>' if member.get('prot') == 'protected' else '') + '</li>'
        extra = ''
        if overrides:
            extra += '<section id="rodin-overrides"><h2><a href="#rodin-overrides">Overrides</a></h2><table class="m-table m-fullwidth"><thead><tr><th>Member</th><th>Base declaration</th></tr></thead><tbody>' + ''.join(overrides) + '</tbody></table></section>'
        if inherited:
            rows = []
            for base, members in inherited.items():
                label = self.relative_name(self.definitions[base].findtext('compoundname'), definition)
                rows.append('<tr><td><a href="' + base + '.html"><code>' + escape(self.overview_argument(label)) + '</code></a></td><td><details><summary>' + str(len(members)) + ' members</summary><ul class="rodin-inherited-members">' + ''.join(members.values()) + '</ul></details></td></tr>')
            extra += '<section id="rodin-inherited"><h2><a href="#rodin-inherited">Inherited interface</a></h2><table class="m-table m-fullwidth rodin-relations"><thead><tr><th>Declared in</th><th>Members</th></tr></thead><tbody>' + ''.join(rows) + '</tbody></table></section>'
        if extra:
            end = max(section.end for section in sections)
            edits.append((end, end, extra))
        return edits

    def group_entries(self, source):
        markup = Markup(source)
        edits = []
        for listing in [n for n in markup.nodes if n.tag == 'dl' and 'm-doc' in n.attrs.get('class', '').split()]:
            ancestor = listing.parent
            while ancestor and ancestor.tag != 'dl':
                ancestor = ancestor.parent
            if ancestor is not None:
                continue
            pairs = []
            children = listing.children
            for index, node in enumerate(children[:-1]):
                if node.tag != 'dt' or children[index + 1].tag != 'dd':
                    continue
                links = re.findall(r'\b(?:class|struct|union)\s+<a\s+href="([^"#]+)\.html"', markup.raw(node))
                if links and links[-1] in self.definitions:
                    pairs.append((links[-1], node, children[index + 1]))
            represented = {entry[0]: entry for entry in pairs}
            for primary in self.families:
                if primary not in represented:
                    continue
                members = list(dict.fromkeys(entry[0] for entry in pairs if self.primaries.get(entry[0]) == primary))
                if not members:
                    continue
                _, _, description = represented[primary]
                original = markup.raw(description)
                edits.append((description.start, description.end,
                              original[:-5] + '<details class="rodin-specialization-disclosure">'
                              '<summary>See specializations</summary>' + self.family(primary, members) + '</details></dd>'))
                for member in members:
                    _, declaration, description = represented[member]
                    edits.append((declaration.start, description.end, ''))
            # Some extraction configurations provide specializations without a
            # primary page. Group these by fully qualified family name without
            # inventing a declaration or linking to a nonexistent primary.
            missing = {}
            for identity, declaration, description in pairs:
                scope, name = scope_name(self.definitions[identity].findtext('compoundname'))
                if '<' not in name or identity in self.primaries:
                    continue
                family = (scope + '::' if scope else '') + name.split('<', 1)[0]
                missing.setdefault(family, []).append((identity, declaration, description))
            for family, entries in missing.items():
                identity, declaration, description = entries[0]
                label = scope_name(family)[1]
                content = ('<dt>' + escape(label) + ' specializations</dt><dd>'
                           '<details class="rodin-specialization-disclosure"><summary>See specializations</summary>' +
                           self.family(None, [entry[0] for entry in entries]) + '</details></dd>')
                edits.append((declaration.start, description.end, content))
                for _, declaration, description in entries[1:]:
                    edits.append((declaration.start, description.end, ''))
        return replace_ranges(source, edits)

    def collapse_declarations(self, source):
        # A stable textual limit makes the overview independent of viewport
        # width. Full declarations, including all links, stay in the disclosure.
        maximum_preview_length = 120
        markup, edits = Markup(source), []
        for listing in [node for node in markup.nodes if node.tag == 'dl' and 'm-doc' in node.attrs.get('class', '').split()]:
            ancestor = listing.parent
            while ancestor and ancestor.tag not in ('dl', 'details'):
                ancestor = ancestor.parent
            if ancestor is not None:
                continue
            for index, declaration in enumerate(listing.children[:-1]):
                description = listing.children[index + 1]
                if declaration.tag != 'dt' or description.tag != 'dd':
                    continue
                plain = re.sub(r'\s+', ' ', markup.content(declaration))
                if len(plain) <= maximum_preview_length:
                    continue
                raw = markup.raw(declaration)
                symbol = re.search(r'\b(class|struct|union)\s+(<a\b[^>]*>.*?</a>)', raw, re.S)
                if symbol:
                    short = symbol[1] + ' ' + symbol[2]
                else:
                    function = re.search(r'(<a\b[^>]*>[^<]*(?:<wbr\s*/?>[^<]*)*</a>)\s*\(', raw)
                    if not function:
                        continue
                    short = function[1]
                opening = raw[:raw.index('>') + 1]
                complete = raw[len(opening):-5]
                brief = markup.raw(description)
                content = (opening + short + '</dt>' + brief[:-5] +
                           '<details class="rodin-declaration-disclosure"><summary>See declaration</summary>'
                           '<div class="rodin-overview-declaration">' + complete + '</div></details></dd>')
                edits.append((declaration.start, description.end, content))
        return replace_ranges(source, edits)

    def deduction_guides(self, identity, source):
        markup, edits = Markup(source), []
        for section in [node for node in markup.nodes if node.tag == 'section' and 'm-doc-details' in node.attrs.get('class', '').split()]:
            if any(node.tag == 'details' and 'rodin-declaration-disclosure' in node.attrs.get('class', '').split()
                   for node in markup.nodes if section.start < node.start < section.end):
                continue
            member = self.members.get(identity + '_1' + section.attrs.get('id', ''))
            # Doxygen represents deduction guides as functions with no return
            # type and a trailing deduced class type. Ordinary trailing-return
            # functions have a nonempty type (normally auto).
            if member is None or member.get('kind') != 'function' or text(member.find('type')):
                continue
            name = member.findtext('name') or ''
            arguments = member.findtext('argsstring') or ''
            result = re.search(r'\)\s*->\s*(.+)$', arguments, re.S)
            if not result or scope_name(result[1].strip())[1].split('<', 1)[0].strip() != name:
                continue
            heading = next((node for node in markup.nodes if node.tag == 'h3' and section.start < node.start < section.end), None)
            if heading is None:
                continue
            raw = markup.raw(heading)
            link = re.search(r'<a\b[^>]*class="m-doc-self"[^>]*>.*?</a>', raw, re.S)
            if not link:
                continue
            declaration = raw[raw.index('>') + 1:-5]
            disclosure = ('<details class="rodin-declaration-disclosure"><summary>See declaration</summary>'
                          '<div class="rodin-overview-declaration">' + declaration + '</div></details>')
            # Keep the source-header annotation immediately below the heading.
            following = next((node for node in heading.parent.children if node.start >= heading.end), None)
            position = following.end if following and 'rodin-reference-include' in following.attrs.get('class', '').split() else heading.end
            edits.append((heading.start, heading.end, '<h3>' + link[0] + '</h3>'))
            edits.append((position, position, disclosure))
        return replace_ranges(source, edits)

    def source_headers(self, source):
        markup, edits = Markup(source), []
        for node in markup.nodes:
            if 'm-doc-include' not in node.attrs.get('class', '').split():
                continue
            anchors = [child for child in markup.nodes if child.tag == 'a' and node.start < child.start < node.end]
            if anchors:
                anchor = anchors[0]
                header = '<a href="' + escape(anchor.attrs.get('href', ''), quote=True) + '"><code>' + escape(markup.content(anchor).removeprefix('<').removesuffix('>')) + '</code></a>'
            else:
                header = '<code>' + escape(markup.content(node).removeprefix('#include').strip().removeprefix('<').removesuffix('>')) + '</code>'
            replacement = '<p class="rodin-reference-include">Defined in ' + header + '</p>'
            ancestor = node.parent
            while ancestor and ancestor.tag not in ('h1', 'h2', 'h3'):
                ancestor = ancestor.parent
            if ancestor:
                # Source location is supporting prose, not part of the signature.
                edits.append((node.start, node.end, ''))
                edits.append((ancestor.end, ancestor.end, replacement))
            else:
                edits.append((node.start, node.end, replacement))
        return replace_ranges(source, edits)

    def remove_related_see_also(self, identity, source):
        definition = self.definitions[identity]
        if definition.get('kind') == 'namespace':
            return source
        # Follow ancestors and descendants separately: sharing a base does not
        # make sibling operations redundant cross-references.
        related = {identity, self.primaries.get(identity, identity)}
        related.update(self.families.get(identity, []))
        seeds = list(related)
        for direction in ['basecompoundref', 'derivedcompoundref']:
            visited, pending = set(), list(seeds)
            while pending:
                current = pending.pop()
                if current in visited or current not in self.definitions:
                    continue
                visited.add(current)
                for target in self.inheritance_edges[direction].get(current, set()):
                    related.add(target)
                    primary = self.primaries.get(target, target)
                    related.add(primary)
                    pending.extend([target, primary])
        covered = {target + '.html' for target in related}
        markup, edits = Markup(source), []
        for aside in [node for node in markup.nodes if node.tag == 'aside']:
            heading = next((child for child in aside.children if child.tag in ['h3', 'h4', 'h5']), None)
            if heading is None or markup.content(heading).lower().rstrip(':') != 'see also':
                continue
            ancestor = aside.parent
            while ancestor and 'm-doc-details' not in ancestor.attrs.get('class', '').split():
                ancestor = ancestor.parent
            if ancestor is not None:
                continue
            removals = []
            for paragraph in aside.children:
                if paragraph.tag != 'p':
                    continue
                raw = markup.raw(paragraph)
                links = re.findall(r'<a\s[^>]*href="([^"]+)"[^>]*>.*?</a>', raw, re.S)
                leftover = re.sub(r'<a\b[^>]*>.*?</a>', '', raw, flags=re.S)
                leftover = unescape(re.sub(r'<[^>]+>', '', leftover))
                if links and all(link in covered for link in links) and re.fullmatch(r'[\s,;|]*', leftover):
                    removals.append((paragraph.start - aside.start, paragraph.end - aside.start, ''))
            if not removals:
                continue
            updated = replace_ranges(markup.raw(aside), removals)
            remainder = Markup(updated)
            bodies = [node for node in remainder.nodes if node.tag == 'aside'][0].children
            # Preserve prose, unrelated links, technical groups and member-level
            # contracts. Only remove the container if its heading is all left.
            if all(child.tag in ['h3', 'h4', 'h5'] or not remainder.content(child) for child in bodies):
                updated = ''
            edits.append((aside.start, aside.end, updated))
        return replace_ranges(source, edits)

    def render(self, identity, source):
        definition = self.definitions[identity]
        markup = Markup(source)
        edits = self.member_signatures(definition, markup)
        edits += self.inheritance(identity, markup)
        if definition.get('kind') != 'namespace':
            edits += self.class_header(identity, markup)
        source = replace_ranges(source, edits)
        if identity in self.families:
            markup = Markup(source)
            tables = [n for n in markup.nodes if n.tag == 'table' and
                      [markup.content(c) for c in markup.nodes if c.tag == 'th' and n.start < c.start < n.end] == ['Specialization', 'Description']]
            family = self.family(identity)
            # Replace only an existing table explicitly belonging to this family.
            matching = [t for t in tables if identity + '.html' in markup.raw(t) or any(m + '.html' in markup.raw(t) for m in self.families[identity])]
            if matching:
                source = replace_ranges(source, [(matching[0].start, matching[0].end, family)])
            else:
                header = next(n for n in markup.nodes if n.tag == 'h1')
                # Keep the new section outside the declaration's preformatted block.
                declaration = next(n for n in markup.nodes if n.tag == 'pre' and n.start >= header.end)
                source = replace_ranges(source, [(declaration.end, declaration.end,
                    '<section id="rodin-specializations"><h2><a href="#rodin-specializations">Specializations</a></h2>' + family + '</section>')])
        source = source.replace('>Base classes</a>', '>Inherits from</a>')
        source = self.source_headers(self.remove_related_see_also(identity, self.collapse_declarations(self.group_entries(source))))
        return self.deduction_guides(identity, source)


def prepare_reference_headers(xml):
    """Give included public API headers real file pages, even without @file prose."""
    xml = Path(xml)
    targets = set()
    for path in xml.glob('*.xml'):
        definition = ET.parse(path).find('compounddef')
        if definition is not None and definition.get('kind') in ('class', 'struct', 'union'):
            includes = definition.find('includes')
            if includes is not None and includes.get('refid'):
                targets.add(includes.get('refid'))
    for identity in targets:
        path = xml / (identity + '.xml')
        if not path.is_file():
            continue
        tree = ET.parse(path)
        definition = tree.find('compounddef')
        brief = definition.find('briefdescription')
        if definition.get('kind') == 'file' and brief is not None and not text(brief) and not text(definition.find('detaileddescription')):
            ET.SubElement(brief, 'para').text = 'Declarations in ' + definition.findtext('compoundname') + '.'
            tree.write(path, encoding='utf-8', xml_declaration=True)


def render_warning_notes(site):
    """Use the native m.css warning component for Doxygen @warning notes."""
    count = 0
    for path in Path(site).rglob('*.html'):
        source = path.read_text()
        if 'm-danger' not in source or 'Warning' not in source:
            continue
        markup = Markup(source)
        edits = []
        for note in markup.nodes:
            classes = note.attrs.get('class', '').split()
            if note.tag != 'aside' or 'm-note' not in classes or 'm-danger' not in classes:
                continue
            heading = next((child for child in note.children if child.tag == 'h4'), None)
            if heading is None or markup.content(heading) != 'Warning':
                continue
            end = source.index('>', note.start) + 1
            opening = source[note.start:end]
            classes = ['m-warning' if value == 'm-danger' else value for value in classes]
            opening = re.sub(r"""(?<![\w:-])class\s*=\s*(["'])(.*?)\1""",
                             lambda match: 'class=' + match.group(1) + ' '.join(classes) + match.group(1),
                             opening, count=1)
            edits.append((note.start, end, opening))
        if edits:
            path.write_text(replace_ranges(source, edits))
            count += len(edits)
    return count


def render_references(xml, site):
    layout = ReferenceLayout(xml, site)
    count = 0
    for identity in layout.definitions:
        path = Path(site) / (identity + '.html')
        if not path.is_file():
            continue
        source = path.read_text()
        rendered = layout.render(identity, source)
        path.write_text(rendered)
        count += 1
    # Header annotations also occur on group and file pages, which do not use
    # the class/namespace declaration layouts above.
    for path in Path(site).glob('*.html'):
        source = path.read_text()
        if 'm-doc-include' not in source:
            continue
        rendered = layout.source_headers(source)
        if rendered != source:
            path.write_text(rendered)
            count += 1
    return count

"""Keep reference links consistent with the pages actually published by m.css."""
from html import escape
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET

from documentation_reference import replace_ranges, text


class Links(HTMLParser):
    """Read destinations and anchor spans without building a full DOM per page."""
    def __init__(self, source):
        super().__init__(convert_charrefs=False)
        self.source = source
        self.lines = [0] + [match.end() for match in re.finditer('\n', source)]
        self.references, self.anchors, self.identifiers = [], [], set()
        self.anchor = None
        self.feed(source)

    def source_offset(self):
        line, column = self.getpos()
        return self.lines[line - 1] + column

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get('id'):
            self.identifiers.add(attrs['id'])
        if tag == 'a' and attrs.get('name'):
            self.identifiers.add(attrs['name'])
        for attribute in ['href', 'src', 'poster']:
            if attrs.get(attribute):
                self.references.append(attrs[attribute])
        if tag == 'a' and attrs.get('href'):
            start = self.source_offset()
            self.anchor = (start, self.source.index('>', start) + 1, attrs['href'])

    def handle_endtag(self, tag):
        if tag == 'a' and self.anchor is not None:
            start, content_start, href = self.anchor
            end = self.source_offset()
            self.anchors.append((start, self.source.index('>', end) + 1,
                                 self.source[content_start:end], href))
            self.anchor = None


def local_target(page, reference):
    url = urlsplit(reference)
    if url.scheme or url.netloc or reference.startswith('/'):
        return None
    if not url.path:
        return page if url.fragment else None
    target = (page.parent / unquote(url.path)).resolve()
    return target / 'index.html' if target.is_dir() else target


def unlink_unpublished_references(xml, site):
    """Retain scope/type text without linking undocumented or private records.

    Doxygen can link a forward declaration or enclosing private scope even when
    m.css intentionally emits no page for it. Only those known XML records are
    unlinked; missing public documentation and stale hand-written URLs remain
    errors for the complete-site check.
    """
    xml, site = Path(xml), Path(site).resolve()
    unpublished = set()
    for path in xml.glob('*.xml'):
        definition = ET.parse(path).find('compounddef')
        if definition is None or definition.get('kind') not in ['class', 'struct', 'union', 'namespace']:
            continue
        if definition.get('prot') == 'private' or not (
                text(definition.find('briefdescription')) or text(definition.find('detaileddescription'))):
            unpublished.add(path.stem + '.html')
    count = 0
    for page in site.rglob('*.html'):
        source = page.read_text()
        markup = Links(source)
        edits = []
        for start, end, content, reference in markup.anchors:
            target = local_target(page, reference)
            if target is None or target.exists() or target.parent != site or target.name not in unpublished:
                continue
            edits.append((start, end, content))
        if edits:
            page.write_text(replace_ranges(source, edits))
            count += len(edits)
    return count


def publish_enum_anchors(xml, site):
    """Give enum values without detailed sections anchors in their declarations."""
    values = set()
    for path in Path(xml).glob('*.xml'):
        for value in ET.parse(path).getroot().iter('enumvalue'):
            identity = value.get('id', '')
            if '_1' in identity:
                values.add(identity.rsplit('_1', 1)[1])
    count = 0
    for page in Path(site).rglob('*.html'):
        source = page.read_text()
        markup = Links(source)
        edits = []
        for start, _, _, reference in markup.anchors:
            if not reference.startswith('#'):
                continue
            fragment = unquote(reference[1:])
            if fragment not in values or fragment in markup.identifiers:
                continue
            edits.append((start + 2, start + 2, ' id="' + escape(fragment, quote=True) + '"'))
            markup.identifiers.add(fragment)
        if edits:
            page.write_text(replace_ranges(source, edits))
            count += len(edits)
    return count


def verify_local_links(site):
    """Reject missing local page, asset and fragment destinations throughout the site."""
    site = Path(site).resolve()
    pages = {}
    for page in site.rglob('*.html'):
        markup = Links(page.read_text())
        pages[page] = (markup.references, markup.identifiers)
    missing = []
    for page, (references, _) in pages.items():
        for reference in references:
            target = local_target(page, reference)
            if target is not None and not target.exists():
                missing.append(f'{page.relative_to(site)} -> {reference}')
            elif target in pages:
                fragment = unquote(urlsplit(reference).fragment)
                if fragment and fragment != '!' and fragment not in pages[target][1]:
                    missing.append(f'{page.relative_to(site)} -> {reference}')
    if missing:
        raise ValueError('Missing documentation destinations:\n' + '\n'.join(sorted(set(missing))))

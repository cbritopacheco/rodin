"""Generate benchmark-only search data for the pinned m.css search engine."""
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import runpy

ROOT = Path(__file__).resolve().parents[1]


class PlainText(HTMLParser):
    def __init__(self, markup):
        super().__init__()
        self.parts = []
        self.feed(markup)

    def handle_data(self, value):
        self.parts.append(value)

    def text(self):
        return ' '.join(' '.join(self.parts).split())


def build_search(pages, modules, assets):
    vendor = ROOT / 'third-party/m.css/documentation'
    api = runpy.run_path(str(vendor / '_search.py'))
    trie, mapping = api['Trie'](), api['ResultMap']()
    entries = []
    for path, page in sorted(pages.items()):
        if not path or path == 'Rodin':
            continue
        module = modules.get(path)
        title = module['title'] if module else Path(path).name
        source = page['description'] or (module['description'] if module else '')
        first = re.search(r'<p(?:\s[^>]*)?>.*?</p>', source, re.S)
        brief = PlainText(module['brief'] if module else first[0] if first else source).text()
        if len(brief) > 280:
            brief = brief[:280].rsplit(' ', 1)[0]
            for opening, closing in [('\\(', '\\)'), ('\\[', '\\]')]:
                if brief.rfind(opening) > brief.rfind(closing):
                    brief = brief[:brief.rfind(opening)].rstrip()
            brief += '…'
        display = path.replace('Rodin/', 'Modules/', 1)
        entry = {'path': path, 'title': title, 'url': path + '/index.html',
                 'display': display, 'brief': brief,
                 'text': PlainText(source).text(),
                 'target': page['record']['target'] if page['record'] else '',
                 'kind': 'Module' if module else 'Benchmark' if page['record'] else 'Collection'}
        entries.append(entry)
        kind = ['Module', 'Benchmark', 'Collection'].index(entry['kind']) + 1
        result = mapping.add(display, entry['url'], flags=api['ResultFlag'](kind << 4))
        # Exact scope lookups support ranked description/multiword matches in the
        # extension. Suffixes retain the documentation's prefix completion.
        parts = path.split('/')
        for depth in range(len(parts)):
            trie.insert('/'.join(parts[depth:]).lower(), result)
        trie.insert(display.lower(), result)
        trie.insert(path.replace('/', '::').lower(), result)
    trie.sort(mapping)
    data = api['serialize_search_data'](
        api['Serializer'](file_offset_bytes=4, result_id_bytes=2, name_size_bytes=2),
        trie, mapping, [(api['CssClass'].INFO, name) for name in ['Module', 'Benchmark', 'Collection']], len(entries))
    engine = (vendor / 'search.js').read_text() + '\n' + (ROOT / 'dev/benchmarks/search.js').read_text()
    engine_name = 'search-' + hashlib.sha256(engine.encode()).hexdigest()[:12] + '.js'
    payload = ('window.RodinBenchmarkSearch=' + json.dumps(entries, ensure_ascii=True) + ';\n').encode() + api['base85encode_search_data'](data)
    data_name = 'searchdata-' + hashlib.sha256(payload).hexdigest()[:12] + '.js'
    (assets / engine_name).write_text(engine)
    (assets / data_name).write_bytes(payload)
    return {'engine': engine_name, 'data': data_name}

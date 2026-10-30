"""General class ranking is global, capped after promotion, and scope-safe."""
from pathlib import Path
import json
import subprocess
import tempfile
import unittest

from documentation_search import assemble_search, prepare_class_search, primary_classes


class SearchRankingTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.xml, self.site = self.root / 'xml', self.root / 'html'
        self.xml.mkdir()
        self.site.mkdir()
        for identity, name, kind in [('first', 'Rodin::First', 'class'),
                                      ('second', 'Rodin::Nested::Second', 'struct'),
                                      ('special', 'Rodin::First&lt;Real&gt;', 'class'),
                                      ('nested', 'Rodin::First&lt;Real&gt;::Nested', 'class'),
                                      ('space', 'Rodin::Space', 'namespace')]:
            (self.xml / (identity + '.xml')).write_text('<doxygen><compounddef id="' + identity + '" kind="' + kind + '"><compoundname>' + name + '</compoundname><briefdescription/><detaileddescription/></compounddef></doxygen>')
            (self.site / (identity + '.html')).touch()

    def test_catalogue_covers_classes_and_structs_in_all_scopes(self):
        records = primary_classes(self.xml, self.site)
        self.assertEqual([entry['url'] for entry in records], ['first.html', 'second.html'])
        self.assertEqual(records[1]['keys'], ['rodin::nested::second', 'nested::second', 'second'])

    def test_undocumented_family_landing_pages_enter_the_index(self):
        import xml.etree.ElementTree as ET
        self.assertEqual(prepare_class_search(self.xml), 1)
        brief = ET.parse(self.xml / 'first.xml').findtext('./compounddef/briefdescription/para')
        self.assertEqual(brief, 'Primary template for First.')
        self.assertIsNone(ET.parse(self.xml / 'second.xml').find('./compounddef/briefdescription/para'))
        self.assertEqual(prepare_class_search(self.xml), 0)

    def test_non_rendered_types_are_not_promoted(self):
        (self.site / 'second.html').unlink()
        self.assertEqual(len(primary_classes(self.xml, self.site)), 1)

    def test_bundle_keeps_upstream_and_fills_catalogue(self):
        (self.site / 'first.html').write_text('<script src="search-v2.js"></script>')
        upstream = self.root / 'upstream.js'
        upstream.write_text('/* upstream licence */\nvar Search = {};')
        extension = self.root / 'extension.js'
        extension.write_text('const classes = __RODIN_SEARCH_CLASSES__;')
        self.assertEqual(assemble_search(self.xml, self.site, upstream, extension, 2), 2)
        result = (self.site / 'search-v2.js').read_text()
        self.assertTrue(result.startswith(upstream.read_text()))
        self.assertIn('Rodin::Nested::Second', result)
        self.assertNotIn('__RODIN_SEARCH_CLASSES__', result)
        self.assertRegex((self.site / 'first.html').read_text(), r'search-v2-[0-9a-f]{12}\.js')

    def test_browser_ranking_and_specific_queries(self):
        extension = (Path(__file__).parents[1] / 'doc' / 'rodin-search.js').read_text()
        catalogue = [{'name': 'Rodin::Field', 'url': 'field.html', 'keys': ['rodin::field', 'field']},
                     {'name': 'Rodin::FieldLoader', 'url': 'loader.html', 'keys': ['rodin::fieldloader', 'fieldloader']}]
        harness = r'''
const assert = require('node:assert/strict');
const primary = {name: 'Rodin::Field', url: 'field.html', flags: 0, suffixLength: 0};
const loader = {name: 'Rodin::FieldLoader', url: 'loader.html', flags: 0, suffixLength: 0};
const special = {name: 'Rodin::Field<Real>', url: 'special.html', flags: 0, suffixLength: 6};
const member = {name: 'Rodin::Field::get', url: 'field.html#get', flags: 0, suffixLength: 0};
var Search = {
  trie: true, maxResults: 2, searchString: 'old', searchStack: [19, 21],
  toUtf8: value => Buffer.from(value).toString('latin1'),
  fromUtf8: value => Buffer.from(value, 'latin1').toString(),
  search(query) {
    this.searchString = query;
    this.searchStack = [19, query.length];
    let results = query === 'Rodin::Field' ? [primary, member] :
      query === 'Rodin::FieldLoader' ? [loader] :
      query.includes('<') ? [special] :
      query.endsWith('::') ? [member] : [special, member];
    return [results.map(r => ({...r})), 'completion'];
  }
};
'''
        assertions = r'''
let [results, completion] = Search.search('Fi');
assert.deepEqual(results.map(r => r.url), ['field.html', 'loader.html']);
assert.equal(results[0].suffixLength, 3);
assert.equal(Search.searchString, 'Fi');
assert.deepEqual(Search.searchStack, [19, 2]);
assert.equal(completion, 'completion');
assert.deepEqual(Search.search('Field<Real>')[0].map(r => r.url), ['special.html']);
assert.deepEqual(Search.search('Field::')[0].map(r => r.url), ['field.html#get']);
Search.maxResults = 5;
results = Search.search('Field')[0];
assert.equal(results.filter(r => r.url === 'field.html').length, 1);
assert.equal(results[0].url, 'field.html');
assert.equal(results[1].url, 'loader.html');
assert.ok(results.some(r => r.url === 'special.html'));
console.log('Search ranking, full-result promotion, highlights, query restoration and explicit queries pass.');
'''
        script = self.root / 'ranking.cjs'
        script.write_text(harness + extension.replace('__RODIN_SEARCH_CLASSES__', json.dumps(catalogue)) + assertions)
        subprocess.run(['node', str(script)], check=True, capture_output=True, text=True)


if __name__ == '__main__':
    unittest.main()

"""Missing private scopes may be plain text; missing public pages must fail CI."""
from pathlib import Path
import tempfile
import unittest

from documentation_links import publish_enum_anchors, unlink_unpublished_references, verify_local_links


class DocumentationLinksTest(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        self.xml, self.site = self.root / 'xml', self.root / 'site'
        self.xml.mkdir()
        self.site.mkdir()

    def definition(self, name, brief='', protection='public'):
        (self.xml / (name + '.xml')).write_text(
            f'<doxygen><compounddef kind="class" prot="{protection}">'
            f'<briefdescription>{brief}</briefdescription><detaileddescription/>'
            '</compounddef></doxygen>')

    def test_unpublished_scope_keeps_text_and_existing_public_links(self):
        self.definition('private', 'Implementation state', 'private')
        self.definition('forward')
        self.definition('public', 'Documented API')
        page = self.site / 'page.html'
        page.write_text('<div class="m-breadcrumb"><a href="private.html">Parent</a></div>'
                        '<a href="forward.html"><code>Forward</code></a>'
                        '<a href="public.html">Public</a>')
        (self.site / 'public.html').write_text('Public API')
        self.assertEqual(unlink_unpublished_references(self.xml, self.site), 2)
        self.assertIn('<div class="m-breadcrumb">Parent</div>', page.read_text())
        self.assertIn('<code>Forward</code>', page.read_text())
        self.assertIn('href="public.html"', page.read_text())
        verify_local_links(self.site)
        self.assertEqual(unlink_unpublished_references(self.xml, self.site), 0)

    def test_missing_public_page_or_unknown_link_cannot_be_silently_removed(self):
        self.definition('public', 'Documented API')
        page = self.site / 'page.html'
        page.write_text('<a href="public.html">Public</a><a href="stale.html">Old</a>')
        self.assertEqual(unlink_unpublished_references(self.xml, self.site), 0)
        with self.assertRaisesRegex(ValueError, 'Missing documentation destinations'):
            verify_local_links(self.site)

    def test_nested_benchmark_paths_and_assets_are_checked(self):
        page = self.site / 'benchmarks' / 'module' / 'index.html'
        page.parent.mkdir(parents=True)
        (self.site / 'index.html').write_text('Docs')
        page.write_text('<a href="../../index.html">Docs</a><script src="../../missing.js"></script>')
        with self.assertRaisesRegex(ValueError, 'missing.js'):
            verify_local_links(self.site)
        (self.site / 'missing.js').write_text('')
        verify_local_links(self.site)

    def test_enum_declaration_gets_real_value_anchor_without_duplicate_ids(self):
        (self.xml / 'enum.xml').write_text('<doxygen><compounddef><sectiondef>'
            '<memberdef kind="enum"><enumvalue id="enum_1avalue"><name>Ready</name>'
            '</enumvalue></memberdef></sectiondef></compounddef></doxygen>')
        page = self.site / 'enum.html'
        page.write_text('<a href="#avalue">Ready</a><a href="#avalue">Ready again</a>')
        self.assertEqual(publish_enum_anchors(self.xml, self.site), 1)
        self.assertEqual(page.read_text().count('id="avalue"'), 1)
        self.assertEqual(publish_enum_anchors(self.xml, self.site), 0)
        verify_local_links(self.site)

    def test_missing_fragment_cannot_pass_when_its_page_exists(self):
        (self.site / 'index.html').write_text('<a href="target.html#missing">Target</a>')
        (self.site / 'target.html').write_text('<h1 id="present">Present</h1>')
        with self.assertRaisesRegex(ValueError, 'target.html#missing'):
            verify_local_links(self.site)

    def test_external_links_and_interactive_fragments_are_not_local_files(self):
        (self.site / 'index.html').write_text('<a href="https://example.test/no.html">External</a>'
            '<a href="#!">Menu</a><a href="#search">Search</a><div id="search"></div><img src="data:image/png;base64,AA">')
        verify_local_links(self.site)


if __name__ == '__main__':
    unittest.main()

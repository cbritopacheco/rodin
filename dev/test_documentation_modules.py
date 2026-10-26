"""Concise module overviews preserve rendered formulas, examples and anchors."""
from types import SimpleNamespace
import unittest
from pathlib import Path
import tempfile

from build_documentation import Representative, scope_html, overview, clear_generated_pages, version_stylesheets, prepare_inline_images, assemble_stylesheets


class ModuleOverviewTest(unittest.TestCase):
    def test_rebuild_removes_stale_pages_but_preserves_formula_cache_and_assets(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            xml, html = root / 'xml', root / 'html'
            xml.mkdir()
            html.mkdir()
            (xml / 'source-headers.xml').write_text('obsolete page')
            (html / 'source-headers.html').write_text('obsolete page')
            (html / 'rodin-navigation.js').write_text('sidebar asset')
            (root / 'm.math.cache').write_text('cached formulas')
            clear_generated_pages(xml, html)
            self.assertFalse((xml / 'source-headers.xml').exists())
            self.assertFalse((html / 'source-headers.html').exists())
            self.assertEqual((html / 'rodin-navigation.js').read_text(), 'sidebar asset')
            self.assertEqual((root / 'm.math.cache').read_text(), 'cached formulas')

    def test_stylesheet_url_changes_with_content_but_is_stable_for_identical_builds(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / 'theme.css'
            source.write_text('code { color: blue; }')
            remote = 'https://fonts.example/style.css'
            first = version_stylesheets([str(source), remote], root / 'output')
            self.assertEqual(first, version_stylesheets([str(source), remote], root / 'output'))
            self.assertEqual(first[1], remote)
            source.write_text('code { color: red; }')
            second = version_stylesheets([str(source), remote], root / 'output')
            self.assertNotEqual(first[0], second[0])
            self.assertEqual(Path(second[0]).read_text(), source.read_text())

    def test_upstream_stylesheets_keep_variables_and_source_order(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            theme, layout = root / 'theme.css', root / 'layout.css'
            theme.write_text('/* Upstream license */ :root { --color: blue; }')
            layout.write_text('body { color: var(--color); }')
            bundle = Path(assemble_stylesheets([theme, layout], root / 'output'))
            self.assertEqual(bundle.read_text(), theme.read_text() + '\n' + layout.read_text())

    def test_inline_avatar_is_staged_without_losing_its_alt_text(self):
        import xml.etree.ElementTree as ET
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            xml, images = root / 'xml', root / 'doc/Images'
            xml.mkdir()
            images.mkdir(parents=True)
            (images / 'Maintainer.png').write_bytes(b'profile image')
            page = xml / 'page.xml'
            page.write_text('<doxygen><compounddef><image type="html" inline="yes" '
                'name="Images/Maintainer.png" alt="Maintainer profile"/></compounddef></doxygen>')
            prepare_inline_images(xml, root / 'doc')
            image = ET.parse(page).find('.//image')
            self.assertEqual(image.get('name'), 'Maintainer.png')
            self.assertEqual(image.get('alt'), 'Maintainer profile')
            self.assertEqual((xml / 'Maintainer.png').read_bytes(), b'profile image')

    def test_display_formula_keeps_balanced_svg_markup(self):
        formula = '<div class="m-math"><svg viewBox="0 0 3 2"><defs><path id="g"/></defs><use xlink:href="#g"/></svg></div>'
        self.assertEqual(Representative('<p>Introduction.</p>' + formula + '<p>After.</p>').html(), formula)

    def test_code_usage_and_void_tags_are_preserved(self):
        code = '<pre class="m-code"><span>Mesh mesh;</span><br/>next line</pre>'
        self.assertEqual(Representative('<p>Intro.</p>' + code).html(), code)

    def test_missing_representative_is_an_error(self):
        with self.assertRaisesRegex(ValueError, 'formula or highlighted usage'):
            Representative('<p>Text alone.</p>').html()

    def test_scoping_preserves_svg_attributes_and_api_links(self):
        source = '<h3 id="usage">Use</h3><svg viewBox="0 0 1 1"><path id="glyph"/><use xlink:href="#glyph"/></svg>'
        source += '<a href="#usage">local</a><a href="#functions">API</a><a href="elsewhere.html#x">other</a>'
        result = scope_html(source, 'geometry-', 'geometry.html')
        self.assertIn('id="geometry-usage"', result)
        self.assertIn('xlink:href="#geometry-glyph"', result)
        self.assertIn('viewBox="0 0 1 1"', result)
        self.assertIn('href="#geometry-usage"', result)
        self.assertIn('href="geometry.html#functions"', result)
        self.assertIn('href="elsewhere.html#x"', result)

    def test_brief_formula_ids_and_api_anchors_are_scoped(self):
        brief = '<svg><path id="g"/><use href="#g"/></svg><a href="#api">API</a>'
        result = overview(SimpleNamespace(brief=brief, description='<p>Numerical operations.</p><pre class="m-code">use();</pre>'), 'math')
        self.assertIn('id="math-brief-g"', result['brief'])
        self.assertIn('href="#math-brief-g"', result['brief'])
        self.assertIn('href="math.html#api"', result['brief'])

    def test_overview_contains_only_summary_and_representative(self):
        description = '<p>Detailed explanation.</p><h2 id="method">Method</h2><div class="m-math"><svg><path id="g"/></svg></div>'
        result = overview(SimpleNamespace(brief='Geometry summary.', description=description), 'geometry')
        self.assertEqual(result['brief'], 'Geometry summary.')
        self.assertEqual(set(result), {'brief', 'representative', 'description'})
        self.assertIn('id="geometry-representative-g"', result['representative'])
        self.assertEqual(result['description'], '<p>Detailed explanation.</p>')
        self.assertNotIn('Detailed explanation.', result['representative'])
        self.assertNotIn('Method', result['description'])


if __name__ == '__main__':
    unittest.main()

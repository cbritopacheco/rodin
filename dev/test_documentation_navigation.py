"""Regression checks for public modules and discoverable specialization references."""
from pathlib import Path
import tempfile
import unittest

from check_documentation_navigation import verify


class DocumentationNavigationTest(unittest.TestCase):
    def fixture(self, root):
        site, xml = root / 'html', root / 'xml'
        site.mkdir()
        xml.mkdir()
        config = root / 'conf.py'
        config.write_text("RODIN_MODULES = [('Core', [('Math', 'namespace_math')])]\n"
                          "RODIN_MODULE_CATEGORIES = {'Core': 'Numerical foundations.'}\n")
        (site / 'modules.html').write_text('<h2>Core</h2><p data-role="category-description" data-category="Core">'
            'Numerical foundations.</p><section id="namespace_math" data-module="namespace_math">'
            '<h3><a href="namespace_math.html">Math</a></h3>'
            '<p data-role="brief"><em>Mathematics.</em></p>'
            '<div data-role="representative"><pre class="m-code">dot(a,b)</pre></div>'
            '<div data-role="description"><p>Numerical operations.</p></div>'
            '</section><a href="configuration.html">Configuration</a>')
        (site / 'namespace_math.html').write_text('<h1>Math <span class="m-thin">module</span></h1>')
        (xml / 'namespace_math.xml').write_text('<doxygen><compounddef><briefdescription>Mathematics</briefdescription>'
            '<detaileddescription><para>Scalar algebra.</para><formula>x+y</formula></detaileddescription></compounddef></doxygen>')
        (site / 'configuration.html').write_text('<a href="_configure_8h_8in.html">Configure.h</a>')
        (site / '_configure_8h_8in.html').write_text('Build configuration macros')
        (xml / '_configure_8h_8in.xml').write_text('<doxygen><compounddef><sectiondef><memberdef kind="define">'
            '<name>RODIN_USE_MPI</name></memberdef></sectiondef></compounddef></doxygen>')
        (xml / 'index.xml').write_text('<doxygenindex><compound kind="class" refid="class_foo"><name>Rodin::Foo</name>'
            '</compound><compound kind="group" refid="group_foo"><name>FooSpecializations</name></compound></doxygenindex>')
        (xml / 'group_foo.xml').write_text('<doxygen><compounddef><title>Foo Template Specializations</title></compounddef></doxygen>')
        (xml / 'class_foo.xml').write_text('<doxygen><compounddef><templateparamlist><param>T</param></templateparamlist>'
            '<detaileddescription><simplesect kind="see"><ref refid="group_foo">Foo Template Specializations</ref>'
            '</simplesect></detaileddescription></compounddef></doxygen>')
        (site / 'class_foo.html').write_text('<a href="group_foo.html">Foo Template Specializations</a>')
        (site / 'group_foo.html').write_text('Specialization reference')
        (site / 'files.html').write_text('<a href="file_header.html">Header.h</a>')
        (site / 'file_header.html').write_text('Header documentation')
        (xml / 'file_header.xml').write_text('<doxygen><compounddef><briefdescription><para>Header.</para>'
            '</briefdescription><detaileddescription/></compounddef></doxygen>')
        index = xml / 'index.xml'
        index.write_text(index.read_text().replace('</doxygenindex>',
            '<compound kind="file" refid="file_header"><name>Header.h</name></compound></doxygenindex>'))
        (site / 'annotated.html').write_text('<a href="class_foo.html">Foo</a>')
        (site / 'rodin-navigation.js').write_text('// Sidebar script')
        for page in site.glob('*.html'):
            page.write_text('<script src="rodin-navigation.js" defer></script>' + page.read_text())
        return site, xml, config

    def test_sidebar_script_must_be_published(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            (args[0] / 'rodin-navigation.js').unlink()
            with self.assertRaisesRegex(ValueError, 'Missing sidebar navigation asset'):
                verify(*args)

    def test_module_pages_must_load_sidebar_script(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            page = args[0] / 'namespace_math.html'
            page.write_text(page.read_text().replace('<script src="rodin-navigation.js" defer></script>', ''))
            with self.assertRaisesRegex(ValueError, 'Page does not load sidebar navigation'):
                verify(*args)

    def test_documented_headers_must_appear_in_file_index(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            (args[0] / 'files.html').write_text('An incomplete file list')
            with self.assertRaisesRegex(ValueError, 'Documented declaration missing from reference index: Header.h'):
                verify(*args)

    def test_rodin_types_must_stay_inside_the_namespace_tree(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            (args[0] / 'annotated.html').write_text('<main><ul class="m-doc"><li>'
                '<a href="class_foo.html">Rodin::Foo</a></li></ul></main>')
            with self.assertRaisesRegex(ValueError, 'Type escaped its namespace hierarchy'):
                verify(*args)

    def test_external_types_must_stay_inside_their_namespace_tree(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            index = args[1] / 'index.xml'
            index.write_text(index.read_text().replace('Rodin::Foo', 'Eigen::internal::Foo'))
            (args[0] / 'annotated.html').write_text('<main><ul class="m-doc"><li>'
                '<a href="class_foo.html">Eigen::internal::Foo</a></li></ul></main>')
            with self.assertRaisesRegex(ValueError, 'Type escaped its namespace hierarchy'):
                verify(*args)

    def test_documented_types_must_appear_in_class_index(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            (args[0] / 'annotated.html').write_text('An incomplete class list')
            with self.assertRaisesRegex(ValueError, 'Documented declaration missing from reference index: Rodin::Foo'):
                verify(*args)

    def test_complete_public_navigation(self):
        with tempfile.TemporaryDirectory() as folder:
            verify(*self.fixture(Path(folder)))

    def test_overview_order_matches_brief_representative_description(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            p = args[0] / 'modules.html'
            p.write_text(p.read_text().replace('data-role="representative"', 'data-role="brief"', 1))
            with self.assertRaisesRegex(ValueError, 'overview order'):
                verify(*args)

    def test_overview_rejects_repeated_svg_ids(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            p = args[0] / 'modules.html'
            p.write_text(p.read_text() + '<svg><path id="g"/><path id="g"/></svg>')
            with self.assertRaisesRegex(ValueError, 'Duplicate HTML or SVG IDs'):
                verify(*args)

    def test_old_technical_group_index_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            with (args[0] / 'modules.html').open('a') as page:
                page.write('<a href="group___rodin_directives.html">Preprocessor directives</a>')
            with self.assertRaisesRegex(ValueError, 'Technical.*groups'):
                verify(*args)

    def test_specializations_must_be_linked_from_the_generic_class(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            (args[0] / 'class_foo.html').write_text('A class with no specialization reference')
            with self.assertRaisesRegex(ValueError, 'Specialization reference link'):
                verify(*args)

    def test_identically_named_traits_and_operators_use_the_owning_namespace(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            index = args[1] / 'index.xml'
            index.write_text(index.read_text().replace('</doxygenindex>',
                '<compound kind="struct" refid="trait_foo"><name>Rodin::FormLanguage::Foo</name></compound></doxygenindex>'))
            (args[1] / 'trait_foo.xml').write_text('<doxygen><compounddef><templateparamlist><param>T</param>'
                '</templateparamlist></compounddef></doxygen>')
            group = args[1] / 'group_foo.xml'
            group.write_text(group.read_text().replace('<title>', '<innerclass>Rodin::Foo&lt; T &gt;</innerclass><title>'))
            verify(*args)

    def test_plain_text_see_also_is_not_a_specialization_link(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            p = args[1] / 'class_foo.xml'
            p.write_text(p.read_text().replace('<ref refid="group_foo">Foo Template Specializations</ref>',
                                              'FooSpecializations'))
            with self.assertRaisesRegex(ValueError, 'Primary template has no See also link'):
                verify(*args)

    def test_configuration_must_link_to_the_rodin_header_reference(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            (args[0] / 'configuration.html').write_text('<a href="_m_m_g_2_configure_8h_8in.html">Configure.h.in</a>')
            with self.assertRaisesRegex(ValueError, 'Configuration does not link'):
                verify(*args)

    def test_modules_need_a_representative_formula_or_usage(self):
        with tempfile.TemporaryDirectory() as folder:
            args = self.fixture(Path(folder))
            p = args[1] / 'namespace_math.xml'
            p.write_text(p.read_text().replace('<formula>x+y</formula>', ''))
            with self.assertRaisesRegex(ValueError, 'representative formula or usage'):
                verify(*args)


if __name__ == '__main__':
    unittest.main()

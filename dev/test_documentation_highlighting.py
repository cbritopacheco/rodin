"""Rodin identifiers are highlighted without changing C++ snippet contents."""
from pathlib import Path
import tempfile
import unittest

from documentation_highlighting import cpp_highlighter, cpp_symbols


class DocumentationHighlightingTest(unittest.TestCase):
    highlight = staticmethod(cpp_highlighter())

    def test_class_and_enum_tokens(self):
        source = '<span class="n">Mesh</span> <span class="n">P1</span> <span class="n">Triangle</span>'
        self.assertEqual(self.highlight(source), '<span class="nc">Mesh</span> <span class="nc">P1</span> <span class="no">Triangle</span>')

    def test_variables_and_literal_text_are_preserved(self):
        source = '<span class="n">mesh</span> <span class="nf">vh</span> <span class="s">"Triangle"</span> <span class="c1">// Mesh</span>'
        self.assertEqual(self.highlight(source), source)

    def test_punctuation_and_existing_highlighting_are_preserved(self):
        source = '<span class="nc">Mesh</span><span class="p">{</span><span class="mi">16</span><span class="p">}</span>'
        self.assertEqual(self.highlight(source), source)

    def test_all_declared_enum_values_are_classified(self):
        highlight = cpp_highlighter(enum_values={'Boundary': {'Mode'}, 'Interior': {'Mode'}})
        self.assertEqual(highlight('<span class="n">Boundary</span> <span class="n">Interior</span>'),
                         '<span class="no">Boundary</span> <span class="no">Interior</span>')

    def test_qualified_enum_value_does_not_recolor_a_class_with_the_same_name(self):
        highlight = cpp_highlighter(type_names={'Point'}, enum_values={'Point': {'Rodin::Polytope::Type'}})
        source = '<span class="n">Point</span> <span class="n">Type</span><span class="o">::</span><span class="n">Point</span>'
        result = highlight(source)
        self.assertTrue(result.startswith('<span class="nc">Point</span>'))
        self.assertTrue(result.endswith('<span class="no">Point</span>'))

    def test_xml_supplies_types_and_scoped_enum_values(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'type.xml').write_text('<doxygen><compounddef kind="class"><compoundname>Rodin::Field&lt; Scalar &gt;</compoundname>'
                '<sectiondef><memberdef kind="enum" strong="yes"><name>Mode</name><qualifiedname>Rodin::Field::Mode</qualifiedname>'
                '<enumvalue><name>Boundary</name></enumvalue></memberdef></sectiondef></compounddef></doxygen>')
            types, enums = cpp_symbols(root)
            self.assertIn('Field', types)
            self.assertIn('Mode', types)
            self.assertEqual(enums['Boundary'], {'Rodin::Field::Mode'})


if __name__ == '__main__':
    unittest.main()

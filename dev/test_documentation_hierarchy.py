"""Regressions for missing Doxygen ownership and hidden directory ancestors."""
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

from documentation_hierarchy import prepare_hierarchy, scope_prefixes, symbol_index


class DocumentationHierarchyTest(unittest.TestCase):
    def fixture(self, root, entries):
        index = ET.Element('doxygenindex')
        for identity, kind, name, children in entries:
            entry = ET.SubElement(index, 'compound', refid=identity, kind=kind)
            ET.SubElement(entry, 'name').text = name
            tree = ET.Element('doxygen')
            definition = ET.SubElement(tree, 'compounddef', id=identity, kind=kind)
            ET.SubElement(definition, 'compoundname').text = name
            ET.SubElement(definition, 'briefdescription')
            ET.SubElement(definition, 'detaileddescription')
            for tag, child in children:
                ET.SubElement(definition, tag, refid=child)
            ET.ElementTree(tree).write(root / (identity + '.xml'))
        ET.ElementTree(index).write(root / 'index.xml')

    def test_scope_resolution_skips_template_argument_namespaces(self):
        self.assertEqual(scope_prefixes('Rodin::Variational::GridFunction<FES,::Vec>::ArrayRead'),
                         ['Rodin::Variational::GridFunction<FES,::Vec>', 'Rodin::Variational', 'Rodin'])

    def test_nested_specialization_retains_its_actual_owner(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, [('namespace', 'namespace', 'Rodin::Variational', []),
                ('owner', 'class', 'Rodin::Variational::GridFunction< FES, ::Vec >', []),
                ('child', 'struct', 'Rodin::Variational::GridFunction< FES, ::Vec >::ArrayRead', [])])
            self.assertEqual(prepare_hierarchy(root)['ownership'], 2)
            owner = ET.parse(root / 'owner.xml').find('compounddef')
            self.assertEqual(owner.find('innerclass').get('refid'), 'child')
            namespace = ET.parse(root / 'namespace.xml').find('compounddef')
            self.assertEqual([c.get('refid') for c in namespace.findall('innerclass')], ['owner'])
            self.assertEqual(prepare_hierarchy(root), {'directories': 0, 'ownership': 0})

    def test_existing_ownership_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, [('owner', 'class', 'Rodin::Owner', [('innerclass', 'child')]),
                ('child', 'struct', 'Rodin::Owner::Child', [])])
            self.assertEqual(prepare_hierarchy(root)['ownership'], 0)
            self.assertEqual(len(ET.parse(root / 'owner.xml').findall('./compounddef/innerclass')), 1)

    def test_missing_specialization_falls_back_to_enclosing_namespace(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, [('namespace', 'namespace', 'Rodin::Math', []),
                ('child', 'class', 'Rodin::Math::Vector< Scalar >::Nested', [])])
            self.assertEqual(prepare_hierarchy(root)['ownership'], 1)
            self.assertEqual(ET.parse(root / 'namespace.xml').find('./compounddef/innerclass').get('refid'), 'child')

    def test_documented_child_survives_undocumented_enclosing_specialization(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            definitions = [('namespace', 'namespace', 'Rodin::Geometry', []),
                ('owner', 'class', 'Rodin::Geometry::Mesh<Context::Local>', []),
                ('child', 'class', 'Rodin::Geometry::Mesh<Context::Local>::Builder', [])]
            self.fixture(root, definitions)
            compounds = {}
            for identity, kind, name, _ in definitions:
                compounds[identity] = SimpleNamespace(kind=kind, name=name,
                    has_details=identity != 'owner', parent={'namespace': None, 'owner': 'namespace', 'child': 'owner'}[identity],
                    url=identity + '.html', brief='', deprecated=None, since=None)
            result = symbol_index(compounds, root / 'index.xml')
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0].name, 'Rodin::Geometry')
            self.assertEqual(result[0].children[0].name, 'Mesh<Context::Local>::Builder')
            self.assertEqual(result[0].children[0].url, 'child.html')

    def test_anonymous_namespace_alias_does_not_capture_public_types(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, [('namespace', 'namespace', 'Rodin::Geometry', []),
                ('alias', 'namespace', 'Rodin::Geometry::@123', []),
                ('child', 'class', 'Rodin::Geometry::Mesh<Local>::Builder', [])])
            tree = ET.parse(root / 'alias.xml')
            tree.find('compounddef/compoundname').text = 'Rodin::Geometry'
            tree.write(root / 'alias.xml')
            prepare_hierarchy(root)
            self.assertEqual(ET.parse(root / 'namespace.xml').find('./compounddef/innerclass').get('refid'), 'child')
            self.assertIsNone(ET.parse(root / 'alias.xml').find('./compounddef/innerclass'))

    def test_library_namespace_precedes_external_customizations(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            names = ['Eigen', 'Rodin', 'std']
            self.fixture(root, [(name, 'namespace', name, []) for name in names])
            compounds = {name: SimpleNamespace(kind='namespace', name=name,
                has_details=True, parent=None, url=name + '.html', brief='',
                deprecated=None, since=None) for name in names}
            result = symbol_index(compounds, root / 'index.xml')
            self.assertEqual([entry.name for entry in result], ['Rodin', 'Eigen', 'std'])

    def test_external_namespace_trees_start_collapsed_and_remain_toggleable(self):
        from types import SimpleNamespace
        from jinja2 import ChoiceLoader, DictLoader, Environment, FileSystemLoader
        from html.parser import HTMLParser

        class Trees(HTMLParser):
            def __init__(self):
                super().__init__()
                self.classes = []

            def handle_starttag(self, tag, attrs):
                if tag == 'li':
                    self.classes.append(dict(attrs).get('class', ''))

        environment = Environment(loader=ChoiceLoader([
            DictLoader({'base-index.html': '{% block main %}{% endblock %}'}),
            FileSystemLoader(Path(__file__).resolve().parents[1] / 'doc/templates')]))
        environment.filters['rtrim'] = str.rstrip
        child = SimpleNamespace(kind='class', name='Type', url='type.html', brief='',
            children=[], is_inline=False, is_final=False, deprecated=None, since=None)
        roots = [SimpleNamespace(kind='namespace', name=name, url=name + '.html',
            brief='', children=[child], is_inline=False, is_final=False,
            deprecated=None, since=None) for name in ['Rodin', 'Eigen', 'std']]
        html = environment.get_template('annotated.html').render(
            index=SimpleNamespace(symbols=roots), CLASS_INDEX_EXPAND_LEVELS=10,
            CLASS_INDEX_EXPAND_INNER=True)
        parser = Trees()
        parser.feed(html)
        self.assertEqual(parser.classes[::2], ['m-doc-collapsible',
            'm-doc-collapsible collapsed', 'm-doc-collapsible collapsed'])
        self.assertEqual(html.count('onclick="return toggle(this)"'), 3)

    def test_private_types_do_not_create_links_to_unrendered_pages(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, [('owner', 'class', 'Hidden', [])])
            tree = ET.parse(root / 'owner.xml')
            tree.find('compounddef').set('prot', 'private')
            tree.write(root / 'owner.xml')
            compound = SimpleNamespace(kind='class', name='Hidden', has_details=True,
                parent=None, url='owner.html', brief='', deprecated=None, since=None)
            self.assertEqual(symbol_index({'owner': compound}, root / 'index.xml'), [])

    def test_documentation_only_directories_are_not_added_to_file_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, [('directory', 'dir', '/temporary/checkout/doc', [('innerfile', 'config')]),
                ('config', 'file', 'conf.py.in', [])])
            self.assertEqual(prepare_hierarchy(root)['directories'], 0)
            self.assertFalse(len(ET.parse(root / 'directory.xml').find('./compounddef/briefdescription')))

    def test_header_directory_is_visible_without_source_boilerplate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root, [('directory', 'dir', 'Rodin/Math', [('innerfile', 'header')]),
                ('header', 'file', 'Header.h', [])])
            tree = ET.parse(root / 'header.xml')
            ET.SubElement(tree.find('compounddef/briefdescription'), 'para').text = 'Header documentation.'
            tree.write(root / 'header.xml')
            self.assertEqual(prepare_hierarchy(root)['directories'], 1)
            description = ET.parse(root / 'directory.xml').find('./compounddef/briefdescription')
            self.assertEqual(''.join(description.itertext()), 'Headers in Rodin/Math.')
            self.assertEqual(prepare_hierarchy(root)['directories'], 0)


if __name__ == '__main__':
    unittest.main()

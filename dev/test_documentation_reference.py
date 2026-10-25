"""Reference presentation preserves C++ contracts, links and member templates."""
from pathlib import Path
import tempfile
import unittest

from documentation_reference import (ReferenceLayout, Markup, parameters,
                                     prepare_reference_headers, render_references, scope_name,
                                     split_arguments, text)
import xml.etree.ElementTree as ET


class ReferenceTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.xml, self.site = self.root / 'xml', self.root / 'html'
        self.xml.mkdir()
        self.site.mkdir()
        self.write('primary', 'Rodin::IO::Loader', '<param><type>FileFormat</type><declname>Fmt</declname></param><param><type>class FES</type></param><param><type>class Data</type></param>')
        self.write('special', 'Rodin::IO::Loader&lt;FileFormat::MEDIT, P1&lt;Range, Mesh&lt;Local&gt;&gt;, Vector&lt;Range&gt;&gt;', '<param><type>class Range</type></param>')
        self.write('other', 'Rodin::Math::Loader')
        self.layout = ReferenceLayout(self.xml, self.site)

    def write(self, identity, name, params='', kind='class'):
        source = ('<doxygen><compounddef id="' + identity + '" kind="' + kind + '"><compoundname>' + name + '</compoundname>' +
                  ('<templateparamlist>' + params + '</templateparamlist>' if params else '') +
                  '<briefdescription><para>Original brief.</para></briefdescription></compounddef></doxygen>')
        (self.xml / (identity + '.xml')).write_text(source)
        (self.site / (identity + '.html')).touch()

    def test_nested_arguments_and_scope(self):
        self.assertEqual(split_arguments('A<B,C>, D<(1,2)>, "a,b"'), ['A<B,C>', 'D<(1,2)>', '"a,b"'])
        self.assertEqual(scope_name('Rodin::Outer<A::B>::Nested<C::D>'), ('Rodin::Outer<A::B>', 'Nested<C::D>'))

    def test_specialization_bindings_do_not_invent_instantiations(self):
        result = self.layout.bindings('special')
        self.assertIn('<code>Fmt</code>', result)
        self.assertIn('<code>FES</code>', result)
        self.assertIn('P1&lt;Range, Mesh&lt;Local&gt;&gt;', result)
        self.assertNotIn('Real', result)
        self.assertEqual(self.layout.primaries['special'], 'primary')
        self.assertNotIn('other', self.layout.families)

    def test_dependent_member_is_not_linked_to_a_same_named_local_alias(self):
        definition = ET.fromstring('<compounddef><ref kindref="member" refid="primary_1a' + '0' * 32 + '">ScalarType</ref></compounddef>')
        result = self.layout.code('typename Traits<T>::ScalarType', definition)
        self.assertNotIn('<a ', result)

    def test_parameter_pack_and_unnamed_parameters(self):
        self.write('pack', 'Rodin::Pack', '<param><type>class... Args</type></param>')
        self.write('packed', 'Rodin::Pack&lt;A, B&lt;C,D&gt;&gt;')
        layout = ReferenceLayout(self.xml, self.site)
        self.assertIn('Args...', layout.bindings('packed'))
        self.assertIn('A, B&lt;C,D&gt;', layout.bindings('packed'))
        definition = ET.fromstring('<compounddef><templateparamlist><param><type>int</type></param></templateparamlist></compounddef>')
        self.assertEqual(parameters(definition)[0][0], '_1')

    def test_short_title_retains_brief_and_full_declaration(self):
        source = '<main><h1>Rodin::IO::Loader&lt;...&gt;</h1><p>Keep <a href="other.html">this link</a>.</p><p>Details.</p></main>'
        result = self.layout.render('special', source)
        self.assertIn('<h1>Loader</h1>', result)
        self.assertIn('href="other.html"', result)
        self.assertIn('<p>Details.</p>', result)
        self.assertIn('Specialization arguments', result)
        self.assertNotIn('Direct inheritance', result)
        plain = Markup(result)
        self.assertIn('class Loader<', plain.content(next(n for n in plain.nodes if n.tag == 'pre')))

    def test_forward_declaration_without_brief_gets_family_table(self):
        source = '<main><h1>Loader</h1><table><tr><td>Details.</td></tr></table></main>'
        result = self.layout.render('primary', source)
        self.assertIn('Full specialization', result)
        self.assertIn('href="special.html"', result)

    def test_declaration_preserves_constraints_and_defaults(self):
        definition = ET.fromstring('<compounddef kind="class"><compoundname>Rodin::Constrained</compoundname><templateparamlist><param><type>class T</type><defval>Real</defval></param></templateparamlist><requiresclause>Concept&lt;T&gt;</requiresclause></compounddef>')
        declaration = self.layout.declaration(definition)
        self.assertIn('class T = Real', declaration)
        self.assertIn('requires Concept<T>', declaration)

    def test_member_prefix_removal_preserves_return_and_method_template(self):
        source = '''<section class="m-doc-details" id="member"><div><h3>
<div class="m-doc-template">template&lt;class Range&gt; template&lt;class U&gt;</div>
<span class="m-doc-wrap-bumper">const <a href="other.html">Math::Vector</a>&amp; Rodin::<wbr />IO::<wbr />Loader&lt;FileFormat::MEDIT, P1&lt;Range, Mesh&lt;Local&gt;&gt;, Vector&lt;Range&gt;&gt;&lt;Range&gt;::<wbr /></span><span class="m-doc-wrap"><span class="m-doc-wrap-bumper"><a href="#member" class="m-doc-self">load</a>(</span>U value) const noexcept</span></h3><p>Original contract.</p></div></section>'''
        result = self.layout.render('special', source)
        self.assertNotIn('Rodin::', result)
        self.assertIn('href="other.html">Math::Vector</a>&amp;', result)
        self.assertIn('template&lt;class U&gt;', result)
        self.assertIn('const noexcept', result)
        self.assertIn('id="member"', result)
        self.assertIn('Original contract.', result)

    def test_namespace_context_keeps_external_qualification(self):
        self.write('namespace', 'Rodin::IO', kind='namespace')
        layout = ReferenceLayout(self.xml, self.site)
        source = '<section class="m-doc-details"><div><h3><span class="m-doc-wrap-bumper">Geometry::Point Rodin::<wbr />IO::<wbr /></span><span class="m-doc-wrap"><a href="#member">read</a>()</span></h3></div></section>'
        result = layout.render('namespace', source)
        self.assertNotIn('Rodin::', result)
        self.assertIn('Geometry::Point', result)

    def test_class_list_groups_only_present_members_preserving_links(self):
        source = '<dl class="m-doc"><dt>class <a href="primary.html">Loader</a></dt><dd>Primary contract.</dd><dt>class <a href="special.html">Loader&lt;...&gt;</a></dt><dd>Special brief.</dd><dt><a href="other.html">Other</a></dt><dd>Other contract.</dd></dl>'
        result = self.layout.group_entries(source)
        self.assertIn('Primary contract.', result)
        self.assertIn('Other contract.', result)
        self.assertEqual(result.count('href="special.html"'), 1)
        self.assertIn('Full specialization', result)
        self.assertIn('href="other.html"', result)
        markup = Markup(result)
        disclosure = next(node for node in markup.nodes if node.attrs.get('class') == 'rodin-specialization-disclosure')
        self.assertNotIn('open', disclosure.attrs)
        self.assertIn('<summary>See specializations</summary>', markup.raw(disclosure))
        self.assertIn('href="special.html"', markup.raw(disclosure))
        self.assertNotIn('Primary contract.', markup.raw(disclosure))

    def test_missing_primary_groups_by_qualified_family_without_invention(self):
        self.write('a', 'Rodin::IO::Missing&lt;A&lt;B,C&gt;, D&gt;')
        self.write('b', 'Rodin::IO::Missing&lt;E&gt;')
        self.write('c', 'Rodin::Math::Missing&lt;A&gt;')
        layout = ReferenceLayout(self.xml, self.site)
        source = '<dl class="m-doc">' + ''.join('<dt>class <a href="' + key + '.html">Missing&lt;...&gt;</a></dt><dd>Brief.</dd>' for key in ['a', 'b', 'c']) + '</dl>'
        result = layout.group_entries(source)
        self.assertEqual(result.count('class="rodin-specialization-disclosure"'), 2)
        self.assertEqual(result.count('Missing specializations</dt>'), 2)
        self.assertNotIn('template&lt;', result.split('Full specialization')[0])
        for key in ['a', 'b', 'c']:
            self.assertEqual(result.count('href="' + key + '.html"'), 1)

    def test_long_overview_declarations_preserve_links_briefs_and_anchors(self):
        parameters = ', '.join('class Parameter' + str(i) for i in range(12))
        raw = '<div class="m-doc-template">template&lt;' + parameters + '&gt;</div>class <a href="primary.html">Loader</a>'
        source = '<dl class="m-doc"><dt id="class-anchor">' + raw + '</dt><dd>Primary brief.</dd></dl>'
        result = self.layout.collapse_declarations(source)
        self.assertIn('<dt id="class-anchor">class <a href="primary.html">Loader</a></dt>', result)
        self.assertIn('Primary brief.', result)
        self.assertIn('<summary>See declaration</summary>', result)
        self.assertIn(raw, result)
        self.assertNotIn(' open', result)
        function = '<span class="m-doc-wrap-bumper">auto <a href="#function">assemble</a>(</span><span class="m-doc-wrap">' + parameters + ') -&gt; <a href="other.html">Result</a></span>'
        result = self.layout.collapse_declarations('<dl class="m-doc"><dt>' + function + '</dt><dd>Function brief.</dd></dl>')
        self.assertIn('<dt><a href="#function">assemble</a></dt>', result)
        self.assertIn(function, result)
        self.assertEqual(result.count('href="other.html"'), 1)
        short = '<dl class="m-doc"><dt>class <a href="primary.html">Loader</a></dt><dd>Short.</dd></dl>'
        self.assertEqual(self.layout.collapse_declarations(short), short)

    def test_deduction_guide_details_preserve_declaration_and_contract(self):
        self.write('space', 'Rodin::IO', kind='namespace')
        path = self.xml / 'space.xml'
        path.write_text(path.read_text().replace('</compounddef>', '<sectiondef><memberdef kind="function" id="space_1guide"><type/><name>Loader</name><argsstring>(const T&amp; value) -&gt; Loader&lt;T&gt;</argsstring></memberdef><memberdef kind="function" id="space_1ordinary"><type>auto</type><name>Loader</name><argsstring>(const T&amp; value) -&gt; Loader&lt;T&gt;</argsstring></memberdef></sectiondef></compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        signature = '<div class="m-doc-template">template&lt;class T&gt;</div><a href="#guide" class="m-doc-self">Loader</a>(const <a href="other.html">T</a>&amp; value) -&gt; Loader&lt;T&gt;'
        header = '<p class="rodin-reference-include">Defined in <a href="header.html">Loader.h</a></p>'
        contract = '<p>Deduction guide.</p><table><tr><td>value</td><td>Input.</td></tr></table>'
        source = '<section class="m-doc-details" id="guide"><div><h3>' + signature + '</h3>' + header + contract + '</div></section>'
        result = layout.deduction_guides('space', source)
        self.assertIn('<h3><a href="#guide" class="m-doc-self">Loader</a></h3>' + header + '<details', result)
        self.assertIn('<summary>See declaration</summary>', result)
        self.assertIn(signature, result)
        self.assertIn('</details>' + contract, result)
        self.assertNotIn(' open', result)
        self.assertEqual(layout.deduction_guides('space', result), result)
        ordinary = source.replace('id="guide"', 'id="ordinary"')
        self.assertEqual(layout.deduction_guides('space', ordinary), ordinary)

    def test_source_headers_are_normalized_on_group_and_file_pages(self):
        xml = self.root / 'empty-xml'
        xml.mkdir()
        annotation = '<h3><div class="m-doc-include m-code"><span>#include</span> <a href="_polytope_8h.html">&lt;Rodin/Geometry/Polytope.h&gt;</a></div>Type</h3>'
        example = '<pre class="m-code"><span>#include</span> &lt;Rodin/Header.h&gt;</pre>'
        for name in ['group_types.html', '_header_8h.html', 'guide.html']:
            (self.site / name).write_text(annotation + example)
        self.assertEqual(render_references(xml, self.site), 3)
        for name in ['group_types.html', '_header_8h.html', 'guide.html']:
            result = (self.site / name).read_text()
            self.assertNotIn('m-doc-include', result)
            self.assertIn('Defined in <a href="_polytope_8h.html"><code>Rodin/Geometry/Polytope.h</code></a>', result)
            self.assertIn(example, result)
        self.assertEqual(render_references(xml, self.site), 0)

    def test_source_header_uses_defined_in_outside_signature(self):
        include = '<div class="m-doc-include m-code m-text-right"><span>#include</span> <a class="cpf" href="header.html">&lt;Rodin/Assembly/OpenMP.h&gt;</a></div>'
        example = '<pre><span>#include &lt;Example.h&gt;</span></pre>'
        source = '<h3><a href="#function">assemble</a>()' + include + '</h3>' + example
        result = self.layout.source_headers(source)
        self.assertIn('</h3><p class="rodin-reference-include">Defined in <a href="header.html"><code>Rodin/Assembly/OpenMP.h</code></a></p>', result)
        self.assertNotIn('m-doc-include', result)
        self.assertIn(example, result)
        self.assertEqual(self.layout.source_headers(result), result)
        unlinked = self.layout.source_headers('<div class="m-doc-include">#include &lt;Local.h&gt;</div>')
        self.assertIn('Defined in <code>Local.h</code>', unlinked)

    def test_typedef_headings_drop_only_the_current_scope(self):
        self.write('space', 'Rodin::IO', kind='namespace')
        layout = ReferenceLayout(self.xml, self.site)
        actual = '<section class="m-doc-details" id="alias"><div><h3>using Rodin::' + '<wbr />IO::<wbr /><a href="#alias" class="m-doc-self">MeshType</a> = <a href="external.html">Geometry::Mesh</a></h3></div></section>'
        result = layout.render('space', actual)
        self.assertNotIn('Rodin::', result)
        self.assertIn('Geometry::Mesh', result)
        self.assertIn('id="alias"', result)
        self.assertIn('href="#alias"', result)
        self.assertIn('href="external.html"', result)
        literal = actual.replace('Geometry::Mesh</a>', 'Geometry::Mesh</a>&lt;&quot;Rodin::IO::&quot;&gt;')
        self.assertIn('Rodin::IO::&quot;', layout.render('space', literal))

    def test_inheritance_keeps_full_types_and_rich_briefs(self):
        self.write('base', 'Rodin::IO::Base')
        self.write('worker', 'Rodin::IO::Worker')
        path = self.xml / 'worker.xml'
        path.write_text(path.read_text().replace('</compounddef>', '<basecompoundref refid="base" prot="public" virt="virtual">Rodin::IO::Base&lt;Math::Vector&lt;Real&gt;&gt;</basecompoundref></compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        source = '<main><h1>Worker</h1><section id="base-classes"><h2>Base classes</h2><dl class="m-doc"><dt>class <a href="base.html">Base</a></dt><dd>Contract with <a href="other.html">details</a>.</dd></dl></section></main>'
        result = layout.render('worker', source)
        self.assertIn('Direct inheritance', result)
        self.assertIn('Inherits from', result)
        self.assertIn('public, virtual', result)
        self.assertIn('href="other.html">details</a>', result)
        self.assertIn('Math::Vector', result)
        self.assertIn('id="base-classes"', result)
        self.assertIn('Full type', result)

    def test_shortened_cross_module_relation_keeps_module_qualification(self):
        self.write('mesh', 'Rodin::Geometry::Mesh')
        self.write('mmgmesh', 'Rodin::MMG::Mesh')
        self.write('submesh', 'Rodin::Geometry::SubMesh')
        path = self.xml / 'mesh.xml'
        path.write_text(path.read_text().replace('</compounddef>',
            '<derivedcompoundref refid="mmgmesh" prot="public">Mesh</derivedcompoundref>'
            '<derivedcompoundref refid="submesh" prot="public">SubMesh&lt;Context::Local&gt;</derivedcompoundref></compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        result = layout.render('mesh', '<main><h1>Mesh</h1><section id="derived-classes"><h2>Derived classes</h2></section></main>')
        self.assertIn('<code>MMG::Mesh</code>', result)
        self.assertIn('title="MMG::Mesh">MMG::Mesh</a>', result)
        self.assertIn('<code>SubMesh&lt;Context::Local&gt;</code>', result)
        self.assertNotIn('<code>Geometry::SubMesh', result)

    def test_diagram_is_bounded_but_every_child_is_listed(self):
        self.write('worker', 'Rodin::IO::Worker')
        relations = ''
        for index in range(5):
            self.write('child' + str(index), 'Rodin::IO::Child' + str(index))
            relations += '<derivedcompoundref refid="child' + str(index) + '" prot="public">Child' + str(index) + '</derivedcompoundref>'
        path = self.xml / 'worker.xml'
        path.write_text(path.read_text().replace('</compounddef>', relations + '</compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        source = '<main><h1>Worker</h1><section id="derived-classes"><h2>Derived classes</h2></section></main>'
        result = layout.render('worker', source)
        self.assertIn('+2 more', result)
        for index in range(5):
            self.assertIn('href="child' + str(index) + '.html"', result)
        diagram = result[result.index('<figure'):result.index('</figure>')]
        self.assertEqual(diagram.count('title="Child'), 3)
        self.assertNotIn('child4.html', diagram)

    def test_distinct_families_precede_repeated_specializations(self):
        self.write('worker', 'Rodin::IO::Worker')
        self.write('third', 'Rodin::IO::Third')
        relations = ''.join('<derivedcompoundref refid="special" prot="public">Loader&lt;' + value + '&gt;</derivedcompoundref>' for value in ['A', 'B', 'C', 'D'])
        relations += '<derivedcompoundref refid="other" prot="public">Rodin::Math::Loader</derivedcompoundref><derivedcompoundref refid="third" prot="public">Third</derivedcompoundref>'
        path = self.xml / 'worker.xml'
        path.write_text(path.read_text().replace('</compounddef>', relations + '</compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        result = layout.render('worker', '<main><h1>Worker</h1><section id="derived-classes"><h2>Derived classes</h2></section></main>')
        page = Markup(result)
        diagram = next(n for n in page.nodes if n.tag == 'figure')
        self.assertEqual(page.raw(diagram).count('<li>'), 3)
        self.assertIn('href="primary.html"', page.raw(diagram))
        self.assertIn('href="other.html"', page.raw(diagram))
        self.assertIn('href="third.html"', page.raw(diagram))
        table = next(n for n in page.nodes if n.tag == 'table')
        tbody = next(n for n in table.children if n.tag == 'tbody')
        self.assertEqual(len(tbody.children), 3)
        self.assertIn('href="other.html"', page.raw(table))
        self.assertIn('href="third.html"', page.raw(table))
        disclosure = next(n for n in page.nodes if n.tag == 'details' and n.attrs.get('class') == 'rodin-derived-disclosure')
        self.assertNotIn('open', disclosure.attrs)
        self.assertIn('See more', page.raw(disclosure))
        self.assertIn('See less', page.raw(disclosure))
        self.assertIn('Loader&lt;B&gt;', page.raw(disclosure))
        self.assertIn('Loader&lt;D&gt;', page.raw(disclosure))

    def test_more_links_count_grouped_entries_for_both_directions(self):
        relations = [ET.fromstring('<basecompoundref refid="special">Loader&lt;' + name + '&gt;</basecompoundref>') for name in ['A', 'B', 'C', 'D']]
        diagram = self.layout.inheritance_diagram('primary', [(relations, 'base-classes'), (relations, 'derived-classes')])
        self.assertIn('href="#base-classes">+3 more</a>', diagram)
        self.assertIn('href="#derived-classes">+3 more</a>', diagram)
        small = self.layout.inheritance_diagram('primary', [(relations[:3], 'base-classes'), ([], 'derived-classes')])
        self.assertNotIn(' more</a>', small)

    def test_overrides_and_inherited_members_link_to_declaring_base(self):
        self.write('base', 'Rodin::IO::Base')
        self.write('worker', 'Rodin::IO::Worker')
        base_member = 'base_1a' + '0' * 32
        own_member = 'worker_1a' + '1' * 32
        path = self.xml / 'base.xml'
        path.write_text(path.read_text().replace('</compounddef>', '<sectiondef><memberdef id="' + base_member + '" kind="function" prot="public"><name>read</name></memberdef></sectiondef></compounddef>'))
        path = self.xml / 'worker.xml'
        extra = '<basecompoundref refid="base" prot="public">Base</basecompoundref>'
        extra += '<sectiondef><memberdef id="' + own_member + '" kind="function" prot="public"><name>read</name><reimplements refid="' + base_member + '">read</reimplements></memberdef></sectiondef>'
        extra += '<listofallmembers><member refid="' + base_member + '" prot="public"><name>read</name></member><member refid="' + base_member + '" prot="private"><name>privateLeak</name></member><member refid="' + base_member + '" prot="public"><name>Base</name></member></listofallmembers>'
        path.write_text(path.read_text().replace('</compounddef>', extra + '</compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        source = '<main><h1>Worker</h1><section id="base-classes"><h2>Base classes</h2></section></main>'
        result = layout.render('worker', source)
        self.assertIn('Overrides', result)
        self.assertIn('Inherited interface', result)
        self.assertNotIn('privateLeak', result)
        self.assertNotIn('<code>Base()</code>', result)
        self.assertIn('href="base.html#a' + '0' * 32 + '"', result)
        self.assertIn('href="worker.html#a' + '1' * 32 + '"', result)

    def test_see_also_omits_hierarchy_and_specialization_duplicates(self):
        self.write('indirect', 'Rodin::IO::MeshLoader')
        path = self.xml / 'primary.xml'
        path.write_text(path.read_text().replace('</compounddef>', '<derivedcompoundref refid="indirect" prot="public">MeshLoader</derivedcompoundref></compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        source = '<aside><h4>See also</h4><p><a href="indirect.html">MeshLoader</a></p><p><a href="special.html">Loader specialization</a></p><p><a href="primary.html">Loader</a></p></aside>'
        self.assertEqual(layout.remove_related_see_also('primary', source), '')

    def test_cleanup_recovers_missing_reverse_inheritance_edges(self):
        self.write('indirect', 'Rodin::IO::MeshLoader')
        path = self.xml / 'indirect.xml'
        path.write_text(path.read_text().replace('</compounddef>', '<basecompoundref refid="primary" prot="public">Loader</basecompoundref></compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        source = '<aside><h4>See also</h4><p><a href="indirect.html">MeshLoader</a></p></aside>'
        self.assertEqual(layout.remove_related_see_also('primary', source), '')

    def test_see_also_keeps_unrelated_links_prose_and_member_contracts(self):
        source = '<aside><h4>See also</h4><p><a href="special.html">Specialization</a></p><p><a href="other.html">Other module</a></p><p><a href="group_types.html">Specialization reference</a></p><p>Read <a href="primary.html">Loader</a> for the contract.</p></aside>'
        result = self.layout.remove_related_see_also('primary', source)
        self.assertNotIn('href="special.html"', result)
        self.assertIn('href="other.html"', result)
        self.assertIn('href="group_types.html"', result)
        self.assertIn('for the contract.', result)
        mixed = '<aside><h4>See also</h4><p><a href="primary.html">Loader</a>, <a href="primary.html#member">load</a></p></aside>'
        self.assertEqual(self.layout.remove_related_see_also('primary', mixed), mixed)
        member = '<section class="m-doc-details"><div><aside><h4>See also</h4><p><a href="primary.html">Loader</a></p></aside></div></section>'
        self.assertEqual(self.layout.remove_related_see_also('primary', member), member)

    def test_hierarchy_cleanup_does_not_remove_sibling_references(self):
        self.write('base', 'Rodin::Base')
        for identity in ['primary', 'other']:
            path = self.xml / (identity + '.xml')
            path.write_text(path.read_text().replace('</compounddef>', '<basecompoundref refid="base" prot="public">Base</basecompoundref></compounddef>'))
        path = self.xml / 'base.xml'
        path.write_text(path.read_text().replace('</compounddef>', '<derivedcompoundref refid="primary">Loader</derivedcompoundref><derivedcompoundref refid="other">Other</derivedcompoundref></compounddef>'))
        layout = ReferenceLayout(self.xml, self.site)
        source = '<aside><h4>See also</h4><p><a href="base.html">Base</a></p><p><a href="other.html">Sibling</a></p></aside>'
        result = layout.remove_related_see_also('primary', source)
        self.assertNotIn('href="base.html"', result)
        self.assertIn('href="other.html"', result)

    def test_headers_get_real_pages_without_overwriting_authored_prose(self):
        self.write('header-user', 'Rodin::Type')
        path = self.xml / 'header-user.xml'
        source = path.read_text().replace('</compounddef>', '<includes refid="header">Header.h</includes></compounddef>')
        path.write_text(source)
        header = self.xml / 'header.xml'
        header.write_text('<doxygen><compounddef kind="file"><compoundname>Header.h</compoundname><briefdescription/><detaileddescription/></compounddef></doxygen>')
        prepare_reference_headers(self.xml)
        self.assertEqual(text(ET.parse(header).find('./compounddef/briefdescription')), 'Declarations in Header.h.')
        prepare_reference_headers(self.xml)
        self.assertEqual(len(ET.parse(header).findall('.//briefdescription/para')), 1)


if __name__ == '__main__':
    unittest.main()

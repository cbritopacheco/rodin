"""Benchmark navigation preserves ownership, metadata and canonical API overviews."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

from benchmark_site import build_site, inventory, publication_path
from build_documentation import export_module_catalogue


class BenchmarkSiteTest(unittest.TestCase):
    def fixture(self, root):
        docs = root / 'docs'
        docs.mkdir()
        (docs / 'theme.css').write_text(':root { --color: black; }')
        (docs / 'index.html').write_text('<header><nav id="navigation"><a id="m-navbar-brand" href="index.html">Rodin</a><a href="pages.html">Pages</a><a href="#search" class="m-doc-search-icon" title="Search" onclick="return showSearch()"><svg viewBox="0 0 16 16"></svg></a></nav></header><footer><nav>Rodin copyright · Created with Doxygen and m.css</nav></footer>')
        (docs / 'Doxygen.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
        overview = {'brief': 'Matrix operations.', 'representative': '<code>Ax</code>',
                    'description': '<p><a href="class_matrix.html">Matrices</a>.</p>'}
        export_module_catalogue([('Core', [('Math', 'namespace_math')])], {'Core': 'Algebra.'},
                                {'namespace_math': overview}, {'namespace_math': 'Rodin::Math'},
                                ['theme.css'], docs)
        description = root / 'matrix.html'
        description.write_text('<p>Matrix multiplication.</p><h2>Timing boundary</h2><p>Inputs excluded.</p>')
        record = {'target': 'MatrixBenchmark', 'owner': 'Rodin/Math', 'executable': '/build/matrix',
                  'state': 'registered', 'path': 'Rodin/Math/SpatialMatrix',
                  'description': str(description), 'format': 'googlecpp', 'profile': 'slow',
                  'metric': 'cpu_time', 'ranks': '1', 'args': ''}
        manifest = root / 'manifest.tsv'
        self.write_manifest(manifest, [record])
        return docs, manifest, record

    def write_manifest(self, path, records):
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(records)

    def test_canonical_module_export_contains_rendered_content_and_no_build_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            docs, _, _ = self.fixture(Path(folder))
            catalogue = json.loads((docs / 'rodin-module-catalogue.json').read_text())
            module = catalogue['groups'][0]['modules'][0]
            self.assertEqual(module['path'], 'Rodin/Math')
            self.assertEqual(module['representative'], '<code>Ax</code>')
            self.assertNotIn(folder, json.dumps(catalogue))

    def test_documentation_group_retains_module_path_without_physical_context(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            docs, manifest, _ = self.fixture(root)
            catalogue = json.loads((docs / 'rodin-module-catalogue.json').read_text())
            module = catalogue['groups'][0]['modules'][0]
            overview = {key: module[key] for key in ('brief', 'representative', 'description')}
            export_module_catalogue(
                [('Infrastructure', [('Serialization', 'group_serialization')])],
                {'Infrastructure': 'Shared services.'}, {'group_serialization': overview},
                {'group_serialization': 'RodinSerialization'}, ['theme.css'], docs,
                {'group_serialization': 'Rodin/Serialization'})
            out = root / 'site'
            build_site(manifest, docs, out, 'Release', '../')
            index = (out / 'index.html').read_text()
            self.assertIn('href="Rodin/Serialization/index.html"', index)
            contexts = index.split('data-category="Physical contexts"', 1)[1]
            self.assertNotIn('Serialization', contexts)
            self.assertFalse((out / 'RodinSerialization').exists())

    def test_documentation_group_requires_explicit_component_path(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, 'Rodin component path'):
                export_module_catalogue(
                    [('Infrastructure', [('Serialization', 'group_serialization')])],
                    {'Infrastructure': 'Shared services.'}, {'group_serialization': {}},
                    {'group_serialization': 'RodinSerialization'}, [], Path(folder))

    def test_site_has_complete_ancestors_relative_assets_and_correct_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            docs, manifest, _ = self.fixture(root)
            out = root / 'site'
            self.assertEqual(build_site(manifest, docs, out, 'Release / real', 'https://example.org/docs/'), 4)
            leaf = (out / 'Rodin/Math/SpatialMatrix/index.html').read_text()
            self.assertIn('../../../assets/site.css', leaf)
            self.assertIn('cpu_time', leaf)
            self.assertIn('Timing boundary', leaf)
            self.assertIn('data-rodin-content', leaf)
            self.assertIn('<main><article>', leaf)
            self.assertIn('class="m-container m-container-inflatable rodin-benchmark-layout"', leaf)
            self.assertIn('class="m-col-l-10 m-push-l-1" data-rodin-content', leaf)
            index = (out / 'index.html').read_text()
            self.assertIn('https://example.org/docs/class_matrix.html', index)
            self.assertNotIn('Convergence', index)
            self.assertNotIn('class="benchmark-breadcrumb"', index)
            self.assertNotIn('class="benchmark-home-link"', index)
            self.assertIn('</svg><svg class="rodin-benchmark-icon"', index)
            self.assertLess(index.index('class="benchmark-magnifier-icon"'),
                            index.index('class="rodin-benchmark-icon"'))
            self.assertIn('class="benchmark-magnifier-icon"', index)
            self.assertIn('class="benchmark-breadcrumb"', leaf)
            self.assertIn('>Rodin</a>', index)
            self.assertIn('href="https://example.org/docs/pages.html"', index)
            self.assertNotIn('Go back to Rodin docs', index)
            self.assertNotIn('data-benchmark-history', index)
            self.assertIn('Rodin/Math/index.html', index)
            footer = '<footer>' + (docs / 'index.html').read_text().split('<footer>', 1)[1]
            self.assertIn(footer, index)
            self.assertIn(footer, leaf)
            self.assertEqual((out / 'assets/Doxygen.svg').read_bytes(),
                             (docs / 'Doxygen.svg').read_bytes())

    def test_documentation_embeds_have_depth_correct_api_links(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            docs,manifest,_=self.fixture(root)
            build_site(manifest,docs,root / 'site','Published configurations','../')
            self.assertIn('href="../class_matrix.html"',(root / 'site/index.html').read_text())
            self.assertIn('href="../../../class_matrix.html"',(root / 'site/Rodin/Math/index.html').read_text())
            self.assertTrue((root / 'site/assets/charts.js').is_file())

    def test_search_uses_generated_benchmark_metadata_and_pinned_engine(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            docs,manifest,_=self.fixture(root)
            out=root / 'site'
            build_site(manifest,docs,out,'Release','../')
            payload=next((out / 'assets').glob('searchdata-*.js')).read_text()
            catalogue=json.loads(payload.splitlines()[0].removeprefix('window.RodinBenchmarkSearch=').removesuffix(';'))
            self.assertEqual({entry['path'] for entry in catalogue},{'Rodin/Math','Rodin/Math/SpatialMatrix'})
            entry=next(entry for entry in catalogue if entry['path'].endswith('SpatialMatrix'))
            self.assertEqual(entry['target'],'MatrixBenchmark')
            self.assertEqual(entry['brief'],'Matrix multiplication.')
            self.assertNotIn(folder,payload)
            engine=next((out / 'assets').glob('search-*.js')).read_text()
            self.assertIn('base85decode',engine)
            self.assertIn('Benchmark metadata',engine)
            self.assertIn('Find benchmarks',(out / 'index.html').read_text())

    def test_unregistered_targets_remain_in_inventory_without_fake_chart_pages(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            docs, manifest, record = self.fixture(root)
            other = {**record, 'target': 'LegacyBenchmark', 'state': 'unregistered', 'path': ''}
            self.write_manifest(manifest, [record, other])
            build_site(manifest, docs, root / 'site', 'Release', '/docs/')
            rows = json.loads((root / 'site/inventory.json').read_text())
            self.assertEqual(len(rows), 2)
            self.assertFalse((root / 'site/LegacyBenchmark').exists())

    def test_duplicate_publication_paths_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            _, manifest, record = self.fixture(Path(folder))
            self.write_manifest(manifest, [record, {**record, 'target': 'OtherBenchmark'}])
            with self.assertRaisesRegex(ValueError, 'Duplicate benchmark publication path'):
                inventory(manifest)

    def test_invalid_metric_ranks_or_missing_description_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            _, manifest, record = self.fixture(Path(folder))
            for field, value in [('metric', 'iterations'), ('ranks', '0'), ('description', '/absent/file.html')]:
                with self.subTest(field=field):
                    self.write_manifest(manifest, [{**record, field: value}])
                    with self.assertRaises(ValueError):
                        inventory(manifest)

    def test_publication_paths_cannot_escape_or_alias_the_output_tree(self):
        for path in ('../escape', '/absolute', 'Rodin/../escape', 'Rodin//Math', 'Rodin/./Math', ''):
            with self.subTest(path=path), self.assertRaises(ValueError):
                publication_path(path)


if __name__ == '__main__':
    unittest.main()

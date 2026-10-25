"""Configured benchmark discovery distinguishes optional and compatibility targets."""
import csv
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class BenchmarkManifestTest(unittest.TestCase):
    def test_configured_graph_and_test_inventory_are_separate(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            bench = root / 'tests/benchmarks'
            bench.mkdir(parents=True)
            (root / 'main.cpp').write_text('int main() {}\n')
            (bench / 'description.html').write_text('<p>Measured operation.</p>')
            (root / 'CMakeLists.txt').write_text('''cmake_minimum_required(VERSION 3.16)
project(ManifestFixture LANGUAGES CXX)
set(CMAKE_BUILD_TYPE Release)
list(APPEND CMAKE_MODULE_PATH "''' + str(ROOT / 'cmake') + '''")
include(RodinBenchmarkManifest)
add_subdirectory(tests/benchmarks)
include(RodinTestManifest)
''')
            (bench / 'CMakeLists.txt').write_text('''add_executable(Registered EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/main.cpp)
rodin_register_benchmark(Registered PATH Rodin/Math/Vector DESCRIPTION description.html
  FORMAT googlecpp PROFILE regular METRIC cpu_time RANKS 1)
if(OPTIONAL_BACKEND)
  add_executable(Backend EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/main.cpp)
  rodin_register_benchmark(Backend PATH Rodin/MPI/Vector DESCRIPTION description.html
    FORMAT googlecpp PROFILE slow METRIC real_time RANKS 2 4)
endif()
rodin_export_benchmarks()
''')
            build = root / 'build'
            for enabled, expected in [('OFF', {'Registered'}),
                                      ('ON', {'Registered', 'Backend'})]:
                subprocess.run(['cmake', '-S', str(root), '-B', str(build),
                                '-DOPTIONAL_BACKEND=' + enabled], check=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                with (build / 'RodinBenchmarkTargets-Release.tsv').open() as stream:
                    rows = {row['target']: row for row in csv.DictReader(stream, delimiter='\t')}
                self.assertEqual(set(rows), expected)
                self.assertEqual(rows['Registered']['path'], 'Rodin/Math/Vector')
                self.assertEqual(rows['Registered']['metric'], 'cpu_time')
                self.assertNotIn('Registered', (build / 'RodinTestTargets-Release.tsv').read_text())

    def test_all_current_benchmark_targets_survive_configuration(self):
        import re
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'tests').symlink_to(ROOT / 'tests', target_is_directory=True)
            (root / 'src').symlink_to(ROOT / 'src', target_is_directory=True)
            sources = '\n'.join(p.read_text() for p in (ROOT / 'tests/benchmarks').rglob('CMakeLists.txt'))
            aliases = sorted(set(re.findall(r'\b(?:Rodin|GTest|benchmark|Boost|Threads)::\w+', sources)))
            declarations = set(re.findall(r'add_executable\(\s*(\w+)', sources))
            cmake = ['cmake_minimum_required(VERSION 3.16)', 'project(CompleteManifest LANGUAGES CXX)',
                     'set(CMAKE_BUILD_TYPE Release)', 'set(RODIN_USE_MPI ON)', 'set(RODIN_USE_PETSC ON)',
                     'list(APPEND CMAKE_MODULE_PATH "' + str(ROOT / 'cmake') + '")']
            for index, alias in enumerate(aliases):
                cmake += ['add_library(dependency' + str(index) + ' INTERFACE)',
                          'add_library(' + alias + ' ALIAS dependency' + str(index) + ')']
            cmake += ['add_subdirectory(tests/benchmarks)']
            (root / 'CMakeLists.txt').write_text('\n'.join(cmake))
            subprocess.run(['cmake', '-S', str(root), '-B', str(root / 'build')], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            with (root / 'build/RodinBenchmarkTargets-Release.tsv').open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual({row['target'] for row in rows}, declarations)
            self.assertEqual(len(rows), len(declarations))
            self.assertNotIn('RodinBenchmarks', declarations)
            self.assertEqual({row['target'] for row in rows if row['state'] == 'registered'},
                             declarations)


    def test_documentation_catalogue_needs_no_executable_or_cpp_compiler(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'tests').symlink_to(ROOT / 'tests', target_is_directory=True)
            (root / 'CMakeLists.txt').write_text('cmake_minimum_required(VERSION 3.16)\n'
                'project(Descriptions LANGUAGES NONE)\n'
                'list(APPEND CMAKE_MODULE_PATH "' + str(ROOT / 'cmake') + '")\n'
                'include(RodinBenchmarkManifest)\n'
                'rodin_export_benchmark_descriptions("${CMAKE_BINARY_DIR}/descriptions.tsv")\n')
            subprocess.run(['cmake', '-S', str(root), '-B', str(root / 'build')], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            with (root / 'build/descriptions.tsv').open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertGreater(len(rows), 70)
            self.assertTrue(all(row['description'] and Path(row['description']).is_file() for row in rows))
            self.assertTrue(all(row['metric'] in ('cpu_time', 'real_time', 'auto') for row in rows))
            self.assertEqual(len({row['path'] for row in rows}), len(rows))


if __name__ == '__main__':
    unittest.main()

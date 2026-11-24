"""Configure the real module target graph independently of optional dependencies."""
import itertools
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class KokkosModuleTest(unittest.TestCase):
    def test_mpi_geometry_never_requires_petsc(self):
        for mpi, petsc in itertools.product([False, True], repeat=2):
            with self.subTest(mpi=mpi, petsc=petsc), tempfile.TemporaryDirectory() as directory:
                source = Path(directory)
                dependencies = ['Rodin::Geometry', 'Rodin::Context', 'Rodin::Variational',
                                'Rodin::Assembly', 'Rodin::Solver', 'Rodin::MPI',
                                'Rodin::PETSc', 'Kokkos::kokkos', 'Kokkos::kokkoskernels']
                (source / 'CMakeLists.txt').write_text(
                    'cmake_minimum_required(VERSION 3.16)\nproject(Modules LANGUAGES CXX)\n'
                    + ''.join(f'add_library({name} INTERFACE IMPORTED)\n' for name in dependencies)
                    + 'function(rodin_install_headers)\nendfunction()\n'
                    + 'set(RODIN_USE_KOKKOS ON)\n'
                    + f'set(RODIN_USE_MPI {"ON" if mpi else "OFF"})\n'
                    + f'set(RODIN_USE_PETSC {"ON" if petsc else "OFF"})\n'
                    + f'add_subdirectory("{ROOT / "src/Rodin/Kokkos"}" kokkos)\n'
                    + 'foreach(name RodinKokkos RodinKokkosMPI RodinKokkosPETSc)\n'
                    + '  if(TARGET ${name})\n'
                    + '    get_target_property(links ${name} INTERFACE_LINK_LIBRARIES)\n'
                    + '    file(APPEND "${CMAKE_BINARY_DIR}/links.txt" "${name}=${links}\\n")\n'
                    + '  endif()\nendforeach()\n')
                build = source / 'build'
                subprocess.run(['cmake', '-S', str(source), '-B', str(build)],
                               check=True, capture_output=True, text=True)
                graph = dict(line.split('=', 1) for line in (build / 'links.txt').read_text().splitlines())
                self.assertEqual('RodinKokkosMPI' in graph, mpi)
                self.assertEqual('RodinKokkosPETSc' in graph, petsc)
                self.assertNotIn('PETSc', graph['RodinKokkos'])
                if mpi:
                    self.assertEqual(graph['RodinKokkosMPI'].split(';'),
                                     ['Rodin::Kokkos', 'Rodin::MPI'])
                if petsc:
                    self.assertIn('Rodin::PETSc', graph['RodinKokkosPETSc'])
                    self.assertEqual('Rodin::Kokkos::MPI' in graph['RodinKokkosPETSc'], mpi)


if __name__ == '__main__':
    unittest.main()

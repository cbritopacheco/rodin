"""Exercise component discovery through installed CMake exports and package templates."""
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class InstalledPackageTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='rodin-package-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.prefix = self.root / 'install'
        source = self.root / 'producer'
        source.mkdir()
        (source / 'CMakeLists.txt').write_text('''
cmake_minimum_required(VERSION 3.16)
project(Export LANGUAGES NONE)
set(RODIN_USE_MPI ON)
set(RODIN_USE_KOKKOS ON)
set(RODIN_USE_PETSC ON)
set(RODIN_USE_SCOTCH ON)
foreach(name Rodin Geometry MPI PETSc Scotch Kokkos KokkosMPI KokkosPETSc)
  add_library(${name} INTERFACE)
  set_target_properties(${name} PROPERTIES EXPORT_NAME ${name})
endforeach()
set_target_properties(KokkosMPI PROPERTIES EXPORT_NAME Kokkos::MPI)
set_target_properties(KokkosPETSc PROPERTIES EXPORT_NAME Kokkos::PETSc)
target_link_libraries(Geometry INTERFACE Rodin)
target_link_libraries(MPI INTERFACE Geometry MPI::MPI_CXX Boost::mpi)
target_link_libraries(PETSc INTERFACE Geometry MPI PETSc::PETSc)
target_link_libraries(Scotch INTERFACE Geometry Scotch::Scotch)
target_link_libraries(Kokkos INTERFACE Geometry Kokkos::kokkos Kokkos::kokkoskernels)
# Exercise wrapped links and a cycle in the imported target closure.
target_link_libraries(KokkosMPI INTERFACE "$<LINK_ONLY:MPI>" Kokkos)
target_link_libraries(Geometry INTERFACE "$<LINK_ONLY:Rodin>")
target_link_libraries(KokkosPETSc INTERFACE KokkosMPI PETSc)
install(TARGETS Rodin Geometry MPI PETSc Scotch Kokkos KokkosMPI KokkosPETSc EXPORT RodinTargets)
install(EXPORT RodinTargets NAMESPACE Rodin:: DESTINATION lib/cmake/Rodin)
configure_file("''' + str(ROOT / 'cmake/RodinConfig.cmake.in') + '''"
  "${CMAKE_CURRENT_BINARY_DIR}/RodinConfig.cmake" @ONLY)
configure_file("''' + str(ROOT / 'cmake/RodinDependencies.cmake.in') + '''"
  "${CMAKE_CURRENT_BINARY_DIR}/RodinDependencies.cmake" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/RodinConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/RodinDependencies.cmake" DESTINATION lib/cmake/Rodin)
''')
        self.run_command(['cmake', '-S', str(source), '-B', str(self.root / 'producer-build')])
        self.run_command(['cmake', '--install', str(self.root / 'producer-build'),
                          '--prefix', str(self.prefix)])
        modules = self.prefix / 'lib/cmake/Rodin'
        # Dependency doubles make discovery observable and do not provide Rodin
        # targets; those come from the actual install(EXPORT) output above.
        for package, targets in {
            'Boost': ['Boost::filesystem', 'Boost::serialization', 'Boost::mpi'],
            'HDF5': ['HDF5::HDF5'], 'Threads': ['Threads::Threads'],
            'Eigen3': ['Eigen3::Eigen'], 'MPI': ['MPI::MPI_CXX'],
            'PETSc': ['PETSc::PETSc'], 'Scotch': ['Scotch::Scotch']
        }.items():
            text = ('file(APPEND "${CMAKE_BINARY_DIR}/discovered.txt" "' + package + '\\n")\n'
                    + f'set({package}_FOUND TRUE)\n')
            if package == 'MPI':
                text += ('if(NOT MPI_FIND_COMPONENTS STREQUAL "CXX")\n'
                         '  message(FATAL_ERROR "Only CXX MPI is required")\nendif()\n')
            for target in targets:
                text += f'if(NOT TARGET {target})\n  add_library({target} INTERFACE IMPORTED)\nendif()\n'
            (modules / f'Find{package}.cmake').write_text(text)
        for package, targets in {'Kokkos': ['Kokkos::kokkos'],
                                 'KokkosKernels': ['Kokkos::kokkoskernels']}.items():
            text = f'file(APPEND "${{CMAKE_BINARY_DIR}}/discovered.txt" "{package}\\n")\n'
            if package == 'Kokkos':
                text += 'set(Kokkos_ENABLE_SERIAL ON)\n'
            for target in targets:
                text += f'if(NOT TARGET {target})\n  add_library({target} INTERFACE IMPORTED)\nendif()\n'
            (modules / f'{package}Config.cmake').write_text(text)
            (modules / f'{package}ConfigVersion.cmake').write_text(
                'set(PACKAGE_VERSION 5.2.2)\nset(PACKAGE_VERSION_COMPATIBLE TRUE)\n')

    def run_command(self, command, success=True):
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def consume(self, body, *, success=True, blocked=False):
        source = self.root / 'consumer'
        source.mkdir(exist_ok=True)
        (source / 'CMakeLists.txt').write_text(
            'cmake_minimum_required(VERSION 3.16)\nproject(Consumer LANGUAGES NONE)\n'
            'set(_original_module_path "${CMAKE_MODULE_PATH}")\n' + body
            + '\nif(NOT "${CMAKE_MODULE_PATH}" STREQUAL "${_original_module_path}")\n'
              '  message(FATAL_ERROR "Package discovery changed CMAKE_MODULE_PATH")\nendif()\n')
        build = self.root / 'consumer-build'
        command = ['cmake', '-S', str(source), '-B', str(build),
                   '-DCMAKE_PREFIX_PATH=' + str(self.prefix),
                   '-DKokkos_DIR=' + str(self.prefix / 'lib/cmake/Rodin'),
                   '-DKokkosKernels_DIR=' + str(self.prefix / 'lib/cmake/Rodin')]
        if blocked:
            command.append('-DCMAKE_DISABLE_FIND_PACKAGE_PETSc=TRUE')
        result = self.run_command(command, success=success)
        log = build / 'discovered.txt'
        return result, set(log.read_text().splitlines()) if log.exists() else set()

    def test_installed_mpi_component_does_not_discover_petsc(self):
        _, dependencies = self.consume('find_package(Rodin REQUIRED COMPONENTS Kokkos::MPI)', blocked=True)
        self.assertIn('MPI', dependencies)
        self.assertIn('KokkosKernels', dependencies)
        self.assertNotIn('PETSc', dependencies)
        self.assertNotIn('Scotch', dependencies)

    def test_scotch_is_discovered_only_for_its_component(self):
        _, dependencies = self.consume('find_package(Rodin REQUIRED COMPONENTS Scotch)', blocked=True)
        self.assertIn('Scotch', dependencies)
        self.assertFalse({'PETSc', 'MPI', 'Kokkos', 'KokkosKernels'} & dependencies)

    def test_core_component_does_not_discover_optional_backends(self):
        _, dependencies = self.consume('find_package(Rodin REQUIRED COMPONENTS Geometry)', blocked=True)
        self.assertFalse({'PETSc', 'MPI', 'Kokkos', 'KokkosKernels'} & dependencies)

    def test_petsc_component_and_legacy_call_discover_the_full_dependencies(self):
        for arguments in ['COMPONENTS Kokkos::PETSc', '']:
            with self.subTest(arguments=arguments):
                _, dependencies = self.consume('find_package(Rodin REQUIRED ' + arguments + ')')
                self.assertTrue({'PETSc', 'MPI', 'Kokkos', 'KokkosKernels'} <= dependencies)

    def test_required_petsc_dependency_is_not_silently_ignored(self):
        result, _ = self.consume('find_package(Rodin REQUIRED COMPONENTS Kokkos::PETSc)',
                                 success=False, blocked=True)
        self.assertIn('PETSc', result.stdout + result.stderr)

    def test_legacy_dependency_variables_remain_visible(self):
        self.consume('''find_package(Rodin REQUIRED)
if(NOT Kokkos_ENABLE_SERIAL OR NOT MPI_FOUND OR NOT PETSc_FOUND)
  message(FATAL_ERROR "Legacy dependency variables were hidden")
endif()''')

    def test_optional_missing_backend_does_not_break_required_components(self):
        self.consume('''find_package(Rodin REQUIRED COMPONENTS Kokkos::MPI
  OPTIONAL_COMPONENTS Kokkos::PETSc NotInstalled)
if(NOT Rodin_FOUND OR NOT Rodin_Kokkos::MPI_FOUND OR Rodin_Kokkos::PETSc_FOUND OR Rodin_NotInstalled_FOUND)
  message(FATAL_ERROR "Incorrect component result")
endif()''', blocked=True)

    def test_unknown_required_component_is_reported(self):
        result, _ = self.consume('find_package(Rodin REQUIRED COMPONENTS NotInstalled)', success=False)
        self.assertIn("component 'NotInstalled' is not installed", result.stdout + result.stderr)

    def test_repeated_discovery_can_add_a_backend_and_recover_after_failure(self):
        self.consume('''find_package(Rodin REQUIRED COMPONENTS Kokkos::MPI)
find_package(Rodin QUIET COMPONENTS NotInstalled)
if(Rodin_FOUND)
  message(FATAL_ERROR "Missing component accepted")
endif()
find_package(Rodin REQUIRED COMPONENTS Kokkos::PETSc)
if(NOT Rodin_FOUND OR NOT Rodin_Kokkos::PETSc_FOUND OR NOT TARGET PETSc::PETSc)
  message(FATAL_ERROR "Repeated discovery did not load the adapter")
endif()''')


if __name__ == '__main__':
    unittest.main()

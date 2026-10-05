add_library(RodinConvergence INTERFACE)
target_include_directories(RodinConvergence
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(RodinConvergence
  INTERFACE
    Rodin::Geometry
    Rodin::QF
    Rodin::Variational
)

add_executable(RodinConvergenceUtilitiesTest ConvergenceTest.cpp ErrorNormTest.cpp)
target_link_libraries(RodinConvergenceUtilitiesTest
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
)
gtest_discover_tests(RodinConvergenceUtilitiesTest
  PROPERTIES
    LABELS convergence
)

if (RODIN_USE_PETSC)
  include(CheckCXXSourceCompiles)
  include(CMakePushCheckState)
  cmake_push_check_state(RESET)
  set(CMAKE_REQUIRED_LIBRARIES PETSc::PETSc)
  check_cxx_source_compiles("#include <petsc.h>
    #ifndef PETSC_USE_COMPLEX
    #error Complex PETSc required
    #endif
    int main() { return 0; }" RODIN_CONVERGENCE_PETSC_COMPLEX)
  if (NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
    check_cxx_source_compiles("#include <petsc.h>
      #ifndef PETSC_HAVE_MUMPS
      #error MUMPS required for Stokes factorization
      #endif
      int main() { return 0; }" RODIN_CONVERGENCE_PETSC_MUMPS)
    if (NOT RODIN_CONVERGENCE_PETSC_MUMPS)
      message(STATUS "PETSc Stokes convergence excluded: PETSc lacks MUMPS")
    endif()
  endif()
  cmake_pop_check_state()
endif()

add_subdirectory(p)
add_subdirectory(h)
add_subdirectory(hp)
add_subdirectory(isoparametric)

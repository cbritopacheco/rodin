add_executable(RodinConvergenceHPETScPoisson Poisson.cpp)
target_link_libraries(RodinConvergenceHPETScPoisson
  PRIVATE
    GTest::gtest
    RodinConvergenceH
    Rodin::PETSc
)
add_test(NAME RodinConvergenceHPETScPoisson
  COMMAND $<TARGET_FILE:RodinConvergenceHPETScPoisson>)
set_tests_properties(RodinConvergenceHPETScPoisson PROPERTIES
  LABELS "convergence;petsc")
rodin_suppress_external_mpi_lsan_for_test(RodinConvergenceHPETScPoisson)

include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPETScPoissonBoundary BoundaryConditions.cpp)
if (RODIN_CONVERGENCE_PETSC_MUMPS)
  target_compile_definitions(RodinConvergenceHPETScPoissonBoundary
    PRIVATE RODIN_POISSON_BOUNDARY_MUMPS)
else()
  message(STATUS "Pure Neumann Poisson convergence excluded: PETSc lacks MUMPS")
endif()

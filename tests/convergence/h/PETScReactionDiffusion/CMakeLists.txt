add_executable(RodinConvergenceHPETScReactionDiffusion ReactionDiffusion.cpp)
include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPETScReactionDiffusionBoundary BoundaryConditions.cpp)
target_link_libraries(RodinConvergenceHPETScReactionDiffusion
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
set(_RODIN_REACTIONDIFFUSION_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_REACTIONDIFFUSION_GEOMETRIES)
  add_test(NAME RodinConvergenceHPETScReactionDiffusion_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPETScReactionDiffusion>
      "--gtest_filter=AllGeometries/PETScReactionDiffusionLocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPETScReactionDiffusion_${geometry} PROPERTIES
    LABELS "convergence;petsc;slow" TIMEOUT 600)
  rodin_suppress_external_mpi_lsan_for_test(
    RodinConvergenceHPETScReactionDiffusion_${geometry})
endforeach()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergenceHPETScReactionDiffusion PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    foreach(geometry IN LISTS _RODIN_REACTIONDIFFUSION_GEOMETRIES)
      add_test(NAME RodinConvergenceHPETScMPIReactionDiffusion_np${np}_${geometry}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergenceHPETScReactionDiffusion>
          "--gtest_filter=AllGeometries/PETScReactionDiffusionMPITest.*/${geometry}")
      set_tests_properties(RodinConvergenceHPETScMPIReactionDiffusion_np${np}_${geometry}
        PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 600 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(
        RodinConvergenceHPETScMPIReactionDiffusion_np${np}_${geometry})
    endforeach()
  endforeach()
endif()

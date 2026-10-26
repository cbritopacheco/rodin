include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPETScLinearElasticityBoundary Traction.cpp
  RATE_CASES MixedTractionP1Rates MixedTractionP2Rates NearlyIncompressibleP2Rates)
# A single tetrahedral P1 history exceeded thirty minutes in CI. Give each
# three-dimensional rate study its own three-hour scheduling allowance.
foreach(geometry Tetrahedron Pyramid Hexahedron Wedge)
  foreach(study MixedTractionP1Rates MixedTractionP2Rates NearlyIncompressibleP2Rates)
    set_tests_properties(RodinConvergenceHPETScLinearElasticityBoundary_${study}_${geometry}
      PROPERTIES TIMEOUT 10800)
    if(RODIN_USE_MPI)
      foreach(np 1 2 3 4)
        set_tests_properties(RodinConvergenceHPETScLinearElasticityBoundary_MPI_np${np}_${study}_${geometry}
          PROPERTIES TIMEOUT 10800)
      endforeach()
    endif()
  endforeach()
endforeach()

add_executable(RodinConvergenceHPETScLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinConvergenceHPETScLinearElasticity
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
set(_RODIN_ELASTICITY_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_ELASTICITY_GEOMETRIES)
  add_test(NAME RodinConvergenceHPETScLinearElasticity_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPETScLinearElasticity>
      "--gtest_filter=AllGeometries/PETScLinearElasticityLocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPETScLinearElasticity_${geometry} PROPERTIES
    LABELS "convergence;petsc;slow" TIMEOUT 600)
  rodin_suppress_external_mpi_lsan_for_test(
    RodinConvergenceHPETScLinearElasticity_${geometry})
endforeach()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergenceHPETScLinearElasticity PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    foreach(geometry IN LISTS _RODIN_ELASTICITY_GEOMETRIES)
      add_test(NAME RodinConvergenceHPETScMPILinearElasticity_np${np}_${geometry}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergenceHPETScLinearElasticity>
          "--gtest_filter=AllGeometries/PETScLinearElasticityMPITest.*/${geometry}")
      set_tests_properties(RodinConvergenceHPETScMPILinearElasticity_np${np}_${geometry}
        PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 600 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(
        RodinConvergenceHPETScMPILinearElasticity_np${np}_${geometry})
    endforeach()
  endforeach()
endif()

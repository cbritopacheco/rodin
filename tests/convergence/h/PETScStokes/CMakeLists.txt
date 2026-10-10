include(${CMAKE_CURRENT_LIST_DIR}/../../PETScRefinement.cmake)
add_executable(RodinConvergenceHPETScStokesStability Stability.cpp)
target_link_libraries(RodinConvergenceHPETScStokesStability
  PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
add_test(NAME RodinConvergenceHPETScStokesStability_Utility
  COMMAND $<TARGET_FILE:RodinConvergenceHPETScStokesStability>
    "--gtest_filter=PETScMixedStabilityTest.*")
set_tests_properties(RodinConvergenceHPETScStokesStability_Utility
  PROPERTIES LABELS "convergence;petsc" TIMEOUT 60)
rodin_suppress_external_mpi_lsan_for_test(RodinConvergenceHPETScStokesStability_Utility)
foreach(geometry Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  set(_stability_test RodinConvergenceHPETScStokesStability_${geometry})
  add_test(NAME ${_stability_test}
    COMMAND $<TARGET_FILE:RodinConvergenceHPETScStokesStability>
      "--gtest_filter=AllGeometries/LocalStabilityTest.*/${geometry}-*P4P3*")
  set_tests_properties(${_stability_test}
    PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
  rodin_suppress_external_mpi_lsan_for_test(${_stability_test})
  foreach(group Rates Controls)
    set(_case P4P3PressureSpectrumAcrossRefinementLevels)
    if (group STREQUAL "Controls")
      set(_case P4P3MissingDivergenceRejected)
    endif()
    set(_test RodinConvergenceHPETScStokesStability_P4P3_${group}_${geometry})
    add_test(NAME ${_test}
      COMMAND $<TARGET_FILE:RodinConvergenceHPETScStokesStability>
        "--gtest_filter=AllGeometries/LocalStabilityTest.${_case}/${geometry}")
    set(_timeout 1800)
    # Rank-zero and local verification use the same independent spectral oracle.
    if (geometry STREQUAL "Pyramid" AND group STREQUAL "Rates")
      set(_timeout 3600)
    endif()
    set_tests_properties(${_test} PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT ${_timeout})
    rodin_suppress_external_mpi_lsan_for_test(${_test})
  endforeach()
endforeach()
rodin_add_petsc_refinement(RodinConvergenceHPETScStokesBoundary Traction.cpp
  GEOMETRIES Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge
  RATE_CASES P2P1Rates P3P2Rates)

add_executable(RodinConvergenceHPETScStokes Stokes.cpp)
target_link_libraries(RodinConvergenceHPETScStokes
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
set(_RODIN_STOKES_GEOMETRIES
  Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_STOKES_GEOMETRIES)
  add_test(NAME RodinConvergenceHPETScStokes_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPETScStokes>
      "--gtest_filter=AllGeometries/PETScStokesLocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPETScStokes_${geometry} PROPERTIES
    LABELS "convergence;petsc;slow" TIMEOUT 600)
  rodin_suppress_external_mpi_lsan_for_test(
    RodinConvergenceHPETScStokes_${geometry})
endforeach()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergenceHPETScStokesStability PRIVATE Rodin::MPI)
  target_link_libraries(RodinConvergenceHPETScStokes PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    set(_stability_utility RodinConvergenceHPETScMPIStokesStability_Utility_np${np})
    add_test(NAME ${_stability_utility}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
        ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
        $<TARGET_FILE:RodinConvergenceHPETScStokesStability>
        "--gtest_filter=PETScMixedStabilityTest.*")
    set_tests_properties(${_stability_utility} PROPERTIES
      LABELS "convergence;petsc;distributed" TIMEOUT 60 PROCESSORS ${np})
    rodin_suppress_external_mpi_lsan_for_test(${_stability_utility})
    foreach(geometry IN LISTS _RODIN_STOKES_GEOMETRIES)
      set(_stability_test RodinConvergenceHPETScMPIStokesStability_np${np}_${geometry})
      add_test(NAME ${_stability_test}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergenceHPETScStokesStability>
          "--gtest_filter=AllGeometries/MPIStabilityTest.*/${geometry}-*P4P3*")
      set_tests_properties(${_stability_test} PROPERTIES
        LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(${_stability_test})
      foreach(group Rates Controls)
        set(_case P4P3PressureSpectrumAcrossRefinementLevels)
        if (group STREQUAL "Controls")
          set(_case P4P3MissingDivergenceRejected)
        endif()
        set(_test RodinConvergenceHPETScMPIStokesStability_P4P3_np${np}_${group}_${geometry})
        add_test(NAME ${_test}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceHPETScStokesStability>
            "--gtest_filter=AllGeometries/MPIStabilityTest.${_case}/${geometry}")
        set(_timeout 1800)
        if (geometry STREQUAL "Pyramid" AND group STREQUAL "Rates")
          set(_timeout 3600)
        endif()
        set_tests_properties(${_test} PROPERTIES
          LABELS "convergence;petsc;distributed;slow" TIMEOUT ${_timeout} PROCESSORS ${np})
        rodin_suppress_external_mpi_lsan_for_test(${_test})
      endforeach()
      add_test(NAME RodinConvergenceHPETScMPIStokes_np${np}_${geometry}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergenceHPETScStokes>
          "--gtest_filter=AllGeometries/PETScStokesMPITest.*/${geometry}")
      set_tests_properties(RodinConvergenceHPETScMPIStokes_np${np}_${geometry}
        PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 600 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(
        RodinConvergenceHPETScMPIStokes_np${np}_${geometry})
    endforeach()
  endforeach()
endif()

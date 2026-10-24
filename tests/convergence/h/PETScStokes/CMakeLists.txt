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

# Independent hierarchy points retain both quadrature orders and all spectral
# diagnostics. The same partition is used by local and distributed contexts.
function(rodin_add_petsc_stokes_pressure_cases geometry rank)
  set(_prefix RodinConvergenceHPETScStokesStability_P4P3)
  set(_suite LocalStabilityTest)
  set(_labels "convergence;petsc;slow")
  set(_command $<TARGET_FILE:RodinConvergenceHPETScStokesStability>)
  if(rank)
    set(_prefix RodinConvergenceHPETScMPIStokesStability_P4P3_np${rank})
    set(_suite MPIStabilityTest)
    set(_labels "convergence;petsc;distributed;slow")
    list(PREPEND _command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank}
      ${_RODIN_MPI_OVERSUBSCRIBE_FLAG})
  endif()
  foreach(point AffineN3 AffineN4 AffineN5 CurvedN3 CurvedN4 CurvedN5 Controls)
    set(_group Rates_${point})
    set(_case P4P3PressureSpectrum${point})
    set(_timeout 3600)
    if(point STREQUAL "Controls")
      set(_group Controls)
      set(_case P4P3MissingDivergenceRejected)
      set(_timeout 1800)
    elseif(geometry STREQUAL "Pyramid" AND point MATCHES "N5$")
      # Two independent order evaluations exceed the former one-hour total.
      set(_timeout 10800)
    endif()
    set(_test ${_prefix}_${_group}_${geometry})
    add_test(NAME ${_test} COMMAND ${_command}
      "--gtest_filter=AllGeometries/${_suite}.${_case}/${geometry}")
    set_tests_properties(${_test} PROPERTIES LABELS "${_labels}" TIMEOUT ${_timeout})
    if(rank)
      set_tests_properties(${_test} PROPERTIES PROCESSORS ${rank})
    endif()
    rodin_suppress_external_mpi_lsan_for_test(${_test})
  endforeach()
endfunction()

foreach(geometry Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  set(_stability_test RodinConvergenceHPETScStokesStability_${geometry})
  add_test(NAME ${_stability_test}
    COMMAND $<TARGET_FILE:RodinConvergenceHPETScStokesStability>
      "--gtest_filter=AllGeometries/LocalStabilityTest.*/${geometry}-*P4P3*")
  set_tests_properties(${_stability_test}
    PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
  rodin_suppress_external_mpi_lsan_for_test(${_stability_test})
  rodin_add_petsc_stokes_pressure_cases(${geometry} 0)
endforeach()
rodin_add_petsc_refinement(RodinConvergenceHPETScStokesBoundary Traction.cpp
  GEOMETRIES Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge
  RATE_CASES P2P1Rates P3P2Rates)


# Coverage CI exceeds 30 minutes for the complete degree-three Pyramid
# hierarchy, and for threaded Tetrahedron/Hexahedron histories. These are
# process budgets only: all three grids, norms and MUMPS checks are retained.
foreach(geometry Tetrahedron Pyramid Hexahedron Wedge)
  set(_stokes_boundary_timeout 3600)
  if (geometry STREQUAL "Pyramid")
    set(_stokes_boundary_timeout 7200)
  endif()
  set_tests_properties(RodinConvergenceHPETScStokesBoundary_P3P2Rates_${geometry}
    PROPERTIES TIMEOUT ${_stokes_boundary_timeout})
  if (RODIN_USE_MPI)
    foreach(np 1 2 3 4)
      set_tests_properties(RodinConvergenceHPETScStokesBoundary_MPI_np${np}_P3P2Rates_${geometry}
        PROPERTIES TIMEOUT ${_stokes_boundary_timeout})
    endforeach()
  endif()
endforeach()

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
      rodin_add_petsc_stokes_pressure_cases(${geometry} ${np})
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

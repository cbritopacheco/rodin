add_executable(RodinConvergenceIsoparametricNonlinearPoisson NonlinearPoisson.cpp)
target_link_libraries(RodinConvergenceIsoparametricNonlinearPoisson
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_NONLINEAR_POISSON_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_CURVED_NONLINEAR_POISSON_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricNonlinearPoisson_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricNonlinearPoisson>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceIsoparametricNonlinearPoisson_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricNonlinearPoisson_${geometry}
      PROPERTIES RESOURCE_LOCK curved_nonlinear_poisson_pyramid)
  endif()
endforeach()

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_executable(RodinConvergenceIsoparametricPETScNonlinearPoisson NonlinearPoisson.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScNonlinearPoisson
    PRIVATE RODIN_CURVED_NONLINEAR_POISSON_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScNonlinearPoisson
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(geometry IN LISTS _RODIN_CURVED_NONLINEAR_POISSON_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScNonlinearPoisson_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScNonlinearPoisson>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}")
    set_tests_properties(RodinConvergenceIsoparametricPETScNonlinearPoisson_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScNonlinearPoisson_${geometry}
        PROPERTIES RESOURCE_LOCK curved_nonlinear_poisson_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScNonlinearPoisson_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScNonlinearPoisson PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_CURVED_NONLINEAR_POISSON_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPINonlinearPoisson_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScNonlinearPoisson>
            "--gtest_filter=AllGeometries/MPITest.*/${geometry}")
        set_tests_properties(
          RodinConvergenceIsoparametricPETScMPINonlinearPoisson_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(
            RodinConvergenceIsoparametricPETScMPINonlinearPoisson_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK curved_nonlinear_poisson_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPINonlinearPoisson_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()

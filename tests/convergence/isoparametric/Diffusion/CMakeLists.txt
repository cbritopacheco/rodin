add_executable(RodinConvergenceIsoparametricDiffusion Diffusion.cpp)
target_link_libraries(RodinConvergenceIsoparametricDiffusion
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_DIFFUSION_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_CURVED_DIFFUSION_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricDiffusion_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricDiffusion>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*Sine*:*Lifted*")
  add_test(NAME RodinConvergenceIsoparametricApproximatedDiffusion_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricDiffusion>
      "--gtest_filter=AllGeometries/LocalTest.Sine*/${geometry}")
  add_test(NAME RodinConvergenceIsoparametricLiftedDiffusion_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricDiffusion>
      "--gtest_filter=AllGeometries/LocalTest.Lifted*/${geometry}-*LiftedSmooth*")
  add_test(NAME RodinConvergenceIsoparametricSmoothLiftedDiffusion_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricDiffusion>
      "--gtest_filter=AllGeometries/LocalTest.LiftedSmooth*/${geometry}")
  set_tests_properties(RodinConvergenceIsoparametricDiffusion_${geometry}
    RodinConvergenceIsoparametricApproximatedDiffusion_${geometry}
    RodinConvergenceIsoparametricLiftedDiffusion_${geometry}
    RodinConvergenceIsoparametricSmoothLiftedDiffusion_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricDiffusion_${geometry}
      RodinConvergenceIsoparametricApproximatedDiffusion_${geometry}
      RodinConvergenceIsoparametricLiftedDiffusion_${geometry}
      RodinConvergenceIsoparametricSmoothLiftedDiffusion_${geometry}
      PROPERTIES RESOURCE_LOCK curved_diffusion_pyramid)
  endif()
endforeach()

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_executable(RodinConvergenceIsoparametricPETScDiffusion Diffusion.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScDiffusion
    PRIVATE RODIN_CURVED_DIFFUSION_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScDiffusion
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(geometry IN LISTS _RODIN_CURVED_DIFFUSION_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScDiffusion_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScDiffusion>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*Sine*:*Lifted*")
    add_test(NAME RodinConvergenceIsoparametricPETScApproximatedDiffusion_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScDiffusion>
        "--gtest_filter=AllGeometries/LocalTest.Sine*/${geometry}")
    add_test(NAME RodinConvergenceIsoparametricPETScLiftedDiffusion_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScDiffusion>
        "--gtest_filter=AllGeometries/LocalTest.Lifted*/${geometry}-*LiftedSmooth*")
    add_test(NAME RodinConvergenceIsoparametricPETScSmoothLiftedDiffusion_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScDiffusion>
        "--gtest_filter=AllGeometries/LocalTest.LiftedSmooth*/${geometry}")
    set_tests_properties(RodinConvergenceIsoparametricPETScDiffusion_${geometry}
      RodinConvergenceIsoparametricPETScApproximatedDiffusion_${geometry}
      RodinConvergenceIsoparametricPETScLiftedDiffusion_${geometry}
      RodinConvergenceIsoparametricPETScSmoothLiftedDiffusion_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScDiffusion_${geometry}
        RodinConvergenceIsoparametricPETScApproximatedDiffusion_${geometry}
        RodinConvergenceIsoparametricPETScLiftedDiffusion_${geometry}
        RodinConvergenceIsoparametricPETScSmoothLiftedDiffusion_${geometry}
        PROPERTIES RESOURCE_LOCK curved_diffusion_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScDiffusion_${geometry})
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScApproximatedDiffusion_${geometry})
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScLiftedDiffusion_${geometry})
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScSmoothLiftedDiffusion_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScDiffusion PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_CURVED_DIFFUSION_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPIDiffusion_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScDiffusion>
            "--gtest_filter=AllGeometries/MPITest.*/${geometry}-*Sine*:*Lifted*")
        add_test(NAME RodinConvergenceIsoparametricPETScMPIApproximatedDiffusion_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScDiffusion>
            "--gtest_filter=AllGeometries/MPITest.Sine*/${geometry}")
        add_test(NAME RodinConvergenceIsoparametricPETScMPILiftedDiffusion_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScDiffusion>
            "--gtest_filter=AllGeometries/MPITest.Lifted*/${geometry}-*LiftedSmooth*")
        add_test(NAME RodinConvergenceIsoparametricPETScMPISmoothLiftedDiffusion_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScDiffusion>
            "--gtest_filter=AllGeometries/MPITest.LiftedSmooth*/${geometry}")
        set_tests_properties(
          RodinConvergenceIsoparametricPETScMPIDiffusion_np${np}_${geometry}
          RodinConvergenceIsoparametricPETScMPIApproximatedDiffusion_np${np}_${geometry}
          RodinConvergenceIsoparametricPETScMPILiftedDiffusion_np${np}_${geometry}
          RodinConvergenceIsoparametricPETScMPISmoothLiftedDiffusion_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(
            RodinConvergenceIsoparametricPETScMPIDiffusion_np${np}_${geometry}
            RodinConvergenceIsoparametricPETScMPIApproximatedDiffusion_np${np}_${geometry}
            RodinConvergenceIsoparametricPETScMPILiftedDiffusion_np${np}_${geometry}
            RodinConvergenceIsoparametricPETScMPISmoothLiftedDiffusion_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK curved_diffusion_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIDiffusion_np${np}_${geometry})
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIApproximatedDiffusion_np${np}_${geometry})
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPILiftedDiffusion_np${np}_${geometry})
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPISmoothLiftedDiffusion_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()

add_executable(RodinConvergenceIsoparametricReactionDiffusion ReactionDiffusion.cpp)
target_link_libraries(RodinConvergenceIsoparametricReactionDiffusion
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_REACTION_DIFFUSION_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_CURVED_REACTION_DIFFUSION_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricReactionDiffusion_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricReactionDiffusion>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
  add_test(NAME RodinConvergenceIsoparametricReactionDiffusion_Approximated_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricReactionDiffusion>
      "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}")
  set_tests_properties(RodinConvergenceIsoparametricReactionDiffusion_${geometry}
    RodinConvergenceIsoparametricReactionDiffusion_Approximated_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricReactionDiffusion_${geometry}
      RodinConvergenceIsoparametricReactionDiffusion_Approximated_${geometry}
      PROPERTIES RESOURCE_LOCK curved_reaction_diffusion_pyramid)
  endif()
endforeach()

foreach(q 1 3)
  foreach(geometry IN LISTS _RODIN_CURVED_REACTION_DIFFUSION_GEOMETRIES)
    set(_native_test RodinConvergenceIsoparametricReactionDiffusion_Q${q}_${geometry})
    add_test(NAME ${_native_test}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricReactionDiffusion>
        "--gtest_filter=AllGeometries/LocalQ${q}Test.*/${geometry}")
    set_tests_properties(${_native_test}
      PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(${_native_test} PROPERTIES RESOURCE_LOCK curved_reaction_diffusion_pyramid)
    endif()
  endforeach()
endforeach()

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_executable(RodinConvergenceIsoparametricPETScReactionDiffusion ReactionDiffusion.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScReactionDiffusion
    PRIVATE RODIN_CURVED_REACTION_DIFFUSION_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScReactionDiffusion
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(q 1 3)
    foreach(geometry IN LISTS _RODIN_CURVED_REACTION_DIFFUSION_GEOMETRIES)
      set(_local_test RodinConvergenceIsoparametricPETScReactionDiffusion_Q${q}_${geometry})
      add_test(NAME ${_local_test}
        COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScReactionDiffusion>
          "--gtest_filter=AllGeometries/LocalQ${q}Test.*/${geometry}")
      set_tests_properties(${_local_test}
        PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
      rodin_suppress_external_mpi_lsan_for_test(${_local_test})
      if (geometry STREQUAL "Pyramid")
        set_tests_properties(${_local_test} PROPERTIES RESOURCE_LOCK curved_reaction_diffusion_pyramid)
      endif()
    endforeach()
  endforeach()
  foreach(geometry IN LISTS _RODIN_CURVED_REACTION_DIFFUSION_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScReactionDiffusion_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScReactionDiffusion>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
    add_test(NAME RodinConvergenceIsoparametricPETScReactionDiffusion_Approximated_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScReactionDiffusion>
        "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}")
    set_tests_properties(RodinConvergenceIsoparametricPETScReactionDiffusion_${geometry}
      RodinConvergenceIsoparametricPETScReactionDiffusion_Approximated_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScReactionDiffusion_${geometry}
        RodinConvergenceIsoparametricPETScReactionDiffusion_Approximated_${geometry}
        PROPERTIES RESOURCE_LOCK curved_reaction_diffusion_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScReactionDiffusion_${geometry})
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScReactionDiffusion_Approximated_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScReactionDiffusion PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(q 1 3)
      foreach(np 1 2 3 4)
        foreach(geometry IN LISTS _RODIN_CURVED_REACTION_DIFFUSION_GEOMETRIES)
          set(_mpi_test RodinConvergenceIsoparametricPETScMPIReactionDiffusion_Q${q}_np${np}_${geometry})
          add_test(NAME ${_mpi_test}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinConvergenceIsoparametricPETScReactionDiffusion>
              "--gtest_filter=AllGeometries/MPIQ${q}Test.*/${geometry}")
          set_tests_properties(${_mpi_test} PROPERTIES
            LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
          rodin_suppress_external_mpi_lsan_for_test(${_mpi_test})
          if (geometry STREQUAL "Pyramid")
            set_tests_properties(${_mpi_test} PROPERTIES RESOURCE_LOCK curved_reaction_diffusion_pyramid)
          endif()
        endforeach()
      endforeach()
    endforeach()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_CURVED_REACTION_DIFFUSION_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPIReactionDiffusion_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScReactionDiffusion>
            "--gtest_filter=AllGeometries/MPITest.*/${geometry}-*MPITest.Approximated*/${geometry}")
        add_test(NAME RodinConvergenceIsoparametricPETScMPIReactionDiffusion_Approximated_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScReactionDiffusion>
            "--gtest_filter=AllGeometries/MPITest.Approximated*/${geometry}")
        set_tests_properties(
          RodinConvergenceIsoparametricPETScMPIReactionDiffusion_np${np}_${geometry}
          RodinConvergenceIsoparametricPETScMPIReactionDiffusion_Approximated_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(
            RodinConvergenceIsoparametricPETScMPIReactionDiffusion_np${np}_${geometry}
            RodinConvergenceIsoparametricPETScMPIReactionDiffusion_Approximated_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK curved_reaction_diffusion_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIReactionDiffusion_np${np}_${geometry})
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIReactionDiffusion_Approximated_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()

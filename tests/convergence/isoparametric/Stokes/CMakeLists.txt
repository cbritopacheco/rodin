add_executable(RodinConvergenceIsoparametricStokes Stokes.cpp)
target_link_libraries(RodinConvergenceIsoparametricStokes
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_STOKES_GEOMETRIES
  Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
add_test(NAME RodinConvergenceIsoparametricStokes_ForwardAccuracy
  COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricStokes>
    "--gtest_filter=Rodin_Convergence_CurvedStokes.NativeCubicWedgePressureForwardAccuracy")
set_tests_properties(RodinConvergenceIsoparametricStokes_ForwardAccuracy
  PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricStokes_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricStokes>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
  add_test(NAME RodinConvergenceIsoparametricStokes_Approximated_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricStokes>
      "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}-*P3Q2*")
  set_tests_properties(RodinConvergenceIsoparametricStokes_${geometry}
    RodinConvergenceIsoparametricStokes_Approximated_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricStokes_${geometry}
      RodinConvergenceIsoparametricStokes_Approximated_${geometry}
      PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
  endif()
endforeach()

# The native cubic-pyramid rate hierarchy exceeded the default 30-minute
# watchdog; the PETSc MPI-one rate case alone took about 25 minutes before
# its independent quadrature study. This one-hour execution policy changes
# neither refinement levels nor numerical acceptance.
set(_RODIN_CURVED_STOKES_Q3_PYRAMID_TIMEOUT 3600)
foreach(q 1 3)
  foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
    set(_native_test RodinConvergenceIsoparametricStokes_Q${q}_${geometry})
    add_test(NAME ${_native_test}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricStokes>
        "--gtest_filter=AllGeometries/LocalQ${q}Test.*/${geometry}")
    set_tests_properties(${_native_test}
      PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(${_native_test} PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
      if (q EQUAL 3)
        set_tests_properties(${_native_test}
          PROPERTIES TIMEOUT ${_RODIN_CURVED_STOKES_Q3_PYRAMID_TIMEOUT})
      endif()
    endif()
  endforeach()
endforeach()

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_executable(RodinConvergenceIsoparametricPETScStokes Stokes.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScStokes
    PRIVATE RODIN_CURVED_STOKES_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScStokes
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(q 1 3)
    foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
      set(_local_test RodinConvergenceIsoparametricPETScStokes_Q${q}_${geometry})
      add_test(NAME ${_local_test}
        COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
          "--gtest_filter=AllGeometries/LocalQ${q}Test.*/${geometry}")
      set_tests_properties(${_local_test}
        PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
      rodin_suppress_external_mpi_lsan_for_test(${_local_test})
      if (geometry STREQUAL "Pyramid")
        set_tests_properties(${_local_test} PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
        if (q EQUAL 3)
          set_tests_properties(${_local_test}
            PROPERTIES TIMEOUT ${_RODIN_CURVED_STOKES_Q3_PYRAMID_TIMEOUT})
        endif()
      endif()
    endforeach()
  endforeach()
  foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScStokes_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
    add_test(NAME RodinConvergenceIsoparametricPETScStokes_Approximated_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
        "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}-*P3Q2*")
    set_tests_properties(RodinConvergenceIsoparametricPETScStokes_${geometry}
      RodinConvergenceIsoparametricPETScStokes_Approximated_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScStokes_${geometry}
        RodinConvergenceIsoparametricPETScStokes_Approximated_${geometry}
        PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScStokes_${geometry})
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScStokes_Approximated_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScStokes PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(q 1 3)
      foreach(np 1 2 3 4)
        foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
          set(_mpi_test RodinConvergenceIsoparametricPETScMPIStokes_Q${q}_np${np}_${geometry})
          add_test(NAME ${_mpi_test}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
              "--gtest_filter=AllGeometries/MPIQ${q}Test.*/${geometry}")
          set_tests_properties(${_mpi_test} PROPERTIES
            LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
          rodin_suppress_external_mpi_lsan_for_test(${_mpi_test})
          if (geometry STREQUAL "Pyramid")
            set_tests_properties(${_mpi_test} PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
            if (q EQUAL 3)
              set_tests_properties(${_mpi_test}
                PROPERTIES TIMEOUT ${_RODIN_CURVED_STOKES_Q3_PYRAMID_TIMEOUT})
            endif()
          endif()
        endforeach()
      endforeach()
    endforeach()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPIStokes_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
            "--gtest_filter=AllGeometries/MPITest.*/${geometry}-*MPITest.Approximated*/${geometry}")
        add_test(NAME RodinConvergenceIsoparametricPETScMPIStokes_Approximated_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
            "--gtest_filter=AllGeometries/MPITest.Approximated*/${geometry}-*P3Q2*")
        set_tests_properties(RodinConvergenceIsoparametricPETScMPIStokes_np${np}_${geometry}
          RodinConvergenceIsoparametricPETScMPIStokes_Approximated_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(RodinConvergenceIsoparametricPETScMPIStokes_np${np}_${geometry}
            RodinConvergenceIsoparametricPETScMPIStokes_Approximated_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIStokes_np${np}_${geometry})
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIStokes_Approximated_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()

# Isolate each complete higher-degree hierarchy from its independent controls.
foreach(group Rates SeparatedRates Controls)
  set(_higher_case "ApproximatedP3Q2*")
  set(_higher_exclusion "")
  if (group STREQUAL "Rates")
    set(_higher_case "ApproximatedP3Q2VelocityPressureRates")
  elseif (group STREQUAL "SeparatedRates")
    set(_higher_case "ApproximatedP3Q2SeparatedPressureRates")
  else()
    set(_higher_exclusion "-*ApproximatedP3Q2*Rates*")
  endif()
  foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
    set(_higher_test RodinConvergenceIsoparametricStokes_P3Q2_${group}_${geometry})
    add_test(NAME ${_higher_test}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricStokes>
        "--gtest_filter=AllGeometries/LocalTest.${_higher_case}/${geometry}${_higher_exclusion}")
    set_tests_properties(${_higher_test} PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(${_higher_test} PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
    endif()
    if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
      set(_higher_test RodinConvergenceIsoparametricPETScStokes_P3Q2_${group}_${geometry})
      add_test(NAME ${_higher_test}
        COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
          "--gtest_filter=AllGeometries/LocalTest.${_higher_case}/${geometry}${_higher_exclusion}")
      set_tests_properties(${_higher_test} PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
      rodin_suppress_external_mpi_lsan_for_test(${_higher_test})
      if (geometry STREQUAL "Pyramid")
        set_tests_properties(${_higher_test} PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
      endif()
      if (RODIN_USE_MPI)
        foreach(np 1 2 3 4)
          set(_higher_test RodinConvergenceIsoparametricPETScMPIStokes_P3Q2_np${np}_${group}_${geometry})
          add_test(NAME ${_higher_test}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
              "--gtest_filter=AllGeometries/MPITest.${_higher_case}/${geometry}${_higher_exclusion}")
          set_tests_properties(${_higher_test} PROPERTIES
            LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
          rodin_suppress_external_mpi_lsan_for_test(${_higher_test})
          if (geometry STREQUAL "Pyramid")
            set_tests_properties(${_higher_test} PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
          endif()
        endforeach()
      endif()
    endif()
  endforeach()
endforeach()

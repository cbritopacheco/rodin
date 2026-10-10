add_executable(RodinAdaptationDeformationMapTest DeformationMapTest.cpp)
target_link_libraries(RodinAdaptationDeformationMapTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationDeformationMapTest PROPERTIES LABELS unit)

add_executable(RodinAdaptationCellDeformationTest CellDeformationTest.cpp)
target_link_libraries(RodinAdaptationCellDeformationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinAdaptationCellDeformationTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationSWIFTHingeStateTest
  SWIFTHingeStateTest.cpp)
target_link_libraries(RodinAdaptationSWIFTHingeStateTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTHingeStateTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationSWIFTLossTest SWIFTLossTest.cpp)
target_link_libraries(RodinAdaptationSWIFTLossTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTLossTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationSWIFTSolverTest SWIFTSolverTest.cpp)
add_executable(RodinAdaptationSWIFTAdaptTest SWIFTAdaptTest.cpp)
target_link_libraries(RodinAdaptationSWIFTAdaptTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTAdaptTest PROPERTIES LABELS unit)

if(TARGET Rodin::MMG)
  add_executable(RodinAdaptationSWIFTExampleParametersTest SWIFTExampleParametersTest.cpp)
  target_link_libraries(RodinAdaptationSWIFTExampleParametersTest
    PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation Rodin::MMG)
  gtest_discover_tests(RodinAdaptationSWIFTExampleParametersTest PROPERTIES LABELS unit)
endif()

add_executable(RodinAdaptationSWIFTRegularityMetricTest SWIFTRegularityMetricTest.cpp)
add_executable(RodinAdaptationSWIFTDirectionalNewtonTest SWIFTDirectionalNewtonTest.cpp)
target_link_libraries(RodinAdaptationSWIFTDirectionalNewtonTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTDirectionalNewtonTest PROPERTIES LABELS unit)
target_link_libraries(RodinAdaptationSWIFTRegularityMetricTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTRegularityMetricTest PROPERTIES LABELS unit)
target_link_libraries(RodinAdaptationSWIFTSolverTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTSolverTest
  PROPERTIES
    LABELS unit
)

if(TARGET SWIFT_ReconstructionP1)
  add_executable(RodinAdaptationSWIFTOutputTest SWIFTOutputTest.cpp)
  target_link_libraries(RodinAdaptationSWIFTOutputTest
    PRIVATE GTest::gtest GTest::gtest_main Rodin::IO)
  foreach(degree IN ITEMS 1 2 3)
    add_dependencies(RodinAdaptationSWIFTOutputTest SWIFT_ReconstructionP${degree})
    foreach(dimension IN ITEMS 2 3)
      set(output_dir "${CMAKE_CURRENT_BINARY_DIR}/swift-output-p${degree}-d${dimension}")
      file(MAKE_DIRECTORY "${output_dir}")
      set(fixture "SWIFT.Output.P${degree}.D${dimension}")
      add_test(NAME ${fixture}.Generate
        COMMAND $<TARGET_FILE:SWIFT_ReconstructionP${degree}>
          --n=4 --dimension=${dimension} --convergence-iterations-outer=1
          --output=${output_dir}/reconstruction)
      set_tests_properties(${fixture}.Generate PROPERTIES
        LABELS unit FIXTURES_SETUP ${fixture} TIMEOUT 120
        WORKING_DIRECTORY "${output_dir}"
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
      add_test(NAME ${fixture}.Verify COMMAND RodinAdaptationSWIFTOutputTest)
      set_tests_properties(${fixture}.Verify PROPERTIES
        LABELS unit FIXTURES_REQUIRED ${fixture}
        ENVIRONMENT "SWIFT_OUTPUT=${output_dir}/reconstruction;SWIFT_DIMENSION=${dimension}")
    endforeach()
  endforeach()
endif()

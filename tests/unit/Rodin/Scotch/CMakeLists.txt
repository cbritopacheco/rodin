add_executable(RodinScotchPartitionerTest PartitionerTest.cpp)
target_link_libraries(RodinScotchPartitionerTest
  PRIVATE GTest::gtest_main Rodin::Scotch)
gtest_discover_tests(RodinScotchPartitionerTest PROPERTIES LABELS unit)

# Link-time wrappers intercept backend calls in static objects.
if(NOT BUILD_SHARED_LIBS AND CMAKE_SYSTEM_NAME STREQUAL "Linux"
    AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_executable(RodinScotchFailureTest FailureTest.cpp)
  target_link_libraries(RodinScotchFailureTest PRIVATE GTest::gtest_main Rodin::Scotch)
  foreach(operation IN ITEMS SCOTCH_graphInit SCOTCH_graphBuild SCOTCH_graphPart
      SCOTCH_graphExit SCOTCH_stratInit SCOTCH_stratExit SCOTCH_stratGraphMap
      SCOTCH_stratSave tmpfile fclose fseek fread ferror)
    target_link_options(RodinScotchFailureTest PRIVATE "LINKER:--wrap=${operation}")
  endforeach()
  gtest_discover_tests(RodinScotchFailureTest PROPERTIES LABELS unit)
endif()

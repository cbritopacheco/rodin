add_executable(RodinPairTest PairTest.cpp)
target_link_libraries(RodinPairTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinPairTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinTupleTest TupleTest.cpp)
target_link_libraries(RodinTupleTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinTupleTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinArrayTest ArrayTest.cpp)
target_link_libraries(RodinArrayTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinArrayTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinCopyableTest CopyableTest.cpp)
target_link_libraries(RodinCopyableTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinCopyableTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMoveableTest MoveableTest.cpp)
target_link_libraries(RodinMoveableTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMoveableTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinFlatSetTest FlatSetTest.cpp)
target_link_libraries(RodinFlatSetTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinFlatSetTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinTypesTest TypesTest.cpp)
target_link_libraries(RodinTypesTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinTypesTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinCastTest CastTest.cpp)
target_link_libraries(RodinCastTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinCastTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinUtilityTest UtilityTest.cpp)
target_link_libraries(RodinUtilityTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinUtilityTest
  PROPERTIES
    LABELS unit
)

add_subdirectory(Alert)
add_subdirectory(IO)
add_subdirectory(QF)
add_subdirectory(Context)
add_subdirectory(Threads)
add_subdirectory(Geometry)
add_subdirectory(Variational)
add_subdirectory(Serialization)
add_subdirectory(Math)
add_subdirectory(Test)
add_subdirectory(FormLanguage)
add_subdirectory(Models)
add_subdirectory(Solver)
add_subdirectory(Assembly)
add_subdirectory(MMG)

if (RODIN_USE_MPI)
  add_subdirectory(MPI)
endif()

if (RODIN_USE_PETSC)
  add_subdirectory(PETSc)
endif()

add_executable(RodinManufacturedNewtonSolverDenseMatrix DenseMatrix.cpp)
target_link_libraries(RodinManufacturedNewtonSolverDenseMatrix
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Alert
  Rodin::Rodin
  Rodin::Solver
)
rodin_discover_tests(RodinManufacturedNewtonSolverDenseMatrix
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedNewtonSolverSparseMatrix SparseMatrix.cpp)
target_link_libraries(RodinManufacturedNewtonSolverSparseMatrix
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Alert
  Rodin::Rodin
  Rodin::Solver
)
rodin_discover_tests(RodinManufacturedNewtonSolverSparseMatrix
  PROPERTIES
    LABELS manufactured
)

foreach(algorithm CG BiCGSTAB GMRES DGMRES LDLT HouseholderQR PartialPivLU
    SimplicialLDLT SimplicialLLT SparseQR)
  set(target Rodin${algorithm}OwnershipTest)
  add_executable(${target} OwnershipTest.cpp)
  target_link_libraries(${target} PRIVATE
    GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::Solver)
  target_compile_definitions(${target} PRIVATE
    RODIN_TEST_SOLVER=${algorithm}
    RODIN_TEST_SOLVER_HEADER="Rodin/Solver/${algorithm}.h")
  if(algorithm MATCHES "^(CG|GMRES|DGMRES|LDLT|HouseholderQR|PartialPivLU)$")
    target_compile_definitions(${target} PRIVATE RODIN_TEST_DENSE)
  endif()
  if(algorithm MATCHES "^(CG|BiCGSTAB|GMRES|DGMRES|SimplicialLDLT|SimplicialLLT|SparseQR)$")
    target_compile_definitions(${target} PRIVATE RODIN_TEST_SPARSE)
  endif()
  gtest_discover_tests(${target} TEST_PREFIX "${algorithm}."
    DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
endforeach()

add_executable(RodinIterativeSolverTest IterativeSolverTest.cpp)
target_link_libraries(RodinIterativeSolverTest
  GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::Solver)
gtest_discover_tests(RodinIterativeSolverTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS unit)

add_executable(RodinNewtonSolverTest NewtonSolverTest.cpp)
target_link_libraries(RodinNewtonSolverTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Alert
  Rodin::Rodin
  Rodin::Solver
)
gtest_discover_tests(RodinNewtonSolverTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinInitialGuessTest InitialGuessTest.cpp)
target_link_libraries(RodinInitialGuessTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinInitialGuessTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinSolverConceptTest SolverConceptTest.cpp)
target_link_libraries(RodinSolverConceptTest
  GTest::gtest
  GTest::gtest_main
  Rodin::Alert
  Rodin::Rodin
  Rodin::Solver
)
gtest_discover_tests(RodinSolverConceptTest
   PROPERTIES
     LABELS unit
 )

if (RODIN_USE_PARU)
  add_executable(RodinParUContractTest ParUContractTest.cpp)
  target_link_libraries(RodinParUContractTest PRIVATE
    GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::Solver)
  gtest_discover_tests(RodinParUContractTest
    DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

  add_executable(RodinParUTest ParUTest.cpp)
  target_link_libraries(RodinParUTest
    GTest::gtest
    GTest::gtest_main
    Rodin::Geometry
    Rodin::Variational
    Rodin::Solver
  )
  gtest_discover_tests(RodinParUTest
    PROPERTIES
      LABELS unit
  )
endif()

if (RODIN_USE_MUMPS)
  add_executable(RodinMUMPSTest MUMPSTest.cpp)
  target_link_libraries(RodinMUMPSTest
    GTest::gtest
    GTest::gtest_main
    Rodin::Geometry
    Rodin::Variational
    Rodin::Solver
  )
  gtest_discover_tests(RodinMUMPSTest
    PROPERTIES
      LABELS unit
  )
endif()

add_executable(RodinDenseSolverTest DenseSolverTest.cpp)
target_link_libraries(RodinDenseSolverTest
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinDenseSolverTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinInfoTest InfoTest.cpp)
target_link_libraries(RodinInfoTest
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinInfoTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinSparseLUTest SparseLUTest.cpp)
target_link_libraries(RodinSparseLUTest
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinSparseLUTest
  PROPERTIES
    LABELS unit
)

if (RODIN_USE_UMFPACK)
  add_executable(RodinUMFPackTest UMFPackTest.cpp)
  target_link_libraries(RodinUMFPackTest
    GTest::gtest
    GTest::gtest_main
    Rodin::Geometry
    Rodin::Variational
    Rodin::Solver
  )
  gtest_discover_tests(RodinUMFPackTest
    PROPERTIES
      LABELS unit
  )
endif()

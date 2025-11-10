set(RodinSolver_HEADERS
  Solver.h
  CG.h
  )

set(RodinSolver_SRCS
  CG.cpp)

# ---- Set targets -----------------------------------------------------------
# add_library(RodinSolver ${RodinSolver_SRCS} ${RodinSolver_HEADERS})
add_library(RodinSolver INTERFACE)
add_library(Rodin::Solver ALIAS RodinSolver)
set_target_properties(RodinSolver PROPERTIES EXPORT_NAME Solver)

# ---- OpenMP ----------------------------------------------------------------
if (RODIN_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(RodinSolver INTERFACE OpenMP::OpenMP_CXX)
endif()

# ---- SuiteSparse -----------------------------------------------------------
if (RODIN_USE_UMFPACK)
  find_package(SuiteSparse REQUIRED COMPONENTS UMFPACK)
  target_link_libraries(RodinSolver INTERFACE SuiteSparse::UMFPACK)
endif()

if (RODIN_USE_SPQR)
  find_package(SuiteSparse REQUIRED COMPONENTS SPQR)
  target_link_libraries(RodinSolver INTERFACE SuiteSparse::SPQR)
endif()

if (RODIN_USE_CHOLMOD)
  find_package(SuiteSparse REQUIRED COMPONENTS CHOLMOD)
  target_link_libraries(RodinSolver INTERFACE SuiteSparse::CHOLMOD)
endif()

if (RODIN_USE_APPLE_ACCELERATE)
  target_sources(RodinSolver INTERFACE AppleAccelerate.h)
endif()

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinSolver
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinSolver
  INTERFACE
  Eigen3::Eigen)

install(TARGETS RodinSolver
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Solver
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Solver
)

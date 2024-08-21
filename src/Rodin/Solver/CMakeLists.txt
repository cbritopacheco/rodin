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

# ---- OpenMP ----------------------------------------------------------------
if (RODIN_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(RodinSolver INTERFACE OpenMP::OpenMP_CXX)
endif()

# ---- SuiteSparse -----------------------------------------------------------
if (RODIN_USE_UMFPACK)
  find_package(SuiteSparse REQUIRED COMPONENTS UMFPACK)
  target_sources(RodinSolver INTERFACE UMFPack.h)
  target_link_libraries(RodinSolver INTERFACE SuiteSparse::UMFPACK)
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


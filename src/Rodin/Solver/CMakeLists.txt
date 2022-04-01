set(RodinSolver_HEADERS
  Solver.h
  CG.h
  )

set(RodinSolver_SRCS
  CG.cpp
  )

# ---- Set targets -----------------------------------------------------------
add_library(RodinSolver ${RodinSolver_SRCS} ${RodinSolver_HEADERS})
add_library(Rodin::Solver ALIAS RodinSolver)

# ---- SuiteSparse -----------------------------------------------------------
if (RODIN_USE_SUITESPARSE)
  target_compile_definitions(RodinSolver PUBLIC RODIN_SOLVER_SUITESPARSE)

  target_sources(RodinSolver PRIVATE
    UMFPack.h
    UMFPack.cpp
    )
endif()

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinSolver
  PUBLIC
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  PRIVATE
  ${MFEM_INCLUDE_DIRS}
  )

target_link_libraries(RodinSolver
  PUBLIC
  RodinVariational
  PRIVATE
  mfem
  )


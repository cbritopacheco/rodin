set(RodinSolver_HEADERS
  Solver.h
  PCG.h
  CG.h
  )

set(RodinSolver_SRCS
  PCG.cpp
  CG.cpp
  )

add_library(RodinSolver ${RodinSolver_SRCS} ${RodinSolver_HEADERS})
add_library(Rodin::Solver ALIAS RodinSolver)


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


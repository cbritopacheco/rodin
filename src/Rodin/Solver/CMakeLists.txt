#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinSolver_HEADERS
  Solver.h
  PCG.h
  )

set(RodinSolver_SRCS
  PCG.cpp
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


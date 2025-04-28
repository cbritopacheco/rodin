#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinPETSc_HEADERS
  ForwardDecls.h
  PETSc.h
  Assembly/Sequential.h
  Assembly/MPI.h
  )

set(RodinPETSc_SRCS
  PETSc.cpp
  )

add_library(RodinPETSc ${RodinPETSc_HEADERS} ${RodinPETSc_SRCS})
add_library(Rodin::PETSc ALIAS RodinPETSc)

target_include_directories(RodinPETSc
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

find_package(PETSc REQUIRED)

target_link_libraries(RodinPETSc
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  Rodin::MPI
  PETSc::PETSc)


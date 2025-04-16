#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinPETSc_HEADERS
  )

set(RodinPETSc_SRCS
  )

add_library(RodinPETSc
  ${RodinPETSc_HEADERS} ${RodinPETSc_SRCS})
add_library(Rodin::PETSc ALIAS RodinPETSc)

find_package(PETSc REQUIRED)

target_link_libraries(RodinPETSc
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  PETSc::PETSc)


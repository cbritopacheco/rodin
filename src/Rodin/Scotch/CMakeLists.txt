#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinScotch_HEADERS
  MeshPartitioner.h
  )

set(RodinScotch_SRCS
  MeshPartitioner.cpp
  )

add_library(RodinScotch
  ${RodinScotch_HEADERS} ${RodinScotch_SRCS})
add_library(Rodin::Scotch ALIAS RodinScotch)

find_package(Scotch REQUIRED)
target_link_libraries(RodinScotch
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  Scotch::Scotch)

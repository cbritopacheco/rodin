#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinExternalScotch_HEADERS
  MeshPartitioner.h
  )

set(RodinExternalScotch_SRCS
  MeshPartitioner.cpp
  )

add_library(RodinExternalScotch
  ${RodinExternalScotch_HEADERS} ${RodinExternalScotch_SRCS})
add_library(Rodin::External::Scotch ALIAS RodinExternalScotch)

find_package(scotch REQUIRED)
target_link_libraries(RodinExternalScotch
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  scotch)

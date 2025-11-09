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
set_target_properties(RodinScotch PROPERTIES EXPORT_NAME Scotch)

target_include_directories(RodinScotch
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

find_package(Scotch REQUIRED)
target_link_libraries(RodinScotch
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  Scotch::Scotch)

install(TARGETS RodinScotch
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Scotch
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Scotch
)

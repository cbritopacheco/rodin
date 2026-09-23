#[[
         Copyright Carlos BRITO PACHECO 2026.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinLocation_HEADERS
  ForwardDecls.h
  AABB.h
)

add_library(RodinLocation INTERFACE)
add_library(Rodin::Location ALIAS RodinLocation)
set_target_properties(RodinLocation PROPERTIES EXPORT_NAME Location)

target_include_directories(RodinLocation
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinLocation
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

target_link_libraries(RodinLocation
  INTERFACE
    Rodin::Geometry
)

foreach(_h IN LISTS RodinLocation_HEADERS)
  target_sources(RodinLocation INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

install(
  TARGETS RodinLocation
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Location
)

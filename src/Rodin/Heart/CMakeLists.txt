#[[
         Copyright Carlos BRITO PACHECO 2021 - 2026.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinHeart_HEADERS
  ForwardDecls.h
  CCMLC2014.h
  CCMLC2014/CCMLC2014.h
  CCMLC2014/PassiveLaw.h
)

add_library(RodinHeart INTERFACE)
add_library(Rodin::Heart ALIAS RodinHeart)
set_target_properties(RodinHeart PROPERTIES EXPORT_NAME Heart)

target_include_directories(RodinHeart
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinHeart
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

target_link_libraries(RodinHeart
  INTERFACE
    Rodin::Solver
)

foreach(_h IN LISTS RodinHeart_HEADERS)
  target_sources(RodinHeart INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

install(
  TARGETS RodinHeart
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Heart
)

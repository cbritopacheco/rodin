#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinAdvection_HEADERS
  Lagrangian.h
)

# Header-only / interface library
add_library(RodinAdvection INTERFACE)
add_library(Rodin::Advection ALIAS RodinAdvection)
set_target_properties(RodinAdvection PROPERTIES EXPORT_NAME Advection)

# ---- Includes / defs -------------------------------------------------------
target_include_directories(RodinAdvection
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinAdvection
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinAdvection
  INTERFACE
    Rodin::Geometry
    Rodin::Variational
)

# Expose headers for IDEs in build tree
foreach(_h IN LISTS RodinAdvection_HEADERS)
  target_sources(RodinAdvection INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinAdvection
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin//Advection
)

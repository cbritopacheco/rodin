#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinModelsAdvection_HEADERS
  Lagrangian.h
)

# Header-only / interface library
add_library(RodinModelsAdvection INTERFACE)
add_library(Rodin::Models::Advection ALIAS RodinModelsAdvection)
set_target_properties(RodinModelsAdvection PROPERTIES EXPORT_NAME ModelsAdvection)

# ---- Includes / defs -------------------------------------------------------
target_include_directories(RodinModelsAdvection
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinModelsAdvection
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinModelsAdvection
  INTERFACE
    Rodin::Geometry
    Rodin::Variational
)

# Expose headers for IDEs in build tree
foreach(_h IN LISTS RodinModelsAdvection_HEADERS)
  target_sources(RodinModelsAdvection INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinModelsAdvection
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Models/Advection
)

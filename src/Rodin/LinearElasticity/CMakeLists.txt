#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinLinearElasticity_HEADERS
  FMM.h
)

# Header-only / interface library
add_library(RodinLinearElasticity INTERFACE)
add_library(Rodin::LinearElasticity ALIAS RodinLinearElasticity)
set_target_properties(RodinLinearElasticity PROPERTIES EXPORT_NAME LinearElasticity)

# ---- Includes / defs -------------------------------------------------------
target_include_directories(RodinLinearElasticity
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinLinearElasticity
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinLinearElasticity
  INTERFACE
    Rodin::Geometry
    Rodin::Variational
)

# Expose headers for IDEs in build tree
foreach(_h IN LISTS RodinLinearElasticity_HEADERS)
  target_sources(RodinLinearElasticity INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinLinearElasticity
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/LinearElasticity
)


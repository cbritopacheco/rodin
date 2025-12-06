#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinModelsEikonal_HEADERS
  FMM.h
)

# Header-only / interface library
add_library(RodinModelsEikonal INTERFACE)
add_library(Rodin::Models::Eikonal ALIAS RodinModelsEikonal)
set_target_properties(RodinModelsEikonal PROPERTIES EXPORT_NAME ModelsEikonal)

# ---- Includes / defs -------------------------------------------------------
target_include_directories(RodinModelsEikonal
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinModelsEikonal
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinModelsEikonal
  INTERFACE
    Rodin::Geometry
    Rodin::Variational
)

# Expose headers for IDEs in build tree
foreach(_h IN LISTS RodinModelsEikonal_HEADERS)
  target_sources(RodinModelsEikonal INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinModelsEikonal
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Models/Eikonal
)

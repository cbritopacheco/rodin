#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinUtility_HEADERS
  Cast.h
  IsSpecialization.h
  OptionalReference.h
  Overloaded.h
)

# Header-only utility helpers
add_library(RodinUtility INTERFACE)
add_library(Rodin::Utility ALIAS RodinUtility)
set_target_properties(RodinUtility PROPERTIES EXPORT_NAME Utility)

target_include_directories(RodinUtility
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinUtility
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# Expose headers for IDEs when building from source
foreach(_h IN LISTS RodinUtility_HEADERS)
  target_sources(RodinUtility INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

install(TARGETS RodinUtility
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Utility
)

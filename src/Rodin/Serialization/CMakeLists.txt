#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinSerialization_HEADERS
  EigenMatrix.h
)

add_library(RodinSerialization INTERFACE)
add_library(Rodin::Serialization ALIAS RodinSerialization)
set_target_properties(RodinSerialization PROPERTIES EXPORT_NAME Serialization)

# Inherit global Rodin interface (includes, defs)
target_include_directories(RodinSerialization
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinSerialization
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# Dependencies
target_link_libraries(RodinSerialization
  INTERFACE
    Boost::serialization
)

# Expose headers for IDE / build tree
foreach(_h IN LISTS RodinSerialization_HEADERS)
  target_sources(RodinSerialization INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# Install target (INTERFACE lib: no RUNTIME/LIBRARY/ARCHIVE destinations)
install(
  TARGETS RodinSerialization
  EXPORT RodinTargets
)

# Install headers on disk
rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Serialization
)

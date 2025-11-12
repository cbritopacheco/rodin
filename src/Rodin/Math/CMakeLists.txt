#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

# Header interface for Math; sources (if any) are built elsewhere or unused here
set(RodinMath_HEADERS
  Constants.h
  Common.h
  Vector.h
  Matrix.h
)

add_library(RodinMath INTERFACE)
add_library(Rodin::Math ALIAS RodinMath)
set_target_properties(RodinMath PROPERTIES EXPORT_NAME Math)

# Inherit the umbrella include interface (root include dir only)
target_include_directories(RodinMath
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

# Inherit Rodin's global compile definitions (e.g. libc++ compatibility)
target_compile_definitions(RodinMath
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# Header-only dependency
target_link_libraries(RodinMath
  INTERFACE
    Eigen3::Eigen
)

# Expose headers for IDEs when building from source
foreach(_h IN LISTS RodinMath_HEADERS)
  target_sources(RodinMath INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# Install target — do NOT set INCLUDES DESTINATION to a subfolder
install(
  TARGETS RodinMath
  EXPORT RodinTargets
)

# Install headers to include/Rodin/Math on disk
rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Math
)

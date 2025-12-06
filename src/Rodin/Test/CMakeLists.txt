#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinTest_HEADERS
  Random.h
  Utility.h
  Random/RandomInteger.h
  Utility/ExpectedResultTable.h
)

set(RodinTest_SRCS
  Utility/ExpectedResultTable.hpp
)

add_library(RodinTest INTERFACE)
add_library(Rodin::Test ALIAS RodinTest)
set_target_properties(RodinTest PROPERTIES EXPORT_NAME Test)

# Expose headers and header-impls to consumers/IDEs
foreach(_file IN LISTS RodinTest_SRCS RodinTest_HEADERS)
  target_sources(RodinTest INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_file}>
  )
endforeach()

target_include_directories(RodinTest
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinTest
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# Install INTERFACE target (no runtime/library/archive needed)
install(
  TARGETS RodinTest
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Test
)

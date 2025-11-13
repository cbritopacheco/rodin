#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinContext_HEADERS
  Local.h
)

set(RodinContext_SRCS
  Local.cpp
)

add_library(RodinContext
  ${RodinContext_SRCS}
  ${RodinContext_HEADERS}
)
add_library(Rodin::Context ALIAS RodinContext)
set_target_properties(RodinContext PROPERTIES EXPORT_NAME Context)

# Reuse Rodin's include dirs (root include path) and propagate them
target_include_directories(RodinContext
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

# If you want Rodin’s global compile definitions (like the libc++ one) too:
target_compile_definitions(RodinContext
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

install(
  TARGETS RodinContext
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Context
)

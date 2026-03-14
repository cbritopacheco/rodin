#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinFormLanguage_HEADERS
  Base.h
  List.h
  ForwardDecls.h
)

set(RodinFormLanguage_SRCS
  Base.cpp
)

add_library(RodinFormLanguage
  ${RodinFormLanguage_SRCS}
  ${RodinFormLanguage_HEADERS}
)
add_library(Rodin::FormLanguage ALIAS RodinFormLanguage)
set_target_properties(RodinFormLanguage PROPERTIES EXPORT_NAME FormLanguage)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinFormLanguage
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinFormLanguage
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinFormLanguage
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/FormLanguage
)

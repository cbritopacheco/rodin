#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinAlert_HEADERS
  Info.h
  Message.h
  Warning.h
  Success.h
  Exception.h
  Notation.h
)

set(RodinAlert_SRCS
  Info.cpp
  Warning.cpp
  Success.cpp
  Notation.cpp
  Exception.cpp
)

add_library(RodinAlert ${RodinAlert_SRCS} ${RodinAlert_HEADERS})
add_library(Rodin::Alert ALIAS RodinAlert)
set_target_properties(RodinAlert PROPERTIES EXPORT_NAME Alert)

# ---- Link targets ----------------------------------------------------------
# Publicly depend on the umbrella Rodin target
target_link_libraries(RodinAlert PUBLIC Rodin::Rodin)

# Termcolor during build + root include dir after install
target_include_directories(RodinAlert
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third-party/termcolor/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinAlert
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

# Install Rodin alert headers under include/Rodin/Alert
rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Alert
)

# Install termcolor so that #include <termcolor/termcolor.hpp> works
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/third-party/termcolor/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

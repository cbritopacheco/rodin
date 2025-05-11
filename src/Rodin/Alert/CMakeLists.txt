#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinAlert_HEADERS
  Info.h
  Message.h
  Warning.h
  Success.h
  Exception.h)

set(RodinAlert_SRCS
  Info.cpp
  Warning.cpp
  Success.cpp
  Notation.cpp
  Exception.cpp)

add_library(RodinAlert ${RodinAlert_SRCS} ${RodinAlert_HEADERS})
add_library(Rodin::Alert ALIAS RodinAlert)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinAlert
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third-party/termcolor/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Alert/termcolor/include>
  )

target_link_libraries(RodinAlert PUBLIC Rodin::Rodin)

install(TARGETS RodinAlert
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Alert
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Alert
)

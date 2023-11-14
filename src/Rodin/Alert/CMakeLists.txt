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
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  ${PROJECT_SOURCE_DIR}/third-party/termcolor/include)


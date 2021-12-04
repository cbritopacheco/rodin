#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinUtility_HEADERS
  Overloaded.h
  IsSpecialization.h
  )

# target_sources() does not resolve relative paths before 3.13
if (${CMAKE_VERSION} VERSION_LESS 3.13)
  list(
    TRANSFORM RodinUtility_HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
endif()

add_library(RodinUtility INTERFACE)
add_library(Rodin::Utility ALIAS RodinUtility)
target_sources(RodinUtility INTERFACE ${RodinUtility_HEADERS})

target_include_directories(RodinUtility
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  )


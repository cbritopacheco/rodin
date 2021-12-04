set(RodinCore_HEADERS
  Common.h
  )

set(RodinCore_SRCS
  Common.hpp
  )

# target_sources() does not resolve relative paths before 3.13
if (${CMAKE_VERSION} VERSION_LESS 3.13)
  list(
    TRANSFORM RodinCore_SRCS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
  list(
    TRANSFORM RodinCore_HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
endif()

add_library(RodinCore INTERFACE)
add_library(Rodin::Core ALIAS RodinCore)
target_sources(RodinCore INTERFACE ${RodinCore_SRCS} ${RodinCore_HEADERS})

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinCore
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  )
target_link_libraries(RodinCore
  INTERFACE
  Eigen3::Eigen
  )



set(RodinMath_HEADERS
  Common.h)

set(RodinMath_SRCS )

# target_sources() does not resolve relative paths before 3.13
if (${CMAKE_VERSION} VERSION_LESS 3.13)
  list(
    TRANSFORM RodinMath_SRCS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
  list(
    TRANSFORM RodinMath_HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
endif()

add_library(RodinMath INTERFACE)
add_library(Rodin::Math ALIAS RodinMath)
target_sources(RodinMath INTERFACE ${RodinMath_SRCS} ${RodinMath_HEADERS})

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinMath
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  )


set(RodinTest_HEADERS
  Random.h
  Utility.h
  Random/RandomInteger.h
  Utility/ExpectedResultTable.h
  )

set(RodinTest_SRCS
  Utility/ExpectedResultTable.hpp
  )

# target_sources() does not resolve relative paths before 3.13
if (${CMAKE_VERSION} VERSION_LESS 3.13)
  list(
    TRANSFORM RodinTest_SRCS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
  list(
    TRANSFORM RodinTest_HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
endif()

add_library(RodinTest INTERFACE)
add_library(Rodin::Test ALIAS RodinTest)
target_sources(RodinTest INTERFACE ${RodinTest_SRCS} ${RodinTest_HEADERS})

target_include_directories(RodinTest
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)



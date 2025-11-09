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

foreach(_file IN LISTS RodinTest_SRCS RodinTest_HEADERS)
  target_sources(RodinTest INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_file}>
  )
endforeach()

target_include_directories(RodinTest
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

install(TARGETS RodinTest
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Test
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Test
)


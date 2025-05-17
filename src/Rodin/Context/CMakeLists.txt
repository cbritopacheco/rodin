set(RodinContext_HEADERS Sequential.h)

set(RodinContext_SRCS Sequential.cpp)

# add_library(RodinContext ${RodinContext_SRCS} ${RodinContext_HEADERS})
add_library(RodinContext INTERFACE)
add_library(Rodin::Context ALIAS RodinContext)

target_include_directories(RodinContext
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

install(TARGETS RodinContext
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Context
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Context
)

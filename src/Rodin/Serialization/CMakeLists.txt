set(RodinSerialization_HEADERS EigenMatrix.h)

set(RodinSerialization_SRCS )
add_library(RodinSerialization INTERFACE)
add_library(Rodin::Serialization ALIAS RodinSerialization)

target_link_libraries(RodinSerialization INTERFACE Boost::serialization)

install(TARGETS RodinSerialization
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Serialization
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Serialization
)

set(RodinUtility_HEADERS
  Cast.h
  IsSpecialization.h
  OptionalReference.h
  Overloaded.h
  )

set(RodinUtility_SRCS Cast.cpp)

# add_library(RodinUtility ${RodinUtility_SRCS} ${RodinUtility_HEADERS})
add_library(RodinUtility INTERFACE)
add_library(Rodin::Utility ALIAS RodinUtility)
set_target_properties(RodinUtility PROPERTIES EXPORT_NAME Utility)

target_include_directories(RodinUtility
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

install(TARGETS RodinUtility
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Utility
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Utility
)

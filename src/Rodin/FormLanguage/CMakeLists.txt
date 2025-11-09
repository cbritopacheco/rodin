set(RodinFormLanguage_HEADERS
  Base.h
  List.h
  ForwardDecls.h)

set(RodinFormLanguage_SRCS Base.cpp)

add_library(RodinFormLanguage ${RodinFormLanguage_SRCS} ${RodinFormLanguage_HEADERS})
add_library(Rodin::FormLanguage ALIAS RodinFormLanguage)
set_target_properties(RodinFormLanguage PROPERTIES EXPORT_NAME FormLanguage)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinFormLanguage
  PUBLIC $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinFormLanguage PUBLIC Rodin::Threads)

install(TARGETS RodinFormLanguage
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/FormLanguage
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/FormLanguage
)

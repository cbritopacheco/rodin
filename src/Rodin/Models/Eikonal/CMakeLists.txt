set(RodinModelsEikonal_HEADERS
  FMM.h
  )

set(RodinModelsEikonal_SRCS
  )

add_library(RodinModelsEikonal INTERFACE)
add_library(Rodin::Models::Eikonal ALIAS RodinModelsEikonal)
set_target_properties(RodinModelsEikonal PROPERTIES EXPORT_NAME Models::Eikonal)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinModelsEikonal
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinModelsEikonal
  INTERFACE
  Rodin::Geometry
  Rodin::Variational)

install(TARGETS RodinModelsEikonal
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Models/Eikonal
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Models/Eikonal
)


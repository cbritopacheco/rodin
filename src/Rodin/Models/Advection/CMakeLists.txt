set(RodinModelsEikonal_HEADERS
  Lagrangian.h
  )

set(RodinModelsAdvection_SRCS
  )

add_library(RodinModelsAdvection INTERFACE)
add_library(Rodin::Models::Advection ALIAS RodinModelsAdvection)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinModelsAdvection
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinModelsAdvection
  INTERFACE
  Rodin::Geometry
  Rodin::Variational)

install(TARGETS RodinModelsAdvection
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Models/Advection
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Models/Advection
)



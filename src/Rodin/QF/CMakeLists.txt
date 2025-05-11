set(RodinQF_HEADERS
  QF1P1.h
  GaussLegendre.h
  GrundmannMoller.h
  GenericPolytopeQuadrature.h
  )

set(RodinQF_SRCS
  QF1P1.cpp
  GaussLegendre.cpp
  GrundmannMoller.cpp
  GenericPolytopeQuadrature.cpp
  )

add_library(RodinQF
  ${RodinQF_SRCS} ${RodinQF_HEADERS})
add_library(Rodin::QF ALIAS RodinQF)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinQF
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinQF
  PUBLIC
  Rodin::Math
  Rodin::Geometry)

install(TARGETS RodinQF
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/QF
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/QF
)

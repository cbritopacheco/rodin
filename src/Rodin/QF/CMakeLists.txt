#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinQF_HEADERS
  QF1P1.h
  QF2P1.h
  GaussLegendre.h
  GrundmannMoller.h
  GenericPolytopeQuadrature.h
)

set(RodinQF_SRCS
  QF1P1.cpp
  QF2P1.cpp
  GaussLegendre.cpp
  GrundmannMoller.cpp
  GenericPolytopeQuadrature.cpp
)

add_library(RodinQF
  ${RodinQF_SRCS}
  ${RodinQF_HEADERS}
)
add_library(Rodin::QF ALIAS RodinQF)
set_target_properties(RodinQF PROPERTIES EXPORT_NAME QF)

# ---- Link / includes -------------------------------------------------------
target_include_directories(RodinQF
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinQF
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

target_link_libraries(RodinQF
  PUBLIC
    Rodin::Math
    Rodin::Geometry
)

# ---- Install ---------------------------------------------------------------
install(
  TARGETS RodinQF
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/QF
)

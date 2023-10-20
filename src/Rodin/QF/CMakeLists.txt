set(RodinQF_HEADERS
  QF1P1.h
  GrundmannMoller.h
  GenericPolytopeQuadrature.h
  )

set(RodinQF_SRCS
  QF1P1.cpp
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


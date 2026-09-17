add_executable(RodinSolidCantileverBeam CantileverBeam.cpp)
target_link_libraries(RodinSolidCantileverBeam
  PUBLIC
  Rodin::Solid
  Rodin::Variational
  Rodin::IO)

add_executable(RodinSolidBlockGravity BlockGravity.cpp)
target_link_libraries(RodinSolidBlockGravity
  PUBLIC
  Rodin::Solid
  Rodin::Variational
  Rodin::IO)

add_executable(RodinSolidActiveContractionPlaneWave ActiveContractionPlaneWave.cpp)
target_link_libraries(RodinSolidActiveContractionPlaneWave
  PUBLIC
  Rodin::Solid
  Rodin::Variational
  Rodin::IO)

add_executable(RodinSolidCantileverBeam CantileverBeam.cpp)
target_link_libraries(RodinSolidCantileverBeam
  PRIVATE
  Rodin::Solid
  Rodin::Variational
  Rodin::IO)

add_executable(RodinSolidBlockGravity BlockGravity.cpp)
target_link_libraries(RodinSolidBlockGravity
  PRIVATE
  Rodin::Solid
  Rodin::Variational
  Rodin::IO)

add_executable(RodinSolidActiveContractionPlaneWave ActiveContractionPlaneWave.cpp)
target_link_libraries(RodinSolidActiveContractionPlaneWave
  PRIVATE
  Rodin::Solid
  Rodin::Variational
  Rodin::IO)

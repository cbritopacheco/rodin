add_executable(MeshIO MeshIO.cpp)
target_link_libraries(MeshIO
  PUBLIC
  Rodin::Geometry
  Rodin::External::MMG)

add_executable(ImplicitDomain ImplicitDomain.cpp)
target_link_libraries(ImplicitDomain
  PUBLIC
  Rodin::Geometry
  Rodin::Variational
  Rodin::External::MMG)

add_executable(DomainDistancing DomainDistancing.cpp)
target_link_libraries(DomainDistancing
  PUBLIC
  Rodin::Geometry
  Rodin::Variational
  Rodin::External::MMG)

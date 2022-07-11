add_executable(ImplicitDomain ImplicitDomain.cpp)
target_link_libraries(ImplicitDomain
  PUBLIC
  Rodin::Mesh
  Rodin::Variational
  Rodin::External::MMG)

add_executable(DomainDistancing DomainDistancing.cpp)
target_link_libraries(DomainDistancing
  PUBLIC
  Rodin::Mesh
  Rodin::Variational
  Rodin::External::MMG)

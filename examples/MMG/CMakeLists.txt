add_executable(MeshIO MeshIO.cpp)
target_link_libraries(MeshIO
  PUBLIC
  Rodin::Geometry
  Rodin::External::MMG)

add_executable(GridFunctionIO GridFunctionIO.cpp)
target_link_libraries(GridFunctionIO
  PUBLIC
  Rodin::Geometry
  Rodin::External::MMG)

add_executable(DomainDiscretization DomainDiscretization.cpp)
target_link_libraries(DomainDiscretization
  PUBLIC
  Rodin::Geometry
  Rodin::External::MMG)

add_executable(MeshOptimization MeshOptimization.cpp)
target_link_libraries(MeshOptimization
  PUBLIC
  Rodin::Geometry
  Rodin::External::MMG)

add_executable(DomainDistancing DomainDistancing.cpp)
target_link_libraries(DomainDistancing
  PUBLIC
  Rodin::Geometry
  Rodin::Variational
  Rodin::External::MMG)

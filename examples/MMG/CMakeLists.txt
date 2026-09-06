add_executable(MeshIO MeshIO.cpp)
target_link_libraries(MeshIO
  PRIVATE
  Rodin::Geometry
  Rodin::MMG)

add_executable(GridFunctionIO GridFunctionIO.cpp)
target_link_libraries(GridFunctionIO
  PRIVATE
  Rodin::Geometry
  Rodin::MMG)

add_executable(DomainDiscretization DomainDiscretization.cpp)
target_link_libraries(DomainDiscretization
  PRIVATE
  Rodin::Geometry
  Rodin::MMG)

add_executable(MeshOptimization MeshOptimization.cpp)
target_link_libraries(MeshOptimization
  PRIVATE
  Rodin::Geometry
  Rodin::MMG)

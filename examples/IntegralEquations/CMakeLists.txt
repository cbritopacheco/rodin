add_executable(NewtonianPotential NewtonianPotential.cpp)
target_link_libraries(NewtonianPotential
  PRIVATE
  Rodin::MMG
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(ElasticDistribution ElasticDistribution.cpp)
target_link_libraries(ElasticDistribution
  PRIVATE
  Rodin::MMG
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )


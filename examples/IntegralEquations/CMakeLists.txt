add_executable(NewtonianPotential NewtonianPotential.cpp)
target_link_libraries(NewtonianPotential
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(ElasticDistribution ElasticDistribution.cpp)
target_link_libraries(ElasticDistribution
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )


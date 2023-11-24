add_executable(EquilibriumDistribution EquilibriumDistribution.cpp)
target_link_libraries(EquilibriumDistribution
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )


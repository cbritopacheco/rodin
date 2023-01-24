add_executable(RodinExampleDGPoisson Poisson.cpp)
add_executable(Rodin::Example::DG::Poisson ALIAS RodinExampleDGPoisson)
target_link_libraries(RodinExampleDGPoisson
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )


add_executable(L2GridFunction L2GridFunction.cpp)
target_link_libraries(L2GridFunction
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

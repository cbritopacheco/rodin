add_executable(AtriumContraction AtriumContraction.cpp)
target_link_libraries(AtriumContraction
  PRIVATE Rodin::Solid Rodin::Variational Rodin::IO Rodin::Solver)

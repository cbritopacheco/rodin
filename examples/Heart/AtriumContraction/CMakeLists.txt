add_executable(AtriumContraction
  AtriumContraction.cpp
  ProblemData.cpp
  Output.cpp)
target_link_libraries(AtriumContraction
  PRIVATE Rodin::Solid Rodin::Variational Rodin::IO Rodin::Solver)

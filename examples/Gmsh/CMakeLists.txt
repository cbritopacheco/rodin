add_executable(IO IO.cpp)
target_link_libraries(IO
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  )


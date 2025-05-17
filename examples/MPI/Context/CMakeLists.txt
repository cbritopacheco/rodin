add_executable(MPI MPI.cpp)
target_link_libraries(MPI
  PUBLIC
  Rodin::MPI
  )


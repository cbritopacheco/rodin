add_executable(MPI MPI.cpp)
target_link_libraries(MPI
  PRIVATE
  Rodin::MPI
  )


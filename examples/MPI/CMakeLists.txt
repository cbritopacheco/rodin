add_executable(MPIMesh MPIMesh.cpp)
target_link_libraries(MPIMesh
  PUBLIC
  Rodin::MPI
  )

add_executable(MPITest MPITest.cpp)
target_link_libraries(MPITest
  PUBLIC
  Rodin::MPI
  )

add_executable(BoostMPITest BoostMPITest.cpp)
target_link_libraries(BoostMPITest
  PUBLIC
  Rodin::MPI
  )

add_executable(BoostMPIBroadcast BoostMPIBroadcast.cpp)
target_link_libraries(BoostMPIBroadcast
  PUBLIC
  Rodin::MPI
  )

add_executable(MPI MPI.cpp)
target_link_libraries(MPI
  PUBLIC
  Rodin::MPI
  )


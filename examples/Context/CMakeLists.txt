if (RODIN_USE_MPI)
  add_executable(MPITest MPITest.cpp)
  target_link_libraries(MPITest
    PUBLIC
    Rodin::Context)

  add_executable(BoostMPITest BoostMPITest.cpp)
  target_link_libraries(BoostMPITest
    PUBLIC
    Rodin::Context)

  add_executable(BoostMPIBroadcast BoostMPIBroadcast.cpp)
  target_link_libraries(BoostMPIBroadcast
    PUBLIC
    Rodin::Context)
endif()

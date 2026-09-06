add_executable(MPI_Boost_Test MPITest.cpp)
target_link_libraries(MPI_Boost_Test PRIVATE Rodin::MPI)

add_executable(MPI_Boost_Broadcast Broadcast.cpp)
target_link_libraries(MPI_Boost_Broadcast PRIVATE Rodin::MPI)


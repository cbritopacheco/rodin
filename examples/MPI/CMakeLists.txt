add_executable(MPITest MPITest.cpp)
target_link_libraries(MPITest PRIVATE Rodin::MPI)

add_subdirectory(Boost)
add_subdirectory(Context)
add_subdirectory(Geometry)
add_subdirectory(Variational)

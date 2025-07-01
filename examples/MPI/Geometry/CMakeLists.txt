add_executable(MPI_Geometry_Mesh Mesh.cpp)
target_link_libraries(MPI_Geometry_Mesh PUBLIC Rodin::MPI)

add_executable(MPI_Geometry_Sharder Sharder.cpp)
target_link_libraries(MPI_Geometry_Sharder PUBLIC Rodin::MPI)


add_executable(MPI_Geometry_Mesh Mesh.cpp)
target_link_libraries(MPI_Geometry_Mesh PRIVATE Rodin::MPI)

add_executable(MPI_Geometry_Sharder Sharder.cpp)
target_link_libraries(MPI_Geometry_Sharder PRIVATE Rodin::MPI)

add_executable(MPI_Geometry_UniformGrid UniformGrid.cpp)
target_link_libraries(MPI_Geometry_UniformGrid PRIVATE Rodin::MPI)

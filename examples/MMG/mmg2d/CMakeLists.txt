add_executable(MeshOptimization MeshOptimization.cpp)
target_link_libraries(MeshOptimization PUBLIC Rodin::External::MMG)

add_executable(MeshGradation MeshGradation.cpp)
target_link_libraries(MeshGradation PUBLIC Rodin::External::MMG)


add_executable(CastMMGMesh2D CastMMGMesh2D.cpp)
target_link_libraries(CastMMGMesh2D PUBLIC Rodin::External::MMG)

add_executable(CastRodinMesh2D CastRodinMesh2D.cpp)
target_link_libraries(CastRodinMesh2D PUBLIC Rodin::External::MMG)

add_executable(CastScalarSolution2D CastScalarSolution2D.cpp)
target_link_libraries(CastScalarSolution2D PUBLIC Rodin::External::MMG)

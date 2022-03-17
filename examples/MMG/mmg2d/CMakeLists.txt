add_executable(MeshOptimization2D MeshOptimization2D.cpp)
target_link_libraries(MeshOptimization2D PUBLIC Rodin::External::MMG)

add_executable(MeshGradation2D MeshGradation2D.cpp)
target_link_libraries(MeshGradation2D PUBLIC Rodin::External::MMG)

add_executable(CastMMGMesh2D CastMMGMesh2D.cpp)
target_link_libraries(CastMMGMesh2D PUBLIC Rodin::External::MMG)

add_executable(CastScalarSolution2D CastScalarSolution2D.cpp)
target_link_libraries(CastScalarSolution2D PUBLIC Rodin::External::MMG)

add_executable(Distancing2D Distancing2D.cpp)
target_link_libraries(Distancing2D PUBLIC Rodin::External::MMG)

add_executable(Redistancing2D Redistancing2D.cpp)
target_link_libraries(Redistancing2D PUBLIC Rodin::External::MMG)

add_executable(ImplicitDomainMeshing2D ImplicitDomainMeshing2D.cpp)
target_link_libraries(ImplicitDomainMeshing2D PUBLIC Rodin::External::MMG)

add_executable(Advection2D Advection2D.cpp)
target_link_libraries(Advection2D PUBLIC Rodin::External::MMG)

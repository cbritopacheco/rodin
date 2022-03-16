add_executable(CastMMGSurfaceMesh CastMMGSurfaceMesh.cpp)
target_link_libraries(CastMMGSurfaceMesh PUBLIC Rodin::External::MMG)

add_executable(CastRodinMeshToMMGSurface CastRodinMeshToMMGSurface.cpp)
target_link_libraries(CastRodinMeshToMMGSurface PUBLIC Rodin::External::MMG)

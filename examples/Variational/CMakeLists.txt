add_executable(Poisson Poisson.cpp)
target_link_libraries(Poisson
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

# add_executable(Elasticity Elasticity.cpp)
# target_link_libraries(Elasticity
#   PUBLIC
#   Rodin::Geometry
#   Rodin::Solver
#   Rodin::Variational
#   )
# 
# add_executable(SurfaceMesh SurfaceMesh.cpp)
# target_link_libraries(SurfaceMesh
#   PUBLIC
#   Rodin::Geometry
#   Rodin::Solver
#   Rodin::Variational
#   )

# add_subdirectory(DG)

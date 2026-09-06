add_executable(Eikonal_FMM FMM.cpp)
target_link_libraries(Eikonal_FMM
  PRIVATE
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational)

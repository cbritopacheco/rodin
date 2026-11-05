add_executable(Eikonal_FMM FMM.cpp)
target_link_libraries(Eikonal_FMM
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational)

add_executable(Eikonal_ReactionEikonal ReactionEikonal.cpp)
target_link_libraries(Eikonal_ReactionEikonal
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational)

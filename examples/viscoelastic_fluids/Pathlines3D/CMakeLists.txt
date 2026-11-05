# Lagrangian pathlines through the XDMF/HDF5 velocity output of
# ArterialLesion3D (sequential: it rebuilds the mesh from the rank pieces).
add_executable(Pathlines3D Pathlines3D.cpp)
target_link_libraries(Pathlines3D
  PUBLIC
  Rodin::Variational
  Rodin::Geometry
  Rodin::Location
  Rodin::IO
  Rodin::Alert)

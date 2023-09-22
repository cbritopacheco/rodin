add_executable(RodinIOMFEM MFEM.cpp)
target_link_libraries(RodinIOMFEM
  PUBLIC
  Rodin::IO
  Rodin::Variational)


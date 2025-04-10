add_executable(Rodin_IO_MEDIT MEDIT.cpp)
target_link_libraries(Rodin_IO_MEDIT
  PUBLIC
  Rodin::IO)

add_executable(RodinIOMFEM MFEM.cpp)
target_link_libraries(RodinIOMFEM
  PUBLIC
  Rodin::IO
  Rodin::Variational)


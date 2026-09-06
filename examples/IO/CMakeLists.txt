add_executable(Rodin_IO_MEDIT MEDIT.cpp)
target_link_libraries(Rodin_IO_MEDIT
  PRIVATE
  Rodin::IO)

add_executable(RodinIOMFEM MFEM.cpp)
target_link_libraries(RodinIOMFEM
  PRIVATE
  Rodin::IO
  Rodin::Variational)

add_executable(RodinIOH1MFEMGridFunction H1MFEMGridFunction.cpp)
target_link_libraries(RodinIOH1MFEMGridFunction
  PRIVATE
  Rodin::IO
  Rodin::Variational)

add_executable(RodinIOXDMF XDMF.cpp)
target_link_libraries(RodinIOXDMF
  PRIVATE
  Rodin::IO
  Rodin::Variational)

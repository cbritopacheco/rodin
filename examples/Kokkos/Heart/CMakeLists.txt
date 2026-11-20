add_executable(KokkosAtriumFlow AtriumFlow.cpp)
target_link_libraries(KokkosAtriumFlow
  PUBLIC
  Rodin::Kokkos::Kernels
  Rodin::Heart
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  Rodin::IO
  )

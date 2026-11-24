add_executable(KokkosBlockGravity BlockGravity.cpp)
target_link_libraries(KokkosBlockGravity
  PUBLIC
  Rodin::Kokkos::Kernels
  Rodin::Solid
  Rodin::Solver
  Rodin::Variational
  Rodin::IO
  )

add_executable(KokkosCantileverBeam CantileverBeam.cpp)
target_link_libraries(KokkosCantileverBeam
  PUBLIC
  Rodin::Kokkos::Kernels
  Rodin::Solid
  Rodin::Solver
  Rodin::Variational
  Rodin::IO
  )

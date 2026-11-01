add_executable(KokkosPoisson Poisson.cpp)
target_link_libraries(KokkosPoisson
  PUBLIC
  Rodin::Kokkos
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  Rodin::IO
  )

add_executable(KokkosMixedBoundary MixedBoundary.cpp)
target_link_libraries(KokkosMixedBoundary
  PUBLIC
  Rodin::Kokkos
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  Rodin::IO
  )

add_executable(KokkosElasticity Elasticity.cpp)
target_link_libraries(KokkosElasticity
  PUBLIC
  Rodin::Kokkos
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  Rodin::IO
  )

add_executable(KokkosStokes Stokes.cpp)
target_link_libraries(KokkosStokes
  PUBLIC
  Rodin::Kokkos
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  Rodin::IO
  )

if (TARGET Rodin::KokkosPETSc)
  add_executable(KokkosSeq_Poisson Seq_Poisson.cpp)
  target_link_libraries(KokkosSeq_Poisson
    PUBLIC
    Rodin::KokkosPETSc
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly
    Rodin::Variational
    Rodin::IO
    )
endif()

add_executable(KokkosNavierStokes NavierStokes.cpp)
target_link_libraries(KokkosNavierStokes
  PUBLIC
  Rodin::Kokkos
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  Rodin::IO
  )

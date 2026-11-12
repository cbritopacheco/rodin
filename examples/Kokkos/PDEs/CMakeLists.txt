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

add_executable(KokkosH1_Assembly H1_Assembly.cpp)
target_link_libraries(KokkosH1_Assembly
  PUBLIC
  Rodin::Kokkos
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  Rodin::Solid
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

if (TARGET Rodin::Kokkos::PETSc)
  add_executable(KokkosSeq_Poisson Seq_Poisson.cpp)
  target_link_libraries(KokkosSeq_Poisson
    PUBLIC
    Rodin::Kokkos::PETSc
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly
    Rodin::Variational
    Rodin::IO
    )

  add_executable(KokkosSeq_NavierStokes Seq_NavierStokes.cpp)
  target_link_libraries(KokkosSeq_NavierStokes
    PUBLIC
    Rodin::Kokkos::PETSc
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly
    Rodin::Variational
    Rodin::IO
    )

  add_executable(KokkosSeq_Elasticity Seq_Elasticity.cpp)
  target_link_libraries(KokkosSeq_Elasticity
    PUBLIC
    Rodin::Kokkos::PETSc
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

if (TARGET Rodin::Kokkos::Kernels)
  add_executable(KokkosKrylov_Poisson Krylov_Poisson.cpp)
  target_link_libraries(KokkosKrylov_Poisson
    PUBLIC
    Rodin::Kokkos::Kernels
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly
    Rodin::Variational
    Rodin::IO
    )
endif()

if (TARGET Rodin::Kokkos::Hypre)
  add_executable(KokkosBoomerAMG_Poisson BoomerAMG_Poisson.cpp)
  target_link_libraries(KokkosBoomerAMG_Poisson
    PUBLIC
    Rodin::Kokkos::Hypre
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly
    Rodin::Variational
    Rodin::IO
    )
endif()

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

# The same example on the PETSc Kokkos path: COO assembly into a PETSc matrix
# of the Kokkos type and an iterative KSP on the device.
if (TARGET Rodin::Kokkos::PETSc)
  add_executable(KokkosAtriumFlowPETSc AtriumFlow.cpp)
  target_compile_definitions(KokkosAtriumFlowPETSc PRIVATE RODIN_ATRIUM_PETSC)
  target_link_libraries(KokkosAtriumFlowPETSc
    PUBLIC
    Rodin::Kokkos::PETSc
    Rodin::Kokkos::Kernels
    Rodin::Heart
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly
    Rodin::Variational
    Rodin::IO
    )

  # Newtonian, Carreau-Yasuda and Owens (log-conformation) blood, PETSc only.
  add_executable(KokkosAtriumRheologyPETSc AtriumRheology.cpp)
  target_link_libraries(KokkosAtriumRheologyPETSc
    PUBLIC
    Rodin::Kokkos::PETSc
    Rodin::Kokkos::Kernels
    Rodin::Heart
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly
    Rodin::Variational
    Rodin::IO
    )
endif()

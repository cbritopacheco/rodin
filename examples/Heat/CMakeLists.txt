if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  # Forced convection (Navier-Stokes + energy, BDF1) in the diverging
  # microchannel of resources/examples/Heat.
  add_executable(ForcedConvection ForcedConvection.cpp)
  target_link_libraries(ForcedConvection
    PUBLIC
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(ForcedConvection PUBLIC Rodin::Scotch)
  endif()
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  # Monolithic (u, p, T) Newton with split orthogonal subgrid scales.
  add_executable(ForcedConvectionNewton ForcedConvectionNewton.cpp)
  target_link_libraries(ForcedConvectionNewton
    PUBLIC
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(ForcedConvectionNewton PUBLIC Rodin::Scotch)
  endif()
endif()

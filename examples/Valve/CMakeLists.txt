if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  # Knife gate valve: stabilised P1/P1 Navier-Stokes (split OSS, Picard,
  # power-law viscosity) and the hydrodynamic force on the gate.
  add_executable(KnifeGateValve KnifeGateValve.cpp)
  target_link_libraries(KnifeGateValve
    PUBLIC
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(KnifeGateValve PUBLIC Rodin::Scotch)
  endif()
endif()

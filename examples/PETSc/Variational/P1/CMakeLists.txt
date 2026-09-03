if (RODIN_USE_MPI)
  add_executable(PETSc_MPI_P1_GridFunction GridFunction.cpp)
  target_link_libraries(PETSc_MPI_P1_GridFunction
    PRIVATE
    Rodin::Math
    Rodin::Solver
    Rodin::Geometry
    Rodin::Variational
    Rodin::PETSc
    )
endif()

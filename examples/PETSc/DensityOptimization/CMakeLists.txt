if (RODIN_USE_MPI)
  add_executable(PETSc_TemperatureMinimization TemperatureMinimization.cpp)
  target_link_libraries(PETSc_TemperatureMinimization
    PRIVATE
    Rodin::Solver
    Rodin::Geometry
    Rodin::Variational
    Rodin::MPI
    Rodin::PETSc
  )

  add_executable(PETSc_TemperatureMinimization3D TemperatureMinimization3D.cpp)
  target_compile_definitions(PETSc_TemperatureMinimization3D
    PRIVATE
    RODIN_EXAMPLE_OUTPUT_DIR="${CMAKE_BINARY_DIR}/tmp"
  )
  target_link_libraries(PETSc_TemperatureMinimization3D
    PRIVATE
    Rodin::Solver
    Rodin::Geometry
    Rodin::Variational
    Rodin::MPI
    Rodin::PETSc
  )
endif()

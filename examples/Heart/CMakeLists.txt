add_executable(CCMLC2014_0D CCMLC2014_0D.cpp)
target_link_libraries(CCMLC2014_0D
  PUBLIC
  Rodin::Heart
  Rodin::Alert
  )

add_executable(CoronarySolid CoronarySolid.cpp)
target_link_libraries(CoronarySolid
  PUBLIC
  Rodin::Heart
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(CoronaryArtery
    CoronaryArtery.cpp
    CoronaryArtery/CoupledLV0DCoronary3D.cpp
    )
  target_link_libraries(CoronaryArtery
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    )

  if (RODIN_USE_SCOTCH)
    target_link_libraries(CoronaryArtery PUBLIC Rodin::Scotch)
  endif()
endif()

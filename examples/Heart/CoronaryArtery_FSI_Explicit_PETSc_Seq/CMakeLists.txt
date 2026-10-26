add_executable(CoronaryArtery_FSI_Explicit_PETSc_Seq CoronaryArtery_FSI_Explicit_PETSc_Seq.cpp)
target_link_libraries(CoronaryArtery_FSI_Explicit_PETSc_Seq
  PUBLIC
  Rodin::Heart
  Rodin::Variational
  Rodin::PETSc
  Rodin::Alert
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly)

target_include_directories(CoronaryArtery_FSI_Explicit_PETSc_Seq
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

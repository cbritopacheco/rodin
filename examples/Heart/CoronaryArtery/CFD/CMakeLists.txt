add_executable(CoronaryArtery
  CoronaryArtery.cpp
  CoupledLV0DCoronary3D.cpp
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

target_include_directories(CoronaryArtery
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

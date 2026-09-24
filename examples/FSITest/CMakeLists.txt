# Manufactured-solution convergence test for the explicit Robin-Robin ALE
# FSI scheme of examples/Heart (2D, meshes generated in the code).
# Sequential, PETSc only.  Regenerate ManufacturedSolution.h with
#   python3 examples/FSITest/generate_mms.py
if (RODIN_USE_PETSC)
  add_executable(FSITest FSITest.cpp)
  target_link_libraries(FSITest
    PUBLIC
    Rodin::Solid
    Rodin::Variational
    Rodin::PETSc
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  target_include_directories(FSITest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../Heart)
endif()

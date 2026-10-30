add_library(AtriumFluid STATIC Fluid.cpp Output.cpp)
target_include_directories(AtriumFluid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(AtriumFluid PUBLIC Rodin::Heart Rodin::Variational Rodin::PETSc
  Rodin::MPI Rodin::Scotch Rodin::Alert Rodin::Solver Rodin::Geometry Rodin::Assembly)

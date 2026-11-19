# Viscoelastic and thixotropic blood in the left atrium (fixed fluid domain):
# the log-conformation solver of examples/viscoelastic_fluids ported to the
# fluid region of the Atrium mesh, with its pressure-driven openings.
add_executable(AtriumViscoelastic AtriumViscoelastic.cpp)
target_include_directories(AtriumViscoelastic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(AtriumViscoelastic PRIVATE Rodin::Heart Rodin::Variational
  Rodin::PETSc Rodin::MPI Rodin::Scotch Rodin::Alert Rodin::Solver Rodin::Geometry
  Rodin::Assembly)

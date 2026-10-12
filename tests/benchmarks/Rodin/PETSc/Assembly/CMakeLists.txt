add_executable(RodinPETScConstrainedAssemblyBenchmarks EXCLUDE_FROM_ALL ConstrainedAssembly.cpp)
target_link_libraries(RodinPETScConstrainedAssemblyBenchmarks PRIVATE
  Rodin::PETSc Rodin::MPI benchmark::benchmark)

add_executable(RodinPETScMassAssemblyBenchmarks EXCLUDE_FROM_ALL MassAssembly.cpp)
target_link_libraries(RodinPETScMassAssemblyBenchmarks PRIVATE
  Rodin::PETSc Rodin::MPI benchmark::benchmark)

add_executable(RodinPETScPhysicsAssemblyBenchmarks EXCLUDE_FROM_ALL PhysicsAssembly.cpp)
target_link_libraries(RodinPETScPhysicsAssemblyBenchmarks PRIVATE
  Rodin::PETSc Rodin::MPI benchmark::benchmark)

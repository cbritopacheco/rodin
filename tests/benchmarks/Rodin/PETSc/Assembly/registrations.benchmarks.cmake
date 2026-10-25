rodin_register_benchmark(RodinPETScConstrainedAssemblyBenchmarks
  PATH Rodin/PETSc/Assembly/Constrained
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinPETScConstrainedAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1 2 4
  CASE_BATCH 256
  MPI_ARGS --assembly_mpi)

rodin_register_benchmark(RodinPETScMassAssemblyBenchmarks
  PATH Rodin/PETSc/Assembly/Mass
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinPETScMassAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1 2 4
  CASE_BATCH 256
  MPI_ARGS --assembly_mpi)

rodin_register_benchmark(RodinPETScPhysicsAssemblyBenchmarks
  PATH Rodin/PETSc/Assembly/Physics
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinPETScPhysicsAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1 2 4
  CASE_BATCH 256
  MPI_ARGS --assembly_mpi)

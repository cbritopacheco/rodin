add_executable(LevelSetSWIFTReconstruction LevelSetSWIFTReconstruction.cpp)
target_link_libraries(LevelSetSWIFTReconstruction
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTReconstructionP2 LevelSetSWIFTReconstruction.cpp)
target_compile_definitions(LevelSetSWIFTReconstructionP2
  PRIVATE
  RODIN_SWIFT_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetSWIFTReconstructionP2
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTReconstruction3D LevelSetSWIFTReconstruction3D.cpp)
target_link_libraries(LevelSetSWIFTReconstruction3D
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTSweep LevelSetSWIFTSweep.cpp)
target_link_libraries(LevelSetSWIFTSweep
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTSweepP2 LevelSetSWIFTSweep.cpp)
target_compile_definitions(LevelSetSWIFTSweepP2
  PRIVATE
  RODIN_SWIFT_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetSWIFTSweepP2
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTReconstruction3DP2 LevelSetSWIFTReconstruction3D.cpp)
target_compile_definitions(LevelSetSWIFTReconstruction3DP2
  PRIVATE
  RODIN_SWIFT_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetSWIFTReconstruction3DP2
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTSweep3D LevelSetSWIFTSweep3D.cpp)
target_link_libraries(LevelSetSWIFTSweep3D
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTSweep3DP2 LevelSetSWIFTSweep3D.cpp)
target_compile_definitions(LevelSetSWIFTSweep3DP2
  PRIVATE
  RODIN_SWIFT_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetSWIFTSweep3DP2
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTAdvection LevelSetSWIFTAdvection.cpp)
target_link_libraries(LevelSetSWIFTAdvection
  PUBLIC
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetSWIFTAdvectionP2 LevelSetSWIFTAdvection.cpp)
target_compile_definitions(LevelSetSWIFTAdvectionP2
  PRIVATE
  RODIN_SWIFT_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetSWIFTAdvectionP2
  PUBLIC
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

foreach(target IN ITEMS
    LevelSetSWIFTReconstruction LevelSetSWIFTReconstructionP2
    LevelSetSWIFTReconstruction3D LevelSetSWIFTReconstruction3DP2
    LevelSetSWIFTSweep LevelSetSWIFTSweepP2
    LevelSetSWIFTSweep3D LevelSetSWIFTSweep3DP2
    LevelSetSWIFTAdvection LevelSetSWIFTAdvectionP2)
  target_link_libraries(${target} PRIVATE Rodin::MMG)
endforeach()


add_executable(LevelSetSWIFTCantilever2D LevelSetSWIFTCantilever2D.cpp)
target_link_libraries(LevelSetSWIFTCantilever2D
  PUBLIC
  Rodin::Advection
  Rodin::Geometry
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational
  )

add_executable(LevelSetSWIFTCantilever2DP2 LevelSetSWIFTCantilever2D.cpp)
target_compile_definitions(LevelSetSWIFTCantilever2DP2
  PRIVATE
  RODIN_SWIFT_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetSWIFTCantilever2DP2
  PRIVATE
  Rodin::Advection
  Rodin::Geometry
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational
  )

if (RODIN_USE_PETSC)
  add_executable(LevelSetSWIFTCantilever3D LevelSetSWIFTCantilever3D.cpp)
  target_link_libraries(LevelSetSWIFTCantilever3D
    PUBLIC
    Rodin::Advection
    Rodin::Geometry
    Rodin::MMG
    Rodin::PETSc
    Rodin::Solver
    Rodin::Variational
    )

  add_executable(LevelSetSWIFTCantilever3DP2 LevelSetSWIFTCantilever3D.cpp)
  target_compile_definitions(LevelSetSWIFTCantilever3DP2
    PRIVATE
    RODIN_SWIFT_P2_DISPLACEMENT=1)
  target_link_libraries(LevelSetSWIFTCantilever3DP2
    PRIVATE
    Rodin::Advection
    Rodin::Geometry
    Rodin::MMG
    Rodin::PETSc
    Rodin::Solver
    Rodin::Variational
    )
endif()

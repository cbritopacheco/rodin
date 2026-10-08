add_executable(Trace Trace.cpp)
target_link_libraries(Trace
  PUBLIC
  Rodin::Geometry)

add_executable(Volume Volume.cpp)
target_link_libraries(Volume
  PUBLIC
  Rodin::Geometry)

add_executable(Displace Displace.cpp)
target_link_libraries(Displace
  PUBLIC
  Rodin::Geometry
  Rodin::Variational)

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

add_executable(Skin Skin.cpp)
target_link_libraries(Skin
  PUBLIC
  Rodin::Geometry)

add_executable(Trim Trim.cpp)
target_link_libraries(Trim
  PUBLIC
  Rodin::Geometry)

add_executable(CCL CCL.cpp)
target_link_libraries(CCL
  PUBLIC
  Rodin::Geometry)

add_executable(Connectivity Connectivity.cpp)
target_link_libraries(Connectivity
  PUBLIC
  Rodin::Geometry)

add_executable(Square Square.cpp)
target_link_libraries(Square
  PUBLIC
  Rodin::Geometry)

add_executable(UniformGrid UniformGrid.cpp)
target_link_libraries(UniformGrid
  PUBLIC
  Rodin::Geometry)

add_executable(ExtrudeTriSurface ExtrudeTriSurface.cpp)
target_link_libraries(ExtrudeTriSurface
  PUBLIC
  Rodin::Geometry)

add_executable(Box Box.cpp)
target_link_libraries(Box
  PUBLIC
  Rodin::Geometry)

add_executable(GreedyPartitioner GreedyPartitioner.cpp)
target_link_libraries(GreedyPartitioner
  PUBLIC
  Rodin::Geometry)

add_executable(BalancedCompactPartitioner BalancedCompactPartitioner.cpp)
target_link_libraries(BalancedCompactPartitioner
  PUBLIC
  Rodin::Geometry)

add_executable(Scale Scale.cpp)
target_link_libraries(Scale
  PUBLIC
  Rodin::Geometry)

add_executable(AnnulusDetection AnnulusDetection.cpp)
target_link_libraries(AnnulusDetection PUBLIC Rodin::Geometry)

add_executable(AnnulusClosing AnnulusClosing.cpp)
target_link_libraries(AnnulusClosing PUBLIC Rodin::Geometry)

add_executable(AttributeDetection AttributeDetection.cpp)
target_link_libraries(AttributeDetection PUBLIC Rodin::Geometry)

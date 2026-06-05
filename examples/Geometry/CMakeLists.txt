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

add_executable(LevelSetLSRReconstruction LevelSetLSRReconstruction.cpp)
target_link_libraries(LevelSetLSRReconstruction
  PUBLIC
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetLSRReconstructionP2 LevelSetLSRReconstruction.cpp)
target_compile_definitions(LevelSetLSRReconstructionP2
  PRIVATE
  RODIN_LSR_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetLSRReconstructionP2
  PUBLIC
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetLSRReconstructionPhiP2 LevelSetLSRReconstruction.cpp)
target_compile_definitions(LevelSetLSRReconstructionPhiP2
  PRIVATE
  RODIN_LSR_PHI_P2=1)
target_link_libraries(LevelSetLSRReconstructionPhiP2
  PUBLIC
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetLSRReconstructionPhiP3 LevelSetLSRReconstruction.cpp)
target_compile_definitions(LevelSetLSRReconstructionPhiP3
  PRIVATE
  RODIN_LSR_PHI_P3=1)
target_link_libraries(LevelSetLSRReconstructionPhiP3
  PUBLIC
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetLSRWavyCircleSweep LevelSetLSRWavyCircleSweep.cpp)
target_link_libraries(LevelSetLSRWavyCircleSweep
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetLSRWavyCircleSweepP2 LevelSetLSRWavyCircleSweep.cpp)
target_compile_definitions(LevelSetLSRWavyCircleSweepP2
  PRIVATE
  RODIN_LSR_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetLSRWavyCircleSweepP2
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetLSRWavyCircleAdvectedSweep
  LevelSetLSRWavyCircleAdvectedSweep.cpp)
target_link_libraries(LevelSetLSRWavyCircleAdvectedSweep
  PUBLIC
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetLSRTwoWavyCirclesMerge LevelSetLSRTwoWavyCirclesMerge.cpp)
target_link_libraries(LevelSetLSRTwoWavyCirclesMerge
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetLSRTwoWavyCirclesMergeP2 LevelSetLSRTwoWavyCirclesMerge.cpp)
target_compile_definitions(LevelSetLSRTwoWavyCirclesMergeP2
  PRIVATE
  RODIN_LSR_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetLSRTwoWavyCirclesMergeP2
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

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

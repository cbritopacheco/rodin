add_executable(Trace Trace.cpp)
target_link_libraries(Trace
  PRIVATE
  Rodin::Geometry)

add_executable(Volume Volume.cpp)
target_link_libraries(Volume
  PRIVATE
  Rodin::Geometry)

add_executable(Displace Displace.cpp)
target_link_libraries(Displace
  PRIVATE
  Rodin::Geometry
  Rodin::Variational)

add_executable(LevelSetWNGIRReconstruction LevelSetWNGIRReconstruction.cpp)
target_link_libraries(LevelSetWNGIRReconstruction
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRReconstructionP2 LevelSetWNGIRReconstruction.cpp)
target_compile_definitions(LevelSetWNGIRReconstructionP2
  PRIVATE
  RODIN_WNGIR_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetWNGIRReconstructionP2
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRReconstruction3D LevelSetWNGIRReconstruction3D.cpp)
target_link_libraries(LevelSetWNGIRReconstruction3D
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRSweep LevelSetWNGIRSweep.cpp)
target_link_libraries(LevelSetWNGIRSweep
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRSweepP2 LevelSetWNGIRSweep.cpp)
target_compile_definitions(LevelSetWNGIRSweepP2
  PRIVATE
  RODIN_WNGIR_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetWNGIRSweepP2
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRReconstruction3DP2 LevelSetWNGIRReconstruction3D.cpp)
target_compile_definitions(LevelSetWNGIRReconstruction3DP2
  PRIVATE
  RODIN_WNGIR_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetWNGIRReconstruction3DP2
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRSweep3D LevelSetWNGIRSweep3D.cpp)
target_link_libraries(LevelSetWNGIRSweep3D
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRSweep3DP2 LevelSetWNGIRSweep3D.cpp)
target_compile_definitions(LevelSetWNGIRSweep3DP2
  PRIVATE
  RODIN_WNGIR_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetWNGIRSweep3DP2
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRAdvection LevelSetWNGIRAdvection.cpp)
target_link_libraries(LevelSetWNGIRAdvection
  PRIVATE
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

add_executable(LevelSetWNGIRAdvectionP2 LevelSetWNGIRAdvection.cpp)
target_compile_definitions(LevelSetWNGIRAdvectionP2
  PRIVATE
  RODIN_WNGIR_P2_DISPLACEMENT=1)
target_link_libraries(LevelSetWNGIRAdvectionP2
  PRIVATE
  Rodin::Advection
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)

foreach(target IN ITEMS
    LevelSetWNGIRReconstruction LevelSetWNGIRReconstructionP2
    LevelSetWNGIRReconstruction3D LevelSetWNGIRReconstruction3DP2
    LevelSetWNGIRSweep LevelSetWNGIRSweepP2
    LevelSetWNGIRSweep3D LevelSetWNGIRSweep3DP2
    LevelSetWNGIRAdvection LevelSetWNGIRAdvectionP2)
  target_link_libraries(${target} PRIVATE Rodin::MMG)
endforeach()

add_executable(Skin Skin.cpp)
target_link_libraries(Skin
  PRIVATE
  Rodin::Geometry)

add_executable(Trim Trim.cpp)
target_link_libraries(Trim
  PRIVATE
  Rodin::Geometry)

add_executable(CCL CCL.cpp)
target_link_libraries(CCL
  PRIVATE
  Rodin::Geometry)

add_executable(Connectivity Connectivity.cpp)
target_link_libraries(Connectivity
  PRIVATE
  Rodin::Geometry)

add_executable(Square Square.cpp)
target_link_libraries(Square
  PRIVATE
  Rodin::Geometry)

add_executable(UniformGrid UniformGrid.cpp)
target_link_libraries(UniformGrid
  PRIVATE
  Rodin::Geometry)

add_executable(Box Box.cpp)
target_link_libraries(Box
  PRIVATE
  Rodin::Geometry)

add_executable(GreedyPartitioner GreedyPartitioner.cpp)
target_link_libraries(GreedyPartitioner
  PRIVATE
  Rodin::Geometry)

add_executable(BalancedCompactPartitioner BalancedCompactPartitioner.cpp)
target_link_libraries(BalancedCompactPartitioner
  PRIVATE
  Rodin::Geometry)

add_executable(Scale Scale.cpp)
target_link_libraries(Scale
  PRIVATE
  Rodin::Geometry)

add_executable(AnnulusDetection AnnulusDetection.cpp)
target_link_libraries(AnnulusDetection PRIVATE Rodin::Geometry)

add_executable(AnnulusClosing AnnulusClosing.cpp)
target_link_libraries(AnnulusClosing PRIVATE Rodin::Geometry)

add_executable(AttributeDetection AttributeDetection.cpp)
target_link_libraries(AttributeDetection PRIVATE Rodin::Geometry)

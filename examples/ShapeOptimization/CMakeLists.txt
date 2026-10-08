add_executable(SimpleCantilever2D SimpleCantilever2D.cpp)
target_link_libraries(SimpleCantilever2D
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::MMG
  )

add_executable(LevelSetCantilever2D LevelSetCantilever2D.cpp)
target_link_libraries(LevelSetCantilever2D
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::MMG
  )

add_executable(LevelSetCantilever3D LevelSetCantilever3D.cpp)
target_link_libraries(LevelSetCantilever3D
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::MMG
  )

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

# add_executable(LevelSetMast2D LevelSetMast2D.cpp)
# target_link_libraries(LevelSetMast2D
#   PUBLIC
#   Rodin::Geometry
#   Rodin::Solver
#   Rodin::Variational
#   Rodin::MMG
#   )
# 
# add_executable(LevelSetArch2D LevelSetArch2D.cpp)
# target_link_libraries(LevelSetArch2D
#   PUBLIC
#   Rodin::Geometry
#   Rodin::Solver
#   Rodin::Variational
#   Rodin::MMG
#   )

# Commented out because it breaks CI
# add_executable(LevelSetEigenvalue2D LevelSetEigenvalue2D.cpp)
# target_link_libraries(LevelSetEigenvalue2D
#   PUBLIC
#   Rodin::Geometry
#   Rodin::Solver
#   Rodin::Variational
#   Rodin::MMG
#   Eigen3::Eigen
#   Spectra::Spectra
#   )

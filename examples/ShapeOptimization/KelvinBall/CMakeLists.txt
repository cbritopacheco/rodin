add_library(KelvinBallSupport STATIC EXCLUDE_FROM_ALL
  Common.cpp
  Configuration.cpp
  Metrics.cpp
  RotatedNitscheIntegrator.cpp
  SewedOutput.cpp
  Sphere.cpp)
target_link_libraries(KelvinBallSupport PUBLIC
  Rodin::Adaptation
  Rodin::Assembly
  Rodin::Geometry
  Rodin::IO
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational)

add_executable(KelvinBall EXCLUDE_FROM_ALL
  KelvinBall.cpp
  KelvinBallOptimization.cpp)
target_link_libraries(KelvinBall PRIVATE KelvinBallSupport Rodin::Advection)

add_executable(KelvinBallSphere EXCLUDE_FROM_ALL KelvinBallSphere.cpp)
target_link_libraries(KelvinBallSphere PRIVATE KelvinBallSupport)

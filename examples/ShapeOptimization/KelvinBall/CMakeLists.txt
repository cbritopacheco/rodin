add_executable(KelvinBall EXCLUDE_FROM_ALL
  KelvinBall.cpp
  KelvinBallOptimization.cpp
  Common.cpp
  Configuration.cpp
  Metrics.cpp
  RotatedNitscheIntegrator.cpp
  SewedOutput.cpp
  Sphere.cpp)
target_link_libraries(KelvinBall
  PRIVATE
  Rodin::Advection
  Rodin::Assembly
  Rodin::Geometry
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational)

add_executable(KelvinBallSphere EXCLUDE_FROM_ALL
  KelvinBallSphere.cpp
  Common.cpp
  Configuration.cpp
  Metrics.cpp
  RotatedNitscheIntegrator.cpp
  SewedOutput.cpp
  Sphere.cpp)
target_link_libraries(KelvinBallSphere
  PRIVATE
  Rodin::Assembly
  Rodin::Geometry
  Rodin::IO
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational)

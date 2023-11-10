add_executable(FMM FMM.cpp)
target_link_libraries(FMM
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(SpaldingTucker SpaldingTucker.cpp)
target_link_libraries(SpaldingTucker
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

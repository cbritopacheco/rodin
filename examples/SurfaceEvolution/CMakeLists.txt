add_executable(FirePropagation FirePropagation.cpp)
target_link_libraries(FirePropagation
  PUBLIC
  Rodin::Mesh
  Rodin::Variational
  Rodin::External::MMG
  )


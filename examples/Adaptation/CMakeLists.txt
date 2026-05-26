add_executable(AdvectedGraphCutCurvedInterface AdvectedGraphCutCurvedInterface.cpp)
target_link_libraries(AdvectedGraphCutCurvedInterface
  PUBLIC
  Rodin::Adaptation
  Rodin::Geometry
  Rodin::IO
  Rodin::Variational)

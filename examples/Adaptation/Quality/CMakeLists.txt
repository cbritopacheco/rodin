add_executable(Adaptation_Quality Quality.cpp)
target_link_libraries(Adaptation_Quality
  PRIVATE Rodin::Geometry Rodin::Variational)

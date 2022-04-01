add_executable(LinePlot LinePlot.cpp)
target_link_libraries(LinePlot
  PUBLIC
  Rodin::Core
  Rodin::Plot
  )

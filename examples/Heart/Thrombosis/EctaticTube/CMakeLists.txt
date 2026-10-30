add_executable(EctaticTubeThrombosis EctaticTubeThrombosis.cpp)
target_link_libraries(EctaticTubeThrombosis
  PUBLIC
  Rodin::Heart
  Rodin::Variational
  Rodin::PETSc
  Rodin::Alert
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly)

target_include_directories(EctaticTubeThrombosis
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../..)

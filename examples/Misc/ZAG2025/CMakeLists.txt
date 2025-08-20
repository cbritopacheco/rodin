# add_executable(MinOneState MinOneState.cpp Update.cpp)
# target_link_libraries(MinOneState
#   PUBLIC
#   Rodin::Geometry
#   Rodin::Solver
#   Rodin::Variational
#   Rodin::External::MMG
#   )

add_executable(AmstutzMethod AmstutzMethod.cpp)
target_link_libraries(AmstutzMethod
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

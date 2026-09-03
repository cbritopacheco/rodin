# add_executable(MinOneState MinOneState.cpp Update.cpp)
# target_link_libraries(MinOneState
#   PRIVATE
#   Rodin::Geometry
#   Rodin::Solver
#   Rodin::Variational
#   Rodin::External::MMG
#   )

add_executable(AmstutzMethod AmstutzMethod.cpp)
target_link_libraries(AmstutzMethod
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(P1Potential P1Potential.cpp)
target_link_libraries(P1Potential
  PUBLIC Rodin::Variational Rodin::External::MMG)

add_executable(P1GridFunctionIO P1GridFunctionIO.cpp)
target_link_libraries(P1GridFunctionIO PUBLIC Rodin::Variational)

add_executable(ScalarP1Build ScalarP1Build.cpp)
target_link_libraries(ScalarP1Build PUBLIC Rodin::Variational)

add_executable(ScalarP1Projection ScalarP1Projection.cpp)
target_link_libraries(ScalarP1Projection PUBLIC Rodin::Variational)

add_executable(VectorP1Build VectorP1Build.cpp)
target_link_libraries(VectorP1Build PUBLIC Rodin::Variational)

add_executable(VectorP1Projection VectorP1Projection.cpp)
target_link_libraries(VectorP1Projection PUBLIC Rodin::Variational)

add_executable(P1GridFunctionGradient P1GridFunctionGradient.cpp)
target_link_libraries(P1GridFunctionGradient PUBLIC Rodin::Variational)

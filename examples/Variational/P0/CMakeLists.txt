add_executable(ScalarP0Build ScalarP0Build.cpp)
target_link_libraries(ScalarP0Build PUBLIC Rodin::Variational)

add_executable(ScalarP0Projection ScalarP0Projection.cpp)
target_link_libraries(ScalarP0Projection PUBLIC Rodin::Variational)

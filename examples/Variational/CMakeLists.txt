add_executable(Normal Normal.cpp)
target_link_libraries(Normal PUBLIC Rodin::Variational)

add_executable(BooleanFunction BooleanFunction.cpp)
target_link_libraries(BooleanFunction PUBLIC Rodin::Variational)

add_subdirectory(P1)

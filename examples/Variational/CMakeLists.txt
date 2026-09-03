add_executable(Normal Normal.cpp)
target_link_libraries(Normal PRIVATE Rodin::Variational)

add_executable(BooleanFunction BooleanFunction.cpp)
target_link_libraries(BooleanFunction PRIVATE Rodin::Variational)

add_subdirectory(P0)
add_subdirectory(P1)
add_subdirectory(H1)

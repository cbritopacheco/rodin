add_executable(KokkosMPIPoissonScaling PoissonScaling.cpp)
target_link_libraries(KokkosMPIPoissonScaling PUBLIC Rodin::Kokkos::MPI)

add_executable(KokkosMPIOseenScaling OseenScaling.cpp)
target_link_libraries(KokkosMPIOseenScaling PUBLIC Rodin::Kokkos::MPI)

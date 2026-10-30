add_executable(InterfaceTraces InterfaceTraces.cpp)
target_link_libraries(InterfaceTraces PUBLIC Rodin::Geometry Rodin::Variational)

add_executable(DGPoisson Poisson.cpp)
target_link_libraries(DGPoisson PUBLIC Rodin::Geometry Rodin::QF Rodin::Variational Rodin::IO)

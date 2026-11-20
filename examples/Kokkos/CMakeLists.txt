if (TARGET Rodin::Kokkos::PETSc)
  add_executable(KokkosOseenAssembly OseenAssembly.cpp)
  target_link_libraries(KokkosOseenAssembly PUBLIC Rodin::Kokkos::PETSc)
endif()

add_subdirectory(PDEs)
add_subdirectory(Solid)
add_subdirectory(Heart)

if (TARGET Rodin::Kokkos::MPI AND TARGET Rodin::Kokkos::PETSc)
  add_subdirectory(MPI)
endif()

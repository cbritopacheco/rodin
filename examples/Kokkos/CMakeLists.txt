if (TARGET Rodin::Kokkos::PETSc)
  add_executable(KokkosOseenAssembly OseenAssembly.cpp)
  target_link_libraries(KokkosOseenAssembly PUBLIC Rodin::Kokkos::PETSc)
endif()

add_subdirectory(PDEs)
add_subdirectory(Solid)

if (TARGET Rodin::Kokkos::MPI)
  add_subdirectory(MPI)
endif()

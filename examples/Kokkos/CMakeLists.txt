if (TARGET Rodin::KokkosPETSc)
  add_executable(KokkosOseenAssembly OseenAssembly.cpp)
  target_link_libraries(KokkosOseenAssembly PUBLIC Rodin::KokkosPETSc)
endif()

add_subdirectory(PDEs)
add_subdirectory(Solid)

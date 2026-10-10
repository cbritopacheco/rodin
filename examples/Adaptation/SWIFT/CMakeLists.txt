foreach(degree IN ITEMS 1 2 3)
  add_executable(SWIFT_ReconstructionP${degree} Reconstruction.cpp)
  target_compile_definitions(SWIFT_ReconstructionP${degree}
    PRIVATE RODIN_SWIFT_RECONSTRUCTION_DEGREE=${degree})
  target_link_libraries(SWIFT_ReconstructionP${degree}
    PRIVATE Rodin::Geometry Rodin::Solver Rodin::Variational)
endforeach()

add_subdirectory(Witness)

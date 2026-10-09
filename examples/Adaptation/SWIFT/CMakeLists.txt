foreach(degree IN ITEMS 1 2 3)
  add_executable(SWIFT_ReconstructionP${degree} ReconstructionP${degree}.cpp)
  target_link_libraries(SWIFT_ReconstructionP${degree}
    PRIVATE Rodin::Geometry Rodin::Solver Rodin::Variational)
endforeach()

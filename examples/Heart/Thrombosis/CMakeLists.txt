# Cada ejemplo se agrega solo si su carpeta esta presente: asi un clon que
# tenga parte del arbol (p. ej. el que corre en pc9) configura igual.
function(add_example_if_present dir)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/CMakeLists.txt")
    add_subdirectory(${dir})
  else()
    message(STATUS "Thrombosis: se saltea ${dir} (no esta en este clon)")
  endif()
endfunction()

if (RODIN_USE_PETSC)
  add_example_if_present(EctaticTube)
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_example_if_present(Aneurysm)
  add_example_if_present(PTAGP)
  add_example_if_present(LeftAtrium2DPTAGX)
  add_example_if_present(PTAGXBench)
endif()

# ---- Add main targets ------------------------------------------------------
set(Rodin_HEADERS
  Core.h
  Alert.h
  Mesh.h
  # Plot module is broken. See below.
  # Plot.h
  Solver.h
  Variational.h
  )
add_library(Rodin INTERFACE)
add_library(Rodin::Rodin ALIAS Rodin)
target_sources(Rodin INTERFACE ${Rodin_HEADERS})


if (RODIN_WITH_PLOT)
# The Plot module is currently broken and barely has any code. So we don't add
# it to the build list.  It is a dependency issue with Magnum + mfem which
# causes the build to fail on some systems. Also I might have modified a little
# the includes so that you will get a file not found error during compilation,
# which should be easy to fix.

# add_subdirectory(Plot)
endif()
add_subdirectory(Core)
add_subdirectory(Mesh)
add_subdirectory(Test)
add_subdirectory(Alert)
add_subdirectory(Solver)
add_subdirectory(Utility)
add_subdirectory(Variational)

target_include_directories(Rodin
  INTERFACE
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/src
  )
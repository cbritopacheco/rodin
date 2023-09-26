# ---- Add main targets ------------------------------------------------------
set(Rodin_HEADERS
  Array.h
  Math.h
  Alert.h
  Geometry.h
  Plot.h
  Solver.h
  Variational.h
  )

# ---- Set compile options for all the submodules ----------------------------
rodin_add_compile_options(LANG C CXX OPTIONS -Wno-old-style-cast -Wno-unused-parameter)

add_library(Rodin INTERFACE)
add_library(Rodin::Rodin ALIAS Rodin)
target_sources(Rodin INTERFACE ${Rodin_HEADERS})
target_link_libraries(Rodin INTERFACE Eigen3::Eigen)

configure_file(Configure.h.in Configure.h @ONLY)

target_include_directories(Rodin
  INTERFACE ${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)

add_subdirectory(QF)
add_subdirectory(Math)
add_subdirectory(IO)
add_subdirectory(Test)
add_subdirectory(Alert)
add_subdirectory(Solver)
add_subdirectory(Utility)
add_subdirectory(Geometry)
add_subdirectory(Variational)
add_subdirectory(FormLanguage)

if (RODIN_WITH_PLOT)
  add_subdirectory(Plot)
endif()


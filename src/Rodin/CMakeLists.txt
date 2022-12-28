# ---- Add main targets ------------------------------------------------------
set(Rodin_HEADERS
  Math.h
  Alert.h
  Geometry.h
  Plot.h
  Solver.h
  Variational.h
  )

# ---- Set compile options for all the submodules ----------------------------
if (${CMAKE_BUILD_TYPE} MATCHES Release)
  rodin_add_compile_options(
    LANG CXX
    OPTIONS
    -w
    -O3
    -march=native)
elseif (${CMAKE_BUILD_TYPE} MATCHES Debug)
  rodin_add_compile_options(
    LANG CXX
    OPTIONS
    -g
    -O0
    -Wall
    -Werror
    -Wextra
    -Wpedantic
    -Wno-old-style-cast
    -Wno-unused-parameter
    -Wno-zero-as-null-pointer-constant
    -Wno-error=extra
    -Wno-error=pedantic
    -Wno-error=cast-function-type
    -Wno-error=deprecated-declarations
    -Wno-error=zero-as-null-pointer-constant)
endif()

add_library(Rodin INTERFACE)
add_library(Rodin::Rodin ALIAS Rodin)
target_sources(Rodin INTERFACE ${Rodin_HEADERS})

configure_file(Configure.h.in Configure.h @ONLY)

target_include_directories(Rodin
  INTERFACE
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/src
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

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


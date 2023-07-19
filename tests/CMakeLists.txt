enable_testing()

include(GoogleTest)

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
    -fPIC
    -Wall
    -Werror
    -Wextra
    -Wpedantic
    -Wno-old-style-cast
    -Wno-cast-function-type
    -Wno-unused-parameter
    -Wno-zero-as-null-pointer-constant
    -Wno-error=extra
    -Wno-error=pedantic
    -Wno-error=sign-compare
    -Wno-error=cast-function-type
    -Wno-error=deprecated-declarations
    -Wno-error=zero-as-null-pointer-constant)
else()
  message(FATAL_ERROR "Unrecognized build: ${CMAKE_BUILD_TYPE}")
endif()

if (RODIN_BUILD_UNIT_TESTS)
  add_subdirectory(unit)
endif()

if (RODIN_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()


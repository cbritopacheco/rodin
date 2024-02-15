enable_testing()
include(GoogleTest)

rodin_add_compile_options(
  LANG CXX
  OPTIONS
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
  -Wno-error=unused-variable
  -Wno-error=unused-but-set-variable
  -Wno-error=cast-function-type
  -Wno-error=deprecated-declarations
  -Wno-error=zero-as-null-pointer-constant)

if (RODIN_BUILD_UNIT_TESTS)
  add_subdirectory(unit)
endif()

if (RODIN_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()


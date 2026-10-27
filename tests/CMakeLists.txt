enable_testing()
include(GoogleTest)

# Record which execution profiles can use each GoogleTest executable before
# building it. Forward its registration unchanged, including filters/properties.
function(rodin_discover_tests target)
  cmake_parse_arguments(PARSE_ARGV 1 _registration
    "NO_PRETTY_TYPES;NO_PRETTY_VALUES"
    "TEST_PREFIX;TEST_SUFFIX;WORKING_DIRECTORY;TEST_LIST;DISCOVERY_TIMEOUT;XML_OUTPUT_DIR;DISCOVERY_MODE"
    "EXTRA_ARGS;DISCOVERY_EXTRA_ARGS;PROPERTIES;TEST_FILTER")
  set(_labels)
  list(FIND _registration_PROPERTIES LABELS _label_index)
  if(NOT _label_index EQUAL -1)
    math(EXPR _label_index "${_label_index} + 1")
    list(GET _registration_PROPERTIES ${_label_index} _labels)
  endif()
  set(_profiles)
  if("unit" IN_LIST _labels AND NOT "slow" IN_LIST _labels)
    list(APPEND _profiles unit)
  endif()
  if("manufactured" IN_LIST _labels AND NOT "slow" IN_LIST _labels)
    list(APPEND _profiles manufactured)
  endif()
  if("convergence" IN_LIST _labels)
    list(APPEND _profiles convergence-baseline)
  endif()
  if("slow" IN_LIST _labels)
    list(APPEND _profiles slow)
  endif()
  if(("unit" IN_LIST _labels AND NOT "distributed" IN_LIST _labels) OR
      ("distributed" IN_LIST _labels AND NOT "slow" IN_LIST _labels) OR
      ("manufactured" IN_LIST _labels AND NOT "slow" IN_LIST _labels))
    list(APPEND _profiles coverage)
  endif()
  set_property(TARGET ${target} APPEND PROPERTY RODIN_TEST_PROFILES ${_profiles})
  if(CMAKE_VERSION VERSION_LESS 3.18)
    # CMake 3.16's GoogleTest module also expands ARGN as a list.
    gtest_discover_tests(${target} ${ARGN})
  else()
    # Preserve argument boundaries, including semicolons and empty arguments.
    # Expanding ARGN directly would alter modern GoogleTest registrations.
    set(_call "gtest_discover_tests(")
    math(EXPR _last "${ARGC} - 1")
    foreach(_index RANGE 0 ${_last})
      set(_value "${ARGV${_index}}")
      set(_bracket "=")
      string(FIND "${_value}" "]${_bracket}]" _position)
      while(NOT _position EQUAL -1)
        string(APPEND _bracket "=")
        string(FIND "${_value}" "]${_bracket}]" _position)
      endwhile()
      string(APPEND _call " [${_bracket}[${_value}]${_bracket}]")
    endforeach()
    cmake_language(EVAL CODE "${_call})")
  endif()
endfunction()

foreach(_cache_mode OFF ON)
  add_test(NAME CMake.ThirdPartyIsolation.${_cache_mode}
    COMMAND ${CMAKE_COMMAND}
      -S ${PROJECT_SOURCE_DIR}/tests/cmake/ThirdPartyIsolation
      -B ${CMAKE_CURRENT_BINARY_DIR}/ThirdPartyIsolation-${_cache_mode}
      -DISOLATION_EXISTING_CACHE=${_cache_mode}
      -DISOLATION_WITH_PY=${RODIN_WITH_PY}
      -DISOLATION_WITH_MCSS=${RODIN_USE_MCSS})
  set_tests_properties(CMake.ThirdPartyIsolation.${_cache_mode}
    PROPERTIES LABELS unit)
endforeach()

rodin_add_compile_options(
  LANG C CXX
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
  -Wno-error=old-style-cast
  -Wno-error=unused-variable
  -Wno-error=unused-but-set-variable
  -Wno-error=cast-function-type
  -Wno-error=deprecated-declarations
  -Wno-error=zero-as-null-pointer-constant
  -Wno-error=deprecated-anon-enum-enum-conversion
  -Wno-error=deprecated-builtins
  -Wno-error=array-bounds
  )

if (RODIN_BUILD_UNIT_TESTS)
  add_subdirectory(unit)
endif()

if (RODIN_BUILD_MANUFACTURED_TESTS)
  add_subdirectory(manufactured)
endif()

if (RODIN_BUILD_CONVERGENCE_TESTS)
  add_subdirectory(convergence)
endif()

if (RODIN_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

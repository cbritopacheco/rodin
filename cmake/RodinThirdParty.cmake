# Configure bundled projects without exporting their private choices to Rodin.
include_guard(GLOBAL)

function(rodin_add_mmg)
  # MMG writes these shared cache entries, including compiler flags with FORCE.
  # Preserve existing user settings and remove entries created only by MMG.
  set(_protected USE_SCOTCH BUILD_TESTING BUILD_DOC CMAKE_BUILD_TYPE)
  foreach(_configuration MAINTAINER RELWITHASSERT)
    foreach(_kind C CXX EXE_LINKER SHARED_LINKER STATIC_LINKER)
      list(APPEND _protected CMAKE_${_kind}_FLAGS_${_configuration})
    endforeach()
  endforeach()
  foreach(_variable IN LISTS _protected)
    get_property(_exists_${_variable} CACHE ${_variable} PROPERTY TYPE SET)
    if(_exists_${_variable})
      get_property(_value_${_variable} CACHE ${_variable} PROPERTY VALUE)
      get_property(_type_${_variable} CACHE ${_variable} PROPERTY TYPE)
      get_property(_help_${_variable} CACHE ${_variable} PROPERTY HELPSTRING)
      get_property(_advanced_${_variable} CACHE ${_variable} PROPERTY ADVANCED)
      get_property(_strings_${_variable} CACHE ${_variable} PROPERTY STRINGS)
    endif()
  endforeach()

  # Typed temporary entries also support MMG's old policies on CMake 3.16.
  # These are restored below, not persistent changes to the user's cache.
  set(USE_SCOTCH OFF CACHE STRING "Private MMG configuration" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "Private MMG configuration" FORCE)
  set(BUILD_DOC OFF CACHE BOOL "Private MMG configuration" FORCE)
  set(USE_SCOTCH OFF)
  set(BUILD_TESTING OFF)
  set(BUILD_DOC OFF)
  add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/mmg
    ${PROJECT_BINARY_DIR}/third-party/mmg EXCLUDE_FROM_ALL)

  foreach(_variable IN LISTS _protected)
    if(_exists_${_variable})
      set(${_variable} "${_value_${_variable}}" CACHE
        ${_type_${_variable}} "${_help_${_variable}}" FORCE)
      set_property(CACHE ${_variable} PROPERTY ADVANCED "${_advanced_${_variable}}")
      set_property(CACHE ${_variable} PROPERTY STRINGS "${_strings_${_variable}}")
    else()
      unset(${_variable} CACHE)
    endif()
  endforeach()
endfunction()

function(rodin_add_googletest)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
  set(gtest_force_shared_crt ON)
  add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/googletest
    ${PROJECT_BINARY_DIR}/third-party/googletest EXCLUDE_FROM_ALL)
  # GoogleTest enables -Wall; silence a warning in its vendored GCC code.
  foreach(_target gtest gtest_main gmock gmock_main)
    if(TARGET ${_target})
      target_compile_options(${_target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>)
    endif()
  endforeach()
endfunction()

function(rodin_add_pybind11)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
  set(PYBIND11_TEST OFF)
  set(PYBIND11_INSTALL OFF)
  set(PYBIND11_FINDPYTHON ON)
  set(Python_FIND_REGISTRY LAST)
  set(Python_FIND_VIRTUALENV FIRST)
  find_package(Python 3.6 COMPONENTS Interpreter Development REQUIRED)
  add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/pybind11
    ${PROJECT_BINARY_DIR}/third-party/pybind11 EXCLUDE_FROM_ALL)
endfunction()

function(rodin_add_termcolor)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
  set(TERMCOLOR_TESTS OFF)
  add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/termcolor
    ${PROJECT_BINARY_DIR}/third-party/termcolor EXCLUDE_FROM_ALL)
  # Alert's public headers include termcolor. Preserve Rodin's existing
  # header installation without installing the vendor's package metadata.
  install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/third-party/termcolor/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
endfunction()

# m.css is a Python documentation tool, not a CMake subproject.
function(rodin_add_mcss configuration)
  set(Python3_FIND_REGISTRY LAST)
  set(Python3_FIND_VIRTUALENV FIRST)
  find_package(Python3 3.6 REQUIRED COMPONENTS Interpreter)
  find_package(LATEX REQUIRED COMPONENTS PDFLATEX)
  add_custom_target(RodinDoxygen
    COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/third-party/m.css/documentation/doxygen.py
      "${configuration}"
    COMMENT "Generate API documentation with m.css"
    VERBATIM)
endfunction()

function(rodin_add_benchmark)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
  if(RODIN_LTO)
    set(BENCHMARK_ENABLE_LTO ON)
  endif()
  set(BENCHMARK_ENABLE_TESTING OFF)
  set(BENCHMARK_ENABLE_WERROR OFF)
  add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/benchmark
    ${PROJECT_BINARY_DIR}/third-party/benchmark EXCLUDE_FROM_ALL)
endfunction()

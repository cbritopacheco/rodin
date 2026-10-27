cmake_minimum_required(VERSION 3.16)
project(ThirdPartyIsolation LANGUAGES C CXX)
include(GNUInstallDirs)

get_filename_component(RODIN_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
include(${RODIN_SOURCE_DIR}/cmake/RodinThirdParty.cmake)
set(PROJECT_SOURCE_DIR "${RODIN_SOURCE_DIR}")

# Contradictory parent settings must remain intact, but must not enable MMG's
# Scotch, tests or documentation. Test both new and already-populated caches.
set(RODIN_USE_SCOTCH ON)
if(ISOLATION_EXISTING_CACHE)
  set(USE_SCOTCH ON CACHE STRING "Parent Scotch choice" FORCE)
  set_property(CACHE USE_SCOTCH PROPERTY STRINGS ON OFF)
  set(BUILD_TESTING ON CACHE BOOL "Parent testing choice" FORCE)
  set(BUILD_DOC ON CACHE BOOL "Parent documentation choice" FORCE)
  set(CMAKE_C_FLAGS_MAINTAINER "-DPARENT_MAINTAINER" CACHE STRING
    "Parent compiler flags" FORCE)
endif()

get_cmake_property(_before CACHE_VARIABLES)
list(FILTER _before INCLUDE REGEX
  "^(USE_SCOTCH|BUILD_TESTING|BUILD_DOC|CMAKE_BUILD_TYPE|CMAKE_.*_FLAGS_.*)$")
foreach(_variable IN LISTS _before)
  get_property(_saved_${_variable} CACHE ${_variable} PROPERTY VALUE)
endforeach()
rodin_add_mmg()
foreach(_variable IN LISTS _before)
  get_property(_value CACHE ${_variable} PROPERTY VALUE)
  if(NOT "${_value}" STREQUAL "${_saved_${_variable}}")
    message(FATAL_ERROR "MMG changed parent cache entry ${_variable}")
  endif()
endforeach()
if(NOT ISOLATION_EXISTING_CACHE)
  foreach(_variable USE_SCOTCH BUILD_TESTING BUILD_DOC CMAKE_C_FLAGS_MAINTAINER)
    if(DEFINED CACHE{${_variable}})
      message(FATAL_ERROR "MMG leaked ${_variable} into the parent cache")
    endif()
  endforeach()
endif()
if(NOT RODIN_USE_SCOTCH)
  message(FATAL_ERROR "Rodin Scotch support was disabled")
endif()
get_directory_property(_definitions DIRECTORY
  "${RODIN_SOURCE_DIR}/third-party/mmg" COMPILE_DEFINITIONS)
if("USE_SCOTCH" IN_LIST _definitions)
  message(FATAL_ERROR "MMG was compiled with Scotch")
endif()
get_target_property(_links libmmg3d_a LINK_LIBRARIES)
if("${_links}" MATCHES "[Ss][Cc][Oo][Tt][Cc][Hh]")
  message(FATAL_ERROR "MMG links Scotch: ${_links}")
endif()

set(gtest_force_shared_crt OFF CACHE BOOL "Parent CRT choice" FORCE)
set(BENCHMARK_ENABLE_TESTING ON CACHE BOOL "Parent benchmark tests" FORCE)
set(BENCHMARK_ENABLE_WERROR ON CACHE BOOL "Parent benchmark warnings" FORCE)
set(BENCHMARK_ENABLE_LTO OFF CACHE BOOL "Parent benchmark LTO" FORCE)
set(RODIN_LTO ON)
rodin_add_googletest()
rodin_add_benchmark()
if(NOT "$CACHE{gtest_force_shared_crt}" STREQUAL "OFF" OR
   NOT "$CACHE{BENCHMARK_ENABLE_TESTING}" STREQUAL "ON" OR
   NOT "$CACHE{BENCHMARK_ENABLE_WERROR}" STREQUAL "ON" OR
   NOT "$CACHE{BENCHMARK_ENABLE_LTO}" STREQUAL "OFF")
  message(FATAL_ERROR "Bundled test libraries overwrote parent cache choices")
endif()
get_directory_property(_benchmark_tests DIRECTORY
  "${RODIN_SOURCE_DIR}/third-party/benchmark" DEFINITION BENCHMARK_ENABLE_TESTING)
if(_benchmark_tests)
  message(FATAL_ERROR "Private benchmark testing override was ignored")
endif()

set(TERMCOLOR_TESTS ON CACHE BOOL "Parent termcolor tests" FORCE)
rodin_add_termcolor()
if(NOT "$CACHE{TERMCOLOR_TESTS}" STREQUAL "ON" OR TARGET test_termcolor)
  message(FATAL_ERROR "Private termcolor test override was ignored")
endif()
if(NOT TARGET termcolor::termcolor)
  message(FATAL_ERROR "Bundled termcolor target was not published")
endif()
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/termcolor.cpp"
  "#include <termcolor/termcolor.hpp>\n#include <iostream>\nint main() { std::cout << termcolor::reset; }\n")
add_executable(IsolationTermcolor "${CMAKE_CURRENT_BINARY_DIR}/termcolor.cpp")
target_link_libraries(IsolationTermcolor PRIVATE termcolor::termcolor)

# Optional bindings/documentation must not impose Python or LaTeX requirements
# on Rodin's C++-only test configurations.
if(ISOLATION_WITH_PY)
  set(PYBIND11_TEST ON CACHE BOOL "Parent pybind11 tests" FORCE)
  set(PYBIND11_INSTALL ON CACHE BOOL "Parent pybind11 install" FORCE)
  set(PYBIND11_FINDPYTHON OFF CACHE STRING "Parent Python discovery mode" FORCE)
  set(CMAKE_PREFIX_PATH "parent-prefix")
  rodin_add_pybind11()
  if(NOT "$CACHE{PYBIND11_TEST}" STREQUAL "ON" OR
     NOT "$CACHE{PYBIND11_INSTALL}" STREQUAL "ON" OR
     NOT "$CACHE{PYBIND11_FINDPYTHON}" STREQUAL "OFF" OR
     NOT CMAKE_PREFIX_PATH STREQUAL "parent-prefix")
    message(FATAL_ERROR "Bundled header libraries changed parent configuration")
  endif()
  if(TARGET pybind11_tests)
    message(FATAL_ERROR "Private vendor tests were enabled")
  endif()
  if(NOT TARGET pybind11::module)
    message(FATAL_ERROR "Bundled pybind11 target was not published")
  endif()

  # Exercise targets and pybind11_extension outside their function scopes,
  # as Alert and the Python bindings do in the real build.
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/headers.cpp"
    "#include <termcolor/termcolor.hpp>\n#include <pybind11/pybind11.h>\nPYBIND11_MODULE(isolation, m) { m.def(\"answer\", [] { return 42; }); }\n")
  add_library(IsolationHeaders MODULE "${CMAKE_CURRENT_BINARY_DIR}/headers.cpp")
  set_target_properties(IsolationHeaders PROPERTIES OUTPUT_NAME isolation)
  target_link_libraries(IsolationHeaders PRIVATE termcolor::termcolor pybind11::module)
  pybind11_extension(IsolationHeaders)
  get_target_property(_suffix IsolationHeaders SUFFIX)
  if(NOT _suffix)
    message(FATAL_ERROR "Python module suffix was lost at the wrapper boundary")
  endif()
endif()
if(ISOLATION_WITH_MCSS)
  rodin_add_mcss("${CMAKE_CURRENT_BINARY_DIR}/conf.py")
  if(NOT TARGET RodinDoxygen)
    message(FATAL_ERROR "m.css documentation target was not published")
  endif()
  if(DEFINED Python3_FIND_REGISTRY OR DEFINED Python3_FIND_VIRTUALENV)
    message(FATAL_ERROR "m.css discovery settings escaped the wrapper")
  endif()
endif()
message(STATUS "Third-party isolation checks passed")

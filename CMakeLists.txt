cmake_minimum_required (VERSION 3.16)

cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0076 NEW)

# Prohibit in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited.")
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

project(Rodin VERSION 0.0.1 LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(RodinUtilities)
include(RodinInstall)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Compute version components
set(RODIN_VERSION         "${PROJECT_VERSION}")
set(RODIN_VERSION_MAJOR   "${PROJECT_VERSION_MAJOR}")
set(RODIN_VERSION_MINOR   "${PROJECT_VERSION_MINOR}")
set(RODIN_VERSION_PATCH   "${PROJECT_VERSION_PATCH}")

# ---- Define install paths --------------------------------------------------
set(RODIN_INSTALL_BINARY_PATH   ${CMAKE_INSTALL_BINDIR}   CACHE PATH
    "Installation directory for executables")
set(RODIN_INSTALL_LIBRARY_PATH  ${CMAKE_INSTALL_LIBDIR}   CACHE PATH
    "Installation directory for libraries")
set(RODIN_INSTALL_INCLUDE_PATH  ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH
    "Installation directory for headers")
set(RODIN_INSTALL_SHARE_PATH    ${CMAKE_INSTALL_DATADIR}  CACHE PATH
    "Installation directory for read-only architecture-independent data")

# Compute version code: major * 10000 + minor * 100 + patch
math(EXPR RODIN_VERSION_CODE "${PROJECT_VERSION_MAJOR} * 10000 + ${PROJECT_VERSION_MINOR} * 100 + ${PROJECT_VERSION_PATCH}")
# Pad to 6 digits manually
if(RODIN_VERSION_CODE LESS 10)
  set(RODIN_VERSION_CODE "00000${RODIN_VERSION_CODE}")
elseif(RODIN_VERSION_CODE LESS 100)
  set(RODIN_VERSION_CODE "0000${RODIN_VERSION_CODE}")
elseif(RODIN_VERSION_CODE LESS 1000)
  set(RODIN_VERSION_CODE "000${RODIN_VERSION_CODE}")
elseif(RODIN_VERSION_CODE LESS 10000)
  set(RODIN_VERSION_CODE "00${RODIN_VERSION_CODE}")
elseif(RODIN_VERSION_CODE LESS 100000)
  set(RODIN_VERSION_CODE "0${RODIN_VERSION_CODE}")
endif()


# ---- Options ---------------------------------------------------------------
option(RODIN_BUILD_SRC                "Build the Rodin source code"               ON)
option(RODIN_BUILD_EXAMPLES           "Build the Rodin examples"                  ON)
option(RODIN_BUILD_UNIT_TESTS         "Build Rodin unit tests"                    OFF)
option(RODIN_BUILD_MANUFACTURED_TESTS "Build Rodin manufactured solution tests"   OFF)
option(RODIN_BUILD_BENCHMARKS         "Build Rodin benchmarks"                    OFF)
option(RODIN_BUILD_DOC                "Build the Rodin documentation"             OFF)
option(RODIN_WITH_PLOT                "Build the Rodin::Plot module"              OFF)
option(RODIN_WITH_PY                  "Build Python bindings"                     OFF)
option(RODIN_USE_MCSS                 "Build documentation with m.css"            OFF)
option(RODIN_USE_MPI                  "Use Rodin with MPI support"                OFF)
option(RODIN_USE_OPENMP               "Use Rodin with OpenMP support"             OFF)
option(RODIN_USE_UMFPACK              "Use Rodin with UMFPACK support"            OFF)
option(RODIN_USE_PASTIX               "Use Rodin with PaStiX support"             OFF)
option(RODIN_USE_KLU                  "Use Rodin with KLU support"                OFF)
option(RODIN_USE_SUPERLU              "Use Rodin with SuperLU support"            OFF)
option(RODIN_USE_SPQR                 "Use Rodin with SPQR support"               OFF)
option(RODIN_USE_PARDISO              "Use Rodin with Pardiso support"            OFF)
option(RODIN_USE_APPLE_ACCELERATE     "Use Rodin with Apple Accelerate support"   OFF)
option(RODIN_USE_SCOTCH               "Use Rodin with Scotch support"             OFF)
option(RODIN_USE_PETSC                "Use Rodin with PETSc support"              OFF)
option(RODIN_SILENCE_WARNINGS         "Silence warnings outputted by Rodin"       OFF)
option(RODIN_SILENCE_EXCEPTIONS       "Silence exceptions thrown by Rodin"        ON)
option(RODIN_CODE_COVERAGE            "Compile with code coverage flags"          OFF)
option(RODIN_LTO                      "Compile with link time optimization"       OFF)
option(RODIN_FEATURE_SUMMARY          "Print Rodin feature summary"               ON)
option(RODIN_THREAD_SAFE              "Enforces thread-safety."                   OFF)
option(RODIN_MULTITHREADED            "Utilizes multithreading capabalities."     ON)
option(RODIN_USE_ASAN                 "Enable AddressSanitizer"                   OFF)
option(RODIN_USE_UBSAN                "Enable UndefinedBehaviorSanitizer"         OFF)
option(RODIN_USE_TSAN                 "Enable ThreadSanitizer"                    OFF)
option(RODIN_USE_MSAN                 "Enable MemorySanitizer (Linux)"            OFF)
option(BUILD_SHARED_LIBS              "Build using shared libraries"              OFF)

set(RODIN_RESOURCES_DIR "${PROJECT_SOURCE_DIR}/resources/")
set(RODIN_THIRD_PARTY_DIR "${PROJECT_SOURCE_DIR}/third-party/")

# ---- Configure compiler ----------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "")

# ---- Set compile options for all the submodules ----------------------------

# Incompatibilities
if(RODIN_WITH_TSAN AND (RODIN_WITH_ASAN OR RODIN_WITH_UBSAN OR RODIN_WITH_MSAN))
  message(FATAL_ERROR "TSan must be used alone.")
endif()
if(RODIN_WITH_MSAN AND (RODIN_WITH_ASAN OR RODIN_WITH_UBSAN OR RODIN_WITH_TSAN))
  message(FATAL_ERROR "MSan must be used alone.")
endif()

if (RODIN_MULTITHREADED)
  set(RODIN_THREAD_SAFE ON CACHE BOOL "" FORCE)
  set(RODIN_USE_OPENMP ON CACHE BOOL "" FORCE)
  set(Boost_USE_MULTITHREADED ON)
else()
  set(RODIN_THREAD_SAFE OFF CACHE BOOL "" FORCE)
  set(RODIN_USE_OPENMP OFF CACHE BOOL "" FORCE)
  # set(Boost_USE_MULTITHREADED OFF)
endif()

# Compose flags once, then apply to all targets defined after this point
set(_SAN_COMPILE "")
set(_SAN_LINK    "")

if(RODIN_USE_ASAN)
  list(APPEND _SAN_COMPILE -fsanitize=address -fno-omit-frame-pointer
                           -fsanitize-address-use-after-scope)
  list(APPEND _SAN_LINK    -fsanitize=address)
endif()

if(RODIN_USE_UBSAN)
  list(APPEND _SAN_COMPILE -fsanitize=undefined)
  list(APPEND _SAN_LINK    -fsanitize=undefined)
  # Optional hard fail:
  # list(APPEND _SAN_COMPILE -fno-sanitize-recover=undefined)
endif()

if(RODIN_USE_TSAN)
  list(APPEND _SAN_COMPILE -fsanitize=thread)
  list(APPEND _SAN_LINK    -fsanitize=thread)
endif()

if(RODIN_USE_MSAN)
  list(APPEND _SAN_COMPILE -fsanitize=memory -fsanitize-memory-track-origins=2
                           -fno-omit-frame-pointer)
  list(APPEND _SAN_LINK    -fsanitize=memory)
endif()

# Apply globally (affects all following targets)
if(_SAN_COMPILE)
  list(APPEND _SAN_COMPILE -fno-inline -fno-inline-functions
                           -fno-omit-frame-pointer -fno-optimize-sibling-calls)
  message(STATUS "Enabling sanitizers with compile options: ${_SAN_COMPILE}")
  add_compile_options(${_SAN_COMPILE})
endif()
if(_SAN_LINK)
  list(APPEND _SAN_LINK -fno-inline -fno-inline-functions
                        -fno-omit-frame-pointer -fno-optimize-sibling-calls)
  add_link_options(${_SAN_LINK})
endif()

rodin_add_compile_options(LANG CXX OPTIONS -fvisibility=hidden)
rodin_add_compile_options(LANG C CXX OPTIONS -fcolor-diagnostics)
rodin_add_compile_options(LANG C CXX OPTIONS -Wno-error=double-promotion -Wno-error=maybe-uninitialized)
rodin_add_compile_options(LANG C CXX OPTIONS -Wno-old-style-cast -Wno-deprecated-anon-enum-enum-conversion)

add_compile_definitions(_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION)
add_compile_definitions(__MAC_OS_X_VERSION_MIN_REQUIRED=__MAC_10_8)

if (RODIN_LTO)
  rodin_add_compile_options(LANG CXX C OPTIONS -flto)
endif()

if (RODIN_CODE_COVERAGE)
  add_link_options(--coverage)
  add_compile_options(--coverage)
  rodin_add_compile_options(LANG CXX C OPTIONS
    -fprofile-update=atomic
    -fprofile-arcs
    -ftest-coverage)
  find_program(LCOV lcov REQUIRED)
  if(DEFINED ENV{GCOV} AND NOT "$ENV{GCOV}" STREQUAL "")
    # Try finding the tool using the environment variable name
    find_program(GCOV "$ENV{GCOV}" REQUIRED)
  else()
    find_program(GCOV gcov REQUIRED)
  endif()
  string(REPLACE " " ";" RODIN_LCOV_FLAGS_LIST "${RODIN_LCOV_FLAGS}")
  add_custom_target(RodinCoverage
    COMMAND ${LCOV} ${RODIN_LCOV_FLAGS_LIST} --gcov-tool ${GCOV} --capture --directory src --directory tests --output-file RodinCoverage.base.info
    COMMAND ${LCOV} ${RODIN_LCOV_FLAGS_LIST} --gcov-tool ${GCOV} --remove RodinCoverage.base.info '/usr/*' '/opt/*' '*/tests/*' '*/third-party/*' '/Library/*' '/Applications/*' --output-file RodinCoverage.info
    COMMAND rm -f RodinCoverage.base.info
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

if (RODIN_BUILD_SRC)
  # ---- Find and build dependencies -----------------------------------------
  # ---- Boost ----
  add_compile_definitions(PRIVATE BOOST_ASIO_HAS_STD_INVOKE_RESULT=1)

  find_package(Boost 1.74 REQUIRED
    COMPONENTS
    filesystem
    serialization
    headers
    )
  include_directories(${Boost_INCLUDE_DIRS})

  # ---- Corrade ----
  message(STATUS "Configuring Corrade...")
  add_subdirectory(third-party/corrade EXCLUDE_FROM_ALL)
  find_package(Corrade REQUIRED Utility)
  set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

  # ---- Eigen ----
  message(STATUS "Configuring Eigen3...")
  find_package(Eigen3 CONFIG REQUIRED NO_MODULE)

  if (RODIN_WITH_PLOT)
    # ---- Magnum ----
    message(STATUS "Configuring Magnum...")
    set(WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)
    add_subdirectory(third-party/magnum EXCLUDE_FROM_ALL)
    find_package(Magnum REQUIRED GL Shaders Sdl2Application)
    set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

    # ---- Magnum-Eigen integration ----
    message(STATUS "Configuring Magnum-Eigen integration...")
    set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third-party/eigen)
    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/third-party/eigen/cmake)
    set(MAGNUM_WITH_EIGEN ON CACHE BOOL "" FORCE)
    add_subdirectory(third-party/magnum-integration EXCLUDE_FROM_ALL)
    find_package(MagnumIntegration REQUIRED Eigen)
  endif()

  add_subdirectory(src)

  if (RODIN_BUILD_EXAMPLES)
    add_subdirectory(examples)
  endif()

  if (RODIN_BUILD_UNIT_TESTS OR RODIN_BUILD_BENCHMARKS OR
      RODIN_BUILD_MANUFACTURED_TESTS)
    # ---- GoogleTest ----
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    add_subdirectory(third-party/googletest EXCLUDE_FROM_ALL)
    if (RODIN_BUILD_BENCHMARKS)
      if (RODIN_LTO)
        set(BENCHMARK_ENABLE_LTO ON CACHE BOOL "" FORCE)
      endif()
      set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
      set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "" FORCE)
      add_subdirectory(third-party/benchmark EXCLUDE_FROM_ALL)
    endif()

    add_subdirectory(tests)
  endif()

  if (RODIN_WITH_PY)
    set(Python3_FIND_REGISTRY LAST)
    set(Python3_FIND_VIRTUALENV FIRST)
    find_package(Python 3.6 COMPONENTS Interpreter Development REQUIRED)
    execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -c
              "import pybind11; print(pybind11.get_cmake_dir())"
      OUTPUT_VARIABLE _tmp_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
    list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
    add_subdirectory(third-party/pybind11)
    add_subdirectory(py)
  endif()
endif()

if (RODIN_BUILD_DOC)
  add_subdirectory(doc)
endif()

if (RODIN_FEATURE_SUMMARY)
  message(STATUS "Rodin feature summary:")
  include(FeatureSummary)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif()

if (RODIN_BUILD_SRC)
  install(
    EXPORT     RodinTargets
    NAMESPACE  Rodin::
    FILE       RodinTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Rodin
  )
endif()

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/RodinConfigVersion.cmake"
  VERSION       ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/RodinConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/RodinConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Rodin
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/RodinConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/RodinConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Rodin
)

install(
  DIRECTORY "${RODIN_RESOURCES_DIR}"
  DESTINATION "${RODIN_INSTALL_SHARE_PATH}/Rodin/resources"
  FILES_MATCHING
    PATTERN "*"
    PATTERN ".*"      EXCLUDE
    PATTERN "*/.*"    EXCLUDE
)

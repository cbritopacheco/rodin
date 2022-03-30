cmake_minimum_required (VERSION 3.12)
cmake_policy(SET CMP0054 NEW)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13)
  cmake_policy(SET CMP0076 NEW)
endif()

# Prohibit in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited.")
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

project(Rodin CXX)

# ---- Options ---------------------------------------------------------------
option(RODIN_BUILD_SRC          "Build the Rodin source code"     ON)
option(RODIN_BUILD_EXAMPLES     "Build the Rodin examples"        ON)
option(RODIN_BUILD_DOC          "Build the Rodin documentation"   OFF)
option(RODIN_USE_MCSS           "Use m.css style documentation"   OFF)
option(RODIN_WITH_PLOT          "Build the Rodin::Plot module"    OFF)
option(RODIN_USE_SUITESPARSE    "Build SuiteSparse integrations"  ON)

# ---- Configuration ---------------------------
include("${PROJECT_SOURCE_DIR}/cmake/configuration.cmake")

# ---- Set the module path ---------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# ---- Configure compiler ----------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT ${CMAKE_CXX_FLAGS_RELEASE})
  add_compile_options(
    -O3
    -w
    -fvisibility=hidden
    -fvisibility-inlines-hidden)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    add_compile_options(-Wno-error=narrowing)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  endif()
endif()

if (NOT ${CMAKE_CXX_FLAGS_DEBUG})
  add_compile_options(
    -O0
    -g
    -ggdb
    -Wno-old-style-cast
    -Wall
    -Wextra
    -pedantic
    -fvisibility=hidden
    -fvisibility-inlines-hidden
    )
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  endif()
endif()

if (RODIN_BUILD_SRC)
  # ---- Find and build dependencies -----------------------------------------
  # ---- MFEM ----
  if (RODIN_USE_SUITESPARSE)
    set(MFEM_USE_SUITESPARSE ON CACHE BOOL TRUE FORCE)
  endif()

  set(BUILD_DOC OFF CACHE BOOL "" FORCE)
  message(STATUS "Configuring mfem ...")
  add_subdirectory(third-party/mfem EXCLUDE_FROM_ALL)
  find_package(MFEM
    REQUIRED PATHS "${PROJECT_BINARY_DIR}/third-party/mfem"
    NO_DEFAULT_PATH)
  message(INFO ${MFEM_INCLUDE_DIRS})

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "")

  # ---- Corrade ----
  message(STATUS "Configuring Corrade...")
  add_subdirectory(third-party/corrade EXCLUDE_FROM_ALL)
  find_package(Corrade REQUIRED Utility)
  set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

  if (RODIN_WITH_PLOT)
    # ---- Magnum ----
    message(STATUS "Configuring Magnum...")
    set(WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)
    add_subdirectory(third-party/magnum EXCLUDE_FROM_ALL)
    find_package(Magnum REQUIRED GL Shaders Sdl2Application)
    set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

    # ---- Eigen ----
    # message(STATUS "Configuring Eigen...")
    # set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    # set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    # add_subdirectory(third-party/eigen)
    # unset(BUILD_TESTING)

    # find_package (Eigen3 3.3 REQUIRED NO_MODULE
    #   HINTS
    #   "${PROJECT_SOURCE_DIR}/third-party/eigen")

    # ---- Magnum-Eigen integration ----
    message(STATUS "Configuring Magnum-Eigen integration...")
    set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third-party/eigen)
    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/third-party/eigen/cmake)
    set(WITH_EIGEN ON CACHE BOOL "" FORCE)
    add_subdirectory(third-party/magnum-integration EXCLUDE_FROM_ALL)
    find_package(MagnumIntegration REQUIRED Eigen)

    # ---- Boost ----
    find_package(Boost 1.65 COMPONENTS System REQUIRED)
  endif()

  add_subdirectory(src)

  if (RODIN_BUILD_EXAMPLES)
    add_subdirectory(examples)
  endif()
endif()

if (RODIN_BUILD_DOC)
  add_subdirectory(doc)
endif()


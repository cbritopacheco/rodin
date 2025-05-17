# @file FindPETSc.cmake
# @brief CMake module to locate the PETSc library and headers.
#
# This module finds the PETSc libraries and headers and defines an imported
# interface target `PETSc::PETSc` for linking against PETSc.
#
# Variables defined by this module:
#
#   PETSc_FOUND             - True if PETSc was found
#   PETSc_INCLUDE_DIRS      - Include directories for PETSc
#   PETSc_LIBRARIES         - Libraries needed to link PETSc
#
# Imported target:
#
#   PETSc::PETSc            - Interface target for linking with PETSc
#
# Usage:
#
#   find_package(PETSc REQUIRED)
#   target_link_libraries(MyTarget PRIVATE PETSc::PETSc)
#
# User may set PETSc_DIR or environment variable PETSC_DIR to aid discovery.
#

include(FindPackageHandleStandardArgs)

# Allow hints via CMake or environment variable
if (DEFINED PETSc_DIR)
  set(_PETSC_HINTS "${PETSc_DIR}")
elseif (DEFINED ENV{PETSC_DIR})
  set(_PETSC_HINTS "$ENV{PETSC_DIR}")
endif()

# Try to find headers and libraries
find_path(PETSc_INCLUDE_DIR
  NAMES petsc.h
  HINTS ${_PETSC_HINTS}
  PATH_SUFFIXES include include/petsc
  DOC "Path to PETSc headers"
)

find_library(PETSc_LIBRARY
  NAMES petsc
  HINTS ${_PETSC_HINTS}
  PATH_SUFFIXES lib
  DOC "Path to PETSc library"
)

# Validate
find_package_handle_standard_args(PETSc
  REQUIRED_VARS PETSc_INCLUDE_DIR PETSc_LIBRARY
)

# Set output variables
if (PETSc_FOUND)
  set(PETSc_INCLUDE_DIRS "${PETSc_INCLUDE_DIR}")
  set(PETSc_LIBRARIES "${PETSc_LIBRARY}")

  if (NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc IMPORTED INTERFACE)
    set_target_properties(PETSc::PETSc PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${PETSc_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${PETSc_LIBRARIES}"
    )
  endif()
endif()

# Hide internals in GUIs
mark_as_advanced(PETSc_INCLUDE_DIR PETSc_LIBRARY)

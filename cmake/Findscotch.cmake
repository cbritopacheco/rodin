# @file FindScotch.cmake
# @brief CMake module to locate the Scotch partitioning library.
#
# This module finds the Scotch libraries and headers and defines an imported
# interface target `Scotch::Scotch` for linking against Scotch.
#
# Variables defined by this module:
#
#   SCOTCH_FOUND           - True if Scotch was found
#   SCOTCH_INCLUDE_DIRS    - Include directories for Scotch
#   SCOTCH_LIBRARIES       - Libraries required to link Scotch
#
# Imported target:
#
#   Scotch::Scotch         - Interface target to link with Scotch
#
# Usage:
#
#   find_package(Scotch REQUIRED)
#   target_link_libraries(MyTarget PRIVATE Scotch::Scotch)
#
# The user may set SCOTCH_DIR or the SCOTCH_DIR environment variable
# to help locate the installation.
#

include(FindPackageHandleStandardArgs)

# Allow user-supplied hints
if (DEFINED SCOTCH_DIR)
  set(_SCOTCH_HINTS "${SCOTCH_DIR}")
elseif (DEFINED ENV{SCOTCH_DIR})
  set(_SCOTCH_HINTS "$ENV{SCOTCH_DIR}")
endif()

# Locate header and libraries
find_path(SCOTCH_INCLUDE_DIR
  NAMES scotch.h
  HINTS ${_SCOTCH_HINTS}
  PATH_SUFFIXES include include/scotch
  DOC "Path to Scotch header files"
)

find_library(SCOTCH_LIBRARY
  NAMES scotch
  HINTS ${_SCOTCH_HINTS}
  PATH_SUFFIXES lib
  DOC "Path to Scotch core library"
)

find_library(SCOTCHERR_LIBRARY
  NAMES scotcherr
  HINTS ${_SCOTCH_HINTS}
  PATH_SUFFIXES lib
  DOC "Path to Scotch error-handling library"
)

# Final values
if (SCOTCH_INCLUDE_DIR)
  set(SCOTCH_INCLUDE_DIRS "${SCOTCH_INCLUDE_DIR}")
endif()

if (SCOTCH_LIBRARY AND SCOTCHERR_LIBRARY)
  set(SCOTCH_LIBRARIES "${SCOTCH_LIBRARY}" "${SCOTCHERR_LIBRARY}")
endif()

# Detection logic
find_package_handle_standard_args(Scotch
  REQUIRED_VARS SCOTCH_INCLUDE_DIRS SCOTCH_LIBRARIES
)

# Define the target
if (SCOTCH_FOUND AND NOT TARGET Scotch::Scotch)
  add_library(Scotch::Scotch IMPORTED INTERFACE)
  set_target_properties(Scotch::Scotch PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SCOTCH_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${SCOTCH_LIBRARIES};m"
  )
endif()

# Hide internals from GUIs
mark_as_advanced(SCOTCH_INCLUDE_DIR SCOTCH_LIBRARY SCOTCHERR_LIBRARY)

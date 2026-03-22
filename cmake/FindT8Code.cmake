# @file FindT8Code.cmake
# @brief CMake module to locate the t8code library and headers.
#
# This module finds the t8code libraries and headers and defines an imported
# interface target `T8Code::T8Code` for linking against t8code.
#
# t8code is a C/C++ library to manage parallel adaptive meshes with various
# element types, using a forest-of-octrees approach for nonconforming AMR.
# See: https://github.com/DLR-AMR/t8code
#
# Variables defined by this module:
#
#   T8Code_FOUND             - True if t8code was found
#   T8Code_INCLUDE_DIRS      - Include directories for t8code
#   T8Code_LIBRARIES         - Libraries needed to link t8code
#
# Imported target:
#
#   T8Code::T8Code           - Interface target for linking with t8code
#
# Usage:
#
#   find_package(T8Code REQUIRED)
#   target_link_libraries(MyTarget PRIVATE T8Code::T8Code)
#
# User may set T8Code_DIR or environment variable T8CODE_DIR to aid discovery.
#

include(FindPackageHandleStandardArgs)

# Allow hints via CMake or environment variable
if (DEFINED T8Code_DIR)
  set(_T8CODE_HINTS "${T8Code_DIR}")
elseif (DEFINED ENV{T8CODE_DIR})
  set(_T8CODE_HINTS "$ENV{T8CODE_DIR}")
endif()

# Try to find headers and libraries
find_path(T8Code_INCLUDE_DIR
  NAMES t8.h
  HINTS ${_T8CODE_HINTS}
  PATH_SUFFIXES include include/t8code
  DOC "Path to t8code headers"
)

find_library(T8Code_LIBRARY
  NAMES t8
  HINTS ${_T8CODE_HINTS}
  PATH_SUFFIXES lib lib64
  DOC "Path to t8code library"
)

# Validate
find_package_handle_standard_args(T8Code
  REQUIRED_VARS T8Code_INCLUDE_DIR T8Code_LIBRARY
)

# Set output variables
if (T8Code_FOUND)
  set(T8Code_INCLUDE_DIRS "${T8Code_INCLUDE_DIR}")
  set(T8Code_LIBRARIES "${T8Code_LIBRARY}")

  if (NOT TARGET T8Code::T8Code)
    add_library(T8Code::T8Code UNKNOWN IMPORTED)
    set_target_properties(T8Code::T8Code PROPERTIES
      IMPORTED_LOCATION "${T8Code_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${T8Code_INCLUDE_DIRS}"
    )
  endif()
endif()

# Hide internals in GUIs
mark_as_advanced(T8Code_INCLUDE_DIR T8Code_LIBRARY)

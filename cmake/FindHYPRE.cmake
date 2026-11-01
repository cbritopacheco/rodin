#[[
         Copyright Carlos BRITO PACHECO 2021 - 2026.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

#[=======================================================================[.rst:
FindHYPRE
=========

Module for locating the hypre library of preconditioners and solvers.

This module defines the following variables:

``HYPRE_FOUND``
   ``TRUE`` iff hypre and MPI have been found.

``HYPRE_VERSION``
   Extracted from ``HYPRE_config.h``.

``HYPRE_INCLUDE_DIRS``
   Include directories holding ``HYPRE.h``.

``HYPRE_LIBRARIES``
   Libraries required to link against hypre.

The following variables may be set to guide the search:

``HYPRE_ROOT`` / ``HYPRE_DIR``
  Installation prefix searched before the default locations. The environment
  variables of the same names are honoured as well.

This module defines the imported target:

``HYPRE::HYPRE``
  The hypre library, with MPI. A hypre built for a GPU carries its CUDA or HIP
  libraries as dependencies of the shared library.
]=======================================================================]

find_path(HYPRE_INCLUDE_DIR
  NAMES HYPRE.h
  HINTS
    ${HYPRE_ROOT}
    ${HYPRE_DIR}
    ENV HYPRE_ROOT
    ENV HYPRE_DIR
  PATH_SUFFIXES include include/hypre
)

find_library(HYPRE_LIBRARY
  NAMES HYPRE
  HINTS
    ${HYPRE_ROOT}
    ${HYPRE_DIR}
    ENV HYPRE_ROOT
    ENV HYPRE_DIR
  PATH_SUFFIXES lib lib64
)

if (HYPRE_INCLUDE_DIR AND EXISTS "${HYPRE_INCLUDE_DIR}/HYPRE_config.h")
  file(STRINGS "${HYPRE_INCLUDE_DIR}/HYPRE_config.h" _HYPRE_VERSION_LINE
    REGEX "^#define[ \t]+HYPRE_RELEASE_VERSION[ \t]+\"[^\"]+\"")
  if (_HYPRE_VERSION_LINE)
    string(REGEX REPLACE "^#define[ \t]+HYPRE_RELEASE_VERSION[ \t]+\"([^\"]+)\".*" "\\1"
      HYPRE_VERSION "${_HYPRE_VERSION_LINE}")
  endif()
  unset(_HYPRE_VERSION_LINE)
endif()

# hypre includes mpi.h in its headers, so MPI is required by the consumer.
find_package(MPI QUIET COMPONENTS C)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE
  REQUIRED_VARS HYPRE_LIBRARY HYPRE_INCLUDE_DIR MPI_C_FOUND
  VERSION_VAR HYPRE_VERSION
)

if (HYPRE_FOUND)
  set(HYPRE_INCLUDE_DIRS "${HYPRE_INCLUDE_DIR}")
  set(HYPRE_LIBRARIES "${HYPRE_LIBRARY}")

  if (NOT TARGET HYPRE::HYPRE)
    add_library(HYPRE::HYPRE UNKNOWN IMPORTED)
    set_target_properties(HYPRE::HYPRE PROPERTIES
      IMPORTED_LOCATION "${HYPRE_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES MPI::MPI_C
    )
  endif()
endif()

mark_as_advanced(
  HYPRE_INCLUDE_DIR
  HYPRE_LIBRARY
)

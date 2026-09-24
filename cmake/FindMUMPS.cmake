#[[
         Copyright Carlos BRITO PACHECO 2021 - 2026.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

#[=======================================================================[.rst:
FindMUMPS
=========

Module for locating the MUMPS multifrontal sparse direct solver.

This module defines the following variables:

``MUMPS_FOUND``
   ``TRUE`` iff MUMPS and its dependencies have been found.

``MUMPS_VERSION``
   Extracted from ``dmumps_c.h``.

``MUMPS_INCLUDE_DIRS``
   Include directories holding ``dmumps_c.h``.

``MUMPS_LIBRARIES``
   Libraries required to link against MUMPS.

The following variables may be set to guide the search:

``MUMPS_ROOT`` / ``MUMPS_DIR``
  Installation prefix searched before the default locations. The environment
  variables of the same names are honoured as well.

This module defines the imported target:

``MUMPS::MUMPS``
  The double-precision MUMPS solver and its dependencies.

MUMPS is built either against a real MPI implementation or against its
sequential ``libmpiseq`` stub. When no stub is present, MPI is required and
linked into the imported target.
]=======================================================================]

find_path(MUMPS_INCLUDE_DIR
  NAMES dmumps_c.h
  HINTS
    ${MUMPS_ROOT}
    ${MUMPS_DIR}
    ENV MUMPS_ROOT
    ENV MUMPS_DIR
  PATH_SUFFIXES include include/mumps include/MUMPS mumps
)

find_library(MUMPS_DMUMPS_LIBRARY
  NAMES dmumps dmumps_seq
  HINTS
    ${MUMPS_ROOT}
    ${MUMPS_DIR}
    ENV MUMPS_ROOT
    ENV MUMPS_DIR
  PATH_SUFFIXES lib lib64
)

find_library(MUMPS_COMMON_LIBRARY
  NAMES mumps_common mumps_common_seq
  HINTS
    ${MUMPS_ROOT}
    ${MUMPS_DIR}
    ENV MUMPS_ROOT
    ENV MUMPS_DIR
  PATH_SUFFIXES lib lib64
)

# Optional: the ordering library shipped with MUMPS, and the sequential MPI
# stub used by builds made without a real MPI implementation.
find_library(MUMPS_PORD_LIBRARY
  NAMES pord pord_seq
  HINTS
    ${MUMPS_ROOT}
    ${MUMPS_DIR}
    ENV MUMPS_ROOT
    ENV MUMPS_DIR
  PATH_SUFFIXES lib lib64
)

find_library(MUMPS_MPISEQ_LIBRARY
  NAMES mpiseq mpiseq_seq
  HINTS
    ${MUMPS_ROOT}
    ${MUMPS_DIR}
    ENV MUMPS_ROOT
    ENV MUMPS_DIR
  PATH_SUFFIXES lib lib64
)

if (MUMPS_INCLUDE_DIR AND EXISTS "${MUMPS_INCLUDE_DIR}/dmumps_c.h")
  file(STRINGS "${MUMPS_INCLUDE_DIR}/dmumps_c.h" _MUMPS_VERSION_LINE
    REGEX "^#define[ \t]+MUMPS_VERSION[ \t]+\"[^\"]+\"")
  if (_MUMPS_VERSION_LINE)
    string(REGEX REPLACE "^#define[ \t]+MUMPS_VERSION[ \t]+\"([^\"]+)\".*" "\\1"
      MUMPS_VERSION "${_MUMPS_VERSION_LINE}")
  endif()
  unset(_MUMPS_VERSION_LINE)
endif()

set(_MUMPS_REQUIRED_VARS MUMPS_DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY MUMPS_INCLUDE_DIR)

# A build without the sequential stub calls into a real MPI implementation, so
# MPI must be available to the consumer as well.
if (NOT MUMPS_MPISEQ_LIBRARY)
  find_package(MPI QUIET COMPONENTS C)
  list(APPEND _MUMPS_REQUIRED_VARS MPI_C_FOUND)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
  REQUIRED_VARS ${_MUMPS_REQUIRED_VARS}
  VERSION_VAR MUMPS_VERSION
)

if (MUMPS_FOUND)
  set(MUMPS_INCLUDE_DIRS "${MUMPS_INCLUDE_DIR}")
  set(MUMPS_LIBRARIES "${MUMPS_DMUMPS_LIBRARY}" "${MUMPS_COMMON_LIBRARY}")

  if (MUMPS_PORD_LIBRARY)
    list(APPEND MUMPS_LIBRARIES "${MUMPS_PORD_LIBRARY}")
  endif()

  if (MUMPS_MPISEQ_LIBRARY)
    list(APPEND MUMPS_LIBRARIES "${MUMPS_MPISEQ_LIBRARY}")
  endif()

  if (NOT TARGET MUMPS::MUMPS)
    add_library(MUMPS::MUMPS UNKNOWN IMPORTED)
    set_target_properties(MUMPS::MUMPS PROPERTIES
      IMPORTED_LOCATION "${MUMPS_DMUMPS_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIRS}"
    )

    set(_MUMPS_LINK_LIBRARIES "${MUMPS_COMMON_LIBRARY}")
    if (MUMPS_PORD_LIBRARY)
      list(APPEND _MUMPS_LINK_LIBRARIES "${MUMPS_PORD_LIBRARY}")
    endif()
    if (MUMPS_MPISEQ_LIBRARY)
      list(APPEND _MUMPS_LINK_LIBRARIES "${MUMPS_MPISEQ_LIBRARY}")
    else()
      list(APPEND _MUMPS_LINK_LIBRARIES MPI::MPI_C)
    endif()

    set_target_properties(MUMPS::MUMPS PROPERTIES
      INTERFACE_LINK_LIBRARIES "${_MUMPS_LINK_LIBRARIES}")
    unset(_MUMPS_LINK_LIBRARIES)
  endif()
endif()

unset(_MUMPS_REQUIRED_VARS)

mark_as_advanced(
  MUMPS_INCLUDE_DIR
  MUMPS_DMUMPS_LIBRARY
  MUMPS_COMMON_LIBRARY
  MUMPS_PORD_LIBRARY
  MUMPS_MPISEQ_LIBRARY
)

#[[
         Copyright Carlos BRITO PACHECO 2021 - 2026.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

#[=======================================================================[.rst:
FindMETIS
=========

Finds the METIS graph partitioning library.

Imported targets
----------------

``METIS::METIS``
  METIS library target.

Result variables
----------------

``METIS_FOUND``
  True if METIS was found.

``METIS_INCLUDE_DIR``
  Directory containing ``metis.h``.

``METIS_LIBRARY``
  METIS library path.
#]=======================================================================]

find_path(METIS_INCLUDE_DIR
  NAMES metis.h)

find_library(METIS_LIBRARY
  NAMES metis)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
  REQUIRED_VARS METIS_INCLUDE_DIR METIS_LIBRARY)

if (METIS_FOUND AND NOT TARGET METIS::METIS)
  add_library(METIS::METIS UNKNOWN IMPORTED)
  set_target_properties(METIS::METIS PROPERTIES
    IMPORTED_LOCATION "${METIS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}")
endif()

mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARY)

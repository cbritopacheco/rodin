# Findscotch.cmake
#
# This module locates the Scotch library and its headers, and sets up the imported
# target "scotch". If a target named "scotch" already exists, its INTERFACE properties
# are overridden with the values found here.
#
# Variables set:
#   SCOTCH_FOUND, SCOTCH_INCLUDE_DIRS, SCOTCH_LIBRARIES
#
# Usage:
#   find_package(scotch REQUIRED)
#   target_link_libraries(your_target PRIVATE scotch)

# Locate the Scotch header file (e.g., scotch.h)
find_path(SCOTCH_INCLUDE_DIR
  NAMES scotch.h
  PATH_SUFFIXES scotch
  DOC "Path to Scotch header files"
)

# Locate the main Scotch library
find_library(SCOTCH_LIBRARY
  NAMES scotch
  DOC "Path to Scotch library"
)

# Locate the Scotch error library
find_library(SCOTCHERR_LIBRARY
  NAMES scotcherr
  DOC "Path to Scotch error library"
)

# Check if all required components were found
if(SCOTCH_INCLUDE_DIR AND SCOTCH_LIBRARY AND SCOTCHERR_LIBRARY)
  set(SCOTCH_FOUND TRUE)
  set(SCOTCH_INCLUDE_DIRS ${SCOTCH_INCLUDE_DIR})
  set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY})
else()
  set(SCOTCH_FOUND FALSE)
endif()

# Either override or create the target "scotch"
if(TARGET scotch)
  # Override the properties of the already defined scotch target
  message(STATUS "Overriding properties of existing target 'scotch'")
  target_include_directories(scotch INTERFACE "${SCOTCH_INCLUDE_DIRS}")
  target_link_libraries(scotch INTERFACE "${SCOTCH_LIBRARY}" "${SCOTCHERR_LIBRARY}" m)
else()
  # Create an INTERFACE target named "scotch"
  add_library(scotch INTERFACE)
  set_target_properties(scotch PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SCOTCH_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${SCOTCH_LIBRARY};${SCOTCHERR_LIBRARY};m"
  )
endif()

# Mark variables as advanced so they don't show in GUIs by default
mark_as_advanced(SCOTCH_INCLUDE_DIR SCOTCH_LIBRARY SCOTCHERR_LIBRARY)

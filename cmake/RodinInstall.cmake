# @brief Installs the files specified by the FILES keyword argument to the
# specified destination directory.
#
# @param FILES        List of files to be installed. These should be relative paths to
#                     the source directory.
# @param DESTINATION  Destination directory where the files will be installed.
#                     This should be an absolute path.
#
# Example usage:
# ```
# rodin_install_files(
#     FILES
#         file1.txt
#         file2.txt
#         subdir/file3.txt
#     DESTINATION
#         ${CMAKE_INSTALL_PREFIX}/my_project
# )
# ```
# This will install `file1.txt` and `file2.txt` to `my_project/` and
# `file3.txt` to `my_project/file3.txt`.
function(rodin_install_files)
    set(options)
    set(oneValueArgs DESTINATION)
    set(multiValueArgs FILES)
    cmake_parse_arguments(RODIN_INSTALL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    if(NOT RODIN_INSTALL_DESTINATION)
        message(FATAL_ERROR "rodin_install_files: The DESTINATION keyword is required for rodin_install_files function.")
    endif()
    foreach(RODIN_INSTALL_FILES_FILE ${RODIN_INSTALL_FILES})
        set(RODIN_INSTALL_SOURCE_FILE
          "${CMAKE_CURRENT_SOURCE_DIR}/${RODIN_INSTALL_FILES_FILE}")
        set(RODIN_INSTALL_DESTINATION_FILE
          "${RODIN_INSTALL_DESTINATION}/${RODIN_INSTALL_FILES_FILE}")
        get_filename_component(RODIN_INSTALL_DESTINATION_DIR
          "${RODIN_INSTALL_DESTINATION_FILE}" DIRECTORY)
        install(
          DIRECTORY "${RODIN_INSTALL_DESTINATION_DIR}"
          DESTINATION "${RODIN_INSTALL_DESTINATION}"
          COMPONENT Development
          OPTIONAL)
        install(
          FILES "${RODIN_INSTALL_SOURCE_FILE}"
          DESTINATION "${RODIN_INSTALL_DESTINATION_DIR}"
          COMPONENT Development)
    endforeach()
endfunction()

# @brief Installs public header files and/or directories into the specified include directory.
#
# @param HEADERS      (optional) List of individual header file paths (relative to source dir).  
# @param DIRECTORIES  (optional) List of directories whose “*.h”/“*.hpp” should be installed recursively.  
# @param DESTINATION  The base include directory under ${CMAKE_INSTALL_INCLUDEDIR} where headers will be placed.  
#
# Example usage (install selected headers):
# ```
# rodin_install_headers(
#   HEADERS
#     FiniteElement.h
#     detail/Integrator.hpp
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Variational
# )
# ```
#
# Example usage (install whole directories):
# ```
# rodin_install_headers(
#   DIRECTORIES
#     Variational
#     Geometry
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin
# )
# ```
function(rodin_install_headers)
  set(oneValueArgs DESTINATION)
  set(multiValueArgs HEADERS DIRECTORIES)
  cmake_parse_arguments(RIH "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT RIH_DESTINATION)
    message(FATAL_ERROR "rodin_install_headers: you must specify DESTINATION")
  endif()

  # individual header files
  foreach(_file IN LISTS RIH_HEADERS)
    get_filename_component(_subdir ${_file} PATH)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
      DESTINATION "${RIH_DESTINATION}/${_subdir}"
    )
  endforeach()

  # entire header directories
  foreach(_dir IN LISTS RIH_DIRECTORIES)
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${_dir}"
      DESTINATION "${RIH_DESTINATION}/${_dir}"
      FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
    )
  endforeach()
endfunction()


include(CheckLinkerFlag)
include (CheckCCompilerFlag)
include (CheckCXXCompilerFlag)


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

# @brief Calculate the relative path for a path given a base directory.
#
# This function calculates the relative path for a path given a base directory.
#
# @param[in] PATH       The path for which the relative path should be calculated.
# @param[in] DIRECTORY  The base directory against which the relative path
#                       should be calculated.
#
# @code
# rodin_get_relpath(REL_PATH PATH "/path/to/some/directory" DIRECTORY "/path/to/base/directory")
# message("Relative Path: ${REL_PATH}")
# @endcode
function(rodin_get_relpath RESULT_VARIABLE)
    cmake_parse_arguments(RODIN_GET_RELPATH "" "PATH;DIRECTORY" "" ${ARGN})
    if(NOT RODIN_GET_RELPATH_PATH)
        message(FATAL_ERROR "rodin_get_relpath: Missing required argument PATH")
    endif()
    # Calculate the relative path
    file(RELATIVE_PATH RODIN_GET_RELPATH_RELATIVE_PATH "${RODIN_GET_RELPATH_DIRECTORY}" "${RODIN_GET_RELPATH_PATH}")
    set(${RESULT_VARIABLE} "${RODIN_GET_RELPATH_RELATIVE_PATH}" PARENT_SCOPE)
endfunction()

# @brief Retrieve a list of files in a specified directory, with filtering options.
#
# The `rodin_get_files_in_dir` function retrieves a list of files in a specified
# directory and stores the result in a CMake variable. This function supports
# filtering by file attributes such as hidden files and directories, allowing
# you to customize the output according to your requirements.
#
# @param RESULT_VARIABLE (Output) The name of the CMake variable that will store
#        the list of filtered files found in the specified directory. The result
#        variable will be available in the caller's scope.
#
# @param PATH (One-Value Argument, Required) The path to the directory where you
#        want to retrieve files.
#
# @param EXCLUDE (Multi-Value Argument, Optional) A list of attributes that you
#        want to exclude from the result. You can specify one or more of the
#        following attributes:
#        - HIDDEN: Excludes hidden files (files starting with a dot).
#        - DIRECTORY: Excludes directories from the result.
#
# @note To use this function, include it in your CMakeLists.txt file and call
#       it with the appropriate parameters. The result variable will contain the
#       filtered list of files and can be accessed for further CMake operations.
#
# @code{.cmake}
# # Include the function in your CMakeLists.txt file.
# include(rodin_get_files_in_dir.cmake)
#
# # Call the function to retrieve files in a directory, excluding hidden files
# # and directories.
# rodin_get_files_in_dir(RESULT_VARIABLE PATH "/path/to/directory" EXCLUDE HIDDEN DIRECTORY)
#
# # Access the result variable in your CMake code.
# message("Files in the directory: ${RESULT_VARIABLE}")
# @endcode
#
# @warning This function uses CMake's `file(GLOB ...)` or `file(GLOB_RECURSE ...)`
#          commands to obtain the list of files and directories based on the
#          provided path.
function(rodin_get_files_in_dir RESULT_VARIABLE)
  set(RODIN_GET_FILES_IN_DIR_OPTIONS GLOB_RECURSE)
  set(RODIN_GET_FILES_IN_DIR_ONE_VALUE_ARGS PATH)
  set(RODIN_GET_FILES_IN_DIR_MULTI_VALUE_ARGS EXCLUDE)
  cmake_parse_arguments(RODIN_GET_FILES_IN_DIR
    "${RODIN_GET_FILES_IN_DIR_OPTIONS}"
    "${RODIN_GET_FILES_IN_DIR_ONE_VALUE_ARGS}"
    "${RODIN_GET_FILES_IN_DIR_MULTI_VALUE_ARGS}"
    ${ARGN})
  if(NOT DEFINED RODIN_GET_FILES_IN_DIR_PATH)
    message(FATAL_ERROR "rodin_get_files_in_dir: You must specify the PATH argument for RODIN_GET_FILES_IN_DIR.")
  endif()
  if(RODIN_GET_FILES_IN_DIR_GLOB_RECURSE)
    file(GLOB_RECURSE RODIN_GET_FILES_IN_DIR_FILES_IN_DIRECTORY "${RODIN_GET_FILES_IN_DIR_PATH}/*")
  else()
    file(GLOB RODIN_GET_FILES_IN_DIR_FILES_IN_DIRECTORY "${RODIN_GET_FILES_IN_DIR_PATH}/*")
  endif()
  set(RODIN_GET_FILES_IN_DIR_FILTERED_FILES "")
  set(RODIN_GET_FILES_IN_DIR_EXCLUDE_OPTIONS ${RODIN_GET_FILES_IN_DIR_EXCLUDE})
  foreach(RODIN_GET_FILES_IN_DIR_FILE ${RODIN_GET_FILES_IN_DIR_FILES_IN_DIRECTORY})
    set(RODIN_GET_FILES_IN_DIR_IS_HIDDEN_FILE FALSE)
    if(RODIN_GET_FILES_IN_DIR_FILE MATCHES "/\\.[^/]*$")
      set(RODIN_GET_FILES_IN_DIR_IS_HIDDEN_FILE TRUE)
    endif()
    set(RODIN_GET_FILES_IN_DIR_IS_DIRECTORY_FILE FALSE)
    if(IS_DIRECTORY ${RODIN_GET_FILES_IN_DIR_FILE})
      set(RODIN_GET_FILES_IN_DIR_IS_DIRECTORY_FILE TRUE)
    endif()
    if(NOT RODIN_GET_FILES_IN_DIR_IS_HIDDEN_FILE OR NOT "HIDDEN" IN_LIST RODIN_GET_FILES_IN_DIR_EXCLUDE_OPTIONS)
      if(NOT RODIN_GET_FILES_IN_DIR_IS_DIRECTORY_FILE OR NOT "DIRECTORY" IN_LIST RODIN_GET_FILES_IN_DIR_EXCLUDE_OPTIONS)
        list(APPEND RODIN_GET_FILES_IN_DIR_FILTERED_FILES ${RODIN_GET_FILES_IN_DIR_FILE})
      endif()
    endif()
  endforeach()
  set(${RESULT_VARIABLE} ${RODIN_GET_FILES_IN_DIR_FILTERED_FILES} PARENT_SCOPE)
endfunction()


#! @brief Add options to the compilation of source files.
#
# Adds options to the COMPILE_OPTIONS directory property. These options are
# used when compiling targets from the current directory and below. Each option
# is added only if the compiler supports it.
#
# @param LANG List of languages for which the options will be added
# @param OPTIONS Compiler options
function(rodin_add_compile_options)
  cmake_parse_arguments(RODIN_ADD_COMPILE_OPTIONS "" "" "LANG;OPTIONS" ${ARGN})
  if (NOT RODIN_ADD_COMPILE_OPTIONS_LANG)
    message (FATAL_ERROR "rodin_add_compile_options: You must specify one or more languages.")
  endif()
  foreach (LANG IN LISTS RODIN_ADD_COMPILE_OPTIONS_LANG)
    if (${LANG} MATCHES CXX)
      foreach (FLAG IN LISTS RODIN_ADD_COMPILE_OPTIONS_OPTIONS)
        string(REGEX REPLACE "[-=+]" "" FLAG_NO_SIGNS ${FLAG})
        check_cxx_compiler_flag(${FLAG} ${LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
        if(${LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
          add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${FLAG}>)
        else()
          message(WARNING "rodin_add_compile_options: ${CMAKE_CXX_COMPILER_ID} does not support flag \"${FLAG}\".")
        endif()
      endforeach()
    elseif (${LANG} MATCHES C)
      foreach (FLAG IN LISTS RODIN_ADD_COMPILE_OPTIONS_OPTIONS)
        string(REGEX REPLACE "[-=+]" "" FLAG_NO_SIGNS ${FLAG})
        check_c_compiler_flag(${FLAG} ${LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
        if(${LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
          add_compile_options($<$<COMPILE_LANGUAGE:C>:${FLAG}>)
        else()
          message(WARNING "rodin_add_compile_options: ${CMAKE_C_COMPILER_ID} does not support flag \"${FLAG}\".")
        endif()
      endforeach()
    elseif (${LANG} MATCHES "")
      # Do nothing
    else ()
      message(FATAL_ERROR "rodin_add_compile_options: Language ${language} not supported.")
    endif()
  endforeach()
endfunction()

#! @brief Add options to the linkage of source files.
#
# Adds options to the LINK_OPTIONS directory property. These options are used
# when compiling targets from the current directory and below. Each option is
# added only if the linker supports it.
#
# @param LANG List of languages for which the options will be added
# @param OPTIONS Compiler options
function(rodin_add_link_options)
  cmake_parse_arguments(RODIN_ADD_LINK_OPTIONS "" "" "LANG;OPTIONS" ${ARGN})
  if (NOT RODIN_ADD_LINK_OPTIONS_LANG)
    message (FATAL_ERROR "rodin_add_link_options: You must specify one or more languages.")
  endif()
  foreach (LANG IN LISTS RODIN_ADD_LINK_OPTIONS_LANG)
    if (${LANG} MATCHES CXX)
      foreach (FLAG IN LISTS RODIN_ADD_LINK_OPTIONS_OPTIONS)
        string(REGEX REPLACE "[-=+]" "" FLAG_NO_SIGNS ${FLAG})
        check_linker_flag(CXX ${FLAG} ${LANG}_LINKER_SUPPORTS_${FLAG_NO_SIGNS})
        if(${LANG}_LINKER_SUPPORTS_${FLAG_NO_SIGNS})
          add_link_options(${FLAG})
        else()
          message(WARNING "rodin_add_link_options: ${CMAKE_CXX_COMPILER_ID} does not support link flag \"${FLAG}\".")
        endif()
      endforeach()
    elseif (${LANG} MATCHES C)
      foreach (FLAG IN LISTS RODIN_ADD_LINK_OPTIONS_OPTIONS)
        string(REGEX REPLACE "[-=+]" "" FLAG_NO_SIGNS ${FLAG})
        check_linker_flag(C ${FLAG} ${LANG}_LINKER_SUPPORTS_${FLAG_NO_SIGNS})
        if(${LANG}_LINKER_SUPPORTS_${FLAG_NO_SIGNS})
          add_link_options(${FLAG})
        else()
          message(WARNING "rodin_add_link_options: ${CMAKE_C_COMPILER_ID} does not support link flag \"${FLAG}\".")
        endif()
      endforeach()
    elseif (${LANG} MATCHES "")
      # Do nothing
    else ()
      message(FATAL_ERROR "rodin_add_link_options: Language ${language} not supported.")
    endif()
  endforeach()
endfunction()


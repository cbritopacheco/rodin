include (CheckCCompilerFlag)
include (CheckCXXCompilerFlag)

#! @brief Add options to the compilation of source files.
#
# Adds options to the COMPILE_OPTIONS directory property. These options are
# used when compiling targets from the current directory and below. Each option
# is added only if the compiler supports it.
#
# @param LANG List of languages for which the options will be added
# @param OPTIONS Compiler options
function(rodin_add_compile_options)
    cmake_parse_arguments(PARSED "" "LANG" "OPTIONS" ${ARGN})
    if (NOT PARSED_LANG)
        message (FATAL_ERROR "You must specify one or more languages.")
    endif()
    foreach (LANG IN LISTS PARSED_LANG)
        if (${LANG} MATCHES CXX)
            foreach (FLAG IN LISTS PARSED_OPTIONS)
                string(REGEX REPLACE "[-=+]" "" FLAG_NO_SIGNS ${FLAG})
                check_cxx_compiler_flag(${FLAG} ${PARSED_LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
                if(${PARSED_LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
                    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${FLAG}>)
                else()
                    message(WARNING "Flag ${FLAG} is not supported in \
                    ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}.")
                endif()
            endforeach()
        elseif(${LANG} MATCHES C)
            foreach (FLAG IN LISTS PARSED_OPTIONS)
                string(REGEX REPLACE "[-=+]" "" flag_no_signs ${FLAG})
                check_c_compiler_flag(${FLAG} ${PARSED_LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
                if(${PARSED_LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
                    add_compile_options($<$<COMPILE_LANGUAGE:C>:${FLAG}>)
                else()
                    message(WARNING "Flag ${FLAG} is not supported in \
                    ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}.")
                endif()
            endforeach()
        else()
            message(FATAL_ERROR "Language ${language} not supported.")
        endif()
    endforeach()
endfunction()

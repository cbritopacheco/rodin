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
    cmake_parse_arguments(PARSED "" "" "LANG;OPTIONS" ${ARGN})
    if (NOT PARSED_LANG)
        message (FATAL_ERROR "You must specify one or more languages.")
    endif()
    foreach (LANG IN LISTS PARSED_LANG)
        if (${LANG} MATCHES CXX)
            foreach (FLAG IN LISTS PARSED_OPTIONS)
                string(REGEX REPLACE "[-=+]" "" FLAG_NO_SIGNS ${FLAG})
                check_cxx_compiler_flag(${FLAG} ${LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
                if(${LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
                    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${FLAG}>)
                endif()
            endforeach()
        elseif (${LANG} MATCHES C)
            foreach (FLAG IN LISTS PARSED_OPTIONS)
                string(REGEX REPLACE "[-=+]" "" FLAG_NO_SIGNS ${FLAG})
                check_c_compiler_flag(${FLAG} ${LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
                if(${LANG}_COMPILER_SUPPORTS_${FLAG_NO_SIGNS})
                    add_compile_options($<$<COMPILE_LANGUAGE:C>:${FLAG}>)
                endif()
            endforeach()
        elseif (${LANG} MATCHES "")
            # Do nothing
        else ()
            message(FATAL_ERROR "Language ${language} not supported.")
        endif()
    endforeach()
endfunction()

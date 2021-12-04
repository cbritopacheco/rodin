# ---- Set colors for displaying messages ------------------------------------
string(ASCII 27 Esc )

set(reset   "${Esc}[m"  )
set(red     "${Esc}[31m")
set(blue    "${Esc}[34m")
set(green   "${Esc}[32m")
set(yellow  "${Esc}[33m")
set(gray    "${Esc}[0;37m")

# ---- Overwrite message function --------------------------------------------
function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif()
endfunction()

set(MESSAGE_QUIET OFF)

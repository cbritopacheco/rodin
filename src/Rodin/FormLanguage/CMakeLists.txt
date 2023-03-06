set(RodinFormLanguage_HEADERS
  Base.h
  List.h
  ForwardDecls.h)

set(RodinFormLanguage_SRCS Base.cpp)

add_library(RodinFormLanguage ${RodinFormLanguage_SRCS} ${RodinFormLanguage_HEADERS})
add_library(Rodin::FormLanguage ALIAS RodinFormLanguage)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinFormLanguage PUBLIC $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(RodinFormLanguage PUBLIC Rodin)




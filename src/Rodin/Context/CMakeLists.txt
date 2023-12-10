set(RodinContext_HEADERS
  Sequential.h)

set(RodinContext_SRCS Cast.cpp)

add_library(RodinContext ${RodinContext_SRCS} ${RodinContext_HEADERS})
add_library(Rodin::Context ALIAS RodinContext)

target_include_directories(RodinContext
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)



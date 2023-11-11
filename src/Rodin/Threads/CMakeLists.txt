set(RodinThreads_HEADERS )

set(RodinThreads_SRCS )

add_library(RodinThreads INTERFACE)
add_library(Rodin::Threads ALIAS RodinThreads)

target_include_directories(RodinThreads
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinThreads INTERFACE Rodin::Alert)



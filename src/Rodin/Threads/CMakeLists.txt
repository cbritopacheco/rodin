set(RodinThreads_HEADERS ThreadPool.h)

set(RodinThreads_SRCS ThreadPool.cpp)

add_library(RodinThreads ${RodinThreads_SRCS} ${RodinThreads_HEADERS})
add_library(Rodin::Threads ALIAS RodinThreads)

target_include_directories(RodinThreads
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinThreads PUBLIC Rodin::Alert)
target_link_libraries(RodinThreads PUBLIC BS_thread_pool)



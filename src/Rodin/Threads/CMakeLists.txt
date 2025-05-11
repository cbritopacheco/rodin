set(RodinThreads_HEADERS ThreadPool.h)

set(RodinThreads_SRCS ThreadPool.cpp)

add_library(RodinThreads ${RodinThreads_SRCS} ${RodinThreads_HEADERS})
add_library(Rodin::Threads ALIAS RodinThreads)

target_include_directories(RodinThreads
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinThreads PUBLIC Rodin::Alert)
target_link_libraries(RodinThreads PUBLIC BS_thread_pool)

install(TARGETS RodinThreads
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Threads
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Threads
)

#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

add_library(RodinThreads INTERFACE)
add_library(Rodin::Threads ALIAS RodinThreads)
set_target_properties(RodinThreads PROPERTIES EXPORT_NAME Threads)

target_include_directories(RodinThreads
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinThreads
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

target_link_libraries(RodinThreads
  INTERFACE
    Rodin::Alert
    Threads::Threads
)

# INTERFACE library: no runtime/library/archive destinations needed
install(
  TARGETS RodinThreads
  EXPORT  RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Threads
)

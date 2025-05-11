set(RodinAssembly_HEADERS
  AssemblyBase.h
  Sequential.h
  Multithreaded.h)

set(RodinAssembly_SRCS
  Sequential.cpp
  Multithreaded.cpp)

add_library(RodinAssembly
  ${RodinAssembly_SRCS} ${RodinAssembly_HEADERS})
add_library(Rodin::Assembly ALIAS RodinAssembly)

if (RODIN_USE_OPENMP)
  target_sources(RodinAssembly PUBLIC OpenMP.h)
  target_sources(RodinAssembly PUBLIC OpenMP.cpp)
  target_link_libraries(RodinAssembly PUBLIC OpenMP::OpenMP_CXX)
endif()

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinAssembly
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinAssembly
  PUBLIC
  Threads::Threads
  Rodin::Utility
  Rodin::FormLanguage
  Rodin::Variational)

install(TARGETS RodinAssembly
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Assembly
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Assembly
)

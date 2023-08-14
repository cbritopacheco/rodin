set(RodinVariational_HEADERS
  AssemblyBase.h
  Native.h)

set(RodinVariational_SRCS
  AssemblyBase.cpp
  Native.cpp)

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
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  PUBLIC
  ${MFEM_INCLUDE_DIRS})

target_link_libraries(RodinAssembly
  PUBLIC
  Rodin::IO
  Rodin::Alert
  Rodin::Utility
  Rodin::Geometry
  Rodin::FormLanguage
  Rodin::Variational)

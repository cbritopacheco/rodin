#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinAssembly_HEADERS
  AssemblyBase.h
  Sequential.h
)

set(RodinAssembly_SRCS
  Sequential.cpp
)

add_library(RodinAssembly
  ${RodinAssembly_SRCS}
  ${RodinAssembly_HEADERS}
)
add_library(Rodin::Assembly ALIAS RodinAssembly)
set_target_properties(RodinAssembly PROPERTIES EXPORT_NAME Assembly)

# ---- Link targets ----------------------------------------------------------
# Reuse Rodin's exported include dirs; expose them PUBLIC
target_include_directories(RodinAssembly
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_link_libraries(RodinAssembly
  PUBLIC
    Rodin::Variational
    Rodin::FormLanguage
)

if (RODIN_USE_OPENMP)
  target_sources(RodinAssembly PRIVATE OpenMP.cpp)
  target_link_libraries(RodinAssembly PUBLIC OpenMP::OpenMP_CXX)
endif()

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinAssembly
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

# Install public headers under include/Rodin/Assembly
rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Assembly
)

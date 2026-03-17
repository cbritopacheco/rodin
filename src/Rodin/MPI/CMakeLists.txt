#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinMPI_HEADERS
  Geometry/Mesh.h
  Geometry/Sharder.h
  Assembly/MPI.h
  IO/HDF5.h
  IO/MFEM.h
  IO/MEDIT.h
)

set(RodinMPI_SRCS
  Geometry/Mesh.cpp
  Geometry/Sharder.cpp
  Context/MPI.cpp
  Assembly/MPI.cpp
)

add_library(RodinMPI
  ${RodinMPI_SRCS}
  ${RodinMPI_HEADERS}
)
add_library(Rodin::MPI ALIAS RodinMPI)
set_target_properties(RodinMPI PROPERTIES EXPORT_NAME MPI)

# ---- Includes / defs -------------------------------------------------------
target_include_directories(RodinMPI
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinMPI
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Dependencies ----------------------------------------------------------
target_link_libraries(RodinMPI
  PUBLIC
    Rodin::Geometry
    Rodin::Alert
    Rodin::IO
    Boost::mpi
    MPI::MPI_CXX
)

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinMPI
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/MPI
)

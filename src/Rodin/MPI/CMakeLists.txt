#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinMPI_HEADERS
  Geometry/Mesh.h
  Geometry/Sharder.h
  Assembly/MPI.h
  )

set(RodinMPI_SRCS
  Geometry/Mesh.cpp
  Geometry/Sharder.cpp
  Context/MPI.cpp
  Assembly/MPI.cpp
  )

add_library(RodinMPI ${RodinMPI_HEADERS} ${RodinMPI_SRCS})
add_library(Rodin::MPI ALIAS RodinMPI)
set_target_properties(RodinMPI PROPERTIES EXPORT_NAME MPI)

target_include_directories(RodinMPI
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

find_package(Boost 1.74 REQUIRED COMPONENTS mpi)

target_link_libraries(RodinMPI
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  Boost::mpi
  MPI::MPI_CXX)

install(TARGETS RodinMPI
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/MPI
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/MPI
)


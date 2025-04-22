#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinMPI_HEADERS
  MPIMesh.h
  MPISharder.h)

set(RodinMPI_SRCS
  MPIMesh.cpp
  MPISharder.cpp
  MPIContext.cpp
  )

add_library(RodinMPI ${RodinMPI_HEADERS} ${RodinMPI_SRCS})
add_library(Rodin::MPI ALIAS RodinMPI)

find_package(Boost 1.74 REQUIRED COMPONENTS mpi)

target_link_libraries(RodinMPI
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  Boost::mpi
  MPI::MPI_CXX)



set(RodinParallel_HEADERS
  MPIParallelizable.h)

set(RodinParallel_SRCS
  MPIParallelizable.cpp)
add_library(RodinParallel ${RodinParallel_SRCS} ${RodinParallel_HEADERS})
add_library(Rodin::Parallel ALIAS RodinParallel)

target_include_directories(RodinParallel
  PRIVATE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

if (RODIN_USE_MPI)
target_link_libraries(RodinParallel
  PUBLIC
  Boost::mpi
  Boost::serialization
  MPI::MPI_CXX)
endif()




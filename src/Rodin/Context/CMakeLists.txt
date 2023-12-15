set(RodinContext_HEADERS Sequential.h)

set(RodinContext_SRCS Sequential.cpp)

add_library(RodinContext ${RodinContext_SRCS} ${RodinContext_HEADERS})
add_library(Rodin::Context ALIAS RodinContext)

target_include_directories(RodinContext
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

if (RODIN_USE_MPI)
  target_sources(RodinContext PUBLIC MPI.h PRIVATE MPI.cpp)
  target_link_libraries(RodinContext PUBLIC Boost::mpi MPI::MPI_CXX)
endif()


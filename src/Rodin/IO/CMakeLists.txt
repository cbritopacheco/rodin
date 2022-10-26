set(RodinIO_HEADERS
  Helpers.h
  Loader.h
  Printer.h
  MeshLoader.h
  MeshPrinter.h
  GridFunctionLoader.h
  # GridFunctionPrinter.h
  )

set(RodinIO_SRCS
  Helpers.cpp
  MeshLoader.cpp
  MeshPrinter.cpp
  GridFunctionLoader.cpp
  # GridFunctionPrinter.cpp
  )
add_library(RodinIO ${RodinIO_SRCS} ${RodinIO_HEADERS})
add_library(Rodin::IO ALIAS RodinIO)

target_include_directories(RodinIO
  PUBLIC
  ${MFEM_INCLUDE_DIRS})

target_link_libraries(RodinIO
  PRIVATE
  mfem
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  Rodin::Variational
  Boost::filesystem)

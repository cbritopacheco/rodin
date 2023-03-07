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

target_link_libraries(RodinIO
  PRIVATE
  mfem
  PUBLIC
  Boost::filesystem
  Rodin::Alert
  Rodin::Geometry
  Rodin::Variational)

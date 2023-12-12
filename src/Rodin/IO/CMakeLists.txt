set(RodinIO_HEADERS
  MFEM.h
  MEDIT.h
  Loader.h
  Printer.h
  MeshLoader.h
  MeshPrinter.h
  GridFunctionLoader.h
  GridFunctionPrinter.h)

set(RodinIO_SRCS
  MFEM.cpp
  MEDIT.cpp
  MeshLoader.cpp
  MeshPrinter.cpp
  GridFunctionLoader.cpp)
add_library(RodinIO ${RodinIO_SRCS} ${RodinIO_HEADERS})
add_library(Rodin::IO ALIAS RodinIO)

target_link_libraries(RodinIO
  PUBLIC
  Rodin::Alert
  Rodin::Geometry
  Rodin::Variational
  Boost::filesystem)
